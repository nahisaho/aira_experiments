"""Evaluation metrics for cloth folding state comparison."""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np
import torch

ArrayLike = np.ndarray | torch.Tensor | Sequence[Sequence[float]]


def _to_tensor(points: ArrayLike) -> torch.Tensor:
    array = np.asarray(points, dtype=np.float32)
    tensor = torch.as_tensor(array, dtype=torch.float32)
    if tensor.ndim != 2:
        raise ValueError("points must have shape (N, D)")
    return tensor


def _to_numpy(points: ArrayLike) -> np.ndarray:
    array = np.asarray(points, dtype=np.float32)
    if array.ndim != 2:
        raise ValueError("points must have shape (N, D)")
    return array


def chamfer_distance(
    points_a: ArrayLike,
    points_b: ArrayLike,
    bidirectional: bool = True,
    squared: bool = False,
) -> float:
    """Compute Chamfer distance between two point sets."""

    a = _to_tensor(points_a)
    b = _to_tensor(points_b)
    distances = torch.cdist(a, b, p=2)
    if squared:
        distances = distances.square()
    a_to_b = distances.min(dim=1).values.mean()
    if not bidirectional:
        return float(a_to_b)
    b_to_a = distances.min(dim=0).values.mean()
    return float(0.5 * (a_to_b + b_to_a))


def earth_movers_distance(
    points_a: ArrayLike,
    points_b: ArrayLike,
    epsilon: float = 0.01,
    max_iterations: int = 100,
) -> float:
    """Approximate Earth Mover's Distance with Sinkhorn iterations."""

    a = _to_tensor(points_a)
    b = _to_tensor(points_b)
    n, m = a.shape[0], b.shape[0]
    if n == 0 or m == 0:
        raise ValueError("point sets must be non-empty")
    mu = torch.full((n,), 1.0 / n, dtype=torch.float32)
    nu = torch.full((m,), 1.0 / m, dtype=torch.float32)
    cost = torch.cdist(a, b, p=2)
    kernel = torch.exp(-cost / max(epsilon, 1e-6)) + 1e-9
    u = torch.ones_like(mu)
    v = torch.ones_like(nu)
    for _ in range(max_iterations):
        u = mu / (kernel @ v + 1e-9)
        v = nu / (kernel.transpose(0, 1) @ u + 1e-9)
    transport = torch.outer(u, v) * kernel
    return float(torch.sum(transport * cost))


def _project_to_grid(points: np.ndarray, bounds: np.ndarray, grid_size: int) -> np.ndarray:
    mins = bounds[0]
    maxs = bounds[1]
    span = np.maximum(maxs - mins, 1e-6)
    normalized = (points[:, :2] - mins) / span
    indices = np.clip((normalized * (grid_size - 1)).astype(np.int32), 0, grid_size - 1)
    occupancy = np.zeros((grid_size, grid_size), dtype=bool)
    for x_idx, y_idx in indices:
        occupancy[y_idx, x_idx] = True
    return occupancy


def intersection_over_union(
    points_a: ArrayLike,
    points_b: ArrayLike,
    grid_size: int = 128,
    padding: float = 0.05,
) -> float:
    """Compute IoU between two 2D cloth projections."""

    a = _to_numpy(points_a)
    b = _to_numpy(points_b)
    stacked = np.vstack([a[:, :2], b[:, :2]])
    mins = stacked.min(axis=0)
    maxs = stacked.max(axis=0)
    span = np.maximum(maxs - mins, 1e-6)
    bounds = np.stack([mins - padding * span, maxs + padding * span], axis=0)
    occ_a = _project_to_grid(a, bounds, grid_size)
    occ_b = _project_to_grid(b, bounds, grid_size)
    union = np.logical_or(occ_a, occ_b).sum()
    if union == 0:
        return 0.0
    intersection = np.logical_and(occ_a, occ_b).sum()
    return float(intersection / union)


def _projected_area(points: np.ndarray, grid_size: int = 128) -> float:
    mins = points[:, :2].min(axis=0)
    maxs = points[:, :2].max(axis=0)
    bounds = np.stack([mins, maxs], axis=0)
    occupancy = _project_to_grid(points, bounds, grid_size)
    total_area = float(np.prod(np.maximum(maxs - mins, 1e-6)))
    occupied_fraction = occupancy.mean()
    return occupied_fraction * total_area


def normalized_coverage(points: ArrayLike, target_area: float, max_area: Optional[float] = None) -> float:
    """Measure normalized 2D area coverage of the cloth."""

    if target_area <= 0.0:
        raise ValueError("target_area must be positive")
    area = _projected_area(_to_numpy(points))
    if max_area is not None:
        area = min(area, max_area)
    return float(np.clip(area / target_area, 0.0, 1.0))


def success_rate(successes: Sequence[bool | float], threshold: float = 0.5) -> float:
    """Compute success rate from booleans or scalar scores."""

    values = np.asarray(successes)
    if values.size == 0:
        return 0.0
    if values.dtype == np.bool_:
        success_mask = values
    else:
        success_mask = values >= threshold
    return float(np.mean(success_mask.astype(np.float32)))


def _mirror_across_line(points: torch.Tensor, point_on_line: torch.Tensor, direction: torch.Tensor) -> torch.Tensor:
    direction = direction / (torch.linalg.norm(direction) + 1e-8)
    relative = points - point_on_line
    parallel = torch.sum(relative * direction, dim=-1, keepdim=True) * direction
    perpendicular = relative - parallel
    return points - 2.0 * perpendicular


def fold_quality_metric(
    current_particles: ArrayLike,
    goal_particles: Optional[ArrayLike] = None,
    fold_line: Optional[ArrayLike] = None,
) -> Mapping[str, float]:
    """Compute fold quality from symmetry and alignment."""

    current = _to_tensor(current_particles)
    if fold_line is None:
        mins, maxs = current[:, :2].min(dim=0).values, current[:, :2].max(dim=0).values
        center = 0.5 * (mins + maxs)
        point_on_line = torch.tensor([center[0], center[1], 0.0], dtype=torch.float32)
        direction = torch.tensor([0.0, 1.0, 0.0], dtype=torch.float32)
    else:
        line = _to_tensor(fold_line)
        point_on_line = line[0]
        direction = line[1] - line[0]
    mirrored = _mirror_across_line(current, point_on_line, direction)
    symmetry_error = torch.cdist(mirrored, current, p=2).min(dim=1).values.mean()
    symmetry = float(torch.exp(-5.0 * symmetry_error).clamp(0.0, 1.0))
    if goal_particles is not None:
        alignment_error = chamfer_distance(current, goal_particles)
        alignment = float(np.exp(-5.0 * alignment_error))
    else:
        line_direction = direction / (torch.linalg.norm(direction) + 1e-8)
        relative = current - point_on_line
        parallel = torch.sum(relative * line_direction, dim=-1, keepdim=True) * line_direction
        perpendicular = relative - parallel
        alignment = float(torch.exp(-10.0 * torch.linalg.norm(perpendicular, dim=-1).mean()).clamp(0.0, 1.0))
    combined = 0.5 * (symmetry + alignment)
    return {"symmetry": symmetry, "alignment": alignment, "combined": combined}


__all__ = [
    "chamfer_distance",
    "earth_movers_distance",
    "fold_quality_metric",
    "intersection_over_union",
    "normalized_coverage",
    "success_rate",
]
