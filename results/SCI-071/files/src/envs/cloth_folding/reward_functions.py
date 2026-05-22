"""Reward functions for cloth folding tasks."""

from __future__ import annotations

from typing import Callable, Mapping, Sequence

import numpy as np
import torch

ArrayLike = np.ndarray | torch.Tensor | Sequence[Sequence[float]]


def _to_tensor(points: ArrayLike) -> torch.Tensor:
    array = np.asarray(points, dtype=np.float32)
    tensor = torch.as_tensor(array, dtype=torch.float32)
    if tensor.ndim != 2:
        raise ValueError("expected an array of shape (N, D)")
    return tensor


def _projected_area(points: torch.Tensor, grid_size: int = 128) -> torch.Tensor:
    xy = points[:, :2]
    mins = xy.min(dim=0).values
    maxs = xy.max(dim=0).values
    span = torch.clamp(maxs - mins, min=1e-6)
    normalized = (xy - mins) / span
    indices = torch.clamp((normalized * (grid_size - 1)).long(), 0, grid_size - 1)
    occupancy = torch.zeros((grid_size, grid_size), dtype=torch.float32)
    occupancy[indices[:, 1], indices[:, 0]] = 1.0
    return occupancy.mean() * torch.prod(span)


def chamfer_distance_reward(current_particles: ArrayLike, goal_particles: ArrayLike) -> torch.Tensor:
    """Reward closeness to a goal particle configuration via negative Chamfer distance."""

    current = _to_tensor(current_particles)
    goal = _to_tensor(goal_particles)
    distances = torch.cdist(current, goal, p=2)
    chamfer = 0.5 * (distances.min(dim=1).values.mean() + distances.min(dim=0).values.mean())
    return -chamfer


def coverage_reward(current_particles: ArrayLike, goal_area: float | Mapping[str, float]) -> torch.Tensor:
    """Reward matching a target 2D coverage area."""

    current = _to_tensor(current_particles)
    if isinstance(goal_area, Mapping):
        target_area = float(goal_area.get("target_area", goal_area.get("area", 0.0)))
    else:
        target_area = float(goal_area)
    if target_area <= 0.0:
        raise ValueError("goal_area must be positive")
    current_area = _projected_area(current)
    score = torch.clamp(current_area / target_area, min=0.0, max=1.0)
    return score


def fold_line_reward(particles: ArrayLike, fold_line: ArrayLike) -> torch.Tensor:
    """Reward alignment of particles with a target fold line."""

    particle_tensor = _to_tensor(particles)
    line = _to_tensor(fold_line)
    if line.shape[0] < 2:
        raise ValueError("fold_line must provide at least two points")
    point_on_line = line[0]
    direction = line[1] - line[0]
    direction = direction / (torch.linalg.norm(direction) + 1e-8)
    relative = particle_tensor - point_on_line
    parallel = torch.sum(relative * direction, dim=-1, keepdim=True) * direction
    perpendicular = relative - parallel
    mean_distance = torch.linalg.norm(perpendicular, dim=-1).mean()
    return torch.exp(-10.0 * mean_distance)


def smoothness_penalty(action_sequence: ArrayLike) -> torch.Tensor:
    """Penalize high curvature in an action sequence."""

    actions = _to_tensor(action_sequence)
    if actions.shape[0] < 3:
        return torch.tensor(0.0, dtype=torch.float32)
    acceleration = actions[2:] - 2.0 * actions[1:-1] + actions[:-2]
    return torch.linalg.norm(acceleration, dim=-1).mean()


def composite_reward(weights_dict: Mapping[str, float]) -> Callable[..., torch.Tensor]:
    """Create a weighted reward combiner.

    The returned callable accepts either a mapping as its first argument or named reward
    components via keyword arguments.
    """

    def _combine(*args: Mapping[str, float | torch.Tensor], **kwargs: float | torch.Tensor) -> torch.Tensor:
        if args:
            components = dict(args[0])
            components.update(kwargs)
        else:
            components = kwargs
        total = torch.tensor(0.0, dtype=torch.float32)
        for name, weight in weights_dict.items():
            if name in components:
                total = total + float(weight) * torch.as_tensor(components[name], dtype=torch.float32)
        return total

    return _combine


__all__ = [
    "chamfer_distance_reward",
    "composite_reward",
    "coverage_reward",
    "fold_line_reward",
    "smoothness_penalty",
]
