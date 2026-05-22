"""Primitive folding actions and simple motion planning utilities."""

from __future__ import annotations

from typing import Iterable, List, Sequence

import numpy as np


ArrayLike = np.ndarray | Sequence[float]


def interpolate_trajectory(
    waypoints: Sequence[ArrayLike],
    steps_per_segment: int = 10,
) -> np.ndarray:
    """Linearly interpolate a multi-segment Cartesian trajectory."""

    if len(waypoints) < 2:
        raise ValueError("at least two waypoints are required")
    trajectory: List[np.ndarray] = []
    waypoint_arrays = [np.asarray(point, dtype=np.float32) for point in waypoints]
    for start, end in zip(waypoint_arrays[:-1], waypoint_arrays[1:]):
        for alpha in np.linspace(0.0, 1.0, steps_per_segment, endpoint=False, dtype=np.float32):
            trajectory.append((1.0 - alpha) * start + alpha * end)
    trajectory.append(waypoint_arrays[-1])
    return np.stack(trajectory, axis=0)


def _collision_free_point(env: object, point: np.ndarray) -> np.ndarray:
    if hasattr(env, "is_position_collision_free"):
        if bool(env.is_position_collision_free(point)):
            return point
    adjusted = point.copy()
    min_height = float(getattr(env, "minimum_gripper_height", 0.02))
    adjusted[2] = max(adjusted[2], min_height)
    return adjusted


def plan_collision_aware_trajectory(
    env: object,
    pick_pos: ArrayLike,
    place_pos: ArrayLike,
    height: float = 0.15,
    via_points: Sequence[ArrayLike] | None = None,
    steps_per_segment: int = 10,
) -> np.ndarray:
    """Plan a simple collision-aware trajectory using lifted via points."""

    pick = np.asarray(pick_pos, dtype=np.float32)
    place = np.asarray(place_pos, dtype=np.float32)
    lift = np.array([0.0, 0.0, height], dtype=np.float32)
    waypoints = [_collision_free_point(env, pick + lift), pick, pick + lift]
    if via_points is not None:
        waypoints.extend(_collision_free_point(env, np.asarray(point, dtype=np.float32)) for point in via_points)
    midpoint = 0.5 * (pick + place)
    midpoint[2] = max(midpoint[2], height)
    waypoints.extend([_collision_free_point(env, midpoint), _collision_free_point(env, place + lift), place, _collision_free_point(env, place + lift)])
    return interpolate_trajectory(waypoints, steps_per_segment=steps_per_segment)


def pick_and_place(
    env: object,
    pick_pos: ArrayLike,
    place_pos: ArrayLike,
    height: float = 0.15,
) -> np.ndarray:
    """Execute a basic pick-and-place motion."""

    trajectory = plan_collision_aware_trajectory(env, pick_pos, place_pos, height=height)
    if hasattr(env, "execute_trajectory"):
        env.execute_trajectory(trajectory, grasp_at=np.asarray(pick_pos, dtype=np.float32), release_at=np.asarray(place_pos, dtype=np.float32))
    elif hasattr(env, "backend") and hasattr(env.backend, "apply_pick_and_place"):
        env.backend.apply_pick_and_place(np.asarray(pick_pos, dtype=np.float32), np.asarray(place_pos, dtype=np.float32), height=height)
    elif hasattr(env, "apply_pick_and_place"):
        env.apply_pick_and_place(np.asarray(pick_pos, dtype=np.float32), np.asarray(place_pos, dtype=np.float32), height=height)
    else:
        raise AttributeError("environment does not expose a pick-and-place execution interface")
    return trajectory


def fold_edge(env: object, edge_id: int, fold_angle: float = np.pi) -> np.ndarray:
    """Fold along an edge by moving the edge center toward the cloth center."""

    if not hasattr(env, "mesh"):
        raise AttributeError("environment must expose a mesh attribute")
    mesh = env.mesh
    if edge_id < 0 or edge_id >= len(mesh.edges):
        raise IndexError("edge_id out of range")
    edge = mesh.edges[edge_id]
    edge_vertices = mesh.vertices[edge]
    pick = edge_vertices.mean(axis=0)
    cloth_center = mesh.vertices.mean(axis=0)
    fold_fraction = float(np.clip(fold_angle / np.pi, 0.0, 1.0))
    place = pick + fold_fraction * (cloth_center - pick)
    place[2] = max(place[2], float(getattr(env, "minimum_gripper_height", 0.02)))
    return pick_and_place(env, pick, place, height=float(getattr(env, "grasp_height", 0.15)))


def flatten(env: object) -> None:
    """Flatten the cloth to its nominal rest configuration."""

    if hasattr(env, "backend") and hasattr(env.backend, "flatten"):
        env.backend.flatten()
    elif hasattr(env, "reset_to_rest_state"):
        env.reset_to_rest_state()
    else:
        raise AttributeError("environment does not expose a flatten operation")


def shake(env: object, amplitude: float = 0.03, frequency: float = 2.0) -> None:
    """Shake the cloth to help unfold wrinkles."""

    if hasattr(env, "backend") and hasattr(env.backend, "shake"):
        env.backend.shake(amplitude=amplitude, frequency=frequency)
        return
    if hasattr(env, "particles"):
        particles = np.asarray(env.particles, dtype=np.float32)
        timesteps = np.linspace(0.0, 1.0, 8, dtype=np.float32)
        for t in timesteps:
            particles[:, 0] += amplitude * np.sin(2.0 * np.pi * frequency * t)
        env.particles = particles
        return
    raise AttributeError("environment does not expose a shake operation")


__all__ = [
    "flatten",
    "fold_edge",
    "interpolate_trajectory",
    "pick_and_place",
    "plan_collision_aware_trajectory",
    "shake",
]
