"""Gym-compatible cloth folding environment with configurable simulation backends."""

from __future__ import annotations

import warnings
from dataclasses import replace
from typing import Any, Dict, Optional, Sequence

import numpy as np
import torch

try:
    import gymnasium as gym
    from gymnasium import spaces
except ImportError:  # pragma: no cover - compatibility fallback
    class _FallbackEnv:
        metadata: Dict[str, Any] = {}

        def reset(self, *, seed: Optional[int] = None, options: Optional[Dict[str, Any]] = None):
            self.np_random = np.random.default_rng(seed)
            return None

    class _FallbackBox:
        def __init__(self, low, high, shape=None, dtype=np.float32):
            self.low = np.array(low if shape is None else np.full(shape, low), dtype=dtype)
            self.high = np.array(high if shape is None else np.full(shape, high), dtype=dtype)
            self.shape = tuple(self.low.shape if shape is None else shape)
            self.dtype = dtype

    class _FallbackDict(dict):
        def __init__(self, spaces_dict):
            super().__init__(spaces_dict)
            self.spaces = spaces_dict

    class _FallbackSpaces:
        Box = _FallbackBox
        Dict = _FallbackDict

    class _FallbackGym:
        Env = _FallbackEnv

    gym = _FallbackGym()
    spaces = _FallbackSpaces()

from src.envs.cloth_folding.cloth_mesh import ClothMesh, generate_rectangular_cloth_mesh
from src.envs.cloth_folding.reward_functions import (
    chamfer_distance_reward,
    composite_reward,
    coverage_reward,
    fold_line_reward,
    smoothness_penalty,
)
from src.utils.config import ClothFoldingEnvConfig, default_cloth_folding_config
from src.utils.metrics import chamfer_distance, fold_quality_metric, intersection_over_union, normalized_coverage


class NumpyClothBackend:
    """Lightweight fallback cloth backend used when SoftGym or Isaac Gym are unavailable."""

    def __init__(self, mesh: ClothMesh, grasp_radius: float, rng: np.random.Generator):
        self.mesh = mesh.copy()
        self.rest_positions = mesh.vertices.copy()
        self.particles = mesh.vertices.copy()
        self.grasp_radius = grasp_radius
        self.rng = rng
        self.adjacency = self._build_adjacency(mesh.edges, mesh.vertices.shape[0])

    @staticmethod
    def _build_adjacency(edges: np.ndarray, num_vertices: int) -> list[list[int]]:
        adjacency = [[] for _ in range(num_vertices)]
        for start, end in edges:
            adjacency[int(start)].append(int(end))
            adjacency[int(end)].append(int(start))
        return adjacency

    def reset(self, randomize_pose: bool, randomize_height: float) -> np.ndarray:
        particles = self.rest_positions.copy()
        if randomize_pose:
            theta = self.rng.uniform(-np.pi / 6.0, np.pi / 6.0)
            rotation = np.array(
                [[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]],
                dtype=np.float32,
            )
            centered = particles - particles.mean(axis=0, keepdims=True)
            particles = centered @ rotation.T
            translation = np.array(
                [self.rng.uniform(-0.05, 0.05), self.rng.uniform(-0.05, 0.05), 0.0], dtype=np.float32
            )
            particles = particles + translation
            particles[:, 2] += self.rng.normal(0.0, randomize_height, size=particles.shape[0]).astype(np.float32)
        self.particles = particles
        self.relax(iterations=10, strength=0.2)
        return self.particles.copy()

    def copy_particles(self) -> np.ndarray:
        return self.particles.copy()

    def set_particles(self, particles: np.ndarray) -> None:
        self.particles = np.asarray(particles, dtype=np.float32).copy()

    def flatten(self) -> None:
        self.particles[:, :2] = 0.6 * self.particles[:, :2] + 0.4 * self.rest_positions[:, :2]
        self.particles[:, 2] = 0.0
        self.relax(iterations=8, strength=0.35)

    def shake(self, amplitude: float, frequency: float) -> None:
        timesteps = np.linspace(0.0, 1.0, 12, dtype=np.float32)
        for t in timesteps:
            offset = amplitude * np.sin(2.0 * np.pi * frequency * t)
            self.particles[:, 0] += offset
            self.relax(iterations=1, strength=0.1)
        self.particles[:, 0] *= 0.98

    def apply_pick_and_place(self, pick_pos: np.ndarray, place_pos: np.ndarray, height: float = 0.15) -> None:
        displacement = np.asarray(place_pos, dtype=np.float32) - np.asarray(pick_pos, dtype=np.float32)
        distances = np.linalg.norm(self.particles[:, :2] - np.asarray(pick_pos, dtype=np.float32)[:2], axis=1)
        weights = np.exp(-0.5 * (distances / max(self.grasp_radius, 1e-6)) ** 2).astype(np.float32)
        if weights.max() < 1e-6:
            nearest = np.argmin(np.linalg.norm(self.particles - pick_pos[None, :], axis=1))
            weights[nearest] = 1.0
        weights = weights / np.maximum(weights.max(), 1e-6)
        lift_profile = height * weights
        self.particles[:, :2] += weights[:, None] * displacement[:2][None, :]
        self.particles[:, 2] += lift_profile
        self.relax(iterations=8, strength=0.25)
        self.particles[:, 2] *= 0.55

    def relax(self, iterations: int = 4, strength: float = 0.2) -> None:
        for _ in range(iterations):
            updated = self.particles.copy()
            for idx, neighbors in enumerate(self.adjacency):
                if not neighbors:
                    continue
                neighbor_mean = self.particles[neighbors].mean(axis=0)
                updated[idx] = (1.0 - strength) * self.particles[idx] + strength * neighbor_mean
            updated[:, 2] = np.maximum(updated[:, 2], 0.0)
            self.particles = updated


class ClothFoldingEnv(gym.Env):
    """A Gym-compatible cloth folding environment with image and particle observations."""

    metadata = {"render_modes": ["rgb_array", "human"], "render_fps": 30}

    def __init__(
        self,
        config: Optional[ClothFoldingEnvConfig] = None,
        backend: Optional[str] = None,
        render_mode: str = "rgb_array",
    ) -> None:
        super().__init__()
        self.config = replace(config, backend=replace(config.backend, backend=backend)) if (config and backend) else (config or default_cloth_folding_config())
        self.config.validate()
        self.render_mode = render_mode
        self.minimum_gripper_height = 0.02
        self.grasp_height = self.config.grasp_height
        self.rng = np.random.default_rng(self.config.seed)
        self.mesh = generate_rectangular_cloth_mesh(
            width=self.config.mesh.width,
            height=self.config.mesh.height,
            resolution_x=self.config.mesh.resolution_x,
            resolution_y=self.config.mesh.resolution_y,
            origin=np.asarray(self.config.mesh.origin, dtype=np.float32),
        )
        self.backend_name, self.backend = self._create_backend(self.config.backend.backend)
        num_particles = self.mesh.vertices.shape[0]
        self.observation_space = spaces.Dict(
            {
                "rgb": spaces.Box(low=0, high=255, shape=(self.config.camera.height, self.config.camera.width, 3), dtype=np.uint8),
                "depth": spaces.Box(low=0.0, high=1.0, shape=(self.config.camera.height, self.config.camera.width), dtype=np.float32),
                "particles": spaces.Box(low=-np.inf, high=np.inf, shape=(num_particles, 3), dtype=np.float32),
            }
        )
        self.action_space = spaces.Box(low=-1.0, high=1.0, shape=(6,), dtype=np.float32)
        self.current_step = 0
        self.fold_type = self.config.supported_fold_types[0]
        self.goal_particles = self.mesh.vertices.copy()
        self.action_history: list[np.ndarray] = []
        self.reward_aggregator = composite_reward(
            {
                "chamfer": self.config.reward.chamfer_weight,
                "coverage": self.config.reward.coverage_weight,
                "fold_line": self.config.reward.fold_line_weight,
                "smoothness": -self.config.reward.smoothness_weight,
            }
        )

    @property
    def particles(self) -> np.ndarray:
        return self.backend.copy_particles()

    @particles.setter
    def particles(self, value: np.ndarray) -> None:
        self.backend.set_particles(value)

    def _create_backend(self, requested_backend: str) -> tuple[str, NumpyClothBackend]:
        if requested_backend == "softgym":
            try:
                __import__("softgym")
            except Exception:
                warnings.warn("SoftGym backend unavailable; falling back to numpy backend.", RuntimeWarning)
        elif requested_backend == "isaac_gym":
            try:
                __import__("isaacgym")
            except Exception:
                warnings.warn("Isaac Gym backend unavailable; falling back to numpy backend.", RuntimeWarning)
        return "numpy" if requested_backend != "numpy" else requested_backend, NumpyClothBackend(self.mesh, self.config.grasp_radius, self.rng)

    def is_position_collision_free(self, position: Sequence[float]) -> bool:
        point = np.asarray(position, dtype=np.float32)
        return bool(point[2] >= self.minimum_gripper_height)

    def reset(self, *, seed: Optional[int] = None, options: Optional[Dict[str, Any]] = None):
        if seed is not None:
            self.rng = np.random.default_rng(seed)
            self.backend.rng = self.rng
        try:
            super().reset(seed=seed)
        except TypeError:
            pass
        self.current_step = 0
        self.action_history.clear()
        options = options or {}
        self.fold_type = str(options.get("fold_type", self.rng.choice(self.config.supported_fold_types)))
        self.backend.reset(self.config.randomize_pose, self.config.randomize_height)
        self.goal_particles = self._sample_goal_particles(self.fold_type)
        observation = self._get_observation()
        info = self._build_info()
        return observation, info

    def step(self, action: np.ndarray):
        action = np.asarray(action, dtype=np.float32)
        if action.shape != (6,):
            raise ValueError("action must have shape (6,)")
        self.current_step += 1
        self.action_history.append(action.copy())
        world_action = self._denormalize_action(action)
        self.backend.apply_pick_and_place(world_action[:3], world_action[3:], height=self.grasp_height)
        observation = self._get_observation()
        reward = float(self._compute_reward())
        success = self._compute_success()
        terminated = success
        truncated = self.current_step >= self.config.max_episode_steps
        info = self._build_info()
        info.update({"success": success, "requested_backend": self.config.backend.backend, "active_backend": self.backend_name})
        return observation, reward, terminated, truncated, info

    def _denormalize_action(self, action: np.ndarray) -> np.ndarray:
        mins = self.mesh.vertices.min(axis=0)
        maxs = self.mesh.vertices.max(axis=0)
        span = np.maximum(maxs - mins, 1e-6)
        pick = mins + 0.5 * (action[:3] + 1.0) * span
        place = mins + 0.5 * (action[3:] + 1.0) * span
        pick[2] = max(pick[2], self.minimum_gripper_height)
        place[2] = max(place[2], self.minimum_gripper_height)
        return np.concatenate([pick, place], axis=0).astype(np.float32)

    def _sample_goal_particles(self, fold_type: str) -> np.ndarray:
        particles = self.mesh.vertices.copy()
        mins = particles.min(axis=0)
        maxs = particles.max(axis=0)
        center = 0.5 * (mins + maxs)
        if fold_type == "half_fold":
            mask = particles[:, 0] > center[0]
            particles[mask, 0] = 2.0 * center[0] - particles[mask, 0]
        elif fold_type == "quarter_fold":
            mask_x = particles[:, 0] > center[0]
            particles[mask_x, 0] = 2.0 * center[0] - particles[mask_x, 0]
            mask_y = particles[:, 1] > center[1]
            particles[mask_y, 1] = 2.0 * center[1] - particles[mask_y, 1]
        elif fold_type == "diagonal_fold":
            relative = particles[:, :2] - center[:2]
            swap_mask = relative[:, 0] > relative[:, 1]
            swapped = relative.copy()
            swapped[swap_mask] = swapped[swap_mask][:, ::-1]
            particles[:, :2] = swapped + center[:2]
        elif fold_type == "sleeve_fold":
            left = particles[:, 0] < mins[0] + 0.25 * (maxs[0] - mins[0])
            right = particles[:, 0] > mins[0] + 0.75 * (maxs[0] - mins[0])
            particles[left, 0] += 0.5 * (center[0] - particles[left, 0])
            particles[right, 0] -= 0.5 * (particles[right, 0] - center[0])
        particles[:, 2] = np.maximum(particles[:, 2], 0.0)
        return particles.astype(np.float32)

    def _goal_area(self) -> float:
        return float(self.config.mesh.width * self.config.mesh.height)

    def _fold_line(self) -> np.ndarray:
        mins = self.mesh.vertices.min(axis=0)
        maxs = self.mesh.vertices.max(axis=0)
        center = 0.5 * (mins + maxs)
        if self.fold_type in {"half_fold", "quarter_fold", "sleeve_fold"}:
            return np.asarray([[center[0], mins[1], 0.0], [center[0], maxs[1], 0.0]], dtype=np.float32)
        return np.asarray([[mins[0], mins[1], 0.0], [maxs[0], maxs[1], 0.0]], dtype=np.float32)

    def _compute_reward(self) -> torch.Tensor:
        current = self.particles
        rewards = {
            "chamfer": chamfer_distance_reward(current, self.goal_particles),
            "coverage": coverage_reward(current, self._goal_area()),
            "fold_line": fold_line_reward(current, self._fold_line()),
            "smoothness": smoothness_penalty(self.action_history) if self.action_history else torch.tensor(0.0),
        }
        return self.reward_aggregator(rewards)

    def _compute_success(self) -> bool:
        quality = fold_quality_metric(self.particles, goal_particles=self.goal_particles, fold_line=self._fold_line())
        return bool(quality["combined"] >= self.config.reward.success_threshold)

    def _build_info(self) -> Dict[str, Any]:
        particles = self.particles
        quality = fold_quality_metric(particles, goal_particles=self.goal_particles, fold_line=self._fold_line())
        return {
            "fold_type": self.fold_type,
            "chamfer_distance": chamfer_distance(particles, self.goal_particles),
            "iou": intersection_over_union(particles, self.goal_particles),
            "coverage": normalized_coverage(particles, self._goal_area()),
            "fold_quality": quality,
            "step": self.current_step,
        }

    def _get_observation(self) -> Dict[str, np.ndarray]:
        rgb, depth = self._render_arrays()
        return {"rgb": rgb, "depth": depth, "particles": self.particles.astype(np.float32)}

    def _render_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        height = self.config.camera.height
        width = self.config.camera.width
        rgb = np.full((height, width, 3), 255, dtype=np.uint8)
        depth = np.ones((height, width), dtype=np.float32)
        particles = self.particles
        mins = self.mesh.vertices.min(axis=0)[:2]
        maxs = self.mesh.vertices.max(axis=0)[:2]
        span = np.maximum(maxs - mins, 1e-6)
        normalized = (particles[:, :2] - mins) / span
        xs = np.clip((normalized[:, 0] * (width - 1)).astype(np.int32), 0, width - 1)
        ys = np.clip((normalized[:, 1] * (height - 1)).astype(np.int32), 0, height - 1)
        z_vals = particles[:, 2]
        if np.ptp(z_vals) < 1e-6:
            normalized_z = np.zeros_like(z_vals)
        else:
            normalized_z = (z_vals - z_vals.min()) / (z_vals.max() - z_vals.min())
        for x_idx, y_idx, z_norm in zip(xs, ys, normalized_z):
            y_slice = slice(max(y_idx - 1, 0), min(y_idx + 2, height))
            x_slice = slice(max(x_idx - 1, 0), min(x_idx + 2, width))
            color = np.array([50, 100 + int(80 * z_norm), 220 - int(60 * z_norm)], dtype=np.uint8)
            rgb[y_slice, x_slice] = color
            depth[y_slice, x_slice] = np.minimum(depth[y_slice, x_slice], 1.0 - float(z_norm))
        return rgb, depth

    def render(self):
        rgb, depth = self._render_arrays()
        if self.render_mode == "human":
            try:
                import matplotlib.pyplot as plt

                plt.figure("cloth-folding")
                plt.clf()
                plt.subplot(1, 2, 1)
                plt.imshow(rgb)
                plt.title("RGB")
                plt.axis("off")
                plt.subplot(1, 2, 2)
                plt.imshow(depth, cmap="viridis")
                plt.title("Depth")
                plt.axis("off")
                plt.pause(0.001)
            except Exception:
                warnings.warn("Human rendering requires matplotlib; returning RGB array instead.", RuntimeWarning)
        return rgb

    def execute_trajectory(self, trajectory: np.ndarray, grasp_at: np.ndarray, release_at: np.ndarray) -> None:
        del trajectory
        self.backend.apply_pick_and_place(grasp_at, release_at, height=self.grasp_height)

    def reset_to_rest_state(self) -> None:
        self.backend.set_particles(self.mesh.vertices.copy())

    def close(self) -> None:
        return None


__all__ = ["ClothFoldingEnv"]
