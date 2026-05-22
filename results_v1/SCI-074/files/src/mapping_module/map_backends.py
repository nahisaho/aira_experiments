from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import math
import numpy as np
import yaml


@dataclass
class TsdfVoxel:
    distance: float = 0.0
    weight: float = 0.0


class VdbFusionBackend:
    """Sparse Python prototype approximating tile-based TSDF behavior."""

    def __init__(self, voxel_size: float, truncation_distance: float, max_weight: float) -> None:
        self.voxel_size = voxel_size
        self.truncation_distance = truncation_distance
        self.max_weight = max_weight
        self._grid: Dict[Tuple[int, int, int], TsdfVoxel] = {}

    def _index(self, point: np.ndarray) -> Tuple[int, int, int]:
        return tuple(np.floor(point / self.voxel_size).astype(int).tolist())

    def integrate_points(self, points: np.ndarray, sensor_origin: np.ndarray, pose_confidence: float = 1.0) -> None:
        for point in points:
            ray = point - sensor_origin
            distance = float(np.linalg.norm(ray))
            if distance < 1e-6:
                continue
            sdf = min(self.truncation_distance, distance)
            weight = min(self.max_weight, pose_confidence * math.exp(-0.08 * distance * distance))
            idx = self._index(point)
            voxel = self._grid.get(idx, TsdfVoxel())
            voxel.distance = (voxel.distance * voxel.weight + sdf * weight) / (voxel.weight + weight)
            voxel.weight = min(self.max_weight, voxel.weight + weight)
            self._grid[idx] = voxel

    def occupancy_points(self, threshold: float = 0.7) -> np.ndarray:
        occupied = []
        for idx, voxel in self._grid.items():
            if voxel.weight >= threshold:
                occupied.append(np.array(idx, dtype=float) * self.voxel_size)
        return np.array(occupied) if occupied else np.zeros((0, 3))

    def serialize(self, path: str) -> None:
        payload = {
            "voxel_size": self.voxel_size,
            "truncation_distance": self.truncation_distance,
            "max_weight": self.max_weight,
            "voxels": [
                {"index": list(idx), "distance": voxel.distance, "weight": voxel.weight}
                for idx, voxel in self._grid.items()
            ],
        }
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with Path(path).open("w", encoding="utf-8") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False)


class OctomapBackend:
    """Sparse occupancy prototype using log-odds updates."""

    def __init__(self, resolution: float, prob_hit: float, prob_miss: float) -> None:
        self.resolution = resolution
        self.log_hit = math.log(prob_hit / (1.0 - prob_hit))
        self.log_miss = math.log(prob_miss / (1.0 - prob_miss))
        self.occupancy: Dict[Tuple[int, int, int], float] = {}

    def _index(self, point: np.ndarray) -> Tuple[int, int, int]:
        return tuple(np.floor(point / self.resolution).astype(int).tolist())

    def insert_point_cloud(self, points: Iterable[np.ndarray], sensor_origin: np.ndarray) -> None:
        for point in points:
            idx = self._index(np.asarray(point, dtype=float))
            self.occupancy[idx] = self.occupancy.get(idx, 0.0) + self.log_hit
            free_mid = sensor_origin + 0.5 * (np.asarray(point, dtype=float) - sensor_origin)
            free_idx = self._index(free_mid)
            self.occupancy[free_idx] = self.occupancy.get(free_idx, 0.0) - abs(self.log_miss)

    def occupied_points(self, threshold: float = 0.0) -> np.ndarray:
        occupied = [np.array(idx, dtype=float) * self.resolution for idx, value in self.occupancy.items() if value > threshold]
        return np.array(occupied) if occupied else np.zeros((0, 3))

    def serialize(self, path: str) -> None:
        payload = {
            "resolution": self.resolution,
            "occupancy": [{"index": list(idx), "log_odds": value} for idx, value in self.occupancy.items()],
        }
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with Path(path).open("w", encoding="utf-8") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False)
