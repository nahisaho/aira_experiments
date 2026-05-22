from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Tuple

import math
import numpy as np


@dataclass
class ImuSample:
    timestamp: float
    accel: np.ndarray
    gyro: np.ndarray


@dataclass
class PreintegratedImu:
    delta_p: np.ndarray = field(default_factory=lambda: np.zeros(3))
    delta_v: np.ndarray = field(default_factory=lambda: np.zeros(3))
    delta_r: np.ndarray = field(default_factory=lambda: np.eye(3))
    dt: float = 0.0


class OnManifoldPreintegrator:
    """Minimal on-manifold IMU preintegration utility for architecture prototyping."""

    def __init__(self, gravity: np.ndarray | None = None) -> None:
        self.gravity = gravity if gravity is not None else np.array([0.0, 0.0, -9.80665])

    @staticmethod
    def _skew(v: np.ndarray) -> np.ndarray:
        return np.array(
            [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]],
            dtype=float,
        )

    def exp_so3(self, omega: np.ndarray) -> np.ndarray:
        theta = float(np.linalg.norm(omega))
        if theta < 1e-9:
            return np.eye(3) + self._skew(omega)
        axis = omega / theta
        axis_hat = self._skew(axis)
        return (
            np.eye(3)
            + math.sin(theta) * axis_hat
            + (1.0 - math.cos(theta)) * (axis_hat @ axis_hat)
        )

    def integrate(
        self,
        samples: Iterable[ImuSample],
        accel_bias: np.ndarray,
        gyro_bias: np.ndarray,
    ) -> PreintegratedImu:
        result = PreintegratedImu()
        previous: ImuSample | None = None
        for sample in samples:
            if previous is None:
                previous = sample
                continue
            dt = sample.timestamp - previous.timestamp
            if dt <= 0.0:
                previous = sample
                continue
            omega = 0.5 * (previous.gyro + sample.gyro) - gyro_bias
            acc = 0.5 * (previous.accel + sample.accel) - accel_bias
            d_r = self.exp_so3(omega * dt)
            result.delta_p += result.delta_v * dt + 0.5 * (result.delta_r @ acc) * dt * dt
            result.delta_v += (result.delta_r @ acc) * dt
            result.delta_r = result.delta_r @ d_r
            result.dt += dt
            previous = sample
        return result


def adaptive_feature_budget(
    texture_entropy: float,
    track_ratio: float,
    blur_score: float,
    min_features: int = 120,
    max_features: int = 260,
) -> int:
    entropy_term = 110.0 * texture_entropy
    track_term = 80.0 * track_ratio
    blur_penalty = 0.8 * max(0.0, 100.0 - blur_score)
    budget = min_features + entropy_term + track_term - blur_penalty
    return int(np.clip(budget, min_features, max_features))


def detect_degenerate_motion(
    mean_parallax_px: float,
    angular_rate_rad_s: float,
    linear_speed_m_s: float,
    hessian_condition_number: float,
) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    if mean_parallax_px < 1.5:
        reasons.append("low_parallax")
    if hessian_condition_number > 1e5:
        reasons.append("ill_conditioned_triangulation")
    ratio = angular_rate_rad_s / (linear_speed_m_s + 1e-3)
    if ratio > 5.0:
        reasons.append("pure_rotation_like")
    return (len(reasons) > 0, reasons)


def estimate_time_offset(
    image_timestamps: np.ndarray,
    imu_timestamps: np.ndarray,
    search_radius_s: float = 0.02,
) -> float:
    if image_timestamps.size == 0 or imu_timestamps.size == 0:
        return 0.0
    image_period = np.median(np.diff(image_timestamps)) if image_timestamps.size > 1 else 0.0
    imu_period = np.median(np.diff(imu_timestamps)) if imu_timestamps.size > 1 else 0.0
    raw_offset = float(image_timestamps[0] - imu_timestamps[0] - 0.5 * (image_period - imu_period))
    return float(np.clip(raw_offset, -search_radius_s, search_radius_s))
