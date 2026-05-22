"""Reality-gap analysis and randomization curricula."""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass, replace
import math
from typing import Any, Optional, Sequence

import torch

from .domain_randomization import DomainRandomizationConfig, IntRange, NumericRange
from .system_identification import Trajectory


@dataclass(frozen=True)
class GapSummary:
    """Aggregated reality-gap metrics across visual and dynamics spaces."""

    mmd: float
    fid: float
    trajectory_distance: float
    aggregate_score: float


@dataclass(frozen=True)
class ProgressiveRandomizationSchedule:
    """Progressively increases randomization strength during training."""

    start_scale: float = 0.25
    end_scale: float = 1.0
    total_steps: int = 100_000
    warmup_steps: int = 0
    schedule_type: str = "linear"

    def value(self, step: int) -> float:
        if step <= self.warmup_steps:
            return self.start_scale
        progress = min(max(step - self.warmup_steps, 0), self.total_steps)
        alpha = progress / max(self.total_steps, 1)
        if self.schedule_type == "linear":
            shaped = alpha
        elif self.schedule_type == "cosine":
            shaped = 0.5 - 0.5 * math.cos(math.pi * alpha)
        else:
            raise ValueError(f"Unsupported schedule_type: {self.schedule_type}")
        return self.start_scale + (self.end_scale - self.start_scale) * shaped


@dataclass(frozen=True)
class CurriculumStage:
    """Performance thresholded curriculum stage."""

    name: str
    randomization_scale: float
    promotion_threshold: float


class RealityGapAnalyzer:
    """Computes reality-gap metrics for visual and dynamical observations."""

    def __init__(self, mmd_bandwidths: Sequence[float] = (0.1, 1.0, 5.0)) -> None:
        self.mmd_bandwidths = tuple(mmd_bandwidths)

    def maximum_mean_discrepancy(self, real: torch.Tensor, simulated: torch.Tensor) -> float:
        real = _flatten_features(real)
        simulated = _flatten_features(simulated)
        xx = self._kernel(real, real)
        yy = self._kernel(simulated, simulated)
        xy = self._kernel(real, simulated)
        return float((xx.mean() + yy.mean() - 2.0 * xy.mean()).item())

    def frechet_inception_distance(self, real: torch.Tensor, simulated: torch.Tensor) -> float:
        real = _flatten_features(real)
        simulated = _flatten_features(simulated)
        mu_real = real.mean(dim=0)
        mu_sim = simulated.mean(dim=0)
        cov_real = _covariance(real)
        cov_sim = _covariance(simulated)
        cov_prod_sqrt = _sqrtm_psd(cov_real @ cov_sim)
        diff = mu_real - mu_sim
        fid = diff.dot(diff) + torch.trace(cov_real + cov_sim - 2.0 * cov_prod_sqrt)
        return float(fid.clamp_min(0.0).item())

    def trajectory_distance(self, real: Trajectory, simulated: Trajectory) -> float:
        steps = min(real.positions.shape[0], simulated.positions.shape[0])
        real_positions = real.positions[:steps]
        sim_positions = simulated.positions[:steps]
        position_distance = torch.linalg.norm(real_positions - sim_positions, dim=-1).mean()
        real_velocities = real.inferred_velocities()[:steps]
        sim_velocities = simulated.inferred_velocities()[:steps]
        velocity_distance = torch.linalg.norm(real_velocities - sim_velocities, dim=-1).mean()
        return float((position_distance + velocity_distance).item())

    def analyze(
        self,
        real_visual_features: torch.Tensor,
        simulated_visual_features: torch.Tensor,
        real_trajectory: Trajectory,
        simulated_trajectory: Trajectory,
        *,
        visual_weight: float = 0.5,
        dynamics_weight: float = 0.5,
    ) -> GapSummary:
        mmd = self.maximum_mean_discrepancy(real_visual_features, simulated_visual_features)
        fid = self.frechet_inception_distance(real_visual_features, simulated_visual_features)
        trajectory_distance = self.trajectory_distance(real_trajectory, simulated_trajectory)
        aggregate = visual_weight * (mmd + fid) + dynamics_weight * trajectory_distance
        return GapSummary(
            mmd=mmd,
            fid=fid,
            trajectory_distance=trajectory_distance,
            aggregate_score=aggregate,
        )

    def _kernel(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        distances = torch.cdist(x, y, p=2) ** 2
        kernels = [torch.exp(-0.5 * distances / (bandwidth**2)) for bandwidth in self.mmd_bandwidths]
        return sum(kernels) / len(kernels)


class RealityGapScheduler:
    """Combines progressive randomization and curriculum scheduling."""

    def __init__(
        self,
        base_config: DomainRandomizationConfig,
        schedule: ProgressiveRandomizationSchedule,
        curriculum: Optional[Sequence[CurriculumStage]] = None,
    ) -> None:
        self.base_config = base_config
        self.schedule = schedule
        self.curriculum = sorted(curriculum or [], key=lambda stage: stage.randomization_scale)
        self._stage_index = 0

    @property
    def current_stage(self) -> Optional[CurriculumStage]:
        if not self.curriculum:
            return None
        return self.curriculum[self._stage_index]

    def update_curriculum(self, performance: float) -> Optional[CurriculumStage]:
        if not self.curriculum:
            return None
        while self._stage_index + 1 < len(self.curriculum):
            next_stage = self.curriculum[self._stage_index + 1]
            if performance < next_stage.promotion_threshold:
                break
            self._stage_index += 1
        return self.curriculum[self._stage_index]

    def config_for_step(self, step: int, performance: Optional[float] = None) -> DomainRandomizationConfig:
        if performance is not None:
            self.update_curriculum(performance)
        scale = self.schedule.value(step)
        if self.current_stage is not None:
            scale = max(scale, self.current_stage.randomization_scale)
        return _scale_config(self.base_config, scale)


def _scale_config(config: DomainRandomizationConfig, factor: float) -> DomainRandomizationConfig:
    return DomainRandomizationConfig(
        material=_scale_dataclass(config.material, factor),
        visual=_scale_dataclass(config.visual, factor),
        dynamics=_scale_dataclass(config.dynamics, factor),
        noise=config.noise,
        adr=config.adr,
    )


def _scale_dataclass(instance: Any, factor: float) -> Any:
    if not is_dataclass(instance):
        return instance
    updates: dict[str, Any] = {}
    for item in fields(instance):
        value = getattr(instance, item.name)
        if isinstance(value, NumericRange):
            updates[item.name] = value.scale(factor)
        elif isinstance(value, IntRange):
            updates[item.name] = value.scale(factor)
        elif is_dataclass(value):
            updates[item.name] = _scale_dataclass(value, factor)
        else:
            updates[item.name] = value
    return replace(instance, **updates)


def _flatten_features(features: torch.Tensor) -> torch.Tensor:
    if features.ndim == 1:
        features = features.unsqueeze(0)
    if features.ndim > 2:
        features = features.reshape(features.shape[0], -1)
    return features.to(dtype=torch.float32)


def _covariance(features: torch.Tensor) -> torch.Tensor:
    centered = features - features.mean(dim=0, keepdim=True)
    denom = max(features.shape[0] - 1, 1)
    return centered.transpose(0, 1) @ centered / denom


def _sqrtm_psd(matrix: torch.Tensor) -> torch.Tensor:
    eigenvalues, eigenvectors = torch.linalg.eig(matrix)
    eigenvalues = eigenvalues.real.clamp_min(0.0)
    eigenvectors = eigenvectors.real
    sqrt_matrix = eigenvectors @ torch.diag(eigenvalues.sqrt()) @ torch.linalg.pinv(eigenvectors)
    return 0.5 * (sqrt_matrix + sqrt_matrix.transpose(0, 1))
