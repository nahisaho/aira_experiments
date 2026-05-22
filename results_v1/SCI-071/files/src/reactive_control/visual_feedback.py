"""Visual feedback control using RGB-D observations and point-cloud alignment."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import torch


@dataclass(frozen=True)
class CameraIntrinsics:
    """Pinhole camera intrinsics."""

    fx: float
    fy: float
    cx: float
    cy: float


@dataclass(frozen=True)
class CameraPose:
    """Extrinsic pose mapping camera points to world coordinates."""

    rotation: torch.Tensor = field(default_factory=lambda: torch.eye(3, dtype=torch.float32))
    translation: torch.Tensor = field(default_factory=lambda: torch.zeros(3, dtype=torch.float32))

    def __post_init__(self) -> None:
        if self.rotation.shape != (3, 3):
            raise ValueError("rotation must have shape [3, 3].")
        if self.translation.shape != (3,):
            raise ValueError("translation must have shape [3].")


@dataclass
class CameraObservation:
    """Single RGB-D observation from one camera."""

    depth: torch.Tensor
    intrinsics: CameraIntrinsics
    rgb: Optional[torch.Tensor] = None
    pose: CameraPose = field(default_factory=CameraPose)
    valid_mask: Optional[torch.Tensor] = None


@dataclass
class PointCloud:
    """Point cloud with optional RGB colors."""

    points: torch.Tensor
    colors: Optional[torch.Tensor] = None

    def __post_init__(self) -> None:
        if self.points.ndim != 2 or self.points.shape[-1] != 3:
            raise ValueError("points must have shape [N, 3].")
        if self.colors is not None and self.colors.shape != self.points.shape:
            raise ValueError("colors must have shape [N, 3].")


@dataclass(frozen=True)
class ICPResult:
    """Result of ICP-based registration."""

    transform: torch.Tensor
    aligned_source: torch.Tensor
    rmse: float
    iterations: int


@dataclass(frozen=True)
class VisualServoingConfig:
    """Controller configuration for visual feedback."""

    jacobian_feature_dim: int = 6
    action_dim: int = 7
    control_gain: float = 0.5
    max_action_norm: float = 0.05
    icp_iterations: int = 15
    icp_tolerance: float = 1e-4
    downsample_points: int = 2048
    jacobian_regularization: float = 1e-3


class PointCloudProcessor:
    """Depth-to-point-cloud conversion and multi-camera fusion."""

    @staticmethod
    def from_depth(observation: CameraObservation) -> PointCloud:
        depth = observation.depth.to(dtype=torch.float32)
        if depth.ndim != 2:
            raise ValueError("depth must have shape [H, W].")
        height, width = depth.shape
        yy, xx = torch.meshgrid(
            torch.arange(height, device=depth.device),
            torch.arange(width, device=depth.device),
            indexing="ij",
        )
        mask = depth > 0.0
        if observation.valid_mask is not None:
            mask = mask & observation.valid_mask.to(dtype=torch.bool, device=depth.device)
        z = depth[mask]
        x = (xx[mask].to(dtype=depth.dtype) - observation.intrinsics.cx) * z / observation.intrinsics.fx
        y = (yy[mask].to(dtype=depth.dtype) - observation.intrinsics.cy) * z / observation.intrinsics.fy
        camera_points = torch.stack([x, y, z], dim=-1)
        world_points = camera_points @ observation.pose.rotation.transpose(0, 1) + observation.pose.translation

        colors = None
        if observation.rgb is not None:
            rgb = observation.rgb.to(dtype=torch.float32)
            if rgb.ndim != 3 or rgb.shape[0] != 3:
                raise ValueError("rgb must have shape [3, H, W].")
            colors = rgb.permute(1, 2, 0)[mask]
        return PointCloud(points=world_points, colors=colors)

    @staticmethod
    def fuse(observations: Sequence[CameraObservation], max_points: Optional[int] = None) -> PointCloud:
        clouds = [PointCloudProcessor.from_depth(observation) for observation in observations]
        if not clouds:
            raise ValueError("At least one observation is required for fusion.")
        points = torch.cat([cloud.points for cloud in clouds], dim=0)
        colors = None
        if all(cloud.colors is not None for cloud in clouds):
            colors = torch.cat([cloud.colors for cloud in clouds if cloud.colors is not None], dim=0)
        if max_points is not None and points.shape[0] > max_points:
            indices = torch.linspace(0, points.shape[0] - 1, max_points, device=points.device).long()
            points = points.index_select(0, indices)
            if colors is not None:
                colors = colors.index_select(0, indices)
        return PointCloud(points=points, colors=colors)


class ICPTracker:
    """Point-to-point ICP for deformable object tracking."""

    def __init__(self, max_iterations: int = 20, tolerance: float = 1e-4) -> None:
        self.max_iterations = max_iterations
        self.tolerance = tolerance

    def register(self, source: torch.Tensor, target: torch.Tensor) -> ICPResult:
        source = _ensure_point_cloud(source)
        target = _ensure_point_cloud(target)
        transform = torch.eye(4, device=source.device, dtype=source.dtype)
        transformed = source.clone()
        previous_error = float("inf")
        iteration = 0
        for iteration in range(1, self.max_iterations + 1):
            distances = torch.cdist(transformed, target)
            correspondences = target[distances.argmin(dim=1)]
            delta = _estimate_rigid_transform(transformed, correspondences)
            transformed = _apply_transform(transformed, delta)
            transform = delta @ transform
            rmse = torch.sqrt(((transformed - correspondences) ** 2).sum(dim=-1).mean()).item()
            if abs(previous_error - rmse) < self.tolerance:
                break
            previous_error = rmse
        return ICPResult(transform=transform, aligned_source=transformed, rmse=previous_error, iterations=iteration)


class DeformableJacobianEstimator:
    """Recursive least-squares estimator for a deformable interaction Jacobian."""

    def __init__(
        self,
        feature_dim: int,
        action_dim: int,
        *,
        regularization: float = 1e-3,
        forgetting_factor: float = 0.98,
        device: Optional[torch.device] = None,
    ) -> None:
        self.feature_dim = feature_dim
        self.action_dim = action_dim
        self.forgetting_factor = forgetting_factor
        self.device = device or torch.device("cpu")
        self.jacobian = torch.zeros(feature_dim, action_dim, device=self.device)
        self.covariance = torch.eye(action_dim, device=self.device) / regularization

    def update(self, action_delta: torch.Tensor, feature_delta: torch.Tensor) -> None:
        action = action_delta.to(device=self.device, dtype=torch.float32).view(-1, 1)
        feature = feature_delta.to(device=self.device, dtype=torch.float32).view(-1, 1)
        gain_denominator = self.forgetting_factor + (action.transpose(0, 1) @ self.covariance @ action).squeeze()
        gain = (self.covariance @ action) / gain_denominator.clamp_min(1e-8)
        prediction_error = feature - self.jacobian @ action
        self.jacobian = self.jacobian + prediction_error @ gain.transpose(0, 1)
        identity = torch.eye(self.action_dim, device=self.device)
        self.covariance = (identity - gain @ action.transpose(0, 1)) @ self.covariance / self.forgetting_factor

    def solve(self, feature_error: torch.Tensor, gain: float = 1.0, damping: float = 1e-3) -> torch.Tensor:
        error = feature_error.to(device=self.device, dtype=torch.float32).view(-1, 1)
        gram = self.jacobian.transpose(0, 1) @ self.jacobian
        damped_inverse = torch.linalg.solve(
            gram + damping * torch.eye(self.action_dim, device=self.device),
            self.jacobian.transpose(0, 1),
        )
        action = -gain * (damped_inverse @ error).squeeze(-1)
        return action


class VisualFeedbackController:
    """Visual servoing controller built on RGB-D point-cloud alignment."""

    def __init__(self, config: Optional[VisualServoingConfig] = None, *, device: Optional[torch.device] = None) -> None:
        self.config = config or VisualServoingConfig()
        self.device = device or torch.device("cpu")
        self.processor = PointCloudProcessor()
        self.tracker = ICPTracker(self.config.icp_iterations, self.config.icp_tolerance)
        self.jacobian_estimator = DeformableJacobianEstimator(
            self.config.jacobian_feature_dim,
            self.config.action_dim,
            regularization=self.config.jacobian_regularization,
            device=self.device,
        )

    def point_cloud_from_observations(self, observations: Sequence[CameraObservation]) -> PointCloud:
        return self.processor.fuse(observations, max_points=self.config.downsample_points)

    def compute_error(self, current: PointCloud, goal: PointCloud) -> tuple[torch.Tensor, ICPResult]:
        registration = self.tracker.register(current.points, goal.points)
        current_centroid = registration.aligned_source.mean(dim=0)
        goal_centroid = goal.points.mean(dim=0)
        translation_error = current_centroid - goal_centroid
        rotation_error = _rotation_vector_from_transform(registration.transform)
        feature_error = torch.cat([translation_error, rotation_error], dim=0)
        return feature_error, registration

    def compute_action(
        self,
        observations: Sequence[CameraObservation],
        goal: PointCloud,
        *,
        previous_action_delta: Optional[torch.Tensor] = None,
        previous_feature_error: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor, ICPResult]:
        current = self.point_cloud_from_observations(observations)
        feature_error, registration = self.compute_error(current, goal)
        if previous_action_delta is not None and previous_feature_error is not None:
            self.jacobian_estimator.update(previous_action_delta, feature_error - previous_feature_error)
        action = self.jacobian_estimator.solve(feature_error, gain=self.config.control_gain)
        action = self._clip_action(action)
        return action, feature_error, registration

    def _clip_action(self, action: torch.Tensor) -> torch.Tensor:
        norm = torch.linalg.norm(action)
        if norm <= self.config.max_action_norm:
            return action
        return action * (self.config.max_action_norm / norm.clamp_min(1e-8))


def _ensure_point_cloud(points: torch.Tensor) -> torch.Tensor:
    if points.ndim != 2 or points.shape[-1] != 3:
        raise ValueError("Point clouds must have shape [N, 3].")
    return points.to(dtype=torch.float32)


def _apply_transform(points: torch.Tensor, transform: torch.Tensor) -> torch.Tensor:
    rotation = transform[:3, :3]
    translation = transform[:3, 3]
    return points @ rotation.transpose(0, 1) + translation


def _estimate_rigid_transform(source: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    source_centroid = source.mean(dim=0)
    target_centroid = target.mean(dim=0)
    centered_source = source - source_centroid
    centered_target = target - target_centroid
    covariance = centered_source.transpose(0, 1) @ centered_target
    u, _, v_t = torch.linalg.svd(covariance)
    rotation = v_t.transpose(0, 1) @ u.transpose(0, 1)
    if torch.linalg.det(rotation) < 0:
        v_t[-1] *= -1
        rotation = v_t.transpose(0, 1) @ u.transpose(0, 1)
    translation = target_centroid - rotation @ source_centroid
    transform = torch.eye(4, device=source.device, dtype=source.dtype)
    transform[:3, :3] = rotation
    transform[:3, 3] = translation
    return transform


def _rotation_vector_from_transform(transform: torch.Tensor) -> torch.Tensor:
    rotation = transform[:3, :3]
    skew = 0.5 * torch.tensor(
        [
            rotation[2, 1] - rotation[1, 2],
            rotation[0, 2] - rotation[2, 0],
            rotation[1, 0] - rotation[0, 1],
        ],
        device=rotation.device,
        dtype=rotation.dtype,
    )
    return skew
