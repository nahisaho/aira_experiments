"""PyTorch actor-critic policy network for reactive deformable manipulation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
from torch import nn
from torch.distributions import Categorical, Normal


@dataclass
class ObservationBatch:
    """Multi-modal observation bundle for policy evaluation."""

    rgb: Optional[torch.Tensor] = None
    depth: Optional[torch.Tensor] = None
    point_cloud: Optional[torch.Tensor] = None
    proprio: Optional[torch.Tensor] = None


@dataclass(frozen=True)
class PolicyNetworkConfig:
    """Configuration for the actor-critic reactive policy."""

    visual_feature_dim: int = 256
    point_feature_dim: int = 256
    proprio_feature_dim: int = 128
    fused_feature_dim: int = 256
    continuous_action_dim: int = 7
    discrete_action_dim: int = 0
    point_cloud_channels: int = 3
    proprio_dim: int = 14
    min_log_std: float = -5.0
    max_log_std: float = 2.0

    @property
    def hybrid_action(self) -> bool:
        return self.discrete_action_dim > 0


@dataclass
class PolicyOutput:
    """Forward-pass outputs for actor and critic heads."""

    value: torch.Tensor
    continuous_mean: torch.Tensor
    continuous_log_std: torch.Tensor
    discrete_logits: Optional[torch.Tensor]
    latent: torch.Tensor


@dataclass
class ActionSample:
    """Sampled action and associated policy statistics."""

    continuous_action: torch.Tensor
    discrete_action: Optional[torch.Tensor]
    log_prob: torch.Tensor
    entropy: torch.Tensor


class CNNEncoder(nn.Module):
    """CNN encoder for stacked RGB-D observations."""

    def __init__(self, in_channels: int, feature_dim: int) -> None:
        super().__init__()
        self.backbone = nn.Sequential(
            nn.Conv2d(in_channels, 32, kernel_size=5, stride=2, padding=2),
            nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, kernel_size=3, stride=2, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(64, 128, kernel_size=3, stride=2, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(128, 256, kernel_size=3, stride=2, padding=1),
            nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool2d((1, 1)),
            nn.Flatten(),
            nn.Linear(256, feature_dim),
            nn.LayerNorm(feature_dim),
            nn.ReLU(inplace=True),
        )

    def forward(self, images: torch.Tensor) -> torch.Tensor:
        return self.backbone(images)


class PointNetEncoder(nn.Module):
    """PointNet-style encoder for point-cloud observations."""

    def __init__(self, input_dim: int, feature_dim: int) -> None:
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 128),
            nn.ReLU(inplace=True),
            nn.Linear(128, feature_dim),
            nn.ReLU(inplace=True),
        )
        self.projection = nn.Sequential(
            nn.Linear(feature_dim, feature_dim),
            nn.LayerNorm(feature_dim),
            nn.ReLU(inplace=True),
        )

    def forward(self, point_cloud: torch.Tensor) -> torch.Tensor:
        if point_cloud.ndim != 3:
            raise ValueError("point_cloud must have shape [B, N, C].")
        point_features = self.mlp(point_cloud)
        pooled = point_features.max(dim=1).values
        return self.projection(pooled)


class MLPEncoder(nn.Module):
    """Lightweight MLP encoder for proprioception."""

    def __init__(self, input_dim: int, feature_dim: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, feature_dim),
            nn.ReLU(inplace=True),
            nn.Linear(feature_dim, feature_dim),
            nn.LayerNorm(feature_dim),
            nn.ReLU(inplace=True),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class ReactivePolicyNetwork(nn.Module):
    """Actor-critic network with separate visual and proprioceptive streams."""

    def __init__(self, config: Optional[PolicyNetworkConfig] = None) -> None:
        super().__init__()
        self.config = config or PolicyNetworkConfig()
        self.visual_encoder = CNNEncoder(in_channels=4, feature_dim=self.config.visual_feature_dim)
        self.pointnet_encoder = PointNetEncoder(
            input_dim=self.config.point_cloud_channels,
            feature_dim=self.config.point_feature_dim,
        )
        self.proprio_encoder = MLPEncoder(
            input_dim=self.config.proprio_dim,
            feature_dim=self.config.proprio_feature_dim,
        )

        fused_input_dim = self.config.fused_feature_dim + self.config.proprio_feature_dim
        self.visual_fusion = nn.Sequential(
            nn.Linear(self.config.visual_feature_dim + self.config.point_feature_dim, self.config.fused_feature_dim),
            nn.ReLU(inplace=True),
            nn.Linear(self.config.fused_feature_dim, self.config.fused_feature_dim),
            nn.LayerNorm(self.config.fused_feature_dim),
            nn.ReLU(inplace=True),
        )
        self.actor_trunk = nn.Sequential(
            nn.Linear(fused_input_dim, self.config.fused_feature_dim),
            nn.ReLU(inplace=True),
            nn.Linear(self.config.fused_feature_dim, self.config.fused_feature_dim),
            nn.ReLU(inplace=True),
        )
        self.critic_trunk = nn.Sequential(
            nn.Linear(fused_input_dim, self.config.fused_feature_dim),
            nn.ReLU(inplace=True),
            nn.Linear(self.config.fused_feature_dim, self.config.fused_feature_dim),
            nn.ReLU(inplace=True),
        )
        self.continuous_mean = nn.Linear(self.config.fused_feature_dim, self.config.continuous_action_dim)
        self.continuous_log_std = nn.Linear(self.config.fused_feature_dim, self.config.continuous_action_dim)
        self.discrete_head = (
            nn.Linear(self.config.fused_feature_dim, self.config.discrete_action_dim)
            if self.config.hybrid_action
            else None
        )
        self.value_head = nn.Linear(self.config.fused_feature_dim, 1)

    def forward(self, observations: ObservationBatch) -> PolicyOutput:
        visual_latent = self._encode_visual(observations)
        proprio_latent = self._encode_proprio(observations)
        latent = torch.cat([visual_latent, proprio_latent], dim=-1)
        actor_latent = self.actor_trunk(latent)
        critic_latent = self.critic_trunk(latent)

        continuous_mean = self.continuous_mean(actor_latent)
        continuous_log_std = self.continuous_log_std(actor_latent).clamp(
            min=self.config.min_log_std,
            max=self.config.max_log_std,
        )
        discrete_logits = self.discrete_head(actor_latent) if self.discrete_head is not None else None
        value = self.value_head(critic_latent).squeeze(-1)
        return PolicyOutput(
            value=value,
            continuous_mean=continuous_mean,
            continuous_log_std=continuous_log_std,
            discrete_logits=discrete_logits,
            latent=actor_latent,
        )

    def sample_action(self, observations: ObservationBatch, deterministic: bool = False) -> ActionSample:
        output = self.forward(observations)
        continuous_dist = Normal(output.continuous_mean, output.continuous_log_std.exp())
        if deterministic:
            continuous_raw = output.continuous_mean
        else:
            continuous_raw = continuous_dist.rsample()
        continuous_action = torch.tanh(continuous_raw)
        continuous_log_prob = continuous_dist.log_prob(continuous_raw).sum(dim=-1)
        squash_correction = torch.log1p(-continuous_action.pow(2) + 1e-6).sum(dim=-1)
        log_prob = continuous_log_prob - squash_correction
        entropy = continuous_dist.entropy().sum(dim=-1)

        discrete_action = None
        if output.discrete_logits is not None:
            discrete_dist = Categorical(logits=output.discrete_logits)
            discrete_action = torch.argmax(output.discrete_logits, dim=-1) if deterministic else discrete_dist.sample()
            log_prob = log_prob + discrete_dist.log_prob(discrete_action)
            entropy = entropy + discrete_dist.entropy()
        return ActionSample(
            continuous_action=continuous_action,
            discrete_action=discrete_action,
            log_prob=log_prob,
            entropy=entropy,
        )

    def evaluate_actions(
        self,
        observations: ObservationBatch,
        continuous_actions: torch.Tensor,
        discrete_actions: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        output = self.forward(observations)
        continuous_actions = continuous_actions.clamp(-0.999, 0.999)
        unsquashed = 0.5 * (torch.log1p(continuous_actions) - torch.log1p(-continuous_actions))
        continuous_dist = Normal(output.continuous_mean, output.continuous_log_std.exp())
        log_prob = continuous_dist.log_prob(unsquashed).sum(dim=-1)
        log_prob = log_prob - torch.log1p(-continuous_actions.pow(2) + 1e-6).sum(dim=-1)
        entropy = continuous_dist.entropy().sum(dim=-1)
        if output.discrete_logits is not None:
            if discrete_actions is None:
                raise ValueError("discrete_actions must be provided for hybrid action policies.")
            discrete_dist = Categorical(logits=output.discrete_logits)
            log_prob = log_prob + discrete_dist.log_prob(discrete_actions)
            entropy = entropy + discrete_dist.entropy()
        return log_prob, entropy, output.value

    def _encode_visual(self, observations: ObservationBatch) -> torch.Tensor:
        if observations.rgb is None or observations.depth is None:
            raise ValueError("Both rgb and depth observations are required.")
        rgb = observations.rgb.to(dtype=torch.float32)
        depth = observations.depth.to(dtype=torch.float32)
        if rgb.ndim != 4 or rgb.shape[1] != 3:
            raise ValueError("rgb must have shape [B, 3, H, W].")
        if depth.ndim != 4 or depth.shape[1] != 1:
            raise ValueError("depth must have shape [B, 1, H, W].")
        stacked = torch.cat([rgb, depth], dim=1)
        visual_latent = self.visual_encoder(stacked)
        if observations.point_cloud is None:
            raise ValueError("point_cloud observations are required.")
        point_latent = self.pointnet_encoder(observations.point_cloud.to(dtype=torch.float32))
        return self.visual_fusion(torch.cat([visual_latent, point_latent], dim=-1))

    def _encode_proprio(self, observations: ObservationBatch) -> torch.Tensor:
        if observations.proprio is None:
            raise ValueError("proprio observations are required.")
        return self.proprio_encoder(observations.proprio.to(dtype=torch.float32))
