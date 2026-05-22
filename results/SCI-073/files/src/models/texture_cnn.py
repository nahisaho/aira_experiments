"""
Module 2: TextureCNN - Texture Classification from Tactile Images
Classifies surface textures using learned tactile features with
multi-scale spatial attention.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision import models


class SpatialAttention(nn.Module):
    """Channel-wise spatial attention for texture-discriminative regions."""

    def __init__(self, channels: int, reduction: int = 16):
        super().__init__()
        self.channel_att = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(),
            nn.Linear(channels, channels // reduction),
            nn.ReLU(inplace=True),
            nn.Linear(channels // reduction, channels),
            nn.Sigmoid(),
        )
        self.spatial_att = nn.Sequential(
            nn.Conv2d(2, 1, kernel_size=7, padding=3),
            nn.Sigmoid(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # Channel attention
        ca = self.channel_att(x).unsqueeze(-1).unsqueeze(-1)
        x = x * ca
        # Spatial attention
        avg_pool = x.mean(dim=1, keepdim=True)
        max_pool = x.max(dim=1, keepdim=True)[0]
        sa = self.spatial_att(torch.cat([avg_pool, max_pool], dim=1))
        return x * sa


class MultiScaleFeatureExtractor(nn.Module):
    """Extract tactile features at multiple spatial scales."""

    def __init__(self, in_channels: int = 3):
        super().__init__()
        self.scales = nn.ModuleList([
            nn.Sequential(
                nn.Conv2d(in_channels, 32, kernel_size=k, padding=k // 2),
                nn.BatchNorm2d(32),
                nn.ReLU(inplace=True),
            )
            for k in [3, 5, 7]
        ])
        self.fuse = nn.Conv2d(32 * 3, 64, 1)
        self.attention = SpatialAttention(64)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        features = [scale(x) for scale in self.scales]
        multi_scale = torch.cat(features, dim=1)
        fused = self.fuse(multi_scale)
        return self.attention(fused)


class GaborFilterBank(nn.Module):
    """Fixed Gabor filter bank for initial texture feature extraction."""

    def __init__(self, num_orientations: int = 8, num_frequencies: int = 4):
        super().__init__()
        self.num_filters = num_orientations * num_frequencies
        filters = self._create_gabor_filters(num_orientations, num_frequencies)
        self.register_buffer('filters', filters)

    def _create_gabor_filters(self, num_orient: int, num_freq: int,
                               ksize: int = 15) -> torch.Tensor:
        import numpy as np
        filters = []
        for theta_idx in range(num_orient):
            theta = theta_idx * np.pi / num_orient
            for freq_idx in range(num_freq):
                frequency = 0.1 + freq_idx * 0.1
                sigma = 3.0 + freq_idx
                half = ksize // 2
                y, x = np.mgrid[-half:half + 1, -half:half + 1].astype(np.float32)
                x_theta = x * np.cos(theta) + y * np.sin(theta)
                y_theta = -x * np.sin(theta) + y * np.cos(theta)
                gb = np.exp(-(x_theta ** 2 + y_theta ** 2) / (2 * sigma ** 2))
                gb *= np.cos(2 * np.pi * frequency * x_theta)
                gb = gb / (np.linalg.norm(gb) + 1e-7)
                filters.append(gb)

        filters = torch.tensor(np.stack(filters), dtype=torch.float32)
        # (N, 1, ksize, ksize) for depthwise-like conv
        return filters.unsqueeze(1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Apply Gabor filters to grayscale input."""
        gray = x.mean(dim=1, keepdim=True)
        # Apply each filter
        responses = F.conv2d(gray, self.filters, padding=self.filters.shape[-1] // 2)
        return responses  # (B, num_filters, H, W)


class TextureCNN(nn.Module):
    """
    Texture classification network combining:
    - Gabor filter bank (hand-crafted features)
    - Multi-scale learned features
    - EfficientNet backbone
    - Spatial attention fusion
    """

    def __init__(self, num_classes: int = 20, backbone: str = "efficientnet_b0",
                 input_size: tuple = (224, 224)):
        super().__init__()
        self.input_size = input_size

        # Branch 1: Gabor filters
        self.gabor = GaborFilterBank(num_orientations=8, num_frequencies=4)
        self.gabor_reducer = nn.Sequential(
            nn.Conv2d(32, 64, 3, 2, 1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool2d(7),
        )

        # Branch 2: Multi-scale extractor
        self.multi_scale = MultiScaleFeatureExtractor(in_channels=3)
        self.ms_pool = nn.AdaptiveAvgPool2d(7)

        # Branch 3: EfficientNet backbone
        if backbone == "efficientnet_b0":
            self.backbone = models.efficientnet_b0(pretrained=False)
            backbone_out = 1280
        else:
            self.backbone = models.resnet18(pretrained=False)
            backbone_out = 512
        self.backbone_features = nn.Sequential(
            *list(self.backbone.features if hasattr(self.backbone, 'features')
                  else list(self.backbone.children())[:-2])
        )

        # Fusion
        total_features = 64 + 64 + backbone_out
        self.classifier = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(),
            nn.Linear(total_features, 512),
            nn.ReLU(inplace=True),
            nn.Dropout(0.4),
            nn.Linear(512, 256),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3),
            nn.Linear(256, num_classes),
        )

        # Temperature scaling for calibration
        self.temperature = nn.Parameter(torch.ones(1) * 1.5)

    def forward(self, tactile_img: torch.Tensor) -> dict:
        """
        Args:
            tactile_img: (B, 3, H, W) tactile image
        Returns:
            dict with logits, probabilities, features
        """
        x = F.interpolate(tactile_img, self.input_size,
                          mode='bilinear', align_corners=False)

        # Gabor features
        gabor_feat = self.gabor(x)
        gabor_feat = self.gabor_reducer(gabor_feat)

        # Multi-scale features
        ms_feat = self.multi_scale(x)
        ms_feat = self.ms_pool(ms_feat)

        # Backbone features
        bb_feat = self.backbone_features(x)
        bb_feat = F.adaptive_avg_pool2d(bb_feat, 7)

        # Concatenate all branches
        combined = torch.cat([gabor_feat, ms_feat, bb_feat], dim=1)
        logits = self.classifier(combined)

        # Calibrated probabilities
        calibrated_logits = logits / self.temperature
        probs = F.softmax(calibrated_logits, dim=1)

        return {
            "logits": logits,
            "probabilities": probs,
            "features": combined.mean(dim=[2, 3]),
        }


class TextureLoss(nn.Module):
    """Label smoothing cross-entropy with optional contrastive regularization."""

    def __init__(self, num_classes: int = 20, label_smoothing: float = 0.1,
                 contrastive_weight: float = 0.1):
        super().__init__()
        self.ce_loss = nn.CrossEntropyLoss(label_smoothing=label_smoothing)
        self.contrastive_weight = contrastive_weight

    def forward(self, pred: dict, target: torch.Tensor) -> dict:
        ce = self.ce_loss(pred["logits"], target)

        # Supervised contrastive loss on feature space
        if self.contrastive_weight > 0 and pred["features"].shape[0] > 1:
            features = F.normalize(pred["features"], dim=1)
            sim_matrix = torch.matmul(features, features.T)
            labels_eq = (target.unsqueeze(0) == target.unsqueeze(1)).float()
            mask = 1.0 - torch.eye(target.shape[0], device=target.device)
            pos_sim = (sim_matrix * labels_eq * mask).sum(1)
            pos_count = (labels_eq * mask).sum(1).clamp(min=1)
            neg_sim = (sim_matrix * (1 - labels_eq) * mask)
            contrastive = -torch.log(
                torch.exp(pos_sim / pos_count) /
                (torch.exp(neg_sim).sum(1) + 1e-8)
            ).mean()
        else:
            contrastive = torch.tensor(0.0, device=target.device)

        total = ce + self.contrastive_weight * contrastive
        return {"total": total, "ce_loss": ce, "contrastive_loss": contrastive}
