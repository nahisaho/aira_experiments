"""
Module 4: GraspStabilityNet - Real-time Grasp Stability Evaluation
Temporal model that evaluates grasp stability from sequential
tactile observations, contact features, and force distributions.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class ContactFeatureExtractor(nn.Module):
    """Extract hand-crafted contact features from tactile data."""

    def __init__(self):
        super().__init__()
        # Learnable contact area segmentation
        self.contact_segmentor = nn.Sequential(
            nn.Conv2d(3, 32, 3, 1, 1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, 3, 2, 1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.Conv2d(64, 32, 3, 2, 1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.Conv2d(32, 1, 1),
            nn.Sigmoid(),
        )

        # Feature summarization
        self.feature_mlp = nn.Sequential(
            nn.Linear(7, 32),  # 7 hand-crafted features
            nn.ReLU(inplace=True),
            nn.Linear(32, 32),
        )

    def forward(self, tactile_img: torch.Tensor,
                force: torch.Tensor) -> torch.Tensor:
        """
        Args:
            tactile_img: (B, 3, H, W)
            force: (B, 6) wrench
        Returns:
            features: (B, 32) contact features
        """
        contact_mask = self.contact_segmentor(tactile_img)
        B = tactile_img.shape[0]

        # Compute features from mask
        contact_area = contact_mask.view(B, -1).mean(dim=1, keepdim=True)
        mask_flat = contact_mask.view(B, -1)
        # Center of pressure (normalized coordinates)
        H, W = contact_mask.shape[2], contact_mask.shape[3]
        y_coords = torch.linspace(0, 1, H, device=tactile_img.device)
        x_coords = torch.linspace(0, 1, W, device=tactile_img.device)
        yy, xx = torch.meshgrid(y_coords, x_coords, indexing='ij')
        yy_flat = yy.flatten().unsqueeze(0).expand(B, -1)
        xx_flat = xx.flatten().unsqueeze(0).expand(B, -1)
        total_pressure = mask_flat.sum(dim=1, keepdim=True).clamp(min=1e-6)
        cop_x = (mask_flat * xx_flat).sum(dim=1, keepdim=True) / total_pressure
        cop_y = (mask_flat * yy_flat).sum(dim=1, keepdim=True) / total_pressure

        # Force magnitude and normal force ratio
        force_mag = force[:, :3].norm(dim=1, keepdim=True)
        normal_ratio = force[:, 2:3].abs() / force_mag.clamp(min=1e-6)
        torque_mag = force[:, 3:].norm(dim=1, keepdim=True)

        features = torch.cat([
            contact_area, cop_x, cop_y,
            force_mag, normal_ratio, torque_mag,
            total_pressure / (H * W),
        ], dim=1)

        return self.feature_mlp(features)


class TemporalEncoder(nn.Module):
    """Bi-directional LSTM for temporal grasp signal encoding."""

    def __init__(self, input_dim: int, hidden_dim: int = 128, num_layers: int = 2):
        super().__init__()
        self.lstm = nn.LSTM(
            input_dim, hidden_dim, num_layers,
            batch_first=True, bidirectional=True, dropout=0.2
        )
        self.temporal_attn = nn.Sequential(
            nn.Linear(hidden_dim * 2, 64),
            nn.Tanh(),
            nn.Linear(64, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, T, D) temporal feature sequence
        Returns:
            output: (B, hidden_dim * 2) attended temporal features
        """
        lstm_out, _ = self.lstm(x)  # (B, T, 2*H)
        # Attention over timesteps
        attn_weights = self.temporal_attn(lstm_out).squeeze(-1)
        attn_weights = F.softmax(attn_weights, dim=1).unsqueeze(-1)
        attended = (lstm_out * attn_weights).sum(dim=1)
        return attended


class GraspStabilityNet(nn.Module):
    """
    Real-time grasp stability predictor.

    Processes temporal sequences of tactile observations to output:
    - Stability score (0-1)
    - Stability class (stable/marginal/unstable)
    - Time-to-failure estimate
    - Corrective action suggestion (increase/decrease/shift)
    """

    def __init__(
        self,
        tactile_feat_dim: int = 128,
        temporal_window: int = 10,
        lstm_hidden: int = 128,
    ):
        super().__init__()
        self.temporal_window = temporal_window

        # Tactile image encoder (lightweight)
        self.img_encoder = nn.Sequential(
            nn.Conv2d(3, 32, 5, 2, 2), nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, 3, 2, 1), nn.ReLU(inplace=True),
            nn.Conv2d(64, 128, 3, 2, 1), nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool2d(1), nn.Flatten(),
        )

        # Contact feature extractor
        self.contact_extractor = ContactFeatureExtractor()

        # Per-frame feature fusion
        frame_dim = 128 + 32 + 6  # img + contact + wrench
        self.frame_fuser = nn.Sequential(
            nn.Linear(frame_dim, 128),
            nn.ReLU(inplace=True),
            nn.Linear(128, 96),
        )

        # Temporal encoder
        self.temporal = TemporalEncoder(96, lstm_hidden, num_layers=2)

        # Delta features (frame-to-frame changes)
        self.delta_encoder = nn.Sequential(
            nn.Linear(96, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 48),
        )
        self.delta_temporal = TemporalEncoder(48, 64, num_layers=1)

        combined_dim = lstm_hidden * 2 + 64 * 2

        # Stability score head
        self.stability_head = nn.Sequential(
            nn.Linear(combined_dim, 128),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3),
            nn.Linear(128, 1),
            nn.Sigmoid(),
        )

        # Stability class head (stable/marginal/unstable)
        self.class_head = nn.Sequential(
            nn.Linear(combined_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 3),
        )

        # Time-to-failure head (regression, in seconds)
        self.ttf_head = nn.Sequential(
            nn.Linear(combined_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 1),
            nn.Softplus(),
        )

        # Corrective action head
        self.action_head = nn.Sequential(
            nn.Linear(combined_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 4),  # increase_force, decrease_force, shift_left, shift_right
        )

    def encode_frame(self, tactile_img: torch.Tensor,
                     wrench: torch.Tensor) -> torch.Tensor:
        img_feat = self.img_encoder(tactile_img)
        contact_feat = self.contact_extractor(tactile_img, wrench)
        return self.frame_fuser(torch.cat([img_feat, contact_feat, wrench], dim=1))

    def forward(self, tactile_sequence: torch.Tensor,
                wrench_sequence: torch.Tensor) -> dict:
        """
        Args:
            tactile_sequence: (B, T, 3, H, W) temporal tactile images
            wrench_sequence: (B, T, 6) temporal wrench readings
        Returns:
            dict with stability_score, stability_class, time_to_failure,
                 corrective_action
        """
        B, T, C, H, W = tactile_sequence.shape

        # Encode each frame
        frame_features = []
        for t in range(T):
            feat = self.encode_frame(
                tactile_sequence[:, t], wrench_sequence[:, t]
            )
            frame_features.append(feat)
        frame_features = torch.stack(frame_features, dim=1)  # (B, T, D)

        # Temporal encoding
        temporal_feat = self.temporal(frame_features)

        # Delta features (temporal differences)
        deltas = frame_features[:, 1:] - frame_features[:, :-1]
        delta_encoded = self.delta_encoder(deltas)
        delta_temporal = self.delta_temporal(delta_encoded)

        combined = torch.cat([temporal_feat, delta_temporal], dim=1)

        return {
            "stability_score": self.stability_head(combined).squeeze(-1),
            "stability_class": self.class_head(combined),
            "time_to_failure": self.ttf_head(combined).squeeze(-1),
            "corrective_action": self.action_head(combined),
        }


class StabilityLoss(nn.Module):
    """Multi-task loss for grasp stability prediction."""

    def __init__(self):
        super().__init__()
        # Learnable task weights (uncertainty-based)
        self.log_vars = nn.Parameter(torch.zeros(4))

    def forward(self, pred: dict, target: dict) -> dict:
        score_loss = F.binary_cross_entropy(
            pred["stability_score"], target["stability_score"]
        )
        class_loss = F.cross_entropy(
            pred["stability_class"], target["stability_class"]
        )
        ttf_loss = F.smooth_l1_loss(
            pred["time_to_failure"], target["time_to_failure"]
        )
        action_loss = F.cross_entropy(
            pred["corrective_action"], target["corrective_action"]
        )

        # Uncertainty-weighted multi-task loss
        losses = [score_loss, class_loss, ttf_loss, action_loss]
        total = sum(
            torch.exp(-self.log_vars[i]) * losses[i] + self.log_vars[i]
            for i in range(4)
        )

        return {
            "total": total,
            "score_loss": score_loss,
            "class_loss": class_loss,
            "ttf_loss": ttf_loss,
            "action_loss": action_loss,
        }
