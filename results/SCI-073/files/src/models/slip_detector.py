"""
Module 5: SlipDetector - Slip Detection & Force Control Feedback
Detects incipient and gross slip events from tactile image sequences
and provides reactive force control adjustments.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import math


class OpticalFlowEstimator(nn.Module):
    """Lightweight optical flow estimation between consecutive tactile frames."""

    def __init__(self, in_channels: int = 6):
        super().__init__()
        # FlowNetS-inspired architecture (simplified)
        self.encoder = nn.Sequential(
            nn.Conv2d(in_channels, 32, 7, 2, 3), nn.LeakyReLU(0.1, inplace=True),
            nn.Conv2d(32, 64, 5, 2, 2), nn.LeakyReLU(0.1, inplace=True),
            nn.Conv2d(64, 128, 5, 2, 2), nn.LeakyReLU(0.1, inplace=True),
            nn.Conv2d(128, 256, 3, 2, 1), nn.LeakyReLU(0.1, inplace=True),
        )
        self.decoder = nn.Sequential(
            nn.ConvTranspose2d(256, 128, 4, 2, 1), nn.LeakyReLU(0.1, inplace=True),
            nn.ConvTranspose2d(128, 64, 4, 2, 1), nn.LeakyReLU(0.1, inplace=True),
            nn.ConvTranspose2d(64, 32, 4, 2, 1), nn.LeakyReLU(0.1, inplace=True),
            nn.ConvTranspose2d(32, 2, 4, 2, 1),  # 2-channel flow (dx, dy)
        )

    def forward(self, frame1: torch.Tensor, frame2: torch.Tensor) -> torch.Tensor:
        """
        Args:
            frame1, frame2: (B, 3, H, W)
        Returns:
            flow: (B, 2, H, W) optical flow field
        """
        x = torch.cat([frame1, frame2], dim=1)
        encoded = self.encoder(x)
        flow = self.decoder(encoded)
        flow = F.interpolate(flow, frame1.shape[2:], mode='bilinear',
                             align_corners=False)
        return flow


class MarkerTracker(nn.Module):
    """Track gel surface markers for precise displacement measurement."""

    def __init__(self, num_markers: int = 100, feature_dim: int = 64):
        super().__init__()
        self.num_markers = num_markers

        # Marker detection head
        self.detector = nn.Sequential(
            nn.Conv2d(3, 32, 3, 1, 1), nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, 3, 2, 1), nn.ReLU(inplace=True),
            nn.Conv2d(64, 64, 3, 1, 1), nn.ReLU(inplace=True),
            nn.Conv2d(64, 1, 1), nn.Sigmoid(),
        )

        # Displacement regression per marker
        self.displacement_head = nn.Sequential(
            nn.Conv2d(64, 32, 3, 1, 1), nn.ReLU(inplace=True),
            nn.Conv2d(32, 2, 1),  # dx, dy per pixel
        )

        self.feature_extractor = nn.Sequential(
            nn.Conv2d(3, 32, 3, 1, 1), nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, 3, 2, 1), nn.ReLU(inplace=True),
        )

    def forward(self, frame1: torch.Tensor, frame2: torch.Tensor) -> dict:
        # Detect markers in frame1
        heatmap = self.detector(frame1)
        # Compute features for displacement
        feat1 = self.feature_extractor(frame1)
        feat2 = self.feature_extractor(frame2)
        diff = feat2 - feat1
        displacements = self.displacement_head(diff)
        displacements = F.interpolate(displacements, frame1.shape[2:],
                                       mode='bilinear', align_corners=False)
        return {"heatmap": heatmap, "displacements": displacements}


class SlipClassifier(nn.Module):
    """Classify slip state from flow and marker displacement features."""

    def __init__(self, flow_feature_dim: int = 128):
        super().__init__()
        # Flow statistics encoder
        self.flow_encoder = nn.Sequential(
            nn.Conv2d(2, 32, 3, 2, 1), nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, 3, 2, 1), nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool2d(4),
            nn.Flatten(),
            nn.Linear(64 * 16, flow_feature_dim),
        )

        # Temporal slip classifier (over flow history)
        self.temporal = nn.GRU(
            flow_feature_dim + 8, 128,
            num_layers=2, batch_first=True, dropout=0.2
        )

        self.classifier = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3),
            nn.Linear(64, 4),  # no_slip, incipient_slip, translational_slip, rotational_slip
        )

        self.magnitude_head = nn.Sequential(
            nn.Linear(128, 32),
            nn.ReLU(inplace=True),
            nn.Linear(32, 1),
            nn.Sigmoid(),
        )

        self.direction_head = nn.Sequential(
            nn.Linear(128, 32),
            nn.ReLU(inplace=True),
            nn.Linear(32, 2),  # slip direction (dx, dy)
        )

    def compute_flow_statistics(self, flow: torch.Tensor) -> torch.Tensor:
        """Compute statistical features from optical flow field."""
        B = flow.shape[0]
        mag = flow.norm(dim=1)  # (B, H, W)
        angle = torch.atan2(flow[:, 1], flow[:, 0])

        stats = torch.stack([
            mag.mean(dim=[1, 2]),
            mag.std(dim=[1, 2]),
            mag.max(dim=2)[0].max(dim=1)[0],
            (mag > 0.1).float().mean(dim=[1, 2]),  # slip area ratio
            angle.mean(dim=[1, 2]),
            angle.std(dim=[1, 2]),
            flow[:, 0].mean(dim=[1, 2]),  # mean dx
            flow[:, 1].mean(dim=[1, 2]),  # mean dy
        ], dim=1)  # (B, 8)
        return stats

    def forward(self, flow_sequence: torch.Tensor) -> dict:
        """
        Args:
            flow_sequence: (B, T, 2, H, W) temporal flow fields
        Returns:
            dict with slip_class, slip_magnitude, slip_direction
        """
        B, T, _, H, W = flow_sequence.shape

        flow_features = []
        for t in range(T):
            feat = self.flow_encoder(flow_sequence[:, t])
            stats = self.compute_flow_statistics(flow_sequence[:, t])
            flow_features.append(torch.cat([feat, stats], dim=1))
        flow_features = torch.stack(flow_features, dim=1)

        temporal_out, _ = self.temporal(flow_features)
        last_hidden = temporal_out[:, -1]

        return {
            "slip_class": self.classifier(last_hidden),
            "slip_magnitude": self.magnitude_head(last_hidden).squeeze(-1),
            "slip_direction": self.direction_head(last_hidden),
        }


class ImpedanceForceController(nn.Module):
    """
    Neural impedance controller for reactive force adjustment.
    Combines classical impedance control with learned residuals.
    """

    def __init__(self, kp: float = 50.0, kd: float = 5.0, ki: float = 0.1,
                 max_force: float = 20.0, min_force: float = 0.5):
        super().__init__()
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.max_force = max_force
        self.min_force = min_force

        # Learned residual controller
        self.residual_net = nn.Sequential(
            nn.Linear(10, 64),  # slip_info(4) + force(3) + target_force(3)
            nn.ReLU(inplace=True),
            nn.Linear(64, 32),
            nn.ReLU(inplace=True),
            nn.Linear(32, 3),  # force correction (dx, dy, dz)
            nn.Tanh(),
        )
        self.residual_scale = nn.Parameter(torch.tensor(2.0))

    def forward(self, current_force: torch.Tensor, target_force: torch.Tensor,
                force_error_integral: torch.Tensor,
                slip_info: torch.Tensor) -> dict:
        """
        Args:
            current_force: (B, 3) current contact force
            target_force: (B, 3) desired force
            force_error_integral: (B, 3) integrated force error
            slip_info: (B, 4) [slip_magnitude, slip_dx, slip_dy, slip_class_idx]
        Returns:
            dict with force_command, impedance_params
        """
        error = target_force - current_force

        # Classical PID
        pid_output = self.kp * error + self.ki * force_error_integral

        # Learned residual
        residual_input = torch.cat([slip_info, current_force, target_force], dim=1)
        residual = self.residual_net(residual_input) * self.residual_scale

        # Slip-reactive force increase
        slip_magnitude = slip_info[:, 0:1]
        slip_gain = 1.0 + 2.0 * slip_magnitude  # Increase force on slip
        force_command = pid_output * slip_gain + residual

        # Clamp to safety limits
        force_magnitude = force_command.norm(dim=1, keepdim=True)
        scale = torch.clamp(
            self.max_force / force_magnitude.clamp(min=1e-6), max=1.0
        )
        force_command = force_command * scale

        return {
            "force_command": force_command,
            "pid_component": pid_output,
            "residual_component": residual,
            "slip_gain": slip_gain.squeeze(-1),
        }


class SlipDetector(nn.Module):
    """
    Complete slip detection and force control module.
    Integrates optical flow estimation, marker tracking,
    slip classification, and impedance force control.
    """

    def __init__(self, **controller_kwargs):
        super().__init__()
        self.flow_estimator = OpticalFlowEstimator()
        self.marker_tracker = MarkerTracker()
        self.slip_classifier = SlipClassifier()
        self.force_controller = ImpedanceForceController(**controller_kwargs)

    def detect_slip(self, tactile_sequence: torch.Tensor) -> dict:
        """
        Args:
            tactile_sequence: (B, T, 3, H, W) temporal tactile images
        Returns:
            dict with flow_sequence, slip_class, slip_magnitude, slip_direction
        """
        B, T, C, H, W = tactile_sequence.shape
        flows = []
        for t in range(T - 1):
            flow = self.flow_estimator(
                tactile_sequence[:, t], tactile_sequence[:, t + 1]
            )
            flows.append(flow)
        flow_seq = torch.stack(flows, dim=1)

        slip_result = self.slip_classifier(flow_seq)
        slip_result["flow_sequence"] = flow_seq
        return slip_result

    def compute_force_response(self, slip_result: dict,
                                current_force: torch.Tensor,
                                target_force: torch.Tensor,
                                error_integral: torch.Tensor) -> dict:
        """Compute force control response based on slip detection."""
        slip_info = torch.cat([
            slip_result["slip_magnitude"].unsqueeze(1),
            slip_result["slip_direction"],
            slip_result["slip_class"].argmax(dim=1, keepdim=True).float(),
        ], dim=1)

        return self.force_controller(
            current_force, target_force, error_integral, slip_info
        )

    def forward(self, tactile_sequence: torch.Tensor,
                current_force: torch.Tensor,
                target_force: torch.Tensor,
                error_integral: torch.Tensor) -> dict:
        slip = self.detect_slip(tactile_sequence)
        control = self.compute_force_response(
            slip, current_force, target_force, error_integral
        )
        return {**slip, **control}
