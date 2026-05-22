"""
Module 1: ContactNet - Contact Shape & Force Distribution Estimation
Estimates 3D contact geometry (depth map, surface normals) and 6-axis
wrench from GelSight/DIGIT tactile RGB images.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision import models


class PhotometricStereoLayer(nn.Module):
    """Differentiable photometric stereo for surface normal estimation."""

    def __init__(self, num_leds: int = 3, image_size: tuple = (64, 64)):
        super().__init__()
        self.num_leds = num_leds
        # Learnable light source directions (initialized to 120° apart)
        angles = torch.linspace(0, 2 * torch.pi, num_leds + 1)[:-1]
        light_dirs = torch.stack([
            torch.cos(angles),
            torch.sin(angles),
            torch.ones(num_leds) * 0.5,
        ], dim=1)
        light_dirs = F.normalize(light_dirs, dim=1)
        self.light_directions = nn.Parameter(light_dirs, requires_grad=True)

    def forward(self, tactile_rgb: torch.Tensor) -> torch.Tensor:
        """
        Args:
            tactile_rgb: (B, 3, H, W) RGB tactile image
        Returns:
            normals: (B, 3, H, W) surface normal map
        """
        B, C, H, W = tactile_rgb.shape
        # Reshape for least-squares: (B, 3, H*W)
        intensities = tactile_rgb.view(B, C, H * W).float()
        # Solve L^T * n = I for normals via pseudo-inverse
        L = self.light_directions  # (3, 3)
        L_pinv = torch.linalg.pinv(L.T)  # (3, 3)
        normals = torch.matmul(L_pinv, intensities)  # (B, 3, H*W)
        normals = F.normalize(normals, dim=1)
        return normals.view(B, 3, H, W)


class DepthReconstructor(nn.Module):
    """Poisson-based depth reconstruction from surface normals using DCT."""

    def __init__(self, height: int = 64, width: int = 64):
        super().__init__()
        self.height = height
        self.width = width

    def forward(self, normals: torch.Tensor) -> torch.Tensor:
        """
        Args:
            normals: (B, 3, H, W) surface normals
        Returns:
            depth: (B, 1, H, W) reconstructed depth map
        """
        nx = normals[:, 0:1]
        ny = normals[:, 1:2]
        nz = normals[:, 2:3].clamp(min=1e-6)

        # Gradient fields
        p = -nx / nz
        q = -ny / nz

        # Divergence
        dp_dx = F.pad(p[:, :, :, 1:] - p[:, :, :, :-1], (0, 1, 0, 0))
        dq_dy = F.pad(q[:, :, 1:, :] - q[:, :, :-1, :], (0, 0, 0, 1))
        divergence = dp_dx + dq_dy

        # Poisson solver via DCT (approximated with FFT for differentiability)
        div_fft = torch.fft.rfft2(divergence.squeeze(1))
        H, W = divergence.shape[2], divergence.shape[3]
        u = torch.arange(H, device=normals.device).float()
        v = torch.arange(W // 2 + 1, device=normals.device).float()
        u, v = torch.meshgrid(u, v, indexing='ij')
        denom = (2 * torch.cos(torch.pi * u / H) +
                 2 * torch.cos(torch.pi * v / W) - 4)
        denom[0, 0] = 1.0  # Avoid division by zero

        depth_fft = div_fft / denom.unsqueeze(0)
        depth_fft[:, 0, 0] = 0  # Zero-mean
        depth = torch.fft.irfft2(depth_fft, s=(H, W))
        return depth.unsqueeze(1)


class ContactNet(nn.Module):
    """
    End-to-end contact estimation network.
    Input:  RGB tactile image (B, 3, H, W)
    Output: depth map (B, 1, 64, 64), normals (B, 3, 64, 64), wrench (B, 6)
    """

    def __init__(self, backbone: str = "resnet18", depth_size: int = 64,
                 force_components: int = 6):
        super().__init__()
        self.depth_size = depth_size

        # Photometric stereo branch (physics-informed)
        self.ps_layer = PhotometricStereoLayer(num_leds=3,
                                                image_size=(depth_size, depth_size))
        self.depth_reconstructor = DepthReconstructor(depth_size, depth_size)

        # Learned encoder
        encoder = getattr(models, backbone)(pretrained=False)
        self.encoder = nn.Sequential(*list(encoder.children())[:-2])  # Remove FC
        enc_channels = 512 if backbone == "resnet18" else 2048

        # Depth refinement decoder
        self.depth_decoder = nn.Sequential(
            nn.ConvTranspose2d(enc_channels, 256, 4, 2, 1),
            nn.BatchNorm2d(256),
            nn.ReLU(inplace=True),
            nn.ConvTranspose2d(256, 128, 4, 2, 1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
            nn.ConvTranspose2d(128, 64, 4, 2, 1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.Conv2d(64, 1, 3, 1, 1),
        )

        # Normal refinement decoder
        self.normal_decoder = nn.Sequential(
            nn.ConvTranspose2d(enc_channels, 256, 4, 2, 1),
            nn.BatchNorm2d(256),
            nn.ReLU(inplace=True),
            nn.ConvTranspose2d(256, 128, 4, 2, 1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
            nn.ConvTranspose2d(128, 64, 4, 2, 1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.Conv2d(64, 3, 3, 1, 1),
        )

        # Force/torque estimation head
        self.force_head = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(),
            nn.Linear(enc_channels, 256),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3),
            nn.Linear(256, force_components),
        )

        # Fusion gate for physics + learned depth
        self.fusion_gate = nn.Sequential(
            nn.Conv2d(2, 16, 3, 1, 1),
            nn.ReLU(inplace=True),
            nn.Conv2d(16, 1, 1),
            nn.Sigmoid(),
        )

    def forward(self, tactile_img: torch.Tensor):
        """
        Args:
            tactile_img: (B, 3, H, W) raw tactile RGB image
        Returns:
            dict with keys: depth, normals, wrench
        """
        B = tactile_img.shape[0]

        # Physics-informed branch
        img_resized = F.interpolate(tactile_img, (self.depth_size, self.depth_size),
                                     mode='bilinear', align_corners=False)
        ps_normals = self.ps_layer(img_resized)
        ps_depth = self.depth_reconstructor(ps_normals)

        # Learned branch
        features = self.encoder(tactile_img)
        learned_depth = self.depth_decoder(features)
        learned_depth = F.interpolate(learned_depth, (self.depth_size, self.depth_size),
                                       mode='bilinear', align_corners=False)
        learned_normals = self.normal_decoder(features)
        learned_normals = F.interpolate(learned_normals, (self.depth_size, self.depth_size),
                                         mode='bilinear', align_corners=False)
        learned_normals = F.normalize(learned_normals, dim=1)

        # Gated fusion of physics + learned depth
        gate = self.fusion_gate(torch.cat([ps_depth, learned_depth], dim=1))
        fused_depth = gate * ps_depth + (1 - gate) * learned_depth

        # Fuse normals (weighted average + renormalize)
        fused_normals = 0.4 * ps_normals + 0.6 * learned_normals
        fused_normals = F.normalize(fused_normals, dim=1)

        # Force estimation
        wrench = self.force_head(features)

        return {
            "depth": fused_depth,
            "normals": fused_normals,
            "wrench": wrench,
            "ps_depth": ps_depth,
            "ps_normals": ps_normals,
        }


class ContactLoss(nn.Module):
    """Combined loss for contact estimation training."""

    def __init__(self, depth_weight: float = 1.0, normal_weight: float = 0.5,
                 force_weight: float = 2.0):
        super().__init__()
        self.w_d = depth_weight
        self.w_n = normal_weight
        self.w_f = force_weight

    def forward(self, pred: dict, target: dict) -> dict:
        # Depth loss: L1 + gradient matching
        depth_l1 = F.l1_loss(pred["depth"], target["depth"])
        pred_grad_x = pred["depth"][:, :, :, 1:] - pred["depth"][:, :, :, :-1]
        tgt_grad_x = target["depth"][:, :, :, 1:] - target["depth"][:, :, :, :-1]
        pred_grad_y = pred["depth"][:, :, 1:, :] - pred["depth"][:, :, :-1, :]
        tgt_grad_y = target["depth"][:, :, 1:, :] - target["depth"][:, :, :-1, :]
        grad_loss = F.l1_loss(pred_grad_x, tgt_grad_x) + F.l1_loss(pred_grad_y, tgt_grad_y)
        depth_loss = depth_l1 + 0.5 * grad_loss

        # Normal loss: cosine similarity
        cos_sim = F.cosine_similarity(pred["normals"], target["normals"], dim=1)
        normal_loss = (1.0 - cos_sim).mean()

        # Force loss: MSE + direction alignment
        force_mse = F.mse_loss(pred["wrench"], target["wrench"])
        force_dir = 1.0 - F.cosine_similarity(
            pred["wrench"][:, :3], target["wrench"][:, :3], dim=1
        ).mean()
        force_loss = force_mse + 0.3 * force_dir

        total = self.w_d * depth_loss + self.w_n * normal_loss + self.w_f * force_loss

        return {
            "total": total,
            "depth_loss": depth_loss,
            "normal_loss": normal_loss,
            "force_loss": force_loss,
        }
