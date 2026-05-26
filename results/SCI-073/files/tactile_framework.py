"""
Tactile Sensing Framework for Object Recognition and Manipulation
=================================================================
PyTorch/IsaacSim-based simulation and learning framework for:
1. Contact shape & force distribution estimation from tactile images
2. Texture classification via deep learning
3. Tactile-visual multimodal fusion
4. Real-time grasp stability evaluation
5. Slip detection and force control feedback
6. Exploratory grasping strategy for unknown objects
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import math
from typing import Tuple, Dict, List, Optional


# ============================================================
# 1. Tactile Image Simulator (GelSight/DIGIT-style)
# ============================================================

class TactileSimulator:
    """Simulates GelSight/DIGIT tactile sensor images.
    Generates synthetic tactile images with contact geometry,
    force distribution, and photometric rendering."""

    def __init__(self, resolution: int = 64, gel_thickness: float = 2.0,
                 elasticity: float = 0.5, num_leds: int = 3):
        self.resolution = resolution
        self.gel_thickness = gel_thickness
        self.elasticity = elasticity
        self.num_leds = num_leds
        # LED directions (RGB illumination from 3 angles)
        self.led_dirs = np.array([
            [-1, 0, 0.5],   # Red from left
            [1, 0, 0.5],    # Green from right
            [0, -1, 0.5],   # Blue from bottom
        ], dtype=np.float32)
        for i in range(self.num_leds):
            self.led_dirs[i] /= np.linalg.norm(self.led_dirs[i])

    def generate_contact_geometry(self, shape: str = 'sphere',
                                  params: Optional[Dict] = None) -> np.ndarray:
        """Generate depth map for various contact shapes."""
        x = np.linspace(-1, 1, self.resolution)
        y = np.linspace(-1, 1, self.resolution)
        X, Y = np.meshgrid(x, y)

        if params is None:
            params = {}
        cx, cy = params.get('center', (0.0, 0.0))
        radius = params.get('radius', 0.5)
        force = params.get('force', 1.0)

        if shape == 'sphere':
            R2 = (X - cx)**2 + (Y - cy)**2
            depth = np.maximum(0, radius**2 - R2)
            depth = np.sqrt(depth) * force * self.elasticity
        elif shape == 'cylinder':
            angle = params.get('angle', 0.0)
            Xr = (X - cx) * np.cos(angle) + (Y - cy) * np.sin(angle)
            depth = np.maximum(0, radius**2 - Xr**2)
            depth = np.sqrt(depth) * force * self.elasticity
        elif shape == 'edge':
            angle = params.get('angle', 0.0)
            Xr = (X - cx) * np.cos(angle) + (Y - cy) * np.sin(angle)
            depth = np.maximum(0, -Xr) * force * self.elasticity * radius
        elif shape == 'flat':
            R2 = (X - cx)**2 + (Y - cy)**2
            mask = R2 < radius**2
            depth = mask.astype(float) * force * self.elasticity * 0.3
        else:
            depth = np.zeros((self.resolution, self.resolution))

        return depth.astype(np.float32)

    def depth_to_normal(self, depth: np.ndarray) -> np.ndarray:
        """Convert depth map to surface normal map."""
        dy, dx = np.gradient(depth)
        normal = np.stack([-dx, -dy, np.ones_like(depth)], axis=-1)
        norm = np.linalg.norm(normal, axis=-1, keepdims=True) + 1e-8
        return (normal / norm).astype(np.float32)

    def render_tactile_image(self, depth: np.ndarray) -> np.ndarray:
        """Photometric rendering under multi-directional LED illumination."""
        normal = self.depth_to_normal(depth)
        image = np.zeros((self.resolution, self.resolution, 3), dtype=np.float32)
        for i, led_dir in enumerate(self.led_dirs):
            intensity = np.maximum(0, np.sum(normal * led_dir, axis=-1))
            image[:, :, i] = intensity
        # Add ambient and noise
        image = image * 0.7 + 0.15
        image += np.random.normal(0, 0.02, image.shape)
        return np.clip(image, 0, 1).astype(np.float32)

    def compute_force_distribution(self, depth: np.ndarray,
                                   youngs_modulus: float = 1e5) -> np.ndarray:
        """Estimate normal force distribution from depth using linear elasticity."""
        pixel_area = (2.0 / self.resolution) ** 2
        force = youngs_modulus * depth * pixel_area / self.gel_thickness
        return force.astype(np.float32)

    def generate_texture(self, texture_type: str = 'smooth',
                         frequency: float = 10.0) -> np.ndarray:
        """Generate texture patterns on contact surface."""
        x = np.linspace(0, 2*np.pi*frequency, self.resolution)
        y = np.linspace(0, 2*np.pi*frequency, self.resolution)
        X, Y = np.meshgrid(x, y)

        if texture_type == 'smooth':
            tex = np.zeros((self.resolution, self.resolution))
        elif texture_type == 'rough':
            tex = np.random.normal(0, 0.05, (self.resolution, self.resolution))
        elif texture_type == 'striped':
            tex = 0.03 * np.sin(X)
        elif texture_type == 'dotted':
            tex = 0.03 * np.sin(X) * np.sin(Y)
        elif texture_type == 'crosshatch':
            tex = 0.02 * (np.sin(X) + np.sin(Y))
        elif texture_type == 'wavy':
            tex = 0.03 * np.sin(X + 0.5 * np.sin(Y))
        elif texture_type == 'grid':
            tex = 0.03 * np.maximum(np.abs(np.sin(X)), np.abs(np.sin(Y)))
        elif texture_type == 'random_bumps':
            tex = np.zeros((self.resolution, self.resolution))
            n_bumps = np.random.randint(5, 15)
            for _ in range(n_bumps):
                bx, by = np.random.uniform(0, self.resolution, 2).astype(int)
                br = np.random.uniform(2, 6)
                xx, yy = np.meshgrid(range(self.resolution), range(self.resolution))
                bump = np.exp(-((xx-bx)**2 + (yy-by)**2) / (2*br**2))
                tex += 0.03 * bump
        else:
            tex = np.zeros((self.resolution, self.resolution))

        return tex.astype(np.float32)


# ============================================================
# 2. Synthetic Dataset
# ============================================================

class TactileDataset(Dataset):
    """Synthetic dataset for tactile sensing experiments."""

    SHAPES = ['sphere', 'cylinder', 'edge', 'flat']
    TEXTURES = ['smooth', 'rough', 'striped', 'dotted', 'crosshatch',
                'wavy', 'grid', 'random_bumps']

    def __init__(self, num_samples: int = 2000, resolution: int = 64,
                 include_visual: bool = True):
        self.simulator = TactileSimulator(resolution=resolution)
        self.num_samples = num_samples
        self.resolution = resolution
        self.include_visual = include_visual
        self.data = self._generate_data()

    def _generate_data(self) -> List[Dict]:
        data = []
        for i in range(self.num_samples):
            shape_idx = i % len(self.SHAPES)
            texture_idx = i % len(self.TEXTURES)
            shape = self.SHAPES[shape_idx]
            texture_type = self.TEXTURES[texture_idx]
            force = np.random.uniform(0.3, 2.0)
            radius = np.random.uniform(0.2, 0.6)
            cx = np.random.uniform(-0.3, 0.3)
            cy = np.random.uniform(-0.3, 0.3)
            angle = np.random.uniform(0, np.pi)

            params = {'center': (cx, cy), 'radius': radius,
                      'force': force, 'angle': angle}

            depth = self.simulator.generate_contact_geometry(shape, params)
            texture = self.simulator.generate_texture(texture_type,
                                                       frequency=np.random.uniform(5, 15))
            depth_with_tex = depth + texture * (depth > 0.01).astype(float)
            tactile_img = self.simulator.render_tactile_image(depth_with_tex)
            force_dist = self.simulator.compute_force_distribution(depth)

            # Slip label: high force + small contact = likely slip
            contact_area = np.sum(depth > 0.01)
            total_force = np.sum(force_dist)
            slip_risk = total_force / (contact_area + 1e-6)
            is_slipping = 1 if slip_risk > 0.15 else 0

            # Grasp stability: function of contact area and force uniformity
            force_std = np.std(force_dist[depth > 0.01]) if contact_area > 0 else 1.0
            stability = np.clip(1.0 - force_std / (np.mean(force_dist[depth > 0.01]) + 1e-6), 0, 1) \
                if contact_area > 10 else 0.0

            sample = {
                'tactile_image': tactile_img.transpose(2, 0, 1),  # CHW
                'depth_map': depth[np.newaxis],
                'force_distribution': force_dist[np.newaxis],
                'shape_label': shape_idx,
                'texture_label': texture_idx,
                'slip_label': is_slipping,
                'stability_score': np.float32(stability),
                'force_magnitude': np.float32(total_force),
                'contact_area': np.float32(contact_area / (self.resolution**2)),
            }

            if self.include_visual:
                visual = self._generate_visual_image(shape, texture_type, params)
                sample['visual_image'] = visual.transpose(2, 0, 1)

            data.append(sample)
        return data

    def _generate_visual_image(self, shape: str, texture: str,
                                params: Dict) -> np.ndarray:
        """Generate a synthetic RGB visual image of the object."""
        img = np.ones((self.resolution, self.resolution, 3), dtype=np.float32) * 0.3
        x = np.linspace(-1, 1, self.resolution)
        y = np.linspace(-1, 1, self.resolution)
        X, Y = np.meshgrid(x, y)
        cx, cy = params['center']
        r = params['radius']
        mask = ((X-cx)**2 + (Y-cy)**2) < r**2

        color_map = {
            'sphere': [0.8, 0.3, 0.2],
            'cylinder': [0.2, 0.7, 0.3],
            'edge': [0.3, 0.3, 0.8],
            'flat': [0.7, 0.7, 0.2],
        }
        color = color_map.get(shape, [0.5, 0.5, 0.5])
        for c in range(3):
            img[:, :, c] = np.where(mask, color[c], img[:, :, c])

        img += np.random.normal(0, 0.03, img.shape)
        return np.clip(img, 0, 1).astype(np.float32)

    def __len__(self):
        return self.num_samples

    def __getitem__(self, idx):
        sample = self.data[idx]
        return {k: torch.tensor(v) if isinstance(v, np.ndarray)
                else torch.tensor(v) for k, v in sample.items()}


# ============================================================
# 3. Neural Network Models
# ============================================================

class ResidualBlock(nn.Module):
    def __init__(self, channels: int):
        super().__init__()
        self.conv1 = nn.Conv2d(channels, channels, 3, padding=1)
        self.bn1 = nn.BatchNorm2d(channels)
        self.conv2 = nn.Conv2d(channels, channels, 3, padding=1)
        self.bn2 = nn.BatchNorm2d(channels)

    def forward(self, x):
        residual = x
        out = F.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        return F.relu(out + residual)


class ContactEstimationNet(nn.Module):
    """Task 1: Estimate contact shape (depth) and force distribution
    from tactile images using an encoder-decoder architecture."""

    def __init__(self, in_channels: int = 3):
        super().__init__()
        # Encoder
        self.enc1 = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(),
            ResidualBlock(32))
        self.enc2 = nn.Sequential(
            nn.Conv2d(32, 64, 3, stride=2, padding=1), nn.BatchNorm2d(64), nn.ReLU(),
            ResidualBlock(64))
        self.enc3 = nn.Sequential(
            nn.Conv2d(64, 128, 3, stride=2, padding=1), nn.BatchNorm2d(128), nn.ReLU(),
            ResidualBlock(128))

        # Decoder for depth
        self.dec_depth3 = nn.Sequential(
            nn.ConvTranspose2d(128, 64, 4, stride=2, padding=1), nn.BatchNorm2d(64), nn.ReLU())
        self.dec_depth2 = nn.Sequential(
            nn.ConvTranspose2d(128, 32, 4, stride=2, padding=1), nn.BatchNorm2d(32), nn.ReLU())
        self.dec_depth1 = nn.Conv2d(64, 1, 3, padding=1)

        # Decoder for force
        self.dec_force3 = nn.Sequential(
            nn.ConvTranspose2d(128, 64, 4, stride=2, padding=1), nn.BatchNorm2d(64), nn.ReLU())
        self.dec_force2 = nn.Sequential(
            nn.ConvTranspose2d(128, 32, 4, stride=2, padding=1), nn.BatchNorm2d(32), nn.ReLU())
        self.dec_force1 = nn.Conv2d(64, 1, 3, padding=1)

    def forward(self, tactile_img: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        e1 = self.enc1(tactile_img)
        e2 = self.enc2(e1)
        e3 = self.enc3(e2)

        # Depth decoder with skip connections
        d3 = self.dec_depth3(e3)
        d2 = self.dec_depth2(torch.cat([d3, e2], dim=1))
        depth = F.relu(self.dec_depth1(torch.cat([d2, e1], dim=1)))

        # Force decoder with skip connections
        f3 = self.dec_force3(e3)
        f2 = self.dec_force2(torch.cat([f3, e2], dim=1))
        force = F.relu(self.dec_force1(torch.cat([f2, e1], dim=1)))

        return depth, force


class TextureClassifier(nn.Module):
    """Task 2: Classify texture from tactile images using CNN + attention."""

    def __init__(self, num_classes: int = 8, in_channels: int = 3):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(64, 128, 3, padding=1), nn.BatchNorm2d(128), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(128, 256, 3, padding=1), nn.BatchNorm2d(256), nn.ReLU(),
        )
        # Channel attention (SE block)
        self.se = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(),
            nn.Linear(256, 64), nn.ReLU(),
            nn.Linear(64, 256), nn.Sigmoid()
        )
        self.classifier = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(),
            nn.Linear(256, 128), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(128, num_classes)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        feat = self.features(x)
        # Apply SE attention
        att = self.se(feat).unsqueeze(-1).unsqueeze(-1)
        feat = feat * att
        return self.classifier(feat)


class MultimodalFusionNet(nn.Module):
    """Task 3: Fuse tactile and visual modalities for enhanced recognition."""

    def __init__(self, num_shape_classes: int = 4, num_texture_classes: int = 8):
        super().__init__()
        # Tactile encoder
        self.tactile_enc = nn.Sequential(
            nn.Conv2d(3, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(64, 128, 3, padding=1), nn.BatchNorm2d(128), nn.ReLU(),
            nn.AdaptiveAvgPool2d(4),
        )
        # Visual encoder
        self.visual_enc = nn.Sequential(
            nn.Conv2d(3, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(64, 128, 3, padding=1), nn.BatchNorm2d(128), nn.ReLU(),
            nn.AdaptiveAvgPool2d(4),
        )
        # Cross-attention fusion
        self.cross_attention = nn.MultiheadAttention(embed_dim=128, num_heads=4, batch_first=True)
        # Classifiers
        self.shape_head = nn.Sequential(
            nn.Linear(256, 64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, num_shape_classes))
        self.texture_head = nn.Sequential(
            nn.Linear(256, 64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, num_texture_classes))

    def forward(self, tactile: torch.Tensor, visual: torch.Tensor
                ) -> Tuple[torch.Tensor, torch.Tensor]:
        t_feat = self.tactile_enc(tactile)  # B, 128, 4, 4
        v_feat = self.visual_enc(visual)
        B = t_feat.size(0)
        t_seq = t_feat.view(B, 128, -1).permute(0, 2, 1)  # B, 16, 128
        v_seq = v_feat.view(B, 128, -1).permute(0, 2, 1)
        # Cross-attention: tactile queries, visual keys/values
        fused, _ = self.cross_attention(t_seq, v_seq, v_seq)
        # Combine
        t_pool = t_seq.mean(dim=1)
        f_pool = fused.mean(dim=1)
        combined = torch.cat([t_pool, f_pool], dim=1)  # B, 256
        return self.shape_head(combined), self.texture_head(combined)


class GraspStabilityNet(nn.Module):
    """Task 4: Predict grasp stability score from tactile images."""

    def __init__(self, in_channels: int = 3):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(64, 128, 3, padding=1), nn.ReLU(),
            nn.AdaptiveAvgPool2d(4),
        )
        self.temporal = nn.LSTM(128 * 16, 64, batch_first=True, num_layers=2)
        self.regressor = nn.Sequential(
            nn.Linear(64, 32), nn.ReLU(),
            nn.Linear(32, 1), nn.Sigmoid()
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, T, C, H, W) or (B, C, H, W)
        if x.dim() == 4:
            x = x.unsqueeze(1)  # add time dim
        B, T, C, H, W = x.shape
        feat = self.features(x.view(B*T, C, H, W))
        feat = feat.view(B, T, -1)
        _, (h_n, _) = self.temporal(feat)
        return self.regressor(h_n[-1]).squeeze(-1)


class SlipDetectionNet(nn.Module):
    """Task 5: Detect slip events from tactile image sequences."""

    def __init__(self, in_channels: int = 3):
        super().__init__()
        self.spatial = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(64, 128, 3, padding=1), nn.ReLU(),
            nn.AdaptiveAvgPool2d(4),
        )
        # Temporal difference module
        self.diff_conv = nn.Sequential(
            nn.Conv2d(128, 64, 1), nn.ReLU(),
            nn.AdaptiveAvgPool2d(1), nn.Flatten()
        )
        self.classifier = nn.Sequential(
            nn.Linear(64, 32), nn.ReLU(),
            nn.Linear(32, 2)  # slip / no-slip
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # For single frame, use spatial features directly
        feat = self.spatial(x)
        diff_feat = self.diff_conv(feat)
        return self.classifier(diff_feat)


class ExploratoryGraspPolicy(nn.Module):
    """Task 6: Policy network for exploratory grasping of unknown objects.
    Uses tactile feedback to determine next grasp action."""

    def __init__(self, tactile_dim: int = 3, action_dim: int = 6):
        super().__init__()
        self.tactile_encoder = nn.Sequential(
            nn.Conv2d(tactile_dim, 32, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(64, 128, 3, padding=1), nn.ReLU(),
            nn.AdaptiveAvgPool2d(4),
            nn.Flatten()
        )
        # State includes tactile features + previous action + force info
        state_dim = 128 * 16 + action_dim + 2
        self.policy = nn.Sequential(
            nn.Linear(state_dim, 256), nn.ReLU(),
            nn.Linear(256, 128), nn.ReLU(),
        )
        self.action_mean = nn.Linear(128, action_dim)
        self.action_log_std = nn.Parameter(torch.zeros(action_dim))
        self.value_head = nn.Linear(128, 1)

    def forward(self, tactile: torch.Tensor, prev_action: torch.Tensor,
                force_info: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        t_feat = self.tactile_encoder(tactile)
        state = torch.cat([t_feat, prev_action, force_info], dim=1)
        hidden = self.policy(state)
        action_mean = self.action_mean(hidden)
        action_std = self.action_log_std.exp().expand_as(action_mean)
        value = self.value_head(hidden)
        return action_mean, action_std, value


# ============================================================
# 4. Force Control Module
# ============================================================

class ForceController:
    """PID-based force controller with slip compensation."""

    def __init__(self, kp: float = 1.0, ki: float = 0.1, kd: float = 0.05,
                 force_target: float = 5.0, slip_gain: float = 2.0):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.force_target = force_target
        self.slip_gain = slip_gain
        self.integral = 0.0
        self.prev_error = 0.0
        self.dt = 0.01

    def step(self, current_force: float, slip_detected: bool) -> float:
        target = self.force_target
        if slip_detected:
            target *= (1 + self.slip_gain)
        error = target - current_force
        self.integral += error * self.dt
        derivative = (error - self.prev_error) / self.dt
        self.prev_error = error
        control = self.kp * error + self.ki * self.integral + self.kd * derivative
        return np.clip(control, -10, 10)


# ============================================================
# 5. IsaacSim-style Environment Interface
# ============================================================

class TactileGraspEnv:
    """Simulated grasping environment (IsaacSim-style interface)."""

    def __init__(self, resolution: int = 64):
        self.simulator = TactileSimulator(resolution=resolution)
        self.resolution = resolution
        self.controller = ForceController()
        self.current_force = 0.0
        self.gripper_pos = np.zeros(6)
        self.object_shapes = TactileDataset.SHAPES
        self.object_textures = TactileDataset.TEXTURES
        self.reset()

    def reset(self) -> Dict[str, np.ndarray]:
        self.current_force = 0.0
        self.gripper_pos = np.zeros(6)
        self.object_shape = np.random.choice(self.object_shapes)
        self.object_texture = np.random.choice(self.object_textures)
        self.contact_params = {
            'center': (np.random.uniform(-0.2, 0.2), np.random.uniform(-0.2, 0.2)),
            'radius': np.random.uniform(0.2, 0.5),
            'force': 0.5,
            'angle': np.random.uniform(0, np.pi)
        }
        return self._get_obs()

    def step(self, action: np.ndarray) -> Tuple[Dict, float, bool, Dict]:
        self.gripper_pos += action[:6] * 0.1
        force_delta = action[0] if len(action) > 0 else 0
        self.contact_params['force'] = np.clip(
            self.contact_params['force'] + force_delta * 0.1, 0, 3.0)

        obs = self._get_obs()
        depth = self.simulator.generate_contact_geometry(
            self.object_shape, self.contact_params)
        force_dist = self.simulator.compute_force_distribution(depth)
        self.current_force = float(np.sum(force_dist))

        contact_area = np.sum(depth > 0.01) / (self.resolution ** 2)
        stability = min(1.0, contact_area * 5) * min(1.0, self.current_force * 0.5)
        reward = stability - 0.1 * np.abs(force_delta)

        done = self.current_force > 50 or stability > 0.95
        info = {'stability': stability, 'force': self.current_force,
                'contact_area': contact_area}
        return obs, reward, done, info

    def _get_obs(self) -> Dict[str, np.ndarray]:
        depth = self.simulator.generate_contact_geometry(
            self.object_shape, self.contact_params)
        texture = self.simulator.generate_texture(self.object_texture)
        depth_with_tex = depth + texture * (depth > 0.01).astype(float)
        tactile_img = self.simulator.render_tactile_image(depth_with_tex)
        return {
            'tactile_image': tactile_img,
            'force': np.array([self.current_force]),
            'gripper_pos': self.gripper_pos.copy()
        }


# ============================================================
# 6. Training & Evaluation Functions
# ============================================================

def train_contact_estimation(model, dataloader, epochs=20, lr=1e-3, device='cpu'):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    history = {'depth_loss': [], 'force_loss': [], 'total_loss': []}

    for epoch in range(epochs):
        model.train()
        epoch_losses = {'depth': 0, 'force': 0, 'total': 0}
        for batch in dataloader:
            tactile = batch['tactile_image'].to(device)
            gt_depth = batch['depth_map'].to(device)
            gt_force = batch['force_distribution'].to(device)

            pred_depth, pred_force = model(tactile)
            depth_loss = F.mse_loss(pred_depth, gt_depth)
            force_loss = F.mse_loss(pred_force, gt_force)
            loss = depth_loss + 0.5 * force_loss

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_losses['depth'] += depth_loss.item()
            epoch_losses['force'] += force_loss.item()
            epoch_losses['total'] += loss.item()

        n = len(dataloader)
        history['depth_loss'].append(epoch_losses['depth'] / n)
        history['force_loss'].append(epoch_losses['force'] / n)
        history['total_loss'].append(epoch_losses['total'] / n)

        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}/{epochs}: depth_loss={history['depth_loss'][-1]:.6f}, "
                  f"force_loss={history['force_loss'][-1]:.6f}")

    return history


def train_texture_classifier(model, dataloader, epochs=20, lr=1e-3, device='cpu'):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()
    history = {'loss': [], 'accuracy': []}

    for epoch in range(epochs):
        model.train()
        total_loss, correct, total = 0, 0, 0
        for batch in dataloader:
            tactile = batch['tactile_image'].to(device)
            labels = batch['texture_label'].to(device)

            logits = model(tactile)
            loss = criterion(logits, labels)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            correct += (logits.argmax(1) == labels).sum().item()
            total += labels.size(0)

        history['loss'].append(total_loss / len(dataloader))
        history['accuracy'].append(correct / total)

        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}/{epochs}: loss={history['loss'][-1]:.4f}, "
                  f"acc={history['accuracy'][-1]:.4f}")

    return history


def train_multimodal(model, dataloader, epochs=20, lr=1e-3, device='cpu'):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()
    history = {'loss': [], 'shape_acc': [], 'texture_acc': []}

    for epoch in range(epochs):
        model.train()
        total_loss = 0
        shape_correct, tex_correct, total = 0, 0, 0
        for batch in dataloader:
            tactile = batch['tactile_image'].to(device)
            visual = batch['visual_image'].to(device)
            shape_labels = batch['shape_label'].to(device)
            tex_labels = batch['texture_label'].to(device)

            shape_logits, tex_logits = model(tactile, visual)
            loss = criterion(shape_logits, shape_labels) + criterion(tex_logits, tex_labels)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            shape_correct += (shape_logits.argmax(1) == shape_labels).sum().item()
            tex_correct += (tex_logits.argmax(1) == tex_labels).sum().item()
            total += shape_labels.size(0)

        history['loss'].append(total_loss / len(dataloader))
        history['shape_acc'].append(shape_correct / total)
        history['texture_acc'].append(tex_correct / total)

        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}/{epochs}: loss={history['loss'][-1]:.4f}, "
                  f"shape_acc={history['shape_acc'][-1]:.4f}, "
                  f"tex_acc={history['texture_acc'][-1]:.4f}")

    return history


def train_slip_detector(model, dataloader, epochs=20, lr=1e-3, device='cpu'):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()
    history = {'loss': [], 'accuracy': [], 'precision': [], 'recall': []}

    for epoch in range(epochs):
        model.train()
        total_loss = 0
        all_preds, all_labels = [], []
        for batch in dataloader:
            tactile = batch['tactile_image'].to(device)
            labels = batch['slip_label'].to(device)

            logits = model(tactile)
            loss = criterion(logits, labels)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            preds = logits.argmax(1)
            all_preds.extend(preds.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

        all_preds = np.array(all_preds)
        all_labels = np.array(all_labels)
        acc = (all_preds == all_labels).mean()
        tp = ((all_preds == 1) & (all_labels == 1)).sum()
        fp = ((all_preds == 1) & (all_labels == 0)).sum()
        fn = ((all_preds == 0) & (all_labels == 1)).sum()
        precision = tp / (tp + fp + 1e-8)
        recall = tp / (tp + fn + 1e-8)

        history['loss'].append(total_loss / len(dataloader))
        history['accuracy'].append(acc)
        history['precision'].append(precision)
        history['recall'].append(recall)

        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}/{epochs}: loss={history['loss'][-1]:.4f}, "
                  f"acc={acc:.4f}, prec={precision:.4f}, rec={recall:.4f}")

    return history


def train_grasp_stability(model, dataloader, epochs=20, lr=1e-3, device='cpu'):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    history = {'loss': [], 'mae': []}

    for epoch in range(epochs):
        model.train()
        total_loss, total_mae = 0, 0
        n_batches = 0
        for batch in dataloader:
            tactile = batch['tactile_image'].to(device)
            gt_stability = batch['stability_score'].to(device)

            pred = model(tactile)
            loss = F.mse_loss(pred, gt_stability)
            mae = F.l1_loss(pred, gt_stability)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            total_mae += mae.item()
            n_batches += 1

        history['loss'].append(total_loss / n_batches)
        history['mae'].append(total_mae / n_batches)

        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}/{epochs}: loss={history['loss'][-1]:.6f}, "
                  f"mae={history['mae'][-1]:.4f}")

    return history


def run_force_control_simulation(n_steps: int = 200) -> Dict:
    """Run force control with slip compensation simulation."""
    controller = ForceController(kp=2.0, ki=0.3, kd=0.1, force_target=5.0)
    forces, targets, controls, slips = [], [], [], []
    current_force = 0.5
    np.random.seed(42)

    for step in range(n_steps):
        # Simulate slip events
        slip = step in range(50, 70) or step in range(120, 140) or step in range(170, 185)
        disturbance = np.random.normal(0, 0.2)
        if slip:
            disturbance -= 2.0

        control = controller.step(current_force, slip)
        current_force += control * 0.01 + disturbance * 0.1
        current_force = max(0, current_force)

        forces.append(current_force)
        targets.append(controller.force_target * (1 + controller.slip_gain) if slip
                       else controller.force_target)
        controls.append(control)
        slips.append(slip)

    return {'forces': forces, 'targets': targets, 'controls': controls, 'slips': slips}


def run_exploratory_grasp(env, policy_model, n_episodes: int = 50, device='cpu'):
    """Run exploratory grasping episodes."""
    policy_model.to(device)
    policy_model.eval()
    results = {'rewards': [], 'stabilities': [], 'forces': [], 'steps': []}

    for ep in range(n_episodes):
        obs = env.reset()
        episode_reward = 0
        prev_action = torch.zeros(1, 6).to(device)
        max_stability = 0

        for step in range(50):
            tactile = torch.tensor(obs['tactile_image'].transpose(2, 0, 1)).unsqueeze(0).to(device)
            force_info = torch.tensor([[obs['force'][0], step / 50.0]]).float().to(device)

            with torch.no_grad():
                action_mean, action_std, _ = policy_model(tactile, prev_action, force_info)
                action = action_mean + action_std * torch.randn_like(action_mean) * 0.1

            obs, reward, done, info = env.step(action.cpu().numpy().flatten())
            episode_reward += reward
            max_stability = max(max_stability, info['stability'])
            prev_action = action

            if done:
                break

        results['rewards'].append(episode_reward)
        results['stabilities'].append(max_stability)
        results['forces'].append(info['force'])
        results['steps'].append(step + 1)

    return results
