"""
Isaac Sim Tactile Simulation Environment
Provides tactile sensor simulation, object interaction, and
domain randomization for training the tactile manipulation system.
"""

import torch
import torch.nn.functional as F
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
import math


@dataclass
class TactileSensorConfig:
    """Configuration for simulated tactile sensor."""
    resolution: Tuple[int, int] = (320, 240)
    gel_mesh_resolution: Tuple[int, int] = (64, 64)
    gel_thickness_mm: float = 3.0
    elastomer_modulus_kpa: float = 100.0
    friction_coefficient: float = 0.8
    contact_stiffness: float = 1000.0
    led_colors: List[str] = field(default_factory=lambda: ["red", "green", "blue"])
    led_angles_deg: List[float] = field(default_factory=lambda: [0, 120, 240])
    noise_std: float = 0.01


@dataclass
class SimConfig:
    """Isaac Sim environment configuration."""
    gravity: Tuple[float, ...] = (0, 0, -9.81)
    timestep: float = 0.001
    render_interval: int = 4
    robot_type: str = "franka_panda"
    num_envs: int = 64
    object_dataset: str = "ycb"
    domain_randomization: bool = True


class GelDeformationModel:
    """
    Finite-element-inspired gel deformation simulator.
    Models the elastomer gel as a neo-Hookean material under contact.
    """

    def __init__(self, config: TactileSensorConfig):
        self.config = config
        self.mesh_h, self.mesh_w = config.gel_mesh_resolution
        self.E = config.elastomer_modulus_kpa * 1000  # Convert to Pa
        self.thickness = config.gel_thickness_mm / 1000  # Convert to m
        self.nu = 0.49  # Poisson's ratio (nearly incompressible)
        self._init_mesh()

    def _init_mesh(self):
        """Initialize rest-state gel mesh."""
        y = np.linspace(-0.01, 0.01, self.mesh_h)
        x = np.linspace(-0.01, 0.01, self.mesh_w)
        self.rest_x, self.rest_y = np.meshgrid(x, y)
        self.rest_z = np.zeros_like(self.rest_x)
        # Stiffness matrix (simplified)
        self.k = self.E * self.thickness / (1 - self.nu ** 2)

    def compute_deformation(self, contact_points: np.ndarray,
                             contact_forces: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Compute gel deformation from contact forces.

        Args:
            contact_points: (N, 3) contact point positions
            contact_forces: (N, 3) forces at each contact point

        Returns:
            dict with depth_map, normal_map, stress_map
        """
        depth_map = np.zeros((self.mesh_h, self.mesh_w))
        stress_map = np.zeros((self.mesh_h, self.mesh_w, 3))

        if len(contact_points) == 0:
            normal_map = np.zeros((self.mesh_h, self.mesh_w, 3))
            normal_map[:, :, 2] = 1.0
            return {
                "depth_map": depth_map,
                "normal_map": normal_map,
                "stress_map": stress_map,
            }

        for i, (pt, force) in enumerate(zip(contact_points, contact_forces)):
            # Boussinesq solution for point load on elastic half-space
            dx = self.rest_x - pt[0]
            dy = self.rest_y - pt[1]
            r = np.sqrt(dx ** 2 + dy ** 2 + 1e-10)

            fz = abs(force[2])
            # Vertical displacement
            uz = fz / (np.pi * self.k) * (1 / r)
            uz = np.clip(uz, 0, self.thickness)

            depth_map += uz

            # Stress components
            sigma_z = -fz / (2 * np.pi) * (3 * self.thickness ** 3) / (r ** 2 + self.thickness ** 2) ** 2.5
            stress_map[:, :, 2] += sigma_z

        # Compute normals from depth gradients
        grad_y, grad_x = np.gradient(depth_map)
        normal_map = np.stack([-grad_x, -grad_y, np.ones_like(grad_x)], axis=-1)
        norms = np.linalg.norm(normal_map, axis=-1, keepdims=True)
        normal_map = normal_map / (norms + 1e-8)

        return {
            "depth_map": depth_map,
            "normal_map": normal_map,
            "stress_map": stress_map,
        }


class TactileImageRenderer:
    """
    Render synthetic tactile images from gel deformation.
    Simulates the photometric stereo imaging of GelSight/DIGIT sensors.
    """

    def __init__(self, config: TactileSensorConfig):
        self.config = config
        self.resolution = config.resolution
        self._init_lighting()

    def _init_lighting(self):
        """Initialize LED lighting model."""
        self.light_dirs = []
        for angle in self.config.led_angles_deg:
            rad = math.radians(angle)
            elevation = math.radians(30)
            ldir = np.array([
                math.cos(rad) * math.cos(elevation),
                math.sin(rad) * math.cos(elevation),
                math.sin(elevation),
            ])
            self.light_dirs.append(ldir / np.linalg.norm(ldir))

        # LED color intensities (R, G, B)
        self.light_colors = np.array([
            [1.0, 0.1, 0.1],  # Red LED
            [0.1, 1.0, 0.1],  # Green LED
            [0.1, 0.1, 1.0],  # Blue LED
        ])

    def render(self, deformation: Dict[str, np.ndarray],
               albedo: float = 0.8) -> np.ndarray:
        """
        Render tactile RGB image from gel deformation.

        Args:
            deformation: output from GelDeformationModel
            albedo: gel surface albedo

        Returns:
            image: (H, W, 3) RGB tactile image in [0, 1]
        """
        normal_map = deformation["normal_map"]
        H, W = self.resolution[1], self.resolution[0]

        # Resize normal map to sensor resolution
        from scipy.ndimage import zoom
        scale_h = H / normal_map.shape[0]
        scale_w = W / normal_map.shape[1]
        normal_hr = zoom(normal_map, (scale_h, scale_w, 1), order=1)
        normal_hr = normal_hr / (np.linalg.norm(normal_hr, axis=-1, keepdims=True) + 1e-8)

        # Phong shading with each LED
        image = np.zeros((H, W, 3))
        for led_idx, (ldir, lcolor) in enumerate(
            zip(self.light_dirs, self.light_colors)
        ):
            # Lambertian diffuse
            dot = np.sum(normal_hr * ldir, axis=-1)
            dot = np.clip(dot, 0, 1)
            for c in range(3):
                image[:, :, c] += albedo * lcolor[c] * dot

        # Add background illumination
        image += 0.05

        # Gamma correction
        image = np.clip(image, 0, 1) ** (1 / 2.2)

        # Add sensor noise
        noise = np.random.randn(*image.shape) * self.config.noise_std
        image = np.clip(image + noise, 0, 1)

        return image.astype(np.float32)


class DomainRandomizer:
    """Domain randomization for sim-to-real transfer."""

    def __init__(self, config: SimConfig):
        self.config = config
        self.rng = np.random.default_rng(42)

    def randomize_object(self) -> Dict:
        """Randomize object physical properties."""
        return {
            "mass": self.rng.uniform(0.01, 2.0),
            "friction": self.rng.uniform(0.3, 1.2),
            "restitution": self.rng.uniform(0.0, 0.5),
            "texture_id": self.rng.integers(0, 100),
        }

    def randomize_lighting(self) -> Dict:
        """Randomize lighting conditions."""
        return {
            "intensity_scale": self.rng.uniform(0.7, 1.3),
            "color_shift": self.rng.uniform(-0.05, 0.05, size=3),
            "ambient": self.rng.uniform(0.02, 0.1),
        }

    def randomize_sensor(self, base_config: TactileSensorConfig) -> TactileSensorConfig:
        """Randomize sensor parameters."""
        import copy
        config = copy.deepcopy(base_config)
        config.noise_std = self.rng.uniform(0.005, 0.02)
        config.elastomer_modulus_kpa *= self.rng.uniform(0.8, 1.2)
        config.friction_coefficient *= self.rng.uniform(0.9, 1.1)
        return config

    def augment_tactile_image(self, image: np.ndarray) -> np.ndarray:
        """Apply image-level augmentations."""
        # Brightness/contrast
        alpha = self.rng.uniform(0.8, 1.2)
        beta = self.rng.uniform(-0.05, 0.05)
        image = alpha * image + beta

        # Color jitter
        color_scale = self.rng.uniform(0.9, 1.1, size=3)
        image = image * color_scale

        # Gaussian blur (occasional)
        if self.rng.random() < 0.3:
            from scipy.ndimage import gaussian_filter
            sigma = self.rng.uniform(0.5, 1.5)
            image = gaussian_filter(image, sigma=(sigma, sigma, 0))

        return np.clip(image, 0, 1).astype(np.float32)


class IsaacSimTactileEnv:
    """
    Isaac Sim environment wrapper for tactile manipulation.

    NOTE: This is a simulation-ready design. Full Isaac Sim integration
    requires the NVIDIA Isaac Sim SDK. This class provides the interface
    and can run with synthetic data for development/testing.
    """

    def __init__(self, sim_config: SimConfig,
                 sensor_config: TactileSensorConfig):
        self.sim_config = sim_config
        self.sensor_config = sensor_config
        self.gel_model = GelDeformationModel(sensor_config)
        self.renderer = TactileImageRenderer(sensor_config)
        self.randomizer = DomainRandomizer(sim_config)
        self.num_envs = sim_config.num_envs
        self.step_count = 0
        self._obs_cache = {}

    def reset(self, env_ids: Optional[List[int]] = None) -> Dict[str, torch.Tensor]:
        """Reset specified environments."""
        if env_ids is None:
            env_ids = list(range(self.num_envs))

        obs = {
            "tactile_img": torch.zeros(
                len(env_ids), 3,
                self.sensor_config.resolution[1],
                self.sensor_config.resolution[0]
            ),
            "visual_img": torch.zeros(len(env_ids), 3, 224, 224),
            "wrench": torch.zeros(len(env_ids), 6),
            "joint_pos": torch.zeros(len(env_ids), 7),
            "joint_vel": torch.zeros(len(env_ids), 7),
            "gripper_pos": torch.zeros(len(env_ids), 1),
        }

        self.step_count = 0
        self._obs_cache = obs
        return obs

    def step(self, actions: torch.Tensor) -> Tuple[
        Dict[str, torch.Tensor], torch.Tensor, torch.Tensor, Dict
    ]:
        """
        Step the simulation.

        Args:
            actions: (num_envs, action_dim) gripper commands

        Returns:
            observations, rewards, dones, info
        """
        B = actions.shape[0]
        self.step_count += 1

        # Simulate contact (synthetic for development)
        contact_points = np.random.randn(3, 3) * 0.005
        contact_forces = np.random.randn(3, 3) * 0.5
        contact_forces[:, 2] = abs(contact_forces[:, 2])

        # Generate observations with domain randomization
        tactile_imgs = []
        for i in range(B):
            if self.sim_config.domain_randomization:
                sensor_cfg = self.randomizer.randomize_sensor(self.sensor_config)
                gel_model = GelDeformationModel(sensor_cfg)
            else:
                gel_model = self.gel_model

            deformation = gel_model.compute_deformation(contact_points, contact_forces)
            img = self.renderer.render(deformation)

            if self.sim_config.domain_randomization:
                img = self.randomizer.augment_tactile_image(img)

            tactile_imgs.append(torch.from_numpy(img).permute(2, 0, 1))

        obs = {
            "tactile_img": torch.stack(tactile_imgs),
            "visual_img": torch.randn(B, 3, 224, 224) * 0.1 + 0.5,
            "wrench": torch.randn(B, 6) * 0.5,
            "joint_pos": torch.randn(B, 7) * 0.1,
            "joint_vel": torch.randn(B, 7) * 0.01,
            "gripper_pos": torch.rand(B, 1) * 0.08,
        }

        # Reward: stability-based
        rewards = torch.randn(B) * 0.1 + 0.5
        dones = torch.zeros(B, dtype=torch.bool)

        if self.step_count >= 1000:
            dones[:] = True

        info = {
            "contact_area": torch.rand(B),
            "slip_detected": torch.zeros(B, dtype=torch.bool),
            "grasp_success": torch.ones(B, dtype=torch.bool),
        }

        self._obs_cache = obs
        return obs, rewards, dones, info

    def get_tactile_data(self) -> Dict[str, torch.Tensor]:
        """Get current tactile sensor data."""
        return {
            "tactile_img": self._obs_cache.get(
                "tactile_img",
                torch.zeros(self.num_envs, 3, 240, 320)
            ),
            "wrench": self._obs_cache.get(
                "wrench",
                torch.zeros(self.num_envs, 6)
            ),
        }


class TactileDataGenerator:
    """Generate synthetic training data using the simulation environment."""

    def __init__(self, env: IsaacSimTactileEnv, save_dir: str = "data/"):
        self.env = env
        self.save_dir = save_dir

    def generate_contact_dataset(self, num_samples: int = 10000) -> Dict:
        """Generate contact estimation training data."""
        import os
        os.makedirs(self.save_dir, exist_ok=True)

        tactile_images = []
        depth_maps = []
        normal_maps = []
        wrenches = []

        for i in range(num_samples):
            # Random contact
            n_contacts = np.random.randint(1, 6)
            points = np.random.randn(n_contacts, 3) * 0.003
            forces = np.random.randn(n_contacts, 3)
            forces[:, 2] = abs(forces[:, 2]) * 2

            deformation = self.env.gel_model.compute_deformation(points, forces)
            image = self.env.renderer.render(deformation)

            if self.env.sim_config.domain_randomization:
                image = self.env.randomizer.augment_tactile_image(image)

            tactile_images.append(image)
            depth_maps.append(deformation["depth_map"])
            normal_maps.append(deformation["normal_map"])
            wrenches.append(forces.sum(axis=0).tolist() + [0, 0, 0])

        dataset = {
            "tactile_images": np.stack(tactile_images),
            "depth_maps": np.stack(depth_maps),
            "normal_maps": np.stack(normal_maps),
            "wrenches": np.array(wrenches),
        }

        np.savez_compressed(
            os.path.join(self.save_dir, "contact_dataset.npz"),
            **dataset
        )
        return dataset

    def generate_texture_dataset(self, num_textures: int = 20,
                                  samples_per_texture: int = 500) -> Dict:
        """Generate texture classification training data."""
        import os
        os.makedirs(self.save_dir, exist_ok=True)

        images = []
        labels = []

        for tex_id in range(num_textures):
            for _ in range(samples_per_texture):
                # Create texture-specific contact pattern
                freq = 0.5 + tex_id * 0.2
                n_contacts = 10 + tex_id * 2
                points = np.random.randn(n_contacts, 3) * (0.002 + tex_id * 0.0005)
                forces = np.ones((n_contacts, 3)) * 0.3
                forces[:, 2] = 0.5 + np.sin(np.arange(n_contacts) * freq) * 0.2

                deformation = self.env.gel_model.compute_deformation(points, forces)
                image = self.env.renderer.render(deformation)

                if self.env.sim_config.domain_randomization:
                    image = self.env.randomizer.augment_tactile_image(image)

                images.append(image)
                labels.append(tex_id)

        dataset = {
            "images": np.stack(images),
            "labels": np.array(labels),
        }

        np.savez_compressed(
            os.path.join(self.save_dir, "texture_dataset.npz"),
            **dataset
        )
        return dataset
