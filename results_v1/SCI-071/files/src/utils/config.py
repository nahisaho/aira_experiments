"""Configuration management for cloth folding experiments."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

try:
    import yaml
except ImportError:  # pragma: no cover - optional dependency guard
    yaml = None


@dataclass
class MaterialConfig:
    """Material parameters for a cloth model."""

    density: float = 0.2
    stretching_stiffness: float = 1.0
    bending_stiffness: float = 0.05
    damping: float = 0.01
    friction: float = 0.6
    thickness: float = 0.002

    def validate(self) -> None:
        if self.density <= 0.0:
            raise ValueError("density must be positive")
        if self.stretching_stiffness <= 0.0:
            raise ValueError("stretching_stiffness must be positive")
        if self.bending_stiffness < 0.0:
            raise ValueError("bending_stiffness must be non-negative")
        if not 0.0 <= self.friction <= 2.0:
            raise ValueError("friction must be in [0, 2]")
        if self.thickness <= 0.0:
            raise ValueError("thickness must be positive")


@dataclass
class MeshConfig:
    """Geometry settings for cloth mesh construction."""

    width: float = 0.6
    height: float = 0.6
    resolution_x: int = 21
    resolution_y: int = 21
    origin: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])

    def validate(self) -> None:
        if self.width <= 0.0 or self.height <= 0.0:
            raise ValueError("cloth width and height must be positive")
        if self.resolution_x < 2 or self.resolution_y < 2:
            raise ValueError("mesh resolution must be at least 2x2")
        if len(self.origin) != 3:
            raise ValueError("origin must contain exactly three values")


@dataclass
class CameraConfig:
    """Top-down observation camera configuration."""

    width: int = 224
    height: int = 224
    near: float = 0.01
    far: float = 2.0
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, 1.0])
    target: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])

    def validate(self) -> None:
        if self.width <= 0 or self.height <= 0:
            raise ValueError("camera image size must be positive")
        if self.near <= 0.0 or self.far <= self.near:
            raise ValueError("camera clipping planes are invalid")
        if len(self.position) != 3 or len(self.target) != 3:
            raise ValueError("camera position and target must be 3D vectors")


@dataclass
class RewardConfig:
    """Weights and thresholds for reward computation."""

    chamfer_weight: float = 1.0
    coverage_weight: float = 0.25
    fold_line_weight: float = 0.2
    smoothness_weight: float = 0.05
    success_threshold: float = 0.92

    def validate(self) -> None:
        if self.success_threshold <= 0.0:
            raise ValueError("success_threshold must be positive")


@dataclass
class BackendConfig:
    """Backend selection and simulation execution options."""

    backend: str = "softgym"
    headless: bool = True
    device: str = "cpu"
    substeps: int = 8
    solver_iterations: int = 10

    def validate(self) -> None:
        valid_backends = {"softgym", "isaac_gym", "numpy"}
        if self.backend not in valid_backends:
            raise ValueError(f"backend must be one of {sorted(valid_backends)}")
        if self.substeps <= 0 or self.solver_iterations <= 0:
            raise ValueError("substeps and solver_iterations must be positive")


@dataclass
class ClothFoldingEnvConfig:
    """Complete configuration for the cloth folding environment."""

    mesh: MeshConfig = field(default_factory=MeshConfig)
    material: MaterialConfig = field(default_factory=MaterialConfig)
    camera: CameraConfig = field(default_factory=CameraConfig)
    reward: RewardConfig = field(default_factory=RewardConfig)
    backend: BackendConfig = field(default_factory=BackendConfig)
    max_episode_steps: int = 50
    grasp_height: float = 0.15
    grasp_radius: float = 0.08
    randomize_pose: bool = True
    randomize_height: float = 0.03
    supported_fold_types: List[str] = field(
        default_factory=lambda: ["half_fold", "quarter_fold", "diagonal_fold", "sleeve_fold"]
    )
    seed: Optional[int] = None

    def validate(self) -> None:
        self.mesh.validate()
        self.material.validate()
        self.camera.validate()
        self.reward.validate()
        self.backend.validate()
        if self.max_episode_steps <= 0:
            raise ValueError("max_episode_steps must be positive")
        if self.grasp_height <= 0.0 or self.grasp_radius <= 0.0:
            raise ValueError("grasp_height and grasp_radius must be positive")
        valid_fold_types = {"half_fold", "quarter_fold", "diagonal_fold", "sleeve_fold"}
        unknown = set(self.supported_fold_types) - valid_fold_types
        if unknown:
            raise ValueError(f"unsupported fold types: {sorted(unknown)}")

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the configuration into a dictionary."""

        return asdict(self)

    def to_yaml(self, path: str | Path) -> None:
        """Save the configuration to a YAML file."""

        if yaml is None:
            raise ImportError("PyYAML is required to save YAML configuration files")
        self.validate()
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(self.to_dict(), handle, sort_keys=False)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ClothFoldingEnvConfig":
        """Construct a configuration from a nested mapping."""

        payload = dict(data)
        mesh = MeshConfig(**payload.pop("mesh", {}))
        material = MaterialConfig(**payload.pop("material", {}))
        camera = CameraConfig(**payload.pop("camera", {}))
        reward = RewardConfig(**payload.pop("reward", {}))
        backend = BackendConfig(**payload.pop("backend", {}))
        config = cls(mesh=mesh, material=material, camera=camera, reward=reward, backend=backend, **payload)
        config.validate()
        return config

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ClothFoldingEnvConfig":
        """Load configuration parameters from a YAML file."""

        if yaml is None:
            raise ImportError("PyYAML is required to load YAML configuration files")
        with Path(path).open("r", encoding="utf-8") as handle:
            payload = yaml.safe_load(handle) or {}
        return cls.from_dict(payload)


def default_cloth_folding_config() -> ClothFoldingEnvConfig:
    """Return a validated default configuration for cloth folding experiments."""

    config = ClothFoldingEnvConfig()
    config.validate()
    return config


__all__ = [
    "BackendConfig",
    "CameraConfig",
    "ClothFoldingEnvConfig",
    "MaterialConfig",
    "MeshConfig",
    "RewardConfig",
    "default_cloth_folding_config",
]
