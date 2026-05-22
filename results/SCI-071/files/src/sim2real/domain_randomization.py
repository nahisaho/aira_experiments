"""Domain randomization utilities for deformable object manipulation."""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass, replace
import math
import random
from typing import Any, Mapping, MutableMapping, Optional

import torch


@dataclass(frozen=True)
class NumericRange:
    """Closed interval for continuous parameter sampling."""

    low: float
    high: float

    def __post_init__(self) -> None:
        if self.low > self.high:
            raise ValueError("Range lower bound must not exceed upper bound.")

    @property
    def width(self) -> float:
        return self.high - self.low

    @property
    def center(self) -> float:
        return 0.5 * (self.low + self.high)

    def sample(self, generator: Optional[torch.Generator] = None) -> float:
        if math.isclose(self.low, self.high):
            return float(self.low)
        if generator is None:
            return random.uniform(self.low, self.high)
        return float(torch.empty(1).uniform_(self.low, self.high, generator=generator).item())

    def scale(self, factor: float) -> "NumericRange":
        half_width = 0.5 * self.width * max(factor, 0.0)
        center = self.center
        return NumericRange(center - half_width, center + half_width)

    def expand(self, factor: float) -> "NumericRange":
        return self.scale(1.0 + max(factor, -0.99))


@dataclass(frozen=True)
class IntRange:
    """Closed interval for integer parameter sampling."""

    low: int
    high: int

    def __post_init__(self) -> None:
        if self.low > self.high:
            raise ValueError("Range lower bound must not exceed upper bound.")

    @property
    def width(self) -> int:
        return self.high - self.low

    @property
    def center(self) -> float:
        return 0.5 * (self.low + self.high)

    def sample(self, generator: Optional[torch.Generator] = None) -> int:
        if self.low == self.high:
            return self.low
        if generator is None:
            return random.randint(self.low, self.high)
        value = torch.randint(self.low, self.high + 1, (1,), generator=generator)
        return int(value.item())

    def scale(self, factor: float) -> "IntRange":
        half_width = max(int(round(0.5 * self.width * max(factor, 0.0))), 0)
        center = self.center
        low = int(math.floor(center - half_width))
        high = int(math.ceil(center + half_width))
        return IntRange(low, max(low, high))

    def expand(self, factor: float) -> "IntRange":
        return self.scale(1.0 + max(factor, -0.99))


@dataclass(frozen=True)
class MaterialRandomizationConfig:
    """Physical parameter randomization for deformable bodies."""

    stiffness: NumericRange = NumericRange(50.0, 500.0)
    damping: NumericRange = NumericRange(0.01, 2.0)
    friction: NumericRange = NumericRange(0.2, 1.2)
    mass: NumericRange = NumericRange(0.05, 5.0)


@dataclass(frozen=True)
class VisualRandomizationConfig:
    """Visual appearance randomization."""

    texture_blend: NumericRange = NumericRange(0.0, 1.0)
    light_intensity: NumericRange = NumericRange(300.0, 1800.0)
    light_color_temperature: NumericRange = NumericRange(3200.0, 7200.0)
    camera_translation_jitter: NumericRange = NumericRange(-0.03, 0.03)
    camera_rotation_jitter_deg: NumericRange = NumericRange(-8.0, 8.0)


@dataclass(frozen=True)
class DynamicsRandomizationConfig:
    """Simulation dynamics randomization."""

    timestep: NumericRange = NumericRange(1.0 / 480.0, 1.0 / 120.0)
    solver_iterations: IntRange = IntRange(4, 32)
    contact_stiffness: NumericRange = NumericRange(100.0, 2000.0)
    contact_damping: NumericRange = NumericRange(0.01, 0.5)
    restitution: NumericRange = NumericRange(0.0, 0.2)


@dataclass(frozen=True)
class NoiseRandomizationConfig:
    """Action and observation perturbation settings."""

    action_std: float = 0.01
    action_bias: float = 0.0
    observation_std: float = 0.005
    observation_dropout: float = 0.0

    def __post_init__(self) -> None:
        if not 0.0 <= self.observation_dropout < 1.0:
            raise ValueError("Observation dropout must be in [0, 1).")


@dataclass(frozen=True)
class ADRConfig:
    """Automatic Domain Randomization parameters."""

    enabled: bool = True
    success_threshold: float = 0.8
    failure_threshold: float = 0.3
    expansion_rate: float = 0.1
    contraction_rate: float = 0.05
    window_size: int = 20

    def __post_init__(self) -> None:
        if self.window_size <= 0:
            raise ValueError("window_size must be positive.")
        if self.failure_threshold > self.success_threshold:
            raise ValueError("failure_threshold must not exceed success_threshold.")


@dataclass(frozen=True)
class DomainRandomizationConfig:
    """Configuration bundle for domain randomization."""

    material: MaterialRandomizationConfig = field(default_factory=MaterialRandomizationConfig)
    visual: VisualRandomizationConfig = field(default_factory=VisualRandomizationConfig)
    dynamics: DynamicsRandomizationConfig = field(default_factory=DynamicsRandomizationConfig)
    noise: NoiseRandomizationConfig = field(default_factory=NoiseRandomizationConfig)
    adr: ADRConfig = field(default_factory=ADRConfig)


@dataclass(frozen=True)
class DomainParameters:
    """Concrete randomized parameters sampled for one rollout."""

    material: Mapping[str, float]
    visual: Mapping[str, float]
    dynamics: Mapping[str, float | int]

    def as_flat_dict(self) -> dict[str, float | int]:
        flat: dict[str, float | int] = {}
        for prefix, values in (
            ("material", self.material),
            ("visual", self.visual),
            ("dynamics", self.dynamics),
        ):
            for key, value in values.items():
                flat[f"{prefix}.{key}"] = value
        return flat


class AutomaticDomainRandomization:
    """Performance-driven randomization boundary adaptation."""

    def __init__(self, config: DomainRandomizationConfig) -> None:
        self._base_config = config
        self._current_config = config
        self._performance_window: list[float] = []

    @property
    def current_config(self) -> DomainRandomizationConfig:
        return self._current_config

    def update(self, performance: float) -> DomainRandomizationConfig:
        adr = self._current_config.adr
        if not adr.enabled:
            return self._current_config

        self._performance_window.append(float(performance))
        if len(self._performance_window) > adr.window_size:
            self._performance_window.pop(0)
        if len(self._performance_window) < adr.window_size:
            return self._current_config

        mean_performance = sum(self._performance_window) / len(self._performance_window)
        if mean_performance >= adr.success_threshold:
            self._current_config = self._apply_scale(
                self._current_config,
                1.0 + adr.expansion_rate,
            )
        elif mean_performance <= adr.failure_threshold:
            shrink = max(1.0 - adr.contraction_rate, 1e-3)
            contracted = self._apply_scale(self._current_config, shrink)
            self._current_config = self._blend_towards_base(contracted)
        return self._current_config

    def _blend_towards_base(self, config: DomainRandomizationConfig) -> DomainRandomizationConfig:
        return DomainRandomizationConfig(
            material=self._blend_dataclass(config.material, self._base_config.material),
            visual=self._blend_dataclass(config.visual, self._base_config.visual),
            dynamics=self._blend_dataclass(config.dynamics, self._base_config.dynamics),
            noise=config.noise,
            adr=config.adr,
        )

    def _blend_dataclass(self, current: Any, base: Any) -> Any:
        updates: MutableMapping[str, Any] = {}
        for item in fields(current):
            current_value = getattr(current, item.name)
            base_value = getattr(base, item.name)
            if isinstance(current_value, NumericRange) and isinstance(base_value, NumericRange):
                low = min(current_value.low, base_value.low)
                high = max(current_value.high, base_value.high)
                updates[item.name] = NumericRange(low, high)
            elif isinstance(current_value, IntRange) and isinstance(base_value, IntRange):
                low = min(current_value.low, base_value.low)
                high = max(current_value.high, base_value.high)
                updates[item.name] = IntRange(low, high)
            else:
                updates[item.name] = current_value
        return replace(current, **updates)

    def _apply_scale(self, config: DomainRandomizationConfig, factor: float) -> DomainRandomizationConfig:
        return DomainRandomizationConfig(
            material=_scale_dataclass_ranges(config.material, factor),
            visual=_scale_dataclass_ranges(config.visual, factor),
            dynamics=_scale_dataclass_ranges(config.dynamics, factor),
            noise=config.noise,
            adr=config.adr,
        )


class DomainRandomizer:
    """Sampler and perturbation helper for sim-to-real training."""

    def __init__(
        self,
        config: DomainRandomizationConfig,
        *,
        seed: Optional[int] = None,
        device: Optional[torch.device] = None,
    ) -> None:
        self._generator = torch.Generator(device=device)
        if seed is not None:
            self._generator.manual_seed(seed)
            random.seed(seed)
        self._adr = AutomaticDomainRandomization(config)

    @property
    def config(self) -> DomainRandomizationConfig:
        return self._adr.current_config

    def sample_domain(self) -> DomainParameters:
        config = self.config
        return DomainParameters(
            material={
                "stiffness": config.material.stiffness.sample(self._generator),
                "damping": config.material.damping.sample(self._generator),
                "friction": config.material.friction.sample(self._generator),
                "mass": config.material.mass.sample(self._generator),
            },
            visual={
                "texture_blend": config.visual.texture_blend.sample(self._generator),
                "light_intensity": config.visual.light_intensity.sample(self._generator),
                "light_color_temperature": config.visual.light_color_temperature.sample(self._generator),
                "camera_translation_jitter": config.visual.camera_translation_jitter.sample(self._generator),
                "camera_rotation_jitter_deg": config.visual.camera_rotation_jitter_deg.sample(self._generator),
            },
            dynamics={
                "timestep": config.dynamics.timestep.sample(self._generator),
                "solver_iterations": config.dynamics.solver_iterations.sample(self._generator),
                "contact_stiffness": config.dynamics.contact_stiffness.sample(self._generator),
                "contact_damping": config.dynamics.contact_damping.sample(self._generator),
                "restitution": config.dynamics.restitution.sample(self._generator),
            },
        )

    def inject_action_noise(self, action: torch.Tensor) -> torch.Tensor:
        noise_cfg = self.config.noise
        noise = torch.randn_like(action, generator=self._generator) * noise_cfg.action_std
        return action + noise + noise_cfg.action_bias

    def inject_observation_noise(self, observation: torch.Tensor) -> torch.Tensor:
        noise_cfg = self.config.noise
        noisy = observation + torch.randn_like(observation, generator=self._generator) * noise_cfg.observation_std
        if noise_cfg.observation_dropout <= 0.0:
            return noisy
        mask = torch.rand_like(noisy, generator=self._generator) >= noise_cfg.observation_dropout
        return noisy * mask

    def update_adr(self, performance: float) -> DomainRandomizationConfig:
        return self._adr.update(performance)


def _scale_dataclass_ranges(instance: Any, factor: float) -> Any:
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
            updates[item.name] = _scale_dataclass_ranges(value, factor)
        else:
            updates[item.name] = value
    return replace(instance, **updates)
