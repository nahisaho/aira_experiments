from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Mapping, Optional, Protocol, Sequence, Tuple, Union, Literal

import numpy as np

from .base_simulator import (
    BaseSimulator,
    SimConfig,
    State,
    primary_positions,
    rasterize_points,
)


@dataclass(slots=True)
class ParticleMaterialConfig:
    material_type: Literal['elastic', 'fluid', 'sand', 'snow'] = 'elastic'
    density: float = 1.0e3
    youngs_modulus: float = 5.0e4
    poisson_ratio: float = 0.25
    viscosity: float = 0.0
    hardening: float = 1.0

    def __post_init__(self) -> None:
        if self.material_type not in {'elastic', 'fluid', 'sand', 'snow'}:
            raise ValueError(f'Unsupported material_type: {self.material_type}')
        if self.density <= 0:
            raise ValueError('density must be positive.')
        if self.youngs_modulus < 0:
            raise ValueError('youngs_modulus cannot be negative.')
        if not -0.999 < self.poisson_ratio < 0.5:
            raise ValueError('poisson_ratio must lie in (-0.999, 0.5).')
        if self.viscosity < 0:
            raise ValueError('viscosity cannot be negative.')
        if self.hardening <= 0:
            raise ValueError('hardening must be positive.')


class MPMBackend(Protocol):
    def reset(self, initial_state: Any = None) -> State: ...
    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]: ...
    def get_state(self) -> State: ...
    def set_state(self, state: State) -> None: ...
    def render(self, mode: str = 'rgb_array') -> np.ndarray: ...
    def get_jacobian(self, state: State, action: Any) -> np.ndarray: ...


@dataclass(slots=True)
class MPMConfig(SimConfig):
    grid_resolution: Tuple[int, int, int] = (32, 32, 32)
    grid_bounds: Tuple[Tuple[float, float, float], Tuple[float, float, float]] = ((-1.0, 0.0, -1.0), (1.0, 2.0, 1.0))
    materials: Dict[str, ParticleMaterialConfig] = field(default_factory=lambda: {'default': ParticleMaterialConfig()})
    apic_enabled: bool = True
    large_deformation_mode: bool = True
    cfl_number: float = 0.4
    boundary_friction: float = 0.0
    backend_factory: Optional[Callable[['MPMConfig'], MPMBackend]] = None
    backend_kwargs: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        SimConfig.__post_init__(self)
        if len(self.grid_resolution) != 3 or any(value <= 0 for value in self.grid_resolution):
            raise ValueError('grid_resolution must contain three positive integers.')
        if len(self.grid_bounds) != 2 or any(len(bound) != 3 for bound in self.grid_bounds):
            raise ValueError('grid_bounds must be ((xmin, ymin, zmin), (xmax, ymax, zmax)).')
        lower = np.asarray(self.grid_bounds[0], dtype=float)
        upper = np.asarray(self.grid_bounds[1], dtype=float)
        if np.any(upper <= lower):
            raise ValueError('grid_bounds upper corner must exceed lower corner in every axis.')
        if not self.materials:
            raise ValueError('At least one material configuration must be provided.')
        if not self.apic_enabled:
            raise ValueError('MPMConfig requires apic_enabled=True for APIC transfers.')
        if self.cfl_number <= 0:
            raise ValueError('cfl_number must be positive.')
        if self.boundary_friction < 0:
            raise ValueError('boundary_friction cannot be negative.')
        self.backend_kwargs = dict(self.backend_kwargs)


class _ReferenceMPMBackend:
    """Compact APIC-style MPM backend for development and integration tests."""

    def __init__(self, config: MPMConfig) -> None:
        self.config = config
        self._state = State()
        self._material_ids = np.zeros(0, dtype=int)
        self._material_names = list(config.materials.keys())
        self._affine = np.zeros((0, 3, 3), dtype=float)

    def reset(self, initial_state: Any = None) -> State:
        state = self._normalize_state(initial_state)
        particles = primary_positions(state)
        velocities = np.zeros_like(particles) if state.particle_velocities is None else np.array(state.particle_velocities, copy=True)
        if state.particle_positions is None:
            state.particle_positions = np.array(particles, copy=True)
        state.positions = np.array(state.particle_positions, copy=True)
        state.particle_velocities = velocities
        state.velocities = np.array(velocities, copy=True)
        state.step_count = 0
        state.time = 0.0
        state.done = False

        material_ids = state.user_data.get('material_ids')
        if material_ids is None:
            self._material_ids = np.zeros(len(state.particle_positions), dtype=int)
        else:
            self._material_ids = np.asarray(material_ids, dtype=int)
        self._affine = np.zeros((len(state.particle_positions), 3, 3), dtype=float)
        state.user_data['material_ids'] = np.array(self._material_ids, copy=True)
        state.info.update(
            {
                'grid_resolution': self.config.grid_resolution,
                'apic_enabled': self.config.apic_enabled,
                'large_deformation_mode': self.config.large_deformation_mode,
            }
        )
        self._state = state
        return self.get_state()

    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]:
        state = self.get_state()
        particles = np.array(state.particle_positions, copy=True)
        velocities = np.array(state.particle_velocities, copy=True)
        external = self._action_to_array(action, particles.shape)
        lower = np.asarray(self.config.grid_bounds[0], dtype=float)
        upper = np.asarray(self.config.grid_bounds[1], dtype=float)
        dt = self.config.time_step / max(1, self.config.substeps)

        for _ in range(self.config.substeps):
            grid_mass, grid_velocity, dominant_material = self._particle_to_grid(particles, velocities + external)
            grid_velocity = self._apply_grid_dynamics(grid_velocity, grid_mass, dominant_material)
            velocities = self._grid_to_particle(particles, grid_velocity) + dt * external
            particles = particles + dt * velocities
            particles, velocities = self._apply_boundaries(particles, velocities, lower, upper)

        state.particle_positions = np.array(particles, copy=True)
        state.positions = np.array(particles, copy=True)
        state.particle_velocities = np.array(velocities, copy=True)
        state.velocities = np.array(velocities, copy=True)
        state.step_count += 1
        state.time += self.config.time_step
        state.user_data['material_ids'] = np.array(self._material_ids, copy=True)
        state.info.update(
            {
                'mean_speed': float(np.linalg.norm(velocities, axis=1).mean()) if velocities.size else 0.0,
                'grid_resolution': self.config.grid_resolution,
            }
        )
        self._state = state
        reward = -float(np.linalg.norm(velocities, axis=1).mean()) if velocities.size else 0.0
        return self.get_state(), reward, state.done, dict(state.info)

    def get_state(self) -> State:
        return self._state.copy()

    def set_state(self, state: State) -> None:
        coerced = state.copy()
        particles = primary_positions(coerced)
        if coerced.particle_positions is None:
            coerced.particle_positions = np.array(particles, copy=True)
        if coerced.positions is None:
            coerced.positions = np.array(coerced.particle_positions, copy=True)
        if coerced.particle_velocities is None:
            coerced.particle_velocities = np.zeros_like(coerced.particle_positions)
        if coerced.velocities is None:
            coerced.velocities = np.array(coerced.particle_velocities, copy=True)
        material_ids = coerced.user_data.get('material_ids')
        if material_ids is None:
            self._material_ids = np.zeros(len(coerced.particle_positions), dtype=int)
        else:
            self._material_ids = np.asarray(material_ids, dtype=int)
        self._affine = np.zeros((len(coerced.particle_positions), 3, 3), dtype=float)
        self._state = coerced

    def render(self, mode: str = 'rgb_array') -> np.ndarray:
        if mode != 'rgb_array':
            raise ValueError(f'Unsupported render mode: {mode}')
        return rasterize_points(self._state.particle_positions, self.config.render, color=(44, 160, 44))

    def get_jacobian(self, state: State, action: Any) -> np.ndarray:
        base_state = state.copy()
        action_array = self._action_to_array(action, primary_positions(base_state).shape)
        if action_array.size == 0:
            return np.zeros((primary_positions(base_state).size, 0), dtype=float)

        epsilon = 1.0e-4
        reference = self.get_state()
        try:
            self.set_state(base_state)
            baseline, _, _, _ = self.step(action_array)
            baseline_vector = primary_positions(baseline).reshape(-1)
            jacobian = np.zeros((baseline_vector.size, action_array.size), dtype=float)
            for column in range(action_array.size):
                perturbed = action_array.reshape(-1).copy()
                perturbed[column] += epsilon
                self.set_state(base_state)
                next_state, _, _, _ = self.step(perturbed.reshape(action_array.shape))
                jacobian[:, column] = (primary_positions(next_state).reshape(-1) - baseline_vector) / epsilon
            return jacobian
        finally:
            self.set_state(reference)

    def _normalize_state(self, initial_state: Any) -> State:
        if initial_state is None:
            particles = np.array(
                [
                    [-0.2, 0.8, -0.2],
                    [-0.2, 0.8, 0.2],
                    [0.2, 0.8, -0.2],
                    [0.2, 0.8, 0.2],
                ],
                dtype=float,
            )
            return State(particle_positions=particles, positions=particles)
        if isinstance(initial_state, Mapping):
            return State(**dict(initial_state))
        if isinstance(initial_state, State):
            return initial_state.copy()
        particles = np.asarray(initial_state, dtype=float)
        return State(particle_positions=particles, positions=particles)

    def _particle_to_grid(self, particles: np.ndarray, velocities: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        resolution = np.asarray(self.config.grid_resolution, dtype=int)
        lower = np.asarray(self.config.grid_bounds[0], dtype=float)
        upper = np.asarray(self.config.grid_bounds[1], dtype=float)
        grid_mass = np.zeros(tuple(resolution), dtype=float)
        grid_velocity = np.zeros(tuple(resolution) + (3,), dtype=float)
        material_votes = np.zeros(tuple(resolution) + (len(self._material_names),), dtype=float)
        if particles.size == 0:
            return grid_mass, grid_velocity, np.zeros(tuple(resolution), dtype=int)

        normalized = np.clip((particles - lower) / (upper - lower), 0.0, 0.999999)
        indices = np.floor(normalized * resolution).astype(int)
        material_ids = np.clip(self._material_ids, 0, len(self._material_names) - 1)
        np.add.at(grid_mass, tuple(indices.T), 1.0)
        for axis in range(3):
            np.add.at(grid_velocity[..., axis], tuple(indices.T), velocities[:, axis])
        for material_index in range(len(self._material_names)):
            mask = material_ids == material_index
            if np.any(mask):
                np.add.at(material_votes[..., material_index], tuple(indices[mask].T), 1.0)
        nonzero = grid_mass > 0
        grid_velocity[nonzero] /= grid_mass[nonzero][..., None]
        dominant_material = np.argmax(material_votes, axis=-1)
        return grid_mass, grid_velocity, dominant_material

    def _apply_grid_dynamics(
        self,
        grid_velocity: np.ndarray,
        grid_mass: np.ndarray,
        dominant_material: np.ndarray,
    ) -> np.ndarray:
        updated = np.array(grid_velocity, copy=True)
        updated += self.config.time_step * np.asarray(self.config.gravity, dtype=float)
        nonzero = grid_mass > 0
        if not np.any(nonzero):
            return updated

        for material_index, name in enumerate(self._material_names):
            material = self.config.materials[name]
            material_mask = nonzero & (dominant_material == material_index)
            if not np.any(material_mask):
                continue
            if material.material_type == 'fluid':
                updated[material_mask] *= max(0.0, 1.0 - material.viscosity * self.config.time_step)
            elif material.material_type == 'sand':
                updated[material_mask, 0] *= 0.8
                updated[material_mask, 2] *= 0.8
            elif material.material_type == 'snow':
                updated[material_mask] *= min(1.5, material.hardening)
            else:
                updated[material_mask] *= 1.0 + min(0.25, material.youngs_modulus / 1.0e6)
        return updated

    def _grid_to_particle(self, particles: np.ndarray, grid_velocity: np.ndarray) -> np.ndarray:
        if particles.size == 0:
            return np.zeros_like(particles)
        resolution = np.asarray(self.config.grid_resolution, dtype=int)
        lower = np.asarray(self.config.grid_bounds[0], dtype=float)
        upper = np.asarray(self.config.grid_bounds[1], dtype=float)
        normalized = np.clip((particles - lower) / (upper - lower), 0.0, 0.999999)
        indices = np.floor(normalized * resolution).astype(int)
        particle_velocity = grid_velocity[tuple(indices.T)]
        if self.config.apic_enabled:
            self._affine = 0.5 * self._affine
            particle_velocity = particle_velocity + np.einsum('nij,nj->ni', self._affine, particles - particles.mean(axis=0, keepdims=True))
        return particle_velocity

    def _apply_boundaries(
        self,
        particles: np.ndarray,
        velocities: np.ndarray,
        lower: np.ndarray,
        upper: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        bounded_particles = np.array(particles, copy=True)
        bounded_velocities = np.array(velocities, copy=True)
        for axis in range(3):
            lower_mask = bounded_particles[:, axis] < lower[axis]
            upper_mask = bounded_particles[:, axis] > upper[axis]
            bounded_particles[lower_mask, axis] = lower[axis]
            bounded_particles[upper_mask, axis] = upper[axis]
            bounded_velocities[lower_mask | upper_mask, axis] *= -(1.0 - self.config.boundary_friction)
        return bounded_particles, bounded_velocities

    @staticmethod
    def _action_to_array(action: Any, shape: Tuple[int, ...]) -> np.ndarray:
        if action is None:
            return np.zeros(shape, dtype=float)
        if isinstance(action, Mapping):
            action = action.get('particle_forces', action.get('force', action))
        action_array = np.asarray(action, dtype=float)
        if action_array.shape == shape:
            return action_array
        if action_array.ndim == 1 and action_array.size == int(np.prod(shape)):
            return action_array.reshape(shape)
        if action_array.size == shape[-1]:
            return np.broadcast_to(action_array.reshape(1, -1), shape).copy()
        raise ValueError(f'Unable to broadcast action with shape {action_array.shape} to {shape}.')


class MPMSimulator(BaseSimulator):
    """Material Point Method simulator wrapper for multi-material deformables."""

    def __init__(self, config: MPMConfig) -> None:
        super().__init__(config)
        self.config: MPMConfig = config
        self._backend = self._build_backend()

    def reset(self, initial_state: Any = None) -> State:
        return self._commit_state(self._backend.reset(initial_state))

    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]:
        state, reward, done, info = self._backend.step(action)
        committed = self._commit_state(state)
        return committed, reward, done, info

    def get_state(self) -> State:
        return self._commit_state(self._backend.get_state())

    def set_state(self, state: Union[State, Mapping[str, Any]]) -> None:
        coerced = self._coerce_state(state)
        self._backend.set_state(coerced)
        self._state = coerced

    def render(self, mode: str = 'rgb_array') -> np.ndarray:
        return self._backend.render(mode=mode)

    def get_jacobian(self, state: State, action: Any) -> np.ndarray:
        return self._backend.get_jacobian(state, action)

    def _build_backend(self) -> MPMBackend:
        if self.config.backend_factory is not None:
            return self.config.backend_factory(self.config)
        return _ReferenceMPMBackend(self.config)


__all__ = [
    'MPMBackend',
    'MPMConfig',
    'MPMSimulator',
    'ParticleMaterialConfig',
]
