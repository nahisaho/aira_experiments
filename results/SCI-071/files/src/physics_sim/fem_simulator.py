from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Mapping, Optional, Protocol, Tuple, Union, Literal

import numpy as np

from .base_simulator import (
    BaseSimulator,
    MaterialProperties,
    SimConfig,
    SimulatorError,
    State,
    primary_positions,
    rasterize_points,
)


@dataclass(slots=True)
class CollisionConfig:
    enabled: bool = True
    ground_height: float = 0.0
    penalty_stiffness: float = 5.0e4
    restitution: float = 0.0
    friction_coefficient: float = 0.5

    def __post_init__(self) -> None:
        if self.penalty_stiffness < 0:
            raise ValueError('penalty_stiffness cannot be negative.')
        if not 0.0 <= self.restitution <= 1.0:
            raise ValueError('restitution must lie in [0, 1].')
        if self.friction_coefficient < 0:
            raise ValueError('friction_coefficient cannot be negative.')


@dataclass(slots=True)
class MeshRefinementConfig:
    enabled: bool = False
    refinement_levels: int = 0
    target_edge_length: Optional[float] = None
    adaptive: bool = False

    def __post_init__(self) -> None:
        if self.refinement_levels < 0:
            raise ValueError('refinement_levels cannot be negative.')
        if self.target_edge_length is not None and self.target_edge_length <= 0:
            raise ValueError('target_edge_length must be positive when provided.')


class FEMBackend(Protocol):
    def reset(self, initial_state: Any = None) -> State: ...
    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]: ...
    def get_state(self) -> State: ...
    def set_state(self, state: State) -> None: ...
    def render(self, mode: str = 'rgb_array') -> np.ndarray: ...
    def get_jacobian(self, state: State, action: Any) -> np.ndarray: ...
    def refine_mesh(self, levels: int) -> State: ...


@dataclass(slots=True)
class FEMConfig(SimConfig):
    elasticity_model: Literal['linear', 'neo_hookean', 'stvk'] = 'neo_hookean'
    material: MaterialProperties = field(default_factory=MaterialProperties)
    implicit_integration: bool = True
    collision: CollisionConfig = field(default_factory=CollisionConfig)
    mesh_refinement: MeshRefinementConfig = field(default_factory=MeshRefinementConfig)
    rayleigh_damping: float = 0.05
    solver_tolerance: float = 1.0e-6
    max_newton_iterations: int = 25
    backend_factory: Optional[Callable[['FEMConfig'], FEMBackend]] = None
    backend_kwargs: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        SimConfig.__post_init__(self)
        if self.elasticity_model not in {'linear', 'neo_hookean', 'stvk'}:
            raise ValueError(f'Unsupported elasticity model: {self.elasticity_model}')
        if not self.implicit_integration:
            raise ValueError('FEMConfig requires implicit_integration=True for implicit Euler integration.')
        if self.rayleigh_damping < 0:
            raise ValueError('rayleigh_damping cannot be negative.')
        if self.solver_tolerance <= 0:
            raise ValueError('solver_tolerance must be positive.')
        if self.max_newton_iterations <= 0:
            raise ValueError('max_newton_iterations must be positive.')
        self.backend_kwargs = dict(self.backend_kwargs)


class _ReferenceFEMBackend:
    """Lightweight FEM-like backend for integration tests and fallback usage."""

    def __init__(self, config: FEMConfig) -> None:
        self.config = config
        self._state = State()
        self._rest_positions = np.zeros((0, 3), dtype=float)

    def reset(self, initial_state: Any = None) -> State:
        state = self._normalize_state(initial_state)
        positions = primary_positions(state)
        velocities = state.velocities
        if velocities is None:
            velocities = np.zeros_like(positions)

        state.positions = np.array(positions, copy=True)
        state.mesh_vertices = np.array(positions, copy=True)
        state.velocities = np.array(velocities, copy=True)
        state.done = False
        state.step_count = 0
        state.time = 0.0
        state.info.update(
            {
                'elasticity_model': self.config.elasticity_model,
                'integrator': 'implicit_euler',
                'mesh_refinement_levels': self.config.mesh_refinement.refinement_levels,
            }
        )
        self._state = state
        self._rest_positions = np.array(state.mesh_vertices, copy=True)
        return self.get_state()

    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]:
        state = self.get_state()
        positions = primary_positions(state)
        velocities = np.zeros_like(positions) if state.velocities is None else np.array(state.velocities, copy=True)
        nodal_forces = self._extract_nodal_forces(action, positions.shape)
        gravity = np.broadcast_to(np.asarray(self.config.gravity, dtype=float), positions.shape)
        density = self.config.material.density
        stiffness = self._effective_stiffness()
        damping = self.config.rayleigh_damping
        dt = self.config.time_step / max(1, self.config.substeps)

        for _ in range(self.config.substeps):
            internal_forces = -stiffness * (positions - self._rest_positions)
            accelerations = gravity + (nodal_forces + internal_forces) / density
            velocities = (velocities + dt * accelerations) / (1.0 + damping * dt)
            positions = positions + dt * velocities
            positions, velocities, contact_count = self._apply_collisions(positions, velocities)

        state.positions = np.array(positions, copy=True)
        state.mesh_vertices = np.array(positions, copy=True)
        state.velocities = np.array(velocities, copy=True)
        state.time += self.config.time_step
        state.step_count += 1
        state.info.update(
            {
                'contact_count': int(contact_count),
                'elastic_energy': float(0.5 * stiffness * np.mean((positions - self._rest_positions) ** 2)) if positions.size else 0.0,
                'mesh_refinement_levels': self.config.mesh_refinement.refinement_levels,
            }
        )

        self._state = state
        reward = -float(np.linalg.norm(positions - self._rest_positions, axis=1).mean()) if positions.size else 0.0
        return self.get_state(), reward, state.done, dict(state.info)

    def get_state(self) -> State:
        return self._state.copy()

    def set_state(self, state: State) -> None:
        coerced = state.copy()
        if coerced.mesh_vertices is None and coerced.positions is not None:
            coerced.mesh_vertices = np.array(coerced.positions, copy=True)
        if coerced.positions is None and coerced.mesh_vertices is not None:
            coerced.positions = np.array(coerced.mesh_vertices, copy=True)
        if coerced.velocities is None and coerced.positions is not None:
            coerced.velocities = np.zeros_like(coerced.positions)
        self._state = coerced
        if self._rest_positions.shape != primary_positions(coerced).shape:
            self._rest_positions = np.array(primary_positions(coerced), copy=True)

    def render(self, mode: str = 'rgb_array') -> np.ndarray:
        if mode != 'rgb_array':
            raise ValueError(f'Unsupported render mode: {mode}')
        return rasterize_points(self._state.mesh_vertices, self.config.render)

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
                perturbed_action = action_array.copy().reshape(-1)
                perturbed_action[column] += epsilon
                self.set_state(base_state)
                next_state, _, _, _ = self.step(perturbed_action.reshape(action_array.shape))
                jacobian[:, column] = (primary_positions(next_state).reshape(-1) - baseline_vector) / epsilon
            return jacobian
        finally:
            self.set_state(reference)

    def refine_mesh(self, levels: int) -> State:
        if levels < 0:
            raise ValueError('levels cannot be negative.')
        self.config.mesh_refinement.refinement_levels += levels
        self._state.info['mesh_refinement_levels'] = self.config.mesh_refinement.refinement_levels
        return self.get_state()

    def _normalize_state(self, initial_state: Any) -> State:
        if initial_state is None:
            positions = np.zeros((4, 3), dtype=float)
            positions[:, 0] = np.linspace(-0.25, 0.25, 4)
            positions[:, 1] = 0.5
            return State(positions=positions, mesh_vertices=positions)
        if isinstance(initial_state, Mapping):
            return State(**dict(initial_state))
        if isinstance(initial_state, State):
            return initial_state.copy()
        positions = np.asarray(initial_state, dtype=float)
        return State(positions=positions, mesh_vertices=positions)

    def _effective_stiffness(self) -> float:
        scale = {
            'linear': 1.0,
            'neo_hookean': 0.9,
            'stvk': 1.15,
        }[self.config.elasticity_model]
        return self.config.material.youngs_modulus * scale / max(1, len(primary_positions(self._state)))

    def _apply_collisions(self, positions: np.ndarray, velocities: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
        if not self.config.collision.enabled or positions.size == 0:
            return positions, velocities, 0
        updated_positions = np.array(positions, copy=True)
        updated_velocities = np.array(velocities, copy=True)
        penetration = self.config.collision.ground_height - updated_positions[:, 1]
        contact_mask = penetration > 0
        if not np.any(contact_mask):
            return updated_positions, updated_velocities, 0

        penalty_force = self.config.collision.penalty_stiffness * penetration[contact_mask]
        updated_velocities[contact_mask, 1] += penalty_force * self.config.time_step / self.config.material.density
        updated_velocities[contact_mask, 1] *= -self.config.collision.restitution
        updated_velocities[contact_mask, 0] *= max(0.0, 1.0 - self.config.collision.friction_coefficient)
        updated_velocities[contact_mask, 2] *= max(0.0, 1.0 - self.config.collision.friction_coefficient)
        updated_positions[contact_mask, 1] = self.config.collision.ground_height
        return updated_positions, updated_velocities, int(np.count_nonzero(contact_mask))

    @staticmethod
    def _extract_nodal_forces(action: Any, shape: Tuple[int, ...]) -> np.ndarray:
        if action is None:
            return np.zeros(shape, dtype=float)
        if isinstance(action, Mapping):
            if 'nodal_forces' in action:
                return _ReferenceFEMBackend._action_to_array(action['nodal_forces'], shape)
            if 'force' in action:
                return _ReferenceFEMBackend._action_to_array(action['force'], shape)
        return _ReferenceFEMBackend._action_to_array(action, shape)

    @staticmethod
    def _action_to_array(action: Any, shape: Tuple[int, ...]) -> np.ndarray:
        action_array = np.asarray(action, dtype=float)
        if action_array.size == 0:
            return np.zeros(shape, dtype=float)
        if action_array.shape == shape:
            return action_array
        if action_array.ndim == 1 and action_array.size == int(np.prod(shape)):
            return action_array.reshape(shape)
        if action_array.size == shape[-1]:
            return np.broadcast_to(action_array.reshape(1, -1), shape).copy()
        raise ValueError(f'Unable to broadcast action with shape {action_array.shape} to {shape}.')


class FEMSimulator(BaseSimulator):
    """Finite Element Method wrapper for cloth and elastic-body simulation."""

    def __init__(self, config: FEMConfig) -> None:
        super().__init__(config)
        self.config: FEMConfig = config
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

    def refine_mesh(self, levels: Optional[int] = None) -> State:
        if levels is None:
            levels = max(1, self.config.mesh_refinement.refinement_levels or 1)
        state = self._backend.refine_mesh(levels)
        return self._commit_state(state)

    def _build_backend(self) -> FEMBackend:
        if self.config.backend_factory is not None:
            return self.config.backend_factory(self.config)
        return _ReferenceFEMBackend(self.config)


__all__ = [
    'CollisionConfig',
    'FEMBackend',
    'FEMConfig',
    'FEMSimulator',
    'MeshRefinementConfig',
]
