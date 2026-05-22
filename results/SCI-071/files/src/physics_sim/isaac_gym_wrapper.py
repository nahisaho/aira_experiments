from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Mapping, Optional, Protocol, Tuple, Union

import numpy as np

from .base_simulator import (
    BaseSimulator,
    DependencyUnavailableError,
    SimConfig,
    State,
    primary_positions,
    rasterize_points,
)


class IsaacGymBackend(Protocol):
    def reset(self, initial_state: Any = None) -> State: ...
    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]: ...
    def get_state(self) -> State: ...
    def set_state(self, state: State) -> None: ...
    def render(self, mode: str = 'rgb_array') -> np.ndarray: ...
    def get_jacobian(self, state: State, action: Any) -> np.ndarray: ...


@dataclass(slots=True)
class IsaacGymConfig(SimConfig):
    num_envs: int = 1
    sim_device: str = 'cuda:0'
    graphics_device_id: int = 0
    use_gpu_pipeline: bool = True
    parallel_env_spacing: float = 1.5
    soft_body_solver_iterations: int = 20
    enable_stress_strain: bool = True
    backend_factory: Optional[Callable[['IsaacGymConfig'], IsaacGymBackend]] = None
    backend_kwargs: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        SimConfig.__post_init__(self)
        if self.num_envs <= 0:
            raise ValueError('num_envs must be positive.')
        if self.graphics_device_id < 0:
            raise ValueError('graphics_device_id cannot be negative.')
        if self.parallel_env_spacing <= 0:
            raise ValueError('parallel_env_spacing must be positive.')
        if self.soft_body_solver_iterations <= 0:
            raise ValueError('soft_body_solver_iterations must be positive.')
        self.backend_kwargs = dict(self.backend_kwargs)


class IsaacGymWrapper(BaseSimulator):
    """Wrapper for Isaac Gym FEM soft-body simulations with parallel GPU environments."""

    def __init__(self, config: IsaacGymConfig, backend: Optional[IsaacGymBackend] = None) -> None:
        super().__init__(config)
        self.config: IsaacGymConfig = config
        self._backend = backend or self._build_backend()

    def reset(self, initial_state: Any = None) -> State:
        return self._commit_state(self._backend.reset(initial_state))

    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]:
        state, reward, done, info = self._backend.step(action)
        committed = self._commit_state(state)
        return committed, float(reward), bool(done), dict(info)

    def get_state(self) -> State:
        return self._commit_state(self._backend.get_state())

    def set_state(self, state: Union[State, Mapping[str, Any]]) -> None:
        coerced = self._coerce_state(state)
        self._backend.set_state(coerced)
        self._state = coerced

    def render(self, mode: str = 'rgb_array') -> np.ndarray:
        try:
            return self._backend.render(mode=mode)
        except NotImplementedError:
            return rasterize_points(primary_positions(self.get_state()), self.config.render, color=(148, 103, 189))

    def get_jacobian(self, state: State, action: Any) -> np.ndarray:
        return self._backend.get_jacobian(state, action)

    def get_stress_tensors(self) -> Optional[np.ndarray]:
        state = self.get_state()
        return None if state.stress is None else np.array(state.stress, copy=True)

    def get_strain_tensors(self) -> Optional[np.ndarray]:
        state = self.get_state()
        return None if state.strain is None else np.array(state.strain, copy=True)

    def _build_backend(self) -> IsaacGymBackend:
        if self.config.backend_factory is not None:
            return self.config.backend_factory(self.config)

        try:
            import isaacgym  # noqa: F401
        except ImportError as exc:
            raise DependencyUnavailableError(
                'Isaac Gym is not installed. Install isaacgym or provide backend_factory/backend to IsaacGymWrapper.'
            ) from exc

        raise DependencyUnavailableError(
            'Isaac Gym was detected, but no backend_factory was supplied. Provide a backend that constructs '
            'the FEM soft-body assets, manages GPU-parallel environments, and returns stress/strain tensors.'
        )


__all__ = [
    'IsaacGymBackend',
    'IsaacGymConfig',
    'IsaacGymWrapper',
]
