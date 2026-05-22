from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union, Literal
import copy

import numpy as np

from .base_simulator import (
    BaseSimulator,
    DependencyUnavailableError,
    SimConfig,
    SimulatorError,
    State,
    primary_positions,
    rasterize_points,
)


@dataclass(slots=True)
class DomainRandomizationConfig:
    enabled: bool = False
    continuous_ranges: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    categorical_choices: Dict[str, Sequence[Any]] = field(default_factory=dict)

    def sample(self, rng: np.random.Generator) -> Dict[str, Any]:
        if not self.enabled:
            return {}
        sampled: Dict[str, Any] = {}
        for key, bounds in self.continuous_ranges.items():
            if len(bounds) != 2:
                raise ValueError(f'continuous range for {key!r} must have two bounds.')
            low, high = bounds
            if high < low:
                raise ValueError(f'continuous range for {key!r} is invalid: {bounds!r}')
            sampled[key] = float(rng.uniform(low, high))
        for key, choices in self.categorical_choices.items():
            if not choices:
                raise ValueError(f'categorical choices for {key!r} cannot be empty.')
            sampled[key] = copy.deepcopy(choices[int(rng.integers(0, len(choices)))])
        return sampled


@dataclass(slots=True)
class SoftGymConfig(SimConfig):
    env_name: Literal['ClothFold', 'ClothFlatten', 'RopeStraighten'] = 'ClothFold'
    observation_mode: str = 'cam_rgb'
    action_mode: str = 'pickerpickplace'
    render_mode: str = 'rgb_array'
    env_kwargs: Dict[str, Any] = field(default_factory=dict)
    domain_randomization: DomainRandomizationConfig = field(default_factory=DomainRandomizationConfig)

    def __post_init__(self) -> None:
        SimConfig.__post_init__(self)
        if self.env_name not in {'ClothFold', 'ClothFlatten', 'RopeStraighten'}:
            raise ValueError(f'Unsupported SoftGym environment: {self.env_name}')
        self.env_kwargs = dict(self.env_kwargs)


class SoftGymWrapper(BaseSimulator):
    """Wrapper around SoftGym cloth and rope manipulation environments."""

    def __init__(self, config: SoftGymConfig, env: Optional[Any] = None) -> None:
        super().__init__(config)
        self.config: SoftGymConfig = config
        self._env = env or self._build_env()
        self._last_observation: Any = None

    def reset(self, initial_state: Any = None) -> State:
        reset_kwargs = dict(self.config.env_kwargs)
        reset_kwargs.update(self.config.domain_randomization.sample(self._rng))
        if isinstance(initial_state, Mapping):
            reset_kwargs.update(dict(initial_state))

        try:
            observation = self._env.reset(**reset_kwargs)
        except TypeError:
            observation = self._env.reset()
            if initial_state is not None and hasattr(self._env, 'set_state'):
                self._env.set_state(initial_state)
                observation = self._observe()

        self._last_observation = observation
        state = self._extract_state(observation)
        return self._commit_state(state)

    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]:
        env_action = self._extract_action(action)
        observation, reward, done, info = self._env.step(env_action)
        self._last_observation = observation
        state = self._extract_state(observation, info=info, done=done)
        committed = self._commit_state(state)
        return committed, float(reward), bool(done), dict(info)

    def get_state(self) -> State:
        if hasattr(self._env, 'get_state'):
            return self._commit_state(self._coerce_state(self._env.get_state()))
        return super().get_state()

    def set_state(self, state: Union[State, Mapping[str, Any]]) -> None:
        coerced = self._coerce_state(state)
        if not hasattr(self._env, 'set_state'):
            raise SimulatorError('SoftGym environment does not expose set_state().')
        self._env.set_state(self._state_to_backend(coerced))
        self._state = coerced

    def render(self, mode: str = 'rgb_array') -> np.ndarray:
        if hasattr(self._env, 'render'):
            rendered = self._env.render(mode=mode)
            if isinstance(rendered, np.ndarray):
                return rendered
        if hasattr(self._env, 'get_image'):
            return np.asarray(self._env.get_image(self.config.render.width, self.config.render.height))
        return rasterize_points(primary_positions(self.get_state()), self.config.render, color=(214, 39, 40))

    def get_particle_positions(self) -> np.ndarray:
        return np.array(self._particle_positions_from_env(), copy=True)

    def get_picker_positions(self) -> Optional[np.ndarray]:
        picker_positions = self._picker_positions_from_env()
        if picker_positions is None:
            return None
        return np.array(picker_positions, copy=True)

    def get_jacobian(self, state: State, action: Any) -> np.ndarray:
        base_state = state.copy()
        action_array = np.asarray(self._extract_action(action), dtype=float)
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

    def _build_env(self) -> Any:
        try:
            from softgym.registered_env import SOFTGYM_ENVS, env_arg_dict
        except ImportError as exc:
            raise DependencyUnavailableError(
                'SoftGym is not installed. Install softgym and its PyFlex dependency to use SoftGymWrapper.'
            ) from exc

        env_kwargs = dict(env_arg_dict.get(self.config.env_name, {}))
        env_kwargs.update(self.config.env_kwargs)
        env_kwargs.setdefault('render_mode', self.config.render_mode)
        env_kwargs.setdefault('observation_mode', self.config.observation_mode)
        env_kwargs.setdefault('action_mode', self.config.action_mode)
        env_ctor = SOFTGYM_ENVS.get(self.config.env_name)
        if env_ctor is None:
            raise SimulatorError(f'SoftGym environment {self.config.env_name!r} is not registered.')
        return env_ctor(**env_kwargs)

    def _observe(self) -> Any:
        if hasattr(self._env, '_get_obs'):
            return self._env._get_obs()
        if self._last_observation is None:
            raise SimulatorError('No observation is available for the current SoftGym environment.')
        return self._last_observation

    def _extract_state(
        self,
        observation: Any,
        *,
        info: Optional[Mapping[str, Any]] = None,
        done: bool = False,
    ) -> State:
        particle_positions = self._particle_positions_from_env()
        picker_positions = self._picker_positions_from_env()
        state = State(
            particle_positions=particle_positions,
            positions=particle_positions,
            done=done,
            info={
                'env_name': self.config.env_name,
                'picker_positions': picker_positions,
                'domain_randomization': self.config.domain_randomization.enabled,
                'raw_observation': observation,
            },
            user_data={'picker_positions': picker_positions},
        )
        if info is not None:
            state.info.update(dict(info))
        return state

    def _particle_positions_from_env(self) -> np.ndarray:
        if hasattr(self._env, 'get_particle_positions'):
            positions = self._env.get_particle_positions()
        else:
            try:
                import pyflex
            except ImportError as exc:
                raise DependencyUnavailableError(
                    'SoftGym particle positions require PyFlex; install pyflex or provide get_particle_positions().' 
                ) from exc
            positions = pyflex.get_positions().reshape(-1, 4)[:, :3]
        return np.asarray(positions, dtype=float)

    def _picker_positions_from_env(self) -> Optional[np.ndarray]:
        action_tool = getattr(self._env, 'action_tool', None)
        if action_tool is None:
            return None
        for attribute in ('picker_pos', 'picker_positions', 'particle_positions'):
            if hasattr(action_tool, attribute):
                return np.asarray(getattr(action_tool, attribute), dtype=float)
        return None

    @staticmethod
    def _extract_action(action: Any) -> np.ndarray:
        if isinstance(action, Mapping):
            action = action.get('picker_action', action.get('action', action))
        return np.asarray(action, dtype=float)

    @staticmethod
    def _state_to_backend(state: State) -> Any:
        payload: Dict[str, Any] = {}
        if state.particle_positions is not None:
            payload['particle_pos'] = np.array(state.particle_positions, copy=True)
        payload.update(copy.deepcopy(state.user_data))
        return payload


__all__ = [
    'DomainRandomizationConfig',
    'SoftGymConfig',
    'SoftGymWrapper',
]
