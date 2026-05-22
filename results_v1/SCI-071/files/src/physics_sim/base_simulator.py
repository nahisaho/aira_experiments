from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, fields as dataclass_fields
from typing import Any, Dict, Mapping, Optional, Tuple, Union
import copy

import numpy as np


class SimulatorError(RuntimeError):
    """Base error for simulator integration failures."""


class DependencyUnavailableError(ImportError, SimulatorError):
    """Raised when an optional simulator dependency is unavailable."""


@dataclass(slots=True)
class RenderConfig:
    width: int = 640
    height: int = 480
    point_radius: int = 2
    background_color: Tuple[int, int, int] = (255, 255, 255)

    def __post_init__(self) -> None:
        if self.width <= 0 or self.height <= 0:
            raise ValueError('Render dimensions must be positive.')
        if self.point_radius < 1:
            raise ValueError('point_radius must be at least 1.')
        if len(self.background_color) != 3 or any(not 0 <= c <= 255 for c in self.background_color):
            raise ValueError('background_color must contain three 8-bit RGB values.')


@dataclass(slots=True)
class MaterialProperties:
    youngs_modulus: float = 1.0e5
    poisson_ratio: float = 0.3
    density: float = 1.0e3

    def __post_init__(self) -> None:
        if self.youngs_modulus <= 0:
            raise ValueError("youngs_modulus must be positive.")
        if not -0.999 < self.poisson_ratio < 0.5:
            raise ValueError('poisson_ratio must lie in (-0.999, 0.5).')
        if self.density <= 0:
            raise ValueError('density must be positive.')


@dataclass(slots=True)
class SimConfig:
    time_step: float = 1.0 / 60.0
    substeps: int = 1
    gravity: Tuple[float, float, float] = (0.0, -9.81, 0.0)
    device: str = 'cpu'
    seed: Optional[int] = None
    render: RenderConfig = field(default_factory=RenderConfig)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.time_step <= 0:
            raise ValueError('time_step must be positive.')
        if self.substeps <= 0:
            raise ValueError('substeps must be a positive integer.')
        if len(self.gravity) != 3:
            raise ValueError('gravity must be a 3-tuple.')
        self.gravity = tuple(float(component) for component in self.gravity)
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class SimulatorState:
    time: float = 0.0
    step_count: int = 0
    positions: Optional[np.ndarray] = None
    velocities: Optional[np.ndarray] = None
    particle_positions: Optional[np.ndarray] = None
    particle_velocities: Optional[np.ndarray] = None
    mesh_vertices: Optional[np.ndarray] = None
    mesh_elements: Optional[np.ndarray] = None
    stress: Optional[np.ndarray] = None
    strain: Optional[np.ndarray] = None
    done: bool = False
    info: Dict[str, Any] = field(default_factory=dict)
    user_data: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        float_arrays = (
            'positions',
            'velocities',
            'particle_positions',
            'particle_velocities',
            'mesh_vertices',
            'stress',
            'strain',
        )
        for name in float_arrays:
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, np.asarray(value, dtype=float))

        if self.mesh_elements is not None:
            self.mesh_elements = np.asarray(self.mesh_elements, dtype=int)

        self.time = float(self.time)
        self.step_count = int(self.step_count)
        self.done = bool(self.done)
        self.info = copy.deepcopy(dict(self.info))
        self.user_data = copy.deepcopy(dict(self.user_data))

    def copy(self) -> 'SimulatorState':
        return SimulatorState(
            time=self.time,
            step_count=self.step_count,
            positions=_copy_array(self.positions),
            velocities=_copy_array(self.velocities),
            particle_positions=_copy_array(self.particle_positions),
            particle_velocities=_copy_array(self.particle_velocities),
            mesh_vertices=_copy_array(self.mesh_vertices),
            mesh_elements=_copy_array(self.mesh_elements),
            stress=_copy_array(self.stress),
            strain=_copy_array(self.strain),
            done=self.done,
            info=copy.deepcopy(self.info),
            user_data=copy.deepcopy(self.user_data),
        )


State = SimulatorState


def _copy_array(value: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if value is None:
        return None
    return np.array(value, copy=True)


def primary_positions(state: State) -> np.ndarray:
    for candidate in (
        state.particle_positions,
        state.mesh_vertices,
        state.positions,
    ):
        if candidate is not None:
            array = np.asarray(candidate, dtype=float)
            if array.ndim == 1:
                return array.reshape(-1, 1)
            return array
    return np.zeros((0, 3), dtype=float)


def rasterize_points(
    points: Optional[np.ndarray],
    render_config: RenderConfig,
    color: Tuple[int, int, int] = (31, 119, 180),
) -> np.ndarray:
    image = np.empty((render_config.height, render_config.width, 3), dtype=np.uint8)
    image[:] = np.asarray(render_config.background_color, dtype=np.uint8)

    if points is None:
        return image

    pts = np.asarray(points, dtype=float)
    if pts.size == 0:
        return image
    if pts.ndim != 2:
        pts = pts.reshape(-1, pts.shape[-1] if pts.ndim > 1 else 1)
    if pts.shape[1] == 1:
        pts = np.concatenate([pts, np.zeros((pts.shape[0], 1))], axis=1)

    planar = pts[:, :2]
    mins = planar.min(axis=0)
    spans = np.maximum(planar.max(axis=0) - mins, 1.0e-8)
    normalized = (planar - mins) / spans
    xs = np.clip((normalized[:, 0] * (render_config.width - 1)).astype(int), 0, render_config.width - 1)
    ys = np.clip((normalized[:, 1] * (render_config.height - 1)).astype(int), 0, render_config.height - 1)

    radius = render_config.point_radius
    pixel_color = np.asarray(color, dtype=np.uint8)
    for x, y in zip(xs, ys):
        y0 = max(0, render_config.height - 1 - y - radius)
        y1 = min(render_config.height, render_config.height - y + radius)
        x0 = max(0, x - radius)
        x1 = min(render_config.width, x + radius + 1)
        image[y0:y1, x0:x1] = pixel_color
    return image


class BaseSimulator(ABC):
    """Abstract interface shared by all physics simulator integrations."""

    def __init__(self, config: SimConfig) -> None:
        self.config = config
        self._state: Optional[State] = None
        self._rng = np.random.default_rng(config.seed)

    @abstractmethod
    def reset(self, initial_state: Any = None) -> State:
        raise NotImplementedError

    @abstractmethod
    def step(self, action: Any) -> Tuple[State, float, bool, Dict[str, Any]]:
        raise NotImplementedError

    def get_state(self) -> State:
        if self._state is None:
            raise SimulatorError('Simulator has not been reset yet.')
        return self._state.copy()

    def set_state(self, state: Union[State, Mapping[str, Any]]) -> None:
        self._state = self._coerce_state(state)

    @abstractmethod
    def render(self, mode: str = 'rgb_array') -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def get_jacobian(self, state: State, action: Any) -> np.ndarray:
        raise NotImplementedError

    def _coerce_state(self, state: Union[State, Mapping[str, Any]]) -> State:
        if isinstance(state, SimulatorState):
            return state.copy()
        if isinstance(state, Mapping):
            state_fields = {item.name for item in dataclass_fields(SimulatorState)}
            direct_values = {key: copy.deepcopy(value) for key, value in state.items() if key in state_fields}
            extra_values = {key: copy.deepcopy(value) for key, value in state.items() if key not in state_fields}
            coerced = SimulatorState(**direct_values)
            coerced.user_data.update(extra_values)
            return coerced
        raise TypeError(f'Unsupported state type: {type(state)!r}')

    def _commit_state(self, state: Union[State, Mapping[str, Any]]) -> State:
        self._state = self._coerce_state(state)
        return self._state.copy()
