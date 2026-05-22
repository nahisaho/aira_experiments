from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Mapping, Protocol, TypeAlias

import torch
from torch import Tensor

StateType: TypeAlias = Tensor | Mapping[str, Tensor]
ActionBounds: TypeAlias = tuple[float, float] | tuple[Tensor, Tensor]


class DynamicsModelProtocol(Protocol):
    """Protocol for learned dynamics models used by planners."""

    def predict_next_state(self, current_state: StateType, action: Tensor) -> StateType:
        """Predict the next state for a batch or single state/action pair."""


@dataclass(slots=True)
class ActionSequence:
    """Container for planned actions and auxiliary planning metadata."""

    actions: Tensor
    predicted_states: list[StateType] = field(default_factory=list)
    cost: float | Tensor | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def to(self, device: torch.device | str) -> "ActionSequence":
        """Move all tensor payloads to the requested device."""
        predicted_states = [move_state_to_device(state, device) for state in self.predicted_states]
        metadata = {
            key: move_state_to_device(value, device) if isinstance(value, (Tensor, dict)) else value
            for key, value in self.metadata.items()
        }
        return ActionSequence(
            actions=self.actions.to(device),
            predicted_states=predicted_states,
            cost=self.cost.to(device) if isinstance(self.cost, Tensor) else self.cost,
            metadata=metadata,
        )


class BasePlanner(ABC):
    """Abstract base class for manipulation sequence planners."""

    def __init__(
        self,
        action_dim: int,
        *,
        device: torch.device | str | None = None,
        action_bounds: ActionBounds | None = None,
        default_horizon: int = 10,
    ) -> None:
        self.action_dim = action_dim
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.action_bounds = _normalize_action_bounds(action_bounds, action_dim, self.device)
        self.default_horizon = default_horizon

    @abstractmethod
    def plan(self, current_state: StateType, goal_state: StateType, horizon: int) -> ActionSequence:
        """Generate an action sequence from the current state toward the goal state."""

    @abstractmethod
    def replan(
        self,
        current_state: StateType,
        goal_state: StateType,
        previous_plan: ActionSequence | None,
    ) -> ActionSequence:
        """Update a previous plan after observing a new state estimate."""

    def clamp_actions(self, actions: Tensor) -> Tensor:
        """Clamp actions to planner bounds when bounds are configured."""
        if self.action_bounds is None:
            return actions
        lower, upper = self.action_bounds
        return torch.max(torch.min(actions, upper), lower)

    def action_smoothness_cost(self, actions: Tensor) -> Tensor:
        """Quadratic penalty on action changes across time."""
        if actions.shape[0] < 2:
            return torch.zeros((), device=actions.device, dtype=actions.dtype)
        deltas = actions[1:] - actions[:-1]
        return deltas.square().mean()

    def constraint_cost(self, actions: Tensor) -> Tensor:
        """Soft penalty for bound violations; zero when actions remain feasible."""
        if self.action_bounds is None:
            return torch.zeros((), device=actions.device, dtype=actions.dtype)
        lower, upper = self.action_bounds
        lower_violation = torch.relu(lower - actions)
        upper_violation = torch.relu(actions - upper)
        return (lower_violation.square() + upper_violation.square()).mean()


def _normalize_action_bounds(
    action_bounds: ActionBounds | None,
    action_dim: int,
    device: torch.device,
) -> tuple[Tensor, Tensor] | None:
    if action_bounds is None:
        return None
    lower, upper = action_bounds
    if isinstance(lower, Tensor):
        lower_tensor = lower.to(device=device, dtype=torch.float32).reshape(1, -1)
    else:
        lower_tensor = torch.full((1, action_dim), float(lower), device=device)
    if isinstance(upper, Tensor):
        upper_tensor = upper.to(device=device, dtype=torch.float32).reshape(1, -1)
    else:
        upper_tensor = torch.full((1, action_dim), float(upper), device=device)
    if lower_tensor.shape[-1] != action_dim or upper_tensor.shape[-1] != action_dim:
        raise ValueError("Action bounds must match action_dim.")
    return lower_tensor, upper_tensor


def move_state_to_device(value: Any, device: torch.device | str) -> Any:
    """Recursively move tensors contained in a state-like object to a device."""
    if isinstance(value, Tensor):
        return value.to(device)
    if isinstance(value, Mapping):
        return {key: move_state_to_device(item, device) for key, item in value.items()}
    return value


def clone_state(state: StateType) -> StateType:
    """Clone a tensor or mapping-based state representation."""
    if isinstance(state, Tensor):
        return state.clone()
    return {key: value.clone() for key, value in state.items()}


def flatten_state(state: StateType) -> Tensor:
    """Flatten a tensor or dictionary state into a single feature vector."""
    if isinstance(state, Tensor):
        return state.reshape(-1)
    flattened = [value.reshape(-1) for _, value in sorted(state.items(), key=lambda item: item[0])]
    if not flattened:
        raise ValueError("State mapping must contain at least one tensor.")
    return torch.cat(flattened, dim=0)


def batched_flatten_state(state: StateType) -> Tensor:
    """Flatten a batched or unbatched state into shape [B, D]."""
    if isinstance(state, Tensor):
        return state.reshape(1, -1) if state.ndim == 1 else state.reshape(state.shape[0], -1)
    components = []
    batch_size: int | None = None
    for _, value in sorted(state.items(), key=lambda item: item[0]):
        if value.ndim == 1:
            value = value.unsqueeze(0)
        current_batch = value.shape[0]
        if batch_size is None:
            batch_size = current_batch
        elif current_batch != batch_size:
            raise ValueError("All state tensors must share the same batch dimension.")
        components.append(value.reshape(current_batch, -1))
    if not components:
        raise ValueError("State mapping must contain at least one tensor.")
    return torch.cat(components, dim=-1)


def state_distance(current_state: StateType, goal_state: StateType, reduction: str = "mean") -> Tensor:
    """Squared Euclidean state-matching cost."""
    current = batched_flatten_state(current_state)
    goal = batched_flatten_state(goal_state)
    if goal.shape[0] == 1 and current.shape[0] > 1:
        goal = goal.expand(current.shape[0], -1)
    if current.shape != goal.shape:
        raise ValueError("Current and goal states must flatten to the same shape.")
    distances = (current - goal).square().mean(dim=-1)
    if reduction == "none":
        return distances
    if reduction == "sum":
        return distances.sum()
    return distances.mean()


def interpolate_state(start_state: StateType, end_state: StateType, alpha: float) -> StateType:
    """Linear interpolation between two compatible states."""
    if isinstance(start_state, Tensor) and isinstance(end_state, Tensor):
        return torch.lerp(start_state, end_state, alpha)
    if isinstance(start_state, Mapping) and isinstance(end_state, Mapping):
        shared_keys = set(start_state) & set(end_state)
        if shared_keys != set(start_state) or shared_keys != set(end_state):
            raise ValueError("State mappings must share the same keys for interpolation.")
        return {
            key: torch.lerp(start_state[key], end_state[key], alpha)
            for key in sorted(shared_keys)
        }
    raise TypeError("State types must match for interpolation.")


__all__ = [
    "ActionSequence",
    "ActionBounds",
    "BasePlanner",
    "DynamicsModelProtocol",
    "StateType",
    "batched_flatten_state",
    "clone_state",
    "flatten_state",
    "interpolate_state",
    "move_state_to_device",
    "state_distance",
]
