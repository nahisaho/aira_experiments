from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import torch
from torch import Tensor

from .base_planner import ActionSequence, BasePlanner, DynamicsModelProtocol, StateType, clone_state, state_distance


ConstraintFn = Callable[[list[StateType], Tensor], Tensor]


@dataclass(slots=True)
class MPCConfig:
    """Configuration for gradient-based receding-horizon MPC."""

    horizon: int = 10
    optimization_steps: int = 100
    learning_rate: float = 5e-2
    state_cost_weight: float = 1.0
    terminal_cost_weight: float = 5.0
    smoothness_cost_weight: float = 0.05
    constraint_cost_weight: float = 10.0


class MPCPlanner(BasePlanner):
    """Gradient-based MPC planner using a learned differentiable dynamics model."""

    def __init__(
        self,
        dynamics_model: DynamicsModelProtocol,
        action_dim: int,
        *,
        config: MPCConfig | None = None,
        constraint_fn: ConstraintFn | None = None,
        device: torch.device | str | None = None,
        action_bounds: tuple[float, float] | tuple[Tensor, Tensor] | None = None,
    ) -> None:
        self.config = config or MPCConfig()
        super().__init__(
            action_dim,
            device=device,
            action_bounds=action_bounds,
            default_horizon=self.config.horizon,
        )
        self.dynamics_model = dynamics_model
        self.constraint_fn = constraint_fn
        self._previous_plan: ActionSequence | None = None

    def plan(self, current_state: StateType, goal_state: StateType, horizon: int) -> ActionSequence:
        horizon = horizon or self.default_horizon
        initial_actions = self._initialize_action_guess(horizon, self._previous_plan)
        action_parameters = torch.nn.Parameter(initial_actions)
        optimizer = torch.optim.Adam([action_parameters], lr=self.config.learning_rate)
        best_cost = float("inf")
        best_actions = initial_actions.detach().clone()
        optimization_history: list[float] = []

        for _ in range(self.config.optimization_steps):
            optimizer.zero_grad(set_to_none=True)
            clamped_actions = self.clamp_actions(action_parameters)
            predicted_states = self._rollout(current_state, clamped_actions)
            cost = self._compute_cost(predicted_states, goal_state, clamped_actions)
            cost.backward()
            optimizer.step()
            optimization_history.append(float(cost.detach().item()))
            with torch.no_grad():
                action_parameters.copy_(self.clamp_actions(action_parameters))
                if cost.item() < best_cost:
                    best_cost = cost.item()
                    best_actions = clamped_actions.detach().clone()

        predicted_states = self._rollout(current_state, best_actions)
        plan = ActionSequence(
            actions=best_actions.detach(),
            predicted_states=[clone_state(state) for state in predicted_states],
            cost=best_cost,
            metadata={"optimizer": "adam", "optimization_history": optimization_history},
        )
        self._previous_plan = plan
        return plan

    def replan(
        self,
        current_state: StateType,
        goal_state: StateType,
        previous_plan: ActionSequence | None,
    ) -> ActionSequence:
        self._previous_plan = previous_plan
        return self.plan(current_state, goal_state, self.default_horizon)

    def _initialize_action_guess(self, horizon: int, previous_plan: ActionSequence | None) -> Tensor:
        if previous_plan is not None and previous_plan.actions.numel() > 0:
            previous_actions = previous_plan.actions.to(self.device)
            shifted = torch.cat([previous_actions[1:], previous_actions[-1:]], dim=0)
            if shifted.shape[0] >= horizon:
                return shifted[:horizon].clone()
            pad = shifted[-1:].expand(horizon - shifted.shape[0], -1)
            return torch.cat([shifted, pad], dim=0)
        return torch.zeros(horizon, self.action_dim, device=self.device)

    def _rollout(self, initial_state: StateType, actions: Tensor) -> list[StateType]:
        state = clone_state(initial_state)
        trajectory: list[StateType] = []
        for action in actions:
            state = self.dynamics_model.predict_next_state(state, action)
            trajectory.append(state)
        return trajectory

    def _compute_cost(self, predicted_states: list[StateType], goal_state: StateType, actions: Tensor) -> Tensor:
        state_cost = torch.stack([state_distance(state, goal_state) for state in predicted_states]).mean()
        terminal_cost = state_distance(predicted_states[-1], goal_state)
        smoothness_cost = self.action_smoothness_cost(actions)
        constraint_cost = self.constraint_cost(actions)
        if self.constraint_fn is not None:
            constraint_cost = constraint_cost + self.constraint_fn(predicted_states, actions)
        return (
            self.config.state_cost_weight * state_cost
            + self.config.terminal_cost_weight * terminal_cost
            + self.config.smoothness_cost_weight * smoothness_cost
            + self.config.constraint_cost_weight * constraint_cost
        )


__all__ = ["MPCConfig", "MPCPlanner"]
