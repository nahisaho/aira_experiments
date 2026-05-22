from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import Tensor

from .base_planner import ActionSequence, BasePlanner, DynamicsModelProtocol, StateType, state_distance


@dataclass(slots=True)
class CEMConfig:
    """Configuration for cross-entropy method planning."""

    horizon: int = 10
    population_size: int = 512
    elite_fraction: float = 0.1
    num_iterations: int = 6
    init_std: float = 1.0
    min_std: float = 1e-2
    action_smoothing: float = 0.25
    smoothness_cost_weight: float = 0.05


@dataclass(slots=True)
class MPPIConfig:
    """Configuration for Model Predictive Path Integral planning."""

    horizon: int = 10
    population_size: int = 1024
    num_iterations: int = 5
    temperature: float = 1.0
    noise_covariance: float | Tensor = 0.5
    action_smoothing: float = 0.25
    smoothness_cost_weight: float = 0.05


class _SamplingPlanner(BasePlanner):
    """Base class for batched GPU trajectory sampling planners."""

    def __init__(
        self,
        dynamics_model: DynamicsModelProtocol,
        action_dim: int,
        *,
        default_horizon: int,
        device: torch.device | str | None = None,
        action_bounds: tuple[float, float] | tuple[Tensor, Tensor] | None = None,
    ) -> None:
        super().__init__(
            action_dim,
            device=device,
            action_bounds=action_bounds,
            default_horizon=default_horizon,
        )
        self.dynamics_model = dynamics_model
        self._previous_plan: ActionSequence | None = None

    def replan(
        self,
        current_state: StateType,
        goal_state: StateType,
        previous_plan: ActionSequence | None,
    ) -> ActionSequence:
        self._previous_plan = previous_plan
        return self.plan(current_state, goal_state, self.default_horizon)

    def _evaluate_trajectories(
        self,
        initial_state: StateType,
        goal_state: StateType,
        actions: Tensor,
        smoothness_weight: float,
    ) -> tuple[Tensor, list[StateType]]:
        population_size, horizon, _ = actions.shape
        state = initial_state
        if isinstance(state, Tensor):
            state = state.to(self.device)
            if state.ndim == 1:
                state = state.unsqueeze(0)
            if state.shape[0] == 1:
                state = state.expand(population_size, *state.shape[1:])
        else:
            state = {
                key: value.to(self.device).unsqueeze(0).expand(population_size, *value.shape)
                if value.ndim == 1
                else value.to(self.device).expand(population_size, *value.shape[1:])
                for key, value in state.items()
            }
        trajectory: list[StateType] = []
        for t in range(horizon):
            state = self.dynamics_model.predict_next_state(state, actions[:, t])
            trajectory.append(state)
        terminal_cost = state_distance(trajectory[-1], goal_state, reduction="none")
        rollout_cost = torch.stack(
            [state_distance(state_t, goal_state, reduction="none") for state_t in trajectory],
            dim=0,
        ).mean(dim=0)
        smoothness_cost = (
            (actions[:, 1:] - actions[:, :-1]).square().mean(dim=(1, 2))
            if horizon > 1
            else torch.zeros(population_size, device=self.device)
        )
        total_cost = rollout_cost + terminal_cost + smoothness_weight * smoothness_cost
        return total_cost, trajectory

    def _apply_smoothing(self, actions: Tensor, smoothing: float) -> Tensor:
        if smoothing <= 0.0 or actions.shape[1] < 3:
            return actions
        padded = torch.cat([actions[:, :1], actions, actions[:, -1:]], dim=1)
        return (1.0 - smoothing) * actions + 0.5 * smoothing * (padded[:, :-2] + padded[:, 2:])

    def _warm_start_mean(self, horizon: int) -> Tensor:
        if self._previous_plan is None or self._previous_plan.actions.numel() == 0:
            return torch.zeros(horizon, self.action_dim, device=self.device)
        previous_actions = self._previous_plan.actions.to(self.device)
        shifted = torch.cat([previous_actions[1:], previous_actions[-1:]], dim=0)
        if shifted.shape[0] >= horizon:
            return shifted[:horizon]
        pad = shifted[-1:].expand(horizon - shifted.shape[0], -1)
        return torch.cat([shifted, pad], dim=0)

    @staticmethod
    def _trajectory_first_state(trajectory: list[StateType]) -> list[StateType]:
        return [
            state[0] if isinstance(state, Tensor) else {key: value[0] for key, value in state.items()}
            for state in trajectory
        ]


class CEMPlanner(_SamplingPlanner):
    """Cross-Entropy Method planner with parallel trajectory evaluation."""

    def __init__(
        self,
        dynamics_model: DynamicsModelProtocol,
        action_dim: int,
        *,
        config: CEMConfig | None = None,
        device: torch.device | str | None = None,
        action_bounds: tuple[float, float] | tuple[Tensor, Tensor] | None = None,
    ) -> None:
        self.config = config or CEMConfig()
        super().__init__(
            dynamics_model,
            action_dim,
            default_horizon=self.config.horizon,
            device=device,
            action_bounds=action_bounds,
        )

    def plan(self, current_state: StateType, goal_state: StateType, horizon: int) -> ActionSequence:
        horizon = horizon or self.default_horizon
        mean = self._warm_start_mean(horizon)
        std = torch.full_like(mean, self.config.init_std)
        best_actions = mean.clone()
        best_cost = float("inf")
        elite_count = max(1, int(self.config.population_size * self.config.elite_fraction))
        iteration_costs: list[float] = []

        for _ in range(self.config.num_iterations):
            samples = mean.unsqueeze(0) + std.unsqueeze(0) * torch.randn(
                self.config.population_size,
                horizon,
                self.action_dim,
                device=self.device,
            )
            samples = self.clamp_actions(self._apply_smoothing(samples, self.config.action_smoothing))
            costs, _ = self._evaluate_trajectories(
                current_state,
                goal_state,
                samples,
                smoothness_weight=self.config.smoothness_cost_weight,
            )
            elite_indices = torch.topk(costs, k=elite_count, largest=False).indices
            elites = samples[elite_indices]
            mean = elites.mean(dim=0)
            std = elites.std(dim=0, unbiased=False).clamp_min(self.config.min_std)
            iteration_costs.append(float(costs.min().item()))
            if costs.min().item() < best_cost:
                best_cost = float(costs.min().item())
                best_actions = samples[costs.argmin()].detach().clone()

        _, trajectory = self._evaluate_trajectories(
            current_state,
            goal_state,
            best_actions.unsqueeze(0),
            smoothness_weight=self.config.smoothness_cost_weight,
        )
        plan = ActionSequence(
            actions=best_actions,
            predicted_states=self._trajectory_first_state(trajectory),
            cost=best_cost,
            metadata={"planner": "cem", "iteration_costs": iteration_costs},
        )
        self._previous_plan = plan
        return plan


class MPPIPlanner(_SamplingPlanner):
    """Model Predictive Path Integral planner with parallel rollout evaluation."""

    def __init__(
        self,
        dynamics_model: DynamicsModelProtocol,
        action_dim: int,
        *,
        config: MPPIConfig | None = None,
        device: torch.device | str | None = None,
        action_bounds: tuple[float, float] | tuple[Tensor, Tensor] | None = None,
    ) -> None:
        self.config = config or MPPIConfig()
        super().__init__(
            dynamics_model,
            action_dim,
            default_horizon=self.config.horizon,
            device=device,
            action_bounds=action_bounds,
        )

    def plan(self, current_state: StateType, goal_state: StateType, horizon: int) -> ActionSequence:
        horizon = horizon or self.default_horizon
        mean = self._warm_start_mean(horizon)
        best_actions = mean.clone()
        best_cost = float("inf")
        iteration_costs: list[float] = []

        covariance = self._covariance_tensor().to(self.device)
        chol = torch.linalg.cholesky(covariance)

        for _ in range(self.config.num_iterations):
            base_noise = torch.randn(self.config.population_size, horizon, self.action_dim, device=self.device)
            noise = torch.einsum("pha,ab->phb", base_noise, chol)
            samples = self.clamp_actions(self._apply_smoothing(mean.unsqueeze(0) + noise, self.config.action_smoothing))
            costs, _ = self._evaluate_trajectories(
                current_state,
                goal_state,
                samples,
                smoothness_weight=self.config.smoothness_cost_weight,
            )
            normalized = costs - costs.min()
            weights = torch.softmax(-normalized / max(self.config.temperature, 1e-6), dim=0)
            mean = torch.sum(weights[:, None, None] * samples, dim=0)
            iteration_costs.append(float(costs.min().item()))
            if costs.min().item() < best_cost:
                best_cost = float(costs.min().item())
                best_actions = samples[costs.argmin()].detach().clone()

        _, trajectory = self._evaluate_trajectories(
            current_state,
            goal_state,
            best_actions.unsqueeze(0),
            smoothness_weight=self.config.smoothness_cost_weight,
        )
        plan = ActionSequence(
            actions=best_actions,
            predicted_states=self._trajectory_first_state(trajectory),
            cost=best_cost,
            metadata={"planner": "mppi", "iteration_costs": iteration_costs},
        )
        self._previous_plan = plan
        return plan

    def _covariance_tensor(self) -> Tensor:
        covariance = self.config.noise_covariance
        if isinstance(covariance, Tensor):
            if covariance.ndim == 1:
                return torch.diag(covariance)
            return covariance
        return torch.eye(self.action_dim) * float(covariance)


__all__ = ["CEMConfig", "CEMPlanner", "MPPIConfig", "MPPIPlanner"]
