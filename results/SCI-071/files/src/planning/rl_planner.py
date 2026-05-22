from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import torch
from torch import Tensor, nn
from torch.distributions import Normal
from torch.nn import functional as F

from .base_planner import ActionSequence, BasePlanner, DynamicsModelProtocol, StateType, batched_flatten_state


LOG_STD_MIN = -5.0
LOG_STD_MAX = 2.0


@dataclass(slots=True)
class RewardShapingConfig:
    """Weights for deformable-object reward shaping terms."""

    chamfer_weight: float = 1.0
    iou_weight: float = 1.0
    smoothness_weight: float = 0.05


@dataclass(slots=True)
class SACConfig:
    """Soft Actor-Critic training configuration."""

    gamma: float = 0.99
    tau: float = 0.005
    actor_lr: float = 3e-4
    critic_lr: float = 3e-4
    alpha_lr: float = 3e-4
    automatic_entropy_tuning: bool = True
    target_entropy: float | None = None


@dataclass(slots=True)
class PPOConfig:
    """Proximal Policy Optimization training configuration."""

    gamma: float = 0.99
    gae_lambda: float = 0.95
    clip_ratio: float = 0.2
    actor_lr: float = 3e-4
    critic_lr: float = 3e-4
    value_coef: float = 0.5
    entropy_coef: float = 0.0
    update_epochs: int = 10
    minibatch_size: int = 128


@dataclass(slots=True)
class CurriculumScheduler:
    """Progressive-difficulty scheduler for curriculum learning."""

    difficulty: float = 0.0
    increment: float = 0.1
    success_threshold: float = 0.8
    max_difficulty: float = 1.0

    def update(self, success_rate: float) -> float:
        if success_rate >= self.success_threshold:
            self.difficulty = min(self.max_difficulty, self.difficulty + self.increment)
        return self.difficulty


class RLAgentProtocol(Protocol):
    """Common interface for inference-time RL planners."""

    def act(self, state: Tensor, deterministic: bool = True) -> Tensor:
        """Select an action for the provided state."""


class GaussianActor(nn.Module):
    """Gaussian policy with tanh squashing."""

    def __init__(self, state_dim: int, action_dim: int, hidden_dim: int = 256) -> None:
        super().__init__()
        self.backbone = nn.Sequential(
            nn.Linear(state_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )
        self.mu = nn.Linear(hidden_dim, action_dim)
        self.log_std = nn.Linear(hidden_dim, action_dim)

    def forward(self, state: Tensor) -> tuple[Tensor, Tensor]:
        features = self.backbone(state)
        mu = self.mu(features)
        log_std = self.log_std(features).clamp(LOG_STD_MIN, LOG_STD_MAX)
        return mu, log_std

    def sample(self, state: Tensor) -> tuple[Tensor, Tensor]:
        mu, log_std = self(state)
        std = log_std.exp()
        distribution = Normal(mu, std)
        raw_action = distribution.rsample()
        squashed_action = torch.tanh(raw_action)
        log_prob = distribution.log_prob(raw_action) - torch.log1p(-squashed_action.square() + 1e-6)
        return squashed_action, log_prob.sum(dim=-1, keepdim=True)

    def deterministic(self, state: Tensor) -> Tensor:
        mu, _ = self(state)
        return torch.tanh(mu)


class QNetwork(nn.Module):
    """Action-value network for actor-critic methods."""

    def __init__(self, state_dim: int, action_dim: int, hidden_dim: int = 256) -> None:
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(state_dim + action_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, state: Tensor, action: Tensor) -> Tensor:
        return self.network(torch.cat([state, action], dim=-1))


class ValueNetwork(nn.Module):
    """State-value network used by PPO."""

    def __init__(self, state_dim: int, hidden_dim: int = 256) -> None:
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(state_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, state: Tensor) -> Tensor:
        return self.network(state)


class SACAgent(RLAgentProtocol):
    """Soft Actor-Critic agent with automatic entropy tuning."""

    def __init__(self, state_dim: int, action_dim: int, *, config: SACConfig | None = None, device: torch.device | str | None = None) -> None:
        self.config = config or SACConfig()
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.actor = GaussianActor(state_dim, action_dim).to(self.device)
        self.q1 = QNetwork(state_dim, action_dim).to(self.device)
        self.q2 = QNetwork(state_dim, action_dim).to(self.device)
        self.q1_target = QNetwork(state_dim, action_dim).to(self.device)
        self.q2_target = QNetwork(state_dim, action_dim).to(self.device)
        self.q1_target.load_state_dict(self.q1.state_dict())
        self.q2_target.load_state_dict(self.q2.state_dict())

        self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=self.config.actor_lr)
        self.q1_optimizer = torch.optim.Adam(self.q1.parameters(), lr=self.config.critic_lr)
        self.q2_optimizer = torch.optim.Adam(self.q2.parameters(), lr=self.config.critic_lr)

        target_entropy = self.config.target_entropy or -float(action_dim)
        self.target_entropy = target_entropy
        self.log_alpha = torch.zeros(1, device=self.device, requires_grad=True)
        self.alpha_optimizer = torch.optim.Adam([self.log_alpha], lr=self.config.alpha_lr)

    @property
    def alpha(self) -> Tensor:
        return self.log_alpha.exp()

    def act(self, state: Tensor, deterministic: bool = True) -> Tensor:
        with torch.no_grad():
            state = state.to(self.device)
            if deterministic:
                return self.actor.deterministic(state)
            action, _ = self.actor.sample(state)
            return action

    def update(self, batch: dict[str, Tensor]) -> dict[str, float]:
        states = batch["states"].to(self.device)
        actions = batch["actions"].to(self.device)
        rewards = batch["rewards"].to(self.device)
        next_states = batch["next_states"].to(self.device)
        dones = batch["dones"].to(self.device)

        with torch.no_grad():
            next_actions, next_log_probs = self.actor.sample(next_states)
            next_q = torch.min(
                self.q1_target(next_states, next_actions),
                self.q2_target(next_states, next_actions),
            ) - self.alpha.detach() * next_log_probs
            targets = rewards + self.config.gamma * (1.0 - dones) * next_q

        q1_loss = F.mse_loss(self.q1(states, actions), targets)
        q2_loss = F.mse_loss(self.q2(states, actions), targets)
        self.q1_optimizer.zero_grad(set_to_none=True)
        q1_loss.backward()
        self.q1_optimizer.step()
        self.q2_optimizer.zero_grad(set_to_none=True)
        q2_loss.backward()
        self.q2_optimizer.step()

        sampled_actions, log_probs = self.actor.sample(states)
        q_estimate = torch.min(self.q1(states, sampled_actions), self.q2(states, sampled_actions))
        actor_loss = (self.alpha.detach() * log_probs - q_estimate).mean()
        self.actor_optimizer.zero_grad(set_to_none=True)
        actor_loss.backward()
        self.actor_optimizer.step()

        alpha_loss = torch.tensor(0.0, device=self.device)
        if self.config.automatic_entropy_tuning:
            alpha_loss = -(self.log_alpha * (log_probs.detach() + self.target_entropy)).mean()
            self.alpha_optimizer.zero_grad(set_to_none=True)
            alpha_loss.backward()
            self.alpha_optimizer.step()

        self._soft_update(self.q1_target, self.q1)
        self._soft_update(self.q2_target, self.q2)
        return {
            "actor_loss": float(actor_loss.item()),
            "q1_loss": float(q1_loss.item()),
            "q2_loss": float(q2_loss.item()),
            "alpha_loss": float(alpha_loss.item()),
            "alpha": float(self.alpha.item()),
        }

    def _soft_update(self, target: nn.Module, source: nn.Module) -> None:
        for target_param, source_param in zip(target.parameters(), source.parameters()):
            target_param.data.mul_(1.0 - self.config.tau).add_(source_param.data, alpha=self.config.tau)


class PPOAgent(RLAgentProtocol):
    """PPO agent with Generalized Advantage Estimation."""

    def __init__(self, state_dim: int, action_dim: int, *, config: PPOConfig | None = None, device: torch.device | str | None = None) -> None:
        self.config = config or PPOConfig()
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.actor = GaussianActor(state_dim, action_dim).to(self.device)
        self.critic = ValueNetwork(state_dim).to(self.device)
        self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=self.config.actor_lr)
        self.critic_optimizer = torch.optim.Adam(self.critic.parameters(), lr=self.config.critic_lr)

    def act(self, state: Tensor, deterministic: bool = True) -> Tensor:
        with torch.no_grad():
            state = state.to(self.device)
            if deterministic:
                return self.actor.deterministic(state)
            action, _ = self.actor.sample(state)
            return action

    def compute_gae(self, rewards: Tensor, values: Tensor, dones: Tensor, next_value: Tensor) -> tuple[Tensor, Tensor]:
        advantages = torch.zeros_like(rewards)
        gae = torch.zeros(1, device=rewards.device)
        for t in reversed(range(rewards.shape[0])):
            mask = 1.0 - dones[t]
            delta = rewards[t] + self.config.gamma * next_value * mask - values[t]
            gae = delta + self.config.gamma * self.config.gae_lambda * mask * gae
            advantages[t] = gae
            next_value = values[t]
        returns = advantages + values
        return advantages, returns

    def update(self, batch: dict[str, Tensor]) -> dict[str, float]:
        states = batch["states"].to(self.device)
        actions = batch["actions"].to(self.device)
        old_log_probs = batch["log_probs"].to(self.device)
        advantages = batch["advantages"].to(self.device)
        returns = batch["returns"].to(self.device)

        advantages = (advantages - advantages.mean()) / (advantages.std(unbiased=False) + 1e-6)
        dataset_size = states.shape[0]
        actor_loss_value = 0.0
        critic_loss_value = 0.0

        for _ in range(self.config.update_epochs):
            permutation = torch.randperm(dataset_size, device=self.device)
            for start in range(0, dataset_size, self.config.minibatch_size):
                indices = permutation[start : start + self.config.minibatch_size]
                minibatch_states = states[indices]
                minibatch_actions = actions[indices]
                minibatch_old_log_probs = old_log_probs[indices]
                minibatch_advantages = advantages[indices]
                minibatch_returns = returns[indices]

                mu, log_std = self.actor(minibatch_states)
                distribution = Normal(mu, log_std.exp())
                raw_actions = torch.atanh(minibatch_actions.clamp(-0.999, 0.999))
                log_probs = distribution.log_prob(raw_actions).sum(dim=-1, keepdim=True)
                ratio = (log_probs - minibatch_old_log_probs).exp()
                unclipped = ratio * minibatch_advantages
                clipped = torch.clamp(ratio, 1.0 - self.config.clip_ratio, 1.0 + self.config.clip_ratio) * minibatch_advantages
                actor_loss = -(torch.min(unclipped, clipped)).mean()

                values = self.critic(minibatch_states)
                critic_loss = F.mse_loss(values, minibatch_returns)
                entropy = distribution.entropy().sum(dim=-1).mean()
                total_loss = actor_loss + self.config.value_coef * critic_loss - self.config.entropy_coef * entropy

                self.actor_optimizer.zero_grad(set_to_none=True)
                self.critic_optimizer.zero_grad(set_to_none=True)
                total_loss.backward()
                self.actor_optimizer.step()
                self.critic_optimizer.step()

                actor_loss_value = float(actor_loss.item())
                critic_loss_value = float(critic_loss.item())

        return {"actor_loss": actor_loss_value, "critic_loss": critic_loss_value}


class RLPlanner(BasePlanner):
    """Planner wrapper around trained SAC or PPO policies with reward shaping."""

    def __init__(
        self,
        agent: RLAgentProtocol,
        action_dim: int,
        *,
        world_model: DynamicsModelProtocol | None = None,
        reward_shaping: RewardShapingConfig | None = None,
        curriculum: CurriculumScheduler | None = None,
        device: torch.device | str | None = None,
        action_bounds: tuple[float, float] | tuple[Tensor, Tensor] | None = None,
        default_horizon: int = 10,
    ) -> None:
        super().__init__(action_dim, device=device, action_bounds=action_bounds, default_horizon=default_horizon)
        self.agent = agent
        self.world_model = world_model
        self.reward_shaping = reward_shaping or RewardShapingConfig()
        self.curriculum = curriculum or CurriculumScheduler()

    def plan(self, current_state: StateType, goal_state: StateType, horizon: int) -> ActionSequence:
        horizon = horizon or self.default_horizon
        state = current_state
        actions = []
        predicted_states: list[StateType] = []
        rewards = []
        previous_action = torch.zeros(1, self.action_dim, device=self.device)

        for _ in range(horizon):
            encoded_state = batched_flatten_state(state).to(self.device)
            action = self.agent.act(encoded_state, deterministic=True)
            action = self.clamp_actions(action)
            action_to_apply = action.squeeze(0) if action.ndim == 2 else action
            actions.append(action_to_apply)
            if self.world_model is not None:
                state = self.world_model.predict_next_state(state, action_to_apply)
            predicted_states.append(state)
            rewards.append(self.compute_shaped_reward(state, goal_state, action, previous_action))
            previous_action = action

        stacked_actions = torch.stack(actions, dim=0) if actions else torch.zeros(0, self.action_dim, device=self.device)
        total_reward = torch.stack(rewards).sum().item() if rewards else 0.0
        return ActionSequence(
            actions=stacked_actions,
            predicted_states=predicted_states,
            cost=-total_reward,
            metadata={"planner": self.agent.__class__.__name__.lower(), "curriculum_difficulty": self.curriculum.difficulty},
        )

    def replan(
        self,
        current_state: StateType,
        goal_state: StateType,
        previous_plan: ActionSequence | None,
    ) -> ActionSequence:
        return self.plan(current_state, goal_state, self.default_horizon)

    def compute_shaped_reward(
        self,
        state: StateType,
        goal_state: StateType,
        action: Tensor,
        previous_action: Tensor | None = None,
    ) -> Tensor:
        current_points = _state_as_points(state).to(self.device)
        goal_points = _state_as_points(goal_state).to(self.device)
        chamfer = chamfer_distance(current_points, goal_points)
        iou = intersection_over_union(current_points, goal_points)
        smoothness_penalty = torch.zeros((), device=self.device)
        if previous_action is not None:
            smoothness_penalty = (action - previous_action).square().mean()
        return (
            -self.reward_shaping.chamfer_weight * chamfer
            + self.reward_shaping.iou_weight * iou
            - self.reward_shaping.smoothness_weight * smoothness_penalty
        )


def chamfer_distance(points_a: Tensor, points_b: Tensor) -> Tensor:
    """Chamfer distance between two point clouds."""
    distances = torch.cdist(points_a, points_b)
    forward = distances.min(dim=-1).values.mean()
    backward = distances.min(dim=-2).values.mean()
    return forward + backward


def intersection_over_union(points_a: Tensor, points_b: Tensor, threshold: float = 0.05) -> Tensor:
    """Approximate IoU by thresholded point correspondences."""
    distances = torch.cdist(points_a, points_b)
    matches_a = (distances.min(dim=-1).values < threshold).float()
    matches_b = (distances.min(dim=-2).values < threshold).float()
    intersection = 0.5 * (matches_a.sum() + matches_b.sum())
    union = points_a.shape[0] + points_b.shape[0] - intersection
    return intersection / union.clamp_min(1.0)


def _state_as_points(state: StateType) -> Tensor:
    if isinstance(state, Tensor):
        tensor_state = state
    else:
        tensor_state = state.get("positions")
        if tensor_state is None:
            tensor_state = batched_flatten_state(state).squeeze(0)
    if tensor_state.ndim == 1:
        tensor_state = tensor_state.reshape(-1, 1)
    if tensor_state.ndim > 2:
        tensor_state = tensor_state.reshape(-1, tensor_state.shape[-1])
    return tensor_state


__all__ = [
    "CurriculumScheduler",
    "PPOAgent",
    "PPOConfig",
    "RLAgentProtocol",
    "RLPlanner",
    "RewardShapingConfig",
    "SACAgent",
    "SACConfig",
    "chamfer_distance",
    "intersection_over_union",
]
