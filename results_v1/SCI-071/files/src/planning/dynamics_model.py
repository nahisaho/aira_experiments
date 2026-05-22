from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Protocol, Sequence

import torch
from torch import Tensor, nn
from torch.nn import functional as F
from torch.utils.data import DataLoader, Dataset

from .base_planner import StateType, batched_flatten_state


class SimulatorProtocol(Protocol):
    """Minimal simulator interface for training data collection."""

    def reset(self) -> Tensor:
        """Reset the simulator and return the initial state."""

    def sample_action(self, state: Tensor) -> Tensor:
        """Sample an exploratory action for the provided state."""

    def step(self, action: Tensor) -> tuple[Tensor, float, bool, dict[str, Any]]:
        """Advance the simulator with the provided action."""


class BaseDynamicsModel(nn.Module):
    """Base class for learned dynamics models."""

    def predict_next_state(self, current_state: StateType, action: Tensor) -> Tensor:
        """Predict the next state tensor."""
        state_tensor = batched_flatten_state(current_state)
        action_tensor = action.reshape(1, -1) if action.ndim == 1 else action
        next_state = self.forward(state_tensor, action_tensor)
        return next_state.squeeze(0) if state_tensor.shape[0] == 1 else next_state

    def rollout(self, initial_state: StateType, actions: Tensor) -> Tensor:
        """Roll out the model over an action horizon."""
        state = batched_flatten_state(initial_state)
        if actions.ndim == 2:
            actions = actions.unsqueeze(0)
        if state.shape[0] == 1 and actions.shape[0] > 1:
            state = state.expand(actions.shape[0], -1)
        trajectory = []
        for t in range(actions.shape[1]):
            state = self.forward(state, actions[:, t])
            trajectory.append(state)
        return torch.stack(trajectory, dim=1)


class MLPDynamicsModel(BaseDynamicsModel):
    """MLP dynamics model for latent-space transition prediction."""

    def __init__(
        self,
        state_dim: int,
        action_dim: int,
        *,
        hidden_dims: Sequence[int] = (256, 256),
        predict_delta: bool = True,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        self.state_dim = state_dim
        self.action_dim = action_dim
        self.predict_delta = predict_delta

        layers: list[nn.Module] = []
        input_dim = state_dim + action_dim
        for hidden_dim in hidden_dims:
            layers.extend(
                [
                    nn.Linear(input_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.SiLU(),
                    nn.Dropout(dropout),
                ]
            )
            input_dim = hidden_dim
        layers.append(nn.Linear(input_dim, state_dim))
        self.network = nn.Sequential(*layers)

    def forward(self, current_state: Tensor, action: Tensor) -> Tensor:
        state = current_state.reshape(current_state.shape[0], -1)
        action = action.reshape(action.shape[0], -1)
        model_input = torch.cat([state, action], dim=-1)
        prediction = self.network(model_input)
        return state + prediction if self.predict_delta else prediction


class MessagePassingLayer(nn.Module):
    """Simple message-passing block for graph-based dynamics."""

    def __init__(self, node_dim: int, edge_dim: int, action_dim: int) -> None:
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * node_dim + edge_dim + action_dim, node_dim),
            nn.SiLU(),
            nn.Linear(node_dim, node_dim),
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * node_dim, node_dim),
            nn.SiLU(),
            nn.Linear(node_dim, node_dim),
        )

    def forward(self, node_features: Tensor, edge_index: Tensor, edge_attr: Tensor, action_embed: Tensor) -> Tensor:
        src, dst = edge_index
        source_features = node_features[:, src]
        target_features = node_features[:, dst]
        repeated_action = action_embed.unsqueeze(1).expand(-1, source_features.shape[1], -1)
        messages = self.edge_mlp(torch.cat([source_features, target_features, edge_attr, repeated_action], dim=-1))

        aggregated = torch.zeros_like(node_features)
        for batch_idx in range(node_features.shape[0]):
            aggregated[batch_idx].index_add_(0, dst, messages[batch_idx])
        return self.node_mlp(torch.cat([node_features, aggregated], dim=-1)) + node_features


class GNNDynamicsModel(BaseDynamicsModel):
    """Message-passing dynamics model for particle or mesh states."""

    def __init__(
        self,
        node_dim: int,
        action_dim: int,
        *,
        latent_dim: int = 128,
        num_layers: int = 3,
        edge_dim: int | None = None,
    ) -> None:
        super().__init__()
        self.node_dim = node_dim
        self.action_dim = action_dim
        self.latent_dim = latent_dim
        self.edge_dim = edge_dim or node_dim

        self.node_encoder = nn.Linear(node_dim, latent_dim)
        self.action_encoder = nn.Linear(action_dim, latent_dim)
        self.layers = nn.ModuleList(
            [MessagePassingLayer(latent_dim, self.edge_dim, latent_dim) for _ in range(num_layers)]
        )
        self.node_decoder = nn.Sequential(
            nn.Linear(latent_dim, latent_dim),
            nn.SiLU(),
            nn.Linear(latent_dim, node_dim),
        )

    def forward(self, current_state: Tensor, action: Tensor) -> Tensor:
        if current_state.ndim != 3:
            raise ValueError("GNNDynamicsModel expects state shape [batch, num_nodes, node_dim].")
        batch_size, num_nodes, _ = current_state.shape
        edge_index = _fully_connected_edge_index(num_nodes, current_state.device)
        edge_attr = _edge_features(current_state, edge_index)

        node_features = self.node_encoder(current_state)
        action_embed = self.action_encoder(action.reshape(batch_size, -1))
        for layer in self.layers:
            node_features = layer(node_features, edge_index, edge_attr, action_embed)
        delta = self.node_decoder(node_features)
        return current_state + delta

    def predict_next_state(self, current_state: StateType, action: Tensor) -> Tensor:
        if not isinstance(current_state, Tensor):
            raise TypeError("GNNDynamicsModel requires tensor states.")
        state = current_state.unsqueeze(0) if current_state.ndim == 2 else current_state
        action_tensor = action.unsqueeze(0) if action.ndim == 1 else action
        next_state = self.forward(state, action_tensor)
        return next_state.squeeze(0) if current_state.ndim == 2 else next_state


class EnsembleDynamicsModel(BaseDynamicsModel):
    """Ensemble model providing mean predictions and uncertainty estimates."""

    def __init__(self, members: Sequence[BaseDynamicsModel]) -> None:
        super().__init__()
        if not members:
            raise ValueError("EnsembleDynamicsModel requires at least one member.")
        self.members = nn.ModuleList(members)

    def forward(self, current_state: Tensor, action: Tensor) -> Tensor:
        mean, _ = self.predict_with_uncertainty(current_state, action)
        return mean

    def predict_with_uncertainty(self, current_state: Tensor, action: Tensor) -> tuple[Tensor, Tensor]:
        predictions = torch.stack([member(current_state, action) for member in self.members], dim=0)
        return predictions.mean(dim=0), predictions.std(dim=0, unbiased=False)


class DynamicsDataset(Dataset[tuple[Tensor, Tensor, Tensor]]):
    """Dataset of state-action-next_state transitions."""

    def __init__(self, states: Tensor, actions: Tensor, next_states: Tensor) -> None:
        if not (len(states) == len(actions) == len(next_states)):
            raise ValueError("States, actions, and next_states must have matching lengths.")
        self.states = states
        self.actions = actions
        self.next_states = next_states

    def __len__(self) -> int:
        return self.states.shape[0]

    def __getitem__(self, index: int) -> tuple[Tensor, Tensor, Tensor]:
        return self.states[index], self.actions[index], self.next_states[index]


@dataclass(slots=True)
class TrainingConfig:
    """Configuration for learned dynamics model fitting."""

    batch_size: int = 256
    epochs: int = 50
    learning_rate: float = 3e-4
    weight_decay: float = 1e-5
    grad_clip_norm: float | None = 1.0


class SimulationDataCollector:
    """Collect transition data from a simulator for dynamics learning."""

    def __init__(self, device: torch.device | str | None = None) -> None:
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))

    def collect(
        self,
        simulator: SimulatorProtocol,
        *,
        num_trajectories: int,
        horizon: int,
        policy: Callable[[Tensor], Tensor] | None = None,
    ) -> DynamicsDataset:
        states: list[Tensor] = []
        actions: list[Tensor] = []
        next_states: list[Tensor] = []

        for _ in range(num_trajectories):
            state = simulator.reset().to(self.device)
            for _ in range(horizon):
                action = policy(state) if policy is not None else simulator.sample_action(state)
                action = action.to(self.device)
                successor, _, done, _ = simulator.step(action)
                successor = successor.to(self.device)
                states.append(state.reshape(-1))
                actions.append(action.reshape(-1))
                next_states.append(successor.reshape(-1))
                state = successor
                if done:
                    break

        return DynamicsDataset(
            states=torch.stack(states),
            actions=torch.stack(actions),
            next_states=torch.stack(next_states),
        )


class DynamicsTrainer:
    """Trainer for single-model and ensemble dynamics learning."""

    def __init__(self, config: TrainingConfig | None = None, *, device: torch.device | str | None = None) -> None:
        self.config = config or TrainingConfig()
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))

    def fit(
        self,
        model: BaseDynamicsModel,
        dataset: DynamicsDataset,
        *,
        validation_dataset: DynamicsDataset | None = None,
    ) -> dict[str, list[float]]:
        model.to(self.device)
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=self.config.learning_rate,
            weight_decay=self.config.weight_decay,
        )
        train_loader = DataLoader(dataset, batch_size=self.config.batch_size, shuffle=True)
        validation_loader = (
            DataLoader(validation_dataset, batch_size=self.config.batch_size, shuffle=False)
            if validation_dataset is not None
            else None
        )

        history = {"train_loss": [], "val_loss": []}
        for _ in range(self.config.epochs):
            model.train()
            epoch_loss = 0.0
            for states, actions, next_states in train_loader:
                states = states.to(self.device)
                actions = actions.to(self.device)
                next_states = next_states.to(self.device)
                predictions = model(states, actions)
                loss = F.mse_loss(predictions, next_states)
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                if self.config.grad_clip_norm is not None:
                    nn.utils.clip_grad_norm_(model.parameters(), self.config.grad_clip_norm)
                optimizer.step()
                epoch_loss += loss.item() * states.shape[0]
            history["train_loss"].append(epoch_loss / len(dataset))

            if validation_loader is None:
                continue
            model.eval()
            validation_loss = 0.0
            with torch.no_grad():
                for states, actions, next_states in validation_loader:
                    states = states.to(self.device)
                    actions = actions.to(self.device)
                    next_states = next_states.to(self.device)
                    predictions = model(states, actions)
                    validation_loss += F.mse_loss(predictions, next_states, reduction="sum").item()
            history["val_loss"].append(validation_loss / len(validation_dataset))
        return history

    def fit_ensemble(self, ensemble: EnsembleDynamicsModel, dataset: DynamicsDataset) -> list[dict[str, list[float]]]:
        histories = []
        dataset_size = len(dataset)
        for member in ensemble.members:
            bootstrap_indices = torch.randint(0, dataset_size, (dataset_size,))
            bootstrap_dataset = DynamicsDataset(
                dataset.states[bootstrap_indices],
                dataset.actions[bootstrap_indices],
                dataset.next_states[bootstrap_indices],
            )
            histories.append(self.fit(member, bootstrap_dataset))
        return histories

    def collect_and_fit(
        self,
        model: BaseDynamicsModel,
        simulator: SimulatorProtocol,
        *,
        num_trajectories: int,
        horizon: int,
        policy: Callable[[Tensor], Tensor] | None = None,
    ) -> tuple[DynamicsDataset, dict[str, list[float]]]:
        collector = SimulationDataCollector(device=self.device)
        dataset = collector.collect(
            simulator,
            num_trajectories=num_trajectories,
            horizon=horizon,
            policy=policy,
        )
        history = self.fit(model, dataset)
        return dataset, history


def _fully_connected_edge_index(num_nodes: int, device: torch.device) -> Tensor:
    indices = torch.arange(num_nodes, device=device)
    src = indices.repeat_interleave(num_nodes)
    dst = indices.repeat(num_nodes)
    mask = src != dst
    return torch.stack([src[mask], dst[mask]], dim=0)


def _edge_features(current_state: Tensor, edge_index: Tensor) -> Tensor:
    src, dst = edge_index
    return current_state[:, dst] - current_state[:, src]


__all__ = [
    "BaseDynamicsModel",
    "DynamicsDataset",
    "DynamicsTrainer",
    "EnsembleDynamicsModel",
    "GNNDynamicsModel",
    "MLPDynamicsModel",
    "SimulationDataCollector",
    "SimulatorProtocol",
    "TrainingConfig",
]
