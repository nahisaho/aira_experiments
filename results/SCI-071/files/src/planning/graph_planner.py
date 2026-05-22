from __future__ import annotations

import heapq
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import torch
from torch import Tensor

from .base_planner import ActionSequence, BasePlanner, StateType, batched_flatten_state, interpolate_state


class LatentModelProtocol(Protocol):
    """Interface for latent encoders used by the graph planner."""

    def encode(self, state: Tensor) -> Tensor:
        """Encode a state tensor into latent space."""


@runtime_checkable
class DecodableLatentModelProtocol(LatentModelProtocol, Protocol):
    """Extension of latent encoders that can decode latent states."""

    def decode(self, latent: Tensor) -> Tensor:
        """Decode a latent vector back into state space."""


class FeasibilitySimulatorProtocol(Protocol):
    """Simulator hook for graph-edge feasibility validation."""

    def is_transition_feasible(self, start_state: StateType, goal_state: StateType) -> bool:
        """Return whether the transition between two subgoals is feasible."""


@dataclass(slots=True)
class GraphPlannerConfig:
    """Configuration for graph-based subgoal planning."""

    horizon: int = 10
    num_keypoints: int = 8
    num_subgoals: int = 6
    neighbor_count: int = 3


@dataclass(slots=True)
class GraphNode:
    """Node in the subgoal search graph."""

    index: int
    state: StateType
    latent: Tensor


class GraphPlanner(BasePlanner):
    """Graph planner using latent subgoals and A* search."""

    def __init__(
        self,
        latent_model: LatentModelProtocol | None,
        *,
        local_planner: BasePlanner | None = None,
        simulator: FeasibilitySimulatorProtocol | None = None,
        config: GraphPlannerConfig | None = None,
        action_dim: int = 0,
        device: torch.device | str | None = None,
    ) -> None:
        self.config = config or GraphPlannerConfig()
        self.latent_model = latent_model
        self.local_planner = local_planner
        self.simulator = simulator
        if local_planner is not None:
            action_dim = local_planner.action_dim
        super().__init__(action_dim, device=device, default_horizon=self.config.horizon)

    def plan(self, current_state: StateType, goal_state: StateType, horizon: int) -> ActionSequence:
        nodes = self._build_nodes(current_state, goal_state)
        graph = self._build_graph(nodes)
        path = self._a_star(nodes, graph)
        actions = []
        segment_plans = []
        for start_idx, goal_idx in zip(path[:-1], path[1:]):
            start_state = nodes[start_idx].state
            target_state = nodes[goal_idx].state
            if self.local_planner is not None:
                segment_plan = self.local_planner.plan(start_state, target_state, max(1, horizon // max(1, len(path) - 1)))
                actions.append(segment_plan.actions)
                segment_plans.append(segment_plan)
        action_tensor = (
            torch.cat(actions, dim=0)
            if actions
            else torch.zeros(max(0, len(path) - 1), self.action_dim, device=self.device)
        )
        return ActionSequence(
            actions=action_tensor,
            predicted_states=[nodes[index].state for index in path[1:]],
            cost=float(len(path) - 1),
            metadata={
                "planner": "graph",
                "path": path,
                "keypoints": self.extract_keypoints(current_state, self.config.num_keypoints),
                "segment_plans": segment_plans,
            },
        )

    def replan(
        self,
        current_state: StateType,
        goal_state: StateType,
        previous_plan: ActionSequence | None,
    ) -> ActionSequence:
        return self.plan(current_state, goal_state, self.default_horizon)

    def extract_keypoints(self, state: StateType, num_keypoints: int) -> Tensor:
        """Extract representative keypoints from a deformable object state."""
        points = self._state_to_points(state)
        if points.shape[0] <= num_keypoints:
            return points
        selected_indices = [0]
        distances = torch.full((points.shape[0],), float("inf"), device=points.device)
        for _ in range(1, num_keypoints):
            last_point = points[selected_indices[-1]].unsqueeze(0)
            distances = torch.minimum(distances, torch.cdist(points, last_point).squeeze(-1))
            selected_indices.append(int(distances.argmax().item()))
        return points[selected_indices]

    def generate_subgoals(self, current_state: StateType, goal_state: StateType) -> list[StateType]:
        """Generate latent-space subgoals between the current and goal states."""
        if self.latent_model is not None and isinstance(self.latent_model, DecodableLatentModelProtocol):
            current_latent = self._encode_state(current_state)
            goal_latent = self._encode_state(goal_state)
            return [
                self.latent_model.decode(torch.lerp(current_latent, goal_latent, alpha))
                for alpha in torch.linspace(0.0, 1.0, self.config.num_subgoals + 2, device=current_latent.device)[1:-1]
            ]
        return [
            interpolate_state(current_state, goal_state, float(alpha))
            for alpha in torch.linspace(0.0, 1.0, self.config.num_subgoals + 2)[1:-1]
        ]

    def _build_nodes(self, current_state: StateType, goal_state: StateType) -> list[GraphNode]:
        states = [current_state, *self.generate_subgoals(current_state, goal_state), goal_state]
        return [GraphNode(index=i, state=state, latent=self._encode_state(state)) for i, state in enumerate(states)]

    def _build_graph(self, nodes: list[GraphNode]) -> dict[int, list[tuple[int, float]]]:
        graph = {node.index: [] for node in nodes}
        latent_matrix = torch.stack([node.latent for node in nodes], dim=0)
        pairwise = torch.cdist(latent_matrix, latent_matrix)
        for node in nodes:
            neighbor_indices = torch.topk(
                pairwise[node.index],
                k=min(self.config.neighbor_count + 1, len(nodes)),
                largest=False,
            ).indices.tolist()
            for neighbor_index in neighbor_indices:
                if neighbor_index == node.index:
                    continue
                if self.simulator is not None and not self.simulator.is_transition_feasible(node.state, nodes[neighbor_index].state):
                    continue
                graph[node.index].append((neighbor_index, float(pairwise[node.index, neighbor_index].item())))
        return graph

    def _a_star(self, nodes: list[GraphNode], graph: dict[int, list[tuple[int, float]]]) -> list[int]:
        start = 0
        goal = len(nodes) - 1
        frontier: list[tuple[float, int]] = [(0.0, start)]
        came_from = {start: None}
        g_costs = {start: 0.0}

        while frontier:
            _, current = heapq.heappop(frontier)
            if current == goal:
                break
            for neighbor, edge_cost in graph[current]:
                candidate_cost = g_costs[current] + edge_cost
                if neighbor not in g_costs or candidate_cost < g_costs[neighbor]:
                    g_costs[neighbor] = candidate_cost
                    heuristic = torch.norm(nodes[neighbor].latent - nodes[goal].latent).item()
                    heapq.heappush(frontier, (candidate_cost + heuristic, neighbor))
                    came_from[neighbor] = current

        if goal not in came_from:
            raise RuntimeError("No feasible subgoal path found by A* search.")

        path = [goal]
        cursor = goal
        while came_from[cursor] is not None:
            cursor = came_from[cursor]
            path.append(cursor)
        return list(reversed(path))

    def _encode_state(self, state: StateType) -> Tensor:
        flattened = batched_flatten_state(state).squeeze(0).to(self.device)
        if self.latent_model is None:
            return flattened
        return self.latent_model.encode(flattened.unsqueeze(0)).squeeze(0)

    def _state_to_points(self, state: StateType) -> Tensor:
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
        return tensor_state.to(self.device)


__all__ = [
    "DecodableLatentModelProtocol",
    "FeasibilitySimulatorProtocol",
    "GraphPlanner",
    "GraphPlannerConfig",
    "GraphNode",
    "LatentModelProtocol",
]
