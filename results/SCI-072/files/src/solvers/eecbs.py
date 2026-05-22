import heapq
import time
import math
import random
import numpy as np
from typing import List, Dict, Optional, Tuple, Set
from collections import defaultdict, deque
import sys
from dataclasses import dataclass

sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')
from src.core.environment import GridEnvironment
from src.core.agent import Agent
from src.core.conflict import ConflictDetector, Conflict, ConflictType
from src.core.constraint import Constraint, ConstraintSet
from src.core.solution import Solution, Path
from src.solvers.base import MAPFSolver
from src.solvers.cbs import CBS


@dataclass
class EECBSNode:
    """High-level node for explicit-estimation CBS."""

    node_id: int
    cost: int
    constraints: ConstraintSet
    paths: Dict[int, List[Tuple[int, int]]]
    conflicts: List[Conflict]
    h_hat: float
    depth: int = 0


class EECBS(CBS):
    """Bounded-suboptimal CBS with OPEN/FOCAL management.

    This implementation keeps the exact CBS low-level search but uses a focal list
    to prefer nodes with fewer expected remaining conflicts while maintaining a
    user-provided suboptimality bound.
    """

    def __init__(
        self,
        env: GridEnvironment,
        agents: List[Agent],
        timeout: int = 300,
        suboptimality_bound: float = 1.5,
        max_expansions: int = 5000,
        seed: int = 0,
        **kwargs,
    ) -> None:
        super().__init__(env, agents, timeout, **kwargs)
        self.suboptimality_bound = max(1.0, float(suboptimality_bound))
        self.max_expansions = max_expansions
        self._rng = random.Random(seed)
        self._conflict_estimates: Dict[Tuple[object, ...], List[float]] = defaultdict(lambda: [0.0, 0.0])
        self.stats.update(
            {
                "suboptimality_bound": self.suboptimality_bound,
                "actual_suboptimality": math.inf,
                "nodes_expanded": 0,
                "runtime": 0.0,
            }
        )

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        try:
            root_paths = self._initial_paths()
            if root_paths is None:
                return None
            root = self._make_eecbs_node(ConstraintSet(), root_paths, 0)
            open_heap: List[Tuple[int, int, EECBSNode]] = [(root.cost, root.node_id, root)]
            active: Dict[int, EECBSNode] = {root.node_id: root}
            next_node_id = 1

            while active and self.stats["nodes_expanded"] < self.max_expansions:
                self._check_timeout()
                while open_heap and open_heap[0][1] not in active:
                    heapq.heappop(open_heap)
                if not open_heap:
                    break

                best_open_cost = open_heap[0][0]
                focal = [
                    node
                    for node in active.values()
                    if node.cost <= self.suboptimality_bound * best_open_cost
                ]
                current = min(
                    focal,
                    key=lambda node: (len(node.conflicts), node.h_hat, node.cost, node.depth, node.node_id),
                )
                active.pop(current.node_id, None)
                self.stats["nodes_expanded"] += 1
                self.stats["expanded_nodes"] = self.stats["nodes_expanded"]

                if not current.conflicts:
                    solution = Solution({aid: Path(aid, path) for aid, path in current.paths.items()})
                    lower_bound = max(1, best_open_cost)
                    self.stats["runtime"] = time.time() - self.start_time
                    self.stats["solution_cost"] = solution.cost
                    self.stats["actual_suboptimality"] = solution.cost / lower_bound
                    return solution

                conflict = self._choose_conflict(current.conflicts)
                parent_cost = current.cost
                for agent_id in (conflict.agent1, conflict.agent2):
                    child_constraints = current.constraints.copy()
                    child_constraints.add(self._constraint_from_conflict(conflict, agent_id))
                    child_paths = dict(current.paths)
                    agent = next(agent for agent in self.agents if agent.id == agent_id)
                    replanned = self._replan_agent_like_eecbs(current, agent, child_constraints)
                    if replanned is None:
                        continue
                    child_paths[agent_id] = replanned
                    child = self._make_eecbs_node(child_constraints, child_paths, next_node_id, current.depth + 1)
                    active[child.node_id] = child
                    heapq.heappush(open_heap, (child.cost, child.node_id, child))
                    next_node_id += 1
                    self.stats["generated_nodes"] += 1
                    self._update_estimate(conflict, max(0.0, child.cost - parent_cost))
        except TimeoutError:
            pass

        self.stats["runtime"] = time.time() - self.start_time
        return None

    def _make_eecbs_node(
        self,
        constraints: ConstraintSet,
        paths: Dict[int, List[Tuple[int, int]]],
        node_id: int,
        depth: int = 0,
    ) -> EECBSNode:
        conflicts = self._detect_conflicts(paths, constraints)
        h_hat = self._estimate_remaining_cost(conflicts)
        return EECBSNode(
            node_id=node_id,
            cost=self._compute_cost(paths),
            constraints=constraints,
            paths=paths,
            conflicts=conflicts,
            h_hat=h_hat,
            depth=depth,
        )

    def _replan_agent_like_eecbs(
        self,
        node: EECBSNode,
        agent: Agent,
        constraints: ConstraintSet,
    ) -> Optional[List[Tuple[int, int]]]:
        current_cost = len(node.paths.get(agent.id, [])) - 1
        latest = constraints.latest_timestep(agent.id)
        horizon = max(current_cost, latest) + getattr(self.env, "width", 1) * getattr(self.env, "height", 1) + 4
        return self.a_star(agent.start, agent.goal, constraints, agent=agent, max_timestep=horizon)

    def _choose_conflict(self, conflicts: List[Conflict]) -> Conflict:
        return min(
            conflicts,
            key=lambda conflict: (
                conflict.timestep,
                self._estimated_conflict_penalty(conflict),
                self._rng.random(),
            ),
        )

    def _estimate_remaining_cost(self, conflicts: List[Conflict]) -> float:
        if not conflicts:
            return 0.0
        penalties = [self._estimated_conflict_penalty(conflict) for conflict in conflicts]
        return float(len(conflicts) + np.mean(penalties))

    def _estimated_conflict_penalty(self, conflict: Conflict) -> float:
        total, count = self._conflict_estimates[self._conflict_signature(conflict)]
        if count <= 0:
            return 1.0
        return max(1.0, total / count)

    def _update_estimate(self, conflict: Conflict, delta: float) -> None:
        stats = self._conflict_estimates[self._conflict_signature(conflict)]
        stats[0] += delta
        stats[1] += 1.0

    def _conflict_signature(self, conflict: Conflict) -> Tuple[object, ...]:
        if conflict.conflict_type == ConflictType.EDGE:
            return (conflict.conflict_type.value, conflict.position, conflict.position2)
        return (conflict.conflict_type.value, conflict.position, conflict.timestep)
