import heapq
import time
import math
import random
import numpy as np
from typing import List, Dict, Optional, Tuple, Set
from collections import defaultdict, deque
import sys
from dataclasses import dataclass, field

sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')
from src.core.environment import GridEnvironment
from src.core.agent import Agent
from src.core.conflict import ConflictDetector, Conflict, ConflictType
from src.core.constraint import Constraint, ConstraintSet
from src.core.solution import Solution, Path
from src.solvers.base import MAPFSolver


@dataclass
class Configuration:
    """Joint agent configuration at a particular timestep."""

    timestep: int
    positions: Dict[int, Tuple[int, int]]


@dataclass
class ConfigNode:
    """Configuration search node with lazily accumulated constraints."""

    configuration: Configuration
    parent: Optional["ConfigNode"] = None
    constraints: Dict[int, Set[Tuple[int, Tuple[int, int]]]] = field(default_factory=lambda: defaultdict(set))


class Generator:
    """PIBT-style successor generator with temporary priority inheritance."""

    def __init__(self, env: GridEnvironment, rng: random.Random) -> None:
        self.env = env
        self.rng = rng

    def generate(
        self,
        positions: Dict[int, Tuple[int, int]],
        goals: Dict[int, Tuple[int, int]],
        priorities: Dict[int, float],
        guides: Dict[int, List[Tuple[int, int]]],
        timestep: int,
        constraints: Dict[int, Set[Tuple[int, Tuple[int, int]]]],
    ) -> Dict[int, Tuple[int, int]]:
        occupancy = {loc: aid for aid, loc in positions.items()}
        reserved: Dict[Tuple[int, int], int] = {}
        chosen: Dict[int, Tuple[int, int]] = {}

        def assign(agent_id: int, chain: Tuple[int, ...] = ()) -> Tuple[int, int]:
            if agent_id in chosen:
                return chosen[agent_id]
            for target in self._candidate_moves(agent_id, positions, goals, guides, timestep):
                if (timestep + 1, target) in constraints.get(agent_id, set()):
                    continue
                if target in reserved and reserved[target] != agent_id:
                    continue
                blocker = occupancy.get(target)
                if blocker is not None and blocker != agent_id:
                    if blocker in chain:
                        continue
                    inherited = priorities[blocker]
                    priorities[blocker] = max(priorities[blocker], priorities[agent_id] + 1e-3)
                    blocker_target = assign(blocker, chain + (agent_id,))
                    priorities[blocker] = inherited
                    if blocker_target == target:
                        continue
                chosen[agent_id] = target
                reserved[target] = agent_id
                return target
            chosen[agent_id] = positions[agent_id]
            reserved[positions[agent_id]] = agent_id
            return chosen[agent_id]

        for agent_id in sorted(priorities, key=lambda aid: (-priorities[aid], aid)):
            assign(agent_id)
        return chosen

    def _candidate_moves(
        self,
        agent_id: int,
        positions: Dict[int, Tuple[int, int]],
        goals: Dict[int, Tuple[int, int]],
        guides: Dict[int, List[Tuple[int, int]]],
        timestep: int,
    ) -> List[Tuple[int, int]]:
        current = positions[agent_id]
        goal = goals[agent_id]
        candidates: List[Tuple[int, int]] = []
        guide = guides.get(agent_id, [])
        if timestep + 1 < len(guide):
            candidates.append(guide[timestep + 1])
        candidates.extend(self.env.neighbors(current, include_wait=False))
        candidates.append(current)
        unique = []
        seen: Set[Tuple[int, int]] = set()
        for candidate in sorted(candidates, key=lambda loc: (self.env.heuristic(loc, goal), loc != current, self.rng.random())):
            if candidate in seen:
                continue
            seen.add(candidate)
            unique.append(candidate)
        return unique


class LaCAM(MAPFSolver):
    """Lazy Constraints Addition for MAPF.

    The solver starts from independent plans and incrementally repairs only the
    local conflicts that are actually observed while simulating joint execution.
    """

    def __init__(
        self,
        env: GridEnvironment,
        agents: List[Agent],
        timeout: int = 300,
        max_timesteps: Optional[int] = None,
        seed: int = 0,
        **kwargs,
    ) -> None:
        super().__init__(env, agents, timeout)
        horizon = max((self._heuristic(agent.start, agent.goal) for agent in agents), default=0)
        self.max_timesteps = max_timesteps or max(16, horizon * 4 + len(agents) * 2)
        self._rng = random.Random(seed)
        np.random.seed(seed)
        self.generator = Generator(env, self._rng)
        self.stats = {
            "configurations_explored": 0,
            "runtime": 0.0,
            "solution_cost": None,
        }

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        guides = self._independent_paths({agent.id: agent.start for agent in self.agents})
        if any(agent.id not in guides for agent in self.agents):
            return None

        current_positions = {agent.id: agent.start for agent in self.agents}
        goals = {agent.id: agent.goal for agent in self.agents}
        traces: Dict[int, List[Tuple[int, int]]] = {agent.id: [agent.start] for agent in self.agents}
        priorities = {agent.id: self._rng.random() for agent in self.agents}
        lazy_constraints: Dict[int, Set[Tuple[int, Tuple[int, int]]]] = defaultdict(set)
        current_node = ConfigNode(Configuration(0, dict(current_positions)))

        try:
            for timestep in range(self.max_timesteps):
                self._check_timeout()
                if all(current_positions[agent.id] == agent.goal for agent in self.agents):
                    break
                guides = self._refresh_guides(guides, current_positions, goals)
                next_positions = self.generator.generate(
                    current_positions,
                    goals,
                    priorities,
                    guides,
                    timestep,
                    lazy_constraints,
                )
                self.stats["configurations_explored"] += 1
                conflicts = self._step_conflicts(current_positions, next_positions)
                repairs = 0
                while conflicts and repairs < len(self.agents) * 2:
                    for conflict in conflicts:
                        loser = min((conflict.agent1, conflict.agent2), key=lambda aid: (priorities[aid], aid))
                        blocked = self._blocked_location(conflict, loser)
                        lazy_constraints[loser].add((timestep + 1, blocked))
                    next_positions = self.generator.generate(
                        current_positions,
                        goals,
                        priorities,
                        guides,
                        timestep,
                        lazy_constraints,
                    )
                    self.stats["configurations_explored"] += 1
                    conflicts = self._step_conflicts(current_positions, next_positions)
                    repairs += 1

                current_positions = next_positions
                current_node = ConfigNode(Configuration(timestep + 1, dict(current_positions)), current_node, lazy_constraints)
                for agent in self.agents:
                    traces[agent.id].append(current_positions[agent.id])
        except TimeoutError:
            pass

        solution = Solution({aid: Path(aid, path) for aid, path in traces.items()})
        self.stats["runtime"] = time.time() - self.start_time
        self.stats["solution_cost"] = solution.cost
        return solution

    def _independent_paths(self, starts: Dict[int, Tuple[int, int]]) -> Dict[int, List[Tuple[int, int]]]:
        paths: Dict[int, List[Tuple[int, int]]] = {}
        for agent in self.agents:
            path = self.a_star(starts[agent.id], agent.goal, ConstraintSet(), agent=agent)
            if path is not None:
                paths[agent.id] = path
        return paths

    def _refresh_guides(
        self,
        guides: Dict[int, List[Tuple[int, int]]],
        current_positions: Dict[int, Tuple[int, int]],
        goals: Dict[int, Tuple[int, int]],
    ) -> Dict[int, List[Tuple[int, int]]]:
        refreshed = dict(guides)
        for agent in self.agents:
            guide = refreshed.get(agent.id)
            if guide is None or not guide or guide[0] != current_positions[agent.id]:
                refreshed[agent.id] = self.a_star(current_positions[agent.id], goals[agent.id], ConstraintSet(), agent=agent) or [current_positions[agent.id]]
        return refreshed

    def _step_conflicts(
        self,
        current_positions: Dict[int, Tuple[int, int]],
        next_positions: Dict[int, Tuple[int, int]],
    ) -> List[Conflict]:
        step_paths = {aid: [current_positions[aid], next_positions[aid]] for aid in current_positions}
        return ConflictDetector.detect_conflicts(step_paths)

    def _blocked_location(self, conflict: Conflict, loser: int) -> Tuple[int, int]:
        if conflict.conflict_type == ConflictType.VERTEX:
            return conflict.position
        return conflict.position2 if loser == conflict.agent1 else conflict.position
