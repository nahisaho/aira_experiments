from abc import ABC, abstractmethod
import heapq
import time
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')
from src.core.environment import GridEnvironment
from src.core.agent import Agent
from src.core.conflict import ConflictDetector, Conflict, ConflictType
from src.core.constraint import Constraint, ConstraintSet
from src.core.solution import Solution, Path


class MAPFSolver(ABC):
    """Base class shared by all MAPF solvers.

    The low-level search is space-time A* over states (x, y, timestep).
    Complexity: O(T * |V| log(T * |V|)) in the explored time horizon T.
    """

    def __init__(self, env, agents, timeout=300):
        self.env = env
        self.agents = agents
        self.timeout = timeout
        self.stats = {}
        self.start_time = None

    @abstractmethod
    def solve(self) -> Optional[Solution]:
        pass

    def _check_timeout(self) -> None:
        if self.start_time is not None and (time.time() - self.start_time) > self.timeout:
            raise TimeoutError(f"{self.__class__.__name__} timed out after {self.timeout}s")

    def _agent_id(self, agent) -> int:
        return getattr(agent, "id", getattr(agent, "name", agent))

    def _heuristic(self, start: Tuple[int, int], goal: Tuple[int, int]) -> int:
        if hasattr(self.env, "heuristic"):
            return int(self.env.heuristic(start, goal))
        return abs(start[0] - goal[0]) + abs(start[1] - goal[1])

    def _neighbors(self, loc: Tuple[int, int]) -> List[Tuple[int, int]]:
        if hasattr(self.env, "neighbors"):
            try:
                return list(self.env.neighbors(loc, include_wait=False))
            except TypeError:
                return list(self.env.neighbors(loc))
        x, y = loc
        candidates = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]
        valid = []
        for nxt in candidates:
            if hasattr(self.env, "valid") and self.env.valid(nxt):
                valid.append(nxt)
            elif hasattr(self.env, "in_bounds") and hasattr(self.env, "passable") and self.env.in_bounds(nxt) and self.env.passable(nxt):
                valid.append(nxt)
        return valid

    def a_star(
        self,
        start,
        goal,
        constraints=None,
        heuristic=None,
        max_timestep=None,
        agent=None,
    ) -> Optional[List]:
        """Space-time A* with wait actions and constraint checking.

        Tie-breaking prefers lower estimated conflict pressure, then smaller heuristic,
        then deeper g-value. This keeps paths optimal while biasing the search toward
        cleaner plans under dense constraints.
        """
        constraint_set = constraints or ConstraintSet()
        heuristic_fn = heuristic or self._heuristic
        agent_id = self._agent_id(agent) if agent is not None else getattr(start, "agent", 0)
        latest_constraint = constraint_set.latest_timestep(agent_id)
        if max_timestep is None:
            width = getattr(self.env, "width", 1)
            height = getattr(self.env, "height", 1)
            max_timestep = max(width * height * 2, heuristic_fn(start, goal) * 2 + 4, latest_constraint + 2)

        if constraint_set.is_forbidden(agent_id, start, 0):
            return None

        open_list = []
        counter = 0
        start_h = heuristic_fn(start, goal)
        heapq.heappush(open_list, (start_h, 0, start_h, 0, counter, start, 0))
        parents: Dict[Tuple[Tuple[int, int], int], Optional[Tuple[Tuple[int, int], int]]] = {(start, 0): None}
        g_scores = defaultdict(lambda: float("inf"))
        g_scores[(start, 0)] = 0

        while open_list:
            self._check_timeout()
            _, conflict_bias, h, neg_g, _, loc, timestep = heapq.heappop(open_list)
            g = -neg_g
            state = (loc, timestep)
            if g != g_scores[state]:
                continue

            if loc == goal and timestep >= latest_constraint and constraint_set.is_goal_safe(agent_id, goal, timestep, max_timestep):
                path = []
                current = state
                while current is not None:
                    path.append(current[0])
                    current = parents[current]
                return list(reversed(path))

            if timestep >= max_timestep:
                continue

            for nxt in self._neighbors(loc) + [loc]:
                next_timestep = timestep + 1
                if constraint_set.is_forbidden(agent_id, nxt, next_timestep, loc):
                    continue
                tentative_g = g + 1
                next_state = (nxt, next_timestep)
                if tentative_g >= g_scores[next_state]:
                    continue
                g_scores[next_state] = tentative_g
                parents[next_state] = state
                next_h = heuristic_fn(nxt, goal)
                future_conflicts = constraint_set.count_future_conflicts(agent_id, nxt, next_timestep)
                counter += 1
                heapq.heappush(
                    open_list,
                    (tentative_g + next_h, future_conflicts, next_h, -tentative_g, counter, nxt, next_timestep),
                )
        return None
