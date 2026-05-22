import heapq
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict
import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')
from src.core.environment import GridEnvironment
from src.core.agent import Agent
from src.core.conflict import ConflictDetector, Conflict, ConflictType
from src.core.constraint import Constraint, ConstraintSet
from src.core.solution import Solution, Path
from src.solvers.base import MAPFSolver


@dataclass(order=True)
class CTNode:
    priority: Tuple[int, int, int]
    cost: int = field(compare=False)
    constraints: ConstraintSet = field(compare=False)
    paths: Dict[int, List[Tuple[int, int]]] = field(compare=False)
    conflicts: List[Conflict] = field(compare=False, default_factory=list)


class CBS(MAPFSolver):
    """Optimal Conflict-Based Search.

    High level complexity is exponential in the number of conflicts in the worst case;
    low-level replans use optimal space-time A*.
    """

    def __init__(self, env, agents, timeout=300, **kwargs):
        super().__init__(env, agents, timeout)
        self.conflict_detector = kwargs.get("conflict_detector", ConflictDetector())
        self.stats = {
            "expanded_nodes": 0,
            "generated_nodes": 0,
            "runtime": 0.0,
            "solution_cost": None,
        }

    def _detect_conflicts(self, paths: Dict[int, List[Tuple[int, int]]], constraints: ConstraintSet) -> List[Conflict]:
        conflicts = self.conflict_detector.detect_conflicts(paths)
        scored = []
        for conflict in conflicts:
            cardinality = self._classify_conflict(conflict, paths, constraints)
            scored.append(Conflict(
                conflict.agent1,
                conflict.agent2,
                conflict.timestep,
                conflict.conflict_type,
                conflict.position,
                conflict.position2,
                cardinality,
            ))
        priority = {"cardinal": 0, "semi-cardinal": 1, "non-cardinal": 2, "unknown": 3}
        scored.sort(key=lambda c: (priority.get(c.cardinality, 3), c.timestep))
        return scored

    def _build_exact_mdd(self, agent: Agent, cost: int, constraints: ConstraintSet) -> Dict[int, Set[Tuple[int, int]]]:
        levels: Dict[int, Set[Tuple[int, int]]] = defaultdict(set)
        levels[0].add(agent.start)
        for t in range(cost):
            for loc in list(levels[t]):
                for nxt in self._neighbors(loc) + [loc]:
                    nt = t + 1
                    if constraints.is_forbidden(agent.id, nxt, nt, loc):
                        continue
                    if nt + self._heuristic(nxt, agent.goal) > cost:
                        continue
                    levels[nt].add(nxt)
        valid: Dict[int, Set[Tuple[int, int]]] = defaultdict(set)
        valid[cost].add(agent.goal)
        for t in range(cost - 1, -1, -1):
            for loc in levels[t]:
                for nxt in self._neighbors(loc) + [loc]:
                    if nxt not in valid[t + 1]:
                        continue
                    if constraints.is_forbidden(agent.id, nxt, t + 1, loc):
                        continue
                    if t + self._heuristic(loc, agent.goal) <= cost:
                        valid[t].add(loc)
                        break
        return valid

    def _classify_conflict(self, conflict: Conflict, paths: Dict[int, List[Tuple[int, int]]], constraints: ConstraintSet) -> str:
        agent1 = next(agent for agent in self.agents if agent.id == conflict.agent1)
        agent2 = next(agent for agent in self.agents if agent.id == conflict.agent2)
        cost1 = len(paths[agent1.id]) - 1
        cost2 = len(paths[agent2.id]) - 1
        mdd1 = self._build_exact_mdd(agent1, cost1, constraints)
        mdd2 = self._build_exact_mdd(agent2, cost2, constraints)
        width1 = len(mdd1.get(conflict.timestep, set()))
        width2 = len(mdd2.get(conflict.timestep, set()))
        if width1 == 1 and width2 == 1:
            return "cardinal"
        if width1 == 1 or width2 == 1:
            return "semi-cardinal"
        return "non-cardinal"

    def _compute_cost(self, paths: Dict[int, List[Tuple[int, int]]]) -> int:
        return sum(max(0, len(path) - 1) for path in paths.values())

    def _make_node(self, constraints: ConstraintSet, paths: Dict[int, List[Tuple[int, int]]]) -> CTNode:
        conflicts = self._detect_conflicts(paths, constraints)
        cost = self._compute_cost(paths)
        priority = (cost, len(conflicts), self.stats["generated_nodes"])
        return CTNode(priority, cost, constraints, paths, conflicts)

    def _initial_paths(self) -> Optional[Dict[int, List[Tuple[int, int]]]]:
        paths = {}
        empty_constraints = ConstraintSet()
        for agent in self.agents:
            path = self.a_star(agent.start, agent.goal, empty_constraints, agent=agent)
            if path is None:
                return None
            paths[agent.id] = path
        return paths

    def _constraint_from_conflict(self, conflict: Conflict, agent_id: int) -> Constraint:
        if conflict.conflict_type == ConflictType.VERTEX:
            return Constraint(agent=agent_id, timestep=conflict.timestep, position=conflict.position)
        if agent_id == conflict.agent1:
            edge = (conflict.position, conflict.position2)
        else:
            edge = (conflict.position2, conflict.position)
        return Constraint(agent=agent_id, timestep=conflict.timestep, edge=edge)

    def _replan_agent(self, node: CTNode, agent: Agent, new_constraints: ConstraintSet) -> Optional[List[Tuple[int, int]]]:
        current_cost = len(node.paths[agent.id]) - 1 if agent.id in node.paths else None
        max_timestep = max(current_cost or 0, new_constraints.latest_timestep(agent.id)) + getattr(self.env, "width", 1) * getattr(self.env, "height", 1) + 2
        return self.a_star(agent.start, agent.goal, new_constraints, agent=agent, max_timestep=max_timestep)

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        try:
            root_paths = self._initial_paths()
            if root_paths is None:
                return None
            root = self._make_node(ConstraintSet(), root_paths)
            open_list = [root]
            self.stats["generated_nodes"] = 1

            while open_list:
                self._check_timeout()
                node = heapq.heappop(open_list)
                self.stats["expanded_nodes"] += 1
                if not node.conflicts:
                    solution = Solution({aid: Path(aid, path) for aid, path in node.paths.items()})
                    self.stats["runtime"] = time.time() - self.start_time
                    self.stats["solution_cost"] = solution.cost
                    return solution

                conflict = node.conflicts[0]
                children: List[CTNode] = []
                for agent_id in (conflict.agent1, conflict.agent2):
                    child_constraints = node.constraints.copy()
                    child_constraints.add(self._constraint_from_conflict(conflict, agent_id))
                    child_paths = dict(node.paths)
                    agent = next(agent for agent in self.agents if agent.id == agent_id)
                    replanned = self._replan_agent(node, agent, child_constraints)
                    if replanned is None:
                        continue
                    child_paths[agent_id] = replanned
                    child = self._make_node(child_constraints, child_paths)
                    self.stats["generated_nodes"] += 1
                    children.append(child)

                bypassed = False
                for child in children:
                    if child.cost == node.cost and len(child.conflicts) < len(node.conflicts):
                        heapq.heappush(open_list, child)
                        bypassed = True
                        break
                if bypassed:
                    continue
                for child in children:
                    heapq.heappush(open_list, child)
        except TimeoutError:
            pass
        self.stats["runtime"] = time.time() - self.start_time
        return None
