import heapq
import time
from dataclasses import dataclass, field
from functools import lru_cache
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
class ICTNode:
    priority: Tuple[int, int]
    total_cost: int = field(compare=False)
    cost_vector: Tuple[int, ...] = field(compare=False)


class MDD:
    """Directed acyclic graph of all paths with exact cost c.

    Construction cost is O(c * |V| * degree).
    """

    def __init__(self, solver: MAPFSolver, agent: Agent, cost: int):
        self.solver = solver
        self.agent = agent
        self.cost = cost
        self.levels: Dict[int, Set[Tuple[int, int]]] = defaultdict(set)
        self.edges: Dict[Tuple[int, Tuple[int, int]], Set[Tuple[int, int]]] = defaultdict(set)
        self._build()

    def _build(self) -> None:
        start, goal = self.agent.start, self.agent.goal
        self.levels[0].add(start)
        for t in range(self.cost):
            for loc in list(self.levels[t]):
                for nxt in self.solver._neighbors(loc) + [loc]:
                    if t + 1 + self.solver._heuristic(nxt, goal) > self.cost:
                        continue
                    self.levels[t + 1].add(nxt)
                    self.edges[(t, loc)].add(nxt)
        valid = {self.cost: {goal}}
        for t in range(self.cost - 1, -1, -1):
            valid[t] = set()
            for loc in self.levels[t]:
                for nxt in self.edges.get((t, loc), set()):
                    if nxt in valid[t + 1]:
                        valid[t].add(loc)
                        break
        self.levels = defaultdict(set, valid)
        filtered_edges: Dict[Tuple[int, Tuple[int, int]], Set[Tuple[int, int]]] = defaultdict(set)
        for (t, loc), successors in self.edges.items():
            if loc not in self.levels[t]:
                continue
            for nxt in successors:
                if nxt in self.levels[t + 1]:
                    filtered_edges[(t, loc)].add(nxt)
        self.edges = filtered_edges


class ICTS(MAPFSolver):
    """Optimal Increasing Cost Tree Search.

    Worst-case complexity is exponential in the number of agents due to joint MDD search,
    but pairwise pruning often removes most cost vectors in practice.
    """

    def __init__(self, env, agents, timeout=300, **kwargs):
        super().__init__(env, agents, timeout)
        self.stats = {
            "expanded_nodes": 0,
            "generated_nodes": 0,
            "runtime": 0.0,
            "solution_cost": None,
            "mdd_computations": 0,
        }
        self._mdd_cache: Dict[Tuple[int, int], MDD] = {}

    def _get_mdd(self, agent: Agent, cost: int) -> Optional[MDD]:
        shortest = self.env.shortest_path_length(agent.start, agent.goal)
        if shortest is None or cost < shortest:
            return None
        key = (agent.id, cost)
        if key not in self._mdd_cache:
            self.stats["mdd_computations"] += 1
            self._mdd_cache[key] = MDD(self, agent, cost)
        return self._mdd_cache[key]

    def _pairwise_feasible(self, mdd1: MDD, mdd2: MDD) -> bool:
        horizon = max(mdd1.cost, mdd2.cost)

        @lru_cache(maxsize=None)
        def dfs(t: int, loc1: Tuple[int, int], loc2: Tuple[int, int]) -> bool:
            if loc1 == loc2:
                return False
            if t == horizon:
                return True
            succ1 = mdd1.edges.get((min(t, mdd1.cost - 1), loc1), {loc1} if t >= mdd1.cost else set())
            succ2 = mdd2.edges.get((min(t, mdd2.cost - 1), loc2), {loc2} if t >= mdd2.cost else set())
            if not succ1:
                succ1 = {loc1}
            if not succ2:
                succ2 = {loc2}
            for nxt1 in succ1:
                for nxt2 in succ2:
                    if nxt1 == nxt2:
                        continue
                    if loc1 == nxt2 and loc2 == nxt1:
                        continue
                    if dfs(t + 1, nxt1, nxt2):
                        return True
            return False

        return dfs(0, mdd1.agent.start, mdd2.agent.start)

    def _joint_solution(self, mdds: List[MDD]) -> Optional[Dict[int, List[Tuple[int, int]]]]:
        horizon = max(mdd.cost for mdd in mdds)

        def successors(mdd: MDD, t: int, loc: Tuple[int, int]) -> Set[Tuple[int, int]]:
            if t >= mdd.cost:
                return {loc}
            return mdd.edges.get((t, loc), {loc})

        @lru_cache(maxsize=None)
        def dfs(t: int, state: Tuple[Tuple[int, int], ...]) -> Optional[Tuple[Tuple[int, int], ...]]:
            if len(set(state)) < len(state):
                return None
            if t == horizon:
                return tuple()
            successor_lists = [sorted(successors(mdd, t, state[i])) for i, mdd in enumerate(mdds)]

            def backtrack(idx: int, chosen: List[Tuple[int, int]]) -> Optional[Tuple[Tuple[int, int], ...]]:
                if idx == len(mdds):
                    for i in range(len(chosen)):
                        for j in range(i + 1, len(chosen)):
                            if chosen[i] == chosen[j]:
                                return None
                            if state[i] == chosen[j] and state[j] == chosen[i]:
                                return None
                    suffix = dfs(t + 1, tuple(chosen))
                    if suffix is None:
                        return None
                    return tuple(chosen) + suffix
                for loc in successor_lists[idx]:
                    chosen.append(loc)
                    result = backtrack(idx + 1, chosen)
                    chosen.pop()
                    if result is not None:
                        return result
                return None

            return backtrack(0, [])

        start_state = tuple(mdd.agent.start for mdd in mdds)
        flattened = dfs(0, start_state)
        if flattened is None:
            return None

        per_agent = {mdd.agent.id: [mdd.agent.start] for mdd in mdds}
        width = len(mdds)
        for offset in range(0, len(flattened), width):
            joint = flattened[offset: offset + width]
            for i, mdd in enumerate(mdds):
                per_agent[mdd.agent.id].append(joint[i])
        return per_agent

    def _trim_path(self, path: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        trimmed = list(path)
        while len(trimmed) > 1 and trimmed[-1] == trimmed[-2]:
            trimmed.pop()
        return trimmed

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        try:
            root_costs = []
            for agent in self.agents:
                dist = self.env.shortest_path_length(agent.start, agent.goal)
                if dist is None:
                    return None
                root_costs.append(dist)
            root_vector = tuple(root_costs)
            open_list = [ICTNode((sum(root_vector), 0), sum(root_vector), root_vector)]
            closed = {root_vector}
            self.stats["generated_nodes"] = 1

            while open_list:
                self._check_timeout()
                node = heapq.heappop(open_list)
                self.stats["expanded_nodes"] += 1
                mdds = []
                feasible = True
                for agent, cost in zip(self.agents, node.cost_vector):
                    mdd = self._get_mdd(agent, cost)
                    if mdd is None or not mdd.levels.get(0):
                        feasible = False
                        break
                    mdds.append(mdd)
                if not feasible:
                    continue
                for i in range(len(mdds)):
                    for j in range(i + 1, len(mdds)):
                        if not self._pairwise_feasible(mdds[i], mdds[j]):
                            feasible = False
                            break
                    if not feasible:
                        break
                if feasible:
                    joint = self._joint_solution(mdds)
                    if joint is not None:
                        solution = Solution({aid: Path(aid, self._trim_path(path)) for aid, path in joint.items()})
                        self.stats["runtime"] = time.time() - self.start_time
                        self.stats["solution_cost"] = solution.cost
                        return solution
                for idx in range(len(node.cost_vector)):
                    child_vector = list(node.cost_vector)
                    child_vector[idx] += 1
                    child_vector = tuple(child_vector)
                    if child_vector in closed:
                        continue
                    closed.add(child_vector)
                    child = ICTNode((sum(child_vector), self.stats["generated_nodes"]), sum(child_vector), child_vector)
                    heapq.heappush(open_list, child)
                    self.stats["generated_nodes"] += 1
        except TimeoutError:
            pass
        self.stats["runtime"] = time.time() - self.start_time
        return None
