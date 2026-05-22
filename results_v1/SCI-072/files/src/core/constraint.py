from dataclasses import dataclass
from typing import DefaultDict, Dict, Iterable, Optional, Set, Tuple
from collections import defaultdict


@dataclass(frozen=True)
class Constraint:
    """Single-agent vertex or edge prohibition."""

    agent: int
    timestep: int
    position: Optional[Tuple[int, int]] = None
    edge: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None

    @property
    def is_vertex(self) -> bool:
        return self.position is not None

    @property
    def is_edge(self) -> bool:
        return self.edge is not None


class ConstraintSet:
    """Lookup-friendly container for CBS/ICTS style constraints."""

    def __init__(self, constraints: Optional[Iterable[Constraint]] = None):
        self.vertex_constraints: DefaultDict[int, Dict[int, Set[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(set))
        self.edge_constraints: DefaultDict[int, Dict[int, Set[Tuple[Tuple[int, int], Tuple[int, int]]]]] = defaultdict(lambda: defaultdict(set))
        self.constraints = []
        if constraints:
            for constraint in constraints:
                self.add(constraint)

    def copy(self) -> "ConstraintSet":
        return ConstraintSet(self.constraints)

    def add(self, constraint: Constraint) -> None:
        self.constraints.append(constraint)
        if constraint.is_vertex:
            self.vertex_constraints[constraint.agent][constraint.timestep].add(constraint.position)
        if constraint.is_edge:
            self.edge_constraints[constraint.agent][constraint.timestep].add(constraint.edge)

    def extend(self, constraints: Iterable[Constraint]) -> "ConstraintSet":
        result = self.copy()
        for constraint in constraints:
            result.add(constraint)
        return result

    def is_forbidden(
        self,
        agent: int,
        position: Tuple[int, int],
        timestep: int,
        prev_position: Optional[Tuple[int, int]] = None,
    ) -> bool:
        if position in self.vertex_constraints.get(agent, {}).get(timestep, set()):
            return True
        if prev_position is not None:
            edge = (prev_position, position)
            if edge in self.edge_constraints.get(agent, {}).get(timestep, set()):
                return True
        return False

    def count_future_conflicts(self, agent: int, position: Tuple[int, int], timestep: int) -> int:
        count = 0
        for t, blocked in self.vertex_constraints.get(agent, {}).items():
            if t >= timestep and position in blocked:
                count += 1
        return count

    def latest_timestep(self, agent: int) -> int:
        vertex_max = max(self.vertex_constraints.get(agent, {0: set()}).keys(), default=0)
        edge_max = max(self.edge_constraints.get(agent, {0: set()}).keys(), default=0)
        return max(vertex_max, edge_max)

    def is_goal_safe(self, agent: int, goal: Tuple[int, int], arrival_timestep: int, horizon: int) -> bool:
        for t in range(arrival_timestep, horizon + 1):
            if self.is_forbidden(agent, goal, t, goal if t > arrival_timestep else None):
                return False
        return True
