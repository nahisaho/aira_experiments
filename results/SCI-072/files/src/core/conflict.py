from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple


class ConflictType(str, Enum):
    VERTEX = "vertex"
    EDGE = "edge"


@dataclass(frozen=True)
class Conflict:
    agent1: int
    agent2: int
    timestep: int
    conflict_type: ConflictType
    position: Tuple[int, int]
    position2: Optional[Tuple[int, int]] = None
    cardinality: str = "unknown"


class ConflictDetector:
    """Detect first or all pairwise conflicts in synchronized paths."""

    @staticmethod
    def _loc_at(path: List[Tuple[int, int]], timestep: int) -> Tuple[int, int]:
        if timestep < len(path):
            return path[timestep]
        return path[-1]

    @classmethod
    def detect_conflicts(cls, paths: Dict[int, List[Tuple[int, int]]]) -> List[Conflict]:
        conflicts: List[Conflict] = []
        agent_ids = sorted(paths)
        makespan = max((len(path) for path in paths.values()), default=0)
        for i, agent1 in enumerate(agent_ids):
            for agent2 in agent_ids[i + 1:]:
                for t in range(makespan):
                    loc1 = cls._loc_at(paths[agent1], t)
                    loc2 = cls._loc_at(paths[agent2], t)
                    if loc1 == loc2:
                        conflicts.append(Conflict(agent1, agent2, t, ConflictType.VERTEX, loc1))
                        break
                    if t == 0:
                        continue
                    prev1 = cls._loc_at(paths[agent1], t - 1)
                    prev2 = cls._loc_at(paths[agent2], t - 1)
                    if prev1 == loc2 and prev2 == loc1:
                        conflicts.append(Conflict(agent1, agent2, t, ConflictType.EDGE, prev1, loc1))
                        break
        return conflicts

    @classmethod
    def first_conflict(cls, paths: Dict[int, List[Tuple[int, int]]]) -> Optional[Conflict]:
        conflicts = cls.detect_conflicts(paths)
        if not conflicts:
            return None
        conflicts.sort(key=lambda c: c.timestep)
        return conflicts[0]
