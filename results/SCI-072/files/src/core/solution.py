from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass
class Path:
    agent: int
    states: List[Tuple[int, int]]

    @property
    def cost(self) -> int:
        last = len(self.states) - 1
        while last > 0 and self.states[last] == self.states[last - 1]:
            last -= 1
        return max(0, last)

    def at(self, timestep: int) -> Tuple[int, int]:
        if timestep < len(self.states):
            return self.states[timestep]
        return self.states[-1]


@dataclass
class Solution:
    paths: Dict[int, Path] = field(default_factory=dict)

    @property
    def cost(self) -> int:
        return sum(path.cost for path in self.paths.values())

    @property
    def makespan(self) -> int:
        return max((path.cost for path in self.paths.values()), default=0)
