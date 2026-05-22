from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class Agent:
    """MAPF agent description shared by discrete and continuous solvers."""

    id: int
    start: Tuple[int, int]
    goal: Tuple[int, int]
    name: str = ""
    radius: float = 0.35
    priority: float = 0.0
