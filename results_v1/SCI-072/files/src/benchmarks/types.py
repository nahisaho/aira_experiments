from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

from dataclasses import dataclass, field
import random
from typing import Any, Iterable, Iterator

from src.core.environment import GridEnvironment as CoreGridEnvironment

Position = tuple[int, int]


@dataclass(slots=True)
class Agent:
    """Benchmark-friendly agent model compatible with the solver package."""

    start: Position
    goal: Position
    id: int | str
    name: str = ''
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def agent_id(self) -> int | str:
        """Alias used by some benchmark utilities."""
        return self.id


class GridEnvironment(CoreGridEnvironment):
    """Grid environment extended with benchmark metadata."""

    def __init__(
        self,
        width: int,
        height: int,
        obstacles: Iterable[Position] | None = None,
        stations: Iterable[Position] | None = None,
        shelves: Iterable[Position] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(width=width, height=height, obstacles=obstacles)
        self.stations: list[Position] = [tuple(cell) for cell in stations or []]
        self.shelves: set[Position] = {tuple(cell) for cell in shelves or []}
        self.metadata: dict[str, Any] = dict(metadata or {})

    def is_walkable(self, position: Position) -> bool:
        """Return True when a position can be occupied by an agent."""
        return self.valid(position)

    def iter_free_cells(self) -> Iterator[Position]:
        """Iterate over all walkable cells."""
        for y in range(self.height):
            for x in range(self.width):
                cell = (x, y)
                if self.is_walkable(cell):
                    yield cell

    def free_cells(self) -> list[Position]:
        """Return all walkable cells as a list."""
        return list(self.iter_free_cells())

    def random_free_cell(self, rng: random.Random | None = None) -> Position:
        """Sample a walkable cell uniformly at random."""
        cells = self.free_cells()
        if not cells:
            raise ValueError('Environment contains no walkable cells.')
        return (rng or random).choice(cells)

    def clone(self) -> 'GridEnvironment':
        """Return a deep copy of the environment state."""
        return GridEnvironment(
            width=self.width,
            height=self.height,
            obstacles=set(self.obstacles),
            stations=list(self.stations),
            shelves=set(self.shelves),
            metadata=dict(self.metadata),
        )
