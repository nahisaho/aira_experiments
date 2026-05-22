from collections import deque
from typing import Iterable, List, Optional, Set, Tuple


class GridEnvironment:
    """2-D 4-connected grid with wait action support."""

    def __init__(self, width: int, height: int, obstacles: Optional[Iterable[Tuple[int, int]]] = None):
        self.width = width
        self.height = height
        self.obstacles: Set[Tuple[int, int]] = set(obstacles or [])

    def in_bounds(self, loc: Tuple[int, int]) -> bool:
        x, y = loc
        return 0 <= x < self.width and 0 <= y < self.height

    def passable(self, loc: Tuple[int, int]) -> bool:
        return loc not in self.obstacles

    def valid(self, loc: Tuple[int, int]) -> bool:
        return self.in_bounds(loc) and self.passable(loc)

    def neighbors(self, loc: Tuple[int, int], include_wait: bool = False) -> List[Tuple[int, int]]:
        x, y = loc
        nbrs = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]
        if include_wait:
            nbrs.append(loc)
        return [n for n in nbrs if self.valid(n)]

    def heuristic(self, src: Tuple[int, int], dst: Tuple[int, int]) -> int:
        return abs(src[0] - dst[0]) + abs(src[1] - dst[1])

    def shortest_path_length(self, start: Tuple[int, int], goal: Tuple[int, int]) -> Optional[int]:
        if not self.valid(start) or not self.valid(goal):
            return None
        if start == goal:
            return 0
        queue = deque([(start, 0)])
        seen = {start}
        while queue:
            loc, dist = queue.popleft()
            for nxt in self.neighbors(loc, include_wait=False):
                if nxt in seen:
                    continue
                if nxt == goal:
                    return dist + 1
                seen.add(nxt)
                queue.append((nxt, dist + 1))
        return None
