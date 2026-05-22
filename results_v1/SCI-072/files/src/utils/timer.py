from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

import time


class Timer:
    """Wall-clock timer context manager."""

    def __init__(self) -> None:
        self.start_time: float | None = None
        self.elapsed = 0.0

    def __enter__(self) -> 'Timer':
        self.start_time = time.perf_counter()
        return self

    def __exit__(self, *args: object) -> None:
        if self.start_time is not None:
            self.elapsed = time.perf_counter() - self.start_time

    @property
    def seconds(self) -> float:
        if self.start_time and self.elapsed == 0:
            return time.perf_counter() - self.start_time
        return self.elapsed
