from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

import math
from pathlib import Path
import random

from .types import GridEnvironment, Position


class MapGenerator:
    """Factory methods for benchmark map generation."""

    @staticmethod
    def generate_empty(width: int, height: int) -> GridEnvironment:
        """Create an empty obstacle-free map."""
        return GridEnvironment(width=width, height=height)

    @staticmethod
    def generate_random(
        width: int,
        height: int,
        obstacle_ratio: float,
        seed: int | None = None,
    ) -> GridEnvironment:
        """Create a random obstacle map."""
        if not 0.0 <= obstacle_ratio < 1.0:
            raise ValueError('obstacle_ratio must be in [0, 1).')

        rng = random.Random(seed)
        cells = [(x, y) for y in range(height) for x in range(width)]
        obstacle_count = min(int(width * height * obstacle_ratio), max(0, len(cells) - 1))
        obstacles = set(rng.sample(cells, obstacle_count))
        env = GridEnvironment(width=width, height=height, obstacles=obstacles, metadata={'generator': 'random', 'seed': seed})
        if not env.free_cells() and obstacles:
            obstacles.pop()
            env = GridEnvironment(width=width, height=height, obstacles=obstacles, metadata={'generator': 'random', 'seed': seed})
        return env

    @staticmethod
    def generate_warehouse(
        width: int,
        height: int,
        aisle_width: int = 3,
        shelf_depth: int = 4,
        num_stations: int = 4,
    ) -> GridEnvironment:
        """Create a warehouse layout with shelves, aisles, and stations."""
        if width < 8 or height < 8:
            raise ValueError('Warehouse maps require both width and height to be at least 8.')
        if aisle_width < 1 or shelf_depth < 1 or num_stations < 1:
            raise ValueError('aisle_width, shelf_depth, and num_stations must be positive.')

        obstacles: set[Position] = set()
        shelves: set[Position] = set()
        cross_aisles = {1, max(2, height // 2), height - 2}
        x = 1
        while x + shelf_depth <= width - 1:
            for shelf_x in range(x, min(x + shelf_depth, width - 1)):
                for y in range(1, height - 1):
                    if y in cross_aisles:
                        continue
                    cell = (shelf_x, y)
                    shelves.add(cell)
                    obstacles.add(cell)
            x += shelf_depth + aisle_width

        station_y = height - 2
        if num_stations == 1:
            station_xs = [width // 2]
        else:
            gap = (width - 3) / (num_stations - 1)
            station_xs = [max(1, min(width - 2, round(1 + index * gap))) for index in range(num_stations)]
        stations = sorted({(station_x, station_y) for station_x in station_xs})

        env = GridEnvironment(
            width=width,
            height=height,
            obstacles=obstacles,
            stations=stations,
            shelves=shelves,
            metadata={'generator': 'warehouse'},
        )

        pickup_points: set[Position] = set()
        for shelf_x, shelf_y in shelves:
            for neighbor in ((shelf_x + 1, shelf_y), (shelf_x - 1, shelf_y), (shelf_x, shelf_y + 1), (shelf_x, shelf_y - 1)):
                if env.is_walkable(neighbor):
                    pickup_points.add(neighbor)
        env.metadata['warehouse_pickups'] = sorted(pickup_points)
        env.metadata['cross_aisles'] = sorted(cross_aisles)
        env.metadata['stations'] = stations
        env.metadata['bottlenecks'] = [cell for cell in pickup_points if cell[1] == height // 2]
        return env

    @staticmethod
    def generate_maze(width: int, height: int, seed: int | None = None) -> GridEnvironment:
        """Create a perfect maze with recursive backtracking."""
        rng = random.Random(seed)
        obstacles = {(x, y) for y in range(height) for x in range(width)}
        env = GridEnvironment(width=width, height=height, obstacles=obstacles, metadata={'generator': 'maze', 'seed': seed})
        if width < 3 or height < 3:
            env.obstacles.clear()
            return env

        start = (1, 1)
        stack = [start]
        visited = {start}
        env.obstacles.discard(start)

        def step_neighbors(cell: Position) -> list[Position]:
            x, y = cell
            candidates = [(x + 2, y), (x - 2, y), (x, y + 2), (x, y - 2)]
            return [candidate for candidate in candidates if 1 <= candidate[0] < width - 1 and 1 <= candidate[1] < height - 1]

        while stack:
            current = stack[-1]
            candidates = [candidate for candidate in step_neighbors(current) if candidate not in visited]
            if not candidates:
                stack.pop()
                continue
            nxt = rng.choice(candidates)
            wall = ((current[0] + nxt[0]) // 2, (current[1] + nxt[1]) // 2)
            env.obstacles.discard(wall)
            env.obstacles.discard(nxt)
            visited.add(nxt)
            stack.append(nxt)

        return env

    @staticmethod
    def generate_room(
        width: int,
        height: int,
        num_rooms: int = 4,
        door_width: int = 1,
    ) -> GridEnvironment:
        """Create a room map separated by walls with door openings."""
        if num_rooms < 1:
            raise ValueError('num_rooms must be at least 1.')
        if door_width < 1:
            raise ValueError('door_width must be at least 1.')

        obstacles: set[Position] = set()
        for x in range(width):
            obstacles.add((x, 0))
            obstacles.add((x, height - 1))
        for y in range(height):
            obstacles.add((0, y))
            obstacles.add((width - 1, y))

        rows = max(1, int(math.sqrt(num_rooms)))
        cols = max(1, math.ceil(num_rooms / rows))

        for row in range(1, rows):
            wall_y = row * height // rows
            doorway_x = max(1, (width // 2) - door_width // 2)
            doorway = set(range(doorway_x, min(width - 1, doorway_x + door_width)))
            for x in range(1, width - 1):
                if x not in doorway:
                    obstacles.add((x, wall_y))

        for col in range(1, cols):
            wall_x = col * width // cols
            doorway_y = max(1, (height // 2) - door_width // 2)
            doorway = set(range(doorway_y, min(height - 1, doorway_y + door_width)))
            for y in range(1, height - 1):
                if y not in doorway:
                    obstacles.add((wall_x, y))

        return GridEnvironment(
            width=width,
            height=height,
            obstacles=obstacles,
            metadata={'generator': 'room', 'num_rooms': num_rooms, 'door_width': door_width},
        )

    @staticmethod
    def save_map(env: GridEnvironment, filepath: str | Path) -> None:
        """Save a map in MovingAI `.map` format."""
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = ['type octile', f'height {env.height}', f'width {env.width}', 'map']
        for y in range(env.height):
            row: list[str] = []
            for x in range(env.width):
                cell = (x, y)
                if cell in env.shelves:
                    row.append('T')
                elif cell in env.obstacles:
                    row.append('@')
                elif cell in env.stations:
                    row.append('S')
                else:
                    row.append('.')
            lines.append(''.join(row))
        path.write_text('\\n'.join(lines) + '\\n', encoding='utf-8')

    @staticmethod
    def load_map(filepath: str | Path) -> GridEnvironment:
        """Load a MovingAI `.map` file."""
        path = Path(filepath)
        raw_lines = path.read_text(encoding='utf-8').splitlines()
        if len(raw_lines) < 4 or raw_lines[0].strip().lower() != 'type octile':
            raise ValueError('Unsupported or malformed MovingAI map file.')

        height = int(raw_lines[1].split()[1])
        width = int(raw_lines[2].split()[1])
        grid_lines = raw_lines[4:4 + height]
        if len(grid_lines) != height:
            raise ValueError('Map height does not match file contents.')

        obstacles: set[Position] = set()
        shelves: set[Position] = set()
        stations: list[Position] = []
        for y, row in enumerate(grid_lines):
            if len(row) != width:
                raise ValueError('Map width does not match file contents.')
            for x, char in enumerate(row):
                cell = (x, y)
                if char in {'@', 'T', '#'}:
                    obstacles.add(cell)
                if char == 'T':
                    shelves.add(cell)
                if char == 'S':
                    stations.append(cell)

        return GridEnvironment(
            width=width,
            height=height,
            obstacles=obstacles,
            stations=stations,
            shelves=shelves,
            metadata={'generator': 'loaded_map', 'source': str(path)},
        )
