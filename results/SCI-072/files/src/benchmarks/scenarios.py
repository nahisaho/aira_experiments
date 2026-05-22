from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

import random

from .types import Agent, GridEnvironment, Position


class ScenarioGenerator:
    """Scenario generation helpers for benchmark experiments."""

    @staticmethod
    def generate_random(
        env: GridEnvironment,
        num_agents: int,
        seed: int | None = None,
    ) -> list[Agent]:
        """Generate random non-colliding start/goal pairs."""
        if num_agents < 1:
            raise ValueError('num_agents must be positive.')

        free_cells = env.free_cells()
        if len(free_cells) < num_agents * 2:
            raise ValueError('Not enough walkable cells for unique starts and goals.')

        rng = random.Random(seed)
        starts = rng.sample(free_cells, num_agents)
        remaining = [cell for cell in free_cells if cell not in starts]
        goals = rng.sample(remaining, num_agents)
        return [
            Agent(start=start, goal=goal, id=index, metadata={'scenario': 'random'})
            for index, (start, goal) in enumerate(zip(starts, goals))
        ]

    @staticmethod
    def generate_warehouse_tasks(
        env: GridEnvironment,
        num_agents: int,
        num_tasks_per_agent: int = 5,
    ) -> tuple[list[Agent], list[list[tuple[Position, Position]]]]:
        """Create pickup/drop-off task queues for warehouse agents."""
        if num_agents < 1 or num_tasks_per_agent < 1:
            raise ValueError('num_agents and num_tasks_per_agent must be positive.')

        stations = env.stations or env.metadata.get('stations') or env.free_cells()[:num_agents]
        pickup_points = env.metadata.get('warehouse_pickups') or ScenarioGenerator._adjacent_walkable_cells(env, env.shelves)
        if not stations or not pickup_points:
            raise ValueError('Warehouse tasks require stations and walkable shelf pickup points.')

        agents: list[Agent] = []
        tasks: list[list[tuple[Position, Position]]] = []
        for index in range(num_agents):
            home_station = stations[index % len(stations)]
            agents.append(
                Agent(
                    start=home_station,
                    goal=home_station,
                    id=index,
                    metadata={'scenario': 'warehouse', 'home_station': home_station},
                )
            )
            task_queue: list[tuple[Position, Position]] = []
            for task_index in range(num_tasks_per_agent):
                pickup = pickup_points[(index * num_tasks_per_agent + task_index) % len(pickup_points)]
                dropoff = stations[(index + task_index) % len(stations)]
                task_queue.append((pickup, dropoff))
            tasks.append(task_queue)
        return agents, tasks

    @staticmethod
    def generate_congested(
        env: GridEnvironment,
        num_agents: int,
        bottleneck_ratio: float = 0.3,
    ) -> list[Agent]:
        """Generate agents that must traverse opposite sides of a bottleneck region."""
        if num_agents < 1:
            raise ValueError('num_agents must be positive.')
        if not 0.0 < bottleneck_ratio < 1.0:
            raise ValueError('bottleneck_ratio must be in (0, 1).')

        free_cells = env.free_cells()
        if len(free_cells) < num_agents * 2:
            raise ValueError('Not enough walkable cells for a congested scenario.')

        if env.width >= env.height:
            band = max(1, int(env.width * bottleneck_ratio / 2))
            center = env.width // 2
            left = [cell for cell in free_cells if cell[0] < center - band]
            right = [cell for cell in free_cells if cell[0] > center + band]
            primary, secondary = left, right
        else:
            band = max(1, int(env.height * bottleneck_ratio / 2))
            center = env.height // 2
            top = [cell for cell in free_cells if cell[1] < center - band]
            bottom = [cell for cell in free_cells if cell[1] > center + band]
            primary, secondary = top, bottom

        if not primary or not secondary:
            return ScenarioGenerator.generate_random(env, num_agents, seed=int(bottleneck_ratio * 10_000))

        rng = random.Random(int(env.width * env.height * bottleneck_ratio))
        rng.shuffle(primary)
        rng.shuffle(secondary)

        agents: list[Agent] = []
        for index in range(num_agents):
            if index % 2 == 0:
                start = primary[(index // 2) % len(primary)]
                goal = secondary[(index // 2) % len(secondary)]
            else:
                start = secondary[(index // 2) % len(secondary)]
                goal = primary[(index // 2) % len(primary)]
            agents.append(Agent(start=start, goal=goal, id=index, metadata={'scenario': 'congested'}))
        return agents

    @staticmethod
    def _adjacent_walkable_cells(env: GridEnvironment, blocked_cells: set[Position]) -> list[Position]:
        """Return walkable cells adjacent to blocked cells."""
        cells: set[Position] = set()
        for x, y in blocked_cells:
            for candidate in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
                if env.is_walkable(candidate):
                    cells.add(candidate)
        return sorted(cells)
