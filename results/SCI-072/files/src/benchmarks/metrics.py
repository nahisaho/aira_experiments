from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

import resource
from typing import Any, Iterable, Sequence

from src.core.solution import Path as SolverPath
from src.core.solution import Solution as SolverSolution

from .types import Agent, GridEnvironment, Position


class MetricsCollector:
    """Compute summary metrics for benchmark runs."""

    @staticmethod
    def compute_all(
        solution: Any,
        env: GridEnvironment,
        agents: Sequence[Agent],
        solver_stats: dict[str, Any] | None,
    ) -> dict[str, Any]:
        """Return all benchmark metrics for a solver output."""
        stats = dict(solver_stats or {})
        paths = MetricsCollector._extract_paths(solution, agents)
        costs = [MetricsCollector._path_cost(path) for path in paths]
        conflicts_remaining = MetricsCollector._count_conflicts(paths)
        runtime = float(stats.get('runtime_seconds', 0.0) or 0.0)
        tasks_completed = int(stats.get('tasks_completed', 0) or 0)
        throughput = float(stats.get('throughput', 0.0) or 0.0)
        if throughput == 0.0 and runtime > 0.0 and tasks_completed > 0:
            throughput = tasks_completed / runtime

        assignment_times = stats.get('task_assignment_times', []) or []
        completion_times = stats.get('task_completion_times', []) or []
        service_time = MetricsCollector._compute_service_time(assignment_times, completion_times)

        memory_usage_mb = float(stats.get('memory_usage_mb', 0.0) or 0.0)
        if memory_usage_mb == 0.0:
            memory_usage_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

        optimal_cost = stats.get('optimal_cost')
        sum_of_costs = float(sum(costs))
        suboptimality_ratio = None
        if optimal_cost not in (None, 0):
            suboptimality_ratio = sum_of_costs / float(optimal_cost)

        success_flag = bool(stats.get('success', bool(paths)))
        valid = MetricsCollector._validate_paths(paths, env, agents)

        return {
            'sum_of_costs': sum_of_costs,
            'makespan': float(max(costs, default=0.0)),
            'runtime_seconds': runtime,
            'success': success_flag and valid and conflicts_remaining == 0,
            'nodes_expanded': int(stats.get('nodes_expanded', stats.get('expanded', 0)) or 0),
            'suboptimality_ratio': suboptimality_ratio,
            'throughput': throughput,
            'conflicts_remaining': conflicts_remaining,
            'memory_usage_mb': memory_usage_mb,
            'flowtime': sum_of_costs,
            'service_time': service_time,
        }

    @staticmethod
    def _extract_paths(solution: Any, agents: Sequence[Agent]) -> list[list[Position]]:
        if solution is None:
            return []
        if isinstance(solution, SolverSolution):
            return [MetricsCollector._normalize_path(path.states) for path in solution.paths.values()]
        if isinstance(solution, dict):
            if 'solution' in solution:
                return MetricsCollector._extract_paths(solution['solution'], agents)
            if 'paths' in solution:
                return MetricsCollector._normalize_paths(solution['paths'], agents)
            return MetricsCollector._normalize_paths(solution, agents)
        return MetricsCollector._normalize_paths(solution, agents)

    @staticmethod
    def _normalize_paths(raw_paths: Any, agents: Sequence[Agent]) -> list[list[Position]]:
        if isinstance(raw_paths, dict):
            paths: list[list[Position]] = []
            for index, agent in enumerate(agents):
                key_candidates = [agent.id, agent.agent_id, index, str(agent.id), str(index)]
                for key in key_candidates:
                    if key in raw_paths:
                        paths.append(MetricsCollector._normalize_path(raw_paths[key]))
                        break
            if paths:
                return paths
            return [MetricsCollector._normalize_path(path) for path in raw_paths.values()]
        return [MetricsCollector._normalize_path(path) for path in list(raw_paths)]

    @staticmethod
    def _normalize_path(path: Iterable[Any]) -> list[Position]:
        if isinstance(path, SolverPath):
            return [tuple(step) for step in path.states]
        normalized: list[Position] = []
        for step in path:
            if isinstance(step, dict):
                normalized.append((int(step['x']), int(step['y'])))
            else:
                x, y = step
                normalized.append((int(x), int(y)))
        return normalized

    @staticmethod
    def _path_cost(path: Sequence[Position]) -> int:
        return max(0, len(path) - 1)

    @staticmethod
    def _validate_paths(paths: Sequence[Sequence[Position]], env: GridEnvironment, agents: Sequence[Agent]) -> bool:
        if len(paths) != len(agents):
            return False
        for agent, path in zip(agents, paths):
            if not path:
                return False
            if tuple(path[0]) != tuple(agent.start):
                return False
            if tuple(path[-1]) != tuple(agent.goal):
                return False
            for cell in path:
                if not env.is_walkable(tuple(cell)):
                    return False
        return True

    @staticmethod
    def _count_conflicts(paths: Sequence[Sequence[Position]]) -> int:
        if not paths:
            return 0
        makespan = max(len(path) for path in paths)
        conflicts = 0
        for timestep in range(makespan):
            positions: dict[Position, int] = {}
            for path in paths:
                position = path[min(timestep, len(path) - 1)]
                positions[position] = positions.get(position, 0) + 1
            conflicts += sum(count - 1 for count in positions.values() if count > 1)

            if timestep == 0:
                continue
            seen_pairs: set[frozenset[Position]] = set()
            for i, path_i in enumerate(paths):
                prev_i = path_i[min(timestep - 1, len(path_i) - 1)]
                curr_i = path_i[min(timestep, len(path_i) - 1)]
                for j in range(i + 1, len(paths)):
                    path_j = paths[j]
                    prev_j = path_j[min(timestep - 1, len(path_j) - 1)]
                    curr_j = path_j[min(timestep, len(path_j) - 1)]
                    pair = frozenset({(prev_i, curr_i), (prev_j, curr_j)})
                    if pair in seen_pairs:
                        continue
                    if prev_i == curr_j and curr_i == prev_j and prev_i != curr_i:
                        conflicts += 1
                    seen_pairs.add(pair)
        return conflicts

    @staticmethod
    def _compute_service_time(
        assignment_times: Sequence[float],
        completion_times: Sequence[float],
    ) -> float | None:
        if not assignment_times or not completion_times:
            return None
        durations = [
            float(completion) - float(assignment)
            for assignment, completion in zip(assignment_times, completion_times)
            if float(completion) >= float(assignment)
        ]
        if not durations:
            return None
        return sum(durations) / len(durations)
