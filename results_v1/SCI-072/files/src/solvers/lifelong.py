import heapq
import time
import math
import random
import numpy as np
from typing import List, Dict, Optional, Tuple, Set
from collections import defaultdict, deque
import sys
from dataclasses import dataclass
from itertools import permutations

sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')
from src.core.environment import GridEnvironment
from src.core.agent import Agent
from src.core.conflict import ConflictDetector, Conflict, ConflictType
from src.core.constraint import Constraint, ConstraintSet
from src.core.solution import Solution, Path
from src.solvers.base import MAPFSolver


class TaskAssigner:
    """Task assignment utilities for lifelong MAPF."""

    def __init__(self, env: GridEnvironment) -> None:
        self.env = env

    def assign_nearest(self, agents: List[Agent], tasks: List[Tuple[int, int]]) -> Dict[int, Tuple[int, int]]:
        remaining = list(tasks)
        assignments: Dict[int, Tuple[int, int]] = {}
        for agent in sorted(agents, key=lambda item: item.id):
            if not remaining:
                break
            best = min(remaining, key=lambda task: self.env.heuristic(agent.start, task))
            assignments[agent.id] = best
            remaining.remove(best)
        return assignments

    def assign_hungarian(self, agents: List[Agent], tasks: List[Tuple[int, int]]) -> Dict[int, Tuple[int, int]]:
        if not agents or not tasks:
            return {}
        size = min(len(agents), len(tasks))
        agents = agents[:size]
        tasks = tasks[:size]
        try:
            from scipy.optimize import linear_sum_assignment

            cost_matrix = np.array([[self.env.heuristic(agent.start, task) for task in tasks] for agent in agents])
            rows, cols = linear_sum_assignment(cost_matrix)
            return {agents[row].id: tasks[col] for row, col in zip(rows, cols)}
        except Exception:
            best_cost = math.inf
            best_perm: Optional[Tuple[Tuple[int, int], ...]] = None
            for perm in permutations(tasks, len(agents)):
                cost = sum(self.env.heuristic(agent.start, task) for agent, task in zip(agents, perm))
                if cost < best_cost:
                    best_cost = cost
                    best_perm = perm
            if best_perm is None:
                return {}
            return {agent.id: task for agent, task in zip(agents, best_perm)}


class LifelongMAPF(MAPFSolver):
    """Rolling-horizon MAPF for continual task execution."""

    def __init__(
        self,
        env: GridEnvironment,
        agents: List[Agent],
        timeout: int = 300,
        tasks: Optional[List[Tuple[int, int]]] = None,
        window_size: int = 8,
        replanning_strategy: str = "FULL",
        assignment_strategy: str = "nearest",
        max_steps: int = 128,
        **kwargs,
    ) -> None:
        super().__init__(env, agents, timeout)
        self.tasks = list(tasks or [])
        self.window_size = max(1, window_size)
        self.replanning_strategy = replanning_strategy.upper()
        self.assignment_strategy = assignment_strategy.lower()
        self.max_steps = max_steps
        self.task_assigner = TaskAssigner(env)
        self.stats = {
            "total_tasks_completed": 0,
            "average_task_time": 0.0,
            "replanning_count": 0,
            "throughput": 0.0,
            "runtime": 0.0,
        }

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        pending_tasks = deque(self.tasks)
        current_positions = {agent.id: agent.start for agent in self.agents}
        current_goals = {agent.id: agent.goal for agent in self.agents}
        active_plans = {agent.id: [agent.start] for agent in self.agents}
        traces = {agent.id: [agent.start] for agent in self.agents}
        completed_flags = {agent.id: False for agent in self.agents}
        task_started = {agent.id: 0 for agent in self.agents}
        completion_times: List[int] = []

        if pending_tasks:
            initial_assignments = self._assign_tasks(self._agent_views(current_positions, current_goals), list(pending_tasks))
            for agent_id, task in initial_assignments.items():
                current_goals[agent_id] = task
                pending_tasks.remove(task)
                completed_flags[agent_id] = False
                task_started[agent_id] = 0

        try:
            for step in range(self.max_steps):
                self._check_timeout()
                window = self._current_window(current_positions)
                to_replan = self._agents_to_replan(current_positions, current_goals, active_plans, completed_flags)
                if to_replan:
                    replanned = self._plan_window(current_positions, current_goals, window, to_replan)
                    active_plans.update(replanned)
                    self.stats["replanning_count"] += 1

                for agent in self.agents:
                    path = active_plans.get(agent.id, [current_positions[agent.id]])
                    next_pos = path[1] if len(path) > 1 else path[0]
                    current_positions[agent.id] = next_pos
                    traces[agent.id].append(next_pos)
                    active_plans[agent.id] = path[1:] if len(path) > 1 else [next_pos]

                for agent in self.agents:
                    if current_positions[agent.id] != current_goals[agent.id] or completed_flags[agent.id]:
                        continue
                    completed_flags[agent.id] = True
                    self.stats["total_tasks_completed"] += 1
                    completion_times.append(step + 1 - task_started.get(agent.id, 0))
                    if pending_tasks:
                        next_task = self._assign_single(agent.id, current_positions[agent.id], pending_tasks)
                        if next_task is not None:
                            current_goals[agent.id] = next_task
                            pending_tasks.remove(next_task)
                            completed_flags[agent.id] = False
                            task_started[agent.id] = step + 1
                            active_plans[agent.id] = [current_positions[agent.id]]

                if not pending_tasks and all(completed_flags.values()):
                    break
        except TimeoutError:
            pass

        elapsed_steps = max((len(path) - 1 for path in traces.values()), default=1)
        self.stats["average_task_time"] = float(sum(completion_times) / max(1, len(completion_times)))
        self.stats["throughput"] = float(self.stats["total_tasks_completed"] / max(1, elapsed_steps))
        self.stats["runtime"] = time.time() - self.start_time
        return Solution({aid: Path(aid, path) for aid, path in traces.items()})

    def _assign_tasks(self, agents: List[Agent], tasks: List[Tuple[int, int]]) -> Dict[int, Tuple[int, int]]:
        if self.assignment_strategy == "hungarian":
            return self.task_assigner.assign_hungarian(agents, tasks)
        return self.task_assigner.assign_nearest(agents, tasks)

    def _assign_single(self, agent_id: int, position: Tuple[int, int], tasks: deque) -> Optional[Tuple[int, int]]:
        if not tasks:
            return None
        return min(tasks, key=lambda task: self.env.heuristic(position, task))

    def _current_window(self, current_positions: Dict[int, Tuple[int, int]]) -> int:
        if self.replanning_strategy != "ADAPTIVE":
            return self.window_size
        congestion = 0
        positions = list(current_positions.values())
        for index, position in enumerate(positions):
            for other in positions[index + 1:]:
                if self.env.heuristic(position, other) <= 2:
                    congestion += 1
        if congestion >= max(1, len(positions) // 2):
            return min(self.window_size + 4, self.env.width * self.env.height + self.window_size)
        return max(2, self.window_size - 1)

    def _agents_to_replan(
        self,
        current_positions: Dict[int, Tuple[int, int]],
        current_goals: Dict[int, Tuple[int, int]],
        active_plans: Dict[int, List[Tuple[int, int]]],
        completed_flags: Dict[int, bool],
    ) -> List[int]:
        if self.replanning_strategy == "FULL":
            return [agent.id for agent in self.agents]
        if self.replanning_strategy == "PARTIAL":
            selected = [agent.id for agent in self.agents if not completed_flags[agent.id] and len(active_plans.get(agent.id, [])) <= 2]
            future_paths = {aid: active_plans.get(aid, [current_positions[aid]])[:3] for aid in current_positions}
            for conflict in ConflictDetector.detect_conflicts(future_paths):
                selected.extend([conflict.agent1, conflict.agent2])
            return sorted(set(selected))
        return [agent.id for agent in self.agents]

    def _plan_window(
        self,
        current_positions: Dict[int, Tuple[int, int]],
        current_goals: Dict[int, Tuple[int, int]],
        window: int,
        selected_agents: List[int],
    ) -> Dict[int, List[Tuple[int, int]]]:
        reservations: Dict[int, List[Tuple[int, int]]] = {}
        plans: Dict[int, List[Tuple[int, int]]] = {}
        ordered = sorted(selected_agents, key=lambda aid: self.env.heuristic(current_positions[aid], current_goals[aid]), reverse=True)
        for agent_id in ordered:
            agent = next(agent for agent in self.agents if agent.id == agent_id)
            constraints = ConstraintSet()
            for reserved_path in reservations.values():
                for timestep in range(1, min(len(reserved_path), window + 1)):
                    constraints.add(Constraint(agent=agent_id, timestep=timestep, position=reserved_path[timestep]))
                    constraints.add(Constraint(agent=agent_id, timestep=timestep, edge=(reserved_path[timestep - 1], reserved_path[timestep])))
            horizon = window + self.env.width * self.env.height + 2
            path = self.a_star(current_positions[agent_id], current_goals[agent_id], constraints, agent=agent, max_timestep=horizon)
            if path is None:
                path = [current_positions[agent_id]] * (window + 1)
            elif len(path) < window + 1:
                path = path + [path[-1]] * (window + 1 - len(path))
            reservations[agent_id] = path
            plans[agent_id] = path
        return plans

    def _agent_views(self, current_positions: Dict[int, Tuple[int, int]], current_goals: Dict[int, Tuple[int, int]]) -> List[Agent]:
        return [
            Agent(id=agent.id, start=current_positions[agent.id], goal=current_goals[agent.id], name=agent.name)
            for agent in self.agents
        ]
