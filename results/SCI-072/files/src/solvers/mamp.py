import heapq
import time
import math
import random
import numpy as np
from typing import List, Dict, Optional, Tuple, Set
from collections import defaultdict, deque
import sys
from dataclasses import dataclass

sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')
from src.core.environment import GridEnvironment
from src.core.agent import Agent
from src.core.conflict import ConflictDetector, Conflict, ConflictType
from src.core.constraint import Constraint, ConstraintSet
from src.core.solution import Solution, Path
from src.solvers.base import MAPFSolver


@dataclass
class ContinuousState:
    """Continuous 2D pose with translational velocity and heading."""

    x: float
    y: float
    vx: float = 0.0
    vy: float = 0.0
    theta: float = 0.0

    def point(self) -> Tuple[float, float]:
        return (self.x, self.y)


@dataclass
class KinematicModel:
    """Simple bounded-velocity point-mass approximation."""

    max_velocity: float = 1.0
    max_acceleration: float = 1.0
    max_omega: float = math.pi
    radius: float = 0.35


@dataclass
class ContinuousConstraint:
    """Continuous CBS-style spatio-temporal forbidden disk."""

    agent_id: int
    time_index: int
    center: Tuple[float, float]
    radius: float


@dataclass
class RRTNode:
    state: ContinuousState
    parent: Optional[int]
    cost: float


class RRTStar:
    """Single-agent continuous-space planner."""

    def __init__(
        self,
        env: GridEnvironment,
        model: KinematicModel,
        step_size: float = 0.8,
        max_iterations: int = 600,
        goal_bias: float = 0.15,
        seed: int = 0,
    ) -> None:
        self.env = env
        self.model = model
        self.step_size = step_size
        self.max_iterations = max_iterations
        self.goal_bias = goal_bias
        self.rng = random.Random(seed)
        self.obstacles = [((x + 0.5), (y + 0.5), 0.5) for x, y in env.obstacles]

    def plan(
        self,
        start: ContinuousState,
        goal: ContinuousState,
        constraints: Optional[List[ContinuousConstraint]] = None,
    ) -> Optional[List[ContinuousState]]:
        constraints = constraints or []
        nodes = [RRTNode(start, None, 0.0)]
        goal_index: Optional[int] = None
        for _ in range(self.max_iterations):
            sample = goal if self.rng.random() < self.goal_bias else self._sample()
            nearest_index = min(range(len(nodes)), key=lambda idx: self._distance(nodes[idx].state, sample))
            new_state = self._steer(nodes[nearest_index].state, sample)
            new_cost = nodes[nearest_index].cost + self._distance(nodes[nearest_index].state, new_state)
            if not self._collision_free(nodes[nearest_index].state, new_state, new_cost, constraints):
                continue
            near = self._near(nodes, new_state)
            parent_index = nearest_index
            parent_cost = new_cost
            for idx in near:
                candidate_cost = nodes[idx].cost + self._distance(nodes[idx].state, new_state)
                if candidate_cost < parent_cost and self._collision_free(nodes[idx].state, new_state, candidate_cost, constraints):
                    parent_index = idx
                    parent_cost = candidate_cost
            nodes.append(RRTNode(new_state, parent_index, parent_cost))
            new_index = len(nodes) - 1
            for idx in near:
                rewired_cost = nodes[new_index].cost + self._distance(nodes[new_index].state, nodes[idx].state)
                if rewired_cost < nodes[idx].cost and self._collision_free(nodes[new_index].state, nodes[idx].state, rewired_cost, constraints):
                    nodes[idx].parent = new_index
                    nodes[idx].cost = rewired_cost
            if self._distance(new_state, goal) <= self.step_size and self._collision_free(new_state, goal, nodes[new_index].cost + self._distance(new_state, goal), constraints):
                nodes.append(RRTNode(goal, new_index, nodes[new_index].cost + self._distance(new_state, goal)))
                goal_index = len(nodes) - 1
                break
        if goal_index is None:
            if self._collision_free(start, goal, self._distance(start, goal), constraints):
                return [start, goal]
            return None
        return self._extract(nodes, goal_index)

    def _sample(self) -> ContinuousState:
        return ContinuousState(
            x=self.rng.uniform(0.5, max(0.5, self.env.width - 0.5)),
            y=self.rng.uniform(0.5, max(0.5, self.env.height - 0.5)),
        )

    def _steer(self, source: ContinuousState, target: ContinuousState) -> ContinuousState:
        dx = target.x - source.x
        dy = target.y - source.y
        dist = max(math.hypot(dx, dy), 1e-9)
        scale = min(self.step_size, dist) / dist
        theta = math.atan2(dy, dx)
        return ContinuousState(source.x + dx * scale, source.y + dy * scale, theta=theta)

    def _near(self, nodes: List[RRTNode], target: ContinuousState) -> List[int]:
        radius = self.step_size * 2.5
        return [idx for idx, node in enumerate(nodes) if self._distance(node.state, target) <= radius]

    def _collision_free(
        self,
        source: ContinuousState,
        target: ContinuousState,
        target_cost: float,
        constraints: List[ContinuousConstraint],
    ) -> bool:
        samples = max(2, int(math.ceil(self._distance(source, target) / max(self.step_size / 4.0, 0.1))))
        for step in range(samples + 1):
            alpha = step / max(1, samples)
            x = source.x + (target.x - source.x) * alpha
            y = source.y + (target.y - source.y) * alpha
            if not (0.0 <= x <= self.env.width and 0.0 <= y <= self.env.height):
                return False
            for ox, oy, radius in self.obstacles:
                if math.dist((x, y), (ox, oy)) <= radius + self.model.radius:
                    return False
            approx_time = int(round((target_cost * alpha) / max(self.model.max_velocity, 1e-6)))
            for constraint in constraints:
                if abs(approx_time - constraint.time_index) <= 1 and math.dist((x, y), constraint.center) <= constraint.radius:
                    return False
        return True

    def _extract(self, nodes: List[RRTNode], index: int) -> List[ContinuousState]:
        path: List[ContinuousState] = []
        current = index
        while current is not None:
            path.append(nodes[current].state)
            current = nodes[current].parent
        path.reverse()
        return path

    def _distance(self, src: ContinuousState, dst: ContinuousState) -> float:
        return math.dist(src.point(), dst.point())


class TrajectoryOptimizer:
    """Shortcut smoother plus simple velocity profiling."""

    def __init__(self, planner: RRTStar, dt: float = 0.5) -> None:
        self.planner = planner
        self.dt = dt

    def optimize(
        self,
        path: List[ContinuousState],
        constraints: Optional[List[ContinuousConstraint]] = None,
    ) -> List[ContinuousState]:
        constraints = constraints or []
        if len(path) <= 2:
            return self.time_parameterize(path)
        shortened = [path[0]]
        index = 0
        while index < len(path) - 1:
            candidate = len(path) - 1
            while candidate > index + 1:
                if self.planner._collision_free(path[index], path[candidate], candidate - index, constraints):
                    break
                candidate -= 1
            shortened.append(path[candidate])
            index = candidate
        return self.time_parameterize(shortened)

    def time_parameterize(self, path: List[ContinuousState]) -> List[ContinuousState]:
        if not path:
            return []
        trajectory = [path[0]]
        for source, target in zip(path, path[1:]):
            distance = math.dist(source.point(), target.point())
            steps = max(1, int(math.ceil(distance / max(self.planner.model.max_velocity * self.dt, 1e-6))))
            for step in range(1, steps + 1):
                alpha = step / steps
                x = source.x + (target.x - source.x) * alpha
                y = source.y + (target.y - source.y) * alpha
                theta = math.atan2(target.y - source.y, target.x - source.x)
                vx = (target.x - source.x) / max(steps * self.dt, 1e-6)
                vy = (target.y - source.y) / max(steps * self.dt, 1e-6)
                trajectory.append(ContinuousState(x, y, vx, vy, theta))
        return trajectory

    def smoothness(self, trajectory: List[ContinuousState]) -> float:
        if len(trajectory) < 3:
            return 0.0
        changes = []
        for prev_state, state, next_state in zip(trajectory, trajectory[1:], trajectory[2:]):
            heading1 = math.atan2(state.y - prev_state.y, state.x - prev_state.x)
            heading2 = math.atan2(next_state.y - state.y, next_state.x - state.x)
            changes.append(abs(heading2 - heading1))
        return float(sum(changes) / len(changes))


@dataclass
class MotionNode:
    constraints: Dict[int, List[ContinuousConstraint]]
    trajectories: Dict[int, List[ContinuousState]]
    cost: float
    node_id: int


class MAMP(MAPFSolver):
    """Continuous-space multi-agent motion planning via RRT* + CBS."""

    def __init__(
        self,
        env: GridEnvironment,
        agents: List[Agent],
        timeout: int = 300,
        kinematics: Optional[KinematicModel] = None,
        dt: float = 0.5,
        max_cbs_expansions: int = 128,
        seed: int = 0,
        **kwargs,
    ) -> None:
        super().__init__(env, agents, timeout)
        self.kinematics = kinematics or KinematicModel()
        self.planner = RRTStar(env, self.kinematics, seed=seed)
        self.optimizer = TrajectoryOptimizer(self.planner, dt)
        self.max_cbs_expansions = max_cbs_expansions
        self.stats = {
            "planning_time": 0.0,
            "trajectory_smoothness": 0.0,
            "min_clearance": math.inf,
        }

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        agent_lookup = {agent.id: agent for agent in self.agents}
        root_trajectories = self._plan_trajectories(defaultdict(list))
        if root_trajectories is None:
            self.stats["planning_time"] = time.time() - self.start_time
            return None
        next_id = 0
        open_heap: List[Tuple[float, int, MotionNode]] = []
        root = MotionNode(defaultdict(list), root_trajectories, self._trajectory_cost(root_trajectories), next_id)
        heapq.heappush(open_heap, (root.cost, root.node_id, root))

        try:
            while open_heap and next_id < self.max_cbs_expansions:
                self._check_timeout()
                _, _, node = heapq.heappop(open_heap)
                collision = self._first_collision(node.trajectories, agent_lookup)
                if collision is None:
                    solution = Solution({aid: self._trajectory_to_path(aid, traj) for aid, traj in node.trajectories.items()})
                    self.stats["planning_time"] = time.time() - self.start_time
                    smoothness = [self.optimizer.smoothness(traj) for traj in node.trajectories.values()]
                    self.stats["trajectory_smoothness"] = float(np.mean(smoothness)) if smoothness else 0.0
                    self.stats["min_clearance"] = self._minimum_clearance(node.trajectories, agent_lookup)
                    return solution
                for agent_id in (collision.agent1, collision.agent2):
                    constraints = defaultdict(list, {aid: list(values) for aid, values in node.constraints.items()})
                    constraints[agent_id].append(
                        ContinuousConstraint(
                            agent_id=agent_id,
                            time_index=collision.timestep,
                            center=((collision.position[0] + 0.5), (collision.position[1] + 0.5)),
                            radius=agent_lookup[collision.agent1].radius + agent_lookup[collision.agent2].radius,
                        )
                    )
                    replanned = self._plan_trajectories(constraints, node.trajectories, only_agent=agent_id)
                    if replanned is None:
                        continue
                    candidates = [replanned]
                    if self._first_collision(replanned, agent_lookup) is not None:
                        delayed = dict(replanned)
                        delayed[agent_id] = self._delay_trajectory(replanned[agent_id], max(1, collision.timestep // 2 + 1))
                        candidates.append(delayed)
                        detoured = dict(replanned)
                        detour_path = self._detour_trajectory(agent_lookup[agent_id], collision.position)
                        if detour_path is not None:
                            detoured[agent_id] = detour_path
                            candidates.append(detoured)
                    for candidate in candidates:
                        if self._first_collision(candidate, agent_lookup) is None:
                            solution = Solution({aid: self._trajectory_to_path(aid, traj) for aid, traj in candidate.items()})
                            self.stats["planning_time"] = time.time() - self.start_time
                            smoothness = [self.optimizer.smoothness(traj) for traj in candidate.values()]
                            self.stats["trajectory_smoothness"] = float(np.mean(smoothness)) if smoothness else 0.0
                            self.stats["min_clearance"] = self._minimum_clearance(candidate, agent_lookup)
                            return solution
                        next_id += 1
                        child = MotionNode(constraints, candidate, self._trajectory_cost(candidate), next_id)
                        heapq.heappush(open_heap, (child.cost, child.node_id, child))
        except TimeoutError:
            pass

        self.stats["planning_time"] = time.time() - self.start_time
        return None

    def _plan_trajectories(
        self,
        constraints: Dict[int, List[ContinuousConstraint]],
        current: Optional[Dict[int, List[ContinuousState]]] = None,
        only_agent: Optional[int] = None,
    ) -> Optional[Dict[int, List[ContinuousState]]]:
        trajectories = dict(current or {})
        for agent in self.agents:
            if only_agent is not None and agent.id != only_agent:
                continue
            start = self._cell_to_state(agent.start)
            goal = self._cell_to_state(agent.goal)
            geometric = self.planner.plan(start, goal, constraints.get(agent.id, []))
            if geometric is None:
                return None
            trajectories[agent.id] = self.optimizer.optimize(geometric, constraints.get(agent.id, []))
        return trajectories

    def _cell_to_state(self, cell: Tuple[int, int]) -> ContinuousState:
        return ContinuousState(cell[0] + 0.5, cell[1] + 0.5)

    def _trajectory_cost(self, trajectories: Dict[int, List[ContinuousState]]) -> float:
        return float(sum(max(0, len(traj) - 1) for traj in trajectories.values()))

    def _first_collision(
        self,
        trajectories: Dict[int, List[ContinuousState]],
        agents: Dict[int, Agent],
    ) -> Optional[Conflict]:
        agent_ids = sorted(trajectories)
        horizon = max((len(path) for path in trajectories.values()), default=0)
        for index, agent1 in enumerate(agent_ids):
            for agent2 in agent_ids[index + 1:]:
                for timestep in range(horizon):
                    state1 = trajectories[agent1][min(timestep, len(trajectories[agent1]) - 1)]
                    state2 = trajectories[agent2][min(timestep, len(trajectories[agent2]) - 1)]
                    if math.dist(state1.point(), state2.point()) <= agents[agent1].radius + agents[agent2].radius:
                        midpoint = ((state1.x + state2.x) / 2.0, (state1.y + state2.y) / 2.0)
                        return Conflict(
                            agent1=agent1,
                            agent2=agent2,
                            timestep=timestep,
                            conflict_type=ConflictType.VERTEX,
                            position=(int(midpoint[0]), int(midpoint[1])),
                        )
        return None

    def _minimum_clearance(self, trajectories: Dict[int, List[ContinuousState]], agents: Dict[int, Agent]) -> float:
        minimum = math.inf
        agent_ids = sorted(trajectories)
        horizon = max((len(path) for path in trajectories.values()), default=0)
        for index, agent1 in enumerate(agent_ids):
            for agent2 in agent_ids[index + 1:]:
                for timestep in range(horizon):
                    state1 = trajectories[agent1][min(timestep, len(trajectories[agent1]) - 1)]
                    state2 = trajectories[agent2][min(timestep, len(trajectories[agent2]) - 1)]
                    clearance = math.dist(state1.point(), state2.point()) - (agents[agent1].radius + agents[agent2].radius)
                    minimum = min(minimum, clearance)
        return 0.0 if math.isinf(minimum) else minimum

    def _delay_trajectory(self, trajectory: List[ContinuousState], wait_steps: int) -> List[ContinuousState]:
        if not trajectory:
            return []
        return [trajectory[0]] * wait_steps + list(trajectory)

    def _detour_trajectory(
        self,
        agent: Agent,
        conflict_cell: Tuple[int, int],
    ) -> Optional[List[ContinuousState]]:
        start = self._cell_to_state(agent.start)
        goal = self._cell_to_state(agent.goal)
        waypoint_candidates = [
            (conflict_cell[0], conflict_cell[1] + 1),
            (conflict_cell[0], conflict_cell[1] - 1),
            (conflict_cell[0] + 1, conflict_cell[1] + 1),
            (conflict_cell[0] - 1, conflict_cell[1] + 1),
            (conflict_cell[0] + 1, conflict_cell[1] - 1),
            (conflict_cell[0] - 1, conflict_cell[1] - 1),
            (conflict_cell[0], conflict_cell[1] + 2),
            (conflict_cell[0], conflict_cell[1] - 2),
            (conflict_cell[0] + 2, conflict_cell[1]),
            (conflict_cell[0] - 2, conflict_cell[1]),
        ]
        for waypoint_cell in waypoint_candidates:
            if not self.env.valid(waypoint_cell):
                continue
            waypoint = self._cell_to_state(waypoint_cell)
            first = self.planner.plan(start, waypoint, [])
            second = self.planner.plan(waypoint, goal, [])
            if first is None or second is None:
                continue
            stitched = first[:-1] + second
            return self.optimizer.time_parameterize(stitched)
        return None

    def _trajectory_to_path(self, agent_id: int, trajectory: List[ContinuousState]) -> Path:
        states: List[Tuple[int, int]] = []
        for state in trajectory:
            cell = (
                min(max(int(state.x), 0), self.env.width - 1),
                min(max(int(state.y), 0), self.env.height - 1),
            )
            if not states or states[-1] != cell:
                states.append(cell)
        return Path(agent_id, states or [(0, 0)])
