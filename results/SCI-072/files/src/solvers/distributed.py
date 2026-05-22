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
class PlanMessage:
    """Shared local plan subject to delay and packet loss."""

    sender: int
    receiver: int
    deliver_at: int
    path: List[Tuple[int, int]]


class DistributedMAPF(MAPFSolver):
    """Decentralized coordination under communication constraints."""

    def __init__(
        self,
        env: GridEnvironment,
        agents: List[Agent],
        timeout: int = 300,
        communication_radius: float = 5.0,
        message_delay: int = 1,
        packet_loss_rate: float = 0.0,
        max_iterations: int = 20,
        seed: int = 0,
        **kwargs,
    ) -> None:
        super().__init__(env, agents, timeout)
        self.communication_radius = float(communication_radius)
        self.message_delay = max(0, int(message_delay))
        self.packet_loss_rate = min(max(float(packet_loss_rate), 0.0), 1.0)
        self.max_iterations = max_iterations
        self._rng = random.Random(seed)
        np.random.seed(seed)
        self.stats = {
            "messages_sent": 0,
            "conflicts_resolved": 0,
            "iterations": 0,
            "convergence_time": 0.0,
        }

    def solve(self) -> Optional[Solution]:
        self.start_time = time.time()
        priorities = {agent.id: self._rng.random() for agent in self.agents}
        local_plans: Dict[int, List[Tuple[int, int]]] = {}
        for agent in self.agents:
            plan = self.a_star(agent.start, agent.goal, ConstraintSet(), agent=agent)
            if plan is None:
                return None
            local_plans[agent.id] = plan

        transit: List[PlanMessage] = []
        inbox: Dict[int, Dict[int, List[Tuple[int, int]]]] = defaultdict(dict)
        stable_rounds = 0
        agent_lookup = {agent.id: agent for agent in self.agents}

        try:
            for iteration in range(self.max_iterations):
                self._check_timeout()
                self.stats["iterations"] = iteration + 1
                changed = False

                for sender in self.agents:
                    for receiver in self.agents:
                        if sender.id == receiver.id:
                            continue
                        if self._distance(sender.start, receiver.start) > self.communication_radius:
                            continue
                        if self._rng.random() < self.packet_loss_rate:
                            continue
                        transit.append(
                            PlanMessage(
                                sender=sender.id,
                                receiver=receiver.id,
                                deliver_at=iteration + self.message_delay,
                                path=list(local_plans[sender.id]),
                            )
                        )
                        self.stats["messages_sent"] += 1

                delivered = [message for message in transit if message.deliver_at <= iteration]
                transit = [message for message in transit if message.deliver_at > iteration]
                for message in delivered:
                    inbox[message.receiver][message.sender] = message.path

                for agent in self.agents:
                    local_view = {agent.id: local_plans[agent.id], **inbox.get(agent.id, {})}
                    conflicts = ConflictDetector.detect_conflicts(local_view)
                    relevant = [conflict for conflict in conflicts if agent.id in (conflict.agent1, conflict.agent2)]
                    if not relevant:
                        continue
                    constraints = ConstraintSet()
                    for conflict in relevant:
                        winner = max((conflict.agent1, conflict.agent2), key=lambda aid: (priorities[aid], aid))
                        loser = conflict.agent2 if winner == conflict.agent1 else conflict.agent1
                        if loser != agent.id:
                            continue
                        self._reserve_path(agent.id, local_view[winner], constraints)
                    if not constraints.constraints:
                        continue
                    horizon = max(len(local_plans[agent.id]) + self.env.width * self.env.height, constraints.latest_timestep(agent.id) + 2)
                    replanned = self.a_star(agent.start, agent.goal, constraints, agent=agent_lookup[agent.id], max_timestep=horizon)
                    if replanned is not None and replanned != local_plans[agent.id]:
                        local_plans[agent.id] = replanned
                        self.stats["conflicts_resolved"] += 1
                        changed = True

                if changed:
                    stable_rounds = 0
                else:
                    stable_rounds += 1
                    if stable_rounds >= max(1, self.message_delay + 1):
                        break
        except TimeoutError:
            pass

        self.stats["convergence_time"] = time.time() - self.start_time
        return Solution({aid: Path(aid, path) for aid, path in local_plans.items()})

    def _reserve_path(self, agent_id: int, path: List[Tuple[int, int]], constraints: ConstraintSet) -> None:
        for timestep, loc in enumerate(path):
            constraints.add(Constraint(agent=agent_id, timestep=timestep, position=loc))
            if timestep > 0:
                constraints.add(Constraint(agent=agent_id, timestep=timestep, edge=(path[timestep - 1], loc)))

    def _distance(self, src: Tuple[int, int], dst: Tuple[int, int]) -> float:
        return math.dist(src, dst)
