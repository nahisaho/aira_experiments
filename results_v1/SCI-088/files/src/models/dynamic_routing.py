"""
Dynamic Rerouting Engine
========================
Incident detection, traffic state estimation, and dynamic rerouting
for accident/construction scenarios.

Methods:
  - Speed anomaly detection for incident identification
  - A* with real-time travel time costs
  - Compliance-weighted route assignment
  - K-shortest paths for route diversity

References:
- Dia, H. (2001). An object-oriented neural network approach to
  short-term traffic forecasting. EJOR.
"""

import numpy as np
import heapq
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from enum import Enum
import logging
import time as time_module

logger = logging.getLogger(__name__)


class IncidentType(Enum):
    ACCIDENT = "accident"
    CONSTRUCTION = "construction"
    BREAKDOWN = "breakdown"
    WEATHER = "weather"


@dataclass
class Incident:
    """Represents a traffic incident."""
    id: str
    type: IncidentType
    link_id: str
    start_time: float
    estimated_duration: float      # seconds
    capacity_reduction: float      # 0.0 = full blockage, 1.0 = no effect
    detected_time: float = 0.0
    confirmed: bool = False
    cleared: bool = False


@dataclass
class LinkState:
    """Real-time state of a road link."""
    id: str
    length: float                  # meters
    free_flow_speed: float         # m/s
    current_speed: float           # m/s
    density: float                 # vehicles/km
    flow: float                    # vehicles/hour
    capacity: float                # vehicles/hour
    num_lanes: int = 2
    incident: Optional[Incident] = None

    @property
    def travel_time(self) -> float:
        """Current travel time in seconds."""
        speed = max(self.current_speed, 0.5)  # min 0.5 m/s
        return self.length / speed

    @property
    def free_flow_time(self) -> float:
        return self.length / self.free_flow_speed

    @property
    def congestion_index(self) -> float:
        """Travel time index: current_tt / free_flow_tt."""
        return self.travel_time / self.free_flow_time


class IncidentDetector:
    """Detect incidents from speed anomalies in real-time traffic data."""

    def __init__(
        self,
        speed_threshold_ratio: float = 0.3,
        confirmation_time: float = 120,
    ):
        self.speed_threshold_ratio = speed_threshold_ratio
        self.confirmation_time = confirmation_time
        self._candidates: Dict[str, Tuple[float, float]] = {}  # link_id -> (detect_time, speed)

    def detect(
        self, link_states: Dict[str, LinkState], current_time: float
    ) -> List[Incident]:
        """Check for speed anomalies indicating incidents.

        Detection logic:
        1. If speed < threshold_ratio * free_flow_speed → candidate
        2. If candidate persists for confirmation_time → confirmed incident
        """
        new_incidents = []

        for lid, state in link_states.items():
            threshold = state.free_flow_speed * self.speed_threshold_ratio

            if state.current_speed < threshold and state.incident is None:
                if lid not in self._candidates:
                    self._candidates[lid] = (current_time, state.current_speed)
                    logger.info(f"Incident candidate on {lid}: speed={state.current_speed:.1f} m/s")
                else:
                    detect_time, _ = self._candidates[lid]
                    if current_time - detect_time >= self.confirmation_time:
                        incident = Incident(
                            id=f"inc_{lid}_{int(current_time)}",
                            type=IncidentType.ACCIDENT,
                            link_id=lid,
                            start_time=detect_time,
                            estimated_duration=1800,  # default 30 min
                            capacity_reduction=0.3,
                            detected_time=current_time,
                            confirmed=True,
                        )
                        new_incidents.append(incident)
                        del self._candidates[lid]
                        logger.warning(f"Incident confirmed on {lid}")
            elif lid in self._candidates and state.current_speed >= threshold:
                del self._candidates[lid]

        return new_incidents


class DynamicRouter:
    """A* router with real-time travel time costs and K-shortest paths."""

    def __init__(self, compliance_rate: float = 0.7):
        self.compliance_rate = compliance_rate
        self.graph: Dict[str, Dict[str, str]] = {}     # node -> {node: link_id}
        self.link_states: Dict[str, LinkState] = {}

    def build_graph(
        self,
        links: List[Tuple[str, str, str, float, float]],
    ):
        """Build graph from link definitions.

        Args:
            links: [(link_id, from_node, to_node, length, free_flow_speed)]
        """
        for lid, from_n, to_n, length, ffs in links:
            if from_n not in self.graph:
                self.graph[from_n] = {}
            self.graph[from_n][to_n] = lid
            self.link_states[lid] = LinkState(
                id=lid, length=length, free_flow_speed=ffs,
                current_speed=ffs, density=0, flow=0,
                capacity=1800,
            )

    def update_link_state(self, link_id: str, speed: float, density: float):
        """Update real-time state of a link."""
        if link_id in self.link_states:
            state = self.link_states[link_id]
            state.current_speed = speed
            state.density = density
            state.flow = speed * density

    def apply_incident(self, incident: Incident):
        """Reduce capacity on incident link."""
        if incident.link_id in self.link_states:
            state = self.link_states[incident.link_id]
            state.incident = incident
            state.current_speed *= incident.capacity_reduction
            state.capacity *= incident.capacity_reduction
            logger.info(f"Applied incident on {incident.link_id}: "
                       f"capacity → {state.capacity:.0f} veh/h")

    def clear_incident(self, link_id: str):
        """Restore link to normal conditions."""
        if link_id in self.link_states:
            state = self.link_states[link_id]
            if state.incident:
                state.current_speed = state.free_flow_speed
                state.capacity = 1800
                state.incident = None
                logger.info(f"Cleared incident on {link_id}")

    def find_shortest_path(
        self, origin: str, destination: str
    ) -> Tuple[List[str], float]:
        """A* shortest path with real-time travel times as costs.

        Returns:
            (path_nodes, total_travel_time)
        """
        if origin not in self.graph:
            return [], float('inf')

        # Priority queue: (cost, node, path)
        pq = [(0.0, origin, [origin])]
        visited: Set[str] = set()

        while pq:
            cost, node, path = heapq.heappop(pq)

            if node == destination:
                return path, cost

            if node in visited:
                continue
            visited.add(node)

            for neighbor, link_id in self.graph.get(node, {}).items():
                if neighbor not in visited:
                    link = self.link_states[link_id]
                    edge_cost = link.travel_time
                    heapq.heappush(pq, (cost + edge_cost, neighbor, path + [neighbor]))

        return [], float('inf')

    def find_k_shortest_paths(
        self, origin: str, destination: str, k: int = 3
    ) -> List[Tuple[List[str], float]]:
        """Find K shortest paths using Yen's algorithm.

        Provides route diversity for distributed rerouting.
        """
        best_path, best_cost = self.find_shortest_path(origin, destination)
        if not best_path:
            return []

        A = [(best_path, best_cost)]
        B: List[Tuple[List[str], float]] = []

        for i in range(1, k):
            for j in range(len(A[-1][0]) - 1):
                spur_node = A[-1][0][j]
                root_path = A[-1][0][:j + 1]

                # Temporarily remove edges used by existing paths
                removed_edges = []
                for path, _ in A:
                    if path[:j + 1] == root_path and j + 1 < len(path):
                        next_node = path[j + 1]
                        if spur_node in self.graph and next_node in self.graph[spur_node]:
                            lid = self.graph[spur_node].pop(next_node)
                            removed_edges.append((spur_node, next_node, lid))

                spur_path, spur_cost = self.find_shortest_path(spur_node, destination)

                # Restore removed edges
                for fn, tn, lid in removed_edges:
                    self.graph[fn][tn] = lid

                if spur_path:
                    total_path = root_path[:-1] + spur_path
                    root_cost = sum(
                        self.link_states[self.graph[root_path[m]][root_path[m + 1]]].travel_time
                        for m in range(len(root_path) - 1)
                        if root_path[m] in self.graph and root_path[m + 1] in self.graph[root_path[m]]
                    )
                    total_cost = root_cost + spur_cost
                    candidate = (total_path, total_cost)
                    if candidate not in B and candidate not in A:
                        B.append(candidate)

            if not B:
                break
            B.sort(key=lambda x: x[1])
            A.append(B.pop(0))

        return A

    def reroute_vehicles(
        self,
        affected_od_pairs: List[Tuple[str, str]],
        num_vehicles_per_od: Dict[Tuple[str, str], int],
    ) -> Dict[Tuple[str, str], List[Tuple[List[str], int]]]:
        """Reroute vehicles affected by incidents.

        Distributes compliant vehicles across K-shortest paths.
        Non-compliant vehicles keep original route.
        """
        rerouting_plan = {}

        for od in affected_od_pairs:
            origin, dest = od
            total = num_vehicles_per_od.get(od, 0)
            compliant = int(total * self.compliance_rate)

            paths = self.find_k_shortest_paths(origin, dest, k=3)
            if not paths:
                continue

            # Distribute compliant vehicles across alternative paths
            # Weight inversely by travel time
            costs = [1.0 / max(c, 1.0) for _, c in paths]
            total_weight = sum(costs)
            allocations = []
            remaining = compliant
            for i, (path, cost) in enumerate(paths):
                if i == len(paths) - 1:
                    alloc = remaining
                else:
                    alloc = int(compliant * costs[i] / total_weight)
                    remaining -= alloc
                allocations.append((path, alloc))

            rerouting_plan[od] = allocations

        return rerouting_plan


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Build simple grid
    router = DynamicRouter(compliance_rate=0.7)
    links = []
    for r in range(4):
        for c in range(4):
            n = f"n_{r}_{c}"
            if c < 3:
                lid = f"link_{r}_{c}_E"
                links.append((lid, n, f"n_{r}_{c+1}", 400, 13.89))
            if r < 3:
                lid = f"link_{r}_{c}_S"
                links.append((lid, n, f"n_{r+1}_{c}", 375, 13.89))

    router.build_graph(links)
    print(f"Graph: {len(router.graph)} nodes, {len(router.link_states)} links")

    # Find path
    path, cost = router.find_shortest_path("n_0_0", "n_3_3")
    print(f"Shortest path: {' → '.join(path)}, cost={cost:.1f}s")

    # Simulate incident
    incident = Incident(
        id="inc_1", type=IncidentType.ACCIDENT,
        link_id="link_1_1_E", start_time=0,
        estimated_duration=1800, capacity_reduction=0.2,
        confirmed=True,
    )
    router.apply_incident(incident)

    path2, cost2 = router.find_shortest_path("n_0_0", "n_3_3")
    print(f"After incident: {' → '.join(path2)}, cost={cost2:.1f}s")

    # K-shortest paths
    kpaths = router.find_k_shortest_paths("n_0_0", "n_3_3", k=3)
    for i, (p, c) in enumerate(kpaths):
        print(f"  Path {i+1}: cost={c:.1f}s, hops={len(p)}")
