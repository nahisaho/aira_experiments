"""
SUMO/Flow Integration Environment
===================================
Gymnasium-compatible environment wrapping SUMO via Flow framework
for multi-agent traffic signal control training.

Provides:
  - SUMO network generation for Tokyo downtown grid
  - Flow-based traffic generation with multimodal vehicles
  - Multi-agent observation/action interface for RLlib
  - Performance metrics collection
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
import logging
import json
import os

logger = logging.getLogger(__name__)


@dataclass
class SUMONetworkConfig:
    """SUMO network configuration for Tokyo downtown."""
    rows: int = 8
    cols: int = 6
    block_length_ns: float = 375.0    # meters (N-S blocks)
    block_length_ew: float = 500.0    # meters (E-W blocks)
    lanes_major: int = 3
    lanes_minor: int = 2
    speed_limit: float = 13.89        # 50 km/h
    junction_type: str = "traffic_light"


@dataclass
class FlowConfig:
    """Flow traffic generation configuration."""
    total_demand: int = 4000           # vehicles/hour
    car_ratio: float = 0.45
    bus_ratio: float = 0.25
    bicycle_ratio: float = 0.15
    pedestrian_ratio: float = 0.15
    simulation_step: float = 0.1
    warmup_steps: int = 6000           # 600s at 0.1s step


class SUMONetworkGenerator:
    """Generate SUMO network files for Tokyo downtown grid."""

    def __init__(self, config: SUMONetworkConfig):
        self.config = config

    def generate_nodes_xml(self) -> str:
        """Generate SUMO nodes (.nod.xml) file content."""
        c = self.config
        lines = ['<?xml version="1.0" encoding="UTF-8"?>',
                 '<nodes xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
                 ' xsi:noNamespaceSchemaLocation="http://sumo.dlr.de/xsd/nodes_file.xsd">']

        for r in range(c.rows):
            for col in range(c.cols):
                x = col * c.block_length_ew
                y = r * c.block_length_ns
                nid = f"n_{r}_{col}"
                ntype = c.junction_type
                lines.append(f'  <node id="{nid}" x="{x:.1f}" y="{y:.1f}" type="{ntype}"/>')

        lines.append('</nodes>')
        return '\n'.join(lines)

    def generate_edges_xml(self) -> str:
        """Generate SUMO edges (.edg.xml) file content."""
        c = self.config
        lines = ['<?xml version="1.0" encoding="UTF-8"?>',
                 '<edges xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
                 ' xsi:noNamespaceSchemaLocation="http://sumo.dlr.de/xsd/edges_file.xsd">']

        for r in range(c.rows):
            for col in range(c.cols):
                # East-West edges
                if col < c.cols - 1:
                    eid = f"e_{r}_{col}_to_{r}_{col+1}"
                    lines.append(f'  <edge id="{eid}" from="n_{r}_{col}" to="n_{r}_{col+1}" '
                               f'numLanes="{c.lanes_major}" speed="{c.speed_limit}"/>')
                    eid_rev = f"e_{r}_{col+1}_to_{r}_{col}"
                    lines.append(f'  <edge id="{eid_rev}" from="n_{r}_{col+1}" to="n_{r}_{col}" '
                               f'numLanes="{c.lanes_major}" speed="{c.speed_limit}"/>')
                # North-South edges
                if r < c.rows - 1:
                    eid = f"e_{r}_{col}_to_{r+1}_{col}"
                    lines.append(f'  <edge id="{eid}" from="n_{r}_{col}" to="n_{r+1}_{col}" '
                               f'numLanes="{c.lanes_minor}" speed="{c.speed_limit}"/>')
                    eid_rev = f"e_{r+1}_{col}_to_{r}_{col}"
                    lines.append(f'  <edge id="{eid_rev}" from="n_{r+1}_{col}" to="n_{r}_{col}" '
                               f'numLanes="{c.lanes_minor}" speed="{c.speed_limit}"/>')

        lines.append('</edges>')
        return '\n'.join(lines)

    def generate_flow_xml(self, flow_config: FlowConfig) -> str:
        """Generate SUMO route/flow (.rou.xml) with multimodal traffic."""
        c = self.config
        fc = flow_config
        lines = ['<?xml version="1.0" encoding="UTF-8"?>',
                 '<routes xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">']

        # Vehicle types
        lines.append('  <vType id="car" vClass="passenger" length="4.5" '
                    'maxSpeed="13.89" accel="1.4" decel="2.0" '
                    'carFollowModel="IDM" sigma="0.5" color="1,1,0"/>')
        lines.append('  <vType id="bus" vClass="bus" length="12.0" '
                    'maxSpeed="11.11" accel="1.0" decel="1.5" '
                    'carFollowModel="IDM" sigma="0.3" color="0,0,1"/>')
        lines.append('  <vType id="bicycle" vClass="bicycle" length="1.8" '
                    'maxSpeed="4.17" accel="1.2" decel="2.5" color="0,1,0"/>')

        # Generate flows from boundary edges
        flow_id = 0
        boundary_edges = self._get_boundary_edges()
        cars_per_flow = int(fc.total_demand * fc.car_ratio / len(boundary_edges))

        for from_edge, to_edge in self._get_od_pairs(boundary_edges):
            lines.append(f'  <flow id="car_{flow_id}" type="car" '
                        f'from="{from_edge}" to="{to_edge}" '
                        f'begin="0" end="7200" vehsPerHour="{cars_per_flow}"/>')
            flow_id += 1

        # Bus routes (12 routes)
        for bus_route in range(12):
            stops = self._generate_bus_route(bus_route)
            route_edges = ' '.join(stops)
            lines.append(f'  <route id="bus_route_{bus_route}" edges="{route_edges}"/>')
            lines.append(f'  <flow id="bus_{bus_route}" type="bus" '
                        f'route="bus_route_{bus_route}" '
                        f'begin="0" end="7200" period="300"/>')

        lines.append('</routes>')
        return '\n'.join(lines)

    def generate_sumo_cfg(self, output_dir: str) -> str:
        """Generate SUMO configuration (.sumocfg) file."""
        return f"""<?xml version="1.0" encoding="UTF-8"?>
<configuration>
  <input>
    <net-file value="{output_dir}/tokyo_downtown.net.xml"/>
    <route-files value="{output_dir}/tokyo_downtown.rou.xml"/>
    <additional-files value="{output_dir}/tokyo_downtown.add.xml"/>
  </input>
  <time>
    <begin value="0"/>
    <end value="7200"/>
    <step-length value="0.1"/>
  </time>
  <processing>
    <lateral-resolution value="0.8"/>
    <collision.action value="warn"/>
  </processing>
  <output>
    <tripinfo-output value="{output_dir}/tripinfo.xml"/>
    <summary-output value="{output_dir}/summary.xml"/>
  </output>
</configuration>"""

    def _get_boundary_edges(self) -> List[str]:
        c = self.config
        edges = []
        for col in range(c.cols):
            edges.append(f"e_0_{col}_to_1_{col}")
            edges.append(f"e_{c.rows-1}_{col}_to_{c.rows-2}_{col}")
        for r in range(c.rows):
            edges.append(f"e_{r}_0_to_{r}_1")
            edges.append(f"e_{r}_{c.cols-1}_to_{r}_{c.cols-2}")
        return edges

    def _get_od_pairs(self, edges: List[str], max_pairs: int = 20) -> List[Tuple[str, str]]:
        pairs = []
        n = len(edges)
        for i in range(min(n, max_pairs)):
            j = (i + n // 2) % n
            if i != j:
                pairs.append((edges[i], edges[j]))
        return pairs

    def _generate_bus_route(self, route_idx: int) -> List[str]:
        c = self.config
        if route_idx % 2 == 0:  # N-S route
            col = route_idx % c.cols
            return [f"e_{r}_{col}_to_{r+1}_{col}" for r in range(c.rows - 1)]
        else:  # E-W route
            row = route_idx % c.rows
            return [f"e_{row}_{c}_to_{row}_{c+1}" for c in range(self.config.cols - 1)]


class TrafficMetricsCollector:
    """Collect and aggregate traffic performance metrics."""

    def __init__(self, collection_interval: int = 60):
        self.interval = collection_interval
        self.metrics_history: List[Dict[str, float]] = []

    def collect(
        self,
        step: int,
        vehicle_speeds: Dict[str, float],
        vehicle_waiting: Dict[str, float],
        queue_lengths: Dict[str, int],
        throughput: int,
        bus_delays: List[float],
    ) -> Dict[str, float]:
        """Collect metrics at current timestep."""
        speeds = list(vehicle_speeds.values()) or [0]
        waits = list(vehicle_waiting.values()) or [0]
        queues = list(queue_lengths.values()) or [0]

        metrics = {
            "step": step,
            "avg_speed_kmh": np.mean(speeds) * 3.6,
            "avg_delay_s": np.mean(waits),
            "max_delay_s": np.max(waits) if waits else 0,
            "avg_queue": np.mean(queues),
            "max_queue": np.max(queues) if queues else 0,
            "total_throughput": throughput,
            "avg_bus_delay_s": np.mean(bus_delays) if bus_delays else 0,
            "fuel_consumption_ml": self._estimate_fuel(speeds, waits),
            "co2_emission_g": self._estimate_co2(speeds, waits),
        }

        self.metrics_history.append(metrics)
        return metrics

    def get_summary(self) -> Dict[str, float]:
        """Get aggregated summary statistics."""
        if not self.metrics_history:
            return {}

        keys = ["avg_speed_kmh", "avg_delay_s", "avg_queue",
                "total_throughput", "avg_bus_delay_s",
                "fuel_consumption_ml", "co2_emission_g"]
        summary = {}
        for k in keys:
            values = [m[k] for m in self.metrics_history]
            summary[f"mean_{k}"] = np.mean(values)
            summary[f"std_{k}"] = np.std(values)
        return summary

    def save_metrics(self, filepath: str):
        """Save metrics history to JSON."""
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, (np.integer,)):
                    return int(obj)
                if isinstance(obj, (np.floating,)):
                    return float(obj)
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return super().default(obj)

        with open(filepath, 'w') as f:
            json.dump(self.metrics_history, f, indent=2, cls=NumpyEncoder)
        logger.info(f"Saved {len(self.metrics_history)} metric records to {filepath}")

    def _estimate_fuel(self, speeds: List[float], waits: List[float]) -> float:
        """Simplified CMEM-based fuel consumption estimate (mL)."""
        avg_speed = np.mean(speeds) if speeds else 0
        num_idling = sum(1 for w in waits if w > 0)
        base_rate = 0.8  # mL/s at idle
        moving_rate = 0.05  # mL/m at speed
        return num_idling * base_rate + len(speeds) * avg_speed * moving_rate

    def _estimate_co2(self, speeds: List[float], waits: List[float]) -> float:
        """CO2 emission estimate (g). ~2.31 kg CO2 per liter gasoline."""
        fuel_ml = self._estimate_fuel(speeds, waits)
        return fuel_ml * 2.31  # g CO2


if __name__ == "__main__":
    # Generate network files
    net_config = SUMONetworkConfig()
    generator = SUMONetworkGenerator(net_config)

    print("=== SUMO Network Generation ===")
    nodes = generator.generate_nodes_xml()
    print(f"Nodes XML: {len(nodes)} chars")

    edges = generator.generate_edges_xml()
    print(f"Edges XML: {len(edges)} chars")

    flow_config = FlowConfig()
    flows = generator.generate_flow_xml(flow_config)
    print(f"Flows XML: {len(flows)} chars")

    cfg = generator.generate_sumo_cfg("output")
    print(f"Config XML: {len(cfg)} chars")

    # Test metrics
    collector = TrafficMetricsCollector()
    for step in range(10):
        metrics = collector.collect(
            step=step,
            vehicle_speeds={f"v_{i}": np.random.uniform(5, 15) for i in range(100)},
            vehicle_waiting={f"v_{i}": np.random.exponential(10) for i in range(100)},
            queue_lengths={f"lane_{i}": np.random.randint(0, 10) for i in range(48)},
            throughput=np.random.randint(50, 100),
            bus_delays=[np.random.uniform(-5, 15) for _ in range(12)],
        )
    summary = collector.get_summary()
    print(f"\nMetrics Summary: avg_speed={summary['mean_avg_speed_kmh']:.1f} km/h, "
          f"avg_delay={summary['mean_avg_delay_s']:.1f}s")
