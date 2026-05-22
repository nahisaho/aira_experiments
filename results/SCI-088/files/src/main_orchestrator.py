"""
Main Simulation Orchestrator
=============================
Integrates all components: SUMO network, IDM models, MARL signal control,
demand estimation, and dynamic rerouting into a unified pipeline.
"""

import numpy as np
import json
import os
import logging
from typing import Dict, Any

# Internal modules
from src.models.idm_model import MultiModalIDM, generate_sumo_vtype_xml
from src.agents.marl_signal_control import (
    create_tokyo_grid, build_mappo_config, IntersectionNetwork
)
from src.models.demand_estimation import (
    KalmanDemandEstimator, HistoricalDemandProfile, MultiSourceFusion
)
from src.models.dynamic_routing import (
    DynamicRouter, IncidentDetector, Incident, IncidentType
)
from src.network.sumo_environment import (
    SUMONetworkGenerator, SUMONetworkConfig, FlowConfig, TrafficMetricsCollector
)

logger = logging.getLogger(__name__)


class SimulationOrchestrator:
    """Main orchestrator for the integrated traffic simulation system.

    Pipeline:
    1. Initialize SUMO network + vehicle types
    2. Initialize MARL signal controllers
    3. Run simulation loop:
       a. Update demand estimation from probe data
       b. Compute RL observations
       c. Select signal actions (MARL)
       d. Step SUMO simulation
       e. Detect incidents → dynamic rerouting
       f. Collect metrics
    4. Output results
    """

    def __init__(self, config_path: str = "configs/simulation_config.yaml"):
        self.config = self._load_config(config_path)
        self.step = 0
        self.total_steps = int(
            self.config.get("simulation", {}).get("time", {}).get("end", 7200)
            / self.config.get("simulation", {}).get("time", {}).get("step_length", 0.1)
        )

        # Components
        self.idm = MultiModalIDM()
        self.signal_network: IntersectionNetwork = None
        self.demand_estimator: KalmanDemandEstimator = None
        self.router: DynamicRouter = None
        self.incident_detector: IncidentDetector = None
        self.metrics: TrafficMetricsCollector = None
        self.demand_profile = HistoricalDemandProfile()
        self.fusion = MultiSourceFusion()

    def _load_config(self, path: str) -> Dict[str, Any]:
        """Load YAML config (with fallback to defaults)."""
        try:
            import yaml
            with open(path, 'r') as f:
                return yaml.safe_load(f)
        except (ImportError, FileNotFoundError):
            logger.warning(f"Config not loaded from {path}, using defaults")
            return {}

    def initialize(self):
        """Initialize all simulation components."""
        logger.info("=" * 60)
        logger.info("Initializing Urban Traffic Simulation System")
        logger.info("=" * 60)

        # 1. SUMO Network
        net_config = SUMONetworkConfig(rows=8, cols=6)
        net_gen = SUMONetworkGenerator(net_config)
        os.makedirs("results/network", exist_ok=True)
        with open("results/network/tokyo_nodes.nod.xml", 'w') as f:
            f.write(net_gen.generate_nodes_xml())
        with open("results/network/tokyo_edges.edg.xml", 'w') as f:
            f.write(net_gen.generate_edges_xml())
        logger.info("SUMO network files generated")

        # 2. Vehicle type XML
        vtype_xml = generate_sumo_vtype_xml(MultiModalIDM.DEFAULT_PARAMS)
        with open("results/network/vtypes.add.xml", 'w') as f:
            f.write(vtype_xml)
        logger.info("Vehicle type definitions generated")

        # 3. MARL Signal Network
        self.signal_network = create_tokyo_grid(rows=8, cols=6)
        logger.info(f"Signal network: {len(self.signal_network.agents)} agents")

        # 4. Demand Estimation
        de_config = self.config.get("demand_estimation", {})
        self.demand_estimator = KalmanDemandEstimator(
            num_zones=de_config.get("od_zones", 25),
            probe_penetration=de_config.get("probe_penetration_rate", 0.15),
        )
        num_links = 100
        H = np.random.dirichlet(
            np.ones(self.demand_estimator.num_od),
            size=num_links
        )
        self.demand_estimator.initialize_assignment_matrix(num_links, H)
        logger.info("Demand estimator initialized")

        # 5. Dynamic Router
        self.router = DynamicRouter(compliance_rate=0.7)
        links = []
        for r in range(8):
            for c in range(6):
                if c < 5:
                    links.append(
                        (f"link_{r}_{c}_E", f"n_{r}_{c}", f"n_{r}_{c+1}", 500, 13.89)
                    )
                if r < 7:
                    links.append(
                        (f"link_{r}_{c}_S", f"n_{r}_{c}", f"n_{r+1}_{c}", 375, 13.89)
                    )
        self.router.build_graph(links)
        logger.info(f"Router: {len(self.router.link_states)} links")

        # 6. Incident Detector
        self.incident_detector = IncidentDetector(
            speed_threshold_ratio=0.3,
            confirmation_time=120,
        )

        # 7. Metrics Collector
        self.metrics = TrafficMetricsCollector(collection_interval=60)

        logger.info("All components initialized successfully")

    def run_synthetic_evaluation(self, num_episodes: int = 3):
        """Run synthetic evaluation without SUMO (for design validation).

        Simulates the full pipeline with synthetic data to validate
        component integration and metric collection.
        """
        logger.info("Starting synthetic evaluation")
        all_results = []

        for episode in range(num_episodes):
            logger.info(f"\n--- Episode {episode + 1}/{num_episodes} ---")
            episode_metrics = []

            for step in range(0, 7200, 60):  # 60-second intervals
                sim_time = step
                hour = 7.0 + step / 3600  # start at 7 AM

                # Demand estimation
                demand_factor = self.demand_profile.get_factor(hour)
                self.demand_estimator.predict(historical_factor=demand_factor)

                # Generate synthetic link counts
                link_counts = np.random.poisson(
                    15 * demand_factor,
                    size=self.demand_estimator.num_links if hasattr(self.demand_estimator, 'num_links') else 100
                ).astype(float)
                self.demand_estimator.update(link_counts)

                # Signal control step
                obs = self.signal_network.get_all_observations()
                actions = {
                    iid: np.random.randint(0, 4)
                    for iid in self.signal_network.agents
                }
                rewards = self.signal_network.apply_all_actions(actions)

                # Collect metrics
                num_vehicles = int(200 * demand_factor)
                metrics = self.metrics.collect(
                    step=step,
                    vehicle_speeds={
                        f"v_{i}": np.random.uniform(3, 14) * (1.0 - 0.3 * demand_factor)
                        for i in range(num_vehicles)
                    },
                    vehicle_waiting={
                        f"v_{i}": np.random.exponential(15 * demand_factor)
                        for i in range(num_vehicles)
                    },
                    queue_lengths={
                        f"lane_{i}": int(np.random.poisson(5 * demand_factor))
                        for i in range(48)
                    },
                    throughput=int(np.random.poisson(80 * demand_factor)),
                    bus_delays=[
                        np.random.normal(10 * demand_factor, 5)
                        for _ in range(12)
                    ],
                )
                episode_metrics.append(metrics)

                # Incident simulation (at step 2400 in episode 1)
                if episode == 1 and step == 2400:
                    incident = Incident(
                        id="inc_test",
                        type=IncidentType.ACCIDENT,
                        link_id="link_3_2_E",
                        start_time=sim_time,
                        estimated_duration=1200,
                        capacity_reduction=0.2,
                        confirmed=True,
                    )
                    self.router.apply_incident(incident)
                    logger.info("Incident injected at link_3_2_E")

                if episode == 1 and step == 3600:
                    self.router.clear_incident("link_3_2_E")
                    logger.info("Incident cleared")

            # Episode summary
            summary = self.metrics.get_summary()
            all_results.append({
                "episode": episode + 1,
                "summary": {k: round(v, 3) for k, v in summary.items()},
            })

        return all_results

    def save_results(self, results: list):
        """Save all results to files."""
        os.makedirs("results", exist_ok=True)

        # Metrics history
        self.metrics.save_metrics("results/metrics_history.json")

        # Summary results
        with open("results/evaluation_summary.json", 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        logger.info("Results saved to results/")

        # RLlib config
        mappo_config = build_mappo_config()
        with open("results/mappo_config.json", 'w') as f:
            # Remove non-serializable items
            serializable = {k: v for k, v in mappo_config.items()
                          if k != "multiagent"}
            serializable["multiagent"] = {
                "num_policies": 1,
                "policy_type": "shared_policy",
                "parameter_sharing": True,
            }
            json.dump(serializable, f, indent=2)
        logger.info("MAPPO config saved")


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    orchestrator = SimulationOrchestrator()
    orchestrator.initialize()
    results = orchestrator.run_synthetic_evaluation(num_episodes=3)
    orchestrator.save_results(results)

    # Print summary
    print("\n" + "=" * 60)
    print("SIMULATION COMPLETE")
    print("=" * 60)
    for r in results:
        ep = r["episode"]
        s = r["summary"]
        print(f"\nEpisode {ep}:")
        print(f"  Avg Speed:      {s.get('mean_avg_speed_kmh', 0):.1f} km/h")
        print(f"  Avg Delay:      {s.get('mean_avg_delay_s', 0):.1f} s")
        print(f"  Avg Queue:      {s.get('mean_avg_queue', 0):.1f} veh")
        print(f"  Throughput:     {s.get('mean_total_throughput', 0):.0f} veh/interval")
        print(f"  Bus Delay:      {s.get('mean_avg_bus_delay_s', 0):.1f} s")
        print(f"  CO2 Emission:   {s.get('mean_co2_emission_g', 0):.1f} g/step")


if __name__ == "__main__":
    main()
