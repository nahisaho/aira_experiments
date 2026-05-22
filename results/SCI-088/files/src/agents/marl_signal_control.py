"""
Multi-Agent Reinforcement Learning for Traffic Signal Control
=============================================================
MAPPO (Multi-Agent PPO) based signal controller using RLlib.

Architecture:
  - Each intersection is an independent agent with local observations
  - Observations include queue lengths, waiting times, current phase,
    neighbor states (within 2-hop radius)
  - Centralized training with decentralized execution (CTDE)
  - Transit Signal Priority (TSP) integrated as reward shaping

References:
- Yu, C., et al. (2022). MAPPO. NeurIPS.
- Wei, H., et al. (2019). PressLight. KDD.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
from enum import IntEnum
import logging

logger = logging.getLogger(__name__)


# =============================================================================
# Signal Phase Definitions
# =============================================================================

class Phase(IntEnum):
    NS_GREEN = 0       # North-South through + right
    NS_LEFT = 1        # North-South left turn
    EW_GREEN = 2       # East-West through + right
    EW_LEFT = 3        # East-West left turn


@dataclass
class SignalTiming:
    """Signal timing for one intersection."""
    cycle_length: int = 90
    phases: List[int] = field(default_factory=lambda: [30, 15, 30, 15])
    yellow: int = 3
    all_red: int = 2
    min_green: int = 8
    max_green: int = 60


# =============================================================================
# Observation & Action Spaces
# =============================================================================

@dataclass
class IntersectionObservation:
    """Observation vector for a single intersection agent."""
    queue_lengths: np.ndarray        # (num_approaches,) e.g., (4,) for NSEW
    waiting_times: np.ndarray        # (num_approaches,)
    current_phase: int
    phase_elapsed: float             # seconds in current phase
    neighbor_queues: np.ndarray      # aggregated neighbor info
    time_of_day: float               # normalized [0, 1]
    demand_estimate: float           # estimated demand level
    bus_approaching: np.ndarray      # (num_approaches,) binary flags

    def to_vector(self) -> np.ndarray:
        """Flatten to fixed-size observation vector."""
        return np.concatenate([
            self.queue_lengths,
            self.waiting_times,
            [self.current_phase / 3.0],
            [self.phase_elapsed / 60.0],
            self.neighbor_queues,
            [self.time_of_day],
            [self.demand_estimate],
            self.bus_approaching,
        ]).astype(np.float32)


# =============================================================================
# Reward Function
# =============================================================================

class RewardCalculator:
    """Composite reward function for traffic signal control."""

    def __init__(self, weights: Optional[Dict[str, float]] = None):
        self.weights = weights or {
            "waiting_time": -0.4,
            "queue_length": -0.3,
            "throughput": 0.2,
            "transit_priority": 0.1,
        }

    def calculate(
        self,
        waiting_time: float,
        queue_length: float,
        throughput: int,
        bus_delay: float,
        prev_waiting: float,
        prev_queue: float,
    ) -> float:
        """Calculate composite reward.

        Uses change-based rewards to reduce variance:
        r = w1 * Δwaiting + w2 * Δqueue + w3 * throughput + w4 * bus_bonus
        """
        w = self.weights

        delta_wait = prev_waiting - waiting_time  # positive = improvement
        delta_queue = prev_queue - queue_length

        # Bus priority bonus: negative delay means bus was prioritized
        bus_bonus = max(0, -bus_delay) if bus_delay != 0 else 0

        reward = (
            w["waiting_time"] * delta_wait
            + w["queue_length"] * delta_queue
            + w["throughput"] * throughput
            + w["transit_priority"] * bus_bonus
        )
        return float(reward)


# =============================================================================
# MAPPO Agent Configuration for RLlib
# =============================================================================

def build_mappo_config(
    num_agents: int = 48,
    obs_dim: int = 52,
    num_actions: int = 4,
    config_overrides: Optional[Dict] = None,
) -> Dict[str, Any]:
    """Build RLlib multi-agent PPO configuration.

    Returns a config dict compatible with ray.rllib.algorithms.ppo.PPOConfig.
    """
    agent_ids = [f"intersection_{i}" for i in range(num_agents)]

    def policy_mapping_fn(agent_id, episode=None, worker=None, **kwargs):
        """All agents share the same policy (parameter sharing)."""
        return "shared_policy"

    config = {
        "framework": "torch",
        "env": "TrafficSignalEnv",

        # Multi-agent setup
        "multiagent": {
            "policies": {
                "shared_policy": (
                    None,                        # default policy class
                    None,                        # obs space (set by env)
                    None,                        # action space (set by env)
                    {
                        "model": {
                            "fcnet_hiddens": [256, 256, 128],
                            "fcnet_activation": "relu",
                            "vf_share_layers": False,
                        },
                    },
                ),
            },
            "policy_mapping_fn": policy_mapping_fn,
            "policies_to_train": ["shared_policy"],
        },

        # PPO hyperparameters
        "lr": 3e-4,
        "gamma": 0.99,
        "lambda": 0.95,
        "clip_param": 0.2,
        "entropy_coeff": 0.01,
        "vf_loss_coeff": 0.5,
        "num_sgd_iter": 10,
        "sgd_minibatch_size": 256,
        "train_batch_size": 4096,

        # Workers
        "num_workers": 8,
        "num_envs_per_worker": 2,
        "rollout_fragment_length": 200,

        # Evaluation
        "evaluation_interval": 10,
        "evaluation_num_episodes": 5,
        "evaluation_config": {
            "explore": False,
        },

        # Stopping criteria
        "stop": {
            "episodes_total": 5000,
        },
    }

    if config_overrides:
        config.update(config_overrides)

    return config


# =============================================================================
# Traffic Signal Environment (SUMO/Flow Interface)
# =============================================================================

class TrafficSignalController:
    """Per-intersection signal controller with RL action interface.

    Integrates with SUMO via TraCI for:
    - Reading detector data (queue, waiting time)
    - Setting signal phases
    - Transit Signal Priority (TSP) logic
    """

    def __init__(
        self,
        intersection_id: str,
        timing: SignalTiming,
        neighbor_ids: List[str],
        tsp_enabled: bool = True,
    ):
        self.id = intersection_id
        self.timing = timing
        self.neighbor_ids = neighbor_ids
        self.tsp_enabled = tsp_enabled

        self.current_phase = Phase.NS_GREEN
        self.phase_timer = 0
        self.yellow_active = False
        self.reward_calc = RewardCalculator()

        # State tracking
        self.prev_waiting = 0.0
        self.prev_queue = 0.0
        self._total_throughput = 0

    def get_observation(self, traci_conn=None, neighbors=None) -> np.ndarray:
        """Collect observation from SUMO via TraCI or from cached state."""
        if traci_conn is not None:
            return self._observe_from_traci(traci_conn, neighbors)
        return self._synthetic_observation()

    def apply_action(self, action: int, traci_conn=None) -> float:
        """Apply RL action (phase selection) and return reward.

        Action semantics:
          0: Keep current phase
          1-3: Switch to phase 1, 2, 3

        With TSP override: if bus is approaching and current phase
        serves the bus, extend green regardless of RL action.
        """
        target_phase = Phase(action)

        # TSP override logic
        if self.tsp_enabled and self._bus_approaching(target_phase):
            target_phase = self.current_phase  # extend current green
            logger.debug(f"{self.id}: TSP override - extending green for bus")

        if target_phase != self.current_phase:
            if self.phase_timer >= self.timing.min_green:
                self._initiate_phase_change(target_phase)
        elif self.phase_timer >= self.timing.max_green:
            next_phase = Phase((self.current_phase + 1) % len(Phase))
            self._initiate_phase_change(next_phase)

        self.phase_timer += 1

        # Calculate reward
        obs = self._synthetic_observation()
        queue = float(np.mean(obs[:4]))
        waiting = float(np.mean(obs[4:8]))
        reward = self.reward_calc.calculate(
            waiting, queue, self._total_throughput, 0.0,
            self.prev_waiting, self.prev_queue
        )
        self.prev_waiting = waiting
        self.prev_queue = queue
        return reward

    def _initiate_phase_change(self, target: Phase):
        """Yellow → All-Red → New Green sequence."""
        self.current_phase = target
        self.phase_timer = 0

    def _bus_approaching(self, target_phase: Phase) -> bool:
        """Check if bus is on approach served by current phase."""
        # Placeholder: in real system, query TraCI for bus positions
        return False

    def _observe_from_traci(self, traci_conn, neighbors) -> np.ndarray:
        """Collect real observation from SUMO TraCI connection."""
        # In production, this queries traci_conn.lane.getLastStepHaltingNumber etc.
        return self._synthetic_observation()

    def _synthetic_observation(self) -> np.ndarray:
        """Generate synthetic observation for testing."""
        rng = np.random.default_rng()
        obs = IntersectionObservation(
            queue_lengths=rng.integers(0, 15, size=4).astype(np.float32),
            waiting_times=rng.uniform(0, 60, size=4).astype(np.float32),
            current_phase=int(self.current_phase),
            phase_elapsed=float(self.phase_timer),
            neighbor_queues=rng.integers(0, 15, size=16).astype(np.float32),
            time_of_day=0.5,
            demand_estimate=1.0,
            bus_approaching=rng.integers(0, 2, size=4).astype(np.float32),
        )
        return obs.to_vector()


# =============================================================================
# Multi-Agent Coordination Graph
# =============================================================================

class IntersectionNetwork:
    """Graph of intersection agents for MARL coordination."""

    def __init__(self, adjacency: Dict[str, List[str]]):
        """
        Args:
            adjacency: {intersection_id: [neighbor_ids]}
        """
        self.adjacency = adjacency
        self.agents: Dict[str, TrafficSignalController] = {}

    def initialize_agents(self, tsp_enabled: bool = True):
        for iid, neighbors in self.adjacency.items():
            self.agents[iid] = TrafficSignalController(
                intersection_id=iid,
                timing=SignalTiming(),
                neighbor_ids=neighbors,
                tsp_enabled=tsp_enabled,
            )
        logger.info(f"Initialized {len(self.agents)} intersection agents")

    def get_all_observations(self) -> Dict[str, np.ndarray]:
        return {iid: agent.get_observation() for iid, agent in self.agents.items()}

    def apply_all_actions(self, actions: Dict[str, int]) -> Dict[str, float]:
        rewards = {}
        for iid, action in actions.items():
            if iid in self.agents:
                rewards[iid] = self.agents[iid].apply_action(action)
        return rewards

    def get_neighbor_info(self, agent_id: str) -> np.ndarray:
        """Aggregate neighbor queue information for agent's observation."""
        neighbors = self.adjacency.get(agent_id, [])
        info = []
        for nid in neighbors[:4]:  # max 4 neighbors
            if nid in self.agents:
                obs = self.agents[nid].get_observation()
                info.append(obs[:4])  # queue lengths only
        while len(info) < 4:
            info.append(np.zeros(4))
        return np.concatenate(info)


# =============================================================================
# Tokyo Downtown Grid Generator
# =============================================================================

def create_tokyo_grid(rows: int = 8, cols: int = 6) -> IntersectionNetwork:
    """Create a grid network approximating Tokyo downtown intersections.

    48 intersections (8 x 6) covering ~3km x 3km area.
    Block size: ~375m x 500m (typical Tokyo downtown block).
    """
    adjacency = {}
    for r in range(rows):
        for c in range(cols):
            iid = f"intersection_{r}_{c}"
            neighbors = []
            if r > 0:
                neighbors.append(f"intersection_{r-1}_{c}")
            if r < rows - 1:
                neighbors.append(f"intersection_{r+1}_{c}")
            if c > 0:
                neighbors.append(f"intersection_{r}_{c-1}")
            if c < cols - 1:
                neighbors.append(f"intersection_{r}_{c+1}")
            adjacency[iid] = neighbors

    network = IntersectionNetwork(adjacency)
    network.initialize_agents(tsp_enabled=True)
    logger.info(f"Created Tokyo grid: {rows}x{cols} = {rows*cols} intersections")
    return network


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Create network
    network = create_tokyo_grid()
    print(f"Network: {len(network.agents)} agents")

    # Test step
    obs = network.get_all_observations()
    print(f"Observation dim per agent: {len(list(obs.values())[0])}")

    actions = {iid: np.random.randint(0, 4) for iid in network.agents}
    rewards = network.apply_all_actions(actions)
    print(f"Mean reward: {np.mean(list(rewards.values())):.4f}")

    # Print RLlib config
    config = build_mappo_config()
    print(f"\nRLlib config keys: {list(config.keys())}")
