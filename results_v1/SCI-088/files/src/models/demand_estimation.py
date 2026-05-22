"""
Real-Time Traffic Demand Estimation
====================================
Kalman Filter + probe data fusion for OD demand estimation.

Methods:
  1. Extended Kalman Filter for state-space demand estimation
  2. Probe vehicle trajectory → link travel times → path flow estimation
  3. Multi-source fusion (probe + loop detector + Bluetooth)

References:
- Antoniou, C., et al. (2007). Dynamic traffic demand estimation. TRC.
- Cascetta, E. (2009). Transportation Systems Analysis. Springer.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


@dataclass
class ODZone:
    """Origin-Destination zone definition."""
    id: int
    centroid: Tuple[float, float]  # (lat, lon)
    area_km2: float
    base_generation: float         # vehicles/hour


@dataclass
class ProbeData:
    """Processed probe vehicle data."""
    timestamp: float
    vehicle_id: str
    link_id: str
    speed: float        # m/s
    travel_time: float  # seconds on link
    position: Tuple[float, float]


class KalmanDemandEstimator:
    """Extended Kalman Filter for real-time OD demand estimation.

    State vector: OD flows [q_11, q_12, ..., q_nn] (n² elements for n zones)
    Measurement: link counts / travel times from detectors + probe data
    """

    def __init__(
        self,
        num_zones: int = 25,
        probe_penetration: float = 0.15,
        update_interval: float = 300,
    ):
        self.num_zones = num_zones
        self.num_od = num_zones * num_zones
        self.probe_penetration = probe_penetration
        self.update_interval = update_interval

        # State: OD demand vector
        self.x = np.ones(self.num_od) * 10.0  # initial demand guess
        # State covariance
        self.P = np.eye(self.num_od) * 100.0
        # Process noise
        self.Q = np.eye(self.num_od) * 5.0
        # Measurement noise (set per update based on data quality)
        self.R_base = 10.0

        # Assignment matrix: links x OD pairs
        # H[l, od] = fraction of OD pair 'od' using link 'l'
        self.H: Optional[np.ndarray] = None
        self.num_links = 0

        self._history: List[np.ndarray] = []

    def initialize_assignment_matrix(
        self, num_links: int, assignment_fractions: np.ndarray
    ):
        """Set the assignment matrix from route choice model.

        Args:
            num_links: number of monitored links
            assignment_fractions: (num_links, num_od) matrix
        """
        self.num_links = num_links
        self.H = assignment_fractions
        self.R = np.eye(num_links) * self.R_base
        logger.info(f"Assignment matrix: {num_links} links × {self.num_od} OD pairs")

    def predict(self, historical_factor: float = 1.0):
        """Prediction step: propagate state with transition model.

        x_{k|k-1} = F * x_{k-1} + process noise
        We use a simple persistence model: F = I (demand stays same)
        with optional historical correction factor.
        """
        F = np.eye(self.num_od) * historical_factor
        self.x = F @ self.x
        self.P = F @ self.P @ F.T + self.Q

    def update(
        self,
        link_counts: np.ndarray,
        probe_speeds: Optional[Dict[str, float]] = None,
    ):
        """Update step: correct state with measurements.

        Args:
            link_counts: observed vehicle counts on monitored links
            probe_speeds: {link_id: average_speed} from probe data
        """
        if self.H is None:
            logger.warning("Assignment matrix not set, skipping update")
            return

        # Predicted measurements
        z_pred = self.H @ self.x

        # Innovation
        y = link_counts - z_pred

        # Innovation covariance
        S = self.H @ self.P @ self.H.T + self.R

        # Kalman gain
        K = self.P @ self.H.T @ np.linalg.inv(S)

        # State update
        self.x = self.x + K @ y
        self.x = np.maximum(self.x, 0)  # demand must be non-negative

        # Covariance update (Joseph form for numerical stability)
        I = np.eye(self.num_od)
        IKH = I - K @ self.H
        self.P = IKH @ self.P @ IKH.T + K @ self.R @ K.T

        self._history.append(self.x.copy())
        logger.debug(f"Demand update: total={self.x.sum():.0f} veh/h")

    def fuse_probe_data(
        self, probe_records: List[ProbeData]
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """Process probe data into link-level measurements.

        1. Scale up probe counts by 1/penetration_rate
        2. Calculate link travel times
        3. Estimate link flows

        Returns:
            (estimated_link_counts, link_speeds)
        """
        link_counts: Dict[str, int] = {}
        link_speeds: Dict[str, List[float]] = {}

        for record in probe_records:
            link_counts[record.link_id] = link_counts.get(record.link_id, 0) + 1
            link_speeds.setdefault(record.link_id, []).append(record.speed)

        # Scale up by penetration rate
        scaled_counts = {
            lid: count / self.probe_penetration
            for lid, count in link_counts.items()
        }

        avg_speeds = {
            lid: np.mean(speeds) for lid, speeds in link_speeds.items()
        }

        return scaled_counts, avg_speeds

    def get_od_matrix(self) -> np.ndarray:
        """Return current OD demand estimate as (zones x zones) matrix."""
        return self.x.reshape(self.num_zones, self.num_zones)

    def get_zone_generation(self) -> np.ndarray:
        """Return total trip generation per zone."""
        od = self.get_od_matrix()
        return od.sum(axis=1)

    def get_zone_attraction(self) -> np.ndarray:
        """Return total trip attraction per zone."""
        od = self.get_od_matrix()
        return od.sum(axis=0)


class MultiSourceFusion:
    """Fusion of multiple traffic data sources.

    Weighted combination of:
    - Loop detector counts (high accuracy, fixed locations)
    - Probe vehicle data (wide coverage, sampled)
    - Bluetooth/WiFi sensors (travel time estimation)
    """

    def __init__(self, source_weights: Optional[Dict[str, float]] = None):
        self.weights = source_weights or {
            "loop_detector": 0.5,
            "probe": 0.3,
            "bluetooth": 0.2,
        }

    def fuse(
        self,
        loop_data: Optional[np.ndarray] = None,
        probe_data: Optional[np.ndarray] = None,
        bt_data: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Weighted fusion of available data sources.

        Missing sources are excluded and weights renormalized.
        """
        sources = {}
        if loop_data is not None:
            sources["loop_detector"] = loop_data
        if probe_data is not None:
            sources["probe"] = probe_data
        if bt_data is not None:
            sources["bluetooth"] = bt_data

        if not sources:
            raise ValueError("No data sources available")

        total_weight = sum(self.weights[k] for k in sources)
        result = sum(
            self.weights[k] / total_weight * v for k, v in sources.items()
        )
        return result


class HistoricalDemandProfile:
    """Time-of-day demand profile from historical data."""

    def __init__(self):
        # Hourly multipliers for typical Tokyo weekday
        self.hourly_factors = np.array([
            0.15, 0.10, 0.08, 0.08, 0.10, 0.30,  # 0-5h
            0.65, 1.00, 0.95, 0.70, 0.60, 0.65,  # 6-11h (morning peak)
            0.75, 0.70, 0.65, 0.70, 0.80, 1.00,  # 12-17h (evening peak)
            0.90, 0.70, 0.55, 0.40, 0.30, 0.20,  # 18-23h
        ])

    def get_factor(self, hour: float) -> float:
        """Get demand scaling factor for given hour of day."""
        idx = int(hour) % 24
        frac = hour - int(hour)
        next_idx = (idx + 1) % 24
        return self.hourly_factors[idx] * (1 - frac) + self.hourly_factors[next_idx] * frac


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Initialize estimator
    estimator = KalmanDemandEstimator(num_zones=5, probe_penetration=0.15)

    # Create simple assignment matrix
    num_links = 20
    H = np.random.dirichlet(np.ones(estimator.num_od), size=num_links)
    estimator.initialize_assignment_matrix(num_links, H)

    # Simulate estimation cycle
    for step in range(10):
        estimator.predict()
        link_counts = H @ (np.ones(estimator.num_od) * 15 + np.random.randn(estimator.num_od) * 3)
        link_counts = np.maximum(link_counts, 0)
        estimator.update(link_counts)

    od = estimator.get_od_matrix()
    print(f"OD Matrix shape: {od.shape}")
    print(f"Total demand: {od.sum():.0f} veh/h")
    print(f"Zone generation: {estimator.get_zone_generation()}")

    # Historical profile
    profile = HistoricalDemandProfile()
    print(f"\nDemand factor at 8:30 AM: {profile.get_factor(8.5):.2f}")
    print(f"Demand factor at 3:00 PM: {profile.get_factor(15.0):.2f}")
