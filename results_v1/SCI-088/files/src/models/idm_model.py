"""
Intelligent Driver Model (IDM) with MOBIL Lane-Change Model
============================================================
Car-following and lane-changing behavior models for SUMO integration.

References:
- Treiber, M., Hennecke, A., & Helbing, D. (2000). PRE, 62(2), 1805.
- Kesting, A., Treiber, M., & Helbing, D. (2007). TRR, 1999(1), 86-94.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, Tuple
import logging

logger = logging.getLogger(__name__)


@dataclass
class IDMParams:
    """Parameters for the Intelligent Driver Model."""
    v0: float = 13.89   # desired speed (m/s)
    T: float = 1.5      # safe time headway (s)
    a: float = 1.4      # max acceleration (m/s²)
    b: float = 2.0      # comfortable deceleration (m/s²)
    delta: int = 4       # acceleration exponent
    s0: float = 2.0      # minimum gap (m)
    s1: float = 0.0      # jam distance parameter (m)
    length: float = 4.5  # vehicle length (m)


@dataclass
class MOBILParams:
    """Parameters for the MOBIL lane-change model."""
    politeness: float = 0.3
    threshold: float = 0.2     # acceleration gain threshold (m/s²)
    max_safe_decel: float = 4.0  # safety criterion (m/s²)


@dataclass
class VehicleState:
    """State of a single vehicle."""
    id: str
    x: float = 0.0       # position (m)
    v: float = 0.0       # speed (m/s)
    lane: int = 0
    mode: str = "car"    # car, bus, bicycle


class IDMModel:
    """Intelligent Driver Model for longitudinal dynamics."""

    def __init__(self, params: IDMParams):
        self.params = params

    def desired_gap(self, v: float, delta_v: float) -> float:
        """Calculate desired minimum gap s*(v, Δv).

        s*(v, Δv) = s0 + s1*sqrt(v/v0) + v*T + v*Δv / (2*sqrt(a*b))
        """
        p = self.params
        interaction_term = (v * delta_v) / (2.0 * np.sqrt(p.a * p.b))
        s_star = p.s0 + p.s1 * np.sqrt(v / p.v0) + v * p.T + max(0, interaction_term)
        return max(s_star, p.s0)

    def acceleration(self, v: float, s: float, delta_v: float) -> float:
        """Calculate IDM acceleration.

        dv/dt = a * [1 - (v/v0)^δ - (s*(v,Δv)/s)²]

        Args:
            v: current speed (m/s)
            s: gap to leader (m)
            delta_v: speed difference v - v_leader (m/s)

        Returns:
            acceleration (m/s²)
        """
        p = self.params
        v_ratio = (v / p.v0) ** p.delta if p.v0 > 0 else 0
        s_star = self.desired_gap(v, delta_v)
        s = max(s, 0.1)  # prevent division by zero
        gap_ratio = (s_star / s) ** 2
        acc = p.a * (1.0 - v_ratio - gap_ratio)
        return np.clip(acc, -p.b * 2, p.a)

    def free_flow_acceleration(self, v: float) -> float:
        """Acceleration on empty road."""
        p = self.params
        return p.a * (1.0 - (v / p.v0) ** p.delta)


class MOBILModel:
    """MOBIL lane-change decision model."""

    def __init__(self, params: MOBILParams, idm: IDMModel):
        self.params = params
        self.idm = idm

    def evaluate_lane_change(
        self,
        ego: VehicleState,
        current_leader: Optional[VehicleState],
        current_follower: Optional[VehicleState],
        target_leader: Optional[VehicleState],
        target_follower: Optional[VehicleState],
    ) -> Tuple[bool, float]:
        """Evaluate whether a lane change is beneficial and safe.

        Returns:
            (should_change, incentive_value)
        """
        p = self.params

        # Current lane accelerations
        acc_ego_curr = self._get_acc(ego, current_leader)
        acc_follower_curr = self._get_acc(current_follower, ego) if current_follower else 0

        # Target lane accelerations (after hypothetical lane change)
        acc_ego_target = self._get_acc(ego, target_leader)
        acc_follower_target = self._get_acc(target_follower, ego) if target_follower else 0

        # Safety criterion: new follower must not brake too hard
        if target_follower and acc_follower_target < -p.max_safe_decel:
            return False, 0.0

        # MOBIL incentive criterion
        # Δa_ego + p * (Δa_follower_new + Δa_follower_old) > threshold
        delta_ego = acc_ego_target - acc_ego_curr

        # follower on old lane benefits when ego leaves
        acc_follower_old_after = self._get_acc(
            current_follower, current_leader
        ) if current_follower else 0
        delta_follower_old = acc_follower_old_after - acc_follower_curr

        delta_follower_new = acc_follower_target - (
            self._get_acc(target_follower, target_leader) if target_follower else 0
        )

        incentive = delta_ego + p.politeness * (delta_follower_new + delta_follower_old)
        should_change = incentive > p.threshold

        return should_change, incentive

    def _get_acc(self, follower: Optional[VehicleState],
                 leader: Optional[VehicleState]) -> float:
        if follower is None:
            return 0.0
        if leader is None:
            return self.idm.free_flow_acceleration(follower.v)
        gap = leader.x - follower.x - self.idm.params.length
        delta_v = follower.v - leader.v
        return self.idm.acceleration(follower.v, gap, delta_v)


class MultiModalIDM:
    """Vehicle-type-specific IDM parameterization."""

    DEFAULT_PARAMS = {
        "car": IDMParams(v0=13.89, T=1.5, a=1.4, b=2.0, s0=2.0, length=4.5),
        "bus": IDMParams(v0=11.11, T=2.0, a=1.0, b=1.5, s0=3.0, length=12.0),
        "bicycle": IDMParams(v0=4.17, T=1.0, a=1.2, b=2.5, s0=1.0, length=1.8),
    }

    def __init__(self, custom_params: Optional[Dict[str, IDMParams]] = None):
        self.models: Dict[str, IDMModel] = {}
        params = custom_params or self.DEFAULT_PARAMS
        for mode, p in params.items():
            self.models[mode] = IDMModel(p)
            logger.info(f"Initialized IDM for mode={mode}: v0={p.v0:.2f} m/s")

    def get_model(self, mode: str) -> IDMModel:
        if mode not in self.models:
            logger.warning(f"Unknown mode '{mode}', falling back to 'car'")
            return self.models["car"]
        return self.models[mode]

    def calibrate_from_trajectory(
        self, mode: str, trajectories: np.ndarray
    ) -> IDMParams:
        """Calibrate IDM parameters from observed trajectory data.

        Uses nonlinear least squares (Levenberg-Marquardt) to fit IDM
        parameters to observed speed-gap-acceleration triples.

        Args:
            mode: vehicle type
            trajectories: (N, 4) array of [time, gap, speed, acceleration]

        Returns:
            Calibrated IDMParams
        """
        from scipy.optimize import least_squares

        current = self.DEFAULT_PARAMS.get(mode, self.DEFAULT_PARAMS["car"])

        def residuals(x):
            v0, T, a_max, b, s0 = x
            params = IDMParams(v0=v0, T=T, a=a_max, b=b, s0=s0)
            model = IDMModel(params)
            predicted = np.array([
                model.acceleration(row[2], row[1], 0.0) for row in trajectories
            ])
            return predicted - trajectories[:, 3]

        x0 = [current.v0, current.T, current.a, current.b, current.s0]
        bounds = ([5, 0.5, 0.5, 0.5, 0.5], [30, 5.0, 3.0, 5.0, 5.0])
        result = least_squares(residuals, x0, bounds=bounds, method='trf')

        calibrated = IDMParams(
            v0=result.x[0], T=result.x[1], a=result.x[2],
            b=result.x[3], s0=result.x[4], length=current.length
        )
        self.models[mode] = IDMModel(calibrated)
        logger.info(f"Calibrated {mode}: v0={calibrated.v0:.2f}, T={calibrated.T:.2f}")
        return calibrated


# --- SUMO Integration Helper ---

def generate_sumo_vtype_xml(params_dict: Dict[str, IDMParams]) -> str:
    """Generate SUMO vType XML from IDM parameters."""
    lines = ['<additional>']
    for mode, p in params_dict.items():
        lines.append(f'  <vType id="{mode}" '
                     f'length="{p.length}" '
                     f'maxSpeed="{p.v0}" '
                     f'accel="{p.a}" '
                     f'decel="{p.b}" '
                     f'sigma="0.5" '
                     f'minGap="{p.s0}" '
                     f'tau="{p.T}" '
                     f'carFollowModel="IDM" '
                     f'lcModel="LC2013">')
        lines.append(f'    <param key="delta" value="{p.delta}"/>')
        lines.append(f'  </vType>')
    lines.append('</additional>')
    return '\n'.join(lines)


if __name__ == "__main__":
    # Quick validation
    mm = MultiModalIDM()
    car_model = mm.get_model("car")

    # Test: following scenario
    v, s, dv = 10.0, 20.0, 2.0
    acc = car_model.acceleration(v, s, dv)
    print(f"Car IDM acceleration: v={v}, s={s}, Δv={dv} → a={acc:.3f} m/s²")

    # Test: free-flow
    acc_free = car_model.free_flow_acceleration(5.0)
    print(f"Free-flow acceleration at 5 m/s: {acc_free:.3f} m/s²")

    # Generate SUMO XML
    xml = generate_sumo_vtype_xml(MultiModalIDM.DEFAULT_PARAMS)
    print("\nSUMO vType XML:")
    print(xml)
