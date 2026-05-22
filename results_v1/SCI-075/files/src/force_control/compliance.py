"""
Force Sensing and Compliance Control Module
=============================================
Impedance/admittance control for safe tissue interaction.
Includes force estimation, filtering, and adaptive compliance.
"""

import numpy as np
from typing import Optional, Tuple
from dataclasses import dataclass
from enum import Enum


class ControlMode(Enum):
    IMPEDANCE = "impedance"
    ADMITTANCE = "admittance"
    HYBRID = "hybrid"


@dataclass
class ImpedanceParams:
    """Impedance controller parameters: M*ddx + D*dx + K*x = F_ext."""
    mass: np.ndarray        # (6,) virtual inertia [kg, kg*m^2]
    damping: np.ndarray     # (6,) virtual damping [Ns/m, Nms/rad]
    stiffness: np.ndarray   # (6,) virtual stiffness [N/m, Nm/rad]

    @classmethod
    def default_suturing(cls) -> 'ImpedanceParams':
        """Default parameters tuned for suturing tasks."""
        return cls(
            mass=np.array([0.5, 0.5, 0.5, 0.05, 0.05, 0.05]),
            damping=np.array([10.0, 10.0, 15.0, 1.0, 1.0, 1.0]),
            stiffness=np.array([200.0, 200.0, 300.0, 20.0, 20.0, 20.0])
        )

    @classmethod
    def for_phase(cls, phase: str) -> 'ImpedanceParams':
        """Phase-specific impedance parameters."""
        configs = {
            'approach': cls(
                mass=np.array([0.5, 0.5, 0.5, 0.05, 0.05, 0.05]),
                damping=np.array([8.0, 8.0, 8.0, 0.8, 0.8, 0.8]),
                stiffness=np.array([300.0, 300.0, 300.0, 30.0, 30.0, 30.0])
            ),
            'insert': cls(
                mass=np.array([0.3, 0.3, 0.3, 0.03, 0.03, 0.03]),
                damping=np.array([15.0, 15.0, 20.0, 1.5, 1.5, 1.5]),
                stiffness=np.array([150.0, 150.0, 200.0, 15.0, 15.0, 15.0])
            ),
            'pull_through': cls(
                mass=np.array([0.5, 0.5, 0.5, 0.05, 0.05, 0.05]),
                damping=np.array([12.0, 12.0, 12.0, 1.2, 1.2, 1.2]),
                stiffness=np.array([250.0, 250.0, 250.0, 25.0, 25.0, 25.0])
            ),
            'knot_tying': cls(
                mass=np.array([0.4, 0.4, 0.4, 0.04, 0.04, 0.04]),
                damping=np.array([10.0, 10.0, 10.0, 1.0, 1.0, 1.0]),
                stiffness=np.array([350.0, 350.0, 350.0, 35.0, 35.0, 35.0])
            ),
        }
        return configs.get(phase, cls.default_suturing())


class ForceFilter:
    """
    Multi-stage force signal processing pipeline.
    Low-pass filter + moving average + spike removal.
    """

    def __init__(self, cutoff_freq: float = 30.0, sample_rate: float = 1000.0,
                 window_size: int = 5, spike_threshold: float = 5.0):
        self.cutoff_freq = cutoff_freq
        self.sample_rate = sample_rate
        self.window_size = window_size
        self.spike_threshold = spike_threshold

        # Low-pass filter coefficients (first-order IIR)
        rc = 1.0 / (2 * np.pi * cutoff_freq)
        dt = 1.0 / sample_rate
        self.alpha = dt / (rc + dt)

        # State
        self.prev_filtered = np.zeros(6)
        self.history: list = []
        self.initialized = False

    def filter(self, raw_force: np.ndarray) -> np.ndarray:
        """
        Apply filtering pipeline to raw force/torque reading.

        Parameters
        ----------
        raw_force : (6,) [Fx, Fy, Fz, Tx, Ty, Tz]

        Returns
        -------
        filtered : (6,) filtered force/torque
        """
        if not self.initialized:
            self.prev_filtered = raw_force.copy()
            self.initialized = True
            return raw_force.copy()

        # Spike detection and removal
        diff = np.abs(raw_force - self.prev_filtered)
        std = np.std(self.history[-min(50, len(self.history)):], axis=0) if len(self.history) > 5 else np.ones(6)
        std = np.maximum(std, 0.01)
        spike_mask = diff > self.spike_threshold * std
        corrected = np.where(spike_mask, self.prev_filtered, raw_force)

        # IIR low-pass filter
        filtered = self.alpha * corrected + (1 - self.alpha) * self.prev_filtered

        # Moving average
        self.history.append(filtered.copy())
        if len(self.history) > self.window_size:
            self.history = self.history[-self.window_size:]
        smoothed = np.mean(self.history, axis=0)

        self.prev_filtered = smoothed
        return smoothed


class ForceEstimator:
    """
    Model-based force estimation for dVRK (no direct F/T sensor).
    Uses motor current → joint torque → external force mapping.
    """

    def __init__(self, n_joints: int = 7):
        self.n_joints = n_joints
        self.gravity_comp = np.zeros(n_joints)
        self.friction_comp = np.zeros(n_joints)
        self._calibrated = False

    def calibrate(self, joint_positions: np.ndarray,
                  motor_currents: np.ndarray,
                  n_samples: int = 100):
        """
        Calibrate gravity and friction compensation using static measurements.
        Should be called with the robot in free-space (no contact).
        """
        self.gravity_comp = motor_currents.mean(axis=0) if len(motor_currents.shape) > 1 else motor_currents
        self._calibrated = True

    def estimate_wrench(self, joint_positions: np.ndarray,
                        joint_velocities: np.ndarray,
                        motor_currents: np.ndarray,
                        jacobian: np.ndarray) -> np.ndarray:
        """
        Estimate external wrench from motor currents.

        Returns
        -------
        wrench : (6,) [Fx, Fy, Fz, Tx, Ty, Tz] in end-effector frame
        """
        # Joint torques from motor currents (simplified)
        torque_motor = motor_currents * 0.1  # current-to-torque ratio
        torque_external = torque_motor - self.gravity_comp

        # Friction compensation (Coulomb + viscous)
        torque_friction = (
            0.1 * np.sign(joint_velocities) +
            0.05 * joint_velocities
        )
        torque_external -= torque_friction

        # Map to Cartesian wrench via Jacobian transpose
        J = jacobian
        if J.shape[0] == 6 and J.shape[1] >= self.n_joints:
            J_pinv = np.linalg.pinv(J)
            # tau = J^T * F => F = (J^T)^{-1} * tau = J^{-T} * tau
            try:
                wrench = np.linalg.lstsq(J.T, torque_external[:J.shape[1]], rcond=None)[0]
            except Exception:
                wrench = np.zeros(6)
        else:
            wrench = np.zeros(6)

        return wrench


class ComplianceController:
    """
    Impedance/Admittance control for compliant tissue interaction.

    Impedance mode: Position input → Force output (for position-controlled robots)
    Admittance mode: Force input → Position output (for torque-controlled robots)
    """

    def __init__(self, params: Optional[ImpedanceParams] = None,
                 mode: ControlMode = ControlMode.IMPEDANCE,
                 dt: float = 0.001):
        self.params = params or ImpedanceParams.default_suturing()
        self.mode = mode
        self.dt = dt

        # State
        self.pos_error = np.zeros(6)
        self.vel_error = np.zeros(6)
        self.pos_desired = np.zeros(6)
        self.vel_desired = np.zeros(6)
        self.force_filter = ForceFilter()
        self.force_estimator = ForceEstimator()

    def set_desired_pose(self, position: np.ndarray, orientation: np.ndarray):
        """Set desired end-effector pose."""
        self.pos_desired[:3] = position
        self.pos_desired[3:] = orientation[:3]  # Use rotation vector or Euler

    def compute_impedance(self, pos_current: np.ndarray,
                          vel_current: np.ndarray,
                          force_measured: np.ndarray) -> np.ndarray:
        """
        Impedance control: compute corrected position/velocity command.

        M * ddx + D * dx + K * (x - x_d) = F_ext
        => x_corrected = x_d + compliance_displacement

        Returns
        -------
        corrected_pose : (6,) corrected position + orientation command
        """
        M = self.params.mass
        D = self.params.damping
        K = self.params.stiffness

        force_filtered = self.force_filter.filter(force_measured)

        # Position and velocity errors
        self.pos_error = pos_current - self.pos_desired

        # Compute compliance displacement
        # In steady state: K * delta_x = F_ext => delta_x = F_ext / K
        # Dynamic: M * ddx + D * dx + K * x = F_ext
        acc = (force_filtered - D * self.vel_error - K * self.pos_error) / M

        self.vel_error += acc * self.dt
        displacement = self.vel_error * self.dt

        corrected = self.pos_desired.copy()
        corrected += displacement

        return corrected

    def compute_admittance(self, force_measured: np.ndarray) -> np.ndarray:
        """
        Admittance control: compute velocity command from measured force.

        F_ext = M * ddx + D * dx + K * x
        => dx = (F_ext - K * x) / D  (simplified, quasi-static)

        Returns
        -------
        velocity_command : (6,) Cartesian velocity command
        """
        M = self.params.mass
        D = self.params.damping
        K = self.params.stiffness

        force_filtered = self.force_filter.filter(force_measured)

        acc = (force_filtered - D * self.vel_error - K * self.pos_error) / M
        self.vel_error += acc * self.dt
        self.pos_error += self.vel_error * self.dt

        return self.vel_error.copy()

    def compute(self, pos_current: np.ndarray, vel_current: np.ndarray,
                force_measured: np.ndarray) -> np.ndarray:
        """Dispatch to active control mode."""
        if self.mode == ControlMode.IMPEDANCE:
            return self.compute_impedance(pos_current, vel_current, force_measured)
        elif self.mode == ControlMode.ADMITTANCE:
            return self.compute_admittance(force_measured)
        else:
            # Hybrid: impedance for translation, admittance for rotation
            result = np.zeros(6)
            result[:3] = self.compute_impedance(
                pos_current, vel_current, force_measured
            )[:3]
            result[3:] = self.compute_admittance(force_measured)[3:]
            return result

    def update_params_for_phase(self, phase: str):
        """Switch impedance parameters for current suturing phase."""
        self.params = ImpedanceParams.for_phase(phase)
        self.pos_error = np.zeros(6)
        self.vel_error = np.zeros(6)


class AdaptiveCompliance:
    """
    Adaptive impedance controller that adjusts parameters based on
    tissue stiffness estimation and interaction forces.
    """

    def __init__(self, base_params: Optional[ImpedanceParams] = None):
        self.base_params = base_params or ImpedanceParams.default_suturing()
        self.estimated_tissue_stiffness = np.zeros(3)
        self.force_history: list = []
        self.displacement_history: list = []

    def update_tissue_estimate(self, force: np.ndarray,
                                displacement: np.ndarray):
        """Update tissue stiffness estimate using force-displacement data."""
        self.force_history.append(force[:3].copy())
        self.displacement_history.append(displacement[:3].copy())

        if len(self.force_history) < 10:
            return

        # Least-squares stiffness estimation: F = K * x
        F = np.array(self.force_history[-50:])
        X = np.array(self.displacement_history[-50:])

        for d in range(3):
            x_vals = X[:, d]
            f_vals = F[:, d]
            if np.std(x_vals) > 1e-6:
                self.estimated_tissue_stiffness[d] = (
                    np.sum(f_vals * x_vals) / np.sum(x_vals ** 2)
                )

    def get_adapted_params(self) -> ImpedanceParams:
        """Get impedance parameters adapted to estimated tissue stiffness."""
        adapted = ImpedanceParams(
            mass=self.base_params.mass.copy(),
            damping=self.base_params.damping.copy(),
            stiffness=self.base_params.stiffness.copy()
        )

        # Reduce controller stiffness for stiffer tissue (more compliant)
        for d in range(3):
            k_tissue = self.estimated_tissue_stiffness[d]
            if k_tissue > 100:
                ratio = min(k_tissue / 500.0, 3.0)
                adapted.stiffness[d] /= ratio
                adapted.damping[d] *= np.sqrt(ratio)

        return adapted
