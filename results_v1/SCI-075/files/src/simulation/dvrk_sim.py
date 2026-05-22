"""
dVRK Simulation Module
========================
SurRoL / PyBullet-based simulation environment for suturing task verification.
Integrates all subsystems: LfD, tissue model, force control, visual servo, safety.
"""

import numpy as np
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass, field
import time
import json
import os

# Internal module imports
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lfd.gmm_gmr import (
    SuturingLfDPipeline, Demonstration, GaussianMixtureRegression
)
from tissue_model.deformation import (
    TissueModelManager, TissueModelType, TissueProperties, MassSpringModel
)
from force_control.compliance import (
    ComplianceController, ImpedanceParams, ControlMode, AdaptiveCompliance
)
from visual_servo.visual_servo import (
    VisualServoController, VisualServoMode, StereoReconstructor, NeedleTracker
)
from safety.constraints import (
    SafetyMonitor, SafetyLimits, SafetyLevel, SafetyAction,
    VelocityLimiter, ControlBarrierFunction
)


@dataclass
class DVRKState:
    """State of a single dVRK PSM (Patient Side Manipulator)."""
    joint_positions: np.ndarray = field(
        default_factory=lambda: np.zeros(7)
    )
    joint_velocities: np.ndarray = field(
        default_factory=lambda: np.zeros(7)
    )
    cartesian_position: np.ndarray = field(
        default_factory=lambda: np.array([0.0, 0.0, -0.10])
    )
    cartesian_orientation: np.ndarray = field(
        default_factory=lambda: np.array([1.0, 0.0, 0.0, 0.0])
    )
    gripper_angle: float = 0.0
    motor_currents: np.ndarray = field(
        default_factory=lambda: np.zeros(7)
    )
    jacobian: np.ndarray = field(
        default_factory=lambda: np.eye(6, 7)
    )


@dataclass
class SimulationConfig:
    """Simulation configuration."""
    dt: float = 0.001                    # 1 kHz control loop
    sim_duration: float = 60.0           # seconds
    render: bool = False
    gravity: np.ndarray = field(
        default_factory=lambda: np.array([0, 0, -9.81])
    )
    tissue_model_type: str = "mass_spring"
    lfd_method: str = "gmr"
    visual_servo_mode: str = "pbvs"
    enable_safety: bool = True
    log_dir: str = "logs"
    results_dir: str = "results"


class DVRKSimulator:
    """
    Simplified dVRK kinematic simulator.
    In production, this would interface with SurRoL/PyBullet.
    """

    # DH parameters for dVRK PSM (simplified)
    DH_PARAMS = [
        # [a, alpha, d, theta_offset]
        [0.0, np.pi/2, 0.0, 0.0],
        [0.0, -np.pi/2, 0.0, 0.0],
        [0.0, np.pi/2, 0.0, 0.0],  # prismatic
        [0.0, 0.0, 0.4318, 0.0],
        [0.0, -np.pi/2, 0.0, 0.0],
        [0.0, np.pi/2, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]

    RCM_HEIGHT = 0.4  # Remote Center of Motion height

    def __init__(self):
        self.state = DVRKState()
        self.psm2_state = DVRKState()  # Second arm
        self._sim_time = 0.0

    def forward_kinematics(self, joint_pos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute forward kinematics for dVRK PSM.

        Returns (position, orientation_quaternion) in base frame.
        """
        # Simplified FK using DH convention
        T = np.eye(4)
        for i, (a, alpha, d, theta_off) in enumerate(self.DH_PARAMS):
            if i >= len(joint_pos):
                break
            theta = joint_pos[i] + theta_off
            if i == 2:  # Prismatic joint
                d = joint_pos[i]
                theta = 0

            ct, st = np.cos(theta), np.sin(theta)
            ca, sa = np.cos(alpha), np.sin(alpha)

            Ti = np.array([
                [ct, -st*ca, st*sa, a*ct],
                [st, ct*ca, -ct*sa, a*st],
                [0, sa, ca, d],
                [0, 0, 0, 1]
            ])
            T = T @ Ti

        position = T[:3, 3]
        # Extract quaternion from rotation matrix
        R = T[:3, :3]
        quat = self._rotation_to_quaternion(R)

        return position, quat

    def compute_jacobian(self, joint_pos: np.ndarray) -> np.ndarray:
        """Compute geometric Jacobian via numerical differentiation."""
        n_joints = len(joint_pos)
        J = np.zeros((6, n_joints))
        eps = 1e-6

        pos0, quat0 = self.forward_kinematics(joint_pos)

        for i in range(n_joints):
            q_pert = joint_pos.copy()
            q_pert[i] += eps
            pos_pert, quat_pert = self.forward_kinematics(q_pert)

            J[:3, i] = (pos_pert - pos0) / eps
            # Simplified angular Jacobian
            J[3:, i] = 2 * (quat_pert[1:] - quat0[1:]) / eps

        return J

    def inverse_kinematics_velocity(self, cartesian_vel: np.ndarray,
                                     joint_pos: np.ndarray) -> np.ndarray:
        """Resolved-rate IK: dq = J^+ * dx."""
        J = self.compute_jacobian(joint_pos)
        J_pinv = np.linalg.pinv(J)
        return J_pinv @ cartesian_vel

    def step(self, joint_velocity_cmd: np.ndarray, dt: float = 0.001):
        """Step simulation forward."""
        self.state.joint_velocities = joint_velocity_cmd
        self.state.joint_positions += joint_velocity_cmd * dt

        pos, quat = self.forward_kinematics(self.state.joint_positions)
        self.state.cartesian_position = pos
        self.state.cartesian_orientation = quat
        self.state.jacobian = self.compute_jacobian(self.state.joint_positions)

        self._sim_time += dt

    @staticmethod
    def _rotation_to_quaternion(R: np.ndarray) -> np.ndarray:
        """Convert 3x3 rotation matrix to quaternion [w, x, y, z]."""
        tr = np.trace(R)
        if tr > 0:
            s = 2.0 * np.sqrt(tr + 1.0)
            w = 0.25 * s
            x = (R[2, 1] - R[1, 2]) / s
            y = (R[0, 2] - R[2, 0]) / s
            z = (R[1, 0] - R[0, 1]) / s
        elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
            s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
            w = (R[2, 1] - R[1, 2]) / s
            x = 0.25 * s
            y = (R[0, 1] + R[1, 0]) / s
            z = (R[0, 2] + R[2, 0]) / s
        elif R[1, 1] > R[2, 2]:
            s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
            w = (R[0, 2] - R[2, 0]) / s
            x = (R[0, 1] + R[1, 0]) / s
            y = 0.25 * s
            z = (R[1, 2] + R[2, 1]) / s
        else:
            s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
            w = (R[1, 0] - R[0, 1]) / s
            x = (R[0, 2] + R[2, 0]) / s
            y = (R[1, 2] + R[2, 1]) / s
            z = 0.25 * s
        return np.array([w, x, y, z])


class SuturingSimulation:
    """
    Complete suturing simulation integrating all subsystems.
    Manages the full semi-autonomous suturing workflow.
    """

    def __init__(self, config: Optional[SimulationConfig] = None):
        self.config = config or SimulationConfig()

        # Initialize subsystems
        self.robot = DVRKSimulator()
        self.tissue = TissueModelManager(TissueProperties())
        self.lfd = SuturingLfDPipeline(
            method=self.config.lfd_method, n_components=5
        )
        self.compliance = ComplianceController(
            mode=ControlMode.IMPEDANCE, dt=self.config.dt
        )
        self.adaptive_compliance = AdaptiveCompliance()
        self.visual_servo = VisualServoController(
            mode=VisualServoMode.PBVS
        )
        self.safety = SafetyMonitor(SafetyLimits())
        self.velocity_limiter = VelocityLimiter(dt=self.config.dt)
        self.cbf = ControlBarrierFunction(alpha=1.0)

        # Setup CBF barriers
        self.cbf.add_force_barrier(self.safety.limits.force_max_critical)
        self.cbf.add_workspace_barrier(np.zeros(3), self.safety.limits.workspace_radius)
        self.cbf.add_strain_barrier(self.safety.limits.max_tissue_strain)

        # Initialize tissue model
        model_type = TissueModelType.MASS_SPRING
        if self.config.tissue_model_type == "fem":
            model_type = TissueModelType.FEM_LINEAR
        self.tissue.initialize(model_type)

        # Logging
        self.logs: List[Dict] = []
        self.metrics: Dict = {
            'total_time': 0.0,
            'phase_times': {},
            'max_forces': {},
            'tracking_errors': [],
            'safety_violations': 0,
            'tissue_max_strain': 0.0,
            'success': False
        }

    def generate_synthetic_demonstrations(self, n_demos: int = 5):
        """Generate synthetic expert demonstrations for testing."""
        np.random.seed(42)

        for phase in SuturingLfDPipeline.PHASES:
            for i in range(n_demos):
                T = 200
                t = np.linspace(0, 2.0, T)

                if phase == 'approach':
                    # Linear approach to tissue surface
                    start = np.array([0.0, 0.0, -0.05]) + np.random.randn(3) * 0.002
                    end = np.array([0.02, 0.01, -0.10]) + np.random.randn(3) * 0.001
                    positions = np.outer(1 - t/t[-1], start) + np.outer(t/t[-1], end)
                    gripper = np.ones(T) * 0.5

                elif phase == 'insert':
                    # Circular needle insertion arc
                    theta = np.linspace(0, np.pi * 0.7, T)
                    r = 0.008 + np.random.randn() * 0.0005
                    positions = np.column_stack([
                        0.02 + r * np.sin(theta),
                        0.01 * np.ones(T),
                        -0.10 - r * (1 - np.cos(theta))
                    ])
                    gripper = np.ones(T) * 0.8

                elif phase == 'pull_through':
                    # Pull suture thread through
                    start = np.array([0.02, 0.01, -0.11])
                    end = np.array([0.0, 0.0, -0.06])
                    positions = np.outer(1 - t/t[-1], start) + np.outer(t/t[-1], end)
                    gripper = np.ones(T) * 0.9

                else:  # knot_tying
                    # Figure-8 knot-tying motion
                    positions = np.column_stack([
                        0.01 * np.sin(2 * np.pi * t / t[-1]),
                        0.01 * np.cos(2 * np.pi * t / t[-1]),
                        -0.08 * np.ones(T)
                    ])
                    gripper = 0.5 + 0.4 * np.sin(np.pi * t / t[-1])

                # Add noise
                positions += np.random.randn(T, 3) * 0.0003

                # Generate orientations (identity + small perturbations)
                orientations = np.tile([1, 0, 0, 0], (T, 1)).astype(float)
                orientations[:, 1:] += np.random.randn(T, 3) * 0.01
                orientations /= np.linalg.norm(orientations, axis=1, keepdims=True)

                # Simulated forces
                forces = np.random.randn(T, 3) * 0.1
                if phase == 'insert':
                    forces[:, 2] += np.linspace(0, 2.0, T)

                demo = Demonstration(
                    positions=positions,
                    orientations=orientations,
                    gripper_angles=gripper,
                    forces=forces,
                    timestamps=t,
                    metadata={'demo_id': i, 'phase': phase}
                )
                self.lfd.add_demonstration(phase, demo)

    def learn_from_demonstrations(self):
        """Train LfD models on collected demonstrations."""
        self.lfd.learn_all_phases()
        stats = self.lfd.get_phase_statistics()
        self.logs.append({
            'event': 'lfd_training_complete',
            'stats': {k: {sk: float(sv) for sk, sv in v.items()} for k, v in stats.items()},
            'timestamp': time.time()
        })
        return stats

    def run_phase(self, phase: str) -> Dict:
        """
        Execute a single suturing phase with all subsystems active.

        Returns metrics for the phase.
        """
        # Generate reference trajectory from LfD
        t_query = np.linspace(0, 1, 200)
        try:
            ref_trajectory = self.lfd.generate_trajectory(phase, t_query)
        except RuntimeError:
            return {'error': f'No trajectory for phase {phase}'}

        # Setup compliance for this phase
        self.compliance.update_params_for_phase(phase)

        phase_metrics = {
            'phase': phase,
            'duration': 0.0,
            'max_force': 0.0,
            'mean_tracking_error': 0.0,
            'max_tracking_error': 0.0,
            'safety_violations': 0,
            'max_strain': 0.0,
            'steps': 0
        }

        tracking_errors = []
        forces_log = []
        positions_log = []

        for step_idx in range(len(ref_trajectory)):
            # Reference pose
            ref_pos = ref_trajectory[step_idx]

            # Current robot state
            cur_pos = self.robot.state.cartesian_position

            # Simulate force measurement (tissue reaction)
            displacement = cur_pos - ref_pos
            if self.tissue.msd_model:
                nearest = np.argmin(np.linalg.norm(
                    self.tissue.msd_model.rest_positions - cur_pos, axis=1
                ))
                self.tissue.msd_model.apply_displacement(
                    nearest, displacement * 0.1
                )
                self.tissue.step(self.config.dt)

            tissue_deform = self.tissue.get_deformation_at(cur_pos)
            simulated_force = np.zeros(6)
            simulated_force[:3] = -self.compliance.params.stiffness[:3] * displacement * 0.01
            if phase == 'insert':
                simulated_force[2] += np.random.uniform(0.5, 2.0)

            # Visual servo (PBVS)
            target_pose = np.zeros(6)
            target_pose[:3] = ref_pos
            vs_velocity = self.visual_servo.compute_pbvs(
                np.concatenate([cur_pos, np.zeros(3)]), target_pose
            )

            # Compliance control
            self.compliance.set_desired_pose(ref_pos, np.zeros(3))
            corrected_pose = self.compliance.compute(
                np.concatenate([cur_pos, np.zeros(3)]),
                np.zeros(6),
                simulated_force
            )

            # Blend LfD reference with visual servo and compliance
            velocity_cmd = 0.6 * vs_velocity + 0.4 * (corrected_pose - np.concatenate([cur_pos, np.zeros(3)])) / self.config.dt

            # Safety check
            strain = self.tissue.msd_model.get_max_strain() if self.tissue.msd_model else 0.0
            safety_state = self.safety.check(
                simulated_force[:3], velocity_cmd[:3],
                cur_pos, strain, phase
            )

            # CBF safety modification
            cbf_state = {
                'force': simulated_force[:3],
                'position': cur_pos,
                'strain': strain
            }
            velocity_cmd = self.cbf.modify_control(velocity_cmd, cbf_state)

            # Velocity limiting
            velocity_cmd = self.velocity_limiter.limit(
                velocity_cmd, safety_state.velocity_scale
            )

            # Apply to robot
            joint_vel = self.robot.inverse_kinematics_velocity(
                velocity_cmd, self.robot.state.joint_positions
            )
            self.robot.step(joint_vel, self.config.dt)

            # Logging
            tracking_error = np.linalg.norm(cur_pos - ref_pos)
            tracking_errors.append(tracking_error)
            force_mag = np.linalg.norm(simulated_force[:3])
            forces_log.append(force_mag)
            positions_log.append(cur_pos.copy())

            if safety_state.violations:
                phase_metrics['safety_violations'] += len(safety_state.violations)

            phase_metrics['steps'] = step_idx + 1

        # Compute phase metrics
        phase_metrics['duration'] = len(ref_trajectory) * self.config.dt
        phase_metrics['max_force'] = float(max(forces_log)) if forces_log else 0.0
        phase_metrics['mean_tracking_error'] = float(np.mean(tracking_errors)) * 1000  # mm
        phase_metrics['max_tracking_error'] = float(np.max(tracking_errors)) * 1000  # mm
        phase_metrics['max_strain'] = float(
            self.tissue.msd_model.get_max_strain() if self.tissue.msd_model else 0.0
        )

        self.metrics['phase_times'][phase] = phase_metrics['duration']
        self.metrics['max_forces'][phase] = phase_metrics['max_force']

        self.logs.append({
            'event': f'phase_{phase}_complete',
            'metrics': phase_metrics,
            'timestamp': time.time()
        })

        return phase_metrics

    def run_full_suturing(self) -> Dict:
        """
        Execute complete suturing sequence through all phases.
        """
        start_time = time.time()
        all_metrics = {}

        for phase in SuturingLfDPipeline.PHASES:
            phase_result = self.run_phase(phase)
            all_metrics[phase] = phase_result

        total_time = time.time() - start_time
        self.metrics['total_time'] = total_time
        self.metrics['success'] = all(
            m.get('max_force', 0) < self.safety.limits.force_max_critical
            for m in all_metrics.values()
        )
        self.metrics['tissue_max_strain'] = max(
            m.get('max_strain', 0) for m in all_metrics.values()
        )

        # Aggregate safety violations
        self.metrics['safety_violations'] = sum(
            m.get('safety_violations', 0) for m in all_metrics.values()
        )

        self.logs.append({
            'event': 'suturing_complete',
            'total_time': total_time,
            'success': self.metrics['success'],
            'timestamp': time.time()
        })

        return {
            'phase_metrics': all_metrics,
            'overall': self.metrics
        }

    def save_results(self, output_dir: str = '.'):
        """Save simulation results and logs."""
        results_dir = os.path.join(output_dir, 'results')
        logs_dir = os.path.join(output_dir, 'logs')
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)

        # Save metrics
        metrics_path = os.path.join(results_dir, 'simulation_metrics.json')
        with open(metrics_path, 'w') as f:
            json.dump(self.metrics, f, indent=2, default=str)

        # Save logs
        log_path = os.path.join(logs_dir, 'process-log.jsonl')
        with open(log_path, 'a') as f:
            for log_entry in self.logs:
                f.write(json.dumps(log_entry, default=str) + '\n')

        return metrics_path, log_path
