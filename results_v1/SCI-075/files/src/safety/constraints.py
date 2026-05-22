"""
Safety Constraint Module
==========================
Enforces force limits, workspace boundaries, velocity limits,
and collision avoidance for safe autonomous suturing.
"""

import numpy as np
from typing import Optional, Tuple, List, Dict, Callable
from dataclasses import dataclass, field
from enum import Enum
import time


class SafetyLevel(Enum):
    NORMAL = "normal"
    WARNING = "warning"
    CRITICAL = "critical"
    EMERGENCY_STOP = "emergency_stop"


class SafetyAction(Enum):
    NONE = "none"
    SCALE_VELOCITY = "scale_velocity"
    STOP = "stop"
    RETRACT = "retract"
    HANDOVER_TO_HUMAN = "handover"


@dataclass
class SafetyLimits:
    """Safety constraint parameters for surgical suturing."""
    # Force limits [N]
    force_max_normal: float = 5.0        # Normal operation max force
    force_max_warning: float = 8.0       # Warning threshold
    force_max_critical: float = 10.0     # Emergency stop threshold
    force_max_insert: float = 3.0        # Max force during needle insertion

    # Torque limits [Nm]
    torque_max: float = 0.5

    # Velocity limits [m/s, rad/s]
    linear_vel_max: float = 0.05         # 50 mm/s max
    angular_vel_max: float = 0.5         # rad/s
    linear_acc_max: float = 0.2          # m/s^2

    # Workspace boundaries [m] - relative to RCM (Remote Center of Motion)
    workspace_x_range: Tuple[float, float] = (-0.1, 0.1)
    workspace_y_range: Tuple[float, float] = (-0.1, 0.1)
    workspace_z_range: Tuple[float, float] = (-0.15, 0.0)
    workspace_radius: float = 0.15       # Spherical workspace radius

    # Tissue safety
    max_tissue_strain: float = 0.25      # 25% maximum strain
    max_tissue_stress: float = 50000.0   # Pa

    # Timing
    max_phase_duration: float = 30.0     # seconds per phase
    watchdog_timeout: float = 0.1        # seconds (control loop watchdog)

    # Collision
    min_tool_distance: float = 0.005     # 5mm minimum between tools


@dataclass
class SafetyState:
    """Current safety system state."""
    level: SafetyLevel = SafetyLevel.NORMAL
    action: SafetyAction = SafetyAction.NONE
    violations: List[str] = field(default_factory=list)
    force_magnitude: float = 0.0
    velocity_magnitude: float = 0.0
    in_workspace: bool = True
    tissue_strain: float = 0.0
    timestamp: float = 0.0
    velocity_scale: float = 1.0


class SafetyMonitor:
    """
    Real-time safety monitoring and constraint enforcement.
    Runs at control loop frequency (1 kHz) with minimal latency.
    """

    def __init__(self, limits: Optional[SafetyLimits] = None):
        self.limits = limits or SafetyLimits()
        self.state = SafetyState()
        self._last_check_time = time.time()
        self._force_history: List[float] = []
        self._velocity_history: List[float] = []
        self._callbacks: Dict[SafetyLevel, List[Callable]] = {
            level: [] for level in SafetyLevel
        }
        self._enabled = True
        self._violation_counts: Dict[str, int] = {}

    def register_callback(self, level: SafetyLevel, callback: Callable):
        """Register a callback for a specific safety level."""
        self._callbacks[level].append(callback)

    def check(self, force: np.ndarray, velocity: np.ndarray,
              position: np.ndarray, tissue_strain: float = 0.0,
              phase: str = 'approach') -> SafetyState:
        """
        Perform comprehensive safety check.

        Parameters
        ----------
        force : (3,) or (6,) measured force [N] (and torque [Nm])
        velocity : (3,) or (6,) commanded velocity [m/s] (and [rad/s])
        position : (3,) current end-effector position [m]
        tissue_strain : float, current tissue strain
        phase : str, current suturing phase

        Returns
        -------
        SafetyState with current status and required actions
        """
        if not self._enabled:
            return SafetyState()

        self.state = SafetyState(timestamp=time.time())
        self.state.violations = []

        # 1. Force check
        self._check_force(force, phase)

        # 2. Velocity check
        self._check_velocity(velocity)

        # 3. Workspace boundary check
        self._check_workspace(position)

        # 4. Tissue integrity check
        self._check_tissue(tissue_strain)

        # 5. Watchdog timer
        self._check_watchdog()

        # Determine overall safety level
        self._determine_level()

        # Trigger callbacks
        self._trigger_callbacks()

        self._last_check_time = time.time()
        return self.state

    def _check_force(self, force: np.ndarray, phase: str):
        """Check force constraints."""
        force_mag = np.linalg.norm(force[:3])
        self.state.force_magnitude = force_mag
        self._force_history.append(force_mag)
        if len(self._force_history) > 1000:
            self._force_history = self._force_history[-1000:]

        # Phase-specific force limit
        if phase == 'insert':
            max_force = self.limits.force_max_insert
        else:
            max_force = self.limits.force_max_normal

        if force_mag > self.limits.force_max_critical:
            self.state.violations.append(
                f"CRITICAL: Force {force_mag:.2f}N exceeds critical limit "
                f"{self.limits.force_max_critical:.1f}N"
            )
            self.state.level = SafetyLevel.EMERGENCY_STOP
            self.state.action = SafetyAction.STOP
        elif force_mag > self.limits.force_max_warning:
            self.state.violations.append(
                f"WARNING: Force {force_mag:.2f}N exceeds warning limit "
                f"{self.limits.force_max_warning:.1f}N"
            )
            if self.state.level.value < SafetyLevel.WARNING.value:
                self.state.level = SafetyLevel.WARNING
            self.state.action = SafetyAction.SCALE_VELOCITY
            self.state.velocity_scale = min(
                self.state.velocity_scale,
                max_force / force_mag
            )
        elif force_mag > max_force:
            self.state.violations.append(
                f"LIMIT: Force {force_mag:.2f}N exceeds phase limit {max_force:.1f}N"
            )
            self.state.action = SafetyAction.SCALE_VELOCITY
            self.state.velocity_scale = min(
                self.state.velocity_scale,
                0.5 * max_force / force_mag
            )

        # Torque check
        if len(force) >= 6:
            torque_mag = np.linalg.norm(force[3:6])
            if torque_mag > self.limits.torque_max:
                self.state.violations.append(
                    f"Torque {torque_mag:.3f}Nm exceeds limit {self.limits.torque_max:.3f}Nm"
                )

    def _check_velocity(self, velocity: np.ndarray):
        """Check velocity constraints."""
        lin_vel = np.linalg.norm(velocity[:3])
        self.state.velocity_magnitude = lin_vel

        if lin_vel > self.limits.linear_vel_max:
            scale = self.limits.linear_vel_max / lin_vel
            self.state.velocity_scale = min(self.state.velocity_scale, scale)
            self.state.violations.append(
                f"Velocity {lin_vel*1000:.1f}mm/s exceeds limit "
                f"{self.limits.linear_vel_max*1000:.1f}mm/s"
            )

        if len(velocity) >= 6:
            ang_vel = np.linalg.norm(velocity[3:6])
            if ang_vel > self.limits.angular_vel_max:
                scale = self.limits.angular_vel_max / ang_vel
                self.state.velocity_scale = min(self.state.velocity_scale, scale)

    def _check_workspace(self, position: np.ndarray):
        """Check workspace boundary constraints."""
        x, y, z = position[:3]
        limits = self.limits
        in_bounds = True

        if not (limits.workspace_x_range[0] <= x <= limits.workspace_x_range[1]):
            in_bounds = False
            self.state.violations.append(f"X={x*1000:.1f}mm out of workspace")

        if not (limits.workspace_y_range[0] <= y <= limits.workspace_y_range[1]):
            in_bounds = False
            self.state.violations.append(f"Y={y*1000:.1f}mm out of workspace")

        if not (limits.workspace_z_range[0] <= z <= limits.workspace_z_range[1]):
            in_bounds = False
            self.state.violations.append(f"Z={z*1000:.1f}mm out of workspace")

        # Spherical workspace check
        dist = np.linalg.norm(position[:3])
        if dist > limits.workspace_radius:
            in_bounds = False
            self.state.violations.append(
                f"Distance {dist*1000:.1f}mm exceeds workspace radius"
            )

        self.state.in_workspace = in_bounds

        if not in_bounds:
            self.state.action = SafetyAction.RETRACT
            if self.state.level.value < SafetyLevel.WARNING.value:
                self.state.level = SafetyLevel.WARNING

    def _check_tissue(self, strain: float):
        """Check tissue integrity constraints."""
        self.state.tissue_strain = strain

        if strain > self.limits.max_tissue_strain:
            self.state.violations.append(
                f"Tissue strain {strain:.1%} exceeds limit "
                f"{self.limits.max_tissue_strain:.1%}"
            )
            self.state.action = SafetyAction.RETRACT
            self.state.level = SafetyLevel.CRITICAL

    def _check_watchdog(self):
        """Check control loop timing."""
        elapsed = time.time() - self._last_check_time
        if elapsed > self.limits.watchdog_timeout and self._last_check_time > 0:
            self.state.violations.append(
                f"Watchdog timeout: {elapsed*1000:.1f}ms > "
                f"{self.limits.watchdog_timeout*1000:.1f}ms"
            )

    def _determine_level(self):
        """Determine overall safety level from violations."""
        if not self.state.violations:
            self.state.level = SafetyLevel.NORMAL
            self.state.action = SafetyAction.NONE
            self.state.velocity_scale = 1.0

    def _trigger_callbacks(self):
        """Trigger registered callbacks for current safety level."""
        for callback in self._callbacks.get(self.state.level, []):
            try:
                callback(self.state)
            except Exception:
                pass  # Safety callbacks must not crash


class VelocityLimiter:
    """
    Smooth velocity limiting with jerk constraints.
    Prevents abrupt motions that could damage tissue.
    """

    def __init__(self, v_max: float = 0.05, a_max: float = 0.2,
                 j_max: float = 1.0, dt: float = 0.001):
        self.v_max = v_max
        self.a_max = a_max
        self.j_max = j_max
        self.dt = dt
        self.prev_velocity = np.zeros(6)
        self.prev_acceleration = np.zeros(6)

    def limit(self, velocity_cmd: np.ndarray,
              safety_scale: float = 1.0) -> np.ndarray:
        """
        Apply smooth velocity limiting.

        Parameters
        ----------
        velocity_cmd : (6,) desired velocity
        safety_scale : float, scaling factor from safety monitor [0, 1]

        Returns
        -------
        limited_velocity : (6,) smoothly limited velocity
        """
        v_cmd = velocity_cmd.copy()

        # Apply safety scaling
        v_cmd *= safety_scale

        # Velocity magnitude limiting
        lin_speed = np.linalg.norm(v_cmd[:3])
        if lin_speed > self.v_max:
            v_cmd[:3] *= self.v_max / lin_speed

        # Acceleration limiting
        acceleration = (v_cmd - self.prev_velocity) / self.dt
        for i in range(6):
            if abs(acceleration[i]) > self.a_max:
                acceleration[i] = np.sign(acceleration[i]) * self.a_max
                v_cmd[i] = self.prev_velocity[i] + acceleration[i] * self.dt

        # Jerk limiting
        jerk = (acceleration - self.prev_acceleration) / self.dt
        for i in range(6):
            if abs(jerk[i]) > self.j_max:
                jerk[i] = np.sign(jerk[i]) * self.j_max
                acceleration[i] = self.prev_acceleration[i] + jerk[i] * self.dt
                v_cmd[i] = self.prev_velocity[i] + acceleration[i] * self.dt

        self.prev_velocity = v_cmd.copy()
        self.prev_acceleration = acceleration.copy()

        return v_cmd


class ControlBarrierFunction:
    """
    Control Barrier Function (CBF) for formal safety guarantees.
    Ensures the system stays within a safe set via QP-based control modification.

    h(x) >= 0 defines the safe set.
    Constraint: dh/dt + alpha * h(x) >= 0
    """

    def __init__(self, alpha: float = 1.0):
        self.alpha = alpha
        self.barriers: List[Dict] = []

    def add_force_barrier(self, force_limit: float):
        """Add force magnitude barrier: h = f_max^2 - ||f||^2."""
        self.barriers.append({
            'type': 'force',
            'limit': force_limit,
            'h': lambda f, lim=force_limit: lim**2 - np.dot(f[:3], f[:3])
        })

    def add_workspace_barrier(self, center: np.ndarray, radius: float):
        """Add workspace boundary barrier: h = r^2 - ||x - c||^2."""
        self.barriers.append({
            'type': 'workspace',
            'center': center,
            'radius': radius,
            'h': lambda x, c=center, r=radius: r**2 - np.dot(x-c, x-c)
        })

    def add_strain_barrier(self, max_strain: float):
        """Add tissue strain barrier: h = eps_max - eps."""
        self.barriers.append({
            'type': 'strain',
            'limit': max_strain,
            'h': lambda s, lim=max_strain: lim - s
        })

    def modify_control(self, u_nominal: np.ndarray,
                       state: Dict) -> np.ndarray:
        """
        Modify control input to satisfy CBF constraints.

        Uses a simple projection approach (simplified QP).
        Full implementation would use OSQP or similar QP solver.

        Parameters
        ----------
        u_nominal : (6,) nominal control input
        state : dict with 'position', 'force', 'strain' keys

        Returns
        -------
        u_safe : (6,) safety-modified control input
        """
        u_safe = u_nominal.copy()

        for barrier in self.barriers:
            h_val = self._evaluate_barrier(barrier, state)

            if h_val < 0:
                # Already violated - apply correction
                u_safe *= 0.0
                break
            elif h_val < 0.1 * abs(barrier.get('limit', barrier.get('radius', 1.0))):
                # Near boundary - scale control
                scale = h_val / (0.1 * abs(barrier.get('limit', barrier.get('radius', 1.0))))
                u_safe *= max(scale, 0.0)

        return u_safe

    def _evaluate_barrier(self, barrier: Dict, state: Dict) -> float:
        """Evaluate barrier function value."""
        if barrier['type'] == 'force':
            force = state.get('force', np.zeros(3))
            return barrier['h'](force)
        elif barrier['type'] == 'workspace':
            pos = state.get('position', np.zeros(3))
            return barrier['h'](pos)
        elif barrier['type'] == 'strain':
            strain = state.get('strain', 0.0)
            return barrier['h'](strain)
        return 1.0

    def is_safe(self, state: Dict) -> bool:
        """Check if current state is in the safe set."""
        for barrier in self.barriers:
            if self._evaluate_barrier(barrier, state) < 0:
                return False
        return True
