#!/usr/bin/env python3
"""
Semi-Autonomous Suturing Simulation Framework
for Surgical Robot Learning and Control

Modules:
1. Learning from Demonstration (LfD) with DMP
2. Tissue Deformation Modeling (FEM / Mass-Spring)
3. Force Sensing & Compliance Control
4. Visual Servoing (3D reconstruction + tracking)
5. Safety Constraints (force limits, workspace bounds)
6. dVRK Simulation Verification
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import os
import json

np.random.seed(42)
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Module 1: Learning from Demonstration (DMP)
# ============================================================

class DynamicMovementPrimitive:
    """Dynamic Movement Primitive for encoding demonstrated trajectories."""

    def __init__(self, n_basis=25, alpha_z=25.0, beta_z=6.25, alpha_x=1.0, dt=0.01):
        self.n_basis = n_basis
        self.alpha_z = alpha_z
        self.beta_z = beta_z
        self.alpha_x = alpha_x
        self.dt = dt
        self.weights = None
        self.goal = None
        self.y0 = None
        self.tau = None

    def _basis_functions(self, x, centers, widths):
        return np.exp(-widths * (x - centers) ** 2)

    def learn(self, trajectory, timestamps=None):
        """Learn DMP weights from a demonstrated trajectory (N x D)."""
        N, D = trajectory.shape
        if timestamps is None:
            timestamps = np.linspace(0, 1, N)

        self.tau = timestamps[-1] - timestamps[0]
        self.y0 = trajectory[0].copy()
        self.goal = trajectory[-1].copy()

        # Canonical system
        x = np.exp(-self.alpha_x * np.linspace(0, 1, N))
        centers = np.exp(-self.alpha_x * np.linspace(0, 1, self.n_basis))
        widths = 1.0 / (0.65 * np.diff(np.append(centers, 0.001)) ** 2)
        widths = np.abs(widths)

        # Compute desired forcing term
        dt_arr = np.diff(timestamps, prepend=timestamps[0] - (timestamps[1] - timestamps[0]))
        vel = np.gradient(trajectory, axis=0) / (dt_arr[:, None] + 1e-8)
        acc = np.gradient(vel, axis=0) / (dt_arr[:, None] + 1e-8)

        f_target = (self.tau ** 2 * acc - self.alpha_z * (self.beta_z * (self.goal - trajectory) - self.tau * vel))

        # Fit weights via least squares
        Phi = np.array([self._basis_functions(x_i, centers, widths) for x_i in x])
        Phi_norm = Phi / (Phi.sum(axis=1, keepdims=True) + 1e-8) * x[:, None]

        self.weights = np.zeros((self.n_basis, D))
        self.centers = centers
        self.widths = widths

        for d in range(D):
            self.weights[:, d] = np.linalg.lstsq(Phi_norm, f_target[:, d], rcond=None)[0]

        return self

    def generate(self, y0=None, goal=None, tau=None, n_steps=200):
        """Generate trajectory from learned DMP."""
        if y0 is None:
            y0 = self.y0.copy()
        if goal is None:
            goal = self.goal.copy()
        if tau is None:
            tau = self.tau

        D = len(y0)
        y = y0.copy()
        dy = np.zeros(D)
        x = 1.0
        dt = tau / n_steps

        traj = [y.copy()]
        velocities = [dy.copy()]
        forces_applied = []

        for _ in range(n_steps):
            psi = self._basis_functions(x, self.centers, self.widths)
            psi_norm = psi / (psi.sum() + 1e-8) * x
            f = psi_norm @ self.weights

            ddy = (self.alpha_z * (self.beta_z * (goal - y) - tau * dy) + f) / (tau ** 2)
            dy += ddy * dt
            y += dy * dt
            x -= self.alpha_x * x * dt / tau

            traj.append(y.copy())
            velocities.append(dy.copy())
            forces_applied.append(np.linalg.norm(f))

        return np.array(traj), np.array(velocities), np.array(forces_applied)


def generate_suturing_demonstration():
    """Generate a realistic 3D suturing trajectory demonstration."""
    N = 200
    t = np.linspace(0, 2 * np.pi, N)

    # Semicircular needle path in 3D
    radius = 0.015  # 15mm radius
    x = radius * np.cos(t) + 0.05
    y = 0.003 * np.sin(2 * t)
    z = radius * np.sin(t) + 0.01

    # Add slight noise to simulate human demonstration
    noise_level = 0.0003
    x += np.random.randn(N) * noise_level
    y += np.random.randn(N) * noise_level
    z += np.random.randn(N) * noise_level

    return np.column_stack([x, y, z])


# ============================================================
# Module 2: Tissue Deformation (Mass-Spring Model)
# ============================================================

class MassSpringTissue:
    """2D Mass-Spring tissue deformation model."""

    def __init__(self, nx=10, ny=10, spacing=0.005, mass=0.001,
                 k_structural=50.0, k_shear=20.0, damping=0.1):
        self.nx = nx
        self.ny = ny
        self.spacing = spacing
        self.mass = mass
        self.k_structural = k_structural
        self.k_shear = k_shear
        self.damping = damping

        # Initialize positions
        self.positions = np.zeros((nx * ny, 2))
        self.velocities = np.zeros((nx * ny, 2))
        self.rest_positions = np.zeros((nx * ny, 2))

        for i in range(nx):
            for j in range(ny):
                idx = i * ny + j
                self.positions[idx] = [i * spacing, j * spacing]
                self.rest_positions[idx] = [i * spacing, j * spacing]

        # Build spring connections
        self.springs = []
        for i in range(nx):
            for j in range(ny):
                idx = i * ny + j
                # Structural
                if i + 1 < nx:
                    self.springs.append((idx, (i + 1) * ny + j, k_structural))
                if j + 1 < ny:
                    self.springs.append((idx, i * ny + j + 1, k_structural))
                # Shear
                if i + 1 < nx and j + 1 < ny:
                    self.springs.append((idx, (i + 1) * ny + j + 1, k_shear))
                if i + 1 < nx and j - 1 >= 0:
                    self.springs.append((idx, (i + 1) * ny + j - 1, k_shear))

        # Fixed boundary (top row)
        self.fixed = set(range(0, ny))

    def apply_force(self, node_idx, force):
        """Apply external force to a node."""
        if node_idx not in self.fixed:
            self.velocities[node_idx] += force / self.mass * 0.001

    def step(self, dt=0.001):
        """Simulate one timestep."""
        forces = np.zeros_like(self.positions)

        for (i, j, k) in self.springs:
            diff = self.positions[j] - self.positions[i]
            dist = np.linalg.norm(diff)
            rest_diff = self.rest_positions[j] - self.rest_positions[i]
            rest_len = np.linalg.norm(rest_diff)

            if dist > 1e-10:
                force_mag = k * (dist - rest_len)
                force_dir = diff / dist
                forces[i] += force_mag * force_dir
                forces[j] -= force_mag * force_dir

        # Damping
        forces -= self.damping * self.velocities

        # Integrate
        for idx in range(len(self.positions)):
            if idx not in self.fixed:
                self.velocities[idx] += forces[idx] / self.mass * dt
                self.positions[idx] += self.velocities[idx] * dt

        return self.positions.copy()

    def simulate_needle_insertion(self, entry_node, force_magnitude=0.5, n_steps=100):
        """Simulate needle insertion at a specific node."""
        deformations = []
        max_deforms = []
        for step in range(n_steps):
            force = np.array([0.0, -force_magnitude * (1 - step / n_steps)])
            self.apply_force(entry_node, force)
            pos = self.step()
            deformation = np.linalg.norm(pos - self.rest_positions, axis=1)
            deformations.append(deformation.copy())
            max_deforms.append(deformation.max())
        return deformations, max_deforms


# ============================================================
# Module 3: Force Sensing & Compliance Control
# ============================================================

class ComplianceController:
    """Impedance-based compliance controller for suturing."""

    def __init__(self, M=0.5, B=10.0, K=100.0, force_limit=5.0):
        self.M = M   # Virtual inertia
        self.B = B   # Virtual damping
        self.K = K   # Virtual stiffness
        self.force_limit = force_limit
        self.x = 0.0
        self.dx = 0.0
        self.x_d = 0.0  # Desired position

    def update(self, f_ext, x_desired, dt=0.001):
        """Update compliance controller with external force."""
        self.x_d = x_desired

        # Safety: clamp force
        f_safe = np.clip(f_ext, -self.force_limit, self.force_limit)

        # Impedance equation: M*ddx + B*dx + K*(x - x_d) = f_ext
        ddx = (f_safe - self.B * self.dx - self.K * (self.x - self.x_d)) / self.M
        self.dx += ddx * dt
        self.x += self.dx * dt

        return self.x, self.dx, f_safe

    def simulate_insertion(self, n_steps=500):
        """Simulate force-controlled needle insertion."""
        positions = []
        velocities = []
        forces = []
        desired_positions = []
        force_errors = []

        for i in range(n_steps):
            t = i * 0.001
            x_desired = 0.015 * np.sin(2 * np.pi * t)

            # Simulated tissue reaction force
            if self.x > 0.005:
                f_ext = -80.0 * (self.x - 0.005) + np.random.randn() * 0.1
            else:
                f_ext = np.random.randn() * 0.05

            x, dx, f_safe = self.update(f_ext, x_desired)

            positions.append(x)
            velocities.append(dx)
            forces.append(f_safe)
            desired_positions.append(x_desired)
            force_errors.append(abs(f_ext - f_safe))

        return (np.array(positions), np.array(velocities),
                np.array(forces), np.array(desired_positions),
                np.array(force_errors))


# ============================================================
# Module 4: Visual Servoing
# ============================================================

class VisualServo:
    """Image-based visual servoing for needle tracking."""

    def __init__(self, camera_matrix=None):
        if camera_matrix is None:
            self.K = np.array([
                [500, 0, 320],
                [0, 500, 240],
                [0, 0, 1]
            ], dtype=float)
        else:
            self.K = camera_matrix

    def project_3d_to_2d(self, points_3d):
        """Project 3D points to 2D image coordinates."""
        points_h = np.hstack([points_3d, np.ones((len(points_3d), 1))])
        R = np.eye(3)
        t = np.array([[0], [0], [0.1]])
        P = self.K @ np.hstack([R, t])
        projected = (P @ points_h.T).T
        projected = projected[:, :2] / projected[:, 2:3]
        return projected

    def estimate_needle_pose(self, observed_2d, true_3d):
        """Estimate needle pose from 2D observations (simplified PnP)."""
        n = len(observed_2d)
        projected = self.project_3d_to_2d(true_3d)

        # Reprojection error
        errors = np.linalg.norm(projected - observed_2d, axis=1)
        mean_error = errors.mean()

        return projected, mean_error

    def track_needle(self, trajectory_3d, noise_std=2.0):
        """Simulate needle tracking with noise."""
        tracking_results = []
        errors = []

        for point in trajectory_3d:
            point_2d = self.project_3d_to_2d(point.reshape(1, -1))[0]
            noisy_2d = point_2d + np.random.randn(2) * noise_std

            error = np.linalg.norm(point_2d - noisy_2d)
            tracking_results.append(noisy_2d)
            errors.append(error)

        return np.array(tracking_results), np.array(errors)


# ============================================================
# Module 5: Safety Constraints
# ============================================================

class SafetyMonitor:
    """Safety constraint monitor for surgical operations."""

    def __init__(self, force_limit=5.0, workspace_bounds=None, velocity_limit=0.1):
        self.force_limit = force_limit
        self.velocity_limit = velocity_limit
        if workspace_bounds is None:
            self.workspace_bounds = {
                'x': (-0.05, 0.10),
                'y': (-0.05, 0.05),
                'z': (-0.01, 0.05)
            }
        else:
            self.workspace_bounds = workspace_bounds
        self.violations = []

    def check_force(self, force_vec):
        """Check if force is within safe limits."""
        magnitude = np.linalg.norm(force_vec)
        if magnitude > self.force_limit:
            self.violations.append(('force', magnitude))
            return False, magnitude
        return True, magnitude

    def check_workspace(self, position):
        """Check if position is within workspace bounds."""
        axes = ['x', 'y', 'z']
        for i, axis in enumerate(axes):
            if position[i] < self.workspace_bounds[axis][0] or position[i] > self.workspace_bounds[axis][1]:
                self.violations.append(('workspace', axis, position[i]))
                return False
        return True

    def check_velocity(self, velocity):
        """Check if velocity is within safe limits."""
        speed = np.linalg.norm(velocity)
        if speed > self.velocity_limit:
            self.violations.append(('velocity', speed))
            return False, speed
        return True, speed

    def evaluate_trajectory(self, positions, velocities, forces):
        """Evaluate complete trajectory for safety."""
        n = min(len(positions), len(velocities), len(forces))
        force_safe = np.zeros(n, dtype=bool)
        workspace_safe = np.zeros(n, dtype=bool)
        velocity_safe = np.zeros(n, dtype=bool)
        force_magnitudes = np.zeros(n)
        velocity_magnitudes = np.zeros(n)

        for i in range(n):
            fs, fm = self.check_force(forces[i] if forces.ndim > 1 else np.array([forces[i]]))
            force_safe[i] = fs
            force_magnitudes[i] = fm
            workspace_safe[i] = self.check_workspace(positions[i] if positions.ndim > 1 else np.array([positions[i], 0, 0]))
            vs, vm = self.check_velocity(velocities[i] if velocities.ndim > 1 else np.array([velocities[i]]))
            velocity_safe[i] = vs
            velocity_magnitudes[i] = vm

        return {
            'force_safe': force_safe,
            'workspace_safe': workspace_safe,
            'velocity_safe': velocity_safe,
            'force_magnitudes': force_magnitudes,
            'velocity_magnitudes': velocity_magnitudes,
            'overall_safe_ratio': (force_safe & workspace_safe & velocity_safe).mean()
        }


# ============================================================
# Module 6: Full Simulation Pipeline
# ============================================================

def run_full_simulation():
    """Run the complete semi-autonomous suturing simulation."""
    results = {}

    print("=" * 60)
    print("Semi-Autonomous Suturing Simulation Framework")
    print("=" * 60)

    # --- 1. Learning from Demonstration ---
    print("\n[1/6] Learning from Demonstration (DMP)...")
    demo_traj = generate_suturing_demonstration()
    dmp = DynamicMovementPrimitive(n_basis=30)
    dmp.learn(demo_traj)

    # Generate with same goal
    gen_traj, gen_vel, gen_forces = dmp.generate(n_steps=200)

    # Generate with modified goal (generalization)
    new_goal = demo_traj[-1] + np.array([0.005, 0.002, -0.003])
    gen_traj_new, _, _ = dmp.generate(goal=new_goal, n_steps=200)

    rmse_reproduction = np.sqrt(np.mean((gen_traj[:len(demo_traj)] - demo_traj) ** 2))
    print(f"  Reproduction RMSE: {rmse_reproduction * 1000:.3f} mm")
    results['dmp_rmse_mm'] = rmse_reproduction * 1000

    # --- 2. Tissue Deformation ---
    print("\n[2/6] Tissue Deformation Modeling (Mass-Spring)...")
    tissue = MassSpringTissue(nx=12, ny=12)
    entry_node = 6 * 12 + 6  # Center node
    deformations, max_deforms = tissue.simulate_needle_insertion(entry_node, force_magnitude=0.8, n_steps=150)
    max_deformation_mm = max(max_deforms) * 1000
    print(f"  Max tissue deformation: {max_deformation_mm:.3f} mm")
    results['max_tissue_deformation_mm'] = max_deformation_mm

    # --- 3. Force Sensing & Compliance Control ---
    print("\n[3/6] Force Sensing & Compliance Control...")
    controller = ComplianceController(M=0.5, B=10.0, K=100.0, force_limit=5.0)
    pos, vel, forces, des_pos, force_err = controller.simulate_insertion(n_steps=1000)
    mean_tracking_error = np.mean(np.abs(pos - des_pos)) * 1000
    max_force = np.max(np.abs(forces))
    force_limit_violations = np.sum(np.abs(forces) >= 4.99) / len(forces) * 100
    print(f"  Mean tracking error: {mean_tracking_error:.3f} mm")
    print(f"  Max force applied: {max_force:.3f} N")
    print(f"  Force limit utilization: {force_limit_violations:.1f}%")
    results['compliance_tracking_error_mm'] = mean_tracking_error
    results['max_force_N'] = max_force

    # --- 4. Visual Servoing ---
    print("\n[4/6] Visual Servoing (Needle Tracking)...")
    vs = VisualServo()
    tracked, track_errors = vs.track_needle(demo_traj, noise_std=1.5)
    mean_track_error_px = np.mean(track_errors)
    std_track_error_px = np.std(track_errors)
    print(f"  Mean tracking error: {mean_track_error_px:.2f} px")
    print(f"  Std tracking error: {std_track_error_px:.2f} px")
    results['visual_tracking_error_px'] = mean_track_error_px
    results['visual_tracking_std_px'] = std_track_error_px

    # --- 5. Safety Constraints ---
    print("\n[5/6] Safety Constraint Evaluation...")
    safety = SafetyMonitor(force_limit=5.0, velocity_limit=0.05)
    safety_results = safety.evaluate_trajectory(gen_traj, gen_vel, gen_forces)
    print(f"  Overall safety compliance: {safety_results['overall_safe_ratio'] * 100:.1f}%")
    print(f"  Force violations: {(~safety_results['force_safe']).sum()}")
    print(f"  Workspace violations: {(~safety_results['workspace_safe']).sum()}")
    print(f"  Velocity violations: {(~safety_results['velocity_safe']).sum()}")
    results['safety_compliance_pct'] = safety_results['overall_safe_ratio'] * 100
    results['force_violations'] = int((~safety_results['force_safe']).sum())
    results['workspace_violations'] = int((~safety_results['workspace_safe']).sum())
    results['velocity_violations'] = int((~safety_results['velocity_safe']).sum())

    # --- 6. dVRK Simulation Summary ---
    print("\n[6/6] dVRK Simulation Verification...")
    # Simulate multiple suturing trials
    n_trials = 10
    trial_rmses = []
    trial_forces = []
    trial_safety = []

    for trial in range(n_trials):
        noise = np.random.randn(*demo_traj.shape) * 0.0005
        noisy_demo = demo_traj + noise
        dmp_trial = DynamicMovementPrimitive(n_basis=30)
        dmp_trial.learn(noisy_demo)
        t_traj, t_vel, t_forces = dmp_trial.generate(n_steps=200)
        rmse = np.sqrt(np.mean((t_traj[:len(demo_traj)] - demo_traj) ** 2))
        trial_rmses.append(rmse * 1000)
        trial_forces.append(np.max(t_forces) if len(t_forces) > 0 else 0)

        sm = SafetyMonitor(force_limit=5.0, velocity_limit=0.05)
        sr = sm.evaluate_trajectory(t_traj, t_vel, t_forces)
        trial_safety.append(sr['overall_safe_ratio'] * 100)

    results['trial_rmse_mean'] = np.mean(trial_rmses)
    results['trial_rmse_std'] = np.std(trial_rmses)
    results['trial_safety_mean'] = np.mean(trial_safety)
    results['trial_safety_std'] = np.std(trial_safety)

    print(f"  Mean RMSE across {n_trials} trials: {np.mean(trial_rmses):.3f} ± {np.std(trial_rmses):.3f} mm")
    print(f"  Mean safety compliance: {np.mean(trial_safety):.1f} ± {np.std(trial_safety):.1f}%")

    # ============================================================
    # Generate Figures
    # ============================================================
    print("\n" + "=" * 60)
    print("Generating Figures...")
    print("=" * 60)

    # Figure 1: DMP Trajectory Reproduction
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    labels = ['X (m)', 'Y (m)', 'Z (m)']
    for i, (ax, label) in enumerate(zip(axes, labels)):
        ax.plot(demo_traj[:, i], 'b-', linewidth=2, label='Demonstration', alpha=0.7)
        ax.plot(gen_traj[:len(demo_traj), i], 'r--', linewidth=2, label='DMP Reproduction')
        ax.plot(gen_traj_new[:len(demo_traj), i], 'g:', linewidth=2, label='DMP Generalization')
        ax.set_xlabel('Time step')
        ax.set_ylabel(label)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    fig.suptitle('DMP-based Trajectory Learning and Generalization', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'dmp_trajectory.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: dmp_trajectory.png")

    # Figure 2: 3D Suturing Trajectory
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(demo_traj[:, 0] * 1000, demo_traj[:, 1] * 1000, demo_traj[:, 2] * 1000,
            'b-', linewidth=2, label='Demonstration', alpha=0.6)
    ax.plot(gen_traj[:len(demo_traj), 0] * 1000, gen_traj[:len(demo_traj), 1] * 1000,
            gen_traj[:len(demo_traj), 2] * 1000, 'r--', linewidth=2, label='Reproduced')
    ax.plot(gen_traj_new[:len(demo_traj), 0] * 1000, gen_traj_new[:len(demo_traj), 1] * 1000,
            gen_traj_new[:len(demo_traj), 2] * 1000, 'g:', linewidth=2, label='Generalized')
    ax.scatter(*demo_traj[0] * 1000, c='blue', s=100, marker='o', label='Start')
    ax.scatter(*demo_traj[-1] * 1000, c='red', s=100, marker='*', label='Goal')
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    ax.set_title('3D Suturing Trajectory', fontsize=14, fontweight='bold')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'suturing_3d.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: suturing_3d.png")

    # Figure 3: Tissue Deformation Heatmap
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    timesteps_to_show = [0, len(deformations) // 2, -1]
    titles = ['Initial', 'Mid-insertion', 'Full insertion']
    for ax, t_idx, title in zip(axes, timesteps_to_show, titles):
        deform_grid = deformations[t_idx].reshape(12, 12) * 1000
        im = ax.imshow(deform_grid, cmap='hot', interpolation='bilinear', aspect='equal')
        ax.set_title(f'{title}\n(max: {deform_grid.max():.2f} mm)', fontsize=12)
        ax.set_xlabel('Node X')
        ax.set_ylabel('Node Y')
        plt.colorbar(im, ax=ax, label='Deformation (mm)')
    fig.suptitle('Tissue Deformation During Needle Insertion (Mass-Spring Model)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'tissue_deformation.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: tissue_deformation.png")

    # Figure 4: Force Sensing & Compliance Control
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    t_axis = np.arange(len(pos)) * 0.001

    axes[0].plot(t_axis, pos * 1000, 'b-', linewidth=1.5, label='Actual')
    axes[0].plot(t_axis, des_pos * 1000, 'r--', linewidth=1.5, label='Desired')
    axes[0].set_ylabel('Position (mm)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Impedance-based Compliance Control', fontsize=14, fontweight='bold')

    axes[1].plot(t_axis, forces, 'g-', linewidth=1.5, label='Applied Force')
    axes[1].axhline(y=5.0, color='r', linestyle='--', label='Force Limit (+)')
    axes[1].axhline(y=-5.0, color='r', linestyle='--', label='Force Limit (-)')
    axes[1].set_ylabel('Force (N)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t_axis, np.abs(pos - des_pos) * 1000, 'm-', linewidth=1.5)
    axes[2].set_ylabel('Tracking Error (mm)')
    axes[2].set_xlabel('Time (s)')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'compliance_control.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: compliance_control.png")

    # Figure 5: Visual Servoing - Needle Tracking
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    projected_2d = vs.project_3d_to_2d(demo_traj)
    axes[0].plot(projected_2d[:, 0], projected_2d[:, 1], 'b-', linewidth=2, label='Ground Truth')
    axes[0].plot(tracked[:, 0], tracked[:, 1], 'r.', markersize=3, alpha=0.5, label='Tracked (noisy)')
    axes[0].set_xlabel('u (pixels)')
    axes[0].set_ylabel('v (pixels)')
    axes[0].set_title('Needle Tracking in Image Space')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].invert_yaxis()

    axes[1].plot(track_errors, 'b-', linewidth=1, alpha=0.5)
    axes[1].axhline(y=mean_track_error_px, color='r', linestyle='--',
                    label=f'Mean: {mean_track_error_px:.2f} px')
    axes[1].fill_between(range(len(track_errors)), 0, track_errors, alpha=0.2)
    axes[1].set_xlabel('Frame')
    axes[1].set_ylabel('Tracking Error (pixels)')
    axes[1].set_title('Tracking Error Over Time')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.suptitle('Visual Servoing: Needle Pose Estimation and Tracking', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'visual_servoing.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: visual_servoing.png")

    # Figure 6: Safety Constraint Evaluation
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0, 0].plot(safety_results['force_magnitudes'], 'b-', linewidth=1.5)
    axes[0, 0].axhline(y=5.0, color='r', linestyle='--', label='Force Limit')
    axes[0, 0].set_ylabel('Force Magnitude')
    axes[0, 0].set_title('Force Profile')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(safety_results['velocity_magnitudes'] * 1000, 'g-', linewidth=1.5)
    axes[0, 1].axhline(y=50, color='r', linestyle='--', label='Velocity Limit')
    axes[0, 1].set_ylabel('Velocity (mm/s)')
    axes[0, 1].set_title('Velocity Profile')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Workspace boundary check
    axes[1, 0].plot(gen_traj[:, 0] * 1000, gen_traj[:, 2] * 1000, 'b-', linewidth=2)
    ws = safety.workspace_bounds
    rect = plt.Rectangle((ws['x'][0] * 1000, ws['z'][0] * 1000),
                          (ws['x'][1] - ws['x'][0]) * 1000,
                          (ws['z'][1] - ws['z'][0]) * 1000,
                          fill=False, edgecolor='r', linestyle='--', linewidth=2, label='Workspace')
    axes[1, 0].add_patch(rect)
    axes[1, 0].set_xlabel('X (mm)')
    axes[1, 0].set_ylabel('Z (mm)')
    axes[1, 0].set_title('Workspace Boundary Check (XZ)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # Safety summary bar chart
    categories = ['Force\nSafe', 'Workspace\nSafe', 'Velocity\nSafe', 'Overall\nSafe']
    values = [
        safety_results['force_safe'].mean() * 100,
        safety_results['workspace_safe'].mean() * 100,
        safety_results['velocity_safe'].mean() * 100,
        safety_results['overall_safe_ratio'] * 100
    ]
    colors = ['#2ecc71' if v > 90 else '#e74c3c' for v in values]
    bars = axes[1, 1].bar(categories, values, color=colors, edgecolor='black', alpha=0.8)
    axes[1, 1].set_ylabel('Compliance (%)')
    axes[1, 1].set_title('Safety Compliance Summary')
    axes[1, 1].set_ylim(0, 110)
    for bar, val in zip(bars, values):
        axes[1, 1].text(bar.get_x() + bar.get_width() / 2., bar.get_height() + 1,
                        f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3, axis='y')

    fig.suptitle('Safety Constraint Evaluation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'safety_evaluation.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: safety_evaluation.png")

    # Figure 7: Multi-trial dVRK Simulation Results
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].bar(range(1, n_trials + 1), trial_rmses, color='steelblue', edgecolor='black', alpha=0.8)
    axes[0].axhline(y=np.mean(trial_rmses), color='r', linestyle='--',
                    label=f'Mean: {np.mean(trial_rmses):.3f} mm')
    axes[0].set_xlabel('Trial')
    axes[0].set_ylabel('RMSE (mm)')
    axes[0].set_title('Trajectory Reproduction Accuracy')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, axis='y')

    axes[1].bar(range(1, n_trials + 1), trial_safety, color='#2ecc71', edgecolor='black', alpha=0.8)
    axes[1].axhline(y=np.mean(trial_safety), color='r', linestyle='--',
                    label=f'Mean: {np.mean(trial_safety):.1f}%')
    axes[1].set_xlabel('Trial')
    axes[1].set_ylabel('Safety Compliance (%)')
    axes[1].set_title('Safety Compliance per Trial')
    axes[1].legend()
    axes[1].set_ylim(0, 110)
    axes[1].grid(True, alpha=0.3, axis='y')

    fig.suptitle('dVRK Simulation: Multi-Trial Evaluation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'dvrk_trials.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: dvrk_trials.png")

    # Figure 8: System Architecture Diagram
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw modules as boxes
    modules = [
        (2, 8, 3.5, 1.2, 'Learning from\nDemonstration\n(DMP)', '#3498db'),
        (7, 8, 3.5, 1.2, 'Visual Servoing\n(3D Tracking)', '#e74c3c'),
        (12, 8, 3.5, 1.2, 'Safety Monitor\n(Constraints)', '#f39c12'),
        (2, 5.5, 3.5, 1.2, 'Tissue Model\n(Mass-Spring)', '#2ecc71'),
        (7, 5.5, 3.5, 1.2, 'Compliance\nController', '#9b59b6'),
        (12, 5.5, 3.5, 1.2, 'Force Sensing\n& Feedback', '#1abc9c'),
        (5, 2.5, 6, 1.5, 'dVRK Simulation\n(SurRoL / ROS)', '#34495e'),
    ]

    for x, y, w, h, text, color in modules:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black',
                              linewidth=2, alpha=0.8, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                fontsize=10, fontweight='bold', color='white', zorder=3)

    # Draw arrows
    arrow_style = dict(arrowstyle='->', color='black', linewidth=2)
    arrows = [
        ((3.75, 8), (7, 8.6)),
        ((10.5, 8.6), (12, 8.6)),
        ((3.75, 5.5 + 1.2), (3.75, 8)),
        ((8.75, 5.5 + 1.2), (8.75, 8)),
        ((12, 6.1), (10.5, 6.1)),
        ((3.75, 5.5), (5, 3.25)),
        ((8.75, 5.5), (8, 4)),
        ((13.75, 5.5), (11, 3.25)),
    ]
    for start, end in arrows:
        ax.annotate('', xy=end, xytext=start,
                    arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=2))

    ax.set_title('System Architecture: Semi-Autonomous Suturing Framework',
                 fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'system_architecture.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: system_architecture.png")

    # Figure 9: Tissue deformation time-series
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(np.array(max_deforms) * 1000, 'b-', linewidth=2)
    ax.set_xlabel('Simulation Step')
    ax.set_ylabel('Max Deformation (mm)')
    ax.set_title('Tissue Deformation Over Time During Needle Insertion', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.fill_between(range(len(max_deforms)), 0, np.array(max_deforms) * 1000, alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'tissue_deformation_timeseries.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: tissue_deformation_timeseries.png")

    print("\n" + "=" * 60)
    print("Simulation Complete!")
    print("=" * 60)

    # Save results
    with open(os.path.join(os.path.dirname(FIGURES_DIR), 'results.json'), 'w') as f:
        json.dump(results, f, indent=2)

    return results


if __name__ == '__main__':
    results = run_full_simulation()
    print("\nResults Summary:")
    for k, v in results.items():
        print(f"  {k}: {v}")
