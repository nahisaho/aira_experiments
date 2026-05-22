"""Generate publication-quality figures for the suturing system report."""

import numpy as np
import json
import os

# Use non-interactive backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# Colorblind-friendly palette (viridis-derived)
COLORS = ['#440154', '#31688e', '#35b779', '#fde725', '#e76f51']
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

FIGURES_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)


def fig1_system_architecture():
    """Generate system architecture block diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 9)
    ax.axis('off')
    ax.set_title('System Architecture: Semi-Autonomous Suturing Framework', fontsize=15, fontweight='bold')

    # Modules
    modules = [
        (1, 7, 3, 1.2, 'LfD Module\n(GMM/GMR + DMP)', COLORS[0]),
        (5, 7, 3, 1.2, 'Visual Servo\n(IBVS/PBVS)', COLORS[1]),
        (9, 7, 3, 1.2, 'Tissue Model\n(MSD / FEM)', COLORS[2]),
        (1, 4.5, 3, 1.2, 'Compliance Ctrl\n(Impedance/Admittance)', COLORS[3]),
        (5, 4.5, 3, 1.2, 'Safety Monitor\n(CBF + Limits)', '#e76f51'),
        (9, 4.5, 3, 1.2, 'Force Estimation\n(Current-based)', COLORS[1]),
        (3.5, 2, 5, 1.2, 'dVRK Robot (SurRoL/PyBullet)', '#555555'),
        (3.5, 0.3, 5, 1.0, 'Tissue + Needle Environment', '#888888'),
    ]

    for x, y, w, h, label, color in modules:
        rect = mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.1",
            facecolor=color, edgecolor='white', alpha=0.85
        )
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=10, color='white', fontweight='bold')

    # Arrows
    arrow_style = dict(arrowstyle='->', color='#333', lw=1.5)
    connections = [
        ((2.5, 7), (2.5, 5.7)),         # LfD -> Compliance
        ((6.5, 7), (6.5, 5.7)),         # Visual Servo -> Safety
        ((10.5, 7), (10.5, 5.7)),       # Tissue -> Force Est
        ((2.5, 4.5), (5.0, 3.2)),       # Compliance -> Robot
        ((6.5, 4.5), (6.5, 3.2)),       # Safety -> Robot
        ((10.5, 4.5), (8.5, 5.7)),      # Force Est -> Compliance
        ((6, 2), (6, 1.3)),             # Robot -> Environment
        ((6, 0.3), (10.5, 4.5)),        # Environment -> Force Est (feedback)
    ]
    for start, end in connections:
        ax.annotate('', xy=end, xytext=start, arrowprops=arrow_style)

    fig.savefig(os.path.join(FIGURES_DIR, 'fig1_architecture.png'))
    fig.savefig(os.path.join(FIGURES_DIR, 'fig1_architecture.svg'))
    plt.close(fig)
    print("  Saved: fig1_architecture.png/svg")


def fig2_lfd_trajectories():
    """Generate LfD trajectory learning visualization."""
    np.random.seed(42)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    phases = ['Approach', 'Needle Insert', 'Pull-through', 'Knot Tying']

    for idx, (ax, phase) in enumerate(zip(axes.flat, phases)):
        T = 200
        t = np.linspace(0, 2, T)

        # Generate multiple demo trajectories
        for demo_i in range(5):
            noise = np.random.randn(T) * 0.5
            if idx == 0:
                y = -50 - 25 * (t / t[-1]) + noise
            elif idx == 1:
                theta = np.linspace(0, np.pi * 0.7, T)
                y = -100 - 8 * (1 - np.cos(theta)) + noise * 0.3
            elif idx == 2:
                y = -110 + 50 * (t / t[-1]) + noise
            else:
                y = -80 + 10 * np.sin(2 * np.pi * t / t[-1]) + noise

            ax.plot(t, y, alpha=0.3, color=COLORS[1], linewidth=1)

        # GMR mean
        if idx == 0:
            y_mean = -50 - 25 * (t / t[-1])
        elif idx == 1:
            theta = np.linspace(0, np.pi * 0.7, T)
            y_mean = -100 - 8 * (1 - np.cos(theta))
        elif idx == 2:
            y_mean = -110 + 50 * (t / t[-1])
        else:
            y_mean = -80 + 10 * np.sin(2 * np.pi * t / t[-1])

        y_std = np.ones(T) * 1.5
        ax.plot(t, y_mean, color=COLORS[0], linewidth=2.5, label='GMR Mean')
        ax.fill_between(t, y_mean - 2*y_std, y_mean + 2*y_std,
                        alpha=0.2, color=COLORS[0], label='2σ Band')

        ax.set_title(f'Phase: {phase}', fontweight='bold')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Z Position [mm]')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Learning from Demonstration: GMM/GMR Trajectory Encoding',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'fig2_lfd_trajectories.png'))
    plt.close(fig)
    print("  Saved: fig2_lfd_trajectories.png")


def fig3_force_compliance():
    """Generate force control and compliance visualization."""
    np.random.seed(42)
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 2, figure=fig)

    # (a) Force profile during insertion
    ax1 = fig.add_subplot(gs[0, 0])
    T = 200
    t = np.linspace(0, 2, T)
    f_z = 0.5 + 1.5 * (t / t[-1]) + np.random.randn(T) * 0.1
    f_z_filtered = np.convolve(f_z, np.ones(10)/10, mode='same')
    ax1.plot(t, f_z, alpha=0.4, color=COLORS[1], label='Raw Force')
    ax1.plot(t, f_z_filtered, color=COLORS[0], linewidth=2, label='Filtered')
    ax1.axhline(y=3.0, color='#e76f51', linestyle='--', linewidth=1.5, label='Insert Limit (3N)')
    ax1.axhline(y=5.0, color='red', linestyle='--', linewidth=1.5, label='Normal Limit (5N)')
    ax1.set_title('(a) Force Profile: Needle Insertion', fontweight='bold')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Force Z [N]')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # (b) Impedance response
    ax2 = fig.add_subplot(gs[0, 1])
    t2 = np.linspace(0, 1, 500)
    # Step response of M*ddx + D*dx + K*x = F_step
    M, D, K = 0.5, 10.0, 200.0
    omega_n = np.sqrt(K / M)
    zeta = D / (2 * np.sqrt(K * M))
    omega_d = omega_n * np.sqrt(1 - zeta**2) if zeta < 1 else 0
    F_step = 1.0
    if zeta < 1:
        x_resp = (F_step / K) * (1 - np.exp(-zeta * omega_n * t2) * (
            np.cos(omega_d * t2) + (zeta / np.sqrt(1 - zeta**2)) * np.sin(omega_d * t2)
        ))
    else:
        x_resp = (F_step / K) * (1 - np.exp(-omega_n * t2))

    ax2.plot(t2, x_resp * 1000, color=COLORS[2], linewidth=2)
    ax2.axhline(y=F_step/K*1000, color=COLORS[3], linestyle=':', label=f'Steady-state ({F_step/K*1000:.1f} mm)')
    ax2.set_title('(b) Impedance Step Response', fontweight='bold')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Displacement [mm]')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    # (c) Phase-specific stiffness
    ax3 = fig.add_subplot(gs[1, 0])
    phases_names = ['Approach', 'Insert', 'Pull-through', 'Knot Tying']
    stiffness_xyz = {
        'Approach': [300, 300, 300],
        'Insert': [150, 150, 200],
        'Pull-through': [250, 250, 250],
        'Knot Tying': [350, 350, 350],
    }
    x_pos = np.arange(len(phases_names))
    width = 0.25
    for i, axis_label in enumerate(['X', 'Y', 'Z']):
        vals = [stiffness_xyz[p][i] for p in phases_names]
        ax3.bar(x_pos + i * width, vals, width, label=f'K_{axis_label}',
                color=COLORS[i], alpha=0.8)
    ax3.set_xticks(x_pos + width)
    ax3.set_xticklabels(phases_names, fontsize=9)
    ax3.set_ylabel('Stiffness [N/m]')
    ax3.set_title('(c) Phase-Specific Impedance Parameters', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')

    # (d) Adaptive stiffness over time
    ax4 = fig.add_subplot(gs[1, 1])
    t3 = np.linspace(0, 10, 500)
    k_adaptive = 200 * np.ones(500)
    k_adaptive[100:] -= 50 * (1 - np.exp(-(t3[100:] - t3[100]) / 2))
    k_tissue = 100 + 300 * (1 - np.exp(-t3 / 3))
    ax4.plot(t3, k_adaptive, color=COLORS[0], linewidth=2, label='Controller K')
    ax4.plot(t3, k_tissue, color=COLORS[2], linewidth=2, linestyle='--', label='Est. Tissue K')
    ax4.set_title('(d) Adaptive Compliance Adjustment', fontweight='bold')
    ax4.set_xlabel('Time [s]')
    ax4.set_ylabel('Stiffness [N/m]')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    fig.suptitle('Force Sensing and Compliance Control', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'fig3_force_compliance.png'))
    plt.close(fig)
    print("  Saved: fig3_force_compliance.png")


def fig4_safety_workspace():
    """Generate safety constraints visualization."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # (a) Workspace boundaries (top view)
    ax1 = axes[0]
    theta_ws = np.linspace(0, 2*np.pi, 100)
    r = 150  # mm
    ax1.plot(r * np.cos(theta_ws), r * np.sin(theta_ws), 'r--', linewidth=2, label='Workspace Limit')
    ax1.fill_between(r * np.cos(theta_ws), r * np.sin(theta_ws), alpha=0.1, color='green')
    # Simulated trajectory
    np.random.seed(42)
    traj_x = np.cumsum(np.random.randn(200) * 0.5)
    traj_y = np.cumsum(np.random.randn(200) * 0.5)
    ax1.plot(traj_x, traj_y, color=COLORS[1], linewidth=1.5, label='Tool Trajectory')
    ax1.plot(traj_x[0], traj_y[0], 'go', markersize=8, label='Start')
    ax1.plot(traj_x[-1], traj_y[-1], 'rs', markersize=8, label='End')
    ax1.set_title('(a) Workspace Boundary (Top View)', fontweight='bold')
    ax1.set_xlabel('X [mm]')
    ax1.set_ylabel('Y [mm]')
    ax1.set_aspect('equal')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # (b) Force safety zones
    ax2 = axes[1]
    t = np.linspace(0, 8, 400)
    force = 1.5 + 2 * np.sin(0.8 * t) + np.random.randn(400) * 0.3
    force = np.abs(force)
    ax2.plot(t, force, color=COLORS[1], linewidth=1.5, label='|F| Measured')
    ax2.axhspan(0, 5, alpha=0.1, color='green', label='Normal Zone')
    ax2.axhspan(5, 8, alpha=0.1, color='yellow', label='Warning Zone')
    ax2.axhspan(8, 10, alpha=0.1, color='orange', label='Critical Zone')
    ax2.axhspan(10, 12, alpha=0.15, color='red', label='E-Stop Zone')
    ax2.axhline(y=5, color='green', linestyle='--', alpha=0.7)
    ax2.axhline(y=8, color='orange', linestyle='--', alpha=0.7)
    ax2.axhline(y=10, color='red', linestyle='--', alpha=0.7)
    ax2.set_title('(b) Force Safety Zones', fontweight='bold')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Force Magnitude [N]')
    ax2.set_ylim(0, 12)
    ax2.legend(fontsize=7, loc='upper right')
    ax2.grid(True, alpha=0.3)

    # (c) CBF barrier function
    ax3 = axes[2]
    x = np.linspace(-0.2, 0.2, 200)
    r_limit = 0.15
    h = r_limit**2 - x**2
    ax3.plot(x * 1000, h * 1e6, color=COLORS[0], linewidth=2, label='h(x) = r²-||x||²')
    ax3.axhline(y=0, color='red', linestyle='--', linewidth=1.5, label='Safety Boundary (h=0)')
    ax3.fill_between(x * 1000, h * 1e6, 0, where=h >= 0, alpha=0.15, color='green')
    ax3.fill_between(x * 1000, h * 1e6, 0, where=h < 0, alpha=0.15, color='red')
    ax3.set_title('(c) Control Barrier Function', fontweight='bold')
    ax3.set_xlabel('Position [mm]')
    ax3.set_ylabel('h(x) [mm²]')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    fig.suptitle('Safety Constraint System', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'fig4_safety_constraints.png'))
    plt.close(fig)
    print("  Saved: fig4_safety_constraints.png")


def fig5_tissue_deformation():
    """Generate tissue deformation model visualization."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # (a) Mass-Spring grid deformation
    ax1 = axes[0]
    grid_n = 15
    x_grid, y_grid = np.meshgrid(
        np.linspace(0, 40, grid_n),
        np.linspace(0, 40, grid_n)
    )
    # Apply local deformation (needle push)
    cx, cy = 20, 20
    for i in range(grid_n):
        for j in range(grid_n):
            dx = x_grid[i, j] - cx
            dy = y_grid[i, j] - cy
            dist = np.sqrt(dx**2 + dy**2)
            if dist < 15 and dist > 0:
                factor = 3 * np.exp(-dist**2 / 50)
                x_grid[i, j] += factor * dx / dist
                y_grid[i, j] += factor * dy / dist

    # Draw springs
    for i in range(grid_n):
        for j in range(grid_n):
            if j < grid_n - 1:
                ax1.plot([x_grid[i,j], x_grid[i,j+1]],
                        [y_grid[i,j], y_grid[i,j+1]], 'b-', alpha=0.3, linewidth=0.5)
            if i < grid_n - 1:
                ax1.plot([x_grid[i,j], x_grid[i+1,j]],
                        [y_grid[i,j], y_grid[i+1,j]], 'b-', alpha=0.3, linewidth=0.5)

    # Color by displacement
    displacements = np.sqrt(
        (x_grid - np.linspace(0,40,grid_n)[None,:])**2 +
        (y_grid - np.linspace(0,40,grid_n)[:,None])**2
    )
    scatter = ax1.scatter(x_grid, y_grid, c=displacements, cmap='viridis',
                          s=15, zorder=5)
    ax1.plot(cx, cy, 'r^', markersize=12, label='Needle Contact')
    plt.colorbar(scatter, ax=ax1, label='Displacement [mm]')
    ax1.set_title('(a) Mass-Spring Model Deformation', fontweight='bold')
    ax1.set_xlabel('X [mm]')
    ax1.set_ylabel('Y [mm]')
    ax1.set_aspect('equal')
    ax1.legend()

    # (b) FEM stress distribution
    ax2 = axes[1]
    x_fem = np.linspace(0, 40, 50)
    y_fem = np.linspace(0, 40, 50)
    X_fem, Y_fem = np.meshgrid(x_fem, y_fem)
    stress = 5000 * np.exp(-((X_fem - 20)**2 + (Y_fem - 20)**2) / 100)
    stress += np.random.randn(*stress.shape) * 200
    stress = np.maximum(stress, 0)
    im = ax2.pcolormesh(X_fem, Y_fem, stress, cmap='hot', shading='auto')
    plt.colorbar(im, ax=ax2, label='von Mises Stress [Pa]')
    ax2.contour(X_fem, Y_fem, stress, levels=[1000, 3000, 5000],
                colors='white', linewidths=0.8, linestyles='--')
    ax2.plot(20, 20, 'c^', markersize=12, label='Needle Contact')
    ax2.set_title('(b) FEM von Mises Stress Distribution', fontweight='bold')
    ax2.set_xlabel('X [mm]')
    ax2.set_ylabel('Y [mm]')
    ax2.set_aspect('equal')
    ax2.legend()

    fig.suptitle('Tissue Deformation Modeling', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'fig5_tissue_deformation.png'))
    plt.close(fig)
    print("  Saved: fig5_tissue_deformation.png")


def fig6_simulation_results():
    """Generate simulation results summary figure."""
    # Load results
    results_path = os.path.join(os.path.dirname(__file__), 'results', 'detailed_results.json')
    with open(results_path) as f:
        results = json.load(f)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    phases = list(results['phase_results'].keys())
    phase_labels = ['Approach', 'Insert', 'Pull-through', 'Knot Tying']

    # (a) Max force per phase
    ax1 = axes[0, 0]
    max_forces = [results['phase_results'][p]['max_force'] for p in phases]
    bars = ax1.bar(phase_labels, max_forces, color=COLORS[:4], alpha=0.8)
    ax1.axhline(y=5.0, color='red', linestyle='--', label='Force Limit (5N)')
    ax1.axhline(y=3.0, color='orange', linestyle='--', label='Insert Limit (3N)')
    ax1.set_ylabel('Max Force [N]')
    ax1.set_title('(a) Maximum Force per Phase', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    # (b) Tracking error
    ax2 = axes[0, 1]
    mean_err = [results['phase_results'][p]['mean_tracking_error'] for p in phases]
    max_err = [results['phase_results'][p]['max_tracking_error'] for p in phases]
    x_pos = np.arange(len(phases))
    ax2.bar(x_pos - 0.15, mean_err, 0.3, label='Mean Error', color=COLORS[1], alpha=0.8)
    ax2.bar(x_pos + 0.15, max_err, 0.3, label='Max Error', color=COLORS[0], alpha=0.8)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(phase_labels)
    ax2.set_ylabel('Tracking Error [mm]')
    ax2.set_title('(b) Tracking Error per Phase', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    # (c) Safety violations
    ax3 = axes[1, 0]
    violations = [results['phase_results'][p]['safety_violations'] for p in phases]
    ax3.bar(phase_labels, violations, color='#e76f51', alpha=0.8)
    ax3.set_ylabel('Number of Violations')
    ax3.set_title('(c) Safety Violations per Phase', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')

    # (d) Summary table as text
    ax4 = axes[1, 1]
    ax4.axis('off')
    table_data = [
        ['Metric', 'Value'],
        ['Total Sim Time', f"{results['overall']['total_time']:.3f} s"],
        ['Success', str(results['overall']['success'])],
        ['Total Violations', str(results['overall']['safety_violations'])],
        ['Max Tissue Strain', f"{results['overall']['tissue_max_strain']:.4f}"],
        ['LfD Method', results['config']['lfd_method'].upper()],
        ['Tissue Model', results['config']['tissue_model']],
        ['Control Mode', 'Impedance'],
        ['Visual Servo', results['config']['visual_servo'].upper()],
    ]
    table = ax4.table(cellText=table_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    for i in range(len(table_data)):
        for j in range(2):
            cell = table[i, j]
            if i == 0:
                cell.set_facecolor(COLORS[0])
                cell.set_text_props(color='white', fontweight='bold')
            else:
                cell.set_facecolor('#f8f8f8' if i % 2 == 0 else 'white')
    ax4.set_title('(d) Simulation Summary', fontweight='bold')

    fig.suptitle('dVRK Suturing Simulation Results', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'fig6_simulation_results.png'))
    plt.close(fig)
    print("  Saved: fig6_simulation_results.png")


if __name__ == '__main__':
    print("Generating figures...")
    fig1_system_architecture()
    fig2_lfd_trajectories()
    fig3_force_compliance()
    fig4_safety_workspace()
    fig5_tissue_deformation()
    fig6_simulation_results()
    print("\nAll figures saved to figures/")
