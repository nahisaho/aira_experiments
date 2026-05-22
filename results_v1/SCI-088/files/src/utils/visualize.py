"""
Visualization scripts for traffic simulation results.
Generates publication-quality figures.
"""

import numpy as np
import json
import os

def generate_figures():
    """Generate all figures for the report."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib not available, generating text-based summaries")
        generate_text_summaries()
        return

    os.makedirs("figures", exist_ok=True)

    # Load metrics
    with open("results/metrics_history.json", 'r') as f:
        metrics = json.load(f)

    # --- Figure 1: System Architecture Diagram (text-based) ---
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Integrated Urban Traffic Simulation System Architecture', fontsize=14, fontweight='bold')

    boxes = [
        (1, 6, 3, 1.2, 'SUMO\nMicrosimulation', '#4ECDC4'),
        (5, 6, 3, 1.2, 'IDM/MOBIL\nVehicle Models', '#45B7D1'),
        (9, 6, 3, 1.2, 'Multimodal\nTraffic Gen', '#96CEB4'),
        (1, 3.5, 3, 1.2, 'MAPPO\nSignal Control', '#FF6B6B'),
        (5, 3.5, 3, 1.2, 'Kalman Filter\nDemand Est.', '#FFEAA7'),
        (9, 3.5, 3, 1.2, 'Dynamic\nRerouting', '#DDA0DD'),
        (3, 1, 6, 1.2, 'Performance Metrics & Evaluation', '#B8B8B8'),
    ]

    for x, y, w, h, label, color in boxes:
        rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                                        facecolor=color, edgecolor='#333', linewidth=1.5, alpha=0.85)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=10, fontweight='bold', color='#333')

    # Arrows
    for sx, sy, ex, ey in [(2.5, 6, 2.5, 4.7), (6.5, 6, 6.5, 4.7),
                            (10.5, 6, 10.5, 4.7), (4, 4.1, 5, 4.1),
                            (8, 4.1, 9, 4.1), (6, 3.5, 6, 2.2)]:
        ax.annotate('', xy=(ex, ey), xytext=(sx, sy),
                   arrowprops=dict(arrowstyle='->', color='#555', lw=1.5))

    plt.tight_layout()
    plt.savefig('figures/fig1_system_architecture.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig1_system_architecture.svg', bbox_inches='tight')
    plt.close()

    # --- Figure 2: Time-series metrics ---
    ep1_metrics = metrics[:120]  # Episode 1
    steps = [m['step'] for m in ep1_metrics]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Speed
    axes[0, 0].plot(steps, [m['avg_speed_kmh'] for m in ep1_metrics], color='#2196F3', linewidth=1.5)
    axes[0, 0].set_ylabel('Average Speed (km/h)')
    axes[0, 0].set_xlabel('Simulation Time (s)')
    axes[0, 0].set_title('(a) Network Average Speed')
    axes[0, 0].grid(True, alpha=0.3)

    # Delay
    axes[0, 1].plot(steps, [m['avg_delay_s'] for m in ep1_metrics], color='#F44336', linewidth=1.5)
    axes[0, 1].set_ylabel('Average Delay (s)')
    axes[0, 1].set_xlabel('Simulation Time (s)')
    axes[0, 1].set_title('(b) Vehicle Average Delay')
    axes[0, 1].grid(True, alpha=0.3)

    # Queue
    axes[1, 0].plot(steps, [m['avg_queue'] for m in ep1_metrics], color='#FF9800', linewidth=1.5)
    axes[1, 0].set_ylabel('Average Queue Length (veh)')
    axes[1, 0].set_xlabel('Simulation Time (s)')
    axes[1, 0].set_title('(c) Intersection Queue Length')
    axes[1, 0].grid(True, alpha=0.3)

    # Throughput
    axes[1, 1].plot(steps, [m['total_throughput'] for m in ep1_metrics], color='#4CAF50', linewidth=1.5)
    axes[1, 1].set_ylabel('Throughput (veh/interval)')
    axes[1, 1].set_xlabel('Simulation Time (s)')
    axes[1, 1].set_title('(d) Network Throughput')
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle('Traffic Performance Metrics - Episode 1 (Morning Peak)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/fig2_performance_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig2_performance_metrics.svg', bbox_inches='tight')
    plt.close()

    # --- Figure 3: Demand Profile ---
    hours = np.linspace(0, 24, 100)
    from src.models.demand_estimation import HistoricalDemandProfile
    profile = HistoricalDemandProfile()
    factors = [profile.get_factor(h) for h in hours]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.fill_between(hours, factors, alpha=0.3, color='#2196F3')
    ax.plot(hours, factors, color='#2196F3', linewidth=2)
    ax.axvspan(7, 9, alpha=0.1, color='red', label='Morning Peak')
    ax.axvspan(17, 19, alpha=0.1, color='orange', label='Evening Peak')
    ax.set_xlabel('Hour of Day')
    ax.set_ylabel('Demand Factor')
    ax.set_title('Tokyo Downtown Traffic Demand Profile (Weekday)')
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 2))
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('figures/fig3_demand_profile.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig3_demand_profile.svg', bbox_inches='tight')
    plt.close()

    # --- Figure 4: Network Topology ---
    fig, ax = plt.subplots(figsize=(8, 10))
    rows, cols = 8, 6
    for r in range(rows):
        for c in range(cols):
            x = c * 500
            y = r * 375
            ax.plot(x, y, 'o', markersize=12, color='#F44336', zorder=5)
            ax.text(x, y + 50, f'({r},{c})', ha='center', fontsize=6, color='#666')
            if c < cols - 1:
                ax.plot([x, x + 500], [y, y], '-', color='#333', linewidth=2)
            if r < rows - 1:
                ax.plot([x, x], [y, y + 375], '-', color='#333', linewidth=1.5)

    ax.set_xlabel('East-West Distance (m)')
    ax.set_ylabel('North-South Distance (m)')
    ax.set_title('Tokyo Downtown Grid Network (8×6 = 48 Intersections)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig('figures/fig4_network_topology.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig4_network_topology.svg', bbox_inches='tight')
    plt.close()

    # --- Figure 5: IDM Acceleration Function ---
    from src.models.idm_model import IDMModel, IDMParams
    model = IDMModel(IDMParams())

    gaps = np.linspace(2, 100, 200)
    speeds = [5.0, 10.0, 13.0]
    fig, ax = plt.subplots(figsize=(8, 5))
    for v in speeds:
        accs = [model.acceleration(v, s, 2.0) for s in gaps]
        ax.plot(gaps, accs, linewidth=2, label=f'v = {v:.0f} m/s ({v*3.6:.0f} km/h)')
    ax.axhline(y=0, color='k', linewidth=0.5, linestyle='--')
    ax.set_xlabel('Gap to Leader (m)')
    ax.set_ylabel('Acceleration (m/s²)')
    ax.set_title('IDM Acceleration vs. Gap (Δv = 2.0 m/s)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('figures/fig5_idm_acceleration.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig5_idm_acceleration.svg', bbox_inches='tight')
    plt.close()

    print(f"Generated 5 figures in figures/")


def generate_text_summaries():
    """Fallback: generate text-based summaries when matplotlib unavailable."""
    os.makedirs("figures", exist_ok=True)
    with open("figures/README.md", 'w') as f:
        f.write("# Figures\n\nGenerate with: `python3 -m src.utils.visualize`\n")


if __name__ == "__main__":
    generate_figures()
