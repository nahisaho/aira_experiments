"""
Visualization Generator
========================
Creates publication-quality figures for all modules.
"""

import numpy as np
import json, os

def create_matplotlib_figures():
    """Generate all figures using matplotlib."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        print("matplotlib not available; skipping figure generation.")
        return False

    os.makedirs("figures", exist_ok=True)
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'figure.dpi': 150,
    })

    # --- Figure 1: CFD Results ---
    with open("results/cfd_results.json") as f:
        cfd = json.load(f)
    with open("results/velocity_profile_data.json") as f:
        vp = json.load(f)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    conditions = cfd["flow_conditions"]
    Qs = [c["flow_rate_mL_min"] for c in conditions]
    Res = [c["reynolds_number"] for c in conditions]
    dPs = [c["pressure_drop_kPa"] for c in conditions]
    mix_effs = [c["mixing_efficiency"] for c in conditions]
    taus = [c["residence_time_s"] for c in conditions]

    axes[0, 0].plot(Qs, Res, 'o-', color='#2196F3', linewidth=2, markersize=8)
    axes[0, 0].set_xlabel("Flow Rate (mL/min)")
    axes[0, 0].set_ylabel("Reynolds Number")
    axes[0, 0].set_title("Reynolds Number vs Flow Rate")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(Qs, dPs, 's-', color='#F44336', linewidth=2, markersize=8)
    axes[0, 1].set_xlabel("Flow Rate (mL/min)")
    axes[0, 1].set_ylabel("Pressure Drop (kPa)")
    axes[0, 1].set_title("Pressure Drop vs Flow Rate")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(vp["y_um"], vp["u_at_midplane_m_s"], '-', color='#4CAF50', linewidth=2)
    axes[1, 0].set_xlabel("y Position (μm)")
    axes[1, 0].set_ylabel("Velocity (m/s)")
    axes[1, 0].set_title("Velocity Profile at Midplane (Q=1 mL/min)")
    axes[1, 0].grid(True, alpha=0.3)

    ax_twin = axes[1, 1]
    color1 = '#FF9800'
    color2 = '#9C27B0'
    ax_twin.plot(Qs, mix_effs, 'D-', color=color1, linewidth=2, markersize=8, label='Mixing Efficiency')
    ax_twin.set_xlabel("Flow Rate (mL/min)")
    ax_twin.set_ylabel("Mixing Efficiency", color=color1)
    ax_twin.tick_params(axis='y', labelcolor=color1)
    ax2 = ax_twin.twinx()
    ax2.plot(Qs, taus, '^-', color=color2, linewidth=2, markersize=8, label='Residence Time')
    ax2.set_ylabel("Residence Time (s)", color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    axes[1, 1].set_title("Mixing & Residence Time vs Flow Rate")
    ax_twin.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig1_cfd_analysis.png", dpi=300, bbox_inches='tight')
    plt.savefig("figures/fig1_cfd_analysis.svg", bbox_inches='tight')
    plt.close()

    # --- Figure 2: RTD Analysis ---
    with open("data/rtd_curves.json") as f:
        rtd = json.load(f)
    with open("results/rtd_results.json") as f:
        rtd_res = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    t = np.array(rtd["time_s"])

    axes[0].plot(t, rtd["E_experimental"], '.', color='gray', alpha=0.3, markersize=2, label='Experimental')
    axes[0].plot(t, rtd["E_true"], '-', color='black', linewidth=2, label='True RTD')
    axes[0].plot(t, rtd["E_axial_dispersion"], '--', color='#2196F3', linewidth=2, label='Axial Dispersion')
    axes[0].plot(t, rtd["E_tanks_in_series"], '-.', color='#F44336', linewidth=1.5, label='Tanks-in-Series')
    axes[0].plot(t, rtd["E_cstr"], ':', color='#4CAF50', linewidth=1.5, label='Single CSTR')
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("E(t) (1/s)")
    axes[0].set_title("Residence Time Distribution")
    axes[0].legend(fontsize=9)
    axes[0].set_xlim(0, 90)
    axes[0].grid(True, alpha=0.3)

    models = list(rtd_res["model_comparison"].keys())
    r2_vals = [rtd_res["model_comparison"][m]["R_squared"] for m in models]
    colors = ['#4CAF50', '#F44336', '#2196F3']
    bars = axes[1].bar(models, r2_vals, color=colors, alpha=0.8, edgecolor='black')
    axes[1].set_ylabel("R² Score")
    axes[1].set_title("Model Fit Comparison")
    axes[1].set_ylim(0.5, 1.05)
    for bar, val in zip(bars, r2_vals):
        axes[1].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                     f'{val:.4f}', ha='center', fontsize=10)
    axes[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig("figures/fig2_rtd_analysis.png", dpi=300, bbox_inches='tight')
    plt.savefig("figures/fig2_rtd_analysis.svg", bbox_inches='tight')
    plt.close()

    # --- Figure 3: Bayesian Optimization ---
    with open("results/bayesian_optimization_results.json") as f:
        bo = json.load(f)

    history = bo["optimization_history"]
    iters = [h["iteration"] for h in history]
    yields = [h["observed_yield"] for h in history]
    best_yields = [h["best_yield_so_far"] for h in history]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].scatter(iters, yields, c='#2196F3', alpha=0.6, s=40, label='Observed')
    axes[0].plot(iters, best_yields, '-', color='#F44336', linewidth=2, label='Best so far')
    axes[0].set_xlabel("BO Iteration")
    axes[0].set_ylabel("Yield (%)")
    axes[0].set_title("Bayesian Optimization Convergence")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Parameter evolution
    temps = [h["parameters"]["temperature_C"] for h in history]
    flows = [h["parameters"]["flow_rate_mL_min"] for h in history]
    ax_t = axes[1]
    ax_t.scatter(temps, flows, c=yields, cmap='viridis', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    cbar = plt.colorbar(ax_t.collections[0], ax=ax_t)
    cbar.set_label("Yield (%)")
    ax_t.set_xlabel("Temperature (°C)")
    ax_t.set_ylabel("Flow Rate (mL/min)")
    ax_t.set_title("Parameter Space Exploration")
    ax_t.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig3_bayesian_optimization.png", dpi=300, bbox_inches='tight')
    plt.savefig("figures/fig3_bayesian_optimization.svg", bbox_inches='tight')
    plt.close()

    # --- Figure 4: Feedback Control ---
    with open("data/control_timeseries.json") as f:
        ts = json.load(f)

    times = [d["time_s"] for d in ts]
    conversions = [d["ir_reading_pct"] for d in ts]
    temps_ctrl = [d["temperature_C"] for d in ts]
    flows_ctrl = [d["flow_rate_mL_min"] for d in ts]

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

    axes[0].plot(times, conversions, '-', color='#2196F3', linewidth=1)
    axes[0].axhline(y=85, color='red', linestyle='--', linewidth=1, label='Setpoint (85%)')
    axes[0].fill_between(times, 83, 87, alpha=0.1, color='green', label='±2% band')
    axes[0].set_ylabel("Conversion (%)")
    axes[0].set_title("Process Control Performance")
    axes[0].legend(loc='lower right')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(times, temps_ctrl, '-', color='#F44336', linewidth=1)
    axes[1].set_ylabel("Temperature (°C)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(times, flows_ctrl, '-', color='#4CAF50', linewidth=1)
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Flow Rate (mL/min)")
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig4_feedback_control.png", dpi=300, bbox_inches='tight')
    plt.savefig("figures/fig4_feedback_control.svg", bbox_inches='tight')
    plt.close()

    # --- Figure 5: Scale-Up Comparison ---
    with open("results/scaleup_results.json") as f:
        su = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    nu = su["numbering_up"]
    n_vals = [r["n_parallel"] for r in nu]
    tp_nu = [r["throughput_g_h"] for r in nu]
    cost_nu = [r["cost_per_g_h"] for r in nu]

    axes[0].bar(range(len(n_vals)), tp_nu, color='#2196F3', alpha=0.8, edgecolor='black')
    axes[0].set_xticks(range(len(n_vals)))
    axes[0].set_xticklabels([str(n) for n in n_vals])
    axes[0].set_xlabel("Number of Parallel Reactors")
    axes[0].set_ylabel("Throughput (g/h)")
    axes[0].set_title("Numbering Up: Throughput")
    axes[0].grid(True, alpha=0.3, axis='y')

    sc = su["scaling_up"]
    sf_vals = [r["scale_factor"] for r in sc]
    ht_vals = [r["heat_transfer_W_m2K"] for r in sc]
    mix_vals = [r["mixing_time_ms"] for r in sc]

    ax1 = axes[1]
    ax1.plot(sf_vals, ht_vals, 'o-', color='#F44336', linewidth=2, markersize=8, label='Heat Transfer')
    ax1.set_xlabel("Scale Factor")
    ax1.set_ylabel("Heat Transfer (W/m²K)", color='#F44336')
    ax1.tick_params(axis='y', labelcolor='#F44336')
    ax3 = ax1.twinx()
    ax3.plot(sf_vals, mix_vals, 's-', color='#4CAF50', linewidth=2, markersize=8, label='Mixing Time')
    ax3.set_ylabel("Mixing Time (ms)", color='#4CAF50')
    ax3.tick_params(axis='y', labelcolor='#4CAF50')
    axes[1].set_title("Scaling Up: Trade-offs")
    ax1.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig5_scaleup_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig("figures/fig5_scaleup_comparison.svg", bbox_inches='tight')
    plt.close()

    print("All figures saved to figures/")
    return True

if __name__ == "__main__":
    create_matplotlib_figures()
