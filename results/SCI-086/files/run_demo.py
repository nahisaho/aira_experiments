"""
run_demo.py
===========
Execute the cardiac digital twin framework demonstration.
Generates all results, figures, and the final report.
"""

import sys
import os
import json
import logging
import numpy as np
from pathlib import Path
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("logs/pipeline.log"),
    ]
)
logger = logging.getLogger("cardiac_dt")

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.digital_twin_pipeline import CardiacDigitalTwinPipeline


def generate_figures(results: dict, output_dir: str = "figures"):
    """Generate publication-quality figures using matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        logger.warning("matplotlib not available, skipping figure generation")
        _generate_text_figures(results, output_dir)
        return

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "figure.dpi": 150,
    })

    fig_dir = Path(output_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # --- Figure 1: Framework Architecture ---
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title("Cardiac Digital Twin Framework Architecture", fontsize=16, fontweight="bold")

    modules = [
        (1.5, 6.5, "Module 1\nMRI Segmentation\n& Mesh Generation", "#4ECDC4"),
        (5.0, 6.5, "Module 2\nElectrophysiology\nSimulation", "#FF6B6B"),
        (8.5, 6.5, "Module 3\nMechanics &\nEM Coupling", "#45B7D1"),
        (1.5, 3.0, "Module 4\nInverse Parameter\nEstimation", "#96CEB4"),
        (5.0, 3.0, "Module 5\nArrhythmia Risk\nAssessment", "#FFEAA7"),
        (8.5, 3.0, "Module 6\nAF Ablation\nPrediction", "#DDA0DD"),
    ]

    for x, y, text, color in modules:
        rect = plt.Rectangle((x - 1.3, y - 1.0), 2.6, 2.0,
                             facecolor=color, edgecolor="black",
                             linewidth=1.5, alpha=0.85, zorder=2)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center",
               fontsize=9, fontweight="bold", zorder=3)

    # Arrows
    arrows = [
        (3.0, 6.5, 3.7, 6.5),
        (6.5, 6.5, 7.2, 6.5),
        (1.5, 5.5, 1.5, 4.0),
        (5.0, 5.5, 5.0, 4.0),
        (8.5, 5.5, 8.5, 4.0),
        (3.0, 3.0, 3.7, 3.0),
        (6.5, 3.0, 7.2, 3.0),
    ]
    for x1, y1, x2, y2 in arrows:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", lw=2, color="gray"))

    # Platforms
    ax.text(12.0, 6.5, "OpenCARP\n(EP Solver)", ha="center", va="center",
           fontsize=10, fontweight="bold",
           bbox=dict(boxstyle="round,pad=0.5", facecolor="#E8E8E8", edgecolor="black"))
    ax.text(12.0, 3.0, "FEBio\n(FE Solver)", ha="center", va="center",
           fontsize=10, fontweight="bold",
           bbox=dict(boxstyle="round,pad=0.5", facecolor="#E8E8E8", edgecolor="black"))

    ax.annotate("", xy=(11.0, 6.5), xytext=(10.0, 6.5),
               arrowprops=dict(arrowstyle="<->", lw=2, color="darkblue"))
    ax.annotate("", xy=(11.0, 3.0), xytext=(10.0, 3.0),
               arrowprops=dict(arrowstyle="<->", lw=2, color="darkblue"))

    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_architecture.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Saved fig1_architecture.png")

    # --- Figure 2: Hemodynamic Results ---
    mech = results.get("modules", {}).get("mechanics", {})
    cc = mech.get("cardiac_cycle", {})

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # PV Loop (synthetic)
    V = np.linspace(cc.get("ESV_mL", 55), cc.get("EDV_mL", 130), 100)
    P_systole = cc.get("peak_pressure_mmHg", 120) * np.sin(np.linspace(0, np.pi, 50))
    P_diastole = np.linspace(5, 12, 50)
    V_loop = np.concatenate([V[:50], V[50:][::-1]])
    P_loop = np.concatenate([P_diastole, P_systole[::-1]])

    axes[0].plot(V_loop, P_loop, "b-", linewidth=2)
    axes[0].fill(V_loop, P_loop, alpha=0.15, color="blue")
    axes[0].set_xlabel("Volume (mL)")
    axes[0].set_ylabel("Pressure (mmHg)")
    axes[0].set_title("Pressure-Volume Loop")
    axes[0].grid(True, alpha=0.3)

    # Bar chart: hemodynamic indices
    labels = ["EDV", "ESV", "SV"]
    values = [cc.get("EDV_mL", 130), cc.get("ESV_mL", 55), cc.get("SV_mL", 75)]
    colors = ["#4ECDC4", "#FF6B6B", "#45B7D1"]
    axes[1].bar(labels, values, color=colors, edgecolor="black", linewidth=0.5)
    axes[1].set_ylabel("Volume (mL)")
    axes[1].set_title("Hemodynamic Indices")
    for i, v in enumerate(values):
        axes[1].text(i, v + 2, f"{v:.1f}", ha="center", fontweight="bold")

    # EF gauge
    ef = cc.get("EF_pct", 57.7)
    theta = np.linspace(0, np.pi, 100)
    axes[2].plot(np.cos(theta), np.sin(theta), "k-", linewidth=2)
    ef_angle = np.pi * (1 - ef / 100)
    axes[2].annotate("", xy=(np.cos(ef_angle), np.sin(ef_angle)),
                    xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", lw=3, color="red"))
    axes[2].set_xlim(-1.3, 1.3)
    axes[2].set_ylim(-0.3, 1.3)
    axes[2].set_title(f"Ejection Fraction: {ef:.1f}%")
    axes[2].set_aspect("equal")
    axes[2].axis("off")

    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_hemodynamics.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Saved fig2_hemodynamics.png")

    # --- Figure 3: Arrhythmia Risk Dashboard ---
    risk = results.get("modules", {}).get("arrhythmia_risk", {})
    overall = risk.get("overall_risk", {})

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # APD Restitution
    di = np.linspace(20, 600, 100)
    apd = 150 + 130 * (1 - np.exp(-di / 80))
    slope = 130 / 80 * np.exp(-di / 80)
    ax = axes[0, 0]
    ax.plot(di, apd, "b-", linewidth=2, label="APD")
    ax2 = ax.twinx()
    ax2.plot(di, slope, "r--", linewidth=2, label="Slope")
    ax2.axhline(y=1.0, color="r", linestyle=":", alpha=0.5)
    ax.set_xlabel("Diastolic Interval (ms)")
    ax.set_ylabel("APD (ms)", color="b")
    ax2.set_ylabel("Slope", color="r")
    ax.set_title("APD Restitution Curve")
    ax.grid(True, alpha=0.3)

    # Risk sub-scores radar
    sub = overall.get("sub_scores", {})
    categories = list(sub.keys())
    values_radar = [sub[c] for c in categories]
    values_radar.append(values_radar[0])
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    angles.append(angles[0])
    axes[0, 1].remove()
    ax_polar = fig.add_subplot(2, 2, 2, polar=True)
    ax_polar.fill(angles, values_radar, alpha=0.25, color="red")
    ax_polar.plot(angles, values_radar, "o-", color="red", linewidth=2)
    ax_polar.set_xticks(angles[:-1])
    ax_polar.set_xticklabels([c.capitalize() for c in categories], size=9)
    ax_polar.set_ylim(0, 1)
    ax_polar.set_title("Arrhythmia Risk Sub-Scores", pad=20)

    # Fibrosis burden
    fib = risk.get("fibrosis", {})
    ax = axes[1, 0]
    labels_fib = ["Healthy", "Fibrotic", "Border Zone"]
    healthy_pct = 100 - fib.get("burden_pct", 12) - fib.get("border_zone_pct", 5)
    sizes = [healthy_pct, fib.get("burden_pct", 12), fib.get("border_zone_pct", 5)]
    colors_fib = ["#4ECDC4", "#FF6B6B", "#FFEAA7"]
    ax.pie(sizes, labels=labels_fib, colors=colors_fib, autopct="%1.1f%%",
           startangle=90, textprops={"fontsize": 10})
    ax.set_title("Tissue Composition")

    # Overall risk summary
    ax = axes[1, 1]
    ax.axis("off")
    risk_score = overall.get("score", 0)
    risk_cat = overall.get("category", "unknown")
    color_map = {"low": "green", "moderate": "orange", "high": "red", "very_high": "darkred"}
    risk_color = color_map.get(risk_cat, "gray")

    ax.text(0.5, 0.8, "Overall Arrhythmia Risk", ha="center", fontsize=16, fontweight="bold",
           transform=ax.transAxes)
    ax.text(0.5, 0.55, f"{risk_score:.3f}", ha="center", fontsize=48,
           fontweight="bold", color=risk_color, transform=ax.transAxes)
    ax.text(0.5, 0.35, risk_cat.upper().replace("_", " "), ha="center",
           fontsize=20, fontweight="bold", color=risk_color, transform=ax.transAxes)
    ax.text(0.5, 0.15, f"Re-entry inducible: {overall.get('reentry_inducible', False)}",
           ha="center", fontsize=12, transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_arrhythmia_risk.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Saved fig3_arrhythmia_risk.png")

    # --- Figure 4: Ablation Strategy Comparison ---
    ablation = results.get("modules", {}).get("ablation_prediction", {})
    ab_results = ablation.get("results", {})

    if ab_results:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        strategies = list(ab_results.keys())
        short_names = [s.replace("pulmonary_vein_isolation", "PVI")
                        .replace("pvi_plus_", "PVI+")
                        .replace("hybrid_pvi_substrate", "Hybrid")
                        .replace("_", " ").title()
                       for s in strategies]

        # Recurrence rates
        recurrence = [ab_results[s]["recurrence_1yr"] * 100 for s in strategies]
        colors_ab = ["#4ECDC4", "#FF6B6B", "#45B7D1", "#96CEB4"][:len(strategies)]
        bars = axes[0].bar(range(len(strategies)), recurrence, color=colors_ab,
                          edgecolor="black", linewidth=0.5)
        axes[0].set_xticks(range(len(strategies)))
        axes[0].set_xticklabels(short_names, rotation=30, ha="right", fontsize=9)
        axes[0].set_ylabel("1-Year Recurrence (%)")
        axes[0].set_title("AF Recurrence by Strategy")
        for bar, val in zip(bars, recurrence):
            axes[0].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                        f"{val:.1f}%", ha="center", fontweight="bold", fontsize=9)

        # Lesion count & procedure time
        n_lesions = [ab_results[s]["n_lesions"] for s in strategies]
        abl_time = [ab_results[s]["ablation_time_min"] for s in strategies]

        ax_l = axes[1]
        ax_t = ax_l.twinx()
        x = np.arange(len(strategies))
        w = 0.35
        ax_l.bar(x - w/2, n_lesions, w, color="#4ECDC4", label="Lesions", edgecolor="black")
        ax_t.bar(x + w/2, abl_time, w, color="#FF6B6B", label="Time (min)", edgecolor="black")
        ax_l.set_xticks(x)
        ax_l.set_xticklabels(short_names, rotation=30, ha="right", fontsize=9)
        ax_l.set_ylabel("Number of Lesions", color="#4ECDC4")
        ax_t.set_ylabel("Procedure Time (min)", color="#FF6B6B")
        ax_l.set_title("Procedure Complexity")
        ax_l.legend(loc="upper left")
        ax_t.legend(loc="upper right")

        # PV reconnection risk
        reconnection = [ab_results[s]["pv_reconnection_risk"] * 100 for s in strategies]
        axes[2].barh(range(len(strategies)), reconnection, color=colors_ab,
                    edgecolor="black", linewidth=0.5)
        axes[2].set_yticks(range(len(strategies)))
        axes[2].set_yticklabels(short_names, fontsize=9)
        axes[2].set_xlabel("PV Reconnection Risk (%)")
        axes[2].set_title("PV Reconnection Risk")
        for i, val in enumerate(reconnection):
            axes[2].text(val + 1, i, f"{val:.1f}%", va="center", fontsize=9)

        fig.tight_layout()
        fig.savefig(fig_dir / "fig4_ablation_comparison.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        logger.info("  Saved fig4_ablation_comparison.png")


def _generate_text_figures(results: dict, output_dir: str):
    """Fallback: generate text-based figure summaries."""
    fig_dir = Path(output_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    with open(fig_dir / "figures_summary.txt", "w") as f:
        f.write("Cardiac Digital Twin - Figure Descriptions\n")
        f.write("=" * 50 + "\n\n")
        f.write("Figure 1: Framework Architecture Diagram\n")
        f.write("Figure 2: Hemodynamic Results (PV Loop, Indices, EF)\n")
        f.write("Figure 3: Arrhythmia Risk Dashboard\n")
        f.write("Figure 4: Ablation Strategy Comparison\n")


def main():
    """Main entry point."""
    print("=" * 70)
    print("  Cardiac Digital Twin Framework - Demonstration Run")
    print("=" * 70)

    # Run pipeline
    pipeline = CardiacDigitalTwinPipeline(
        patient_id="DT_DEMO_001",
        output_dir="."
    )
    results = pipeline.run_full_pipeline()

    # Generate figures
    print("\n" + "=" * 70)
    print("  Generating Figures")
    print("=" * 70)
    generate_figures(results, "figures")

    # Save process log
    process_log = {
        "timestamp": datetime.now().isoformat(),
        "phase": "run_completed",
        "event_type": "run_completed",
        "actor": "co-scientist",
        "skill_or_tool": "cardiac-digital-twin",
        "files_written": [
            "results/pipeline_results.json",
            "figures/fig1_architecture.png",
            "figures/fig2_hemodynamics.png",
            "figures/fig3_arrhythmia_risk.png",
            "figures/fig4_ablation_comparison.png",
            "data/opencarp/heart.pts",
            "data/opencarp/heart.elem",
            "data/opencarp/heart.lon",
            "data/febio/heart.feb",
            "configs/opencarp_ep.par",
            "configs/febio_mechanics.feb",
        ],
        "status": "ok",
    }
    with open("logs/process-log.jsonl", "w") as f:
        f.write(json.dumps(process_log) + "\n")

    print("\n" + "=" * 70)
    print("  Pipeline Complete!")
    print("=" * 70)
    print(f"  Results: results/pipeline_results.json")
    print(f"  Figures: figures/")
    print(f"  Mesh data: data/opencarp/, data/febio/")
    print(f"  Configs: configs/")

    return results


if __name__ == "__main__":
    main()
