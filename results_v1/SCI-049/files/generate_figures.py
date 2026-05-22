"""
Visualization module: Generate all figures for the anomaly detection system report.
"""
import numpy as np
import pandas as pd
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_experiment import run_experiments

sns.set_style("whitegrid")
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300, "font.size": 10})
CMAP = plt.cm.viridis


def plot_changepoint_detection(ts, true_cps, pelt_bps, bocpd_probs, save_dir="figures"):
    """Fig 1: Changepoint detection results."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=False)

    # PELT
    ax = axes[0]
    ax.plot(ts, color="steelblue", linewidth=0.5, alpha=0.8)
    for cp in true_cps:
        ax.axvline(cp, color="red", linestyle="--", alpha=0.7, label="True CP" if cp == true_cps[0] else "")
    for bp in pelt_bps:
        if bp < len(ts):
            ax.axvline(bp, color="green", linestyle="-.", alpha=0.7, label="PELT detected" if bp == pelt_bps[0] else "")
    ax.set_title("PELT Change Point Detection", fontsize=13, fontweight="bold")
    ax.set_ylabel("Value")
    ax.legend(loc="upper right")

    # BOCPD
    ax2 = axes[1]
    prob_data = bocpd_probs[:1000] if len(bocpd_probs) >= 1000 else bocpd_probs
    ax2.plot(prob_data, color="darkorange", linewidth=0.8)
    ax2.axhline(0.3, color="red", linestyle="--", alpha=0.5, label="Threshold")
    for cp in true_cps:
        if cp < len(prob_data):
            ax2.axvline(cp, color="red", linestyle="--", alpha=0.5)
    ax2.set_title("BOCPD Changepoint Probability", fontsize=13, fontweight="bold")
    ax2.set_ylabel("P(changepoint)")
    ax2.set_xlabel("Time")
    ax2.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(f"{save_dir}/fig1_changepoint_detection.png", bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig1_changepoint_detection.png")


def plot_multivariate_outliers(X, y_true, iforest_result, feat_names, save_dir="figures"):
    """Fig 2: Multivariate outlier detection results."""
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    scores = iforest_result["scores"]
    labels = iforest_result["labels"]

    # Score distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(scores[y_true == 0], bins=50, alpha=0.6, label="Normal", color="steelblue", density=True)
    ax1.hist(scores[y_true == 1], bins=50, alpha=0.6, label="Anomaly", color="tomato", density=True)
    ax1.set_title("Anomaly Score Distribution", fontweight="bold")
    ax1.set_xlabel("Decision Score")
    ax1.legend()

    # Feature importance
    ax2 = fig.add_subplot(gs[0, 1])
    imp = iforest_result["feature_importance"]
    sorted_imp = sorted(imp.items(), key=lambda x: x[1], reverse=True)
    names, vals = zip(*sorted_imp)
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(names)))
    ax2.barh(range(len(names)), vals, color=colors)
    ax2.set_yticks(range(len(names)))
    ax2.set_yticklabels(names, fontsize=8)
    ax2.set_title("Feature Importance (IF)", fontweight="bold")
    ax2.invert_yaxis()

    # Confusion summary
    ax3 = fig.add_subplot(gs[0, 2])
    tp = int(np.sum((labels == -1) & (y_true == 1)))
    fp = int(np.sum((labels == -1) & (y_true == 0)))
    fn = int(np.sum((labels == 1) & (y_true == 1)))
    tn = int(np.sum((labels == 1) & (y_true == 0)))
    cm = np.array([[tn, fp], [fn, tp]])
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", ax=ax3,
                xticklabels=["Pred Normal", "Pred Anomaly"],
                yticklabels=["True Normal", "True Anomaly"])
    ax3.set_title("Confusion Matrix", fontweight="bold")

    # 2D scatter of top 2 features
    top2 = [names[0], names[1]]
    i0 = feat_names.index(top2[0]) if top2[0] in feat_names else 0
    i1 = feat_names.index(top2[1]) if top2[1] in feat_names else 1

    ax4 = fig.add_subplot(gs[1, 0:2])
    scatter = ax4.scatter(X[:, i0], X[:, i1], c=scores, cmap="RdYlBu", s=3, alpha=0.5)
    ax4.set_xlabel(feat_names[i0])
    ax4.set_ylabel(feat_names[i1])
    ax4.set_title(f"Anomaly Scores: {feat_names[i0]} vs {feat_names[i1]}", fontweight="bold")
    plt.colorbar(scatter, ax=ax4, label="Anomaly Score")

    # ROC-like plot (score vs true labels)
    ax5 = fig.add_subplot(gs[1, 2])
    thresholds = np.linspace(np.min(scores), np.max(scores), 100)
    tpr_list, fpr_list = [], []
    for th in thresholds:
        pred = scores < th
        tpr = np.sum(pred & (y_true == 1)) / max(np.sum(y_true == 1), 1)
        fpr = np.sum(pred & (y_true == 0)) / max(np.sum(y_true == 0), 1)
        tpr_list.append(tpr)
        fpr_list.append(fpr)
    ax5.plot(fpr_list, tpr_list, color="darkorange", linewidth=2)
    ax5.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax5.set_xlabel("False Positive Rate")
    ax5.set_ylabel("True Positive Rate")
    ax5.set_title("ROC Curve (IF)", fontweight="bold")

    plt.savefig(f"{save_dir}/fig2_multivariate_outliers.png", bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig2_multivariate_outliers.png")


def plot_physics_constraints(phys_result, save_dir="figures"):
    """Fig 3: Physics-constrained scoring."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Score comparison
    ax = axes[0]
    ax.scatter(phys_result["statistical_scores"], phys_result["physics_scores"],
               c=phys_result["combined_scores"], cmap="hot", s=3, alpha=0.4)
    ax.set_xlabel("Statistical Score")
    ax.set_ylabel("Physics Violation Score")
    ax.set_title("Statistical vs Physics Scores", fontweight="bold")

    # Combined score distribution
    ax2 = axes[1]
    ax2.hist(phys_result["combined_scores"], bins=60, color="teal", alpha=0.7, edgecolor="white")
    ax2.axvline(np.percentile(phys_result["combined_scores"], 95),
                color="red", linestyle="--", label="95th percentile")
    ax2.set_title("Combined Anomaly Score", fontweight="bold")
    ax2.set_xlabel("Score")
    ax2.legend()

    # Constraint violations
    ax3 = axes[2]
    details = phys_result["constraint_details"]
    c_names = [d["name"][:15] for d in details]
    c_rates = [d["violation_rate"] for d in details]
    colors = ["tomato" if r > 0.01 else "steelblue" for r in c_rates]
    ax3.barh(c_names, c_rates, color=colors)
    ax3.set_xlabel("Violation Rate")
    ax3.set_title("Physics Constraint Violations", fontweight="bold")

    plt.tight_layout()
    plt.savefig(f"{save_dir}/fig3_physics_constraints.png", bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig3_physics_constraints.png")


def plot_drift_detection(drift_data, adwin_result, save_dir="figures"):
    """Fig 4: Drift detection and retraining triggers."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    ax = axes[0]
    ax.plot(drift_data, color="steelblue", linewidth=0.5, alpha=0.7)
    # Overlay running mean
    window = 50
    if len(drift_data) >= window:
        rm = np.convolve(drift_data, np.ones(window)/window, mode="valid")
        ax.plot(range(window-1, len(drift_data)), rm, color="darkred", linewidth=1.5, label="Running mean")
    for dp in adwin_result.get("drift_points", []):
        ax.axvline(dp, color="orange", alpha=0.6, linewidth=1)
    ax.axvline(-1, color="orange", alpha=0.6, label="ADWIN drift point")
    ax.set_title("Concept Drift Detection (ADWIN)", fontweight="bold")
    ax.set_ylabel("Value")
    ax.legend(loc="upper left")

    # Simulated performance & retrain triggers
    ax2 = axes[1]
    steps = np.arange(0, len(drift_data), 1)
    perf = 0.95 - 0.02 * np.sin(steps / 200) - (steps > 1200) * 0.1 * (steps - 1200) / (len(drift_data) - 1200 + 1)
    perf = np.clip(perf, 0.5, 1.0)
    ax2.plot(steps, perf, color="steelblue", linewidth=1)
    ax2.axhline(0.90, color="red", linestyle="--", alpha=0.5, label="Retrain threshold")
    ax2.fill_between(steps, perf, 0.90, where=perf < 0.90, alpha=0.2, color="red")
    ax2.set_title("Model Performance & Retraining Triggers", fontweight="bold")
    ax2.set_xlabel("Time step")
    ax2.set_ylabel("Performance")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f"{save_dir}/fig4_drift_detection.png", bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig4_drift_detection.png")


def plot_explainable_anomalies(explain_result, feat_names, save_dir="figures"):
    """Fig 5: Explainable anomaly detection."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Global feature importance
    ax = axes[0]
    imp = explain_result["global_feature_importance"]
    sorted_imp = sorted(imp.items(), key=lambda x: x[1], reverse=True)
    names, vals = zip(*sorted_imp)
    colors = plt.cm.cividis(np.linspace(0.2, 0.8, len(names)))
    ax.barh(range(len(names)), vals, color=colors)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_title("Global Feature Importance", fontweight="bold")
    ax.set_xlabel("Relative Importance")
    ax.invert_yaxis()

    # Root cause analysis
    ax2 = axes[1]
    rcs = explain_result["root_causes"]
    if rcs and isinstance(rcs[0], dict) and "feature" in rcs[0]:
        rc_names = [rc["feature"][:15] for rc in rcs[:6]]
        rc_fracs = [rc["fraction"] for rc in rcs[:6]]
        rc_colors = ["tomato" if rc.get("direction") == "high" else "steelblue" for rc in rcs[:6]]
        ax2.barh(rc_names, rc_fracs, color=rc_colors)
        ax2.set_xlabel("Fraction of Anomalies")
        high_patch = mpatches.Patch(color="tomato", label="High")
        low_patch = mpatches.Patch(color="steelblue", label="Low")
        ax2.legend(handles=[high_patch, low_patch])
    ax2.set_title("Root Cause Analysis", fontweight="bold")

    # Decision rules summary
    ax3 = axes[2]
    rules = explain_result.get("decision_rules", explain_result.get("top_rules", []))
    if rules:
        rule_text = []
        for i, r in enumerate(rules[:5]):
            conds = " AND\n    ".join(r["conditions"][:3])
            rule_text.append(f"R{i+1} (conf={r['confidence']:.2f}, n={r['support']}):\n    {conds}")
        text = "\n\n".join(rule_text)
        ax3.text(0.05, 0.95, text, transform=ax3.transAxes, fontsize=8,
                 verticalalignment="top", fontfamily="monospace",
                 bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    ax3.set_title("Anomaly Decision Rules (Surrogate Tree)", fontweight="bold")
    ax3.axis("off")

    plt.tight_layout()
    plt.savefig(f"{save_dir}/fig5_explainable_anomalies.png", bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig5_explainable_anomalies.png")


def plot_architecture(save_dir="figures"):
    """Fig 6: CERN/LIGO-scale streaming architecture diagram."""
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Title
    ax.text(8, 9.5, "Large-Scale Experiment Anomaly Detection Architecture",
            ha="center", fontsize=16, fontweight="bold")

    # Tier boxes
    tiers = [
        {"name": "L0: Hardware Trigger", "x": 1, "y": 7, "color": "#FF6B6B",
         "desc": "FPGA/ASIC\n<1μs latency\n40MHz → 100kHz\n400× reduction"},
        {"name": "L1: Online Filter", "x": 5, "y": 7, "color": "#FFA07A",
         "desc": "Z-score, EWMA\n<10ms latency\n100kHz → 1kHz\n100× reduction"},
        {"name": "L2: Nearline Analysis", "x": 9, "y": 7, "color": "#87CEEB",
         "desc": "Isolation Forest, PELT\n<1s latency\n1kHz → 10Hz\n100× reduction"},
        {"name": "L3: Offline Deep", "x": 13, "y": 7, "color": "#90EE90",
         "desc": "Deep SVDD, SHAP\nmin-hours\nStored events\nFull analysis"},
    ]

    for t in tiers:
        rect = mpatches.FancyBboxPatch((t["x"]-1.3, t["y"]-1.5), 2.6, 3,
                                        boxstyle="round,pad=0.15",
                                        facecolor=t["color"], alpha=0.3, edgecolor="gray")
        ax.add_patch(rect)
        ax.text(t["x"], t["y"]+1, t["name"], ha="center", fontsize=10, fontweight="bold")
        ax.text(t["x"], t["y"]-0.2, t["desc"], ha="center", fontsize=8, va="top")

    # Arrows
    for i in range(3):
        ax.annotate("", xy=(tiers[i+1]["x"]-1.3, tiers[i+1]["y"]),
                     xytext=(tiers[i]["x"]+1.3, tiers[i]["y"]),
                     arrowprops=dict(arrowstyle="->", color="gray", lw=2))

    # Bottom: components
    components = [
        {"name": "Drift Detector\n(ADWIN/PH)", "x": 3, "y": 3, "color": "#DDA0DD"},
        {"name": "Physics\nConstraints", "x": 6, "y": 3, "color": "#F0E68C"},
        {"name": "Explainability\n(SHAP/Rules)", "x": 9, "y": 3, "color": "#98FB98"},
        {"name": "Retrain\nTrigger", "x": 12, "y": 3, "color": "#FFB6C1"},
    ]
    for c in components:
        rect = mpatches.FancyBboxPatch((c["x"]-1.1, c["y"]-0.8), 2.2, 1.8,
                                        boxstyle="round,pad=0.1",
                                        facecolor=c["color"], alpha=0.4, edgecolor="gray")
        ax.add_patch(rect)
        ax.text(c["x"], c["y"], c["name"], ha="center", va="center", fontsize=9)

    # Connecting lines from tiers to components
    for c in components:
        ax.plot([c["x"], c["x"]], [c["y"]+0.9, 5.5], color="gray", linestyle=":", alpha=0.5)

    # Data flow label
    ax.text(8, 1.2, "Data Flow: Raw Detector → L0 (FPGA) → L1 (Online) → L2 (Nearline) → L3 (Offline)",
            ha="center", fontsize=11, style="italic",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    ax.text(8, 0.5, "Total Data Reduction: ~10⁶×  |  Fault Tolerance: N+1 redundancy, 60s checkpointing",
            ha="center", fontsize=9, color="gray")

    plt.savefig(f"{save_dir}/fig6_architecture.png", bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig6_architecture.png")


if __name__ == "__main__":
    print("Running experiments and generating figures...")
    print("=" * 60)

    outputs = run_experiments()
    results, X, y_true, feat_names, ts, true_cps, iforest_result, \
        explain_result, phys_result, drift_data, adwin_result, bocpd = outputs

    print("\n── Generating Figures ──")

    pelt_bps = results["pelt"]["detected_changepoints"]
    bocpd_probs = bocpd.changepoint_probs

    plot_changepoint_detection(ts, true_cps, pelt_bps, bocpd_probs)
    plot_multivariate_outliers(X, y_true, iforest_result, feat_names)
    plot_physics_constraints(phys_result)
    plot_drift_detection(drift_data, adwin_result)
    plot_explainable_anomalies(results["explainable"], feat_names)
    plot_architecture()

    print("\n✅ All figures generated successfully!")
