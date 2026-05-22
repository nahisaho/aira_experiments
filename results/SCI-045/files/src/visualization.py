"""
Visualization module for epigenetic clock results.
All figure text (axes, legends, titles) in English.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

FIG_DIR = Path("figures")
FIG_DIR.mkdir(exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300, "font.size": 10})


def plot_prediction_scatter(y_true, y_pred, label, filename):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_true, y_pred, alpha=0.4, s=15, c="steelblue")
    mn, mx = min(y_true.min(), y_pred.min()) - 5, max(y_true.max(), y_pred.max()) + 5
    ax.plot([mn, mx], [mn, mx], "k--", lw=1, label="y=x")
    from sklearn.metrics import mean_absolute_error, r2_score
    from scipy.stats import pearsonr
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    r, _ = pearsonr(y_true, y_pred)
    ax.set_xlabel("Chronological Age (years)")
    ax.set_ylabel("Predicted DNAm Age (years)")
    ax.set_title(f"{label}\nMAE={mae:.2f}, R²={r2:.3f}, r={r:.3f}")
    ax.legend()
    ax.set_xlim(mn, mx)
    ax.set_ylim(mn, mx)
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")


def plot_model_comparison(metrics_list, filename="model_comparison.png"):
    df = pd.DataFrame(metrics_list)
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    colors = sns.color_palette("viridis", len(df))
    for i, (metric, ax) in enumerate(zip(["MAE", "R2", "Pearson_r"], axes)):
        bars = ax.bar(df["model"], df[metric], color=colors)
        ax.set_ylabel(metric)
        ax.set_title(f"Model Comparison: {metric}")
        ax.tick_params(axis="x", rotation=45)
        for bar, val in zip(bars, df[metric]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f"{val:.3f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")


def plot_tissue_performance(tissue_results, filename="tissue_performance.png"):
    tissues = list(tissue_results.keys())
    maes = [tissue_results[t]["MAE"] for t in tissues]
    r2s = [tissue_results[t]["R2"] for t in tissues]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    colors = sns.color_palette("cividis", len(tissues))
    ax1.bar(tissues, maes, color=colors)
    ax1.set_ylabel("MAE (years)")
    ax1.set_title("Tissue-Specific Clock Performance: MAE")
    ax1.tick_params(axis="x", rotation=45)
    ax2.bar(tissues, r2s, color=colors)
    ax2.set_ylabel("R²")
    ax2.set_title("Tissue-Specific Clock Performance: R²")
    ax2.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")


def plot_age_acceleration(accel_values, groups=None, filename="age_acceleration.png"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].hist(accel_values, bins=30, color="steelblue", edgecolor="white", alpha=0.8)
    axes[0].axvline(0, color="red", linestyle="--", lw=1.5)
    axes[0].set_xlabel("Age Acceleration (years)")
    axes[0].set_ylabel("Frequency")
    axes[0].set_title("Distribution of Age Acceleration")
    if groups is not None:
        unique_groups = sorted(groups.unique())
        group_data = [accel_values[groups == g] for g in unique_groups]
        bp = axes[1].boxplot(group_data, labels=unique_groups, patch_artist=True)
        colors = sns.color_palette("viridis", len(unique_groups))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
        axes[1].axhline(0, color="red", linestyle="--", lw=1)
        axes[1].set_ylabel("Age Acceleration (years)")
        axes[1].set_title("Age Acceleration by Group")
        axes[1].tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")


def plot_intervention_effects(intv_results, filename="intervention_effects.png"):
    interventions = [k for k in intv_results.keys() if k != "control_mean_accel"]
    deltas = [intv_results[k]["delta_vs_control"] for k in interventions]
    pvals = [intv_results[k]["p_value"] for k in interventions]
    cohens = [intv_results[k]["cohens_d"] for k in interventions]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    colors = ["#2ecc71" if p < 0.05 else "#e74c3c" for p in pvals]
    bars = ax1.bar(interventions, deltas, color=colors, edgecolor="white")
    ax1.axhline(0, color="gray", linestyle="--")
    ax1.set_ylabel("Δ Age Acceleration vs Control (years)")
    ax1.set_title("Intervention Effects on Biological Age")
    ax1.tick_params(axis="x", rotation=45)
    for bar, p in zip(bars, pvals):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                 sig, ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax2.bar(interventions, cohens, color=sns.color_palette("viridis", len(interventions)))
    ax2.axhline(0.2, color="gray", linestyle=":", label="Small effect (0.2)")
    ax2.axhline(0.5, color="gray", linestyle="--", label="Medium effect (0.5)")
    ax2.axhline(0.8, color="gray", linestyle="-", label="Large effect (0.8)")
    ax2.set_ylabel("Cohen's d")
    ax2.set_title("Effect Size of Interventions")
    ax2.tick_params(axis="x", rotation=45)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")


def plot_training_history(history, filename="training_history.png"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(history["train_loss"], label="Training Loss", color="steelblue")
    if history.get("val_loss"):
        ax.plot(history["val_loss"], label="Validation Loss", color="coral")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss (Huber)")
    ax.set_title("Deep Clock Training History")
    ax.legend()
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")


def plot_longevity_comparison(accel_long, accel_norm, filename="longevity_validation.png"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(accel_norm, bins=25, alpha=0.6, label="Normal Aging", color="steelblue", edgecolor="white")
    ax.hist(accel_long, bins=25, alpha=0.6, label="Longevity Cohort", color="coral", edgecolor="white")
    ax.axvline(0, color="black", linestyle="--", lw=1)
    ax.set_xlabel("Age Acceleration (years)")
    ax.set_ylabel("Frequency")
    ax.set_title("Age Acceleration: Longevity vs Normal Cohort")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {FIG_DIR / filename}")
