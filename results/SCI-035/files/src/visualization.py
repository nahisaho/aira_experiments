"""
Visualization Module: Publication-quality figures for the evaluation framework.
All text (labels, titles, legends) in English.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional

FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(exist_ok=True)

PALETTE = {
    "SA": "#1f77b4",
    "SQA": "#ff7f0e",
    "QAOA(p=2)": "#2ca02c",
    "Greedy": "#d62728",
    "ReverseAnnealing": "#9467bd",
    "GreedyLocalSearch": "#d62728",
    "BruteForce": "#8c564b",
}

sns.set_theme(style="whitegrid", palette="colorblind")
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "figure.dpi": 150,
})


# ------------------------------------------------------------------ #
#  Figure 1: Solver Comparison (energy distribution)                 #
# ------------------------------------------------------------------ #
def plot_solver_comparison(df: pd.DataFrame, filename: str = "fig1_solver_comparison.pdf"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # --- Bar chart: best energy ---
    df_valid = df.dropna(subset=["best_energy"])
    solvers = df_valid["solver"].tolist()
    colors = [PALETTE.get(s.split("(")[0], "#aec7e8") for s in solvers]

    axes[0].barh(solvers, df_valid["best_energy"], color=colors, edgecolor="k", linewidth=0.5)
    axes[0].set_xlabel("Best QUBO Energy")
    axes[0].set_title("Best Energy Achieved by Solver")
    axes[0].axvline(0, color="k", linewidth=0.8, linestyle="--")

    # --- Bar chart: elapsed time ---
    if "elapsed_sec" in df_valid.columns:
        times = df_valid["elapsed_sec"].fillna(0)
        axes[1].barh(solvers, times, color=colors, edgecolor="k", linewidth=0.5)
        axes[1].set_xlabel("Wall-clock Time (s)")
        axes[1].set_title("Computation Time by Solver")
        axes[1].set_xscale("log")

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


# ------------------------------------------------------------------ #
#  Figure 2: Scaling Analysis                                         #
# ------------------------------------------------------------------ #
def plot_scaling_analysis(df: pd.DataFrame, filename: str = "fig2_scaling_analysis.pdf"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    df_valid = df.dropna(subset=["elapsed_sec", "n"])
    solver_groups = df_valid.groupby("solver")

    for solver_name, group in solver_groups:
        color = PALETTE.get(solver_name, None)
        x = group.groupby("n")["elapsed_sec"].mean()
        axes[0].plot(x.index, x.values, marker="o", label=solver_name, color=color)
        x_energy = group.groupby("n")["best_energy"].mean()
        axes[1].plot(x_energy.index, x_energy.values, marker="s", label=solver_name, color=color)

    axes[0].set_xlabel("Problem Size (n variables)")
    axes[0].set_ylabel("Wall-clock Time (s)")
    axes[0].set_title("Time-to-Solution vs Problem Size")
    axes[0].legend()
    axes[0].set_yscale("log")

    axes[1].set_xlabel("Problem Size (n variables)")
    axes[1].set_ylabel("Best QUBO Energy")
    axes[1].set_title("Solution Quality vs Problem Size")
    axes[1].legend()

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


# ------------------------------------------------------------------ #
#  Figure 3: Annealing Schedule Comparison                            #
# ------------------------------------------------------------------ #
def plot_schedule_comparison(df: pd.DataFrame, filename: str = "fig3_schedule_comparison.pdf"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Energy by schedule
    x = range(len(df))
    bars = axes[0].bar(x, df["best_energy"], color=sns.color_palette("viridis", len(df)),
                       edgecolor="k", linewidth=0.5)
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(df["schedule"].tolist(), rotation=45, ha="right")
    axes[0].set_ylabel("Best QUBO Energy")
    axes[0].set_title("Best Energy by Annealing Schedule")

    # Error bars for mean ± std
    axes[1].bar(x, df["mean_energy"], yerr=df.get("std_energy", None),
                color=sns.color_palette("cividis", len(df)),
                edgecolor="k", linewidth=0.5, capsize=4)
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(df["schedule"].tolist(), rotation=45, ha="right")
    axes[1].set_ylabel("Mean Energy ± Std")
    axes[1].set_title("Mean Solution Quality by Schedule")

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


# ------------------------------------------------------------------ #
#  Figure 4: Embedding Analysis                                       #
# ------------------------------------------------------------------ #
def plot_embedding_analysis(df: pd.DataFrame, filename: str = "fig4_embedding_analysis.pdf"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for strategy, group in df.groupby("strategy"):
        axes[0].plot(group["problem_size"], group["avg_chain_length"],
                     marker="o", label=strategy)
        axes[1].plot(group["problem_size"], group["overhead"],
                     marker="s", label=strategy)

    axes[0].set_xlabel("Problem Size (logical qubits)")
    axes[0].set_ylabel("Average Chain Length")
    axes[0].set_title("Embedding Chain Length vs Problem Size")
    axes[0].legend()

    axes[1].set_xlabel("Problem Size (logical qubits)")
    axes[1].set_ylabel("Physical/Logical Qubit Ratio")
    axes[1].set_title("Embedding Overhead vs Problem Size")
    axes[1].legend()

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


# ------------------------------------------------------------------ #
#  Figure 5: VRP Route Visualization                                  #
# ------------------------------------------------------------------ #
def plot_vrp_routes(vrp_data: dict, N_key: int, filename: str = "fig5_vrp_routes.pdf"):
    if N_key not in vrp_data:
        return None
    data = vrp_data[N_key]
    coords = np.array(data["coords"])

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(coords[:, 0], coords[:, 1], s=100, zorder=5, c="steelblue", edgecolors="k")
    for i, (x, y) in enumerate(coords):
        ax.annotate(f"C{i}" if i > 0 else "Depot",
                    (x + 1.5, y + 1.5), fontsize=9, zorder=6)

    # Draw depot
    ax.scatter(coords[0, 0], coords[0, 1], s=200, marker="*", c="gold", edgecolors="k", zorder=6)

    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.set_title(f"VRP Instance: {N_key} cities, 2 vehicles")

    # Solver performance annotation
    solver_rows = data.get("solver_results", [])
    if solver_rows:
        best = min(solver_rows, key=lambda r: r.get("best_energy", float("inf")))
        ax.annotate(
            f"Best: {best.get('solver','?')} | E={best.get('best_energy',0):.2f}",
            xy=(0.02, 0.02), xycoords="axes fraction", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
        )

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


# ------------------------------------------------------------------ #
#  Figure 6: Reverse Annealing Schedule                               #
# ------------------------------------------------------------------ #
def plot_reverse_annealing_schedule(filename: str = "fig6_reverse_annealing.pdf"):
    """Illustrative plot of the reverse annealing protocol."""
    fig, ax = plt.subplots(figsize=(8, 4))

    steps_down = 20
    hold_time_steps = 5
    steps_up = 30
    s_target = 0.3

    s_values = (
        [1.0 - (1.0 - s_target) * (i + 1) / steps_down for i in range(steps_down)]
        + [s_target] * hold_time_steps
        + [s_target + (1.0 - s_target) * (i + 1) / steps_up for i in range(steps_up)]
    )

    t = list(range(len(s_values)))
    ax.plot(t, s_values, color="steelblue", linewidth=2)
    ax.axhline(s_target, color="tomato", linestyle="--", linewidth=1.2, label=f"s_target = {s_target}")
    ax.axvspan(steps_down, steps_down + hold_time_steps, alpha=0.15, color="orange", label="Hold phase")
    ax.set_xlabel("Annealing Step")
    ax.set_ylabel("Annealing Parameter s")
    ax.set_title("Reverse Annealing Schedule: s vs Time")
    ax.legend()
    ax.set_ylim(0, 1.05)

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


# ------------------------------------------------------------------ #
#  Figure 7: QUBO Coefficient Distribution                            #
# ------------------------------------------------------------------ #
def plot_qubo_distribution(Q: dict, title: str = "VRP QUBO", filename: str = "fig7_qubo_distribution.pdf"):
    vals = list(Q.values())
    linear = [v for (i, j), v in Q.items() if i == j]
    quadratic = [v for (i, j), v in Q.items() if i != j]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].hist(linear, bins=40, color="steelblue", edgecolor="k", linewidth=0.3, alpha=0.8)
    axes[0].set_xlabel("Coefficient Value")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"{title}: Linear (Diagonal) Terms")

    axes[1].hist(quadratic, bins=40, color="coral", edgecolor="k", linewidth=0.3, alpha=0.8)
    axes[1].set_xlabel("Coefficient Value")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"{title}: Quadratic (Off-diagonal) Terms")

    plt.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return str(out)


if __name__ == "__main__":
    # Test with dummy data
    import pandas as pd
    df_dummy = pd.DataFrame({
        "solver": ["SA", "SQA", "QAOA(p=2)", "Greedy"],
        "best_energy": [-12.3, -11.8, -10.2, -10.9],
        "mean_energy": [-11.0, -10.5, -9.8, -10.2],
        "elapsed_sec": [0.5, 1.2, 5.3, 0.3],
        "success_rate": [0.8, 0.7, 0.5, 0.6],
    })
    plot_solver_comparison(df_dummy)
    print("Visualization test passed.")
