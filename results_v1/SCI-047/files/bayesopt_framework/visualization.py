"""
Visualization module for Bayesian Optimization framework.
Generates publication-quality figures for all experiments.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]


def plot_kernel_comparison(kernel_results, filename="kernel_comparison.png"):
    """Bar chart of kernel NLPD scores."""
    fig, ax = plt.subplots(figsize=(8, 5))

    names = [r["kernel"] for r in kernel_results]
    means = [r["mean_nlpd"] for r in kernel_results]
    stds = [r["std_nlpd"] for r in kernel_results]

    bars = ax.bar(names, means, yerr=stds, capsize=5, color=COLORS[:len(names)],
                  edgecolor="black", linewidth=0.5, alpha=0.85)

    ax.set_ylabel("Mean NLPD (lower is better)")
    ax.set_title("Kernel Comparison via Cross-Validation")
    ax.set_xlabel("Kernel")

    # Highlight best
    best_idx = np.argmin(means)
    bars[best_idx].set_edgecolor("red")
    bars[best_idx].set_linewidth(2.5)

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_lengthscales(lengthscales, param_names=None, filename="lengthscales.png"):
    """Bar chart of ARD lengthscales (feature importance proxy)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    dim = len(lengthscales)

    if param_names is None:
        param_names = [f"x{i+1}" for i in range(dim)]

    inv_ls = 1.0 / np.array(lengthscales)
    inv_ls_norm = inv_ls / inv_ls.sum()

    ax.barh(param_names, inv_ls_norm, color=COLORS[0], edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Relative Importance (1/lengthscale, normalized)")
    ax.set_title("Feature Importance from ARD Lengthscales")
    ax.invert_yaxis()

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_acquisition_comparison(acq_results, filename="acquisition_comparison.png"):
    """Convergence curves for different acquisition functions."""
    fig, ax = plt.subplots(figsize=(9, 5))

    for i, (name, data) in enumerate(acq_results.items()):
        mean_traj = np.array(data["mean_trajectory"])
        std_traj = np.array(data["std_trajectory"])
        iters = np.arange(len(mean_traj))

        ax.plot(iters, mean_traj, color=COLORS[i], label=name, linewidth=2)
        ax.fill_between(iters, mean_traj - std_traj, mean_traj + std_traj,
                         color=COLORS[i], alpha=0.15)

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best Observed Value")
    ax.set_title("Acquisition Function Comparison (Hartmann-6)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_batch_comparison(batch_results, filename="batch_comparison.png"):
    """Compare batch methods: best value vs. total evaluations."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    for i, (method, data) in enumerate(batch_results.items()):
        n_evals = data["n_evals"]
        best_vals = data["best_values"]
        ax1.plot(n_evals, best_vals, color=COLORS[i], label=method,
                 linewidth=2, marker="o", markersize=4)

    ax1.set_xlabel("Total Evaluations")
    ax1.set_ylabel("Best Observed Value")
    ax1.set_title("Batch BO: Convergence by Evaluations")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Bar chart of final values and timing
    methods = list(batch_results.keys())
    finals = [batch_results[m]["final_best"] for m in methods]
    times = [batch_results[m]["total_time"] for m in methods]

    x = np.arange(len(methods))
    ax2.bar(x - 0.2, finals, 0.35, color=COLORS[0], label="Final Best Value")
    ax2_twin = ax2.twinx()
    ax2_twin.bar(x + 0.2, times, 0.35, color=COLORS[1], label="Total Time (s)")

    ax2.set_xticks(x)
    ax2.set_xticklabels(methods)
    ax2.set_ylabel("Best Value")
    ax2_twin.set_ylabel("Time (s)")
    ax2.set_title("Batch Methods: Quality vs. Cost")
    ax2.legend(loc="upper left")
    ax2_twin.legend(loc="upper right")

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_pareto_front(pareto_Y, hv_history, filename="pareto_front.png"):
    """Plot Pareto front and hypervolume convergence."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    pareto_arr = np.array(pareto_Y)
    ax1.scatter(pareto_arr[:, 0], pareto_arr[:, 1], c=COLORS[0], s=60,
                edgecolors="black", linewidth=0.5, zorder=5)

    sorted_idx = np.argsort(pareto_arr[:, 0])
    ax1.plot(pareto_arr[sorted_idx, 0], pareto_arr[sorted_idx, 1],
             color=COLORS[0], alpha=0.5, linewidth=1.5)

    ax1.set_xlabel("Objective 1 (Yield %)")
    ax1.set_ylabel("Objective 2 (Selectivity %)")
    ax1.set_title("Pareto Front")
    ax1.grid(True, alpha=0.3)

    ax2.plot(range(len(hv_history)), hv_history, color=COLORS[2], linewidth=2)
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Hypervolume")
    ax2.set_title("Hypervolume Convergence")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_highdim_comparison(hd_results, filename="highdim_comparison.png"):
    """Compare high-dimensional methods."""
    fig, ax = plt.subplots(figsize=(9, 5))

    for i, (method, data) in enumerate(hd_results.items()):
        bv = data["best_values"]
        ax.plot(range(len(bv)), bv, color=COLORS[i], label=method, linewidth=2)

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best Observed Value")
    ax.set_title("High-Dimensional BO Comparison (D=25, d_eff=6)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_chemical_convergence(convergence, filename="chemical_convergence.png"):
    """Plot convergence for chemical optimization."""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(range(len(convergence)), convergence, color=COLORS[0], linewidth=2)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best Yield (%)")
    ax.set_title("Chemical Reaction Yield Optimization Convergence")
    ax.grid(True, alpha=0.3)
    ax.axhline(y=convergence[-1], color=COLORS[3], linestyle="--", alpha=0.5,
               label=f"Final: {convergence[-1]:.1f}%")
    ax.legend()

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_chemical_pareto(pareto_front, filename="chemical_pareto.png"):
    """Plot Pareto front for chemical reaction optimization."""
    fig, ax = plt.subplots(figsize=(8, 5))

    yields = [p["yield"] for p in pareto_front]
    sels = [p["selectivity"] for p in pareto_front]

    ax.scatter(yields, sels, c=COLORS[0], s=80, edgecolors="black",
               linewidth=0.5, zorder=5)

    sorted_idx = np.argsort(yields)
    ax.plot([yields[i] for i in sorted_idx], [sels[i] for i in sorted_idx],
            color=COLORS[0], alpha=0.4, linewidth=1.5)

    ax.set_xlabel("Yield (%)")
    ax.set_ylabel("Selectivity (%)")
    ax.set_title("Chemical Reaction: Yield-Selectivity Pareto Front")
    ax.grid(True, alpha=0.3)

    # Annotate extreme points
    max_yield_idx = np.argmax(yields)
    max_sel_idx = np.argmax(sels)
    ax.annotate(f"Max Yield: {yields[max_yield_idx]:.1f}%",
                xy=(yields[max_yield_idx], sels[max_yield_idx]),
                xytext=(5, 10), textcoords="offset points", fontsize=9)
    ax.annotate(f"Max Selectivity: {sels[max_sel_idx]:.1f}%",
                xy=(yields[max_sel_idx], sels[max_sel_idx]),
                xytext=(5, -15), textcoords="offset points", fontsize=9)

    plt.tight_layout()
    path = os.path.join(FIGDIR, filename)
    fig.savefig(path)
    plt.close(fig)
    return path
