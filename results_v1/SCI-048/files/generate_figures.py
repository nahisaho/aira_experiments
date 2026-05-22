"""
Generate all figures for the PINN benchmark report.
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 9,
    "figure.dpi": 150,
})

FIGURES_DIR = "figures"
RESULTS_DIR = "results"


def load_json(name):
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None


def plot_fourier_features():
    data = load_json("fourier_features")
    if not data:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # L2 error comparison
    names = list(data.keys())
    errors = [data[n]["l2_relative_error"] for n in names]
    short_names = [n.replace("Multi-scale ", "MS\n") for n in names]
    colors = ["#1f77b4", "#2ca02c", "#d62728"]

    axes[0].bar(range(len(names)), errors, color=colors[:len(names)])
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(short_names, fontsize=9)
    axes[0].set_ylabel("Relative L2 Error")
    axes[0].set_title("Multi-scale Fourier Feature: L2 Error Comparison")
    axes[0].set_yscale("log")

    # Loss convergence
    for i, name in enumerate(names):
        hist = data[name]["loss_history"]
        epochs = [h["epoch"] for h in hist]
        losses = [h["loss"] for h in hist]
        axes[1].semilogy(epochs, losses, label=name.split("(")[0].strip(), color=colors[i])

    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Total Loss")
    axes[1].set_title("Training Loss Convergence")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig1_fourier_features.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIGURES_DIR, "fig1_fourier_features.svg"), bbox_inches="tight")
    plt.close()
    print("  Saved fig1_fourier_features")


def plot_inverse_problem():
    data = load_json("inverse_problem")
    if not data:
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # D convergence
    hist = data["history"]
    axes[0].plot(hist["epoch"], hist["D_est"], "b-", linewidth=2, label="Estimated D")
    axes[0].axhline(y=data["D_true"], color="r", linestyle="--", linewidth=2, label=f"True D = {data['D_true']}")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("D")
    axes[0].set_title("Parameter Estimation Convergence")
    axes[0].legend()

    # Prediction with uncertainty
    x = np.array(data["test_x"]).flatten()
    mean = np.array(data["test_mean"]).flatten()
    std = np.array(data["test_std"]).flatten()
    exact = np.array(data["test_exact"]).flatten()

    axes[1].plot(x, exact, "r-", linewidth=2, label="Exact")
    axes[1].plot(x, mean, "b--", linewidth=2, label="PINN Mean")
    axes[1].fill_between(x, mean - 2 * std, mean + 2 * std, alpha=0.3, color="blue", label="±2σ (MC Dropout)")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("u(x, t=0.25)")
    axes[1].set_title("Prediction with Uncertainty")
    axes[1].legend()

    # Loss curves
    axes[2].semilogy(hist["epoch"], hist["data_loss"], label="Data Loss")
    axes[2].semilogy(hist["epoch"], hist["pde_loss"], label="PDE Loss")
    axes[2].semilogy(hist["epoch"], hist["loss"], label="Total Loss")
    axes[2].set_xlabel("Epoch")
    axes[2].set_ylabel("Loss")
    axes[2].set_title("Loss Components")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig2_inverse_problem.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIGURES_DIR, "fig2_inverse_problem.svg"), bbox_inches="tight")
    plt.close()
    print("  Saved fig2_inverse_problem")


def plot_causal_training():
    data = load_json("causal_training")
    if not data:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Error by time slice
    for mode in ["standard", "causal"]:
        if mode not in data:
            continue
        errors = data[mode]["errors_by_time"]
        t_vals = [e["t"] for e in errors]
        l2_vals = [e["l2_error"] for e in errors]
        marker = "o-" if mode == "standard" else "s-"
        axes[0].plot(t_vals, l2_vals, marker, linewidth=2, markersize=8, label=f"{mode.capitalize()}")

    axes[0].set_xlabel("Time t")
    axes[0].set_ylabel("Relative L2 Error")
    axes[0].set_title("Error vs Time (Causal vs Standard)")
    axes[0].legend()
    axes[0].set_yscale("log")

    # Loss convergence
    for mode in ["standard", "causal"]:
        if mode not in data:
            continue
        hist = data[mode]["loss_history"]
        epochs = [h["epoch"] for h in hist]
        losses = [h["loss"] for h in hist]
        axes[1].semilogy(epochs, losses, linewidth=2, label=f"{mode.capitalize()}")

    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Total Loss")
    axes[1].set_title("Training Loss Convergence")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig3_causal_training.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIGURES_DIR, "fig3_causal_training.svg"), bbox_inches="tight")
    plt.close()
    print("  Saved fig3_causal_training")


def plot_adaptive_collocation():
    data = load_json("adaptive_collocation")
    if not data:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = list(data.keys())
    errors = [data[n]["l2_relative_error"] for n in names]
    n_pts = [data[n]["n_collocation_points"] for n in names]

    colors = ["#1f77b4", "#2ca02c", "#d62728"]
    bars = axes[0].bar(range(len(names)), errors, color=colors[:len(names)])
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(names, fontsize=10)
    axes[0].set_ylabel("Relative L2 Error")
    axes[0].set_title("Adaptive Collocation: Error Comparison")
    axes[0].set_yscale("log")

    for i, (bar, n) in enumerate(zip(bars, n_pts)):
        axes[0].text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                     f"n={n}", ha="center", va="bottom", fontsize=9)

    for i, name in enumerate(names):
        hist = data[name]["loss_history"]
        epochs = [h["epoch"] for h in hist]
        losses = [h["loss"] for h in hist]
        axes[1].semilogy(epochs, losses, label=name, color=colors[i])

    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Total Loss")
    axes[1].set_title("Training Loss Convergence")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig4_adaptive_collocation.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIGURES_DIR, "fig4_adaptive_collocation.svg"), bbox_inches="tight")
    plt.close()
    print("  Saved fig4_adaptive_collocation")


def plot_operator_learning():
    data = load_json("operator_learning")
    if not data:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = list(data.keys())
    errors = [data[n]["l2_relative_error"] for n in names]
    params = [data[n]["params"] for n in names]

    colors = ["#1f77b4", "#ff7f0e"]
    axes[0].bar(range(len(names)), errors, color=colors[:len(names)])
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(names)
    axes[0].set_ylabel("Relative L2 Error")
    axes[0].set_title("DeepONet vs FNO: Error Comparison")

    for i, name in enumerate(names):
        hist = data[name]["loss_history"]
        epochs = [h["epoch"] for h in hist]
        losses = [h["loss"] for h in hist]
        axes[1].semilogy(epochs, losses, label=f"{name} ({params[i]:,} params)", color=colors[i])

    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("MSE Loss")
    axes[1].set_title("Training Loss Convergence")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig5_operator_learning.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIGURES_DIR, "fig5_operator_learning.svg"), bbox_inches="tight")
    plt.close()
    print("  Saved fig5_operator_learning")


def plot_navier_stokes():
    data = load_json("navier_stokes")
    if not data:
        return

    re_keys = [k for k in data.keys() if k.startswith("Re")]
    n_re = len(re_keys)

    fig, axes = plt.subplots(1, max(n_re, 2), figsize=(6 * max(n_re, 2), 5))
    if n_re == 1:
        axes = [axes, plt.subplot(122)]

    for i, key in enumerate(re_keys):
        d = data[key]
        y = np.array(d["y_coords"]).flatten()
        u = np.array(d["u_centerline_vert"]).flatten()
        axes[i].plot(u, y, "b-", linewidth=2, label="PINN")
        axes[i].set_xlabel("u velocity")
        axes[i].set_ylabel("y")
        axes[i].set_title(f"Vertical Centerline (Re={d['Re']})")
        axes[i].legend()
        axes[i].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig6_navier_stokes.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIGURES_DIR, "fig6_navier_stokes.svg"), bbox_inches="tight")
    plt.close()
    print("  Saved fig6_navier_stokes")


def plot_summary_table():
    """Create a summary comparison figure."""
    all_data = load_json("all_results")
    if not all_data:
        return

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis("off")

    headers = ["Module", "Method", "Key Metric", "Value"]
    rows = []

    if "fourier_features" in all_data and "error" not in all_data["fourier_features"]:
        d = all_data["fourier_features"]
        for name, r in d.items():
            rows.append(["1. Fourier Features", name, "L2 Error", f"{r['l2_relative_error']:.6f}"])

    if "inverse_problem" in all_data and "error" not in all_data["inverse_problem"]:
        d = all_data["inverse_problem"]
        rows.append(["2. Inverse Problem", "MC Dropout PINN", "D rel. error", f"{d['D_relative_error']:.4%}"])
        rows.append(["", "", "Mean σ", f"{d['mean_uncertainty']:.6f}"])

    if "causal_training" in all_data and "error" not in all_data["causal_training"]:
        d = all_data["causal_training"]
        if "standard" in d:
            rows.append(["3. Causal Training", "Standard", "Avg L2 Error", f"{d['standard']['avg_l2_error']:.6f}"])
        if "causal" in d:
            rows.append(["", "Causal", "Avg L2 Error", f"{d['causal']['avg_l2_error']:.6f}"])

    if "adaptive_collocation" in all_data and "error" not in all_data["adaptive_collocation"]:
        d = all_data["adaptive_collocation"]
        for name, r in d.items():
            rows.append(["4. Adaptive Colloc.", name, "L2 Error", f"{r['l2_relative_error']:.6f}"])

    if "operator_learning" in all_data and "error" not in all_data["operator_learning"]:
        d = all_data["operator_learning"]
        for name, r in d.items():
            rows.append(["5. Operator Learning", name, "L2 Error", f"{r['l2_relative_error']:.6f}"])

    if "navier_stokes" in all_data and "error" not in all_data["navier_stokes"]:
        d = all_data["navier_stokes"]
        for key, r in d.items():
            if isinstance(r, dict) and "Re" in r:
                rows.append(["6. Navier-Stokes", f"Re={r['Re']}", "Cont. MSE", f"{r['continuity_mse']:.2e}"])

    if rows:
        table = ax.table(cellText=rows, colLabels=headers, loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.0, 1.5)
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor("#4472C4")
                cell.set_text_props(color="white", fontweight="bold")
            elif row % 2 == 0:
                cell.set_facecolor("#D6E4F0")

    ax.set_title("Summary of All PINN Benchmark Results", fontsize=14, fontweight="bold", pad=20)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig7_summary_table.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved fig7_summary_table")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    print("Generating figures...")
    plot_fourier_features()
    plot_inverse_problem()
    plot_causal_training()
    plot_adaptive_collocation()
    plot_operator_learning()
    plot_navier_stokes()
    plot_summary_table()
    print("All figures generated.")


if __name__ == "__main__":
    main()
