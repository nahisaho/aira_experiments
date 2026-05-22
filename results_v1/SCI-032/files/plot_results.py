"""
Publication-quality figures for surface code simulation results.
All figure text in English. Uses colorblind-friendly palettes.
"""

import numpy as np
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import LogLocator, LogFormatter
from typing import List, Dict

RESULTS_DIR = "/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/results"
FIGURES_DIR = "/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# Colorblind-friendly palette (Okabe-Ito)
COLORS = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]
DISTANCE_COLORS = {3: "#56B4E9", 5: "#009E73", 7: "#E69F00", 9: "#D55E00"}
LINESTYLES = ["-", "--", "-.", ":"]


def load(filename):
    with open(f"{RESULTS_DIR}/{filename}") as f:
        return json.load(f)


def save_fig(fig, name, dpi=300):
    path_svg = f"{FIGURES_DIR}/{name}.svg"
    path_png = f"{FIGURES_DIR}/{name}.png"
    fig.savefig(path_svg, bbox_inches="tight", format="svg")
    fig.savefig(path_png, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path_png}, {path_svg}", flush=True)
    return path_png, path_svg


def plot_noise_model_comparison():
    """Figure 1: Logical error rates for different noise models."""
    data = load("noise_model_comparison.json")
    p = np.array(data["error_rates"])

    fig, ax = plt.subplots(figsize=(7, 5))

    models = [
        ("depolarizing", "Depolarizing", COLORS[0], "o", "-"),
        ("amplitude_damping", "Amplitude Damping (T1)", COLORS[1], "s", "--"),
        ("phase_damping", "Phase Damping (T2*)", COLORS[2], "^", "-."),
        ("combined", "Combined (T1+T2*+Dep)", COLORS[5], "D", ":"),
    ]

    for key, label, color, marker, ls in models:
        rates = np.array(data[key])
        ax.semilogy(p * 100, np.clip(rates, 1e-5, 1), marker=marker,
                    label=label, color=color, linestyle=ls,
                    markersize=6, linewidth=1.8)

    ax.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax.set_ylabel("Logical Error Rate per Shot", fontsize=12)
    ax.set_title(f"Noise Model Comparison (d={data['distance']}, {data['rounds']} rounds)",
                 fontsize=13)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_xlim(0, max(p) * 100 * 1.05)

    fig.tight_layout()
    return save_fig(fig, "fig1_noise_model_comparison")


def plot_threshold_analysis():
    """Figure 2: Threshold analysis - logical error rate vs physical error rate."""
    data = load("mwpm_threshold.json")
    p = np.array(data["error_rates"])
    distances = data["distances"]
    logical_rates = np.array(data["logical_rates"])
    ci_half = np.array(data["ci_half"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: Linear scale for threshold crossing
    for i, d in enumerate(distances):
        rates = logical_rates[i]
        ci = ci_half[i]
        color = DISTANCE_COLORS.get(d, COLORS[i])
        ax1.plot(p * 100, rates, "o-", label=f"d={d}", color=color,
                 linewidth=2, markersize=5)
        ax1.fill_between(p * 100,
                         np.clip(rates - ci, 0, 1),
                         np.clip(rates + ci, 0, 1),
                         alpha=0.2, color=color)

    th = data.get("threshold_estimate", 0.006)
    th_unc = data.get("threshold_uncertainty", 0.001)
    ax1.axvline(th * 100, color="black", linestyle="--", linewidth=1.5,
                label=f"p_th ≈ {th*100:.2f}%")
    ax1.axvspan((th - th_unc) * 100, (th + th_unc) * 100, alpha=0.15, color="gray")

    ax1.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax1.set_ylabel("Logical Error Rate per Shot", fontsize=12)
    ax1.set_title("Surface Code Threshold (MWPM, Circuit-Level)", fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(logical_rates.max(), 0.05))

    # Right: Log scale
    for i, d in enumerate(distances):
        rates = logical_rates[i]
        color = DISTANCE_COLORS.get(d, COLORS[i])
        valid = np.array(rates) > 1e-6
        if valid.any():
            ax2.semilogy(p[valid] * 100, np.array(rates)[valid], "o-",
                         label=f"d={d}", color=color, linewidth=2, markersize=5)

    ax2.axvline(th * 100, color="black", linestyle="--", linewidth=1.5,
                label=f"p_th ≈ {th*100:.2f}%")
    ax2.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax2.set_ylabel("Logical Error Rate per Shot (log)", fontsize=12)
    ax2.set_title("Threshold — Log Scale", fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return save_fig(fig, "fig2_threshold_analysis")


def plot_decoder_comparison():
    """Figure 3: MWPM vs Union-Find decoder comparison."""
    data = load("decoder_comparison.json")
    p = np.array(data["error_rates"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Logical error rates
    p_mwpm = np.array(data["mwpm_logical_rates"])
    p_uf = np.array(data["uf_logical_rates"])
    ci_mwpm = np.array(data["mwpm_ci"])
    ci_uf = np.array(data["uf_ci"])

    ax1.plot(p * 100, p_mwpm, "o-", color=COLORS[0], label="MWPM", linewidth=2, markersize=6)
    ax1.fill_between(p * 100, np.clip(p_mwpm - ci_mwpm, 0, 1),
                     np.clip(p_mwpm + ci_mwpm, 0, 1), alpha=0.2, color=COLORS[0])
    ax1.plot(p * 100, p_uf, "s--", color=COLORS[1], label="Union-Find", linewidth=2, markersize=6)
    ax1.fill_between(p * 100, np.clip(p_uf - ci_uf, 0, 1),
                     np.clip(p_uf + ci_uf, 0, 1), alpha=0.2, color=COLORS[1])

    ax1.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax1.set_ylabel("Logical Error Rate per Shot", fontsize=12)
    ax1.set_title(f"MWPM vs Union-Find Logical Error Rate\n(d={data['distance']}, {data['rounds']} rounds)",
                  fontsize=11)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Decoding time comparison
    t_mwpm = np.array(data["mwpm_decode_time_ms"]) * 1000  # convert to μs
    t_uf = np.array(data["uf_decode_time_ms"]) * 1000

    x = np.arange(len(p))
    width = 0.35
    bars1 = ax2.bar(x - width/2, t_mwpm, width, label="MWPM", color=COLORS[0], alpha=0.85)
    bars2 = ax2.bar(x + width/2, t_uf, width, label="Union-Find", color=COLORS[1], alpha=0.85)

    ax2.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax2.set_ylabel("Decode Time per Shot (μs)", fontsize=12)
    ax2.set_title("Decoding Time Comparison", fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{v*100:.1f}" for v in p], fontsize=9)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    return save_fig(fig, "fig3_decoder_comparison")


def plot_non_pauli_noise():
    """Figure 4: Non-Pauli noise impact."""
    data = load("non_pauli_noise.json")
    p = np.array(data["error_rates"])

    fig, ax = plt.subplots(figsize=(7, 5))

    noise_models = [
        ("depolarizing", "Depolarizing only", COLORS[0], "o", "-"),
        ("leakage", "Dep. + Leakage (10%)", COLORS[1], "s", "--"),
        ("meas_error", "Dep. + Meas. Error (2x)", COLORS[2], "^", "-."),
        ("combined", "Combined Non-Pauli", COLORS[5], "D", ":"),
    ]

    for key, label, color, marker, ls in noise_models:
        rates = np.array(data[key])
        valid = ~np.isnan(rates)
        if valid.any():
            ax.semilogy(p[valid] * 100, np.clip(rates[valid], 1e-5, 1),
                        marker=marker, label=label, color=color,
                        linestyle=ls, markersize=6, linewidth=1.8)

    ax.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax.set_ylabel("Logical Error Rate per Shot", fontsize=12)
    ax.set_title("Non-Pauli Noise Impact on Logical Error Rate", fontsize=13)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return save_fig(fig, "fig4_non_pauli_noise")


def plot_lattice_surgery():
    """Figure 5: Lattice surgery (logical CNOT) error rates."""
    data = load("lattice_surgery.json")
    p = np.array(data["error_rates"])
    distances = data["distances"]
    cnot_rates = np.array(data["cnot_rates"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for i, d in enumerate(distances):
        rates = cnot_rates[i]
        color = DISTANCE_COLORS.get(d, COLORS[i])
        ax1.semilogy(p * 100, np.clip(rates, 1e-6, 1), "o-",
                     label=f"d={d}", color=color, linewidth=2, markersize=6)

    # Reference line: p_cnot = p (no protection)
    ax1.semilogy(p * 100, p, "k:", linewidth=1.5, label="p_L = p (no coding)", alpha=0.7)

    ax1.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax1.set_ylabel("Logical CNOT Error Rate", fontsize=12)
    ax1.set_title("Lattice Surgery CNOT Error Rate", fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, which="both")

    # T-gate rates via magic state distillation
    t_gate_rates = data.get("t_gate_rates", {})
    for i, d in enumerate(distances):
        d_str = str(d)
        if d_str in t_gate_rates:
            t_rates = np.array(t_gate_rates[d_str])
            color = DISTANCE_COLORS.get(d, COLORS[i])
            valid = t_rates > 1e-10
            if valid.any():
                ax2.semilogy(p[valid] * 100, np.clip(t_rates[valid], 1e-10, 1),
                             "s--", label=f"d={d}", color=color, linewidth=2, markersize=6)

    ax2.set_xlabel("Physical Error Rate p (%)", fontsize=12)
    ax2.set_ylabel("Logical T-gate Error Rate (15-to-1 Distillation)", fontsize=12)
    ax2.set_title("Magic State Distillation (T-Gate)", fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return save_fig(fig, "fig5_lattice_surgery")


def plot_decoder_scaling():
    """Figure 6: Decoder performance vs code distance."""
    data = load("decoder_scaling.json")
    distances = data["distances"]
    n_qubits = data["num_qubits"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Logical rates vs distance
    p_mwpm = np.array(data["mwpm_logical_rates"])
    p_uf = np.array(data["uf_logical_rates"])

    ax1.semilogy(distances, np.clip(p_mwpm, 1e-5, 1), "o-", color=COLORS[0],
                 label="MWPM", linewidth=2, markersize=8)
    ax1.semilogy(distances, np.clip(p_uf, 1e-5, 1), "s--", color=COLORS[1],
                 label="Union-Find", linewidth=2, markersize=8)

    ax1.set_xlabel("Code Distance d", fontsize=12)
    ax1.set_ylabel("Logical Error Rate per Shot", fontsize=12)
    ax1.set_title(f"Logical Error Rate vs Distance (p={data['p_physical']*100:.1f}%)",
                  fontsize=12)
    ax1.set_xticks(distances)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, which="both")

    # Decode time vs num_qubits
    t_mwpm = np.array(data["mwpm_time_ms"]) * 1000
    t_uf = np.array(data["uf_time_ms"]) * 1000

    ax2.loglog(n_qubits, np.clip(t_mwpm, 1e-3, 1e6), "o-", color=COLORS[0],
               label="MWPM", linewidth=2, markersize=8)
    ax2.loglog(n_qubits, np.clip(t_uf, 1e-3, 1e6), "s--", color=COLORS[1],
               label="Union-Find", linewidth=2, markersize=8)

    # Annotate points with distance
    for i, (n, t_m, t_u, d) in enumerate(zip(n_qubits, t_mwpm, t_uf, distances)):
        ax2.annotate(f"d={d}", (n, t_m), textcoords="offset points",
                     xytext=(5, 5), fontsize=9)

    ax2.set_xlabel("Number of Data Qubits", fontsize=12)
    ax2.set_ylabel("Decode Time per Shot (μs)", fontsize=12)
    ax2.set_title("Decoder Time Scaling with System Size", fontsize=12)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return save_fig(fig, "fig6_decoder_scaling")


def plot_summary_panel():
    """Figure 7: Summary panel combining key results."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle("Surface Code Error Correction: Summary", fontsize=14, fontweight="bold")

    # Panel A: Threshold
    data_th = load("mwpm_threshold.json")
    ax = axes[0, 0]
    p = np.array(data_th["error_rates"])
    for i, d in enumerate(data_th["distances"]):
        rates = np.array(data_th["logical_rates"][i])
        ax.plot(p * 100, rates, "o-", label=f"d={d}",
                color=DISTANCE_COLORS.get(d, COLORS[i]), linewidth=1.8, markersize=4)
    th = data_th.get("threshold_estimate", 0.006)
    ax.axvline(th * 100, color="k", linestyle="--", linewidth=1.5, label=f"p_th={th*100:.2f}%")
    ax.set_xlabel("Physical Error Rate p (%)", fontsize=10)
    ax.set_ylabel("Logical Error Rate", fontsize=10)
    ax.set_title("(A) Threshold Analysis (MWPM)", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel B: Noise models
    data_nm = load("noise_model_comparison.json")
    ax = axes[0, 1]
    p_nm = np.array(data_nm["error_rates"])
    for key, label, color in [
        ("depolarizing", "Depolarizing", COLORS[0]),
        ("amplitude_damping", "Amplitude Damping", COLORS[1]),
        ("phase_damping", "Phase Damping", COLORS[2]),
    ]:
        rates = np.array(data_nm[key])
        ax.semilogy(p_nm * 100, np.clip(rates, 1e-5, 1), "o-",
                    label=label, color=color, linewidth=1.8, markersize=4)
    ax.set_xlabel("Physical Error Rate p (%)", fontsize=10)
    ax.set_ylabel("Logical Error Rate", fontsize=10)
    ax.set_title("(B) Noise Model Comparison", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")

    # Panel C: Non-Pauli effects
    data_np = load("non_pauli_noise.json")
    ax = axes[1, 0]
    p_np = np.array(data_np["error_rates"])
    for key, label, color in [
        ("depolarizing", "Depolarizing only", COLORS[0]),
        ("leakage", "With Leakage", COLORS[1]),
        ("meas_error", "With Meas. Error", COLORS[2]),
    ]:
        rates = np.array(data_np[key])
        valid = ~np.isnan(rates) & (np.array(rates) > 0)
        if valid.any():
            ax.semilogy(p_np[valid] * 100, np.clip(np.array(rates)[valid], 1e-5, 1),
                        "o-", label=label, color=color, linewidth=1.8, markersize=4)
    ax.set_xlabel("Physical Error Rate p (%)", fontsize=10)
    ax.set_ylabel("Logical Error Rate", fontsize=10)
    ax.set_title("(C) Non-Pauli Noise Impact", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")

    # Panel D: Lattice surgery
    data_ls = load("lattice_surgery.json")
    ax = axes[1, 1]
    p_ls = np.array(data_ls["error_rates"])
    for i, d in enumerate(data_ls["distances"]):
        rates = np.array(data_ls["cnot_rates"][i])
        ax.semilogy(p_ls * 100, np.clip(rates, 1e-6, 1), "o-",
                    label=f"CNOT d={d}", color=DISTANCE_COLORS.get(d, COLORS[i]),
                    linewidth=1.8, markersize=4)
    ax.set_xlabel("Physical Error Rate p (%)", fontsize=10)
    ax.set_ylabel("Logical CNOT Error Rate", fontsize=10)
    ax.set_title("(D) Lattice Surgery CNOT", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return save_fig(fig, "fig7_summary_panel")


if __name__ == "__main__":
    print("Generating figures...", flush=True)

    figs = []
    figs.extend(plot_noise_model_comparison())
    figs.extend(plot_threshold_analysis())
    figs.extend(plot_decoder_comparison())
    figs.extend(plot_non_pauli_noise())
    figs.extend(plot_lattice_surgery())
    figs.extend(plot_decoder_scaling())
    figs.extend(plot_summary_panel())

    print(f"\nAll figures saved to {FIGURES_DIR}", flush=True)
