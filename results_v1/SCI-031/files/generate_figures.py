"""
Figure generation: All plots for the VQE noise-resilience study.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json
import os

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'figure.dpi': 150,
})
COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00', '#56B4E9']

os.makedirs("figures", exist_ok=True)

# ── Fig 1: VQE Convergence ─────────────────────────────────────────────────

def plot_convergence(data_path="results/ansatz_comparison.json"):
    with open(data_path) as f:
        data = json.load(f)

    conv = data["vqe_convergence"]
    he_e = conv["hardware_efficient"]["energies"]
    uccsd_e = conv["uccsd"]["energies"]
    exact_e = conv["exact_energy_h2"]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    steps = range(len(he_e))
    ax.plot(steps, he_e, color=COLORS[0], linewidth=2, label="Hardware-Efficient Ansatz")
    ax.plot(range(len(uccsd_e)), uccsd_e, color=COLORS[1], linewidth=2, label="UCCSD-Inspired Ansatz")
    ax.axhline(exact_e, color='k', linestyle='--', linewidth=1.5, label=f"FCI Exact ({exact_e} Ha)")
    ax.fill_between(steps, exact_e - 0.0016, exact_e + 0.0016,
                    alpha=0.15, color='gray', label="Chemical Accuracy (±1.6 mHa)")
    ax.set_xlabel("Optimization Steps")
    ax.set_ylabel("Energy (Ha)")
    ax.set_title("VQE Convergence: Hardware-Efficient vs UCCSD-Inspired (H₂, STO-3G)")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/fig1_vqe_convergence.png", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig1_vqe_convergence.pdf", bbox_inches='tight')
    plt.close()
    print("Saved: figures/fig1_vqe_convergence.{png,pdf}")

# ── Fig 2: Measurement Cost ────────────────────────────────────────────────

def plot_measurement_cost(data_path="results/measurement_cost.json"):
    with open(data_path) as f:
        data = json.load(f)

    shot_data = data["shot_scaling"]
    n_terms = [d["n_terms"] for d in shot_data]
    naive = [d["naive_shots"] for d in shot_data]
    grouped = [d["grouped_shots"] for d in shot_data]
    shadow = [d["shadow_shots"] for d in shot_data]

    shadow_est = data["classical_shadow"]["shadow_estimates"]
    n_shadows = [s["n_shadows"] for s in shadow_est]
    abs_err = [s["abs_error"] for s in shadow_est]
    exact_e = data["classical_shadow"]["exact_energy"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: Shot scaling
    ax = axes[0]
    ax.loglog(n_terms, naive, 'o-', color=COLORS[0], linewidth=2, label="Naive (per term)")
    ax.loglog(n_terms, grouped, 's-', color=COLORS[1], linewidth=2, label="QWC Grouping")
    ax.loglog(n_terms, shadow, '^-', color=COLORS[2], linewidth=2, label="Classical Shadow")
    ax.set_xlabel("Number of Pauli Terms")
    ax.set_ylabel("Required Shots")
    ax.set_title("Measurement Shot Scaling (ε = 0.01)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: Shadow convergence
    ax = axes[1]
    ax.plot(n_shadows, abs_err, 'o-', color=COLORS[3], linewidth=2, markersize=8)
    ref_line = [0.01 / np.sqrt(n) * 50 for n in n_shadows]
    ax.plot(n_shadows, ref_line, '--', color='gray', linewidth=1.5, label="1/√N reference")
    ax.axhline(0.001, color='red', linestyle=':', linewidth=1.5, label="1 mHa threshold")
    ax.set_xlabel("Number of Shadow Samples")
    ax.set_ylabel("|Error| (Ha)")
    ax.set_title(f"Classical Shadow Energy Convergence\n(Exact: {exact_e:.4f} Ha)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig2_measurement_cost.png", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig2_measurement_cost.pdf", bbox_inches='tight')
    plt.close()
    print("Saved: figures/fig2_measurement_cost.{png,pdf}")

# ── Fig 3: Barren Plateau ─────────────────────────────────────────────────

def plot_barren_plateau(data_path="results/barren_plateau.json"):
    with open(data_path) as f:
        data = json.load(f)

    depth_deep = data["depth_analysis"]["deep_circuit"]
    depth_shallow = data["depth_analysis"]["shallow_local"]
    qubit_data = data["qubit_size_analysis"]
    theory = data["theoretical_scaling"]
    init_rand = data["init_comparison"]["random_init"]
    init_ident = data["init_comparison"]["identity_block_init"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel A: Gradient variance vs depth
    ax = axes[0]
    depths_d = [d["depth"] for d in depth_deep]
    vars_d = [d["grad_variance"] for d in depth_deep]
    vars_s = [d["grad_variance"] for d in depth_shallow]
    ax.semilogy(depths_d, vars_d, 'o-', color=COLORS[0], linewidth=2, label="Deep Random Circuit")
    ax.semilogy(depths_d, vars_s, 's-', color=COLORS[1], linewidth=2, label="Shallow Local Circuit")
    ax.set_xlabel("Circuit Depth (# Layers)")
    ax.set_ylabel("Gradient Variance")
    ax.set_title("Barren Plateau: Gradient Variance vs Depth")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel B: Gradient variance vs system size
    ax = axes[1]
    n_qubits_list = [d["n_qubits"] for d in qubit_data]
    vars_q = [d["grad_variance"] for d in qubit_data]
    theory_n = theory["n_qubits"]
    theory_v = theory["variance_scaling"]
    ax.semilogy(n_qubits_list, vars_q, 'o-', color=COLORS[2], linewidth=2, label="Observed (depth=4)")
    ax.semilogy(theory_n[:len(n_qubits_list)], theory_v[:len(n_qubits_list)],
                '--', color='gray', linewidth=2, label="Theory: 2⁻ⁿ")
    ax.set_xlabel("Number of Qubits")
    ax.set_ylabel("Gradient Variance")
    ax.set_title("Barren Plateau: Gradient Variance vs System Size")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel C: Init comparison
    ax = axes[2]
    steps_r = range(len(init_rand))
    steps_i = range(len(init_ident))
    ax.plot(steps_r, init_rand, color=COLORS[0], linewidth=2, label="Random Init")
    ax.plot(steps_i, init_ident, color=COLORS[1], linewidth=2, label="Identity-Block Init")
    ax.set_xlabel("Optimization Steps")
    ax.set_ylabel("Energy (arb. units)")
    ax.set_title("Initialization Strategies Comparison")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig3_barren_plateau.png", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig3_barren_plateau.pdf", bbox_inches='tight')
    plt.close()
    print("Saved: figures/fig3_barren_plateau.{png,pdf}")

# ── Fig 4: Error Mitigation ────────────────────────────────────────────────

def plot_error_mitigation(data_path="results/error_mitigation.json"):
    with open(data_path) as f:
        data = json.load(f)

    results = data["results"]
    noise_levels = [r["noise_level"] for r in results]
    exact = results[0]["exact"]
    errors = {
        "Noisy (baseline)": [r["error_noisy"] * 1000 for r in results],
        "ZNE": [r["error_zne"] * 1000 for r in results],
        "PEC": [r["error_pec"] * 1000 for r in results],
        "CDR": [r["error_cdr"] * 1000 for r in results],
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel A: Error vs noise level
    ax = axes[0]
    for (label, err_list), color in zip(errors.items(), COLORS):
        ax.semilogy([n * 100 for n in noise_levels], err_list,
                    'o-', color=color, linewidth=2, label=label)
    ax.axhline(1.6, color='red', linestyle='--', linewidth=1.5,
               label="Chemical Accuracy (1.6 mHa)")
    ax.set_xlabel("Noise Rate (%)")
    ax.set_ylabel("Absolute Error (mHa)")
    ax.set_title("Error Mitigation Methods: Accuracy vs Noise")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel B: Bar chart at moderate noise (noise_levels[2])
    ax = axes[1]
    if len(results) > 2:
        mid = results[2]  # moderate noise
        methods = ['Noisy', 'ZNE', 'PEC', 'CDR']
        values = [mid[f'error_{m.lower()}'] * 1000 for m in methods]
        bars = ax.bar(methods, values, color=COLORS[:4], edgecolor='black', linewidth=0.8)
        ax.axhline(1.6, color='red', linestyle='--', linewidth=1.5, label="Chemical Accuracy")
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    f'{val:.2f}', ha='center', va='bottom', fontsize=10)
        ax.set_ylabel("Absolute Error (mHa)")
        ax.set_title(f"Error Comparison at Noise = {mid['noise_level']*100:.1f}%")
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig("figures/fig4_error_mitigation.png", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig4_error_mitigation.pdf", bbox_inches='tight')
    plt.close()
    print("Saved: figures/fig4_error_mitigation.{png,pdf}")

# ── Fig 5: Fermion Mapping ────────────────────────────────────────────────

def plot_fermion_mapping(data_path="results/fermion_mapping.json"):
    with open(data_path) as f:
        data = json.load(f)

    scaling = data["scaling_analysis"]
    n_modes = [s["n_modes"] for s in scaling]
    jw_w = [s["jw_weight"] for s in scaling]
    bk_w = [s["bk_weight"] for s in scaling]
    par_w = [s["parity_weight"] for s in scaling]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel A: Pauli weight scaling
    ax = axes[0]
    ax.plot(n_modes, jw_w, 'o-', color=COLORS[0], linewidth=2, label="Jordan-Wigner O(n)")
    ax.plot(n_modes, bk_w, 's-', color=COLORS[1], linewidth=2, label="Bravyi-Kitaev O(log n)")
    ax.plot(n_modes, par_w, '^-', color=COLORS[2], linewidth=2, label="Parity O(log n)")
    n_arr = np.array(n_modes)
    ax.plot(n_modes, np.log2(n_arr + 1), '--', color='gray', linewidth=1.5, label="log₂(n) reference")
    ax.set_xlabel("Number of Modes (Spin-Orbitals)")
    ax.set_ylabel("Max Pauli Weight")
    ax.set_title("Fermion-Qubit Mapping: Pauli Weight Scaling")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel B: Molecule comparison
    ax = axes[1]
    molecules = ['H₂', 'LiH', 'H₂O']
    jw_qubits = [data.get("h2", {}).get("n_qubits_jw", 4),
                 data.get("lih", {}).get("n_qubits_jw", 4),
                 4]
    bk_qubits = [data.get("h2", {}).get("n_qubits_bk", 4),
                 data.get("lih", {}).get("n_qubits_bk", 4),
                 4]
    par_qubits = [data.get("h2", {}).get("n_qubits_parity", 2),
                  data.get("lih", {}).get("n_qubits_parity", 2),
                  2]

    x = np.arange(len(molecules))
    width = 0.25
    ax.bar(x - width, jw_qubits, width, label="Jordan-Wigner", color=COLORS[0], edgecolor='black')
    ax.bar(x, bk_qubits, width, label="Bravyi-Kitaev", color=COLORS[1], edgecolor='black')
    ax.bar(x + width, par_qubits, width, label="Parity (2-qubit reduction)", color=COLORS[2], edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(molecules)
    ax.set_ylabel("Number of Qubits Required")
    ax.set_title("Qubit Requirements by Mapping and Molecule")
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig("figures/fig5_fermion_mapping.png", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig5_fermion_mapping.pdf", bbox_inches='tight')
    plt.close()
    print("Saved: figures/fig5_fermion_mapping.{png,pdf}")

# ── Fig 6: Benchmark Summary ──────────────────────────────────────────────

def plot_benchmark(data_path="results/benchmark.json"):
    with open(data_path) as f:
        data = json.load(f)

    summary = data["summary"]
    molecules = ['H2', 'LiH', 'H2O']
    ansatze = ['hardware_efficient', 'uccsd']
    labels_ans = ['Hardware-Efficient', 'UCCSD-Inspired']

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: Energy errors (noiseless)
    ax = axes[0]
    x = np.arange(len(molecules))
    width = 0.35
    for ai, (ansatz, label_ans) in enumerate(zip(ansatze, labels_ans)):
        errors = []
        for mol in molecules:
            key = f"{mol.lower()}_{ansatz}"
            if key in summary:
                errors.append(summary[key]["error_mHa"])
            else:
                errors.append(0.0)
        offset = (ai - 0.5) * width
        ax.bar(x + offset, errors, width, label=label_ans,
               color=COLORS[ai], edgecolor='black', alpha=0.85)

    ax.axhline(1.6, color='red', linestyle='--', linewidth=2, label="Chemical Accuracy (1.6 mHa)")
    ax.set_xticks(x)
    ax.set_xticklabels(molecules)
    ax.set_ylabel("Energy Error (mHa)")
    ax.set_title("VQE Energy Error vs Classical Reference (Noiseless)")
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Panel B: Noise sensitivity for H2 UCCSD
    ax = axes[1]
    bench_h2 = data["benchmark_results"].get("h2", {})
    if "uccsd" in bench_h2:
        noise_res = bench_h2["uccsd"]["noise_results"]
        noise_lvls = [r["noise_level"] * 100 for r in noise_res]
        energies = [r["final_energy"] for r in noise_res]
        errors_mha = [r["error_hartree"] * 1000 for r in noise_res]

        ax2 = ax.twinx()
        line1 = ax.plot(noise_lvls, energies, 'o-', color=COLORS[0], linewidth=2, label="VQE Energy")
        ax.axhline(bench_h2["uccsd"]["ref_energy"], color='k', linestyle='--',
                   linewidth=1.5, label="FCI Reference")
        line2 = ax2.plot(noise_lvls, errors_mha, 's--', color=COLORS[3], linewidth=2, label="Error (mHa)")
        ax2.axhline(1.6, color='red', linestyle=':', linewidth=1.5)
        ax.set_xlabel("Noise Rate (%)")
        ax.set_ylabel("Energy (Ha)")
        ax2.set_ylabel("Absolute Error (mHa)", color=COLORS[3])
        ax.set_title("H₂ UCCSD-VQE: Noise Sensitivity")
        lines = line1 + line2
        labs = [l.get_label() for l in lines]
        ax.legend(lines, labs, loc='upper left')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig6_benchmark.png", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig6_benchmark.pdf", bbox_inches='tight')
    plt.close()
    print("Saved: figures/fig6_benchmark.{png,pdf}")

# ── Main ─────────────────────────────────────────────────────────────────

def main():
    print("Generating all figures...")
    plot_convergence()
    plot_measurement_cost()
    plot_barren_plateau()
    plot_error_mitigation()
    plot_fermion_mapping()
    plot_benchmark()
    print("All figures generated.")

if __name__ == "__main__":
    main()
