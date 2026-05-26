#!/usr/bin/env python3
"""
Main experiment script — runs all experiments and generates figures.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import json

from src.parts_catalog import PartsCatalog, PartType
from src.circuit_spec import (
    make_toggle_switch, make_repressilator, CircuitSpec,
    LogicGate, GateType, Signal
)
from src.stochastic_sim import (
    build_toggle_switch_model, build_repressilator_model,
    gillespie_ssa, tau_leaping
)
from src.robust_design import (
    robustness_score, optimize_circuit_params,
    toggle_switch_bistability_score, repressilator_oscillation_score,
    latin_hypercube_sample
)
from src.context_effects import ContextPredictor
from src.design_pipeline import DesignPipeline

FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'figure.dpi': 150,
})

print("=" * 60)
print("Synthetic Gene Circuit Design Framework — Experiments")
print("=" * 60)

# ================================================================
# Experiment 1: Toggle Switch — Gillespie vs Tau-Leaping
# ================================================================
print("\n[Exp 1] Toggle Switch: Gillespie vs Tau-Leaping...")

ts_model = build_toggle_switch_model({"IPTG": 0.0, "aTc": 0.0})

# Gillespie SSA
t_gill, s_gill = gillespie_ssa(ts_model, t_end=300.0, seed=42)
# Tau-leaping
t_tau, s_tau = tau_leaping(ts_model, t_end=300.0, tau=0.5, seed=42)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Toggle Switch: Stochastic Simulation Comparison", fontweight='bold')

# Gillespie
axes[0, 0].plot(t_gill, s_gill[:, 0], 'b-', alpha=0.8, label='LacI', linewidth=0.8)
axes[0, 0].plot(t_gill, s_gill[:, 1], 'r-', alpha=0.8, label='TetR', linewidth=0.8)
axes[0, 0].set_title('Gillespie SSA')
axes[0, 0].set_xlabel('Time (min)')
axes[0, 0].set_ylabel('Molecule count')
axes[0, 0].legend()

# Tau-leaping
axes[0, 1].plot(t_tau, s_tau[:, 0], 'b-', alpha=0.8, label='LacI', linewidth=0.8)
axes[0, 1].plot(t_tau, s_tau[:, 1], 'r-', alpha=0.8, label='TetR', linewidth=0.8)
axes[0, 1].set_title('Tau-Leaping (τ=0.5)')
axes[0, 1].set_xlabel('Time (min)')
axes[0, 1].set_ylabel('Molecule count')
axes[0, 1].legend()

# Multiple Gillespie trajectories for bistability visualization
lacI_final = []
tetR_final = []
for seed in range(20):
    t_g, s_g = tau_leaping(ts_model, t_end=500.0, tau=0.5, seed=seed)
    n = len(t_g)
    start = int(0.8 * n)
    lacI_final.append(np.mean(s_g[start:, 0]))
    tetR_final.append(np.mean(s_g[start:, 1]))

axes[1, 0].scatter(lacI_final, tetR_final, c='purple', alpha=0.7, s=60)
axes[1, 0].set_xlabel('LacI steady-state')
axes[1, 0].set_ylabel('TetR steady-state')
axes[1, 0].set_title('Bistability Phase Space (20 runs)')
axes[1, 0].axline((0, 0), slope=1, color='gray', linestyle='--', alpha=0.5)

# Inducer response
iptg_vals = np.linspace(0, 1, 15)
lacI_means = []
tetR_means = []
for iptg in iptg_vals:
    m = build_toggle_switch_model({"IPTG": iptg, "aTc": 0.0})
    t_s, s_s = tau_leaping(m, t_end=400.0, tau=0.5, seed=42)
    n = len(t_s)
    start = int(0.7 * n)
    lacI_means.append(np.mean(s_s[start:, 0]))
    tetR_means.append(np.mean(s_s[start:, 1]))

axes[1, 1].plot(iptg_vals, lacI_means, 'bo-', label='LacI', markersize=5)
axes[1, 1].plot(iptg_vals, tetR_means, 'rs-', label='TetR', markersize=5)
axes[1, 1].set_xlabel('IPTG concentration')
axes[1, 1].set_ylabel('Steady-state level')
axes[1, 1].set_title('Inducer Dose-Response')
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig1_toggle_switch_simulation.png', bbox_inches='tight')
plt.close()
print("  -> Figure 1 saved.")

# ================================================================
# Experiment 2: Repressilator Oscillation
# ================================================================
print("\n[Exp 2] Repressilator Oscillation...")

rep_model = build_repressilator_model()
t_rep_g, s_rep_g = gillespie_ssa(rep_model, t_end=500.0, seed=42, max_steps=5_000_000)
t_rep_t, s_rep_t = tau_leaping(rep_model, t_end=500.0, tau=0.2, seed=42)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Repressilator: Stochastic Oscillation Dynamics", fontweight='bold')

# Gillespie proteins
for i, (name, color) in enumerate([("LacI", "blue"), ("TetR", "red"), ("cI", "green")]):
    idx = rep_model.species.index(name)
    axes[0, 0].plot(t_rep_g, s_rep_g[:, idx], color=color, alpha=0.8,
                    label=name, linewidth=0.8)
axes[0, 0].set_title('Gillespie SSA — Proteins')
axes[0, 0].set_xlabel('Time (min)')
axes[0, 0].set_ylabel('Molecule count')
axes[0, 0].legend()

# Tau-leaping proteins
for i, (name, color) in enumerate([("LacI", "blue"), ("TetR", "red"), ("cI", "green")]):
    idx = rep_model.species.index(name)
    axes[0, 1].plot(t_rep_t, s_rep_t[:, idx], color=color, alpha=0.8,
                    label=name, linewidth=0.8)
axes[0, 1].set_title('Tau-Leaping (τ=0.2) — Proteins')
axes[0, 1].set_xlabel('Time (min)')
axes[0, 1].set_ylabel('Molecule count')
axes[0, 1].legend()

# mRNA dynamics
for i, (name, color) in enumerate([("mRNA_lacI", "blue"), ("mRNA_tetR", "red"), ("mRNA_cI", "green")]):
    idx = rep_model.species.index(name)
    axes[1, 0].plot(t_rep_t, s_rep_t[:, idx], color=color, alpha=0.8,
                    label=name, linewidth=0.8)
axes[1, 0].set_title('mRNA Dynamics (Tau-Leaping)')
axes[1, 0].set_xlabel('Time (min)')
axes[1, 0].set_ylabel('mRNA count')
axes[1, 0].legend()

# Phase portrait LacI vs TetR
lacI_idx = rep_model.species.index("LacI")
tetR_idx = rep_model.species.index("TetR")
axes[1, 1].plot(s_rep_t[:, lacI_idx], s_rep_t[:, tetR_idx],
                'k-', alpha=0.3, linewidth=0.5)
axes[1, 1].plot(s_rep_t[0, lacI_idx], s_rep_t[0, tetR_idx],
                'go', markersize=10, label='Start')
axes[1, 1].plot(s_rep_t[-1, lacI_idx], s_rep_t[-1, tetR_idx],
                'rs', markersize=10, label='End')
axes[1, 1].set_xlabel('LacI')
axes[1, 1].set_ylabel('TetR')
axes[1, 1].set_title('Phase Portrait')
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig2_repressilator_oscillation.png', bbox_inches='tight')
plt.close()
print("  -> Figure 2 saved.")

# ================================================================
# Experiment 3: Robustness Analysis
# ================================================================
print("\n[Exp 3] Robustness Analysis under Parameter Uncertainty...")

# Toggle switch robustness
ts_param_ranges = {
    "alpha1": (1.5, 5.0),
    "alpha2": (1.5, 5.0),
    "K1": (20.0, 80.0),
    "K2": (20.0, 80.0),
    "n1": (1.5, 3.5),
    "n2": (1.5, 3.5),
    "delta1": (0.02, 0.1),
    "delta2": (0.02, 0.1),
}

ts_mean, ts_std, ts_scores = robustness_score(
    build_toggle_switch_model, ts_param_ranges,
    toggle_switch_bistability_score,
    n_samples=60, t_end=500.0, tau=0.5
)
print(f"  Toggle switch robustness: {ts_mean:.4f} ± {ts_std:.4f}")

# Repressilator robustness
rep_param_ranges = {
    "alpha": (1.5, 5.0),
    "alpha0": (0.001, 0.05),
    "K": (20.0, 80.0),
    "n": (1.5, 3.5),
    "delta_m": (0.05, 0.2),
    "delta_p": (0.01, 0.05),
    "beta": (0.2, 1.0),
}

rep_mean, rep_std, rep_scores = robustness_score(
    build_repressilator_model, rep_param_ranges,
    repressilator_oscillation_score,
    n_samples=60, t_end=500.0, tau=0.3
)
print(f"  Repressilator robustness: {rep_mean:.4f} ± {rep_std:.4f}")

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle("Robustness Analysis under Parameter Uncertainty", fontweight='bold')

axes[0].hist(ts_scores, bins=15, color='steelblue', edgecolor='black', alpha=0.8)
axes[0].axvline(ts_mean, color='red', linestyle='--', label=f'Mean={ts_mean:.3f}')
axes[0].set_xlabel('Bistability Score')
axes[0].set_ylabel('Count')
axes[0].set_title('Toggle Switch Robustness')
axes[0].legend()

axes[1].hist(rep_scores, bins=15, color='coral', edgecolor='black', alpha=0.8)
axes[1].axvline(rep_mean, color='red', linestyle='--', label=f'Mean={rep_mean:.3f}')
axes[1].set_xlabel('Oscillation Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Repressilator Robustness')
axes[1].legend()

# Sensitivity analysis: vary one parameter at a time
param_names = list(ts_param_ranges.keys())
sensitivities = []
for pname in param_names:
    lo, hi = ts_param_ranges[pname]
    scores_lo = []
    scores_hi = []
    for _ in range(10):
        seed = np.random.randint(1000)
        model_lo = build_toggle_switch_model({pname: lo})
        t_l, s_l = tau_leaping(model_lo, 400.0, 0.5, seed)
        scores_lo.append(toggle_switch_bistability_score(t_l, s_l, model_lo.species))

        model_hi = build_toggle_switch_model({pname: hi})
        t_h, s_h = tau_leaping(model_hi, 400.0, 0.5, seed)
        scores_hi.append(toggle_switch_bistability_score(t_h, s_h, model_hi.species))

    sensitivity = abs(np.mean(scores_hi) - np.mean(scores_lo))
    sensitivities.append(sensitivity)

axes[2].barh(param_names, sensitivities, color='teal', edgecolor='black')
axes[2].set_xlabel('Sensitivity (Δ score)')
axes[2].set_title('Parameter Sensitivity (Toggle)')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig3_robustness_analysis.png', bbox_inches='tight')
plt.close()
print("  -> Figure 3 saved.")

# ================================================================
# Experiment 4: Parameter Optimization
# ================================================================
print("\n[Exp 4] Evolutionary Parameter Optimization...")

# Toggle switch optimization
best_ts_params, best_ts_score, ts_opt_history = optimize_circuit_params(
    build_toggle_switch_model, ts_param_ranges,
    toggle_switch_bistability_score,
    n_iterations=25, n_samples_per_iter=15,
    t_end=400.0, tau=0.5, seed=42
)
print(f"  Toggle switch best score: {best_ts_score:.4f}")
print(f"  Best params: { {k: round(v, 3) for k, v in best_ts_params.items()} }")

# Repressilator optimization
best_rep_params, best_rep_score, rep_opt_history = optimize_circuit_params(
    build_repressilator_model, rep_param_ranges,
    repressilator_oscillation_score,
    n_iterations=25, n_samples_per_iter=15,
    t_end=500.0, tau=0.3, seed=42
)
print(f"  Repressilator best score: {best_rep_score:.4f}")
print(f"  Best params: { {k: round(v, 3) for k, v in best_rep_params.items()} }")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
fig.suptitle("Evolutionary Optimization Convergence", fontweight='bold')

axes[0].plot(range(1, len(ts_opt_history) + 1), ts_opt_history, 'bo-', markersize=4)
axes[0].set_xlabel('Generation')
axes[0].set_ylabel('Best Bistability Score')
axes[0].set_title('Toggle Switch Optimization')
axes[0].grid(True, alpha=0.3)

axes[1].plot(range(1, len(rep_opt_history) + 1), rep_opt_history, 'rs-', markersize=4)
axes[1].set_xlabel('Generation')
axes[1].set_ylabel('Best Oscillation Score')
axes[1].set_title('Repressilator Optimization')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig4_optimization_convergence.png', bbox_inches='tight')
plt.close()
print("  -> Figure 4 saved.")

# ================================================================
# Experiment 5: Optimized vs Default Circuit Comparison
# ================================================================
print("\n[Exp 5] Optimized vs Default Circuit Comparison...")

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Default vs Optimized Circuit Performance", fontweight='bold')

# Toggle switch: default vs optimized
ts_default = build_toggle_switch_model()
t_d, s_d = tau_leaping(ts_default, 500.0, 0.5, seed=42)
ts_opt = build_toggle_switch_model(best_ts_params)
t_o, s_o = tau_leaping(ts_opt, 500.0, 0.5, seed=42)

axes[0, 0].plot(t_d, s_d[:, 0], 'b-', alpha=0.7, label='LacI')
axes[0, 0].plot(t_d, s_d[:, 1], 'r-', alpha=0.7, label='TetR')
axes[0, 0].set_title('Toggle Switch — Default Parameters')
axes[0, 0].set_xlabel('Time (min)')
axes[0, 0].set_ylabel('Molecules')
axes[0, 0].legend()

axes[0, 1].plot(t_o, s_o[:, 0], 'b-', alpha=0.7, label='LacI')
axes[0, 1].plot(t_o, s_o[:, 1], 'r-', alpha=0.7, label='TetR')
axes[0, 1].set_title('Toggle Switch — Optimized Parameters')
axes[0, 1].set_xlabel('Time (min)')
axes[0, 1].set_ylabel('Molecules')
axes[0, 1].legend()

# Repressilator: default vs optimized
rep_default = build_repressilator_model()
t_rd, s_rd = tau_leaping(rep_default, 500.0, 0.2, seed=42)
rep_opt_model = build_repressilator_model(best_rep_params)
t_ro, s_ro = tau_leaping(rep_opt_model, 500.0, 0.2, seed=42)

for name, color in [("LacI", "blue"), ("TetR", "red"), ("cI", "green")]:
    idx = rep_default.species.index(name)
    axes[1, 0].plot(t_rd, s_rd[:, idx], color=color, alpha=0.7, label=name)
axes[1, 0].set_title('Repressilator — Default Parameters')
axes[1, 0].set_xlabel('Time (min)')
axes[1, 0].set_ylabel('Molecules')
axes[1, 0].legend()

for name, color in [("LacI", "blue"), ("TetR", "red"), ("cI", "green")]:
    idx = rep_opt_model.species.index(name)
    axes[1, 1].plot(t_ro, s_ro[:, idx], color=color, alpha=0.7, label=name)
axes[1, 1].set_title('Repressilator — Optimized Parameters')
axes[1, 1].set_xlabel('Time (min)')
axes[1, 1].set_ylabel('Molecules')
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig5_default_vs_optimized.png', bbox_inches='tight')
plt.close()
print("  -> Figure 5 saved.")

# ================================================================
# Experiment 6: Context Effects Analysis
# ================================================================
print("\n[Exp 6] Genetic Context Effects Analysis...")

predictor = ContextPredictor()

# Toggle switch assembly
ts_spec = make_toggle_switch()
ts_assembly = ["pTet", "B0034", "LacI", "B0015", "pLac", "B0034", "TetR", "B0015"]
ts_effects = predictor.compute_circuit_context_effects(ts_assembly)
ts_insulation = predictor.insulation_recommendation(ts_assembly)

# Repressilator assembly
rep_spec = make_repressilator()
rep_assembly = ["pLambda", "B0034", "LacI", "B0015",
                "pLac", "B0034", "TetR", "B0015",
                "pTet", "B0034", "cI", "B0015"]
rep_effects = predictor.compute_circuit_context_effects(rep_assembly)
rep_insulation = predictor.insulation_recommendation(rep_assembly)

# With B0010 (weaker terminator) to show context effect
ts_assembly_weak = ["pTet", "B0034", "LacI", "B0010", "pLac", "B0034", "TetR", "B0010"]
ts_effects_weak = predictor.compute_circuit_context_effects(ts_assembly_weak)
ts_insulation_weak = predictor.insulation_recommendation(ts_assembly_weak)

print(f"  Toggle (B0015): {len(ts_insulation)} issues found")
print(f"  Toggle (B0010): {len(ts_insulation_weak)} issues found")
print(f"  Repressilator: {len(rep_insulation)} issues found")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Genetic Context Effects Analysis", fontweight='bold')

# Context effects for toggle with B0015
labels_ts = [f"{u}→{d}" for u, d, _ in ts_effects]
fcs_ts = [fc for _, _, fc in ts_effects]
colors_ts = ['green' if fc >= 0.9 else 'orange' if fc >= 0.8 else 'red' for fc in fcs_ts]
axes[0].barh(labels_ts, fcs_ts, color=colors_ts, edgecolor='black')
axes[0].axvline(1.0, color='gray', linestyle='--')
axes[0].axvline(0.9, color='red', linestyle=':', alpha=0.5)
axes[0].set_xlabel('Expression Fold Change')
axes[0].set_title('Toggle Switch (B0015)')
axes[0].set_xlim(0.6, 1.2)

# Context effects for toggle with B0010 (weaker terminator)
labels_tw = [f"{u}→{d}" for u, d, _ in ts_effects_weak]
fcs_tw = [fc for _, _, fc in ts_effects_weak]
colors_tw = ['green' if fc >= 0.9 else 'orange' if fc >= 0.8 else 'red' for fc in fcs_tw]
axes[1].barh(labels_tw, fcs_tw, color=colors_tw, edgecolor='black')
axes[1].axvline(1.0, color='gray', linestyle='--')
axes[1].axvline(0.9, color='red', linestyle=':', alpha=0.5)
axes[1].set_xlabel('Expression Fold Change')
axes[1].set_title('Toggle Switch (B0010 — weak term.)')
axes[1].set_xlim(0.6, 1.2)

# Context effects for repressilator
labels_rep = [f"{u}→{d}" for u, d, _ in rep_effects]
fcs_rep = [fc for _, _, fc in rep_effects]
colors_rep = ['green' if fc >= 0.9 else 'orange' if fc >= 0.8 else 'red' for fc in fcs_rep]
axes[2].barh(labels_rep, fcs_rep, color=colors_rep, edgecolor='black')
axes[2].axvline(1.0, color='gray', linestyle='--')
axes[2].axvline(0.9, color='red', linestyle=':', alpha=0.5)
axes[2].set_xlabel('Expression Fold Change')
axes[2].set_title('Repressilator (B0015)')
axes[2].set_xlim(0.6, 1.2)

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig6_context_effects.png', bbox_inches='tight')
plt.close()
print("  -> Figure 6 saved.")

# ================================================================
# Experiment 7: Design Pipeline End-to-End
# ================================================================
print("\n[Exp 7] Full Design Pipeline Execution...")

pipeline = DesignPipeline()

# Toggle switch pipeline
ts_result = pipeline.design_circuit(ts_spec, t_end=500.0)
print(f"  Toggle switch: {ts_result['n_gates']} gates, {ts_result['n_feedbacks']} feedbacks")
print(f"  Context issues: {len(ts_result['insulation_recommendations'])}")

# Repressilator pipeline
rep_result = pipeline.design_circuit(rep_spec, t_end=500.0)
print(f"  Repressilator: {rep_result['n_gates']} gates, {rep_result['n_feedbacks']} feedbacks")
print(f"  Context issues: {len(rep_result['insulation_recommendations'])}")

# NAND gate example
nand_spec = CircuitSpec(name="nand_gate")
nand_spec.inputs = [Signal("A", is_input=True), Signal("B", is_input=True)]
nand_spec.outputs = [Signal("Y", is_output=True)]
nand_spec.gates = [
    LogicGate("G1", GateType.NAND, ["A", "B"], "Y",
              promoter="pTac", rbs="B0034", cds="GFP", terminator="B0015"),
]
nand_result = pipeline.design_circuit(nand_spec, t_end=300.0)
nand_tt = nand_spec.get_truth_table()

print(f"\n  NAND Gate Truth Table:")
for inputs, outputs in nand_tt.items():
    print(f"    {inputs} -> {outputs}")

print(f"\n  NAND Verilog:\n{nand_result.get('verilog', 'N/A')}")

# ================================================================
# Experiment 8: Tau-leaping step size comparison
# ================================================================
print("\n[Exp 8] Tau-Leaping Step Size Comparison...")

tau_values = [0.05, 0.1, 0.5, 1.0, 2.0, 5.0]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle("Tau-Leaping: Effect of Step Size (τ) on Toggle Switch", fontweight='bold')

for idx, tau_val in enumerate(tau_values):
    ax = axes[idx // 3, idx % 3]
    m = build_toggle_switch_model()
    t_s, s_s = tau_leaping(m, t_end=300.0, tau=tau_val, seed=42)
    ax.plot(t_s, s_s[:, 0], 'b-', alpha=0.8, label='LacI', linewidth=0.8)
    ax.plot(t_s, s_s[:, 1], 'r-', alpha=0.8, label='TetR', linewidth=0.8)
    ax.set_title(f'τ = {tau_val}')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Molecules')
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig7_tau_step_comparison.png', bbox_inches='tight')
plt.close()
print("  -> Figure 7 saved.")

# ================================================================
# Experiment 9: Robustness heatmap for key parameter pairs
# ================================================================
print("\n[Exp 9] Parameter Robustness Heatmap...")

alpha_range = np.linspace(1.5, 5.0, 12)
K_range = np.linspace(20.0, 80.0, 12)
heatmap = np.zeros((len(alpha_range), len(K_range)))

for i, alpha in enumerate(alpha_range):
    for j, K in enumerate(K_range):
        m = build_toggle_switch_model({"alpha1": alpha, "alpha2": alpha, "K1": K, "K2": K})
        t_s, s_s = tau_leaping(m, t_end=400.0, tau=0.5, seed=42)
        heatmap[i, j] = toggle_switch_bistability_score(t_s, s_s, m.species)

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
im = ax.imshow(heatmap, origin='lower', aspect='auto',
               extent=[K_range[0], K_range[-1], alpha_range[0], alpha_range[-1]],
               cmap='viridis')
ax.set_xlabel('K (half-saturation constant)')
ax.set_ylabel('α (max expression rate)')
ax.set_title('Toggle Switch Bistability Score\nacross α and K', fontweight='bold')
plt.colorbar(im, ax=ax, label='Bistability Score')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig8_robustness_heatmap.png', bbox_inches='tight')
plt.close()
print("  -> Figure 8 saved.")

# ================================================================
# Summary Statistics
# ================================================================
print("\n" + "=" * 60)
print("EXPERIMENT SUMMARY")
print("=" * 60)

summary = {
    "toggle_switch": {
        "default_bistability_score": round(toggle_switch_bistability_score(
            t_d, s_d, ts_default.species), 4),
        "optimized_bistability_score": round(toggle_switch_bistability_score(
            t_o, s_o, ts_opt.species), 4),
        "robustness_mean": round(ts_mean, 4),
        "robustness_std": round(ts_std, 4),
        "best_optimized_params": {k: round(v, 4) for k, v in best_ts_params.items()},
        "context_issues_B0015": len(ts_insulation),
        "context_issues_B0010": len(ts_insulation_weak),
    },
    "repressilator": {
        "default_oscillation_score": round(repressilator_oscillation_score(
            t_rd, s_rd, rep_default.species), 4),
        "optimized_oscillation_score": round(repressilator_oscillation_score(
            t_ro, s_ro, rep_opt_model.species), 4),
        "robustness_mean": round(rep_mean, 4),
        "robustness_std": round(rep_std, 4),
        "best_optimized_params": {k: round(v, 4) for k, v in best_rep_params.items()},
        "context_issues": len(rep_insulation),
    },
    "nand_gate_truth_table": {str(k): v for k, v in nand_tt.items()},
}

with open("experiment_results.json", "w") as f:
    json.dump(summary, f, indent=2)

print(json.dumps(summary, indent=2))
print(f"\nAll figures saved to {FIGDIR}/")
print("Results saved to experiment_results.json")
print("Done!")
