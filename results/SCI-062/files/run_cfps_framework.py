#!/usr/bin/env python3
"""
CFPS Productivity Optimization Framework — Master Runner
Executes all 6 modules and generates comprehensive results.
"""

import os, sys, json
from datetime import datetime

os.chdir(os.path.dirname(os.path.abspath(__file__)))

print("=" * 70)
print("  CFPS Productivity Optimization Framework")
print(f"  Execution started: {datetime.now().isoformat()}")
print("=" * 70)

# Ensure directories
for d in ['figures', 'results', 'data', 'logs']:
    os.makedirs(d, exist_ok=True)

# --- Module 1 ---
print("\n[1/6] Transcription-Translation Coupled Model...")
from cfps_transcription_translation import run_simulation, plot_dynamics, resource_competition_analysis
sol = run_simulation()
plot_dynamics(sol)
m1_results = resource_competition_analysis()
with open('results/m1_txn_tln_results.json', 'w') as f:
    json.dump(m1_results, f, indent=2)
print(f"  ✓ Optimal DNA: {m1_results['optimal_DNA_nM']:.1f} nM, Max yield: {m1_results['max_protein_yield_nM']:.1f} nM")

# --- Module 2 ---
print("\n[2/6] Energy Regeneration System Comparison...")
from cfps_energy_regeneration import run_energy_comparison, plot_energy_comparison, compute_metrics, plot_radar_comparison
e_results = run_energy_comparison()
plot_energy_comparison(e_results)
m2_metrics = compute_metrics(e_results)
plot_radar_comparison(m2_metrics)
with open('results/m2_energy_metrics.json', 'w') as f:
    json.dump(m2_metrics, f, indent=2)
for name, m in m2_metrics.items():
    print(f"  ✓ {name}: sustain={m['sustain_time_min']:.0f} min, peak ATP={m['peak_ATP_mM']:.2f} mM")

# --- Module 3 ---
print("\n[3/6] Ion Concentration Optimization...")
from cfps_ion_optimization import plot_optimization_maps, bayesian_optimization_1d_demo, plot_bayesian_opt, grid_search_optimum
plot_optimization_maps()
X_obs, y_obs, best_vals = bayesian_optimization_1d_demo()
plot_bayesian_opt(X_obs, y_obs, best_vals)
m3_opt = grid_search_optimum()
with open('results/m3_ion_optimization.json', 'w') as f:
    json.dump(m3_opt, f, indent=2)
print(f"  ✓ Optimal: Mg={m3_opt['Mg_mM']:.1f}, K={m3_opt['K_mM']:.0f}, Spd={m3_opt['Spermidine_mM']:.1f}, yield={m3_opt['max_yield']:.4f}")

# --- Module 4 ---
print("\n[4/6] mRNA Stability & Ribosome Loading...")
from cfps_mrna_ribosome import plot_mrna_stability, plot_codon_sensitivity
m4_hl, m4_ribo = plot_mrna_stability()
plot_codon_sensitivity()
m4_results = {}
for name in m4_hl:
    m4_results[name] = {'halflife_min': m4_hl[name], **m4_ribo[name]}
with open('results/m4_mrna_ribosome.json', 'w') as f:
    json.dump(m4_results, f, indent=2)
for name, r in m4_results.items():
    print(f"  ✓ {name}: t½={r['halflife_min']:.1f} min, ribo/mRNA={r['ribosomes_per_mRNA']:.1f}")

# --- Module 5 ---
print("\n[5/6] Scale-Up Design...")
from cfps_scaleup import run_all_modes, plot_scaleup_comparison, plot_scaleup_volume_design, compute_scaleup_metrics
sol_b, sol_s, sol_c = run_all_modes()
plot_scaleup_comparison(sol_b, sol_s, sol_c)
plot_scaleup_volume_design()
m5_metrics = compute_scaleup_metrics(sol_b, sol_s, sol_c)
with open('results/m5_scaleup_metrics.json', 'w') as f:
    json.dump(m5_metrics, f, indent=2)
for name, m in m5_metrics.items():
    print(f"  ✓ {name}: peak={m['peak_protein_nM']:.1f} nM, productivity={m['volumetric_productivity_nM_per_hr']:.1f} nM/hr")

# --- Module 6 ---
print("\n[6/6] Membrane Protein / Nanodisc Case Study...")
from cfps_membrane_protein import run_nanodisc_titration, plot_nanodisc_case_study, compute_mp_metrics
mp_results = run_nanodisc_titration()
plot_nanodisc_case_study(mp_results)
m6_metrics = compute_mp_metrics(mp_results)
with open('results/m6_nanodisc_metrics.json', 'w') as f:
    json.dump(m6_metrics, f, indent=2)
for key, m in m6_metrics.items():
    print(f"  ✓ {key}: inserted={m['final_inserted_uM']:.2f} µM, efficiency={m['insertion_efficiency_pct']:.1f}%")

# --- Summary ---
all_results = {
    'module_1_txn_tln': m1_results,
    'module_2_energy': m2_metrics,
    'module_3_ion_opt': m3_opt,
    'module_4_mrna': m4_results,
    'module_5_scaleup': m5_metrics,
    'module_6_nanodisc': m6_metrics,
    'timestamp': datetime.now().isoformat(),
}
with open('results/all_results_summary.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\n" + "=" * 70)
print("  All modules completed successfully!")
print(f"  Results saved to: results/all_results_summary.json")
print(f"  Figures: figures/fig1-fig11")
print("=" * 70)
