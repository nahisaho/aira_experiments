#!/usr/bin/env python3
"""
Main simulation runner for Fischer-Tropsch Microkinetic Modeling
================================================================
Runs the complete FT case study and generates all outputs.
"""

import sys
import os
import json
import numpy as np
from datetime import datetime

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from microkinetic_framework.ft_synthesis import (
    run_ft_case_study, FT_TRANSITION_STATES, FT_SURFACE_SPECIES,
    FT_GAS_SPECIES, build_stoichiometric_matrices, ft_rate_expressions
)
from microkinetic_framework.rate_constants import calculate_tst_rate, arrhenius_parameters
from microkinetic_framework.adsorption import (
    langmuir_isotherm, temkin_isotherm, fractal_isotherm, AdsorptionParameters
)
from microkinetic_framework.lateral import (
    LateralInteractionParams, solve_coverage_self_consistent,
    mean_field_interaction_energy
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

CMAP = plt.cm.viridis


def save_json(data, filepath):
    """Save dict to JSON, handling numpy types."""
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2, cls=NpEncoder, ensure_ascii=False)


def plot_energy_diagram(results, filepath):
    """Plot reaction energy diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))

    intermediate_energies = [0.0, -1.30, -0.50, -0.60, -0.70, -1.10, -1.40, -2.00, -1.85, -2.35]
    labels = ['CO(g)+H2(g)', 'CO*', 'H*', 'HCO*', 'CH2O*', 'CH2*', 'CH3*', 'O*', 'OH*', 'Products']

    rc_list = results['rate_constants']
    ts_labels = list(rc_list.keys())

    # Plot intermediates
    x_pos = np.arange(len(labels)) * 2
    for i, (x, E, lbl) in enumerate(zip(x_pos, intermediate_energies, labels)):
        color = CMAP(i / len(labels))
        ax.plot([x - 0.4, x + 0.4], [E, E], '-', color=color, lw=3)
        ax.text(x, E - 0.08, lbl, ha='center', va='top', fontsize=7, rotation=30)

    # Plot TS barriers
    for i in range(len(intermediate_energies) - 1):
        E_reactant = intermediate_energies[i]
        if i < len(ts_labels):
            E_act = list(rc_list.values())[i]['E_act_forward']
        else:
            E_act = 0.5
        E_ts = E_reactant + E_act
        E_product = intermediate_energies[i + 1]

        x_r = x_pos[i]
        x_p = x_pos[i + 1]
        x_ts = (x_r + x_p) / 2

        # Spline-like path through TS
        x_path = np.linspace(x_r + 0.4, x_p - 0.4, 50)
        # Parabolic path
        t = (x_path - x_r - 0.4) / (x_p - 0.4 - x_r - 0.4)
        E_path = (1 - t)**2 * E_reactant + 2 * (1 - t) * t * E_ts + t**2 * E_product

        ax.plot(x_path, E_path, '--', color='gray', alpha=0.5, lw=1)
        ax.plot(x_ts, E_ts, 'v', color='red', markersize=6)

    ax.set_xlabel('Reaction Coordinate')
    ax.set_ylabel('Energy [eV]')
    ax.set_title('FT Synthesis Energy Diagram on Co(0001)')
    ax.set_xticks([])
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_rate_constants(results, filepath):
    """Plot rate constants bar chart."""
    rc = results['rate_constants']
    labels = list(rc.keys())
    k_fwd = [rc[l]['k_forward'] for l in labels]
    k_rev = [rc[l]['k_reverse'] for l in labels]

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(labels))
    w = 0.35

    colors_fwd = [CMAP(0.3)] * len(labels)
    colors_rev = [CMAP(0.7)] * len(labels)

    ax.bar(x - w/2, np.log10(np.maximum(k_fwd, 1e-30)), w, label='k_forward', color=colors_fwd)
    ax.bar(x + w/2, np.log10(np.maximum(k_rev, 1e-30)), w, label='k_reverse', color=colors_rev)

    ax.set_xlabel('Elementary Step')
    ax.set_ylabel('log10(k) [1/s]')
    ax.set_title('Rate Constants at T=500 K (TST + Wigner Tunneling)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_rds_analysis(results, filepath):
    """Plot degree of rate control."""
    rds = results['rds_analysis']
    labels = rds['step_labels']
    X_RC = rds['X_RC']

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [CMAP(0.8) if abs(x) > 0.1 else CMAP(0.3) for x in X_RC]
    bars = ax.bar(range(len(labels)), X_RC, color=colors, edgecolor='black', lw=0.5)

    ax.set_xlabel('Elementary Step')
    ax.set_ylabel('Degree of Rate Control (X_RC)')
    ax.set_title(f'Campbell\'s Degree of Rate Control — RDS: {rds["rds"]}')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.axhline(y=0, color='black', lw=0.5)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_coverages(results, filepath):
    """Plot steady-state surface coverages."""
    cov = results['steady_state_coverages']
    species = list(cov.keys())
    values = [cov[sp] for sp in species]

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = [CMAP(i / len(species)) for i in range(len(species))]
    ax.bar(range(len(species)), values, color=colors, edgecolor='black', lw=0.5)

    ax.set_xlabel('Surface Species')
    ax.set_ylabel('Coverage (θ)')
    ax.set_title('Steady-State Surface Coverages at T=500 K, P=20 bar')
    ax.set_xticks(range(len(species)))
    ax.set_xticklabels(species, rotation=45, ha='right')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-15)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_adsorption_isotherms(filepath):
    """Plot comparison of adsorption isotherm models."""
    params_CO = AdsorptionParameters('CO', delta_H_ads=-1.30, delta_S_ads=-1.5e-3)
    T = 500.0
    pressures = np.linspace(0.01, 30, 200)

    theta_lang = [langmuir_isotherm(P, T, params_CO) for P in pressures]
    theta_temkin = [temkin_isotherm(P, T, params_CO, alpha=0.5, delta_E=0.3) for P in pressures]
    theta_fractal = [fractal_isotherm(P, T, params_CO, D_f=2.5) for P in pressures]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(pressures, theta_lang, '-', color=CMAP(0.2), lw=2, label='Langmuir')
    ax.plot(pressures, theta_temkin, '--', color=CMAP(0.5), lw=2, label='Temkin (α=0.5)')
    ax.plot(pressures, theta_fractal, '-.', color=CMAP(0.8), lw=2, label='Fractal (D=2.5)')

    ax.set_xlabel('Pressure [bar]')
    ax.set_ylabel('Coverage (θ)')
    ax.set_title('CO Adsorption Isotherms on Co(0001) at T=500 K')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_arrhenius(results, filepath):
    """Plot Arrhenius analysis for C-O scission (key step)."""
    arrh = results['arrhenius']
    temps = np.array(arrh['temperatures'])
    k_fwd = np.array(arrh['k_forward'])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(1000 / temps, np.log10(np.maximum(k_fwd, 1e-30)), 'o-',
            color=CMAP(0.4), markersize=3, lw=1.5)

    ax.set_xlabel('1000/T [1/K]')
    ax.set_ylabel('log10(k_forward) [1/s]')
    ax.set_title(f'Arrhenius Plot — C-O Scission (Ea = {arrh["Ea_forward"]:.3f} eV)')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_reactor_profiles(results, filepath):
    """Plot PFR reactor profiles."""
    rr = results['reactor_result_obj']
    if rr.reactor_type != "PFR":
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Gas composition profiles
    ax1 = axes[0]
    for sp in FT_GAS_SPECIES:
        if sp in rr.gas_compositions:
            data = rr.gas_compositions[sp]
            ax1.plot(rr.positions * 1000, data, lw=2, label=sp)

    ax1.set_xlabel('Catalyst Weight [g]')
    ax1.set_ylabel('Mole Fraction')
    ax1.set_title('Gas Composition along PFR')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Surface coverages
    ax2 = axes[1]
    for sp in FT_SURFACE_SPECIES:
        if sp in rr.surface_coverages:
            data = rr.surface_coverages[sp]
            if np.max(data) > 1e-10:
                ax2.plot(rr.positions * 1000, data, lw=2, label=sp)

    ax2.set_xlabel('Catalyst Weight [g]')
    ax2.set_ylabel('Coverage (θ)')
    ax2.set_title('Surface Coverage along PFR')
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_lateral_interaction_effect(filepath):
    """Plot effect of lateral interactions on coverage."""
    T = 500.0
    P_CO_values = np.linspace(0.1, 15, 50)

    cov_no_lat = []
    cov_with_lat = []

    interactions = [
        LateralInteractionParams(('CO', 'CO'), epsilon_nn=-0.10, z_nn=6),
        LateralInteractionParams(('CO', 'H'), epsilon_nn=-0.02, z_nn=6),
        LateralInteractionParams(('H', 'H'), epsilon_nn=-0.01, z_nn=6),
    ]

    for P_CO in P_CO_values:
        # Without lateral interactions
        cov_clean = solve_coverage_self_consistent(
            {'CO': P_CO, 'H': 10.0},
            T,
            {'CO': -1.30, 'H': -0.50},
            {'CO': -1.5e-3, 'H': -0.8e-3},
            []  # No interactions
        )
        cov_no_lat.append(cov_clean.get('CO', 0))

        # With lateral interactions
        cov_lat = solve_coverage_self_consistent(
            {'CO': P_CO, 'H': 10.0},
            T,
            {'CO': -1.30, 'H': -0.50},
            {'CO': -1.5e-3, 'H': -0.8e-3},
            interactions
        )
        cov_with_lat.append(cov_lat.get('CO', 0))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(P_CO_values, cov_no_lat, '-', color=CMAP(0.3), lw=2, label='Without lateral interactions')
    ax.plot(P_CO_values, cov_with_lat, '--', color=CMAP(0.7), lw=2, label='With CO-CO repulsion')

    ax.set_xlabel('P_CO [bar]')
    ax.set_ylabel('CO Coverage (θ_CO)')
    ax.set_title('Effect of Lateral Interactions on CO Coverage (T=500 K)')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def plot_temperature_sensitivity(filepath):
    """Plot temperature sensitivity of key rates and coverages."""
    temps = np.linspace(400, 650, 30)
    tof_values = []
    conv_values = []

    for T in temps:
        try:
            res = run_ft_case_study(T=T, P_total=20.0, H2_CO_ratio=2.0,
                                     reactor_type="PFR", catalyst_mass=0.1)
            tof_values.append(res['reactor']['TOF'])
            conv_values.append(res['reactor']['conversion'])
        except Exception:
            tof_values.append(0)
            conv_values.append(0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(temps, tof_values, 'o-', color=CMAP(0.4), lw=2, markersize=4)
    ax1.set_xlabel('Temperature [K]')
    ax1.set_ylabel('TOF [1/s]')
    ax1.set_title('Temperature Dependence of TOF')
    ax1.grid(alpha=0.3)

    ax2.plot(temps, [c * 100 for c in conv_values], 's-', color=CMAP(0.7), lw=2, markersize=4)
    ax2.set_xlabel('Temperature [K]')
    ax2.set_ylabel('CO Conversion [%]')
    ax2.set_title('Temperature Dependence of CO Conversion')
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    print(f"  Saved: {filepath}")


def main():
    print("=" * 70)
    print("Microkinetic Modeling Framework — Fischer-Tropsch Case Study")
    print("=" * 70)
    print(f"Start time: {datetime.now().isoformat()}")
    print()

    # --- Run main simulation ---
    print("[1/9] Running FT synthesis case study (T=500K, P=20bar, H2/CO=2)...")
    results = run_ft_case_study(T=500.0, P_total=20.0, H2_CO_ratio=2.0,
                                 tunneling="wigner", reactor_type="PFR",
                                 catalyst_mass=0.1)

    # --- Print summary ---
    print("\n--- Rate Constants ---")
    for label, rc in results['rate_constants'].items():
        print(f"  {label:25s}: k_fwd = {rc['k_forward']:.3e}, "
              f"Ea = {rc['E_act_forward']:.2f} eV, κ = {rc['tunneling_correction']:.3f}")

    print("\n--- Adsorption (Langmuir) ---")
    for sp, cov in results['adsorption']['langmuir_coverages'].items():
        print(f"  θ({sp}) = {cov:.4f}")

    print("\n--- Steady-State Coverages ---")
    for sp, cov in results['steady_state_coverages'].items():
        print(f"  θ({sp}) = {cov:.4e}")

    print("\n--- RDS Analysis ---")
    rds = results['rds_analysis']
    for i, (label, xrc) in enumerate(zip(rds['step_labels'], rds['X_RC'])):
        marker = " ◄ RDS" if i == rds['rds_index'] else ""
        print(f"  {label:25s}: X_RC = {xrc:+.4f}{marker}")

    print(f"\n  Rate-Determining Step: {rds['rds']}")

    print("\n--- Energy Span Analysis ---")
    esp = results['energy_span']
    print(f"  Energy span: {esp['energy_span']:.3f} eV")
    print(f"  TDTS: {esp['TDTS_label']} (E = {esp['TDTS_energy']:.3f} eV)")
    print(f"  TDI energy: {esp['TDI_energy']:.3f} eV")
    print(f"  Estimated TOF: {esp['TOF_estimate']:.3e} s^-1")

    print("\n--- Reactor Results (PFR) ---")
    rx = results['reactor']
    print(f"  CO Conversion: {rx['conversion']*100:.2f}%")
    print(f"  TOF: {rx['TOF']:.3e} s^-1")
    print(f"  STY: {rx['STY']:.3e} mol/(kg_cat·s)")
    print(f"  Selectivities:")
    for sp, sel in rx['selectivities'].items():
        print(f"    {sp}: {sel*100:.2f}%")

    # --- Save numerical results ---
    print("\n[2/9] Saving numerical results...")
    results_export = {k: v for k, v in results.items() if k != 'reactor_result_obj'}
    save_json(results_export, 'results/ft_simulation_results.json')
    print("  Saved: results/ft_simulation_results.json")

    # --- Generate figures ---
    print("\n[3/9] Plotting energy diagram...")
    plot_energy_diagram(results, 'figures/fig1_energy_diagram.png')

    print("[4/9] Plotting rate constants...")
    plot_rate_constants(results, 'figures/fig2_rate_constants.png')

    print("[5/9] Plotting RDS analysis...")
    plot_rds_analysis(results, 'figures/fig3_rds_analysis.png')

    print("[6/9] Plotting surface coverages...")
    plot_coverages(results, 'figures/fig4_coverages.png')

    print("[7/9] Plotting adsorption isotherms...")
    plot_adsorption_isotherms('figures/fig5_adsorption_isotherms.png')

    print("[8/9] Plotting Arrhenius analysis...")
    plot_arrhenius(results, 'figures/fig6_arrhenius.png')

    print("[9/9] Plotting reactor profiles & sensitivity...")
    plot_reactor_profiles(results, 'figures/fig7_reactor_profiles.png')
    plot_lateral_interaction_effect('figures/fig8_lateral_interactions.png')
    plot_temperature_sensitivity('figures/fig9_temperature_sensitivity.png')

    print(f"\nCompleted: {datetime.now().isoformat()}")
    print("=" * 70)


if __name__ == '__main__':
    main()
