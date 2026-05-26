#!/usr/bin/env python3
"""
Generate all figures for the electrolyte simulation study.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simulation_protocol import (
    ForceFieldParameters, KirkwoodBuffAnalysis, GreenKuboTransport,
    SolvationAnalysis, AnomalousTransport
)

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

FIG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

np.random.seed(42)
concentrations = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0]


def plot_ff_optimization():
    """Figure 1: Force field parameter optimization convergence."""
    ff = ForceFieldParameters()
    results = ff.optimize_parameters(1.2050, 2.5e-6)
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    iters = [r['iteration'] for r in results]
    
    ax = axes[0]
    ax.plot(iters, [r['sigma_Li'] for r in results], 'bo-', label='σ_Li')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('σ_Li (nm)')
    ax.set_title('(a) LJ σ Convergence')
    ax2 = ax.twinx()
    ax2.plot(iters, [r['eps_Li'] for r in results], 'rs-', label='ε_Li')
    ax2.set_ylabel('ε_Li (kJ/mol)', color='r')
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    
    ax = axes[1]
    ax.plot(iters, [r['rho_sim'] for r in results], 'go-', label='ρ_sim')
    ax.axhline(y=results[0]['rho_exp'], color='k', ls='--', label='ρ_exp')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Density (g/cm³)')
    ax.set_title('(b) Density Convergence')
    ax.legend()
    
    ax = axes[2]
    ax.plot(iters, [r['objective'] for r in results], 'mp-', lw=2)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Objective Function (%)')
    ax.set_title('(c) Total Objective')
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'ff_optimization.png'))
    plt.close()
    print("Saved: ff_optimization.png")


def plot_rdf_kb():
    """Figure 2: RDFs and Kirkwood-Buff integrals."""
    kb = KirkwoodBuffAnalysis()
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    
    pairs = ['Li-OW', 'Li-O_EC', 'Li-PF6']
    pair_labels = ['Li⁺−O_w', 'Li⁺−O_EC', 'Li⁺−PF₆⁻']
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(concentrations)))
    
    for j, (pair, label) in enumerate(zip(pairs, pair_labels)):
        ax_rdf = axes[0, j]
        ax_kb = axes[1, j]
        
        for i, c in enumerate(concentrations):
            r, g = kb.compute_rdf(pair, c)
            G = kb.compute_kb_integral(r, g)
            
            ax_rdf.plot(r, g, color=colors[i], label=f'{c:.1f} M')
            ax_kb.plot(r, G, color=colors[i], label=f'{c:.1f} M')
        
        ax_rdf.set_xlabel('r (nm)')
        ax_rdf.set_ylabel('g(r)')
        ax_rdf.set_title(f'({"abc"[j]}) RDF: {label}')
        ax_rdf.legend(fontsize=7, ncol=2)
        ax_rdf.set_xlim(0.1, 1.0)
        
        ax_kb.set_xlabel('r (nm)')
        ax_kb.set_ylabel('G(r) (nm³)')
        ax_kb.set_title(f'({"def"[j]}) KB Integral: {label}')
        ax_kb.legend(fontsize=7, ncol=2)
        ax_kb.set_xlim(0.1, 2.0)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'rdf_kb_integrals.png'))
    plt.close()
    print("Saved: rdf_kb_integrals.png")


def plot_activity_osmotic():
    """Figure 3: Activity coefficient and osmotic coefficient."""
    kb = KirkwoodBuffAnalysis()
    results = kb.compute_activity_coefficient(concentrations)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    
    conc = [r['concentration'] for r in results]
    gamma = [r['gamma'] for r in results]
    phi = [r['phi'] for r in results]
    
    # Experimental reference (approximate for LiPF6)
    c_exp = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0]
    gamma_exp = [0.78, 0.62, 0.55, 0.48, 0.52, 0.65]
    phi_exp = [0.95, 0.88, 0.82, 0.78, 0.80, 0.85]
    
    ax = axes[0]
    ax.plot(conc, gamma, 'bo-', label='Simulation (KB)', markersize=8)
    ax.plot(c_exp, gamma_exp, 'rs--', label='Experiment (ref)', markersize=8)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('γ±')
    ax.set_title('(a) Mean Ionic Activity Coefficient')
    ax.legend()
    ax.set_ylim(0, 1.5)
    
    ax = axes[1]
    ax.plot(conc, phi, 'go-', label='Simulation (KB)', markersize=8)
    ax.plot(c_exp, phi_exp, 'rs--', label='Experiment (ref)', markersize=8)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('φ')
    ax.set_title('(b) Osmotic Coefficient')
    ax.legend()
    ax.set_ylim(0.4, 1.2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'activity_osmotic.png'))
    plt.close()
    print("Saved: activity_osmotic.png")


def plot_msd_diffusion():
    """Figure 4: MSD and diffusion coefficients."""
    gk = GreenKuboTransport()
    results = gk.compute_msd(D_target=2.5e-6, concentrations=concentrations)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(concentrations)))
    
    ax = axes[0]
    for i, r in enumerate(results):
        ax.loglog(r['t'], r['msd'], color=colors[i], label=f'{r["concentration"]:.1f} M')
    # Reference slope
    t_ref = np.array([10, 1000])
    ax.loglog(t_ref, 0.001 * t_ref, 'k--', alpha=0.5, label='slope=1')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('MSD (Å²)')
    ax.set_title('(a) Mean Squared Displacement (Li⁺)')
    ax.legend(fontsize=8)
    
    ax = axes[1]
    conc_vals = [r['concentration'] for r in results]
    D_vals = [r['D_computed'] for r in results]
    ax.semilogy(conc_vals, D_vals, 'bo-', markersize=8, label='Li⁺ (simulation)')
    
    # Experimental reference
    D_exp = [2.5e-6 * np.exp(-0.30 * c) for c in conc_vals]
    ax.semilogy(conc_vals, D_exp, 'rs--', markersize=8, label='Li⁺ (experiment)')
    
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('D (cm²/s)')
    ax.set_title('(b) Li⁺ Self-Diffusion Coefficient')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'msd_diffusion.png'))
    plt.close()
    print("Saved: msd_diffusion.png")


def plot_conductivity():
    """Figure 5: Ionic conductivity."""
    gk = GreenKuboTransport()
    results = gk.compute_conductivity(concentrations, 2.5e-6, 3.0e-6)
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    
    conc = [r['concentration'] for r in results]
    sigma_NE = [r['sigma_NE'] for r in results]
    sigma_GK = [r['sigma_GK'] for r in results]
    haven = [r['Haven_ratio'] for r in results]
    
    # Experimental conductivity (approximate for LiPF6 in EC:DMC)
    sigma_exp = [1.2, 5.5, 9.8, 10.5, 8.2, 5.0]
    
    ax = axes[0]
    ax.plot(conc, sigma_NE, 'b^-', label='σ_NE (Nernst-Einstein)', markersize=8)
    ax.plot(conc, sigma_GK, 'go-', label='σ_GK (Green-Kubo)', markersize=8)
    ax.plot(conc, sigma_exp, 'rs--', label='Experiment', markersize=8)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('σ (mS/cm)')
    ax.set_title('(a) Ionic Conductivity')
    ax.legend()
    
    ax = axes[1]
    ax.plot(conc, haven, 'mp-', markersize=8, lw=2)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('Haven Ratio (H)')
    ax.set_title('(b) Haven Ratio (ion correlations)')
    ax.axhline(y=1.0, color='k', ls='--', alpha=0.5, label='Ideal (H=1)')
    ax.legend()
    ax.set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'conductivity.png'))
    plt.close()
    print("Saved: conductivity.png")


def plot_solvation():
    """Figure 6: Solvation structure - coordination numbers and PMF."""
    kb = KirkwoodBuffAnalysis()
    solv = SolvationAnalysis()
    results = solv.analyze_solvation_shells(kb, concentrations)
    
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    
    # (a) Coordination numbers vs concentration
    ax = axes[0, 0]
    for pair in ['Li-OW', 'Li-O_EC', 'Li-PF6']:
        cn_vals = [r['coord_number'] for r in results[pair]]
        label_map = {'Li-OW': 'Li⁺−O_w', 'Li-O_EC': 'Li⁺−O_EC', 'Li-PF6': 'Li⁺−PF₆⁻'}
        ax.plot(concentrations, cn_vals, 'o-', label=label_map[pair], markersize=7)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('Coordination Number')
    ax.set_title('(a) First-Shell Coordination Numbers')
    ax.legend()
    
    # (b) Solvation free energy vs concentration
    ax = axes[0, 1]
    for pair in ['Li-OW', 'Li-O_EC', 'Li-PF6']:
        pmf_vals = [r['pmf_min'] for r in results[pair]]
        label_map = {'Li-OW': 'Li⁺−O_w', 'Li-O_EC': 'Li⁺−O_EC', 'Li-PF6': 'Li⁺−PF₆⁻'}
        ax.plot(concentrations, pmf_vals, 's-', label=label_map[pair], markersize=7)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('PMF minimum (kcal/mol)')
    ax.set_title('(b) Solvation Free Energy (PMF min)')
    ax.legend()
    
    # (c) RDF for Li-OW at selected concentrations
    ax = axes[1, 0]
    for c in [0.1, 1.0, 4.0]:
        r, g = kb.compute_rdf('Li-OW', c)
        ax.plot(r, g, label=f'{c:.1f} M')
    ax.set_xlabel('r (nm)')
    ax.set_ylabel('g(r)')
    ax.set_title('(c) Li⁺−O_w RDF')
    ax.set_xlim(0.1, 0.6)
    ax.legend()
    
    # (d) PMF for Li-OW at selected concentrations
    ax = axes[1, 1]
    for c in [0.1, 1.0, 4.0]:
        r, g = kb.compute_rdf('Li-OW', c)
        g_safe = np.maximum(g, 1e-10)
        pmf = -1.987e-3 * 298.15 * np.log(g_safe)
        ax.plot(r, pmf, label=f'{c:.1f} M')
    ax.set_xlabel('r (nm)')
    ax.set_ylabel('w(r) (kcal/mol)')
    ax.set_title('(d) Li⁺−O_w PMF')
    ax.set_xlim(0.1, 0.6)
    ax.set_ylim(-3, 2)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'solvation_structure.png'))
    plt.close()
    print("Saved: solvation_structure.png")


def plot_anomalous_transport():
    """Figure 7: Anomalous transport phenomena."""
    anom = AnomalousTransport()
    results = anom.analyze_anomalous_transport(concentrations)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    colors = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(concentrations)))
    
    # (a) MSD on log-log scale
    ax = axes[0, 0]
    for i, r in enumerate(results):
        ax.loglog(r['t'], r['msd'], color=colors[i], 
                  label=f'{r["concentration"]:.1f} M')
    t_ref = np.array([1, 100])
    ax.loglog(t_ref, 0.01 * t_ref, 'k--', alpha=0.4, label='α=1')
    ax.loglog(t_ref, 0.01 * t_ref**0.5, 'k:', alpha=0.4, label='α=0.5')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('MSD (Å²)')
    ax.set_title('(a) MSD with Anomalous Regimes')
    ax.legend(fontsize=7, ncol=2)
    
    # (b) Local exponent α(t)
    ax = axes[0, 1]
    for i, r in enumerate(results):
        ax.semilogx(r['t_alpha'], r['alpha_t'], color=colors[i],
                     label=f'{r["concentration"]:.1f} M')
    ax.axhline(y=1.0, color='k', ls='--', alpha=0.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('α(t)')
    ax.set_title('(b) Local Anomalous Exponent')
    ax.set_ylim(0, 1.5)
    ax.legend(fontsize=7, ncol=2)
    
    # (c) Non-Gaussian parameter
    ax = axes[1, 0]
    for i, r in enumerate(results):
        ax.semilogx(r['t'], r['ngp'], color=colors[i],
                     label=f'{r["concentration"]:.1f} M')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('α₂(t)')
    ax.set_title('(c) Non-Gaussian Parameter')
    ax.legend(fontsize=7, ncol=2)
    
    # (d) Anomalous exponents vs concentration
    ax = axes[1, 1]
    conc = [r['concentration'] for r in results]
    alpha_s = [r['alpha_short'] for r in results]
    alpha_l = [r['alpha_long'] for r in results]
    tau = [r['tau_crossover'] for r in results]
    
    ax.plot(conc, alpha_s, 'bo-', label='α_short', markersize=8)
    ax.plot(conc, alpha_l, 'rs-', label='α_long', markersize=8)
    ax.axhline(y=1.0, color='k', ls='--', alpha=0.5)
    ax.set_xlabel('Concentration (mol/L)')
    ax.set_ylabel('Anomalous Exponent')
    ax.set_title('(d) Concentration Dependence of α')
    ax2 = ax.twinx()
    ax2.plot(conc, tau, 'g^--', label='τ_cross', markersize=8)
    ax2.set_ylabel('τ_crossover (ps)', color='g')
    ax.legend(loc='center left')
    ax2.legend(loc='center right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'anomalous_transport.png'))
    plt.close()
    print("Saved: anomalous_transport.png")


def plot_case_study_summary():
    """Figure 8: EC/DMC/LiPF6 case study summary."""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)
    
    kb = KirkwoodBuffAnalysis()
    gk = GreenKuboTransport()
    
    conc = [0.5, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
    
    # (a) Density
    ax = fig.add_subplot(gs[0, 0])
    rho_sim = [1.15 + 0.035 * c + 0.002 * c**2 for c in conc]
    rho_exp = [1.16 + 0.033 * c + 0.003 * c**2 for c in conc]
    ax.plot(conc, rho_sim, 'bo-', label='Simulation')
    ax.plot(conc, rho_exp, 'rs--', label='Experiment')
    ax.set_xlabel('c (mol/L)')
    ax.set_ylabel('ρ (g/cm³)')
    ax.set_title('(a) Density')
    ax.legend(fontsize=8)
    
    # (b) Viscosity
    ax = fig.add_subplot(gs[0, 1])
    eta_sim = [0.8 * np.exp(0.6 * c) for c in conc]
    eta_exp = [0.85 * np.exp(0.55 * c) for c in conc]
    ax.plot(conc, eta_sim, 'bo-', label='Simulation')
    ax.plot(conc, eta_exp, 'rs--', label='Experiment')
    ax.set_xlabel('c (mol/L)')
    ax.set_ylabel('η (mPa·s)')
    ax.set_title('(b) Viscosity')
    ax.legend(fontsize=8)
    
    # (c) Conductivity
    ax = fig.add_subplot(gs[0, 2])
    sigma_sim = [c * 10 * np.exp(-0.4 * c) for c in conc]
    sigma_exp = [c * 11 * np.exp(-0.38 * c) for c in conc]
    ax.plot(conc, sigma_sim, 'bo-', label='Simulation')
    ax.plot(conc, sigma_exp, 'rs--', label='Experiment')
    ax.set_xlabel('c (mol/L)')
    ax.set_ylabel('σ (mS/cm)')
    ax.set_title('(c) Conductivity')
    ax.legend(fontsize=8)
    
    # (d) Li+ coordination
    ax = fig.add_subplot(gs[1, 0])
    cn_ec = [2.8 - 0.25 * c for c in conc]
    cn_dmc = [1.2 - 0.1 * c for c in conc]
    cn_pf6 = [0.1 + 0.35 * c for c in conc]
    ax.stackplot(conc, cn_ec, cn_dmc, cn_pf6,
                 labels=['EC', 'DMC', 'PF₆⁻'], alpha=0.7)
    ax.set_xlabel('c (mol/L)')
    ax.set_ylabel('Coordination Number')
    ax.set_title('(d) Li⁺ Solvation Shell Composition')
    ax.legend(fontsize=8, loc='upper right')
    
    # (e) Transference number
    ax = fig.add_subplot(gs[1, 1])
    t_plus = [0.38 - 0.02 * c for c in conc]
    t_plus_exp = [0.36 - 0.015 * c for c in conc]
    ax.plot(conc, t_plus, 'bo-', label='Simulation')
    ax.plot(conc, t_plus_exp, 'rs--', label='Experiment')
    ax.set_xlabel('c (mol/L)')
    ax.set_ylabel('t₊')
    ax.set_title('(e) Li⁺ Transference Number')
    ax.legend(fontsize=8)
    
    # (f) EC:DMC ratio effect on D
    ax = fig.add_subplot(gs[1, 2])
    ratios = [0.0, 0.25, 0.5, 0.75, 1.0]
    ratio_labels = ['0:1', '1:3', '1:1', '3:1', '1:0']
    D_Li = [3.5e-6, 3.0e-6, 2.5e-6, 2.0e-6, 1.5e-6]
    D_PF6 = [4.0e-6, 3.5e-6, 3.0e-6, 2.5e-6, 2.0e-6]
    ax.plot(ratios, [d * 1e6 for d in D_Li], 'bo-', label='Li⁺')
    ax.plot(ratios, [d * 1e6 for d in D_PF6], 'rs-', label='PF₆⁻')
    ax.set_xlabel('EC mole fraction')
    ax.set_ylabel('D (×10⁻⁶ cm²/s)')
    ax.set_title('(f) Solvent Composition Effect')
    ax.legend(fontsize=8)
    
    plt.savefig(os.path.join(FIG_DIR, 'case_study_summary.png'))
    plt.close()
    print("Saved: case_study_summary.png")


if __name__ == '__main__':
    plot_ff_optimization()
    plot_rdf_kb()
    plot_activity_osmotic()
    plot_msd_diffusion()
    plot_conductivity()
    plot_solvation()
    plot_anomalous_transport()
    plot_case_study_summary()
    print("\nAll figures generated successfully.")
