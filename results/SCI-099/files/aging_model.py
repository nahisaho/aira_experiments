"""
Integrated Mathematical Model of Aging Mechanisms
==================================================
ODE-based framework integrating:
- Hallmarks of Aging interaction network
- Reliability Theory + Antagonistic Pleiotropy
- Senolytics effect prediction
- Caloric Restriction / Rapamycin / NAD+ precursor modeling
- Cross-species lifespan scaling
- Intervention combination optimization
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import json
import os
from datetime import datetime
from itertools import product

# ──────────────────────────────────────────────
# 1. State Variables (Hallmarks of Aging)
# ──────────────────────────────────────────────
STATE_NAMES = [
    'T',    # 0: Telomere length (normalized, 1→0)
    'E',    # 1: Epigenetic integrity (1→0)
    'M',    # 2: Mitochondrial function (1→0)
    'P',    # 3: Proteostasis capacity (1→0)
    'N',    # 4: Nutrient sensing dysregulation (0→1, damage accumulates)
    'S',    # 5: Senescent cell fraction (0→1)
    'I',    # 6: Chronic inflammation (0→1)
    'SC',   # 7: Stem cell function (1→0)
    'G',    # 8: Genomic stability (1→0)
    'IC',   # 9: Intercellular communication quality (1→0)
    'D',    # 10: Overall damage (Reliability Theory aggregate)
    'R',    # 11: Reproductive fitness proxy (Antagonistic Pleiotropy)
]

N_STATES = len(STATE_NAMES)

# ──────────────────────────────────────────────
# 2. Default Parameters
# ──────────────────────────────────────────────
def default_params():
    """Return default kinetic parameters for human-like aging."""
    return dict(
        # Telomere shortening
        k_T=0.012, alpha_TG=0.005,
        # Epigenetic drift
        k_E=0.008, alpha_EM=0.004, alpha_ES=0.003,
        # Mitochondrial decline
        k_M=0.010, alpha_MG=0.006, alpha_ME=0.003,
        # Proteostasis decline
        k_P=0.007, alpha_PM=0.005, alpha_PE=0.002,
        # Nutrient sensing dysregulation
        k_N=0.006, alpha_NM=0.004,
        # Senescent cell accumulation
        k_S=0.005, alpha_ST=0.008, alpha_SG=0.004, alpha_SM=0.003,
        gamma_S=0.002,  # immune clearance rate
        # Inflammation (SASP-driven)
        k_I=0.003, alpha_IS=0.010, alpha_IM=0.004,
        gamma_I=0.005,  # resolution rate
        # Stem cell exhaustion
        k_SC=0.006, alpha_SCT=0.005, alpha_SCI=0.004, alpha_SCE=0.003,
        # Genomic instability
        k_G=0.009, alpha_GT=0.004, alpha_GM=0.003,
        # Intercellular communication
        k_IC=0.005, alpha_ICI=0.006, alpha_ICS=0.003,
        # Reliability theory aggregate damage
        k_D=0.15, beta_D=2.5,  # Gompertz-like acceleration
        # Antagonistic pleiotropy
        k_R_early=0.05, k_R_late=0.02, t_repro=30.0,
        # Mortality hazard (calibrated for ~80y median human lifespan)
        h0=0.004, h_exp=2.8,
    )


# ──────────────────────────────────────────────
# 3. Intervention Parameters
# ──────────────────────────────────────────────
def intervention_params(cr=0.0, rapa=0.0, nad=0.0, senolytic=0.0):
    """
    cr: caloric restriction intensity [0,1]
    rapa: rapamycin dose level [0,1]
    nad: NAD+ precursor level [0,1]
    senolytic: senolytic clearance boost [0,1]
    """
    return dict(
        cr=cr,
        rapa=rapa,
        nad=nad,
        senolytic=senolytic,
    )


# ──────────────────────────────────────────────
# 4. ODE System
# ──────────────────────────────────────────────
def aging_odes(t, y, params, interv):
    T, E, M, P, N, S, I, SC, G, IC, D, R = y
    p = params
    iv = interv

    # Clamp values
    T = np.clip(T, 0, 1)
    E = np.clip(E, 0, 1)
    M = np.clip(M, 0, 1)
    P = np.clip(P, 0, 1)
    N = np.clip(N, 0, 1)
    S = np.clip(S, 0, 1)
    I = np.clip(I, 0, 1)
    SC = np.clip(SC, 0, 1)
    G = np.clip(G, 0, 1)
    IC = np.clip(IC, 0, 1)
    D = np.clip(D, 0, 10)
    R = np.clip(R, 0, 1)

    # Intervention modifiers
    cr_factor = 1.0 - 0.3 * iv['cr']       # CR slows damage by up to 30%
    rapa_factor = 1.0 - 0.25 * iv['rapa']   # Rapamycin: mTOR inhibition
    nad_factor = 1.0 + 0.2 * iv['nad']      # NAD+ boosts mitochondrial function
    seno_clear = iv['senolytic'] * 0.05      # Senolytic clearance boost

    # dT/dt: Telomere shortening (accelerated by genomic instability)
    dT = -(p['k_T'] + p['alpha_TG'] * (1 - G)) * T * cr_factor

    # dE/dt: Epigenetic drift (accelerated by mitochondrial dysfunction & senescence)
    dE = -(p['k_E'] + p['alpha_EM'] * (1 - M) + p['alpha_ES'] * S) * E * cr_factor

    # dM/dt: Mitochondrial function decline (modulated by NAD+, genomic instability, epigenetic)
    dM = -(p['k_M'] + p['alpha_MG'] * (1 - G) + p['alpha_ME'] * (1 - E)) * M / nad_factor

    # dP/dt: Proteostasis decline
    dP = -(p['k_P'] + p['alpha_PM'] * (1 - M) + p['alpha_PE'] * (1 - E)) * P * rapa_factor

    # dN/dt: Nutrient sensing dysregulation (accumulates)
    dN = (p['k_N'] + p['alpha_NM'] * (1 - M)) * (1 - N) * rapa_factor * cr_factor

    # dS/dt: Senescent cell accumulation
    senescence_prod = (p['k_S'] + p['alpha_ST'] * (1 - T) + p['alpha_SG'] * (1 - G)
                       + p['alpha_SM'] * (1 - M)) * (1 - S)
    senescence_clear = (p['gamma_S'] + seno_clear) * S
    dS = senescence_prod - senescence_clear

    # dI/dt: Chronic inflammation (SASP-driven)
    dI = (p['k_I'] + p['alpha_IS'] * S + p['alpha_IM'] * (1 - M)) * (1 - I) - p['gamma_I'] * I

    # dSC/dt: Stem cell exhaustion
    dSC = -(p['k_SC'] + p['alpha_SCT'] * (1 - T) + p['alpha_SCI'] * I
            + p['alpha_SCE'] * (1 - E)) * SC * cr_factor

    # dG/dt: Genomic instability
    dG = -(p['k_G'] + p['alpha_GT'] * (1 - T) + p['alpha_GM'] * (1 - M)) * G

    # dIC/dt: Intercellular communication decline
    dIC = -(p['k_IC'] + p['alpha_ICI'] * I + p['alpha_ICS'] * S) * IC

    # dD/dt: Aggregate damage (Reliability Theory — Gompertz-like acceleration)
    damage_sources = (1 - T) + (1 - E) + (1 - M) + (1 - P) + N + S + I + (1 - SC) + (1 - G) + (1 - IC)
    dD = p['k_D'] * (damage_sources / 10.0) ** p['beta_D']

    # dR/dt: Reproductive fitness (Antagonistic Pleiotropy)
    if t < p['t_repro']:
        dR = p['k_R_early'] * (1 - R) - p['k_R_late'] * D * R
    else:
        dR = -p['k_R_late'] * D * R * (1 + 0.5 * N)

    return [dT, dE, dM, dP, dN, dS, dI, dSC, dG, dIC, dD, dR]


# ──────────────────────────────────────────────
# 5. Simulation Runner
# ──────────────────────────────────────────────
def run_simulation(params=None, interv=None, t_span=(0, 120), y0=None, t_eval=None):
    if params is None:
        params = default_params()
    if interv is None:
        interv = intervention_params()
    if y0 is None:
        y0 = [1.0, 1.0, 1.0, 1.0, 0.0, 0.01, 0.02, 1.0, 1.0, 1.0, 0.0, 0.0]
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 1000)

    sol = solve_ivp(
        aging_odes, t_span, y0,
        args=(params, interv),
        t_eval=t_eval,
        method='RK45',
        max_step=0.5,
        rtol=1e-8, atol=1e-10,
    )
    return sol


def compute_mortality_hazard(sol, params):
    """Compute time-dependent mortality hazard h(t) = h0 * exp(h_exp * D(t))."""
    D = sol.y[10]
    exponent = np.clip(params['h_exp'] * D, 0, 500)
    h = params['h0'] * np.exp(exponent)
    return h


def compute_survival(sol, params):
    """Compute survival curve S(t) = exp(-∫h(τ)dτ)."""
    h = compute_mortality_hazard(sol, params)
    t = sol.t
    cumulative_hazard = np.cumsum(h[:-1] * np.diff(t))
    cumulative_hazard = np.insert(cumulative_hazard, 0, 0)
    survival = np.exp(-cumulative_hazard)
    return survival


def compute_median_lifespan(sol, params):
    survival = compute_survival(sol, params)
    idx = np.searchsorted(-survival, -0.5)
    if idx >= len(sol.t):
        return sol.t[-1]
    return sol.t[idx]


# ──────────────────────────────────────────────
# 6. Cross-Species Scaling
# ──────────────────────────────────────────────
def species_params(body_mass_kg, metabolic_rate_factor=1.0, dna_repair_factor=1.0):
    """
    Scale aging parameters by body mass and DNA repair capacity.
    Kleiber's law: metabolic rate ∝ M^0.75
    Lifespan ∝ M^0.25 (approximate)
    """
    p = default_params()
    # Human reference: ~70 kg
    mass_ratio = body_mass_kg / 70.0
    scaling = mass_ratio ** (-0.25) * metabolic_rate_factor

    # Faster metabolism → faster aging clocks
    for key in ['k_T', 'k_E', 'k_M', 'k_P', 'k_N', 'k_S', 'k_G', 'k_SC', 'k_IC', 'k_D']:
        p[key] *= scaling

    # Better DNA repair → slower genomic instability and telomere attrition
    p['k_G'] /= dna_repair_factor
    p['k_T'] /= (dna_repair_factor ** 0.5)
    p['alpha_SG'] /= dna_repair_factor

    return p


SPECIES_DB = {
    'Mouse':           dict(mass=0.03,  metab=7.0, repair=0.5),
    'Rat':             dict(mass=0.3,   metab=5.0, repair=0.6),
    'Dog':             dict(mass=20.0,  metab=2.0, repair=0.8),
    'Human':           dict(mass=70.0,  metab=1.0, repair=1.0),
    'Elephant':        dict(mass=5000,  metab=0.5, repair=1.3),
    'Bowhead Whale':   dict(mass=80000, metab=0.3, repair=2.0),
    'Naked Mole Rat':  dict(mass=0.035, metab=1.5, repair=3.5),
    'Greenland Shark': dict(mass=1000,  metab=0.2, repair=1.5),
}


# ──────────────────────────────────────────────
# 7. Visualization
# ──────────────────────────────────────────────
def setup_plot_style():
    sns.set_theme(style='whitegrid', font_scale=1.1)
    plt.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })


def plot_hallmarks_trajectory(sol, title='Hallmarks of Aging Trajectories', fname='figures/fig1_hallmarks.png'):
    setup_plot_style()
    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    axes = axes.flatten()
    colors = sns.color_palette('viridis', N_STATES)
    for i, (name, color) in enumerate(zip(STATE_NAMES, colors)):
        ax = axes[i]
        ax.plot(sol.t, sol.y[i], color=color, linewidth=2)
        ax.set_title(name, fontsize=12, fontweight='bold')
        ax.set_xlabel('Age (years)')
        ax.set_ylabel('Level')
        ax.set_ylim(-0.05, max(1.05, np.max(sol.y[i]) * 1.1))
    fig.suptitle(title, fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname


def plot_interaction_network(fname='figures/fig2_network.png'):
    """Visualize the hallmark interaction network as a heatmap."""
    setup_plot_style()
    labels = ['Telomere', 'Epigenetic', 'Mitochond.', 'Proteostasis',
              'Nutr.Sens.', 'Senescence', 'Inflamm.', 'Stem Cell',
              'Genomic', 'Intercomm.']
    n = len(labels)
    # Interaction matrix (positive = promotes damage, negative = protective)
    W = np.zeros((n, n))
    # Telomere → Senescence, Stem Cell, Genomic
    W[0, 5] = 0.8; W[0, 7] = 0.5; W[0, 8] = 0.4
    # Epigenetic → Mitochondria, Proteostasis, Stem Cell
    W[1, 2] = 0.3; W[1, 3] = 0.2; W[1, 7] = 0.3
    # Mitochondria → Epigenetic, Proteostasis, Senescence, Inflammation, Genomic
    W[2, 1] = 0.4; W[2, 3] = 0.5; W[2, 5] = 0.3; W[2, 6] = 0.4; W[2, 8] = 0.3
    # Proteostasis → (downstream effects relatively contained)
    W[3, 5] = 0.2
    # Nutrient sensing → Mitochondria, Proteostasis
    W[4, 2] = 0.4; W[4, 3] = 0.3
    # Senescence → Inflammation, Stem Cell, Intercomm.
    W[5, 6] = 1.0; W[5, 7] = 0.4; W[5, 9] = 0.3
    # Inflammation → Senescence, Stem Cell, Genomic, Intercomm.
    W[6, 5] = 0.3; W[6, 7] = 0.4; W[6, 8] = 0.2; W[6, 9] = 0.6
    # Genomic → Telomere, Senescence, Mitochondria
    W[8, 0] = 0.4; W[8, 5] = 0.4; W[8, 2] = 0.3

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(W, cmap='YlOrRd', vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticks(range(n)); ax.set_yticklabels(labels)
    ax.set_xlabel('Target Hallmark (damaged by)')
    ax.set_ylabel('Source Hallmark (damages)')
    ax.set_title('Hallmarks of Aging Interaction Network', fontweight='bold')
    for i in range(n):
        for j in range(n):
            if W[i, j] > 0:
                ax.text(j, i, f'{W[i,j]:.1f}', ha='center', va='center',
                        color='white' if W[i, j] > 0.5 else 'black', fontsize=9)
    plt.colorbar(im, ax=ax, label='Interaction Strength')
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname


def plot_intervention_comparison(results_dict, params, fname='figures/fig3_interventions.png'):
    """Compare survival curves under different interventions."""
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    colors = sns.color_palette('Set2', len(results_dict))

    for (name, sol), color in zip(results_dict.items(), colors):
        survival = compute_survival(sol, params)
        median = compute_median_lifespan(sol, params)
        ax1.plot(sol.t, survival, label=f'{name} (med={median:.1f}y)', color=color, linewidth=2)

        # Plot damage accumulation
        ax2.plot(sol.t, sol.y[10], label=name, color=color, linewidth=2)

    ax1.set_xlabel('Age (years)'); ax1.set_ylabel('Survival Probability')
    ax1.set_title('Survival Curves Under Interventions', fontweight='bold')
    ax1.legend(fontsize=9); ax1.set_ylim(0, 1.05)

    ax2.set_xlabel('Age (years)'); ax2.set_ylabel('Aggregate Damage (D)')
    ax2.set_title('Damage Accumulation', fontweight='bold')
    ax2.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname


def plot_senolytic_timing(fname='figures/fig4_senolytics.png'):
    """Simulate senolytics started at different ages."""
    setup_plot_style()
    params = default_params()
    start_ages = [30, 40, 50, 60, 70]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    colors = sns.color_palette('viridis', len(start_ages) + 1)

    # Baseline
    sol_base = run_simulation(params=params)
    surv_base = compute_survival(sol_base, params)
    ax1.plot(sol_base.t, surv_base, '--', color='gray', linewidth=2, label='No Treatment')
    ax2.plot(sol_base.t, sol_base.y[5], '--', color='gray', linewidth=2, label='No Treatment')

    medians = {'No Treatment': compute_median_lifespan(sol_base, params)}

    for i, start in enumerate(start_ages):
        # Phase 1: no treatment
        t1 = np.linspace(0, start, 500)
        sol1 = run_simulation(params=params, t_span=(0, start), t_eval=t1)
        y_start = sol1.y[:, -1]

        # Phase 2: with senolytics
        interv = intervention_params(senolytic=0.8)
        t2 = np.linspace(start, 120, 500)
        sol2 = run_simulation(params=params, interv=interv, t_span=(start, 120), y0=y_start.tolist(), t_eval=t2)

        # Stitch together
        t_full = np.concatenate([sol1.t, sol2.t])
        y_full = np.concatenate([sol1.y, sol2.y], axis=1)

        class StitchedSol:
            pass
        sol_s = StitchedSol()
        sol_s.t = t_full
        sol_s.y = y_full

        surv = compute_survival(sol_s, params)
        median = compute_median_lifespan(sol_s, params)
        medians[f'Start@{start}'] = median

        ax1.plot(t_full, surv, color=colors[i], linewidth=2, label=f'Senolytic @{start}y (med={median:.1f})')
        ax2.plot(t_full, y_full[5], color=colors[i], linewidth=2, label=f'Start @{start}y')

    ax1.set_xlabel('Age (years)'); ax1.set_ylabel('Survival')
    ax1.set_title('Senolytic Timing: Survival Curves', fontweight='bold')
    ax1.legend(fontsize=8); ax1.set_ylim(0, 1.05)

    ax2.set_xlabel('Age (years)'); ax2.set_ylabel('Senescent Cell Fraction')
    ax2.set_title('Senescent Cell Dynamics', fontweight='bold')
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname, medians


def plot_species_comparison(fname='figures/fig5_species.png'):
    """Cross-species lifespan comparison."""
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    colors = sns.color_palette('Set2', len(SPECIES_DB))

    species_results = {}
    for (name, spec), color in zip(SPECIES_DB.items(), colors):
        p = species_params(spec['mass'], spec['metab'], spec['repair'])
        sol = run_simulation(params=p, t_span=(0, 300))
        surv = compute_survival(sol, p)
        median = compute_median_lifespan(sol, p)
        species_results[name] = dict(mass=spec['mass'], median=median,
                                     repair=spec['repair'], metab=spec['metab'])
        ax1.plot(sol.t, surv, color=color, linewidth=2, label=f'{name} ({median:.0f}y)')

    ax1.set_xlabel('Age (years)'); ax1.set_ylabel('Survival')
    ax1.set_title('Cross-Species Survival Curves', fontweight='bold')
    ax1.legend(fontsize=8)
    ax1.set_xlim(0, 300)

    # Scaling plot: mass vs lifespan
    masses = [v['mass'] for v in species_results.values()]
    lifespans = [v['median'] for v in species_results.values()]
    names = list(species_results.keys())
    ax2.scatter(masses, lifespans, s=100, c=colors[:len(names)], zorder=5)
    for i, n in enumerate(names):
        ax2.annotate(n, (masses[i], lifespans[i]), fontsize=8,
                     xytext=(5, 5), textcoords='offset points')
    ax2.set_xscale('log'); ax2.set_yscale('log')
    ax2.set_xlabel('Body Mass (kg)'); ax2.set_ylabel('Median Lifespan (years)')
    ax2.set_title('Body Mass vs Lifespan (Log-Log)', fontweight='bold')

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname, species_results


def plot_combination_optimization(opt_results, fname='figures/fig6_optimization.png'):
    """Visualize intervention optimization results."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    df = pd.DataFrame(opt_results)

    # Bar chart of top combinations
    top = df.nlargest(10, 'median_lifespan')
    axes[0].barh(range(len(top)), top['median_lifespan'], color=sns.color_palette('viridis', len(top)))
    axes[0].set_yticks(range(len(top)))
    axes[0].set_yticklabels(top['label'], fontsize=8)
    axes[0].set_xlabel('Median Lifespan (years)')
    axes[0].set_title('Top 10 Intervention Combinations', fontweight='bold')

    # Heatmap: CR vs Rapamycin
    pivot1 = df.pivot_table(values='median_lifespan', index='cr', columns='rapa', aggfunc='mean')
    sns.heatmap(pivot1, ax=axes[1], cmap='viridis', annot=True, fmt='.1f', cbar_kws={'label': 'Lifespan (y)'})
    axes[1].set_title('CR × Rapamycin', fontweight='bold')
    axes[1].set_xlabel('Rapamycin'); axes[1].set_ylabel('Caloric Restriction')

    # Heatmap: NAD+ vs Senolytic
    pivot2 = df.pivot_table(values='median_lifespan', index='nad', columns='senolytic', aggfunc='mean')
    sns.heatmap(pivot2, ax=axes[2], cmap='viridis', annot=True, fmt='.1f', cbar_kws={'label': 'Lifespan (y)'})
    axes[2].set_title('NAD+ × Senolytic', fontweight='bold')
    axes[2].set_xlabel('Senolytic'); axes[2].set_ylabel('NAD+ Precursor')

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname


def plot_reliability_pleiotropy(sol, params, fname='figures/fig7_reliability_pleiotropy.png'):
    """Plot Reliability Theory damage and Antagonistic Pleiotropy trade-off."""
    setup_plot_style()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

    # Damage accumulation (Reliability Theory)
    ax1.plot(sol.t, sol.y[10], color='#d62728', linewidth=2)
    ax1.set_xlabel('Age (years)'); ax1.set_ylabel('Cumulative Damage D(t)')
    ax1.set_title('Reliability Theory: Damage Accumulation', fontweight='bold')
    ax1.fill_between(sol.t, 0, sol.y[10], alpha=0.2, color='#d62728')

    # Mortality hazard (Gompertz)
    h = compute_mortality_hazard(sol, params)
    ax2.semilogy(sol.t, h, color='#2ca02c', linewidth=2)
    ax2.set_xlabel('Age (years)'); ax2.set_ylabel('Mortality Hazard h(t)')
    ax2.set_title('Gompertz-like Mortality Hazard', fontweight='bold')

    # Antagonistic Pleiotropy: R vs D trade-off
    ax3.plot(sol.t, sol.y[11], color='#1f77b4', linewidth=2, label='Reproductive Fitness R(t)')
    ax3_twin = ax3.twinx()
    ax3_twin.plot(sol.t, sol.y[10], color='#d62728', linewidth=2, linestyle='--', label='Damage D(t)')
    ax3.set_xlabel('Age (years)'); ax3.set_ylabel('Reproductive Fitness', color='#1f77b4')
    ax3_twin.set_ylabel('Damage', color='#d62728')
    ax3.set_title('Antagonistic Pleiotropy Trade-off', fontweight='bold')
    ax3.axvline(x=params['t_repro'], color='gray', linestyle=':', alpha=0.7, label=f'Repro. peak (~{params["t_repro"]:.0f}y)')
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, fontsize=8)

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname


def plot_mechanism_models(fname='figures/fig8_mechanism_models.png'):
    """Plot CR, Rapamycin, NAD+ mechanism diagrams as dose-response curves."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    params = default_params()
    doses = np.linspace(0, 1, 20)

    # CR dose-response
    cr_lifespans = []
    for d in doses:
        interv = intervention_params(cr=d)
        sol = run_simulation(params=params, interv=interv)
        cr_lifespans.append(compute_median_lifespan(sol, params))
    axes[0].plot(doses * 40, cr_lifespans, 'o-', color='#2ca02c', linewidth=2)
    axes[0].set_xlabel('Caloric Restriction (%)')
    axes[0].set_ylabel('Median Lifespan (years)')
    axes[0].set_title('Caloric Restriction Dose-Response', fontweight='bold')

    # Rapamycin dose-response
    rapa_lifespans = []
    for d in doses:
        interv = intervention_params(rapa=d)
        sol = run_simulation(params=params, interv=interv)
        rapa_lifespans.append(compute_median_lifespan(sol, params))
    axes[1].plot(doses, rapa_lifespans, 's-', color='#9467bd', linewidth=2)
    axes[1].set_xlabel('Rapamycin Dose (normalized)')
    axes[1].set_ylabel('Median Lifespan (years)')
    axes[1].set_title('Rapamycin Dose-Response', fontweight='bold')

    # NAD+ dose-response
    nad_lifespans = []
    for d in doses:
        interv = intervention_params(nad=d)
        sol = run_simulation(params=params, interv=interv)
        nad_lifespans.append(compute_median_lifespan(sol, params))
    axes[2].plot(doses, nad_lifespans, '^-', color='#e377c2', linewidth=2)
    axes[2].set_xlabel('NAD+ Precursor Dose (normalized)')
    axes[2].set_ylabel('Median Lifespan (years)')
    axes[2].set_title('NAD+ Precursor Dose-Response', fontweight='bold')

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname


# ──────────────────────────────────────────────
# 8. Combination Optimization (Grid Search)
# ──────────────────────────────────────────────
def run_combination_grid(levels=None):
    """Grid search over intervention combinations."""
    if levels is None:
        levels = [0.0, 0.3, 0.6, 1.0]
    params = default_params()
    results = []
    for cr, rapa, nad, seno in product(levels, repeat=4):
        interv = intervention_params(cr=cr, rapa=rapa, nad=nad, senolytic=seno)
        sol = run_simulation(params=params, interv=interv)
        median = compute_median_lifespan(sol, params)
        results.append(dict(
            cr=cr, rapa=rapa, nad=nad, senolytic=seno,
            median_lifespan=median,
            label=f'CR={cr:.1f} R={rapa:.1f} N={nad:.1f} S={seno:.1f}'
        ))
    return results


# ──────────────────────────────────────────────
# 9. Main Execution
# ──────────────────────────────────────────────
def main():
    timestamp = datetime.now().isoformat()
    log_entries = []

    def log(phase, event, **kwargs):
        entry = dict(timestamp=datetime.now().isoformat(), phase=phase,
                     event_type=event, actor='co-scientist',
                     skill_or_tool='aging_model', **kwargs)
        log_entries.append(entry)

    log('init', 'run_started')
    print("=" * 60)
    print("Integrated Mathematical Model of Aging")
    print("=" * 60)

    params = default_params()

    # ── Baseline simulation ──
    print("\n[1/8] Running baseline simulation...")
    sol_base = run_simulation(params=params)
    baseline_median = compute_median_lifespan(sol_base, params)
    print(f"  Baseline median lifespan: {baseline_median:.1f} years")
    log('execute', 'baseline_simulation', median_lifespan=baseline_median)

    # ── Hallmarks trajectory ──
    print("[2/8] Plotting hallmark trajectories...")
    f1 = plot_hallmarks_trajectory(sol_base)
    log('execute', 'file_written', files_written=[f1])

    # ── Interaction network ──
    print("[3/8] Plotting interaction network...")
    f2 = plot_interaction_network()
    log('execute', 'file_written', files_written=[f2])

    # ── Reliability Theory & Pleiotropy ──
    print("[4/8] Plotting reliability theory & antagonistic pleiotropy...")
    f7 = plot_reliability_pleiotropy(sol_base, params)
    log('execute', 'file_written', files_written=[f7])

    # ── Intervention comparisons ──
    print("[5/8] Running intervention simulations...")
    interventions = {
        'No Treatment': intervention_params(),
        'CR 30%': intervention_params(cr=0.75),
        'Rapamycin': intervention_params(rapa=0.8),
        'NAD+ Precursor': intervention_params(nad=0.8),
        'Senolytic': intervention_params(senolytic=0.8),
        'CR+Rapa+NAD+': intervention_params(cr=0.5, rapa=0.6, nad=0.7),
        'All Combined': intervention_params(cr=0.7, rapa=0.8, nad=0.8, senolytic=0.8),
    }
    results = {}
    intervention_medians = {}
    for name, interv in interventions.items():
        sol = run_simulation(params=params, interv=interv)
        results[name] = sol
        median = compute_median_lifespan(sol, params)
        intervention_medians[name] = median
        print(f"  {name}: median = {median:.1f} years")

    f3 = plot_intervention_comparison(results, params)
    log('execute', 'file_written', files_written=[f3])

    # ── Senolytic timing ──
    print("[6/8] Simulating senolytic timing...")
    f4, seno_medians = plot_senolytic_timing()
    for k, v in seno_medians.items():
        print(f"  {k}: {v:.1f} years")
    log('execute', 'file_written', files_written=[f4])

    # ── Cross-species comparison ──
    print("[7/8] Running cross-species comparison...")
    f5, species_results = plot_species_comparison()
    for name, data in species_results.items():
        print(f"  {name}: mass={data['mass']:.2f}kg, lifespan={data['median']:.1f}y")
    log('execute', 'file_written', files_written=[f5])

    # ── Mechanism dose-response ──
    print("[7.5/8] Plotting mechanism dose-response curves...")
    f8 = plot_mechanism_models()
    log('execute', 'file_written', files_written=[f8])

    # ── Combination optimization ──
    print("[8/8] Running combination optimization grid search...")
    combo_results = run_combination_grid(levels=[0.0, 0.25, 0.5, 0.75, 1.0])
    f6 = plot_combination_optimization(combo_results)
    log('execute', 'file_written', files_written=[f6])

    # Find best combination
    df_combo = pd.DataFrame(combo_results)
    best = df_combo.loc[df_combo['median_lifespan'].idxmax()]
    print(f"\n  Best combination: {best['label']}")
    print(f"  Best median lifespan: {best['median_lifespan']:.1f} years")
    print(f"  Lifespan gain: +{best['median_lifespan'] - baseline_median:.1f} years ({(best['median_lifespan']/baseline_median - 1)*100:.1f}%)")

    # ── Save results ──
    results_summary = {
        'baseline_median_lifespan': baseline_median,
        'intervention_medians': intervention_medians,
        'senolytic_timing_medians': seno_medians,
        'species_lifespans': {k: v['median'] for k, v in species_results.items()},
        'best_combination': {
            'cr': float(best['cr']),
            'rapa': float(best['rapa']),
            'nad': float(best['nad']),
            'senolytic': float(best['senolytic']),
            'median_lifespan': float(best['median_lifespan']),
        },
        'timestamp': timestamp,
    }

    with open('results/simulation_results.json', 'w') as f:
        json.dump(results_summary, f, indent=2)

    df_combo.to_csv('results/combination_grid.csv', index=False)

    # Save species data
    df_species = pd.DataFrame(species_results).T
    df_species.to_csv('results/species_comparison.csv')

    log('report', 'report_finalized')
    log('complete', 'run_completed', status='ok')

    with open('logs/process-log.jsonl', 'w') as f:
        for entry in log_entries:
            f.write(json.dumps(entry) + '\n')

    print("\n" + "=" * 60)
    print("All simulations complete. Files saved.")
    print("=" * 60)

    return results_summary


if __name__ == '__main__':
    main()
