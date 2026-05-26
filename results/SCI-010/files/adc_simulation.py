#!/usr/bin/env python3
"""
ADC Payload-Linker Optimization Computational Platform
=======================================================
Comprehensive simulation platform for Antibody-Drug Conjugate optimization:
1. DAR distribution and therapeutic window modeling
2. Linker cleavage mechanism simulation
3. Bystander effect mathematical model
4. Plasma stability vs tumor release optimization
5. PK model integration
6. HER2-targeted ADC (T-DXd analog) case study
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import os
import json
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

sns.set_theme(style="whitegrid", font_scale=1.1)
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. DAR Distribution and Therapeutic Window Model
# ============================================================

@dataclass
class DARParameters:
    max_dar: int = 8
    target_dar: float = 4.0
    dar_std: float = 1.2
    clearance_base: float = 0.3   # L/day for DAR=0
    clearance_slope: float = 0.05  # additional clearance per DAR unit
    potency_per_dar: float = 0.25  # relative potency per drug unit
    toxicity_per_dar: float = 0.06

def generate_dar_distribution(params: DARParameters, n_molecules: int = 100000) -> np.ndarray:
    """Generate DAR distribution using truncated Gaussian centered at target_dar."""
    dar_values = np.random.normal(params.target_dar, params.dar_std, n_molecules)
    dar_values = np.clip(np.round(dar_values), 0, params.max_dar).astype(int)
    return dar_values

def dar_species_fractions(dar_values: np.ndarray, max_dar: int = 8) -> np.ndarray:
    fractions = np.zeros(max_dar + 1)
    for d in range(max_dar + 1):
        fractions[d] = np.sum(dar_values == d) / len(dar_values)
    return fractions

def compute_therapeutic_window(params: DARParameters):
    """Compute efficacy and toxicity as function of average DAR."""
    avg_dars = np.linspace(0, 8, 100)
    efficacy = 1 - np.exp(-params.potency_per_dar * avg_dars)
    toxicity = 1 - np.exp(-params.toxicity_per_dar * avg_dars)
    therapeutic_index = efficacy / (toxicity + 0.01)
    return avg_dars, efficacy, toxicity, therapeutic_index

def plot_dar_distribution(params: DARParameters):
    dar_values = generate_dar_distribution(params)
    fractions = dar_species_fractions(dar_values, params.max_dar)
    avg_dars, efficacy, toxicity, ti = compute_therapeutic_window(params)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # DAR distribution
    dar_range = np.arange(params.max_dar + 1)
    colors = plt.cm.RdYlGn_r(fractions / fractions.max())
    axes[0].bar(dar_range, fractions * 100, color=colors, edgecolor='black', linewidth=0.5)
    axes[0].set_xlabel('DAR Species')
    axes[0].set_ylabel('Fraction (%)')
    axes[0].set_title(f'DAR Distribution (Mean={params.target_dar})')
    axes[0].set_xticks(dar_range)

    # Clearance by DAR
    clearances = params.clearance_base + params.clearance_slope * dar_range
    axes[1].bar(dar_range, clearances, color='steelblue', edgecolor='black', linewidth=0.5)
    axes[1].set_xlabel('DAR Species')
    axes[1].set_ylabel('Clearance (L/day)')
    axes[1].set_title('DAR-Dependent Clearance')
    axes[1].set_xticks(dar_range)

    # Therapeutic window
    axes[2].plot(avg_dars, efficacy, 'g-', lw=2, label='Efficacy')
    axes[2].plot(avg_dars, toxicity, 'r-', lw=2, label='Toxicity')
    axes[2].plot(avg_dars, ti / ti.max(), 'b--', lw=2, label='Therapeutic Index (norm)')
    optimal_dar = avg_dars[np.argmax(ti)]
    axes[2].axvline(optimal_dar, color='purple', ls=':', lw=1.5, label=f'Optimal DAR={optimal_dar:.1f}')
    axes[2].fill_between(avg_dars, efficacy, toxicity, where=efficacy > toxicity,
                         alpha=0.15, color='green', label='Therapeutic Window')
    axes[2].set_xlabel('Average DAR')
    axes[2].set_ylabel('Normalized Response')
    axes[2].set_title('Therapeutic Window Analysis')
    axes[2].legend(fontsize=8, loc='center right')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/dar_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[1] DAR distribution plot saved. Optimal DAR = {optimal_dar:.2f}")
    return fractions, optimal_dar

# ============================================================
# 2. Linker Cleavage Mechanism Simulation
# ============================================================

@dataclass
class LinkerParams:
    # Acid-sensitive (hydrazone) parameters
    acid_k_neutral: float = 0.001   # h^-1 at pH 7.4
    acid_k_acidic: float = 0.5      # h^-1 at pH 5.0
    acid_pH_half: float = 6.0
    acid_hill: float = 3.0
    # Enzyme-cleavable (Val-Cit-cathepsin B) parameters
    enzyme_vmax: float = 0.8        # h^-1
    enzyme_km: float = 5.0          # μM
    enzyme_conc_lyso: float = 20.0  # μM in lysosome
    enzyme_conc_plasma: float = 0.1 # μM in plasma
    # Reducible (disulfide) parameters
    reduce_k_plasma: float = 0.002  # h^-1 (low GSH ~2μM)
    reduce_k_cyto: float = 0.3      # h^-1 (high GSH ~10mM)
    reduce_gsh_half: float = 100.0  # μM
    reduce_hill: float = 2.0

def acid_cleavage_rate(pH: float, params: LinkerParams) -> float:
    """pH-dependent cleavage rate using Hill equation."""
    f = 1 / (1 + (pH / params.acid_pH_half) ** params.acid_hill)
    return params.acid_k_neutral + (params.acid_k_acidic - params.acid_k_neutral) * f

def enzyme_cleavage_rate(enzyme_conc: float, params: LinkerParams) -> float:
    """Michaelis-Menten enzyme cleavage."""
    return params.enzyme_vmax * enzyme_conc / (params.enzyme_km + enzyme_conc)

def reducible_cleavage_rate(gsh_conc: float, params: LinkerParams) -> float:
    """GSH-dependent disulfide reduction."""
    f = gsh_conc ** params.reduce_hill / (params.reduce_gsh_half ** params.reduce_hill + gsh_conc ** params.reduce_hill)
    return params.reduce_k_plasma + (params.reduce_k_cyto - params.reduce_k_plasma) * f

def simulate_linker_cleavage(params: LinkerParams, t_max: float = 72.0):
    t = np.linspace(0, t_max, 500)

    # pH profile: plasma (7.4) -> endosome transition at t=6h -> lysosome at t=12h
    pH_profile = np.where(t < 6, 7.4, np.where(t < 12, 7.4 - (7.4 - 5.0) * (t - 6) / 6, 5.0))
    # Enzyme concentration: low in plasma, ramps up in lysosome
    enzyme_profile = np.where(t < 6, params.enzyme_conc_plasma,
                              np.where(t < 12, params.enzyme_conc_plasma + (params.enzyme_conc_lyso - params.enzyme_conc_plasma) * (t - 6) / 6,
                                       params.enzyme_conc_lyso))
    # GSH concentration: low extracellular, high intracellular
    gsh_profile = np.where(t < 6, 2.0,
                           np.where(t < 12, 2.0 + (10000 - 2.0) * (t - 6) / 6, 10000.0))

    # Integrate cleavage for each linker type
    intact_acid = np.ones_like(t)
    intact_enzyme = np.ones_like(t)
    intact_reduce = np.ones_like(t)

    dt = t[1] - t[0]
    for i in range(1, len(t)):
        k_acid = acid_cleavage_rate(pH_profile[i], params)
        k_enzyme = enzyme_cleavage_rate(enzyme_profile[i], params)
        k_reduce = reducible_cleavage_rate(gsh_profile[i], params)
        intact_acid[i] = intact_acid[i-1] * np.exp(-k_acid * dt)
        intact_enzyme[i] = intact_enzyme[i-1] * np.exp(-k_enzyme * dt)
        intact_reduce[i] = intact_reduce[i-1] * np.exp(-k_reduce * dt)

    return t, intact_acid, intact_enzyme, intact_reduce, pH_profile, enzyme_profile, gsh_profile

def plot_linker_cleavage(params: LinkerParams):
    t, acid, enzyme, reduce, pH, enz_conc, gsh = simulate_linker_cleavage(params)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Environmental conditions
    ax = axes[0, 0]
    ax2 = ax.twinx()
    ax.plot(t, pH, 'b-', lw=2, label='pH')
    ax2.plot(t, enz_conc, 'r--', lw=2, label='Cathepsin B (μM)')
    ax.set_ylabel('pH', color='blue')
    ax2.set_ylabel('Enzyme (μM)', color='red')
    ax.set_xlabel('Time (h)')
    ax.set_title('Microenvironment Transition')
    ax.axvspan(0, 6, alpha=0.1, color='blue', label='Plasma')
    ax.axvspan(6, 12, alpha=0.1, color='yellow', label='Endosome')
    ax.axvspan(12, 72, alpha=0.1, color='red', label='Lysosome')
    ax.legend(loc='center left', fontsize=8)

    # Intact linker fraction
    ax = axes[0, 1]
    ax.plot(t, acid * 100, 'b-', lw=2, label='Acid-sensitive (Hydrazone)')
    ax.plot(t, enzyme * 100, 'r-', lw=2, label='Enzyme-cleavable (Val-Cit)')
    ax.plot(t, reduce * 100, 'g-', lw=2, label='Reducible (Disulfide)')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Intact Linker (%)')
    ax.set_title('Linker Cleavage Kinetics')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 105)

    # Payload release (= 1 - intact)
    ax = axes[1, 0]
    ax.plot(t, (1 - acid) * 100, 'b-', lw=2, label='Acid-sensitive')
    ax.plot(t, (1 - enzyme) * 100, 'r-', lw=2, label='Enzyme-cleavable')
    ax.plot(t, (1 - reduce) * 100, 'g-', lw=2, label='Reducible')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Payload Released (%)')
    ax.set_title('Cumulative Payload Release')
    ax.legend(fontsize=9)

    # Half-life comparison
    ax = axes[1, 1]
    half_lives = {}
    for name, data in [('Acid', acid), ('Enzyme', enzyme), ('Reducible', reduce)]:
        idx = np.argmin(np.abs(data - 0.5))
        half_lives[name] = t[idx]
    bars = ax.bar(half_lives.keys(), half_lives.values(),
                  color=['steelblue', 'indianred', 'seagreen'], edgecolor='black')
    ax.set_ylabel('Half-life (h)')
    ax.set_title('Linker Cleavage Half-life (post-internalization)')
    for bar, val in zip(bars, half_lives.values()):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{val:.1f}h', ha='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/linker_cleavage.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[2] Linker cleavage plot saved. Half-lives: {half_lives}")
    return half_lives

# ============================================================
# 3. Bystander Effect Mathematical Model
# ============================================================

def bystander_diffusion_model(params=None):
    """
    2D reaction-diffusion model for bystander effect in tumor tissue.
    PDE: ∂C/∂t = D·∇²C - k_uptake·C + S(x,y,t)
    """
    if params is None:
        params = {
            'D': 1e-7,          # cm²/s diffusion coefficient
            'k_uptake': 0.01,   # s^-1 cellular uptake rate
            'k_release': 0.005, # s^-1 payload release rate
            'cell_diameter': 15e-4,  # cm (15 μm)
            'tumor_radius': 0.05,    # cm (500 μm)
            'antigen_pos_frac': 0.6,
            'nx': 100,
            'dt': 0.5,         # seconds
            'total_time': 3600 * 24,  # 24 hours in seconds
            'save_interval': 3600,
        }

    nx = params['nx']
    L = params['tumor_radius'] * 2
    dx = L / nx
    dt = params['dt']
    D = params['D']
    k_up = params['k_uptake']
    k_rel = params['k_release']
    n_steps = int(params['total_time'] / dt)
    save_every = int(params['save_interval'] / dt)

    # Create spatial grid
    x = np.linspace(-L/2, L/2, nx)
    y = np.linspace(-L/2, L/2, nx)
    X, Y = np.meshgrid(x, y)

    # Antigen-positive cells: random distribution
    np.random.seed(42)
    antigen_map = np.random.rand(nx, nx) < params['antigen_pos_frac']
    # Only within tumor radius
    tumor_mask = (X**2 + Y**2) < params['tumor_radius']**2
    antigen_map = antigen_map & tumor_mask

    # Source term: payload released from Ag+ cells
    C = np.zeros((nx, nx))  # drug concentration field
    cell_viability = np.ones((nx, nx))  # 1=alive, 0=dead
    snapshots = []
    viability_snapshots = []
    times = []

    r = D * dt / dx**2
    if r > 0.25:
        dt = 0.2 * dx**2 / D
        n_steps = int(params['total_time'] / dt)
        save_every = max(1, int(params['save_interval'] / dt))

    for step in range(n_steps):
        # Source: Ag+ cells release payload
        source = k_rel * antigen_map.astype(float) * cell_viability * dt

        # 2D diffusion (explicit finite difference)
        C_new = C.copy()
        C_new[1:-1, 1:-1] = C[1:-1, 1:-1] + D * dt / dx**2 * (
            C[2:, 1:-1] + C[:-2, 1:-1] + C[1:-1, 2:] + C[1:-1, :-2] - 4 * C[1:-1, 1:-1]
        )
        # Add source and uptake
        C_new += source
        C_new -= k_up * C_new * dt * tumor_mask
        C_new = np.maximum(C_new, 0)

        # Zero-flux boundary
        C_new[0, :] = C_new[1, :]
        C_new[-1, :] = C_new[-2, :]
        C_new[:, 0] = C_new[:, 1]
        C_new[:, -1] = C_new[:, -2]

        C = C_new

        # Cell killing: probability based on drug concentration
        kill_prob = 1 - np.exp(-0.5 * C * dt)
        dying = (np.random.rand(nx, nx) < kill_prob) & tumor_mask
        cell_viability[dying] = 0

        if step % save_every == 0:
            snapshots.append(C.copy())
            viability_snapshots.append(cell_viability.copy())
            times.append(step * dt / 3600)  # hours

    return snapshots, viability_snapshots, times, X, Y, tumor_mask, antigen_map, params

def plot_bystander_effect():
    snapshots, viab_snaps, times, X, Y, tumor_mask, ag_map, params = bystander_diffusion_model()

    # Select 4 timepoints
    indices = [0, len(times)//4, len(times)//2, -1]
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))

    for col, idx in enumerate(indices):
        # Drug concentration
        ax = axes[0, col]
        im = ax.imshow(snapshots[idx], extent=[-params['tumor_radius']*1e4, params['tumor_radius']*1e4,
                        -params['tumor_radius']*1e4, params['tumor_radius']*1e4],
                       cmap='hot', origin='lower')
        ax.set_title(f't = {times[idx]:.0f} h', fontsize=11)
        if col == 0:
            ax.set_ylabel('Drug Concentration\n(μm)')
        plt.colorbar(im, ax=ax, fraction=0.046)

        # Cell viability
        ax = axes[1, col]
        viab = viab_snaps[idx].copy()
        # Color: green=alive Ag+, blue=alive Ag-, red=dead
        rgb = np.zeros((*viab.shape, 3))
        alive_agpos = (viab > 0.5) & ag_map
        alive_agneg = (viab > 0.5) & ~ag_map & tumor_mask
        dead = (viab < 0.5) & tumor_mask
        rgb[alive_agpos] = [0.2, 0.7, 0.2]  # green
        rgb[alive_agneg] = [0.2, 0.4, 0.8]  # blue
        rgb[dead] = [0.8, 0.2, 0.2]         # red
        ax.imshow(rgb, extent=[-params['tumor_radius']*1e4, params['tumor_radius']*1e4,
                  -params['tumor_radius']*1e4, params['tumor_radius']*1e4], origin='lower')
        ax.set_title(f't = {times[idx]:.0f} h', fontsize=11)
        if col == 0:
            ax.set_ylabel('Cell Viability\n(μm)')

    # Add legend for viability
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=[0.2,0.7,0.2], label='Ag+ Alive'),
                       Patch(facecolor=[0.2,0.4,0.8], label='Ag- Alive'),
                       Patch(facecolor=[0.8,0.2,0.2], label='Dead')]
    axes[1, -1].legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.suptitle('Bystander Effect: Payload Diffusion and Cell Killing', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/bystander_effect.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Summary statistics
    initial_alive = np.sum(tumor_mask)
    final_alive = np.sum(viab_snaps[-1][tumor_mask] > 0.5)
    agpos_killed = np.sum((viab_snaps[-1] < 0.5) & ag_map)
    agneg_killed = np.sum((viab_snaps[-1] < 0.5) & ~ag_map & tumor_mask)
    total_agpos = np.sum(ag_map)
    total_agneg = np.sum(tumor_mask & ~ag_map)

    stats = {
        'total_cells': int(initial_alive),
        'final_alive': int(final_alive),
        'kill_rate_pct': float((1 - final_alive/initial_alive) * 100),
        'agpos_kill_pct': float(agpos_killed / max(total_agpos, 1) * 100),
        'agneg_kill_pct': float(agneg_killed / max(total_agneg, 1) * 100),
    }
    print(f"[3] Bystander effect plot saved. Stats: {stats}")
    return stats

# ============================================================
# 4. Plasma Stability vs Tumor Release Optimization
# ============================================================

def stability_release_optimization():
    """
    Optimize the balance between plasma stability and tumor payload release.
    Model: Two-phase kinetics with stability parameter α and release parameter β.
    """
    def plasma_stability(alpha, t):
        """Fraction of intact ADC in plasma."""
        return np.exp(-alpha * t)

    def tumor_release(beta, t, t_internalize=6.0):
        """Fraction of payload released in tumor after internalization."""
        t_eff = np.maximum(t - t_internalize, 0)
        return 1 - np.exp(-beta * t_eff)

    def therapeutic_score(params_opt, t):
        alpha, beta = params_opt
        stability = plasma_stability(alpha, t)
        release = tumor_release(beta, t)
        # We want high stability in plasma (low alpha) and high release in tumor (high beta)
        # But there's a trade-off: more stable linker → slower release
        coupling = 0.3  # coupling factor
        effective_release = release * np.exp(-coupling / (beta + 0.01))
        auc_efficacy = np.trapezoid(effective_release * stability, t)
        auc_toxicity = np.trapezoid((1 - stability) * 0.5, t)
        return -(auc_efficacy - auc_toxicity)

    t = np.linspace(0, 168, 500)  # 1 week

    # Grid search for visualization
    alphas = np.linspace(0.005, 0.1, 50)
    betas = np.linspace(0.01, 1.0, 50)
    A, B = np.meshgrid(alphas, betas)
    scores = np.zeros_like(A)

    for i in range(len(betas)):
        for j in range(len(alphas)):
            scores[i, j] = -therapeutic_score([alphas[j], betas[i]], t)

    # Find optimal
    opt_idx = np.unravel_index(np.argmax(scores), scores.shape)
    opt_alpha = alphas[opt_idx[1]]
    opt_beta = betas[opt_idx[0]]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Heatmap
    ax = axes[0]
    im = ax.contourf(A, B, scores, levels=30, cmap='viridis')
    ax.plot(opt_alpha, opt_beta, 'r*', markersize=15, label=f'Optimal (α={opt_alpha:.3f}, β={opt_beta:.2f})')
    ax.set_xlabel('Plasma Deconjugation Rate α (h⁻¹)')
    ax.set_ylabel('Tumor Release Rate β (h⁻¹)')
    ax.set_title('Therapeutic Score Landscape')
    plt.colorbar(im, ax=ax)
    ax.legend(fontsize=9)

    # Optimal kinetics
    ax = axes[1]
    stab = plasma_stability(opt_alpha, t)
    rel = tumor_release(opt_beta, t)
    ax.plot(t, stab * 100, 'b-', lw=2, label='Plasma Intact ADC')
    ax.plot(t, rel * 100, 'r-', lw=2, label='Tumor Payload Released')
    ax.fill_between(t, stab * 100, alpha=0.1, color='blue')
    ax.fill_between(t, rel * 100, alpha=0.1, color='red')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Fraction (%)')
    ax.set_title(f'Optimal Kinetics (α={opt_alpha:.3f}, β={opt_beta:.2f})')
    ax.legend()

    # Compare different linker strategies
    ax = axes[2]
    strategies = {
        'Conservative': (0.01, 0.1),
        'Balanced': (opt_alpha, opt_beta),
        'Aggressive': (0.08, 0.8),
    }
    for name, (a, b) in strategies.items():
        score = -therapeutic_score([a, b], t)
        stab = plasma_stability(a, t)
        rel = tumor_release(b, t)
        eff = rel * stab
        ax.plot(t, eff * 100, lw=2, label=f'{name} (score={score:.2f})')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Effective Tumor Drug (%)')
    ax.set_title('Linker Strategy Comparison')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/stability_release.png', dpi=150, bbox_inches='tight')
    plt.close()

    results = {
        'optimal_alpha': float(opt_alpha),
        'optimal_beta': float(opt_beta),
        'optimal_score': float(scores[opt_idx]),
    }
    print(f"[4] Stability-release optimization saved. Optimal: α={opt_alpha:.4f}, β={opt_beta:.3f}")
    return results

# ============================================================
# 5. Integrated PK Model
# ============================================================

@dataclass
class PKParams:
    """Two-compartment PK model parameters for ADC (T-DXd-like)."""
    CL: float = 0.41       # L/day clearance
    V1: float = 2.74       # L central volume
    V2: float = 5.93       # L peripheral volume
    Q: float = 0.65        # L/day intercompartmental clearance
    # Payload (DXd) parameters
    CL_dxd: float = 18.5 * 24   # L/day (converted from L/h)
    V_dxd: float = 40.0   # L
    k_release: float = 0.02  # day^-1 in vivo deconjugation rate
    # Target-mediated disposition
    kon: float = 0.1       # (nM·day)^-1
    koff: float = 0.01     # day^-1
    kint: float = 0.2      # day^-1 internalization rate
    R0: float = 10.0       # nM baseline receptor
    ksyn: float = 2.0      # nM/day receptor synthesis

def pk_ode(t, y, params: PKParams, dose_times, dose_amount):
    """
    State variables:
    y[0] = ADC in central (nM)
    y[1] = ADC in peripheral (nM)
    y[2] = Free payload in central (nM)
    y[3] = Free receptor (nM)
    y[4] = ADC-receptor complex (nM)
    y[5] = Tumor ADC concentration (nM)
    """
    C1, C2, Cp, R, AR, Ct = y

    # ADC two-compartment
    dC1 = -(params.CL / params.V1) * C1 - (params.Q / params.V1) * (C1 - C2) \
          - params.kon * C1 * R + params.koff * AR - params.k_release * C1
    dC2 = (params.Q / params.V2) * (C1 - C2)

    # Free payload
    dCp = params.k_release * C1 * params.V1 / params.V_dxd \
          + params.kint * AR * params.V1 / params.V_dxd \
          - (params.CL_dxd / params.V_dxd) * Cp

    # Receptor dynamics (TMDD)
    dR = params.ksyn - params.koff * AR - params.kon * C1 * R - 0.01 * R + 0.01 * params.R0
    dAR = params.kon * C1 * R - params.koff * AR - params.kint * AR

    # Tumor ADC (simplified distribution)
    k_tumor_uptake = 0.05  # day^-1
    k_tumor_elim = 0.02    # day^-1
    dCt = k_tumor_uptake * C1 - k_tumor_elim * Ct

    return [dC1, dC2, dCp, dR, dAR, dCt]

def simulate_pk(params: PKParams, dose_mg_kg: float = 5.4, body_weight: float = 70,
                n_cycles: int = 6, cycle_days: int = 21):
    """Simulate multiple dosing cycles."""
    mw_adc = 150000  # Da
    dose_mg = dose_mg_kg * body_weight
    dose_nmol = dose_mg * 1e6 / mw_adc  # nmol
    dose_nM = dose_nmol / params.V1  # nM (V1 in L, dose in nmol → nmol/L = nM)

    total_days = n_cycles * cycle_days
    dose_times = [i * cycle_days for i in range(n_cycles)]

    t_span = (0, total_days)
    t_eval = np.linspace(0, total_days, 2000)

    y0 = [0, 0, 0, params.R0, 0, 0]

    all_t = []
    all_y = []
    current_y = y0.copy()

    for cycle in range(n_cycles):
        t_start = cycle * cycle_days
        t_end = (cycle + 1) * cycle_days
        current_y[0] += dose_nM  # bolus dose

        t_cycle = np.linspace(t_start, t_end, 500)
        sol = solve_ivp(pk_ode, (t_start, t_end), current_y,
                        t_eval=t_cycle, args=(params, dose_times, dose_nM),
                        method='LSODA', max_step=0.1)
        all_t.append(sol.t)
        all_y.append(sol.y)
        current_y = [sol.y[i, -1] for i in range(6)]

    t_full = np.concatenate(all_t)
    y_full = np.concatenate(all_y, axis=1)

    return t_full, y_full, dose_times, dose_nM

def plot_pk_model(params: PKParams):
    t, y, dose_times, dose_nM = simulate_pk(params)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # ADC plasma concentration
    ax = axes[0, 0]
    ax.semilogy(t, np.maximum(y[0], 1e-3), 'b-', lw=2, label='Central')
    ax.semilogy(t, np.maximum(y[1], 1e-3), 'b--', lw=1.5, label='Peripheral')
    for dt in dose_times:
        ax.axvline(dt, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('ADC Concentration (nM)')
    ax.set_title('ADC Pharmacokinetics')
    ax.legend()

    # Free payload
    ax = axes[0, 1]
    ax.plot(t, y[2], 'r-', lw=2)
    for dt in dose_times:
        ax.axvline(dt, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Free Payload (nM)')
    ax.set_title('Released DXd Concentration')

    # TMDD - receptor dynamics
    ax = axes[1, 0]
    ax.plot(t, y[3], 'g-', lw=2, label='Free HER2')
    ax.plot(t, y[4], 'm-', lw=2, label='ADC-HER2 Complex')
    for dt in dose_times:
        ax.axvline(dt, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Concentration (nM)')
    ax.set_title('Target-Mediated Drug Disposition')
    ax.legend()

    # Tumor ADC concentration
    ax = axes[1, 1]
    ax.plot(t, y[5], 'darkorange', lw=2)
    for dt in dose_times:
        ax.axvline(dt, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Tumor ADC (nM)')
    ax.set_title('Tumor ADC Concentration')

    plt.suptitle('Integrated PK Model (T-DXd Analog, 5.4 mg/kg Q3W)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/pk_model.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Compute PK metrics
    cmax = float(np.max(y[0]))
    # AUC over first cycle
    first_cycle = t < 21
    auc_first = float(np.trapezoid(y[0, first_cycle], t[first_cycle]))
    trough = float(y[0, np.argmin(np.abs(t - 21))])

    metrics = {
        'Cmax_nM': cmax,
        'AUC_0_21_nM_day': auc_first,
        'Ctrough_day21_nM': trough,
        'Max_free_payload_nM': float(np.max(y[2])),
        'Max_tumor_ADC_nM': float(np.max(y[5])),
    }
    print(f"[5] PK model plot saved. Cmax={cmax:.1f} nM, AUC={auc_first:.0f} nM·day")
    return metrics

# ============================================================
# 6. HER2-Targeted ADC Case Study (T-DXd Analog) with Monte Carlo
# ============================================================

def monte_carlo_efficacy_simulation(n_patients: int = 500, n_tumor_cells: int = 10000):
    """
    Monte Carlo simulation of ADC efficacy across patient population.
    Varies: HER2 expression, tumor size, DAR, clearance, body weight.
    """
    np.random.seed(123)

    # Patient parameter distributions
    her2_expression = np.random.lognormal(mean=np.log(50), sigma=0.6, size=n_patients)  # receptors/cell (thousands)
    tumor_volume = np.random.lognormal(mean=np.log(20), sigma=0.5, size=n_patients)  # cm³
    body_weight = np.random.normal(70, 15, n_patients)
    body_weight = np.clip(body_weight, 40, 130)
    dar_mean = np.random.normal(4.0, 0.3, n_patients)
    clearance = np.random.lognormal(mean=np.log(0.41), sigma=0.3, size=n_patients)

    # Dose: 5.4 mg/kg
    dose = 5.4 * body_weight  # mg

    # Simplified efficacy model
    # Exposure ∝ Dose / CL
    exposure = dose / clearance  # AUC proxy (mg·day/L)

    # Efficacy depends on exposure × HER2 expression
    drug_effect = exposure * her2_expression / 1000
    efficacy = 1 - np.exp(-0.02 * drug_effect)

    # Tumor response (RECIST-like)
    tumor_shrinkage = efficacy * 0.8 + np.random.normal(0, 0.1, n_patients)
    tumor_shrinkage = np.clip(tumor_shrinkage, -0.2, 1.0)

    # Response categories
    cr = tumor_shrinkage >= 0.95  # Complete response
    pr = (tumor_shrinkage >= 0.3) & (tumor_shrinkage < 0.95)  # Partial response
    sd = (tumor_shrinkage >= -0.2) & (tumor_shrinkage < 0.3)  # Stable disease
    pd = tumor_shrinkage < -0.2   # Progressive disease

    orr = float((np.sum(cr) + np.sum(pr)) / n_patients * 100)
    dcr = float((np.sum(cr) + np.sum(pr) + np.sum(sd)) / n_patients * 100)

    # Safety: Grade 3+ AE probability
    toxicity_prob = 0.3 * (clearance < np.percentile(clearance, 25)).astype(float) + \
                    0.1 * (dar_mean > 5).astype(float) + \
                    np.random.uniform(0, 0.15, n_patients)
    grade3_ae = toxicity_prob > 0.3

    # Plotting
    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

    # 1. HER2 expression distribution
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(her2_expression, bins=40, color='steelblue', edgecolor='black', alpha=0.8)
    ax.set_xlabel('HER2 Expression (×10³ receptors/cell)')
    ax.set_ylabel('Count')
    ax.set_title('HER2 Expression Distribution')

    # 2. Exposure vs HER2
    ax = fig.add_subplot(gs[0, 1])
    scatter = ax.scatter(her2_expression, exposure, c=tumor_shrinkage, cmap='RdYlGn',
                         s=10, alpha=0.6)
    plt.colorbar(scatter, ax=ax, label='Tumor Shrinkage')
    ax.set_xlabel('HER2 Expression')
    ax.set_ylabel('Exposure (AUC proxy)')
    ax.set_title('Exposure-Response Relationship')

    # 3. Waterfall plot (tumor response)
    ax = fig.add_subplot(gs[0, 2])
    sorted_shrinkage = np.sort(tumor_shrinkage)[::-1] * 100
    colors_wf = ['green' if s >= 30 else 'blue' if s >= 0 else 'red' for s in sorted_shrinkage]
    ax.bar(range(n_patients), sorted_shrinkage, color=colors_wf, width=1.0)
    ax.axhline(30, color='black', ls='--', lw=1, alpha=0.5)
    ax.axhline(0, color='black', ls='-', lw=0.5)
    ax.set_xlabel('Patients')
    ax.set_ylabel('Tumor Change (%)')
    ax.set_title(f'Waterfall Plot (ORR={orr:.1f}%)')

    # 4. Response categories
    ax = fig.add_subplot(gs[1, 0])
    categories = ['CR', 'PR', 'SD', 'PD']
    counts = [np.sum(cr), np.sum(pr), np.sum(sd), np.sum(pd)]
    colors_cat = ['darkgreen', 'limegreen', 'gold', 'red']
    bars = ax.bar(categories, counts, color=colors_cat, edgecolor='black')
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{count}\n({count/n_patients*100:.1f}%)', ha='center', fontsize=9)
    ax.set_ylabel('Number of Patients')
    ax.set_title('Response Categories (RECIST)')

    # 5. DAR impact on efficacy
    ax = fig.add_subplot(gs[1, 1])
    dar_bins = np.linspace(2.5, 5.5, 7)
    dar_groups = np.digitize(dar_mean, dar_bins)
    for g in range(1, len(dar_bins)):
        mask = dar_groups == g
        if np.sum(mask) > 0:
            ax.scatter(dar_mean[mask], tumor_shrinkage[mask] * 100, s=10, alpha=0.5,
                       label=f'DAR {dar_bins[g-1]:.1f}-{dar_bins[g]:.1f}' if g <= 3 else None)
    ax.set_xlabel('Average DAR')
    ax.set_ylabel('Tumor Shrinkage (%)')
    ax.set_title('DAR vs Efficacy')

    # 6. Clearance vs toxicity
    ax = fig.add_subplot(gs[1, 2])
    ax.scatter(clearance[~grade3_ae], exposure[~grade3_ae], s=10, alpha=0.4, c='blue', label='No Grade 3+ AE')
    ax.scatter(clearance[grade3_ae], exposure[grade3_ae], s=15, alpha=0.6, c='red', label='Grade 3+ AE')
    ax.set_xlabel('Clearance (L/day)')
    ax.set_ylabel('Exposure')
    ax.set_title('Safety Profile')
    ax.legend(fontsize=8)

    # 7. PFS Kaplan-Meier-like curve
    ax = fig.add_subplot(gs[2, 0])
    pfs_months = -np.log(1 - efficacy + 0.01) * 8 + np.random.exponential(2, n_patients)
    pfs_months = np.clip(pfs_months, 0.5, 36)
    sorted_pfs = np.sort(pfs_months)
    survival = np.linspace(1, 0, len(sorted_pfs))
    ax.step(sorted_pfs, survival * 100, 'b-', lw=2)
    median_pfs = sorted_pfs[len(sorted_pfs)//2]
    ax.axhline(50, color='gray', ls='--', alpha=0.5)
    ax.axvline(median_pfs, color='red', ls=':', label=f'Median PFS={median_pfs:.1f} mo')
    ax.set_xlabel('Time (months)')
    ax.set_ylabel('PFS (%)')
    ax.set_title('Progression-Free Survival')
    ax.legend()

    # 8. Body weight vs exposure
    ax = fig.add_subplot(gs[2, 1])
    ax.scatter(body_weight, exposure, c=efficacy, cmap='viridis', s=10, alpha=0.5)
    ax.set_xlabel('Body Weight (kg)')
    ax.set_ylabel('Exposure (AUC proxy)')
    ax.set_title('Weight-Based Dosing Impact')

    # 9. Summary table
    ax = fig.add_subplot(gs[2, 2])
    ax.axis('off')
    summary_data = [
        ['Parameter', 'Value'],
        ['N patients', str(n_patients)],
        ['ORR', f'{orr:.1f}%'],
        ['DCR', f'{dcr:.1f}%'],
        ['Median PFS', f'{median_pfs:.1f} months'],
        ['Grade 3+ AE', f'{np.sum(grade3_ae)/n_patients*100:.1f}%'],
        ['Mean DAR', f'{np.mean(dar_mean):.2f}'],
        ['Mean CL', f'{np.mean(clearance):.3f} L/day'],
    ]
    table = ax.table(cellText=summary_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)
    for i in range(len(summary_data)):
        for j in range(2):
            cell = table[i, j]
            if i == 0:
                cell.set_facecolor('#4472C4')
                cell.set_text_props(color='white', fontweight='bold')
            else:
                cell.set_facecolor('#D9E2F3' if i % 2 == 0 else 'white')
    ax.set_title('Clinical Summary', fontsize=12, fontweight='bold', pad=20)

    plt.suptitle('HER2-Targeted ADC (T-DXd Analog) Monte Carlo Case Study', fontsize=14, fontweight='bold')
    plt.savefig(f'{FIGDIR}/case_study_tdxd.png', dpi=150, bbox_inches='tight')
    plt.close()

    results = {
        'ORR_pct': orr,
        'DCR_pct': dcr,
        'median_PFS_months': float(median_pfs),
        'grade3_AE_pct': float(np.sum(grade3_ae) / n_patients * 100),
        'CR': int(np.sum(cr)),
        'PR': int(np.sum(pr)),
        'SD': int(np.sum(sd)),
        'PD': int(np.sum(pd)),
        'mean_tumor_shrinkage_pct': float(np.mean(tumor_shrinkage) * 100),
    }
    print(f"[6] Case study plot saved. ORR={orr:.1f}%, Median PFS={median_pfs:.1f} mo")
    return results

# ============================================================
# 7. DAR Optimization Monte Carlo
# ============================================================

def dar_optimization_monte_carlo(n_simulations: int = 10000):
    """Monte Carlo optimization of DAR for maximizing therapeutic index."""
    np.random.seed(456)

    dar_targets = np.linspace(1, 8, 50)
    results = {'dar': [], 'ti_mean': [], 'ti_std': [], 'efficacy': [], 'toxicity': []}

    for dar_target in dar_targets:
        ti_samples = []
        eff_samples = []
        tox_samples = []
        for _ in range(n_simulations // len(dar_targets)):
            dar = np.random.normal(dar_target, 1.0)
            dar = np.clip(dar, 0, 8)
            clearance = 0.3 + 0.05 * dar + np.random.normal(0, 0.02)
            efficacy = 1 - np.exp(-0.25 * dar) + np.random.normal(0, 0.03)
            toxicity = 1 - np.exp(-0.06 * dar) + np.random.normal(0, 0.02)
            efficacy = np.clip(efficacy, 0, 1)
            toxicity = np.clip(toxicity, 0.01, 1)
            ti = efficacy / toxicity
            ti_samples.append(ti)
            eff_samples.append(efficacy)
            tox_samples.append(toxicity)

        results['dar'].append(dar_target)
        results['ti_mean'].append(np.mean(ti_samples))
        results['ti_std'].append(np.std(ti_samples))
        results['efficacy'].append(np.mean(eff_samples))
        results['toxicity'].append(np.mean(tox_samples))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Therapeutic index vs DAR
    ax = axes[0]
    ti_mean = np.array(results['ti_mean'])
    ti_std = np.array(results['ti_std'])
    dar_arr = np.array(results['dar'])
    ax.plot(dar_arr, ti_mean, 'b-', lw=2)
    ax.fill_between(dar_arr, ti_mean - ti_std, ti_mean + ti_std, alpha=0.2, color='blue')
    opt_dar = dar_arr[np.argmax(ti_mean)]
    ax.axvline(opt_dar, color='red', ls='--', label=f'Optimal DAR={opt_dar:.1f}')
    ax.set_xlabel('Target DAR')
    ax.set_ylabel('Therapeutic Index')
    ax.set_title('Monte Carlo DAR Optimization')
    ax.legend()

    # Efficacy and toxicity
    ax = axes[1]
    ax.plot(dar_arr, results['efficacy'], 'g-', lw=2, label='Efficacy')
    ax.plot(dar_arr, results['toxicity'], 'r-', lw=2, label='Toxicity')
    ax.fill_between(dar_arr, results['efficacy'], results['toxicity'],
                    where=np.array(results['efficacy']) > np.array(results['toxicity']),
                    alpha=0.15, color='green')
    ax.set_xlabel('Target DAR')
    ax.set_ylabel('Probability')
    ax.set_title('Efficacy vs Toxicity Trade-off')
    ax.legend()

    # Distribution of TI at optimal DAR
    ax = axes[2]
    np.random.seed(789)
    dar_opt_samples = np.random.normal(opt_dar, 1.0, 5000)
    dar_opt_samples = np.clip(dar_opt_samples, 0, 8)
    eff = 1 - np.exp(-0.25 * dar_opt_samples)
    tox = np.maximum(1 - np.exp(-0.06 * dar_opt_samples), 0.01)
    ti_dist = eff / tox
    ax.hist(ti_dist, bins=50, color='steelblue', edgecolor='black', alpha=0.8, density=True)
    ax.axvline(np.mean(ti_dist), color='red', ls='--', lw=2, label=f'Mean TI={np.mean(ti_dist):.2f}')
    ax.set_xlabel('Therapeutic Index')
    ax.set_ylabel('Density')
    ax.set_title(f'TI Distribution at Optimal DAR={opt_dar:.1f}')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/dar_optimization_mc.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"[7] DAR optimization MC saved. Optimal DAR = {opt_dar:.1f}, Mean TI = {np.max(ti_mean):.2f}")
    return {'optimal_dar': float(opt_dar), 'max_ti': float(np.max(ti_mean))}

# ============================================================
# Main execution
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("ADC Payload-Linker Optimization Computational Platform")
    print("=" * 70)

    all_results = {}

    print("\n--- Module 1: DAR Distribution & Therapeutic Window ---")
    dar_params = DARParameters()
    fractions, optimal_dar = plot_dar_distribution(dar_params)
    all_results['dar'] = {'fractions': fractions.tolist(), 'optimal_dar': optimal_dar}

    print("\n--- Module 2: Linker Cleavage Mechanisms ---")
    linker_params = LinkerParams()
    half_lives = plot_linker_cleavage(linker_params)
    all_results['linker'] = half_lives

    print("\n--- Module 3: Bystander Effect Model ---")
    bystander_stats = plot_bystander_effect()
    all_results['bystander'] = bystander_stats

    print("\n--- Module 4: Plasma Stability vs Tumor Release ---")
    stability_results = stability_release_optimization()
    all_results['stability'] = stability_results

    print("\n--- Module 5: Integrated PK Model ---")
    pk_params = PKParams()
    pk_metrics = plot_pk_model(pk_params)
    all_results['pk'] = pk_metrics

    print("\n--- Module 6: HER2-Targeted ADC Case Study ---")
    case_results = monte_carlo_efficacy_simulation()
    all_results['case_study'] = case_results

    print("\n--- Module 7: DAR Optimization Monte Carlo ---")
    mc_results = dar_optimization_monte_carlo()
    all_results['dar_mc'] = mc_results

    # Save results
    with open('simulation_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)

    print("\n" + "=" * 70)
    print("All simulations completed successfully!")
    print(f"Results saved to simulation_results.json")
    print(f"Figures saved to {FIGDIR}/")
    print("=" * 70)
