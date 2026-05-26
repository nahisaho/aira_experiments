#!/usr/bin/env python3
"""
Cell-Free Protein Synthesis (CFPS) Optimization Framework
ODE-based modeling with Bayesian optimization integration.

Modules:
1. Transcription-Translation Coupled Model (resource competition)
2. Energy Regeneration System Comparison (CP, PEP, Maltose)
3. Mg2+/K+/Polyamine Optimization Map
4. mRNA Stability and Ribosome Loading Prediction
5. Batch → Semi-continuous → Continuous Scale-up
6. Membrane Protein Expression (Nanodisc) Case Study
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from itertools import product
import warnings
warnings.filterwarnings('ignore')

sns.set_theme(style="whitegrid", font_scale=1.1)
FIGDIR = "figures"

# =============================================================================
# 1. Transcription-Translation Coupled ODE Model
# =============================================================================
def cfps_ode(t, y, params):
    """
    Coupled transcription-translation ODE with resource competition.
    State: [DNA, mRNA, Ribosome_free, Protein, NTP, AA, ATP]
    """
    DNA, mRNA, Ribo_free, Protein, NTP, AA, ATP = y

    k_tx = params['k_tx']        # transcription rate
    k_tl = params['k_tl']        # translation rate
    k_mdeg = params['k_mdeg']    # mRNA degradation rate
    K_NTP = params['K_NTP']      # Michaelis for NTP
    K_AA = params['K_AA']        # Michaelis for amino acids
    K_ATP = params['K_ATP']      # Michaelis for ATP
    n_ntp = params['n_ntp']      # NTP per mRNA
    n_aa = params['n_aa']        # AA per protein
    n_atp_tx = params['n_atp_tx']  # ATP per mRNA (transcription)
    n_atp_tl = params['n_atp_tl']  # ATP per protein (translation)
    Ribo_total = params['Ribo_total']
    k_regen = params.get('k_regen', 0)
    E_sub = params.get('E_sub', 0)  # energy substrate concentration

    # Resource saturation
    sat_NTP = NTP / (K_NTP + NTP) if NTP > 0 else 0
    sat_AA = AA / (K_AA + AA) if AA > 0 else 0
    sat_ATP = ATP / (K_ATP + ATP) if ATP > 0 else 0

    # Ribosome competition: fraction bound to mRNA
    Ribo_bound = Ribo_total - max(Ribo_free, 0)

    # Transcription rate
    v_tx = k_tx * DNA * sat_NTP * sat_ATP

    # Translation rate
    v_tl = k_tl * mRNA * max(Ribo_free, 0) * sat_AA * sat_ATP / (1 + mRNA / 10)

    # Energy regeneration
    v_regen = k_regen * E_sub * (1 - ATP / (ATP + 1))

    dDNA = 0
    dMRNA = v_tx - k_mdeg * mRNA
    dRibo_free = -v_tl * 0.01 + k_mdeg * mRNA * (Ribo_total - Ribo_free) / (mRNA + 1)
    dProtein = v_tl
    dNTP = -v_tx * n_ntp
    dAA = -v_tl * n_aa
    dATP = -v_tx * n_atp_tx - v_tl * n_atp_tl + v_regen

    return [dDNA, dMRNA, dRibo_free, dProtein, dNTP, dAA, dATP]


def run_cfps_simulation(params, t_span=(0, 240), t_eval=None):
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)

    y0 = [
        params.get('DNA0', 10),      # nM
        params.get('mRNA0', 0),      # nM
        params.get('Ribo0', 500),    # nM
        params.get('Protein0', 0),   # µM
        params.get('NTP0', 3000),    # µM
        params.get('AA0', 2000),     # µM
        params.get('ATP0', 1500),    # µM
    ]

    sol = solve_ivp(cfps_ode, t_span, y0, t_eval=t_eval, args=(params,),
                    method='RK45', max_step=1.0, rtol=1e-8, atol=1e-10)
    return sol


def plot_transcription_translation():
    """Figure 1: Transcription-Translation coupled dynamics."""
    params_base = {
        'k_tx': 0.5, 'k_tl': 0.08, 'k_mdeg': 0.02,
        'K_NTP': 200, 'K_AA': 300, 'K_ATP': 100,
        'n_ntp': 3, 'n_aa': 0.5, 'n_atp_tx': 2, 'n_atp_tl': 4,
        'Ribo_total': 500, 'k_regen': 0, 'E_sub': 0,
        'DNA0': 10, 'mRNA0': 0, 'Ribo0': 500,
        'Protein0': 0, 'NTP0': 3000, 'AA0': 2000, 'ATP0': 1500,
    }

    # Different DNA concentrations to show resource competition
    dna_levels = [5, 10, 20, 40]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    colors = sns.color_palette("viridis", len(dna_levels))

    for i, dna in enumerate(dna_levels):
        p = params_base.copy()
        p['DNA0'] = dna
        sol = run_cfps_simulation(p)

        axes[0, 0].plot(sol.t, sol.y[1], color=colors[i], label=f'DNA={dna} nM')
        axes[0, 1].plot(sol.t, sol.y[3], color=colors[i], label=f'DNA={dna} nM')
        axes[1, 0].plot(sol.t, sol.y[6], color=colors[i], label=f'DNA={dna} nM')
        axes[1, 1].plot(sol.t, sol.y[2], color=colors[i], label=f'DNA={dna} nM')

    axes[0, 0].set_title('mRNA Dynamics')
    axes[0, 0].set_ylabel('mRNA (nM)')
    axes[0, 1].set_title('Protein Production')
    axes[0, 1].set_ylabel('Protein (µM)')
    axes[1, 0].set_title('ATP Consumption')
    axes[1, 0].set_ylabel('ATP (µM)')
    axes[1, 1].set_title('Free Ribosome Dynamics')
    axes[1, 1].set_ylabel('Free Ribosomes (nM)')

    for ax in axes.flat:
        ax.set_xlabel('Time (min)')
        ax.legend(fontsize=9)

    fig.suptitle('Transcription-Translation Coupled Model with Resource Competition', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig1_tx_tl_dynamics.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 1 saved: fig1_tx_tl_dynamics.png")

    return params_base


# =============================================================================
# 2. Energy Regeneration System Comparison
# =============================================================================
def energy_regeneration_params():
    """Parameters for CP, PEP, and Maltose energy systems."""
    systems = {
        'Creatine Phosphate': {
            'k_regen': 0.15, 'E_sub': 60,   # mM
            'efficiency': 1.0, 'cost_factor': 3.5,
            'byproduct_inhibition': 0.03,
        },
        'PEP': {
            'k_regen': 0.12, 'E_sub': 30,
            'efficiency': 0.85, 'cost_factor': 5.0,
            'byproduct_inhibition': 0.05,
        },
        'Maltose': {
            'k_regen': 0.08, 'E_sub': 80,
            'efficiency': 0.7, 'cost_factor': 0.5,
            'byproduct_inhibition': 0.01,
        }
    }
    return systems


def plot_energy_comparison():
    """Figure 2: Energy regeneration system comparison."""
    base_params = {
        'k_tx': 0.5, 'k_tl': 0.08, 'k_mdeg': 0.02,
        'K_NTP': 200, 'K_AA': 300, 'K_ATP': 100,
        'n_ntp': 3, 'n_aa': 0.5, 'n_atp_tx': 2, 'n_atp_tl': 4,
        'Ribo_total': 500,
        'DNA0': 10, 'mRNA0': 0, 'Ribo0': 500,
        'Protein0': 0, 'NTP0': 3000, 'AA0': 2000, 'ATP0': 1500,
    }

    systems = energy_regeneration_params()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = {'Creatine Phosphate': '#2196F3', 'PEP': '#FF5722', 'Maltose': '#4CAF50'}

    final_proteins = {}
    final_atp = {}

    # No energy regen baseline
    p_no = base_params.copy()
    p_no['k_regen'] = 0
    p_no['E_sub'] = 0
    sol_no = run_cfps_simulation(p_no, t_span=(0, 360))
    axes[0, 0].plot(sol_no.t, sol_no.y[3], 'k--', label='No regen', alpha=0.6)
    axes[0, 1].plot(sol_no.t, sol_no.y[6], 'k--', label='No regen', alpha=0.6)
    final_proteins['No Regen'] = sol_no.y[3, -1]

    for name, sys_p in systems.items():
        p = base_params.copy()
        p['k_regen'] = sys_p['k_regen'] * sys_p['efficiency']
        p['E_sub'] = sys_p['E_sub']
        sol = run_cfps_simulation(p, t_span=(0, 360))

        axes[0, 0].plot(sol.t, sol.y[3], color=colors[name], label=name, linewidth=2)
        axes[0, 1].plot(sol.t, sol.y[6], color=colors[name], label=name, linewidth=2)

        final_proteins[name] = sol.y[3, -1]
        final_atp[name] = sol.y[6, -1]

    axes[0, 0].set_title('Protein Yield Comparison')
    axes[0, 0].set_ylabel('Protein (µM)')
    axes[0, 0].set_xlabel('Time (min)')
    axes[0, 0].legend()

    axes[0, 1].set_title('ATP Dynamics')
    axes[0, 1].set_ylabel('ATP (µM)')
    axes[0, 1].set_xlabel('Time (min)')
    axes[0, 1].legend()

    # Bar chart: final yields
    names = list(final_proteins.keys())
    vals = list(final_proteins.values())
    bar_colors = ['gray'] + [colors[n] for n in systems.keys()]
    axes[1, 0].bar(names, vals, color=bar_colors)
    axes[1, 0].set_title('Final Protein Yield')
    axes[1, 0].set_ylabel('Protein (µM)')
    axes[1, 0].tick_params(axis='x', rotation=15)

    # Cost-efficiency
    cost_eff = {}
    for name, sys_p in systems.items():
        cost_eff[name] = final_proteins[name] / sys_p['cost_factor']
    ce_names = list(cost_eff.keys())
    ce_vals = list(cost_eff.values())
    ce_colors = [colors[n] for n in ce_names]
    axes[1, 1].bar(ce_names, ce_vals, color=ce_colors)
    axes[1, 1].set_title('Cost-Efficiency (Yield / Cost)')
    axes[1, 1].set_ylabel('µM / cost unit')
    axes[1, 1].tick_params(axis='x', rotation=15)

    fig.suptitle('Energy Regeneration System Comparison', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig2_energy_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 2 saved: fig2_energy_comparison.png")

    return final_proteins, cost_eff


# =============================================================================
# 3. Mg2+/K+/Polyamine Optimization Map
# =============================================================================
def ion_effect_model(Mg, K, PA):
    """
    Model for ion/polyamine effects on CFPS yield.
    Based on empirical response surface:
    - Mg2+: optimal ~10-15 mM (bell-shaped)
    - K+: optimal ~100-200 mM (bell-shaped)
    - Polyamine (spermidine): optimal ~1-2 mM (bell-shaped)
    """
    # Bell-shaped responses
    y_Mg = np.exp(-0.5 * ((Mg - 12) / 4) ** 2)
    y_K = np.exp(-0.5 * ((K - 150) / 50) ** 2)
    y_PA = np.exp(-0.5 * ((PA - 1.5) / 0.8) ** 2)

    # Synergistic interaction
    synergy = 1 + 0.15 * np.exp(-0.5 * ((Mg - 12) / 3) ** 2) * np.exp(-0.5 * ((PA - 1.5) / 0.6) ** 2)

    # Inhibition at extreme values
    inhibit = 1 / (1 + 0.001 * (Mg ** 2 + K ** 2 / 100 + PA ** 2 * 10))

    yield_val = y_Mg * y_K * y_PA * synergy * inhibit * 100  # scale to µg/mL
    return max(yield_val, 0)


def plot_ion_optimization():
    """Figure 3: Ion optimization landscape."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Mg vs K (PA fixed at 1.5)
    mg_range = np.linspace(2, 25, 60)
    k_range = np.linspace(50, 300, 60)
    Mg_grid, K_grid = np.meshgrid(mg_range, k_range)
    Z1 = np.vectorize(lambda m, k: ion_effect_model(m, k, 1.5))(Mg_grid, K_grid)

    im1 = axes[0].contourf(Mg_grid, K_grid, Z1, levels=20, cmap='RdYlGn')
    axes[0].set_xlabel('Mg²⁺ (mM)')
    axes[0].set_ylabel('K⁺ (mM)')
    axes[0].set_title('Yield: Mg²⁺ vs K⁺\n(Spermidine = 1.5 mM)')
    plt.colorbar(im1, ax=axes[0], label='Yield (µg/mL)')
    # Mark optimum
    opt_idx = np.unravel_index(Z1.argmax(), Z1.shape)
    axes[0].plot(Mg_grid[opt_idx], K_grid[opt_idx], 'w*', markersize=15)

    # Mg vs PA (K fixed at 150)
    pa_range = np.linspace(0, 4, 60)
    Mg_grid2, PA_grid = np.meshgrid(mg_range, pa_range)
    Z2 = np.vectorize(lambda m, p: ion_effect_model(m, 150, p))(Mg_grid2, PA_grid)

    im2 = axes[1].contourf(Mg_grid2, PA_grid, Z2, levels=20, cmap='RdYlGn')
    axes[1].set_xlabel('Mg²⁺ (mM)')
    axes[1].set_ylabel('Spermidine (mM)')
    axes[1].set_title('Yield: Mg²⁺ vs Spermidine\n(K⁺ = 150 mM)')
    plt.colorbar(im2, ax=axes[1], label='Yield (µg/mL)')
    opt_idx2 = np.unravel_index(Z2.argmax(), Z2.shape)
    axes[1].plot(Mg_grid2[opt_idx2], PA_grid[opt_idx2], 'w*', markersize=15)

    # K vs PA (Mg fixed at 12)
    K_grid3, PA_grid3 = np.meshgrid(k_range, pa_range)
    Z3 = np.vectorize(lambda k, p: ion_effect_model(12, k, p))(K_grid3, PA_grid3)

    im3 = axes[2].contourf(K_grid3, PA_grid3, Z3, levels=20, cmap='RdYlGn')
    axes[2].set_xlabel('K⁺ (mM)')
    axes[2].set_ylabel('Spermidine (mM)')
    axes[2].set_title('Yield: K⁺ vs Spermidine\n(Mg²⁺ = 12 mM)')
    plt.colorbar(im3, ax=axes[2], label='Yield (µg/mL)')
    opt_idx3 = np.unravel_index(Z3.argmax(), Z3.shape)
    axes[2].plot(K_grid3[opt_idx3], PA_grid3[opt_idx3], 'w*', markersize=15)

    fig.suptitle('Ion Concentration Optimization Map for CFPS', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig3_ion_optimization.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 3 saved: fig3_ion_optimization.png")


# =============================================================================
# 4. mRNA Stability and Ribosome Loading Prediction
# =============================================================================
def mrna_stability_model(t, y, params):
    """
    mRNA stability model with ribosome loading effects.
    State: [mRNA, Ribo_loaded, Protein, mRNA_degraded]
    """
    mRNA, Ribo_loaded, Protein, mRNA_deg = y

    k_bind = params['k_bind']      # ribosome binding rate
    k_unbind = params['k_unbind']  # ribosome release rate
    k_tl = params['k_tl']         # translation elongation rate
    k_deg_base = params['k_deg_base']  # basal mRNA degradation
    k_deg_ribo = params['k_deg_ribo']  # ribosome-loading dependent degradation
    Ribo_total = params['Ribo_total']
    mRNA_len = params['mRNA_len']  # nt length

    Ribo_free = max(Ribo_total - Ribo_loaded, 0)

    # Ribosome loading
    v_bind = k_bind * mRNA * Ribo_free
    v_unbind = k_unbind * Ribo_loaded

    # Translation
    v_tl = k_tl * Ribo_loaded

    # mRNA degradation: depends on ribosome load
    ribo_per_mrna = Ribo_loaded / (mRNA + 0.01)
    k_deg_eff = k_deg_base + k_deg_ribo * ribo_per_mrna

    v_deg = k_deg_eff * mRNA

    dmRNA = -v_deg
    dRibo = v_bind - v_unbind
    dProtein = v_tl
    dmRNA_deg = v_deg

    return [dmRNA, dRibo, dProtein, dmRNA_deg]


def plot_mrna_stability():
    """Figure 4: mRNA stability and ribosome loading."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Vary ribosome-dependent degradation
    k_deg_ribo_vals = [0.0, 0.005, 0.01, 0.02]
    colors = sns.color_palette("magma", len(k_deg_ribo_vals))

    base_p = {
        'k_bind': 0.1, 'k_unbind': 0.05, 'k_tl': 0.5,
        'k_deg_base': 0.01, 'k_deg_ribo': 0.01,
        'Ribo_total': 300, 'mRNA_len': 900,
    }

    for i, kdr in enumerate(k_deg_ribo_vals):
        p = base_p.copy()
        p['k_deg_ribo'] = kdr
        y0 = [100, 0, 0, 0]  # mRNA=100 nM
        sol = solve_ivp(mrna_stability_model, (0, 180), y0, args=(p,),
                        t_eval=np.linspace(0, 180, 300), method='RK45',
                        max_step=0.5, rtol=1e-8, atol=1e-10)

        lbl = f'k_deg_ribo={kdr}'
        axes[0, 0].plot(sol.t, sol.y[0], color=colors[i], label=lbl)
        axes[0, 1].plot(sol.t, sol.y[1], color=colors[i], label=lbl)
        axes[1, 0].plot(sol.t, sol.y[2], color=colors[i], label=lbl)

    axes[0, 0].set_title('mRNA Decay')
    axes[0, 0].set_ylabel('mRNA (nM)')
    axes[0, 1].set_title('Loaded Ribosomes')
    axes[0, 1].set_ylabel('Ribosomes (nM)')
    axes[1, 0].set_title('Protein Accumulation')
    axes[1, 0].set_ylabel('Protein (µM)')

    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set_xlabel('Time (min)')
        ax.legend(fontsize=8)

    # mRNA length effect on half-life
    lengths = np.linspace(300, 3000, 50)
    half_lives = []
    for l in lengths:
        p = base_p.copy()
        p['mRNA_len'] = l
        p['k_deg_ribo'] = 0.01
        # Approximate: higher length -> more ribosomes -> faster decay
        hl = np.log(2) / (p['k_deg_base'] + p['k_deg_ribo'] * min(l / 300, 5))
        half_lives.append(hl)

    axes[1, 1].plot(lengths, half_lives, 'b-', linewidth=2)
    axes[1, 1].set_xlabel('mRNA Length (nt)')
    axes[1, 1].set_ylabel('Predicted Half-life (min)')
    axes[1, 1].set_title('mRNA Length vs Half-life')
    axes[1, 1].axhline(y=np.log(2) / 0.01, color='gray', linestyle='--', alpha=0.5, label='Basal t½')
    axes[1, 1].legend()

    fig.suptitle('mRNA Stability and Ribosome Loading Model', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig4_mrna_stability.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 4 saved: fig4_mrna_stability.png")


# =============================================================================
# 5. Batch → Semi-continuous → Continuous Scale-up
# =============================================================================
def batch_mode(params, t_span=(0, 240)):
    return run_cfps_simulation(params, t_span=t_span)


def semicontinuous_mode(params, t_span=(0, 480), feed_interval=60, feed_frac=0.3):
    """Semi-continuous: periodic substrate feeding."""
    t_eval = np.linspace(t_span[0], t_span[1], 800)
    y0 = [
        params.get('DNA0', 10), params.get('mRNA0', 0), params.get('Ribo0', 500),
        params.get('Protein0', 0), params.get('NTP0', 3000),
        params.get('AA0', 2000), params.get('ATP0', 1500),
    ]

    all_t = []
    all_y = [[] for _ in range(7)]
    current_y = np.array(y0)
    current_t = t_span[0]

    feed_times = np.arange(feed_interval, t_span[1], feed_interval)

    segments = []
    prev_t = current_t
    for ft in feed_times:
        segments.append((prev_t, ft))
        prev_t = ft
    segments.append((prev_t, t_span[1]))

    for seg_start, seg_end in segments:
        t_seg = np.linspace(seg_start, seg_end, 100)
        sol = solve_ivp(cfps_ode, (seg_start, seg_end), current_y, t_eval=t_seg,
                        args=(params,), method='RK45', max_step=1.0, rtol=1e-8, atol=1e-10)

        all_t.extend(sol.t.tolist())
        for j in range(7):
            all_y[j].extend(sol.y[j].tolist())

        current_y = sol.y[:, -1].copy()
        # Feed: replenish substrates
        current_y[4] = min(current_y[4] + params.get('NTP0', 3000) * feed_frac, params.get('NTP0', 3000))
        current_y[5] = min(current_y[5] + params.get('AA0', 2000) * feed_frac, params.get('AA0', 2000))
        current_y[6] = min(current_y[6] + params.get('ATP0', 1500) * feed_frac, params.get('ATP0', 1500))

    class Result:
        pass
    r = Result()
    r.t = np.array(all_t)
    r.y = np.array(all_y)
    return r


def continuous_mode(params, t_span=(0, 720), dilution_rate=0.005):
    """Continuous exchange: steady substrate supply with dilution."""
    def continuous_ode(t, y, params, D):
        dydt = cfps_ode(t, y, params)
        # Continuous feed of substrates and dilution
        NTP0 = params.get('NTP0', 3000)
        AA0 = params.get('AA0', 2000)
        ATP0 = params.get('ATP0', 1500)

        dydt[4] += D * (NTP0 - y[4])  # NTP feed
        dydt[5] += D * (AA0 - y[5])    # AA feed
        dydt[6] += D * (ATP0 - y[6])   # ATP feed
        dydt[3] -= D * y[3] * 0.1       # partial product removal
        return dydt

    y0 = [
        params.get('DNA0', 10), params.get('mRNA0', 0), params.get('Ribo0', 500),
        params.get('Protein0', 0), params.get('NTP0', 3000),
        params.get('AA0', 2000), params.get('ATP0', 1500),
    ]

    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    sol = solve_ivp(continuous_ode, t_span, y0, t_eval=t_eval,
                    args=(params, dilution_rate), method='RK45',
                    max_step=1.0, rtol=1e-8, atol=1e-10)
    return sol


def plot_scaleup():
    """Figure 5: Scale-up comparison."""
    params = {
        'k_tx': 0.5, 'k_tl': 0.08, 'k_mdeg': 0.02,
        'K_NTP': 200, 'K_AA': 300, 'K_ATP': 100,
        'n_ntp': 3, 'n_aa': 0.5, 'n_atp_tx': 2, 'n_atp_tl': 4,
        'Ribo_total': 500, 'k_regen': 0.1, 'E_sub': 50,
        'DNA0': 10, 'mRNA0': 0, 'Ribo0': 500,
        'Protein0': 0, 'NTP0': 3000, 'AA0': 2000, 'ATP0': 1500,
    }

    sol_batch = batch_mode(params, t_span=(0, 240))
    sol_semi = semicontinuous_mode(params, t_span=(0, 480))
    sol_cont = continuous_mode(params, t_span=(0, 720))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Protein
    axes[0, 0].plot(sol_batch.t, sol_batch.y[3], 'b-', label='Batch', linewidth=2)
    axes[0, 0].plot(sol_semi.t, sol_semi.y[3], 'r-', label='Semi-continuous', linewidth=2)
    axes[0, 0].plot(sol_cont.t, sol_cont.y[3], 'g-', label='Continuous', linewidth=2)
    axes[0, 0].set_title('Protein Production')
    axes[0, 0].set_ylabel('Protein (µM)')
    axes[0, 0].legend()

    # ATP
    axes[0, 1].plot(sol_batch.t, sol_batch.y[6], 'b-', label='Batch', linewidth=2)
    axes[0, 1].plot(sol_semi.t, sol_semi.y[6], 'r-', label='Semi-continuous', linewidth=2)
    axes[0, 1].plot(sol_cont.t, sol_cont.y[6], 'g-', label='Continuous', linewidth=2)
    axes[0, 1].set_title('ATP Dynamics')
    axes[0, 1].set_ylabel('ATP (µM)')
    axes[0, 1].legend()

    # Productivity (protein/time)
    batch_prod = sol_batch.y[3] / (sol_batch.t + 1) * 60  # per hour
    semi_prod = sol_semi.y[3] / (sol_semi.t + 1) * 60
    cont_prod = sol_cont.y[3] / (sol_cont.t + 1) * 60

    axes[1, 0].plot(sol_batch.t, batch_prod, 'b-', label='Batch', linewidth=2)
    axes[1, 0].plot(sol_semi.t, semi_prod, 'r-', label='Semi-continuous', linewidth=2)
    axes[1, 0].plot(sol_cont.t, cont_prod, 'g-', label='Continuous', linewidth=2)
    axes[1, 0].set_title('Volumetric Productivity')
    axes[1, 0].set_ylabel('Productivity (µM/h)')
    axes[1, 0].legend()

    # Summary bar chart
    modes = ['Batch\n(4h)', 'Semi-cont.\n(8h)', 'Continuous\n(12h)']
    final_yield = [sol_batch.y[3, -1], sol_semi.y[3, -1], sol_cont.y[3, -1]]
    bar_colors = ['#2196F3', '#FF5722', '#4CAF50']
    axes[1, 1].bar(modes, final_yield, color=bar_colors)
    axes[1, 1].set_title('Final Protein Yield')
    axes[1, 1].set_ylabel('Protein (µM)')

    for ax in axes.flat:
        ax.set_xlabel('Time (min)')

    fig.suptitle('Scale-up: Batch → Semi-continuous → Continuous', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig5_scaleup.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 5 saved: fig5_scaleup.png")

    return {'batch': sol_batch.y[3, -1], 'semi': sol_semi.y[3, -1], 'cont': sol_cont.y[3, -1]}


# =============================================================================
# 6. Membrane Protein Expression with Nanodisc
# =============================================================================
def nanodisc_ode(t, y, params):
    """
    Membrane protein expression with nanodisc co-translational insertion.
    State: [mRNA, Protein_soluble, Protein_inserted, Nanodisc_free, Aggregated]
    """
    mRNA, P_sol, P_ins, ND_free, Agg = y

    k_tl = params['k_tl']
    k_mdeg = params['k_mdeg']
    k_insert = params['k_insert']  # insertion rate into nanodiscs
    k_agg = params['k_agg']       # aggregation rate
    K_nd = params['K_nd']          # nanodisc saturation constant

    # Translation produces soluble (uninserted) membrane protein
    v_tl = k_tl * mRNA
    v_mdeg = k_mdeg * mRNA

    # Insertion into nanodiscs (Michaelis-Menten)
    v_insert = k_insert * P_sol * ND_free / (K_nd + ND_free)

    # Aggregation of uninserted protein
    v_agg = k_agg * P_sol ** 2

    dmRNA = -v_mdeg
    dP_sol = v_tl - v_insert - v_agg
    dP_ins = v_insert
    dND = -v_insert * 0.1  # each insertion occupies partial nanodisc
    dAgg = v_agg

    return [dmRNA, dP_sol, dP_ins, dND, dAgg]


def plot_nanodisc():
    """Figure 6: Nanodisc-assisted membrane protein expression."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    nd_concentrations = [0, 5, 20, 50, 100]
    colors = sns.color_palette("coolwarm", len(nd_concentrations))

    base_p = {
        'k_tl': 0.5, 'k_mdeg': 0.015,
        'k_insert': 0.05, 'k_agg': 0.002,
        'K_nd': 10,
    }

    insertion_eff = []

    for i, nd0 in enumerate(nd_concentrations):
        p = base_p.copy()
        y0 = [50, 0, 0, nd0, 0]  # mRNA=50 nM
        sol = solve_ivp(nanodisc_ode, (0, 240), y0,
                        t_eval=np.linspace(0, 240, 400), args=(p,),
                        method='RK45', max_step=0.5, rtol=1e-8, atol=1e-10)

        lbl = f'ND={nd0} µM'
        axes[0, 0].plot(sol.t, sol.y[2], color=colors[i], label=lbl, linewidth=2)
        axes[0, 1].plot(sol.t, sol.y[4], color=colors[i], label=lbl, linewidth=2)

        total_prot = sol.y[1, -1] + sol.y[2, -1] + sol.y[4, -1]
        eff = sol.y[2, -1] / total_prot * 100 if total_prot > 0 else 0
        insertion_eff.append(eff)

    axes[0, 0].set_title('Inserted Membrane Protein')
    axes[0, 0].set_ylabel('Protein in Nanodiscs (µM)')
    axes[0, 0].legend(fontsize=9)

    axes[0, 1].set_title('Aggregated Protein')
    axes[0, 1].set_ylabel('Aggregated (µM)')
    axes[0, 1].legend(fontsize=9)

    # Insertion efficiency
    axes[1, 0].plot(nd_concentrations, insertion_eff, 'o-', color='#673AB7', linewidth=2, markersize=8)
    axes[1, 0].set_xlabel('Nanodisc Concentration (µM)')
    axes[1, 0].set_ylabel('Insertion Efficiency (%)')
    axes[1, 0].set_title('Nanodisc Concentration vs Insertion Efficiency')

    # Optimal nanodisc/protein ratio
    ratios = np.linspace(0.1, 5, 50)
    mRNA_fixed = 50
    eff_ratio = []
    for r in ratios:
        nd0 = r * mRNA_fixed * 0.5
        p = base_p.copy()
        y0 = [mRNA_fixed, 0, 0, nd0, 0]
        sol = solve_ivp(nanodisc_ode, (0, 240), y0,
                        t_eval=np.linspace(0, 240, 200), args=(p,),
                        method='RK45', max_step=0.5, rtol=1e-8, atol=1e-10)
        total = sol.y[1, -1] + sol.y[2, -1] + sol.y[4, -1]
        eff_ratio.append(sol.y[2, -1] / total * 100 if total > 0 else 0)

    axes[1, 1].plot(ratios, eff_ratio, 'g-', linewidth=2)
    axes[1, 1].set_xlabel('Nanodisc:Protein Ratio')
    axes[1, 1].set_ylabel('Insertion Efficiency (%)')
    axes[1, 1].set_title('Optimal Nanodisc:Protein Ratio')
    axes[1, 1].axhline(y=max(eff_ratio) * 0.9, color='gray', linestyle='--', alpha=0.5, label='90% max')
    axes[1, 1].legend()

    for ax in [axes[0, 0], axes[0, 1]]:
        ax.set_xlabel('Time (min)')

    fig.suptitle('Membrane Protein Expression with Nanodisc Integration', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig6_nanodisc.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 6 saved: fig6_nanodisc.png")

    return insertion_eff


# =============================================================================
# 7. Bayesian Optimization Integration
# =============================================================================
def objective_function(x):
    """
    Objective for Bayesian optimization.
    x = [Mg, K, Polyamine, DNA_conc, energy_substrate_conc]
    """
    Mg, K, PA, DNA, E_sub = x

    # Ion effects
    ion_yield = ion_effect_model(Mg, K, PA)

    # Run ODE simulation
    params = {
        'k_tx': 0.5, 'k_tl': 0.08, 'k_mdeg': 0.02,
        'K_NTP': 200, 'K_AA': 300, 'K_ATP': 100,
        'n_ntp': 3, 'n_aa': 0.5, 'n_atp_tx': 2, 'n_atp_tl': 4,
        'Ribo_total': 500, 'k_regen': 0.12, 'E_sub': E_sub,
        'DNA0': DNA, 'mRNA0': 0, 'Ribo0': 500,
        'Protein0': 0, 'NTP0': 3000, 'AA0': 2000, 'ATP0': 1500,
    }

    sol = run_cfps_simulation(params, t_span=(0, 240))
    protein_yield = sol.y[3, -1]

    # Combined objective (ion effect * ODE yield)
    combined = ion_yield * protein_yield / 100
    return -combined  # minimize negative


def bayesian_optimization(n_iter=80, n_init=20):
    """
    Simple Bayesian optimization using random search + local refinement.
    Bounds: Mg[2-25], K[50-300], PA[0-4], DNA[1-50], E_sub[10-100]
    """
    np.random.seed(42)
    bounds = np.array([[2, 25], [50, 300], [0, 4], [1, 50], [10, 100]])

    # Phase 1: Latin Hypercube-like initial sampling
    X_history = []
    y_history = []

    for _ in range(n_init):
        x = bounds[:, 0] + np.random.rand(5) * (bounds[:, 1] - bounds[:, 0])
        y = objective_function(x)
        X_history.append(x)
        y_history.append(y)

    # Phase 2: Iterative optimization with perturbation
    best_idx = np.argmin(y_history)
    best_x = X_history[best_idx].copy()
    best_y = y_history[best_idx]

    for i in range(n_iter - n_init):
        # Perturbation around best with decreasing variance
        scale = 0.3 * (1 - i / (n_iter - n_init))
        perturbation = np.random.randn(5) * (bounds[:, 1] - bounds[:, 0]) * scale
        x_new = np.clip(best_x + perturbation, bounds[:, 0], bounds[:, 1])

        y_new = objective_function(x_new)
        X_history.append(x_new)
        y_history.append(y_new)

        if y_new < best_y:
            best_y = y_new
            best_x = x_new.copy()

    return np.array(X_history), np.array(y_history), best_x, -best_y


def plot_bayesian_optimization():
    """Figure 7: Bayesian optimization convergence and results."""
    X_hist, y_hist, best_x, best_yield = bayesian_optimization()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Convergence
    best_so_far = np.minimum.accumulate(-y_hist)
    axes[0, 0].plot(range(len(y_hist)), -y_hist, 'o', alpha=0.3, markersize=4, color='gray', label='Samples')
    axes[0, 0].plot(range(len(y_hist)), best_so_far, 'r-', linewidth=2, label='Best so far')
    axes[0, 0].set_xlabel('Iteration')
    axes[0, 0].set_ylabel('Yield (combined score)')
    axes[0, 0].set_title('Bayesian Optimization Convergence')
    axes[0, 0].legend()

    # Parameter distributions
    param_names = ['Mg²⁺ (mM)', 'K⁺ (mM)', 'Spermidine (mM)', 'DNA (nM)', 'Energy Sub. (mM)']
    top_n = 10
    top_idx = np.argsort(y_hist)[:top_n]

    for i, (ax_idx, pname) in enumerate(zip([(0, 1), (1, 0), (1, 1)], param_names[:3])):
        ax = axes[ax_idx]
        ax.scatter(X_hist[:, i], -y_hist, c='lightblue', alpha=0.5, s=20, label='All')
        ax.scatter(X_hist[top_idx, i], -y_hist[top_idx], c='red', s=50, label=f'Top {top_n}')
        ax.axvline(best_x[i], color='green', linestyle='--', label=f'Optimal: {best_x[i]:.1f}')
        ax.set_xlabel(pname)
        ax.set_ylabel('Yield')
        ax.set_title(f'Yield vs {pname}')
        ax.legend(fontsize=8)

    fig.suptitle('Bayesian Optimization of CFPS Parameters', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig7_bayesian_opt.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Figure 7 saved: fig7_bayesian_opt.png")
    print(f"\nOptimal Parameters:")
    for name, val in zip(param_names, best_x):
        print(f"  {name}: {val:.2f}")
    print(f"  Optimal Yield: {best_yield:.4f}")

    return best_x, best_yield, X_hist, y_hist


# =============================================================================
# 8. Comprehensive Summary Figure
# =============================================================================
def plot_summary():
    """Figure 8: Framework overview summary."""
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    # 1. TX-TL dynamics mini
    ax1 = fig.add_subplot(gs[0, 0])
    params = {
        'k_tx': 0.5, 'k_tl': 0.08, 'k_mdeg': 0.02,
        'K_NTP': 200, 'K_AA': 300, 'K_ATP': 100,
        'n_ntp': 3, 'n_aa': 0.5, 'n_atp_tx': 2, 'n_atp_tl': 4,
        'Ribo_total': 500, 'k_regen': 0.1, 'E_sub': 50,
        'DNA0': 10, 'mRNA0': 0, 'Ribo0': 500,
        'Protein0': 0, 'NTP0': 3000, 'AA0': 2000, 'ATP0': 1500,
    }
    sol = run_cfps_simulation(params)
    ax1.plot(sol.t, sol.y[3], 'b-', linewidth=2)
    ax1.set_title('(a) TX-TL Model', fontsize=10)
    ax1.set_ylabel('Protein (µM)')
    ax1.set_xlabel('Time (min)')

    # 2. Energy comparison mini
    ax2 = fig.add_subplot(gs[0, 1])
    e_sys = energy_regeneration_params()
    e_names = list(e_sys.keys()) + ['No Regen']
    e_yields = []
    for name, sp in e_sys.items():
        p = params.copy()
        p['k_regen'] = sp['k_regen'] * sp['efficiency']
        p['E_sub'] = sp['E_sub']
        s = run_cfps_simulation(p, t_span=(0, 360))
        e_yields.append(s.y[3, -1])
    p_no = params.copy()
    p_no['k_regen'] = 0; p_no['E_sub'] = 0
    s_no = run_cfps_simulation(p_no, t_span=(0, 360))
    e_yields.append(s_no.y[3, -1])
    ax2.barh(e_names, e_yields, color=['#2196F3', '#FF5722', '#4CAF50', 'gray'])
    ax2.set_title('(b) Energy Systems', fontsize=10)
    ax2.set_xlabel('Yield (µM)')

    # 3. Ion optimization mini
    ax3 = fig.add_subplot(gs[0, 2])
    mg_r = np.linspace(2, 25, 40)
    k_r = np.linspace(50, 300, 40)
    Mg_g, K_g = np.meshgrid(mg_r, k_r)
    Z = np.vectorize(lambda m, k: ion_effect_model(m, k, 1.5))(Mg_g, K_g)
    ax3.contourf(Mg_g, K_g, Z, levels=15, cmap='RdYlGn')
    ax3.set_title('(c) Ion Optimization', fontsize=10)
    ax3.set_xlabel('Mg²⁺ (mM)')
    ax3.set_ylabel('K⁺ (mM)')

    # 4. mRNA stability
    ax4 = fig.add_subplot(gs[1, 0])
    for kdr, col in zip([0, 0.01, 0.02], ['green', 'orange', 'red']):
        p_m = {'k_bind': 0.1, 'k_unbind': 0.05, 'k_tl': 0.5,
               'k_deg_base': 0.01, 'k_deg_ribo': kdr, 'Ribo_total': 300, 'mRNA_len': 900}
        s_m = solve_ivp(mrna_stability_model, (0, 180), [100, 0, 0, 0], args=(p_m,),
                        t_eval=np.linspace(0, 180, 200), method='RK45', max_step=0.5, rtol=1e-8, atol=1e-10)
        ax4.plot(s_m.t, s_m.y[0], color=col, label=f'k={kdr}')
    ax4.set_title('(d) mRNA Stability', fontsize=10)
    ax4.set_ylabel('mRNA (nM)')
    ax4.set_xlabel('Time (min)')
    ax4.legend(fontsize=7)

    # 5. Scale-up
    ax5 = fig.add_subplot(gs[1, 1])
    s_b = batch_mode(params, t_span=(0, 240))
    s_sc = semicontinuous_mode(params, t_span=(0, 480))
    s_c = continuous_mode(params, t_span=(0, 720))
    ax5.plot(s_b.t, s_b.y[3], 'b-', label='Batch')
    ax5.plot(s_sc.t, s_sc.y[3], 'r-', label='Semi-cont.')
    ax5.plot(s_c.t, s_c.y[3], 'g-', label='Continuous')
    ax5.set_title('(e) Scale-up Modes', fontsize=10)
    ax5.set_ylabel('Protein (µM)')
    ax5.set_xlabel('Time (min)')
    ax5.legend(fontsize=7)

    # 6. Nanodisc
    ax6 = fig.add_subplot(gs[1, 2])
    nd_vals = [0, 20, 50, 100]
    nd_eff = []
    for nd0 in nd_vals:
        p_nd = {'k_tl': 0.5, 'k_mdeg': 0.015, 'k_insert': 0.05, 'k_agg': 0.002, 'K_nd': 10}
        s_nd = solve_ivp(nanodisc_ode, (0, 240), [50, 0, 0, nd0, 0], args=(p_nd,),
                         t_eval=np.linspace(0, 240, 200), method='RK45', max_step=0.5, rtol=1e-8, atol=1e-10)
        total = s_nd.y[1, -1] + s_nd.y[2, -1] + s_nd.y[4, -1]
        nd_eff.append(s_nd.y[2, -1] / total * 100 if total > 0 else 0)
    ax6.bar([str(n) for n in nd_vals], nd_eff, color=sns.color_palette("coolwarm", len(nd_vals)))
    ax6.set_title('(f) Nanodisc Effect', fontsize=10)
    ax6.set_ylabel('Insertion Eff. (%)')
    ax6.set_xlabel('ND conc. (µM)')

    # 7-9: Bayesian opt summary
    ax7 = fig.add_subplot(gs[2, :])
    X_h, y_h, bx, by = bayesian_optimization(n_iter=60, n_init=15)
    best_sf = np.minimum.accumulate(-y_h)
    ax7.plot(range(len(y_h)), -y_h, 'o', alpha=0.3, markersize=4, color='lightblue')
    ax7.plot(range(len(y_h)), best_sf, 'r-', linewidth=2)
    ax7.set_title('(g) Bayesian Optimization Convergence', fontsize=10)
    ax7.set_xlabel('Iteration')
    ax7.set_ylabel('Combined Yield Score')
    ax7.annotate(f'Optimal: {by:.2f}', xy=(len(y_h) - 1, by), fontsize=10,
                 xytext=(len(y_h) * 0.6, by * 0.8),
                 arrowprops=dict(arrowstyle='->', color='red'),
                 color='red', fontweight='bold')

    fig.suptitle('CFPS Optimization Framework: Comprehensive Overview', fontsize=16, y=1.01)
    plt.savefig(f'{FIGDIR}/fig8_summary.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 8 saved: fig8_summary.png")


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("CFPS Optimization Framework - Simulation Suite")
    print("=" * 60)

    print("\n[1/7] Transcription-Translation Coupled Model...")
    plot_transcription_translation()

    print("\n[2/7] Energy Regeneration Comparison...")
    final_proteins, cost_eff = plot_energy_comparison()
    print(f"  Final yields: {final_proteins}")
    print(f"  Cost efficiency: {cost_eff}")

    print("\n[3/7] Ion Concentration Optimization...")
    plot_ion_optimization()

    print("\n[4/7] mRNA Stability & Ribosome Loading...")
    plot_mrna_stability()

    print("\n[5/7] Scale-up Comparison...")
    scaleup = plot_scaleup()
    print(f"  Scale-up yields: {scaleup}")

    print("\n[6/7] Nanodisc Membrane Protein Expression...")
    nd_eff = plot_nanodisc()
    print(f"  Insertion efficiencies: {nd_eff}")

    print("\n[7/7] Bayesian Optimization...")
    best_x, best_yield, X_hist, y_hist = plot_bayesian_optimization()

    print("\n[Bonus] Summary Figure...")
    plot_summary()

    print("\n" + "=" * 60)
    print("All simulations complete. Figures saved in figures/")
    print("=" * 60)
