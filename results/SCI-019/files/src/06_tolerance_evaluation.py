#!/usr/bin/env python3
"""
Module 6: In silico evaluation of immune tolerance restoration strategies.
Models Treg expansion, antigen-specific tolerance, and combination therapies.
Uses logistic growth terms and LSODA stiff solver for stability.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

np.random.seed(42)

def tolerance_ode(t, y, params):
    """
    ODE for immune tolerance dynamics with saturation terms.
    States: [Teff, Treg, APC, AutoAg, IL2, TGFb, IL10_t, Inflammation]
    """
    Teff, Treg, APC, AutoAg, IL2, TGFb, IL10_t, Infl = [max(v, 0) for v in y]
    p = params

    K = p['K_cap']  # Carrying capacity for cell populations

    # Effector T cells: logistic growth, activated by APC+AutoAg, suppressed by Treg
    dTeff = (p['k_Teff_act'] * APC * AutoAg / (1 + p['Ki_Treg_Teff'] * Treg) *
             (1 - Teff / K) - p['d_Teff'] * Teff -
             p['k_Treg_kill'] * Treg * Teff / (Teff + p['Km_kill']))

    # Regulatory T cells: induced by TGF-beta, maintained by IL-2
    dTreg = (p['k_Treg_ind'] * TGFb / (1 + p['Ki_Infl_Treg'] * Infl) *
             (1 - Treg / K) +
             p['k_Treg_IL2'] * IL2 * Treg / (IL2 + p['Km_IL2']) * (1 - Treg / K) -
             p['d_Treg'] * Treg)

    # Antigen presenting cells
    dAPC = p['k_APC_act'] * AutoAg / (AutoAg + 1) * (1 + p['Ka_Infl_APC'] * Infl / (Infl + 1)) * (1 - APC / K) - p['d_APC'] * APC

    # Autoantigen
    dAutoAg = p['k_Ag_prod'] / (1 + 0.1 * AutoAg) - p['d_Ag'] * AutoAg - p['k_Ag_clear'] * Treg * AutoAg / (AutoAg + 1)

    # IL-2
    dIL2 = p['k_IL2_prod'] * Teff / (Teff + 1) - p['k_IL2_cons'] * (Teff + Treg) * IL2 / (IL2 + p['Km_IL2']) - p['d_IL2'] * IL2

    # TGF-beta
    dTGFb = p['k_TGFb_prod'] * Treg / (Treg + 1) + p['k_TGFb_base'] - p['d_TGFb'] * TGFb

    # IL-10 (tolerance-associated)
    dIL10_t = p['k_IL10_prod'] * Treg / (Treg + 1) - p['d_IL10'] * IL10_t

    # Inflammation score
    dInfl = (p['k_Infl_Teff'] * Teff / (Teff + 5) + p['k_Infl_APC'] * APC / (APC + 5) -
             p['k_Infl_Treg'] * Treg * Infl / (Infl + 1) / (Treg + 1) -
             p['k_Infl_IL10'] * IL10_t / (IL10_t + 1) - p['d_Infl'] * Infl)

    return [dTeff, dTreg, dAPC, dAutoAg, dIL2, dTGFb, dIL10_t, dInfl]

def get_autoimmune_params():
    """Parameters for established autoimmune state."""
    return {
        'K_cap': 15.0,
        'k_Teff_act': 2.0, 'Ki_Treg_Teff': 0.3, 'd_Teff': 0.15,
        'k_Treg_kill': 0.15, 'Km_kill': 1.0,
        'k_Treg_ind': 0.5, 'Ki_Infl_Treg': 0.5, 'k_Treg_IL2': 0.15, 'Km_IL2': 0.5, 'd_Treg': 0.12,
        'k_APC_act': 1.5, 'Ka_Infl_APC': 0.4, 'd_APC': 0.25,
        'k_Ag_prod': 1.0, 'd_Ag': 0.3, 'k_Ag_clear': 0.08,
        'k_IL2_prod': 0.8, 'k_IL2_cons': 0.4, 'd_IL2': 0.5,
        'k_TGFb_prod': 0.3, 'k_TGFb_base': 0.15, 'd_TGFb': 0.3,
        'k_IL10_prod': 0.5, 'd_IL10': 0.3,
        'k_Infl_Teff': 1.0, 'k_Infl_APC': 0.4, 'k_Infl_Treg': 0.5,
        'k_Infl_IL10': 0.4, 'd_Infl': 0.2
    }

def apply_strategy(params, strategy):
    """Apply tolerance restoration strategy."""
    p = params.copy()
    if strategy == 'Treg_expansion':
        p['k_Treg_ind'] *= 3.0
        p['k_Treg_IL2'] *= 2.0
    elif strategy == 'low_dose_IL2':
        p['k_Treg_IL2'] *= 2.5
        p['Km_IL2'] *= 0.5  # Higher IL-2 sensitivity for Treg
    elif strategy == 'tolerogenic_DC':
        p['k_APC_act'] *= 0.3
        p['k_Treg_ind'] *= 2.0
        p['k_TGFb_base'] *= 2.0
    elif strategy == 'antigen_specific':
        p['k_Ag_clear'] *= 3.0
        p['k_Treg_ind'] *= 1.5
        p['Ki_Treg_Teff'] *= 2.0
    elif strategy == 'combination':
        p['k_Treg_ind'] *= 2.5
        p['k_Treg_IL2'] *= 2.0
        p['k_APC_act'] *= 0.4
        p['k_Ag_clear'] *= 2.0
        p['k_TGFb_base'] *= 1.5
    elif strategy == 'none':
        pass
    return p

def run_tolerance_simulation():
    """Run and plot tolerance restoration simulations."""
    y0 = [3.0, 1.5, 2.0, 2.0, 0.8, 0.5, 0.5, 2.5]  # Established autoimmune state
    
    strategies = ['none', 'Treg_expansion', 'low_dose_IL2', 'tolerogenic_DC',
                  'antigen_specific', 'combination']
    strategy_labels = ['No Treatment', 'Treg Expansion', 'Low-dose IL-2',
                       'Tolerogenic DC', 'Antigen-specific', 'Combination']
    colors = ['#95A5A6', '#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6']
    
    base_params = get_autoimmune_params()
    
    # First establish disease (30 days), then treat (120 days)
    sol_pre = solve_ivp(tolerance_ode, (0, 30), y0, args=(base_params,),
                        t_eval=np.linspace(0, 30, 300), method='LSODA', max_step=0.5)
    y_disease = sol_pre.y[:, -1]
    
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    var_names = ['Teff', 'Treg', 'APC', 'AutoAg', 'IL-2', 'TGF-β', 'IL-10', 'Inflammation']
    
    steady_states = {}
    
    for strat, label, color in zip(strategies, strategy_labels, colors):
        treated_params = apply_strategy(base_params, strat)
        sol = solve_ivp(tolerance_ode, (30, 150), y_disease, args=(treated_params,),
                       t_eval=np.linspace(30, 150, 1200), method='LSODA', max_step=0.5)
        
        steady_states[label] = sol.y[:, -1]
        
        for i, (ax, vn) in enumerate(zip(axes.flatten(), var_names)):
            # Pre-treatment
            if strat == strategies[0]:
                ax.plot(sol_pre.t, sol_pre.y[i], 'k-', linewidth=1.5, alpha=0.5)
                ax.axvline(x=30, color='gray', linestyle=':', alpha=0.4)
            ax.plot(sol.t, sol.y[i], color=color, linewidth=1.5, label=label, alpha=0.8)
            ax.set_title(vn, fontsize=11, fontweight='bold')
            ax.set_xlabel('Time (days)')
            ax.grid(True, alpha=0.2)
    
    axes[0, 0].legend(fontsize=7, loc='upper right')
    plt.suptitle('In Silico Immune Tolerance Restoration', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/tolerance_restoration.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return steady_states

def plot_strategy_comparison(steady_states):
    """Compare strategies at steady state."""
    var_names = ['Teff', 'Treg', 'APC', 'AutoAg', 'IL-2', 'TGF-β', 'IL-10', 'Inflammation']
    
    # Tolerance score: high Treg, high IL-10, low Teff, low Inflammation
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    strategies = list(steady_states.keys())
    
    # Inflammation reduction
    baseline_infl = steady_states['No Treatment'][7]
    reductions = [(baseline_infl - ss[7]) / baseline_infl * 100 for ss in steady_states.values()]
    
    colors = ['#95A5A6', '#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6']
    axes[0].bar(range(len(strategies)), reductions, color=colors, edgecolor='black')
    axes[0].set_xticks(range(len(strategies)))
    axes[0].set_xticklabels(strategies, rotation=30, ha='right', fontsize=9)
    axes[0].set_ylabel('Inflammation Reduction (%)')
    axes[0].set_title('Inflammation Reduction by Strategy')
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Treg/Teff ratio
    ratios = [ss[1] / (ss[0] + 0.01) for ss in steady_states.values()]
    axes[1].bar(range(len(strategies)), ratios, color=colors, edgecolor='black')
    axes[1].set_xticks(range(len(strategies)))
    axes[1].set_xticklabels(strategies, rotation=30, ha='right', fontsize=9)
    axes[1].set_ylabel('Treg/Teff Ratio')
    axes[1].set_title('Treg/Teff Ratio by Strategy')
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/strategy_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return reductions, ratios

if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    
    print("Running tolerance restoration simulations...")
    steady_states = run_tolerance_simulation()
    
    print("Comparing strategies...")
    reductions, ratios = plot_strategy_comparison(steady_states)
    
    strategies = list(steady_states.keys())
    print("\nStrategy Comparison:")
    print(f"{'Strategy':<25} {'Infl. Reduction':<20} {'Treg/Teff Ratio':<15}")
    print("-" * 60)
    for s, r, ratio in zip(strategies, reductions, ratios):
        print(f"{s:<25} {r:>8.1f}%           {ratio:>8.3f}")
    
    print("\nModule 6 complete.")
