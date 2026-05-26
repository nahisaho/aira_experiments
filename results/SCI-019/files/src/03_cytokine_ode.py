#!/usr/bin/env python3
"""
Module 3: Cytokine network dynamic modeling using ODE systems.
Models TNF-α, IL-6, IL-17, IL-10, IFN-γ, and Treg interactions.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

np.random.seed(42)

def cytokine_ode(t, y, params):
    """
    ODE system for cytokine network dynamics.
    State variables: [TNF, IL6, IL17, IL10, IFNG, Treg, Mact, Th17_pop]
    """
    TNF, IL6, IL17, IL10, IFNG, Treg, Mact, Th17 = y
    p = params
    
    # TNF-alpha dynamics: produced by macrophages, inhibited by IL-10
    dTNF = p['k_TNF_prod'] * Mact / (1 + p['Ki_IL10_TNF'] * IL10) - p['d_TNF'] * TNF
    
    # IL-6: produced by macrophages and fibroblasts, amplified by TNF
    dIL6 = p['k_IL6_prod'] * Mact * (1 + p['Ka_TNF_IL6'] * TNF) / (1 + p['Ki_IL10_IL6'] * IL10) - p['d_IL6'] * IL6
    
    # IL-17: produced by Th17 cells, amplified by IL-6
    dIL17 = p['k_IL17_prod'] * Th17 * (1 + p['Ka_IL6_IL17'] * IL6) / (1 + p['Ki_IL10_IL17'] * IL10) - p['d_IL17'] * IL17
    
    # IL-10: produced by Tregs and M2 macrophages
    dIL10 = p['k_IL10_prod'] * Treg + p['k_IL10_Mact'] * Mact * 0.3 - p['d_IL10'] * IL10
    
    # IFN-gamma: produced by Th1/NK, modulated by IL-10
    dIFNG = p['k_IFNG_prod'] / (1 + p['Ki_IL10_IFNG'] * IL10) - p['d_IFNG'] * IFNG
    
    # Treg dynamics: suppressed in inflammatory milieu
    dTreg = p['k_Treg_prod'] / (1 + p['Ki_IL6_Treg'] * IL6 + p['Ki_TNF_Treg'] * TNF) - p['d_Treg'] * Treg
    
    # Activated macrophages: recruited by TNF and IFN-gamma
    dMact = p['k_Mact_prod'] * (1 + p['Ka_TNF_Mact'] * TNF + p['Ka_IFNG_Mact'] * IFNG) / \
            (1 + p['Ki_IL10_Mact'] * IL10) - p['d_Mact'] * Mact
    
    # Th17 population: expanded by IL-6, IL-23(implicit), suppressed by Treg
    dTh17 = p['k_Th17_prod'] * (1 + p['Ka_IL6_Th17'] * IL6) / (1 + p['Ki_Treg_Th17'] * Treg) - p['d_Th17'] * Th17
    
    return [dTNF, dIL6, dIL17, dIL10, dIFNG, dTreg, dMact, dTh17]

def get_ra_params():
    """Parameters for RA inflammatory state."""
    return {
        'k_TNF_prod': 2.5, 'Ki_IL10_TNF': 0.5, 'd_TNF': 0.8,
        'k_IL6_prod': 3.0, 'Ka_TNF_IL6': 0.8, 'Ki_IL10_IL6': 0.4, 'd_IL6': 0.7,
        'k_IL17_prod': 2.0, 'Ka_IL6_IL17': 1.0, 'Ki_IL10_IL17': 0.3, 'd_IL17': 0.6,
        'k_IL10_prod': 1.0, 'k_IL10_Mact': 0.3, 'd_IL10': 0.5,
        'k_IFNG_prod': 1.8, 'Ki_IL10_IFNG': 0.6, 'd_IFNG': 0.7,
        'k_Treg_prod': 0.8, 'Ki_IL6_Treg': 0.5, 'Ki_TNF_Treg': 0.3, 'd_Treg': 0.3,
        'k_Mact_prod': 1.5, 'Ka_TNF_Mact': 0.6, 'Ka_IFNG_Mact': 0.4, 'Ki_IL10_Mact': 0.3, 'd_Mact': 0.4,
        'k_Th17_prod': 1.2, 'Ka_IL6_Th17': 0.8, 'Ki_Treg_Th17': 0.5, 'd_Th17': 0.3
    }

def get_healthy_params():
    """Parameters for healthy immune homeostasis."""
    params = get_ra_params()
    params['k_TNF_prod'] = 0.8
    params['k_IL6_prod'] = 1.0
    params['k_IL17_prod'] = 0.5
    params['k_IL10_prod'] = 2.0
    params['k_Treg_prod'] = 2.0
    params['k_IFNG_prod'] = 0.8
    params['k_Mact_prod'] = 0.6
    params['k_Th17_prod'] = 0.4
    return params

def simulate_treatment(params, drug):
    """Simulate drug intervention on cytokine network."""
    treated = params.copy()
    if drug == 'anti-TNF':
        treated['d_TNF'] *= 3.0  # Accelerated TNF clearance
    elif drug == 'anti-IL6R':
        treated['Ka_TNF_IL6'] *= 0.1  # Block IL-6 signaling
        treated['d_IL6'] *= 2.0
    elif drug == 'JAK-inhibitor':
        treated['k_IL6_prod'] *= 0.4
        treated['k_IFNG_prod'] *= 0.5
        treated['k_IL17_prod'] *= 0.5
    elif drug == 'CTLA4-Ig':
        treated['k_Th17_prod'] *= 0.4
        treated['k_Treg_prod'] *= 1.5
    return treated

def run_simulation(params, y0, t_span=(0, 50), t_eval=None):
    """Run ODE simulation."""
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)
    sol = solve_ivp(cytokine_ode, t_span, y0, args=(params,),
                    t_eval=t_eval, method='RK45', max_step=0.1)
    return sol

def plot_cytokine_dynamics():
    """Plot cytokine dynamics for RA vs healthy and treatments."""
    y0 = [1.0, 1.0, 0.5, 1.5, 0.8, 2.0, 1.0, 0.5]
    
    # RA vs Healthy
    ra_params = get_ra_params()
    hc_params = get_healthy_params()
    
    sol_ra = run_simulation(ra_params, y0)
    sol_hc = run_simulation(hc_params, y0)
    
    var_names = ['TNF-α', 'IL-6', 'IL-17', 'IL-10', 'IFN-γ', 'Treg', 'M1 Mac', 'Th17']
    colors = ['#E74C3C', '#3498DB', '#F39C12', '#2ECC71', '#9B59B6', '#1ABC9C', '#E67E22', '#E91E63']
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    for i, (name, color) in enumerate(zip(var_names, colors)):
        axes[i].plot(sol_ra.t, sol_ra.y[i], color=color, linewidth=2, label='RA')
        axes[i].plot(sol_hc.t, sol_hc.y[i], color=color, linewidth=2, linestyle='--', label='Healthy')
        axes[i].set_title(name, fontsize=12, fontweight='bold')
        axes[i].set_xlabel('Time (days)')
        axes[i].set_ylabel('Concentration (a.u.)')
        axes[i].legend(fontsize=9)
        axes[i].grid(True, alpha=0.3)
    
    plt.suptitle('Cytokine Network Dynamics: RA vs Healthy', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/cytokine_dynamics.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_treatment_response():
    """Plot treatment effects on cytokine dynamics."""
    y0 = [1.0, 1.0, 0.5, 1.5, 0.8, 2.0, 1.0, 0.5]
    ra_params = get_ra_params()
    
    # First run RA for 20 days, then apply treatment
    sol_pre = run_simulation(ra_params, y0, t_span=(0, 20), t_eval=np.linspace(0, 20, 200))
    y_at_20 = sol_pre.y[:, -1]
    
    drugs = ['anti-TNF', 'anti-IL6R', 'JAK-inhibitor', 'CTLA4-Ig']
    drug_colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    cytokines_to_plot = [0, 1, 2, 3, 5, 7]  # TNF, IL6, IL17, IL10, Treg, Th17
    cytokine_names = ['TNF-α', 'IL-6', 'IL-17', 'IL-10', 'Treg', 'Th17']
    
    for ax_idx, (ci, cn) in enumerate(zip(cytokines_to_plot, cytokine_names)):
        ax = axes.flatten()[ax_idx]
        # Plot pre-treatment
        ax.plot(sol_pre.t, sol_pre.y[ci], 'k-', linewidth=2, label='Pre-treatment')
        ax.axvline(x=20, color='gray', linestyle=':', alpha=0.5)
        
        for drug, dc in zip(drugs, drug_colors):
            treated_params = simulate_treatment(ra_params, drug)
            sol_post = run_simulation(treated_params, y_at_20, t_span=(20, 50),
                                     t_eval=np.linspace(20, 50, 300))
            ax.plot(sol_post.t, sol_post.y[ci], color=dc, linewidth=2, label=drug)
        
        ax.set_title(cn, fontsize=12, fontweight='bold')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Concentration (a.u.)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('Treatment Effects on Cytokine Network', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/treatment_response_ode.png', dpi=150, bbox_inches='tight')
    plt.close()

def compute_steady_states():
    """Compute and compare steady-state values."""
    y0 = [1.0, 1.0, 0.5, 1.5, 0.8, 2.0, 1.0, 0.5]
    var_names = ['TNF-α', 'IL-6', 'IL-17', 'IL-10', 'IFN-γ', 'Treg', 'M1 Mac', 'Th17']
    
    results = {}
    for condition, params_fn in [('RA', get_ra_params), ('Healthy', get_healthy_params)]:
        sol = run_simulation(params_fn(), y0, t_span=(0, 100))
        ss = sol.y[:, -1]
        results[condition] = dict(zip(var_names, ss))
    
    # Treatments
    ra_params = get_ra_params()
    sol_pre = run_simulation(ra_params, y0, t_span=(0, 100))
    y_ss_ra = sol_pre.y[:, -1]
    
    for drug in ['anti-TNF', 'anti-IL6R', 'JAK-inhibitor', 'CTLA4-Ig']:
        treated = simulate_treatment(ra_params, drug)
        sol = run_simulation(treated, y_ss_ra, t_span=(0, 100))
        ss = sol.y[:, -1]
        results[drug] = dict(zip(var_names, ss))
    
    return results

if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    
    print("Running cytokine ODE dynamics...")
    plot_cytokine_dynamics()
    
    print("Running treatment simulations...")
    plot_treatment_response()
    
    print("Computing steady states...")
    ss = compute_steady_states()
    for cond, vals in ss.items():
        print(f"  {cond}:")
        for k, v in vals.items():
            print(f"    {k}: {v:.3f}")
    
    print("Module 3 complete.")
