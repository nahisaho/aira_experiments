"""
Module 2: Energy Regeneration System Comparison
Creatine Phosphate (CP), Phosphoenolpyruvate (PEP), Maltose systems
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

# Energy regeneration kinetics for each system
energy_systems = {
    'Creatine Phosphate': {
        'k_regen': 3.5,       # regeneration rate constant (mM/min)
        'K_substrate': 2.0,   # Km for CP (mM)
        'initial_sub': 30.0,  # initial CP concentration (mM)
        'k_inhibit': 0.01,    # product inhibition constant
        'byproduct_tox': 0.005, # byproduct toxicity factor
        'cost_per_mM': 0.8,   # relative cost
        'color': '#2196F3',
        'label': 'Creatine Phosphate (CP)',
    },
    'PEP': {
        'k_regen': 2.8,
        'K_substrate': 1.5,
        'initial_sub': 30.0,
        'k_inhibit': 0.05,    # pyruvate inhibition is significant
        'byproduct_tox': 0.02,
        'cost_per_mM': 1.2,
        'color': '#FF9800',
        'label': 'Phosphoenolpyruvate (PEP)',
    },
    'Maltose': {
        'k_regen': 1.8,
        'K_substrate': 5.0,
        'initial_sub': 50.0,
        'k_inhibit': 0.002,   # minimal product inhibition
        'byproduct_tox': 0.001,
        'cost_per_mM': 0.3,
        'color': '#4CAF50',
        'label': 'Maltose',
    },
}

def energy_regen_odes(t, y, sys_params):
    ATP, substrate, byproduct = y
    ATP = max(ATP, 0)
    substrate = max(substrate, 0)
    byproduct = max(byproduct, 0)

    sp = sys_params
    # ATP regeneration (Michaelis-Menten with product inhibition)
    v_regen = sp['k_regen'] * (substrate / (sp['K_substrate'] + substrate)) * \
              (1 / (1 + sp['k_inhibit'] * byproduct))

    # ATP consumption (basal + protein synthesis)
    v_consume = 0.5 + 0.3 * ATP / (0.5 + ATP)  # simplified consumption

    # Toxicity reduces overall system activity
    tox_factor = 1 / (1 + sp['byproduct_tox'] * byproduct**2)

    dATP_dt = v_regen * tox_factor - v_consume
    dsubstrate_dt = -v_regen * 0.5  # substrate consumed
    dbyproduct_dt = v_regen * 0.3   # byproduct accumulation

    return [dATP_dt, dsubstrate_dt, dbyproduct_dt]

def run_energy_comparison(t_span=(0, 480)):
    results = {}
    for name, sp in energy_systems.items():
        y0 = [2.0, sp['initial_sub'], 0.0]  # [ATP, substrate, byproduct]
        sol = solve_ivp(energy_regen_odes, t_span, y0, args=(sp,),
                        method='RK45', dense_output=True, max_step=1.0,
                        rtol=1e-8, atol=1e-10)
        results[name] = sol
    return results

def plot_energy_comparison(results, save_path='figures/fig3_energy_systems.png'):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    t = np.linspace(0, 480, 500)

    titles = ['ATP Concentration', 'Substrate Depletion', 'Byproduct Accumulation']
    ylabels = ['[ATP] (mM)', 'Substrate Remaining (mM)', 'Byproduct (mM)']

    for name, sol in results.items():
        sp = energy_systems[name]
        y = sol.sol(t)
        for i in range(3):
            axes[i].plot(t, y[i], color=sp['color'], linewidth=2, label=sp['label'])

    for i in range(3):
        axes[i].set_xlabel('Time (min)', fontsize=11)
        axes[i].set_ylabel(ylabels[i], fontsize=11)
        axes[i].set_title(titles[i], fontsize=13, fontweight='bold')
        axes[i].legend(fontsize=9)
        axes[i].grid(True, alpha=0.3)

    plt.suptitle('Energy Regeneration System Comparison for CFPS',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

def compute_metrics(results):
    t = np.linspace(0, 480, 1000)
    metrics = {}
    for name, sol in results.items():
        y = sol.sol(t)
        sp = energy_systems[name]
        atp = y[0]
        # ATP sustain time: time above 1 mM
        above_threshold = t[atp > 1.0]
        sustain_time = above_threshold[-1] if len(above_threshold) > 0 else 0

        # Total ATP produced (integral)
        total_atp = np.trapz(np.maximum(atp, 0), t)

        # Cost efficiency
        cost_eff = total_atp / (sp['initial_sub'] * sp['cost_per_mM'])

        metrics[name] = {
            'sustain_time_min': float(sustain_time),
            'total_ATP_integral': float(total_atp),
            'peak_ATP_mM': float(np.max(atp)),
            'final_byproduct_mM': float(y[2, -1]),
            'cost_efficiency': float(cost_eff),
        }
    return metrics

def plot_radar_comparison(metrics, save_path='figures/fig4_energy_radar.png'):
    categories = ['Sustain Time', 'Total ATP', 'Peak ATP', 'Low Toxicity', 'Cost Efficiency']

    # Normalize metrics
    max_vals = {
        'sustain_time_min': max(m['sustain_time_min'] for m in metrics.values()),
        'total_ATP_integral': max(m['total_ATP_integral'] for m in metrics.values()),
        'peak_ATP_mM': max(m['peak_ATP_mM'] for m in metrics.values()),
        'final_byproduct_mM': max(m['final_byproduct_mM'] for m in metrics.values()),
        'cost_efficiency': max(m['cost_efficiency'] for m in metrics.values()),
    }

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]

    for name, m in metrics.items():
        vals = [
            m['sustain_time_min'] / max_vals['sustain_time_min'],
            m['total_ATP_integral'] / max_vals['total_ATP_integral'],
            m['peak_ATP_mM'] / max_vals['peak_ATP_mM'],
            1 - m['final_byproduct_mM'] / max_vals['final_byproduct_mM'],
            m['cost_efficiency'] / max_vals['cost_efficiency'],
        ]
        vals += vals[:1]
        ax.plot(angles, vals, 'o-', linewidth=2, label=name,
                color=energy_systems[name]['color'])
        ax.fill(angles, vals, alpha=0.1, color=energy_systems[name]['color'])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_title('Energy System Multi-Criteria Comparison', fontsize=14,
                 fontweight='bold', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=10)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

if __name__ == '__main__':
    print("=== Module 2: Energy Regeneration System Comparison ===")
    results = run_energy_comparison()
    plot_energy_comparison(results)
    metrics = compute_metrics(results)
    plot_radar_comparison(metrics)

    for name, m in metrics.items():
        print(f"  {name}: sustain={m['sustain_time_min']:.0f} min, "
              f"peak ATP={m['peak_ATP_mM']:.2f} mM, cost_eff={m['cost_efficiency']:.2f}")

    with open('results/m2_energy_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    print("  Saved: results/m2_energy_metrics.json")
