"""
Module 5: Batch → Semi-Continuous → Continuous CFPS Scale-Up Design
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

# --- Common parameters ---
base_params = {
    'k_tx': 0.8, 'k_tl': 2.0, 'k_deg_mRNA': 0.03, 'k_deg_protein': 0.001,
    'K_NTP': 0.5, 'K_aa': 0.1,
    'k_NTP_consume': 0.2, 'k_aa_consume': 0.08,
    'k_regen': 2.5, 'K_regen_sub': 2.0,
}

def batch_odes(t, y, p):
    mRNA, protein, NTP, aa, regen_sub = [max(yi, 0) for yi in y]

    v_tx = p['k_tx'] * (NTP / (p['K_NTP'] + NTP))
    v_tl = p['k_tl'] * mRNA * (aa / (p['K_aa'] + aa)) * (NTP / (p['K_NTP'] + NTP))
    v_regen = p['k_regen'] * (regen_sub / (p['K_regen_sub'] + regen_sub))

    return [
        v_tx - p['k_deg_mRNA'] * mRNA,
        v_tl - p['k_deg_protein'] * protein,
        -p['k_NTP_consume'] * (v_tx + v_tl) + v_regen,
        -p['k_aa_consume'] * v_tl,
        -v_regen * 0.3,
    ]

def semicontinuous_odes(t, y, p, feed_interval=60, feed_amount_NTP=2.0, feed_amount_aa=1.5):
    """Semi-continuous: periodic substrate feeding"""
    dydt = batch_odes(t, y, p)

    # Feeding pulses (approximated as smooth periodic function)
    pulse = np.exp(-((t % feed_interval) - 1)**2 / 0.5)  # narrow Gaussian pulse
    dydt[2] += feed_amount_NTP * pulse   # NTP feed
    dydt[3] += feed_amount_aa * pulse    # amino acid feed
    dydt[4] += 1.0 * pulse              # energy substrate feed

    return dydt

def continuous_odes(t, y, p, dilution_rate=0.005, feed_NTP=8.0, feed_aa=5.0):
    """Continuous exchange: constant feed and withdrawal"""
    mRNA, protein, NTP, aa, regen_sub = [max(yi, 0) for yi in y]

    v_tx = p['k_tx'] * (NTP / (p['K_NTP'] + NTP))
    v_tl = p['k_tl'] * mRNA * (aa / (p['K_aa'] + aa)) * (NTP / (p['K_NTP'] + NTP))
    v_regen = p['k_regen'] * (regen_sub / (p['K_regen_sub'] + regen_sub))

    D = dilution_rate

    return [
        v_tx - p['k_deg_mRNA'] * mRNA - D * mRNA,
        v_tl - p['k_deg_protein'] * protein - D * protein,
        -p['k_NTP_consume'] * (v_tx + v_tl) + v_regen + D * (feed_NTP - NTP),
        -p['k_aa_consume'] * v_tl + D * (feed_aa - aa),
        -v_regen * 0.3 + D * (10.0 - regen_sub),
    ]

def run_all_modes():
    y0 = [0.0, 0.0, 8.0, 5.0, 30.0]  # mRNA, protein, NTP, aa, energy_substrate
    t_span = (0, 720)  # 12 hours
    p = base_params.copy()

    sol_batch = solve_ivp(batch_odes, t_span, y0, args=(p,),
                          method='RK45', dense_output=True, max_step=1.0,
                          rtol=1e-8, atol=1e-10)

    sol_semi = solve_ivp(semicontinuous_odes, t_span, y0, args=(p,),
                         method='RK45', dense_output=True, max_step=0.5,
                         rtol=1e-8, atol=1e-10)

    sol_cont = solve_ivp(continuous_odes, t_span, y0, args=(p,),
                         method='RK45', dense_output=True, max_step=1.0,
                         rtol=1e-8, atol=1e-10)

    return sol_batch, sol_semi, sol_cont

def plot_scaleup_comparison(sol_batch, sol_semi, sol_cont,
                            save_path='figures/fig9_scaleup_comparison.png'):
    t = np.linspace(0, 720, 1000)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    labels_row = ['mRNA (nM)', 'Protein (nM)', 'NTP (mM)', 'Amino Acids (mM)', 'Energy Substrate (mM)']
    modes = [
        ('Batch', sol_batch, '#2196F3'),
        ('Semi-Continuous', sol_semi, '#FF9800'),
        ('Continuous (CECF)', sol_cont, '#4CAF50'),
    ]

    for idx in range(5):
        row, col = divmod(idx, 3)
        ax = axes[row][col]
        for label, sol, color in modes:
            y = sol.sol(t)
            ax.plot(t / 60, y[idx], color=color, linewidth=2, label=label)
        ax.set_xlabel('Time (hours)', fontsize=11)
        ax.set_ylabel(labels_row[idx], fontsize=11)
        ax.set_title(labels_row[idx].split('(')[0].strip(), fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Productivity comparison in last panel
    ax = axes[1][2]
    t_hours = t / 60
    for label, sol, color in modes:
        y = sol.sol(t)
        ax.plot(t_hours, y[1], color=color, linewidth=2, label=label)

    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('Cumulative Protein (nM)', fontsize=11)
    ax.set_title('Productivity Comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.suptitle('CFPS Scale-Up: Batch vs Semi-Continuous vs Continuous',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

def compute_scaleup_metrics(sol_batch, sol_semi, sol_cont):
    t = np.linspace(0, 720, 1000)
    metrics = {}
    for name, sol in [('Batch', sol_batch), ('Semi-Continuous', sol_semi), ('Continuous', sol_cont)]:
        y = sol.sol(t)
        protein = y[1]
        metrics[name] = {
            'final_protein_nM': float(protein[-1]),
            'peak_protein_nM': float(np.max(protein)),
            'time_to_50pct_peak_min': float(t[np.argmax(protein > 0.5 * np.max(protein))]),
            'volumetric_productivity_nM_per_hr': float(np.max(protein) / 12.0),
            'NTP_final_mM': float(y[2, -1]),
            'aa_final_mM': float(y[3, -1]),
        }
    return metrics

def plot_scaleup_volume_design(save_path='figures/fig10_volume_scaling.png'):
    """Design chart for volume scaling considerations"""
    volumes = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0])  # mL
    # Scaling factors based on literature
    mixing_efficiency = 1.0 / (1 + 0.01 * volumes)
    O2_transfer = 1.0 / (1 + 0.005 * volumes**0.67)
    heat_dissipation = 1.0 / (1 + 0.002 * volumes**0.5)
    cost_per_rxn = 0.5 + 0.3 * volumes  # relative cost

    relative_yield = mixing_efficiency * O2_transfer * heat_dissipation

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].semilogx(volumes, mixing_efficiency * 100, 'o-', label='Mixing Efficiency', linewidth=2)
    axes[0].semilogx(volumes, O2_transfer * 100, 's-', label='O₂ Transfer', linewidth=2)
    axes[0].semilogx(volumes, heat_dissipation * 100, '^-', label='Heat Dissipation', linewidth=2)
    axes[0].semilogx(volumes, relative_yield * 100, 'D-', label='Overall Yield', linewidth=2,
                     color='red')
    axes[0].set_xlabel('Reaction Volume (mL)', fontsize=11)
    axes[0].set_ylabel('Relative Efficiency (%)', fontsize=11)
    axes[0].set_title('Scale-Up Parameters vs Volume', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    # Cost-yield trade-off
    axes[1].semilogx(volumes, relative_yield / (cost_per_rxn / cost_per_rxn[0]), 'o-',
                     color='#4CAF50', linewidth=2, markersize=8)
    axes[1].set_xlabel('Reaction Volume (mL)', fontsize=11)
    axes[1].set_ylabel('Yield / Cost Ratio (relative)', fontsize=11)
    axes[1].set_title('Cost-Yield Trade-off', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    plt.suptitle('CFPS Volume Scale-Up Design Considerations',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

if __name__ == '__main__':
    print("=== Module 5: Scale-Up Design ===")
    sol_b, sol_s, sol_c = run_all_modes()
    plot_scaleup_comparison(sol_b, sol_s, sol_c)
    plot_scaleup_volume_design()

    metrics = compute_scaleup_metrics(sol_b, sol_s, sol_c)
    for name, m in metrics.items():
        print(f"  {name}: final={m['final_protein_nM']:.1f} nM, "
              f"peak={m['peak_protein_nM']:.1f} nM, "
              f"productivity={m['volumetric_productivity_nM_per_hr']:.1f} nM/hr")

    with open('results/m5_scaleup_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    print("  Saved: results/m5_scaleup_metrics.json")
