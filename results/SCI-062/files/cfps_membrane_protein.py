"""
Module 6: Membrane Protein Expression with Nanodisc Integration — Case Study
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

# --- Parameters for membrane protein expression in nanodiscs ---
params_mp = {
    # Expression
    'k_tx': 0.5,           # lower for MP due to template complexity
    'k_tl': 1.2,           # reduced translation rate
    'k_deg_mRNA': 0.04,    # slightly faster degradation
    'k_deg_protein': 0.002,
    'K_NTP': 0.5,
    'K_aa': 0.15,

    # Nanodisc parameters
    'ND_total': 50.0,      # total nanodisc concentration (µM)
    'K_insertion': 5.0,    # Km for co-translational insertion (µM)
    'k_insert': 0.3,       # insertion rate constant
    'k_misfolded': 0.15,   # misfolding rate (no nanodisc)
    'k_aggregation': 0.05, # aggregation rate for free MP

    # Energy
    'NTP_init': 8.0,
    'aa_init': 5.0,
    'k_NTP_consume': 0.2,
    'k_aa_consume': 0.1,
    'k_regen': 2.0,
    'regen_sub_init': 30.0,
    'K_regen_sub': 2.0,
}

def mp_nanodisc_odes(t, y, p):
    """
    State variables:
    y[0]: mRNA
    y[1]: nascent (unfolded) membrane protein
    y[2]: nanodisc-inserted (correctly folded) MP
    y[3]: misfolded/aggregated MP
    y[4]: free nanodiscs
    y[5]: NTP
    y[6]: amino acids
    y[7]: energy substrate
    """
    mRNA, MP_nascent, MP_inserted, MP_misfolded, ND_free, NTP, aa, regen_sub = \
        [max(yi, 0) for yi in y]

    v_tx = p['k_tx'] * (NTP / (p['K_NTP'] + NTP))
    v_tl = p['k_tl'] * mRNA * (aa / (p['K_aa'] + aa)) * (NTP / (p['K_NTP'] + NTP))
    v_regen = p['k_regen'] * (regen_sub / (p['K_regen_sub'] + regen_sub))

    # Co-translational insertion into nanodiscs
    v_insert = p['k_insert'] * MP_nascent * (ND_free / (p['K_insertion'] + ND_free))

    # Misfolding (competes with insertion)
    v_misfold = p['k_misfolded'] * MP_nascent * (1 - ND_free / (p['K_insertion'] + ND_free))

    # Aggregation of misfolded
    v_aggregate = p['k_aggregation'] * MP_misfolded

    return [
        v_tx - p['k_deg_mRNA'] * mRNA,                    # mRNA
        v_tl - v_insert - v_misfold - v_aggregate * 0.1,  # nascent MP
        v_insert,                                           # inserted MP (stable in ND)
        v_misfold - v_aggregate,                           # misfolded
        -v_insert,                                          # free nanodiscs consumed
        -p['k_NTP_consume'] * (v_tx + v_tl) + v_regen,    # NTP
        -p['k_aa_consume'] * v_tl,                         # amino acids
        -v_regen * 0.3,                                     # energy substrate
    ]

def run_nanodisc_titration():
    """Test different nanodisc concentrations"""
    nd_concs = [0, 10, 25, 50, 100, 200]
    results = {}

    for nd in nd_concs:
        p = params_mp.copy()
        p['ND_total'] = nd
        y0 = [0, 0, 0, 0, nd, p['NTP_init'], p['aa_init'], p['regen_sub_init']]
        sol = solve_ivp(mp_nanodisc_odes, [0, 360], y0, args=(p,),
                        method='RK45', dense_output=True, max_step=0.5,
                        rtol=1e-8, atol=1e-10)
        results[nd] = sol
    return results

def plot_nanodisc_case_study(results, save_path='figures/fig11_nanodisc_mp.png'):
    t = np.linspace(0, 360, 500)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    cmap = plt.cm.viridis
    nd_concs = sorted(results.keys())
    colors = [cmap(i / (len(nd_concs) - 1)) for i in range(len(nd_concs))]

    # Panel 1: Inserted MP vs time for each ND conc
    for i, nd in enumerate(nd_concs):
        y = results[nd].sol(t)
        axes[0, 0].plot(t, y[2], color=colors[i], linewidth=2, label=f'ND={nd} µM')
    axes[0, 0].set_xlabel('Time (min)', fontsize=11)
    axes[0, 0].set_ylabel('Inserted MP (µM)', fontsize=11)
    axes[0, 0].set_title('Correctly Folded MP in Nanodiscs', fontsize=13, fontweight='bold')
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].grid(True, alpha=0.3)

    # Panel 2: Misfolded MP
    for i, nd in enumerate(nd_concs):
        y = results[nd].sol(t)
        axes[0, 1].plot(t, y[3], color=colors[i], linewidth=2, label=f'ND={nd} µM')
    axes[0, 1].set_xlabel('Time (min)', fontsize=11)
    axes[0, 1].set_ylabel('Misfolded MP (µM)', fontsize=11)
    axes[0, 1].set_title('Misfolded/Aggregated MP', fontsize=13, fontweight='bold')
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].grid(True, alpha=0.3)

    # Panel 3: Insertion efficiency (ratio)
    for i, nd in enumerate(nd_concs):
        y = results[nd].sol(t)
        total = y[2] + y[3] + y[1]
        eff = np.where(total > 0.01, np.divide(y[2], total, where=total > 0.01, out=np.zeros_like(total)) * 100, 0)
        axes[0, 2].plot(t, eff, color=colors[i], linewidth=2, label=f'ND={nd} µM')
    axes[0, 2].set_xlabel('Time (min)', fontsize=11)
    axes[0, 2].set_ylabel('Insertion Efficiency (%)', fontsize=11)
    axes[0, 2].set_title('Co-translational Insertion Efficiency', fontsize=13, fontweight='bold')
    axes[0, 2].legend(fontsize=9)
    axes[0, 2].grid(True, alpha=0.3)

    # Panel 4: Final yields bar chart
    final_inserted = []
    final_misfolded = []
    for nd in nd_concs:
        y = results[nd].sol(t)
        final_inserted.append(y[2, -1])
        final_misfolded.append(y[3, -1])

    x_pos = np.arange(len(nd_concs))
    axes[1, 0].bar(x_pos - 0.15, final_inserted, 0.3, label='Inserted (folded)',
                   color='#4CAF50', edgecolor='black')
    axes[1, 0].bar(x_pos + 0.15, final_misfolded, 0.3, label='Misfolded',
                   color='#F44336', edgecolor='black')
    axes[1, 0].set_xticks(x_pos)
    axes[1, 0].set_xticklabels([f'{nd}' for nd in nd_concs])
    axes[1, 0].set_xlabel('[Nanodisc] (µM)', fontsize=11)
    axes[1, 0].set_ylabel('Final Yield (µM)', fontsize=11)
    axes[1, 0].set_title('Final Yield vs Nanodisc Concentration', fontsize=13, fontweight='bold')
    axes[1, 0].legend(fontsize=10)
    axes[1, 0].grid(True, alpha=0.3, axis='y')

    # Panel 5: Free nanodiscs remaining
    for i, nd in enumerate(nd_concs):
        y = results[nd].sol(t)
        axes[1, 1].plot(t, y[4], color=colors[i], linewidth=2, label=f'ND={nd} µM')
    axes[1, 1].set_xlabel('Time (min)', fontsize=11)
    axes[1, 1].set_ylabel('Free Nanodiscs (µM)', fontsize=11)
    axes[1, 1].set_title('Nanodisc Consumption', fontsize=13, fontweight='bold')
    axes[1, 1].legend(fontsize=9)
    axes[1, 1].grid(True, alpha=0.3)

    # Panel 6: Cost-benefit (yield per nanodisc used)
    yield_per_nd = []
    for nd in nd_concs:
        y = results[nd].sol(t)
        used_nd = nd - y[4, -1]
        yld = y[2, -1]
        yield_per_nd.append(yld / max(used_nd, 0.01))

    axes[1, 2].bar(x_pos, yield_per_nd, 0.5, color='#FF9800', edgecolor='black')
    axes[1, 2].set_xticks(x_pos)
    axes[1, 2].set_xticklabels([f'{nd}' for nd in nd_concs])
    axes[1, 2].set_xlabel('[Nanodisc] (µM)', fontsize=11)
    axes[1, 2].set_ylabel('Yield per Nanodisc Used', fontsize=11)
    axes[1, 2].set_title('Nanodisc Utilization Efficiency', fontsize=13, fontweight='bold')
    axes[1, 2].grid(True, alpha=0.3, axis='y')

    plt.suptitle('Membrane Protein Expression in Nanodiscs: Case Study',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

def compute_mp_metrics(results):
    t = np.linspace(0, 360, 1000)
    metrics = {}
    nd_concs = sorted(results.keys())
    for nd in nd_concs:
        y = results[nd].sol(t)
        inserted = y[2, -1]
        misfolded = y[3, -1]
        total = inserted + misfolded + y[1, -1]
        metrics[f'ND_{nd}_uM'] = {
            'final_inserted_uM': float(inserted),
            'final_misfolded_uM': float(misfolded),
            'insertion_efficiency_pct': float(inserted / max(total, 0.01) * 100),
            'nanodisc_used_uM': float(nd - y[4, -1]),
            'yield_per_nanodisc': float(inserted / max(nd - y[4, -1], 0.01)),
        }
    return metrics

if __name__ == '__main__':
    print("=== Module 6: Membrane Protein / Nanodisc Case Study ===")
    results = run_nanodisc_titration()
    plot_nanodisc_case_study(results)
    metrics = compute_mp_metrics(results)

    for key, m in metrics.items():
        print(f"  {key}: inserted={m['final_inserted_uM']:.2f} µM, "
              f"efficiency={m['insertion_efficiency_pct']:.1f}%")

    with open('results/m6_nanodisc_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    print("  Saved: results/m6_nanodisc_metrics.json")
