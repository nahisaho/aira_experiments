"""
Module 1: Transcription-Translation Coupled Model with Resource Competition
ODE-based model for CFPS system dynamics
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os

# --- Parameters ---
params = {
    # Transcription
    'k_tx': 0.8,          # transcription rate (nM/min)
    'K_NTP_tx': 0.5,      # Km for NTP in transcription (mM)
    'k_deg_mRNA': 0.03,   # mRNA degradation rate (1/min)
    # Translation
    'k_tl': 2.0,          # translation rate (nM/min/nM_mRNA)
    'K_aa': 0.1,          # Km for amino acids (mM)
    'K_NTP_tl': 0.3,      # Km for NTP in translation (mM - GTP for elongation)
    'k_deg_protein': 0.001, # protein degradation (1/min)
    # Resource competition
    'K_ribo': 50.0,       # ribosome Km (nM)
    'R_total': 500.0,     # total ribosome concentration (nM)
    'RNAP_total': 100.0,  # total RNAP concentration (nM)
    'K_RNAP': 20.0,       # RNAP Km (nM)
    # Energy
    'k_NTP_consume_tx': 0.05,  # NTP consumption per transcription event
    'k_NTP_consume_tl': 0.15,  # NTP consumption per translation event
    'k_aa_consume': 0.08,      # amino acid consumption rate
    # DNA template
    'DNA': 10.0,          # DNA template concentration (nM)
}

def cfps_odes(t, y, p):
    mRNA, protein, NTP, aa = y
    mRNA = max(mRNA, 0)
    protein = max(protein, 0)
    NTP = max(NTP, 0)
    aa = max(aa, 0)

    # Resource competition: free ribosomes
    ribo_free = p['R_total'] / (1 + mRNA / p['K_ribo'])
    rnap_free = p['RNAP_total'] / (1 + p['DNA'] / p['K_RNAP'])

    # Transcription rate (Michaelis-Menten for NTP)
    v_tx = p['k_tx'] * rnap_free * (p['DNA'] / (p['K_RNAP'] + p['DNA'])) * (NTP / (p['K_NTP_tx'] + NTP))

    # Translation rate (depends on mRNA, ribosomes, amino acids, GTP)
    v_tl = p['k_tl'] * mRNA * (ribo_free / (p['K_ribo'] + ribo_free)) * \
           (aa / (p['K_aa'] + aa)) * (NTP / (p['K_NTP_tl'] + NTP))

    dmRNA_dt = v_tx - p['k_deg_mRNA'] * mRNA
    dprotein_dt = v_tl - p['k_deg_protein'] * protein
    dNTP_dt = -p['k_NTP_consume_tx'] * v_tx - p['k_NTP_consume_tl'] * v_tl
    daa_dt = -p['k_aa_consume'] * v_tl

    return [dmRNA_dt, dprotein_dt, dNTP_dt, daa_dt]

def run_simulation(p=None, t_span=(0, 300), y0=None):
    if p is None:
        p = params.copy()
    if y0 is None:
        y0 = [0.0, 0.0, 8.0, 5.0]  # [mRNA, protein, NTP(mM), aa(mM)]

    sol = solve_ivp(cfps_odes, t_span, y0, args=(p,),
                    method='RK45', dense_output=True,
                    max_step=0.5, rtol=1e-8, atol=1e-10)
    return sol

def plot_dynamics(sol, save_path='figures/fig1_txn_tln_dynamics.png'):
    t = np.linspace(sol.t[0], sol.t[-1], 500)
    y = sol.sol(t)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    labels = ['mRNA (nM)', 'Protein (nM)', 'NTP (mM)', 'Amino Acids (mM)']
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0']

    for i, ax in enumerate(axes.flat):
        ax.plot(t, y[i], color=colors[i], linewidth=2)
        ax.set_xlabel('Time (min)', fontsize=11)
        ax.set_ylabel(labels[i], fontsize=11)
        ax.set_title(labels[i].split('(')[0].strip(), fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(t[0], t[-1])

    fig.suptitle('CFPS Transcription-Translation Coupled Dynamics\nwith Resource Competition',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

def resource_competition_analysis(save_path='figures/fig2_resource_competition.png'):
    """Analyze effect of DNA template concentration on resource competition"""
    dna_range = np.linspace(1, 100, 20)
    final_protein = []
    final_mRNA = []
    peak_rate = []

    for dna in dna_range:
        p = params.copy()
        p['DNA'] = dna
        sol = run_simulation(p)
        t_dense = np.linspace(0, 300, 500)
        y_dense = sol.sol(t_dense)
        final_protein.append(y_dense[1, -1])
        final_mRNA.append(np.max(y_dense[0]))
        # protein production rate
        dt = t_dense[1] - t_dense[0]
        rates = np.diff(y_dense[1]) / dt
        peak_rate.append(np.max(rates))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].plot(dna_range, final_protein, 'o-', color='#4CAF50', linewidth=2)
    axes[0].set_xlabel('[DNA template] (nM)', fontsize=11)
    axes[0].set_ylabel('Final Protein Yield (nM)', fontsize=11)
    axes[0].set_title('Protein Yield vs DNA Loading', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(dna_range, final_mRNA, 's-', color='#2196F3', linewidth=2)
    axes[1].set_xlabel('[DNA template] (nM)', fontsize=11)
    axes[1].set_ylabel('Peak mRNA (nM)', fontsize=11)
    axes[1].set_title('Peak mRNA vs DNA Loading', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(dna_range, peak_rate, '^-', color='#FF5722', linewidth=2)
    axes[2].set_xlabel('[DNA template] (nM)', fontsize=11)
    axes[2].set_ylabel('Peak Translation Rate (nM/min)', fontsize=11)
    axes[2].set_title('Translation Rate Saturation', fontsize=13, fontweight='bold')
    axes[2].grid(True, alpha=0.3)

    plt.suptitle('Resource Competition Analysis: DNA Template Titration',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

    return {
        'optimal_DNA_nM': float(dna_range[np.argmax(final_protein)]),
        'max_protein_yield_nM': float(max(final_protein)),
        'saturation_DNA_nM': float(dna_range[np.argmax(peak_rate)])
    }

if __name__ == '__main__':
    print("=== Module 1: Transcription-Translation Coupled Model ===")
    sol = run_simulation()
    plot_dynamics(sol)
    results = resource_competition_analysis()
    print(f"  Optimal DNA: {results['optimal_DNA_nM']:.1f} nM")
    print(f"  Max protein yield: {results['max_protein_yield_nM']:.1f} nM")

    with open('results/m1_txn_tln_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("  Saved: results/m1_txn_tln_results.json")
