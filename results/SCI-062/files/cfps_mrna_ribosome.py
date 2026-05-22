"""
Module 4: mRNA Stability and Ribosome Loading Prediction Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

# --- mRNA Stability Parameters ---
mrna_features = {
    'GFP': {
        'length_nt': 720, 'gc_content': 0.55, 'has_5UTR_stem': True,
        'kozak_strength': 0.9, 'polyA_length': 50, 'rare_codons_frac': 0.05,
        'rbs_strength': 0.85, 'label': 'GFP (720 nt)',
    },
    'Luciferase': {
        'length_nt': 1653, 'gc_content': 0.48, 'has_5UTR_stem': False,
        'kozak_strength': 0.7, 'polyA_length': 80, 'rare_codons_frac': 0.12,
        'rbs_strength': 0.75, 'label': 'Luciferase (1653 nt)',
    },
    'MP_GPCR': {
        'length_nt': 1050, 'gc_content': 0.42, 'has_5UTR_stem': False,
        'kozak_strength': 0.6, 'polyA_length': 60, 'rare_codons_frac': 0.18,
        'rbs_strength': 0.65, 'label': 'GPCR (1050 nt)',
    },
    'scFv': {
        'length_nt': 750, 'gc_content': 0.52, 'has_5UTR_stem': True,
        'kozak_strength': 0.85, 'polyA_length': 40, 'rare_codons_frac': 0.08,
        'rbs_strength': 0.80, 'label': 'scFv (750 nt)',
    },
}

def predict_mrna_halflife(features):
    """Predict mRNA half-life based on sequence features (minutes)"""
    base_hl = 25.0  # base half-life in minutes

    # GC content effect (moderate GC is best)
    gc_factor = np.exp(-5 * (features['gc_content'] - 0.50)**2)

    # Length penalty (longer mRNAs degrade faster)
    length_factor = np.exp(-0.0003 * features['length_nt'])

    # 5' UTR structure protection
    utr_factor = 1.5 if features['has_5UTR_stem'] else 1.0

    # polyA tail stabilization
    polyA_factor = 1.0 + 0.005 * features['polyA_length']

    halflife = base_hl * gc_factor * length_factor * utr_factor * polyA_factor
    return halflife

def predict_ribosome_loading(features, ribosome_conc=500.0):
    """
    Predict ribosome occupancy and translation initiation rate.
    Returns dict with loading metrics.
    """
    # Initiation rate depends on RBS/Kozak strength and 5'UTR accessibility
    k_init = 0.5 * features['rbs_strength'] * features['kozak_strength']
    if features['has_5UTR_stem']:
        k_init *= 0.7  # structured 5'UTR reduces initiation

    # Elongation rate modulated by rare codons
    k_elong_base = 10.0  # codons per second
    k_elong = k_elong_base * (1 - 0.5 * features['rare_codons_frac'])

    # Ribosome spacing (minimum ~30 nt between ribosomes)
    footprint = 30  # nt
    max_ribosomes = features['length_nt'] / footprint

    # Steady-state ribosome density
    transit_time = features['length_nt'] / (3 * k_elong)  # seconds (3 nt per codon)
    ribo_per_mrna = min(k_init * transit_time, max_ribosomes)

    # Protein production rate
    protein_rate = ribo_per_mrna * k_elong * 3 / features['length_nt']  # proteins/sec/mRNA

    return {
        'k_init_per_sec': float(k_init),
        'k_elong_codons_per_sec': float(k_elong),
        'ribosomes_per_mRNA': float(ribo_per_mrna),
        'max_ribosomes': float(max_ribosomes),
        'transit_time_sec': float(transit_time),
        'protein_rate_per_sec_per_mRNA': float(protein_rate),
    }

def mrna_dynamics_ode(t, y, halflife, k_tx=0.5):
    """mRNA dynamics with synthesis and degradation"""
    mRNA = max(y[0], 0)
    k_deg = np.log(2) / halflife
    dmRNA_dt = k_tx - k_deg * mRNA
    return [dmRNA_dt]

def plot_mrna_stability(save_path='figures/fig7_mrna_stability.png'):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    colors = ['#2196F3', '#FF9800', '#E91E63', '#4CAF50']

    halflife_data = {}
    for i, (name, feat) in enumerate(mrna_features.items()):
        hl = predict_mrna_halflife(feat)
        halflife_data[name] = hl

        sol = solve_ivp(mrna_dynamics_ode, [0, 180], [0.0], args=(hl,),
                        dense_output=True, max_step=0.5, rtol=1e-8, atol=1e-10)
        t = np.linspace(0, 180, 300)
        y = sol.sol(t)
        axes[0].plot(t, y[0], color=colors[i], linewidth=2, label=f"{feat['label']} (t½={hl:.1f} min)")

    axes[0].set_xlabel('Time (min)', fontsize=11)
    axes[0].set_ylabel('[mRNA] (a.u.)', fontsize=11)
    axes[0].set_title('mRNA Accumulation Dynamics', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    # Ribosome loading comparison
    names = list(mrna_features.keys())
    ribo_data = {n: predict_ribosome_loading(mrna_features[n]) for n in names}

    ribo_per = [ribo_data[n]['ribosomes_per_mRNA'] for n in names]
    axes[1].barh(range(len(names)), ribo_per, color=colors, edgecolor='black', alpha=0.8)
    axes[1].set_yticks(range(len(names)))
    axes[1].set_yticklabels([mrna_features[n]['label'] for n in names], fontsize=10)
    axes[1].set_xlabel('Ribosomes per mRNA', fontsize=11)
    axes[1].set_title('Ribosome Loading Density', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='x')

    # Protein production rate
    prot_rate = [ribo_data[n]['protein_rate_per_sec_per_mRNA'] * 60 for n in names]  # per min
    axes[2].barh(range(len(names)), prot_rate, color=colors, edgecolor='black', alpha=0.8)
    axes[2].set_yticks(range(len(names)))
    axes[2].set_yticklabels([mrna_features[n]['label'] for n in names], fontsize=10)
    axes[2].set_xlabel('Protein Production Rate (proteins/min/mRNA)', fontsize=11)
    axes[2].set_title('Translation Efficiency', fontsize=13, fontweight='bold')
    axes[2].grid(True, alpha=0.3, axis='x')

    plt.suptitle('mRNA Stability & Ribosome Loading Prediction',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

    return halflife_data, ribo_data

def plot_codon_sensitivity(save_path='figures/fig8_codon_sensitivity.png'):
    """Sensitivity analysis: rare codon fraction vs productivity"""
    rare_fracs = np.linspace(0, 0.35, 50)
    fig, ax = plt.subplots(figsize=(8, 5))

    for name, feat in mrna_features.items():
        rates = []
        for rf in rare_fracs:
            f = feat.copy()
            f['rare_codons_frac'] = rf
            rd = predict_ribosome_loading(f)
            rates.append(rd['protein_rate_per_sec_per_mRNA'] * 60)
        ax.plot(rare_fracs * 100, rates, linewidth=2, label=feat['label'])

    ax.set_xlabel('Rare Codon Fraction (%)', fontsize=11)
    ax.set_ylabel('Protein Rate (proteins/min/mRNA)', fontsize=11)
    ax.set_title('Effect of Codon Usage on Translation Efficiency', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

if __name__ == '__main__':
    print("=== Module 4: mRNA Stability & Ribosome Loading ===")
    halflife_data, ribo_data = plot_mrna_stability()
    plot_codon_sensitivity()

    results = {}
    for name in mrna_features:
        results[name] = {
            'halflife_min': halflife_data[name],
            **ribo_data[name],
        }

    for name, r in results.items():
        print(f"  {name}: t½={r['halflife_min']:.1f} min, "
              f"ribo/mRNA={r['ribosomes_per_mRNA']:.1f}, "
              f"rate={r['protein_rate_per_sec_per_mRNA']*60:.3f} prot/min/mRNA")

    with open('results/m4_mrna_ribosome.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("  Saved: results/m4_mrna_ribosome.json")
