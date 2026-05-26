#!/usr/bin/env python3
"""
Module 2: Immune cell subset deconvolution (CIBERSORTx-style analysis)
Simulates CIBERSORTx immune deconvolution for RA synovial tissue.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

np.random.seed(42)

CELL_TYPES = [
    'Naive CD4+ T', 'Memory CD4+ T', 'Th1', 'Th17', 'Treg',
    'Naive CD8+ T', 'Cytotoxic CD8+ T', 'NK cells',
    'Naive B', 'Memory B', 'Plasma cells',
    'Monocytes CD14+', 'Monocytes CD16+', 'M1 Macrophages', 'M2 Macrophages',
    'Dendritic cells', 'Neutrophils', 'Mast cells',
    'Fibroblast-like Synoviocytes', 'Endothelial cells'
]

def simulate_deconvolution(n_ra=60, n_hc=60):
    """Simulate CIBERSORTx-like cell fraction estimates."""
    n_total = n_ra + n_hc
    labels = np.array(['RA'] * n_ra + ['HC'] * n_hc)
    
    # Base proportions (healthy)
    base = {
        'Naive CD4+ T': 0.08, 'Memory CD4+ T': 0.06, 'Th1': 0.03, 'Th17': 0.02, 'Treg': 0.04,
        'Naive CD8+ T': 0.05, 'Cytotoxic CD8+ T': 0.04, 'NK cells': 0.06,
        'Naive B': 0.04, 'Memory B': 0.03, 'Plasma cells': 0.02,
        'Monocytes CD14+': 0.10, 'Monocytes CD16+': 0.04, 'M1 Macrophages': 0.05, 'M2 Macrophages': 0.06,
        'Dendritic cells': 0.03, 'Neutrophils': 0.15, 'Mast cells': 0.02,
        'Fibroblast-like Synoviocytes': 0.05, 'Endothelial cells': 0.03
    }
    
    # RA-specific shifts
    ra_shift = {
        'Th1': 0.04, 'Th17': 0.05, 'Treg': -0.02,
        'Memory CD4+ T': 0.03, 'Plasma cells': 0.04,
        'M1 Macrophages': 0.06, 'M2 Macrophages': -0.02,
        'Monocytes CD14+': 0.03, 'Neutrophils': 0.05,
        'Fibroblast-like Synoviocytes': 0.04,
        'Naive CD4+ T': -0.03, 'NK cells': -0.02, 'Naive B': -0.01
    }
    
    fractions = np.zeros((n_total, len(CELL_TYPES)))
    for i in range(n_total):
        for j, ct in enumerate(CELL_TYPES):
            val = base[ct]
            if i < n_ra and ct in ra_shift:
                val += ra_shift[ct]
            val += np.random.normal(0, 0.008)
            fractions[i, j] = max(val, 0.001)
        fractions[i] /= fractions[i].sum()
    
    df = pd.DataFrame(fractions, columns=CELL_TYPES)
    df['Group'] = labels
    return df

def plot_cell_proportions(df):
    """Stacked bar plot of cell type proportions."""
    fig, ax = plt.subplots(figsize=(16, 6))
    
    ra_mean = df[df['Group'] == 'RA'][CELL_TYPES].mean()
    hc_mean = df[df['Group'] == 'HC'][CELL_TYPES].mean()
    
    x = np.arange(len(CELL_TYPES))
    width = 0.35
    
    ax.bar(x - width/2, ra_mean, width, label='RA', color='#E74C3C', alpha=0.8)
    ax.bar(x + width/2, hc_mean, width, label='HC', color='#3498DB', alpha=0.8)
    
    ax.set_ylabel('Mean Fraction')
    ax.set_title('Immune Cell Deconvolution: RA vs Healthy Controls')
    ax.set_xticks(x)
    ax.set_xticklabels(CELL_TYPES, rotation=45, ha='right', fontsize=8)
    ax.legend()
    
    # Mark significant differences
    for j, ct in enumerate(CELL_TYPES):
        ra_vals = df[df['Group'] == 'RA'][ct]
        hc_vals = df[df['Group'] == 'HC'][ct]
        _, p = stats.ttest_ind(ra_vals, hc_vals)
        if p < 0.001:
            ax.text(j, max(ra_mean[ct], hc_mean[ct]) + 0.005, '***', ha='center', fontsize=8)
        elif p < 0.01:
            ax.text(j, max(ra_mean[ct], hc_mean[ct]) + 0.005, '**', ha='center', fontsize=8)
        elif p < 0.05:
            ax.text(j, max(ra_mean[ct], hc_mean[ct]) + 0.005, '*', ha='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('figures/immune_deconvolution.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_deconv_heatmap(df):
    """Heatmap of cell fractions per sample."""
    fig, ax = plt.subplots(figsize=(14, 8))
    
    sorted_df = df.sort_values('Group')
    data = sorted_df[CELL_TYPES].values.T
    groups = sorted_df['Group'].values
    
    sns.heatmap(data, cmap='YlOrRd', ax=ax, xticklabels=False,
                yticklabels=CELL_TYPES)
    ax.set_xlabel('Samples (RA | HC)')
    ax.set_title('Immune Cell Fraction Heatmap')
    
    # Add group annotation
    ra_end = (groups == 'RA').sum()
    ax.axvline(x=ra_end, color='white', linewidth=2)
    ax.text(ra_end/2, -0.5, 'RA', ha='center', fontsize=12, fontweight='bold', color='red')
    ax.text(ra_end + (len(groups)-ra_end)/2, -0.5, 'HC', ha='center', fontsize=12, fontweight='bold', color='blue')
    
    plt.tight_layout()
    plt.savefig('figures/deconvolution_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()

def compute_statistics(df):
    """Compute differential cell fractions."""
    results = []
    for ct in CELL_TYPES:
        ra_vals = df[df['Group'] == 'RA'][ct]
        hc_vals = df[df['Group'] == 'HC'][ct]
        t_stat, p_val = stats.ttest_ind(ra_vals, hc_vals)
        fc = ra_vals.mean() / hc_vals.mean()
        results.append({
            'Cell Type': ct,
            'RA Mean': ra_vals.mean(),
            'HC Mean': hc_vals.mean(),
            'Fold Change': fc,
            'log2FC': np.log2(fc),
            't-statistic': t_stat,
            'p-value': p_val,
            'Significant': p_val < 0.05
        })
    return pd.DataFrame(results)

if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    
    print("Simulating CIBERSORTx deconvolution...")
    df = simulate_deconvolution()
    
    print("Computing statistics...")
    stats_df = compute_statistics(df)
    sig = stats_df[stats_df['Significant']]
    print(f"  Significantly altered cell types: {len(sig)}")
    for _, row in sig.iterrows():
        direction = "↑" if row['log2FC'] > 0 else "↓"
        print(f"    {row['Cell Type']}: {direction} log2FC={row['log2FC']:.2f}, p={row['p-value']:.2e}")
    
    print("Creating cell proportion plot...")
    plot_cell_proportions(df)
    
    print("Creating heatmap...")
    plot_deconv_heatmap(df)
    
    stats_df.to_csv('src/deconvolution_stats.csv', index=False)
    df.to_csv('src/deconvolution_fractions.csv', index=False)
    
    print("Module 2 complete.")
