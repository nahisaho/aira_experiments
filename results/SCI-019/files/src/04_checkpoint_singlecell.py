#!/usr/bin/env python3
"""
Module 4: Immune checkpoint molecule expression - single-cell level analysis.
Simulates scRNA-seq-like data for checkpoint molecules in RA.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from scipy import stats
import os

np.random.seed(42)

CHECKPOINTS = ['PD-1', 'PD-L1', 'CTLA-4', 'LAG-3', 'TIM-3', 'TIGIT', 'VISTA', 'ICOS', 'CD28', 'BTLA']
CELL_TYPES = ['CD4+ T', 'CD8+ T', 'Treg', 'Th17', 'B cell', 'NK', 'Monocyte', 'DC']

def simulate_singlecell_data(n_cells=5000):
    """Simulate single-cell checkpoint expression data."""
    n_ra = n_cells // 2
    n_hc = n_cells - n_ra
    
    cells = []
    
    # Cell type proportions
    ct_props_ra = [0.15, 0.12, 0.05, 0.10, 0.12, 0.08, 0.25, 0.13]
    ct_props_hc = [0.18, 0.15, 0.10, 0.03, 0.10, 0.12, 0.20, 0.12]
    
    for condition, n, props in [('RA', n_ra, ct_props_ra), ('HC', n_hc, ct_props_hc)]:
        ct_assignments = np.random.choice(CELL_TYPES, size=n, p=props)
        
        for i in range(n):
            ct = ct_assignments[i]
            expr = {}
            
            # Base expression varies by cell type
            base = {cp: np.random.exponential(0.5) for cp in CHECKPOINTS}
            
            # Cell-type-specific patterns
            if ct in ['CD4+ T', 'CD8+ T']:
                base['PD-1'] = np.random.exponential(1.5)
                base['CTLA-4'] = np.random.exponential(1.0)
                base['CD28'] = np.random.exponential(2.0)
            if ct == 'Treg':
                base['CTLA-4'] = np.random.exponential(3.0)
                base['PD-1'] = np.random.exponential(1.2)
                base['TIGIT'] = np.random.exponential(2.0)
                base['ICOS'] = np.random.exponential(2.5)
            if ct == 'Th17':
                base['PD-1'] = np.random.exponential(0.8)
                base['ICOS'] = np.random.exponential(1.5)
            if ct == 'NK':
                base['TIGIT'] = np.random.exponential(2.0)
                base['TIM-3'] = np.random.exponential(1.5)
                base['LAG-3'] = np.random.exponential(1.0)
            if ct == 'Monocyte':
                base['PD-L1'] = np.random.exponential(1.5)
                base['VISTA'] = np.random.exponential(2.0)
                base['TIM-3'] = np.random.exponential(1.2)
            if ct == 'DC':
                base['PD-L1'] = np.random.exponential(2.0)
                base['CD28'] = np.random.exponential(0.3)
                base['VISTA'] = np.random.exponential(1.5)
            
            # RA upregulation of inhibitory checkpoints
            if condition == 'RA':
                for cp in ['PD-1', 'LAG-3', 'TIM-3', 'TIGIT', 'VISTA']:
                    base[cp] *= np.random.uniform(1.3, 2.0)
                base['CTLA-4'] *= np.random.uniform(1.2, 1.5)
                base['CD28'] *= np.random.uniform(0.6, 0.9)  # Downregulated
            
            expr.update(base)
            expr['CellType'] = ct
            expr['Condition'] = condition
            cells.append(expr)
    
    return pd.DataFrame(cells)

def plot_tsne_checkpoint(df):
    """t-SNE visualization of checkpoint expression."""
    expr_data = df[CHECKPOINTS].values
    
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    coords = tsne.fit_transform(expr_data)
    
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    
    # By condition
    for cond, color in [('RA', '#E74C3C'), ('HC', '#3498DB')]:
        mask = df['Condition'] == cond
        axes[0].scatter(coords[mask, 0], coords[mask, 1], c=color, alpha=0.3, s=5, label=cond)
    axes[0].set_title('t-SNE by Condition', fontsize=13)
    axes[0].legend(markerscale=5)
    axes[0].set_xlabel('t-SNE 1')
    axes[0].set_ylabel('t-SNE 2')
    
    # By cell type
    ct_colors = plt.cm.Set3(np.linspace(0, 1, len(CELL_TYPES)))
    for ct, color in zip(CELL_TYPES, ct_colors):
        mask = df['CellType'] == ct
        axes[1].scatter(coords[mask, 0], coords[mask, 1], c=[color], alpha=0.4, s=5, label=ct)
    axes[1].set_title('t-SNE by Cell Type', fontsize=13)
    axes[1].legend(markerscale=5, fontsize=8, loc='upper right')
    axes[1].set_xlabel('t-SNE 1')
    axes[1].set_ylabel('t-SNE 2')
    
    # PD-1 expression overlay
    pd1 = df['PD-1'].values
    sc = axes[2].scatter(coords[:, 0], coords[:, 1], c=pd1, cmap='Reds', alpha=0.4, s=5, vmax=np.percentile(pd1, 95))
    axes[2].set_title('PD-1 Expression', fontsize=13)
    plt.colorbar(sc, ax=axes[2], label='Expression')
    axes[2].set_xlabel('t-SNE 1')
    axes[2].set_ylabel('t-SNE 2')
    
    plt.suptitle('Single-Cell Immune Checkpoint Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/checkpoint_tsne.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_checkpoint_dotplot(df):
    """Dot plot of checkpoint expression by cell type and condition."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    for ax, cond in zip(axes, ['RA', 'HC']):
        sub = df[df['Condition'] == cond]
        
        mean_expr = sub.groupby('CellType')[CHECKPOINTS].mean()
        pct_expr = sub.groupby('CellType')[CHECKPOINTS].apply(lambda x: (x > 0.5).mean())
        
        for i, ct in enumerate(CELL_TYPES):
            for j, cp in enumerate(CHECKPOINTS):
                size = pct_expr.loc[ct, cp] * 200
                color = mean_expr.loc[ct, cp]
                ax.scatter(j, i, s=size, c=color, cmap='Reds', vmin=0, vmax=3,
                          edgecolors='black', linewidth=0.5)
        
        ax.set_xticks(range(len(CHECKPOINTS)))
        ax.set_xticklabels(CHECKPOINTS, rotation=45, ha='right')
        ax.set_yticks(range(len(CELL_TYPES)))
        ax.set_yticklabels(CELL_TYPES)
        ax.set_title(f'{cond} - Checkpoint Expression', fontsize=13)
        ax.set_xlim(-0.5, len(CHECKPOINTS) - 0.5)
        ax.set_ylim(-0.5, len(CELL_TYPES) - 0.5)
        ax.grid(True, alpha=0.2)
    
    plt.suptitle('Checkpoint Expression Dot Plot (size=% expressing, color=mean)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/checkpoint_dotplot.png', dpi=150, bbox_inches='tight')
    plt.close()

def compute_checkpoint_stats(df):
    """Compute RA vs HC differential checkpoint expression per cell type."""
    results = []
    for ct in CELL_TYPES:
        for cp in CHECKPOINTS:
            ra_vals = df[(df['Condition'] == 'RA') & (df['CellType'] == ct)][cp]
            hc_vals = df[(df['Condition'] == 'HC') & (df['CellType'] == ct)][cp]
            if len(ra_vals) > 5 and len(hc_vals) > 5:
                t_stat, p_val = stats.mannwhitneyu(ra_vals, hc_vals, alternative='two-sided')
                results.append({
                    'CellType': ct, 'Checkpoint': cp,
                    'RA_mean': ra_vals.mean(), 'HC_mean': hc_vals.mean(),
                    'log2FC': np.log2((ra_vals.mean() + 0.01) / (hc_vals.mean() + 0.01)),
                    'p_value': p_val
                })
    return pd.DataFrame(results)

if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    
    print("Simulating single-cell checkpoint data (5000 cells)...")
    df = simulate_singlecell_data(5000)
    
    print("Creating t-SNE visualization...")
    plot_tsne_checkpoint(df)
    
    print("Creating dot plot...")
    plot_checkpoint_dotplot(df)
    
    print("Computing differential expression stats...")
    stats_df = compute_checkpoint_stats(df)
    sig = stats_df[stats_df['p_value'] < 0.05]
    print(f"  Significant checkpoint-cell type pairs: {len(sig)}/{len(stats_df)}")
    
    stats_df.to_csv('src/checkpoint_stats.csv', index=False)
    print("Module 4 complete.")
