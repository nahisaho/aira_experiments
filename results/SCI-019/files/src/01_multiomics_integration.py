#!/usr/bin/env python3
"""
Module 1: Multi-omics data integration (Transcriptome, Proteome, Metabolome)
Simulates and integrates multi-omics data for autoimmune disease analysis.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import CCA
from scipy import stats
import os

np.random.seed(42)

N_SAMPLES = 120  # 60 RA patients, 60 healthy controls
N_GENES = 500
N_PROTEINS = 200
N_METABOLITES = 150

def generate_multiomics_data():
    """Generate synthetic multi-omics data with disease-associated signals."""
    labels = np.array(['RA'] * 60 + ['HC'] * 60)
    
    # Transcriptome: log2 expression
    base_expr = np.random.normal(8, 2, (N_SAMPLES, N_GENES))
    # Inject RA-associated DE genes (first 50 genes upregulated, next 30 downregulated)
    base_expr[:60, :50] += np.random.normal(2.5, 0.5, (60, 50))
    base_expr[:60, 50:80] -= np.random.normal(1.8, 0.4, (60, 30))
    
    gene_names = [f"Gene_{i}" for i in range(N_GENES)]
    # Key immune genes
    immune_genes = ['TNF', 'IL6', 'IL1B', 'IL17A', 'IL10', 'IFNG', 'TGFB1', 'IL23A',
                    'CTLA4', 'PDCD1', 'LAG3', 'HAVCR2', 'CD274', 'ICOS', 'CD28', 'FOXP3']
    for i, g in enumerate(immune_genes):
        gene_names[i] = g
    
    transcriptome = pd.DataFrame(base_expr, columns=gene_names)
    transcriptome['Group'] = labels
    
    # Proteome
    base_prot = np.random.normal(6, 1.5, (N_SAMPLES, N_PROTEINS))
    base_prot[:60, :30] += np.random.normal(1.8, 0.3, (60, 30))
    base_prot[:60, 30:50] -= np.random.normal(1.2, 0.3, (60, 20))
    
    prot_names = [f"Protein_{i}" for i in range(N_PROTEINS)]
    key_prots = ['CRP', 'SAA', 'MMP3', 'VEGF', 'IL6_prot', 'TNF_prot', 'RF', 'ACPA']
    for i, p in enumerate(key_prots):
        prot_names[i] = p
    
    proteome = pd.DataFrame(base_prot, columns=prot_names)
    proteome['Group'] = labels
    
    # Metabolome
    base_met = np.random.normal(5, 1, (N_SAMPLES, N_METABOLITES))
    base_met[:60, :20] += np.random.normal(1.5, 0.3, (60, 20))
    base_met[:60, 20:35] -= np.random.normal(1.0, 0.2, (60, 15))
    
    met_names = [f"Metabolite_{i}" for i in range(N_METABOLITES)]
    key_mets = ['Tryptophan', 'Kynurenine', 'Lactate', 'Succinate', 'Itaconate', 'PGE2']
    for i, m in enumerate(key_mets):
        met_names[i] = m
    
    metabolome = pd.DataFrame(base_met, columns=met_names)
    metabolome['Group'] = labels
    
    return transcriptome, proteome, metabolome, labels

def perform_pca_integration(transcriptome, proteome, metabolome, labels):
    """PCA-based integration of multi-omics data."""
    scaler = StandardScaler()
    
    t_data = scaler.fit_transform(transcriptome.drop('Group', axis=1))
    p_data = scaler.fit_transform(proteome.drop('Group', axis=1))
    m_data = scaler.fit_transform(metabolome.drop('Group', axis=1))
    
    # Concatenate and perform joint PCA
    combined = np.hstack([t_data, p_data, m_data])
    pca = PCA(n_components=10)
    pca_result = pca.fit_transform(combined)
    
    return pca, pca_result

def plot_multiomics_pca(pca_result, labels, pca):
    """Plot PCA of integrated multi-omics data."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    colors = {'RA': '#E74C3C', 'HC': '#3498DB'}
    
    # PC1 vs PC2
    for group in ['RA', 'HC']:
        mask = labels == group
        axes[0].scatter(pca_result[mask, 0], pca_result[mask, 1],
                       c=colors[group], label=group, alpha=0.7, s=50, edgecolors='white')
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    axes[0].set_title('Multi-omics PCA: PC1 vs PC2')
    axes[0].legend()
    
    # PC2 vs PC3
    for group in ['RA', 'HC']:
        mask = labels == group
        axes[1].scatter(pca_result[mask, 1], pca_result[mask, 2],
                       c=colors[group], label=group, alpha=0.7, s=50, edgecolors='white')
    axes[1].set_xlabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    axes[1].set_ylabel(f'PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)')
    axes[1].set_title('Multi-omics PCA: PC2 vs PC3')
    axes[1].legend()
    
    # Variance explained
    axes[2].bar(range(1, 11), pca.explained_variance_ratio_ * 100, color='#2ECC71', edgecolor='black')
    axes[2].set_xlabel('Principal Component')
    axes[2].set_ylabel('Variance Explained (%)')
    axes[2].set_title('Scree Plot')
    
    plt.tight_layout()
    plt.savefig('figures/multiomics_pca.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_de_volcano(transcriptome, labels):
    """Volcano plot of differentially expressed genes."""
    ra = transcriptome[labels == 'RA'].drop('Group', axis=1)
    hc = transcriptome[labels == 'HC'].drop('Group', axis=1)
    
    log2fc = ra.mean() - hc.mean()
    pvals = []
    for col in ra.columns:
        _, p = stats.ttest_ind(ra[col], hc[col])
        pvals.append(p)
    pvals = np.array(pvals)
    neg_log10p = -np.log10(pvals + 1e-300)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    sig_up = (log2fc > 1) & (pvals < 0.05)
    sig_down = (log2fc < -1) & (pvals < 0.05)
    ns = ~sig_up & ~sig_down
    
    ax.scatter(log2fc[ns], neg_log10p[ns], c='gray', alpha=0.3, s=20, label='NS')
    ax.scatter(log2fc[sig_up], neg_log10p[sig_up], c='#E74C3C', alpha=0.7, s=30, label=f'Up ({sig_up.sum()})')
    ax.scatter(log2fc[sig_down], neg_log10p[sig_down], c='#3498DB', alpha=0.7, s=30, label=f'Down ({sig_down.sum()})')
    
    # Label top genes
    top_genes = log2fc.abs().nlargest(8).index
    for gene in top_genes:
        idx = list(transcriptome.columns).index(gene)
        ax.annotate(gene, (log2fc[gene], neg_log10p[idx]), fontsize=8, fontweight='bold')
    
    ax.axhline(-np.log10(0.05), ls='--', c='gray', alpha=0.5)
    ax.axvline(1, ls='--', c='gray', alpha=0.5)
    ax.axvline(-1, ls='--', c='gray', alpha=0.5)
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Differential Expression: RA vs Healthy Controls')
    ax.legend()
    plt.tight_layout()
    plt.savefig('figures/volcano_plot.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return log2fc, pvals

def plot_omics_correlation_heatmap(transcriptome, proteome, metabolome):
    """Cross-omics correlation heatmap for key features."""
    key_t = ['TNF', 'IL6', 'IL1B', 'IL17A', 'IL10', 'IFNG', 'CTLA4', 'PDCD1']
    key_p = ['CRP', 'SAA', 'MMP3', 'VEGF', 'IL6_prot', 'TNF_prot', 'RF', 'ACPA']
    key_m = ['Tryptophan', 'Kynurenine', 'Lactate', 'Succinate', 'Itaconate', 'PGE2']
    
    combined = pd.concat([
        transcriptome[key_t],
        proteome[key_p],
        metabolome[key_m]
    ], axis=1)
    
    corr = combined.corr()
    
    fig, ax = plt.subplots(figsize=(14, 11))
    mask = np.triu(np.ones_like(corr, dtype=bool))
    sns.heatmap(corr, mask=mask, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                square=True, linewidths=0.5, ax=ax, annot=True, fmt='.2f', annot_kws={'size': 7})
    ax.set_title('Cross-omics Correlation Heatmap\n(Transcriptome | Proteome | Metabolome)', fontsize=14)
    plt.tight_layout()
    plt.savefig('figures/cross_omics_correlation.png', dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    
    print("Generating multi-omics data...")
    transcriptome, proteome, metabolome, labels = generate_multiomics_data()
    
    print("Performing PCA integration...")
    pca, pca_result = perform_pca_integration(transcriptome, proteome, metabolome, labels)
    print(f"  Top 3 PCs explain: {sum(pca.explained_variance_ratio_[:3])*100:.1f}% variance")
    
    print("Creating volcano plot...")
    log2fc, pvals = plot_de_volcano(transcriptome, labels)
    sig_genes = ((log2fc.abs() > 1) & (pvals < 0.05)).sum()
    print(f"  Significant DE genes: {sig_genes}")
    
    print("Creating PCA plot...")
    plot_multiomics_pca(pca_result, labels, pca)
    
    print("Creating correlation heatmap...")
    plot_omics_correlation_heatmap(transcriptome, proteome, metabolome)
    
    # Save data for downstream modules
    transcriptome.to_csv('src/transcriptome_data.csv', index=False)
    proteome.to_csv('src/proteome_data.csv', index=False)
    metabolome.to_csv('src/metabolome_data.csv', index=False)
    
    print("Module 1 complete.")
