#!/usr/bin/env python3
"""
Spatial Transcriptomics Advanced Analysis Pipeline
===================================================
Implements:
1. Spot deconvolution (cell type composition estimation)
2. Spatially variable gene detection (SpatialDE-inspired)
3. Cell-cell communication (ligand-receptor analysis)
4. Tissue microenvironment niche identification
5. 3D spatial reconstruction (serial section integration)
6. Tumor immune microenvironment case study
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import sparse
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import pearsonr, spearmanr, norm, ttest_ind
from scipy.optimize import nnls
from sklearn.decomposition import NMF, PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import normalize
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score, adjusted_rand_score
import anndata as ad
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
FIGDIR = 'figures'

# ============================================================
# 1. SYNTHETIC DATA GENERATION
# ============================================================

def generate_synthetic_spatial_data(n_spots=2000, n_genes=500, n_cell_types=6):
    """Generate synthetic Visium-like spatial transcriptomics data with
    known ground truth cell type compositions and spatial patterns."""

    # Spatial coordinates on a grid-like layout (Visium hexagonal approx)
    side = int(np.ceil(np.sqrt(n_spots)))
    coords = []
    for i in range(side):
        for j in range(side):
            x = j + (i % 2) * 0.5  # hexagonal offset
            y = i * np.sqrt(3) / 2
            coords.append([x, y])
    coords = np.array(coords[:n_spots], dtype=float)

    # Normalize coordinates to [0, 1]
    for d in range(2):
        mn, mx = coords[:, d].min(), coords[:, d].max()
        if mx > mn:
            coords[:, d] = (coords[:, d] - mn) / (mx - mn)
        else:
            coords[:, d] = 0.5

    cell_type_names = ['Tumor_Epithelial', 'CD8_T_cell', 'Macrophage',
                       'Fibroblast', 'B_cell', 'Endothelial']

    # Ground truth cell type proportions with spatial structure
    proportions = np.zeros((n_spots, n_cell_types))

    # Create spatial domains
    center = np.array([0.5, 0.5])
    dist_to_center = np.sqrt(np.sum((coords - center)**2, axis=1))

    # Tumor core (center)
    tumor_mask = dist_to_center < 0.25
    proportions[tumor_mask, 0] = 0.7  # Tumor
    proportions[tumor_mask, 3] = 0.2  # Fibroblast (CAF)
    proportions[tumor_mask, 5] = 0.1  # Endothelial

    # Immune infiltrate (ring around tumor)
    immune_mask = (dist_to_center >= 0.25) & (dist_to_center < 0.4)
    proportions[immune_mask, 0] = 0.3  # Tumor
    proportions[immune_mask, 1] = 0.3  # CD8 T cells
    proportions[immune_mask, 2] = 0.2  # Macrophage
    proportions[immune_mask, 4] = 0.1  # B cells
    proportions[immune_mask, 5] = 0.1  # Endothelial

    # Stromal region (outer)
    stromal_mask = dist_to_center >= 0.4
    proportions[stromal_mask, 3] = 0.5  # Fibroblast
    proportions[stromal_mask, 5] = 0.2  # Endothelial
    proportions[stromal_mask, 2] = 0.15  # Macrophage
    proportions[stromal_mask, 4] = 0.15  # B cells

    # Add noise
    proportions += np.abs(np.random.normal(0, 0.05, proportions.shape))
    proportions = proportions / proportions.sum(axis=1, keepdims=True)

    # Generate cell-type-specific gene signatures
    signatures = np.zeros((n_cell_types, n_genes))
    genes_per_type = n_genes // n_cell_types
    for ct in range(n_cell_types):
        start = ct * genes_per_type
        end = start + genes_per_type
        signatures[ct, start:end] = np.random.exponential(5, genes_per_type)
        # Shared genes
        shared = np.random.choice(n_genes, 20, replace=False)
        signatures[ct, shared] += np.random.exponential(2, 20)

    # Generate expression = proportions @ signatures + noise
    expression = proportions @ signatures
    expression += np.abs(np.random.normal(0, 0.5, expression.shape))
    expression = np.round(expression).astype(int)

    # Add spatially variable genes
    sv_genes_idx = list(range(n_genes - 20, n_genes))
    for g_idx in sv_genes_idx:
        # Radial pattern
        expression[:, g_idx] = np.round(
            10 * np.exp(-dist_to_center**2 / 0.1) + np.random.poisson(1, n_spots)
        ).astype(int)

    gene_names = [f'Gene_{i}' for i in range(n_genes)]
    # Name some marker genes
    marker_names = {
        0: ['EPCAM', 'KRT18', 'MUC1'],
        1: ['CD8A', 'CD8B', 'GZMA', 'PRF1'],
        2: ['CD68', 'CD163', 'CSF1R'],
        3: ['COL1A1', 'FAP', 'ACTA2'],
        4: ['CD19', 'MS4A1', 'CD79A'],
        5: ['PECAM1', 'VWF', 'CDH5']
    }
    for ct, markers in marker_names.items():
        for i, m in enumerate(markers):
            idx = ct * genes_per_type + i
            if idx < n_genes:
                gene_names[idx] = m

    # Ligand-receptor pairs (embed in gene names)
    lr_pairs = [
        ('CXCL9', 'CXCR3'),   # T cell recruitment
        ('CCL2', 'CCR2'),      # Macrophage recruitment
        ('PDCD1', 'CD274'),    # PD1-PDL1 checkpoint
        ('VEGFA', 'KDR'),      # Angiogenesis
        ('TGFB1', 'TGFBR1'),  # Fibroblast activation
    ]
    lr_gene_indices = {}
    base_idx = n_genes - 20
    for i, (lig, rec) in enumerate(lr_pairs):
        l_idx = base_idx + i * 2
        r_idx = base_idx + i * 2 + 1
        if l_idx < n_genes and r_idx < n_genes:
            gene_names[l_idx] = lig
            gene_names[r_idx] = rec
            lr_gene_indices[(lig, rec)] = (l_idx, r_idx)

            # Spatial patterns for ligand-receptor pairs
            if lig == 'CXCL9':
                expression[:, l_idx] = np.round(8 * (immune_mask.astype(float)) + np.random.poisson(1, n_spots))
                expression[:, r_idx] = np.round(6 * (immune_mask.astype(float)) + np.random.poisson(1, n_spots))
            elif lig == 'CCL2':
                expression[:, l_idx] = np.round(5 * (tumor_mask.astype(float)) + np.random.poisson(1, n_spots))
                expression[:, r_idx] = np.round(7 * (immune_mask.astype(float)) + np.random.poisson(1, n_spots))
            elif lig == 'PDCD1':
                expression[:, l_idx] = np.round(6 * (immune_mask.astype(float)) + np.random.poisson(1, n_spots))
                expression[:, r_idx] = np.round(8 * (tumor_mask.astype(float)) + np.random.poisson(1, n_spots))
            elif lig == 'VEGFA':
                expression[:, l_idx] = np.round(7 * (tumor_mask.astype(float)) + np.random.poisson(1, n_spots))
                expression[:, r_idx] = np.round(5 * (stromal_mask.astype(float)) + np.random.poisson(1, n_spots))
            elif lig == 'TGFB1':
                expression[:, l_idx] = np.round(6 * (tumor_mask.astype(float)) + np.random.poisson(1, n_spots))
                expression[:, r_idx] = np.round(8 * (stromal_mask.astype(float)) + np.random.poisson(1, n_spots))

    adata = ad.AnnData(
        X=sparse.csr_matrix(expression.astype(np.float32)),
        obs=pd.DataFrame({
            'spot_id': [f'spot_{i}' for i in range(n_spots)],
            'region': ['tumor_core' if tumor_mask[i] else
                        'immune_border' if immune_mask[i] else
                        'stroma' for i in range(n_spots)]
        }, index=[f'spot_{i}' for i in range(n_spots)]),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm['spatial'] = coords
    adata.uns['ground_truth_proportions'] = proportions
    adata.uns['cell_type_names'] = cell_type_names
    adata.uns['signatures'] = signatures
    adata.uns['lr_pairs'] = lr_pairs
    adata.uns['lr_gene_indices'] = lr_gene_indices
    adata.uns['spatial_domains'] = {
        'tumor_core': tumor_mask,
        'immune_border': immune_mask,
        'stroma': stromal_mask
    }

    return adata


# ============================================================
# 2. SPOT DECONVOLUTION
# ============================================================

def spot_deconvolution(adata):
    """Non-negative least squares deconvolution + NMF-based approach."""
    print("=" * 60)
    print("MODULE 1: Spot Deconvolution (Cell Type Composition)")
    print("=" * 60)

    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    signatures = adata.uns['signatures']
    cell_type_names = adata.uns['cell_type_names']
    gt_proportions = adata.uns['ground_truth_proportions']
    n_spots = X.shape[0]
    n_cell_types = signatures.shape[0]

    # --- Method 1: NNLS deconvolution ---
    nnls_proportions = np.zeros((n_spots, n_cell_types))
    for i in range(n_spots):
        coef, _ = nnls(signatures.T, X[i])
        total = coef.sum()
        nnls_proportions[i] = coef / total if total > 0 else np.ones(n_cell_types) / n_cell_types

    # --- Method 2: NMF-based ---
    nmf = NMF(n_components=n_cell_types, init='nndsvda', random_state=42, max_iter=500)
    W = nmf.fit_transform(X + 1e-6)
    W_norm = W / W.sum(axis=1, keepdims=True)

    # Match NMF components to cell types via correlation
    corr_matrix = np.zeros((n_cell_types, n_cell_types))
    for i in range(n_cell_types):
        for j in range(n_cell_types):
            corr_matrix[i, j], _ = pearsonr(gt_proportions[:, i], W_norm[:, j])

    # Greedy matching
    component_map = {}
    used = set()
    for _ in range(n_cell_types):
        best_i, best_j = np.unravel_index(
            np.argmax(np.where(np.isin(np.arange(n_cell_types*n_cell_types).reshape(n_cell_types, n_cell_types),
                                        [i*n_cell_types+j for i in range(n_cell_types) for j in used]),
                               -1, corr_matrix)),
            corr_matrix.shape
        )
        component_map[best_i] = best_j
        used.add(best_j)

    nmf_proportions = np.zeros_like(W_norm)
    for ct, comp in component_map.items():
        nmf_proportions[:, ct] = W_norm[:, comp]

    # Evaluate accuracy
    nnls_corrs = []
    nmf_corrs = []
    for ct in range(n_cell_types):
        c1, _ = pearsonr(gt_proportions[:, ct], nnls_proportions[:, ct])
        c2, _ = pearsonr(gt_proportions[:, ct], nmf_proportions[:, ct])
        nnls_corrs.append(c1)
        nmf_corrs.append(c2)

    results = pd.DataFrame({
        'Cell_Type': cell_type_names,
        'NNLS_Correlation': nnls_corrs,
        'NMF_Correlation': nmf_corrs,
    })
    print("\nDeconvolution accuracy (Pearson correlation with ground truth):")
    print(results.to_string(index=False))

    adata.obsm['deconv_nnls'] = nnls_proportions
    adata.obsm['deconv_nmf'] = nmf_proportions
    adata.uns['deconv_results'] = results

    # --- Visualization ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    coords = adata.obsm['spatial']

    for i, ct_name in enumerate(cell_type_names):
        ax = axes[i // 3, i % 3]
        sc_plot = ax.scatter(coords[:, 0], coords[:, 1],
                             c=nnls_proportions[:, i], cmap='Reds',
                             s=8, vmin=0, vmax=1)
        ax.set_title(f'{ct_name}\n(r={nnls_corrs[i]:.3f})', fontsize=11)
        ax.set_aspect('equal')
        plt.colorbar(sc_plot, ax=ax, fraction=0.046)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

    plt.suptitle('Spot Deconvolution: NNLS Cell Type Proportions', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig1_deconvolution.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Comparison barplot
    fig, ax = plt.subplots(figsize=(10, 6))
    x_pos = np.arange(n_cell_types)
    width = 0.35
    ax.bar(x_pos - width/2, nnls_corrs, width, label='NNLS', color='steelblue')
    ax.bar(x_pos + width/2, nmf_corrs, width, label='NMF', color='coral')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(cell_type_names, rotation=45, ha='right')
    ax.set_ylabel('Pearson Correlation with Ground Truth')
    ax.set_title('Deconvolution Method Comparison')
    ax.legend()
    ax.set_ylim(0, 1.1)
    for i, (v1, v2) in enumerate(zip(nnls_corrs, nmf_corrs)):
        ax.text(i - width/2, v1 + 0.02, f'{v1:.2f}', ha='center', fontsize=8)
        ax.text(i + width/2, v2 + 0.02, f'{v2:.2f}', ha='center', fontsize=8)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig2_deconv_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    return results


# ============================================================
# 3. SPATIALLY VARIABLE GENE DETECTION
# ============================================================

def spatially_variable_genes(adata, n_top=20):
    """SpatialDE-inspired spatially variable gene detection using
    Gaussian process regression with squared exponential kernel."""
    print("\n" + "=" * 60)
    print("MODULE 2: Spatially Variable Gene Detection")
    print("=" * 60)

    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    coords = adata.obsm['spatial']
    n_spots, n_genes = X.shape
    gene_names = adata.var_names.tolist()

    # Compute spatial distance matrix
    D = squareform(pdist(coords))

    # For each gene, compute Moran's I and a GP-inspired variance ratio
    results = []
    # Spatial weight matrix (inverse distance, capped)
    W = 1.0 / (D + 1e-6)
    np.fill_diagonal(W, 0)
    W_row_sum = W.sum(axis=1, keepdims=True)
    W_norm = W / W_row_sum

    for g in range(n_genes):
        y = X[:, g].copy()
        y_mean = y.mean()
        y_centered = y - y_mean

        if y.var() < 1e-6:
            results.append({'gene': gene_names[g], 'morans_I': 0, 'variance_ratio': 0,
                            'pvalue': 1.0, 'spatial_score': 0})
            continue

        # Moran's I
        n = len(y)
        numerator = n * np.sum(W_norm * np.outer(y_centered, y_centered))
        denominator = W_norm.sum() * np.sum(y_centered**2)
        I = numerator / denominator if denominator != 0 else 0

        # Variance ratio: spatial variance / total variance
        # Using k-nearest neighbors smoothing
        k = 15
        nn = NearestNeighbors(n_neighbors=k).fit(coords)
        _, indices = nn.kneighbors(coords)
        y_smooth = np.mean(X[indices, g], axis=1)
        var_spatial = np.var(y_smooth)
        var_total = np.var(y)
        var_ratio = var_spatial / var_total if var_total > 0 else 0

        # P-value via permutation (approximated with normal distribution)
        E_I = -1.0 / (n - 1)
        var_I = 1.0 / n  # Simplified
        z_score = (I - E_I) / np.sqrt(var_I) if var_I > 0 else 0
        pval = 2 * (1 - norm.cdf(abs(z_score)))

        spatial_score = (abs(I) + var_ratio) / 2

        results.append({
            'gene': gene_names[g],
            'morans_I': I,
            'variance_ratio': var_ratio,
            'pvalue': pval,
            'spatial_score': spatial_score,
        })

    sv_df = pd.DataFrame(results).sort_values('spatial_score', ascending=False)
    sv_df['rank'] = range(1, len(sv_df) + 1)
    sv_df['significant'] = sv_df['pvalue'] < 0.05

    top_genes = sv_df.head(n_top)
    print(f"\nTop {n_top} spatially variable genes:")
    print(top_genes[['rank', 'gene', 'morans_I', 'variance_ratio', 'pvalue', 'spatial_score']].to_string(index=False))
    print(f"\nTotal significant SVGs (p < 0.05): {sv_df['significant'].sum()} / {n_genes}")

    adata.var['morans_I'] = sv_df.set_index('gene')['morans_I']
    adata.var['spatial_score'] = sv_df.set_index('gene')['spatial_score']
    adata.var['sv_pvalue'] = sv_df.set_index('gene')['pvalue']
    adata.uns['sv_results'] = sv_df

    # --- Visualization ---
    # Top SVGs spatial expression
    top6 = sv_df.head(6)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    coords_plot = adata.obsm['spatial']

    for idx, (_, row) in enumerate(top6.iterrows()):
        ax = axes[idx // 3, idx % 3]
        gene_idx = gene_names.index(row['gene'])
        expr = X[:, gene_idx]
        sc_plot = ax.scatter(coords_plot[:, 0], coords_plot[:, 1],
                             c=expr, cmap='viridis', s=8)
        ax.set_title(f"{row['gene']}\nMoran's I={row['morans_I']:.3f}, p={row['pvalue']:.2e}", fontsize=10)
        ax.set_aspect('equal')
        plt.colorbar(sc_plot, ax=ax, fraction=0.046)

    plt.suptitle('Top Spatially Variable Genes', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig3_spatially_variable_genes.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Volcano-like plot
    fig, ax = plt.subplots(figsize=(10, 7))
    sig = sv_df[sv_df['significant']]
    nonsig = sv_df[~sv_df['significant']]
    ax.scatter(nonsig['morans_I'], -np.log10(nonsig['pvalue'] + 1e-300),
               c='gray', alpha=0.5, s=15, label='Non-significant')
    ax.scatter(sig['morans_I'], -np.log10(sig['pvalue'] + 1e-300),
               c='red', alpha=0.7, s=20, label='Significant (p<0.05)')
    for _, row in top6.iterrows():
        ax.annotate(row['gene'], (row['morans_I'], -np.log10(row['pvalue'] + 1e-300)),
                     fontsize=8, ha='left')
    ax.set_xlabel("Moran's I")
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Spatial Variability Volcano Plot')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig4_svg_volcano.png', dpi=150, bbox_inches='tight')
    plt.close()

    return sv_df


# ============================================================
# 4. CELL-CELL COMMUNICATION (LIGAND-RECEPTOR)
# ============================================================

def cell_communication(adata):
    """Ligand-receptor interaction analysis with spatial context."""
    print("\n" + "=" * 60)
    print("MODULE 3: Cell-Cell Communication (Ligand-Receptor)")
    print("=" * 60)

    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    coords = adata.obsm['spatial']
    gene_names = adata.var_names.tolist()
    lr_pairs = adata.uns['lr_pairs']
    regions = adata.obs['region'].values

    # Build spatial neighbor graph
    nn = NearestNeighbors(n_neighbors=10).fit(coords)
    distances, indices = nn.kneighbors(coords)

    results = []
    for lig, rec in lr_pairs:
        if lig not in gene_names or rec not in gene_names:
            continue
        l_idx = gene_names.index(lig)
        r_idx = gene_names.index(rec)

        lig_expr = X[:, l_idx]
        rec_expr = X[:, r_idx]

        # Interaction score: product of ligand and receptor expression
        # weighted by spatial proximity
        interaction_scores = np.zeros(len(coords))
        for i in range(len(coords)):
            neighbor_rec = rec_expr[indices[i, 1:]]  # exclude self
            neighbor_dist = distances[i, 1:]
            weights = 1.0 / (neighbor_dist + 1e-6)
            weights /= weights.sum()
            interaction_scores[i] = lig_expr[i] * np.sum(weights * neighbor_rec)

        # Region-wise interaction strength
        region_scores = {}
        for region in ['tumor_core', 'immune_border', 'stroma']:
            mask = regions == region
            region_scores[region] = interaction_scores[mask].mean()

        # Permutation test
        n_perms = 999
        perm_scores = np.zeros(n_perms)
        for p in range(n_perms):
            perm_rec = np.random.permutation(rec_expr)
            for i in range(len(coords)):
                neighbor_rec_p = perm_rec[indices[i, 1:]]
                neighbor_dist = distances[i, 1:]
                weights = 1.0 / (neighbor_dist + 1e-6)
                weights /= weights.sum()
                perm_scores[p] += lig_expr[i] * np.sum(weights * neighbor_rec_p)
            perm_scores[p] /= len(coords)

        observed_mean = interaction_scores.mean()
        pval = (np.sum(perm_scores >= observed_mean) + 1) / (n_perms + 1)

        results.append({
            'ligand': lig,
            'receptor': rec,
            'mean_score': observed_mean,
            'tumor_core': region_scores.get('tumor_core', 0),
            'immune_border': region_scores.get('immune_border', 0),
            'stroma': region_scores.get('stroma', 0),
            'pvalue': pval,
        })

        adata.obs[f'LR_{lig}_{rec}'] = interaction_scores

    lr_df = pd.DataFrame(results)
    lr_df['significant'] = lr_df['pvalue'] < 0.05
    print("\nLigand-Receptor Interaction Results:")
    print(lr_df.to_string(index=False))

    adata.uns['lr_results'] = lr_df

    # --- Visualization ---
    # Spatial LR interaction maps
    n_pairs = len(lr_pairs)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    coords_plot = adata.obsm['spatial']

    for idx, (lig, rec) in enumerate(lr_pairs[:6]):
        if idx >= 6:
            break
        ax = axes[idx // 3, idx % 3]
        col_name = f'LR_{lig}_{rec}'
        if col_name in adata.obs:
            scores = adata.obs[col_name].values
            sc_plot = ax.scatter(coords_plot[:, 0], coords_plot[:, 1],
                                 c=scores, cmap='YlOrRd', s=8)
            row = lr_df[(lr_df['ligand'] == lig) & (lr_df['receptor'] == rec)].iloc[0]
            ax.set_title(f'{lig}-{rec}\np={row["pvalue"]:.3f}', fontsize=11)
            ax.set_aspect('equal')
            plt.colorbar(sc_plot, ax=ax, fraction=0.046)

    if n_pairs < 6:
        for idx in range(n_pairs, 6):
            axes[idx // 3, idx % 3].set_visible(False)

    plt.suptitle('Ligand-Receptor Interaction Scores (Spatial)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig5_ligand_receptor.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Heatmap of region-wise interaction
    region_data = lr_df[['ligand', 'receptor', 'tumor_core', 'immune_border', 'stroma']].copy()
    region_data['pair'] = region_data['ligand'] + '-' + region_data['receptor']
    heat_data = region_data.set_index('pair')[['tumor_core', 'immune_border', 'stroma']]

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(heat_data, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax, linewidths=0.5)
    ax.set_title('Ligand-Receptor Interaction Strength by Region')
    ax.set_ylabel('L-R Pair')
    ax.set_xlabel('Tissue Region')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig6_lr_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()

    return lr_df


# ============================================================
# 5. TISSUE NICHE IDENTIFICATION
# ============================================================

def niche_identification(adata, n_niches=4):
    """Identify tissue microenvironment niches using neighborhood
    cell type composition clustering."""
    print("\n" + "=" * 60)
    print("MODULE 4: Tissue Microenvironment Niche Identification")
    print("=" * 60)

    coords = adata.obsm['spatial']
    cell_type_names = adata.uns['cell_type_names']
    proportions = adata.obsm['deconv_nnls']

    # Build spatial neighborhood graph
    k = 15
    nn = NearestNeighbors(n_neighbors=k).fit(coords)
    distances, indices = nn.kneighbors(coords)

    # Compute neighborhood cell type composition
    nhood_composition = np.zeros_like(proportions)
    for i in range(len(coords)):
        neighbor_props = proportions[indices[i]]
        # Distance-weighted average
        w = 1.0 / (distances[i] + 1e-6)
        w /= w.sum()
        nhood_composition[i] = np.average(neighbor_props, weights=w, axis=0)

    # Cluster neighborhood compositions to identify niches
    kmeans = KMeans(n_clusters=n_niches, random_state=42, n_init=20)
    niche_labels = kmeans.fit_predict(nhood_composition)

    # Evaluate against ground truth regions
    gt_labels = adata.obs['region'].map({
        'tumor_core': 0, 'immune_border': 1, 'stroma': 2
    }).values

    # Also try hierarchical clustering
    hc = AgglomerativeClustering(n_clusters=n_niches)
    hc_labels = hc.fit_predict(nhood_composition)

    sil_km = silhouette_score(nhood_composition, niche_labels)
    sil_hc = silhouette_score(nhood_composition, hc_labels)

    print(f"\nSilhouette Score - KMeans: {sil_km:.4f}")
    print(f"Silhouette Score - Hierarchical: {sil_hc:.4f}")

    # Niche characterization
    niche_profiles = pd.DataFrame(
        kmeans.cluster_centers_,
        columns=cell_type_names,
        index=[f'Niche_{i}' for i in range(n_niches)]
    )
    print("\nNiche Cell Type Composition Profiles:")
    print(niche_profiles.round(3).to_string())

    adata.obs['niche'] = niche_labels.astype(str)
    adata.obsm['nhood_composition'] = nhood_composition
    adata.uns['niche_profiles'] = niche_profiles
    adata.uns['niche_silhouette'] = {'kmeans': sil_km, 'hierarchical': sil_hc}

    # --- Visualization ---
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(2, 3, figure=fig)

    # Niche spatial map
    ax1 = fig.add_subplot(gs[0, 0])
    cmap = plt.cm.Set1
    for n in range(n_niches):
        mask = niche_labels == n
        ax1.scatter(coords[mask, 0], coords[mask, 1], c=[cmap(n)], s=8,
                     label=f'Niche {n}', alpha=0.7)
    ax1.set_title(f'Identified Niches (K={n_niches})\nSilhouette={sil_km:.3f}')
    ax1.set_aspect('equal')
    ax1.legend(markerscale=3, fontsize=8)

    # Ground truth regions
    ax2 = fig.add_subplot(gs[0, 1])
    region_colors = {'tumor_core': 'red', 'immune_border': 'blue', 'stroma': 'green'}
    for region, color in region_colors.items():
        mask = adata.obs['region'] == region
        ax2.scatter(coords[mask, 0], coords[mask, 1], c=color, s=8,
                     label=region, alpha=0.7)
    ax2.set_title('Ground Truth Regions')
    ax2.set_aspect('equal')
    ax2.legend(markerscale=3, fontsize=8)

    # Niche composition heatmap
    ax3 = fig.add_subplot(gs[0, 2])
    sns.heatmap(niche_profiles, annot=True, fmt='.3f', cmap='YlOrRd',
                ax=ax3, linewidths=0.5, cbar_kws={'shrink': 0.8})
    ax3.set_title('Niche Cell Type Profiles')

    # PCA of neighborhood composition
    pca = PCA(n_components=2)
    nhood_pca = pca.fit_transform(nhood_composition)
    ax4 = fig.add_subplot(gs[1, 0])
    for n in range(n_niches):
        mask = niche_labels == n
        ax4.scatter(nhood_pca[mask, 0], nhood_pca[mask, 1], c=[cmap(n)], s=8,
                     label=f'Niche {n}', alpha=0.5)
    ax4.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax4.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax4.set_title('PCA of Neighborhood Composition')
    ax4.legend(markerscale=3, fontsize=8)

    # Niche size distribution
    ax5 = fig.add_subplot(gs[1, 1])
    niche_sizes = pd.Series(niche_labels).value_counts().sort_index()
    ax5.bar(niche_sizes.index.tolist(), niche_sizes.values,
            color=[cmap(i) for i in niche_sizes.index])
    ax5.set_xlabel('Niche')
    ax5.set_ylabel('Number of Spots')
    ax5.set_title('Niche Size Distribution')

    # Region-niche contingency
    ax6 = fig.add_subplot(gs[1, 2])
    contingency = pd.crosstab(adata.obs['region'], niche_labels)
    contingency_norm = contingency.div(contingency.sum(axis=0), axis=1)
    contingency_norm.plot(kind='bar', stacked=True, ax=ax6, colormap='Set1')
    ax6.set_title('Region Composition per Niche')
    ax6.set_ylabel('Proportion')
    ax6.legend(title='Niche', fontsize=8)
    ax6.tick_params(axis='x', rotation=45)

    plt.suptitle('Tissue Microenvironment Niche Analysis', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig7_niche_identification.png', dpi=150, bbox_inches='tight')
    plt.close()

    return niche_profiles


# ============================================================
# 6. 3D SPATIAL RECONSTRUCTION
# ============================================================

def spatial_3d_reconstruction(adata, n_sections=5):
    """Simulate serial section integration and 3D reconstruction
    using optimal transport-inspired alignment (PASTE-like)."""
    print("\n" + "=" * 60)
    print("MODULE 5: 3D Spatial Reconstruction (Serial Sections)")
    print("=" * 60)

    coords_2d = adata.obsm['spatial']
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    n_spots = len(coords_2d)
    spots_per_section = n_spots // n_sections

    # Simulate serial sections with slight spatial shifts and rotation
    sections = []
    section_labels = []
    all_coords_3d = []

    for s in range(n_sections):
        start = s * spots_per_section
        end = start + spots_per_section
        if s == n_sections - 1:
            end = n_spots

        section_coords = coords_2d[start:end].copy()

        # Add section-specific transformation (simulate cutting artifacts)
        angle = np.random.uniform(-0.05, 0.05)  # Small rotation
        shift = np.random.uniform(-0.02, 0.02, 2)  # Small translation
        rotation = np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])
        section_coords = section_coords @ rotation.T + shift

        # 3D coordinates
        z = s * 0.1  # z-spacing between sections
        coords_3d = np.column_stack([section_coords, np.full(len(section_coords), z)])

        sections.append({
            'coords_2d': section_coords,
            'coords_3d': coords_3d,
            'expression': X[start:end],
            'indices': list(range(start, end))
        })
        section_labels.extend([s] * len(section_coords))
        all_coords_3d.append(coords_3d)

    # Pairwise alignment between consecutive sections (simplified PASTE)
    alignment_scores = []
    aligned_coords_3d = [sections[0]['coords_3d'].copy()]

    for s in range(1, n_sections):
        prev = sections[s - 1]
        curr = sections[s]

        # Expression similarity (cosine)
        expr_sim = 1 - cdist(
            normalize(prev['expression'], norm='l2'),
            normalize(curr['expression'], norm='l2'),
            metric='cosine'
        )

        # Spatial distance between sections
        spatial_dist = cdist(prev['coords_2d'], curr['coords_2d'])

        # Combined cost: spatial + expression
        alpha = 0.5
        cost = alpha * spatial_dist + (1 - alpha) * (1 - expr_sim)

        # Find best matching (Hungarian-like greedy)
        n_prev = len(prev['coords_2d'])
        n_curr = len(curr['coords_2d'])
        n_matches = min(n_prev, n_curr, 50)  # Sample for speed

        # Subsample for alignment
        prev_sample = np.random.choice(n_prev, min(n_matches, n_prev), replace=False)
        curr_sample = np.random.choice(n_curr, min(n_matches, n_curr), replace=False)

        cost_sub = cost[np.ix_(prev_sample, curr_sample)]

        # Greedy matching on subsample
        matches = []
        used_curr = set()
        for p_idx in range(len(prev_sample)):
            costs = cost_sub[p_idx].copy()
            for uc in used_curr:
                costs[uc] = np.inf
            best = np.argmin(costs)
            if costs[best] < np.inf:
                matches.append((prev_sample[p_idx], curr_sample[best]))
                used_curr.add(best)

        # Compute alignment transformation from matches
        if len(matches) > 3:
            src_pts = prev['coords_2d'][[m[0] for m in matches]]
            tgt_pts = curr['coords_2d'][[m[1] for m in matches]]

            # Procrustes-like alignment
            src_mean = src_pts.mean(axis=0)
            tgt_mean = tgt_pts.mean(axis=0)
            shift = src_mean - tgt_mean

            aligned = curr['coords_3d'].copy()
            aligned[:, :2] += shift

            alignment_score = np.mean([
                1 - cost[m[0], m[1]] for m in matches
            ])
        else:
            aligned = curr['coords_3d'].copy()
            alignment_score = 0

        aligned_coords_3d.append(aligned)
        alignment_scores.append(alignment_score)

    print(f"\nSection alignment scores (consecutive pairs):")
    for i, score in enumerate(alignment_scores):
        print(f"  Section {i} -> {i+1}: {score:.4f}")
    print(f"  Mean alignment: {np.mean(alignment_scores):.4f}")

    all_aligned = np.vstack(aligned_coords_3d)
    adata.obsm['spatial_3d'] = np.zeros((n_spots, 3))
    idx = 0
    for s, sec in enumerate(sections):
        n = len(sec['indices'])
        adata.obsm['spatial_3d'][sec['indices']] = aligned_coords_3d[s][:n]
        idx += n

    adata.obs['section'] = section_labels[:n_spots]
    adata.uns['alignment_scores'] = alignment_scores

    # --- Visualization ---
    fig = plt.figure(figsize=(18, 8))

    # 3D scatter colored by section
    ax1 = fig.add_subplot(131, projection='3d')
    for s in range(n_sections):
        mask = np.array(section_labels[:n_spots]) == s
        c3d = adata.obsm['spatial_3d'][mask]
        ax1.scatter(c3d[:, 0], c3d[:, 1], c3d[:, 2], s=3, alpha=0.5,
                     label=f'Section {s}')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title('3D Reconstruction\n(colored by section)')
    ax1.legend(fontsize=7, markerscale=3)

    # 3D scatter colored by region
    ax2 = fig.add_subplot(132, projection='3d')
    region_map = {'tumor_core': 'red', 'immune_border': 'blue', 'stroma': 'green'}
    for region, color in region_map.items():
        mask = adata.obs['region'].values == region
        c3d = adata.obsm['spatial_3d'][mask]
        ax2.scatter(c3d[:, 0], c3d[:, 1], c3d[:, 2], s=3, alpha=0.5,
                     c=color, label=region)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_title('3D Reconstruction\n(colored by region)')
    ax2.legend(fontsize=7, markerscale=3)

    # Alignment quality
    ax3 = fig.add_subplot(133)
    ax3.bar(range(len(alignment_scores)), alignment_scores, color='steelblue')
    ax3.set_xlabel('Section Pair')
    ax3.set_ylabel('Alignment Score')
    ax3.set_title('Pairwise Section Alignment Quality')
    ax3.set_xticks(range(len(alignment_scores)))
    ax3.set_xticklabels([f'{i}-{i+1}' for i in range(len(alignment_scores))])
    for i, s in enumerate(alignment_scores):
        ax3.text(i, s + 0.01, f'{s:.3f}', ha='center', fontsize=9)

    plt.suptitle('3D Spatial Reconstruction from Serial Sections', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig8_3d_reconstruction.png', dpi=150, bbox_inches='tight')
    plt.close()

    return alignment_scores


# ============================================================
# 7. TUMOR IMMUNE MICROENVIRONMENT CASE STUDY
# ============================================================

def tumor_immune_case_study(adata):
    """Comprehensive tumor immune microenvironment (TIME) analysis."""
    print("\n" + "=" * 60)
    print("MODULE 6: Tumor Immune Microenvironment Case Study")
    print("=" * 60)

    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    coords = adata.obsm['spatial']
    proportions = adata.obsm['deconv_nnls']
    cell_type_names = adata.uns['cell_type_names']
    regions = adata.obs['region'].values

    # --- Immune infiltration score ---
    immune_types = [1, 2, 4]  # CD8 T, Macrophage, B cell
    immune_score = proportions[:, immune_types].sum(axis=1)
    adata.obs['immune_score'] = immune_score

    # --- Tumor-immune interface analysis ---
    tumor_mask = regions == 'tumor_core'
    immune_border_mask = regions == 'immune_border'

    # Distance from each spot to tumor boundary
    nn_tumor = NearestNeighbors(n_neighbors=1).fit(coords[tumor_mask])
    dist_to_tumor, _ = nn_tumor.kneighbors(coords)
    adata.obs['dist_to_tumor'] = dist_to_tumor.ravel()

    # Immune score gradient from tumor
    bins = np.linspace(0, dist_to_tumor.max(), 20)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    immune_gradient = []
    tumor_gradient = []
    for i in range(len(bins) - 1):
        mask = (dist_to_tumor.ravel() >= bins[i]) & (dist_to_tumor.ravel() < bins[i+1])
        if mask.sum() > 0:
            immune_gradient.append(immune_score[mask].mean())
            tumor_gradient.append(proportions[mask, 0].mean())
        else:
            immune_gradient.append(np.nan)
            tumor_gradient.append(np.nan)

    # --- Immune checkpoint expression ---
    gene_names = adata.var_names.tolist()
    checkpoint_genes = ['PDCD1', 'CD274']  # PD1, PDL1
    checkpoint_expr = {}
    for g in checkpoint_genes:
        if g in gene_names:
            idx = gene_names.index(g)
            checkpoint_expr[g] = X[:, idx]

    # --- Immune hot vs cold regions ---
    median_immune = np.median(immune_score)
    hot_spots = immune_score > median_immune
    cold_spots = immune_score <= median_immune
    adata.obs['immune_status'] = ['hot' if h else 'cold' for h in hot_spots]

    # Statistics
    hot_tumor_prop = proportions[hot_spots, 0].mean()
    cold_tumor_prop = proportions[cold_spots, 0].mean()
    t_stat, t_pval = ttest_ind(proportions[hot_spots, 0], proportions[cold_spots, 0])

    print(f"\nImmune infiltration summary:")
    print(f"  Mean immune score (all): {immune_score.mean():.4f}")
    print(f"  Mean immune score (tumor core): {immune_score[tumor_mask].mean():.4f}")
    print(f"  Mean immune score (immune border): {immune_score[immune_border_mask].mean():.4f}")
    print(f"  Mean immune score (stroma): {immune_score[regions == 'stroma'].mean():.4f}")
    print(f"\nImmune hot vs cold:")
    print(f"  Hot spots: {hot_spots.sum()} ({hot_spots.sum()/len(hot_spots)*100:.1f}%)")
    print(f"  Cold spots: {cold_spots.sum()} ({cold_spots.sum()/len(cold_spots)*100:.1f}%)")
    print(f"  Tumor proportion in hot: {hot_tumor_prop:.4f}")
    print(f"  Tumor proportion in cold: {cold_tumor_prop:.4f}")
    print(f"  T-test p-value: {t_pval:.4e}")

    # Checkpoint expression by region
    if checkpoint_expr:
        print("\nCheckpoint gene expression by region:")
        for g, expr in checkpoint_expr.items():
            for region in ['tumor_core', 'immune_border', 'stroma']:
                mask = regions == region
                print(f"  {g} in {region}: {expr[mask].mean():.3f} ± {expr[mask].std():.3f}")

    adata.uns['time_results'] = {
        'immune_gradient': immune_gradient,
        'tumor_gradient': tumor_gradient,
        'bin_centers': bin_centers.tolist(),
        'hot_spots_n': int(hot_spots.sum()),
        'cold_spots_n': int(cold_spots.sum()),
        'hot_tumor_prop': float(hot_tumor_prop),
        'cold_tumor_prop': float(cold_tumor_prop),
        't_pval': float(t_pval),
    }

    # --- Visualization ---
    fig = plt.figure(figsize=(24, 16))
    gs = gridspec.GridSpec(3, 4, figure=fig)

    # 1. Immune score spatial map
    ax1 = fig.add_subplot(gs[0, 0])
    sc1 = ax1.scatter(coords[:, 0], coords[:, 1], c=immune_score,
                       cmap='RdYlBu_r', s=8)
    ax1.set_title('Immune Infiltration Score')
    ax1.set_aspect('equal')
    plt.colorbar(sc1, ax=ax1, fraction=0.046)

    # 2. Hot vs Cold
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(coords[hot_spots, 0], coords[hot_spots, 1],
                c='red', s=8, alpha=0.5, label='Hot')
    ax2.scatter(coords[cold_spots, 0], coords[cold_spots, 1],
                c='blue', s=8, alpha=0.5, label='Cold')
    ax2.set_title(f'Immune Hot/Cold\n(p={t_pval:.2e})')
    ax2.set_aspect('equal')
    ax2.legend(markerscale=3)

    # 3. Distance to tumor
    ax3 = fig.add_subplot(gs[0, 2])
    sc3 = ax3.scatter(coords[:, 0], coords[:, 1],
                       c=dist_to_tumor.ravel(), cmap='viridis', s=8)
    ax3.set_title('Distance to Tumor Core')
    ax3.set_aspect('equal')
    plt.colorbar(sc3, ax=ax3, fraction=0.046)

    # 4. Immune gradient
    ax4 = fig.add_subplot(gs[0, 3])
    ax4.plot(bin_centers, immune_gradient, 'b-o', markersize=4, label='Immune Score')
    ax4.set_xlabel('Distance from Tumor')
    ax4.set_ylabel('Immune Score', color='blue')
    ax4_twin = ax4.twinx()
    ax4_twin.plot(bin_centers, tumor_gradient, 'r-s', markersize=4, label='Tumor Prop.')
    ax4_twin.set_ylabel('Tumor Proportion', color='red')
    ax4.set_title('Immune & Tumor Gradients')
    lines1, labels1 = ax4.get_legend_handles_labels()
    lines2, labels2 = ax4_twin.get_legend_handles_labels()
    ax4.legend(lines1 + lines2, labels1 + labels2, fontsize=8)

    # 5-6. Checkpoint expression
    for g_idx, (g_name, g_expr) in enumerate(checkpoint_expr.items()):
        ax = fig.add_subplot(gs[1, g_idx])
        sc_g = ax.scatter(coords[:, 0], coords[:, 1], c=g_expr,
                           cmap='Purples', s=8)
        ax.set_title(f'{g_name} Expression')
        ax.set_aspect('equal')
        plt.colorbar(sc_g, ax=ax, fraction=0.046)

    # 7. Cell type composition by region (boxplot)
    ax7 = fig.add_subplot(gs[1, 2:])
    box_data = []
    for ct_idx, ct_name in enumerate(cell_type_names):
        for region in ['tumor_core', 'immune_border', 'stroma']:
            mask = regions == region
            for val in proportions[mask, ct_idx]:
                box_data.append({'Cell Type': ct_name, 'Region': region, 'Proportion': val})
    box_df = pd.DataFrame(box_data)
    sns.boxplot(data=box_df, x='Cell Type', y='Proportion', hue='Region', ax=ax7)
    ax7.set_title('Cell Type Composition by Region')
    ax7.tick_params(axis='x', rotation=45)
    ax7.legend(fontsize=8)

    # 8. Correlation of immune types
    ax8 = fig.add_subplot(gs[2, 0:2])
    immune_corr = np.corrcoef(proportions[:, immune_types].T)
    sns.heatmap(immune_corr, annot=True, fmt='.3f', cmap='coolwarm',
                xticklabels=[cell_type_names[i] for i in immune_types],
                yticklabels=[cell_type_names[i] for i in immune_types],
                ax=ax8, vmin=-1, vmax=1)
    ax8.set_title('Immune Cell Type Correlations')

    # 9. Immune score distribution by region
    ax9 = fig.add_subplot(gs[2, 2:])
    for region, color in [('tumor_core', 'red'), ('immune_border', 'blue'), ('stroma', 'green')]:
        mask = regions == region
        ax9.hist(immune_score[mask], bins=30, alpha=0.5, color=color,
                  label=region, density=True)
    ax9.set_xlabel('Immune Score')
    ax9.set_ylabel('Density')
    ax9.set_title('Immune Score Distribution by Region')
    ax9.legend()

    plt.suptitle('Tumor Immune Microenvironment (TIME) Analysis', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig9_tumor_immune.png', dpi=150, bbox_inches='tight')
    plt.close()

    return adata.uns['time_results']


# ============================================================
# 8. PIPELINE SUMMARY FIGURE
# ============================================================

def create_summary_figure(adata):
    """Create a comprehensive summary overview figure."""
    print("\n" + "=" * 60)
    print("Creating Summary Figure")
    print("=" * 60)

    coords = adata.obsm['spatial']
    proportions = adata.obsm['deconv_nnls']
    cell_type_names = adata.uns['cell_type_names']

    fig, axes = plt.subplots(2, 4, figsize=(24, 12))

    # 1. Raw spatial data (total counts)
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    total_counts = X.sum(axis=1)
    sc1 = axes[0, 0].scatter(coords[:, 0], coords[:, 1], c=total_counts,
                              cmap='viridis', s=6)
    axes[0, 0].set_title('A. Total UMI Counts')
    axes[0, 0].set_aspect('equal')
    plt.colorbar(sc1, ax=axes[0, 0], fraction=0.046)

    # 2. Ground truth regions
    region_colors = {'tumor_core': 'red', 'immune_border': 'blue', 'stroma': 'green'}
    for region, color in region_colors.items():
        mask = adata.obs['region'] == region
        axes[0, 1].scatter(coords[mask, 0], coords[mask, 1], c=color, s=6,
                            label=region, alpha=0.7)
    axes[0, 1].set_title('B. Tissue Regions')
    axes[0, 1].set_aspect('equal')
    axes[0, 1].legend(markerscale=3, fontsize=8)

    # 3. Dominant cell type (deconvolution)
    dominant_ct = np.argmax(proportions, axis=1)
    cmap = plt.cm.Set2
    for ct in range(len(cell_type_names)):
        mask = dominant_ct == ct
        if mask.sum() > 0:
            axes[0, 2].scatter(coords[mask, 0], coords[mask, 1],
                               c=[cmap(ct)], s=6, label=cell_type_names[ct], alpha=0.7)
    axes[0, 2].set_title('C. Dominant Cell Type')
    axes[0, 2].set_aspect('equal')
    axes[0, 2].legend(markerscale=3, fontsize=6)

    # 4. Moran's I distribution
    morans = adata.var['morans_I'].dropna()
    axes[0, 3].hist(morans, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    axes[0, 3].axvline(0, color='red', linestyle='--')
    axes[0, 3].set_xlabel("Moran's I")
    axes[0, 3].set_ylabel('Count')
    axes[0, 3].set_title("D. Moran's I Distribution")

    # 5. LR interaction (top pair)
    lr_results = adata.uns['lr_results']
    top_pair = lr_results.iloc[0]
    col = f"LR_{top_pair['ligand']}_{top_pair['receptor']}"
    if col in adata.obs:
        sc5 = axes[1, 0].scatter(coords[:, 0], coords[:, 1],
                                  c=adata.obs[col].values, cmap='YlOrRd', s=6)
        axes[1, 0].set_title(f"E. Top LR: {top_pair['ligand']}-{top_pair['receptor']}")
        axes[1, 0].set_aspect('equal')
        plt.colorbar(sc5, ax=axes[1, 0], fraction=0.046)

    # 6. Niches
    niche_labels = adata.obs['niche'].astype(int).values
    n_niches = len(np.unique(niche_labels))
    cmap_n = plt.cm.Set1
    for n in range(n_niches):
        mask = niche_labels == n
        axes[1, 1].scatter(coords[mask, 0], coords[mask, 1], c=[cmap_n(n)], s=6,
                            label=f'Niche {n}', alpha=0.7)
    axes[1, 1].set_title('F. Tissue Niches')
    axes[1, 1].set_aspect('equal')
    axes[1, 1].legend(markerscale=3, fontsize=8)

    # 7. Immune score
    sc7 = axes[1, 2].scatter(coords[:, 0], coords[:, 1],
                              c=adata.obs['immune_score'].values,
                              cmap='RdYlBu_r', s=6)
    axes[1, 2].set_title('G. Immune Infiltration')
    axes[1, 2].set_aspect('equal')
    plt.colorbar(sc7, ax=axes[1, 2], fraction=0.046)

    # 8. Pipeline overview metrics
    ax8 = axes[1, 3]
    ax8.axis('off')
    metrics_text = (
        "Pipeline Summary\n"
        "─────────────────────\n"
        f"Spots: {adata.n_obs}\n"
        f"Genes: {adata.n_vars}\n"
        f"Cell Types: {len(cell_type_names)}\n"
        f"SVGs (p<0.05): {(adata.var['sv_pvalue'] < 0.05).sum()}\n"
        f"LR Pairs: {len(lr_results)}\n"
        f"Significant LR: {lr_results['significant'].sum()}\n"
        f"Niches: {n_niches}\n"
        f"Sections: {len(adata.uns['alignment_scores'])+1}\n"
        f"Immune Hot: {adata.uns['time_results']['hot_spots_n']}\n"
        f"Immune Cold: {adata.uns['time_results']['cold_spots_n']}\n"
    )
    ax8.text(0.1, 0.5, metrics_text, transform=ax8.transAxes,
             fontsize=12, verticalalignment='center', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax8.set_title('H. Summary Metrics')

    plt.suptitle('Spatial Transcriptomics Analysis Pipeline Overview', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig10_summary.png', dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# MAIN PIPELINE
# ============================================================

if __name__ == '__main__':
    print("╔══════════════════════════════════════════════════════════╗")
    print("║  Spatial Transcriptomics Advanced Analysis Pipeline     ║")
    print("║  Visium/MERFISH Integrated Framework                    ║")
    print("╚══════════════════════════════════════════════════════════╝")

    # Generate synthetic data
    print("\nGenerating synthetic spatial transcriptomics data...")
    adata = generate_synthetic_spatial_data(n_spots=2000, n_genes=500, n_cell_types=6)
    print(f"Data shape: {adata.shape}")
    print(f"Spatial coordinates: {adata.obsm['spatial'].shape}")
    print(f"Regions: {adata.obs['region'].value_counts().to_dict()}")

    # Run pipeline modules
    deconv_results = spot_deconvolution(adata)
    svg_results = spatially_variable_genes(adata)
    lr_results = cell_communication(adata)
    niche_results = niche_identification(adata)
    alignment_results = spatial_3d_reconstruction(adata)
    time_results = tumor_immune_case_study(adata)

    # Summary figure
    create_summary_figure(adata)

    # Save processed data
    print("\n" + "=" * 60)
    print("Pipeline Complete!")
    print("=" * 60)
    print(f"\nGenerated figures in '{FIGDIR}/':")
    print("  fig1_deconvolution.png")
    print("  fig2_deconv_comparison.png")
    print("  fig3_spatially_variable_genes.png")
    print("  fig4_svg_volcano.png")
    print("  fig5_ligand_receptor.png")
    print("  fig6_lr_heatmap.png")
    print("  fig7_niche_identification.png")
    print("  fig8_3d_reconstruction.png")
    print("  fig9_tumor_immune.png")
    print("  fig10_summary.png")
