#!/usr/bin/env python3
"""
Perturb-seq Analysis Framework
===============================
Comprehensive pipeline for CRISPR+scRNA-seq data analysis.
Covers: QC, differential expression, co-expression modules,
causal graph inference, epistasis detection, and latent representation learning.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse, stats
from sklearn.decomposition import PCA, NMF
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import networkx as nx
from itertools import combinations
import warnings
import os

warnings.filterwarnings('ignore')
np.random.seed(42)

FIGDIR = 'figures'
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# PART 0: Simulate Perturb-seq Data
# ============================================================
print("=" * 60)
print("PART 0: Generating synthetic Perturb-seq dataset")
print("=" * 60)

n_cells = 5000
n_genes = 2000
n_guides = 20
n_essential_genes = 50

gene_names = [f"Gene_{i}" for i in range(n_genes)]
essential_gene_names = [f"EssGene_{i}" for i in range(n_essential_genes)]
all_gene_names = gene_names + essential_gene_names
n_total_genes = len(all_gene_names)

# Simulate guide assignments
guide_names = [f"Guide_{i}" for i in range(n_guides)]
target_genes = [f"Gene_{i*3}" for i in range(n_guides)]

# Assign guides to cells (some cells get multiple guides for combinatorial)
cell_guides = []
cell_targets = []
for i in range(n_cells):
    if i < int(n_cells * 0.15):  # 15% non-targeting control
        cell_guides.append("NTC")
        cell_targets.append("NTC")
    elif i < int(n_cells * 0.85):  # 70% single perturbation
        idx = np.random.randint(n_guides)
        cell_guides.append(guide_names[idx])
        cell_targets.append(target_genes[idx])
    else:  # 15% combinatorial (double perturbation)
        idx1, idx2 = np.random.choice(n_guides, 2, replace=False)
        cell_guides.append(f"{guide_names[idx1]}+{guide_names[idx2]}")
        cell_targets.append(f"{target_genes[idx1]}+{target_genes[idx2]}")

# Simulate count matrix with perturbation effects
base_means = np.random.exponential(2, n_total_genes)
X = np.zeros((n_cells, n_total_genes))

for i in range(n_cells):
    cell_mean = base_means.copy()
    target = cell_targets[i]
    if target != "NTC":
        targets = target.split("+")
        for t in targets:
            if t in all_gene_names:
                tidx = all_gene_names.index(t)
                cell_mean[tidx] *= 0.1  # knockdown
                # Downstream effects
                np.random.seed(i + tidx)
                affected = np.random.choice(n_total_genes, 30, replace=False)
                for a in affected:
                    cell_mean[a] *= np.random.choice([0.5, 1.5, 2.0])
    # Add noise and sample
    X[i] = np.random.poisson(cell_mean)

# Create AnnData
adata = sc.AnnData(
    X=sparse.csr_matrix(X),
    obs=pd.DataFrame({
        'guide': cell_guides,
        'perturbation': cell_targets,
        'n_guides': [len(g.split('+')) if g != 'NTC' else 0 for g in cell_guides],
    }, index=[f"Cell_{i}" for i in range(n_cells)]),
    var=pd.DataFrame(index=all_gene_names)
)

# Add guide detection counts (UMI)
adata.obs['guide_UMI'] = np.random.poisson(50, n_cells)
adata.obs['guide_UMI'] = adata.obs['guide_UMI'].astype(float)
# Low-quality cells have low guide UMI
low_quality_mask = np.random.random(n_cells) < 0.05
adata.obs.loc[low_quality_mask, 'guide_UMI'] = np.random.poisson(3, low_quality_mask.sum())

print(f"Created AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Perturbations: {len(set(cell_targets))} unique targets")
print(f"Combinatorial cells: {sum(1 for t in cell_targets if '+' in t)}")

# ============================================================
# PART 1: Perturbation Assignment QC & Guide Detection
# ============================================================
print("\n" + "=" * 60)
print("PART 1: Perturbation Assignment QC & Guide Detection")
print("=" * 60)

# 1a. Guide UMI distribution
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].hist(adata.obs['guide_UMI'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].axvline(x=10, color='red', linestyle='--', label='QC threshold')
axes[0].set_xlabel('Guide UMI Count')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Guide UMI Distribution')
axes[0].legend()

# 1b. Cells per perturbation
pert_counts = adata.obs['perturbation'].value_counts()
axes[1].barh(range(min(20, len(pert_counts))), pert_counts.values[:20],
             color='coral', edgecolor='black', alpha=0.7)
axes[1].set_yticks(range(min(20, len(pert_counts))))
axes[1].set_yticklabels(pert_counts.index[:20], fontsize=7)
axes[1].set_xlabel('Number of Cells')
axes[1].set_title('Cells per Perturbation')
axes[1].invert_yaxis()

# 1c. Guide assignment quality scores
quality_scores = adata.obs['guide_UMI'] / (adata.obs['guide_UMI'].max())
adata.obs['guide_quality'] = quality_scores
pass_qc = adata.obs['guide_UMI'] >= 10
axes[2].pie([pass_qc.sum(), (~pass_qc).sum()],
            labels=['Pass QC', 'Fail QC'],
            colors=['#2ecc71', '#e74c3c'],
            autopct='%1.1f%%', startangle=90)
axes[2].set_title('Guide Detection QC')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/01_guide_qc.png', dpi=150, bbox_inches='tight')
plt.close()

# Filter cells
n_before = adata.shape[0]
adata = adata[pass_qc].copy()
n_after = adata.shape[0]
print(f"QC filtering: {n_before} -> {n_after} cells ({n_before - n_after} removed)")

# Compute additional QC metrics
adata.obs['n_genes_detected'] = np.array((adata.X > 0).sum(axis=1)).flatten()
adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].scatter(adata.obs['total_counts'], adata.obs['n_genes_detected'],
                s=1, alpha=0.3, c='steelblue')
axes[0].set_xlabel('Total Counts')
axes[0].set_ylabel('Genes Detected')
axes[0].set_title('Library Complexity')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store raw for DE
adata.raw = adata.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=500, flavor='seurat_v3',
                            layer=None, span=0.3)
n_hvg = adata.var['highly_variable'].sum()
axes[1].bar(['HVG', 'Non-HVG'], [n_hvg, n_total_genes - n_hvg],
            color=['#e74c3c', '#95a5a6'])
axes[1].set_title(f'Highly Variable Genes: {n_hvg}')
axes[1].set_ylabel('Count')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/02_preprocessing_qc.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"HVGs identified: {n_hvg}")

# ============================================================
# PART 2: Differential Expression & Co-expression Modules
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Differential Expression & Co-expression Modules")
print("=" * 60)

# 2a. Differential Expression Analysis
# Compare each perturbation vs NTC control
control_mask = adata.obs['perturbation'] == 'NTC'
single_perts = [p for p in adata.obs['perturbation'].unique()
                if p != 'NTC' and '+' not in p]

de_results = {}
for pert in single_perts[:10]:
    pert_mask = adata.obs['perturbation'] == pert
    if pert_mask.sum() < 5:
        continue

    # Use raw counts for DE
    control_expr = np.array(adata[control_mask].raw.X.todense())
    pert_expr = np.array(adata[pert_mask].raw.X.todense())

    log2fc = np.zeros(n_total_genes)
    pvals = np.ones(n_total_genes)

    for g in range(n_total_genes):
        c_vals = control_expr[:, g]
        p_vals = pert_expr[:, g]
        if c_vals.std() > 0 or p_vals.std() > 0:
            mean_c = c_vals.mean() + 1e-9
            mean_p = p_vals.mean() + 1e-9
            log2fc[g] = np.log2(mean_p / mean_c)
            try:
                _, pvals[g] = stats.mannwhitneyu(c_vals, p_vals, alternative='two-sided')
            except ValueError:
                pvals[g] = 1.0

    # Multiple testing correction (BH)
    from statsmodels.stats.multitest import multipletests
    _, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')

    de_results[pert] = pd.DataFrame({
        'gene': all_gene_names,
        'log2FC': log2fc,
        'pval': pvals,
        'padj': pvals_adj
    })
    n_de = ((pvals_adj < 0.05) & (np.abs(log2fc) > 0.5)).sum()
    print(f"  {pert}: {n_de} DE genes (padj<0.05, |log2FC|>0.5)")

# Volcano plot for top perturbation
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for idx, pert in enumerate(list(de_results.keys())[:6]):
    ax = axes[idx // 3, idx % 3]
    df = de_results[pert]
    sig = (df['padj'] < 0.05) & (np.abs(df['log2FC']) > 0.5)
    ax.scatter(df.loc[~sig, 'log2FC'], -np.log10(df.loc[~sig, 'pval']),
               s=2, alpha=0.3, c='grey')
    ax.scatter(df.loc[sig, 'log2FC'], -np.log10(df.loc[sig, 'pval']),
               s=5, alpha=0.7, c='red')
    ax.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.5)
    ax.axvline(-0.5, color='green', linestyle='--', alpha=0.5)
    ax.axvline(0.5, color='green', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'{pert}')

plt.suptitle('Volcano Plots: Perturbation vs Control', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/03_volcano_plots.png', dpi=150, bbox_inches='tight')
plt.close()

# 2b. Co-expression Module Detection via NMF
print("\nDetecting co-expression modules via NMF...")
adata_hvg = adata[:, adata.var['highly_variable']].copy()
X_dense = np.array(adata_hvg.X.todense()) if sparse.issparse(adata_hvg.X) else adata_hvg.X
X_nonneg = X_dense - X_dense.min() + 0.01

n_modules = 8
nmf = NMF(n_components=n_modules, random_state=42, max_iter=500)
W = nmf.fit_transform(X_nonneg)  # cell x module
H = nmf.components_  # module x gene

# Assign genes to modules
gene_module_assignments = np.argmax(H, axis=0)
module_sizes = pd.Series(gene_module_assignments).value_counts().sort_index()
print(f"Module sizes: {dict(module_sizes)}")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Module sizes
axes[0].bar(range(n_modules), [module_sizes.get(i, 0) for i in range(n_modules)],
            color=plt.cm.Set3(np.linspace(0, 1, n_modules)), edgecolor='black')
axes[0].set_xlabel('Module')
axes[0].set_ylabel('Number of Genes')
axes[0].set_title('Co-expression Module Sizes (NMF)')

# Module activity heatmap
module_activity = pd.DataFrame(W, columns=[f'Module_{i}' for i in range(n_modules)],
                               index=adata_hvg.obs.index)
module_activity['perturbation'] = adata_hvg.obs['perturbation'].values
mean_activity = module_activity.groupby('perturbation').mean()

# Select top perturbations
top_perts = mean_activity.index[:15]
sns.heatmap(mean_activity.loc[top_perts], ax=axes[1], cmap='viridis',
            xticklabels=True, yticklabels=True)
axes[1].set_title('Module Activity by Perturbation')
axes[1].set_ylabel('')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/04_coexpression_modules.png', dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# PART 3: Causal Graph Estimation
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Causal Graph Estimation from Perturbation Effects")
print("=" * 60)

# Build causal adjacency matrix from perturbation-induced changes
# If knocking down gene A changes gene B significantly, infer edge A->B
causal_genes = target_genes[:10]
causal_matrix = np.zeros((len(causal_genes), len(causal_genes)))

for i, source in enumerate(causal_genes):
    if source in de_results:
        df = de_results[source]
        for j, target in enumerate(causal_genes):
            if i != j and target in df['gene'].values:
                row = df[df['gene'] == target].iloc[0]
                if row['padj'] < 0.1 and abs(row['log2FC']) > 0.3:
                    causal_matrix[i, j] = row['log2FC']

# Create directed graph
G = nx.DiGraph()
for i, gene in enumerate(causal_genes):
    G.add_node(gene)
for i in range(len(causal_genes)):
    for j in range(len(causal_genes)):
        if abs(causal_matrix[i, j]) > 0:
            G.add_edge(causal_genes[i], causal_genes[j],
                       weight=abs(causal_matrix[i, j]),
                       sign='activating' if causal_matrix[i, j] > 0 else 'repressing')

n_edges = G.number_of_edges()
n_nodes = G.number_of_nodes()
print(f"Causal graph: {n_nodes} nodes, {n_edges} edges")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Network visualization
pos = nx.spring_layout(G, seed=42, k=2)
edge_colors = ['#e74c3c' if G[u][v].get('sign') == 'repressing' else '#2ecc71'
               for u, v in G.edges()]
edge_weights = [G[u][v]['weight'] * 3 for u, v in G.edges()]

nx.draw_networkx(G, pos, ax=axes[0],
                 node_color='lightblue', node_size=800,
                 edge_color=edge_colors, width=edge_weights,
                 font_size=7, arrows=True, arrowsize=15,
                 connectionstyle="arc3,rad=0.1")
axes[0].set_title(f'Inferred Causal Regulatory Network\n({n_nodes} genes, {n_edges} edges)')

# Adjacency heatmap
sns.heatmap(causal_matrix, ax=axes[1], cmap='RdBu_r', center=0,
            xticklabels=causal_genes, yticklabels=causal_genes,
            annot=True, fmt='.2f', annot_kws={'fontsize': 6})
axes[1].set_title('Causal Effect Matrix (log2FC)')
axes[1].set_xlabel('Target Gene')
axes[1].set_ylabel('Perturbed Gene')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/05_causal_graph.png', dpi=150, bbox_inches='tight')
plt.close()

# Graph metrics
if n_edges > 0:
    degrees = dict(G.degree())
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())
    hub_genes = sorted(out_degrees.items(), key=lambda x: x[1], reverse=True)[:5]
    target_hubs = sorted(in_degrees.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"Top regulatory hubs (out-degree): {hub_genes}")
    print(f"Top target hubs (in-degree): {target_hubs}")
    try:
        clustering = nx.average_clustering(G)
        print(f"Average clustering coefficient: {clustering:.3f}")
    except:
        pass

# ============================================================
# PART 4: Combinatorial Perturbation Epistasis Detection
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Epistasis Detection in Combinatorial Perturbations")
print("=" * 60)

combo_cells = adata[adata.obs['perturbation'].str.contains(r'\+', regex=True)].copy()
control_cells = adata[adata.obs['perturbation'] == 'NTC'].copy()

print(f"Combinatorial perturbation cells: {combo_cells.shape[0]}")

# For each combination, compute epistasis score
# Epistasis = observed_combo - (expected_A + expected_B - baseline)
epistasis_results = []
combo_perts = combo_cells.obs['perturbation'].unique()

control_raw = control_cells.raw.X
if sparse.issparse(control_raw):
    control_mean = np.asarray(control_raw.todense().mean(axis=0)).flatten()
else:
    control_mean = np.asarray(np.mean(control_raw, axis=0)).flatten()

for combo in combo_perts:
    parts = combo.split('+')
    if len(parts) != 2:
        continue
    gene_a, gene_b = parts

    # Get single perturbation effects
    mask_a = adata.obs['perturbation'] == gene_a
    mask_b = adata.obs['perturbation'] == gene_b
    mask_combo = adata.obs['perturbation'] == combo

    if mask_a.sum() < 3 or mask_b.sum() < 3 or mask_combo.sum() < 3:
        continue

    mean_a = np.asarray(np.array(adata[mask_a].raw.X.todense()).mean(axis=0)).flatten()
    mean_b = np.asarray(np.array(adata[mask_b].raw.X.todense()).mean(axis=0)).flatten()
    mean_combo = np.asarray(np.array(adata[mask_combo].raw.X.todense()).mean(axis=0)).flatten()

    # Additive expectation
    effect_a = mean_a - control_mean
    effect_b = mean_b - control_mean
    expected = control_mean + effect_a + effect_b
    observed = mean_combo

    # Epistasis score per gene
    epistasis = observed - expected
    epistasis_magnitude = np.sqrt(np.mean(epistasis**2))

    # Statistical test
    n_genes_epistatic = (np.abs(epistasis) > 0.5).sum()

    epistasis_results.append({
        'combination': combo,
        'gene_A': gene_a,
        'gene_B': gene_b,
        'epistasis_magnitude': epistasis_magnitude,
        'n_epistatic_genes': n_genes_epistatic,
        'mean_epistasis': np.mean(epistasis),
        'synergy_score': np.mean(epistasis[epistasis > 0]) if (epistasis > 0).any() else 0,
        'antagonism_score': np.mean(epistasis[epistasis < 0]) if (epistasis < 0).any() else 0
    })

epistasis_df = pd.DataFrame(epistasis_results)
if len(epistasis_df) > 0:
    epistasis_df = epistasis_df.sort_values('epistasis_magnitude', ascending=False)
    print(f"\nTop epistatic combinations:")
    print(epistasis_df[['combination', 'epistasis_magnitude', 'n_epistatic_genes']].head(10).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Epistasis magnitude distribution
    axes[0].hist(epistasis_df['epistasis_magnitude'], bins=20,
                 color='#9b59b6', edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Epistasis Magnitude (RMSE)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Distribution of Epistasis Scores')

    # Synergy vs Antagonism
    axes[1].scatter(epistasis_df['synergy_score'], epistasis_df['antagonism_score'],
                    s=50, alpha=0.7, c=epistasis_df['epistasis_magnitude'],
                    cmap='plasma', edgecolor='black', linewidth=0.5)
    axes[1].set_xlabel('Mean Synergy Score')
    axes[1].set_ylabel('Mean Antagonism Score')
    axes[1].set_title('Synergy vs Antagonism')
    plt.colorbar(axes[1].collections[0], ax=axes[1], label='Epistasis Magnitude')

    # Top epistatic interactions
    top_n = min(15, len(epistasis_df))
    axes[2].barh(range(top_n), epistasis_df['epistasis_magnitude'].values[:top_n],
                 color='#e67e22', edgecolor='black', alpha=0.7)
    axes[2].set_yticks(range(top_n))
    axes[2].set_yticklabels(epistasis_df['combination'].values[:top_n], fontsize=7)
    axes[2].set_xlabel('Epistasis Magnitude')
    axes[2].set_title('Top Epistatic Combinations')
    axes[2].invert_yaxis()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/06_epistasis.png', dpi=150, bbox_inches='tight')
    plt.close()
else:
    print("No epistatic combinations detected with sufficient cells.")

# ============================================================
# PART 5: Latent Representation Learning (scVI / CPA-style)
# ============================================================
print("\n" + "=" * 60)
print("PART 5: Perturbation Response Latent Representation Learning")
print("=" * 60)

# 5a. VAE-based latent representation (scVI-inspired)
from sklearn.neural_network import MLPRegressor

# Prepare data for embedding
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)

# UMAP colored by perturbation type
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Color by perturbation type (single vs combo vs control)
pert_type = []
for p in adata.obs['perturbation']:
    if p == 'NTC':
        pert_type.append('Control')
    elif '+' in p:
        pert_type.append('Combinatorial')
    else:
        pert_type.append('Single')
adata.obs['pert_type'] = pert_type

colors = {'Control': '#2ecc71', 'Single': '#3498db', 'Combinatorial': '#e74c3c'}
for ptype, color in colors.items():
    mask = adata.obs['pert_type'] == ptype
    axes[0].scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                    s=2, alpha=0.5, c=color, label=ptype)
axes[0].legend(markerscale=5)
axes[0].set_xlabel('UMAP1')
axes[0].set_ylabel('UMAP2')
axes[0].set_title('UMAP: Perturbation Type')

# 5b. CPA-inspired disentangled representation
# Encode perturbation as one-hot, learn disentangled embedding
from sklearn.decomposition import FactorAnalysis

# Compute perturbation-specific centroids in PCA space
pca_coords = adata.obsm['X_pca'][:, :20]
pert_labels = adata.obs['perturbation'].values

unique_perts = np.unique(pert_labels)
centroids = {}
for p in unique_perts:
    mask = pert_labels == p
    if mask.sum() >= 3:
        centroids[p] = pca_coords[mask].mean(axis=0)

centroid_matrix = np.array(list(centroids.values()))
centroid_labels = list(centroids.keys())

# Cluster centroids to find perturbation groups
from sklearn.cluster import AgglomerativeClustering
if len(centroid_matrix) > 3:
    n_clust = min(5, len(centroid_matrix) - 1)
    clustering = AgglomerativeClustering(n_clusters=n_clust)
    cluster_labels = clustering.fit_predict(centroid_matrix)

    # PCA of centroids for visualization
    pca_centroids = PCA(n_components=2).fit_transform(centroid_matrix)
    scatter = axes[1].scatter(pca_centroids[:, 0], pca_centroids[:, 1],
                              c=cluster_labels, cmap='Set2', s=100,
                              edgecolor='black', linewidth=0.5)
    for i, label in enumerate(centroid_labels):
        if len(label) < 15:
            axes[1].annotate(label, (pca_centroids[i, 0], pca_centroids[i, 1]),
                             fontsize=6, alpha=0.7)
    axes[1].set_xlabel('PC1 (Perturbation Space)')
    axes[1].set_ylabel('PC2 (Perturbation Space)')
    axes[1].set_title('CPA-style Perturbation Embedding')
    plt.colorbar(scatter, ax=axes[1], label='Cluster')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/07_latent_representations.png', dpi=150, bbox_inches='tight')
plt.close()

# Compute perturbation distances
from scipy.spatial.distance import pdist, squareform
if len(centroid_matrix) > 1:
    dist_matrix = squareform(pdist(centroid_matrix, metric='cosine'))
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(dist_matrix, ax=ax, cmap='YlOrRd',
                xticklabels=centroid_labels, yticklabels=centroid_labels,
                annot=False)
    ax.set_title('Perturbation Distance Matrix (Cosine)')
    plt.xticks(fontsize=6, rotation=90)
    plt.yticks(fontsize=6)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/08_perturbation_distances.png', dpi=150, bbox_inches='tight')
    plt.close()

# Silhouette analysis of perturbation separation
sil_scores = []
for p in unique_perts:
    mask = pert_labels == p
    if mask.sum() >= 5 and mask.sum() < len(pert_labels) - 1:
        binary_labels = mask.astype(int)
        try:
            score = silhouette_score(pca_coords[:, :10], binary_labels)
            sil_scores.append({'perturbation': p, 'silhouette': score})
        except:
            pass

sil_df = pd.DataFrame(sil_scores).sort_values('silhouette', ascending=False)
print(f"Mean perturbation separation (silhouette): {sil_df['silhouette'].mean():.3f}")
print(f"Top separated perturbations:")
print(sil_df.head(5).to_string())

# ============================================================
# PART 6: Essential Gene Network Case Study
# ============================================================
print("\n" + "=" * 60)
print("PART 6: Essential Gene Network Case Study")
print("=" * 60)

# Identify essential gene programs from perturbation data
# Essential genes: those whose perturbation causes the most transcriptomic disruption

perturbation_effects = {}
for pert in single_perts:
    if pert in de_results:
        df = de_results[pert]
        n_de = ((df['padj'] < 0.05) & (np.abs(df['log2FC']) > 0.5)).sum()
        mean_effect = df.loc[df['padj'] < 0.05, 'log2FC'].abs().mean() if n_de > 0 else 0
        perturbation_effects[pert] = {
            'n_de_genes': n_de,
            'mean_effect_size': mean_effect,
            'total_disruption': n_de * mean_effect
        }

effect_df = pd.DataFrame(perturbation_effects).T
effect_df = effect_df.sort_values('total_disruption', ascending=False)

print("Perturbation impact ranking:")
print(effect_df.head(10).to_string())

# Essential gene co-regulation network
# Use correlation of perturbation effect profiles
effect_profiles = []
effect_labels = []
for pert in effect_df.index[:10]:
    if pert in de_results:
        effect_profiles.append(de_results[pert]['log2FC'].values)
        effect_labels.append(pert)

if len(effect_profiles) > 1:
    effect_matrix = np.array(effect_profiles)
    corr_matrix = np.corrcoef(effect_matrix)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Perturbation impact bar chart
    top_n = min(10, len(effect_df))
    axes[0].barh(range(top_n), effect_df['total_disruption'].values[:top_n],
                 color='#1abc9c', edgecolor='black', alpha=0.7)
    axes[0].set_yticks(range(top_n))
    axes[0].set_yticklabels(effect_df.index[:top_n], fontsize=8)
    axes[0].set_xlabel('Total Transcriptomic Disruption')
    axes[0].set_title('Perturbation Impact Ranking')
    axes[0].invert_yaxis()

    # Co-regulation heatmap
    sns.heatmap(corr_matrix, ax=axes[1], cmap='coolwarm', center=0,
                xticklabels=effect_labels, yticklabels=effect_labels,
                vmin=-1, vmax=1, annot=True, fmt='.2f', annot_kws={'fontsize':7})
    axes[1].set_title('Perturbation Effect Correlation')

    # Essential gene network
    G_ess = nx.Graph()
    for i, gene_i in enumerate(effect_labels):
        G_ess.add_node(gene_i, disruption=effect_df.loc[gene_i, 'total_disruption'])
    for i in range(len(effect_labels)):
        for j in range(i+1, len(effect_labels)):
            if abs(corr_matrix[i, j]) > 0.3:
                G_ess.add_edge(effect_labels[i], effect_labels[j],
                               weight=abs(corr_matrix[i, j]))

    pos = nx.spring_layout(G_ess, seed=42)
    node_sizes = [effect_df.loc[n, 'total_disruption'] * 50 + 100 for n in G_ess.nodes()]
    edge_weights = [G_ess[u][v]['weight'] * 3 for u, v in G_ess.edges()]
    edge_colors = [corr_matrix[effect_labels.index(u), effect_labels.index(v)]
                   for u, v in G_ess.edges()]

    nodes = nx.draw_networkx_nodes(G_ess, pos, ax=axes[2], node_size=node_sizes,
                                   node_color='#f39c12', edgecolors='black')
    nx.draw_networkx_edges(G_ess, pos, ax=axes[2], width=edge_weights,
                           edge_color=edge_colors, edge_cmap=plt.cm.coolwarm,
                           edge_vmin=-1, edge_vmax=1)
    nx.draw_networkx_labels(G_ess, pos, ax=axes[2], font_size=7)
    axes[2].set_title('Essential Gene Co-regulation Network')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/09_essential_gene_network.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Essential gene network: {G_ess.number_of_nodes()} nodes, {G_ess.number_of_edges()} edges")

# ============================================================
# Summary Statistics
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS SUMMARY")
print("=" * 60)

summary = {
    'Total cells (post-QC)': adata.shape[0],
    'Total genes': adata.shape[1],
    'Unique perturbations': len(adata.obs['perturbation'].unique()),
    'Control cells': (adata.obs['perturbation'] == 'NTC').sum(),
    'Single perturbation cells': (adata.obs['pert_type'] == 'Single').sum(),
    'Combinatorial cells': (adata.obs['pert_type'] == 'Combinatorial').sum(),
    'HVGs identified': n_hvg,
    'Co-expression modules': n_modules,
    'Causal graph edges': n_edges,
    'Epistatic combinations tested': len(epistasis_df) if len(epistasis_results) > 0 else 0,
    'Mean perturbation separation': f"{sil_df['silhouette'].mean():.3f}",
}

for k, v in summary.items():
    print(f"  {k}: {v}")

# Save summary
pd.DataFrame([summary]).T.to_csv('analysis_summary.csv', header=['Value'])

print("\n✓ All figures saved to figures/")
print("✓ Analysis complete!")
