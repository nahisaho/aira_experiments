#!/usr/bin/env python3
"""
Multi-omics Single-Cell Integration Pipeline
=============================================
Integrates scRNA-seq, scATAC-seq, and methylation data using:
1. Preprocessing (QC, normalization, dimensionality reduction)
2. Anchor-based cross-modality integration
3. VAE-based latent space integration
4. Cell lineage inference (RNA velocity + pseudotime)
5. Gene regulatory network (GRN) inference comparison
6. Tumor microenvironment immune cell subtype classification
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import (
    silhouette_score, adjusted_rand_score, normalized_mutual_info_score,
    confusion_matrix, classification_report
)
from sklearn.decomposition import PCA
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

np.random.seed(42)
torch.manual_seed(42)

FIGURE_DIR = "figures"
os.makedirs(FIGURE_DIR, exist_ok=True)

sc.settings.figdir = FIGURE_DIR
sc.settings.verbosity = 1

# ============================================================
# 1. Simulate realistic multi-omics single-cell data
# ============================================================

CELL_TYPES = [
    'CD8+ T cell', 'CD4+ T cell', 'Treg', 'NK cell',
    'B cell', 'Macrophage M1', 'Macrophage M2',
    'Dendritic cell', 'Fibroblast', 'Tumor cell'
]

N_CELLS = 1500
N_GENES_RNA = 800
N_PEAKS_ATAC = 600
N_CPG_SITES = 400

def simulate_cell_type_programs(n_cell_types, n_features, sparsity=0.7):
    """Generate cell-type-specific feature programs."""
    programs = np.random.randn(n_cell_types, n_features) * 0.5
    mask = np.random.rand(n_cell_types, n_features) < sparsity
    programs[mask] = 0
    for i in range(n_cell_types):
        marker_start = i * (n_features // n_cell_types)
        marker_end = marker_start + max(n_features // (n_cell_types * 2), 10)
        marker_end = min(marker_end, n_features)
        programs[i, marker_start:marker_end] = np.abs(np.random.randn(marker_end - marker_start)) * 2
    return programs

def simulate_pseudotime_trajectory(n_cells, n_branches=3):
    """Simulate branching pseudotime trajectory."""
    pseudotime = np.zeros(n_cells)
    branch_labels = np.zeros(n_cells, dtype=int)
    cells_per_branch = n_cells // n_branches
    for b in range(n_branches):
        start = b * cells_per_branch
        end = start + cells_per_branch if b < n_branches - 1 else n_cells
        pseudotime[start:end] = np.sort(np.random.beta(2, 5, end - start))
        branch_labels[start:end] = b
    return pseudotime, branch_labels

def generate_multiomics_data():
    """Generate synthetic scRNA-seq, scATAC-seq, and methylation data."""
    print("Generating synthetic multi-omics data...")

    cell_type_assignments = np.random.choice(len(CELL_TYPES), N_CELLS, 
        p=[0.15, 0.12, 0.05, 0.08, 0.10, 0.08, 0.07, 0.05, 0.10, 0.20])
    cell_type_labels = np.array([CELL_TYPES[i] for i in cell_type_assignments])

    pseudotime, branch_labels = simulate_pseudotime_trajectory(N_CELLS)

    # --- scRNA-seq ---
    rna_programs = simulate_cell_type_programs(len(CELL_TYPES), N_GENES_RNA, sparsity=0.6)
    rna_base = rna_programs[cell_type_assignments]
    trajectory_effect = np.outer(pseudotime, np.random.randn(N_GENES_RNA) * 0.3)
    rna_counts = np.exp(rna_base + trajectory_effect + np.random.randn(N_CELLS, N_GENES_RNA) * 0.3)
    rna_counts = np.random.poisson(rna_counts).astype(np.float32)

    gene_names = [f"Gene_{i}" for i in range(N_GENES_RNA)]
    immune_markers = {
        'CD3E': 0, 'CD8A': 1, 'CD4': 2, 'FOXP3': 3, 'NKG7': 4,
        'CD19': 5, 'CD68': 6, 'CD163': 7, 'CLEC9A': 8, 'COL1A1': 9,
        'EPCAM': 10, 'IL2RA': 11, 'GZMB': 12, 'PDCD1': 13, 'CTLA4': 14,
        'MKI67': 15, 'TP53': 16, 'BRCA1': 17, 'VEGFA': 18, 'HIF1A': 19
    }
    for marker, idx in immune_markers.items():
        if idx < N_GENES_RNA:
            gene_names[idx] = marker

    cell_ids = [f"Cell_{i}" for i in range(N_CELLS)]
    adata_rna = ad.AnnData(
        X=sparse.csr_matrix(rna_counts),
        obs=pd.DataFrame({
            'cell_type': cell_type_labels,
            'pseudotime': pseudotime,
            'branch': branch_labels
        }, index=cell_ids),
        var=pd.DataFrame(index=gene_names)
    )

    # --- scATAC-seq ---
    atac_programs = simulate_cell_type_programs(len(CELL_TYPES), N_PEAKS_ATAC, sparsity=0.75)
    atac_base = atac_programs[cell_type_assignments]
    atac_trajectory = np.outer(pseudotime, np.random.randn(N_PEAKS_ATAC) * 0.2)
    atac_signal = 1 / (1 + np.exp(-(atac_base + atac_trajectory + np.random.randn(N_CELLS, N_PEAKS_ATAC) * 0.5)))
    atac_binary = (atac_signal > 0.5).astype(np.float32)

    peak_names = [f"chr{np.random.randint(1,23)}:{np.random.randint(1e6,1e8)}-{np.random.randint(1e6,1e8)}" 
                  for _ in range(N_PEAKS_ATAC)]
    adata_atac = ad.AnnData(
        X=sparse.csr_matrix(atac_binary),
        obs=pd.DataFrame({
            'cell_type': cell_type_labels,
            'pseudotime': pseudotime,
            'branch': branch_labels
        }, index=cell_ids),
        var=pd.DataFrame(index=peak_names)
    )

    # --- Methylation ---
    meth_programs = simulate_cell_type_programs(len(CELL_TYPES), N_CPG_SITES, sparsity=0.65)
    meth_base = meth_programs[cell_type_assignments]
    meth_trajectory = np.outer(pseudotime, np.random.randn(N_CPG_SITES) * 0.15)
    meth_beta = 1 / (1 + np.exp(-(meth_base + meth_trajectory + np.random.randn(N_CELLS, N_CPG_SITES) * 0.4)))

    cpg_names = [f"CpG_{i}" for i in range(N_CPG_SITES)]
    adata_meth = ad.AnnData(
        X=meth_beta.astype(np.float32),
        obs=pd.DataFrame({
            'cell_type': cell_type_labels,
            'pseudotime': pseudotime,
            'branch': branch_labels
        }, index=cell_ids),
        var=pd.DataFrame(index=cpg_names)
    )

    print(f"  RNA-seq: {adata_rna.shape}")
    print(f"  ATAC-seq: {adata_atac.shape}")
    print(f"  Methylation: {adata_meth.shape}")

    return adata_rna, adata_atac, adata_meth

# ============================================================
# 2. Preprocessing
# ============================================================

def preprocess_rna(adata):
    """QC, normalization, dimensionality reduction for scRNA-seq."""
    print("\n=== Preprocessing scRNA-seq ===")
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()

    adata.layers['raw_counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log_normalized'] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=min(1500, adata.n_vars), flavor='seurat_v3',
                                 layer='raw_counts')
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50)
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)

    print(f"  After QC: {adata.shape}")
    print(f"  HVGs: {adata.var['highly_variable'].sum()}")
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")
    return adata

def preprocess_atac(adata):
    """QC, normalization, dimensionality reduction for scATAC-seq."""
    print("\n=== Preprocessing scATAC-seq ===")
    adata.obs['n_peaks'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    sc.pp.filter_cells(adata, min_genes=100)

    # TF-IDF normalization
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X.copy()
    tf = X / (X.sum(axis=1, keepdims=True) + 1e-8)
    idf = np.log1p(X.shape[0] / (X.sum(axis=0) + 1e-8))
    tfidf = tf * idf
    adata.layers['tfidf'] = tfidf

    # LSI (Latent Semantic Indexing) via SVD
    from sklearn.decomposition import TruncatedSVD
    svd = TruncatedSVD(n_components=50, random_state=42)
    lsi = svd.fit_transform(tfidf)
    adata.obsm['X_lsi'] = lsi[:, 1:]  # Remove first component (depth-correlated)

    sc.pp.neighbors(adata, use_rep='X_lsi', n_neighbors=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)

    print(f"  After QC: {adata.shape}")
    print(f"  LSI components: {adata.obsm['X_lsi'].shape[1]}")
    return adata

def preprocess_methylation(adata):
    """QC, normalization, dimensionality reduction for methylation data."""
    print("\n=== Preprocessing Methylation ===")
    X = adata.X if not sparse.issparse(adata.X) else adata.X.toarray()

    # Filter low-variance CpG sites
    site_var = X.var(axis=0)
    high_var_mask = site_var > np.percentile(site_var, 25)
    adata = adata[:, high_var_mask].copy()
    X = adata.X if not sparse.issparse(adata.X) else adata.X.toarray()

    # M-value transformation: M = log2(beta / (1 - beta))
    beta_clipped = np.clip(X, 0.01, 0.99)
    m_values = np.log2(beta_clipped / (1 - beta_clipped))
    adata.layers['m_values'] = m_values

    scaler = StandardScaler()
    m_scaled = scaler.fit_transform(m_values)
    pca = PCA(n_components=50, random_state=42)
    adata.obsm['X_pca'] = pca.fit_transform(m_scaled)

    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)

    print(f"  After QC: {adata.shape}")
    print(f"  PCA components: {adata.obsm['X_pca'].shape[1]}")
    return adata

# ============================================================
# 3. Anchor-based integration
# ============================================================

def find_anchors(X1, X2, k=30):
    """Find mutual nearest neighbors as anchors between two modalities."""
    print("  Finding MNN anchors...")
    nn1 = NearestNeighbors(n_neighbors=k).fit(X1)
    nn2 = NearestNeighbors(n_neighbors=k).fit(X2)

    dist12, idx12 = nn2.kneighbors(X1)
    dist21, idx21 = nn1.kneighbors(X2)

    anchors = []
    for i in range(X1.shape[0]):
        for j_idx in range(min(5, k)):
            j = idx12[i, j_idx]
            if i in idx21[j, :5]:
                anchors.append((i, j))

    print(f"  Found {len(anchors)} MNN anchor pairs")
    return anchors

def anchor_based_integration(adata_rna, adata_atac, adata_meth):
    """Integrate modalities using anchor-based approach (MNN)."""
    print("\n=== Anchor-based Integration ===")

    # Use reduced representations
    X_rna = adata_rna.obsm['X_pca'][:, :30]
    X_atac = adata_atac.obsm['X_lsi'][:, :30]
    X_meth = adata_meth.obsm['X_pca'][:, :30]

    # Find anchors between pairs
    anchors_rna_atac = find_anchors(X_rna, X_atac)
    anchors_rna_meth = find_anchors(X_rna, X_meth)

    # Compute correction vectors from anchors
    def compute_correction(X_ref, X_query, anchors):
        if len(anchors) == 0:
            return np.zeros_like(X_query)
        anchor_arr = np.array(anchors)
        corrections = X_ref[anchor_arr[:, 0]] - X_query[anchor_arr[:, 1]]
        mean_correction = corrections.mean(axis=0)
        return X_query + mean_correction

    X_atac_corrected = compute_correction(X_rna, X_atac, anchors_rna_atac)
    X_meth_corrected = compute_correction(X_rna, X_meth, anchors_rna_meth)

    # Concatenate for integrated space
    integrated = np.concatenate([X_rna, X_atac_corrected, X_meth_corrected], axis=1)
    pca = PCA(n_components=30, random_state=42)
    integrated_pca = pca.fit_transform(StandardScaler().fit_transform(integrated))

    print(f"  Integrated space: {integrated_pca.shape}")
    return integrated_pca, anchors_rna_atac, anchors_rna_meth

# ============================================================
# 4. VAE-based Integration
# ============================================================

class MultiOmicsVAE(nn.Module):
    """Variational Autoencoder for multi-omics integration."""
    def __init__(self, input_dims, latent_dim=20, hidden_dim=128):
        super().__init__()
        self.latent_dim = latent_dim
        total_input = sum(input_dims)

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(total_input, hidden_dim * 2),
            nn.BatchNorm1d(hidden_dim * 2),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
        )
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

        # Modality-specific decoders
        self.decoders = nn.ModuleList()
        for dim in input_dims:
            self.decoders.append(nn.Sequential(
                nn.Linear(latent_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(0.1),
                nn.Linear(hidden_dim, hidden_dim * 2),
                nn.BatchNorm1d(hidden_dim * 2),
                nn.ReLU(),
                nn.Linear(hidden_dim * 2, dim),
            ))

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def decode(self, z):
        return [decoder(z) for decoder in self.decoders]

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recons = self.decode(z)
        return recons, mu, logvar, z

def vae_loss(recons, x_splits, mu, logvar, beta=1.0):
    """ELBO loss with modality-specific reconstruction."""
    recon_loss = 0
    for recon, x in zip(recons, x_splits):
        recon_loss += nn.functional.mse_loss(recon, x, reduction='sum')
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + beta * kl_loss, recon_loss, kl_loss

def train_vae(adata_rna, adata_atac, adata_meth, latent_dim=20, epochs=50, batch_size=128, lr=1e-3):
    """Train multi-omics VAE."""
    print("\n=== Training Multi-Omics VAE ===")

    X_rna = adata_rna.obsm['X_pca'][:, :30].astype(np.float32)
    X_atac = adata_atac.obsm['X_lsi'][:, :30].astype(np.float32)
    X_meth = adata_meth.obsm['X_pca'][:, :30].astype(np.float32)

    scaler_rna = StandardScaler().fit(X_rna)
    scaler_atac = StandardScaler().fit(X_atac)
    scaler_meth = StandardScaler().fit(X_meth)

    X_rna_s = scaler_rna.transform(X_rna)
    X_atac_s = scaler_atac.transform(X_atac)
    X_meth_s = scaler_meth.transform(X_meth)

    X_combined = np.concatenate([X_rna_s, X_atac_s, X_meth_s], axis=1)
    input_dims = [X_rna_s.shape[1], X_atac_s.shape[1], X_meth_s.shape[1]]

    dataset = TensorDataset(torch.FloatTensor(X_combined))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    model = MultiOmicsVAE(input_dims, latent_dim=latent_dim)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    train_losses = []
    recon_losses = []
    kl_losses = []

    model.train()
    for epoch in range(epochs):
        epoch_loss = 0
        epoch_recon = 0
        epoch_kl = 0
        n_batches = 0

        # KL annealing
        beta = min(1.0, epoch / (epochs * 0.3))

        for (batch,) in loader:
            optimizer.zero_grad()
            x_splits = torch.split(batch, input_dims, dim=1)
            recons, mu, logvar, z = model(batch)
            loss, recon, kl = vae_loss(recons, x_splits, mu, logvar, beta=beta)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            epoch_loss += loss.item()
            epoch_recon += recon.item()
            epoch_kl += kl.item()
            n_batches += 1

        scheduler.step()
        train_losses.append(epoch_loss / n_batches)
        recon_losses.append(epoch_recon / n_batches)
        kl_losses.append(epoch_kl / n_batches)

        if (epoch + 1) % 20 == 0:
            print(f"  Epoch {epoch+1}/{epochs}: Loss={train_losses[-1]:.1f}, "
                  f"Recon={recon_losses[-1]:.1f}, KL={kl_losses[-1]:.1f}")

    # Extract latent space
    model.eval()
    with torch.no_grad():
        X_tensor = torch.FloatTensor(X_combined)
        h = model.encoder(X_tensor)
        mu = model.fc_mu(h)
        latent = mu.numpy()

    print(f"  Latent space: {latent.shape}")
    return latent, model, train_losses, recon_losses, kl_losses

# ============================================================
# 5. Cell Lineage Inference
# ============================================================

def simulate_rna_velocity(adata_rna):
    """Simulate RNA velocity using spliced/unspliced ratio model."""
    print("\n=== RNA Velocity & Pseudotime Analysis ===")

    n_cells, n_genes = adata_rna.shape
    pseudotime = adata_rna.obs['pseudotime'].values

    # Simulate spliced/unspliced counts
    X = adata_rna.X.toarray() if sparse.issparse(adata_rna.X) else adata_rna.X
    alpha = np.random.gamma(2, 1, n_genes)  # transcription rates
    beta = np.random.gamma(1, 0.5, n_genes)  # splicing rates
    gamma = np.random.gamma(1, 0.5, n_genes)  # degradation rates

    u_ss = alpha / beta  # unspliced steady state
    s_ss = alpha / gamma  # spliced steady state

    # Model velocity as change in spliced counts along pseudotime
    pt_expanded = pseudotime[:, np.newaxis]
    velocity = np.zeros((n_cells, n_genes))
    for g in range(min(n_genes, 500)):
        k = np.random.choice([-1, 1]) * np.random.gamma(1, 1)
        velocity[:, g] = k * np.exp(-((pt_expanded[:, 0] - 0.5)**2) / 0.2)

    adata_rna.layers['velocity'] = velocity[:, :n_genes]
    adata_rna.layers['spliced'] = X
    unspliced = np.abs(X * 0.3 + np.random.poisson(0.5, X.shape))
    adata_rna.layers['unspliced'] = unspliced.astype(np.float32)

    # Compute velocity graph (transition probabilities)
    X_pca = adata_rna.obsm['X_pca']
    nn = NearestNeighbors(n_neighbors=30).fit(X_pca)
    _, indices = nn.kneighbors(X_pca)

    velocity_pca = PCA(n_components=50, random_state=42).fit_transform(
        velocity[:, :min(n_genes, 500)])

    # Cosine similarity between velocity and cell-cell displacement
    transition_probs = np.zeros((n_cells, n_cells))
    for i in range(n_cells):
        neighbors = indices[i]
        displacements = X_pca[neighbors] - X_pca[i]
        vel_i = velocity_pca[i]
        cos_sim = np.dot(displacements, vel_i) / (
            np.linalg.norm(displacements, axis=1) * np.linalg.norm(vel_i) + 1e-8)
        cos_sim = np.clip(cos_sim, 0, None)
        if cos_sim.sum() > 0:
            transition_probs[i, neighbors] = cos_sim / cos_sim.sum()

    adata_rna.uns['velocity_graph'] = sparse.csr_matrix(transition_probs)

    # Diffusion pseudotime
    from sklearn.manifold import SpectralEmbedding
    se = SpectralEmbedding(n_components=1, affinity='nearest_neighbors',
                           n_neighbors=15, random_state=42)
    dpt = se.fit_transform(X_pca)
    adata_rna.obs['dpt_pseudotime'] = (dpt[:, 0] - dpt[:, 0].min()) / (dpt[:, 0].max() - dpt[:, 0].min())

    # Correlation with ground truth
    corr_pt, _ = spearmanr(adata_rna.obs['pseudotime'], adata_rna.obs['dpt_pseudotime'])
    print(f"  Pseudotime correlation (Spearman): {corr_pt:.3f}")
    print(f"  Velocity vectors computed for {n_cells} cells")

    return adata_rna

# ============================================================
# 6. GRN Inference Methods Comparison
# ============================================================

def grn_correlation_based(X, gene_names, top_k=500):
    """Simple correlation-based GRN inference."""
    corr_matrix = np.corrcoef(X.T)
    np.fill_diagonal(corr_matrix, 0)
    edges = []
    n = corr_matrix.shape[0]
    for i in range(n):
        for j in range(i+1, n):
            if abs(corr_matrix[i, j]) > 0.3:
                edges.append((gene_names[i], gene_names[j], corr_matrix[i, j]))
    edges.sort(key=lambda x: abs(x[2]), reverse=True)
    return edges[:top_k]

def grn_mutual_information(X, gene_names, top_k=500):
    """Mutual information-based GRN inference (ARACNE-like)."""
    from sklearn.feature_selection import mutual_info_regression
    n_genes = X.shape[1]
    mi_matrix = np.zeros((n_genes, n_genes))
    for i in range(min(n_genes, 50)):
        mi = mutual_info_regression(X[:, :min(n_genes, 50)], X[:, i], random_state=42, n_neighbors=5)
        mi_matrix[i, :min(n_genes, 50)] = mi[:min(n_genes, 50)]
        mi_matrix[:min(n_genes, 50), i] = mi[:min(n_genes, 50)]
    np.fill_diagonal(mi_matrix, 0)

    edges = []
    for i in range(min(n_genes, 50)):
        for j in range(i+1, min(n_genes, 50)):
            if mi_matrix[i, j] > 0.05:
                edges.append((gene_names[i], gene_names[j], mi_matrix[i, j]))
    edges.sort(key=lambda x: x[2], reverse=True)
    return edges[:top_k]

def grn_genie3_like(X, gene_names, top_k=500):
    """Random Forest-based GRN inference (GENIE3-like)."""
    from sklearn.ensemble import RandomForestRegressor
    n_genes = min(X.shape[1], 30)
    importance_matrix = np.zeros((n_genes, n_genes))

    for target in range(n_genes):
        predictors = np.delete(np.arange(n_genes), target)
        X_pred = X[:, predictors]
        y = X[:, target]
        rf = RandomForestRegressor(n_estimators=50, max_depth=5, random_state=42, n_jobs=-1)
        rf.fit(X_pred, y)
        imp = rf.feature_importances_
        for k, p in enumerate(predictors):
            importance_matrix[p, target] = imp[k]

    edges = []
    for i in range(n_genes):
        for j in range(n_genes):
            if i != j and importance_matrix[i, j] > 0.02:
                edges.append((gene_names[i], gene_names[j], importance_matrix[i, j]))
    edges.sort(key=lambda x: x[2], reverse=True)
    return edges[:top_k]

def compare_grn_methods(adata_rna):
    """Compare GRN inference methods."""
    print("\n=== GRN Inference Comparison ===")

    X = adata_rna.X.toarray() if sparse.issparse(adata_rna.X) else adata_rna.X
    gene_names = list(adata_rna.var_names)

    # Use top variable genes for efficiency
    gene_var = X.var(axis=0)
    top_idx = np.argsort(gene_var)[-50:]
    X_sub = X[:, top_idx]
    genes_sub = [gene_names[i] for i in top_idx]

    print("  Method 1: Correlation-based...")
    edges_corr = grn_correlation_based(X_sub, genes_sub)

    print("  Method 2: Mutual Information (ARACNE-like)...")
    edges_mi = grn_mutual_information(X_sub, genes_sub)

    print("  Method 3: Random Forest (GENIE3-like)...")
    edges_rf = grn_genie3_like(X_sub, genes_sub)

    results = {
        'Correlation': {'edges': len(edges_corr), 'top_edges': edges_corr[:10]},
        'Mutual Information': {'edges': len(edges_mi), 'top_edges': edges_mi[:10]},
        'Random Forest': {'edges': len(edges_rf), 'top_edges': edges_rf[:10]},
    }

    for method, res in results.items():
        print(f"  {method}: {res['edges']} edges")

    return results, edges_corr, edges_mi, edges_rf

# ============================================================
# 7. Immune Cell Subtype Classification
# ============================================================

def classify_immune_subtypes(adata_rna, latent_space):
    """Classify immune cell subtypes in tumor microenvironment."""
    print("\n=== Immune Cell Subtype Classification ===")

    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.svm import SVC
    from sklearn.model_selection import StratifiedKFold, cross_val_score

    labels = adata_rna.obs['cell_type'].values
    le = LabelEncoder()
    y = le.fit_transform(labels)

    X = latent_space

    # Compare classifiers
    classifiers = {
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'Gradient Boosting': GradientBoostingClassifier(n_estimators=100, random_state=42),
        'SVM (RBF)': SVC(kernel='rbf', random_state=42),
    }

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results = {}

    for name, clf in classifiers.items():
        scores = cross_val_score(clf, X, y, cv=cv, scoring='accuracy')
        results[name] = {
            'mean_accuracy': scores.mean(),
            'std_accuracy': scores.std(),
            'scores': scores
        }
        print(f"  {name}: {scores.mean():.3f} ± {scores.std():.3f}")

    # Train best classifier for confusion matrix
    best_clf = RandomForestClassifier(n_estimators=200, random_state=42)
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    best_clf.fit(X_train, y_train)
    y_pred = best_clf.predict(X_test)

    accuracy = (y_pred == y_test).mean()
    ari = adjusted_rand_score(y_test, y_pred)
    nmi = normalized_mutual_info_score(y_test, y_pred)

    print(f"\n  Final classifier (RF-200):")
    print(f"    Accuracy: {accuracy:.3f}")
    print(f"    ARI: {ari:.3f}")
    print(f"    NMI: {nmi:.3f}")

    cm = confusion_matrix(y_test, y_pred)
    class_names = le.classes_

    return results, cm, class_names, accuracy, ari, nmi

# ============================================================
# 8. Visualization
# ============================================================

def plot_preprocessing_qc(adata_rna, adata_atac, adata_meth):
    """Plot QC metrics for all modalities."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # RNA QC
    axes[0, 0].hist(adata_rna.obs['n_counts'], bins=50, color='#2196F3', alpha=0.7, edgecolor='black')
    axes[0, 0].set_xlabel('Total Counts')
    axes[0, 0].set_ylabel('Number of Cells')
    axes[0, 0].set_title('scRNA-seq: Count Distribution')

    axes[0, 1].hist(adata_atac.obs['n_peaks'], bins=50, color='#4CAF50', alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('Number of Peaks')
    axes[0, 1].set_ylabel('Number of Cells')
    axes[0, 1].set_title('scATAC-seq: Peak Distribution')

    X_meth = adata_meth.X if not sparse.issparse(adata_meth.X) else adata_meth.X.toarray()
    axes[0, 2].hist(X_meth.mean(axis=1), bins=50, color='#FF9800', alpha=0.7, edgecolor='black')
    axes[0, 2].set_xlabel('Mean β-value')
    axes[0, 2].set_ylabel('Number of Cells')
    axes[0, 2].set_title('Methylation: β-value Distribution')

    # UMAPs by cell type
    for idx, (adata, title, cmap) in enumerate([
        (adata_rna, 'scRNA-seq UMAP', 'tab10'),
        (adata_atac, 'scATAC-seq UMAP', 'tab10'),
        (adata_meth, 'Methylation UMAP', 'tab10')
    ]):
        umap = adata.obsm['X_umap']
        cell_types = adata.obs['cell_type']
        unique_types = cell_types.unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_types)))
        for ct_idx, ct in enumerate(unique_types):
            mask = cell_types == ct
            axes[1, idx].scatter(umap[mask, 0], umap[mask, 1], c=[colors[ct_idx]],
                               label=ct, s=3, alpha=0.6)
        axes[1, idx].set_title(title)
        axes[1, idx].set_xlabel('UMAP1')
        axes[1, idx].set_ylabel('UMAP2')

    axes[1, 0].legend(bbox_to_anchor=(0, -0.3), loc='upper left', ncol=3, fontsize=7, markerscale=3)
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig1_preprocessing_qc.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig1_preprocessing_qc.png")

def plot_integration_comparison(adata_rna, integrated_anchor, latent_vae):
    """Compare anchor-based vs VAE integration."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    cell_types = adata_rna.obs['cell_type'].values
    unique_types = np.unique(cell_types)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_types)))
    color_map = {ct: colors[i] for i, ct in enumerate(unique_types)}

    # Original RNA UMAP
    umap_rna = adata_rna.obsm['X_umap']
    for ct in unique_types:
        mask = cell_types == ct
        axes[0].scatter(umap_rna[mask, 0], umap_rna[mask, 1],
                       c=[color_map[ct]], label=ct, s=5, alpha=0.6)
    axes[0].set_title('Original scRNA-seq UMAP')
    axes[0].set_xlabel('UMAP1')
    axes[0].set_ylabel('UMAP2')

    # Anchor-based integration
    import umap as umap_lib
    reducer = umap_lib.UMAP(n_components=2, random_state=42)
    anchor_umap = reducer.fit_transform(integrated_anchor)
    for ct in unique_types:
        mask = cell_types == ct
        axes[1].scatter(anchor_umap[mask, 0], anchor_umap[mask, 1],
                       c=[color_map[ct]], label=ct, s=5, alpha=0.6)
    axes[1].set_title('Anchor-based Integration UMAP')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')

    # VAE integration
    vae_umap = reducer.fit_transform(latent_vae)
    for ct in unique_types:
        mask = cell_types == ct
        axes[2].scatter(vae_umap[mask, 0], vae_umap[mask, 1],
                       c=[color_map[ct]], label=ct, s=5, alpha=0.6)
    axes[2].set_title('VAE Integration UMAP')
    axes[2].set_xlabel('UMAP1')
    axes[2].set_ylabel('UMAP2')

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=8,
              markerscale=3, bbox_to_anchor=(0.5, -0.08))
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig2_integration_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig2_integration_comparison.png")
    return anchor_umap, vae_umap

def plot_vae_training(train_losses, recon_losses, kl_losses):
    """Plot VAE training curves."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    epochs = range(1, len(train_losses) + 1)
    axes[0].plot(epochs, train_losses, 'b-', linewidth=1.5)
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Total Loss')
    axes[0].set_title('Total ELBO Loss')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(epochs, recon_losses, 'r-', linewidth=1.5)
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('Reconstruction Loss')
    axes[1].set_title('Reconstruction Loss')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(epochs, kl_losses, 'g-', linewidth=1.5)
    axes[2].set_xlabel('Epoch')
    axes[2].set_ylabel('KL Divergence')
    axes[2].set_title('KL Divergence')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig3_vae_training.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig3_vae_training.png")

def plot_pseudotime_velocity(adata_rna):
    """Plot pseudotime and velocity analysis."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    umap = adata_rna.obsm['X_umap']

    # Ground truth pseudotime
    sc1 = axes[0].scatter(umap[:, 0], umap[:, 1], c=adata_rna.obs['pseudotime'],
                         cmap='viridis', s=5, alpha=0.7)
    plt.colorbar(sc1, ax=axes[0], shrink=0.8)
    axes[0].set_title('Ground Truth Pseudotime')
    axes[0].set_xlabel('UMAP1')
    axes[0].set_ylabel('UMAP2')

    # Inferred pseudotime
    sc2 = axes[1].scatter(umap[:, 0], umap[:, 1], c=adata_rna.obs['dpt_pseudotime'],
                         cmap='viridis', s=5, alpha=0.7)
    plt.colorbar(sc2, ax=axes[1], shrink=0.8)
    axes[1].set_title('Inferred Pseudotime (DPT)')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')

    # Pseudotime correlation
    axes[2].scatter(adata_rna.obs['pseudotime'], adata_rna.obs['dpt_pseudotime'],
                   s=3, alpha=0.3, c='steelblue')
    corr, _ = spearmanr(adata_rna.obs['pseudotime'], adata_rna.obs['dpt_pseudotime'])
    axes[2].set_xlabel('Ground Truth Pseudotime')
    axes[2].set_ylabel('Inferred Pseudotime')
    axes[2].set_title(f'Pseudotime Correlation (ρ={corr:.3f})')
    axes[2].plot([0, 1], [0, 1], 'r--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig4_pseudotime_velocity.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig4_pseudotime_velocity.png")

def plot_grn_comparison(grn_results):
    """Plot GRN inference comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    methods = list(grn_results.keys())
    n_edges = [grn_results[m]['edges'] for m in methods]
    colors = ['#2196F3', '#4CAF50', '#FF9800']

    # Number of edges
    axes[0].bar(methods, n_edges, color=colors, edgecolor='black', alpha=0.8)
    axes[0].set_ylabel('Number of Edges')
    axes[0].set_title('GRN: Number of Inferred Edges')
    for i, v in enumerate(n_edges):
        axes[0].text(i, v + 5, str(v), ha='center', fontsize=10, fontweight='bold')

    # Top edge weights distribution
    for idx, method in enumerate(methods):
        weights = [abs(e[2]) for e in grn_results[method]['top_edges']]
        axes[1].bar(np.arange(len(weights)) + idx * 0.25, weights, width=0.25,
                   color=colors[idx], label=method, alpha=0.8)
    axes[1].set_xlabel('Edge Rank')
    axes[1].set_ylabel('Edge Weight')
    axes[1].set_title('Top 10 Edge Weights')
    axes[1].legend(fontsize=8)

    # Venn-like overlap (simplified as grouped bars)
    # Compute pairwise overlap
    edge_sets = {}
    for method_key, edges_list_key in [('Correlation', 'Correlation'),
                                        ('Mutual Information', 'Mutual Information'),
                                        ('Random Forest', 'Random Forest')]:
        edge_sets[method_key] = set(
            (e[0], e[1]) if e[0] < e[1] else (e[1], e[0])
            for e in grn_results[method_key]['top_edges']
        )

    overlap_data = {}
    method_pairs = [('Correlation', 'Mutual Information'),
                    ('Correlation', 'Random Forest'),
                    ('Mutual Information', 'Random Forest')]
    pair_labels = ['Corr∩MI', 'Corr∩RF', 'MI∩RF']
    overlaps = [len(edge_sets[a] & edge_sets[b]) for a, b in method_pairs]
    all_overlap = len(edge_sets['Correlation'] & edge_sets['Mutual Information'] & edge_sets['Random Forest'])

    bars = axes[2].bar(pair_labels + ['All three'], overlaps + [all_overlap],
                      color=['#9C27B0', '#E91E63', '#00BCD4', '#FFC107'],
                      edgecolor='black', alpha=0.8)
    axes[2].set_ylabel('Shared Edges')
    axes[2].set_title('GRN Method Agreement')
    for bar, val in zip(bars, overlaps + [all_overlap]):
        axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                    str(val), ha='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig5_grn_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig5_grn_comparison.png")

def plot_immune_classification(cm, class_names, clf_results):
    """Plot immune cell classification results."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Confusion matrix
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    im = axes[0].imshow(cm_norm, interpolation='nearest', cmap='Blues', vmin=0, vmax=1)
    plt.colorbar(im, ax=axes[0], shrink=0.8)
    axes[0].set_xticks(range(len(class_names)))
    axes[0].set_yticks(range(len(class_names)))
    axes[0].set_xticklabels(class_names, rotation=45, ha='right', fontsize=7)
    axes[0].set_yticklabels(class_names, fontsize=7)
    axes[0].set_xlabel('Predicted')
    axes[0].set_ylabel('True')
    axes[0].set_title('Normalized Confusion Matrix')

    # Add text to cells
    for i in range(len(class_names)):
        for j in range(len(class_names)):
            color = 'white' if cm_norm[i, j] > 0.5 else 'black'
            axes[0].text(j, i, f'{cm_norm[i, j]:.2f}', ha='center', va='center',
                        color=color, fontsize=7)

    # Classifier comparison
    methods = list(clf_results.keys())
    means = [clf_results[m]['mean_accuracy'] for m in methods]
    stds = [clf_results[m]['std_accuracy'] for m in methods]
    colors = ['#2196F3', '#4CAF50', '#FF9800']

    bars = axes[1].bar(methods, means, yerr=stds, color=colors, edgecolor='black',
                      alpha=0.8, capsize=5)
    axes[1].set_ylabel('Accuracy')
    axes[1].set_title('Classifier Comparison (5-fold CV)')
    axes[1].set_ylim(0, 1.0)
    for bar, mean in zip(bars, means):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                    f'{mean:.3f}', ha='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig6_immune_classification.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig6_immune_classification.png")

def plot_integration_metrics(adata_rna, integrated_anchor, latent_vae):
    """Plot integration quality metrics."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cell_types = adata_rna.obs['cell_type'].values
    le = LabelEncoder()
    y = le.fit_transform(cell_types)

    # Compute silhouette scores
    methods = ['RNA PCA', 'Anchor-based', 'VAE Latent']
    embeddings = [adata_rna.obsm['X_pca'][:, :30], integrated_anchor, latent_vae]
    sil_scores = []
    for emb in embeddings:
        sil = silhouette_score(emb, y, sample_size=min(1000, len(y)))
        sil_scores.append(sil)

    colors = ['#2196F3', '#4CAF50', '#FF9800']
    bars = axes[0].bar(methods, sil_scores, color=colors, edgecolor='black', alpha=0.8)
    axes[0].set_ylabel('Silhouette Score')
    axes[0].set_title('Integration Quality: Silhouette Score')
    axes[0].set_ylim(0, max(sil_scores) * 1.3)
    for bar, score in zip(bars, sil_scores):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{score:.3f}', ha='center', fontsize=10, fontweight='bold')

    # Leiden clustering ARI/NMI for each method
    from sklearn.cluster import KMeans
    ari_scores = []
    nmi_scores = []
    for emb in embeddings:
        kmeans = KMeans(n_clusters=len(CELL_TYPES), random_state=42, n_init=10)
        pred = kmeans.fit_predict(emb)
        ari_scores.append(adjusted_rand_score(y, pred))
        nmi_scores.append(normalized_mutual_info_score(y, pred))

    x_pos = np.arange(len(methods))
    width = 0.35
    bars1 = axes[1].bar(x_pos - width/2, ari_scores, width, label='ARI', color='#3F51B5',
                       edgecolor='black', alpha=0.8)
    bars2 = axes[1].bar(x_pos + width/2, nmi_scores, width, label='NMI', color='#E91E63',
                       edgecolor='black', alpha=0.8)
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(methods)
    axes[1].set_ylabel('Score')
    axes[1].set_title('Clustering Quality: ARI & NMI')
    axes[1].legend()
    axes[1].set_ylim(0, max(max(ari_scores), max(nmi_scores)) * 1.3)

    for bar, score in zip(bars1, ari_scores):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{score:.3f}', ha='center', fontsize=8)
    for bar, score in zip(bars2, nmi_scores):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{score:.3f}', ha='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig7_integration_metrics.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig7_integration_metrics.png")
    return sil_scores, ari_scores, nmi_scores

def plot_latent_space_analysis(latent_vae, adata_rna):
    """Analyze VAE latent space structure."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    import umap as umap_lib
    reducer = umap_lib.UMAP(n_components=2, random_state=42)
    vae_umap = reducer.fit_transform(latent_vae)

    # By cell type
    cell_types = adata_rna.obs['cell_type'].values
    unique_types = np.unique(cell_types)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_types)))
    for i, ct in enumerate(unique_types):
        mask = cell_types == ct
        axes[0].scatter(vae_umap[mask, 0], vae_umap[mask, 1],
                       c=[colors[i]], label=ct, s=5, alpha=0.6)
    axes[0].set_title('VAE Latent Space by Cell Type')
    axes[0].set_xlabel('UMAP1')
    axes[0].set_ylabel('UMAP2')
    axes[0].legend(fontsize=6, markerscale=3, loc='upper right')

    # By pseudotime
    sc1 = axes[1].scatter(vae_umap[:, 0], vae_umap[:, 1],
                         c=adata_rna.obs['pseudotime'], cmap='viridis', s=5, alpha=0.7)
    plt.colorbar(sc1, ax=axes[1], shrink=0.8)
    axes[1].set_title('VAE Latent Space by Pseudotime')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')

    # Latent dimension variance
    latent_var = latent_vae.var(axis=0)
    sorted_var = np.sort(latent_var)[::-1]
    axes[2].bar(range(len(sorted_var)), sorted_var, color='steelblue', alpha=0.8)
    axes[2].set_xlabel('Latent Dimension')
    axes[2].set_ylabel('Variance')
    axes[2].set_title('VAE Latent Dimension Variance')

    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig8_latent_space.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig8_latent_space.png")

# ============================================================
# Main Pipeline
# ============================================================

def main():
    print("=" * 70)
    print("Multi-Omics Single-Cell Integration Pipeline")
    print("=" * 70)

    # 1. Generate data
    adata_rna, adata_atac, adata_meth = generate_multiomics_data()

    # 2. Preprocess
    adata_rna = preprocess_rna(adata_rna)
    adata_atac = preprocess_atac(adata_atac)
    adata_meth = preprocess_methylation(adata_meth)

    # 3. Visualize preprocessing
    print("\n=== Generating Visualizations ===")
    plot_preprocessing_qc(adata_rna, adata_atac, adata_meth)

    # 4. Anchor-based integration
    integrated_anchor, anchors_ra, anchors_rm = anchor_based_integration(
        adata_rna, adata_atac, adata_meth)

    # 5. VAE integration
    latent_vae, vae_model, train_losses, recon_losses, kl_losses = train_vae(
        adata_rna, adata_atac, adata_meth)

    # 6. Integration visualizations
    anchor_umap, vae_umap = plot_integration_comparison(adata_rna, integrated_anchor, latent_vae)
    plot_vae_training(train_losses, recon_losses, kl_losses)

    # 7. Cell lineage (velocity + pseudotime)
    adata_rna = simulate_rna_velocity(adata_rna)
    plot_pseudotime_velocity(adata_rna)

    # 8. GRN comparison
    grn_results, edges_corr, edges_mi, edges_rf = compare_grn_methods(adata_rna)
    plot_grn_comparison(grn_results)

    # 9. Immune classification
    clf_results, cm, class_names, accuracy, ari, nmi = classify_immune_subtypes(
        adata_rna, latent_vae)
    plot_immune_classification(cm, class_names, clf_results)

    # 10. Integration metrics
    sil_scores, ari_scores, nmi_scores = plot_integration_metrics(
        adata_rna, integrated_anchor, latent_vae)

    # 11. Latent space analysis
    plot_latent_space_analysis(latent_vae, adata_rna)

    # Save summary statistics
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)

    summary = {
        'n_cells': N_CELLS,
        'n_genes_rna': N_GENES_RNA,
        'n_peaks_atac': N_PEAKS_ATAC,
        'n_cpg_sites': N_CPG_SITES,
        'n_cell_types': len(CELL_TYPES),
        'anchors_rna_atac': len(anchors_ra),
        'anchors_rna_meth': len(anchors_rm),
        'vae_latent_dim': latent_vae.shape[1],
        'vae_final_loss': train_losses[-1],
        'silhouette_rna': sil_scores[0],
        'silhouette_anchor': sil_scores[1],
        'silhouette_vae': sil_scores[2],
        'ari_rna': ari_scores[0],
        'ari_anchor': ari_scores[1],
        'ari_vae': ari_scores[2],
        'nmi_rna': nmi_scores[0],
        'nmi_anchor': nmi_scores[1],
        'nmi_vae': nmi_scores[2],
        'grn_edges_corr': grn_results['Correlation']['edges'],
        'grn_edges_mi': grn_results['Mutual Information']['edges'],
        'grn_edges_rf': grn_results['Random Forest']['edges'],
        'classification_accuracy': accuracy,
        'classification_ari': ari,
        'classification_nmi': nmi,
    }

    for k, v in summary.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
        else:
            print(f"  {k}: {v}")

    # Save summary as JSON
    import json
    summary_serializable = {k: float(v) if isinstance(v, (np.floating,)) else v for k, v in summary.items()}
    with open('pipeline_summary.json', 'w') as f:
        json.dump(summary_serializable, f, indent=2)

    print(f"\nAll figures saved to {FIGURE_DIR}/")
    print("Summary saved to pipeline_summary.json")
    return summary

if __name__ == '__main__':
    summary = main()
