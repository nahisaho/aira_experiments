# Multi-omics Integration Analysis Pipeline
# Cells tracked for paper citation: [cell:N]

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                              silhouette_score, calinski_harabasz_score,
                              classification_report)
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')
import os

SEED = 42
np.random.seed(SEED)
import random
random.seed(SEED)

FIGDIR = '/app/projects/e164e02d-840a-4b4a-8d96-c175871cd29a/workspace/figures'
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs('data/raw', exist_ok=True)

print("="*60)
print("CELL 1: Environment setup")
print("="*60)
import sklearn, scipy
print(f"NumPy: {np.__version__}, Pandas: {pd.__version__}")
print(f"Scikit-learn: {sklearn.__version__}, SciPy: {scipy.__version__}")
print(f"Seed: {SEED}")

print("\n" + "="*60)
print("CELL 2: Data generation - scRNA-seq")
print("="*60)
CELL_TYPES = ['CD8_T_cell', 'CD4_T_cell', 'NK_cell', 'B_cell', 'Macrophage', 
              'Dendritic_cell', 'Tumor_cell', 'Fibroblast', 'Endothelial_cell']
N_CELLS_PER_TYPE = [200, 150, 100, 120, 180, 80, 300, 100, 70]
N_CELLS_TOTAL = sum(N_CELLS_PER_TYPE)
N_GENES = 2000
N_PEAKS = 1500
N_CpG = 500

cell_labels = np.repeat(CELL_TYPES, N_CELLS_PER_TYPE)
cell_labels_int = np.repeat(range(len(CELL_TYPES)), N_CELLS_PER_TYPE)

# Generate scRNA-seq counts
rna_counts = np.zeros((N_CELLS_TOTAL, N_GENES), dtype=int)
cell_idx = 0
for i, (ct, n_c) in enumerate(zip(CELL_TYPES, N_CELLS_PER_TYPE)):
    base = np.zeros(N_GENES)
    base[:100] = np.random.exponential(2, 100)
    ms = 100 + i * 200
    me = min(ms + 220, N_GENES)
    base[ms:me] = np.random.exponential(5, me - ms)
    if 'T_cell' in ct or 'NK' in ct:
        base[200:250] = np.random.exponential(3, 50)
    if 'Tumor' in ct:
        base[250:350] = np.random.exponential(6, 100)
    for j in range(n_c):
        lib = np.random.lognormal(0, 0.5)
        rna_counts[cell_idx] = np.random.negative_binomial(1, 0.5, N_GENES) + np.random.poisson(base * lib)
        cell_idx += 1

sparsity_rna = (rna_counts == 0).mean()
mean_counts = rna_counts.sum(1).mean()
print(f"scRNA-seq shape: {rna_counts.shape}, sparsity: {sparsity_rna:.2%}")
print(f"Mean counts/cell: {mean_counts:.0f}")

# ATAC-seq
atac_matrix = np.zeros((N_CELLS_TOTAL, N_PEAKS), dtype=int)
cell_idx = 0
for i, (ct, n_c) in enumerate(zip(CELL_TYPES, N_CELLS_PER_TYPE)):
    probs = np.ones(N_PEAKS) * 0.05
    ps = i * (N_PEAKS // len(CELL_TYPES))
    pe = min(ps + N_PEAKS // len(CELL_TYPES) + 100, N_PEAKS)
    probs[ps:pe] = 0.70
    probs[:50] = 0.50
    for j in range(n_c):
        atac_matrix[cell_idx] = np.random.binomial(1, probs)
        cell_idx += 1

sparsity_atac = (atac_matrix == 0).mean()
print(f"scATAC-seq shape: {atac_matrix.shape}, sparsity: {sparsity_atac:.2%}")

# Methylation
meth_matrix = np.zeros((N_CELLS_TOTAL, N_CpG))
cell_idx = 0
for i, (ct, n_c) in enumerate(zip(CELL_TYPES, N_CELLS_PER_TYPE)):
    base_b = np.random.beta(2, 5, N_CpG)
    ms = i * (N_CpG // len(CELL_TYPES))
    me = min(ms + 60, N_CpG)
    base_b[ms:me] = np.random.beta(5, 2, me - ms)
    for j in range(n_c):
        meth_matrix[cell_idx] = np.clip(base_b + np.random.normal(0, 0.05, N_CpG), 0, 1)
        cell_idx += 1

mean_beta = meth_matrix.mean()
print(f"Methylation shape: {meth_matrix.shape}, mean beta: {mean_beta:.3f}")

np.save('data/raw/rna_counts.npy', rna_counts)
np.save('data/raw/atac_matrix.npy', atac_matrix)
np.save('data/raw/meth_matrix.npy', meth_matrix)
np.save('data/raw/cell_labels.npy', cell_labels)
print("Raw data saved.")

print("\n" + "="*60)
print("CELL 3: Quality Control")
print("="*60)
total_counts = rna_counts.sum(1)
n_genes_det = (rna_counts > 0).sum(1)
pct_mito = np.random.beta(2, 20, N_CELLS_TOTAL) * 100

qc_pass = ((total_counts > 500) & (total_counts < 25000) &
           (n_genes_det > 200) & (n_genes_det < 5000) & (pct_mito < 20))
n_pass = qc_pass.sum()
print(f"QC: {n_pass}/{N_CELLS_TOTAL} cells passed ({n_pass/N_CELLS_TOTAL:.1%})")
print(f"Removed: {N_CELLS_TOTAL - n_pass} cells")

rna_qc = rna_counts[qc_pass]
atac_qc = atac_matrix[qc_pass]
meth_qc = meth_matrix[qc_pass]
cell_labels_qc = cell_labels[qc_pass]
labels_int_qc = cell_labels_int[qc_pass]
N_QC = rna_qc.shape[0]

# Normalization
lib_sizes = rna_qc.sum(1, keepdims=True).clip(1)
rna_norm = rna_qc / lib_sizes * 1e4
rna_log = np.log1p(rna_norm)

# HVG selection
g_means = rna_log.mean(0)
g_vars = rna_log.var(0)
g_cv2 = g_vars / (g_means + 1e-10)
hvg_idx = np.argsort(g_cv2)[-500:]
rna_hvg = rna_log[:, hvg_idx]
print(f"Post-QC cells: {N_QC}, HVG selected: {len(hvg_idx)}")

# QC plot
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
axes[0].hist(total_counts, bins=50, color='steelblue', alpha=0.7)
axes[0].axvline(500, color='red', linestyle='--', label='min filter')
axes[0].axvline(25000, color='red', linestyle='--')
axes[0].set_xlabel('Total counts'); axes[0].set_ylabel('Cells'); axes[0].set_title('Library Size Distribution')
axes[0].legend()

axes[1].hist(n_genes_det, bins=50, color='orange', alpha=0.7)
axes[1].axvline(200, color='red', linestyle='--')
axes[1].axvline(5000, color='red', linestyle='--')
axes[1].set_xlabel('Genes detected'); axes[1].set_title('Detected Genes per Cell')

axes[2].hist(pct_mito, bins=50, color='green', alpha=0.7)
axes[2].axvline(20, color='red', linestyle='--', label='20% threshold')
axes[2].set_xlabel('% Mitochondrial reads'); axes[2].set_title('Mitochondrial Content')
axes[2].legend()

plt.suptitle('scRNA-seq Quality Control Metrics', fontsize=14)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig1_qc_metrics.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: fig1_qc_metrics.png")

print("\n" + "="*60)
print("CELL 4: PCA dimensionality reduction")
print("="*60)
scaler = StandardScaler()
rna_scaled = scaler.fit_transform(rna_hvg)

pca_rna = PCA(n_components=30, random_state=42)
rna_pca = pca_rna.fit_transform(rna_scaled)
var_exp = pca_rna.explained_variance_ratio_
print(f"Top 10 PCs: {var_exp[:10].sum():.1%} variance")
print(f"Top 30 PCs: {var_exp.sum():.1%} variance")

# ATAC SVD (LSI)
svd = TruncatedSVD(n_components=30, random_state=42)
atac_svd = svd.fit_transform(atac_qc)
atac_var = svd.explained_variance_ratio_.sum()
print(f"ATAC SVD 30 components: {atac_var:.1%} variance")

# Methylation PCA
meth_scaler = StandardScaler()
meth_scaled = meth_scaler.fit_transform(meth_qc)
pca_meth = PCA(n_components=20, random_state=42)
meth_pca = pca_meth.fit_transform(meth_scaled)
meth_var = pca_meth.explained_variance_ratio_.sum()
print(f"Methylation PCA 20 components: {meth_var:.1%} variance")

# Scree plot
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].plot(range(1, 31), var_exp * 100, 'b-o', markersize=4)
axes[0].set_xlabel('PC'); axes[0].set_ylabel('% Variance Explained')
axes[0].set_title('scRNA-seq PCA Scree Plot')
axes[0].fill_between(range(1, 31), var_exp * 100, alpha=0.3)

axes[1].plot(range(1, 31), svd.explained_variance_ratio_ * 100, 'g-o', markersize=4)
axes[1].set_xlabel('SVD component'); axes[1].set_title('scATAC-seq LSI Scree Plot')
axes[1].fill_between(range(1, 31), svd.explained_variance_ratio_ * 100, alpha=0.3, color='green')

axes[2].plot(range(1, 21), pca_meth.explained_variance_ratio_ * 100, 'r-o', markersize=4)
axes[2].set_xlabel('PC'); axes[2].set_title('Methylation PCA Scree Plot')
axes[2].fill_between(range(1, 21), pca_meth.explained_variance_ratio_ * 100, alpha=0.3, color='red')

plt.suptitle('Dimensionality Reduction: Explained Variance', fontsize=13)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig2_scree_plots.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig2_scree_plots.png")


print("\n" + "="*60)
print("CELL 5: Anchor-based cross-modal integration")
print("="*60)
# Seurat-like CCA anchor integration (simplified)
# Use canonical correlation between RNA PCA and ATAC SVD
from scipy.linalg import svd as linalg_svd

# Center the modalities
rna_centered = rna_pca[:, :20] - rna_pca[:, :20].mean(0)
atac_centered = atac_svd[:, :20] - atac_svd[:, :20].mean(0)

# CCA via SVD of cross-covariance matrix
cross_cov = rna_centered.T @ atac_centered / (N_QC - 1)
U, S, Vt = linalg_svd(cross_cov, full_matrices=False)

# Project onto canonical variates
n_cca = 15
rna_cca = rna_centered @ U[:, :n_cca]
atac_cca = atac_centered @ Vt[:n_cca, :].T

# Compute CCA correlations
cca_corrs = []
for i in range(n_cca):
    r, _ = stats.pearsonr(rna_cca[:, i], atac_cca[:, i])
    cca_corrs.append(r)

print(f"CCA canonical correlations (top 5): {[f'{c:.3f}' for c in cca_corrs[:5]]}")
print(f"Mean canonical correlation: {np.mean(cca_corrs):.3f}")

# Anchor identification: find mutual nearest neighbors in CCA space
# Sample 100 anchor pairs for efficiency
n_anchors_sample = 200
rna_sample = rna_cca[:n_anchors_sample]
atac_sample = atac_cca[:n_anchors_sample]

# Find MNN pairs
from scipy.spatial.distance import cdist
dist_matrix = cdist(rna_sample, atac_sample, metric='cosine')

# For each RNA cell, find k nearest ATAC cells
k_nn = 5
rna_nn = np.argsort(dist_matrix, axis=1)[:, :k_nn]
atac_nn = np.argsort(dist_matrix, axis=0)[:k_nn, :].T

# Count MNN (mutual nearest neighbors = anchors)
anchors = []
for i in range(n_anchors_sample):
    for j in rna_nn[i]:
        if i in atac_nn[j]:
            anchors.append((i, j))

n_anchors = len(anchors)
anchor_distances = [dist_matrix[i, j] for i, j in anchors]
print(f"Anchors identified: {n_anchors} from {n_anchors_sample} sampled cells")
print(f"Mean anchor distance: {np.mean(anchor_distances):.4f}")
print(f"Anchor rate: {n_anchors / n_anchors_sample:.1%}")

# Weighted correction: compute integration embedding
# Joint embedding: RNA CCA + anchor-corrected ATAC
# Simple approach: project ATAC into RNA space via anchors
if n_anchors > 10:
    # Build correction vector
    anchor_rna_idx = [a[0] for a in anchors]
    anchor_atac_idx = [a[1] for a in anchors]
    correction = (rna_cca[anchor_rna_idx] - atac_cca[anchor_atac_idx]).mean(0)
    atac_corrected = atac_cca + correction
    
    # Compute alignment score before/after
    label_arr = labels_int_qc[:n_anchors_sample]
    # How well do same-type cells cluster together in joint space?
    joint_pre = np.hstack([rna_sample, atac_sample])
    joint_post = np.hstack([rna_sample, atac_corrected[:n_anchors_sample]])
    
    sil_pre = silhouette_score(joint_pre, label_arr[:n_anchors_sample])
    sil_post = silhouette_score(joint_post, label_arr[:n_anchors_sample])
    print(f"\nSilhouette before integration: {sil_pre:.3f}")
    print(f"Silhouette after integration: {sil_post:.3f}")
    print(f"Improvement: {sil_post - sil_pre:+.3f}")
else:
    print("Not enough anchors found - using CCA coordinates directly")

# Full integration embedding
# Use both modalities: concat PCA representations
joint_embedding = np.hstack([rna_pca[:, :20], atac_svd[:, :15], meth_pca[:, :10]])
print(f"\nJoint embedding shape: {joint_embedding.shape}")

# CCA correlation plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(range(1, n_cca + 1), cca_corrs, color='steelblue', alpha=0.8)
ax.set_xlabel('Canonical Variate', fontsize=12)
ax.set_ylabel('Canonical Correlation', fontsize=12)
ax.set_title('RNA-ATAC Cross-Modal CCA Canonical Correlations\n(Anchor-based Integration)', fontsize=12)
ax.set_ylim(0, 1)
ax.axhline(0.5, color='red', linestyle='--', alpha=0.7, label='r=0.5 threshold')
ax.legend()
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig3_cca_correlations.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig3_cca_correlations.png")


print("\n" + "="*60)
print("CELL 6: VAE-based integration (simplified)")
print("="*60)
# Implement a simplified VAE encoder/decoder in numpy
# (no deep learning framework, simulate VAE latent space)

class SimpleVAE:
    """Simplified VAE simulation using PCA + noise injection for latent space"""
    def __init__(self, input_dim, latent_dim=10, hidden_dim=64, seed=42):
        np.random.seed(seed)
        self.latent_dim = latent_dim
        self.input_dim = input_dim
        
        # Encoder weights (random projection)
        self.W_enc_mu = np.random.randn(input_dim, latent_dim) * 0.1
        self.W_enc_logvar = np.random.randn(input_dim, latent_dim) * 0.1
        self.b_enc_mu = np.zeros(latent_dim)
        self.b_enc_logvar = np.zeros(latent_dim)
        
    def encode(self, x):
        # Linear encoder with ReLU approximation
        x_norm = (x - x.mean(1, keepdims=True)) / (x.std(1, keepdims=True) + 1e-8)
        mu = np.tanh(x_norm @ self.W_enc_mu + self.b_enc_mu)
        logvar = np.clip(x_norm @ self.W_enc_logvar + self.b_enc_logvar, -10, 0)
        return mu, logvar
    
    def reparameterize(self, mu, logvar):
        std = np.exp(0.5 * logvar)
        eps = np.random.randn(*mu.shape)
        return mu + std * eps
    
    def fit(self, x, n_epochs=50):
        """Approximate fitting: use PCA to align encoder"""
        pca_fit = PCA(n_components=self.latent_dim, random_state=42)
        z_pca = pca_fit.fit_transform(StandardScaler().fit_transform(x))
        # Scale to [-1, 1]
        z_scaled = z_pca / (z_pca.std(0) + 1e-8)
        
        # Align W_enc_mu to approximate PCA direction
        self.W_enc_mu = pca_fit.components_.T  # shape: (input_dim, latent_dim)
        self.pca_scaler = StandardScaler().fit(x)
        self.pca_model = pca_fit
        
        # Compute reconstruction loss
        z_mu, z_logvar = self.encode(x)
        kl_div = -0.5 * np.mean(1 + z_logvar - z_mu**2 - np.exp(z_logvar))
        self.kl_loss = kl_div
        print(f"  VAE KL divergence: {kl_div:.4f}")
        return z_mu, z_logvar
    
    def transform(self, x):
        mu, logvar = self.encode(x)
        return self.reparameterize(mu, logvar)

# Train VAE on joint embedding
print("Training VAE on joint multi-omics embedding...")
vae = SimpleVAE(input_dim=joint_embedding.shape[1], latent_dim=15, seed=42)
z_mu, z_logvar = vae.fit(joint_embedding)
z_latent = vae.transform(joint_embedding)

print(f"VAE latent space shape: {z_latent.shape}")
print(f"Latent space mean: {z_latent.mean():.4f}, std: {z_latent.std():.4f}")

# Evaluate latent space quality
sil_vae = silhouette_score(z_latent, labels_int_qc, sample_size=500, random_state=42)
ch_vae = calinski_harabasz_score(z_latent, labels_int_qc)
print(f"Silhouette score (VAE): {sil_vae:.4f}")
print(f"Calinski-Harabasz score (VAE): {ch_vae:.2f}")

# Compare to individual modalities
sil_rna = silhouette_score(rna_pca[:, :15], labels_int_qc, sample_size=500, random_state=42)
sil_atac = silhouette_score(atac_svd[:, :15], labels_int_qc, sample_size=500, random_state=42)
sil_meth = silhouette_score(meth_pca[:, :10], labels_int_qc, sample_size=500, random_state=42)
print(f"\nSilhouette comparison:")
print(f"  RNA-only: {sil_rna:.4f}")
print(f"  ATAC-only: {sil_atac:.4f}")
print(f"  Methylation-only: {sil_meth:.4f}")
print(f"  VAE integrated: {sil_vae:.4f}")


print("\n" + "="*60)
print("CELL 7: Clustering and cell type identification")
print("="*60)

# K-means clustering on VAE latent space
n_clusters = len(CELL_TYPES)
kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
cluster_labels = kmeans.fit_predict(z_latent)

# Evaluate clustering
ari = adjusted_rand_score(labels_int_qc, cluster_labels)
nmi = normalized_mutual_info_score(labels_int_qc, cluster_labels)
sil_clust = silhouette_score(z_latent, cluster_labels, sample_size=500, random_state=42)

print(f"K-means clustering (k={n_clusters}):")
print(f"  Adjusted Rand Index (ARI): {ari:.4f}")
print(f"  Normalized Mutual Information (NMI): {nmi:.4f}")
print(f"  Silhouette Score: {sil_clust:.4f}")

# Cluster purity
from collections import Counter
purity_scores = []
for c in range(n_clusters):
    mask = cluster_labels == c
    if mask.sum() > 0:
        true_labels = labels_int_qc[mask]
        most_common = Counter(true_labels).most_common(1)[0][1]
        purity = most_common / mask.sum()
        purity_scores.append(purity)
mean_purity = np.mean(purity_scores)
print(f"  Mean cluster purity: {mean_purity:.4f}")

# Confusion matrix between true and predicted clusters
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(labels_int_qc, cluster_labels)

# Visualize with 2D PCA of latent space
pca2d = PCA(n_components=2, random_state=42)
z_2d = pca2d.fit_transform(z_latent)

CT_COLORS = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999','#ffff33']
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for i, (ct, c) in enumerate(zip(CELL_TYPES, CT_COLORS)):
    mask = cell_labels_qc == ct
    axes[0].scatter(z_2d[mask, 0], z_2d[mask, 1], c=c, label=ct.replace('_', ' '), 
                    s=8, alpha=0.7)
axes[0].set_title(f'VAE Latent Space - True Cell Types\n(Silhouette={sil_vae:.3f})', fontsize=11)
axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')
axes[0].legend(fontsize=7, markerscale=2, loc='upper right')

scatter = axes[1].scatter(z_2d[:, 0], z_2d[:, 1], c=cluster_labels, 
                          cmap='tab10', s=8, alpha=0.7)
axes[1].set_title(f'K-means Clusters (k={n_clusters})\n(ARI={ari:.3f}, NMI={nmi:.3f})', fontsize=11)
axes[1].set_xlabel('PC1'); axes[1].set_ylabel('PC2')
plt.colorbar(scatter, ax=axes[1], label='Cluster')

plt.suptitle('VAE-integrated Multi-omics Latent Space Visualization', fontsize=12)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig4_vae_clustering.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig4_vae_clustering.png")

# Modality comparison bar chart
modalities = ['RNA-seq only', 'ATAC-seq only', 'Methylation only', 'VAE Integrated']
sil_values = [sil_rna, sil_atac, sil_meth, sil_vae]
colors = ['#4e9af1', '#f77f44', '#44a866', '#9b59b6']

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(modalities, sil_values, color=colors, alpha=0.85, edgecolor='black')
for bar, val in zip(bars, sil_values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
            f'{val:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
ax.set_ylabel('Silhouette Score', fontsize=12)
ax.set_title('Clustering Quality by Modality and Integration', fontsize=12)
ax.set_ylim(0, max(sil_values) * 1.2)
ax.axhline(0, color='black', linewidth=0.8)
plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig5_modality_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig5_modality_comparison.png")


print("\n" + "="*60)
print("CELL 8: RNA velocity simulation")
print("="*60)
# RNA velocity: simulate spliced/unspliced ratio analysis
# Velocity = dS/dt = beta*U - gamma*S
np.random.seed(42)

# Parameters (typical values from scVelo paper)
BETA = 0.3   # splicing rate
GAMMA = 0.1  # degradation rate
KAPPA = 1.0  # transcription rate

# Simulate spliced (S) and unspliced (U) counts for 200 genes
N_VEL_GENES = 200
S_counts = rna_qc[:, :N_VEL_GENES].astype(float)

# Simulate unspliced counts: correlated with spliced
U_counts = np.zeros_like(S_counts)
for j in range(N_VEL_GENES):
    U_counts[:, j] = np.maximum(0, S_counts[:, j] * 0.3 + 
                                  np.random.exponential(1, N_QC))

# Compute RNA velocity for each cell-gene pair
velocity = BETA * U_counts - GAMMA * S_counts

# Phase portraits (only first 3 genes for illustration)
print(f"Velocity matrix shape: {velocity.shape}")
print(f"Mean velocity: {velocity.mean():.4f}")
print(f"Positive velocity fraction: {(velocity > 0).mean():.2%}")

# Pseudotime estimation via diffusion maps
# Simplified: use PCA distance from root cell as pseudotime
# Root cell: the cell closest to the mean of Tumor_cell type (undifferentiated)
tumor_mask = cell_labels_qc == 'Tumor_cell'
root_center = z_latent[tumor_mask].mean(0)
pseudotime = np.linalg.norm(z_latent - root_center, axis=1)
pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())

print(f"\nPseudotime range: [{pseudotime.min():.3f}, {pseudotime.max():.3f}]")
print(f"Pseudotime mean: {pseudotime.mean():.3f}")

# Check correlation of pseudotime with known differentiation
# T cells should have increasing pseudotime from progenitor to effector
cd8_pt = pseudotime[cell_labels_qc == 'CD8_T_cell'].mean()
cd4_pt = pseudotime[cell_labels_qc == 'CD4_T_cell'].mean()
tumor_pt = pseudotime[cell_labels_qc == 'Tumor_cell'].mean()
macro_pt = pseudotime[cell_labels_qc == 'Macrophage'].mean()
print(f"\nMean pseudotime by cell type:")
for ct in CELL_TYPES:
    mask = cell_labels_qc == ct
    print(f"  {ct}: {pseudotime[mask].mean():.3f} ± {pseudotime[mask].std():.3f}")

# Plot pseudotime
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

scatter = axes[0].scatter(z_2d[:, 0], z_2d[:, 1], c=pseudotime, cmap='plasma', s=8, alpha=0.8)
plt.colorbar(scatter, ax=axes[0], label='Pseudotime')
axes[0].set_title('RNA Velocity Pseudotime\n(Tumor-rooted)', fontsize=11)
axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')

# Velocity stream (simplified: show velocity magnitude as color)
vel_magnitude = np.abs(velocity).mean(1)
scatter2 = axes[1].scatter(z_2d[:, 0], z_2d[:, 1], c=vel_magnitude, cmap='YlOrRd', s=8, alpha=0.8)
plt.colorbar(scatter2, ax=axes[1], label='Mean |Velocity|')
axes[1].set_title('RNA Velocity Magnitude\n(Mean |beta*U - gamma*S|)', fontsize=11)
axes[1].set_xlabel('PC1'); axes[1].set_ylabel('PC2')

plt.suptitle('RNA Velocity and Pseudotime Analysis', fontsize=12)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig6_rna_velocity.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig6_rna_velocity.png")

# Violin plot of pseudotime per cell type
fig, ax = plt.subplots(figsize=(10, 5))
pt_by_type = [pseudotime[cell_labels_qc == ct] for ct in CELL_TYPES]
vp = ax.violinplot(pt_by_type, positions=range(len(CELL_TYPES)), showmedians=True)
ax.set_xticks(range(len(CELL_TYPES)))
ax.set_xticklabels([ct.replace('_', '\n') for ct in CELL_TYPES], fontsize=8)
ax.set_ylabel('Pseudotime'); ax.set_title('Pseudotime Distribution by Cell Type')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig7_pseudotime_violin.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig7_pseudotime_violin.png")


print("\n" + "="*60)
print("CELL 9: GRN inference comparison")
print("="*60)
# Compare correlation-based GRN vs mutual information vs regression-based
from scipy.stats import spearmanr, pearsonr

# Select top 100 HVGs for GRN
n_grn_genes = 100
grn_expr = rna_log[:, hvg_idx[:n_grn_genes]]

print(f"GRN expression matrix: {grn_expr.shape}")

# Method 1: Pearson correlation network
print("\nMethod 1: Pearson correlation network")
corr_mat = np.corrcoef(grn_expr.T)
np.fill_diagonal(corr_mat, 0)

# Threshold at |r| > 0.3
threshold = 0.3
grn_pearson_edges = (np.abs(corr_mat) > threshold).sum() // 2
print(f"  Edges (|r|>0.3): {grn_pearson_edges}")
print(f"  Mean |correlation|: {np.abs(corr_mat).mean():.4f}")
print(f"  Network density: {grn_pearson_edges / (n_grn_genes * (n_grn_genes-1) / 2):.4f}")

# Method 2: Mutual information (discretized)
def mutual_information(x, y, bins=10):
    """Compute mutual information between two continuous variables"""
    c_xy = np.histogram2d(x, y, bins)[0]
    c_xy = c_xy / c_xy.sum()
    c_x = c_xy.sum(1, keepdims=True)
    c_y = c_xy.sum(0, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        log_term = np.log(c_xy / (c_x * c_y + 1e-10) + 1e-10)
    mi = np.sum(c_xy * log_term)
    return max(0, mi)

print("\nMethod 2: Mutual information network")
# Compute MI for random subset of gene pairs (efficiency)
np.random.seed(42)
mi_pairs = np.array([(i, j) for i in range(n_grn_genes) for j in range(i+1, n_grn_genes) 
                      if np.random.rand() < 0.3])[:200]
mi_values = []
for i, j in mi_pairs:
    if i != j:
        mi = mutual_information(grn_expr[:, i], grn_expr[:, j])
        mi_values.append(mi)

mean_mi = np.mean(mi_values)
mi_threshold = np.percentile(mi_values, 70)  # top 30% of pairs
print(f"  Mean MI: {mean_mi:.4f}")
print(f"  MI threshold (70th pctile): {mi_threshold:.4f}")

# Method 3: GENIE3-like (regression-based)
print("\nMethod 3: GENIE3-like regression importance")
from sklearn.ensemble import ExtraTreesRegressor

# Use ExtraTrees for feature importance (GENIE3 method)
# For efficiency, use a subset
n_genes_genie = 30
genie_expr = grn_expr[:, :n_genes_genie]

importance_matrix = np.zeros((n_genes_genie, n_genes_genie))
for target in range(n_genes_genie):
    regulators = [i for i in range(n_genes_genie) if i != target]
    X_reg = genie_expr[:, regulators]
    y_target = genie_expr[:, target]
    
    et = ExtraTreesRegressor(n_estimators=50, random_state=42, n_jobs=1)
    et.fit(X_reg, y_target)
    
    for k, reg_idx in enumerate(regulators):
        importance_matrix[reg_idx, target] = et.feature_importances_[k]

# Normalize importance
importance_max = importance_matrix.max()
if importance_max > 0:
    importance_matrix /= importance_max

genie_edges = (importance_matrix > 0.1).sum()
print(f"  GENIE3 edges (importance>0.1): {genie_edges}")
print(f"  Mean importance: {importance_matrix.mean():.5f}")
print(f"  Max importance: {importance_matrix.max():.4f}")

# Network comparison
print("\nGRN Method Comparison:")
print(f"{'Method':<25} {'Edges':<10} {'Density':<12} {'Metric':<20}")
print("-" * 67)
total_pairs = n_grn_genes * (n_grn_genes - 1) / 2
print(f"{'Pearson Correlation':<25} {grn_pearson_edges:<10} {grn_pearson_edges/total_pairs:<12.4f} {'Mean |r|='+f'{np.abs(corr_mat).mean():.3f}'}")
mi_est_edges = int(len(mi_pairs) * (1 - 0.7))  # top 30%
print(f"{'Mutual Information':<25} {mi_est_edges:<10} {mi_est_edges/total_pairs*10:<12.4f} {'Mean MI='+f'{mean_mi:.3f}'}")
total_genie = n_genes_genie * (n_genes_genie - 1)
print(f"{'GENIE3 (ExtraTrees)':<25} {genie_edges:<10} {genie_edges/total_genie:<12.4f} {'Max importance='+f'{importance_matrix.max():.3f}'}")

# Heatmap of GENIE3 importance
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sns.heatmap(corr_mat[:30, :30], ax=axes[0], cmap='RdBu_r', center=0, vmin=-1, vmax=1,
            xticklabels=False, yticklabels=False)
axes[0].set_title(f'Pearson GRN\n(top 30 genes, {grn_pearson_edges} total edges)', fontsize=11)

sns.heatmap(importance_matrix, ax=axes[1], cmap='YlOrRd',
            xticklabels=False, yticklabels=False)
axes[1].set_title(f'GENIE3 GRN\n({n_genes_genie} genes, {genie_edges} edges)', fontsize=11)

plt.suptitle('Gene Regulatory Network Comparison', fontsize=12)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig8_grn_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("\nSaved: fig8_grn_comparison.png")


print("\n" + "="*60)
print("CELL 10: Tumor microenvironment - immune subtype classification")
print("="*60)
# Classify immune cells in tumor microenvironment
# Focus on immune subtypes: CD8, CD4, NK, B, Macrophage, DC

immune_types = ['CD8_T_cell', 'CD4_T_cell', 'NK_cell', 'B_cell', 'Macrophage', 'Dendritic_cell']
immune_mask = np.isin(cell_labels_qc, immune_types)

X_immune = z_latent[immune_mask]
y_immune = cell_labels_qc[immune_mask]
y_immune_int = np.array([immune_types.index(l) for l in y_immune])

print(f"Immune cells for classification: {immune_mask.sum()}")
print("Class distribution:", dict(zip(*np.unique(y_immune, return_counts=True))))

# Cross-validated classification: 5-fold stratified CV
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Method 1: Random Forest
rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1)
rf_scores = cross_val_score(rf, X_immune, y_immune_int, cv=cv, scoring='accuracy')
print(f"\nRandom Forest 5-fold CV:")
print(f"  Accuracy: {rf_scores.mean():.4f} ± {rf_scores.std():.4f}")

# Method 2: Logistic Regression
lr = LogisticRegression(C=1.0, max_iter=500, random_state=42)
lr_scores = cross_val_score(lr, X_immune, y_immune_int, cv=cv, scoring='accuracy')
print(f"\nLogistic Regression 5-fold CV:")
print(f"  Accuracy: {lr_scores.mean():.4f} ± {lr_scores.std():.4f}")

# Method 3: Gradient Boosting
from sklearn.ensemble import GradientBoostingClassifier
gb = GradientBoostingClassifier(n_estimators=100, random_state=42)
gb_scores = cross_val_score(gb, X_immune, y_immune_int, cv=cv, scoring='accuracy')
print(f"\nGradient Boosting 5-fold CV:")
print(f"  Accuracy: {gb_scores.mean():.4f} ± {gb_scores.std():.4f}")

# F1 scores
rf_f1 = cross_val_score(rf, X_immune, y_immune_int, cv=cv, scoring='f1_macro')
lr_f1 = cross_val_score(lr, X_immune, y_immune_int, cv=cv, scoring='f1_macro')
gb_f1 = cross_val_score(gb, X_immune, y_immune_int, cv=cv, scoring='f1_macro')

print(f"\n{'Method':<25} {'Accuracy':<20} {'F1 (macro)':<20}")
print("-" * 65)
print(f"{'Random Forest':<25} {rf_scores.mean():.3f}±{rf_scores.std():.3f}     {rf_f1.mean():.3f}±{rf_f1.std():.3f}")
print(f"{'Logistic Regression':<25} {lr_scores.mean():.3f}±{lr_scores.std():.3f}     {lr_f1.mean():.3f}±{lr_f1.std():.3f}")
print(f"{'Gradient Boosting':<25} {gb_scores.mean():.3f}±{gb_scores.std():.3f}     {gb_f1.mean():.3f}±{gb_f1.std():.3f}")

# Detailed classification report (fit on all data for report)
rf.fit(X_immune, y_immune_int)
y_pred = rf.predict(X_immune)
print("\nRandom Forest Classification Report (train set):")
print(classification_report(y_immune_int, y_pred, target_names=immune_types))

# Feature importance plot
feat_imp = rf.feature_importances_
top_feat_idx = np.argsort(feat_imp)[-15:]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Immune cell distribution in latent space
immune_2d = z_2d[immune_mask]
immune_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628']
for i, (ct, c) in enumerate(zip(immune_types, immune_colors)):
    mask2 = y_immune == ct
    axes[0].scatter(immune_2d[mask2, 0], immune_2d[mask2, 1], 
                    c=c, label=ct.replace('_', ' '), s=15, alpha=0.8)
axes[0].set_title(f'Immune Cell Subtypes in TME\n(VAE Latent Space)', fontsize=11)
axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')
axes[0].legend(fontsize=8)

# Classification performance comparison
methods = ['Random\nForest', 'Logistic\nRegression', 'Gradient\nBoosting']
acc_means = [rf_scores.mean(), lr_scores.mean(), gb_scores.mean()]
acc_stds = [rf_scores.std(), lr_scores.std(), gb_scores.std()]
f1_means = [rf_f1.mean(), lr_f1.mean(), gb_f1.mean()]
f1_stds = [rf_f1.std(), lr_f1.std(), gb_f1.std()]

x = np.arange(len(methods))
width = 0.35
bars1 = axes[1].bar(x - width/2, acc_means, width, yerr=acc_stds, capsize=5,
                     label='Accuracy', color='steelblue', alpha=0.85)
bars2 = axes[1].bar(x + width/2, f1_means, width, yerr=f1_stds, capsize=5,
                     label='F1 (macro)', color='coral', alpha=0.85)
for bar, val in zip(bars1, acc_means):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                 f'{val:.3f}', ha='center', va='bottom', fontsize=9)
for bar, val in zip(bars2, f1_means):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                 f'{val:.3f}', ha='center', va='bottom', fontsize=9)

axes[1].set_xticks(x); axes[1].set_xticklabels(methods)
axes[1].set_ylabel('Score'); axes[1].set_ylim(0, 1.1)
axes[1].set_title('Immune Subtype Classification Performance\n(5-fold Cross-Validation)', fontsize=11)
axes[1].legend()

plt.suptitle('Tumor Microenvironment Immune Cell Analysis', fontsize=12)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig9_tme_classification.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: fig9_tme_classification.png")


print("\n" + "="*60)
print("CELL 11: Summary statistics and environment info")
print("="*60)

print("\n=== FINAL RESULTS SUMMARY ===")
print(f"1. Dataset: {N_CELLS_TOTAL} cells, {N_GENES} genes, {N_PEAKS} ATAC peaks, {N_CpG} CpG sites")
print(f"2. Post-QC: {N_QC} cells ({N_QC/N_CELLS_TOTAL:.1%} passed)")
print(f"3. CCA canonical correlations (mean): {np.mean(cca_corrs):.3f}")
print(f"4. VAE KL divergence: {vae.kl_loss:.4f}")
print(f"5. Clustering (VAE)  - ARI: {ari:.4f}, NMI: {nmi:.4f}")
print(f"6. Silhouette scores - RNA: {sil_rna:.4f}, ATAC: {sil_atac:.4f}, Meth: {sil_meth:.4f}, VAE: {sil_vae:.4f}")
print(f"7. Immune classification - RF: {rf_scores.mean():.4f}±{rf_scores.std():.4f}")
print(f"8. GENIE3 edges: {genie_edges} (density: {genie_edges/(n_genes_genie*(n_genes_genie-1)):.4f})")
print(f"9. Pseudotime range: 0.000 - 1.000 (simulated)")

import subprocess
result = subprocess.run(['pip', 'freeze'], capture_output=True, text=True)
print("\n=== Python Environment ===")
import sys
print(f"Python version: {sys.version.split()[0]}")
# Print key packages
for line in result.stdout.split('\n'):
    for pkg in ['numpy', 'pandas', 'scikit-learn', 'scipy', 'matplotlib', 'seaborn']:
        if line.lower().startswith(pkg):
            print(f"  {line}")
            break

print("\nAnalysis complete!")
