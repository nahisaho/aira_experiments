import json
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')
np.random.seed(42)

try:
    import scanpy as sc
except Exception:
    sc = None

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / 'data'
RESULTS_DIR = BASE_DIR / 'results'
FIG_DIR = BASE_DIR / 'figures'
LOG_PATH = BASE_DIR / 'logs' / 'process-log.jsonl'
SKILL_NAME = 'co-scientist-multi-omics'


def log_event(phase, event_type, status='ok', files_written=None, extra=None):
    entry = {
        'timestamp': datetime.utcnow().isoformat(),
        'phase': phase,
        'event_type': event_type,
        'actor': 'co-scientist',
        'skill': SKILL_NAME,
        'status': status,
        'files_written': files_written or [],
    }
    if extra:
        entry.update(extra)
    with open(LOG_PATH, 'a', encoding='utf-8') as fh:
        fh.write(json.dumps(entry, ensure_ascii=False) + '\n')


def build_similarity_graph(embedding, n_neighbors=20):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1, metric='euclidean')
    nbrs.fit(embedding)
    distances, indices = nbrs.kneighbors(embedding)
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    sigma = np.median(distances[distances > 0]) + 1e-8
    similarities = np.exp(-(distances ** 2) / (2 * sigma ** 2))

    rows = np.repeat(np.arange(embedding.shape[0]), n_neighbors)
    cols = indices.ravel()
    sim = sparse.csr_matrix((similarities.ravel(), (rows, cols)), shape=(embedding.shape[0], embedding.shape[0]))
    sim = 0.5 * (sim + sim.T)
    mean_sim = np.asarray(sim.sum(axis=1)).ravel() / max(n_neighbors, 1)
    return sim, mean_sim, distances, indices


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    log_event('integration', 'module_started')

    if sc is None:
        raise ImportError('scanpy is required for WNN integration')

    rna = sc.read_h5ad(DATA_DIR / 'rna_processed.h5ad')
    atac = sc.read_h5ad(DATA_DIR / 'atac_processed.h5ad')
    meth = sc.read_h5ad(DATA_DIR / 'methylation_processed.h5ad')

    rna_emb = np.asarray(rna.obsm['X_pca'][:, :30], dtype=np.float32)
    atac_emb = np.asarray(atac.obsm['X_lsi'][:, :30], dtype=np.float32)
    meth_emb = np.asarray(meth.obsm['X_pca'][:, :30], dtype=np.float32)

    rna_sim, rna_strength, _, _ = build_similarity_graph(rna_emb)
    atac_sim, atac_strength, _, _ = build_similarity_graph(atac_emb)
    meth_sim, meth_strength, _, _ = build_similarity_graph(meth_emb)

    weights = np.vstack([rna_strength, atac_strength, meth_strength]).T
    weights = weights / np.clip(weights.sum(axis=1, keepdims=True), 1e-8, None)

    avg_weights = 0.5 * (weights[:, None, :] + weights[None, :, :])
    joint_conn = (
        rna_sim.multiply(avg_weights[:, :, 0]) +
        atac_sim.multiply(avg_weights[:, :, 1]) +
        meth_sim.multiply(avg_weights[:, :, 2])
    ).tocsr()
    joint_conn = 0.5 * (joint_conn + joint_conn.T)
    joint_conn.setdiag(0)
    joint_conn.eliminate_zeros()

    joint_embedding = np.hstack([
        StandardScaler().fit_transform(rna_emb) * weights[:, [0]],
        StandardScaler().fit_transform(atac_emb) * weights[:, [1]],
        StandardScaler().fit_transform(meth_emb) * weights[:, [2]],
    ]).astype(np.float32)

    joint = AnnData(np.zeros((rna.n_obs, 1), dtype=np.float32), obs=rna.obs.copy())
    joint.obsm['X_joint'] = joint_embedding
    joint.obsp['connectivities'] = joint_conn
    joint.obsp['distances'] = sparse.csr_matrix(1 - joint_conn.toarray() / max(joint_conn.max(), 1e-8))
    joint.uns['neighbors'] = {
        'connectivities_key': 'connectivities',
        'distances_key': 'distances',
        'params': {'method': 'wnn', 'n_neighbors': 20},
    }
    joint.obs['rna_weight'] = weights[:, 0]
    joint.obs['atac_weight'] = weights[:, 1]
    joint.obs['methylation_weight'] = weights[:, 2]

    sc.tl.leiden(joint, resolution=0.6, key_added='wnn_cluster')
    sc.tl.umap(joint)

    metrics = pd.DataFrame([
        {'embedding': 'RNA_PCA', 'silhouette_score': float(silhouette_score(rna_emb, rna.obs['cell_type']))},
        {'embedding': 'ATAC_LSI', 'silhouette_score': float(silhouette_score(atac_emb, atac.obs['cell_type']))},
        {'embedding': 'Methylation_PCA', 'silhouette_score': float(silhouette_score(meth_emb, meth.obs['cell_type']))},
        {'embedding': 'Joint_WNN', 'silhouette_score': float(silhouette_score(joint_embedding, joint.obs['cell_type']))},
    ])
    metrics['mean_weight'] = [weights[:, 0].mean(), weights[:, 1].mean(), weights[:, 2].mean(), np.nan]
    metrics_path = RESULTS_DIR / 'integration_metrics.csv'
    metrics.to_csv(metrics_path, index=False)

    joint.write_h5ad(DATA_DIR / 'wnn_integrated.h5ad')

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.ravel()
    codes = pd.Categorical(joint.obs['cell_type']).codes
    axes[0].scatter(joint.obsm['X_umap'][:, 0], joint.obsm['X_umap'][:, 1], c=codes, cmap='tab10', s=12)
    axes[0].set_title('WNN UMAP by cell type')
    axes[0].set_xlabel('UMAP1')
    axes[0].set_ylabel('UMAP2')

    sc1 = axes[1].scatter(joint.obsm['X_umap'][:, 0], joint.obsm['X_umap'][:, 1], c=joint.obs['rna_weight'], cmap='viridis', s=12)
    axes[1].set_title('RNA modality weight')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')
    plt.colorbar(sc1, ax=axes[1], fraction=0.046)

    sc2 = axes[2].scatter(joint.obsm['X_umap'][:, 0], joint.obsm['X_umap'][:, 1], c=joint.obs['atac_weight'], cmap='cividis', s=12)
    axes[2].set_title('ATAC modality weight')
    axes[2].set_xlabel('UMAP1')
    axes[2].set_ylabel('UMAP2')
    plt.colorbar(sc2, ax=axes[2], fraction=0.046)

    sc3 = axes[3].scatter(joint.obsm['X_umap'][:, 0], joint.obsm['X_umap'][:, 1], c=joint.obs['methylation_weight'], cmap='magma', s=12)
    axes[3].set_title('Methylation modality weight')
    axes[3].set_xlabel('UMAP1')
    axes[3].set_ylabel('UMAP2')
    plt.colorbar(sc3, ax=axes[3], fraction=0.046)

    plt.tight_layout()
    fig_path = FIG_DIR / 'wnn_integration.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    files_written = [str(metrics_path), str(DATA_DIR / 'wnn_integrated.h5ad'), str(fig_path)]
    log_event('integration', 'module_completed', files_written=files_written, extra={'summary': metrics.to_dict(orient='records')})
    return {
        'joint_silhouette': float(metrics.loc[metrics['embedding'] == 'Joint_WNN', 'silhouette_score'].iloc[0]),
        'files_written': files_written,
    }


if __name__ == '__main__':
    main()
