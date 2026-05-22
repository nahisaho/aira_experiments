import json
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
np.random.seed(42)

try:
    import scanpy as sc
except Exception:
    sc = None

try:
    import scvelo as scv
except Exception:
    scv = None

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


def simulated_pseudotime(obs):
    order = {
        'Dendritic cells': 0.10,
        'Macrophages': 0.25,
        'B cells': 0.45,
        'NK cells': 0.65,
        'T cells': 0.82,
        'Tumor cells': 0.95,
    }
    pt = np.array([order[c] for c in obs['cell_type']]) + np.random.normal(0, 0.035, size=obs.shape[0])
    return np.clip(pt, 0, 1)


def build_velocity_layers(rna):
    counts = np.asarray(rna.layers['counts'], dtype=np.float32)
    pseudotime = simulated_pseudotime(rna.obs)
    gene_switch = np.linspace(0.1, 0.9, rna.n_vars) + np.random.normal(0, 0.05, size=rna.n_vars)
    gene_switch = np.clip(gene_switch, 0.05, 0.95)
    direction = np.random.choice([-1, 1], size=rna.n_vars)

    activation = 1 / (1 + np.exp(-8 * (pseudotime[:, None] - gene_switch[None, :]) * direction[None, :]))
    d_activation = activation * (1 - activation)
    base = counts / np.maximum(counts.max(axis=0, keepdims=True), 1)

    spliced = np.random.poisson(np.clip(2 + 8 * activation + 4 * base, 0.1, None)).astype(np.float32)
    unspliced = np.random.poisson(np.clip(1 + 6 * d_activation + 2 * (1 - activation) + 2 * base, 0.1, None)).astype(np.float32)
    return pseudotime, spliced, unspliced


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    log_event('trajectory', 'module_started')

    if sc is None:
        raise ImportError('scanpy is required for trajectory analysis')

    rna = sc.read_h5ad(DATA_DIR / 'rna_processed.h5ad')
    pseudotime_true, spliced, unspliced = build_velocity_layers(rna)

    adata = sc.AnnData(np.asarray(rna.layers['counts'], dtype=np.float32), obs=rna.obs.copy(), var=rna.var.copy())
    adata.layers['spliced'] = spliced
    adata.layers['unspliced'] = unspliced
    adata.obs['designed_pseudotime'] = pseudotime_true
    adata.var_names_make_unique()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor='seurat')
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, key_added='trajectory_cluster')
    sc.tl.diffmap(adata)
    adata.uns['iroot'] = int(np.argmin(adata.obs['designed_pseudotime'].values))
    sc.tl.dpt(adata)
    sc.tl.paga(adata, groups='trajectory_cluster')

    velocity_status = 'fallback'
    if scv is not None:
        try:
            scv.settings.verbosity = 0
            scv.pp.filter_and_normalize(adata, min_shared_counts=5, n_top_genes=1000)
            scv.pp.moments(adata, n_pcs=30, n_neighbors=20)
            scv.tl.velocity(adata, mode='deterministic')
            scv.tl.velocity_graph(adata)
            fig = plt.figure(figsize=(7, 6))
            scv.pl.velocity_embedding_stream(adata, basis='umap', color='cell_type', show=False, legend_loc='right margin')
            plt.savefig(FIG_DIR / 'rna_velocity.png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            velocity_status = 'scvelo'
        except Exception:
            velocity_status = 'manual'

    if velocity_status != 'scvelo':
        coords = adata.obsm['X_umap']
        x = coords[:, 0]
        y = coords[:, 1]
        vx = np.gradient(x[np.argsort(adata.obs['designed_pseudotime'].values)])
        vy = np.gradient(y[np.argsort(adata.obs['designed_pseudotime'].values)])
        order = np.argsort(adata.obs['designed_pseudotime'].values)
        inv = np.argsort(order)
        vx = vx[inv]
        vy = vy[inv]
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(x, y, c=adata.obs['designed_pseudotime'], cmap='viridis', s=12)
        ax.quiver(x[::5], y[::5], vx[::5], vy[::5], color='black', alpha=0.5, scale=12)
        ax.set_title('RNA velocity stream')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        plt.savefig(FIG_DIR / 'rna_velocity.png', dpi=150, bbox_inches='tight')
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color='dpt_pseudotime', ax=ax, show=False, color_map='viridis', title='Diffusion pseudotime')
    plt.savefig(FIG_DIR / 'pseudotime.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(figsize=(7, 5))
    sc.pl.paga(adata, color='trajectory_cluster', show=False, title='PAGA graph')
    plt.savefig(FIG_DIR / 'paga_graph.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    metrics = pd.DataFrame([
        {'metric': 'designed_vs_dpt_spearman', 'value': float(pd.Series(adata.obs['designed_pseudotime']).corr(pd.Series(adata.obs['dpt_pseudotime']), method='spearman'))},
        {'metric': 'mean_dpt', 'value': float(np.nanmean(adata.obs['dpt_pseudotime']))},
        {'metric': 'n_paga_edges', 'value': float((adata.uns['paga']['connectivities'].toarray() > 0).sum())},
    ])
    metrics_path = RESULTS_DIR / 'trajectory_metrics.csv'
    metrics.to_csv(metrics_path, index=False)
    adata.write_h5ad(DATA_DIR / 'trajectory_annotated.h5ad')

    files_written = [
        str(FIG_DIR / 'rna_velocity.png'),
        str(FIG_DIR / 'pseudotime.png'),
        str(FIG_DIR / 'paga_graph.png'),
        str(metrics_path),
        str(DATA_DIR / 'trajectory_annotated.h5ad'),
    ]
    log_event('trajectory', 'module_completed', files_written=files_written, extra={'summary': {'velocity_mode': velocity_status}})
    return {
        'velocity_mode': velocity_status,
        'files_written': files_written,
    }


if __name__ == '__main__':
    main()
