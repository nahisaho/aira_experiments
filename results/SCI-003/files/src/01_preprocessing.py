import json
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')
np.random.seed(42)

try:
    import scanpy as sc
except Exception:
    sc = None

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / 'data'
FIG_DIR = BASE_DIR / 'figures'
LOG_PATH = BASE_DIR / 'logs' / 'process-log.jsonl'
PREPROCESS_LOG = DATA_DIR / 'preprocessing-log.md'
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


MARKERS = {
    'CD8_T': ['CD3D', 'CD8A', 'CD8B', 'GZMB', 'PRF1'],
    'CD4_T': ['CD3D', 'CD4', 'IL7R', 'FOXP3'],
    'B_cell': ['CD19', 'MS4A1', 'CD79A', 'IGHM'],
    'NK': ['NCAM1', 'KLRD1', 'NKG7', 'GNLY'],
    'Macrophage': ['CD68', 'CD163', 'MRC1', 'CSF1R'],
    'DC': ['ITGAX', 'HLA-DRA', 'CLEC9A', 'FCER1A'],
    'Tumor': ['MKI67', 'TOP2A', 'EPCAM', 'CDH1'],
}
MAIN_TYPE_MAP = {
    'CD8_T': 'T cells',
    'CD4_T': 'T cells',
    'B_cell': 'B cells',
    'NK': 'NK cells',
    'Macrophage': 'Macrophages',
    'DC': 'Dendritic cells',
    'Tumor': 'Tumor cells',
}
CELL_COUNTS = {
    'CD8_T': 80,
    'CD4_T': 80,
    'B_cell': 80,
    'NK': 80,
    'Macrophage': 80,
    'DC': 50,
    'Tumor': 50,
}


def make_gene_list(n_genes=2000):
    mito = [
        'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L',
        'MT-ND5', 'MT-ND6', 'MT-ATP6', 'MT-ATP8', 'MT-CYB', 'MT-RNR1', 'MT-RNR2', 'MT-TF',
        'MT-TV', 'MT-TL1', 'MT-TL2', 'MT-TS', 'MT-TY', 'MT-TR', 'MT-TH', 'MT-TQ', 'MT-TP',
    ]
    genes = []
    for gene in sorted({g for geneset in MARKERS.values() for g in geneset}):
        genes.append(gene)
    genes.extend([g for g in mito if g not in genes])
    idx = 1
    while len(genes) < n_genes:
        gene = f'GENE{idx:04d}'
        if gene not in genes:
            genes.append(gene)
        idx += 1
    return genes[:n_genes]


def simulate_rna():
    genes = make_gene_list()
    n_genes = len(genes)
    cells = []
    immune_subtypes = []
    main_types = []
    for subtype, count in CELL_COUNTS.items():
        for _ in range(count):
            immune_subtypes.append(subtype)
            main_types.append(MAIN_TYPE_MAP[subtype])
            cells.append(f'cell_{len(cells):03d}')

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    base_means = np.random.gamma(shape=1.6, scale=1.4, size=n_genes)
    programs = {}
    pool = [g for g in genes if not g.startswith('MT-') and g not in {m for v in MARKERS.values() for m in v}]
    for subtype in CELL_COUNTS:
        programs[subtype] = np.random.choice(pool, size=120, replace=False)

    x = np.zeros((len(cells), n_genes), dtype=np.float32)
    mito_idx = [gene_to_idx[g] for g in genes if g.startswith('MT-')]
    for i, subtype in enumerate(immune_subtypes):
        lam = base_means.copy()
        size_factor = np.random.lognormal(mean=0.3, sigma=0.35)
        lam *= size_factor
        for gene in MARKERS[subtype]:
            lam[gene_to_idx[gene]] += np.random.uniform(8, 14)
        for gene in programs[subtype]:
            lam[gene_to_idx[gene]] += np.random.uniform(1.5, 3.0)
        if subtype == 'Tumor':
            for gene in ['MKI67', 'TOP2A', 'EPCAM', 'CDH1']:
                lam[gene_to_idx[gene]] += np.random.uniform(10, 16)
        if subtype in {'CD4_T', 'CD8_T'}:
            lam[gene_to_idx['CD3D']] += 8
        mito_scale = 0.08 if np.random.rand() > 0.06 else 0.25
        lam[mito_idx] = np.random.gamma(shape=1.2, scale=mito_scale, size=len(mito_idx)) * size_factor
        noise = np.random.normal(0, 0.15, size=n_genes)
        lam = np.clip(lam + noise, 0.01, None)
        x[i] = np.random.poisson(lam).astype(np.float32)

    obs = pd.DataFrame({
        'cell_type': main_types,
        'immune_subtype': immune_subtypes,
    }, index=cells)
    var = pd.DataFrame(index=genes)
    adata = AnnData(x, obs=obs, var=var)
    adata.layers['counts'] = adata.X.copy()
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    return adata


def simulate_atac(obs):
    n_cells = obs.shape[0]
    n_peaks = 5000
    peaks = [f'peak_{i:05d}' for i in range(n_peaks)]
    base_prob = np.random.beta(0.8, 8, size=n_peaks)
    peak_programs = {}
    for subtype in CELL_COUNTS:
        peak_programs[subtype] = np.random.choice(np.arange(n_peaks), size=500, replace=False)

    x = np.zeros((n_cells, n_peaks), dtype=np.float32)
    for i, subtype in enumerate(obs['immune_subtype'].tolist()):
        p = base_prob.copy()
        p[peak_programs[subtype]] += np.random.uniform(0.18, 0.28)
        if subtype in {'CD4_T', 'CD8_T'}:
            p[np.arange(150)] += 0.12
        p = np.clip(p + np.random.normal(0, 0.01, size=n_peaks), 0.001, 0.95)
        depth = np.random.poisson(2.0, size=n_peaks)
        x[i] = np.random.binomial(np.maximum(depth, 1), p).astype(np.float32)
    adata = AnnData(x, obs=obs.copy(), var=pd.DataFrame(index=peaks))
    adata.layers['counts'] = adata.X.copy()
    return adata


def simulate_methylation(rna_counts, obs, genes):
    n_cells = obs.shape[0]
    n_cpg = 1000
    cpgs = [f'CpG_{i:04d}' for i in range(n_cpg)]
    linked_gene_idx = np.random.choice(np.arange(len(genes)), size=n_cpg, replace=True)
    expr = np.asarray(rna_counts, dtype=np.float32)
    expr_norm = expr / np.maximum(expr.max(axis=0, keepdims=True), 1)

    subtype_offsets = {
        'CD8_T': 0.02,
        'CD4_T': 0.05,
        'B_cell': -0.02,
        'NK': 0.00,
        'Macrophage': 0.08,
        'DC': 0.03,
        'Tumor': -0.10,
    }
    x = np.zeros((n_cells, n_cpg), dtype=np.float32)
    for j in range(n_cpg):
        g_idx = linked_gene_idx[j]
        beta = 0.82 - 0.55 * expr_norm[:, g_idx]
        beta += np.array([subtype_offsets[s] for s in obs['immune_subtype']])
        beta += np.random.normal(0, 0.05, size=n_cells)
        x[:, j] = np.clip(beta, 0.01, 0.99)
    adata = AnnData(x, obs=obs.copy(), var=pd.DataFrame({'linked_gene': [genes[i] for i in linked_gene_idx]}, index=cpgs))
    adata.layers['beta'] = adata.X.copy()
    return adata


def preprocess_rna(adata):
    if sc is None:
        raise ImportError('scanpy is required for RNA preprocessing')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[(adata.obs['n_genes_by_counts'] >= 200) & (adata.obs['n_genes_by_counts'] <= 5000) & (adata.obs['pct_counts_mt'] <= 20)].copy()
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    if adata.var['highly_variable'].sum() < 2000:
        adata.var['highly_variable'] = True
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=20)
    sc.tl.umap(adata)
    return adata


def preprocess_atac(adata):
    x = np.asarray(adata.layers['counts']).astype(np.float32)
    cell_sum = np.maximum(x.sum(axis=1, keepdims=True), 1)
    term_freq = x / cell_sum
    doc_freq = np.maximum((x > 0).sum(axis=0), 1)
    idf = np.log(1 + x.shape[0] / doc_freq)
    tfidf = term_freq * idf
    tfidf = np.nan_to_num(tfidf)
    adata.layers['tfidf'] = tfidf.astype(np.float32)
    svd = TruncatedSVD(n_components=50, random_state=42)
    adata.obsm['X_lsi'] = svd.fit_transform(tfidf)
    if sc is not None:
        sc.pp.neighbors(adata, use_rep='X_lsi', n_neighbors=20)
        sc.tl.umap(adata)
    return adata


def preprocess_methylation(adata):
    x = np.asarray(adata.layers['beta']).astype(np.float32)
    variances = x.var(axis=0)
    threshold = np.quantile(variances, 0.01)
    keep = variances > threshold
    if keep.sum() < x.shape[1]:
        keep = np.ones(x.shape[1], dtype=bool)
    adata = adata[:, keep].copy()
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(np.asarray(adata.X, dtype=np.float32))
    adata.X = x_scaled.astype(np.float32)
    pca = PCA(n_components=min(50, adata.n_vars, adata.n_obs - 1), random_state=42)
    adata.obsm['X_pca'] = pca.fit_transform(x_scaled)
    if sc is not None:
        sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=20)
        sc.tl.umap(adata)
    return adata


def save_preprocessing_summary(rna, atac, meth):
    with open(PREPROCESS_LOG, 'w', encoding='utf-8') as fh:
        fh.write('# Preprocessing Log\n\n')
        fh.write(f'- Timestamp: {datetime.utcnow().isoformat()}\n')
        fh.write('- RNA: QC filters (min_genes=200, max_genes=5000, max_mito_pct=20), normalize_total, log1p, HVG, PCA, UMAP\n')
        fh.write('- ATAC: TF-IDF normalization, LSI (TruncatedSVD), UMAP\n')
        fh.write('- Methylation: low-variance CpG check, standardization, PCA, UMAP\n')
        fh.write(f'- RNA cells retained: {rna.n_obs}, genes: {rna.n_vars}\n')
        fh.write(f'- ATAC cells retained: {atac.n_obs}, peaks: {atac.n_vars}\n')
        fh.write(f'- Methylation cells retained: {meth.n_obs}, CpGs: {meth.n_vars}\n')


def plot_modalities(rna, atac, meth):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.ravel()
    axes[0].scatter(rna.obs['total_counts'], rna.obs['pct_counts_mt'], c='steelblue', s=12, alpha=0.7)
    axes[0].set_title('RNA QC')
    axes[0].set_xlabel('Total counts')
    axes[0].set_ylabel('Mito percent')

    axes[1].scatter(rna.obsm['X_umap'][:, 0], rna.obsm['X_umap'][:, 1], c=pd.Categorical(rna.obs['cell_type']).codes, cmap='tab10', s=12)
    axes[1].set_title('RNA UMAP')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')

    axes[2].scatter(atac.obsm['X_umap'][:, 0], atac.obsm['X_umap'][:, 1], c=pd.Categorical(atac.obs['cell_type']).codes, cmap='tab10', s=12)
    axes[2].set_title('ATAC UMAP')
    axes[2].set_xlabel('UMAP1')
    axes[2].set_ylabel('UMAP2')

    axes[3].hist(rna.obs['n_genes_by_counts'], bins=30, color='slateblue', alpha=0.8)
    axes[3].set_title('Detected genes per cell')
    axes[3].set_xlabel('Genes')
    axes[3].set_ylabel('Cell count')

    axes[4].scatter(meth.obsm['X_umap'][:, 0], meth.obsm['X_umap'][:, 1], c=pd.Categorical(meth.obs['cell_type']).codes, cmap='tab10', s=12)
    axes[4].set_title('Methylation UMAP')
    axes[4].set_xlabel('UMAP1')
    axes[4].set_ylabel('UMAP2')

    axes[5].hist(np.asarray(meth.layers['beta']).ravel(), bins=40, color='darkorange', alpha=0.8)
    axes[5].set_title('Methylation beta values')
    axes[5].set_xlabel('Beta value')
    axes[5].set_ylabel('Frequency')

    plt.tight_layout()
    out = FIG_DIR / 'preprocessing_overview.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return [str(out)]


def main():
    DATA_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    LOG_PATH.parent.mkdir(exist_ok=True)
    log_event('preprocessing', 'module_started')

    rna = simulate_rna()
    rna_processed = preprocess_rna(rna)
    atac = simulate_atac(rna_processed.obs)
    atac_processed = preprocess_atac(atac)
    meth = simulate_methylation(rna_processed.layers['counts'], rna_processed.obs, rna_processed.var_names.tolist())
    meth_processed = preprocess_methylation(meth)

    files_written = []
    rna_path = DATA_DIR / 'rna_processed.h5ad'
    atac_path = DATA_DIR / 'atac_processed.h5ad'
    meth_path = DATA_DIR / 'methylation_processed.h5ad'
    meta_path = DATA_DIR / 'cell_metadata.csv'
    rna_processed.write_h5ad(rna_path)
    atac_processed.write_h5ad(atac_path)
    meth_processed.write_h5ad(meth_path)
    rna_processed.obs.to_csv(meta_path)
    files_written.extend([str(rna_path), str(atac_path), str(meth_path), str(meta_path)])
    files_written.extend(plot_modalities(rna_processed, atac_processed, meth_processed))
    save_preprocessing_summary(rna_processed, atac_processed, meth_processed)
    files_written.append(str(PREPROCESS_LOG))

    log_event(
        'preprocessing',
        'module_completed',
        files_written=files_written,
        extra={
            'summary': {
                'rna_shape': [int(rna_processed.n_obs), int(rna_processed.n_vars)],
                'atac_shape': [int(atac_processed.n_obs), int(atac_processed.n_vars)],
                'methylation_shape': [int(meth_processed.n_obs), int(meth_processed.n_vars)],
            }
        },
    )
    return {
        'rna_cells': int(rna_processed.n_obs),
        'rna_genes': int(rna_processed.n_vars),
        'atac_peaks': int(atac_processed.n_vars),
        'methylation_cpgs': int(meth_processed.n_vars),
        'files_written': files_written,
    }


if __name__ == '__main__':
    main()
