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
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import mutual_info_regression

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


MARKER_GENES = [
    'CD3D', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'FOXP3', 'CD19', 'MS4A1', 'CD79A', 'IGHM',
    'NCAM1', 'KLRD1', 'NKG7', 'GNLY', 'CD68', 'CD163', 'MRC1', 'CSF1R', 'ITGAX', 'HLA-DRA',
    'CLEC9A', 'FCER1A', 'MKI67', 'TOP2A', 'EPCAM', 'CDH1'
]


def top_edges_from_matrix(matrix, genes, top_n=500, threshold=0.0, absolute=True):
    tri = np.triu_indices_from(matrix, k=1)
    values = matrix[tri]
    score = np.abs(values) if absolute else values
    keep = score > threshold
    tri_i = tri[0][keep]
    tri_j = tri[1][keep]
    vals = values[keep]
    order = np.argsort(np.abs(vals))[::-1][:top_n]
    edges = [(genes[tri_i[k]], genes[tri_j[k]], float(vals[k])) for k in order]
    return edges


def network_metrics(edges):
    graph = nx.Graph()
    for a, b, w in edges:
        graph.add_edge(a, b, weight=abs(w))
    n = max(graph.number_of_nodes(), 1)
    density = nx.density(graph) if n > 1 else 0.0
    degrees = np.array([deg for _, deg in graph.degree()]) if graph.number_of_nodes() else np.array([0])
    avg_degree = float(degrees.mean()) if len(degrees) else 0.0
    hubs = int((degrees > (degrees.mean() + 2 * degrees.std())).sum()) if len(degrees) else 0
    return density, avg_degree, hubs, graph


def edge_set(edges):
    return {tuple(sorted((a, b))) for a, b, _ in edges}


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    log_event('grn_inference', 'module_started')

    if sc is None:
        raise ImportError('scanpy is required to load AnnData objects')

    rna = sc.read_h5ad(DATA_DIR / 'rna_processed.h5ad')
    var_df = pd.DataFrame({
        'gene': rna.var_names,
        'variance': np.var(np.asarray(rna.X, dtype=np.float32), axis=0),
    }).sort_values('variance', ascending=False)
    selected = list(dict.fromkeys(MARKER_GENES + var_df['gene'].head(124).tolist()))
    selected = [g for g in selected if g in rna.var_names][:150]
    x = np.asarray(rna[:, selected].X, dtype=np.float32)

    corr = np.corrcoef(x.T)
    corr_edges = top_edges_from_matrix(corr, selected, top_n=500, threshold=0.3, absolute=True)

    mi_matrix = np.zeros((len(selected), len(selected)), dtype=np.float32)
    for i in range(len(selected)):
        target = x[:, i]
        mi = mutual_info_regression(x, target, random_state=42)
        mi_matrix[i, :] = mi
    mi_matrix = 0.5 * (mi_matrix + mi_matrix.T)
    np.fill_diagonal(mi_matrix, 0)
    mi_edges = top_edges_from_matrix(mi_matrix, selected, top_n=500, threshold=0.0, absolute=False)

    tf_genes = list(np.random.choice(selected, size=20, replace=False))
    rf_edges = []
    predictors = x[:, [selected.index(g) for g in tf_genes]]
    target_genes = [g for g in selected if g not in tf_genes][:80]
    for target_gene in target_genes:
        y = x[:, selected.index(target_gene)]
        model = RandomForestRegressor(n_estimators=60, max_depth=6, random_state=42, n_jobs=-1)
        model.fit(predictors, y)
        for tf_gene, importance in zip(tf_genes, model.feature_importances_):
            rf_edges.append((tf_gene, target_gene, float(importance)))
    rf_edges = sorted(rf_edges, key=lambda z: abs(z[2]), reverse=True)[:500]

    rows = []
    graphs = {}
    for method, edges in [('Pearson', corr_edges), ('MutualInformation', mi_edges), ('RandomForest', rf_edges)]:
        density, avg_degree, hubs, graph = network_metrics(edges)
        graphs[method] = graph
        rows.append({'method': method, 'metric': 'network_density', 'value': density})
        rows.append({'method': method, 'metric': 'average_degree', 'value': avg_degree})
        rows.append({'method': method, 'metric': 'hub_genes', 'value': hubs})
        rows.append({'method': method, 'metric': 'edge_count', 'value': len(edges)})

    jaccard_pairs = [
        ('Pearson', 'MutualInformation', corr_edges, mi_edges),
        ('Pearson', 'RandomForest', corr_edges, rf_edges),
        ('MutualInformation', 'RandomForest', mi_edges, rf_edges),
    ]
    for a, b, ea, eb in jaccard_pairs:
        sa, sb = edge_set(ea), edge_set(eb)
        jac = len(sa & sb) / max(len(sa | sb), 1)
        rows.append({'method': f'{a}_vs_{b}', 'metric': 'jaccard', 'value': jac})

    metrics = pd.DataFrame(rows)
    metrics_path = RESULTS_DIR / 'grn_metrics.csv'
    metrics.to_csv(metrics_path, index=False)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    layout_seed = 42
    for ax, (method, edges) in zip(axes[0], [('Pearson', corr_edges[:50]), ('MutualInformation', mi_edges[:50]), ('RandomForest', rf_edges[:50])]):
        graph = nx.Graph()
        for a, b, w in edges:
            graph.add_edge(a, b, weight=abs(w))
        pos = nx.spring_layout(graph, seed=layout_seed, k=0.45)
        widths = [0.5 + 2.5 * d['weight'] / max([ed['weight'] for _, _, ed in graph.edges(data=True)] + [1]) for _, _, d in graph.edges(data=True)]
        nx.draw_networkx(graph, pos=pos, ax=ax, with_labels=False, node_size=55, width=widths, edge_color='gray', node_color='skyblue')
        ax.set_title(f'{method} top-50 edges')
        ax.axis('off')

    summary = metrics.pivot_table(index='metric', columns='method', values='value', aggfunc='first')
    summary.loc[['network_density', 'average_degree', 'hub_genes']].T.plot(kind='bar', ax=axes[1, 0], colormap='viridis')
    axes[1, 0].set_title('GRN metric comparison')
    axes[1, 0].set_xlabel('Method')
    axes[1, 0].set_ylabel('Value')
    axes[1, 0].tick_params(axis='x', rotation=25)

    jac_df = metrics[metrics['metric'] == 'jaccard'].copy()
    axes[1, 1].bar(jac_df['method'], jac_df['value'], color=['#4c72b0', '#55a868', '#c44e52'])
    axes[1, 1].set_title('Jaccard similarity')
    axes[1, 1].set_xlabel('Method pair')
    axes[1, 1].set_ylabel('Jaccard')
    axes[1, 1].tick_params(axis='x', rotation=20)

    top_genes = sorted(nx.degree_centrality(graphs['Pearson']).items(), key=lambda x: x[1], reverse=True)[:10]
    axes[1, 2].bar([g for g, _ in top_genes], [v for _, v in top_genes], color='mediumpurple')
    axes[1, 2].set_title('Pearson hub genes')
    axes[1, 2].set_xlabel('Gene')
    axes[1, 2].set_ylabel('Degree centrality')
    axes[1, 2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    fig_path = FIG_DIR / 'grn_comparison.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    files_written = [str(metrics_path), str(fig_path)]
    log_event('grn_inference', 'module_completed', files_written=files_written, extra={'summary': {'tf_genes': tf_genes}})
    return {
        'pearson_edges': len(corr_edges),
        'mi_edges': len(mi_edges),
        'rf_edges': len(rf_edges),
        'files_written': files_written,
    }


if __name__ == '__main__':
    main()
