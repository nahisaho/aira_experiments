import json
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_predict

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

MARKERS = {
    'CD8_T': ['CD3D', 'CD8A', 'CD8B', 'GZMB', 'PRF1'],
    'CD4_T': ['CD3D', 'CD4', 'IL7R', 'FOXP3'],
    'B_cell': ['CD19', 'MS4A1', 'CD79A', 'IGHM'],
    'NK': ['NCAM1', 'KLRD1', 'NKG7', 'GNLY'],
    'Macrophage': ['CD68', 'CD163', 'MRC1', 'CSF1R'],
    'DC': ['ITGAX', 'HLA-DRA', 'CLEC9A', 'FCER1A'],
    'Tumor': ['MKI67', 'TOP2A', 'EPCAM', 'CDH1'],
}


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


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    log_event('tme_classification', 'module_started')

    if sc is None:
        raise ImportError('scanpy is required for TME classification')

    rna = sc.read_h5ad(DATA_DIR / 'rna_processed.h5ad')
    for label, genes in MARKERS.items():
        present = [g for g in genes if g in rna.var_names]
        sc.tl.score_genes(rna, gene_list=present, score_name=f'{label}_score', use_raw=False)

    score_cols = [f'{label}_score' for label in MARKERS]
    score_mat = rna.obs[score_cols].copy()
    rna.obs['marker_prediction'] = score_mat.idxmax(axis=1).str.replace('_score', '', regex=False)

    true_labels = rna.obs['immune_subtype'].astype(str).values
    hvgs = rna.var_names[rna.var.get('highly_variable', pd.Series(True, index=rna.var_names)).values]
    x = np.asarray(rna[:, hvgs].X, dtype=np.float32)
    model = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1, class_weight='balanced')
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    predicted = cross_val_predict(model, x, true_labels, cv=cv, method='predict')
    accuracy = accuracy_score(true_labels, predicted)
    report = classification_report(true_labels, predicted, output_dict=True, zero_division=0)
    cm_labels = list(MARKERS.keys())
    cm = confusion_matrix(true_labels, predicted, labels=cm_labels)

    metrics_rows = []
    for label in cm_labels:
        metrics_rows.append({
            'subtype': label,
            'precision': report.get(label, {}).get('precision', 0.0),
            'recall': report.get(label, {}).get('recall', 0.0),
            'f1_score': report.get(label, {}).get('f1-score', 0.0),
            'support': report.get(label, {}).get('support', 0),
        })
    metrics_rows.append({'subtype': 'overall', 'precision': accuracy, 'recall': accuracy, 'f1_score': accuracy, 'support': len(true_labels)})
    metrics = pd.DataFrame(metrics_rows)
    metrics_path = RESULTS_DIR / 'tme_metrics.csv'
    metrics.to_csv(metrics_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    codes = pd.Categorical(rna.obs['marker_prediction']).codes
    axes[0].scatter(rna.obsm['X_umap'][:, 0], rna.obsm['X_umap'][:, 1], c=codes, cmap='tab10', s=12)
    axes[0].set_title('Marker-based TME classification')
    axes[0].set_xlabel('UMAP1')
    axes[0].set_ylabel('UMAP2')

    axes[1].scatter(rna.obsm['X_umap'][:, 0], rna.obsm['X_umap'][:, 1], c=pd.Categorical(predicted).codes, cmap='tab10', s=12)
    axes[1].set_title('RandomForest TME prediction')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')
    plt.tight_layout()
    fig_path = FIG_DIR / 'tme_classification.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=cm_labels, yticklabels=cm_labels, ax=ax)
    ax.set_title('TME confusion matrix')
    ax.set_xlabel('Predicted label')
    ax.set_ylabel('True label')
    cm_path = FIG_DIR / 'tme_confusion_matrix.png'
    plt.savefig(cm_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    rna.obs['rf_prediction'] = predicted
    rna.write_h5ad(DATA_DIR / 'tme_classified.h5ad')

    files_written = [str(metrics_path), str(fig_path), str(cm_path), str(DATA_DIR / 'tme_classified.h5ad')]
    log_event('tme_classification', 'module_completed', files_written=files_written, extra={'summary': {'accuracy': accuracy}})
    return {
        'accuracy': float(accuracy),
        'files_written': files_written,
    }


if __name__ == '__main__':
    main()
