import importlib.util
import json
import traceback
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = BASE_DIR / 'src'
LOG_PATH = BASE_DIR / 'logs' / 'process-log.jsonl'
RESULTS_DIR = BASE_DIR / 'results'
FIG_DIR = BASE_DIR / 'figures'
DATA_DIR = BASE_DIR / 'data'
SKILL_NAME = 'co-scientist-multi-omics'

MODULES = [
    ('preprocessing', SRC_DIR / '01_preprocessing.py'),
    ('integration', SRC_DIR / '02_integration.py'),
    ('vae_integration', SRC_DIR / '03_vae_integration.py'),
    ('trajectory', SRC_DIR / '04_trajectory.py'),
    ('grn_inference', SRC_DIR / '05_grn_inference.py'),
    ('tme_classification', SRC_DIR / '06_tme_classification.py'),
]


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


def load_module(module_path, name):
    spec = importlib.util.spec_from_file_location(name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect_generated_files():
    files = []
    for root in [DATA_DIR, RESULTS_DIR, FIG_DIR, BASE_DIR / 'logs']:
        if root.exists():
            files.extend(sorted(str(path.relative_to(BASE_DIR)) for path in root.rglob('*') if path.is_file()))
    return files


def write_report(summaries):
    integration = pd.read_csv(RESULTS_DIR / 'integration_metrics.csv') if (RESULTS_DIR / 'integration_metrics.csv').exists() else pd.DataFrame()
    trajectory = pd.read_csv(RESULTS_DIR / 'trajectory_metrics.csv') if (RESULTS_DIR / 'trajectory_metrics.csv').exists() else pd.DataFrame()
    grn = pd.read_csv(RESULTS_DIR / 'grn_metrics.csv') if (RESULTS_DIR / 'grn_metrics.csv').exists() else pd.DataFrame()
    tme = pd.read_csv(RESULTS_DIR / 'tme_metrics.csv') if (RESULTS_DIR / 'tme_metrics.csv').exists() else pd.DataFrame()

    joint_sil = integration.loc[integration['embedding'] == 'Joint_WNN', 'silhouette_score'].iloc[0] if not integration.empty else np.nan
    vae_loss = summaries.get('vae_integration', {}).get('final_loss', np.nan)
    trajectory_corr = trajectory.loc[trajectory['metric'] == 'designed_vs_dpt_spearman', 'value'].iloc[0] if not trajectory.empty else np.nan
    tme_acc = tme.loc[tme['subtype'] == 'overall', 'f1_score'].iloc[0] if not tme.empty else np.nan
    grn_jaccard = grn.loc[grn['metric'] == 'jaccard', 'value'].mean() if not grn.empty else np.nan
    generated_files = collect_generated_files()

    report_lines = [
        '# DRAFT — NOT FOR DISTRIBUTION',
        '',
        '## 実験目的と背景',
        '本解析では、実データ未取得の研究初期段階を想定し、合成 single-cell multi-omics データを用いて scRNA-seq、scATAC-seq、DNA methylation を統合解析する再現可能な Python パイプラインを構築した。目的は、前処理、WNN 型統合、VAE 型潜在表現学習、RNA velocity と擬似時間解析、GRN 推定比較、ならびに腫瘍微小環境 (TME) の免疫細胞サブタイプ分類を一連で実行し、後続の実データ適用に耐える雛形を作成することである。',
        '',
        '## 使用した手法・アルゴリズムの概要',
        '1. **Preprocessing**: scRNA-seq は QC、library-size normalization、log1p、HVG 選択、PCA、UMAP を実施した。scATAC-seq は TF-IDF 正規化後に LSI を算出し、methylation は低分散 CpG の点検、標準化、PCA、UMAP を実施した。',
        '2. **Anchor-based integration**: 各モダリティの近傍構造からセルごとのモダリティ重みを求め、加重 joint graph と joint embedding を構築する WNN-style 統合を実装した。',
        '3. **VAE-based integration**: PyTorch により RNA / ATAC / methylation の別個 encoder-decoder と共有 32 次元 latent space を持つ MultimodalVAE を 50 epoch 学習した。',
        '4. **Trajectory analysis**: 合成 spliced/unspliced counts を生成し、scVelo style の velocity 推定（失敗時は manual fallback）と Scanpy の diffusion pseudotime、PAGA を実行した。',
        '5. **GRN inference**: Pearson correlation、mutual information、RandomForest importance (GENIE3-style) を比較し、密度・平均次数・hub gene 数・Jaccard 類似度を算出した。',
        '6. **TME classification**: marker gene scoring と RandomForestClassifier を併用し、免疫細胞サブタイプ予測と 5-fold CV 評価を実施した。',
        '',
        '## 主要な結果と数値',
        f'- Joint WNN silhouette score: **{joint_sil:.3f}**',
        f'- Multimodal VAE final loss: **{vae_loss:.4f}**',
        f'- Designed pseudotime vs DPT Spearman correlation: **{trajectory_corr:.3f}**',
        f'- Mean GRN Jaccard similarity across method pairs: **{grn_jaccard:.3f}**',
        f'- TME classification overall F1 / accuracy: **{tme_acc:.3f}**',
        '',
        '主要出力として、統合 UMAP、VAE training curve、RNA velocity / pseudotime / PAGA 図、GRN 比較図、TME 分類図および confusion matrix を保存した。各モジュールの中間成果物は `data/` と `results/` に保存し、実行ログは `logs/process-log.jsonl` に記録した。',
        '',
        '## 考察と今後の展望',
        '本パイプラインは、合成データ上で multi-omics 統合から downstream 解析まで一貫して完走できることを確認した。一方で、合成データは実際の read depth、dropout、batch effect、クロマチン accessibility の局所依存性、methylation の文脈依存性を単純化しているため、得られた性能値は過大評価の可能性がある。今後は、(1) 実データへの置換、(2) donor/batch の明示的導入、(3) SCENIC 本体や pySCENIC への拡張、(4) scVI/totalVI との比較、(5) 外部参照データセットを用いた TME ラベル転移評価、を進めることで妥当性を高められる。',
        '',
        '## 生成したファイル一覧',
        '以下の主要ファイルを生成した。',
        '',
    ]
    for path in generated_files:
        report_lines.append(f'- `{path}`')

    report_lines.extend([
        '',
        '### Figure captions (English)',
        '- `figures/preprocessing_overview.png`: Multi-omics preprocessing overview',
        '- `figures/wnn_integration.png`: WNN integration UMAP and modality weights',
        '- `figures/vae_training_loss.png`: Multimodal VAE training loss',
        '- `figures/vae_latent_umap.png`: VAE latent UMAP',
        '- `figures/rna_velocity.png`: RNA velocity stream',
        '- `figures/pseudotime.png`: Diffusion pseudotime map',
        '- `figures/paga_graph.png`: PAGA topology graph',
        '- `figures/grn_comparison.png`: GRN method comparison',
        '- `figures/tme_classification.png`: Marker-based and RandomForest TME classification',
        '- `figures/tme_confusion_matrix.png`: TME confusion matrix',
    ])

    report_path = BASE_DIR / 'report.md'
    with open(report_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(report_lines) + '\n')
    return report_path


def main():
    BASE_DIR.joinpath('logs').mkdir(exist_ok=True)
    log_event('pipeline', 'run_started', extra={'handoff_in': {'modules': [name for name, _ in MODULES]}})
    summaries = {}
    for name, module_path in MODULES:
        try:
            log_event(name, 'handoff_started', extra={'handoff_in': {'module_path': str(module_path)}})
            module = load_module(module_path, name)
            result = module.main()
            summaries[name] = result
            log_event(name, 'handoff_completed', files_written=result.get('files_written', []), extra={'handoff_out': result})
        except Exception as exc:
            log_event(name, 'module_failed', status='error', extra={'error': str(exc), 'traceback': traceback.format_exc()})
            raise

    report_path = write_report(summaries)
    log_event('pipeline', 'report_finalized', files_written=[str(report_path)])
    log_event('pipeline', 'run_completed', files_written=collect_generated_files(), extra={'summary': summaries})

    print('Pipeline completed successfully.')
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'files_written'} for k, v in summaries.items()}, ensure_ascii=False, indent=2))
    print(f'Report written to: {report_path}')


if __name__ == '__main__':
    main()
