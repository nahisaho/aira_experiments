"""
CRISPR-Cas9 Off-Target Prediction — Main Pipeline Runner.

Executes:
  1. Synthetic data generation (stand-in for GUIDE-seq / CIRCLE-seq)
  2. Feature engineering
  3. 5-fold cross-validation training
  4. Evaluation plots (ROC, PR, CV summary)
  5. Lightweight SHAP analysis
  6. Attention map visualisation
  7. Data flow diagram
  8. Results serialisation
"""

import sys
import os
import json
import logging
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import torch

# ── Ensure workspace root is on sys.path ──────────────────────────────────────
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import (
    generate_synthetic_dataset, CRISPRFeatureBuilder,
)
from src.model    import CRISPROffTargetModel
from src.train    import cross_validate, OffTargetDataset
from src.evaluate import (
    compute_metrics, find_optimal_threshold,
    plot_roc_curves, plot_pr_curves, plot_cv_summary,
    plot_mismatch_importance,
)
from src.interpretability import (
    CRISPRSHAPExplainer, extract_attention_maps,
    plot_positional_shap_heatmap, plot_attention_heatmap,
    save_shap_summary, SCALAR_FEATURE_NAMES,
)
from src.dataflow_diagram import create_dataflow_diagram

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"
LOGS_DIR    = ROOT / "logs"
for d in [RESULTS_DIR, FIGURES_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOGS_DIR / "process-log.jsonl"


def log_event(phase: str, event_type: str, skill: str, **kwargs):
    entry = {
        "timestamp":    datetime.utcnow().isoformat() + "Z",
        "phase":        phase,
        "event_type":   event_type,
        "actor":        "co-scientist",
        "skill_or_tool": skill,
        **kwargs,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


# ─────────────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    log_event("pipeline", "run_started", "co-scientist-crispr-design",
              status="ok", note="CRISPR off-target prediction pipeline started")

    # ── 1. Data generation ────────────────────────────────────────────────────
    logger.info("=== Phase 1: Data Generation ===")
    log_event("data", "handoff_started", "preprocessing",
              note="Generating synthetic GUIDE-seq/CIRCLE-seq data")

    df = generate_synthetic_dataset(n_guides=50, sites_per_guide=100, seed=42)
    df.to_csv(RESULTS_DIR / "synthetic_dataset.csv", index=False)

    log_event("data", "file_written", "preprocessing",
              files_written=["results/synthetic_dataset.csv"],
              n_samples=len(df), pos_rate=float(df["label"].mean()))

    # ── 2. Feature engineering ────────────────────────────────────────────────
    logger.info("=== Phase 2: Feature Engineering ===")
    builder = CRISPRFeatureBuilder(include_epigenetics=True)
    X_seq, X_scalar, y = builder.transform(df)

    feature_info = {
        "X_seq_shape":    list(X_seq.shape),
        "X_scalar_shape": list(X_scalar.shape),
        "y_shape":        list(y.shape),
        "pos_rate":       float(y.mean()),
        "seq_channels":   "guide_OH(4) + target_OH(4) + mismatch_type(15) = 23",
        "scalar_channels": "positional_mismatch(23) + epigenetic(8) = 31",
    }
    with open(RESULTS_DIR / "feature_info.json", "w") as f:
        json.dump(feature_info, f, indent=2)
    log_event("features", "file_written", "preprocessing",
              files_written=["results/feature_info.json"], **feature_info)

    # ── 3. Cross-validation training ──────────────────────────────────────────
    logger.info("=== Phase 3: 5-Fold Cross-Validation Training ===")
    log_event("training", "handoff_started", "train", note="Starting CV training")

    fold_metrics = cross_validate(
        X_seq, X_scalar, y,
        n_splits=5, epochs=40, batch_size=64, lr=3e-4,
        seed=42, save_dir=str(RESULTS_DIR),
    )

    aurocs = [m["best_auroc"] for m in fold_metrics]
    auprcs = [m["auprc"]      for m in fold_metrics]
    f1s    = [m["f1"]         for m in fold_metrics]

    cv_summary = {
        "mean_auroc": float(np.mean(aurocs)),
        "std_auroc":  float(np.std(aurocs)),
        "mean_auprc": float(np.mean(auprcs)),
        "std_auprc":  float(np.std(auprcs)),
        "mean_f1":    float(np.mean(f1s)),
        "std_f1":     float(np.std(f1s)),
    }
    logger.info("CV Summary: AUROC=%.4f±%.4f, AUPRC=%.4f±%.4f",
                cv_summary["mean_auroc"], cv_summary["std_auroc"],
                cv_summary["mean_auprc"], cv_summary["std_auprc"])

    log_event("training", "handoff_completed", "train",
              files_written=["results/cv_results.json"], **cv_summary)

    # ── 4. Final model on full data (for evaluation & SHAP) ───────────────────
    logger.info("=== Phase 4: Final Model Evaluation ===")
    device = torch.device("cpu")
    model  = CRISPROffTargetModel(
        seq_in_channels=X_seq.shape[2],
        seq_length=X_seq.shape[1],
        scalar_dim=X_scalar.shape[1],
    ).to(device)

    # Quick train for demo (use best fold checkpoint in production)
    from src.train import train_epoch, evaluate as eval_model, OffTargetDataset
    from torch.utils.data import DataLoader, random_split

    full_ds  = OffTargetDataset(X_seq, X_scalar, y)
    n_val    = int(0.2 * len(full_ds))
    n_train  = len(full_ds) - n_val
    train_ds, val_ds = random_split(full_ds, [n_train, n_val],
                                    generator=torch.Generator().manual_seed(42))
    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
    val_loader   = DataLoader(val_ds,   batch_size=64)

    optim = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optim, T_max=40)
    for epoch in range(40):
        train_epoch(model, train_loader, optim, device)
        scheduler.step()

    # Collect predictions on validation set
    model.eval()
    all_probs, all_labels = [], []
    with torch.no_grad():
        for xs, xsc, yl in val_loader:
            logits, _ = model(xs.to(device), xsc.to(device))
            all_probs.extend(torch.sigmoid(logits).cpu().numpy())
            all_labels.extend(yl.numpy())

    y_prob = np.array(all_probs)
    y_true = np.array(all_labels).astype(int)
    opt_t  = find_optimal_threshold(y_true, y_prob)
    metrics = compute_metrics(y_true, y_prob, threshold=opt_t)
    metrics["optimal_threshold"] = float(opt_t)
    metrics["n_params"] = model.count_parameters()

    with open(RESULTS_DIR / "final_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    logger.info("Final metrics: AUROC=%.4f  AUPRC=%.4f  F1=%.4f",
                metrics["auroc"], metrics["auprc"], metrics["f1"])
    log_event("evaluation", "file_written", "evaluate",
              files_written=["results/final_metrics.json"], **metrics)

    # ── 5. Evaluation plots ───────────────────────────────────────────────────
    logger.info("=== Phase 5: Evaluation Plots ===")
    roc_data = [{"y_true": y_true, "y_prob": y_prob, "label": "CNN+Attention"}]
    plot_roc_curves(roc_data, save_path=str(FIGURES_DIR / "roc_curves.png"))
    plot_pr_curves(roc_data,  save_path=str(FIGURES_DIR / "pr_curves.png"))
    plot_cv_summary(fold_metrics, save_path=str(FIGURES_DIR / "cv_summary.png"))
    log_event("evaluation", "file_written", "evaluate",
              files_written=["figures/roc_curves.png", "figures/pr_curves.png",
                             "figures/cv_summary.png"])

    # ── 6. SHAP interpretability ──────────────────────────────────────────────
    logger.info("=== Phase 6: SHAP Interpretability ===")
    log_event("interpretability", "handoff_started", "shap",
              note="Running KernelSHAP on 30 validation samples")

    val_idx  = val_ds.indices
    bg_idx   = train_ds.indices
    x_seq_v  = X_seq[val_idx[:30]]
    x_scal_v = X_scalar[val_idx[:30]]
    y_v      = y[val_idx[:30]]

    explainer  = CRISPRSHAPExplainer(
        model, X_seq[bg_idx], X_scalar[bg_idx],
        n_background=50, seed=42,
    )
    shap_vals  = explainer.explain(x_seq_v, x_scal_v, n_samples=50)
    pos_shap   = explainer.aggregate_positional(shap_vals)
    scal_shap  = explainer.aggregate_scalar(shap_vals)

    # Save SHAP summary
    scalar_feature_vals = x_scal_v
    save_shap_summary(scal_shap, SCALAR_FEATURE_NAMES,
                      save_path=str(RESULTS_DIR / "shap_summary.json"))

    # Positional SHAP heatmap
    plot_positional_shap_heatmap(
        pos_shap, y_v,
        save_path=str(FIGURES_DIR / "positional_shap_heatmap.png"),
    )

    # Per-position mean |SHAP|
    from src.evaluate import plot_mismatch_importance
    plot_mismatch_importance(
        pos_shap.mean(axis=0),
        save_path=str(FIGURES_DIR / "mismatch_importance.png"),
    )

    log_event("interpretability", "handoff_completed", "shap",
              files_written=["results/shap_summary.json",
                             "figures/positional_shap_heatmap.png",
                             "figures/mismatch_importance.png"])

    # ── 7. Attention map ──────────────────────────────────────────────────────
    logger.info("=== Phase 7: Attention Map Visualisation ===")
    x_seq_t  = torch.from_numpy(x_seq_v[:1]).float()
    x_scal_t = torch.from_numpy(x_scal_v[:1]).float()
    attn = extract_attention_maps(model, x_seq_t, x_scal_t)
    plot_attention_heatmap(attn, save_path=str(FIGURES_DIR / "attention_heatmap.png"))
    log_event("interpretability", "file_written", "attention",
              files_written=["figures/attention_heatmap.png"])

    # ── 8. Data flow diagram ──────────────────────────────────────────────────
    logger.info("=== Phase 8: Data Flow Diagram ===")
    create_dataflow_diagram(save_path=str(FIGURES_DIR / "dataflow_diagram.png"))
    log_event("visualization", "file_written", "dataflow",
              files_written=["figures/dataflow_diagram.png"])

    # ── 9. Benchmark plan ─────────────────────────────────────────────────────
    benchmark = {
        "description": "Performance benchmark plan for CRISPR off-target prediction",
        "datasets": [
            {"name": "GUIDE-seq (Tsai 2015)", "n_sites": 702,   "n_guides": 13,  "assay": "in-vivo"},
            {"name": "CIRCLE-seq (Tsai 2017)", "n_sites": 1574, "n_guides": 10,  "assay": "in-vitro"},
            {"name": "SITE-seq (Cameron 2017)", "n_sites": 840, "n_guides": 8,   "assay": "in-vitro"},
            {"name": "CHANGE-seq (Lazzarotto 2020)", "n_sites": 9340, "n_guides": 110, "assay": "in-vitro"},
        ],
        "baselines": ["Cas-OFFinder", "CRISPOR", "DeepCRISPR", "CRISPR-ML"],
        "metrics": ["AUROC", "AUPRC", "Recall@10%FPR", "Precision@50%Recall", "F1"],
        "cross_validation": "Leave-one-guide-out (LOGO) + 5-fold stratified",
        "statistical_test": "DeLong test for AUROC comparison",
        "hardware": "NVIDIA A100 (40GB) or CPU fallback",
        "estimated_runtime_minutes": {"train_5fold": 30, "shap_100samples": 15},
    }
    with open(RESULTS_DIR / "benchmark_plan.json", "w") as f:
        json.dump(benchmark, f, indent=2)

    # ── 10. Final summary ─────────────────────────────────────────────────────
    elapsed = time.time() - t0
    summary = {
        "pipeline_version": "1.0.0",
        "elapsed_seconds":  round(elapsed, 1),
        "cv_auroc":    f"{cv_summary['mean_auroc']:.4f} ± {cv_summary['std_auroc']:.4f}",
        "cv_auprc":    f"{cv_summary['mean_auprc']:.4f} ± {cv_summary['std_auprc']:.4f}",
        "final_auroc": f"{metrics['auroc']:.4f}",
        "final_auprc": f"{metrics['auprc']:.4f}",
        "n_params":    metrics["n_params"],
        "files": {
            "code":   ["src/preprocessing.py", "src/model.py", "src/train.py",
                       "src/evaluate.py", "src/interpretability.py",
                       "src/dataflow_diagram.py", "run_pipeline.py"],
            "results": ["results/synthetic_dataset.csv", "results/feature_info.json",
                        "results/cv_results.json", "results/final_metrics.json",
                        "results/shap_summary.json", "results/benchmark_plan.json"],
            "figures": ["figures/dataflow_diagram.png", "figures/roc_curves.png",
                        "figures/pr_curves.png", "figures/cv_summary.png",
                        "figures/positional_shap_heatmap.png",
                        "figures/mismatch_importance.png",
                        "figures/attention_heatmap.png"],
        },
    }
    with open(RESULTS_DIR / "pipeline_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log_event("pipeline", "run_completed", "co-scientist-crispr-design",
              status="ok", elapsed_seconds=elapsed, **cv_summary)

    logger.info("Pipeline completed in %.1f seconds.", elapsed)
    logger.info("CV AUROC: %s", summary["cv_auroc"])
    logger.info("CV AUPRC: %s", summary["cv_auprc"])
    return summary, cv_summary, metrics


if __name__ == "__main__":
    summary, cv_summary, metrics = main()
