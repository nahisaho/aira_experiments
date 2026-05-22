"""
CRISPR-Cas9 Off-Target Prediction — Evaluation Metrics and Visualisations.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, auc,
    precision_recall_curve, average_precision_score,
    roc_auc_score, confusion_matrix, ConfusionMatrixDisplay,
)
from typing import List, Dict, Tuple
from pathlib import Path
import json


# ─── Metric Computation ───────────────────────────────────────────────────────

def compute_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    threshold: float = 0.5,
) -> Dict[str, float]:
    """Comprehensive binary classification metrics."""
    y_pred = (y_prob >= threshold).astype(int)
    tp = ((y_pred == 1) & (y_true == 1)).sum()
    fp = ((y_pred == 1) & (y_true == 0)).sum()
    fn = ((y_pred == 0) & (y_true == 1)).sum()
    tn = ((y_pred == 0) & (y_true == 0)).sum()

    precision = tp / (tp + fp + 1e-8)
    recall    = tp / (tp + fn + 1e-8)
    specificity = tn / (tn + fp + 1e-8)
    f1        = 2 * precision * recall / (precision + recall + 1e-8)
    mcc_num   = tp * tn - fp * fn
    mcc_den   = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)
    mcc       = mcc_num / mcc_den

    return {
        "auroc":       float(roc_auc_score(y_true, y_prob)),
        "auprc":       float(average_precision_score(y_true, y_prob)),
        "precision":   float(precision),
        "recall":      float(recall),
        "specificity": float(specificity),
        "f1":          float(f1),
        "mcc":         float(mcc),
        "tp": int(tp), "fp": int(fp), "fn": int(fn), "tn": int(tn),
    }


def find_optimal_threshold(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    """Find threshold maximising F1 score on the validation set."""
    thresholds = np.linspace(0.01, 0.99, 100)
    best_f1, best_t = 0.0, 0.5
    for t in thresholds:
        y_pred = (y_prob >= t).astype(int)
        tp = ((y_pred == 1) & (y_true == 1)).sum()
        fp = ((y_pred == 1) & (y_true == 0)).sum()
        fn = ((y_pred == 0) & (y_true == 1)).sum()
        prec = tp / (tp + fp + 1e-8)
        rec  = tp / (tp + fn + 1e-8)
        f1   = 2 * prec * rec / (prec + rec + 1e-8)
        if f1 > best_f1:
            best_f1, best_t = f1, t
    return best_t


# ─── Plotting ─────────────────────────────────────────────────────────────────

PALETTE = {
    "primary":   "#0072B2",
    "secondary": "#E69F00",
    "accent":    "#009E73",
    "danger":    "#D55E00",
    "gray":      "#999999",
}


def plot_roc_curves(
    results: List[Dict],   # list of {"y_true": ..., "y_prob": ..., "label": ...}
    save_path: str = "figures/roc_curves.png",
) -> None:
    """Plot ROC curves for multiple models or folds."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))

    colors = [PALETTE["primary"], PALETTE["secondary"], PALETTE["accent"],
              PALETTE["danger"], PALETTE["gray"]]

    for i, res in enumerate(results):
        fpr, tpr, _ = roc_curve(res["y_true"], res["y_prob"])
        auroc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=colors[i % len(colors)], lw=2,
                label=f"{res['label']} (AUROC={auroc:.3f})")

    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title("ROC Curves — CRISPR Off-Target Prediction", fontsize=13)
    ax.legend(loc="lower right", fontsize=10)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_pr_curves(
    results: List[Dict],
    save_path: str = "figures/pr_curves.png",
) -> None:
    """Plot Precision-Recall curves."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))

    colors = [PALETTE["primary"], PALETTE["secondary"], PALETTE["accent"],
              PALETTE["danger"], PALETTE["gray"]]

    for i, res in enumerate(results):
        precision, recall, _ = precision_recall_curve(res["y_true"], res["y_prob"])
        auprc = average_precision_score(res["y_true"], res["y_prob"])
        baseline = res["y_true"].mean()
        ax.plot(recall, precision, color=colors[i % len(colors)], lw=2,
                label=f"{res['label']} (AUPRC={auprc:.3f})")

    ax.axhline(baseline, color="k", linestyle="--", lw=1, alpha=0.5,
               label=f"Baseline (prevalence={baseline:.2f})")
    ax.set_xlabel("Recall", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_title("Precision-Recall Curves — CRISPR Off-Target Prediction", fontsize=13)
    ax.legend(loc="upper right", fontsize=10)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_cv_summary(
    fold_metrics: List[Dict],
    save_path: str = "figures/cv_summary.png",
) -> None:
    """Bar chart of per-fold AUROC and AUPRC with mean ± std."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    folds  = [f"Fold {m['fold']}" for m in fold_metrics]
    aurocs = [m["best_auroc"] for m in fold_metrics]
    auprcs = [m["auprc"] for m in fold_metrics]

    x = np.arange(len(folds))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    bars1 = ax.bar(x - width / 2, aurocs, width, label="AUROC",
                   color=PALETTE["primary"], alpha=0.85)
    bars2 = ax.bar(x + width / 2, auprcs, width, label="AUPRC",
                   color=PALETTE["secondary"], alpha=0.85)

    ax.axhline(np.mean(aurocs), color=PALETTE["primary"], linestyle="--", lw=1.5,
               label=f"Mean AUROC={np.mean(aurocs):.3f}±{np.std(aurocs):.3f}")
    ax.axhline(np.mean(auprcs), color=PALETTE["secondary"], linestyle="--", lw=1.5,
               label=f"Mean AUPRC={np.mean(auprcs):.3f}±{np.std(auprcs):.3f}")

    ax.set_xticks(x)
    ax.set_xticklabels(folds)
    ax.set_ylabel("Score", fontsize=12)
    ax.set_title("5-Fold Cross-Validation Performance", fontsize=13)
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_mismatch_importance(
    importances: np.ndarray,
    save_path: str = "figures/mismatch_importance.png",
) -> None:
    """Bar chart of per-position mismatch importance (from SHAP or attention)."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    n = len(importances)
    positions = [f"P{i+1}" for i in range(n)]
    colors = [PALETTE["danger"] if 8 <= i <= 19 else PALETTE["primary"] for i in range(n)]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(positions, importances, color=colors, alpha=0.85)
    ax.set_xlabel("Sequence Position (1=5'-end, 23=PAM)", fontsize=11)
    ax.set_ylabel("Mean |SHAP Value|", fontsize=11)
    ax.set_title("Per-Position Feature Importance (Seed Region Highlighted)", fontsize=12)
    # legend patches
    import matplotlib.patches as mpatches
    seed_patch = mpatches.Patch(color=PALETTE["danger"], label="Seed region (P9–P20)")
    non_patch  = mpatches.Patch(color=PALETTE["primary"], label="Non-seed")
    ax.legend(handles=[seed_patch, non_patch], fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
