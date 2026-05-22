"""
CRISPR-Cas9 Off-Target Prediction — Training Loop with Cross-Validation.
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, Subset
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from typing import Dict, List, Tuple
import logging
import json
from pathlib import Path

logger = logging.getLogger(__name__)


# ─── Dataset ──────────────────────────────────────────────────────────────────

class OffTargetDataset(Dataset):
    def __init__(
        self,
        X_seq:    np.ndarray,   # (N, L, C)
        X_scalar: np.ndarray,   # (N, D)
        y:        np.ndarray,   # (N,)
    ):
        self.X_seq    = torch.from_numpy(X_seq).float()
        self.X_scalar = torch.from_numpy(X_scalar).float()
        self.y        = torch.from_numpy(y).float()

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X_seq[idx], self.X_scalar[idx], self.y[idx]


# ─── Training Utilities ───────────────────────────────────────────────────────

def focal_bce_loss(
    logits: torch.Tensor, targets: torch.Tensor,
    gamma: float = 2.0, pos_weight: float = 3.0
) -> torch.Tensor:
    """
    Focal binary cross-entropy for class imbalance.
    Downweights easy examples (pt → 1) and focuses on hard negatives.
    """
    pw = torch.tensor([pos_weight], device=logits.device)
    bce = nn.functional.binary_cross_entropy_with_logits(
        logits, targets, pos_weight=pw, reduction="none"
    )
    probs = torch.sigmoid(logits)
    pt = torch.where(targets == 1, probs, 1 - probs)
    focal_weight = (1 - pt) ** gamma
    return (focal_weight * bce).mean()


def train_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    pos_weight: float = 3.0,
) -> float:
    model.train()
    total_loss = 0.0
    for x_seq, x_scalar, y in loader:
        x_seq    = x_seq.to(device)
        x_scalar = x_scalar.to(device)
        y        = y.to(device)

        optimizer.zero_grad()
        logits, _ = model(x_seq, x_scalar)
        loss = focal_bce_loss(logits, y, pos_weight=pos_weight)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)


@torch.no_grad()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
) -> Dict[str, float]:
    model.eval()
    all_probs, all_labels = [], []
    for x_seq, x_scalar, y in loader:
        logits, _ = model(x_seq.to(device), x_scalar.to(device))
        probs = torch.sigmoid(logits).cpu().numpy()
        all_probs.extend(probs)
        all_labels.extend(y.numpy())

    probs  = np.array(all_probs)
    labels = np.array(all_labels)
    auroc  = roc_auc_score(labels, probs) if len(np.unique(labels)) > 1 else 0.5
    auprc  = average_precision_score(labels, probs) if len(np.unique(labels)) > 1 else 0.0
    preds  = (probs >= 0.5).astype(int)
    tp = ((preds == 1) & (labels == 1)).sum()
    fp = ((preds == 1) & (labels == 0)).sum()
    fn = ((preds == 0) & (labels == 1)).sum()
    precision = tp / (tp + fp + 1e-8)
    recall    = tp / (tp + fn + 1e-8)
    f1        = 2 * precision * recall / (precision + recall + 1e-8)
    return {
        "auroc":     float(auroc),
        "auprc":     float(auprc),
        "precision": float(precision),
        "recall":    float(recall),
        "f1":        float(f1),
    }


# ─── Cross-Validation Trainer ─────────────────────────────────────────────────

def cross_validate(
    X_seq:      np.ndarray,
    X_scalar:   np.ndarray,
    y:          np.ndarray,
    n_splits:   int = 5,
    epochs:     int = 30,
    batch_size: int = 64,
    lr:         float = 3e-4,
    seed:       int = 42,
    save_dir:   str = "results",
    device_str: str = "cpu",
) -> List[Dict[str, float]]:
    """
    Stratified k-fold cross-validation.
    Saves per-fold metrics and the best model checkpoint.
    """
    np.random.seed(seed)
    torch.manual_seed(seed)
    device    = torch.device(device_str)
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)

    dataset = OffTargetDataset(X_seq, X_scalar, y)
    skf     = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    fold_metrics = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(X_seq, y)):
        logger.info("─── Fold %d/%d ───", fold + 1, n_splits)

        from src.model import CRISPROffTargetModel  # local import for clean reload
        model = CRISPROffTargetModel(
            seq_in_channels=X_seq.shape[2],
            seq_length=X_seq.shape[1],
            scalar_dim=X_scalar.shape[1],
        ).to(device)

        optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

        train_loader = DataLoader(
            Subset(dataset, train_idx), batch_size=batch_size, shuffle=True
        )
        val_loader = DataLoader(
            Subset(dataset, val_idx), batch_size=batch_size
        )

        best_auroc  = 0.0
        best_state  = None
        history: List[Dict] = []

        pos_weight = (y[train_idx] == 0).sum() / max((y[train_idx] == 1).sum(), 1)

        for epoch in range(epochs):
            train_loss = train_epoch(model, train_loader, optimizer, device, pos_weight)
            metrics    = evaluate(model, val_loader, device)
            metrics["train_loss"] = train_loss
            metrics["epoch"]      = epoch + 1
            history.append(metrics)
            scheduler.step()

            if metrics["auroc"] > best_auroc:
                best_auroc = metrics["auroc"]
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

            if (epoch + 1) % 10 == 0:
                logger.info(
                    "  Epoch %d: loss=%.4f AUROC=%.4f AUPRC=%.4f",
                    epoch + 1, train_loss, metrics["auroc"], metrics["auprc"],
                )

        # Save best model for this fold
        torch.save(best_state, save_path / f"model_fold{fold+1}.pt")

        final = history[-1].copy()
        final["best_auroc"] = best_auroc
        final["fold"]       = fold + 1
        fold_metrics.append(final)
        logger.info("  Best AUROC: %.4f", best_auroc)

    # Summary
    aurocs = [m["best_auroc"] for m in fold_metrics]
    logger.info(
        "CV Summary — AUROC: %.4f ± %.4f",
        np.mean(aurocs), np.std(aurocs),
    )

    summary = {
        "fold_metrics": fold_metrics,
        "mean_auroc":   float(np.mean(aurocs)),
        "std_auroc":    float(np.std(aurocs)),
        "mean_auprc":   float(np.mean([m["auprc"] for m in fold_metrics])),
    }
    with open(save_path / "cv_results.json", "w") as f:
        json.dump(summary, f, indent=2)

    return fold_metrics
