"""
CRISPR-Cas9 Off-Target Prediction: Training & Evaluation Pipeline
Implements training loop, evaluation metrics, and cross-validation.
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    precision_recall_curve, roc_curve, f1_score,
    confusion_matrix, accuracy_score
)
from typing import Dict, List, Tuple
import time
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from model import EpiCRISPRNet, BaselineCNN, SequenceOnlyModel
from data_preprocessing import generate_synthetic_dataset, preprocess_dataset, create_cv_splits


class FocalLoss(nn.Module):
    """Focal loss for handling class imbalance in off-target prediction."""

    def __init__(self, alpha: float = 0.25, gamma: float = 2.0):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        bce = nn.functional.binary_cross_entropy_with_logits(logits, targets, reduction='none')
        pt = torch.exp(-bce)
        alpha_t = self.alpha * targets + (1 - self.alpha) * (1 - targets)
        focal = alpha_t * (1 - pt) ** self.gamma * bce
        return focal.mean()


def train_epoch(
    model: nn.Module,
    dataloader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device
) -> float:
    model.train()
    total_loss = 0
    for X_batch, y_batch in dataloader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()
        logits = model(X_batch)
        loss = criterion(logits, y_batch)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        total_loss += loss.item() * X_batch.size(0)
    return total_loss / len(dataloader.dataset)


def evaluate(
    model: nn.Module,
    dataloader: DataLoader,
    criterion: nn.Module,
    device: torch.device
) -> Dict:
    model.eval()
    all_logits, all_labels = [], []
    total_loss = 0

    with torch.no_grad():
        for X_batch, y_batch in dataloader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            logits = model(X_batch)
            loss = criterion(logits, y_batch)
            total_loss += loss.item() * X_batch.size(0)
            all_logits.append(logits.cpu().numpy())
            all_labels.append(y_batch.cpu().numpy())

    all_logits = np.concatenate(all_logits)
    all_labels = np.concatenate(all_labels)
    all_probs = 1 / (1 + np.exp(-all_logits))
    all_preds = (all_probs >= 0.5).astype(int)

    metrics = {
        'loss': total_loss / len(dataloader.dataset),
        'auroc': roc_auc_score(all_labels, all_probs),
        'auprc': average_precision_score(all_labels, all_probs),
        'accuracy': accuracy_score(all_labels, all_preds),
        'f1': f1_score(all_labels, all_preds, zero_division=0),
        'probs': all_probs,
        'labels': all_labels,
        'preds': all_preds
    }
    return metrics


def cross_validate(
    X: np.ndarray,
    y: np.ndarray,
    model_class,
    model_kwargs: Dict,
    n_folds: int = 5,
    n_epochs: int = 30,
    batch_size: int = 64,
    lr: float = 1e-3,
    device: torch.device = None,
    model_name: str = "Model"
) -> Dict:
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    splits = create_cv_splits(X, y, n_folds=n_folds)
    fold_results = []

    print(f"\n{'='*60}")
    print(f"Cross-validation: {model_name}")
    print(f"{'='*60}")

    for fold, (train_idx, val_idx) in enumerate(splits):
        print(f"\nFold {fold+1}/{n_folds}")

        X_train = torch.FloatTensor(X[train_idx])
        y_train = torch.FloatTensor(y[train_idx])
        X_val = torch.FloatTensor(X[val_idx])
        y_val = torch.FloatTensor(y[val_idx])

        train_loader = DataLoader(
            TensorDataset(X_train, y_train),
            batch_size=batch_size, shuffle=True
        )
        val_loader = DataLoader(
            TensorDataset(X_val, y_val),
            batch_size=batch_size
        )

        model = model_class(**model_kwargs).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)
        criterion = FocalLoss(alpha=0.25, gamma=2.0)

        best_auroc = 0
        best_metrics = None

        for epoch in range(n_epochs):
            train_loss = train_epoch(model, train_loader, optimizer, criterion, device)
            val_metrics = evaluate(model, val_loader, criterion, device)
            scheduler.step()

            if val_metrics['auroc'] > best_auroc:
                best_auroc = val_metrics['auroc']
                best_metrics = val_metrics

            if (epoch + 1) % 10 == 0:
                print(f"  Epoch {epoch+1}: train_loss={train_loss:.4f}, "
                      f"val_auroc={val_metrics['auroc']:.4f}, "
                      f"val_auprc={val_metrics['auprc']:.4f}")

        fold_results.append(best_metrics)
        print(f"  Best AUROC: {best_auroc:.4f}")

    # Aggregate results
    avg_metrics = {
        'auroc_mean': np.mean([r['auroc'] for r in fold_results]),
        'auroc_std': np.std([r['auroc'] for r in fold_results]),
        'auprc_mean': np.mean([r['auprc'] for r in fold_results]),
        'auprc_std': np.std([r['auprc'] for r in fold_results]),
        'f1_mean': np.mean([r['f1'] for r in fold_results]),
        'f1_std': np.std([r['f1'] for r in fold_results]),
        'accuracy_mean': np.mean([r['accuracy'] for r in fold_results]),
        'accuracy_std': np.std([r['accuracy'] for r in fold_results]),
    }

    print(f"\n--- {model_name} Summary ---")
    print(f"  AUROC: {avg_metrics['auroc_mean']:.4f} ± {avg_metrics['auroc_std']:.4f}")
    print(f"  AUPRC: {avg_metrics['auprc_mean']:.4f} ± {avg_metrics['auprc_std']:.4f}")
    print(f"  F1:    {avg_metrics['f1_mean']:.4f} ± {avg_metrics['f1_std']:.4f}")

    return {
        'fold_results': fold_results,
        'avg_metrics': avg_metrics,
        'model_name': model_name
    }


def run_full_experiment(
    n_samples: int = 8000,
    n_epochs: int = 30,
    n_folds: int = 5,
    batch_size: int = 64
) -> Dict:
    """Run the full comparative experiment."""
    print("=" * 60)
    print("CRISPR-Cas9 Off-Target Prediction Experiment")
    print("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Generate and preprocess data
    print("\n1. Generating synthetic dataset...")
    data = generate_synthetic_dataset(n_samples=n_samples)
    X, y = preprocess_dataset(data)
    print(f"   X: {X.shape}, y: {y.shape}")
    print(f"   Positive: {int(y.sum())}, Negative: {int(len(y) - y.sum())}")

    results = {}

    # Model 1: EpiCRISPR-Net (proposed)
    results['EpiCRISPR-Net'] = cross_validate(
        X, y,
        model_class=EpiCRISPRNet,
        model_kwargs={'input_channels': 31, 'seq_channels': 27, 'epi_channels': 4},
        n_folds=n_folds, n_epochs=n_epochs, batch_size=batch_size,
        device=device, model_name="EpiCRISPR-Net (Proposed)"
    )

    # Model 2: Baseline CNN
    results['BaselineCNN'] = cross_validate(
        X, y,
        model_class=BaselineCNN,
        model_kwargs={'input_channels': 31},
        n_folds=n_folds, n_epochs=n_epochs, batch_size=batch_size,
        device=device, model_name="Baseline CNN"
    )

    # Model 3: Sequence-only (ablation)
    results['SequenceOnly'] = cross_validate(
        X, y,
        model_class=SequenceOnlyModel,
        model_kwargs={'seq_channels': 27},
        n_folds=n_folds, n_epochs=n_epochs, batch_size=batch_size,
        device=device, model_name="Sequence-Only (Ablation)"
    )

    # Print comparison
    print("\n" + "=" * 60)
    print("COMPARISON TABLE")
    print("=" * 60)
    print(f"{'Model':<30} {'AUROC':<18} {'AUPRC':<18} {'F1':<18}")
    print("-" * 84)
    for name, res in results.items():
        m = res['avg_metrics']
        print(f"{res['model_name']:<30} "
              f"{m['auroc_mean']:.4f}±{m['auroc_std']:.4f}  "
              f"{m['auprc_mean']:.4f}±{m['auprc_std']:.4f}  "
              f"{m['f1_mean']:.4f}±{m['f1_std']:.4f}")

    return results


if __name__ == '__main__':
    results = run_full_experiment(n_samples=8000, n_epochs=30, n_folds=5)

    # Save results
    save_results = {}
    for name, res in results.items():
        save_results[name] = res['avg_metrics']

    os.makedirs('../results', exist_ok=True)
    with open('../results/experiment_results.json', 'w') as f:
        json.dump(save_results, f, indent=2)
    print("\nResults saved to results/experiment_results.json")
