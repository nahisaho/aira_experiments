"""
CRISPR-Cas9 Off-Target Prediction: Interpretability Module
Implements SHAP-based model interpretation for clinical applicability.
"""

import numpy as np
import torch
import shap
from typing import Dict, Optional
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from model import EpiCRISPRNet
from data_preprocessing import generate_synthetic_dataset, preprocess_dataset


FEATURE_NAMES = (
    # gRNA one-hot (4)
    ['gRNA_A', 'gRNA_C', 'gRNA_G', 'gRNA_T'] +
    # Target one-hot (4)
    ['Target_A', 'Target_C', 'Target_G', 'Target_T'] +
    # Mismatch types (16)
    [f'MM_{a}{b}' for a in 'ACGT' for b in 'ACGT'] +
    # Positional (3)
    ['Mismatch_Binary', 'PAM_Distance', 'Consec_Mismatch'] +
    # Epigenetic (4)
    ['Chromatin_Access', 'Methylation', 'H3K4me3', 'H3K27ac']
)

FEATURE_GROUPS = {
    'gRNA Sequence': list(range(0, 4)),
    'Target Sequence': list(range(4, 8)),
    'Mismatch Type': list(range(8, 24)),
    'Positional': list(range(24, 27)),
    'Epigenetic': list(range(27, 31))
}


class _WrapModel(torch.nn.Module):
    """Wrap model to output 2D tensor for SHAP compatibility."""
    def __init__(self, model):
        super().__init__()
        self.model = model
    def forward(self, x):
        return self.model(x).unsqueeze(-1)


def compute_shap_values(
    model: EpiCRISPRNet,
    X_background: np.ndarray,
    X_explain: np.ndarray,
    device: torch.device = None,
    n_background: int = 100
) -> np.ndarray:
    """
    Compute SHAP values using GradientExplainer (more compatible with attention).
    Returns shape: (n_explain, seq_len, n_features)
    """
    if device is None:
        device = torch.device('cpu')

    model.eval()
    model.to(device)
    wrapped = _WrapModel(model).to(device)

    bg = torch.FloatTensor(X_background[:n_background]).to(device)
    explain = torch.FloatTensor(X_explain).to(device)

    explainer = shap.GradientExplainer(wrapped, bg)
    shap_values = explainer.shap_values(explain)

    if isinstance(shap_values, list):
        shap_values = shap_values[0]

    return np.array(shap_values)


def compute_feature_group_importance(shap_values: np.ndarray) -> Dict[str, float]:
    """Aggregate SHAP values by feature group."""
    importance = {}
    for group_name, indices in FEATURE_GROUPS.items():
        group_shap = shap_values[:, :, indices]
        importance[group_name] = float(np.mean(np.abs(group_shap)))
    return importance


def compute_position_importance(shap_values: np.ndarray) -> np.ndarray:
    """Compute importance score per position (averaged over features and samples)."""
    return np.mean(np.abs(shap_values), axis=(0, 2))


def compute_shap_for_clinical(
    model: EpiCRISPRNet,
    X_data: np.ndarray,
    n_background: int = 100,
    n_explain: int = 50,
    device: torch.device = None
) -> Dict:
    """
    Full SHAP analysis pipeline for clinical interpretability.
    """
    if device is None:
        device = torch.device('cpu')

    shap_vals = compute_shap_values(
        model, X_data, X_data[:n_explain],
        device=device, n_background=n_background
    )

    group_importance = compute_feature_group_importance(shap_vals)
    position_importance = compute_position_importance(shap_vals)

    # Per-feature importance (averaged over positions and samples)
    feature_importance = {}
    for i, name in enumerate(FEATURE_NAMES):
        feature_importance[name] = float(np.mean(np.abs(shap_vals[:, :, i])))

    # Top features
    sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)

    return {
        'shap_values': shap_vals,
        'group_importance': group_importance,
        'position_importance': position_importance,
        'feature_importance': feature_importance,
        'top_features': sorted_features[:15],
        'feature_names': FEATURE_NAMES
    }


if __name__ == '__main__':
    print("Generating data for SHAP analysis...")
    data = generate_synthetic_dataset(n_samples=1000)
    X, y = preprocess_dataset(data)

    model = EpiCRISPRNet()

    print("Computing SHAP values...")
    results = compute_shap_for_clinical(
        model, X, n_background=50, n_explain=20
    )

    print("\nFeature Group Importance:")
    for group, imp in sorted(results['group_importance'].items(), key=lambda x: x[1], reverse=True):
        print(f"  {group}: {imp:.6f}")

    print("\nTop 10 Features:")
    for name, imp in results['top_features'][:10]:
        print(f"  {name}: {imp:.6f}")

    print("\nPosition Importance (first 10):")
    for i, imp in enumerate(results['position_importance'][:10]):
        print(f"  Position {i+1}: {imp:.6f}")
