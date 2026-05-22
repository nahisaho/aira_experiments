"""
CRISPR-Cas9 Off-Target Prediction — SHAP-based Interpretability.

Provides:
  1. KernelSHAP wrapper for the CNN+Attention model.
  2. Per-position, per-feature, and per-mismatch-type SHAP analysis.
  3. Summary beeswarm + heatmap plots.
"""

import numpy as np
import torch
import shap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple, Optional
import json
import logging

logger = logging.getLogger(__name__)

# ─── Constants (must match preprocessing.py) ──────────────────────────────────
SEQ_LEN         = 23
SEQ_CHANNELS    = 23   # 4 (guide OH) + 4 (target OH) + 15 (mismatch type)
SCALAR_DIM      = 31   # 23 positional + 8 epigenetic

FEATURE_NAMES_POS = [f"Pos{i+1}_mismatch" for i in range(SEQ_LEN)]
FEATURE_NAMES_EPI = [
    "ATAC_min", "ATAC_p33", "ATAC_p66", "ATAC_max",
    "Meth_mean", "Meth_std", "Meth_hyper_frac", "Meth_unmeth_frac",
]
SCALAR_FEATURE_NAMES = FEATURE_NAMES_POS + FEATURE_NAMES_EPI


# ─── Model Wrapper for SHAP ───────────────────────────────────────────────────

class FlatModelWrapper:
    """
    Wraps CRISPROffTargetModel for SHAP KernelExplainer.
    Input: flat (N, SEQ_LEN*SEQ_CHANNELS + SCALAR_DIM) array.
    Output: predicted probability (N,).
    """

    def __init__(self, model: torch.nn.Module, device: str = "cpu"):
        self.model  = model
        self.device = torch.device(device)
        self.model.eval()
        self.seq_flat = SEQ_LEN * SEQ_CHANNELS  # 529

    def __call__(self, X_flat: np.ndarray) -> np.ndarray:
        X_flat  = np.atleast_2d(X_flat).astype(np.float32)
        x_seq   = X_flat[:, : self.seq_flat].reshape(-1, SEQ_LEN, SEQ_CHANNELS)
        x_scal  = X_flat[:, self.seq_flat :]

        x_seq_t  = torch.from_numpy(x_seq).to(self.device)
        x_scal_t = torch.from_numpy(x_scal).to(self.device)

        with torch.no_grad():
            logits, _ = self.model(x_seq_t, x_scal_t)
            probs = torch.sigmoid(logits).cpu().numpy()
        return probs


def flatten_inputs(x_seq: np.ndarray, x_scalar: np.ndarray) -> np.ndarray:
    """Flatten (N, L, C) + (N, D) → (N, L*C + D)."""
    N = x_seq.shape[0]
    return np.concatenate([x_seq.reshape(N, -1), x_scalar], axis=1)


# ─── SHAP Explainer ───────────────────────────────────────────────────────────

class CRISPRSHAPExplainer:
    """
    KernelSHAP explainer for off-target prediction model.
    Uses a background (reference) set to compute expected model output.
    """

    def __init__(
        self,
        model: torch.nn.Module,
        background_seq:    np.ndarray,
        background_scalar: np.ndarray,
        n_background:      int = 100,
        device:            str = "cpu",
        seed:              int = 42,
    ):
        np.random.seed(seed)
        self.wrapper  = FlatModelWrapper(model, device)
        self.seq_flat = SEQ_LEN * SEQ_CHANNELS

        # Sample background set
        N = background_seq.shape[0]
        idx = np.random.choice(N, min(n_background, N), replace=False)
        bg_flat = flatten_inputs(background_seq[idx], background_scalar[idx])

        logger.info("Initialising KernelSHAP with %d background samples…", len(idx))
        self.explainer = shap.KernelExplainer(self.wrapper, bg_flat)

    def explain(
        self,
        x_seq:    np.ndarray,   # (N, L, C)
        x_scalar: np.ndarray,   # (N, D)
        n_samples: int = 50,
    ) -> np.ndarray:
        """
        Compute SHAP values for N samples.
        Returns (N, SEQ_LEN*SEQ_CHANNELS + SCALAR_DIM) array.
        """
        X_flat = flatten_inputs(x_seq, x_scalar)
        logger.info("Computing SHAP values for %d samples (nsamples=%d)…",
                    len(X_flat), n_samples)
        shap_vals = self.explainer.shap_values(X_flat, nsamples=n_samples, silent=True)
        return np.array(shap_vals)

    def aggregate_positional(self, shap_vals: np.ndarray) -> np.ndarray:
        """
        Aggregate SHAP values over the sequence channel dimension → (N, SEQ_LEN).
        Uses sum of |SHAP| per position.
        """
        seq_shap = shap_vals[:, : self.seq_flat].reshape(-1, SEQ_LEN, SEQ_CHANNELS)
        return np.abs(seq_shap).sum(axis=2)   # (N, SEQ_LEN)

    def aggregate_scalar(self, shap_vals: np.ndarray) -> np.ndarray:
        """Return SHAP values for scalar features → (N, SCALAR_DIM)."""
        return shap_vals[:, self.seq_flat :]


# ─── Attention Map Extraction ─────────────────────────────────────────────────

@torch.no_grad()
def extract_attention_maps(
    model:    torch.nn.Module,
    x_seq:    torch.Tensor,   # (B, L, C)
    x_scalar: torch.Tensor,   # (B, D)
) -> np.ndarray:
    """
    Extract attention weight maps from model.
    Returns (B, num_heads, L, L) numpy array.
    """
    model.eval()
    _, attn_weights = model(x_seq, x_scalar)
    return attn_weights.cpu().numpy()


# ─── Visualisations ───────────────────────────────────────────────────────────

def plot_shap_beeswarm(
    shap_vals:    np.ndarray,   # (N, SCALAR_DIM)
    feature_vals: np.ndarray,   # (N, SCALAR_DIM)
    feature_names: list,
    save_path: str = "figures/shap_beeswarm.png",
    max_display: int = 20,
) -> None:
    """SHAP beeswarm plot for scalar features."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)

    # Sort by mean |SHAP|
    mean_abs = np.abs(shap_vals).mean(axis=0)
    top_idx  = np.argsort(mean_abs)[::-1][:max_display]
    sv_top   = shap_vals[:, top_idx]
    fv_top   = feature_vals[:, top_idx]
    names_top = [feature_names[i] for i in top_idx]

    shap.summary_plot(
        sv_top, fv_top,
        feature_names=names_top,
        show=False, plot_type="dot",
    )
    plt.title("SHAP Feature Importance — Off-Target Prediction", fontsize=12)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("SHAP beeswarm saved: %s", save_path)


def plot_positional_shap_heatmap(
    pos_shap: np.ndarray,   # (N, SEQ_LEN)
    labels:   np.ndarray,   # (N,) binary
    save_path: str = "figures/positional_shap_heatmap.png",
) -> None:
    """
    Heatmap of mean |SHAP| per position, split by positive / negative label.
    """
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    pos_mean = np.abs(pos_shap[labels == 1]).mean(axis=0)
    neg_mean = np.abs(pos_shap[labels == 0]).mean(axis=0)

    data   = np.vstack([pos_mean, neg_mean])
    yticks = ["Off-target (+)", "Non-target (-)"]
    xticks = [f"P{i+1}" for i in range(SEQ_LEN)]

    fig, ax = plt.subplots(figsize=(14, 3))
    sns.heatmap(
        data, ax=ax, cmap="viridis",
        xticklabels=xticks, yticklabels=yticks,
        linewidths=0.3, linecolor="white",
        cbar_kws={"label": "Mean |SHAP|"},
    )
    ax.set_title("Positional SHAP Heatmap by Label", fontsize=12)
    ax.tick_params(axis="x", labelsize=8)
    # Highlight seed region
    for i in range(8, 20):
        ax.add_patch(plt.Rectangle((i, 0), 1, 2, fill=False,
                                   edgecolor="red", lw=1.5))
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Positional SHAP heatmap saved: %s", save_path)


def plot_attention_heatmap(
    attn_weights: np.ndarray,  # (B, L, L) or (B, heads, L, L)
    sample_idx: int = 0,
    save_path: str = "figures/attention_heatmap.png",
) -> None:
    """Plot mean-head attention weight matrix for a single sample."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    a = attn_weights[sample_idx]
    # If heads dimension is present (ndim==3), average over heads
    mean_attn = a.mean(axis=0) if a.ndim == 3 else a   # → (L, L)

    pos_labels = [f"P{i+1}" for i in range(mean_attn.shape[0])]
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(
        mean_attn, ax=ax, cmap="Blues",
        xticklabels=pos_labels, yticklabels=pos_labels,
        cbar_kws={"label": "Attention Weight"},
        linewidths=0.2,
    )
    ax.set_xlabel("Key Position", fontsize=10)
    ax.set_ylabel("Query Position", fontsize=10)
    ax.set_title(f"Mean-Head Self-Attention Map (Sample {sample_idx})", fontsize=11)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Attention heatmap saved: %s", save_path)


def save_shap_summary(
    shap_vals:    np.ndarray,
    feature_names: list,
    save_path: str = "results/shap_summary.json",
) -> None:
    """Save mean |SHAP| per feature to JSON."""
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    mean_abs = np.abs(shap_vals).mean(axis=0)
    summary  = {name: float(v) for name, v in zip(feature_names, mean_abs)}
    summary_sorted = dict(sorted(summary.items(), key=lambda x: -x[1]))
    with open(save_path, "w") as f:
        json.dump(summary_sorted, f, indent=2)
    logger.info("SHAP summary saved: %s", save_path)
