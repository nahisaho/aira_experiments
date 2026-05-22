"""
Transfer learning framework for multi-device plasma disruption prediction.

Strategy: Domain-Adversarial Neural Network (DANN) + Fine-tuning protocol
- Source: JET (Joint European Torus) — large historical disruption database
- Target: ITER (extrapolation) / KSTAR (cross-device validation)

Transfer axes:
1. Feature normalisation: per-device IP / BT / minor-radius normalisation
2. Domain-adversarial training: learn device-invariant representations
3. Physics-guided alignment: Troyon/Greenwald fractions are device-agnostic
4. Progressive fine-tuning: frozen backbone → unfreeze TCN → full fine-tune
"""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset


# ─── Device normalisation parameters ─────────────────────────────────────────

@dataclass
class DeviceParams:
    """
    Physical parameters used to normalise signals to dimensionless units.
    Enables cross-device feature alignment.
    """
    name: str
    ip_max_ma: float    # Maximum plasma current [MA]
    bt_nom_t: float     # Nominal toroidal field [T]
    r_major_m: float    # Major radius [m]
    a_minor_m: float    # Minor radius [m]
    b_volume_m3: float  # Plasma volume [m³]
    n_greenwald_ref: float  # Reference Greenwald density [10^20 m^-3]
    # Normalisation statistics (set from data)
    signal_mean: Dict[str, float] = field(default_factory=dict)
    signal_std:  Dict[str, float] = field(default_factory=dict)


# Pre-configured device parameters
DEVICE_PARAMS: Dict[str, DeviceParams] = {
    "JET": DeviceParams(
        name="JET",
        ip_max_ma=4.0,
        bt_nom_t=3.45,
        r_major_m=2.96,
        a_minor_m=0.96,
        b_volume_m3=100.0,
        n_greenwald_ref=1.4,
    ),
    "KSTAR": DeviceParams(
        name="KSTAR",
        ip_max_ma=2.0,
        bt_nom_t=3.5,
        r_major_m=1.8,
        a_minor_m=0.5,
        b_volume_m3=18.0,
        n_greenwald_ref=2.5,
    ),
    "ITER": DeviceParams(
        name="ITER",
        ip_max_ma=15.0,
        bt_nom_t=5.3,
        r_major_m=6.2,
        a_minor_m=2.0,
        b_volume_m3=840.0,
        n_greenwald_ref=1.2,
    ),
    "ASDEX_UG": DeviceParams(
        name="ASDEX_UG",
        ip_max_ma=1.4,
        bt_nom_t=3.1,
        r_major_m=1.65,
        a_minor_m=0.5,
        b_volume_m3=14.0,
        n_greenwald_ref=1.8,
    ),
}


# ─── Normalisation layer ──────────────────────────────────────────────────────

class DeviceNormaliser(nn.Module):
    """
    Learnable per-device normalisation of input signals.
    Uses a combination of:
    1. Physics-based scaling (Troyon / Greenwald fractions)
    2. Learned affine transform (gamma, beta) per signal channel
    """

    def __init__(self, n_signals: int):
        super().__init__()
        # Learnable scale and offset (initialised to identity)
        self.gamma = nn.Parameter(torch.ones(n_signals))
        self.beta  = nn.Parameter(torch.zeros(n_signals))
        self.register_buffer("running_mean", torch.zeros(n_signals))
        self.register_buffer("running_std",  torch.ones(n_signals))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # Standardise with running statistics, then apply learnable affine
        x_norm = (x - self.running_mean) / (self.running_std + 1e-8)
        return self.gamma * x_norm + self.beta

    @torch.no_grad()
    def update_statistics(self, x: torch.Tensor, momentum: float = 0.1):
        """Update running mean/std from a batch."""
        mean = x.mean(dim=(0, 1))
        std  = x.std(dim=(0, 1))
        self.running_mean = (1 - momentum) * self.running_mean + momentum * mean
        self.running_std  = (1 - momentum) * self.running_std  + momentum * std


# ─── Gradient-reversal layer (for DANN) ──────────────────────────────────────

class GradientReversalFunction(torch.autograd.Function):
    """Reverses gradients during backpropagation (Ganin et al., 2016)."""

    @staticmethod
    def forward(ctx, x: torch.Tensor, alpha: float) -> torch.Tensor:
        ctx.alpha = alpha
        return x.clone()

    @staticmethod
    def backward(ctx, grad_output: torch.Tensor) -> Tuple[torch.Tensor, None]:
        return -ctx.alpha * grad_output, None


class GradientReversal(nn.Module):
    def __init__(self, alpha: float = 1.0):
        super().__init__()
        self.alpha = alpha

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return GradientReversalFunction.apply(x, self.alpha)


# ─── Domain classifier (DANN discriminator) ──────────────────────────────────

class DomainClassifier(nn.Module):
    """
    Classifies which device/domain a sample comes from.
    Coupled with gradient reversal to learn device-invariant features.
    """

    def __init__(self, in_features: int, n_devices: int = 4, alpha: float = 1.0):
        super().__init__()
        self.grl = GradientReversal(alpha)
        self.net = nn.Sequential(
            nn.Linear(in_features, 256), nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(256, 128),         nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(128, n_devices),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(self.grl(x))


# ─── Transfer-aware model wrapper ────────────────────────────────────────────

class TransferMHDModel(nn.Module):
    """
    Wraps HybridMHDPredictor with device-specific normalisation layers
    and a DANN domain classifier for adversarial domain alignment.
    """

    def __init__(
        self,
        backbone: nn.Module,
        device_names: List[str],
        n_signals: int = 32,
        domain_alpha: float = 1.0,
    ):
        super().__init__()
        self.backbone = backbone
        self.device_names = device_names
        self.n_devices = len(device_names)

        # Per-device normalisation layers
        self.normalisers = nn.ModuleDict({
            dev: DeviceNormaliser(n_signals) for dev in device_names
        })

        # Infer feature dimension from backbone
        fused_dim = (
            backbone.cfg.attn_d_model + backbone.cfg.physics_hidden
        )
        self.domain_clf = DomainClassifier(fused_dim, self.n_devices, domain_alpha)

    def forward(
        self,
        signals_seq: torch.Tensor,
        physics_feats: torch.Tensor,
        device_name: str,
    ) -> Dict[str, torch.Tensor]:

        # Apply device-specific normalisation
        if device_name in self.normalisers:
            norm = self.normalisers[device_name]
            signals_seq = norm(signals_seq)

        # Backbone forward pass
        out = self.backbone(signals_seq, physics_feats)

        # Domain adversarial prediction
        # Extract fused representation (re-use backbone internals via hook or direct call)
        # Here we approximate using physics_residuals projection input
        # In production, hook into the fusion layer output
        domain_logits = self.domain_clf(out["physics_residuals"])
        out["domain_logits"] = domain_logits

        return out


# ─── Transfer learning training protocol ─────────────────────────────────────

@dataclass
class TransferConfig:
    """Protocol for staged fine-tuning: JET pretrain → ITER/KSTAR fine-tune."""

    # Stage 0: Pre-train on JET
    pretrain_epochs: int = 100
    pretrain_lr: float = 1e-3

    # Stage 1: DANN domain alignment (JET + target data, unlabelled target ok)
    dann_epochs: int = 30
    dann_lr: float = 5e-4
    dann_alpha: float = 0.5  # GRL strength

    # Stage 2: Progressive unfreezing fine-tune on target device
    finetune_stages: List[Dict] = field(default_factory=lambda: [
        {"unfreeze": "normaliser",    "epochs": 10, "lr": 1e-3},
        {"unfreeze": "tcn_last2",     "epochs": 20, "lr": 5e-4},
        {"unfreeze": "full_backbone", "epochs": 30, "lr": 1e-4},
    ])

    # Regularisation
    l2_lambda: float = 1e-4
    mixup_alpha: float = 0.2  # Mixup between source and target batches


def progressive_unfreeze(
    model: TransferMHDModel,
    stage: str,
) -> None:
    """Freeze/unfreeze backbone layers according to transfer stage."""

    # First freeze everything
    for p in model.backbone.parameters():
        p.requires_grad = False

    if stage == "normaliser":
        # Only normaliser + task heads
        for norm in model.normalisers.values():
            for p in norm.parameters():
                p.requires_grad = True
        for head_name in ["cls_head", "ttd_head", "stability_head"]:
            for p in getattr(model.backbone, head_name).parameters():
                p.requires_grad = True

    elif stage == "tcn_last2":
        # Normaliser + task heads + last 2 TCN blocks
        progressive_unfreeze(model, "normaliser")
        tcn_blocks = list(model.backbone.tcn.net.children())
        for block in tcn_blocks[-2:]:
            for p in block.parameters():
                p.requires_grad = True
        for p in model.backbone.fusion.parameters():
            p.requires_grad = True

    elif stage == "full_backbone":
        # Unfreeze everything
        for p in model.backbone.parameters():
            p.requires_grad = True
        for p in model.normalisers.parameters():
            p.requires_grad = True

    n_trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[Transfer] Stage '{stage}': {n_trainable:,} trainable params")


# ─── Physics-guided domain alignment ─────────────────────────────────────────

def physics_guided_alignment_loss(
    src_physics: torch.Tensor,   # (B_src, n_phys)
    tgt_physics: torch.Tensor,   # (B_tgt, n_phys)
    physics_names: List[str],
) -> torch.Tensor:
    """
    Align distributions of device-agnostic physics features across devices.

    Dimensionless quantities (Troyon fraction, Greenwald fraction, q95/2, etc.)
    should have similar distributions across tokamaks at equivalent operating points.
    Uses Maximum Mean Discrepancy (MMD) with RBF kernel.
    """
    # Only align physics-invariant features
    invariant_indices = [
        i for i, n in enumerate(physics_names)
        if any(key in n for key in ["troyon", "greenwald", "q95_margin", "h98"])
    ]
    if not invariant_indices:
        return torch.tensor(0.0, device=src_physics.device)

    src = src_physics[:, invariant_indices]
    tgt = tgt_physics[:, invariant_indices]

    return _mmd_rbf(src, tgt)


def _mmd_rbf(
    x: torch.Tensor,
    y: torch.Tensor,
    bandwidth: float = 1.0,
) -> torch.Tensor:
    """Maximum Mean Discrepancy with RBF kernel."""
    xx = _rbf_kernel(x, x, bandwidth)
    yy = _rbf_kernel(y, y, bandwidth)
    xy = _rbf_kernel(x, y, bandwidth)
    return xx.mean() + yy.mean() - 2.0 * xy.mean()


def _rbf_kernel(
    a: torch.Tensor, b: torch.Tensor, bandwidth: float
) -> torch.Tensor:
    sq_dist = torch.cdist(a, b, p=2).pow(2)
    return torch.exp(-sq_dist / (2.0 * bandwidth ** 2))


# ─── ITER extrapolation scaling laws ─────────────────────────────────────────

def iter_scaling_correction(
    jet_features: np.ndarray,
    feature_names: List[str],
) -> np.ndarray:
    """
    Apply physics-based scaling corrections when extrapolating JET predictions to ITER.

    Key scaling differences:
    - ITER: R = 6.2 m vs JET: R = 2.96 m  → Alfvén time × √(R/a)
    - ITER: B = 5.3 T, I = 15 MA           → τE ∝ I^0.97 B^0.08 (ITER98)
    - ITER: D-T operation                   → dilution effects
    - ITER: larger machine → slower growth rates (τ_A scales as R)

    This is a simplified linear correction; a full simulation uses JINTRAC/TRANSP.
    """
    jet  = DEVICE_PARAMS["JET"]
    iter_ = DEVICE_PARAMS["ITER"]

    # Alfvén time scaling: τ_A ∝ R / v_A ∝ R √ρ / B
    alfven_scale = (iter_.r_major_m / jet.r_major_m) * (jet.bt_nom_t / iter_.bt_nom_t)

    corrected = features.copy() if hasattr(features := jet_features, 'copy') else np.copy(jet_features)

    # Scale time-derivative features by Alfvén time ratio
    deriv_indices = [i for i, n in enumerate(feature_names) if "_d1_" in n or "_d2_" in n]
    corrected[:, deriv_indices] /= alfven_scale

    # Scale frequency features by inverse Alfvén time
    freq_indices = [i for i, n in enumerate(feature_names) if "peak_freq" in n or "spectral_centroid" in n]
    corrected[:, freq_indices] *= alfven_scale

    return corrected


# ─── Checkpoint / serialisation ──────────────────────────────────────────────

def save_transfer_checkpoint(
    model: TransferMHDModel,
    config: TransferConfig,
    metrics: Dict[str, float],
    path: Path,
) -> None:
    torch.save({
        "model_state":  model.state_dict(),
        "transfer_cfg": asdict(config),
        "metrics":      metrics,
        "device_names": model.device_names,
    }, path)
    print(f"Checkpoint saved → {path}")


def load_transfer_checkpoint(
    model: TransferMHDModel,
    path: Path,
    strict: bool = True,
) -> Dict:
    ckpt = torch.load(path, map_location="cpu")
    model.load_state_dict(ckpt["model_state"], strict=strict)
    return ckpt["metrics"]
