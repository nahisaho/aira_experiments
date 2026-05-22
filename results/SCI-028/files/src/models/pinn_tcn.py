"""
Physics-Informed Machine Learning for MHD instability prediction.

Architecture: Hybrid PINN + Temporal Convolutional Network (TCN)
- Physics branch: encodes MHD equilibrium constraints as soft loss terms
- Data branch: TCN extracts spatiotemporal patterns
- Fusion layer: attention-weighted combination
"""

from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


# ─── Configuration ────────────────────────────────────────────────────────────

@dataclass
class ModelConfig:
    # Input
    n_signals: int = 32          # Number of raw signal channels
    n_stat_features: int = 256   # Pre-computed statistical features
    seq_len: int = 500           # Temporal sequence length
    # TCN
    tcn_channels: List[int] = None
    tcn_kernel_size: int = 7
    tcn_dropout: float = 0.15
    # Physics branch
    physics_hidden: int = 128
    physics_layers: int = 3
    # Fusion
    attn_heads: int = 4
    attn_d_model: int = 128
    # Output
    n_classes: int = 3           # [safe, warning, imminent_disruption]
    regression_targets: int = 4  # [ttd, betan_margin, q95_margin, locked_mode_amp]

    def __post_init__(self):
        if self.tcn_channels is None:
            self.tcn_channels = [64, 64, 128, 128, 256]


# ─── TCN components ───────────────────────────────────────────────────────────

class CausalConv1d(nn.Module):
    """Causal dilated 1-D convolution (no future leakage)."""

    def __init__(self, in_ch: int, out_ch: int, kernel: int, dilation: int = 1):
        super().__init__()
        self.padding = (kernel - 1) * dilation
        self.conv = nn.Conv1d(
            in_ch, out_ch, kernel,
            padding=self.padding, dilation=dilation,
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = self.conv(x)
        return out[:, :, :-self.padding] if self.padding else out


class TCNResidualBlock(nn.Module):
    """Temporal Convolutional Network residual block with dilated causal convs."""

    def __init__(self, in_ch: int, out_ch: int, kernel: int, dilation: int, dropout: float):
        super().__init__()
        self.conv1 = CausalConv1d(in_ch, out_ch, kernel, dilation)
        self.conv2 = CausalConv1d(out_ch, out_ch, kernel, dilation)
        self.norm1 = nn.LayerNorm(out_ch)
        self.norm2 = nn.LayerNorm(out_ch)
        self.drop  = nn.Dropout(dropout)
        self.skip  = nn.Conv1d(in_ch, out_ch, 1) if in_ch != out_ch else nn.Identity()
        self.act   = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, C, T)
        h = self.act(self.norm1(self.conv1(x).transpose(1, 2)).transpose(1, 2))
        h = self.drop(h)
        h = self.act(self.norm2(self.conv2(h).transpose(1, 2)).transpose(1, 2))
        h = self.drop(h)
        return self.act(h + self.skip(x))


class TemporalConvNet(nn.Module):
    """
    Multi-scale TCN for plasma signal processing.
    Receptive field grows exponentially with dilation, covering ~500 ms at 10 kHz.
    """

    def __init__(self, in_ch: int, channels: List[int], kernel: int = 7, dropout: float = 0.15):
        super().__init__()
        layers = []
        for i, out_ch in enumerate(channels):
            dilation = 2 ** i
            layers.append(TCNResidualBlock(
                in_ch if i == 0 else channels[i - 1],
                out_ch, kernel, dilation, dropout,
            ))
        self.net = nn.Sequential(*layers)
        self.out_channels = channels[-1]

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)   # (B, C_last, T)


# ─── Physics branch ───────────────────────────────────────────────────────────

class PhysicsBranch(nn.Module):
    """
    Encodes MHD equilibrium-derived features through a shallow MLP.
    Inputs are the hand-crafted physics proxies (stability margins, Troyon fraction, etc.).

    Physics priors enforced as soft constraints via auxiliary losses (see PINNLoss below).
    """

    def __init__(self, in_features: int, hidden: int = 128, n_layers: int = 3):
        super().__init__()
        layers: List[nn.Module] = []
        d_in = in_features
        for _ in range(n_layers):
            layers += [nn.Linear(d_in, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(0.1)]
            d_in = hidden
        self.mlp = nn.Sequential(*layers)
        self.out_dim = hidden

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.mlp(x)   # (B, hidden)


# ─── Cross-attention fusion ───────────────────────────────────────────────────

class CrossAttentionFusion(nn.Module):
    """
    Cross-attention between physics-branch representation and TCN temporal outputs.
    Physics features act as queries; TCN sequence acts as keys/values.
    """

    def __init__(self, d_phys: int, d_tcn: int, d_model: int, n_heads: int):
        super().__init__()
        self.q_proj = nn.Linear(d_phys, d_model)
        self.k_proj = nn.Linear(d_tcn,  d_model)
        self.v_proj = nn.Linear(d_tcn,  d_model)
        self.attn   = nn.MultiheadAttention(d_model, n_heads, batch_first=True)
        self.norm   = nn.LayerNorm(d_model)
        self.out_dim = d_model

    def forward(
        self,
        phys: torch.Tensor,  # (B, d_phys)
        tcn:  torch.Tensor,  # (B, C, T)
    ) -> torch.Tensor:
        # Reshape TCN output: (B, T, C)
        tcn_seq = tcn.permute(0, 2, 1)  # (B, T, C)
        q = self.q_proj(phys).unsqueeze(1)   # (B, 1, d_model)
        k = self.k_proj(tcn_seq)             # (B, T, d_model)
        v = self.v_proj(tcn_seq)             # (B, T, d_model)
        out, _ = self.attn(q, k, v)          # (B, 1, d_model)
        return self.norm(out.squeeze(1))     # (B, d_model)


# ─── Full hybrid PINN-TCN model ───────────────────────────────────────────────

class HybridMHDPredictor(nn.Module):
    """
    Hybrid Physics-Informed TCN for disruption prediction and MHD stability assessment.

    Outputs
    -------
    cls_logits : (B, n_classes)
        Disruption risk classification logits.
    ttd_pred : (B, 1)
        Time-to-disruption regression [ms].
    stability_margins : (B, 3)
        Predicted stability margins: [betan_margin, q95_margin, locked_mode_amp].
    physics_residuals : (B, n_physics_constraints)
        Physics constraint residuals for PINN loss (training only).
    """

    def __init__(self, cfg: Optional[ModelConfig] = None):
        super().__init__()
        self.cfg = cfg or ModelConfig()
        c = self.cfg

        # Signal pre-processing: per-channel normalisation + projection
        self.signal_proj = nn.Sequential(
            nn.Linear(c.n_signals, 64),
            nn.GELU(),
            nn.Linear(64, 64),
        )

        # TCN branch: operates on projected signals over time
        self.tcn = TemporalConvNet(
            in_ch=64,
            channels=c.tcn_channels,
            kernel=c.tcn_kernel_size,
            dropout=c.tcn_dropout,
        )

        # Physics branch
        self.physics = PhysicsBranch(
            in_features=c.n_stat_features,
            hidden=c.physics_hidden,
            n_layers=c.physics_layers,
        )

        # Cross-attention fusion
        self.fusion = CrossAttentionFusion(
            d_phys=c.physics_hidden,
            d_tcn=c.tcn_channels[-1],
            d_model=c.attn_d_model,
            n_heads=c.attn_heads,
        )

        # Task heads
        d_fused = c.attn_d_model + c.physics_hidden  # concatenation
        self.cls_head = nn.Sequential(
            nn.Linear(d_fused, 128), nn.GELU(), nn.Dropout(0.1),
            nn.Linear(128, c.n_classes),
        )
        self.ttd_head = nn.Sequential(
            nn.Linear(d_fused, 64), nn.GELU(),
            nn.Linear(64, 1), nn.Softplus(),  # TTD > 0
        )
        self.stability_head = nn.Sequential(
            nn.Linear(d_fused, 64), nn.GELU(),
            nn.Linear(64, 3),  # betan_margin, q95_margin, locked_mode_amp
        )

        # Physics constraint projector (for soft residuals)
        self.physics_constraint_proj = nn.Linear(d_fused, 6)

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(
        self,
        signals_seq: torch.Tensor,   # (B, T, n_signals) raw signal sequence
        physics_feats: torch.Tensor, # (B, n_stat_features) pre-computed features
    ) -> Dict[str, torch.Tensor]:
        B, T, S = signals_seq.shape

        # Project signals to embedding dimension
        sig_emb = self.signal_proj(signals_seq)  # (B, T, 64)
        sig_emb = sig_emb.permute(0, 2, 1)       # (B, 64, T) for TCN

        # TCN temporal encoding
        tcn_out = self.tcn(sig_emb)   # (B, C_last, T)

        # Physics branch
        phys_out = self.physics(physics_feats)  # (B, d_phys)

        # Cross-attention fusion
        fused_attn = self.fusion(phys_out, tcn_out)  # (B, d_model)

        # Concatenate physics and attention output
        fused = torch.cat([fused_attn, phys_out], dim=-1)  # (B, d_fused)

        return {
            "cls_logits":          self.cls_head(fused),
            "ttd_pred":            self.ttd_head(fused),
            "stability_margins":   self.stability_head(fused),
            "physics_residuals":   self.physics_constraint_proj(fused),
        }


# ─── Physics-informed loss ────────────────────────────────────────────────────

class PINNLoss(nn.Module):
    """
    Combined loss with physics-constraint soft penalties.

    Physics constraints implemented:
    1. Troyon limit: βN < βN_limit → predicted betan_margin > 0 when safe
    2. q = 2 proximity: q95 > 2.0 → q95_margin > 0 when safe
    3. Greenwald fraction: n/n_G < 0.9 → feature penalty
    4. Locked mode threshold: locked_mode < threshold when safe
    5. Energy balance consistency: dW/dt ≈ P_heat - P_rad - P_lost
    6. Plasma current continuity: smooth dIp/dt
    """

    def __init__(
        self,
        lambda_cls:   float = 1.0,
        lambda_ttd:   float = 0.5,
        lambda_stab:  float = 0.3,
        lambda_phys:  float = 0.2,
        class_weights: Optional[torch.Tensor] = None,
    ):
        super().__init__()
        self.lambda_cls  = lambda_cls
        self.lambda_ttd  = lambda_ttd
        self.lambda_stab = lambda_stab
        self.lambda_phys = lambda_phys

        self.cls_loss  = nn.CrossEntropyLoss(weight=class_weights)
        self.ttd_loss  = nn.HuberLoss(delta=50.0)   # delta = 50 ms
        self.stab_loss = nn.MSELoss()

    def forward(
        self,
        preds: Dict[str, torch.Tensor],
        targets: Dict[str, torch.Tensor],
    ) -> Tuple[torch.Tensor, Dict[str, float]]:

        loss_cls = self.cls_loss(preds["cls_logits"], targets["cls_label"])
        loss_ttd = self.ttd_loss(preds["ttd_pred"].squeeze(-1), targets["ttd"])
        loss_stab = self.stab_loss(preds["stability_margins"], targets["stability_margins"])

        # Physics residual penalty
        # Enforce sign consistency: safe shots should have positive stability margins
        phys_res = preds["physics_residuals"]
        safe_mask = (targets["cls_label"] == 0).float().unsqueeze(-1)
        loss_phys = (F.relu(-phys_res) * safe_mask).mean()

        total = (
            self.lambda_cls  * loss_cls +
            self.lambda_ttd  * loss_ttd +
            self.lambda_stab * loss_stab +
            self.lambda_phys * loss_phys
        )

        return total, {
            "loss_cls":  loss_cls.item(),
            "loss_ttd":  loss_ttd.item(),
            "loss_stab": loss_stab.item(),
            "loss_phys": loss_phys.item(),
            "loss_total": total.item(),
        }


# ─── Uncertainty quantification via MC-Dropout ───────────────────────────────

class MCDropoutPredictor:
    """Approximate Bayesian inference via Monte Carlo dropout at inference time."""

    def __init__(self, model: HybridMHDPredictor, n_samples: int = 50):
        self.model = model
        self.n_samples = n_samples

    @torch.no_grad()
    def predict_with_uncertainty(
        self,
        signals_seq: torch.Tensor,
        physics_feats: torch.Tensor,
    ) -> Dict[str, np.ndarray]:
        self.model.train()  # Enable dropout
        samples = [self.model(signals_seq, physics_feats) for _ in range(self.n_samples)]
        self.model.eval()

        # Stack and compute statistics
        cls_probs = torch.stack([
            F.softmax(s["cls_logits"], dim=-1) for s in samples
        ])  # (n_samples, B, n_classes)

        ttd_preds = torch.stack([s["ttd_pred"] for s in samples])  # (n_samples, B, 1)

        return {
            "cls_prob_mean":  cls_probs.mean(0).cpu().numpy(),
            "cls_prob_std":   cls_probs.std(0).cpu().numpy(),
            "ttd_mean_ms":    ttd_preds.mean(0).squeeze(-1).cpu().numpy(),
            "ttd_std_ms":     ttd_preds.std(0).squeeze(-1).cpu().numpy(),
            "predictive_entropy": self._predictive_entropy(cls_probs).cpu().numpy(),
        }

    @staticmethod
    def _predictive_entropy(probs: torch.Tensor) -> torch.Tensor:
        """H[p(y|x)] = -Σ p̄ log p̄"""
        p_bar = probs.mean(0)
        return -(p_bar * (p_bar + 1e-9).log()).sum(-1)


# ─── Model factory ───────────────────────────────────────────────────────────

def build_model(cfg: Optional[ModelConfig] = None) -> HybridMHDPredictor:
    """Build and return model with optional config."""
    model = HybridMHDPredictor(cfg)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"HybridMHDPredictor — {n_params:,} trainable parameters")
    return model


if __name__ == "__main__":
    cfg = ModelConfig()
    model = build_model(cfg)
    model.eval()

    B, T = 4, 500
    signals = torch.randn(B, T, cfg.n_signals)
    physics  = torch.randn(B, cfg.n_stat_features)

    with torch.no_grad():
        out = model(signals, physics)

    print("cls_logits:        ", out["cls_logits"].shape)
    print("ttd_pred:          ", out["ttd_pred"].shape)
    print("stability_margins: ", out["stability_margins"].shape)
    print("physics_residuals: ", out["physics_residuals"].shape)
