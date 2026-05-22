"""
CRISPR-Cas9 Off-Target Prediction — CNN + Attention Model Architecture.

Architecture overview
─────────────────────
Input (sequence channel)
  └─ Conv1D block ×3 (local motif detection)
       └─ Multi-Head Self-Attention (position-aware context)
            └─ Global Average Pool
                 └─ Concat with scalar features (epigenetics, positional)
                      └─ MLP head → sigmoid (binary classification)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple


# ─── Building Blocks ──────────────────────────────────────────────────────────

class ConvBlock(nn.Module):
    """Conv1D → BatchNorm → GELU → Dropout."""

    def __init__(self, in_channels: int, out_channels: int,
                 kernel_size: int = 3, dropout: float = 0.1):
        super().__init__()
        self.conv = nn.Conv1d(
            in_channels, out_channels, kernel_size,
            padding=kernel_size // 2, bias=False,
        )
        self.bn   = nn.BatchNorm1d(out_channels)
        self.act  = nn.GELU()
        self.drop = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.drop(self.act(self.bn(self.conv(x))))


class MultiHeadSelfAttention(nn.Module):
    """
    Scaled dot-product multi-head self-attention over the sequence length axis.
    Input shape: (B, C, L)  →  Output: (B, C, L)
    """

    def __init__(self, embed_dim: int, num_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        assert embed_dim % num_heads == 0, "embed_dim must be divisible by num_heads"
        self.attn = nn.MultiheadAttention(
            embed_dim, num_heads, dropout=dropout, batch_first=True
        )
        self.norm = nn.LayerNorm(embed_dim)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        # x: (B, C, L) → transpose to (B, L, C) for MHA
        x_t = x.permute(0, 2, 1)
        attn_out, attn_weights = self.attn(x_t, x_t, x_t)
        out = self.norm(x_t + attn_out)          # residual + layer norm
        return out.permute(0, 2, 1), attn_weights  # back to (B, C, L)


class MLPHead(nn.Module):
    """Two-layer MLP classifier with residual dropout."""

    def __init__(self, in_dim: int, hidden_dim: int = 128, dropout: float = 0.2):
        super().__init__()
        self.fc1  = nn.Linear(in_dim, hidden_dim)
        self.act  = nn.GELU()
        self.drop = nn.Dropout(dropout)
        self.fc2  = nn.Linear(hidden_dim, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.drop(self.act(self.fc1(x)))
        return self.fc2(x).squeeze(-1)


# ─── Main Model ───────────────────────────────────────────────────────────────

class CRISPROffTargetModel(nn.Module):
    """
    CNN + Multi-Head Attention model for CRISPR off-target site prediction.

    Parameters
    ──────────
    seq_in_channels : int
        Number of input channels per position (default 23: 4+4+15).
    seq_length      : int
        Sequence length (default 23 for 20-nt guide + PAM).
    scalar_dim      : int
        Dimension of scalar features (positional mismatch + epigenetics).
    conv_channels   : list[int]
        Output channels for each Conv1D block.
    attn_heads      : int
        Number of attention heads.
    mlp_hidden      : int
        Hidden units in MLP head.
    dropout         : float
        Dropout probability applied throughout.
    """

    def __init__(
        self,
        seq_in_channels: int = 23,
        seq_length:      int = 23,
        scalar_dim:      int = 31,
        conv_channels:   Tuple[int, ...] = (64, 128, 256),
        attn_heads:      int = 4,
        mlp_hidden:      int = 128,
        dropout:         float = 0.2,
    ):
        super().__init__()
        self.seq_length = seq_length

        # ── Conv stack ────────────────────────────────────────────────────────
        conv_blocks = []
        in_ch = seq_in_channels
        for out_ch in conv_channels:
            conv_blocks.append(ConvBlock(in_ch, out_ch, kernel_size=3, dropout=dropout))
            in_ch = out_ch
        self.conv_stack = nn.Sequential(*conv_blocks)

        # ── Attention ─────────────────────────────────────────────────────────
        self.attention = MultiHeadSelfAttention(
            embed_dim=conv_channels[-1], num_heads=attn_heads, dropout=dropout
        )

        # ── Positional encoding (learnable) ──────────────────────────────────
        self.pos_emb = nn.Parameter(
            torch.zeros(1, conv_channels[-1], seq_length)
        )
        nn.init.trunc_normal_(self.pos_emb, std=0.02)

        # ── Global context pooling ────────────────────────────────────────────
        self.gap = nn.AdaptiveAvgPool1d(1)   # global average pool
        self.gmp = nn.AdaptiveMaxPool1d(1)   # global max  pool

        # ── Scalar feature encoder ────────────────────────────────────────────
        self.scalar_enc = nn.Sequential(
            nn.Linear(scalar_dim, 64),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(64, 64),
            nn.GELU(),
        )

        # ── Fusion + MLP head ─────────────────────────────────────────────────
        fused_dim = conv_channels[-1] * 2 + 64  # gap + gmp + scalar
        self.head = MLPHead(fused_dim, mlp_hidden, dropout)

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, nonlinearity="relu")
            elif isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(
        self,
        x_seq:    torch.Tensor,   # (B, L, C_in)
        x_scalar: torch.Tensor,   # (B, scalar_dim)
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Returns
        ───────
        logits      : (B,)  — raw pre-sigmoid scores
        attn_weights: (B, heads, L, L)  — for interpretability
        """
        # Conv expects (B, C, L)
        x = x_seq.permute(0, 2, 1)        # (B, C_in, L)
        x = self.conv_stack(x)             # (B, 256, L)
        x = x + self.pos_emb              # learnable positional encoding

        x, attn_weights = self.attention(x)  # (B, 256, L)

        gap_feat = self.gap(x).squeeze(-1)   # (B, 256)
        gmp_feat = self.gmp(x).squeeze(-1)   # (B, 256)
        scalar_feat = self.scalar_enc(x_scalar)  # (B, 64)

        fused  = torch.cat([gap_feat, gmp_feat, scalar_feat], dim=1)  # (B, 576)
        logits = self.head(fused)                                      # (B,)
        return logits, attn_weights

    @torch.no_grad()
    def predict_proba(
        self, x_seq: torch.Tensor, x_scalar: torch.Tensor
    ) -> torch.Tensor:
        """Return predicted probabilities (B,)."""
        self.eval()
        logits, _ = self(x_seq, x_scalar)
        return torch.sigmoid(logits)

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ─── Model Summary ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    model = CRISPROffTargetModel()
    B = 8
    x_seq    = torch.randn(B, 23, 23)
    x_scalar = torch.randn(B, 31)
    logits, attn = model(x_seq, x_scalar)
    print(f"Logits shape     : {logits.shape}")
    print(f"Attention shape  : {attn.shape}")
    print(f"Trainable params : {model.count_parameters():,}")
