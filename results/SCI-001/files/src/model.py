"""
CRISPR-Cas9 Off-Target Prediction Model: CNN + Attention Architecture
EpiCRISPR-Net: Epigenetics-integrated CNN-Attention model for off-target prediction.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import math
from typing import Optional, Tuple


class MultiScaleCNNBlock(nn.Module):
    """Multi-scale CNN block with parallel convolutions of different kernel sizes."""

    def __init__(self, in_channels: int, out_channels: int):
        super().__init__()
        assert out_channels % 3 == 0, "out_channels must be divisible by 3"
        branch_ch = out_channels // 3

        self.conv3 = nn.Sequential(
            nn.Conv1d(in_channels, branch_ch, kernel_size=3, padding=1),
            nn.BatchNorm1d(branch_ch),
            nn.GELU()
        )
        self.conv5 = nn.Sequential(
            nn.Conv1d(in_channels, branch_ch, kernel_size=5, padding=2),
            nn.BatchNorm1d(branch_ch),
            nn.GELU()
        )
        self.conv7 = nn.Sequential(
            nn.Conv1d(in_channels, branch_ch, kernel_size=7, padding=3),
            nn.BatchNorm1d(branch_ch),
            nn.GELU()
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return torch.cat([self.conv3(x), self.conv5(x), self.conv7(x)], dim=1)


class MultiHeadSelfAttention(nn.Module):
    """Multi-head self-attention mechanism for sequence modeling."""

    def __init__(self, d_model: int, n_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        assert d_model % n_heads == 0
        self.d_model = d_model
        self.n_heads = n_heads
        self.d_k = d_model // n_heads

        self.W_q = nn.Linear(d_model, d_model)
        self.W_k = nn.Linear(d_model, d_model)
        self.W_v = nn.Linear(d_model, d_model)
        self.W_o = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(d_model)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        batch_size, seq_len, _ = x.shape
        residual = x

        Q = self.W_q(x).view(batch_size, seq_len, self.n_heads, self.d_k).transpose(1, 2)
        K = self.W_k(x).view(batch_size, seq_len, self.n_heads, self.d_k).transpose(1, 2)
        V = self.W_v(x).view(batch_size, seq_len, self.n_heads, self.d_k).transpose(1, 2)

        scores = torch.matmul(Q, K.transpose(-2, -1)) / math.sqrt(self.d_k)
        attn_weights = F.softmax(scores, dim=-1)
        attn_weights = self.dropout(attn_weights)

        context = torch.matmul(attn_weights, V)
        context = context.transpose(1, 2).contiguous().view(batch_size, seq_len, self.d_model)

        output = self.W_o(context)
        output = self.dropout(output)
        output = self.layer_norm(output + residual)

        return output, attn_weights


class EpigeneticFusionModule(nn.Module):
    """Gated fusion of sequence and epigenetic features."""

    def __init__(self, seq_dim: int, epi_dim: int, out_dim: int):
        super().__init__()
        self.seq_proj = nn.Linear(seq_dim, out_dim)
        self.epi_proj = nn.Linear(epi_dim, out_dim)
        self.gate = nn.Sequential(
            nn.Linear(seq_dim + epi_dim, out_dim),
            nn.Sigmoid()
        )

    def forward(self, seq_features: torch.Tensor, epi_features: torch.Tensor) -> torch.Tensor:
        seq_proj = self.seq_proj(seq_features)
        epi_proj = self.epi_proj(epi_features)
        gate = self.gate(torch.cat([seq_features, epi_features], dim=-1))
        return gate * seq_proj + (1 - gate) * epi_proj


class EpiCRISPRNet(nn.Module):
    """
    EpiCRISPR-Net: CNN + Attention model with epigenetic integration.
    
    Architecture:
    1. Input splitting: sequence features (27 ch) + epigenetic features (4 ch)
    2. Multi-scale CNN for local pattern extraction
    3. Epigenetic fusion via gated mechanism
    4. Multi-head self-attention for long-range dependencies
    5. Classification head with dropout
    """

    def __init__(
        self,
        input_channels: int = 31,
        seq_channels: int = 27,
        epi_channels: int = 4,
        cnn_channels: int = 96,
        attention_dim: int = 96,
        n_heads: int = 4,
        n_attention_layers: int = 2,
        dropout: float = 0.3,
        seq_len: int = 23
    ):
        super().__init__()
        self.seq_channels = seq_channels
        self.epi_channels = epi_channels

        # Multi-scale CNN for sequence features
        self.cnn_block1 = MultiScaleCNNBlock(seq_channels, cnn_channels)
        self.cnn_block2 = MultiScaleCNNBlock(cnn_channels, cnn_channels)
        self.cnn_dropout = nn.Dropout(dropout)

        # Epigenetic feature processing
        self.epi_encoder = nn.Sequential(
            nn.Linear(epi_channels, 32),
            nn.GELU(),
            nn.Linear(32, cnn_channels),
            nn.GELU()
        )

        # Gated fusion
        self.fusion = EpigeneticFusionModule(cnn_channels, cnn_channels, attention_dim)

        # Attention layers
        self.attention_layers = nn.ModuleList([
            MultiHeadSelfAttention(attention_dim, n_heads, dropout)
            for _ in range(n_attention_layers)
        ])

        # Classification head
        self.classifier = nn.Sequential(
            nn.Linear(attention_dim * seq_len, 128),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(128, 64),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1)
        )

        self._attention_weights = None

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        batch_size = x.shape[0]

        # Split input: x shape = (batch, seq_len, total_channels)
        seq_features = x[:, :, :self.seq_channels]  # (B, L, 27)
        epi_features = x[:, :, self.seq_channels:]   # (B, L, 4)

        # CNN expects (B, C, L)
        seq_cnn = seq_features.transpose(1, 2)
        seq_cnn = self.cnn_block1(seq_cnn)
        seq_cnn = self.cnn_block2(seq_cnn)
        seq_cnn = self.cnn_dropout(seq_cnn)
        seq_cnn = seq_cnn.transpose(1, 2)  # (B, L, cnn_channels)

        # Encode epigenetic features
        epi_encoded = self.epi_encoder(epi_features)  # (B, L, cnn_channels)

        # Gated fusion
        fused = self.fusion(seq_cnn, epi_encoded)  # (B, L, attention_dim)

        # Self-attention
        attn_out = fused
        for attn_layer in self.attention_layers:
            attn_out, attn_w = attn_layer(attn_out)
        self._attention_weights = attn_w

        # Classify
        flat = attn_out.reshape(batch_size, -1)
        logits = self.classifier(flat)

        return logits.squeeze(-1)

    def get_attention_weights(self) -> Optional[torch.Tensor]:
        return self._attention_weights


class BaselineCNN(nn.Module):
    """Simple CNN baseline without attention or epigenetic integration."""

    def __init__(self, input_channels: int = 31, seq_len: int = 23):
        super().__init__()
        self.conv1 = nn.Sequential(
            nn.Conv1d(input_channels, 64, kernel_size=3, padding=1),
            nn.BatchNorm1d(64), nn.ReLU()
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(64, 128, kernel_size=3, padding=1),
            nn.BatchNorm1d(128), nn.ReLU()
        )
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            nn.Linear(128, 64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.transpose(1, 2)
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.pool(x).squeeze(-1)
        return self.fc(x).squeeze(-1)


class SequenceOnlyModel(nn.Module):
    """Model using only sequence features (no epigenetics) for ablation."""

    def __init__(self, seq_channels: int = 27, seq_len: int = 23):
        super().__init__()
        self.cnn = nn.Sequential(
            nn.Conv1d(seq_channels, 64, kernel_size=3, padding=1),
            nn.BatchNorm1d(64), nn.GELU(),
            nn.Conv1d(64, 96, kernel_size=5, padding=2),
            nn.BatchNorm1d(96), nn.GELU()
        )
        self.fc = nn.Sequential(
            nn.Linear(96 * seq_len, 128), nn.GELU(), nn.Dropout(0.3),
            nn.Linear(128, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        seq = x[:, :, :27].transpose(1, 2)
        h = self.cnn(seq)
        h = h.reshape(h.size(0), -1)
        return self.fc(h).squeeze(-1)


if __name__ == '__main__':
    model = EpiCRISPRNet()
    print(f"EpiCRISPR-Net parameters: {sum(p.numel() for p in model.parameters()):,}")

    x = torch.randn(8, 23, 31)
    out = model(x)
    print(f"Input: {x.shape} -> Output: {out.shape}")
    print(f"Attention weights: {model.get_attention_weights().shape}")

    baseline = BaselineCNN()
    print(f"\nBaseline CNN parameters: {sum(p.numel() for p in baseline.parameters()):,}")

    seq_only = SequenceOnlyModel()
    print(f"Sequence-only parameters: {sum(p.numel() for p in seq_only.parameters()):,}")
