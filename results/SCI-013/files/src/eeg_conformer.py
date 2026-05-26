"""EEG Conformer architecture for EEG decoding."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import torch
from torch import nn


class PatchEmbedding(nn.Module):
    """Temporal-spatial EEG patch embedding."""

    def __init__(self, n_channels: int, embed_dim: int = 64, kernel_length: int = 25, patch_stride: int = 8) -> None:
        super().__init__()
        self.temporal = nn.Conv2d(1, 16, kernel_size=(1, kernel_length), padding=(0, kernel_length // 2), bias=False)
        self.temporal_bn = nn.BatchNorm2d(16)
        self.spatial = nn.Conv2d(16, embed_dim, kernel_size=(n_channels, 1), bias=False)
        self.spatial_bn = nn.BatchNorm2d(embed_dim)
        self.activation = nn.GELU()
        self.pool = nn.AvgPool1d(kernel_size=patch_stride, stride=patch_stride)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.unsqueeze(1)
        x = self.activation(self.temporal_bn(self.temporal(x)))
        x = self.activation(self.spatial_bn(self.spatial(x))).squeeze(2)
        x = self.pool(x)
        return x.transpose(1, 2)


class ConformerBlock(nn.Module):
    """Conformer block with attention, convolution, and feed-forward modules."""

    def __init__(self, embed_dim: int = 64, num_heads: int = 4, ff_mult: int = 4, conv_kernel: int = 15, dropout: float = 0.1) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(embed_dim)
        self.attn = nn.MultiheadAttention(embed_dim, num_heads, dropout=dropout, batch_first=True)
        self.norm2 = nn.LayerNorm(embed_dim)
        self.ff1 = nn.Sequential(
            nn.Linear(embed_dim, ff_mult * embed_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(ff_mult * embed_dim, embed_dim),
            nn.Dropout(dropout),
        )
        self.norm3 = nn.LayerNorm(embed_dim)
        self.conv_pw1 = nn.Conv1d(embed_dim, 2 * embed_dim, kernel_size=1)
        self.conv_dw = nn.Conv1d(embed_dim, embed_dim, kernel_size=conv_kernel, padding=conv_kernel // 2, groups=embed_dim)
        self.conv_bn = nn.BatchNorm1d(embed_dim)
        self.conv_pw2 = nn.Conv1d(embed_dim, embed_dim, kernel_size=1)
        self.norm4 = nn.LayerNorm(embed_dim)
        self.ff2 = nn.Sequential(
            nn.Linear(embed_dim, ff_mult * embed_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(ff_mult * embed_dim, embed_dim),
            nn.Dropout(dropout),
        )
        self.last_attention: Optional[torch.Tensor] = None

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        attn_input = self.norm1(x)
        attn_out, attn_weights = self.attn(attn_input, attn_input, attn_input, need_weights=True, average_attn_weights=False)
        self.last_attention = attn_weights.detach()
        x = x + attn_out
        x = x + 0.5 * self.ff1(self.norm2(x))
        conv_in = self.norm3(x).transpose(1, 2)
        conv_in = self.conv_pw1(conv_in)
        a, b = conv_in.chunk(2, dim=1)
        conv_in = a * torch.sigmoid(b)
        conv_in = self.conv_bn(self.conv_dw(conv_in))
        conv_in = torch.nn.functional.silu(conv_in)
        conv_in = self.conv_pw2(conv_in).transpose(1, 2)
        x = x + conv_in
        x = x + 0.5 * self.ff2(self.norm4(x))
        return x


class EEGConformer(nn.Module):
    """Convolutional Transformer for EEG classification."""

    def __init__(self, n_channels: int = 22, n_classes: int = 4, embed_dim: int = 64, depth: int = 3, num_heads: int = 4) -> None:
        super().__init__()
        self.patch_embedding = PatchEmbedding(n_channels=n_channels, embed_dim=embed_dim)
        self.blocks = nn.ModuleList([ConformerBlock(embed_dim=embed_dim, num_heads=num_heads) for _ in range(depth)])
        self.norm = nn.LayerNorm(embed_dim)
        self.classifier = nn.Sequential(nn.Linear(embed_dim, embed_dim), nn.GELU(), nn.Linear(embed_dim, n_classes))
        self._last_tokens: Optional[torch.Tensor] = None

    def forward(self, x: torch.Tensor, return_attention: bool = False) -> torch.Tensor | Tuple[torch.Tensor, torch.Tensor]:
        tokens = self.patch_embedding(x)
        for block in self.blocks:
            tokens = block(tokens)
        tokens = self.norm(tokens)
        if tokens.requires_grad:
            tokens.retain_grad()
        self._last_tokens = tokens
        pooled = tokens.mean(dim=1)
        logits = self.classifier(pooled)
        if return_attention:
            attention = self.get_attention_map()
            return logits, attention if attention is not None else torch.empty(0)
        return logits

    def get_attention_map(self) -> Optional[torch.Tensor]:
        if not self.blocks or self.blocks[-1].last_attention is None:
            return None
        return self.blocks[-1].last_attention.mean(dim=1)

    def class_activation_map(self, x: torch.Tensor, target_class: Optional[int] = None) -> np.ndarray:
        self.zero_grad(set_to_none=True)
        logits = self.forward(x)
        if isinstance(logits, tuple):
            logits = logits[0]
        if target_class is None:
            target_class = int(torch.argmax(logits, dim=1)[0].item())
        score = logits[:, target_class].sum()
        score.backward()
        if self._last_tokens is None or self._last_tokens.grad is None:
            raise RuntimeError("Run a forward/backward pass before requesting class activation map.")
        cam = torch.relu((self._last_tokens.grad * self._last_tokens).sum(dim=-1))
        return cam.detach().cpu().numpy()
