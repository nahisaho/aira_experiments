"""Discrete diffusion model for synthetic antibody CDR-H3 sequence generation."""
from __future__ import annotations

import math
from typing import Callable, Optional, Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F

from data_utils import (
    AA_VOCAB_SIZE,
    AMINO_ACIDS,
    MAX_CDR_H3_LEN,
    MIN_CDR_H3_LEN,
    PAD_IDX,
    decode_sequence,
)


class TimeEmbedding(nn.Module):
    """Sinusoidal time embedding followed by a small MLP."""

    def __init__(self, hidden_dim: int) -> None:
        super().__init__()
        self.hidden_dim = hidden_dim
        self.mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        half_dim = self.hidden_dim // 2
        frequencies = torch.exp(
            -math.log(10000.0) * torch.arange(half_dim, device=t.device, dtype=torch.float32) / max(half_dim - 1, 1)
        )
        angles = t.float().unsqueeze(-1) * frequencies.unsqueeze(0)
        embedding = torch.cat([angles.sin(), angles.cos()], dim=-1)
        if embedding.size(-1) < self.hidden_dim:
            embedding = F.pad(embedding, (0, self.hidden_dim - embedding.size(-1)))
        return self.mlp(embedding)


class EquivariantAttention(nn.Module):
    """Attention with distance-aware bias derived from relative coordinates."""

    def __init__(self, hidden_dim: int, num_heads: int) -> None:
        super().__init__()
        if hidden_dim % num_heads != 0:
            raise ValueError("hidden_dim must be divisible by num_heads")
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        self.q_proj = nn.Linear(hidden_dim, hidden_dim)
        self.k_proj = nn.Linear(hidden_dim, hidden_dim)
        self.v_proj = nn.Linear(hidden_dim, hidden_dim)
        self.out_proj = nn.Linear(hidden_dim, hidden_dim)
        self.distance_bias = nn.Sequential(
            nn.Linear(1, num_heads),
            nn.Tanh(),
        )

    def forward(self, x: torch.Tensor, coords: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        batch_size, seq_len, hidden_dim = x.shape
        q = self.q_proj(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(self.head_dim)
        relative = coords[:, :, None, :] - coords[:, None, :, :]
        distances = torch.norm(relative, dim=-1, keepdim=True)
        bias = self.distance_bias(distances).permute(0, 3, 1, 2)
        scores = scores + bias
        if mask is not None:
            pair_mask = mask[:, None, None, :] * mask[:, None, :, None]
            scores = scores.masked_fill(pair_mask == 0, -1e4)
        attention = torch.softmax(scores, dim=-1)
        output = torch.matmul(attention, v).transpose(1, 2).contiguous().view(batch_size, seq_len, hidden_dim)
        return self.out_proj(output)


class DenoisingBlock(nn.Module):
    """Transformer block with geometry-aware attention."""

    def __init__(self, hidden_dim: int, num_heads: int, dropout: float = 0.1) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.attn = EquivariantAttention(hidden_dim, num_heads)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.ff = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
        )

    def forward(self, x: torch.Tensor, coords: torch.Tensor, mask: Optional[torch.Tensor]) -> torch.Tensor:
        x = x + self.attn(self.norm1(x), coords, mask)
        x = x + self.ff(self.norm2(x))
        return x


class DenoisingNetwork(nn.Module):
    """Transformer denoiser for discrete sequence diffusion."""

    def __init__(
        self,
        vocab_size: int,
        hidden_dim: int = 128,
        num_layers: int = 4,
        num_heads: int = 4,
        max_length: int = MAX_CDR_H3_LEN,
        condition_dim: Optional[int] = None,
    ) -> None:
        super().__init__()
        self.token_embedding = nn.Embedding(vocab_size, hidden_dim)
        self.position_embedding = nn.Embedding(max_length, hidden_dim)
        self.time_embedding = TimeEmbedding(hidden_dim)
        self.condition_proj = nn.Linear(condition_dim, hidden_dim) if condition_dim is not None else None
        self.layers = nn.ModuleList([DenoisingBlock(hidden_dim, num_heads) for _ in range(num_layers)])
        self.norm = nn.LayerNorm(hidden_dim)
        self.output = nn.Linear(hidden_dim, vocab_size)

    def forward(
        self,
        x_t: torch.Tensor,
        timesteps: torch.Tensor,
        coords: torch.Tensor,
        condition: Optional[torch.Tensor] = None,
        mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        batch_size, seq_len = x_t.shape
        positions = torch.arange(seq_len, device=x_t.device).unsqueeze(0).expand(batch_size, -1)
        hidden = self.token_embedding(x_t) + self.position_embedding(positions)
        hidden = hidden + self.time_embedding(timesteps).unsqueeze(1)
        if condition is not None and self.condition_proj is not None:
            hidden = hidden + self.condition_proj(condition).unsqueeze(1)
        for layer in self.layers:
            hidden = layer(hidden, coords, mask)
        return self.output(self.norm(hidden))


class DiscreteDiffusionModel(nn.Module):
    """D3PM-style discrete diffusion model for amino-acid sequences."""

    def __init__(
        self,
        vocab_size: int = AA_VOCAB_SIZE,
        hidden_dim: int = 128,
        num_layers: int = 4,
        num_heads: int = 4,
        max_length: int = MAX_CDR_H3_LEN,
        T: int = 100,
        condition_dim: Optional[int] = None,
    ) -> None:
        super().__init__()
        self.vocab_size = vocab_size
        self.max_length = max_length
        self.T = T
        self.standard_vocab_size = len(AMINO_ACIDS)
        self.denoiser = DenoisingNetwork(
            vocab_size=vocab_size,
            hidden_dim=hidden_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            max_length=max_length,
            condition_dim=condition_dim,
        )
        betas = torch.linspace(1e-3, 0.08, T)
        alphas = 1.0 - betas
        alpha_bars = torch.cumprod(alphas, dim=0)
        self.register_buffer("betas", betas)
        self.register_buffer("alphas", alphas)
        self.register_buffer("alpha_bars", alpha_bars)

    def _default_coords(self, lengths: torch.Tensor, device: torch.device) -> torch.Tensor:
        batch_size = lengths.size(0)
        seq_len = int(lengths.max().item()) if lengths.numel() else self.max_length
        seq_len = max(seq_len, self.max_length)
        positions = torch.arange(seq_len, device=device, dtype=torch.float32)
        base = torch.stack([
            positions,
            0.2 * torch.sin(positions / 2.5),
            0.2 * torch.cos(positions / 3.0),
        ], dim=-1)
        coords = base.unsqueeze(0).expand(batch_size, -1, -1).clone()
        return coords[:, : self.max_length]

    def q_sample(self, x0: torch.Tensor, timesteps: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """Apply forward diffusion by replacing tokens with random amino acids."""
        if mask is None:
            mask = (x0 != PAD_IDX).float()
        noise_prob = (1.0 - self.alpha_bars[timesteps]).unsqueeze(-1)
        random_tokens = torch.randint(0, self.standard_vocab_size, x0.shape, device=x0.device)
        replace_mask = (torch.rand_like(x0.float()) < noise_prob) & (mask > 0)
        x_t = torch.where(replace_mask, random_tokens, x0)
        return torch.where(mask > 0, x_t, torch.full_like(x_t, PAD_IDX))

    def predict_start_logits(
        self,
        x_t: torch.Tensor,
        timesteps: torch.Tensor,
        condition: Optional[torch.Tensor] = None,
        coords: Optional[torch.Tensor] = None,
        mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Predict clean-sequence logits from a noisy sample."""
        if mask is None:
            mask = (x_t != PAD_IDX).float()
        if coords is None:
            lengths = mask.sum(dim=-1).long().clamp_min(MIN_CDR_H3_LEN)
            coords = self._default_coords(lengths, x_t.device)
        return self.denoiser(x_t, timesteps, coords, condition=condition, mask=mask)

    def reverse_step(
        self,
        x_t: torch.Tensor,
        timesteps: torch.Tensor,
        condition: Optional[torch.Tensor] = None,
        coords: Optional[torch.Tensor] = None,
        mask: Optional[torch.Tensor] = None,
        temperature: float = 1.0,
        logits_override: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Single reverse diffusion step."""
        if mask is None:
            mask = (x_t != PAD_IDX).float()
        logits = logits_override if logits_override is not None else self.predict_start_logits(x_t, timesteps, condition, coords, mask)
        probs = torch.softmax(logits / max(temperature, 1e-6), dim=-1)
        sampled = torch.multinomial(probs.view(-1, self.vocab_size), 1).view_as(x_t)
        beta_t = self.betas[timesteps].unsqueeze(-1)
        keep_pred = (torch.rand_like(x_t.float()) > beta_t) & (mask > 0)
        random_tokens = torch.randint(0, self.standard_vocab_size, x_t.shape, device=x_t.device)
        x_prev = torch.where(keep_pred, sampled, random_tokens)
        x_prev = torch.where(mask > 0, x_prev, torch.full_like(x_prev, PAD_IDX))
        return x_prev, probs

    def training_step(
        self,
        batch: dict,
        optimizer: Optional[torch.optim.Optimizer] = None,
    ) -> float:
        """Run one denoising objective step and optionally update parameters."""
        self.train()
        x0 = batch["tokens"]
        mask = batch.get("mask", (x0 != PAD_IDX).float())
        coords = batch.get("coords")
        condition = batch.get("antigen_features")
        batch_size = x0.size(0)
        timesteps = torch.randint(0, self.T, (batch_size,), device=x0.device)
        x_t = self.q_sample(x0, timesteps, mask)
        logits = self.predict_start_logits(x_t, timesteps, condition=condition, coords=coords, mask=mask)
        loss = F.cross_entropy(logits.view(-1, self.vocab_size), x0.view(-1), reduction="none")
        loss = (loss * mask.reshape(-1)).sum() / mask.sum().clamp_min(1.0)
        if optimizer is not None:
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        return float(loss.detach().item())

    @torch.no_grad()
    def sample(
        self,
        batch_size: int,
        seq_length: Optional[int | Sequence[int]] = None,
        condition: Optional[torch.Tensor] = None,
        coords: Optional[torch.Tensor] = None,
        device: Optional[torch.device] = None,
        temperature: float = 1.0,
        return_tokens: bool = False,
        logit_modifier: Optional[Callable[[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor], torch.Tensor]] = None,
    ) -> list[str] | torch.Tensor:
        """Generate new CDR-H3 sequences by reverse diffusion."""
        if device is None:
            device = next(self.parameters()).device
        if seq_length is None:
            lengths = torch.randint(MIN_CDR_H3_LEN, self.max_length + 1, (batch_size,), device=device)
        elif isinstance(seq_length, Sequence):
            lengths = torch.tensor(list(seq_length), device=device, dtype=torch.long)
        else:
            lengths = torch.full((batch_size,), int(seq_length), device=device, dtype=torch.long)
        mask = (torch.arange(self.max_length, device=device).unsqueeze(0) < lengths.unsqueeze(1)).float()
        x_t = torch.randint(0, self.standard_vocab_size, (batch_size, self.max_length), device=device)
        x_t = torch.where(mask > 0, x_t, torch.full_like(x_t, PAD_IDX))
        coords = coords if coords is not None else self._default_coords(lengths, device)
        for step in reversed(range(self.T)):
            timestep = torch.full((batch_size,), step, device=device, dtype=torch.long)
            logits = self.predict_start_logits(x_t, timestep, condition=condition, coords=coords, mask=mask)
            if logit_modifier is not None:
                logits = logit_modifier(logits, x_t, timestep, mask)
            x_t, _ = self.reverse_step(
                x_t,
                timestep,
                condition=condition,
                coords=coords,
                mask=mask,
                temperature=temperature,
                logits_override=logits,
            )
        if return_tokens:
            return x_t
        return [decode_sequence(tokens) for tokens in x_t]


__all__ = [
    "DiscreteDiffusionModel",
    "DenoisingNetwork",
    "EquivariantAttention",
    "TimeEmbedding",
]
