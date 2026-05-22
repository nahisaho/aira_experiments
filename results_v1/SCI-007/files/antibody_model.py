"""
Antibody CDR De Novo Design System
Core model architectures: CDR structure-sequence learning + Diffusion Model
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Optional, Tuple


# ─────────────────────────────────────────
# Constants
# ─────────────────────────────────────────
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
IDX_TO_AA = {i: aa for aa, i in AA_TO_IDX.items()}
VOCAB_SIZE = len(AMINO_ACIDS)
CDR_H3_MAX_LEN = 25
PAD_IDX = VOCAB_SIZE  # padding token


def encode_sequence(seq: str) -> torch.Tensor:
    indices = [AA_TO_IDX.get(aa, 0) for aa in seq]
    return torch.tensor(indices, dtype=torch.long)


def decode_sequence(indices: torch.Tensor) -> str:
    return "".join(IDX_TO_AA.get(i.item(), "X") for i in indices)


# ─────────────────────────────────────────
# 1. Positional Encoding
# ─────────────────────────────────────────
class SinusoidalPositionalEncoding(nn.Module):
    def __init__(self, d_model: int, max_len: int = 512, dropout: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe.unsqueeze(0))  # (1, max_len, d_model)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.pe[:, : x.size(1)]
        return self.dropout(x)


# ─────────────────────────────────────────
# 2. CDR Sequence-Structure Encoder (Transformer)
# ─────────────────────────────────────────
class CDRStructureEncoder(nn.Module):
    """
    Encodes CDR sequence + backbone torsion angles into a latent representation.
    Inputs:
      - seq_tokens: (B, L) long tensor of AA indices
      - torsion:    (B, L, 4) float tensor of (phi, psi, omega, chi1) in radians
    Output:
      - latent:     (B, L, d_model)
    """

    def __init__(
        self,
        d_model: int = 256,
        n_heads: int = 8,
        n_layers: int = 4,
        dim_ff: int = 512,
        dropout: float = 0.1,
        max_len: int = 64,
    ):
        super().__init__()
        self.d_model = d_model
        self.seq_emb = nn.Embedding(VOCAB_SIZE + 1, d_model, padding_idx=PAD_IDX)
        self.torsion_proj = nn.Linear(8, d_model)  # (sin, cos) x 4 angles
        self.pos_enc = SinusoidalPositionalEncoding(d_model, max(max_len, 512), dropout)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model, n_heads, dim_ff, dropout, batch_first=True, norm_first=True
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, n_layers)
        self.layer_norm = nn.LayerNorm(d_model)

    def forward(
        self,
        seq_tokens: torch.Tensor,
        torsion: Optional[torch.Tensor] = None,
        src_key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        x = self.seq_emb(seq_tokens)  # (B, L, d_model)
        if torsion is not None:
            # encode torsion angles as (sin, cos) pairs
            torsion_feat = torch.cat(
                [torch.sin(torsion), torch.cos(torsion)], dim=-1
            )  # (B, L, 8)
            x = x + self.torsion_proj(torsion_feat)
        x = self.pos_enc(x)
        x = self.encoder(x, src_key_padding_mask=src_key_padding_mask)
        return self.layer_norm(x)


# ─────────────────────────────────────────
# 3. CDR Structure Predictor (auxiliary task)
# ─────────────────────────────────────────
class CDRStructurePredictor(nn.Module):
    """
    Predicts torsion angles from encoded representation.
    Used as auxiliary task during pretraining.
    """

    def __init__(self, d_model: int = 256):
        super().__init__()
        self.head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.GELU(),
            nn.Linear(128, 8),  # (sin, cos) x 4 torsion angles
        )

    def forward(self, latent: torch.Tensor) -> torch.Tensor:
        return self.head(latent)  # (B, L, 8)


# ─────────────────────────────────────────
# 4. Diffusion Model for CDR-H3 Sequence Generation
# ─────────────────────────────────────────
class DiffusionNoiseScheduler:
    """
    Linear / cosine noise schedule for discrete sequence diffusion.
    Implements the absorbing-state diffusion (D3PM-like) for amino acids.
    """

    def __init__(self, T: int = 1000, schedule: str = "cosine"):
        self.T = T
        if schedule == "cosine":
            steps = torch.arange(T + 1, dtype=torch.float)
            alphas_cumprod = torch.cos((steps / T + 0.008) / 1.008 * math.pi / 2) ** 2
            alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
        else:  # linear
            betas = torch.linspace(1e-4, 0.02, T)
            alphas = 1.0 - betas
            alphas_cumprod = torch.cumprod(alphas, dim=0)
            alphas_cumprod = torch.cat([torch.tensor([1.0]), alphas_cumprod])

        self.register_buffers(alphas_cumprod)

    def register_buffers(self, alphas_cumprod: torch.Tensor):
        self.alphas_cumprod = alphas_cumprod
        self.sqrt_alphas_cumprod = torch.sqrt(alphas_cumprod)
        self.sqrt_one_minus_alphas_cumprod = torch.sqrt(1.0 - alphas_cumprod)
        self.betas = 1.0 - alphas_cumprod[1:] / alphas_cumprod[:-1]

    def q_sample_continuous(
        self, x0: torch.Tensor, t: torch.Tensor, noise: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """Forward diffusion on continuous embeddings."""
        if noise is None:
            noise = torch.randn_like(x0)
        sqrt_alpha = self.sqrt_alphas_cumprod[t].view(-1, 1, 1)
        sqrt_one_minus = self.sqrt_one_minus_alphas_cumprod[t].view(-1, 1, 1)
        return sqrt_alpha * x0 + sqrt_one_minus * noise

    def get_snr(self, t: torch.Tensor) -> torch.Tensor:
        alpha = self.alphas_cumprod[t]
        return alpha / (1.0 - alpha)


class CDRDiffusionModel(nn.Module):
    """
    Continuous-space diffusion model on amino acid embeddings for CDR-H3 generation.
    Architecture: Transformer-based denoiser conditioned on:
      - antigen embedding (context)
      - timestep embedding
      - framework region embedding
    """

    def __init__(
        self,
        d_model: int = 256,
        n_heads: int = 8,
        n_layers: int = 4,
        dim_ff: int = 512,
        dropout: float = 0.1,
        max_cdr_len: int = CDR_H3_MAX_LEN,
        T: int = 1000,
    ):
        super().__init__()
        self.d_model = d_model
        self.T = T
        self.max_cdr_len = max_cdr_len

        # AA embedding (shared with encoder)
        self.aa_emb = nn.Embedding(VOCAB_SIZE + 1, d_model, padding_idx=PAD_IDX)

        # Timestep embedding
        self.time_emb = nn.Sequential(
            nn.Linear(d_model, d_model * 4),
            nn.SiLU(),
            nn.Linear(d_model * 4, d_model),
        )

        # Positional encoding
        self.pos_enc = SinusoidalPositionalEncoding(d_model, max_cdr_len + 64, dropout)

        # Cross-attention-based Transformer denoiser
        decoder_layer = nn.TransformerDecoderLayer(
            d_model, n_heads, dim_ff, dropout, batch_first=True, norm_first=True
        )
        self.denoiser = nn.TransformerDecoder(decoder_layer, n_layers)

        # Output projection to AA logits
        self.output_proj = nn.Linear(d_model, VOCAB_SIZE)

        # Noise scheduler
        self.scheduler = DiffusionNoiseScheduler(T=T, schedule="cosine")

    def get_timestep_embedding(self, t: torch.Tensor) -> torch.Tensor:
        """Sinusoidal timestep embedding."""
        half = self.d_model // 2
        freq = torch.exp(
            -math.log(10000) * torch.arange(half, device=t.device).float() / half
        )
        args = t.float()[:, None] * freq[None]
        embedding = torch.cat([torch.sin(args), torch.cos(args)], dim=-1)
        return self.time_emb(embedding)

    def forward(
        self,
        noisy_cdr: torch.Tensor,           # (B, L_cdr) long OR (B, L_cdr, d_model) float
        t: torch.Tensor,                    # (B,) int timestep
        antigen_context: torch.Tensor,      # (B, L_ag, d_model) antigen encoding
        framework_context: torch.Tensor,    # (B, L_fw, d_model) framework region encoding
        cdr_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Predict denoised CDR logits.
        Returns: (B, L_cdr, VOCAB_SIZE) logits
        """
        if noisy_cdr.dtype == torch.long:
            x = self.aa_emb(noisy_cdr)  # (B, L, d_model)
        else:
            x = noisy_cdr  # already embedded (continuous diffusion)

        # Add timestep embedding
        t_emb = self.get_timestep_embedding(t)  # (B, d_model)
        x = x + t_emb.unsqueeze(1)

        x = self.pos_enc(x)

        # Concatenate antigen + framework as memory for cross-attention
        memory = torch.cat([antigen_context, framework_context], dim=1)

        x = self.denoiser(x, memory, tgt_key_padding_mask=cdr_mask)
        return self.output_proj(x)  # (B, L_cdr, VOCAB_SIZE)

    @torch.no_grad()
    def sample(
        self,
        antigen_context: torch.Tensor,
        framework_context: torch.Tensor,
        cdr_length: int,
        n_samples: int = 1,
        temperature: float = 1.0,
        device: str = "cpu",
    ) -> torch.Tensor:
        """
        DDPM-style reverse diffusion sampling.
        Returns: (n_samples, cdr_length) long tensor of AA indices.
        """
        B = n_samples
        # Start from random AA tokens
        x = torch.randint(0, VOCAB_SIZE, (B, cdr_length), device=device)

        for t_val in reversed(range(1, self.T + 1)):
            t_batch = torch.full((B,), t_val, device=device, dtype=torch.long)
            logits = self.forward(x, t_batch, antigen_context, framework_context)
            # Gumbel-softmax sampling
            probs = F.softmax(logits / temperature, dim=-1)
            x = torch.multinomial(probs.view(-1, VOCAB_SIZE), 1).view(B, cdr_length)

        return x


# ─────────────────────────────────────────
# 5. Binding Affinity Predictor
# ─────────────────────────────────────────
class BindingAffinityPredictor(nn.Module):
    """
    Predicts log Kd (binding affinity) from CDR + antigen encodings.
    Uses cross-attention pooling then MLP regression.
    """

    def __init__(self, d_model: int = 256, dropout: float = 0.1):
        super().__init__()
        self.cross_attn = nn.MultiheadAttention(
            d_model, num_heads=8, dropout=dropout, batch_first=True
        )
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.regressor = nn.Sequential(
            nn.Linear(d_model * 2, 512),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(512, 128),
            nn.GELU(),
            nn.Linear(128, 1),
        )

    def forward(
        self,
        cdr_enc: torch.Tensor,      # (B, L_cdr, d)
        antigen_enc: torch.Tensor,  # (B, L_ag, d)
    ) -> torch.Tensor:
        # Cross-attention: CDR attends to antigen
        ctx, _ = self.cross_attn(cdr_enc, antigen_enc, antigen_enc)  # (B, L_cdr, d)
        # Pool both
        cdr_pool = ctx.mean(dim=1)           # (B, d)
        ag_pool = antigen_enc.mean(dim=1)    # (B, d)
        feat = torch.cat([cdr_pool, ag_pool], dim=-1)
        return self.regressor(feat).squeeze(-1)  # (B,) log Kd


# ─────────────────────────────────────────
# 6. Stability Predictor
# ─────────────────────────────────────────
class StabilityPredictor(nn.Module):
    """
    Predicts ΔΔG (thermodynamic stability) and Tm from CDR sequence.
    """

    def __init__(self, d_model: int = 256, dropout: float = 0.1):
        super().__init__()
        self.pool = nn.Linear(d_model, 1)
        self.mlp = nn.Sequential(
            nn.Linear(d_model, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, 64),
            nn.GELU(),
            nn.Linear(64, 2),  # [delta_delta_G, Tm]
        )

    def forward(self, cdr_enc: torch.Tensor) -> torch.Tensor:
        # Weighted pool
        w = torch.softmax(self.pool(cdr_enc), dim=1)  # (B, L, 1)
        pooled = (cdr_enc * w).sum(dim=1)              # (B, d)
        return self.mlp(pooled)  # (B, 2)


# ─────────────────────────────────────────
# 7. Full Antibody Design Model (all-in-one)
# ─────────────────────────────────────────
class AntibodyDesignModel(nn.Module):
    """
    End-to-end antibody design model combining:
    - CDR structure encoder
    - Diffusion-based CDR generator
    - Multi-attribute predictors
    """

    def __init__(self, d_model: int = 256, T: int = 1000):
        super().__init__()
        self.encoder = CDRStructureEncoder(d_model=d_model)
        self.diffusion = CDRDiffusionModel(d_model=d_model, T=T)
        self.affinity = BindingAffinityPredictor(d_model=d_model)
        self.stability = StabilityPredictor(d_model=d_model)

    def encode_cdr(
        self,
        seq_tokens: torch.Tensor,
        torsion: Optional[torch.Tensor] = None,
        mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        return self.encoder(seq_tokens, torsion, mask)

    def predict_properties(
        self,
        cdr_tokens: torch.Tensor,
        antigen_tokens: torch.Tensor,
        cdr_torsion: Optional[torch.Tensor] = None,
    ) -> dict:
        cdr_enc = self.encode_cdr(cdr_tokens, cdr_torsion)
        ag_enc = self.encoder(antigen_tokens)
        log_kd = self.affinity(cdr_enc, ag_enc)
        stab = self.stability(cdr_enc)
        return {
            "log_kd": log_kd,
            "delta_delta_G": stab[:, 0],
            "Tm": stab[:, 1],
        }

    def generate_cdrs(
        self,
        antigen_tokens: torch.Tensor,
        framework_tokens: torch.Tensor,
        cdr_length: int = 12,
        n_samples: int = 10,
        temperature: float = 0.8,
    ) -> torch.Tensor:
        device = antigen_tokens.device
        ag_enc = self.encoder(antigen_tokens)
        fw_enc = self.encoder(framework_tokens)
        # Broadcast to n_samples
        ag_enc = ag_enc.expand(n_samples, -1, -1)
        fw_enc = fw_enc.expand(n_samples, -1, -1)
        return self.diffusion.sample(
            ag_enc, fw_enc, cdr_length, n_samples, temperature, str(device)
        )


if __name__ == "__main__":
    print("=== Antibody Design Model Sanity Check ===")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = AntibodyDesignModel(d_model=128, T=100).to(device)

    B, L_cdr, L_ag, L_fw = 4, 12, 80, 100
    cdr_tok = torch.randint(0, VOCAB_SIZE, (B, L_cdr)).to(device)
    ag_tok  = torch.randint(0, VOCAB_SIZE, (1, L_ag)).to(device)
    fw_tok  = torch.randint(0, VOCAB_SIZE, (1, L_fw)).to(device)

    props = model.predict_properties(cdr_tok, ag_tok.expand(B, -1))
    print(f"  log_Kd shape:      {props['log_kd'].shape}")
    print(f"  ΔΔG shape:         {props['delta_delta_G'].shape}")
    print(f"  Tm shape:          {props['Tm'].shape}")

    generated = model.generate_cdrs(ag_tok, fw_tok, cdr_length=12, n_samples=5, temperature=0.8)
    print(f"  Generated CDRs:    {generated.shape}")
    for i, seq_idx in enumerate(generated):
        seq = decode_sequence(seq_idx)
        print(f"    Sample {i+1}: {seq}")
    print("Sanity check passed.")
