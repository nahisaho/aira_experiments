"""
Training Pipeline for the Antibody Design System
Synthetic data generation + model training + evaluation.
"""

import math
import time
import json
import os
import random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from typing import List, Tuple, Dict, Optional

from antibody_model import (
    AntibodyDesignModel, CDRStructureEncoder, CDRStructurePredictor,
    VOCAB_SIZE, PAD_IDX, AMINO_ACIDS, encode_sequence, decode_sequence,
    CDR_H3_MAX_LEN
)

# ─────────────────────────────────────────
# Reproducibility
# ─────────────────────────────────────────
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)


# ─────────────────────────────────────────
# 1. Synthetic Dataset Generation
# ─────────────────────────────────────────
def generate_synthetic_cdr_dataset(
    n_samples: int = 5000,
    min_len: int = 6,
    max_len: int = CDR_H3_MAX_LEN,
    seed: int = 42,
) -> List[Dict]:
    """
    Generate synthetic CDR-H3 sequences with pseudo-physical labels.
    In production, replace with real experimental data from SAbDab, PDB, etc.
    """
    rng = np.random.default_rng(seed)
    aa_list = list(AMINO_ACIDS)

    # Amino acid frequency profile similar to real CDR-H3 sequences
    aa_freq = np.array([
        2.0, 0.5, 3.0, 3.0, 0.2, 4.0, 2.0, 5.0, 1.0, 2.0,
        3.0, 1.5, 0.5, 3.0, 0.5, 4.0, 2.0, 0.8, 4.0, 3.0,
    ])
    aa_freq = aa_freq / aa_freq.sum()

    data = []
    for _ in range(n_samples):
        # CDR-H3 length distribution: typically 8–20 AA
        length = int(rng.choice(np.arange(min_len, max_len + 1), p=None))
        # Lengths follow roughly a Poisson distribution centered at 12
        weights = np.array([
            np.exp(-0.5 * ((l - 12) / 4) ** 2)
            for l in range(min_len, max_len + 1)
        ])
        weights /= weights.sum()
        length = int(rng.choice(np.arange(min_len, max_len + 1), p=weights))

        seq = "".join(rng.choice(aa_list, p=aa_freq) for _ in range(length))

        # Pseudo-torsion angles
        torsion = rng.uniform(-math.pi, math.pi, size=(length, 4)).astype(np.float32)

        # Pseudo-labels based on sequence properties
        hydro_vals = {
            "A": 1.8, "C": 2.5, "F": 2.8, "I": 4.5, "L": 3.8,
            "M": 1.9, "V": 4.2, "W": -0.9, "Y": -1.3,
        }
        hydro = np.mean([hydro_vals.get(aa, 0.0) for aa in seq])

        # Log Kd: more hydrophobic CDRs tend to bind better (simplified)
        log_kd = float(-8.0 - hydro * 0.3 + rng.normal(0, 0.5))

        # Tm: stability decreases with length and high flexibility
        tm = float(70.0 - length * 0.3 + hydro * 0.5 + rng.normal(0, 2.0))

        data.append({
            "sequence": seq,
            "torsion": torsion,
            "log_kd": log_kd,
            "tm": max(40.0, min(90.0, tm)),
            "delta_delta_G": float(rng.normal(0, 1.0)),
            "humanization": float(rng.beta(8, 2)),   # skewed towards high
            "immunogenicity": float(rng.beta(2, 8)),  # skewed towards low
        })

    return data


class AntibodyDataset(Dataset):
    """PyTorch dataset for antibody training."""

    def __init__(
        self,
        data: List[Dict],
        max_len: int = CDR_H3_MAX_LEN,
        antigen_len: int = 80,
    ):
        self.data = data
        self.max_len = max_len
        self.antigen_len = antigen_len
        # Fixed random antigen context (placeholder for PD-L1 sequence segment)
        rng = np.random.default_rng(SEED)
        self.ag_tokens = torch.tensor(
            rng.integers(0, VOCAB_SIZE, antigen_len), dtype=torch.long
        )

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx: int) -> Dict:
        item = self.data[idx]
        seq = item["sequence"]
        tok = encode_sequence(seq)
        L = len(tok)

        # Pad sequence
        tok_padded = F.pad(tok, (0, max(0, self.max_len - L)), value=PAD_IDX)
        mask = torch.zeros(self.max_len, dtype=torch.bool)
        mask[L:] = True  # True = padding position

        # Pad torsion
        torsion = torch.tensor(item["torsion"], dtype=torch.float32)
        torsion_padded = F.pad(torsion, (0, 0, 0, max(0, self.max_len - L)))

        return {
            "seq_tokens": tok_padded,
            "torsion": torsion_padded,
            "mask": mask,
            "antigen": self.ag_tokens,
            "log_kd": torch.tensor(item["log_kd"], dtype=torch.float32),
            "tm": torch.tensor(item["tm"], dtype=torch.float32),
            "delta_delta_G": torch.tensor(item["delta_delta_G"], dtype=torch.float32),
            "humanization": torch.tensor(item["humanization"], dtype=torch.float32),
            "length": L,
        }


# ─────────────────────────────────────────
# 2. Training Loop
# ─────────────────────────────────────────
class AntibodyTrainer:
    def __init__(
        self,
        model: AntibodyDesignModel,
        device: str = "cpu",
        lr: float = 3e-4,
        warmup_steps: int = 1000,
    ):
        self.model = model.to(device)
        self.device = device
        self.torsion_predictor = CDRStructurePredictor(model.encoder.d_model).to(device)
        self.optimizer = torch.optim.AdamW(
            list(model.parameters()) + list(self.torsion_predictor.parameters()),
            lr=lr, weight_decay=1e-4
        )

        # Cosine annealing with warmup
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
            self.optimizer, T_0=1000, T_mult=2
        )
        self.global_step = 0

    def compute_losses(self, batch: Dict) -> Dict[str, torch.Tensor]:
        seq_tokens = batch["seq_tokens"].to(self.device)
        torsion = batch["torsion"].to(self.device)
        mask = batch["mask"].to(self.device)
        antigen = batch["antigen"].to(self.device)
        log_kd_gt = batch["log_kd"].to(self.device)
        tm_gt = batch["tm"].to(self.device)
        ddg_gt = batch["delta_delta_G"].to(self.device)

        B = seq_tokens.size(0)

        # CDR encoding
        cdr_enc = self.model.encoder(seq_tokens, torsion, mask)

        # Antigen encoding (expand single antigen for whole batch)
        # antigen shape: (B, L_ag) — already batched from DataLoader
        ag_enc = self.model.encoder(antigen)

        # Torsion reconstruction (auxiliary task)
        torsion_pred = self.torsion_predictor(cdr_enc)
        # target: (sin, cos) pairs
        torsion_target = torch.cat([torch.sin(torsion), torch.cos(torsion)], dim=-1)
        torsion_mask = (~mask).unsqueeze(-1).float()
        torsion_loss = (F.mse_loss(torsion_pred, torsion_target, reduction="none") * torsion_mask).mean()

        # Affinity prediction loss
        log_kd_pred = self.model.affinity(cdr_enc, ag_enc)
        affinity_loss = F.mse_loss(log_kd_pred, log_kd_gt)

        # Stability prediction loss
        stab_pred = self.model.stability(cdr_enc)  # (B, 2): [ΔΔG, Tm]
        ddg_pred = stab_pred[:, 0]
        tm_pred = stab_pred[:, 1]
        stability_loss = F.mse_loss(ddg_pred, ddg_gt) + F.mse_loss(tm_pred, tm_gt / 90.0)

        # Diffusion denoising loss (simplified: predict tokens at t=T/2)
        T_half = self.model.diffusion.T // 2
        t_batch = torch.full((B,), T_half, device=self.device, dtype=torch.long)
        # Noisy tokens (random replacement for discrete diffusion)
        noisy_tokens = seq_tokens.clone()
        noise_mask = torch.rand_like(seq_tokens.float()) < 0.3
        random_tokens = torch.randint(0, VOCAB_SIZE, seq_tokens.shape, device=self.device)
        noisy_tokens = torch.where(noise_mask, random_tokens, noisy_tokens)

        fw_enc = ag_enc  # use antigen as framework proxy in this simplified version
        denoised_logits = self.model.diffusion(noisy_tokens, t_batch, ag_enc, fw_enc)
        denoising_loss = F.cross_entropy(
            denoised_logits.view(-1, VOCAB_SIZE),
            seq_tokens.view(-1),
            ignore_index=PAD_IDX,
        )

        # Total loss
        total = (
            0.4 * denoising_loss
            + 0.25 * affinity_loss
            + 0.20 * stability_loss
            + 0.15 * torsion_loss
        )

        return {
            "total": total,
            "denoising": denoising_loss.detach(),
            "affinity": affinity_loss.detach(),
            "stability": stability_loss.detach(),
            "torsion": torsion_loss.detach(),
        }

    def train_epoch(self, loader: DataLoader) -> Dict[str, float]:
        self.model.train()
        epoch_losses = {k: 0.0 for k in ["total", "denoising", "affinity", "stability", "torsion"]}
        n_batches = 0

        for batch in loader:
            self.optimizer.zero_grad()
            losses = self.compute_losses(batch)
            losses["total"].backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()
            self.scheduler.step()
            self.global_step += 1

            for k, v in losses.items():
                epoch_losses[k] += v.item() if isinstance(v, torch.Tensor) else v
            n_batches += 1

        return {k: v / n_batches for k, v in epoch_losses.items()}

    @torch.no_grad()
    def evaluate(self, loader: DataLoader) -> Dict[str, float]:
        self.model.eval()
        all_log_kd_pred, all_log_kd_true = [], []
        all_tm_pred, all_tm_true = [], []
        total_loss = 0.0
        n_batches = 0

        for batch in loader:
            seq_tokens = batch["seq_tokens"].to(self.device)
            torsion = batch["torsion"].to(self.device)
            mask = batch["mask"].to(self.device)
            antigen = batch["antigen"].to(self.device)
            log_kd_gt = batch["log_kd"]
            tm_gt = batch["tm"]

            B = seq_tokens.size(0)
            cdr_enc = self.model.encoder(seq_tokens, torsion, mask)
            ag_enc = self.model.encoder(antigen)

            log_kd_pred = self.model.affinity(cdr_enc, ag_enc).cpu()
            stab = self.model.stability(cdr_enc).cpu()
            tm_pred = stab[:, 1] * 90.0

            all_log_kd_pred.extend(log_kd_pred.numpy())
            all_log_kd_true.extend(log_kd_gt.numpy())
            all_tm_pred.extend(tm_pred.numpy())
            all_tm_true.extend(tm_gt.numpy())
            n_batches += 1

        from scipy.stats import pearsonr, spearmanr
        kd_preds = np.array(all_log_kd_pred)
        kd_true = np.array(all_log_kd_true)
        tm_preds = np.array(all_tm_pred)
        tm_true = np.array(all_tm_true)

        kd_r, _ = pearsonr(kd_preds, kd_true)
        kd_rho, _ = spearmanr(kd_preds, kd_true)
        kd_rmse = float(np.sqrt(np.mean((kd_preds - kd_true) ** 2)))
        tm_r, _ = pearsonr(tm_preds, tm_true)
        tm_rmse = float(np.sqrt(np.mean((tm_preds - tm_true) ** 2)))

        return {
            "log_kd_pearson_r": round(kd_r, 4),
            "log_kd_spearman_rho": round(kd_rho, 4),
            "log_kd_rmse": round(kd_rmse, 4),
            "tm_pearson_r": round(tm_r, 4),
            "tm_rmse": round(tm_rmse, 4),
        }


# ─────────────────────────────────────────
# 3. Full Training Run
# ─────────────────────────────────────────
def train_model(
    d_model: int = 256,
    n_epochs: int = 30,
    batch_size: int = 64,
    n_train: int = 4000,
    n_val: int = 500,
    lr: float = 3e-4,
    device: str = "cpu",
) -> Tuple[AntibodyDesignModel, Dict]:
    """Full training pipeline with synthetic data."""
    print(f"[Training] d_model={d_model}, epochs={n_epochs}, device={device}")

    # Generate datasets
    print("  Generating synthetic training data...")
    train_data = generate_synthetic_cdr_dataset(n_train, seed=42)
    val_data = generate_synthetic_cdr_dataset(n_val, seed=99)

    train_ds = AntibodyDataset(train_data)
    val_ds = AntibodyDataset(val_data)

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False, num_workers=0)

    # Initialize model
    model = AntibodyDesignModel(d_model=d_model, T=200)
    trainer = AntibodyTrainer(model, device=device, lr=lr)

    history = {
        "train_loss": [], "val_kd_pearson": [], "val_tm_pearson": [],
        "val_kd_rmse": [], "val_tm_rmse": [], "val_kd_spearman": []
    }

    print(f"  Starting training: {len(train_ds)} train / {len(val_ds)} val samples")
    for epoch in range(n_epochs):
        t0 = time.time()
        train_losses = trainer.train_epoch(train_loader)
        val_metrics = trainer.evaluate(val_loader)

        history["train_loss"].append(train_losses["total"])
        history["val_kd_pearson"].append(val_metrics["log_kd_pearson_r"])
        history["val_tm_pearson"].append(val_metrics["tm_pearson_r"])
        history["val_kd_rmse"].append(val_metrics["log_kd_rmse"])
        history["val_tm_rmse"].append(val_metrics["tm_rmse"])
        history["val_kd_spearman"].append(val_metrics["log_kd_spearman_rho"])

        elapsed = time.time() - t0
        if (epoch + 1) % 5 == 0 or epoch == 0:
            print(
                f"  Epoch {epoch+1:3d}/{n_epochs} | "
                f"Loss={train_losses['total']:.4f} | "
                f"Kd_r={val_metrics['log_kd_pearson_r']:.3f} | "
                f"Kd_ρ={val_metrics['log_kd_spearman_rho']:.3f} | "
                f"Kd_RMSE={val_metrics['log_kd_rmse']:.3f} | "
                f"Tm_r={val_metrics['tm_pearson_r']:.3f} | "
                f"{elapsed:.1f}s"
            )

    return model, history


if __name__ == "__main__":
    model, history = train_model(d_model=128, n_epochs=5, n_train=500, n_val=100, batch_size=32)
    print("Final val log_Kd Pearson r:", history["val_kd_pearson"][-1])
