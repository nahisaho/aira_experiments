"""
Training pipeline for the ESM AI Emulator.

Implements:
- Multi-objective loss (MSE + physics constraints)
- Learning rate scheduling with warm-up
- Gradient clipping for stability
- Ensemble training with diverse initialization
- Validation with early stopping
"""

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from typing import Dict, List
import json
import os
from datetime import datetime


class EmulatorTrainer:
    """Trainer for ESM Emulator with physics-constrained learning."""

    def __init__(self, model, config: dict = None, device: str = "cpu"):
        self.model = model.to(device)
        self.device = device
        self.config = config or self.default_config()
        self.history: List[Dict] = []

        self.optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=self.config["lr"],
            weight_decay=self.config["weight_decay"],
        )

        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
            self.optimizer,
            T_0=self.config["scheduler_t0"],
            T_mult=2,
        )

    @staticmethod
    def default_config() -> dict:
        return {
            "lr": 1e-3,
            "weight_decay": 1e-5,
            "epochs": 50,
            "batch_size": 16,
            "grad_clip": 1.0,
            "scheduler_t0": 10,
            "val_split": 0.2,
            "patience": 10,
            "min_delta": 1e-4,
        }

    def train_epoch(self, dataloader: DataLoader) -> Dict[str, float]:
        self.model.train()
        epoch_losses = {"total": 0, "mse": 0}
        n_batches = 0

        for seq, target, scenario_id, forcing in dataloader:
            seq = seq.to(self.device)
            target = target.to(self.device)
            scenario_id = scenario_id.to(self.device)
            forcing = forcing.to(self.device)

            self.optimizer.zero_grad()
            outputs = self.model(seq, scenario_id, forcing)
            losses = self.model.compute_loss(outputs, target)

            losses["total"].backward()
            torch.nn.utils.clip_grad_norm_(
                self.model.parameters(), self.config["grad_clip"]
            )
            self.optimizer.step()

            for k, v in losses.items():
                epoch_losses[k] = epoch_losses.get(k, 0) + v.item()
            n_batches += 1

        return {k: v / max(n_batches, 1) for k, v in epoch_losses.items()}

    @torch.no_grad()
    def validate(self, dataloader: DataLoader) -> Dict[str, float]:
        self.model.eval()
        val_losses = {"total": 0, "mse": 0}
        n_batches = 0

        for seq, target, scenario_id, forcing in dataloader:
            seq = seq.to(self.device)
            target = target.to(self.device)
            scenario_id = scenario_id.to(self.device)
            forcing = forcing.to(self.device)

            outputs = self.model(seq, scenario_id, forcing)
            losses = self.model.compute_loss(outputs, target)

            for k, v in losses.items():
                val_losses[k] = val_losses.get(k, 0) + v.item()
            n_batches += 1

        return {k: v / max(n_batches, 1) for k, v in val_losses.items()}

    def train(self, dataset, log_dir: str = "logs") -> Dict:
        """Full training loop with validation and early stopping."""
        cfg = self.config
        n_val = int(len(dataset) * cfg["val_split"])
        n_train = len(dataset) - n_val

        train_ds, val_ds = random_split(
            dataset, [n_train, n_val],
            generator=torch.Generator().manual_seed(42),
        )

        train_loader = DataLoader(
            train_ds, batch_size=cfg["batch_size"], shuffle=True, drop_last=True,
        )
        val_loader = DataLoader(
            val_ds, batch_size=cfg["batch_size"], shuffle=False,
        )

        best_val_loss = float("inf")
        patience_counter = 0

        for epoch in range(cfg["epochs"]):
            train_losses = self.train_epoch(train_loader)
            val_losses = self.validate(val_loader)
            self.scheduler.step()

            record = {
                "epoch": epoch,
                "train": train_losses,
                "val": val_losses,
                "lr": self.optimizer.param_groups[0]["lr"],
            }
            self.history.append(record)

            # Early stopping
            if val_losses["total"] < best_val_loss - cfg["min_delta"]:
                best_val_loss = val_losses["total"]
                patience_counter = 0
                best_state = {k: v.cpu().clone()
                              for k, v in self.model.state_dict().items()}
            else:
                patience_counter += 1
                if patience_counter >= cfg["patience"]:
                    break

        # Restore best model
        if best_state:
            self.model.load_state_dict(best_state)

        return {
            "best_val_loss": best_val_loss,
            "epochs_trained": len(self.history),
            "history": self.history,
        }


def train_ensemble(model_factory, dataset, n_members: int = 5,
                   config: dict = None, device: str = "cpu") -> List:
    """
    Train an ensemble of emulators with diverse initializations.

    Each member uses a different random seed for parameter initialization
    and data shuffling to maximize prediction diversity.
    """
    ensemble_results = []

    for i in range(n_members):
        torch.manual_seed(i * 1000 + 42)
        model = model_factory()
        trainer = EmulatorTrainer(model, config, device)
        result = trainer.train(dataset)
        result["member_id"] = i
        ensemble_results.append({
            "member_id": i,
            "best_val_loss": result["best_val_loss"],
            "epochs_trained": result["epochs_trained"],
            "model_state": model.state_dict(),
        })

    return ensemble_results
