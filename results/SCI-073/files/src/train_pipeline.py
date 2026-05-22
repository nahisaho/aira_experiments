"""
Unified Training Pipeline for Tactile Manipulation System.
Supports training all 6 modules with configurable hyperparameters,
logging, and evaluation.
"""

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
import numpy as np
import json
import os
import time
from pathlib import Path
from typing import Dict, Optional


class TactileDataset(torch.utils.data.Dataset):
    """Generic tactile dataset loader."""

    def __init__(self, data_path: str, task: str = "contact"):
        self.task = task
        if os.path.exists(data_path):
            data = np.load(data_path, allow_pickle=True)
            self._load_data(data)
        else:
            self._generate_synthetic(1000)

    def _load_data(self, data):
        if self.task == "contact":
            self.images = torch.from_numpy(data["tactile_images"]).permute(0, 3, 1, 2).float()
            self.depth = torch.from_numpy(data["depth_maps"]).unsqueeze(1).float()
            self.normals = torch.from_numpy(data["normal_maps"]).permute(0, 3, 1, 2).float()
            self.wrench = torch.from_numpy(data["wrenches"]).float()
        elif self.task == "texture":
            self.images = torch.from_numpy(data["images"]).permute(0, 3, 1, 2).float()
            self.labels = torch.from_numpy(data["labels"]).long()

    def _generate_synthetic(self, n: int):
        """Generate minimal synthetic data for testing."""
        if self.task == "contact":
            self.images = torch.randn(n, 3, 240, 320)
            self.depth = torch.randn(n, 1, 64, 64)
            self.normals = torch.randn(n, 3, 64, 64)
            self.normals = torch.nn.functional.normalize(self.normals, dim=1)
            self.wrench = torch.randn(n, 6)
        elif self.task == "texture":
            self.images = torch.randn(n, 3, 224, 224)
            self.labels = torch.randint(0, 20, (n,))
        elif self.task == "multimodal":
            self.tactile = torch.randn(n, 3, 224, 224)
            self.visual = torch.randn(n, 3, 224, 224)
            self.object_labels = torch.randint(0, 50, (n,))
            self.material_labels = torch.randint(0, 20, (n,))
            self.grasp_quality = torch.rand(n)
        elif self.task == "stability":
            self.tactile_seq = torch.randn(n, 10, 3, 64, 64)
            self.wrench_seq = torch.randn(n, 10, 6)
            self.stability_score = torch.rand(n)
            self.stability_class = torch.randint(0, 3, (n,))
            self.ttf = torch.rand(n) * 5
            self.action = torch.randint(0, 4, (n,))
        elif self.task == "slip":
            self.tactile_seq = torch.randn(n, 10, 3, 64, 64)
            self.force = torch.randn(n, 3)
            self.target_force = torch.randn(n, 3)
            self.slip_class = torch.randint(0, 4, (n,))
        elif self.task == "exploration":
            self.visual = torch.randn(n, 3, 224, 224)
            self.true_quality = torch.rand(n)
            self.safety_labels = torch.ones(n, 50)

    def __len__(self):
        if self.task == "contact":
            return len(self.images)
        elif self.task == "texture":
            return len(self.images)
        elif self.task == "multimodal":
            return len(self.tactile)
        elif self.task in ("stability", "slip"):
            return len(self.tactile_seq)
        elif self.task == "exploration":
            return len(self.visual)
        return 0

    def __getitem__(self, idx):
        if self.task == "contact":
            return {
                "image": self.images[idx],
                "depth": self.depth[idx],
                "normals": self.normals[idx],
                "wrench": self.wrench[idx],
            }
        elif self.task == "texture":
            return {"image": self.images[idx], "label": self.labels[idx]}
        elif self.task == "multimodal":
            return {
                "tactile": self.tactile[idx],
                "visual": self.visual[idx],
                "object_label": self.object_labels[idx],
                "material_label": self.material_labels[idx],
                "grasp_quality": self.grasp_quality[idx],
            }
        elif self.task == "stability":
            return {
                "tactile_seq": self.tactile_seq[idx],
                "wrench_seq": self.wrench_seq[idx],
                "stability_score": self.stability_score[idx],
                "stability_class": self.stability_class[idx],
                "time_to_failure": self.ttf[idx],
                "corrective_action": self.action[idx],
            }
        elif self.task == "slip":
            return {
                "tactile_seq": self.tactile_seq[idx],
                "force": self.force[idx],
                "target_force": self.target_force[idx],
                "slip_class": self.slip_class[idx],
            }
        elif self.task == "exploration":
            return {
                "visual": self.visual[idx],
                "true_quality": self.true_quality[idx],
                "safety_labels": self.safety_labels[idx],
            }


class MetricsTracker:
    """Track and log training metrics."""

    def __init__(self, log_dir: str = "logs/"):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.history = {}

    def update(self, phase: str, metrics: Dict[str, float], epoch: int):
        key = f"{phase}"
        if key not in self.history:
            self.history[key] = []
        self.history[key].append({"epoch": epoch, **metrics})

    def save(self, filename: str = "training_metrics.json"):
        serializable = {}
        for k, v in self.history.items():
            serializable[k] = []
            for entry in v:
                serializable[k].append(
                    {kk: float(vv) if isinstance(vv, (torch.Tensor, np.floating)) else vv
                     for kk, vv in entry.items()}
                )
        with open(self.log_dir / filename, 'w') as f:
            json.dump(serializable, f, indent=2)


class Trainer:
    """Unified trainer for all modules."""

    def __init__(self, model: nn.Module, loss_fn: nn.Module,
                 optimizer: optim.Optimizer, scheduler=None,
                 device: str = "cpu", log_dir: str = "logs/"):
        self.model = model.to(device)
        self.loss_fn = loss_fn.to(device) if hasattr(loss_fn, 'to') else loss_fn
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.device = device
        self.metrics = MetricsTracker(log_dir)
        self.best_val_loss = float('inf')

    def train_epoch(self, dataloader: DataLoader, task: str) -> Dict[str, float]:
        self.model.train()
        total_losses = {}
        n_batches = 0

        for batch in dataloader:
            batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                     for k, v in batch.items()}

            self.optimizer.zero_grad()

            if task == "contact":
                pred = self.model(batch["image"])
                target = {"depth": batch["depth"], "normals": batch["normals"],
                          "wrench": batch["wrench"]}
                losses = self.loss_fn(pred, target)
            elif task == "texture":
                pred = self.model(batch["image"])
                losses = self.loss_fn(pred, batch["label"])
            elif task == "multimodal":
                pred = self.model(batch["tactile"], batch["visual"])
                target = {
                    "object_label": batch["object_label"],
                    "material_label": batch["material_label"],
                    "grasp_quality": batch["grasp_quality"],
                }
                losses = self.loss_fn(pred, target)
            elif task == "stability":
                pred = self.model(batch["tactile_seq"], batch["wrench_seq"])
                target = {
                    "stability_score": batch["stability_score"],
                    "stability_class": batch["stability_class"],
                    "time_to_failure": batch["time_to_failure"],
                    "corrective_action": batch["corrective_action"],
                }
                losses = self.loss_fn(pred, target)
            elif task == "slip":
                pred = self.model(
                    batch["tactile_seq"], batch["force"],
                    batch["target_force"],
                    torch.zeros_like(batch["force"])
                )
                losses = {"total": nn.functional.cross_entropy(
                    pred["slip_class"], batch["slip_class"]
                )}
            elif task == "exploration":
                pred = self.model(batch["visual"])
                target = {
                    "true_quality": batch["true_quality"],
                    "safety_labels": batch["safety_labels"],
                }
                losses = self.loss_fn(pred, target)
            else:
                raise ValueError(f"Unknown task: {task}")

            losses["total"].backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
            self.optimizer.step()

            for k, v in losses.items():
                if k not in total_losses:
                    total_losses[k] = 0.0
                total_losses[k] += v.item()
            n_batches += 1

        return {k: v / max(n_batches, 1) for k, v in total_losses.items()}

    @torch.no_grad()
    def evaluate(self, dataloader: DataLoader, task: str) -> Dict[str, float]:
        self.model.eval()
        total_losses = {}
        n_batches = 0
        correct = 0
        total_samples = 0

        for batch in dataloader:
            batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                     for k, v in batch.items()}

            if task == "contact":
                pred = self.model(batch["image"])
                target = {"depth": batch["depth"], "normals": batch["normals"],
                          "wrench": batch["wrench"]}
                losses = self.loss_fn(pred, target)
            elif task == "texture":
                pred = self.model(batch["image"])
                losses = self.loss_fn(pred, batch["label"])
                correct += (pred["logits"].argmax(1) == batch["label"]).sum().item()
                total_samples += batch["label"].shape[0]
            elif task == "multimodal":
                pred = self.model(batch["tactile"], batch["visual"])
                target = {
                    "object_label": batch["object_label"],
                    "material_label": batch["material_label"],
                    "grasp_quality": batch["grasp_quality"],
                }
                losses = self.loss_fn(pred, target)
                correct += (pred["object_logits"].argmax(1) == batch["object_label"]).sum().item()
                total_samples += batch["object_label"].shape[0]
            elif task == "stability":
                pred = self.model(batch["tactile_seq"], batch["wrench_seq"])
                target = {
                    "stability_score": batch["stability_score"],
                    "stability_class": batch["stability_class"],
                    "time_to_failure": batch["time_to_failure"],
                    "corrective_action": batch["corrective_action"],
                }
                losses = self.loss_fn(pred, target)
            elif task == "slip":
                pred = self.model(
                    batch["tactile_seq"], batch["force"],
                    batch["target_force"],
                    torch.zeros_like(batch["force"])
                )
                losses = {"total": nn.functional.cross_entropy(
                    pred["slip_class"], batch["slip_class"]
                )}
            elif task == "exploration":
                pred = self.model(batch["visual"])
                target = {
                    "true_quality": batch["true_quality"],
                    "safety_labels": batch["safety_labels"],
                }
                losses = self.loss_fn(pred, target)
            else:
                raise ValueError(f"Unknown task: {task}")

            for k, v in losses.items():
                if k not in total_losses:
                    total_losses[k] = 0.0
                total_losses[k] += v.item()
            n_batches += 1

        result = {k: v / max(n_batches, 1) for k, v in total_losses.items()}
        if total_samples > 0:
            result["accuracy"] = correct / total_samples
        return result

    def fit(self, train_loader: DataLoader, val_loader: DataLoader,
            task: str, epochs: int = 100, save_dir: str = "results/models/"):
        """Full training loop with early stopping."""
        os.makedirs(save_dir, exist_ok=True)
        patience = 15
        no_improve = 0

        for epoch in range(epochs):
            train_metrics = self.train_epoch(train_loader, task)
            val_metrics = self.evaluate(val_loader, task)

            self.metrics.update("train", train_metrics, epoch)
            self.metrics.update("val", val_metrics, epoch)

            if self.scheduler:
                self.scheduler.step()

            # Early stopping
            if val_metrics["total"] < self.best_val_loss:
                self.best_val_loss = val_metrics["total"]
                no_improve = 0
                torch.save(self.model.state_dict(),
                           os.path.join(save_dir, f"{task}_best.pth"))
            else:
                no_improve += 1
                if no_improve >= patience:
                    break

            if epoch % 10 == 0:
                print(f"[{task}] Epoch {epoch}: "
                      f"train_loss={train_metrics['total']:.4f} "
                      f"val_loss={val_metrics['total']:.4f}")

        self.metrics.save(f"{task}_metrics.json")
        return self.metrics.history


def build_and_train_all(device: str = "cpu", quick_test: bool = True):
    """Build and train all modules. Set quick_test=True for validation."""
    from src.models.contact_net import ContactNet, ContactLoss
    from src.models.texture_cnn import TextureCNN, TextureLoss
    from src.models.cross_modal_transformer import CrossModalTransformer, MultiModalLoss
    from src.models.grasp_stability_net import GraspStabilityNet, StabilityLoss
    from src.models.slip_detector import SlipDetector
    from src.models.exploratory_grasping import ExploratoryGraspingPolicy, ExplorationLoss

    results = {}
    n_samples = 64 if quick_test else 5000
    epochs = 3 if quick_test else 100
    batch_size = 8 if quick_test else 32

    modules = [
        ("contact", ContactNet, ContactLoss, {}),
        ("texture", TextureCNN, TextureLoss, {}),
        ("multimodal", CrossModalTransformer, MultiModalLoss, {}),
        ("stability", GraspStabilityNet, StabilityLoss, {}),
        ("slip", SlipDetector, None, {}),
        ("exploration", ExploratoryGraspingPolicy, ExplorationLoss, {}),
    ]

    for task, ModelClass, LossClass, kwargs in modules:
        print(f"\n{'='*60}")
        print(f"Training Module: {task}")
        print(f"{'='*60}")

        dataset = TactileDataset("", task=task)
        dataset._generate_synthetic(n_samples)

        train_size = int(0.8 * len(dataset))
        val_size = len(dataset) - train_size
        train_ds, val_ds = random_split(dataset, [train_size, val_size])

        train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_ds, batch_size=batch_size)

        model = ModelClass(**kwargs)
        loss_fn = LossClass() if LossClass else nn.CrossEntropyLoss()
        optimizer = optim.AdamW(model.parameters(), lr=1e-4, weight_decay=1e-5)
        scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

        trainer = Trainer(model, loss_fn, optimizer, scheduler,
                          device=device, log_dir="logs/")
        history = trainer.fit(train_loader, val_loader, task,
                              epochs=epochs, save_dir="results/models/")

        # Compute parameter count
        n_params = sum(p.numel() for p in model.parameters())
        results[task] = {
            "parameters": n_params,
            "best_val_loss": trainer.best_val_loss,
            "epochs_trained": len(history.get("train", [])),
        }
        print(f"  Parameters: {n_params:,}")
        print(f"  Best Val Loss: {trainer.best_val_loss:.4f}")

    return results


if __name__ == "__main__":
    results = build_and_train_all(device="cpu", quick_test=True)
    print("\n" + "=" * 60)
    print("Training Summary")
    print("=" * 60)
    for task, info in results.items():
        print(f"  {task}: {info['parameters']:,} params, "
              f"val_loss={info['best_val_loss']:.4f}")
