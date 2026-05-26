"""P300 decoding, transfer learning, and adaptive online classification."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch
from sklearn.metrics import accuracy_score, cohen_kappa_score
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


class EEGNetP300(nn.Module):
    """Compact EEGNet-inspired network for P300 detection."""

    def __init__(self, n_channels: int = 22, n_samples: int = 200, dropout: float = 0.3) -> None:
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(1, 8, kernel_size=(1, 31), padding=(0, 15), bias=False),
            nn.BatchNorm2d(8),
            nn.Conv2d(8, 16, kernel_size=(n_channels, 1), groups=8, bias=False),
            nn.BatchNorm2d(16),
            nn.ELU(),
            nn.AvgPool2d((1, 4)),
            nn.Dropout(dropout),
            nn.Conv2d(16, 16, kernel_size=(1, 15), padding=(0, 7), groups=16, bias=False),
            nn.Conv2d(16, 32, kernel_size=(1, 1), bias=False),
            nn.BatchNorm2d(32),
            nn.ELU(),
            nn.AvgPool2d((1, 4)),
            nn.Dropout(dropout),
        )
        reduced = max(n_samples // 16, 1)
        self.classifier = nn.Linear(32 * reduced, 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.unsqueeze(1)
        x = self.features(x)
        x = x.flatten(start_dim=1)
        return self.classifier(x)



def _loader(x: np.ndarray, y: np.ndarray, batch_size: int = 32, shuffle: bool = True) -> DataLoader:
    return DataLoader(TensorDataset(torch.tensor(x, dtype=torch.float32), torch.tensor(y, dtype=torch.long)), batch_size=batch_size, shuffle=shuffle)


@dataclass
class TransferLearningP300:
    """Transfer-learning helper for subject-to-subject P300 adaptation."""

    model: EEGNetP300
    device: Optional[str] = None

    def __post_init__(self) -> None:
        self.device = self.device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)

    def fit_source(self, x: np.ndarray, y: np.ndarray, epochs: int = 6, lr: float = 1e-3) -> None:
        train_network(self.model, x, y, epochs=epochs, lr=lr, device=self.device)

    def adapt_target(
        self,
        x_target: np.ndarray,
        y_target: Optional[np.ndarray] = None,
        pseudo_threshold: float = 0.85,
        epochs: int = 3,
        lr: float = 5e-4,
    ) -> np.ndarray:
        if y_target is None:
            probabilities = predict_proba(self.model, x_target, device=self.device)
            confidence = probabilities.max(axis=1)
            pseudo_y = probabilities.argmax(axis=1)
            mask = confidence >= pseudo_threshold
            if not np.any(mask):
                mask = np.argsort(confidence)[-max(1, len(confidence) // 4) :]
                x_adapt = x_target[mask]
                y_adapt = pseudo_y[mask]
            else:
                x_adapt = x_target[mask]
                y_adapt = pseudo_y[mask]
        else:
            x_adapt = x_target
            y_adapt = y_target
        if len(x_adapt):
            train_network(self.model, x_adapt, y_adapt, epochs=epochs, lr=lr, device=self.device)
        return predict_proba(self.model, x_target, device=self.device)

    def evaluate(self, x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        preds = predict_proba(self.model, x, device=self.device).argmax(axis=1)
        return {"accuracy": float(accuracy_score(y, preds)), "kappa": float(cohen_kappa_score(y, preds))}


@dataclass
class AdaptiveP300Classifier:
    """Online P300 classifier with EMA adaptation and pseudo-labeling."""

    model: EEGNetP300
    ema_decay: float = 0.98
    device: Optional[str] = None

    def __post_init__(self) -> None:
        self.device = self.device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=5e-4)
        self.criterion = nn.CrossEntropyLoss()
        self.shadow: Dict[str, torch.Tensor] = {k: v.detach().clone() for k, v in self.model.state_dict().items()}

    def _update_shadow(self) -> None:
        for key, value in self.model.state_dict().items():
            self.shadow[key] = self.ema_decay * self.shadow[key] + (1.0 - self.ema_decay) * value.detach()

    def adapt(self, x: np.ndarray, y: Optional[np.ndarray] = None, pseudo_threshold: float = 0.9) -> Dict[str, float]:
        x_tensor = torch.tensor(x, dtype=torch.float32, device=self.device)
        if y is None:
            probabilities = torch.softmax(self.model(x_tensor), dim=1)
            confidence, pseudo = torch.max(probabilities, dim=1)
            mask = confidence >= pseudo_threshold
            if not torch.any(mask):
                return {"adapted_samples": 0.0}
            x_tensor = x_tensor[mask]
            y_tensor = pseudo[mask]
        else:
            y_tensor = torch.tensor(y, dtype=torch.long, device=self.device)
        self.model.train()
        self.optimizer.zero_grad()
        logits = self.model(x_tensor)
        loss = self.criterion(logits, y_tensor)
        loss.backward()
        self.optimizer.step()
        self._update_shadow()
        return {"adapted_samples": float(len(x_tensor)), "loss": float(loss.item())}

    @torch.no_grad()
    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        state = self.model.state_dict()
        self.model.load_state_dict(self.shadow, strict=False)
        self.model.eval()
        x_tensor = torch.tensor(x, dtype=torch.float32, device=self.device)
        probabilities = torch.softmax(self.model(x_tensor), dim=1).cpu().numpy()
        self.model.load_state_dict(state, strict=False)
        return probabilities


@dataclass
class P300SpellerSimulation:
    """Synthetic 6x6 P300 speller stimulation and ERP generation."""

    sfreq: float = 250.0
    n_channels: int = 22

    def __post_init__(self) -> None:
        symbols = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890")
        self.matrix = np.array(symbols[:36]).reshape(6, 6)

    def find_symbol(self, symbol: str) -> Tuple[int, int]:
        location = np.argwhere(self.matrix == symbol.upper())
        if not len(location):
            raise ValueError(f"Unknown symbol: {symbol}")
        return tuple(map(int, location[0]))

    def flash_sequence(self, target_symbol: str, repetitions: int = 8) -> List[Tuple[int, int, int]]:
        target_row, target_col = self.find_symbol(target_symbol)
        rng = np.random.default_rng(42)
        sequence: List[Tuple[int, int, int]] = []
        for _ in range(repetitions):
            order = list(range(12))
            rng.shuffle(order)
            for event in order:
                if event < 6:
                    is_target = int(event == target_row)
                    sequence.append((event, -1, is_target))
                else:
                    col = event - 6
                    is_target = int(col == target_col)
                    sequence.append((-1, col, is_target))
        return sequence

    def generate_trials(self, target_symbol: str, repetitions: int = 8, n_samples: int = 200) -> Tuple[np.ndarray, np.ndarray]:
        sequence = self.flash_sequence(target_symbol, repetitions=repetitions)
        rng = np.random.default_rng(1)
        time = np.linspace(0, 0.8, n_samples, endpoint=False)
        p300 = np.exp(-0.5 * ((time - 0.32) / 0.06) ** 2)
        trials = []
        labels = []
        for _, _, label in sequence:
            trial = 0.15 * rng.standard_normal((self.n_channels, n_samples))
            if label:
                scalp_profile = np.linspace(1.1, 0.7, self.n_channels)[:, None]
                trial += 1.8 * scalp_profile * p300
            trials.append(trial)
            labels.append(label)
        return np.asarray(trials), np.asarray(labels)



def train_network(
    model: nn.Module,
    x: np.ndarray,
    y: np.ndarray,
    epochs: int = 6,
    lr: float = 1e-3,
    batch_size: int = 32,
    device: Optional[str] = None,
) -> List[float]:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    model.train()
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    loader = _loader(x, y, batch_size=batch_size, shuffle=True)
    losses: List[float] = []
    for _ in range(epochs):
        epoch_losses = []
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            optimizer.zero_grad()
            loss = criterion(model(batch_x), batch_y)
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.item()))
        losses.append(float(np.mean(epoch_losses)))
    return losses


@torch.no_grad()
def predict_proba(model: nn.Module, x: np.ndarray, device: Optional[str] = None) -> np.ndarray:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model.eval().to(device)
    x_tensor = torch.tensor(x, dtype=torch.float32, device=device)
    return torch.softmax(model(x_tensor), dim=1).cpu().numpy()
