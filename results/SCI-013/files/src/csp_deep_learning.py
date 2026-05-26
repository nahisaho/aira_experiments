"""Filter-bank CSP and deep learning models for motor imagery EEG."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch
from scipy import linalg, signal
from sklearn.metrics import accuracy_score, cohen_kappa_score
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


BANDS: List[Tuple[float, float]] = [(4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 28), (28, 32)]



def bandpass_filter(data: np.ndarray, sfreq: float, band: Tuple[float, float], order: int = 4) -> np.ndarray:
    sos = signal.butter(order, band, btype="bandpass", fs=sfreq, output="sos")
    return signal.sosfiltfilt(sos, data, axis=-1)



def _trial_covariance(trial: np.ndarray) -> np.ndarray:
    cov = trial @ trial.T / max(trial.shape[-1] - 1, 1)
    cov /= np.trace(cov) + 1e-8
    cov += 1e-6 * np.eye(cov.shape[0])
    return cov


@dataclass
class FilterBankCSP:
    """One-vs-rest filter-bank common spatial patterns."""

    sfreq: float = 250.0
    bands: Sequence[Tuple[float, float]] = tuple(BANDS)
    n_components: int = 4

    def __post_init__(self) -> None:
        self.filters_: Dict[Tuple[int, int], np.ndarray] = {}
        self.patterns_: Dict[Tuple[int, int], np.ndarray] = {}
        self.classes_: Optional[np.ndarray] = None

    def fit(self, x: np.ndarray, y: np.ndarray) -> "FilterBankCSP":
        x = np.asarray(x, dtype=float)
        y = np.asarray(y)
        self.classes_ = np.unique(y)
        for band_idx, band in enumerate(self.bands):
            filtered = bandpass_filter(x, self.sfreq, band)
            for cls in self.classes_:
                cov_class = np.mean([_trial_covariance(trial) for trial in filtered[y == cls]], axis=0)
                cov_rest = np.mean([_trial_covariance(trial) for trial in filtered[y != cls]], axis=0)
                eigvals, eigvecs = linalg.eigh(cov_class, cov_class + cov_rest)
                order = np.argsort(eigvals)
                eigvecs = eigvecs[:, order]
                half = max(self.n_components // 2, 1)
                filters = np.concatenate([eigvecs[:, :half].T, eigvecs[:, -half:].T], axis=0)
                self.filters_[(band_idx, int(cls))] = filters
                self.patterns_[(band_idx, int(cls))] = np.linalg.pinv(filters).T
        return self

    def transform(self, x: np.ndarray) -> np.ndarray:
        if self.classes_ is None:
            raise RuntimeError("FilterBankCSP must be fitted before transform.")
        x = np.asarray(x, dtype=float)
        features: List[np.ndarray] = []
        for band_idx, band in enumerate(self.bands):
            filtered = bandpass_filter(x, self.sfreq, band)
            for cls in self.classes_:
                filters = self.filters_[(band_idx, int(cls))]
                projected = np.einsum("fc,nct->nft", filters, filtered)
                variance = np.log(np.var(projected, axis=-1) + 1e-8)
                features.append(variance)
        return np.concatenate(features, axis=1)

    def fit_transform(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return self.fit(x, y).transform(x)

    def get_patterns(self) -> Dict[Tuple[int, int], np.ndarray]:
        return self.patterns_


class CSPNet(nn.Module):
    """Deep motor-imagery classifier with learnable CSP filters and temporal modeling."""

    def __init__(
        self,
        n_channels: int,
        n_classes: int = 4,
        n_csp: int = 8,
        conv_channels: int = 32,
        lstm_hidden: int = 64,
        dropout: float = 0.3,
    ) -> None:
        super().__init__()
        self.spatial_filters = nn.Parameter(torch.randn(n_csp, n_channels) * 0.1)
        self.batch_norm = nn.BatchNorm1d(n_csp)
        self.conv = nn.Conv1d(n_csp, conv_channels, kernel_size=7, padding=3)
        self.dropout = nn.Dropout(dropout)
        self.lstm = nn.LSTM(conv_channels, lstm_hidden, batch_first=True, bidirectional=True)
        self.attention = nn.Linear(lstm_hidden * 2, 1)
        self.fc = nn.Linear(lstm_hidden * 2, n_classes)
        self._attention_weights: Optional[torch.Tensor] = None

    def forward(self, x: torch.Tensor, return_attention: bool = False) -> torch.Tensor | Tuple[torch.Tensor, torch.Tensor]:
        x = torch.einsum("oc,bct->bot", self.spatial_filters, x)
        x = self.batch_norm(x)
        x = torch.nn.functional.gelu(self.conv(x))
        x = self.dropout(x)
        x = x.transpose(1, 2)
        seq, _ = self.lstm(x)
        weights = torch.softmax(self.attention(seq), dim=1)
        self._attention_weights = weights.detach()
        pooled = torch.sum(weights * seq, dim=1)
        logits = self.fc(pooled)
        if return_attention:
            return logits, weights.squeeze(-1)
        return logits

    @property
    def attention_weights(self) -> Optional[torch.Tensor]:
        return self._attention_weights


def _make_loader(x: np.ndarray, y: np.ndarray, batch_size: int = 32, shuffle: bool = True) -> DataLoader:
    dataset = TensorDataset(torch.tensor(x, dtype=torch.float32), torch.tensor(y, dtype=torch.long))
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)



def train_cspnet(
    model: nn.Module,
    train_x: np.ndarray,
    train_y: np.ndarray,
    val_x: Optional[np.ndarray] = None,
    val_y: Optional[np.ndarray] = None,
    epochs: int = 8,
    lr: float = 1e-3,
    batch_size: int = 32,
    device: Optional[str] = None,
) -> Dict[str, List[float]]:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    loader = _make_loader(train_x, train_y, batch_size=batch_size, shuffle=True)
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    history: Dict[str, List[float]] = {"loss": [], "val_accuracy": []}
    for _ in range(epochs):
        model.train()
        losses = []
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            optimizer.zero_grad()
            logits = model(batch_x)
            loss = criterion(logits, batch_y)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.item()))
        history["loss"].append(float(np.mean(losses)))
        if val_x is not None and val_y is not None:
            metrics = evaluate_classifier(model, val_x, val_y, device=device)
            history["val_accuracy"].append(metrics["accuracy"])
    return history


@torch.no_grad()
def evaluate_classifier(model: nn.Module, x: np.ndarray, y: np.ndarray, device: Optional[str] = None) -> Dict[str, np.ndarray | float]:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model.eval().to(device)
    tensor_x = torch.tensor(x, dtype=torch.float32, device=device)
    logits = model(tensor_x)
    predictions = torch.argmax(logits, dim=1).cpu().numpy()
    accuracy = accuracy_score(y, predictions)
    kappa = cohen_kappa_score(y, predictions)
    return {"accuracy": float(accuracy), "kappa": float(kappa), "predictions": predictions}
