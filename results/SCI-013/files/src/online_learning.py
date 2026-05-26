"""Online learning and concept-drift adaptation for BCI systems."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Deque, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score



def _inv_sqrt(matrix: np.ndarray) -> np.ndarray:
    eigvals, eigvecs = np.linalg.eigh(matrix)
    eigvals = np.clip(eigvals, 1e-8, None)
    return eigvecs @ np.diag(eigvals ** -0.5) @ eigvecs.T


@dataclass
class ConceptDriftDetector:
    """Approximate ADWIN detector using adaptive window hypothesis testing."""

    delta: float = 0.002
    min_window: int = 20

    def __post_init__(self) -> None:
        self.window: List[float] = []
        self.drift_points: List[int] = []
        self.t = 0

    def update(self, value: float) -> bool:
        self.t += 1
        self.window.append(float(value))
        if len(self.window) < 2 * self.min_window:
            return False
        values = np.asarray(self.window, dtype=float)
        for cut in range(self.min_window, len(values) - self.min_window + 1, max(self.min_window // 2, 1)):
            left = values[:cut]
            right = values[cut:]
            diff = abs(left.mean() - right.mean())
            var = values.var() + 1e-8
            eps = np.sqrt(2 * var * np.log(2 / self.delta) * (1 / len(left) + 1 / len(right)))
            eps += 2 * np.log(2 / self.delta) * (1 / len(left) + 1 / len(right)) / 3
            if diff > eps:
                self.window = self.window[cut:]
                self.drift_points.append(self.t)
                return True
        return False


@dataclass
class OnlineLearner:
    """Incremental linear learner with replay buffer and Euclidean alignment."""

    classes: Sequence[int]
    replay_size: int = 256
    use_alignment: bool = True

    def __post_init__(self) -> None:
        self.model = SGDClassifier(loss="log_loss", learning_rate="optimal", random_state=7)
        self.replay_buffer: Deque[Tuple[np.ndarray, int]] = deque(maxlen=self.replay_size)
        self.reference_cov: Optional[np.ndarray] = None
        self.is_fitted = False
        self.statistics: Dict[str, List[float]] = {"accuracy": [], "buffer_fill": []}

    def fit_baseline(self, x: np.ndarray, y: np.ndarray) -> "OnlineLearner":
        x_aligned = self.euclidean_alignment(x)
        self.model.partial_fit(x_aligned, y, classes=np.asarray(self.classes))
        self.is_fitted = True
        return self

    def euclidean_alignment(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if x.ndim == 3:
            covs = np.stack([trial @ trial.T / max(trial.shape[-1] - 1, 1) for trial in x], axis=0)
            self.reference_cov = covs.mean(axis=0) if self.reference_cov is None else 0.95 * self.reference_cov + 0.05 * covs.mean(axis=0)
            whitening = _inv_sqrt(self.reference_cov)
            return np.asarray([(whitening @ trial).reshape(-1) for trial in x])
        return x

    def update(self, x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y)
        if x.ndim == 3 and self.use_alignment:
            features = self.euclidean_alignment(x)
        else:
            features = x.reshape(len(x), -1)
        if not self.is_fitted:
            self.fit_baseline(features, y)
        predictions = self.model.predict(features)
        accuracy = accuracy_score(y, predictions)
        for feat, label in zip(features, y):
            self.replay_buffer.append((feat, int(label)))
        replay = list(self.replay_buffer)
        if replay:
            replay_x = np.stack([item[0] for item in replay], axis=0)
            replay_y = np.asarray([item[1] for item in replay])
            self.model.partial_fit(replay_x, replay_y)
        self.statistics["accuracy"].append(float(accuracy))
        self.statistics["buffer_fill"].append(float(len(self.replay_buffer)))
        return {"accuracy": float(accuracy), "buffer_fill": float(len(self.replay_buffer))}

    def predict(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if x.ndim == 3 and self.use_alignment:
            features = self.euclidean_alignment(x)
        else:
            features = x.reshape(len(x), -1)
        return self.model.predict(features)

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if x.ndim == 3 and self.use_alignment:
            features = self.euclidean_alignment(x)
        else:
            features = x.reshape(len(x), -1)
        if hasattr(self.model, "decision_function"):
            scores = self.model.decision_function(features)
            if scores.ndim == 1:
                scores = np.column_stack([-scores, scores])
            scores = scores - scores.max(axis=1, keepdims=True)
            probabilities = np.exp(scores)
            probabilities /= probabilities.sum(axis=1, keepdims=True) + 1e-8
            return probabilities
        probabilities = np.nan_to_num(self.model.predict_proba(features), nan=1.0 / len(self.classes), posinf=1.0, neginf=0.0)
        probabilities /= probabilities.sum(axis=1, keepdims=True) + 1e-8
        return probabilities

    def get_statistics(self) -> Dict[str, List[float]]:
        return self.statistics


@dataclass
class EnsembleAdapter:
    """Dynamic ensemble with adaptive weighted voting."""

    models: List[OnlineLearner] = field(default_factory=list)
    weights: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.weights is None:
            self.weights = np.ones(len(self.models), dtype=float) if self.models else np.empty(0)

    def add_model(self, model: OnlineLearner, weight: float = 1.0) -> None:
        self.models.append(model)
        self.weights = np.append(self.weights, weight)

    def predict(self, x: np.ndarray) -> np.ndarray:
        if not self.models:
            raise RuntimeError("No models available in the ensemble.")
        probabilities = [model.predict_proba(x) for model in self.models]
        weighted = np.tensordot(self.weights / np.sum(self.weights), np.stack(probabilities, axis=0), axes=(0, 0))
        return weighted.argmax(axis=1)

    def update_weights(self, x: np.ndarray, y: np.ndarray, smoothing: float = 0.8) -> np.ndarray:
        scores = []
        for model in self.models:
            preds = model.predict(x)
            scores.append(accuracy_score(y, preds) + 1e-3)
        scores = np.asarray(scores)
        normalized = scores / scores.sum()
        self.weights = smoothing * self.weights + (1.0 - smoothing) * normalized
        return self.weights
