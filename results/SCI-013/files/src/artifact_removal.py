"""Real-time EEG artifact removal utilities."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Deque, List, Optional

import numpy as np
from scipy import linalg, stats


ArrayLike = np.ndarray



def _ensure_2d(data: ArrayLike) -> ArrayLike:
    data = np.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data[None, :]
    if data.ndim != 2:
        raise ValueError("Expected data with shape (channels, samples).")
    return data



def _covariance(data: ArrayLike) -> ArrayLike:
    centered = data - data.mean(axis=1, keepdims=True)
    cov = centered @ centered.T / max(centered.shape[1] - 1, 1)
    cov += 1e-6 * np.eye(centered.shape[0])
    return cov



def _matrix_inv_sqrt(matrix: ArrayLike) -> ArrayLike:
    eigvals, eigvecs = np.linalg.eigh(matrix)
    eigvals = np.clip(eigvals, 1e-8, None)
    return eigvecs @ np.diag(eigvals ** -0.5) @ eigvecs.T


@dataclass
class RealTimeASR:
    """Sliding-window Artifact Subspace Reconstruction for streaming EEG."""

    sfreq: float = 250.0
    window_size: float = 0.5
    cutoff: float = 3.0
    history: Deque[ArrayLike] = field(default_factory=deque)

    def __post_init__(self) -> None:
        self.window_samples = max(8, int(self.window_size * self.sfreq))
        self.buffer: Deque[ArrayLike] = deque(maxlen=8)
        self.reference_cov: Optional[ArrayLike] = None
        self.reference_mean: Optional[ArrayLike] = None
        self.reference_std: Optional[ArrayLike] = None

    def fit_baseline(self, data: ArrayLike) -> "RealTimeASR":
        data = _ensure_2d(data)
        windows = []
        for start in range(0, max(data.shape[1] - self.window_samples + 1, 1), max(self.window_samples // 2, 1)):
            stop = min(start + self.window_samples, data.shape[1])
            window = data[:, start:stop]
            if window.shape[1] >= 8:
                windows.append(_covariance(window))
        if not windows:
            windows = [_covariance(data)]
        eigvals = np.vstack([np.linalg.eigvalsh(cov) for cov in windows])
        self.reference_cov = np.mean(windows, axis=0)
        self.reference_mean = eigvals.mean(axis=0)
        self.reference_std = eigvals.std(axis=0) + 1e-6
        return self

    def _clean(self, chunk: ArrayLike, cov: ArrayLike) -> ArrayLike:
        if self.reference_cov is None or self.reference_mean is None or self.reference_std is None:
            return chunk
        eigvals, eigvecs = np.linalg.eigh(cov)
        threshold = self.reference_mean + self.cutoff * self.reference_std
        clipped = np.minimum(eigvals, threshold)
        scaling = np.sqrt(clipped / np.clip(eigvals, 1e-8, None))
        projector = eigvecs @ np.diag(scaling) @ eigvecs.T
        return projector @ chunk

    def process(self, chunk: ArrayLike) -> ArrayLike:
        chunk = _ensure_2d(chunk)
        self.buffer.append(chunk)
        stacked = np.concatenate(list(self.buffer), axis=1)
        if stacked.shape[1] > self.window_samples:
            stacked = stacked[:, -self.window_samples :]
        cleaned = self._clean(chunk, _covariance(stacked))
        self.history.append(cleaned)
        return cleaned

    def get_clean_signal(self) -> ArrayLike:
        if not self.history:
            return np.empty((0, 0))
        return np.concatenate(list(self.history), axis=1)


@dataclass
class OnlineICA:
    """Incremental ICA using a natural-gradient update in whitened space."""

    n_components: Optional[int] = None
    learning_rate: float = 0.05
    kurtosis_threshold: float = 5.0
    history: List[ArrayLike] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.mean_: Optional[ArrayLike] = None
        self.whitening_: Optional[ArrayLike] = None
        self.dewhitening_: Optional[ArrayLike] = None
        self.unmixing_: Optional[ArrayLike] = None

    def fit_baseline(self, data: ArrayLike) -> "OnlineICA":
        data = _ensure_2d(data)
        n_channels = data.shape[0]
        self.n_components = self.n_components or n_channels
        self.mean_ = data.mean(axis=1, keepdims=True)
        cov = _covariance(data)
        eigvals, eigvecs = np.linalg.eigh(cov)
        idx = np.argsort(eigvals)[::-1][: self.n_components]
        eigvals = np.clip(eigvals[idx], 1e-8, None)
        eigvecs = eigvecs[:, idx]
        self.whitening_ = np.diag(eigvals ** -0.5) @ eigvecs.T
        self.dewhitening_ = eigvecs @ np.diag(eigvals ** 0.5)
        random_matrix = np.random.default_rng(7).standard_normal((self.n_components, self.n_components))
        q, _ = np.linalg.qr(random_matrix)
        self.unmixing_ = q
        return self

    def _decorrelate(self, matrix: ArrayLike) -> ArrayLike:
        s, u = np.linalg.eigh(matrix @ matrix.T)
        s = np.clip(s, 1e-8, None)
        return (u @ np.diag(s ** -0.5) @ u.T) @ matrix

    def process(self, chunk: ArrayLike) -> ArrayLike:
        chunk = _ensure_2d(chunk)
        if self.mean_ is None or self.whitening_ is None or self.dewhitening_ is None or self.unmixing_ is None:
            self.fit_baseline(chunk)
        assert self.mean_ is not None and self.whitening_ is not None
        assert self.dewhitening_ is not None and self.unmixing_ is not None
        self.mean_ = 0.995 * self.mean_ + 0.005 * chunk.mean(axis=1, keepdims=True)
        centered = chunk - self.mean_
        white = self.whitening_ @ centered
        y = self.unmixing_ @ white
        g = np.tanh(y)
        update = (np.eye(self.unmixing_.shape[0]) - (g @ y.T) / max(y.shape[1], 1)) @ self.unmixing_
        self.unmixing_ = self._decorrelate(self.unmixing_ + self.learning_rate * update)
        sources = self.unmixing_ @ white
        kurtosis = np.abs(stats.kurtosis(sources, axis=1, fisher=False, bias=False, nan_policy="omit"))
        kurtosis = np.nan_to_num(kurtosis, nan=0.0)
        variances = sources.var(axis=1)
        variance_limit = variances.mean() + 2.0 * variances.std()
        keep = (kurtosis < self.kurtosis_threshold) & (variances < variance_limit)
        cleaned_sources = sources * keep[:, None]
        reconstruction = self.dewhitening_ @ np.linalg.pinv(self.unmixing_) @ cleaned_sources + self.mean_
        self.history.append(reconstruction)
        return reconstruction

    def get_clean_signal(self) -> ArrayLike:
        if not self.history:
            return np.empty((0, 0))
        return np.concatenate(self.history, axis=1)


@dataclass
class ArtifactRemovalPipeline:
    """Combined ASR + incremental ICA pipeline for real-time EEG cleaning."""

    sfreq: float = 250.0

    def __post_init__(self) -> None:
        self.asr = RealTimeASR(sfreq=self.sfreq)
        self.ica = OnlineICA()
        self.history: List[ArrayLike] = []

    def fit_baseline(self, data: ArrayLike) -> "ArtifactRemovalPipeline":
        data = _ensure_2d(data)
        self.asr.fit_baseline(data)
        baseline_asr = self.asr.process(data)
        self.ica.fit_baseline(baseline_asr)
        self.history.clear()
        return self

    def process(self, chunk: ArrayLike) -> ArrayLike:
        chunk = _ensure_2d(chunk)
        cleaned_asr = self.asr.process(chunk)
        cleaned = self.ica.process(cleaned_asr)
        self.history.append(cleaned)
        return cleaned

    def get_clean_signal(self) -> ArrayLike:
        if not self.history:
            return np.empty((0, 0))
        return np.concatenate(self.history, axis=1)
