"""Main real-time BCI processing pipeline."""

from __future__ import annotations

import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import signal

try:
    from artifact_removal import ArtifactRemovalPipeline
except ImportError:  # pragma: no cover
    from .artifact_removal import ArtifactRemovalPipeline


@dataclass
class RealTimeBCIPipeline:
    """End-to-end real-time BCI pipeline with synthetic data acquisition."""

    paradigm: str = "MI"
    n_channels: int = 22
    sfreq: float = 250.0

    def __post_init__(self) -> None:
        self.artifact_pipeline = ArtifactRemovalPipeline(sfreq=self.sfreq)
        self.latencies: List[float] = []
        self.accuracy_log: List[float] = []
        self.predictions: List[int] = []
        self.targets: List[int] = []
        self.rng = np.random.default_rng(11)
        baseline = 0.05 * self.rng.standard_normal((self.n_channels, int(self.sfreq * 2)))
        self.artifact_pipeline.fit_baseline(baseline)

    def simulate_eeg_chunk(self, chunk_size: int = 250) -> Tuple[np.ndarray, int]:
        t = np.arange(chunk_size) / self.sfreq
        if self.paradigm.upper() == "MI":
            label = int(self.rng.integers(0, 4))
            base_freqs = [10, 12, 18, 22]
            groups = [(0, 6), (6, 12), (12, 17), (17, 22)]
            chunk = 0.15 * self.rng.standard_normal((self.n_channels, chunk_size))
            start, stop = groups[label]
            freq = base_freqs[label]
            chunk[start:stop] += 0.8 * np.sin(2 * np.pi * freq * t)
        elif self.paradigm.upper() == "P300":
            label = int(self.rng.random() > 0.75)
            chunk = 0.12 * self.rng.standard_normal((self.n_channels, chunk_size))
            if label:
                erp = np.exp(-0.5 * ((t - 0.32) / 0.05) ** 2)
                chunk += np.linspace(1.0, 0.6, self.n_channels)[:, None] * erp
        else:
            freqs = [10, 12, 15]
            label = int(self.rng.integers(0, len(freqs)))
            chunk = 0.1 * self.rng.standard_normal((self.n_channels, chunk_size))
            chunk += 0.8 * np.sin(2 * np.pi * freqs[label] * t)
        artifact = np.zeros_like(chunk)
        if self.rng.random() > 0.8:
            artifact[:, chunk_size // 3 : chunk_size // 3 + 10] += 1.5
        return chunk + artifact, label

    def preprocess(self, chunk: np.ndarray) -> np.ndarray:
        detrended = signal.detrend(chunk, axis=-1, type="linear")
        sos = signal.butter(4, [1, 40], btype="bandpass", fs=self.sfreq, output="sos")
        return signal.sosfiltfilt(sos, detrended, axis=-1)

    def classify(self, chunk: np.ndarray) -> int:
        paradigm = self.paradigm.upper()
        if paradigm == "MI":
            band_powers = []
            for band in [(8, 12), (12, 20), (20, 28), (28, 36)]:
                sos = signal.butter(2, band, btype="bandpass", fs=self.sfreq, output="sos")
                filtered = signal.sosfiltfilt(sos, chunk, axis=-1)
                band_powers.append(filtered.var(axis=1).mean())
            return int(np.argmax(band_powers)) % 4
        if paradigm == "P300":
            window = chunk[:, int(0.25 * self.sfreq) : int(0.45 * self.sfreq)].mean()
            return int(window > 0.15)
        freqs = np.fft.rfftfreq(chunk.shape[1], 1 / self.sfreq)
        spectrum = np.abs(np.fft.rfft(chunk.mean(axis=0)))
        target_freqs = np.array([10, 12, 15])
        scores = [spectrum[np.argmin(np.abs(freqs - freq))] for freq in target_freqs]
        return int(np.argmax(scores))

    def process_once(self, chunk_size: int = 250) -> Dict[str, float | int]:
        raw, target = self.simulate_eeg_chunk(chunk_size=chunk_size)
        start = time.perf_counter()
        clean = self.artifact_pipeline.process(self.preprocess(raw))
        pred = self.classify(clean)
        latency_ms = (time.perf_counter() - start) * 1000.0
        self.latencies.append(latency_ms)
        self.targets.append(target)
        self.predictions.append(pred)
        self.accuracy_log.append(float(pred == target))
        return {"target": target, "prediction": pred, "latency_ms": latency_ms, "correct": int(pred == target)}

    def run(self, n_chunks: int = 50, chunk_size: int = 250) -> Dict[str, object]:
        results = [self.process_once(chunk_size=chunk_size) for _ in range(n_chunks)]
        return {
            "results": results,
            "mean_latency_ms": float(np.mean(self.latencies)),
            "accuracy": float(np.mean(self.accuracy_log)),
            "latencies": self.latencies,
        }
