"""
Module 4: Concept Drift Detection & Model Retraining Triggers
- ADWIN (Adaptive Windowing)
- Page-Hinkley Test
- DDM (Drift Detection Method)
- Automated retraining pipeline
"""
import numpy as np
from collections import deque
from dataclasses import dataclass
from typing import Optional


class ADWINDetector:
    """Adaptive Windowing drift detector."""

    def __init__(self, delta=0.002, max_window=5000):
        self.delta = delta
        self.max_window = max_window
        self.window = deque(maxlen=max_window)
        self.drift_points = []
        self.total = 0.0
        self.variance = 0.0
        self.width = 0

    def update(self, value: float) -> bool:
        self.window.append(value)
        self.width = len(self.window)

        if self.width < 10:
            return False

        return self._check_cut()

    def _check_cut(self) -> bool:
        data = np.array(self.window)
        n = len(data)

        for i in range(max(5, n // 10), n - max(5, n // 10)):
            left = data[:i]
            right = data[i:]
            n0, n1 = len(left), len(right)
            mu0, mu1 = np.mean(left), np.mean(right)

            m = 1.0 / (1.0/n0 + 1.0/n1)
            eps = np.sqrt(np.log(4.0 / self.delta) / (2.0 * m))

            if np.abs(mu0 - mu1) >= eps:
                # Remove older half
                for _ in range(i):
                    if self.window:
                        self.window.popleft()
                return True
        return False

    def detect_batch(self, data: np.ndarray) -> dict:
        drift_points = []
        for i, val in enumerate(data):
            if self.update(val):
                drift_points.append(i)

        self.drift_points = drift_points
        return {
            "drift_detected": len(drift_points) > 0,
            "drift_points": drift_points,
            "n_drifts": len(drift_points),
            "final_window_size": len(self.window),
        }


class PageHinkleyDetector:
    """Page-Hinkley test for drift detection."""

    def __init__(self, delta=0.005, threshold=50, alpha=0.9999):
        self.delta = delta
        self.threshold = threshold
        self.alpha = alpha
        self.sum_ = 0.0
        self.x_mean = 0.0
        self.n = 0
        self.m_t = 0.0
        self.M_t = 0.0

    def update(self, x: float) -> bool:
        self.n += 1
        self.x_mean = self.x_mean + (x - self.x_mean) / self.n
        self.m_t = self.alpha * self.m_t + (x - self.x_mean - self.delta)
        self.M_t = max(self.M_t, self.m_t)

        return (self.M_t - self.m_t) > self.threshold

    def detect_batch(self, data: np.ndarray) -> dict:
        drift_points = []
        for i, val in enumerate(data):
            if self.update(val):
                drift_points.append(i)
                self._reset()
        return {
            "drift_detected": len(drift_points) > 0,
            "drift_points": drift_points,
            "n_drifts": len(drift_points),
        }

    def _reset(self):
        self.sum_ = 0.0
        self.x_mean = 0.0
        self.n = 0
        self.m_t = 0.0
        self.M_t = 0.0


@dataclass
class RetrainingDecision:
    should_retrain: bool
    reason: str
    drift_type: str
    severity: float
    recommended_window: Optional[int] = None


class RetrainingTrigger:
    """Automated model retraining trigger based on drift detection."""

    def __init__(self, performance_threshold=0.1, drift_patience=3,
                 min_retrain_interval=100):
        self.perf_threshold = performance_threshold
        self.drift_patience = drift_patience
        self.min_interval = min_retrain_interval
        self.consecutive_drifts = 0
        self.last_retrain = 0
        self.performance_history = []

    def evaluate(self, step: int, current_perf: float, baseline_perf: float,
                 drift_detected: bool) -> RetrainingDecision:
        self.performance_history.append(current_perf)
        perf_drop = baseline_perf - current_perf

        if step - self.last_retrain < self.min_interval:
            return RetrainingDecision(False, "Too soon since last retrain",
                                      "none", 0.0)

        if drift_detected:
            self.consecutive_drifts += 1
        else:
            self.consecutive_drifts = max(0, self.consecutive_drifts - 1)

        # Sudden drift: large performance drop
        if perf_drop > self.perf_threshold * 2:
            self.last_retrain = step
            self.consecutive_drifts = 0
            return RetrainingDecision(
                True, f"Sudden performance drop: {perf_drop:.4f}",
                "sudden", perf_drop / baseline_perf,
                recommended_window=self.min_interval
            )

        # Gradual drift: sustained small drifts
        if self.consecutive_drifts >= self.drift_patience:
            self.last_retrain = step
            self.consecutive_drifts = 0
            return RetrainingDecision(
                True, f"Gradual drift: {self.consecutive_drifts} consecutive detections",
                "gradual", perf_drop / baseline_perf if baseline_perf > 0 else 0,
                recommended_window=self.min_interval * 2
            )

        # Performance below threshold
        if perf_drop > self.perf_threshold:
            self.last_retrain = step
            return RetrainingDecision(
                True, f"Performance below threshold: drop={perf_drop:.4f}",
                "degradation", perf_drop / baseline_perf,
                recommended_window=self.min_interval
            )

        return RetrainingDecision(False, "No retrain needed", "none", 0.0)
