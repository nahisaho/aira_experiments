"""
Module 6: Streaming Anomaly Detection Pipeline
- CERN/LIGO-scale data processing architecture
- Online anomaly detection with River
- Windowed batch processing
- Pipeline orchestration
"""
import numpy as np
import time
import json
from collections import deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Callable
from enum import Enum


class AlertSeverity(Enum):
    INFO = "info"
    WARNING = "warning"
    CRITICAL = "critical"


@dataclass
class Alert:
    timestamp: float
    severity: AlertSeverity
    source: str
    message: str
    score: float
    metadata: Dict = field(default_factory=dict)


@dataclass
class PipelineMetrics:
    total_processed: int = 0
    total_anomalies: int = 0
    total_alerts: int = 0
    processing_time_ms: float = 0.0
    throughput_per_sec: float = 0.0
    window_anomaly_rate: float = 0.0


class StreamingAnomalyPipeline:
    """End-to-end streaming anomaly detection pipeline for large-scale experiments."""

    def __init__(self, window_size=1000, alert_threshold=0.8,
                 batch_interval=100, n_features=10):
        self.window_size = window_size
        self.alert_threshold = alert_threshold
        self.batch_interval = batch_interval
        self.n_features = n_features

        # Sliding windows
        self.data_window = deque(maxlen=window_size)
        self.score_window = deque(maxlen=window_size)
        self.alert_buffer: List[Alert] = []

        # Online statistics
        self.running_mean = np.zeros(n_features)
        self.running_var = np.ones(n_features)
        self.n_seen = 0

        # Metrics
        self.metrics = PipelineMetrics()

        # Pipeline stages
        self.preprocessors: List[Callable] = []
        self.detectors: List[Callable] = []
        self.postprocessors: List[Callable] = []

    def add_preprocessor(self, fn: Callable):
        self.preprocessors.append(fn)

    def add_detector(self, fn: Callable):
        self.detectors.append(fn)

    def add_postprocessor(self, fn: Callable):
        self.postprocessors.append(fn)

    def process_event(self, event: np.ndarray) -> Dict:
        """Process a single event through the full pipeline."""
        start = time.time()

        # 1. Preprocessing
        processed = event.copy()
        for pp in self.preprocessors:
            processed = pp(processed)

        # 2. Update online statistics
        self._update_stats(processed)

        # 3. Anomaly detection
        scores = []
        for det in self.detectors:
            score = det(processed, self.running_mean, self.running_var)
            scores.append(score)

        combined_score = float(np.mean(scores)) if scores else 0.0

        # 4. Add to windows
        self.data_window.append(processed)
        self.score_window.append(combined_score)

        # 5. Check for alerts
        alert = None
        if combined_score > self.alert_threshold:
            alert = Alert(
                timestamp=time.time(),
                severity=self._classify_severity(combined_score),
                source="streaming_pipeline",
                message=f"Anomaly detected: score={combined_score:.4f}",
                score=combined_score,
                metadata={"event_index": self.metrics.total_processed}
            )
            self.alert_buffer.append(alert)
            self.metrics.total_alerts += 1
            self.metrics.total_anomalies += 1

        # 6. Postprocessing
        result = {
            "score": combined_score,
            "is_anomaly": combined_score > self.alert_threshold,
            "alert": alert,
        }
        for pp in self.postprocessors:
            result = pp(result)

        # 7. Update metrics
        elapsed = (time.time() - start) * 1000
        self.metrics.total_processed += 1
        self.metrics.processing_time_ms = elapsed
        self.metrics.throughput_per_sec = 1000.0 / max(elapsed, 0.001)

        if len(self.score_window) > 0:
            self.metrics.window_anomaly_rate = (
                sum(1 for s in self.score_window if s > self.alert_threshold) /
                len(self.score_window)
            )

        return result

    def process_batch(self, batch: np.ndarray) -> Dict:
        """Process a batch of events."""
        results = []
        for event in batch:
            results.append(self.process_event(event))

        scores = [r["score"] for r in results]
        anomalies = [r for r in results if r["is_anomaly"]]

        return {
            "n_processed": len(batch),
            "n_anomalies": len(anomalies),
            "anomaly_rate": len(anomalies) / len(batch) if len(batch) > 0 else 0,
            "mean_score": float(np.mean(scores)),
            "max_score": float(np.max(scores)),
            "metrics": self._get_metrics_dict(),
        }

    def _update_stats(self, x):
        self.n_seen += 1
        delta = x - self.running_mean
        self.running_mean += delta / self.n_seen
        delta2 = x - self.running_mean
        self.running_var += (delta * delta2 - self.running_var) / self.n_seen

    def _classify_severity(self, score):
        if score > 0.95:
            return AlertSeverity.CRITICAL
        elif score > 0.85:
            return AlertSeverity.WARNING
        return AlertSeverity.INFO

    def _get_metrics_dict(self):
        return {
            "total_processed": self.metrics.total_processed,
            "total_anomalies": self.metrics.total_anomalies,
            "total_alerts": self.metrics.total_alerts,
            "processing_time_ms": self.metrics.processing_time_ms,
            "throughput_per_sec": self.metrics.throughput_per_sec,
            "window_anomaly_rate": self.metrics.window_anomaly_rate,
        }

    def get_summary(self) -> Dict:
        return {
            "metrics": self._get_metrics_dict(),
            "n_alerts": len(self.alert_buffer),
            "alert_severity_counts": {
                "critical": sum(1 for a in self.alert_buffer if a.severity == AlertSeverity.CRITICAL),
                "warning": sum(1 for a in self.alert_buffer if a.severity == AlertSeverity.WARNING),
                "info": sum(1 for a in self.alert_buffer if a.severity == AlertSeverity.INFO),
            },
        }


# ── Standard Detectors for Pipeline ──

def mahalanobis_detector(x, mean, var):
    """Mahalanobis distance-based anomaly score (0-1 normalized)."""
    safe_var = np.maximum(var, 1e-10)
    d = np.sqrt(np.sum((x - mean)**2 / safe_var))
    return float(1.0 - np.exp(-d / np.sqrt(len(x))))


def zscore_detector(x, mean, var):
    """Max z-score across features."""
    safe_std = np.sqrt(np.maximum(var, 1e-10))
    z = np.max(np.abs(x - mean) / safe_std)
    return float(min(z / 5.0, 1.0))


def ewma_detector(x, mean, var, lambda_=0.3):
    """EWMA-based detector."""
    safe_std = np.sqrt(np.maximum(var, 1e-10))
    z = np.abs(x - mean) / safe_std
    ewma_score = np.max(z)  # simplified
    return float(min(ewma_score / 4.0, 1.0))


# ── CERN/LIGO Architecture Blueprint ──

ARCHITECTURE_SPEC = {
    "name": "Large-Scale Experiment Anomaly Detection Architecture",
    "tiers": {
        "L0_hardware_trigger": {
            "description": "Hardware-level trigger (FPGA/ASIC)",
            "latency": "<1μs",
            "data_rate": "40MHz → ~100kHz",
            "function": "Basic threshold and pattern matching",
            "reduction_factor": "400x"
        },
        "L1_online_filter": {
            "description": "Software trigger with simple ML models",
            "latency": "<10ms",
            "data_rate": "100kHz → ~1kHz",
            "function": "Online anomaly scoring (z-score, EWMA)",
            "reduction_factor": "100x",
            "implementation": "StreamingAnomalyPipeline with zscore_detector"
        },
        "L2_nearline_analysis": {
            "description": "Near-real-time batch analysis",
            "latency": "<1s",
            "data_rate": "1kHz → ~10Hz",
            "function": "Isolation Forest, changepoint detection",
            "reduction_factor": "100x",
            "implementation": "IsolationForestDetector + PELTDetector"
        },
        "L3_offline_deep": {
            "description": "Full offline analysis with deep models",
            "latency": "minutes-hours",
            "data_rate": "Stored events",
            "function": "Deep SVDD, SHAP explanations, physics constraints",
            "implementation": "DeepSVDDDetector + PhysicsConstrainedScorer + ExplainableAnomalyDetector"
        },
    },
    "data_flow": [
        "Raw Detector → L0 (FPGA) → L1 (Online) → L2 (Nearline) → L3 (Offline)",
        "Each tier reduces data volume by ~100x",
        "Total reduction: O(10^6) from raw to stored",
    ],
    "fault_tolerance": {
        "redundancy": "N+1 detector nodes per tier",
        "checkpointing": "Model state saved every 60s",
        "dead_letter_queue": "Failed events routed to separate analysis",
        "monitoring": "Prometheus metrics + Grafana dashboards",
    },
}
