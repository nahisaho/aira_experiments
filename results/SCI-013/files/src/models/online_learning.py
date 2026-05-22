"""
Online Learning and Concept Drift Detection/Adaptation for EEG-BCI.
Implements:
- ADWIN (Adaptive Windowing) for drift detection
- Page-Hinkley test for sequential drift detection
- DDM (Drift Detection Method) 
- Incremental LDA with forgetting factor
- Reservoir sampling for balanced online learning
"""

import numpy as np
from typing import Optional, List, Tuple, Dict, Deque
from collections import deque
import time


# ---------------------------------------------------------------------------
# Drift Detection: ADWIN
# ---------------------------------------------------------------------------

class ADWIN:
    """
    ADWIN (ADaptive WINdowing) drift detector.
    Maintains a variable-length window and detects distribution change
    by testing if any two sub-windows have significantly different means.

    Reference: Bifet & Gavalda (2007), Learning from Time-Changing Data with
               Adaptive Windowing.
    """

    def __init__(self, delta: float = 0.002):
        self.delta = delta
        self._window: Deque[float] = deque()
        self._total_sum: float = 0.0
        self._total_n: int = 0
        self.drift_detected: bool = False
        self.n_detections: int = 0
        self._last_detection_idx: int = 0
        self._sample_count: int = 0

    def add_element(self, value: float) -> bool:
        """
        Add a new performance value (e.g., accuracy, loss).
        Returns True if drift is detected.
        """
        self._window.append(value)
        self._total_sum += value
        self._total_n += 1
        self._sample_count += 1
        self.drift_detected = False

        if self._total_n < 4:
            return False

        # Test all possible splits
        drift_found = False
        n0 = 0
        sum0 = 0.0
        for i, v in enumerate(self._window):
            n0 += 1
            sum0 += v
            n1 = self._total_n - n0
            sum1 = self._total_sum - sum0
            if n1 <= 0:
                break

            mu0 = sum0 / n0
            mu1 = sum1 / n1
            mu_total = self._total_sum / self._total_n

            epsilon_cut = self._compute_epsilon_cut(n0, n1)
            if abs(mu0 - mu1) >= epsilon_cut:
                # Drift detected: drop older part of window
                drift_found = True
                # Remove elements up to split point
                for _ in range(n0):
                    removed = self._window.popleft()
                    self._total_sum -= removed
                self._total_n -= n0
                break

        if drift_found:
            self.drift_detected = True
            self.n_detections += 1
            self._last_detection_idx = self._sample_count

        return drift_found

    def _compute_epsilon_cut(self, n0: int, n1: int) -> float:
        """Hoeffding-based epsilon cut threshold."""
        n = n0 + n1
        m = 1.0 / (1.0 / n0 + 1.0 / n1)
        return float(np.sqrt(np.log(4 * n / self.delta) / (2 * m)))

    @property
    def mean(self) -> float:
        if self._total_n == 0:
            return 0.0
        return self._total_sum / self._total_n

    @property
    def window_size(self) -> int:
        return self._total_n


class PageHinkley:
    """
    Page-Hinkley sequential drift detection test.
    Detects a persistent shift in the mean of a monitored variable.

    Reference: Page (1954), Continuous Inspection Schemes.
    """

    def __init__(self, min_instances: int = 30, delta: float = 0.005,
                 threshold: float = 50.0, alpha: float = 0.9999):
        self.min_instances = min_instances
        self.delta = delta
        self.threshold = threshold
        self.alpha = alpha
        self._n: int = 0
        self._sum: float = 0.0
        self._ph_sum: float = 0.0
        self._min_ph: float = float("inf")
        self.drift_detected: bool = False
        self.n_detections: int = 0

    def add_element(self, value: float) -> bool:
        """Add new observation. Returns True if drift detected."""
        self._n += 1
        self._sum += value
        mean_t = self._sum / self._n
        # PH cumulative sum
        self._ph_sum = self.alpha * self._ph_sum + (value - mean_t - self.delta)
        self._min_ph = min(self._min_ph, self._ph_sum)

        self.drift_detected = False
        if self._n >= self.min_instances:
            ph_stat = self._ph_sum - self._min_ph
            if ph_stat > self.threshold:
                self.drift_detected = True
                self.n_detections += 1
                # Reset
                self._n = 0
                self._sum = 0.0
                self._ph_sum = 0.0
                self._min_ph = float("inf")

        return self.drift_detected


class DDM:
    """
    DDM (Drift Detection Method) based on binomial error monitoring.

    Reference: Gama et al. (2004), Learning with drift detection.
    """

    def __init__(self, warning_level: float = 2.0, drift_level: float = 3.0,
                 min_instances: int = 30):
        self.warning_level = warning_level
        self.drift_level = drift_level
        self.min_instances = min_instances
        self._n: int = 0
        self._error_sum: float = 0.0
        self._min_pi: float = float("inf")
        self._min_si: float = float("inf")
        self.drift_detected: bool = False
        self.warning_detected: bool = False
        self.n_detections: int = 0

    def add_element(self, is_error: bool) -> Tuple[bool, bool]:
        """
        Add new prediction result.
        Returns (drift_detected, warning_detected).
        """
        self._n += 1
        if is_error:
            self._error_sum += 1.0
        p_i = self._error_sum / self._n
        s_i = float(np.sqrt(p_i * (1 - p_i) / self._n))

        self.drift_detected = False
        self.warning_detected = False

        if self._n >= self.min_instances:
            if p_i + s_i < self._min_pi + self._min_si:
                self._min_pi = p_i
                self._min_si = s_i

            if p_i + s_i > self._min_pi + self.drift_level * self._min_si:
                self.drift_detected = True
                self.n_detections += 1
                self._reset()
            elif p_i + s_i > self._min_pi + self.warning_level * self._min_si:
                self.warning_detected = True

        return self.drift_detected, self.warning_detected

    def _reset(self) -> None:
        self._n = 0
        self._error_sum = 0.0
        self._min_pi = float("inf")
        self._min_si = float("inf")


# ---------------------------------------------------------------------------
# Online Incremental LDA
# ---------------------------------------------------------------------------

class IncrementalLDA:
    """
    Incremental LDA with exponential forgetting for concept drift.
    Maintains running class statistics that decay over time.
    """

    def __init__(self, n_classes: int, n_features: int,
                 forgetting_factor: float = 0.99,
                 shrinkage: float = 0.1):
        self.n_classes = n_classes
        self.n_features = n_features
        self.lam = forgetting_factor
        self.shrinkage = shrinkage
        # Running statistics
        self._class_sum = np.zeros((n_classes, n_features), dtype=np.float64)
        self._class_sq_sum = np.zeros((n_classes, n_features, n_features), dtype=np.float64)
        self._class_count = np.zeros(n_classes, dtype=np.float64)
        # Trained discriminant
        self.W_: Optional[np.ndarray] = None  # (n_features, n_classes-1) or (n_features,)
        self.b_: Optional[np.ndarray] = None
        self._is_fitted = False

    def partial_fit(self, x: np.ndarray, label: int) -> None:
        """
        Update with a single labeled sample using forgetting factor.
        x: (n_features,)
        """
        # Decay all statistics
        self._class_sum *= self.lam
        self._class_sq_sum *= self.lam
        self._class_count *= self.lam
        # Accumulate new sample
        self._class_sum[label] += x
        self._class_sq_sum[label] += np.outer(x, x)
        self._class_count[label] += 1.0
        # Re-compute discriminant if enough data
        if self._class_count.min() >= 3:
            self._refit()

    def partial_fit_batch(self, X: np.ndarray, y: np.ndarray) -> None:
        """Batch update: process each sample sequentially."""
        for xi, yi in zip(X, y):
            self.partial_fit(xi, int(yi))

    def _refit(self) -> None:
        """Recompute LDA from current running statistics."""
        total_count = self._class_count.sum()
        if total_count < self.n_classes * 3:
            return

        # Class means
        means = np.array([
            self._class_sum[c] / (self._class_count[c] + 1e-8)
            for c in range(self.n_classes)
        ])  # (n_classes, n_features)

        global_mean = (self._class_sum.sum(axis=0)) / (total_count + 1e-8)

        # Within-class scatter
        Sw = np.zeros((self.n_features, self.n_features))
        for c in range(self.n_classes):
            if self._class_count[c] < 1:
                continue
            E_xx = self._class_sq_sum[c] / (self._class_count[c] + 1e-8)
            mu_c = means[c]
            Sw += self._class_count[c] * (E_xx - np.outer(mu_c, mu_c))
        Sw /= (total_count + 1e-8)

        # Regularization
        trace_ratio = np.trace(Sw) / self.n_features
        Sw += self.shrinkage * trace_ratio * np.eye(self.n_features)

        # Store for prediction
        self._means = means
        self._global_mean = global_mean
        self._Sw_inv = np.linalg.pinv(Sw)
        self._class_priors = self._class_count / (total_count + 1e-8)
        self._is_fitted = True

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """X: (n_samples x n_features) -> (n_samples x n_classes) probabilities."""
        if not self._is_fitted:
            return np.ones((len(X), self.n_classes)) / self.n_classes

        single = X.ndim == 1
        if single:
            X = X[np.newaxis]

        scores = np.zeros((len(X), self.n_classes))
        for c in range(self.n_classes):
            diff = X - self._means[c]
            mahal = np.sum(diff @ self._Sw_inv * diff, axis=1)
            scores[:, c] = -0.5 * mahal + np.log(self._class_priors[c] + 1e-12)

        scores -= scores.max(axis=1, keepdims=True)
        proba = np.exp(scores)
        proba /= proba.sum(axis=1, keepdims=True)
        return proba[0] if single else proba

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.argmax(self.predict_proba(X), axis=-1)


# ---------------------------------------------------------------------------
# Reservoir Sampling for Class-Balanced Online Buffer
# ---------------------------------------------------------------------------

class BalancedReservoir:
    """
    Class-balanced reservoir sampler for online EEG learning.
    Maintains a fixed-size buffer with balanced class distribution.
    """

    def __init__(self, capacity_per_class: int, n_classes: int, n_features: int):
        self.capacity = capacity_per_class
        self.n_classes = n_classes
        self.n_features = n_features
        self._buffers = [
            np.zeros((capacity_per_class, n_features), dtype=np.float32)
            for _ in range(n_classes)
        ]
        self._counts = np.zeros(n_classes, dtype=int)
        self._rng = np.random.RandomState(99)

    def add(self, x: np.ndarray, label: int) -> None:
        """Add sample to reservoir with random replacement after filling."""
        buf = self._buffers[label]
        n = self._counts[label]
        if n < self.capacity:
            buf[n] = x
        else:
            # Reservoir sampling: replace random existing sample
            idx = self._rng.randint(0, n + 1)
            if idx < self.capacity:
                buf[idx] = x
        self._counts[label] += 1

    def sample(self, n_per_class: int) -> Tuple[np.ndarray, np.ndarray]:
        """Draw balanced batch from reservoir."""
        X_list, y_list = [], []
        for c in range(self.n_classes):
            available = min(self._counts[c], self.capacity)
            if available == 0:
                continue
            k = min(n_per_class, available)
            idx = self._rng.choice(available, size=k, replace=False)
            X_list.append(self._buffers[c][idx])
            y_list.extend([c] * k)
        if not X_list:
            return np.zeros((0, self.n_features)), np.zeros(0, dtype=int)
        return np.concatenate(X_list, axis=0), np.array(y_list)

    @property
    def total_samples(self) -> int:
        return int(self._counts.sum())


# ---------------------------------------------------------------------------
# Full Online Learning Pipeline
# ---------------------------------------------------------------------------

class OnlineLearningPipeline:
    """
    Integrates:
    - Multiple drift detectors (ADWIN + DDM)
    - Incremental LDA classifier
    - Balanced reservoir for model retraining triggers
    - Automatic model reset/retrain on confirmed drift
    """

    def __init__(self, n_classes: int, n_features: int,
                 reservoir_capacity: int = 200,
                 forgetting_factor: float = 0.995,
                 retrain_every_n: int = 50):
        self.n_classes = n_classes
        self.n_features = n_features
        self.retrain_every_n = retrain_every_n

        self.classifier = IncrementalLDA(
            n_classes=n_classes,
            n_features=n_features,
            forgetting_factor=forgetting_factor,
        )
        self.reservoir = BalancedReservoir(
            capacity_per_class=reservoir_capacity,
            n_classes=n_classes,
            n_features=n_features,
        )
        self.adwin = ADWIN(delta=0.002)
        self.ddm = DDM(warning_level=2.0, drift_level=3.0)
        self.page_hinkley = PageHinkley(threshold=50.0)

        # Metrics
        self._n_processed = 0
        self._correct = 0
        self._drift_events: List[Dict] = []
        self._accuracy_history: List[float] = []
        self._window_acc: Deque[float] = deque(maxlen=100)

    def process_sample(self, x: np.ndarray, true_label: int) -> Dict:
        """
        Process one labeled sample online.
        Returns status dictionary.
        """
        t_start = time.perf_counter()
        result: Dict = {"drift_adwin": False, "drift_ddm": False,
                        "drift_ph": False, "retrained": False}

        # Predict
        pred = int(np.asarray(self.classifier.predict(x)).flat[0])
        is_correct = (pred == true_label)
        is_error = not is_correct

        # Update metrics
        self._n_processed += 1
        self._correct += int(is_correct)
        self._window_acc.append(float(is_correct))
        running_acc = np.mean(list(self._window_acc))
        self._accuracy_history.append(running_acc)

        # Update drift detectors
        adwin_drift = self.adwin.add_element(float(is_correct))
        ddm_drift, ddm_warning = self.ddm.add_element(is_error)
        ph_drift = self.page_hinkley.add_element(float(is_error))

        result["drift_adwin"] = adwin_drift
        result["drift_ddm"] = ddm_drift
        result["drift_ph"] = ph_drift
        result["pred"] = pred
        result["true"] = true_label
        result["correct"] = is_correct
        result["running_accuracy"] = running_acc

        # On confirmed drift: log event
        if adwin_drift or ddm_drift:
            event = {
                "sample_idx": self._n_processed,
                "detector": "ADWIN" if adwin_drift else "DDM",
                "running_acc": running_acc,
            }
            self._drift_events.append(event)
            result["drift_event"] = event

        # Online update: incremental LDA
        self.classifier.partial_fit(x, true_label)
        self.reservoir.add(x, true_label)

        # Periodic batch retrain from reservoir
        if self._n_processed % self.retrain_every_n == 0:
            X_batch, y_batch = self.reservoir.sample(n_per_class=20)
            if len(X_batch) >= self.n_classes * 5:
                self.classifier.partial_fit_batch(X_batch, y_batch)
                result["retrained"] = True

        result["latency_ms"] = (time.perf_counter() - t_start) * 1000
        return result

    def get_summary(self) -> Dict:
        return {
            "total_samples": self._n_processed,
            "overall_accuracy": self._correct / max(self._n_processed, 1),
            "n_drift_events_adwin": self.adwin.n_detections,
            "n_drift_events_ddm": self.ddm.n_detections,
            "adwin_window_size": self.adwin.window_size,
            "reservoir_total": self.reservoir.total_samples,
            "drift_events": self._drift_events,
        }


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo_online_learning():
    """Simulate concept drift in motor imagery BCI and test detection."""
    print("=== Online Learning & Concept Drift Detection Demo ===\n")
    rng = np.random.RandomState(7)
    n_classes = 4
    n_features = 16  # CSP features
    n_samples = 2000

    pipeline = OnlineLearningPipeline(
        n_classes=n_classes,
        n_features=n_features,
        reservoir_capacity=200,
        retrain_every_n=50,
    )

    # Simulate: drift occurs at sample 800 and 1500
    drift_points = [800, 1500]
    class_means = rng.randn(n_classes, n_features)  # initial class distribution

    print(f"Simulating {n_samples} samples with drift at {drift_points}...")
    n_drifts_detected = 0

    for i in range(n_samples):
        # Shift class means at drift points
        if i == drift_points[0]:
            class_means += rng.randn(n_classes, n_features) * 2.0  # large shift
            print(f"  [INJECTED DRIFT at sample {i}]")
        elif i == drift_points[1]:
            class_means += rng.randn(n_classes, n_features) * 1.5
            print(f"  [INJECTED DRIFT at sample {i}]")

        true_label = rng.randint(0, n_classes)
        x = class_means[true_label] + rng.randn(n_features) * 0.8

        result = pipeline.process_sample(x, true_label)
        if result.get("drift_adwin") or result.get("drift_ddm"):
            n_drifts_detected += 1

    summary = pipeline.get_summary()
    print(f"\n=== Online Learning Summary ===")
    print(f"  Total samples         : {summary['total_samples']}")
    print(f"  Overall accuracy      : {summary['overall_accuracy']*100:.1f}%")
    print(f"  ADWIN drift events    : {summary['n_drift_events_adwin']}")
    print(f"  DDM drift events      : {summary['n_drift_events_ddm']}")
    print(f"  Drift events log:")
    for ev in summary["drift_events"][:5]:
        print(f"    Sample {ev['sample_idx']:4d} | {ev['detector']:6s} | acc={ev['running_acc']*100:.1f}%")

    # Window-level accuracy
    acc_history = np.array(pipeline._accuracy_history)
    print(f"\n  Final 100-sample rolling accuracy: {acc_history[-100:].mean()*100:.1f}%")
    print(f"  First 200-sample accuracy        : {acc_history[:200].mean()*100:.1f}%")

    return {
        "overall_accuracy": summary["overall_accuracy"],
        "adwin_detections": summary["n_drift_events_adwin"],
        "ddm_detections": summary["n_drift_events_ddm"],
        "final_rolling_acc": float(acc_history[-100:].mean()),
    }


if __name__ == "__main__":
    demo_online_learning()
