"""
CSP (Common Spatial Pattern) + Deep Learning for Motor Imagery Classification.
Implements online CSP, EEGNet, ShallowConvNet, and a hybrid CSP-DNN pipeline.
"""

import numpy as np
from scipy import linalg
from typing import Optional, List, Tuple, Dict
import time


# ---------------------------------------------------------------------------
# CSP Implementation
# ---------------------------------------------------------------------------

class CommonSpatialPattern:
    """
    Multi-class CSP using One-vs-Rest strategy.
    Supports online update via exponential moving average of covariances.
    """

    def __init__(self, n_components: int = 4, reg: float = 1e-4,
                 n_classes: int = 4, online_alpha: float = 0.02):
        self.n_components = n_components
        self.reg = reg
        self.n_classes = n_classes
        self.online_alpha = online_alpha
        # Per-class spatial filters (n_classes, n_components, n_ch)
        self.W_: Optional[np.ndarray] = None
        self.n_channels_: Optional[int] = None
        # Running covariances for online update (n_classes, n_ch, n_ch)
        self._running_cov: Optional[np.ndarray] = None
        self._class_counts: Optional[np.ndarray] = None

    @staticmethod
    def _covariance(X: np.ndarray) -> np.ndarray:
        """Normalized covariance of (n_ch x n_samples)."""
        C = X @ X.T
        return C / (np.trace(C) + 1e-8)

    def _fit_binary(self, cov_a: np.ndarray, cov_b: np.ndarray) -> np.ndarray:
        """
        Solve generalized eigenvalue problem: cov_a @ w = λ (cov_a+cov_b) @ w.
        Returns filters (n_ch x n_components) = best + worst components.
        """
        cov_total = cov_a + cov_b + self.reg * np.eye(cov_a.shape[0])
        eigenvalues, eigenvectors = linalg.eigh(cov_a, cov_total)
        # Sort descending by eigenvalue
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        # Take top and bottom n_components // 2 filters
        half = self.n_components // 2
        filters = np.hstack([
            eigenvectors[:, :half],
            eigenvectors[:, -half:]
        ])  # (n_ch, n_components)
        return filters.T  # (n_components, n_ch)

    def fit(self, X: np.ndarray, y: np.ndarray) -> "CommonSpatialPattern":
        """
        Fit CSP on epochs (n_trials x n_ch x n_times).
        y: integer class labels (0..n_classes-1).
        """
        n_trials, n_ch, _ = X.shape
        self.n_channels_ = n_ch
        classes = np.unique(y)
        self.n_classes = len(classes)

        # Compute per-class mean covariance
        class_covs = {}
        for c in classes:
            trials_c = X[y == c]
            covs = np.array([self._covariance(t) for t in trials_c])
            class_covs[c] = covs.mean(axis=0)

        # Initialize running covariances
        self._running_cov = np.array([class_covs[c] for c in sorted(classes)])
        self._class_counts = np.array([np.sum(y == c) for c in sorted(classes)], dtype=float)

        # One-vs-Rest CSP filters
        W_list = []
        for c in sorted(classes):
            # "Rest" = mean of all other classes
            other_classes = [cc for cc in sorted(classes) if cc != c]
            cov_others = np.mean([class_covs[cc] for cc in other_classes], axis=0)
            W_c = self._fit_binary(class_covs[c], cov_others)  # (n_comp, n_ch)
            W_list.append(W_c)

        self.W_ = np.array(W_list)  # (n_classes, n_components, n_ch)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Extract log-variance CSP features.
        X: (n_trials x n_ch x n_times) or (n_ch x n_times) single trial.
        Returns: (n_trials x n_features) where n_features = n_classes * n_components.
        """
        if self.W_ is None:
            raise RuntimeError("CSP not fitted.")
        single = X.ndim == 2
        if single:
            X = X[np.newaxis]

        features = []
        for trial in X:
            trial_feats = []
            for c in range(self.n_classes):
                projected = self.W_[c] @ trial  # (n_comp, n_times)
                log_var = np.log(np.var(projected, axis=1) + 1e-8)
                trial_feats.append(log_var)
            features.append(np.concatenate(trial_feats))

        result = np.array(features)
        return result[0] if single else result

    def update_online(self, trial: np.ndarray, label: int) -> None:
        """
        Exponential moving average update for online adaptation.
        trial: (n_ch x n_times) single trial.
        """
        if self._running_cov is None:
            return
        new_cov = self._covariance(trial)
        alpha = self.online_alpha
        self._running_cov[label] = (1 - alpha) * self._running_cov[label] + alpha * new_cov
        # Refit filters from updated covariances
        self._refit_from_running()

    def _refit_from_running(self) -> None:
        """Refit spatial filters from current running covariances."""
        W_list = []
        n_classes = len(self._running_cov)
        for c in range(n_classes):
            other = [i for i in range(n_classes) if i != c]
            cov_others = np.mean(self._running_cov[other], axis=0)
            W_c = self._fit_binary(self._running_cov[c], cov_others)
            W_list.append(W_c)
        self.W_ = np.array(W_list)


# ---------------------------------------------------------------------------
# Pure-NumPy Deep Learning Layers (for inference without PyTorch dependency)
# ---------------------------------------------------------------------------

class Conv2dLayer:
    """Lightweight 2D conv layer for inference (no backprop)."""

    def __init__(self, weight: np.ndarray, bias: Optional[np.ndarray] = None,
                 stride: Tuple[int, int] = (1, 1), padding: Tuple[int, int] = (0, 0)):
        self.weight = weight  # (out_ch, in_ch, kH, kW)
        self.bias = bias
        self.stride = stride
        self.padding = padding

    def forward(self, x: np.ndarray) -> np.ndarray:
        """x: (batch, in_ch, H, W) → (batch, out_ch, H', W')"""
        out_ch, in_ch, kH, kW = self.weight.shape
        B, C, H, W = x.shape
        pH, pW = self.padding
        sH, sW = self.stride
        if pH > 0 or pW > 0:
            x = np.pad(x, ((0, 0), (0, 0), (pH, pH), (pW, pW)))
        H_out = (x.shape[2] - kH) // sH + 1
        W_out = (x.shape[3] - kW) // sW + 1
        out = np.zeros((B, out_ch, H_out, W_out), dtype=np.float32)
        for oc in range(out_ch):
            for i in range(H_out):
                for j in range(W_out):
                    patch = x[:, :, i*sH:i*sH+kH, j*sW:j*sW+kW]
                    out[:, oc, i, j] = np.sum(patch * self.weight[oc][None], axis=(1, 2, 3))
            if self.bias is not None:
                out[:, oc] += self.bias[oc]
        return out


class BatchNorm2dLayer:
    """Batch normalization (inference mode)."""

    def __init__(self, weight: np.ndarray, bias: np.ndarray,
                 running_mean: np.ndarray, running_var: np.ndarray, eps: float = 1e-5):
        self.weight = weight
        self.bias = bias
        self.running_mean = running_mean
        self.running_var = running_var
        self.eps = eps

    def forward(self, x: np.ndarray) -> np.ndarray:
        mean = self.running_mean[:, None, None]
        var = self.running_var[:, None, None]
        x_norm = (x - mean[None]) / np.sqrt(var[None] + self.eps)
        return x_norm * self.weight[None, :, None, None] + self.bias[None, :, None, None]


def elu(x: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    return np.where(x > 0, x, alpha * (np.exp(np.clip(x, -88, 0)) - 1))


def log_softmax(x: np.ndarray) -> np.ndarray:
    x_max = x.max(axis=-1, keepdims=True)
    log_sum = np.log(np.sum(np.exp(x - x_max), axis=-1, keepdims=True))
    return x - x_max - log_sum


# ---------------------------------------------------------------------------
# EEGNet (NumPy inference / PyTorch-style architecture description)
# ---------------------------------------------------------------------------

class EEGNetConfig:
    """EEGNet hyperparameters (Lawhern et al., 2018)."""
    n_classes: int = 4
    n_channels: int = 64
    n_times: int = 256          # 1 s at 256 Hz
    F1: int = 8                 # temporal filters
    D: int = 2                  # depth multiplier
    F2: int = 16                # pointwise filters (F1 * D)
    kern_len: int = 64          # temporal kernel = sfreq/2
    dropout: float = 0.5


def describe_eegnet(cfg: EEGNetConfig) -> str:
    """Returns a string description of EEGNet architecture."""
    F2 = cfg.F1 * cfg.D
    desc = f"""
EEGNet Architecture (Lawhern et al., 2018)
==========================================
Input : (1, {cfg.n_channels}, {cfg.n_times})

Block 1 — Temporal Convolution
  Conv2D  : (F1={cfg.F1}, 1, kern_len={cfg.kern_len})  → (F1, C, T)   [same padding, no bias]
  BatchNorm2d(F1)
  DepthwiseConv2D(F1→{cfg.F1*cfg.D}, (C,1), groups=F1)  [depthwise separable]
  BatchNorm2d(F1*D={F2})
  ELU activation
  AvgPool2D(1,4) → ({F2}, 1, T/4)
  Dropout(p={cfg.dropout})

Block 2 — Separable Convolution
  SeparableConv2D({F2}, kern=16) → ({F2}, 1, T/4)   [same padding]
  BatchNorm2d(F2={F2})
  ELU activation
  AvgPool2D(1,8) → ({F2}, 1, T/32)
  Dropout(p={cfg.dropout})

Classification Head
  Flatten → ({F2 * (cfg.n_times // 32)},)
  Linear({F2 * (cfg.n_times // 32)} → {cfg.n_classes})
  LogSoftmax

Total parameters ≈ {_count_eegnet_params(cfg):,}
"""
    return desc.strip()


def _count_eegnet_params(cfg: EEGNetConfig) -> int:
    F2 = cfg.F1 * cfg.D
    p = 0
    p += cfg.F1 * cfg.kern_len          # temporal conv
    p += 2 * cfg.F1                     # BN1
    p += cfg.F1 * cfg.D * cfg.n_channels  # depthwise conv
    p += 2 * F2                         # BN2
    p += F2 * 16                        # separable conv
    p += 2 * F2                         # BN3
    flat = F2 * (cfg.n_times // 32)
    p += flat * cfg.n_classes + cfg.n_classes
    return p


# ---------------------------------------------------------------------------
# Shallow ConvNet (Schirrmeister et al., 2017)
# ---------------------------------------------------------------------------

def describe_shallowconvnet(n_ch: int = 64, n_times: int = 256,
                             n_classes: int = 4) -> str:
    return f"""
ShallowConvNet Architecture (Schirrmeister et al., 2017)
=========================================================
Input : (1, {n_ch}, {n_times})

Layer 1 — Temporal Convolution
  Conv2D(1→40, kernel=(1,25))  → (40, {n_ch}, {n_times-24})

Layer 2 — Spatial Convolution  
  Conv2D(40→40, kernel=({n_ch},1)) → (40, 1, {n_times-24})
  BatchNorm2d
  Square activation

Layer 3 — Average Pooling
  AvgPool2D(kernel=(1,75), stride=(1,15))
  Log activation

Classifier
  Dropout(0.5)
  Linear(40×? → {n_classes})
  LogSoftmax
""".strip()


# ---------------------------------------------------------------------------
# Hybrid CSP + DNN Classifier
# ---------------------------------------------------------------------------

class CSPNumpyClassifier:
    """
    Lightweight online CSP classifier backed by LDA or softmax.
    Suitable for embedded real-time deployment.
    """

    def __init__(self, csp: CommonSpatialPattern, n_classes: int = 4,
                 shrinkage: float = 0.1):
        self.csp = csp
        self.n_classes = n_classes
        self.shrinkage = shrinkage
        # LDA parameters (fitted later)
        self.means_: Optional[np.ndarray] = None      # (n_classes, n_features)
        self.pooled_cov_inv_: Optional[np.ndarray] = None
        self.priors_: Optional[np.ndarray] = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> "CSPNumpyClassifier":
        """
        X: (n_trials x n_ch x n_times), y: (n_trials,)
        """
        self.csp.fit(X, y)
        features = self.csp.transform(X)  # (n_trials, n_features)
        self._fit_lda(features, y)
        return self

    def _fit_lda(self, features: np.ndarray, y: np.ndarray) -> None:
        """Fit regularized LDA (Ledoit-Wolf shrinkage)."""
        classes = np.unique(y)
        n_features = features.shape[1]
        self.means_ = np.zeros((len(classes), n_features))
        self.priors_ = np.zeros(len(classes))

        within_cov = np.zeros((n_features, n_features))
        for i, c in enumerate(classes):
            Xc = features[y == c]
            self.means_[i] = Xc.mean(axis=0)
            self.priors_[i] = len(Xc) / len(y)
            diff = Xc - self.means_[i]
            within_cov += diff.T @ diff

        within_cov /= (len(y) - len(classes))
        # Ledoit-Wolf shrinkage
        trace_ratio = np.trace(within_cov) / n_features
        within_cov = (1 - self.shrinkage) * within_cov + self.shrinkage * trace_ratio * np.eye(n_features)
        self.pooled_cov_inv_ = np.linalg.inv(within_cov)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Returns (n_trials x n_classes) probability estimates."""
        single = X.ndim == 2
        if single:
            X = X[np.newaxis]
        features = self.csp.transform(X)  # (n_trials, n_feat)

        scores = np.zeros((len(features), self.n_classes))
        for i in range(self.n_classes):
            diff = features - self.means_[i]
            mahal = np.sum(diff @ self.pooled_cov_inv_ * diff, axis=1)
            scores[:, i] = -0.5 * mahal + np.log(self.priors_[i] + 1e-12)

        # Softmax
        scores -= scores.max(axis=1, keepdims=True)
        proba = np.exp(scores)
        proba /= proba.sum(axis=1, keepdims=True)
        return proba[0] if single else proba

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.argmax(self.predict_proba(X), axis=-1)

    def update_online(self, trial: np.ndarray, label: int,
                      lr: float = 0.1) -> None:
        """
        Online gradient update for LDA means.
        trial: (n_ch x n_times)
        """
        self.csp.update_online(trial, label)
        feature = self.csp.transform(trial)  # (n_features,)
        # Move class mean toward new observation
        if self.means_ is not None:
            self.means_[label] = (1 - lr) * self.means_[label] + lr * feature


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo_csp_classifier():
    """Synthetic 4-class motor imagery classification demo."""
    print("=== CSP + LDA Motor Imagery Classification Demo ===\n")
    rng = np.random.RandomState(0)
    n_classes = 4
    n_channels = 32
    n_times = 256   # 1 s at 256 Hz
    n_trials_per_class = 50

    # Simulate class-discriminative spatial patterns
    patterns = rng.randn(n_classes, n_channels)  # one pattern per class
    X_list, y_list = [], []
    for c in range(n_classes):
        for _ in range(n_trials_per_class):
            base = rng.randn(n_channels, n_times) * 5.0
            # Add class-specific spatial activity
            base += np.outer(patterns[c], np.sin(2 * np.pi * (10 + c * 5) *
                             np.arange(n_times) / 256)) * 20
            X_list.append(base)
            y_list.append(c)

    X = np.array(X_list)
    y = np.array(y_list)

    # Shuffle
    idx = rng.permutation(len(y))
    X, y = X[idx], y[idx]
    split = int(0.8 * len(y))
    X_train, X_test = X[:split], X[split:]
    y_train, y_test = y[:split], y[split:]

    # Fit
    csp = CommonSpatialPattern(n_components=4, n_classes=n_classes)
    clf = CSPNumpyClassifier(csp, n_classes=n_classes)
    t0 = time.perf_counter()
    clf.fit(X_train, y_train)
    fit_time = time.perf_counter() - t0

    # Evaluate
    y_pred = clf.predict(X_test)
    acc = float(np.mean(y_pred == y_test))

    # Per-class accuracy
    per_class = {}
    for c in range(n_classes):
        mask = y_test == c
        per_class[c] = float(np.mean(y_pred[mask] == y_test[mask]))

    print(f"Training samples : {len(y_train)}")
    print(f"Test samples     : {len(y_test)}")
    print(f"Fit time         : {fit_time*1000:.1f} ms")
    print(f"Overall accuracy : {acc*100:.1f}%")
    print("Per-class accuracy:")
    mi_labels = ["Left Hand", "Right Hand", "Feet", "Rest"]
    for c, label in enumerate(mi_labels):
        print(f"  Class {c} ({label:12s}): {per_class[c]*100:.1f}%")

    # Online update simulation
    print("\nSimulating 20 online updates...")
    for i in range(20):
        trial = X_test[i]
        true_label = int(y_test[i])
        clf.update_online(trial, true_label, lr=0.05)

    y_pred_adapted = clf.predict(X_test)
    acc_adapted = float(np.mean(y_pred_adapted == y_test))
    print(f"Accuracy after online adaptation: {acc_adapted*100:.1f}%")

    # Architecture descriptions
    print("\n" + describe_eegnet(EEGNetConfig()))
    print("\n" + describe_shallowconvnet())

    return {
        "accuracy": acc,
        "accuracy_after_online": acc_adapted,
        "fit_time_ms": fit_time * 1000,
        "per_class_accuracy": per_class,
    }


if __name__ == "__main__":
    results = demo_csp_classifier()
