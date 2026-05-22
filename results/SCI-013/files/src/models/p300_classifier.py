"""
P300 Speller Adaptive Classifier with Transfer Learning.
Implements XDAWN spatial filtering, LDA/Riemannian classifier,
and domain adaptation for cross-subject transfer.
"""

import numpy as np
from scipy import linalg, signal
from typing import Optional, List, Tuple, Dict
import time


# ---------------------------------------------------------------------------
# XDAWN Spatial Filter (Rivet et al., 2009)
# ---------------------------------------------------------------------------

class XDAWNSpatialFilter:
    """
    XDAWN: Xsupervised signal Denoising And discriminant component ANalysis.
    Optimizes SNR of averaged ERPs (P300 paradigm).

    Reference: Rivet et al. (2009), xDAWN Algorithm to Enhance Evoked Potentials:
               Application to Brain–Computer Interface.
    """

    def __init__(self, n_components: int = 4, signal_cov_estimator: str = "epoch"):
        self.n_components = n_components
        self.signal_cov_estimator = signal_cov_estimator
        self.filters_: Optional[np.ndarray] = None    # (n_comp, n_ch)
        self.patterns_: Optional[np.ndarray] = None   # (n_ch, n_comp)
        self.evoked_: Optional[np.ndarray] = None     # (n_ch, n_times)

    def fit(self, X: np.ndarray, y: np.ndarray) -> "XDAWNSpatialFilter":
        """
        X: (n_trials x n_ch x n_times), y: binary (0/1), 1 = P300 target.
        """
        n_trials, n_ch, n_times = X.shape
        # Estimate evoked signal (average over target trials)
        target_mask = y == 1
        self.evoked_ = X[target_mask].mean(axis=0)  # (n_ch, n_times)

        # Signal covariance: Toeplitz structure of the averaged ERP
        A = self.evoked_ @ self.evoked_.T + 1e-8 * np.eye(n_ch)  # (n_ch, n_ch)

        # Noise covariance: total data covariance
        total_cov = np.zeros((n_ch, n_ch))
        for trial in X:
            total_cov += trial @ trial.T
        total_cov /= (n_trials * n_times)
        total_cov += 1e-6 * np.eye(n_ch)

        # Solve generalized eigenvalue: A @ w = lambda B @ w
        eigenvalues, eigenvectors = linalg.eigh(A, total_cov)
        # Sort descending
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        self.filters_ = eigenvectors[:, :self.n_components].T   # (n_comp, n_ch)
        self.patterns_ = linalg.pinv(self.filters_).T            # (n_ch, n_comp)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        X: (n_trials x n_ch x n_times) -> (n_trials x n_comp x n_times)
        or single trial (n_ch x n_times) -> (n_comp x n_times).
        """
        if self.filters_ is None:
            raise RuntimeError("XDAWN not fitted.")
        single = X.ndim == 2
        if single:
            X = X[np.newaxis]
        out = np.array([self.filters_ @ trial for trial in X])
        return out[0] if single else out


# ---------------------------------------------------------------------------
# Riemannian Geometry-based Classifier (SPD Manifold)
# ---------------------------------------------------------------------------

class RiemannianMDM:
    """
    Minimum Distance to Mean (MDM) classifier on SPD manifold.
    Uses Riemannian distance to class centroids.

    Reference: Barachant et al. (2012), Multiclass Brain-Computer Interface
               Classification by Riemannian Geometry.
    """

    def __init__(self, metric: str = "riemann", reg: float = 1e-6):
        self.metric = metric
        self.reg = reg
        self.class_means_: Optional[List[np.ndarray]] = None
        self.classes_: Optional[np.ndarray] = None

    @staticmethod
    def _covariance(X: np.ndarray, reg: float = 1e-6) -> np.ndarray:
        """Regularized covariance of (n_ch x n_times)."""
        C = X @ X.T / X.shape[1]
        C += reg * np.eye(C.shape[0])
        return C

    @staticmethod
    def _matrix_log(C: np.ndarray) -> np.ndarray:
        eigenvalues, eigenvectors = linalg.eigh(C)
        eigenvalues = np.maximum(eigenvalues, 1e-10)
        return eigenvectors @ np.diag(np.log(eigenvalues)) @ eigenvectors.T

    @staticmethod
    def _matrix_exp(S: np.ndarray) -> np.ndarray:
        eigenvalues, eigenvectors = linalg.eigh(S)
        return eigenvectors @ np.diag(np.exp(eigenvalues)) @ eigenvectors.T

    @staticmethod
    def _matrix_sqrt(C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        eigenvalues, eigenvectors = linalg.eigh(C)
        eigenvalues = np.maximum(eigenvalues, 1e-10)
        sqrt_vals = np.sqrt(eigenvalues)
        C_sqrt = eigenvectors @ np.diag(sqrt_vals) @ eigenvectors.T
        C_invsqrt = eigenvectors @ np.diag(1.0 / sqrt_vals) @ eigenvectors.T
        return C_sqrt, C_invsqrt

    def _riemannian_distance(self, C1: np.ndarray, C2: np.ndarray) -> float:
        _, C1_invsqrt = self._matrix_sqrt(C1)
        M = C1_invsqrt @ C2 @ C1_invsqrt
        eigenvalues = linalg.eigvalsh(M)
        eigenvalues = np.maximum(eigenvalues, 1e-10)
        return float(np.sqrt(np.sum(np.log(eigenvalues) ** 2)))

    def _frechet_mean(self, covs: List[np.ndarray],
                       max_iter: int = 50, tol: float = 1e-7) -> np.ndarray:
        """Frechet mean on Riemannian manifold (gradient descent)."""
        M = np.mean(covs, axis=0)
        for _ in range(max_iter):
            M_sqrt, M_invsqrt = self._matrix_sqrt(M)
            S = np.mean([
                M_sqrt @ self._matrix_log(M_invsqrt @ C @ M_invsqrt) @ M_sqrt
                for C in covs
            ], axis=0)
            M_new = M_sqrt @ self._matrix_exp(M_invsqrt @ S @ M_invsqrt) @ M_sqrt
            if np.linalg.norm(M_new - M, 'fro') < tol:
                break
            M = M_new
        return M

    def fit(self, X: np.ndarray, y: np.ndarray) -> "RiemannianMDM":
        """X: (n_trials x n_ch x n_times), y: class labels."""
        self.classes_ = np.unique(y)
        covs = [self._covariance(trial, self.reg) for trial in X]
        self.class_means_ = []
        for c in self.classes_:
            class_covs = [covs[i] for i in range(len(y)) if y[i] == c]
            mean_cov = self._frechet_mean(class_covs)
            self.class_means_.append(mean_cov)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self.class_means_ is None:
            raise RuntimeError("MDM not fitted.")
        single = X.ndim == 2
        if single:
            X = X[np.newaxis]
        predictions = []
        for trial in X:
            cov = self._covariance(trial, self.reg)
            dists = [self._riemannian_distance(cov, m) for m in self.class_means_]
            predictions.append(self.classes_[np.argmin(dists)])
        result = np.array(predictions)
        return result[0] if single else result

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        if self.class_means_ is None:
            raise RuntimeError("MDM not fitted.")
        single = X.ndim == 2
        if single:
            X = X[np.newaxis]
        probas = []
        for trial in X:
            cov = self._covariance(trial, self.reg)
            dists = np.array([self._riemannian_distance(cov, m) for m in self.class_means_])
            scores = -dists
            scores -= scores.max()
            proba = np.exp(scores)
            proba /= proba.sum()
            probas.append(proba)
        result = np.array(probas)
        return result[0] if single else result


# ---------------------------------------------------------------------------
# Euclidean Alignment (EA) for EEG Transfer Learning
# ---------------------------------------------------------------------------

class EuclideanAlignment:
    """
    Euclidean Alignment for rapid cross-subject EEG transfer.

    Reference: He & Wu (2020), Transfer Learning for Brain-Computer Interfaces:
               A Euclidean Space Data Alignment Approach.
    """

    def __init__(self, reg: float = 1e-6):
        self.reg = reg
        self.R_: Optional[np.ndarray] = None

    def compute_reference(self, X: np.ndarray) -> np.ndarray:
        """Compute per-subject reference matrix R = mean(C_i)^{-1/2}."""
        covs = np.array([t @ t.T / t.shape[1] for t in X])
        R_mean = covs.mean(axis=0) + self.reg * np.eye(covs.shape[1])
        eigenvalues, eigenvectors = linalg.eigh(R_mean)
        eigenvalues = np.maximum(eigenvalues, 1e-10)
        R = eigenvectors @ np.diag(1.0 / np.sqrt(eigenvalues)) @ eigenvectors.T
        return R

    def fit(self, X: np.ndarray) -> "EuclideanAlignment":
        self.R_ = self.compute_reference(X)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.R_ is None:
            raise RuntimeError("EA not fitted.")
        return np.array([self.R_ @ trial for trial in X])


# ---------------------------------------------------------------------------
# P300 Speller Adaptive Classifier
# ---------------------------------------------------------------------------

class P300SpellerConfig:
    n_channels: int = 32
    sfreq: float = 256.0
    epoch_tmin: float = 0.0
    epoch_tmax: float = 0.8
    n_repetitions: int = 6
    n_rows: int = 6
    n_cols: int = 6
    baseline_duration: float = 0.2
    n_xdawn_components: int = 4
    use_riemannian: bool = True
    use_ea: bool = True


class P300AdaptiveClassifier:
    """
    Adaptive P300 speller classifier combining:
    - XDAWN spatial filtering
    - Riemannian MDM or LDA classification
    - Euclidean Alignment for cross-session adaptation
    - Online Bayesian updating for drift correction
    """

    def __init__(self, config: P300SpellerConfig):
        self.config = config
        self.xdawn = XDAWNSpatialFilter(n_components=config.n_xdawn_components)
        self.ea = EuclideanAlignment() if config.use_ea else None
        self.clf = RiemannianMDM()
        self._is_fitted = False
        self._online_buffer_X: List[np.ndarray] = []
        self._online_buffer_y: List[int] = []
        self._online_update_count = 0

    def fit(self, X: np.ndarray, y: np.ndarray) -> "P300AdaptiveClassifier":
        if self.ea is not None:
            self.ea.fit(X)
            X_aligned = self.ea.transform(X)
        else:
            X_aligned = X
        self.xdawn.fit(X_aligned, y)
        X_xdawn = self.xdawn.transform(X_aligned)
        self.clf.fit(X_xdawn, y)
        self._is_fitted = True
        return self

    def transfer_fit(self, X_source: np.ndarray, y_source: np.ndarray,
                     X_target_unlabeled: Optional[np.ndarray] = None,
                     X_target_labeled: Optional[np.ndarray] = None,
                     y_target_labeled: Optional[np.ndarray] = None) -> "P300AdaptiveClassifier":
        """Transfer learning from source subject to new target subject."""
        if self.ea is not None:
            self.ea.fit(X_source)
            X_source_aligned = self.ea.transform(X_source)
        else:
            X_source_aligned = X_source

        self.xdawn.fit(X_source_aligned, y_source)
        X_source_xdawn = self.xdawn.transform(X_source_aligned)
        self.clf.fit(X_source_xdawn, y_source)
        self._is_fitted = True

        if X_target_labeled is not None and y_target_labeled is not None and len(X_target_labeled) >= 10:
            if self.ea is not None:
                R_target = self.ea.compute_reference(X_target_labeled)
                X_tgt_aligned = np.array([R_target @ t for t in X_target_labeled])
            else:
                X_tgt_aligned = X_target_labeled
            X_tgt_xdawn = self.xdawn.transform(X_tgt_aligned)
            n_blend = min(len(X_tgt_xdawn) * 3, len(X_source_xdawn))
            X_combined = np.concatenate([X_source_xdawn[:n_blend], X_tgt_xdawn], axis=0)
            y_combined = np.concatenate([y_source[:n_blend], y_target_labeled])
            self.clf.fit(X_combined, y_combined)

        return self

    def predict_character(self, row_epochs: np.ndarray,
                          col_epochs: np.ndarray,
                          n_rows: int = 6, n_cols: int = 6) -> Tuple[int, int, np.ndarray]:
        """Decode character from row/column flash responses."""
        if not self._is_fitted:
            raise RuntimeError("Classifier not fitted.")

        all_epochs = np.concatenate([row_epochs, col_epochs], axis=0)
        if self.ea is not None:
            R = self.ea.compute_reference(all_epochs)
            all_aligned = np.array([R @ t for t in all_epochs])
        else:
            all_aligned = all_epochs

        row_aligned = all_aligned[:n_rows]
        col_aligned = all_aligned[n_rows:]
        row_xdawn = self.xdawn.transform(row_aligned)
        col_xdawn = self.xdawn.transform(col_aligned)
        row_proba = self.clf.predict_proba(row_xdawn)
        col_proba = self.clf.predict_proba(col_xdawn)
        row_scores = row_proba[:, 1]
        col_scores = col_proba[:, 1]
        best_row = int(np.argmax(row_scores))
        best_col = int(np.argmax(col_scores))
        score_matrix = np.outer(row_scores, col_scores)
        return best_row, best_col, score_matrix

    def update_online(self, epoch: np.ndarray, label: int) -> None:
        self._online_buffer_X.append(epoch)
        self._online_buffer_y.append(label)
        self._online_update_count += 1
        if len(self._online_buffer_X) >= 20:
            X_new = np.array(self._online_buffer_X[-100:])
            y_new = np.array(self._online_buffer_y[-100:])
            if np.sum(y_new == 1) >= 5:
                X_new_xdawn = self.xdawn.transform(X_new)
                self.clf.fit(X_new_xdawn, y_new)

    def get_info_transfer_rate(self, n_symbols: int, n_correct: int,
                                trial_duration_s: float) -> Dict[str, float]:
        p = n_correct / max(n_symbols, 1)
        n = self.config.n_rows * self.config.n_cols
        if p <= 0 or p >= 1:
            B = 0.0
        else:
            B = np.log2(n) + p * np.log2(p) + (1 - p) * np.log2((1 - p) / (n - 1))
        B = max(B, 0.0)
        itr = B * 60.0 / max(trial_duration_s, 0.1)
        return {
            "bits_per_selection": B,
            "itr_bits_per_min": itr,
            "accuracy": p,
            "n_choices": float(n),
        }


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo_p300_classifier():
    print("=== P300 Speller Adaptive Classifier Demo ===\n")
    rng = np.random.RandomState(1)
    cfg = P300SpellerConfig()
    n_times = int((cfg.epoch_tmax - cfg.epoch_tmin) * cfg.sfreq)

    def make_p300_data(n_target: int, n_nontarget: int, snr_db: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
        t = np.linspace(cfg.epoch_tmin, cfg.epoch_tmax, n_times)
        p300_peak = np.exp(-0.5 * ((t - 0.3) / 0.05) ** 2)
        channel_weights = rng.rand(cfg.n_channels) * 2 - 0.5
        channel_weights[:8] += 1.5
        signal_amplitude = 10 ** (snr_db / 20)
        X, y = [], []
        for _ in range(n_target):
            noise = rng.randn(cfg.n_channels, n_times) * 10
            sig = np.outer(channel_weights, p300_peak) * signal_amplitude
            X.append(noise + sig)
            y.append(1)
        for _ in range(n_nontarget):
            noise = rng.randn(cfg.n_channels, n_times) * 10
            X.append(noise)
            y.append(0)
        return np.array(X), np.array(y)

    X_src, y_src = make_p300_data(200, 1000, snr_db=6.0)
    X_tgt_labeled, y_tgt_labeled = make_p300_data(30, 150, snr_db=4.0)
    X_tgt_test, y_tgt_test = make_p300_data(50, 250, snr_db=4.0)

    print(f"Source data  : {len(y_src)} trials ({y_src.sum()} targets)")
    print(f"Target train : {len(y_tgt_labeled)} trials ({y_tgt_labeled.sum()} targets)")
    print(f"Target test  : {len(y_tgt_test)} trials ({y_tgt_test.sum()} targets)\n")

    clf = P300AdaptiveClassifier(cfg)
    t0 = time.perf_counter()
    clf.fit(X_src, y_src)
    X_tgt_test_aligned_src = clf.ea.transform(X_tgt_test) if clf.ea else X_tgt_test
    src_only_predictions = clf.clf.predict(clf.xdawn.transform(X_tgt_test_aligned_src))
    src_acc = float(np.mean(src_only_predictions == y_tgt_test))

    clf2 = P300AdaptiveClassifier(cfg)
    clf2.transfer_fit(X_src, y_src, X_target_labeled=X_tgt_labeled, y_target_labeled=y_tgt_labeled)
    X_tgt_test_aligned = clf2.ea.transform(X_tgt_test) if clf2.ea else X_tgt_test
    tgt_predictions = clf2.clf.predict(clf2.xdawn.transform(X_tgt_test_aligned))
    tgt_acc = float(np.mean(tgt_predictions == y_tgt_test))
    fit_time = time.perf_counter() - t0

    print(f"Source-only accuracy (target domain) : {src_acc*100:.1f}%")
    print(f"Transfer learning accuracy           : {tgt_acc*100:.1f}%")
    print(f"Improvement                          : {(tgt_acc-src_acc)*100:+.1f}%")
    print(f"Total fit time                       : {fit_time*1000:.0f} ms")

    row_ep, _ = make_p300_data(6, 0, snr_db=8.0)
    col_ep, _ = make_p300_data(6, 0, snr_db=8.0)
    best_row, best_col, score_matrix = clf2.predict_character(row_ep, col_ep)
    print(f"\nCharacter prediction: row={best_row}, col={best_col}")

    itr_info = clf2.get_info_transfer_rate(n_symbols=36, n_correct=30, trial_duration_s=12.0)
    print("\nInformation Transfer Rate:")
    for k, v in itr_info.items():
        print(f"  {k:30s}: {v:.3f}")

    return {
        "source_only_accuracy": src_acc,
        "transfer_accuracy": tgt_acc,
        "improvement": tgt_acc - src_acc,
        "itr_bits_per_min": itr_info["itr_bits_per_min"],
    }


if __name__ == "__main__":
    demo_p300_classifier()
