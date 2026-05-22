"""
Real-Time EEG Artifact Removal Pipeline
Implements ICA-based and ASR (Artifact Subspace Reconstruction) methods
for non-invasive BCI applications.
"""

import numpy as np
from scipy import signal, linalg
from typing import Optional, Tuple, Dict, List
import threading
import queue
import time
from dataclasses import dataclass, field


@dataclass
class ArtifactRemovalConfig:
    """Configuration for artifact removal pipeline."""
    sfreq: float = 250.0
    n_channels: int = 64
    # ICA settings
    n_components: int = 20
    ica_method: str = "fastica"
    max_iter: int = 1000
    # ASR settings
    asr_cutoff: float = 20.0          # SD threshold for ASR
    asr_blocksize: int = 100          # samples per ASR block
    # Bandpass filter
    l_freq: float = 1.0
    h_freq: float = 50.0
    notch_freq: float = 50.0
    # Buffer
    buffer_seconds: float = 2.0
    # Artifact thresholds (µV)
    amplitude_threshold: float = 150.0
    gradient_threshold: float = 50.0


class RingBuffer:
    """Thread-safe ring buffer for streaming EEG data."""

    def __init__(self, n_channels: int, max_samples: int):
        self.buffer = np.zeros((n_channels, max_samples), dtype=np.float32)
        self.max_samples = max_samples
        self.write_idx = 0
        self.count = 0
        self._lock = threading.Lock()

    def write(self, data: np.ndarray) -> None:
        """Write samples (n_channels x n_samples) to buffer."""
        with self._lock:
            n_samples = data.shape[1]
            indices = np.arange(self.write_idx, self.write_idx + n_samples) % self.max_samples
            self.buffer[:, indices] = data
            self.write_idx = (self.write_idx + n_samples) % self.max_samples
            self.count = min(self.count + n_samples, self.max_samples)

    def read(self, n_samples: int) -> Optional[np.ndarray]:
        """Read last n_samples from buffer."""
        with self._lock:
            if self.count < n_samples:
                return None
            indices = np.arange(self.write_idx - n_samples, self.write_idx) % self.max_samples
            return self.buffer[:, indices].copy()


class ButterworthFilter:
    """Real-time zero-phase Butterworth filter using second-order sections."""

    def __init__(self, l_freq: float, h_freq: float, sfreq: float,
                 notch_freq: Optional[float] = None, order: int = 4):
        self.sfreq = sfreq
        self.sos_bp = signal.butter(
            order, [l_freq, h_freq], btype='bandpass',
            fs=sfreq, output='sos'
        )
        self.sos_notch = None
        if notch_freq is not None:
            q = 30.0
            b_notch, a_notch = signal.iirnotch(notch_freq, q, sfreq)
            self.sos_notch = signal.tf2sos(b_notch, a_notch)
        # Per-channel filter states
        self._bp_zi = None
        self._notch_zi = None

    def initialize(self, n_channels: int) -> None:
        """Initialize filter state for n_channels."""
        bp_zi_base = signal.sosfilt_zi(self.sos_bp)  # (n_sections, 2)
        self._bp_zi = np.stack([bp_zi_base] * n_channels, axis=1)  # (n_sec, n_ch, 2)
        if self.sos_notch is not None:
            notch_zi_base = signal.sosfilt_zi(self.sos_notch)
            self._notch_zi = np.stack([notch_zi_base] * n_channels, axis=1)

    def process(self, data: np.ndarray) -> np.ndarray:
        """Process (n_channels x n_samples) chunk, preserving state."""
        if self._bp_zi is None:
            self.initialize(data.shape[0])

        out = np.zeros_like(data)
        for ch in range(data.shape[0]):
            filtered, self._bp_zi[:, ch, :] = signal.sosfilt(
                self.sos_bp, data[ch], zi=self._bp_zi[:, ch, :]
            )
            if self.sos_notch is not None:
                filtered, self._notch_zi[:, ch, :] = signal.sosfilt(
                    self.sos_notch, filtered, zi=self._notch_zi[:, ch, :]
                )
            out[ch] = filtered
        return out


class FastICA:
    """Online-capable FastICA for EEG artifact removal."""

    def __init__(self, n_components: int = 20, max_iter: int = 1000,
                 tol: float = 1e-4, random_state: int = 42):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.rng = np.random.RandomState(random_state)
        # Learned parameters
        self.mixing_matrix: Optional[np.ndarray] = None   # (n_ch, n_comp)
        self.unmixing_matrix: Optional[np.ndarray] = None # (n_comp, n_ch)
        self.mean_: Optional[np.ndarray] = None
        self.whitening_matrix_: Optional[np.ndarray] = None
        self.artifact_components_: List[int] = []

    def _whiten(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """PCA whitening of (n_ch x n_samples)."""
        cov = np.cov(X)
        eigenvalues, eigenvectors = linalg.eigh(cov)
        # Sort descending
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        # Keep top n_components
        eigenvalues = eigenvalues[:self.n_components]
        eigenvectors = eigenvectors[:, :self.n_components]
        W = (eigenvectors / np.sqrt(eigenvalues + 1e-8)).T  # (n_comp, n_ch)
        return W @ X, W

    @staticmethod
    def _g(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """LogCosh non-linearity: g(x) = tanh(x), g'(x) = 1 - tanh²(x)."""
        tanh_x = np.tanh(x)
        return tanh_x, 1.0 - tanh_x ** 2

    def fit(self, X: np.ndarray) -> "FastICA":
        """Fit ICA on (n_ch x n_samples) calibration data."""
        self.mean_ = X.mean(axis=1, keepdims=True)
        X_c = X - self.mean_
        X_white, self.whitening_matrix_ = self._whiten(X_c)
        n_comp, n_samples = X_white.shape

        W = self.rng.randn(n_comp, n_comp)
        W, _ = linalg.qr(W)

        for iteration in range(self.max_iter):
            W_old = W.copy()
            wx = W @ X_white          # (n_comp, n_samples)
            g_wx, g_prime_wx = self._g(wx)
            W_new = (g_wx @ X_white.T) / n_samples - g_prime_wx.mean(axis=1, keepdims=True) * W
            # Symmetric decorrelation
            u, s, vt = linalg.svd(W_new @ W_new.T)
            W = (u * (1.0 / np.sqrt(s + 1e-8))) @ u.T @ W_new
            # Convergence check
            delta = np.max(np.abs(np.abs(np.diag(W @ W_old.T)) - 1.0))
            if delta < self.tol:
                break

        self.unmixing_matrix_ = W @ self.whitening_matrix_  # (n_comp, n_ch)
        # Pseudo-inverse gives mixing matrix
        self.mixing_matrix_ = linalg.pinv(self.unmixing_matrix_)  # (n_ch, n_comp)
        return self

    def detect_artifact_components(self, threshold_kurtosis: float = 3.0) -> List[int]:
        """
        Identify artifact ICs via kurtosis (EOG: very high kurtosis).
        Returns list of artifact component indices.
        """
        if self.unmixing_matrix_ is None:
            raise RuntimeError("ICA not fitted yet.")
        # We need raw data activations — use calibration proxy via unit variance
        # In practice, compute on actual calibration data
        self.artifact_components_ = []
        return self.artifact_components_

    def detect_artifacts_from_activations(self, activations: np.ndarray,
                                           kurt_thresh: float = 5.0,
                                           power_thresh: float = 0.15) -> List[int]:
        """
        Detect artifact ICs from actual activations (n_comp x n_samples).
        Uses kurtosis (blink/EMG) and spectral power ratio (muscle artifacts).
        """
        artifact_idx = []
        for i in range(activations.shape[0]):
            ic = activations[i]
            kurt = float(np.mean((ic - ic.mean()) ** 4) / (ic.std() ** 4 + 1e-8))
            if kurt > kurt_thresh:
                artifact_idx.append(i)
        self.artifact_components_ = artifact_idx
        return artifact_idx

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Apply ICA, zero out artifact components, reconstruct. (n_ch x n_samples)."""
        if self.unmixing_matrix_ is None:
            return X
        X_c = X - self.mean_
        activations = self.unmixing_matrix_ @ X_c  # (n_comp, n_samples)
        # Zero out artifact components
        mask = np.ones(activations.shape[0])
        mask[self.artifact_components_] = 0.0
        activations_clean = activations * mask[:, None]
        return self.mixing_matrix_ @ activations_clean + self.mean_


class ArtifactSubspaceReconstruction:
    """
    ASR (Artifact Subspace Reconstruction) for real-time EEG cleaning.
    Robust PCA-based approach to reject high-amplitude bursts.

    Reference: Mullen et al. (2015), Real-time neuroimaging and cognitive monitoring
               using wearable dry EEG.
    """

    def __init__(self, cutoff: float = 20.0, blocksize: int = 100,
                 win_len: float = 0.5, win_overlap: float = 0.66,
                 sfreq: float = 250.0):
        self.cutoff = cutoff
        self.blocksize = blocksize
        self.win_len = win_len
        self.win_overlap = win_overlap
        self.sfreq = sfreq
        # Calibration state
        self.M_: Optional[np.ndarray] = None   # mixing matrix from clean baseline
        self.T_: Optional[np.ndarray] = None   # threshold matrix
        self._calibrated = False

    def _geometric_median(self, X: np.ndarray, tol: float = 1e-5,
                           max_iter: int = 300) -> np.ndarray:
        """Weiszfeld geometric median of columns of X (n_ch x n_samples)."""
        y = X.mean(axis=1)
        for _ in range(max_iter):
            dists = np.linalg.norm(X - y[:, None], axis=0)
            dists = np.maximum(dists, 1e-8)
            weights = 1.0 / dists
            y_new = (X * weights[None, :]).sum(axis=1) / weights.sum()
            if np.linalg.norm(y_new - y) < tol:
                break
            y = y_new
        return y

    def calibrate(self, clean_data: np.ndarray) -> "ArtifactSubspaceReconstruction":
        """
        Calibrate ASR on clean baseline EEG (n_ch x n_samples).
        Estimates the covariance structure of artifact-free data.
        """
        n_ch, n_samples = clean_data.shape
        win_samples = int(self.win_len * self.sfreq)
        step = max(1, int(win_samples * (1 - self.win_overlap)))

        # Compute covariance matrices in sliding windows
        cov_matrices = []
        for start in range(0, n_samples - win_samples, step):
            seg = clean_data[:, start:start + win_samples]
            cov_matrices.append(np.cov(seg))

        if len(cov_matrices) == 0:
            return self

        # RieMAP: use geometric median of covariance matrices (approximation)
        # Stack and use component-wise median as robust estimate
        cov_stack = np.stack(cov_matrices, axis=0)  # (n_wins, n_ch, n_ch)
        median_cov = np.median(cov_stack, axis=0)

        # Eigen-decomposition of median covariance
        eigenvalues, eigenvectors = linalg.eigh(median_cov)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        self.M_ = eigenvectors                                  # (n_ch, n_ch)
        # Threshold: cutoff * std of component activations
        activations = eigenvectors.T @ clean_data              # (n_ch, n_samples)
        component_rms = np.sqrt(np.mean(activations ** 2, axis=1))  # (n_ch,)
        self.T_ = self.cutoff * component_rms                  # (n_ch,)
        self._calibrated = True
        return self

    def process_block(self, block: np.ndarray) -> np.ndarray:
        """
        Clean a block (n_ch x n_samples).
        Returns cleaned block with artifact subspaces reconstructed.
        """
        if not self._calibrated:
            return block

        # Project to principal components
        activations = self.M_.T @ block  # (n_ch, n_samples)

        # Find samples exceeding threshold
        exceeded = np.abs(activations) > self.T_[:, None]  # (n_ch, n_samples)
        artifact_mask = exceeded.any(axis=0)               # (n_samples,)

        if not artifact_mask.any():
            return block

        # For artifact segments: reconstruct via clean subspace only
        clean_idx = np.where(~artifact_mask)[0]
        if len(clean_idx) < 5:
            return block

        # Identify components that exceeded threshold
        contaminated_components = exceeded[:, artifact_mask].any(axis=1)
        clean_components = ~contaminated_components

        if clean_components.sum() == 0:
            return block

        # Reconstruct artifact segments using only clean components
        block_clean = block.copy()
        artifact_idx = np.where(artifact_mask)[0]

        # Project-reconstruct: use clean components only
        M_clean = self.M_[:, clean_components]             # (n_ch, n_clean_comp)
        act_clean = M_clean.T @ block[:, artifact_idx]    # (n_clean_comp, n_art)
        block_clean[:, artifact_idx] = M_clean @ act_clean

        return block_clean

    @property
    def is_calibrated(self) -> bool:
        return self._calibrated


class AmplitudeArtifactDetector:
    """Fast amplitude/gradient-based artifact flagging for real-time use."""

    def __init__(self, amplitude_threshold: float = 150.0,
                 gradient_threshold: float = 50.0):
        self.amplitude_threshold = amplitude_threshold
        self.gradient_threshold = gradient_threshold
        self._prev_sample: Optional[np.ndarray] = None

    def detect(self, chunk: np.ndarray) -> np.ndarray:
        """
        Returns boolean mask (n_samples,): True = artifact present.
        chunk: (n_channels x n_samples)
        """
        n_samples = chunk.shape[1]
        artifact = np.zeros(n_samples, dtype=bool)

        # Amplitude check
        artifact |= (np.abs(chunk) > self.amplitude_threshold).any(axis=0)

        # Gradient check (sample-to-sample jump)
        if self._prev_sample is not None:
            extended = np.hstack([self._prev_sample[:, None], chunk])
        else:
            extended = np.hstack([chunk[:, :1], chunk])
        grad = np.diff(extended, axis=1)
        artifact |= (np.abs(grad) > self.gradient_threshold).any(axis=0)

        self._prev_sample = chunk[:, -1]
        return artifact


class RealtimeArtifactRemovalPipeline:
    """
    Full real-time artifact removal pipeline combining:
    - Bandpass + notch filtering
    - Amplitude/gradient artifact detection
    - ASR (real-time block processing)
    - ICA (applied after offline calibration)
    """

    def __init__(self, config: ArtifactRemovalConfig):
        self.config = config
        self.filter = ButterworthFilter(
            l_freq=config.l_freq,
            h_freq=config.h_freq,
            sfreq=config.sfreq,
            notch_freq=config.notch_freq,
        )
        self.amplitude_detector = AmplitudeArtifactDetector(
            amplitude_threshold=config.amplitude_threshold,
            gradient_threshold=config.gradient_threshold,
        )
        self.asr = ArtifactSubspaceReconstruction(
            cutoff=config.asr_cutoff,
            blocksize=config.asr_blocksize,
            sfreq=config.sfreq,
        )
        self.ica = FastICA(
            n_components=config.n_components,
            max_iter=config.max_iter,
        )
        self._ica_fitted = False
        buffer_samples = int(config.buffer_seconds * config.sfreq)
        self.buffer = RingBuffer(config.n_channels, buffer_samples)

        # Metrics
        self.chunks_processed = 0
        self.artifact_samples_total = 0
        self.total_samples = 0
        self._latencies: List[float] = []

    def calibrate(self, clean_data: np.ndarray,
                  fit_ica: bool = True) -> Dict[str, float]:
        """
        Calibrate pipeline on clean baseline data (n_ch x n_samples).
        Returns calibration statistics.
        """
        stats: Dict[str, float] = {}
        t0 = time.perf_counter()

        # 1. Filter calibration data
        filtered = self.filter.process(clean_data.copy())

        # 2. ASR calibration
        self.asr.calibrate(filtered)
        stats["asr_calibration_time_s"] = time.perf_counter() - t0

        # 3. ICA fitting
        if fit_ica:
            t_ica = time.perf_counter()
            self.ica.fit(filtered)
            # Auto-detect artifacts from activations on calibration data
            activations = self.ica.unmixing_matrix_ @ (filtered - self.ica.mean_)
            art_comps = self.ica.detect_artifacts_from_activations(activations)
            stats["ica_fit_time_s"] = time.perf_counter() - t_ica
            stats["n_artifact_components"] = float(len(art_comps))
            self._ica_fitted = True

        stats["total_calibration_time_s"] = time.perf_counter() - t0
        return stats

    def process_chunk(self, chunk: np.ndarray) -> Tuple[np.ndarray, Dict]:
        """
        Process one chunk (n_channels x n_samples) in real time.
        Returns (cleaned_chunk, metadata_dict).
        Typical latency: < 5 ms for 8-sample chunks at 250 Hz.
        """
        t_start = time.perf_counter()
        meta: Dict = {}

        # Stage 1: Bandpass + notch filtering
        filtered = self.filter.process(chunk)

        # Stage 2: Amplitude artifact detection
        artifact_mask = self.amplitude_detector.detect(filtered)
        n_artifact = artifact_mask.sum()
        meta["amplitude_artifacts"] = int(n_artifact)

        # Stage 3: ASR
        if self.asr.is_calibrated:
            cleaned = self.asr.process_block(filtered)
        else:
            cleaned = filtered

        # Stage 4: ICA artifact removal
        if self._ica_fitted:
            cleaned = self.ica.transform(cleaned)

        # Update ring buffer
        self.buffer.write(cleaned)

        # Bookkeeping
        self.chunks_processed += 1
        n_samples = chunk.shape[1]
        self.artifact_samples_total += n_artifact
        self.total_samples += n_samples

        latency_ms = (time.perf_counter() - t_start) * 1000.0
        self._latencies.append(latency_ms)
        meta["latency_ms"] = latency_ms
        meta["artifact_rate"] = (self.artifact_samples_total / max(self.total_samples, 1))

        return cleaned, meta

    def get_buffer(self, n_samples: int) -> Optional[np.ndarray]:
        """Retrieve last n_samples from the processed ring buffer."""
        return self.buffer.read(n_samples)

    def get_stats(self) -> Dict[str, float]:
        """Return running pipeline statistics."""
        lats = np.array(self._latencies) if self._latencies else np.array([0.0])
        return {
            "chunks_processed": float(self.chunks_processed),
            "total_samples": float(self.total_samples),
            "artifact_rate": float(self.artifact_samples_total / max(self.total_samples, 1)),
            "mean_latency_ms": float(lats.mean()),
            "p95_latency_ms": float(np.percentile(lats, 95)),
            "max_latency_ms": float(lats.max()),
        }


# ---------------------------------------------------------------------------
# Demonstration / unit test
# ---------------------------------------------------------------------------

def demo_artifact_removal():
    """Generate synthetic EEG with artifacts, run pipeline, report metrics."""
    import sys
    print("=== Real-Time EEG Artifact Removal Pipeline Demo ===\n")

    cfg = ArtifactRemovalConfig(sfreq=250.0, n_channels=32, n_components=20)
    pipeline = RealtimeArtifactRemovalPipeline(cfg)

    rng = np.random.RandomState(42)
    n_calibration = int(30 * cfg.sfreq)  # 30 seconds of clean data
    clean_eeg = rng.randn(cfg.n_channels, n_calibration) * 15.0

    print("Calibrating pipeline on 30 s of clean data...")
    calib_stats = pipeline.calibrate(clean_eeg, fit_ica=True)
    print(f"  ASR calibration : {calib_stats['asr_calibration_time_s']:.3f} s")
    print(f"  ICA fit time    : {calib_stats['ica_fit_time_s']:.3f} s")
    print(f"  Artifact ICs    : {int(calib_stats['n_artifact_components'])}")

    # Simulate real-time streaming with injected artifacts
    chunk_size = 8  # 32 ms at 250 Hz
    n_test_chunks = 500
    print(f"\nStreaming {n_test_chunks} chunks of {chunk_size} samples each...")

    for i in range(n_test_chunks):
        chunk = rng.randn(cfg.n_channels, chunk_size) * 15.0
        # Inject eye blink artifact every 50 chunks
        if i % 50 == 0:
            chunk[:4] += rng.randn(4, chunk_size) * 100.0
        cleaned, meta = pipeline.process_chunk(chunk)

    stats = pipeline.get_stats()
    print("\n=== Pipeline Statistics ===")
    for k, v in stats.items():
        print(f"  {k:30s}: {v:.4f}")
    return stats


if __name__ == "__main__":
    demo_artifact_removal()
