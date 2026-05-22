"""
Utility functions: entropy, mutual information, and synthetic data generation.
"""
import numpy as np
from itertools import combinations
from scipy.stats import entropy as scipy_entropy
from scipy.signal import welch


def entropy(p: np.ndarray, base: float = 2) -> float:
    """Shannon entropy H(X) in bits."""
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p) / np.log(base)))


def joint_entropy(x: np.ndarray, y: np.ndarray, bins: int = 32) -> float:
    """Estimate joint entropy H(X,Y) via histogram."""
    hist, _, _ = np.histogram2d(x, y, bins=bins)
    hist = hist / hist.sum()
    return entropy(hist.ravel())


def mutual_information(x: np.ndarray, y: np.ndarray, bins: int = 32) -> float:
    """Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y)."""
    hx = entropy(np.histogram(x, bins=bins)[0])
    hy = entropy(np.histogram(y, bins=bins)[0])
    hxy = joint_entropy(x, y, bins=bins)
    return max(0.0, hx + hy - hxy)


def conditional_entropy(x: np.ndarray, y: np.ndarray, bins: int = 32) -> float:
    """H(X|Y) = H(X,Y) - H(Y)."""
    return joint_entropy(x, y, bins) - entropy(np.histogram(y, bins=bins)[0])


def transfer_entropy(x: np.ndarray, y: np.ndarray, lag: int = 1, bins: int = 16) -> float:
    """
    Transfer entropy TE(X→Y): information flow from X to Y.
    TE(X→Y) = H(Y_t|Y_{t-1}) - H(Y_t|Y_{t-1}, X_{t-1})
    """
    y_t = y[lag:]
    y_lag = y[:-lag]
    x_lag = x[:-lag]

    # H(Y_t | Y_{t-1})
    h_y_given_ylag = conditional_entropy(y_t, y_lag, bins)

    # H(Y_t | Y_{t-1}, X_{t-1}) via 3-var joint
    n = len(y_t)
    # Discretize
    def disc(v):
        v_min, v_max = v.min(), v.max()
        if v_max == v_min:
            return np.zeros(len(v), dtype=int)
        return ((v - v_min) / (v_max - v_min) * (bins - 1)).astype(int)

    yt_d = disc(y_t)
    yl_d = disc(y_lag)
    xl_d = disc(x_lag)

    # Joint index
    joint_yx = yt_d * bins + yl_d
    joint_yxy = yt_d * bins**2 + yl_d * bins + xl_d

    # H(Y_t, Y_lag)
    h_yt_yl = entropy(np.bincount(joint_yx, minlength=bins**2))
    # H(Y_lag)
    h_yl = entropy(np.bincount(yl_d, minlength=bins))
    # H(Y_t, Y_lag, X_lag)
    h_yt_yl_xl = entropy(np.bincount(joint_yxy, minlength=bins**3))
    # H(Y_lag, X_lag)
    h_yl_xl = entropy(np.bincount(yl_d * bins + xl_d, minlength=bins**2))

    te = h_yt_yl - h_yl - h_yt_yl_xl + h_yl_xl
    return max(0.0, te)


def lempel_ziv_complexity(binary_sequence: np.ndarray) -> float:
    """
    Lempel-Ziv-76 complexity (normalized).
    Fast O(n log n) implementation via copy phrases.
    Used in PCI calculation.
    """
    s = list(binary_sequence.astype(int))
    n = len(s)
    if n == 0:
        return 0.0

    # LZ76 complexity counting
    c = 1
    l = 1
    i = 0
    k = 1
    k_max = 1
    stop = False

    while not stop:
        if s[i + k - 1] == s[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1
                stop = True
        else:
            if k > k_max:
                k_max = k
            i += 1
            if i == l:
                c += 1
                l += k_max
                if l + 1 > n:
                    stop = True
                else:
                    i = 0
                    k = 1
                    k_max = 1
            else:
                k = 1

    # Normalize: C(n) / (n / log2(n))
    if n > 1:
        norm = n / max(np.log2(n), 1.0)
        return float(c / norm)
    return 0.0


def spectral_entropy(signal: np.ndarray, fs: float = 256.0) -> float:
    """Normalized spectral entropy of a time-series signal."""
    freqs, psd = welch(signal, fs=fs, nperseg=min(256, len(signal)))
    psd = psd / psd.sum()
    return entropy(psd)


def generate_anesthesia_data(
    n_channels: int = 8,
    n_samples: int = 2048,
    fs: float = 256.0,
    consciousness_level: float = 1.0,
    noise_level: float = 0.1,
    seed: int = 42,
) -> np.ndarray:
    """
    Generate synthetic multi-channel EEG/LFP data parameterized by
    consciousness level (0 = deep anesthesia, 1 = fully awake).

    Model:
    - Awake: high-frequency, low-amplitude, spatially decorrelated
    - Anesthesia: slow-wave (delta), high-amplitude, spatially correlated bursts

    Returns shape (n_channels, n_samples)
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n_samples) / fs

    # Frequency band power weights shift with consciousness
    # Awake: gamma/beta dominant; Anesthesia: delta dominant
    delta_power = 1.0 - 0.8 * consciousness_level   # 0.2–1.0
    alpha_power = 0.3 + 0.3 * consciousness_level   # 0.3–0.6
    beta_power  = 0.1 + 0.6 * consciousness_level   # 0.1–0.7
    gamma_power = 0.05 + 0.7 * consciousness_level  # 0.05–0.75

    data = np.zeros((n_channels, n_samples))

    for ch in range(n_channels):
        # Delta (1–4 Hz)
        f_delta = rng.uniform(1, 4)
        sig = delta_power * rng.uniform(2, 4) * np.sin(2 * np.pi * f_delta * t + rng.uniform(0, 2 * np.pi))

        # Alpha (8–12 Hz)
        f_alpha = rng.uniform(8, 12)
        sig += alpha_power * rng.uniform(0.5, 1.5) * np.sin(2 * np.pi * f_alpha * t + rng.uniform(0, 2 * np.pi))

        # Beta (15–30 Hz)
        f_beta = rng.uniform(15, 30)
        sig += beta_power * rng.uniform(0.3, 1.0) * np.sin(2 * np.pi * f_beta * t + rng.uniform(0, 2 * np.pi))

        # Gamma (30–80 Hz)
        f_gamma = rng.uniform(30, 80)
        sig += gamma_power * rng.uniform(0.1, 0.5) * np.sin(2 * np.pi * f_gamma * t + rng.uniform(0, 2 * np.pi))

        # Burst suppression (low consciousness_level only)
        if consciousness_level < 0.3:
            suppression_mask = (rng.random(n_samples) < 0.6).astype(float)
            sig *= suppression_mask

        # Independent noise
        sig += noise_level * rng.standard_normal(n_samples)

        # Inter-channel correlation increases as consciousness decreases
        if ch > 0:
            correlation = (1 - consciousness_level) * 0.6
            sig = (1 - correlation) * sig + correlation * data[0]

        data[ch] = sig

    return data


def discretize_channels(data: np.ndarray, n_states: int = 2) -> np.ndarray:
    """Binarize or discretize multi-channel data for Φ calculation."""
    thresholds = np.median(data, axis=1, keepdims=True)
    return (data > thresholds).astype(int)
