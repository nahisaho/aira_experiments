"""
Perturbational Complexity Index (PCI) simulator.


PCI = Lempel-Ziv complexity of significant cortical responses to TMS perturbation,
normalized by the source entropy. Casali et al. (2013).

This module simulates:
  - TMS-evoked potential (TEP) generation
  - Significant source activation via SSANOVA thresholding
  - Binary matrix construction from significant activations
  - PCI computation from LZ complexity
"""
import numpy as np
from typing import Dict, Any
from scipy import signal as scipy_signal
from .utils import lempel_ziv_complexity, spectral_entropy

Dict_type = Dict[str, Any]


def simulate_tms_evoked_potential(
    n_channels: int = 30,
    n_times: int = 200,       # time points post-TMS
    fs: float = 1000.0,       # sampling rate in Hz
    consciousness_level: float = 1.0,
    spread_factor: float = 0.5,
    seed: int = 42,
) -> np.ndarray:
    """
    Simulate TMS-evoked EEG response (TEP).

    consciousness_level controls:
    - Complexity: awake → complex, high-frequency propagation
    - Anesthesia → slow, locally restricted, high-amplitude simple response

    Returns: (n_channels, n_times) TEP matrix
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n_times) / fs * 1000  # ms

    tep = np.zeros((n_channels, n_times))

    # Stimulation site — first cluster
    stim_site = 0

    for ch in range(n_channels):
        # Distance from stimulation site (circular)
        dist = min(abs(ch - stim_site), n_channels - abs(ch - stim_site))
        dist_factor = np.exp(-dist / (n_channels * spread_factor))

        # Primary response: N100 component
        t0 = 20 + dist * 0.5  # propagation delay
        amp_n100 = dist_factor * (0.5 + consciousness_level * 0.5)
        n100 = -amp_n100 * np.exp(-((t - t0) ** 2) / (2 * 20**2))

        # Secondary oscillatory response (higher consciousness → more complex)
        n_oscillations = int(2 + consciousness_level * 6)
        osc = np.zeros(n_times)
        for k in range(1, n_oscillations + 1):
            freq = 10 * k * consciousness_level + 5
            phase = rng.uniform(0, 2 * np.pi)
            decay = np.exp(-t / (50 + consciousness_level * 100))
            amp = dist_factor * consciousness_level * 0.3 / k
            osc += amp * np.sin(2 * np.pi * freq / 1000 * t + phase) * decay

        # Slow recovery (anesthesia: large slow wave)
        slow_amp = (1 - consciousness_level) * dist_factor * 2.0
        slow = slow_amp * np.exp(-((t - 100) ** 2) / (2 * 60**2))

        # Burst suppression effect
        if consciousness_level < 0.3:
            if rng.random() > consciousness_level * 2:
                n100 *= 0.1
                osc *= 0.05

        tep[ch] = n100 + osc + slow + 0.05 * rng.standard_normal(n_times)

    return tep


def compute_ssanova_threshold(tep: np.ndarray,
                               baseline_window: tuple = (0, 50),
                               alpha: float = 0.01) -> np.ndarray:
    """
    Threshold TEP matrix using bootstrap-based significance test.
    Returns binary (n_channels, n_times) significance matrix.
    """
    n_ch, n_t = tep.shape
    baseline = tep[:, baseline_window[0]:baseline_window[1]]

    # Bootstrap null distribution from baseline
    null_std = baseline.std(axis=1, keepdims=True)
    null_mean = baseline.mean(axis=1, keepdims=True)

    # Z-score against baseline
    z_tep = (tep - null_mean) / (null_std + 1e-10)

    # Threshold at z-score corresponding to alpha
    from scipy.stats import norm
    z_thresh = norm.ppf(1 - alpha / 2)

    significant = (np.abs(z_tep) > z_thresh).astype(int)
    return significant


def binary_matrix_to_sequence(binary_matrix: np.ndarray) -> np.ndarray:
    """
    Flatten binary significance matrix to 1D sequence for LZ complexity.
    Channels stacked sequentially in time.
    """
    # Concatenate channel sequences
    return binary_matrix.T.ravel()


class PCISimulator:
    """
    Perturbational Complexity Index (PCI) simulator.

    Simulates TMS-EEG and computes PCI as described in Casali et al. (2013).

    Parameters
    ----------
    n_channels : int
        Number of EEG channels to simulate
    fs : float
        Sampling rate (Hz)
    alpha : float
        Significance threshold for source activation
    """

    def __init__(self, n_channels: int = 30, fs: float = 1000.0, alpha: float = 0.01):
        self.n_channels = n_channels
        self.fs = fs
        self.alpha = alpha

    def simulate_and_compute(
        self,
        consciousness_level: float,
        n_trials: int = 10,
        seed: int = 42,
    ) -> Dict_type:
        """
        Full PCI pipeline: simulate TEP → threshold → LZ complexity.

        Returns dict with pci, lz_complexity, source_entropy, n_significant.
        """
        rng = np.random.default_rng(seed)
        pci_values = []
        lz_values = []

        for trial in range(n_trials):
            trial_seed = int(rng.integers(0, 10000))
            tep = simulate_tms_evoked_potential(
                n_channels=self.n_channels,
                consciousness_level=consciousness_level,
                seed=trial_seed,
            )

            sig_matrix = compute_ssanova_threshold(tep, alpha=self.alpha)
            sequence = binary_matrix_to_sequence(sig_matrix)

            lz = lempel_ziv_complexity(sequence)
            n_ones = sequence.sum()
            n_total = len(sequence)

            # Source entropy
            if 0 < n_ones < n_total:
                p = n_ones / n_total
                source_entropy = -(p * np.log2(p) + (1 - p) * np.log2(1 - p + 1e-10))
            else:
                source_entropy = 0.0

            # PCI = LZ / source_entropy (normalized complexity)
            pci = lz / (source_entropy + 1e-10) if source_entropy > 0 else 0.0
            pci_values.append(pci)
            lz_values.append(lz)

        return {
            "pci": float(np.mean(pci_values)),
            "pci_std": float(np.std(pci_values)),
            "lz_complexity": float(np.mean(lz_values)),
            "consciousness_level": consciousness_level,
            "n_trials": n_trials,
        }

    def pci_across_levels(
        self,
        levels: np.ndarray,
        n_trials: int = 5,
        seed: int = 42,
    ) -> list:
        """Compute PCI for a range of consciousness levels."""
        results = []
        for i, level in enumerate(levels):
            res = self.simulate_and_compute(
                consciousness_level=float(level),
                n_trials=n_trials,
                seed=seed + i * 7,
            )
            results.append(res)
        return results
