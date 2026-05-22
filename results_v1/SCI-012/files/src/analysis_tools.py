"""
Analysis tools for SNN outputs:
  - Firing rate estimation (PSTH, instantaneous rate)
  - Phase synchronisation (PLV, PPC)
  - Information transfer (mutual information, transfer entropy)
  - Spectral analysis (LFP power spectrum)
"""

import numpy as np
from typing import Optional, Tuple, List, Dict
from scipy import signal, stats


# ---------------------------------------------------------------------------
# Firing rate estimation
# ---------------------------------------------------------------------------

def psth(spike_times: np.ndarray, T_ms: float,
         bin_ms: float = 10.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Peri-Stimulus Time Histogram (PSTH).

    Args:
        spike_times: 1D array of spike times [ms]
        T_ms:        total simulation duration [ms]
        bin_ms:      bin width [ms]

    Returns:
        (counts, bin_edges) where counts is in spikes/bin
    """
    bins = np.arange(0, T_ms + bin_ms, bin_ms)
    counts, edges = np.histogram(spike_times, bins=bins)
    return counts.astype(float), edges


def instantaneous_rate(spike_times: np.ndarray, T_ms: float,
                       sigma_ms: float = 10.0,
                       dt: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Gaussian-kernel smoothed instantaneous firing rate.

    Returns:
        (rate_Hz, time_ms) arrays
    """
    t = np.arange(0, T_ms, dt)
    rate = np.zeros(len(t))
    for sp in spike_times:
        rate += np.exp(-0.5 * ((t - sp) / sigma_ms)**2)
    rate /= (sigma_ms * np.sqrt(2 * np.pi)) * 1e-3  # convert to Hz
    return rate, t


def population_firing_rate(spike_records: List[Tuple[float, int]],
                            N: int, T_ms: float,
                            bin_ms: float = 10.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Population firing rate from spike records [(t_ms, neuron_id), ...].

    Returns:
        (rate_Hz, bin_centers_ms)
    """
    if not spike_records:
        bins = np.arange(0, T_ms, bin_ms)
        return np.zeros(len(bins)-1), bins[:-1] + bin_ms/2
    times = np.array([s[0] for s in spike_records])
    bins  = np.arange(0, T_ms + bin_ms, bin_ms)
    counts, edges = np.histogram(times, bins=bins)
    rate = counts / (N * bin_ms * 1e-3)  # Hz
    centers = 0.5 * (edges[:-1] + edges[1:])
    return rate, centers


# ---------------------------------------------------------------------------
# Phase synchronisation
# ---------------------------------------------------------------------------

def compute_plv(lfp1: np.ndarray, lfp2: np.ndarray,
                fs: float, f_low: float, f_high: float) -> float:
    """
    Phase Locking Value (PLV) between two LFP signals in a frequency band.

    Args:
        lfp1, lfp2: LFP signals (same length)
        fs:         sampling rate [Hz]
        f_low/high: bandpass range [Hz]

    Returns:
        PLV in [0, 1]
    """
    sos = signal.butter(4, [f_low, f_high], btype='bandpass',
                        fs=fs, output='sos')
    x1 = signal.sosfiltfilt(sos, lfp1)
    x2 = signal.sosfiltfilt(sos, lfp2)
    phi1 = np.angle(signal.hilbert(x1))
    phi2 = np.angle(signal.hilbert(x2))
    dphi = phi1 - phi2
    return float(np.abs(np.mean(np.exp(1j * dphi))))


def compute_ppc(spike_phases: np.ndarray) -> float:
    """
    Pairwise Phase Consistency (PPC, Vinck et al. 2010).

    Args:
        spike_phases: phase of each spike relative to oscillation [radians]

    Returns:
        PPC in [0, 1]
    """
    n = len(spike_phases)
    if n < 2:
        return 0.0
    z = np.exp(1j * spike_phases)
    plv = np.abs(z.mean())
    ppc = (n * plv**2 - 1) / (n - 1)
    return float(max(0.0, ppc))


def oscillation_power_spectrum(lfp: np.ndarray, fs: float,
                               f_max: float = 150.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Welch power spectral density of LFP.

    Returns:
        (freqs, psd) up to f_max [Hz]
    """
    nperseg = min(len(lfp), int(fs))
    freqs, psd = signal.welch(lfp, fs=fs, nperseg=nperseg)
    mask = freqs <= f_max
    return freqs[mask], psd[mask]


def band_power(lfp: np.ndarray, fs: float,
               f_low: float, f_high: float) -> float:
    """Integrate PSD within frequency band [f_low, f_high] Hz."""
    freqs, psd = oscillation_power_spectrum(lfp, fs)
    mask = (freqs >= f_low) & (freqs <= f_high)
    return float(np.trapezoid(psd[mask], freqs[mask]) if hasattr(np, 'trapezoid') else np.trapz(psd[mask], freqs[mask]))


# ---------------------------------------------------------------------------
# Information-theoretic measures
# ---------------------------------------------------------------------------

def mutual_information_binned(x: np.ndarray, y: np.ndarray,
                               n_bins: int = 10) -> float:
    """
    Mutual information I(X;Y) via histogram binning.

    Args:
        x, y:   1D arrays (same length)
        n_bins: bins per variable

    Returns:
        MI in bits
    """
    c_xy, _, _ = np.histogram2d(x, y, bins=n_bins)
    c_x = c_xy.sum(axis=1)
    c_y = c_xy.sum(axis=0)
    n   = c_xy.sum()
    if n == 0:
        return 0.0

    p_xy = c_xy / n
    p_x  = c_x  / n
    p_y  = c_y  / n

    mi = 0.0
    for i in range(n_bins):
        for j in range(n_bins):
            if p_xy[i, j] > 0 and p_x[i] > 0 and p_y[j] > 0:
                mi += p_xy[i, j] * np.log2(p_xy[i, j] / (p_x[i] * p_y[j]))
    return float(mi)


def transfer_entropy(source: np.ndarray, target: np.ndarray,
                     lag: int = 1, n_bins: int = 4) -> float:
    """
    Transfer entropy TE(source → target) via binning.

    TE = I(target_future ; source_past | target_past)

    Args:
        source, target: discrete or continuous 1D time series (same length)
        lag:            temporal lag [samples]
        n_bins:         histogram bins

    Returns:
        TE in bits
    """
    n = len(target)
    if n <= lag + 1:
        return 0.0

    # Bin
    def _bin(arr):
        mn, mx = arr.min(), arr.max()
        if mx == mn:
            return np.zeros(len(arr), dtype=int)
        edges = np.linspace(mn, mx, n_bins + 1)
        return np.digitize(arr, edges[:-1]) - 1

    s = _bin(source)
    t = _bin(target)

    t_fut  = t[lag:]
    t_past = t[:-lag]
    s_past = s[:-lag]
    n_t = n - lag

    # Joint / marginal histograms
    def _joint3(a, b, c):
        h = np.zeros((n_bins, n_bins, n_bins))
        for i in range(len(a)):
            h[a[i], b[i], c[i]] += 1
        return h / (h.sum() + 1e-12)

    def _joint2(a, b):
        h = np.zeros((n_bins, n_bins))
        for i in range(len(a)):
            h[a[i], b[i]] += 1
        return h / (h.sum() + 1e-12)

    p_tft_tp_sp = _joint3(t_fut, t_past, s_past)
    p_tp_sp     = _joint2(t_past, s_past)
    p_tft_tp    = _joint2(t_fut, t_past)
    p_tp        = p_tp_sp.sum(axis=1)

    te = 0.0
    for a in range(n_bins):
        for b in range(n_bins):
            for c in range(n_bins):
                p3 = p_tft_tp_sp[a, b, c]
                p2s = p_tp_sp[b, c]
                p2t = p_tft_tp[a, b]
                p1  = p_tp[b]
                if p3 > 0 and p2s > 0 and p2t > 0 and p1 > 0:
                    te += p3 * np.log2(p3 * p1 / (p2s * p2t))
    return float(max(0.0, te))


# ---------------------------------------------------------------------------
# Spike-train similarity
# ---------------------------------------------------------------------------

def van_rossum_distance(spikes_a: np.ndarray, spikes_b: np.ndarray,
                         tau: float = 10.0, T_ms: float = 1000.0,
                         dt: float = 1.0) -> float:
    """
    Van Rossum (2001) spike train distance.

    Args:
        spikes_a/b: spike time arrays [ms]
        tau:        decay time constant [ms]
        T_ms:       total duration [ms]
        dt:         evaluation step [ms]

    Returns:
        Distance (lower = more similar)
    """
    t = np.arange(0, T_ms, dt)
    f_a = np.zeros(len(t))
    f_b = np.zeros(len(t))
    for sp in spikes_a:
        idx = int(sp / dt)
        if idx < len(t):
            f_a[idx:] += np.exp(-(t[idx:] - sp) / tau)
    for sp in spikes_b:
        idx = int(sp / dt)
        if idx < len(t):
            f_b[idx:] += np.exp(-(t[idx:] - sp) / tau)
    _trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    return float(np.sqrt(_trapz((f_a - f_b)**2, t) / tau))


# ---------------------------------------------------------------------------
# LFP approximation from population spike trains
# ---------------------------------------------------------------------------

def spikes_to_lfp(spike_records: List[Tuple[float, int]],
                   T_ms: float, dt: float = 1.0,
                   tau_lfp: float = 5.0) -> np.ndarray:
    """
    Approximate LFP as low-pass filtered population spike density.

    Returns:
        LFP array, shape (T,)
    """
    T = int(T_ms / dt)
    pop_rate = np.zeros(T)
    for (t_ms, _) in spike_records:
        idx = min(int(t_ms / dt), T - 1)
        pop_rate[idx] += 1.0

    # Exponential low-pass filter
    alpha = dt / tau_lfp
    lfp = np.zeros(T)
    for i in range(1, T):
        lfp[i] = lfp[i-1] * (1 - alpha) + pop_rate[i] * alpha
    return lfp


# ---------------------------------------------------------------------------
# Summary statistics for a full simulation
# ---------------------------------------------------------------------------

def analyse_simulation(spike_records: Dict[str, List[Tuple[float, int]]],
                       N_scaled: Dict[str, int],
                       T_ms: float,
                       dt_lfp: float = 1.0) -> Dict:
    """
    Compute comprehensive analysis of a Potjans-Diesmann simulation run.

    Returns dict of per-population and cross-population metrics.
    """
    results = {}
    lfps    = {}

    for pop_name, records in spike_records.items():
        N = N_scaled[pop_name]
        # Firing rate time series
        rate, bins = population_firing_rate(records, N, T_ms)
        # LFP
        lfp = spikes_to_lfp(records, T_ms, dt=dt_lfp)
        lfps[pop_name] = lfp
        # Spectral
        fs = 1000.0 / dt_lfp   # Hz
        freqs, psd = oscillation_power_spectrum(lfp, fs)
        gamma_power = band_power(lfp, fs, 30, 80)
        beta_power  = band_power(lfp, fs, 15, 30)
        theta_power = band_power(lfp, fs, 4,  8)

        results[pop_name] = {
            "mean_rate_Hz":   float(rate.mean()),
            "peak_rate_Hz":   float(rate.max()),
            "n_spikes":       len(records),
            "gamma_power":    gamma_power,
            "beta_power":     beta_power,
            "theta_power":    theta_power,
            "rate_series":    rate,
            "rate_bins_ms":   bins,
        }

    # Cross-population phase synchrony (L23E ↔ L4E if both present)
    if "L23E" in lfps and "L4E" in lfps:
        fs = 1000.0 / dt_lfp
        plv_gamma = compute_plv(lfps["L23E"], lfps["L4E"], fs,
                                f_low=30.0, f_high=80.0)
        results["cross"] = {"L23E_L4E_PLV_gamma": plv_gamma}

    return results
