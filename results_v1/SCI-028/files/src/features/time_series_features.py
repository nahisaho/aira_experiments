"""
Time-series feature engineering for tokamak disruption prediction.
Covers MHD diagnostics, magnetics, Thomson scattering, ECE, and soft X-ray signals.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from scipy import signal, stats
from scipy.fft import rfft, rfftfreq


# ─── Signal groups ───────────────────────────────────────────────────────────

MAGNETIC_SIGNALS = [
    "ip",          # Plasma current [MA]
    "bt",          # Toroidal field [T]
    "li",          # Internal inductance
    "q95",         # Safety factor at 95% flux surface
    "betan",       # Normalised beta
    "kappa",       # Elongation
    "delta_u",     # Upper triangularity
    "delta_l",     # Lower triangularity
    "zmag",        # Magnetic axis height [m]
    "aminor",      # Minor radius [m]
    "drsep",       # D-Rsep (primary → secondary X-point distance)
]

THERMAL_SIGNALS = [
    "te_core",     # Core electron temperature [keV]
    "te_ped",      # Pedestal electron temperature [keV]
    "ne_core",     # Core electron density [10^19 m-3]
    "ne_ped",      # Pedestal electron density [10^19 m-3]
    "wdia",        # Diamagnetic energy [MJ]
    "wmhd",        # MHD energy [MJ]
    "tauE",        # Energy confinement time [s]
    "h98",         # H98(y,2) confinement enhancement factor
]

MHD_SIGNALS = [
    "mirnov_rms",  # Mirnov coil RMS amplitude
    "mirnov_n1",   # n=1 mode amplitude
    "mirnov_n2",   # n=2 mode amplitude
    "mirnov_m2n1", # (m,n)=(2,1) NTM amplitude
    "mirnov_m3n2", # (m,n)=(3,2) NTM amplitude
    "locked_mode", # Locked mode amplitude [a.u.]
    "vloop",       # Loop voltage [V]
    "dip_dt",      # dIp/dt [MA/s]
]

AUXILIARY_SIGNALS = [
    "p_nbi",       # NBI heating power [MW]
    "p_icrh",      # ICRH power [MW]
    "p_ecrh",      # ECRH power [MW]
    "p_rad",       # Total radiated power [MW]
    "p_hfs_rad",   # HFS radiated power [MW]
    "gas_puff",    # Gas puff rate [Pa m^3 s^-1]
]

ALL_SIGNALS = MAGNETIC_SIGNALS + THERMAL_SIGNALS + MHD_SIGNALS + AUXILIARY_SIGNALS


# ─── Feature specification ────────────────────────────────────────────────────

@dataclass
class WindowConfig:
    """Sliding window parameters for feature extraction."""
    short_ms: float = 50.0    # Short window: Mirnov / fast MHD [ms]
    medium_ms: float = 200.0  # Medium window: thermal + beta [ms]
    long_ms: float = 500.0    # Long window: global trend [ms]
    horizon_ms: float = 30.0  # Prediction horizon [ms]
    sample_rate_hz: float = 10_000.0  # Acquisition rate [Hz]

    @property
    def short_n(self) -> int:
        return int(self.short_ms * 1e-3 * self.sample_rate_hz)

    @property
    def medium_n(self) -> int:
        return int(self.medium_ms * 1e-3 * self.sample_rate_hz)

    @property
    def long_n(self) -> int:
        return int(self.long_ms * 1e-3 * self.sample_rate_hz)


@dataclass
class DisruptionFeatureSet:
    """Complete feature vector for a single time step."""
    # Statistical moments
    stat_features: np.ndarray = field(default_factory=lambda: np.array([]))
    # Spectral features
    spectral_features: np.ndarray = field(default_factory=lambda: np.array([]))
    # Rate-of-change features
    derivative_features: np.ndarray = field(default_factory=lambda: np.array([]))
    # MHD stability proxies
    stability_features: np.ndarray = field(default_factory=lambda: np.array([]))
    # Cross-signal correlation features
    correlation_features: np.ndarray = field(default_factory=lambda: np.array([]))

    def to_vector(self) -> np.ndarray:
        return np.concatenate([
            self.stat_features,
            self.spectral_features,
            self.derivative_features,
            self.stability_features,
            self.correlation_features,
        ])


# ─── Statistical features ─────────────────────────────────────────────────────

def statistical_features(x: np.ndarray, prefix: str = "") -> Dict[str, float]:
    """
    Compute statistical moments over a window.
    Returns mean, std, skewness, kurtosis, min, max, trend slope.
    """
    feats: Dict[str, float] = {}
    if len(x) < 4:
        return feats

    feats[f"{prefix}_mean"] = float(np.mean(x))
    feats[f"{prefix}_std"] = float(np.std(x, ddof=1))
    feats[f"{prefix}_skew"] = float(stats.skew(x))
    feats[f"{prefix}_kurt"] = float(stats.kurtosis(x))
    feats[f"{prefix}_min"] = float(np.min(x))
    feats[f"{prefix}_max"] = float(np.max(x))
    feats[f"{prefix}_range"] = float(np.ptp(x))

    # Linear trend slope (normalised by mean to get relative drift)
    t = np.arange(len(x), dtype=float)
    slope, _, _, _, _ = stats.linregress(t, x)
    feats[f"{prefix}_slope"] = float(slope)

    return feats


# ─── Spectral features ────────────────────────────────────────────────────────

# Key MHD frequency bands [Hz]
MHD_BANDS = {
    "sawteeth": (1, 5),       # Sawtooth precursors
    "tearing":  (5, 30),      # Classical tearing modes
    "ntm":      (5, 20),      # Neoclassical tearing modes
    "alf_low":  (50, 200),    # Low-frequency Alfvén eigenmodes
    "alf_high": (200, 1000),  # High-frequency TAE / EAE
    "hf_mhd":  (1000, 4000), # High-frequency MHD
}


def spectral_features(
    x: np.ndarray,
    fs: float = 10_000.0,
    prefix: str = "",
) -> Dict[str, float]:
    """
    Compute power spectral density features in MHD-relevant bands.
    """
    feats: Dict[str, float] = {}
    if len(x) < 64:
        return feats

    n = len(x)
    freqs = rfftfreq(n, d=1.0 / fs)
    psd = (np.abs(rfft(x)) ** 2) / n

    for band_name, (f_lo, f_hi) in MHD_BANDS.items():
        mask = (freqs >= f_lo) & (freqs < f_hi)
        band_power = float(np.sum(psd[mask]))
        feats[f"{prefix}_psd_{band_name}"] = band_power

    # Dominant frequency and its power
    peak_idx = int(np.argmax(psd[1:]) + 1)  # skip DC
    feats[f"{prefix}_peak_freq"] = float(freqs[peak_idx])
    feats[f"{prefix}_peak_power"] = float(psd[peak_idx])

    # Spectral centroid
    total_power = float(np.sum(psd[1:]))
    if total_power > 0:
        feats[f"{prefix}_spectral_centroid"] = float(
            np.sum(freqs[1:] * psd[1:]) / total_power
        )
    else:
        feats[f"{prefix}_spectral_centroid"] = 0.0

    return feats


# ─── Derivative / rate-of-change features ─────────────────────────────────────

def derivative_features(
    x: np.ndarray,
    dt: float,
    prefix: str = "",
    order: int = 2,
) -> Dict[str, float]:
    """
    Compute first and second temporal derivatives using Savitzky-Golay filter.
    """
    feats: Dict[str, float] = {}
    if len(x) < 11:
        return feats

    window = min(11, len(x) if len(x) % 2 != 0 else len(x) - 1)
    window = window if window >= 5 else 5

    dx = signal.savgol_filter(x, window_length=window, polyorder=3, deriv=1, delta=dt)
    feats[f"{prefix}_d1_mean"] = float(np.mean(dx))
    feats[f"{prefix}_d1_max_abs"] = float(np.max(np.abs(dx)))
    feats[f"{prefix}_d1_std"] = float(np.std(dx, ddof=1))

    if order >= 2:
        d2x = signal.savgol_filter(x, window_length=window, polyorder=3, deriv=2, delta=dt)
        feats[f"{prefix}_d2_mean"] = float(np.mean(d2x))
        feats[f"{prefix}_d2_max_abs"] = float(np.max(np.abs(d2x)))

    return feats


# ─── MHD stability proxy features ─────────────────────────────────────────────

def mhd_stability_features(
    signals: Dict[str, np.ndarray],
    cfg: WindowConfig,
) -> Dict[str, float]:
    """
    Compute physics-motivated stability proxy features.

    Key proxies:
    - Normalised beta proximity to Troyon limit
    - q=2 surface proximity (sawtooth / tearing trigger)
    - Greenwald density fraction
    - Locked-mode growth rate
    - Confinement degradation rate
    - Poloidal beta proxy
    """
    feats: Dict[str, float] = {}

    # --- Troyon normalised beta limit proximity ---
    if "betan" in signals and "ip" in signals and "bt" in signals and "aminor" in signals:
        betan = signals["betan"]
        ip_ma = signals["ip"]        # [MA]
        bt_t  = signals["bt"]        # [T]
        a_m   = signals["aminor"]    # [m]
        # Troyon limit: βN_max ≈ 3.5–4.0 (device dependent)
        troyon_limit = 3.5
        troyon_fraction = betan / (troyon_limit + 1e-9)
        feats["betan_troyon_frac"] = float(np.mean(troyon_fraction[-cfg.medium_n:]))

    # --- Greenwald density fraction ---
    if "ne_core" in signals and "ip" in signals and "aminor" in signals:
        ne = signals["ne_core"]   # [10^19 m^-3]
        ip = signals["ip"]        # [MA]
        a  = signals["aminor"]    # [m]
        # n_G = I_p [MA] / (π a² [m²]) × 10^20 m^-3
        # -> fraction = n_e / n_G
        n_greenwald = ip / (np.pi * a**2 + 1e-9)  # [10^19 m^-3 × 10]
        greenwald_frac = ne / (n_greenwald * 10.0 + 1e-9)
        feats["greenwald_fraction"] = float(np.mean(greenwald_frac[-cfg.medium_n:]))

    # --- Locked mode rate of change ---
    if "locked_mode" in signals:
        lm = signals["locked_mode"]
        dt = 1.0 / cfg.sample_rate_hz
        if len(lm) > 11:
            dlm = np.gradient(lm, dt)
            feats["locked_mode_growth"] = float(np.mean(dlm[-cfg.short_n:]))
            feats["locked_mode_rms"] = float(np.sqrt(np.mean(lm[-cfg.medium_n:]**2)))

    # --- H-factor confinement degradation ---
    if "h98" in signals:
        h98 = signals["h98"]
        if len(h98) >= cfg.long_n:
            dh = np.gradient(h98, 1.0 / cfg.sample_rate_hz)
            feats["h98_mean"] = float(np.mean(h98[-cfg.medium_n:]))
            feats["h98_degradation"] = float(np.mean(dh[-cfg.long_n:]))

    # --- q95 proximity to q=2 ---
    if "q95" in signals:
        q95 = signals["q95"]
        # q=2 is the primary tearing surface; lower q95 → greater instability risk
        feats["q95_mean"]   = float(np.mean(q95[-cfg.medium_n:]))
        feats["q95_min"]    = float(np.min(q95[-cfg.long_n:]))
        feats["q95_margin"] = float(np.mean(q95[-cfg.medium_n:]) - 2.0)

    # --- Radiated power fraction ---
    if "p_rad" in signals and "p_nbi" in signals and "p_icrh" in signals and "p_ecrh" in signals:
        p_heat = (signals["p_nbi"] + signals["p_icrh"] + signals["p_ecrh"])
        p_rad  = signals["p_rad"]
        f_rad  = p_rad / (p_heat + 1e-3)
        feats["f_rad_mean"] = float(np.mean(f_rad[-cfg.medium_n:]))
        feats["f_rad_max"]  = float(np.max(f_rad[-cfg.long_n:]))

    return feats


# ─── Cross-signal correlation features ────────────────────────────────────────

def correlation_features(
    signals: Dict[str, np.ndarray],
    pairs: Optional[List[Tuple[str, str]]] = None,
    window_n: int = 500,
) -> Dict[str, float]:
    """
    Compute lagged cross-correlations between key signal pairs.
    Disruption precursors often manifest as changes in inter-signal coherence.
    """
    if pairs is None:
        pairs = [
            ("mirnov_rms", "locked_mode"),
            ("betan",      "mirnov_n1"),
            ("ne_ped",     "te_ped"),
            ("wmhd",       "p_rad"),
            ("dip_dt",     "vloop"),
            ("li",         "q95"),
        ]

    feats: Dict[str, float] = {}
    for sig_a, sig_b in pairs:
        if sig_a not in signals or sig_b not in signals:
            continue
        a = signals[sig_a][-window_n:]
        b = signals[sig_b][-window_n:]
        n = min(len(a), len(b))
        if n < 10:
            continue
        a, b = a[:n], b[:n]
        # Pearson correlation
        rho, pval = stats.pearsonr(a, b)
        key = f"corr_{sig_a}_{sig_b}"
        feats[f"{key}_rho"] = float(rho)
        feats[f"{key}_pval"] = float(pval)

        # Cross-correlation peak lag (in samples)
        xcorr = np.correlate(a - a.mean(), b - b.mean(), mode="full")
        lag_idx = int(np.argmax(np.abs(xcorr))) - (n - 1)
        feats[f"{key}_lag"] = float(lag_idx)

    return feats


# ─── Main feature extractor ───────────────────────────────────────────────────

class TokamakFeatureExtractor:
    """
    Extracts the complete feature vector from a multi-channel signal buffer.

    The extractor operates on a rolling buffer updated at each control cycle
    (nominally every 1 ms for a 1 kHz control loop, or every 0.1 ms at 10 kHz).
    """

    def __init__(self, cfg: Optional[WindowConfig] = None):
        self.cfg = cfg or WindowConfig()
        self._dt = 1.0 / self.cfg.sample_rate_hz

    def extract(
        self,
        signals: Dict[str, np.ndarray],
        return_dict: bool = False,
    ) -> np.ndarray | Dict[str, float]:
        """
        Extract all features from a signal buffer.

        Parameters
        ----------
        signals : dict[str, ndarray]
            Mapping from signal name to 1-D time series (length ≥ long_n).
        return_dict : bool
            If True, return feature dict (for debugging); else return numpy vector.
        """
        feats: Dict[str, float] = {}

        for sig_name, x in signals.items():
            if len(x) < 4:
                continue
            # Short window (fast MHD)
            x_short  = x[-self.cfg.short_n:]
            # Medium window (thermal / global)
            x_medium = x[-self.cfg.medium_n:]
            # Long window (trend)
            x_long   = x[-self.cfg.long_n:]

            feats.update(statistical_features(x_short,  prefix=f"{sig_name}_s"))
            feats.update(statistical_features(x_medium, prefix=f"{sig_name}_m"))
            feats.update(statistical_features(x_long,   prefix=f"{sig_name}_l"))

            feats.update(spectral_features(x_short,  fs=self.cfg.sample_rate_hz, prefix=f"{sig_name}_s"))
            feats.update(spectral_features(x_medium, fs=self.cfg.sample_rate_hz, prefix=f"{sig_name}_m"))

            feats.update(derivative_features(x_short,  dt=self._dt, prefix=f"{sig_name}_s"))
            feats.update(derivative_features(x_medium, dt=self._dt, prefix=f"{sig_name}_m"))

        feats.update(mhd_stability_features(signals, self.cfg))
        feats.update(correlation_features(signals, window_n=self.cfg.medium_n))

        if return_dict:
            return feats
        return np.array(list(feats.values()), dtype=np.float32)

    @property
    def feature_names(self) -> List[str]:
        """Return feature names (requires a dummy signal call)."""
        dummy = {s: np.zeros(self.cfg.long_n) for s in ALL_SIGNALS[:4]}
        d = self.extract(dummy, return_dict=True)
        return list(d.keys())
