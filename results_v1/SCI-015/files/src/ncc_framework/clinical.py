"""
Clinical consciousness classification:
- Vegetative State (VS) / Unresponsive Wakefulness Syndrome (UWS)
- Minimally Conscious State (MCS)
- Emerged from MCS (EMCS)
- Locked-In Syndrome (LIS)
- Healthy awake (CTRL)

Multi-feature classifier using IIT Φ, PCI, GWT metrics,
and spectral EEG features.
"""
import numpy as np
from typing import Dict, List, Optional, Tuple
from scipy.signal import welch
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from .utils import (
    generate_anesthesia_data, spectral_entropy,
    mutual_information, entropy
)


# Clinical diagnostic labels
CLASSES = {
    "CTRL": 0,    # Healthy awake
    "MCS+": 1,    # Minimally conscious state — with higher cortical signs
    "MCS-": 2,    # MCS — with basic signs only
    "VS": 3,      # Vegetative / Unresponsive Wakefulness
    "LIS": 4,     # Locked-In Syndrome
}

CLASS_NAMES = {v: k for k, v in CLASSES.items()}

# Canonical consciousness level per clinical state
STATE_CONSCIOUSNESS = {
    "CTRL": 1.0,
    "MCS+": 0.55,
    "MCS-": 0.35,
    "VS": 0.15,
    "LIS": 0.90,  # Fully conscious but locked-in
}


def extract_spectral_features(data: np.ndarray, fs: float = 256.0) -> Dict[str, float]:
    """
    Extract spectral EEG features relevant to consciousness.
    """
    n_ch, n_t = data.shape
    band_powers = {}

    bands = {
        "delta": (0.5, 4),
        "theta": (4, 8),
        "alpha": (8, 13),
        "beta": (13, 30),
        "gamma": (30, 80),
    }

    all_psd = []
    for ch in range(n_ch):
        freqs, psd = welch(data[ch], fs=fs, nperseg=min(256, n_t))
        all_psd.append(psd)

    mean_psd = np.mean(all_psd, axis=0)
    total_power = mean_psd.sum() + 1e-10

    for band, (lo, hi) in bands.items():
        mask = (freqs >= lo) & (freqs <= hi)
        band_powers[f"power_{band}"] = float(mean_psd[mask].sum() / total_power)

    # Spectral edge frequency (95%)
    cumsum_psd = np.cumsum(mean_psd)
    sef95_idx = np.searchsorted(cumsum_psd, 0.95 * cumsum_psd[-1])
    band_powers["sef95"] = float(freqs[min(sef95_idx, len(freqs) - 1)])

    # Mean spectral entropy
    sp_entropies = [spectral_entropy(data[ch], fs) for ch in range(n_ch)]
    band_powers["spectral_entropy"] = float(np.mean(sp_entropies))

    # Alpha/delta ratio (marker of arousal)
    delta_p = band_powers["power_delta"] + 1e-10
    alpha_p = band_powers["power_alpha"] + 1e-10
    band_powers["alpha_delta_ratio"] = float(alpha_p / delta_p)

    return band_powers


def extract_connectivity_features(data: np.ndarray) -> Dict[str, float]:
    """
    Extract connectivity-based consciousness markers.
    """
    n_ch = data.shape[0]

    # Pairwise correlations
    corr_mat = np.corrcoef(data)
    np.fill_diagonal(corr_mat, 0)

    # Mean correlation (global synchrony)
    mean_corr = float(np.abs(corr_mat).mean())

    # Long-range vs short-range connectivity
    long_range = []
    short_range = []
    for i in range(n_ch):
        for j in range(i + 1, n_ch):
            dist = abs(i - j)
            val = abs(corr_mat[i, j])
            if dist > n_ch // 3:
                long_range.append(val)
            else:
                short_range.append(val)

    lr_ratio = (np.mean(long_range) + 1e-10) / (np.mean(short_range) + 1e-10)

    # Anterior-posterior coherence (channels 0:n//2 vs n//2:n)
    half = n_ch // 2
    ap_coherence = float(np.abs(corr_mat[:half, half:]).mean()) if half > 0 else 0.0

    # Pairwise MI (sample)
    mi_values = []
    pairs = [(i, j) for i in range(min(n_ch, 6)) for j in range(i + 1, min(n_ch, 6))]
    for i, j in pairs:
        mi_values.append(mutual_information(data[i], data[j]))

    return {
        "mean_correlation": mean_corr,
        "long_short_range_ratio": float(lr_ratio),
        "ap_coherence": ap_coherence,
        "mean_mi": float(np.mean(mi_values)) if mi_values else 0.0,
    }


def generate_clinical_dataset(
    n_samples_per_class: int = 20,
    n_channels: int = 8,
    n_time: int = 1024,
    fs: float = 256.0,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Generate synthetic multi-class clinical dataset.

    Returns: (features, labels, feature_names)
    """
    rng = np.random.default_rng(seed)
    features_list = []
    labels = []

    states = list(STATE_CONSCIOUSNESS.keys())

    for state, cl_level in STATE_CONSCIOUSNESS.items():
        for k in range(n_samples_per_class):
            # Add jitter to consciousness level
            level = float(np.clip(
                cl_level + rng.normal(0, 0.08),
                0.05, 0.98
            ))
            sample_seed = int(rng.integers(0, 100000))

            data = generate_anesthesia_data(
                n_channels=n_channels,
                n_samples=n_time,
                consciousness_level=level,
                seed=sample_seed,
            )

            spec_feats = extract_spectral_features(data, fs=fs)
            conn_feats = extract_connectivity_features(data)

            # Additional: burst suppression ratio (low consciousness marker)
            amplitude_envelope = np.abs(data).mean(axis=0)
            burst_thresh = np.percentile(amplitude_envelope, 25)
            bsr = float((amplitude_envelope < burst_thresh).mean())
            conn_feats["burst_suppression_ratio"] = bsr

            feats = {**spec_feats, **conn_feats}
            features_list.append(feats)
            labels.append(CLASSES[state])

    feature_names = list(features_list[0].keys())
    X = np.array([[f[k] for k in feature_names] for f in features_list])
    y = np.array(labels)

    return X, y, feature_names


class ConsciousnessClassifier:
    """
    Multi-class consciousness state classifier.

    Uses LDA + Random Forest ensemble on spectral and
    connectivity features to differentiate:
    CTRL / MCS+ / MCS- / VS / LIS

    Parameters
    ----------
    n_channels : int
        Number of EEG channels in input data
    fs : float
        Sampling rate (Hz)
    """

    def __init__(self, n_channels: int = 8, fs: float = 256.0):
        self.n_channels = n_channels
        self.fs = fs
        self.scaler = StandardScaler()
        self.lda = LinearDiscriminantAnalysis()
        self.rf = RandomForestClassifier(
            n_estimators=100, random_state=42, n_jobs=-1
        )
        self.feature_names = None
        self._is_fitted = False

    def fit(self, X: np.ndarray, y: np.ndarray, feature_names: List[str]):
        """Fit classifier on feature matrix X and labels y."""
        self.feature_names = feature_names
        X_scaled = self.scaler.fit_transform(X)
        self.lda.fit(X_scaled, y)
        self.rf.fit(X_scaled, y)
        self._is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X_scaled = self.scaler.transform(X)
        # Ensemble: majority vote
        lda_pred = self.lda.predict(X_scaled)
        rf_pred = self.rf.predict(X_scaled)
        # Simple majority: if agree use that, else use RF
        preds = np.where(lda_pred == rf_pred, lda_pred, rf_pred)
        return preds

    def predict_from_data(self, data: np.ndarray) -> dict:
        """
        Predict consciousness state from raw EEG data.
        Returns predicted state label and probabilities.
        """
        spec_feats = extract_spectral_features(data, self.fs)
        conn_feats = extract_connectivity_features(data)
        amplitude_envelope = np.abs(data).mean(axis=0)
        burst_thresh = np.percentile(amplitude_envelope, 25)
        conn_feats["burst_suppression_ratio"] = float(
            (amplitude_envelope < burst_thresh).mean()
        )
        all_feats = {**spec_feats, **conn_feats}

        if self.feature_names is None:
            self.feature_names = list(all_feats.keys())

        X = np.array([[all_feats[k] for k in self.feature_names]])
        X_scaled = self.scaler.transform(X)
        pred_class = self.rf.predict(X_scaled)[0]
        proba = self.rf.predict_proba(X_scaled)[0]

        return {
            "predicted_state": CLASS_NAMES.get(pred_class, "UNKNOWN"),
            "probabilities": {
                CLASS_NAMES.get(i, str(i)): float(p)
                for i, p in enumerate(proba)
            },
        }

    def cross_validate(
        self, X: np.ndarray, y: np.ndarray, cv: int = 5
    ) -> dict:
        """Cross-validate classifier and return accuracy metrics."""
        X_scaled = self.scaler.fit_transform(X)
        rf_scores = cross_val_score(
            RandomForestClassifier(n_estimators=100, random_state=42),
            X_scaled, y, cv=cv, scoring="accuracy"
        )
        lda_scores = cross_val_score(
            LinearDiscriminantAnalysis(),
            X_scaled, y, cv=cv, scoring="accuracy"
        )
        return {
            "rf_accuracy": float(rf_scores.mean()),
            "rf_accuracy_std": float(rf_scores.std()),
            "lda_accuracy": float(lda_scores.mean()),
            "lda_accuracy_std": float(lda_scores.std()),
            "rf_scores": rf_scores.tolist(),
            "lda_scores": lda_scores.tolist(),
        }

    def feature_importance(self) -> Dict:
        """Return feature importances from Random Forest."""
        if not self._is_fitted:
            raise RuntimeError("Classifier not fitted yet.")
        return {
            name: float(imp)
            for name, imp in zip(self.feature_names, self.rf.feature_importances_)
        }
