"""
evaluation.py — Evaluation metrics and comparison framework.

Implements standard weather prediction metrics:
  - RMSE (Root Mean Square Error)
  - ACC (Anomaly Correlation Coefficient)
  - MAE (Mean Absolute Error)
  - SSIM-like spatial structure scores
  - Physics consistency scores
  - Comparison against GFS/ECMWF baselines (simulated)
"""

import numpy as np
import torch
from typing import Dict, List, Tuple, Optional
import json


VARIABLE_NAMES = {
    'surface': ['MSLP (Pa)', 'T2m (K)', 'U10 (m/s)', 'V10 (m/s)', 'TP (mm)'],
    'pressure': ['T (K)', 'U (m/s)', 'V (m/s)', 'q (kg/kg)', 'Z (m²/s²)'],
}

PRESSURE_LEVELS = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]


def latitude_weights(lat: np.ndarray) -> np.ndarray:
    """Compute area-weighted factors based on latitude."""
    cos_lat = np.cos(np.deg2rad(lat))
    return cos_lat / cos_lat.mean()


def rmse(pred: np.ndarray, target: np.ndarray,
         weights: Optional[np.ndarray] = None) -> float:
    """Weighted Root Mean Square Error."""
    diff2 = (pred - target) ** 2
    if weights is not None:
        diff2 = diff2 * weights.reshape(-1, 1) if diff2.ndim > 1 else diff2 * weights
    return float(np.sqrt(np.mean(diff2)))


def mae(pred: np.ndarray, target: np.ndarray,
        weights: Optional[np.ndarray] = None) -> float:
    """Weighted Mean Absolute Error."""
    diff = np.abs(pred - target)
    if weights is not None:
        diff = diff * weights.reshape(-1, 1) if diff.ndim > 1 else diff * weights
    return float(np.mean(diff))


def anomaly_correlation(pred: np.ndarray, target: np.ndarray,
                        climatology: np.ndarray,
                        weights: Optional[np.ndarray] = None) -> float:
    """
    Anomaly Correlation Coefficient (ACC).

    ACC = Σ w*(pred' * target') / sqrt(Σ w*pred'² * Σ w*target'²)
    where pred' = pred - climatology, target' = target - climatology
    """
    pred_anom = pred - climatology
    target_anom = target - climatology

    if weights is None:
        weights = np.ones(pred.shape[0])
    w = weights.reshape(-1, 1) if pred_anom.ndim > 1 else weights

    num = np.sum(w * pred_anom * target_anom)
    denom = np.sqrt(np.sum(w * pred_anom**2) * np.sum(w * target_anom**2))
    return float(num / (denom + 1e-10))


def spatial_correlation(pred: np.ndarray, target: np.ndarray) -> float:
    """Pearson correlation of spatial patterns."""
    pred_flat = pred.flatten()
    target_flat = target.flatten()
    corr = np.corrcoef(pred_flat, target_flat)[0, 1]
    return float(corr) if not np.isnan(corr) else 0.0


def energy_spectrum_score(pred: np.ndarray, target: np.ndarray,
                          n_lat: int, n_lon: int) -> float:
    """
    Compare energy spectra (simplified).
    Measures how well the model preserves multi-scale variability.
    """
    pred_2d = pred.reshape(n_lat, n_lon)
    target_2d = target.reshape(n_lat, n_lon)

    pred_fft = np.abs(np.fft.fft2(pred_2d)) ** 2
    target_fft = np.abs(np.fft.fft2(target_2d)) ** 2

    pred_spectrum = np.sort(pred_fft.flatten())[::-1][:20]
    target_spectrum = np.sort(target_fft.flatten())[::-1][:20]

    pred_norm = pred_spectrum / (pred_spectrum.sum() + 1e-10)
    target_norm = target_spectrum / (target_spectrum.sum() + 1e-10)

    return float(1.0 - 0.5 * np.sum(np.abs(pred_norm - target_norm)))


def physics_consistency_score(
    pred_state: np.ndarray,
    n_levels: int = 13,
) -> Dict[str, float]:
    """
    Evaluate physics consistency of predicted state.

    Checks:
      1. Humidity non-negativity
      2. Temperature lapse rate reasonableness
      3. Wind magnitude reasonableness
    """
    n_surface = 5
    n_vars_per_level = 5

    # Extract variables
    pressure_data = pred_state[:, n_surface:]
    temp = pressure_data[:, :n_levels]
    u = pressure_data[:, n_levels:2*n_levels]
    v = pressure_data[:, 2*n_levels:3*n_levels]
    q = pressure_data[:, 3*n_levels:4*n_levels]

    # Humidity non-negativity
    humidity_violation = float(np.mean(q < 0) * 100)

    # Temperature lapse rate (should generally decrease with height in troposphere)
    temp_diff = np.diff(temp[:, 6:], axis=1)  # from 400 hPa downward
    lapse_rate_ok = float(np.mean(temp_diff > 0) * 100)

    # Wind speed reasonableness (< 150 m/s)
    wind_speed = np.sqrt(u**2 + v**2)
    wind_ok = float(np.mean(wind_speed < 150) * 100)

    return {
        'humidity_valid_pct': 100.0 - humidity_violation,
        'lapse_rate_ok_pct': lapse_rate_ok,
        'wind_reasonable_pct': wind_ok,
        'overall_physics_score': (
            (100 - humidity_violation) + lapse_rate_ok + wind_ok
        ) / 3,
    }


class NWPBaseline:
    """
    Simulated NWP baseline scores for comparison.

    Based on published RMSE values from ECMWF IFS and NCEP GFS
    for Z500, T850, and surface variables at various lead times.
    """

    # Published approximate RMSE values (simplified)
    BASELINES = {
        'ECMWF_IFS': {
            6: {'z500_rmse': 15.0, 't850_rmse': 0.5, 'u850_rmse': 1.0, 'acc_z500': 0.998},
            24: {'z500_rmse': 45.0, 't850_rmse': 1.0, 'u850_rmse': 2.0, 'acc_z500': 0.990},
            120: {'z500_rmse': 180.0, 't850_rmse': 2.5, 'u850_rmse': 4.0, 'acc_z500': 0.880},
        },
        'GFS': {
            6: {'z500_rmse': 20.0, 't850_rmse': 0.6, 'u850_rmse': 1.2, 'acc_z500': 0.996},
            24: {'z500_rmse': 55.0, 't850_rmse': 1.3, 'u850_rmse': 2.5, 'acc_z500': 0.985},
            120: {'z500_rmse': 220.0, 't850_rmse': 3.0, 'u850_rmse': 5.0, 'acc_z500': 0.840},
        },
        'GraphCast_published': {
            6: {'z500_rmse': 12.0, 't850_rmse': 0.45, 'u850_rmse': 0.9, 'acc_z500': 0.999},
            24: {'z500_rmse': 38.0, 't850_rmse': 0.9, 'u850_rmse': 1.8, 'acc_z500': 0.993},
            120: {'z500_rmse': 155.0, 't850_rmse': 2.2, 'u850_rmse': 3.5, 'acc_z500': 0.900},
        },
    }

    @classmethod
    def get_baselines(cls, lead_time_hours: int) -> Dict[str, Dict[str, float]]:
        """Get baseline metrics for a given lead time."""
        result = {}
        for model_name, lt_data in cls.BASELINES.items():
            if lead_time_hours in lt_data:
                result[model_name] = lt_data[lead_time_hours]
        return result


class WeatherEvaluator:
    """Full evaluation pipeline for weather prediction models."""

    def __init__(self, lat: np.ndarray, lon: np.ndarray,
                 n_lat: int = 46, n_lon: int = 90):
        self.lat = lat
        self.lon = lon
        self.n_lat = n_lat
        self.n_lon = n_lon
        self.weights = latitude_weights(lat)

    def evaluate_forecast(
        self,
        pred: np.ndarray,
        target: np.ndarray,
        climatology: np.ndarray,
        lead_time_hours: int,
    ) -> Dict:
        """
        Full evaluation of a single forecast.

        Args:
            pred: (N, F) predicted state
            target: (N, F) target state
            climatology: (N, F) climatological mean
            lead_time_hours: forecast lead time in hours

        Returns:
            Dict with all metrics
        """
        n_surface = 5
        n_levels = 13

        metrics = {
            'lead_time_hours': lead_time_hours,
        }

        # Overall RMSE/MAE
        metrics['overall_rmse'] = rmse(pred, target, self.weights)
        metrics['overall_mae'] = mae(pred, target, self.weights)

        # Surface variable metrics
        for i, name in enumerate(VARIABLE_NAMES['surface']):
            metrics[f'rmse_{name}'] = rmse(
                pred[:, i], target[:, i], self.weights
            )

        # Key pressure level metrics (Z500, T850)
        z_start = n_surface + 4 * n_levels  # geopotential offset
        z500_idx = z_start + 7  # 500 hPa is index 7
        t_start = n_surface
        t850_idx = t_start + 10  # 850 hPa is index 10

        metrics['z500_rmse'] = rmse(pred[:, z500_idx], target[:, z500_idx], self.weights)
        metrics['t850_rmse'] = rmse(pred[:, t850_idx], target[:, t850_idx], self.weights)

        u_start = n_surface + n_levels
        u850_idx = u_start + 10
        metrics['u850_rmse'] = rmse(pred[:, u850_idx], target[:, u850_idx], self.weights)

        # ACC for Z500
        metrics['acc_z500'] = anomaly_correlation(
            pred[:, z500_idx], target[:, z500_idx],
            climatology[:, z500_idx], self.weights
        )

        # Spatial correlation
        metrics['spatial_corr_z500'] = spatial_correlation(
            pred[:, z500_idx], target[:, z500_idx]
        )

        # Energy spectrum score
        metrics['spectrum_score_z500'] = energy_spectrum_score(
            pred[:, z500_idx], target[:, z500_idx], self.n_lat, self.n_lon
        )

        # Physics consistency
        physics = physics_consistency_score(pred, n_levels)
        metrics.update({f'physics_{k}': v for k, v in physics.items()})

        # Compare with NWP baselines
        baselines = NWPBaseline.get_baselines(lead_time_hours)
        metrics['nwp_baselines'] = baselines

        return metrics

    def evaluate_multi_step(
        self,
        predictions: np.ndarray,
        targets: np.ndarray,
        climatology: np.ndarray,
        lead_times_hours: List[int],
    ) -> List[Dict]:
        """Evaluate multiple forecast lead times."""
        results = []
        for i, lt in enumerate(lead_times_hours):
            result = self.evaluate_forecast(
                predictions[i], targets[i], climatology, lt
            )
            results.append(result)
        return results

    @staticmethod
    def format_comparison_table(results: List[Dict]) -> str:
        """Format results as a comparison table."""
        lines = []
        lines.append("=" * 90)
        lines.append(f"{'Lead Time':>12s} | {'Z500 RMSE':>10s} | {'T850 RMSE':>10s} | "
                      f"{'U850 RMSE':>10s} | {'ACC Z500':>10s} | {'Physics':>10s}")
        lines.append("-" * 90)

        for r in results:
            lt = r['lead_time_hours']
            lines.append(
                f"{'Our Model':>12s} | "
                f"{r.get('z500_rmse', 0):>10.2f} | "
                f"{r.get('t850_rmse', 0):>10.3f} | "
                f"{r.get('u850_rmse', 0):>10.3f} | "
                f"{r.get('acc_z500', 0):>10.4f} | "
                f"{r.get('physics_overall_physics_score', 0):>9.1f}%"
            )

            # Add baseline rows
            for model_name, baseline in r.get('nwp_baselines', {}).items():
                lines.append(
                    f"{model_name:>12s} | "
                    f"{baseline.get('z500_rmse', 0):>10.2f} | "
                    f"{baseline.get('t850_rmse', 0):>10.3f} | "
                    f"{baseline.get('u850_rmse', 0):>10.3f} | "
                    f"{baseline.get('acc_z500', 0):>10.4f} | "
                    f"{'N/A':>10s}"
                )
            lines.append("-" * 90)

        lines.append("=" * 90)
        return "\n".join(lines)
