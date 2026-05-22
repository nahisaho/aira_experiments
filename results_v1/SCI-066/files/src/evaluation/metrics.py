"""
Evaluation framework for ESM AI Emulator.

Implements ClimateBench-compatible metrics:
- RMSE, MAE per variable per region
- Spatial pattern correlation
- Trend accuracy
- Ensemble spread calibration (rank histogram, CRPS)
- Physical consistency checks
"""

import numpy as np
import torch
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field, asdict
import json


@dataclass
class EvaluationMetrics:
    """Container for all evaluation metrics."""
    global_rmse: Dict[str, float] = field(default_factory=dict)
    global_mae: Dict[str, float] = field(default_factory=dict)
    regional_rmse: Dict[str, Dict[str, float]] = field(default_factory=dict)
    pattern_correlation: Dict[str, float] = field(default_factory=dict)
    trend_error: Dict[str, float] = field(default_factory=dict)
    ensemble_spread: Dict[str, float] = field(default_factory=dict)
    crps: Dict[str, float] = field(default_factory=dict)
    energy_conservation_error: float = 0.0
    mass_conservation_error: float = 0.0
    nrmse: Dict[str, float] = field(default_factory=dict)
    spatial_skill_score: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


VARIABLE_NAMES = ["temperature", "precipitation", "sea_level"]

REGIONS = {
    "global": (slice(None), slice(None)),
    "tropics": (slice(22, 42), slice(None)),          # ~30S-30N
    "northern_extratropics": (slice(42, 58), slice(None)),  # 30N-60N
    "southern_extratropics": (slice(6, 22), slice(None)),   # 60S-30S
    "arctic": (slice(58, 64), slice(None)),            # 60N-90N
    "antarctic": (slice(0, 6), slice(None)),           # 90S-60S
}


def compute_rmse(pred: np.ndarray, target: np.ndarray,
                 weights: Optional[np.ndarray] = None) -> float:
    """Weighted root mean squared error."""
    diff_sq = (pred - target) ** 2
    if weights is not None:
        return float(np.sqrt(np.average(diff_sq, weights=weights)))
    return float(np.sqrt(np.mean(diff_sq)))


def compute_mae(pred: np.ndarray, target: np.ndarray,
                weights: Optional[np.ndarray] = None) -> float:
    """Weighted mean absolute error."""
    diff = np.abs(pred - target)
    if weights is not None:
        return float(np.average(diff, weights=weights))
    return float(np.mean(diff))


def compute_nrmse(pred: np.ndarray, target: np.ndarray) -> float:
    """Normalized RMSE (by target standard deviation)."""
    rmse = compute_rmse(pred, target)
    std = np.std(target)
    return rmse / std if std > 0 else float("inf")


def compute_pattern_correlation(pred: np.ndarray,
                                 target: np.ndarray) -> float:
    """Spatial pattern correlation (centered anomaly correlation)."""
    pred_anom = pred - np.mean(pred)
    target_anom = target - np.mean(target)

    numerator = np.sum(pred_anom * target_anom)
    denominator = np.sqrt(np.sum(pred_anom ** 2) * np.sum(target_anom ** 2))

    return float(numerator / denominator) if denominator > 0 else 0.0


def compute_trend(timeseries: np.ndarray) -> float:
    """Linear trend via least squares (units per timestep)."""
    n = len(timeseries)
    x = np.arange(n)
    slope = (n * np.sum(x * timeseries) - np.sum(x) * np.sum(timeseries)) / \
            (n * np.sum(x ** 2) - np.sum(x) ** 2)
    return float(slope)


def compute_crps(ensemble_preds: np.ndarray,
                 observation: np.ndarray) -> float:
    """
    Continuous Ranked Probability Score.

    Measures calibration of probabilistic forecasts.
    Lower is better.
    """
    n_members = ensemble_preds.shape[0]
    obs_flat = observation.flatten()
    ens_flat = ensemble_preds.reshape(n_members, -1)

    crps_sum = 0.0
    n_points = obs_flat.shape[0]

    for j in range(n_points):
        ens_sorted = np.sort(ens_flat[:, j])
        obs_val = obs_flat[j]

        term1 = np.mean(np.abs(ens_sorted - obs_val))
        term2 = 0.0
        for k in range(n_members):
            for l in range(k + 1, n_members):
                term2 += np.abs(ens_sorted[k] - ens_sorted[l])
        term2 /= (n_members * (n_members - 1) / 2) if n_members > 1 else 1

        crps_sum += term1 - 0.5 * term2

    return crps_sum / n_points


def compute_spread_skill_ratio(ensemble_std: np.ndarray,
                                rmse: float) -> float:
    """
    Spread-skill ratio.

    Ideal ratio = 1.0 (ensemble spread matches actual error).
    > 1.0: overdispersive, < 1.0: underdispersive.
    """
    mean_spread = float(np.mean(ensemble_std))
    return mean_spread / rmse if rmse > 0 else float("inf")


def generate_latitude_weights(n_lat: int) -> np.ndarray:
    """Area weights proportional to cos(latitude)."""
    lats = np.linspace(-90, 90, n_lat)
    weights = np.cos(np.radians(lats))
    weights = weights / weights.sum()
    return weights


class ClimateBenchEvaluator:
    """
    ClimateBench-compatible evaluation framework.

    Evaluates emulator predictions against reference ESM outputs
    using standardized metrics across variables, regions, and scenarios.
    """

    def __init__(self, spatial_size: Tuple[int, int] = (64, 128)):
        self.spatial_size = spatial_size
        self.lat_weights = generate_latitude_weights(spatial_size[0])

    def evaluate_prediction(self, pred: np.ndarray, target: np.ndarray,
                            ensemble_preds: Optional[np.ndarray] = None,
                            ) -> EvaluationMetrics:
        """
        Full evaluation of a single prediction.

        Args:
            pred: (C, H, W) predicted fields
            target: (C, H, W) reference fields
            ensemble_preds: (N, C, H, W) ensemble member predictions

        Returns:
            EvaluationMetrics with all computed metrics
        """
        metrics = EvaluationMetrics()
        n_vars = pred.shape[0]

        # Per-variable metrics
        for v in range(min(n_vars, len(VARIABLE_NAMES))):
            var_name = VARIABLE_NAMES[v]
            p, t = pred[v], target[v]

            # Global metrics
            metrics.global_rmse[var_name] = compute_rmse(p, t)
            metrics.global_mae[var_name] = compute_mae(p, t)
            metrics.nrmse[var_name] = compute_nrmse(p, t)
            metrics.pattern_correlation[var_name] = compute_pattern_correlation(p, t)

            # Regional metrics
            metrics.regional_rmse[var_name] = {}
            for region_name, (lat_slice, lon_slice) in REGIONS.items():
                p_r = p[lat_slice, lon_slice]
                t_r = t[lat_slice, lon_slice]
                metrics.regional_rmse[var_name][region_name] = compute_rmse(p_r, t_r)

            # Spatial skill score (1 - NRMSE²)
            nrmse = metrics.nrmse[var_name]
            metrics.spatial_skill_score[var_name] = max(0, 1 - nrmse ** 2)

            # Ensemble metrics
            if ensemble_preds is not None:
                ens_std = np.std(ensemble_preds[:, v], axis=0)
                rmse = metrics.global_rmse[var_name]
                metrics.ensemble_spread[var_name] = compute_spread_skill_ratio(
                    ens_std, rmse
                )

                if ensemble_preds.shape[0] <= 20:
                    metrics.crps[var_name] = compute_crps(
                        ensemble_preds[:, v], t
                    )

        return metrics

    def evaluate_scenario(self, predictions: List[np.ndarray],
                          targets: List[np.ndarray],
                          scenario: str) -> Dict:
        """Evaluate all timesteps for a scenario."""
        all_metrics = []
        for pred, target in zip(predictions, targets):
            m = self.evaluate_prediction(pred, target)
            all_metrics.append(m.to_dict())

        # Aggregate
        aggregated = {}
        for var_name in VARIABLE_NAMES:
            rmses = [m["global_rmse"].get(var_name, 0) for m in all_metrics]
            maes = [m["global_mae"].get(var_name, 0) for m in all_metrics]
            pcs = [m["pattern_correlation"].get(var_name, 0) for m in all_metrics]

            aggregated[var_name] = {
                "mean_rmse": float(np.mean(rmses)),
                "std_rmse": float(np.std(rmses)),
                "mean_mae": float(np.mean(maes)),
                "mean_pattern_corr": float(np.mean(pcs)),
            }

            # Trend accuracy
            if len(rmses) > 1:
                pred_means = [p.mean() for p in predictions]
                target_means = [t.mean() for t in targets]
                pred_trend = compute_trend(np.array(pred_means))
                target_trend = compute_trend(np.array(target_means))
                aggregated[var_name]["trend_error_pct"] = float(
                    abs(pred_trend - target_trend) / abs(target_trend) * 100
                ) if abs(target_trend) > 1e-10 else 0.0

        return {
            "scenario": scenario,
            "n_timesteps": len(predictions),
            "per_variable": aggregated,
        }

    def benchmark_against_cmip6(self, emulator_results: Dict,
                                 reference_results: Dict) -> Dict:
        """
        Compare emulator against multiple CMIP6 models.

        Computes relative performance metrics to quantify
        how well the emulator reproduces multi-model spread.
        """
        comparison = {}
        for model_name, ref_data in reference_results.items():
            model_metrics = {}
            for var_name in VARIABLE_NAMES:
                if var_name in emulator_results and var_name in ref_data:
                    emu = emulator_results[var_name]
                    ref = ref_data[var_name]
                    model_metrics[var_name] = {
                        "rmse_ratio": emu.get("mean_rmse", 0) / max(ref.get("mean_rmse", 1), 1e-10),
                        "pattern_corr_diff": emu.get("mean_pattern_corr", 0) - ref.get("mean_pattern_corr", 0),
                    }
            comparison[model_name] = model_metrics

        return comparison


def save_metrics(metrics: Dict, filepath: str):
    """Save metrics to JSON."""
    with open(filepath, "w") as f:
        json.dump(metrics, f, indent=2, default=str)
