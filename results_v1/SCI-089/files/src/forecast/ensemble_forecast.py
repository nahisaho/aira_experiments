"""Ensemble calibration and forecast evaluation utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import stats


def _to_float_array(values: ArrayLike) -> NDArray[np.float64]:
    """Convert array-like input to a floating NumPy array."""
    return np.asarray(values, dtype=float)


@dataclass
class ForecastMetrics:
    """Collection of deterministic and probabilistic forecast metrics."""

    def compute_deterministic_metrics(
        self,
        predictions: ArrayLike,
        observations: ArrayLike,
    ) -> dict[str, float]:
        """Compute MAE, RMSE, and bias for deterministic forecasts."""
        pred = _to_float_array(predictions).reshape(-1)
        obs = _to_float_array(observations).reshape(-1)
        error = pred - obs
        return {
            "mae": float(np.mean(np.abs(error))),
            "rmse": float(np.sqrt(np.mean(error ** 2))),
            "bias": float(np.mean(error)),
        }

    def compute_probabilistic_metrics(
        self,
        quantile_predictions: Mapping[str, ArrayLike],
        observations: ArrayLike,
    ) -> dict[str, float]:
        """Compute reliability, sharpness, resolution, and pinball loss."""
        obs = _to_float_array(observations).reshape(-1)
        quantile_levels, quantile_values = self._extract_quantiles(quantile_predictions)
        coverage = []
        pinball_terms = []
        for level, values in zip(quantile_levels, quantile_values.T, strict=False):
            values = values.reshape(-1)
            coverage.append(np.mean(obs <= values) - level)
            diff = obs - values
            pinball_terms.append(np.mean(np.maximum(level * diff, (level - 1.0) * diff)))

        lower = np.asarray(
            [np.interp(0.1, quantile_levels, row) for row in quantile_values],
            dtype=float,
        )
        upper = np.asarray(
            [np.interp(0.9, quantile_levels, row) for row in quantile_values],
            dtype=float,
        )
        sharpness = float(np.mean(upper - lower))
        median = np.asarray(
            [np.interp(0.5, quantile_levels, row) for row in quantile_values],
            dtype=float,
        )
        resolution = float(np.var(median))
        return {
            "reliability": float(np.mean(np.abs(coverage))),
            "sharpness": sharpness,
            "resolution": resolution,
            "pinball_loss": float(np.mean(pinball_terms)),
        }

    def _extract_quantiles(
        self,
        quantile_predictions: Mapping[str, ArrayLike],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Extract sorted quantile levels and a value matrix from a mapping."""
        entries = []
        for key, value in quantile_predictions.items():
            if key.startswith("q"):
                entries.append((float(key[1:]) / 100.0, _to_float_array(value).reshape(-1)))
        if not entries:
            raise KeyError("At least one quantile key such as 'q10' is required.")
        entries.sort(key=lambda item: item[0])
        levels = np.asarray([item[0] for item in entries], dtype=float)
        values = np.vstack([item[1] for item in entries]).T
        return levels, values


class EnsembleForecastSystem:
    """Calibrate ensemble forecasts and evaluate probabilistic skill.

    The calibration stage uses a Gaussian EMOS formulation where the ensemble
    mean and spread are linearly mapped to predictive mean and standard
    deviation parameters.
    """

    def __init__(self) -> None:
        self.calibration_: dict[str, float] | None = None
        self.metrics = ForecastMetrics()

    def calibrate(
        self,
        observations: ArrayLike,
        raw_ensemble: ArrayLike,
    ) -> dict[str, NDArray[np.float64] | dict[str, float]]:
        """Calibrate a raw ensemble using a simple EMOS regression."""
        obs = _to_float_array(observations).reshape(-1)
        ensemble = _to_float_array(raw_ensemble)
        if ensemble.ndim != 2 or ensemble.shape[0] != obs.size:
            raise ValueError("raw_ensemble must have shape (n_samples, n_members).")

        ens_mean = ensemble.mean(axis=1)
        ens_std = ensemble.std(axis=1, ddof=1)
        design = np.column_stack([np.ones(obs.size), ens_mean, ens_std])
        coefficients, *_ = np.linalg.lstsq(design, obs, rcond=None)
        calibrated_mean = design @ coefficients
        residual = obs - calibrated_mean
        variance_design = np.column_stack([np.ones(obs.size), ens_std ** 2])
        variance_coef, *_ = np.linalg.lstsq(variance_design, residual ** 2, rcond=None)
        calibrated_std = np.sqrt(np.maximum(variance_design @ variance_coef, 1e-6))

        self.calibration_ = {
            "a": float(coefficients[0]),
            "b": float(coefficients[1]),
            "c": float(coefficients[2]),
            "v0": float(variance_coef[0]),
            "v1": float(variance_coef[1]),
        }
        return {
            "mean": calibrated_mean,
            "std": calibrated_std,
            "parameters": self.calibration_,
        }

    def evaluate(
        self,
        predictions: Mapping[str, ArrayLike],
        observations: ArrayLike,
    ) -> dict[str, Any]:
        """Evaluate calibrated forecasts with deterministic and probabilistic metrics."""
        obs = _to_float_array(observations).reshape(-1)
        mean_prediction = _to_float_array(predictions["mean"]).reshape(-1)
        std_prediction = np.maximum(_to_float_array(predictions["std"]).reshape(-1), 1e-6)
        deterministic = self.metrics.compute_deterministic_metrics(mean_prediction, obs)

        quantile_levels = np.array([0.1, 0.25, 0.5, 0.75, 0.9], dtype=float)
        quantile_predictions = {
            f"q{int(level * 100):02d}": stats.norm.ppf(level, loc=mean_prediction, scale=std_prediction)
            for level in quantile_levels
        }
        probabilistic = self.metrics.compute_probabilistic_metrics(quantile_predictions, obs)
        reliability = self._compute_reliability_diagram(mean_prediction, std_prediction, obs)
        crps = self._gaussian_crps(mean_prediction, std_prediction, obs)
        skill_scores = self._compute_skill_scores(mean_prediction, obs)
        return {
            "deterministic": deterministic,
            "probabilistic": probabilistic,
            "reliability_diagram": reliability,
            "crps": float(np.mean(crps)),
            "skill_scores": skill_scores,
        }

    def _compute_reliability_diagram(
        self,
        mean: NDArray[np.float64],
        std: NDArray[np.float64],
        observations: NDArray[np.float64],
        thresholds: NDArray[np.float64] | None = None,
    ) -> dict[str, NDArray[np.float64]]:
        """Compute observed frequency versus forecast probability pairs."""
        if thresholds is None:
            thresholds = np.quantile(observations, np.linspace(0.1, 0.9, 9))
        forecast_probabilities = np.column_stack(
            [1.0 - stats.norm.cdf(threshold, loc=mean, scale=std) for threshold in thresholds]
        )
        observed_events = np.column_stack([(observations >= threshold).astype(float) for threshold in thresholds])
        bins = np.linspace(0.0, 1.0, 11)
        binned_forecast = np.full((len(thresholds), len(bins) - 1), np.nan, dtype=float)
        binned_observed = np.full_like(binned_forecast, np.nan)

        for threshold_idx in range(len(thresholds)):
            probs = forecast_probabilities[:, threshold_idx]
            events = observed_events[:, threshold_idx]
            for bin_idx in range(len(bins) - 1):
                mask = (probs >= bins[bin_idx]) & (probs < bins[bin_idx + 1])
                if np.any(mask):
                    binned_forecast[threshold_idx, bin_idx] = probs[mask].mean()
                    binned_observed[threshold_idx, bin_idx] = events[mask].mean()
        return {
            "thresholds": thresholds,
            "forecast_probability": binned_forecast,
            "observed_frequency": binned_observed,
        }

    def _gaussian_crps(
        self,
        mean: NDArray[np.float64],
        std: NDArray[np.float64],
        observations: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Compute the CRPS for Gaussian predictive distributions."""
        z = (observations - mean) / std
        return std * (z * (2.0 * stats.norm.cdf(z) - 1.0) + 2.0 * stats.norm.pdf(z) - 1.0 / np.sqrt(np.pi))

    def _compute_skill_scores(
        self,
        predictions: NDArray[np.float64],
        observations: NDArray[np.float64],
    ) -> dict[str, float]:
        """Compute skill scores against persistence and climatology baselines."""
        persistence = np.roll(observations, 1)
        persistence[0] = observations[0]
        climatology = np.full_like(observations, observations.mean())
        mse_forecast = np.mean((predictions - observations) ** 2)
        mse_persistence = np.mean((persistence - observations) ** 2)
        mse_climatology = np.mean((climatology - observations) ** 2)
        return {
            "against_persistence": float(1.0 - mse_forecast / max(mse_persistence, 1e-6)),
            "against_climatology": float(1.0 - mse_forecast / max(mse_climatology, 1e-6)),
        }
