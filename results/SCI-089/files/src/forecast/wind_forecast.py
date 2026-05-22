"""Probabilistic wind power forecasting utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import stats
from sklearn.neighbors import KernelDensity


def _as_float_array(values: ArrayLike, ndim: int | None = None) -> NDArray[np.float64]:
    """Convert array-like input into a float ``ndarray``."""
    array = np.asarray(values, dtype=float)
    if ndim is not None and array.ndim != ndim:
        raise ValueError(f"Expected {ndim} dimensions, received {array.ndim}.")
    return array.astype(float, copy=False)


@dataclass
class _NumpyLSTMCell:
    """Small PyTorch-style LSTM cell implemented with NumPy only."""

    input_size: int
    hidden_size: int
    rng: np.random.Generator

    def __post_init__(self) -> None:
        scale = 1.0 / np.sqrt(self.input_size + self.hidden_size)
        self.w_ih = self.rng.normal(0.0, scale, size=(4 * self.hidden_size, self.input_size))
        self.w_hh = self.rng.normal(0.0, scale, size=(4 * self.hidden_size, self.hidden_size))
        self.b = np.zeros(4 * self.hidden_size, dtype=float)

    def forward(self, sequence: NDArray[np.float64]) -> NDArray[np.float64]:
        """Run a sequence through the cell and return the final hidden state."""
        h = np.zeros(self.hidden_size, dtype=float)
        c = np.zeros(self.hidden_size, dtype=float)
        for x_t in sequence:
            gates = self.w_ih @ x_t + self.w_hh @ h + self.b
            i, f, g, o = np.split(gates, 4)
            i = 1.0 / (1.0 + np.exp(-i))
            f = 1.0 / (1.0 + np.exp(-f))
            g = np.tanh(g)
            o = 1.0 / (1.0 + np.exp(-o))
            c = f * c + i * g
            h = o * np.tanh(c)
        return h


class NWPWindPreprocessor:
    """Preprocess NWP wind data and extrapolate it to hub height."""

    def __init__(
        self,
        hub_height: float = 100.0,
        reference_height: float = 10.0,
        roughness_length: float = 0.1,
        shear_exponent: float = 0.14,
    ) -> None:
        self.hub_height = float(hub_height)
        self.reference_height = float(reference_height)
        self.roughness_length = float(roughness_length)
        self.shear_exponent = float(shear_exponent)

    def preprocess(self, nwp_data: Mapping[str, ArrayLike]) -> dict[str, NDArray[np.float64]]:
        """Prepare wind fields for machine-learning models."""
        required = ("u", "v")
        missing = [key for key in required if key not in nwp_data]
        if missing:
            raise KeyError(f"Missing required wind fields: {missing}")
        u = _as_float_array(nwp_data["u"])
        v = _as_float_array(nwp_data["v"])
        speed = self.compute_wind_speed(u, v)
        direction = (np.degrees(np.arctan2(v, u)) + 360.0) % 360.0
        hub_speed = self.extrapolate_to_hub_height(speed)
        corrected_speed = self._apply_terrain_roughness_correction(hub_speed)
        return {
            "u": u,
            "v": v,
            "wind_speed": speed,
            "wind_direction": direction,
            "hub_height_speed": corrected_speed,
        }

    def extrapolate_to_hub_height(
        self,
        wind_speed: ArrayLike,
        method: str = "power_law",
    ) -> NDArray[np.float64]:
        """Extrapolate wind speed from reference to hub height."""
        speed = _as_float_array(wind_speed)
        method_normalized = method.lower()
        if method_normalized == "power_law":
            factor = (self.hub_height / self.reference_height) ** self.shear_exponent
            return speed * factor
        if method_normalized == "log_law":
            numerator = np.log(self.hub_height / self.roughness_length)
            denominator = np.log(self.reference_height / self.roughness_length)
            return speed * numerator / denominator
        raise ValueError("method must be either 'power_law' or 'log_law'.")

    def compute_wind_speed(self, u_component: ArrayLike, v_component: ArrayLike) -> NDArray[np.float64]:
        """Compute scalar wind speed from horizontal components."""
        u = _as_float_array(u_component)
        v = _as_float_array(v_component)
        return np.hypot(u, v)

    def _apply_terrain_roughness_correction(
        self,
        wind_speed: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Apply a first-order terrain roughness correction."""
        correction = 1.0 - 0.08 * np.tanh(5.0 * self.roughness_length)
        return np.maximum(wind_speed * correction, 0.0)


class WindMLPredictor:
    """Wind power predictor with regime-dependent modeling and KDE outputs."""

    def __init__(
        self,
        hidden_size: int = 16,
        bandwidth: float = 0.08,
        random_state: int = 42,
    ) -> None:
        self.hidden_size = int(hidden_size)
        self.bandwidth = float(bandwidth)
        self.random_state = int(random_state)
        self.rng = np.random.default_rng(random_state)
        self.encoder: _NumpyLSTMCell | None = None
        self.regime_models: dict[str, NDArray[np.float64]] = {}
        self.kde_models: dict[str, KernelDensity] = {}
        self.cut_in = 3.0
        self.rated = 12.0
        self.cut_out = 25.0

    def train(self, historical_data: Mapping[str, ArrayLike]) -> None:
        """Train regime-dependent linear readouts on top of a NumPy LSTM encoder."""
        sequences = _as_float_array(historical_data["features"], ndim=3)
        target = _as_float_array(historical_data["power"], ndim=1)
        wind_speed = _as_float_array(historical_data["wind_speed"], ndim=1)
        if sequences.shape[0] != target.shape[0] or target.shape[0] != wind_speed.shape[0]:
            raise ValueError("features, power, and wind_speed must share the sample dimension.")

        self.encoder = _NumpyLSTMCell(sequences.shape[-1], self.hidden_size, self.rng)
        latent = np.vstack([self.encoder.forward(sequence) for sequence in sequences])
        design = np.hstack([latent, wind_speed[:, None], np.ones((len(target), 1))])
        regimes = self._classify_regime(wind_speed)

        for regime in ("low", "medium", "high"):
            mask = regimes == regime
            if not np.any(mask):
                continue
            coefficients, *_ = np.linalg.lstsq(design[mask], target[mask], rcond=None)
            prediction = design[mask] @ coefficients
            residuals = target[mask] - prediction
            kde = KernelDensity(kernel="gaussian", bandwidth=self.bandwidth)
            kde.fit(residuals[:, None])
            self.regime_models[regime] = coefficients
            self.kde_models[regime] = kde

    def predict(self, features: Mapping[str, ArrayLike]) -> NDArray[np.float64]:
        """Generate deterministic wind power forecasts with hysteresis-aware power curves."""
        encoded, wind_speed = self._encode_features(features)
        regimes = self._classify_regime(wind_speed)
        predictions = np.zeros(wind_speed.shape[0], dtype=float)
        previous_speed = _as_float_array(features.get("previous_wind_speed", wind_speed), ndim=1)

        for idx, regime in enumerate(regimes):
            coefficients = self.regime_models.get(regime)
            if coefficients is None:
                coefficients = next(iter(self.regime_models.values()))
            design_row = np.concatenate([encoded[idx], [wind_speed[idx], 1.0]])
            base_prediction = float(design_row @ coefficients)
            curve_prediction = self._power_curve_with_hysteresis(wind_speed[idx], previous_speed[idx])
            predictions[idx] = 0.5 * base_prediction + 0.5 * curve_prediction
        return np.clip(predictions, 0.0, 1.2)

    def predict_distribution(self, features: Mapping[str, ArrayLike]) -> dict[str, NDArray[np.float64]]:
        """Estimate predictive densities with regime-specific kernel density models."""
        mean_prediction = self.predict(features)
        encoded, wind_speed = self._encode_features(features)
        regimes = self._classify_regime(wind_speed)
        support = np.linspace(0.0, 1.2, 200)
        densities = np.zeros((len(mean_prediction), support.size), dtype=float)

        for idx, regime in enumerate(regimes):
            kde = self.kde_models.get(regime)
            if kde is None:
                kde = next(iter(self.kde_models.values()))
            centered_support = support - mean_prediction[idx]
            log_density = kde.score_samples(centered_support[:, None])
            density = np.exp(log_density)
            density /= np.trapz(density, support)
            densities[idx] = density
        cdf = np.cumsum(densities, axis=1)
        cdf /= cdf[:, [-1]]
        return {
            "support": support,
            "density": densities,
            "cdf": cdf,
            "mean": mean_prediction,
        }

    def _encode_features(
        self,
        features: Mapping[str, ArrayLike],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Encode sequential features into fixed-width latent vectors."""
        if self.encoder is None:
            raise RuntimeError("The predictor must be trained before inference.")
        sequences = _as_float_array(features["features"], ndim=3)
        wind_speed = _as_float_array(features["wind_speed"], ndim=1)
        latent = np.vstack([self.encoder.forward(sequence) for sequence in sequences])
        return latent, wind_speed

    def _classify_regime(self, wind_speed: NDArray[np.float64]) -> NDArray[np.str_]:
        """Partition wind speed into low, medium, and high regimes."""
        regimes = np.full(wind_speed.shape, "medium", dtype="<U6")
        regimes[wind_speed < self.cut_in] = "low"
        regimes[wind_speed >= self.rated] = "high"
        return regimes

    def _power_curve_with_hysteresis(self, wind_speed: float, previous_speed: float) -> float:
        """Approximate turbine power curve with ramp-up and ramp-down hysteresis."""
        effective_cut_out = self.cut_out - 1.5 if wind_speed < previous_speed else self.cut_out
        if wind_speed < self.cut_in or wind_speed >= effective_cut_out:
            return 0.0
        if wind_speed <= self.rated:
            normalized = (wind_speed - self.cut_in) / max(self.rated - self.cut_in, 1e-6)
            return normalized ** 3
        return 1.0


class WindScenarioGenerator:
    """Generate spatio-temporal wind power scenarios."""

    def __init__(self, ar_coefficient: float = 0.85, random_state: int = 42) -> None:
        self.ar_coefficient = float(ar_coefficient)
        self.rng = np.random.default_rng(random_state)

    def generate_scenarios(
        self,
        forecast: Mapping[str, ArrayLike],
        correlation_matrix: ArrayLike,
        n_scenarios: int = 100,
    ) -> NDArray[np.float64]:
        """Generate Gaussian-copula scenarios with AR(1) temporal smoothing."""
        mean = _as_float_array(forecast["mean"])
        std = np.maximum(_as_float_array(forecast.get("std", np.full_like(mean, 0.1))), 1e-6)
        correlation = _as_float_array(correlation_matrix, ndim=2)
        if correlation.shape != (mean.size, mean.size):
            raise ValueError("correlation_matrix must have shape (n_sites, n_sites).")

        innovations = self.rng.multivariate_normal(
            mean=np.zeros(mean.size),
            cov=correlation,
            size=n_scenarios,
        )
        horizon = int(np.asarray(forecast.get("horizon", 24)).reshape(()))
        scenarios = np.zeros((n_scenarios, horizon, mean.size), dtype=float)
        state = innovations
        for step in range(horizon):
            noise = self.rng.multivariate_normal(
                mean=np.zeros(mean.size),
                cov=correlation,
                size=n_scenarios,
            )
            state = self.ar_coefficient * state + np.sqrt(1.0 - self.ar_coefficient ** 2) * noise
            uniforms = stats.norm.cdf(state)
            scenarios[:, step, :] = np.clip(
                mean + std * stats.norm.ppf(np.clip(uniforms, 1e-5, 1.0 - 1e-5)),
                0.0,
                1.2,
            )
        return scenarios
