"""Probabilistic solar power forecasting utilities.

This module provides a compact implementation of a renewable energy forecasting
stack for solar assets. The implementation focuses on robust numerical
building blocks that can be reused in research prototypes and lightweight
production pipelines without depending on domain-specific packages.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import stats
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.distance import cdist
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.neighbors import NearestNeighbors


def _as_1d_float(values: ArrayLike) -> NDArray[np.float64]:
    """Convert an array-like object to a one-dimensional float array."""
    array = np.asarray(values, dtype=float)
    if array.ndim == 0:
        return array.reshape(1)
    return array.astype(float, copy=False).reshape(-1)


def _as_time_array(values: ArrayLike) -> NDArray[np.datetime64]:
    """Convert inputs to a one-dimensional ``datetime64`` array."""
    array = np.asarray(values, dtype="datetime64[m]")
    if array.ndim == 0:
        return array.reshape(1)
    return array.reshape(-1)


@dataclass
class SolarTrainingBundle:
    """Container for fitted solar forecasting artifacts."""

    feature_matrix: NDArray[np.float64]
    target: NDArray[np.float64]
    residuals: NDArray[np.float64]


class NWPSolarPreprocessor:
    """Preprocess NWP data for solar probabilistic forecasting.

    Parameters
    ----------
    site_latitudes:
        Latitude of target forecast sites in degrees.
    site_longitudes:
        Longitude of target forecast sites in degrees.
    altitude:
        Site altitude in meters above sea level.
    linke_turbidity:
        Effective Linke turbidity used by the simplified Ineichen clear-sky
        approximation.
    """

    def __init__(
        self,
        site_latitudes: ArrayLike,
        site_longitudes: ArrayLike,
        altitude: float | ArrayLike = 0.0,
        linke_turbidity: float = 3.0,
    ) -> None:
        self.site_latitudes = _as_1d_float(site_latitudes)
        self.site_longitudes = _as_1d_float(site_longitudes)
        self.altitude = np.asarray(altitude, dtype=float)
        self.linke_turbidity = float(linke_turbidity)

    def preprocess(self, nwp_data: Mapping[str, ArrayLike]) -> dict[str, NDArray[np.float64]]:
        """Preprocess NWP features and derive clear-sky diagnostics.

        Parameters
        ----------
        nwp_data:
            Dictionary containing at least ``ghi``, ``dni``, ``dhi``,
            ``temperature``, ``cloud_cover`` and ``times``. If the arrays are
            defined on a latitude/longitude grid, ``latitudes`` and
            ``longitudes`` can be supplied and the method interpolates them to
            the configured sites.

        Returns
        -------
        dict
            Processed features keyed by variable name.
        """
        required = ("ghi", "dni", "dhi", "temperature", "cloud_cover", "times")
        missing = [key for key in required if key not in nwp_data]
        if missing:
            raise KeyError(f"Missing required NWP fields: {missing}")

        processed: dict[str, NDArray[np.float64]] = {}
        for key in ("ghi", "dni", "dhi", "temperature", "cloud_cover"):
            processed[key] = np.asarray(nwp_data[key], dtype=float)

        times = _as_time_array(nwp_data["times"])
        solar_position = self._compute_solar_position(times)
        clear_sky = self._ineichen_clear_sky(times, solar_position["cos_zenith"])

        if "latitudes" in nwp_data and "longitudes" in nwp_data:
            interpolated = self.interpolate_to_sites(
                {key: processed[key] for key in processed},
                _as_1d_float(nwp_data["latitudes"]),
                _as_1d_float(nwp_data["longitudes"]),
            )
            processed.update(interpolated)

        processed["clear_sky_ghi"] = clear_sky["ghi"]
        processed["clear_sky_dni"] = clear_sky["dni"]
        processed["clear_sky_dhi"] = clear_sky["dhi"]
        processed["clear_sky_index"] = self.compute_clear_sky_index(
            processed["ghi"], clear_sky["ghi"]
        )
        processed["apparent_zenith"] = solar_position["zenith"]
        processed["solar_elevation"] = 90.0 - solar_position["zenith"]
        processed["times_numeric"] = times.astype("datetime64[m]").astype(float)
        return processed

    def compute_clear_sky_index(
        self,
        ghi: ArrayLike | None = None,
        clear_sky_ghi: ArrayLike | None = None,
    ) -> NDArray[np.float64]:
        """Compute the clear-sky index.

        Parameters
        ----------
        ghi:
            Measured or forecast global horizontal irradiance.
        clear_sky_ghi:
            Clear-sky irradiance computed from the Ineichen approximation.
        """
        if ghi is None or clear_sky_ghi is None:
            raise ValueError("Both ghi and clear_sky_ghi must be provided.")
        ghi_array = np.asarray(ghi, dtype=float)
        clear_array = np.maximum(np.asarray(clear_sky_ghi, dtype=float), 1e-6)
        return np.clip(ghi_array / clear_array, 0.0, 2.0)

    def interpolate_to_sites(
        self,
        fields: Mapping[str, ArrayLike],
        latitudes: ArrayLike,
        longitudes: ArrayLike,
    ) -> dict[str, NDArray[np.float64]]:
        """Interpolate gridded NWP fields to the configured sites.

        Parameters
        ----------
        fields:
            Mapping of variable names to arrays with trailing ``(lat, lon)``
            dimensions.
        latitudes, longitudes:
            Grid coordinates corresponding to the NWP fields.
        """
        latitudes_1d = _as_1d_float(latitudes)
        longitudes_1d = _as_1d_float(longitudes)
        targets = np.column_stack([self.site_latitudes, self.site_longitudes])
        interpolated: dict[str, NDArray[np.float64]] = {}

        for name, values in fields.items():
            array = np.asarray(values, dtype=float)
            if array.ndim < 2 or array.shape[-2:] != (latitudes_1d.size, longitudes_1d.size):
                interpolated[name] = array
                continue

            flattened_shape = (-1, latitudes_1d.size, longitudes_1d.size)
            time_slices = array.reshape(flattened_shape)
            site_values = []
            for slice_2d in time_slices:
                interpolator = RegularGridInterpolator(
                    (latitudes_1d, longitudes_1d),
                    slice_2d,
                    bounds_error=False,
                    fill_value=None,
                )
                site_values.append(interpolator(targets))
            interpolated[name] = np.vstack(site_values)
        return interpolated

    def _compute_solar_position(
        self,
        times: NDArray[np.datetime64],
    ) -> dict[str, NDArray[np.float64]]:
        """Approximate solar position using standard astronomical formulas."""
        timestamps = times.astype("datetime64[m]")
        midnight = timestamps.astype("datetime64[D]")
        minutes = (timestamps - midnight).astype("timedelta64[m]").astype(float)
        doy = (
            midnight.astype("datetime64[D]")
            - midnight.astype("datetime64[Y]")
        ).astype(int) + 1

        hour = minutes / 60.0
        gamma = 2.0 * np.pi * (doy - 1 + (hour - 12.0) / 24.0) / 365.0
        decl = (
            0.006918
            - 0.399912 * np.cos(gamma)
            + 0.070257 * np.sin(gamma)
            - 0.006758 * np.cos(2 * gamma)
            + 0.000907 * np.sin(2 * gamma)
            - 0.002697 * np.cos(3 * gamma)
            + 0.00148 * np.sin(3 * gamma)
        )
        eqtime = 229.18 * (
            0.000075
            + 0.001868 * np.cos(gamma)
            - 0.032077 * np.sin(gamma)
            - 0.014615 * np.cos(2 * gamma)
            - 0.040849 * np.sin(2 * gamma)
        )

        latitude = np.deg2rad(np.mean(self.site_latitudes))
        longitude = np.mean(self.site_longitudes)
        time_offset = eqtime + 4.0 * longitude
        tst = (minutes + time_offset) % 1440.0
        hour_angle = np.deg2rad(tst / 4.0 - 180.0)

        cos_zenith = np.clip(
            np.sin(latitude) * np.sin(decl)
            + np.cos(latitude) * np.cos(decl) * np.cos(hour_angle),
            0.0,
            1.0,
        )
        zenith = np.degrees(np.arccos(np.clip(cos_zenith, -1.0, 1.0)))
        return {"cos_zenith": cos_zenith, "zenith": zenith}

    def _ineichen_clear_sky(
        self,
        times: NDArray[np.datetime64],
        cos_zenith: NDArray[np.float64],
    ) -> dict[str, NDArray[np.float64]]:
        """Simplified Ineichen clear-sky irradiance approximation."""
        doy = (
            times.astype("datetime64[D]") - times.astype("datetime64[Y]")
        ).astype(int) + 1
        eccentricity = 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)
        i0 = 1367.0 * eccentricity
        altitude_km = np.mean(np.asarray(self.altitude, dtype=float)) / 1000.0
        fh1 = np.exp(-altitude_km / 8.0)
        fh2 = np.exp(-altitude_km / 1.25)
        cg1 = 0.868 + 5.09e-5 * np.mean(np.asarray(self.altitude, dtype=float))
        cg2 = 0.0387 + 3.92e-5 * np.mean(np.asarray(self.altitude, dtype=float))
        air_mass = 1.0 / np.maximum(cos_zenith + 0.15 * (93.885 - np.degrees(np.arccos(np.clip(cos_zenith, 0.0, 1.0)))) ** -1.253, 1e-3)
        tl = max(self.linke_turbidity, 1.0)
        dni = cg1 * i0 * cos_zenith * np.exp(-cg2 * air_mass * (fh1 + fh2 * (tl - 1.0)))
        dni = np.maximum(dni, 0.0)
        ghi = np.maximum(dni * np.maximum(cos_zenith, 0.0) + 0.1 * i0 * cos_zenith, 0.0)
        dhi = np.maximum(ghi - dni * cos_zenith, 0.0)
        return {"ghi": ghi, "dni": dni, "dhi": dhi}


class SolarMLPredictor:
    """Machine-learning solar power predictor with quantile support.

    The deterministic forecast is represented by a gradient boosting regressor.
    Additional gradient boosting models trained with quantile loss provide
    probabilistic intervals. An analog ensemble built on nearest neighbours of
    historical feature vectors is used as an extra uncertainty quantification
    mechanism.
    """

    def __init__(
        self,
        quantiles: Sequence[float] = (0.1, 0.25, 0.5, 0.75, 0.9),
        n_analogs: int = 20,
        random_state: int = 42,
    ) -> None:
        self.quantiles = tuple(float(q) for q in quantiles)
        self.n_analogs = int(n_analogs)
        self.random_state = int(random_state)
        self.point_model = GradientBoostingRegressor(random_state=random_state)
        self.quantile_models = {
            q: GradientBoostingRegressor(
                loss="quantile",
                alpha=q,
                random_state=random_state,
            )
            for q in self.quantiles
        }
        self.analog_model = NearestNeighbors(n_neighbors=max(1, self.n_analogs))
        self.training_bundle: SolarTrainingBundle | None = None

    def train(self, historical_data: Mapping[str, ArrayLike]) -> None:
        """Train deterministic, quantile, and analog ensemble models."""
        features, target = self._build_feature_matrix(historical_data, fit_mode=True)
        self.point_model.fit(features, target)
        point_prediction = self.point_model.predict(features)
        residuals = target - point_prediction

        for model in self.quantile_models.values():
            model.fit(features, target)
        self.analog_model = NearestNeighbors(n_neighbors=min(max(1, self.n_analogs), len(target)))
        self.analog_model.fit(features)
        self.training_bundle = SolarTrainingBundle(features, target, residuals)

    def predict(self, features: Mapping[str, ArrayLike]) -> NDArray[np.float64]:
        """Generate a deterministic forecast."""
        if self.training_bundle is None:
            raise RuntimeError("The predictor must be trained before calling predict().")
        feature_matrix, _ = self._build_feature_matrix(features, fit_mode=False)
        return self.point_model.predict(feature_matrix)

    def predict_quantiles(self, features: Mapping[str, ArrayLike]) -> dict[str, NDArray[np.float64]]:
        """Predict quantile forecasts and analog ensemble uncertainty summaries."""
        if self.training_bundle is None:
            raise RuntimeError("The predictor must be trained before calling predict_quantiles().")

        feature_matrix, _ = self._build_feature_matrix(features, fit_mode=False)
        quantile_predictions = {
            f"q{int(q * 100):02d}": model.predict(feature_matrix)
            for q, model in self.quantile_models.items()
        }

        neighbor_idx = self.analog_model.kneighbors(feature_matrix, return_distance=False)
        analog_targets = self.training_bundle.target[neighbor_idx]
        analog_residuals = self.training_bundle.residuals[neighbor_idx]
        quantile_predictions["analog_mean"] = analog_targets.mean(axis=1)
        quantile_predictions["analog_std"] = analog_residuals.std(axis=1)
        return quantile_predictions

    def _build_feature_matrix(
        self,
        data: Mapping[str, ArrayLike],
        fit_mode: bool,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64] | None]:
        """Construct a feature matrix from temporal, NWP, and lagged inputs."""
        times = _as_time_array(data["times"])
        solar_position = self._solar_feature_block(times, data)
        nwp_features = [
            np.asarray(data[key], dtype=float).reshape(len(times), -1)
            for key in ("ghi", "dni", "dhi", "temperature", "cloud_cover", "clear_sky_index")
            if key in data
        ]
        lagged = [
            np.asarray(data[key], dtype=float).reshape(len(times), -1)
            for key in sorted(k for k in data if k.startswith("lag_"))
        ]
        temporal = self._temporal_feature_block(times)
        feature_blocks = [solar_position, temporal, *nwp_features, *lagged]
        feature_matrix = np.hstack(feature_blocks)
        target = None
        if fit_mode:
            if "power" not in data:
                raise KeyError("Historical data must include 'power' for training.")
            target = np.asarray(data["power"], dtype=float).reshape(len(times))
        return feature_matrix, target

    def _solar_feature_block(
        self,
        times: NDArray[np.datetime64],
        data: Mapping[str, ArrayLike],
    ) -> NDArray[np.float64]:
        """Build solar-position-based explanatory variables."""
        latitude = float(np.mean(np.asarray(data.get("latitude", 35.0), dtype=float)))
        longitude = float(np.mean(np.asarray(data.get("longitude", 135.0), dtype=float)))
        preprocessor = NWPSolarPreprocessor([latitude], [longitude])
        solar_position = preprocessor._compute_solar_position(times)
        zenith_rad = np.deg2rad(solar_position["zenith"])
        azimuth_proxy = np.sin(2.0 * np.pi * times.astype("datetime64[m]").astype(float) / 1440.0)
        return np.column_stack(
            [
                np.cos(zenith_rad),
                np.sin(zenith_rad),
                azimuth_proxy,
            ]
        )

    def _temporal_feature_block(self, times: NDArray[np.datetime64]) -> NDArray[np.float64]:
        """Encode calendar and diurnal cycles as cyclical features."""
        minutes = (times - times.astype("datetime64[D]")).astype("timedelta64[m]").astype(float)
        day_fraction = minutes / 1440.0
        doy = (times.astype("datetime64[D]") - times.astype("datetime64[Y]")).astype(int) + 1
        return np.column_stack(
            [
                np.sin(2.0 * np.pi * day_fraction),
                np.cos(2.0 * np.pi * day_fraction),
                np.sin(2.0 * np.pi * doy / 365.25),
                np.cos(2.0 * np.pi * doy / 365.25),
            ]
        )


class SolarScenarioGenerator:
    """Generate correlated solar power scenarios from quantile forecasts."""

    def __init__(self, random_state: int = 42) -> None:
        self.random_state = np.random.default_rng(random_state)
        self._last_scenarios: NDArray[np.float64] | None = None

    def generate_scenarios(
        self,
        forecast: Mapping[str, ArrayLike],
        n_scenarios: int = 100,
    ) -> NDArray[np.float64]:
        """Generate correlated scenarios using a Gaussian copula approximation."""
        quantile_levels, quantile_values = self._extract_quantiles(forecast)
        median = self._interp_quantile(0.5, quantile_levels, quantile_values)
        spread = np.maximum(
            self._interp_quantile(0.9, quantile_levels, quantile_values)
            - self._interp_quantile(0.1, quantile_levels, quantile_values),
            1e-6,
        ) / 2.563
        correlation = self._estimate_correlation(median)
        latent = self.random_state.multivariate_normal(
            mean=np.zeros(median.shape[0]),
            cov=correlation,
            size=n_scenarios,
        )
        uniforms = stats.norm.cdf(latent)
        scenarios = np.empty_like(uniforms)
        for idx in range(median.shape[0]):
            scenarios[:, idx] = np.interp(
                uniforms[:, idx],
                quantile_levels,
                quantile_values[idx],
                left=quantile_values[idx, 0],
                right=quantile_values[idx, -1],
            )
        self._last_scenarios = scenarios
        return scenarios

    def reduce_scenarios(self, n_reduced: int = 10) -> dict[str, NDArray[np.float64]]:
        """Reduce the scenario set using a simple k-medoids procedure."""
        if self._last_scenarios is None:
            raise RuntimeError("Call generate_scenarios() before reduce_scenarios().")
        scenarios = self._last_scenarios
        n_reduced = min(max(1, int(n_reduced)), len(scenarios))
        medoids = self._initialize_medoids(scenarios, n_reduced)
        distances = cdist(scenarios, scenarios[medoids])
        labels = distances.argmin(axis=1)
        reduced = []
        weights = []
        for cluster in range(n_reduced):
            members = np.where(labels == cluster)[0]
            if members.size == 0:
                continue
            cluster_distances = cdist(scenarios[members], scenarios[members])
            medoid_index = members[np.argmin(cluster_distances.sum(axis=1))]
            reduced.append(scenarios[medoid_index])
            weights.append(members.size / len(scenarios))
        return {
            "scenarios": np.asarray(reduced, dtype=float),
            "weights": np.asarray(weights, dtype=float),
        }

    def _extract_quantiles(
        self,
        forecast: Mapping[str, ArrayLike],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Extract sorted quantile levels and values from a forecast mapping."""
        quantile_items = []
        for key, value in forecast.items():
            if not key.startswith("q"):
                continue
            quantile_items.append((float(key[1:]) / 100.0, np.asarray(value, dtype=float).reshape(-1)))
        if not quantile_items:
            raise KeyError("Forecast must contain quantile keys such as 'q10' or 'q90'.")
        quantile_items.sort(key=lambda item: item[0])
        levels = np.asarray([item[0] for item in quantile_items], dtype=float)
        values = np.vstack([item[1] for item in quantile_items]).T
        return levels, values

    def _interp_quantile(
        self,
        level: float,
        quantile_levels: NDArray[np.float64],
        quantile_values: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Interpolate a target quantile level across lead times."""
        return np.asarray(
            [
                np.interp(level, quantile_levels, quantile_values_row)
                for quantile_values_row in quantile_values
            ],
            dtype=float,
        )

    def _estimate_correlation(self, median: NDArray[np.float64]) -> NDArray[np.float64]:
        """Construct a smooth correlation matrix over forecast horizons."""
        n_steps = median.shape[0]
        distances = np.abs(np.subtract.outer(np.arange(n_steps), np.arange(n_steps)))
        scale = max(n_steps / 6.0, 1.0)
        correlation = np.exp(-distances / scale)
        correlation += np.eye(n_steps) * 1e-6
        return correlation

    def _initialize_medoids(self, scenarios: NDArray[np.float64], n_medoids: int) -> NDArray[np.int64]:
        """Initialize medoids by sampling diverse trajectories."""
        if n_medoids == 1:
            return np.array([0], dtype=int)
        norms = np.linalg.norm(scenarios - scenarios.mean(axis=0), axis=1)
        order = np.argsort(norms)
        indices = np.linspace(0, len(order) - 1, n_medoids, dtype=int)
        return order[indices]
