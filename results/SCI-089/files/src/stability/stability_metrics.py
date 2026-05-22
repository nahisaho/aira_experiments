"""Grid stability metrics and reporting helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np


@dataclass(slots=True)
class StabilityMetrics:
    """Collection of transient, frequency, and voltage stability metrics."""

    nominal_frequency: float = 50.0

    def compute_tsi(self, angles: Mapping[str, Sequence[float]] | Sequence[Sequence[float]]) -> float:
        """Compute the transient stability index from machine angle trajectories."""

        angle_matrix = _as_matrix(angles)
        angle_matrix = np.unwrap(angle_matrix, axis=1)
        max_spread = np.max(angle_matrix.max(axis=0) - angle_matrix.min(axis=0))
        max_spread_deg = np.degrees(max_spread)
        tsi = 100.0 * (360.0 - max_spread_deg) / (360.0 + max_spread_deg)
        return float(np.clip(tsi, 0.0, 100.0))

    def compute_fsi(self, frequency: Mapping[str, Sequence[float]] | Sequence[float]) -> float:
        """Compute a normalized frequency stability index."""

        if isinstance(frequency, Mapping):
            values = np.asarray(frequency["frequency"], dtype=float)
            time = np.asarray(
                frequency.get("time", np.arange(values.size, dtype=float)),
                dtype=float,
            )
        else:
            values = np.asarray(frequency, dtype=float)
            time = np.arange(values.size, dtype=float)

        max_deviation = np.max(np.abs(values - self.nominal_frequency))
        rocof = np.max(np.abs(np.gradient(values, time, edge_order=1))) if values.size > 1 else 0.0
        deviation_score = max(0.0, 1.0 - max_deviation / 2.0)
        rocof_score = max(0.0, 1.0 - rocof / 2.0)
        return float(100.0 * (0.7 * deviation_score + 0.3 * rocof_score))

    def compute_voltage_margin(self, network: Any, bus: str | int) -> dict[str, float | str | int]:
        """Compute voltage stability margin using P-V curve data when available."""

        pv_curves = _get_value(network, "pv_curves", {})
        curve = pv_curves.get(bus) if isinstance(pv_curves, Mapping) else None
        if curve is not None:
            powers = np.asarray(curve["power"], dtype=float)
            collapse_power = float(np.max(powers))
        else:
            collapse_power = float(_nested_lookup(network, "collapse_power", bus, 0.0))

        current_power = float(
            _nested_lookup(network, "operating_power", bus, _nested_lookup(network, "bus_power", bus, 0.0))
        )
        margin = collapse_power - current_power
        margin_ratio = margin / max(collapse_power, 1e-9)
        return {
            "bus": bus,
            "current_power": current_power,
            "collapse_power": collapse_power,
            "margin": margin,
            "margin_ratio": margin_ratio,
        }

    def compute_scr(self, network: Any, bus: str | int) -> float:
        """Compute short-circuit ratio at the point of common coupling."""

        short_circuit_mva = float(_nested_lookup(network, "short_circuit_mva", bus, 0.0))
        renewable_rating = float(_nested_lookup(network, "renewable_rating_mva", bus, _nested_lookup(network, "renewable_rating_mw", bus, 0.0)))
        if renewable_rating <= 0.0:
            raise ValueError(f"Renewable plant rating is not available for bus {bus!r}.")
        return short_circuit_mva / renewable_rating


@dataclass(slots=True)
class StabilityReport:
    """Traffic-light stability reporting and recommendation engine."""

    metrics: StabilityMetrics = field(default_factory=StabilityMetrics)
    _last_report: dict[str, Any] = field(default_factory=dict, init=False)

    def generate_report(self, network: Any, results: Mapping[str, Any]) -> dict[str, Any]:
        """Generate a comprehensive stability assessment report."""

        metrics_report: dict[str, Any] = {}
        transient = results.get("transient")
        if transient is not None:
            metrics_report["transient_stability_index"] = self.metrics.compute_tsi(transient.angles)

        frequency = results.get("frequency")
        if frequency is not None:
            metrics_report["frequency_stability_index"] = self.metrics.compute_fsi(frequency)

        monitored_buses = results.get("monitored_buses", _get_value(network, "monitored_buses", []))
        if monitored_buses:
            metrics_report["voltage_margin"] = {
                str(bus): self.metrics.compute_voltage_margin(network, bus) for bus in monitored_buses
            }
            metrics_report["short_circuit_ratio"] = {
                str(bus): self.metrics.compute_scr(network, bus) for bus in monitored_buses
            }

        traffic_lights = self._build_traffic_lights(metrics_report)
        recommendations = self._build_recommendations(metrics_report, traffic_lights)
        self._last_report = {
            "metrics": metrics_report,
            "traffic_lights": traffic_lights,
            "recommendations": recommendations,
        }
        return self._last_report

    def get_recommendations(self) -> list[str]:
        """Return recommendations from the last generated report."""

        return list(self._last_report.get("recommendations", []))

    def _build_traffic_lights(self, metrics_report: Mapping[str, Any]) -> dict[str, Any]:
        traffic_lights: dict[str, Any] = {}
        for key, value in metrics_report.items():
            if key in {"transient_stability_index", "frequency_stability_index"}:
                traffic_lights[key] = _score_to_colour(float(value), green=80.0, yellow=60.0)
            elif key == "voltage_margin":
                traffic_lights[key] = {
                    bus: _score_to_colour(float(entry["margin_ratio"]), green=0.2, yellow=0.1)
                    for bus, entry in value.items()
                }
            elif key == "short_circuit_ratio":
                traffic_lights[key] = {
                    bus: _score_to_colour(float(entry), green=3.0, yellow=2.0)
                    for bus, entry in value.items()
                }
        return traffic_lights

    def _build_recommendations(
        self,
        metrics_report: Mapping[str, Any],
        traffic_lights: Mapping[str, Any],
    ) -> list[str]:
        recommendations: list[str] = []
        if traffic_lights.get("transient_stability_index") in {"yellow", "red"}:
            recommendations.append("Increase damping support or reduce fault clearing time.")
        if traffic_lights.get("frequency_stability_index") in {"yellow", "red"}:
            recommendations.append("Add primary reserve, fast frequency response, or virtual inertia.")

        voltage_flags = traffic_lights.get("voltage_margin", {})
        if isinstance(voltage_flags, Mapping):
            for bus, colour in voltage_flags.items():
                if colour in {"yellow", "red"}:
                    recommendations.append(f"Review reactive support and loading at bus {bus}.")

        scr_flags = traffic_lights.get("short_circuit_ratio", {})
        if isinstance(scr_flags, Mapping):
            for bus, colour in scr_flags.items():
                if colour in {"yellow", "red"}:
                    recommendations.append(f"Strengthen grid connection or reduce inverter penetration at bus {bus}.")

        if not recommendations and metrics_report:
            recommendations.append("No corrective action required under the assessed operating point.")
        return recommendations


def _as_matrix(angles: Mapping[str, Sequence[float]] | Sequence[Sequence[float]]) -> np.ndarray:
    if isinstance(angles, Mapping):
        return np.vstack([np.asarray(series, dtype=float) for series in angles.values()])
    return np.vstack([np.asarray(series, dtype=float) for series in angles])


def _get_value(container: Any, key: str, default: Any = None) -> Any:
    if isinstance(container, Mapping):
        return container.get(key, default)
    return getattr(container, key, default)


def _nested_lookup(container: Any, key: str, index: str | int, default: Any = None) -> Any:
    values = _get_value(container, key)
    if isinstance(values, Mapping):
        return values.get(index, default)
    return default


def _score_to_colour(value: float, green: float, yellow: float) -> str:
    if value >= green:
        return "green"
    if value >= yellow:
        return "yellow"
    return "red"
