"""Frequency response and frequency security models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import signal
from scipy.integrate import solve_ivp


@dataclass(slots=True)
class GovernorModel:
    """Dynamic governor model for steam, gas, or hydro turbines."""

    governor_type: str = "steam"
    droop: float = 0.05
    Tg: float = 0.2
    Tt: float = 0.5
    Tw: float = 1.5
    reheat_factor: float = 0.3
    nominal_frequency: float = 50.0

    def response(self, delta_f: float, t: float) -> float:
        """Return per-unit governor power response to a frequency deviation."""

        steady_state = (-delta_f / max(self.nominal_frequency * self.droop, 1e-9))
        servo = 1.0 - np.exp(-max(t, 0.0) / max(self.Tg, 1e-9))

        if self.governor_type == "steam":
            turbine = 1.0 - np.exp(-max(t, 0.0) / max(self.Tt, 1e-9))
            reheat = self.reheat_factor + (1.0 - self.reheat_factor) * turbine
            return float(steady_state * servo * reheat)
        if self.governor_type == "gas":
            turbine = 1.0 - np.exp(-max(t, 0.0) / max(self.Tt, 1e-9))
            return float(steady_state * servo * turbine)
        if self.governor_type == "hydro":
            water_column = 1.0 - np.exp(-max(t, 0.0) / max(self.Tw, 1e-9))
            water_hammer = 1.0 - 0.5 * np.exp(-max(t, 0.0) / max(self.Tw, 1e-9))
            return float(steady_state * servo * water_column * water_hammer)
        raise ValueError(f"Unsupported governor type: {self.governor_type}")

    def get_transfer_function(self) -> signal.TransferFunction:
        """Return a continuous-time transfer function representation."""

        gain = 1.0 / max(self.droop, 1e-9)
        if self.governor_type == "steam":
            numerator = gain * np.array([self.reheat_factor * self.Tt, 1.0])
            denominator = np.convolve([self.Tg, 1.0], [self.Tt, 1.0])
        elif self.governor_type == "gas":
            numerator = np.array([gain])
            denominator = np.convolve([self.Tg, 1.0], [self.Tt, 1.0])
        elif self.governor_type == "hydro":
            numerator = gain * np.array([-0.5 * self.Tw, 1.0])
            denominator = np.convolve([self.Tg, 1.0], [self.Tw, 1.0])
        else:
            raise ValueError(f"Unsupported governor type: {self.governor_type}")
        return signal.TransferFunction(numerator, denominator)


@dataclass(slots=True)
class InverterModel:
    """Grid-following or grid-forming inverter frequency support model."""

    inverter_type: str = "grid_following"
    rated_power: float = 100.0
    droop: float = 0.05
    virtual_inertia: float = 0.0
    pll_bandwidth: float = 15.0
    fast_frequency_gain: float = 0.0
    voltage_support_gain: float = 0.0
    nominal_frequency: float = 50.0

    def compute_power_injection(self, frequency: float, voltage: float) -> float:
        """Compute active power injection as a function of frequency and voltage."""

        delta_f = frequency - self.nominal_frequency
        voltage_error = 1.0 - voltage
        if self.inverter_type == "grid_forming":
            droop_power = -self.rated_power * delta_f / max(self.nominal_frequency * self.droop, 1e-9)
            voltage_support = self.voltage_support_gain * voltage_error * self.rated_power
            return float(np.clip(droop_power + voltage_support, -self.rated_power, self.rated_power))
        if self.inverter_type == "grid_following":
            pll_factor = self.pll_bandwidth / (self.pll_bandwidth + abs(delta_f) + 1e-9)
            droop_power = -0.5 * self.rated_power * pll_factor * delta_f / max(self.nominal_frequency * self.droop, 1e-9)
            voltage_support = 0.5 * self.voltage_support_gain * voltage_error * self.rated_power
            return float(np.clip(droop_power + voltage_support, -self.rated_power, self.rated_power))
        raise ValueError(f"Unsupported inverter type: {self.inverter_type}")

    def virtual_inertia_response(self, rocof: float) -> float:
        """Compute inertial and fast frequency response power injection."""

        inertial_power = -self.virtual_inertia * self.rated_power * rocof / max(self.nominal_frequency, 1e-9)
        fast_response = self.fast_frequency_gain * max(-rocof, 0.0)
        return float(np.clip(inertial_power + fast_response, -self.rated_power, self.rated_power))


@dataclass(slots=True)
class SystemFrequencyModel:
    """Aggregate system frequency response model."""

    nominal_frequency: float = 50.0
    load_damping: float = 1.0
    base_power: float = 1000.0
    governors: list[GovernorModel] = field(default_factory=list)
    inverters: list[InverterModel] = field(default_factory=list)
    system_inertia: float = 5.0

    def compute_inertia(self, generators: Sequence[Any]) -> float:
        """Estimate aggregate inertia from online synchronous generators and IBRs."""

        weighted_inertia = 0.0
        online_power = 0.0
        for generator in generators:
            is_online = bool(_get_value(generator, "online", _get_value(generator, "is_online", True)))
            if not is_online:
                continue
            rating = float(_get_value(generator, "rating", _get_value(generator, "mva", _get_value(generator, "power", 0.0))))
            inertia = float(_get_value(generator, "inertia_constant", _get_value(generator, "H", 0.0)))
            weighted_inertia += inertia * rating
            online_power += rating

        virtual_inertia = sum(inverter.virtual_inertia * inverter.rated_power for inverter in self.inverters)
        equivalent_base = max(online_power + sum(inverter.rated_power for inverter in self.inverters), 1e-9)
        self.system_inertia = (weighted_inertia + virtual_inertia) / equivalent_base
        return self.system_inertia

    def simulate_frequency_event(
        self,
        delta_p: float,
        duration: float = 30.0,
        voltage: float = 1.0,
    ) -> dict[str, np.ndarray | float]:
        """Simulate system frequency evolution following a power imbalance."""

        equivalent_inertia = max(
            self.system_inertia + sum(inverter.virtual_inertia for inverter in self.inverters),
            1e-6,
        )
        initial_rocof = self.compute_rocof(delta_p)

        def dynamics(t: float, state: np.ndarray) -> np.ndarray:
            delta_f = float(state[0])
            governor_support = sum(governor.response(delta_f, t) for governor in self.governors)
            inverter_droop = sum(
                inverter.compute_power_injection(self.nominal_frequency + delta_f, voltage)
                for inverter in self.inverters
            ) / max(self.base_power, 1e-9)
            inertial_support = sum(
                inverter.virtual_inertia_response(initial_rocof) for inverter in self.inverters
            ) / max(self.base_power, 1e-9)
            damping = self.load_damping * delta_f / max(self.nominal_frequency, 1e-9)
            mismatch = governor_support + inverter_droop + inertial_support - damping - (delta_p / max(self.base_power, 1e-9))
            derivative = self.nominal_frequency * mismatch / (2.0 * equivalent_inertia)
            return np.array([derivative], dtype=float)

        time = np.linspace(0.0, duration, int(duration / 0.01) + 1)
        solution = solve_ivp(dynamics, (0.0, duration), np.array([0.0]), t_eval=time, max_step=0.05)
        frequency = self.nominal_frequency + solution.y[0]
        rocof = np.gradient(frequency, solution.t, edge_order=2)
        nadir = float(np.min(frequency))

        return {
            "time": solution.t,
            "frequency": frequency,
            "rocof": rocof,
            "nadir": nadir,
            "delta_p": float(delta_p),
        }

    def compute_rocof(self, delta_p: float) -> float:
        """Compute initial rate of change of frequency."""

        equivalent_inertia = max(
            self.system_inertia + sum(inverter.virtual_inertia for inverter in self.inverters),
            1e-6,
        )
        return float(
            -self.nominal_frequency * delta_p / (2.0 * equivalent_inertia * max(self.base_power, 1e-9))
        )

    def compute_nadir(self, delta_p: float) -> float:
        """Estimate frequency nadir from a dynamic simulation."""

        trajectory = self.simulate_frequency_event(delta_p)
        return float(trajectory["nadir"])


@dataclass(slots=True)
class FrequencySecurityAssessment:
    """Frequency security screening and UFLS assessment."""

    frequency_model: SystemFrequencyModel
    nadir_limit: float = 49.0
    rocof_limit: float = 1.0
    ufls_stages: list[tuple[float, float]] = field(
        default_factory=lambda: [(49.0, 0.05), (48.8, 0.10), (48.5, 0.15)]
    )
    recovery_gain: float = 1.0

    def screen_contingencies(
        self,
        network: Any,
        contingency_list: Sequence[Mapping[str, Any]],
    ) -> list[dict[str, Any]]:
        """Screen N-1 contingencies for frequency security."""

        generators = _get_value(network, "generators", [])
        if generators:
            self.frequency_model.compute_inertia(generators)

        results: list[dict[str, Any]] = []
        for contingency in contingency_list:
            delta_p = float(contingency.get("delta_p", contingency.get("power_loss", 0.0)))
            trajectory = self.frequency_model.simulate_frequency_event(delta_p)
            ufls = self.simulate_ufls(trajectory)
            adjusted_frequency = np.asarray(ufls["adjusted_frequency"], dtype=float)
            rocof = np.asarray(trajectory["rocof"], dtype=float)
            results.append(
                {
                    "contingency": contingency.get("name", "unnamed"),
                    "delta_p": delta_p,
                    "nadir": float(np.min(adjusted_frequency)),
                    "rocof": float(np.max(np.abs(rocof))),
                    "ufls_activated": ufls["activated_stages"],
                    "is_secure": bool(
                        np.min(adjusted_frequency) >= self.nadir_limit
                        and np.max(np.abs(rocof)) <= self.rocof_limit
                    ),
                }
            )
        return results

    def compute_max_loss(self, network: Any, tolerance: float = 1.0) -> float:
        """Compute the maximum generation loss that satisfies nadir security."""

        generators = _get_value(network, "generators", [])
        if generators:
            self.frequency_model.compute_inertia(generators)
        upper = float(sum(_get_value(generator, "power", _get_value(generator, "rating", 0.0)) for generator in generators) or self.frequency_model.base_power)
        lower = 0.0
        for _ in range(30):
            midpoint = 0.5 * (lower + upper)
            nadir = self.frequency_model.compute_nadir(midpoint)
            if nadir >= self.nadir_limit:
                lower = midpoint
            else:
                upper = midpoint
            if upper - lower <= tolerance:
                break
        return lower

    def simulate_ufls(
        self,
        frequency_trajectory: Mapping[str, Any] | Sequence[float],
    ) -> dict[str, Any]:
        """Simulate staged under-frequency load shedding."""

        if isinstance(frequency_trajectory, Mapping):
            frequency = np.asarray(frequency_trajectory["frequency"], dtype=float).copy()
            time = np.asarray(
                frequency_trajectory.get("time", np.arange(frequency.size, dtype=float)),
                dtype=float,
            )
        else:
            frequency = np.asarray(frequency_trajectory, dtype=float).copy()
            time = np.arange(frequency.size, dtype=float)

        activated: list[dict[str, float]] = []
        total_shed = 0.0
        adjusted = frequency.copy()
        for threshold, fraction in self.ufls_stages:
            below = np.where(adjusted <= threshold)[0]
            if below.size == 0:
                continue
            index = int(below[0])
            total_shed += fraction
            recovery = self.recovery_gain * fraction * self.frequency_model.nominal_frequency
            adjusted[index:] = adjusted[index:] + recovery
            activated.append(
                {
                    "threshold": float(threshold),
                    "fraction": float(fraction),
                    "time": float(time[index]),
                }
            )

        return {
            "time": time,
            "adjusted_frequency": adjusted,
            "activated_stages": activated,
            "total_shed": float(total_shed),
        }


def _get_value(container: Any, key: str, default: Any = None) -> Any:
    if isinstance(container, Mapping):
        return container.get(key, default)
    return getattr(container, key, default)
