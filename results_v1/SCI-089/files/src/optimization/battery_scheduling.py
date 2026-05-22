"""Battery energy storage system scheduling models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog


ArrayLike = Sequence[float] | np.ndarray


@dataclass(slots=True)
class BatteryParameters:
    """Core battery parameters for scheduling and degradation tracking."""

    capacity_kwh: float
    max_power_kw: float
    efficiency: float
    min_soc: float
    max_soc: float
    degradation_rate: float
    replacement_cost_per_kwh: float = 150.0
    temperature_c: float = 25.0
    calendar_aging_rate: float = 0.02
    cycle_life_at_full_dod: float = 4_000.0


@dataclass(slots=True)
class BatteryState:
    """Mutable battery state values."""

    soc: float
    throughput_kwh: float = 0.0
    equivalent_cycles: float = 0.0
    degradation_fraction: float = 0.0
    temperature_c: float = 25.0


@dataclass(slots=True)
class BatteryDispatchResult:
    """Optimal battery dispatch result for a single unit."""

    charge_kw: np.ndarray
    discharge_kw: np.ndarray
    soc: np.ndarray
    grid_import_kw: np.ndarray
    grid_export_kw: np.ndarray
    regulation_kw: np.ndarray
    objective_value: float
    expected_profit: float
    degradation_cost: float
    status: str
    message: str


@dataclass(slots=True)
class FleetDispatchResult:
    """Aggregated fleet scheduling result."""

    per_battery: dict[str, BatteryDispatchResult]
    aggregate_charge_kw: np.ndarray
    aggregate_discharge_kw: np.ndarray
    aggregate_soc_kwh: np.ndarray
    fairness_index: float
    status: str
    message: str


class BatteryModel:
    """Lithium-ion battery model with efficiency and degradation effects."""

    def __init__(
        self,
        capacity_kwh: float,
        max_power_kw: float,
        efficiency: float,
        min_soc: float,
        max_soc: float,
        degradation_rate: float,
        initial_soc: Optional[float] = None,
        temperature_c: float = 25.0,
    ) -> None:
        self.parameters = BatteryParameters(
            capacity_kwh=capacity_kwh,
            max_power_kw=max_power_kw,
            efficiency=efficiency,
            min_soc=min_soc,
            max_soc=max_soc,
            degradation_rate=degradation_rate,
            temperature_c=temperature_c,
        )
        start_soc = initial_soc if initial_soc is not None else 0.5 * (min_soc + max_soc)
        self.state = BatteryState(
            soc=float(np.clip(start_soc, min_soc, max_soc)),
            temperature_c=temperature_c,
        )

    def charge(self, power: float, duration: float) -> float:
        """Charge the battery and return stored energy in kWh."""

        power = float(np.clip(power, 0.0, self.parameters.max_power_kw))
        duration = max(float(duration), 0.0)
        efficiency = self._effective_efficiency(power)
        room_kwh = (self.parameters.max_soc - self.state.soc) * self.parameters.capacity_kwh
        input_energy = min(power * duration, room_kwh / max(efficiency, 1e-9))
        stored_energy = input_energy * efficiency
        self.state.soc += stored_energy / self.parameters.capacity_kwh
        self._apply_aging(abs(input_energy), duration, stored_energy / max(self.parameters.capacity_kwh, 1e-9))
        return stored_energy

    def discharge(self, power: float, duration: float) -> float:
        """Discharge the battery and return delivered energy in kWh."""

        power = float(np.clip(power, 0.0, self.parameters.max_power_kw))
        duration = max(float(duration), 0.0)
        efficiency = self._effective_efficiency(power)
        available_kwh = (self.state.soc - self.parameters.min_soc) * self.parameters.capacity_kwh
        internal_energy = min(power * duration / max(efficiency, 1e-9), available_kwh)
        delivered_energy = internal_energy * efficiency
        self.state.soc -= internal_energy / self.parameters.capacity_kwh
        self._apply_aging(abs(internal_energy), duration, internal_energy / max(self.parameters.capacity_kwh, 1e-9))
        return delivered_energy

    def get_soc(self) -> float:
        """Return the current state of charge in per-unit."""

        return float(self.state.soc)

    def get_degradation_cost(self) -> float:
        """Return the estimated degradation cost in monetary units."""

        degraded_capacity = self.state.degradation_fraction * self.parameters.capacity_kwh
        return degraded_capacity * self.parameters.replacement_cost_per_kwh

    def _effective_efficiency(self, power: float) -> float:
        temperature_factor = max(0.85, 1.0 - 0.002 * abs(self.state.temperature_c - 25.0))
        power_ratio = power / max(self.parameters.max_power_kw, 1e-9)
        curve_factor = max(0.9, 1.0 - 0.05 * power_ratio**2)
        return float(np.clip(self.parameters.efficiency * temperature_factor * curve_factor, 0.7, 0.999))

    def _apply_aging(self, throughput_kwh: float, duration_h: float, dod: float) -> None:
        self.state.throughput_kwh += throughput_kwh
        cycle_stress = 1.0 + max(dod - 0.8, 0.0) * 1.5
        equivalent_cycles = throughput_kwh / max(2.0 * self.parameters.capacity_kwh, 1e-9)
        equivalent_cycles *= cycle_stress
        self.state.equivalent_cycles += equivalent_cycles
        cycle_degradation = equivalent_cycles * self.parameters.degradation_rate
        calendar_degradation = (
            self.parameters.calendar_aging_rate
            * self.parameters.degradation_rate
            * duration_h
            / (24.0 * 365.0)
        )
        self.state.degradation_fraction = float(
            np.clip(self.state.degradation_fraction + cycle_degradation + calendar_degradation, 0.0, 1.0)
        )


class BatteryScheduler:
    """LP-based rolling-horizon scheduler for battery dispatch.

    The formulation is a convex relaxation of mixed-integer battery scheduling,
    keeping the API compatible with downstream MILP-style workflows.
    """

    def __init__(
        self,
        battery: BatteryModel,
        interval_minutes: int = 15,
        peak_weight: float = 5.0,
        degradation_weight: float = 1.0,
    ) -> None:
        self.battery = battery
        self.interval_minutes = interval_minutes
        self.peak_weight = peak_weight
        self.degradation_weight = degradation_weight
        self._rolling_inputs: dict[str, Any] = {}
        self._last_result: BatteryDispatchResult | None = None

    def optimize_schedule(
        self,
        prices: ArrayLike | Mapping[str, ArrayLike],
        load: ArrayLike,
        renewable: ArrayLike,
        horizon_hours: int = 24,
    ) -> BatteryDispatchResult:
        """Optimize arbitrage, peak shaving, and regulation dispatch."""

        dt = self.interval_minutes / 60.0
        load_array = np.asarray(load, dtype=float)
        renewable_array = np.asarray(renewable, dtype=float)
        buy_price, sell_price, regulation_price = self._parse_prices(prices)
        horizon_steps = min(
            int(round(horizon_hours / dt)),
            load_array.size,
            renewable_array.size,
            buy_price.size,
            sell_price.size,
            regulation_price.size,
        )
        load_array = load_array[:horizon_steps]
        renewable_array = renewable_array[:horizon_steps]
        buy_price = buy_price[:horizon_steps]
        sell_price = sell_price[:horizon_steps]
        regulation_price = regulation_price[:horizon_steps]

        params = self.battery.parameters
        efficiency = self.battery._effective_efficiency(0.5 * params.max_power_kw)
        n_t = horizon_steps
        n_vars = 6 * n_t + 2
        grid_exchange_limit = float(
            max(
                np.max(np.abs(load_array - renewable_array)) + params.max_power_kw,
                params.max_power_kw,
            )
        )

        charge_slice = slice(0, n_t)
        discharge_slice = slice(n_t, 2 * n_t)
        import_slice = slice(2 * n_t, 3 * n_t)
        export_slice = slice(3 * n_t, 4 * n_t)
        regulation_slice = slice(4 * n_t, 5 * n_t)
        soc_slice = slice(5 * n_t, 6 * n_t + 1)
        peak_index = 6 * n_t + 1

        c = np.zeros(n_vars)
        c[import_slice] = buy_price * dt
        c[export_slice] = -sell_price * dt
        c[regulation_slice] = -regulation_price * dt
        throughput_cost = (
            params.degradation_rate
            * params.replacement_cost_per_kwh
            * self.degradation_weight
            / max(params.capacity_kwh, 1e-9)
        )
        c[charge_slice] = throughput_cost * dt
        c[discharge_slice] = throughput_cost * dt
        c[peak_index] = self.peak_weight

        bounds: list[tuple[float | None, float | None]] = []
        bounds.extend([(0.0, params.max_power_kw)] * n_t)
        bounds.extend([(0.0, params.max_power_kw)] * n_t)
        bounds.extend([(0.0, grid_exchange_limit)] * n_t)
        bounds.extend([(0.0, grid_exchange_limit)] * n_t)
        bounds.extend([(0.0, params.max_power_kw)] * n_t)
        bounds.extend(
            [
                (params.min_soc * params.capacity_kwh, params.max_soc * params.capacity_kwh)
                for _ in range(n_t + 1)
            ]
        )
        bounds.append((0.0, None))

        a_eq: list[np.ndarray] = []
        b_eq: list[float] = []
        a_ub: list[np.ndarray] = []
        b_ub: list[float] = []

        initial_soc_kwh = self.battery.state.soc * params.capacity_kwh
        row = np.zeros(n_vars)
        row[soc_slice.start] = 1.0
        a_eq.append(row)
        b_eq.append(initial_soc_kwh)

        for t_idx in range(n_t):
            soc_row = np.zeros(n_vars)
            soc_row[soc_slice.start + t_idx + 1] = 1.0
            soc_row[soc_slice.start + t_idx] = -1.0
            soc_row[charge_slice.start + t_idx] = -efficiency * dt
            soc_row[discharge_slice.start + t_idx] = dt / max(efficiency, 1e-9)
            a_eq.append(soc_row)
            b_eq.append(0.0)

            balance_row = np.zeros(n_vars)
            balance_row[import_slice.start + t_idx] = 1.0
            balance_row[export_slice.start + t_idx] = -1.0
            balance_row[charge_slice.start + t_idx] = -1.0
            balance_row[discharge_slice.start + t_idx] = 1.0
            a_eq.append(balance_row)
            b_eq.append(float(load_array[t_idx] - renewable_array[t_idx]))

            power_row = np.zeros(n_vars)
            power_row[charge_slice.start + t_idx] = 1.0
            power_row[discharge_slice.start + t_idx] = 1.0
            power_row[regulation_slice.start + t_idx] = 1.0
            a_ub.append(power_row)
            b_ub.append(params.max_power_kw)

            peak_row = np.zeros(n_vars)
            peak_row[import_slice.start + t_idx] = 1.0
            peak_row[peak_index] = -1.0
            a_ub.append(peak_row)
            b_ub.append(0.0)

        lp_result = linprog(
            c=c,
            A_ub=np.vstack(a_ub) if a_ub else None,
            b_ub=np.asarray(b_ub, dtype=float) if b_ub else None,
            A_eq=np.vstack(a_eq) if a_eq else None,
            b_eq=np.asarray(b_eq, dtype=float) if b_eq else None,
            bounds=bounds,
            method="highs",
        )

        if not lp_result.success:
            result = BatteryDispatchResult(
                charge_kw=np.zeros(n_t),
                discharge_kw=np.zeros(n_t),
                soc=np.full(n_t + 1, self.battery.state.soc),
                grid_import_kw=np.zeros(n_t),
                grid_export_kw=np.zeros(n_t),
                regulation_kw=np.zeros(n_t),
                objective_value=float("inf"),
                expected_profit=float("-inf"),
                degradation_cost=0.0,
                status="failed",
                message=lp_result.message,
            )
            self._last_result = result
            return result

        x = np.asarray(lp_result.x, dtype=float)
        charge_kw = x[charge_slice]
        discharge_kw = x[discharge_slice]
        grid_import_kw = x[import_slice]
        grid_export_kw = x[export_slice]
        regulation_kw = x[regulation_slice]
        soc = x[soc_slice] / max(params.capacity_kwh, 1e-9)
        degradation_cost = float(np.sum((charge_kw + discharge_kw) * throughput_cost * dt))
        gross_profit = float(np.sum((sell_price * grid_export_kw - buy_price * grid_import_kw + regulation_price * regulation_kw) * dt))
        result = BatteryDispatchResult(
            charge_kw=charge_kw,
            discharge_kw=discharge_kw,
            soc=np.clip(soc, params.min_soc, params.max_soc),
            grid_import_kw=grid_import_kw,
            grid_export_kw=grid_export_kw,
            regulation_kw=regulation_kw,
            objective_value=float(lp_result.fun),
            expected_profit=gross_profit - degradation_cost,
            degradation_cost=degradation_cost,
            status="optimal",
            message=lp_result.message,
        )
        self._rolling_inputs = {
            "prices": prices,
            "load": load_array,
            "renewable": renewable_array,
            "horizon_hours": horizon_hours,
        }
        self._last_result = result
        return result

    def update_rolling(self, new_data: Mapping[str, Any]) -> BatteryDispatchResult:
        """Update the rolling-horizon forecast and re-optimize."""

        if not self._rolling_inputs:
            raise RuntimeError("No prior optimization exists. Call optimize_schedule first.")
        merged: dict[str, Any] = dict(self._rolling_inputs)
        merged.update(new_data)
        return self.optimize_schedule(
            prices=merged["prices"],
            load=merged["load"],
            renewable=merged["renewable"],
            horizon_hours=int(merged.get("horizon_hours", 24)),
        )

    @staticmethod
    def _parse_prices(prices: ArrayLike | Mapping[str, ArrayLike]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if isinstance(prices, Mapping):
            buy = np.asarray(prices.get("buy", prices.get("energy", 0.0)), dtype=float)
            sell = np.asarray(prices.get("sell", buy), dtype=float)
            regulation = np.asarray(prices.get("regulation", np.zeros_like(buy)), dtype=float)
            return buy, sell, regulation
        buy = np.asarray(prices, dtype=float)
        sell = buy.copy()
        regulation = np.zeros_like(buy)
        return buy, sell, regulation


class FleetScheduler:
    """Coordinator for multiple battery systems operating as a VPP."""

    def optimize_fleet(
        self,
        batteries: Sequence[BatteryModel],
        market_signals: Mapping[str, Any],
    ) -> FleetDispatchResult:
        """Optimize fleet dispatch with degradation-aware fair allocation."""

        if not batteries:
            raise ValueError("At least one battery is required.")

        prices = market_signals["prices"]
        load = np.asarray(market_signals.get("load", 0.0), dtype=float)
        renewable = np.asarray(market_signals.get("renewable", np.zeros_like(load)), dtype=float)
        horizon_hours = int(market_signals.get("horizon_hours", 24))

        weights = np.asarray(
            [
                battery.parameters.capacity_kwh
                * max(0.1, 1.0 - battery.state.degradation_fraction)
                / (1.0 + battery.state.equivalent_cycles)
                for battery in batteries
            ],
            dtype=float,
        )
        weights = weights / np.sum(weights)

        per_battery: dict[str, BatteryDispatchResult] = {}
        aggregate_charge = None
        aggregate_discharge = None
        aggregate_soc = None
        utilizations = []

        for idx, (battery, weight) in enumerate(zip(batteries, weights)):
            scheduler = BatteryScheduler(battery)
            result = scheduler.optimize_schedule(
                prices=prices,
                load=load * weight,
                renewable=renewable * weight,
                horizon_hours=horizon_hours,
            )
            per_battery[f"battery_{idx}"] = result
            utilizations.append(float(np.sum(result.discharge_kw + result.charge_kw)))
            aggregate_charge = result.charge_kw if aggregate_charge is None else aggregate_charge + result.charge_kw
            aggregate_discharge = result.discharge_kw if aggregate_discharge is None else aggregate_discharge + result.discharge_kw
            soc_kwh = result.soc * battery.parameters.capacity_kwh
            aggregate_soc = soc_kwh if aggregate_soc is None else aggregate_soc + soc_kwh

        utilization_vector = np.asarray(utilizations, dtype=float)
        fairness_index = float(
            (np.sum(utilization_vector) ** 2)
            / max(len(utilization_vector) * np.sum(utilization_vector**2), 1e-9)
        )
        return FleetDispatchResult(
            per_battery=per_battery,
            aggregate_charge_kw=np.asarray(aggregate_charge, dtype=float),
            aggregate_discharge_kw=np.asarray(aggregate_discharge, dtype=float),
            aggregate_soc_kwh=np.asarray(aggregate_soc, dtype=float),
            fairness_index=fairness_index,
            status="optimal" if all(result.status == "optimal" for result in per_battery.values()) else "degraded",
            message="Fleet optimization completed.",
        )
