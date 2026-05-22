"""Demand response models and coordinated BESS/DR scheduling."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog

from .battery_scheduling import BatteryModel, FleetDispatchResult, FleetScheduler


ArrayLike = Sequence[float] | np.ndarray


@dataclass(slots=True)
class LoadFlexibilityParameters:
    """Parameters describing customer-side flexibility."""

    category: str
    max_shift_hours: int = 0
    max_curtailment_fraction: float = 0.2
    comfort_min: Optional[float] = None
    comfort_max: Optional[float] = None
    process_deadline: Optional[int] = None
    rebound_ratio: float = 0.15
    participation_base: float = 0.7
    participation_sensitivity: float = 0.02
    price_elasticity: float = 0.1


@dataclass(slots=True)
class FlexibleLoad:
    """Flexible load with a baseline profile and operating parameters."""

    name: str
    baseline_kw: np.ndarray
    parameters: LoadFlexibilityParameters
    ambient_temperature: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        self.baseline_kw = np.asarray(self.baseline_kw, dtype=float)
        if self.ambient_temperature is not None:
            self.ambient_temperature = np.asarray(self.ambient_temperature, dtype=float)
            if self.ambient_temperature.size != self.baseline_kw.size:
                raise ValueError("ambient_temperature must align with baseline_kw.")


@dataclass(slots=True)
class FlexibilityProfile:
    """Computed DR flexibility for a single load."""

    name: str
    shiftable_kw: np.ndarray
    curtailable_kw: np.ndarray
    interruptible_kw: np.ndarray
    participation_probability: np.ndarray
    rebound_kw: np.ndarray


@dataclass(slots=True)
class DRScheduleResult:
    """Demand response optimization output."""

    reductions_kw: dict[str, np.ndarray]
    incentives: np.ndarray
    adjusted_load_kw: np.ndarray
    expected_savings: float
    objective_value: float
    status: str
    message: str


@dataclass(slots=True)
class CoordinatedDispatchResult:
    """Joint DR and BESS coordination result."""

    dr_result: DRScheduleResult
    fleet_result: FleetDispatchResult
    residual_history: list[float]
    dual_signal: np.ndarray
    status: str
    message: str


class DemandResponseModel:
    """Flexible-load model including participation and rebound effects."""

    def __init__(self, indoor_setpoint_c: float = 22.0) -> None:
        self.indoor_setpoint_c = indoor_setpoint_c

    def model_flexibility(self, load_profile: FlexibleLoad | Mapping[str, Any]) -> FlexibilityProfile:
        """Translate a load profile into DR potentials by category."""

        load = self._coerce_load(load_profile)
        profile = load.baseline_kw
        params = load.parameters
        ambient = load.ambient_temperature

        shift_fraction = 0.3 if params.category == "shiftable" else 0.05
        curtailable_fraction = params.max_curtailment_fraction if params.category in {"curtailable", "interruptible"} else 0.1
        interruptible_fraction = 0.8 if params.category == "interruptible" else 0.0

        thermal_stress = np.zeros_like(profile)
        if ambient is not None and params.comfort_max is not None:
            thermal_stress = np.maximum(ambient - params.comfort_max, 0.0)
        comfort_modifier = np.clip(1.0 - 0.03 * thermal_stress, 0.4, 1.0)
        participation = np.clip(
            params.participation_base + params.participation_sensitivity * comfort_modifier,
            0.05,
            0.99,
        )

        shiftable = profile * shift_fraction * comfort_modifier
        curtailable = profile * curtailable_fraction * comfort_modifier
        interruptible = profile * interruptible_fraction * comfort_modifier
        rebound = (shiftable + 0.5 * curtailable) * params.rebound_ratio
        return FlexibilityProfile(
            name=load.name,
            shiftable_kw=shiftable,
            curtailable_kw=curtailable,
            interruptible_kw=interruptible,
            participation_probability=np.asarray(participation, dtype=float),
            rebound_kw=rebound,
        )

    def compute_dr_potential(
        self,
        load_profiles: Sequence[FlexibleLoad | Mapping[str, Any]],
    ) -> dict[str, np.ndarray]:
        """Aggregate demand response potential across all flexible loads."""

        modeled = [self.model_flexibility(load_profile) for load_profile in load_profiles]
        total_shiftable = np.sum([profile.shiftable_kw for profile in modeled], axis=0)
        total_curtailable = np.sum([profile.curtailable_kw for profile in modeled], axis=0)
        total_interruptible = np.sum([profile.interruptible_kw for profile in modeled], axis=0)
        total_rebound = np.sum([profile.rebound_kw for profile in modeled], axis=0)
        return {
            "shiftable_kw": total_shiftable,
            "curtailable_kw": total_curtailable,
            "interruptible_kw": total_interruptible,
            "rebound_kw": total_rebound,
        }

    @staticmethod
    def _coerce_load(load_profile: FlexibleLoad | Mapping[str, Any]) -> FlexibleLoad:
        if isinstance(load_profile, FlexibleLoad):
            return load_profile
        parameters = load_profile["parameters"]
        if not isinstance(parameters, LoadFlexibilityParameters):
            parameters = LoadFlexibilityParameters(**parameters)
        return FlexibleLoad(
            name=str(load_profile["name"]),
            baseline_kw=np.asarray(load_profile["baseline_kw"], dtype=float),
            parameters=parameters,
            ambient_temperature=np.asarray(load_profile["ambient_temperature"], dtype=float)
            if load_profile.get("ambient_temperature") is not None
            else None,
        )


class DRScheduler:
    """Optimal scheduler for incentive-based and price-responsive DR."""

    def __init__(self, base_incentive: float = 20.0, incentive_cap: float = 200.0) -> None:
        self.base_incentive = base_incentive
        self.incentive_cap = incentive_cap
        self.model = DemandResponseModel()

    def schedule_dr(
        self,
        loads: Sequence[FlexibleLoad | Mapping[str, Any]],
        prices: ArrayLike,
        constraints: Mapping[str, Any],
    ) -> DRScheduleResult:
        """Schedule flexible demand while respecting customer-side constraints."""

        load_objects = [self.model._coerce_load(load) for load in loads]
        price_array = np.asarray(prices, dtype=float)
        if not load_objects:
            raise ValueError("At least one flexible load is required.")
        horizon = min([load.baseline_kw.size for load in load_objects] + [price_array.size])
        price_array = price_array[:horizon]

        flexibility = {load.name: self.model.model_flexibility(load) for load in load_objects}
        incentives = self.compute_incentives(price_array, flexibility.values())
        battery_schedule = constraints.get("battery_schedule")
        battery_support = np.zeros(horizon)
        if isinstance(battery_schedule, FleetDispatchResult):
            raw_support = battery_schedule.aggregate_discharge_kw[:horizon] - battery_schedule.aggregate_charge_kw[:horizon]
            battery_support = np.maximum(raw_support, 0.0)

        n_loads = len(load_objects)
        n_vars = n_loads * horizon
        c = np.zeros(n_vars)
        bounds: list[tuple[float | None, float | None]] = []

        baseline_total = np.sum([load.baseline_kw[:horizon] for load in load_objects], axis=0)
        for l_idx, load in enumerate(load_objects):
            profile = flexibility[load.name]
            availability = (
                profile.shiftable_kw + profile.curtailable_kw + profile.interruptible_kw
            ) * profile.participation_probability
            elasticity_scale = 1.0 + load.parameters.price_elasticity * (price_array / max(np.mean(price_array), 1e-9) - 1.0)
            availability = np.maximum(availability * elasticity_scale, 0.0)
            discomfort_cost = 0.25 * incentives + 0.1 * np.maximum(price_array - np.mean(price_array), 0.0)
            c[l_idx * horizon : (l_idx + 1) * horizon] = -(price_array - discomfort_cost)
            bounds.extend((0.0, float(value)) for value in availability)

        a_ub: list[np.ndarray] = []
        b_ub: list[float] = []
        max_net_load = constraints.get("max_net_load")
        if max_net_load is not None:
            net_limit = np.asarray(max_net_load, dtype=float)
            if net_limit.size == 1:
                net_limit = np.full(horizon, float(net_limit))
            for t_idx in range(horizon):
                row = np.zeros(n_vars)
                for l_idx in range(n_loads):
                    row[l_idx * horizon + t_idx] = -1.0
                a_ub.append(row)
                b_ub.append(float(net_limit[t_idx] - baseline_total[t_idx] + battery_support[t_idx]))

        for l_idx, load in enumerate(load_objects):
            deadline = load.parameters.process_deadline
            if deadline is None:
                continue
            row = np.zeros(n_vars)
            for t_idx in range(min(deadline + 1, horizon)):
                row[l_idx * horizon + t_idx] = 1.0
            max_shiftable = float(np.sum(flexibility[load.name].shiftable_kw[: min(deadline + 1, horizon)]))
            a_ub.append(row)
            b_ub.append(max_shiftable)

        lp_result = linprog(
            c=c,
            A_ub=np.vstack(a_ub) if a_ub else None,
            b_ub=np.asarray(b_ub, dtype=float) if b_ub else None,
            bounds=bounds,
            method="highs",
        )

        if not lp_result.success:
            return DRScheduleResult(
                reductions_kw={load.name: np.zeros(horizon) for load in load_objects},
                incentives=incentives,
                adjusted_load_kw=baseline_total.copy(),
                expected_savings=0.0,
                objective_value=float("inf"),
                status="failed",
                message=lp_result.message,
            )

        x = np.asarray(lp_result.x, dtype=float).reshape(n_loads, horizon)
        reductions = {load.name: x[idx] for idx, load in enumerate(load_objects)}
        rebound = self._build_rebound(load_objects, reductions)
        adjusted_load = baseline_total - np.sum(x, axis=0) + rebound - battery_support
        dr_value = float(np.sum((price_array * np.sum(x, axis=0) - incentives * np.sum(x, axis=0)) * 0.25))
        return DRScheduleResult(
            reductions_kw=reductions,
            incentives=incentives,
            adjusted_load_kw=adjusted_load,
            expected_savings=dr_value,
            objective_value=float(lp_result.fun),
            status="optimal",
            message=lp_result.message,
        )

    def compute_incentives(
        self,
        prices: ArrayLike,
        flexibility_profiles: Optional[Sequence[FlexibilityProfile]] = None,
    ) -> np.ndarray:
        """Compute incentive levels from price scarcity and flexibility depth."""

        price_array = np.asarray(prices, dtype=float)
        normalized_price = (price_array - np.min(price_array)) / max(np.ptp(price_array), 1e-9)
        scarcity = np.ones_like(price_array)
        if flexibility_profiles is not None:
            available = np.sum(
                [
                    profile.shiftable_kw + profile.curtailable_kw + profile.interruptible_kw
                    for profile in flexibility_profiles
                ],
                axis=0,
            )
            scarcity = 1.0 / np.maximum(available / max(np.max(available), 1e-9), 0.1)
        incentives = self.base_incentive * (1.0 + normalized_price) * scarcity
        return np.clip(incentives, 0.0, self.incentive_cap)

    def _build_rebound(
        self,
        loads: Sequence[FlexibleLoad],
        reductions: Mapping[str, np.ndarray],
    ) -> np.ndarray:
        horizon = next(iter(reductions.values())).size
        rebound = np.zeros(horizon)
        for load in loads:
            reduction = reductions[load.name]
            shift_window = max(1, load.parameters.max_shift_hours * 4)
            rebound_energy = reduction * load.parameters.rebound_ratio
            for t_idx, value in enumerate(rebound_energy):
                if value <= 0.0:
                    continue
                target = min(horizon - 1, t_idx + shift_window)
                rebound[target] += value
        return rebound


class CoordinatedScheduler:
    """Hierarchical DR + BESS scheduler with ADMM-style coordination."""

    def __init__(self, max_iterations: int = 8, rho: float = 1.0) -> None:
        self.max_iterations = max_iterations
        self.rho = rho
        self.dr_scheduler = DRScheduler()
        self.fleet_scheduler = FleetScheduler()

    def optimize_joint(
        self,
        batteries: Sequence[BatteryModel],
        dr_resources: Sequence[FlexibleLoad | Mapping[str, Any]],
        system_constraints: Mapping[str, Any],
    ) -> CoordinatedDispatchResult:
        """Jointly optimize batteries and demand response resources."""

        prices = np.asarray(system_constraints["prices"], dtype=float)
        load = np.asarray(system_constraints.get("load"), dtype=float)
        renewable = np.asarray(system_constraints.get("renewable", np.zeros_like(load)), dtype=float)
        net_limit = np.asarray(system_constraints.get("max_net_load", np.full_like(load, np.max(load))), dtype=float)
        if net_limit.size == 1:
            net_limit = np.full(load.size, float(net_limit))

        dual_signal = np.zeros_like(prices, dtype=float)
        residual_history: list[float] = []
        dr_result: DRScheduleResult | None = None
        fleet_result: FleetDispatchResult | None = None

        for _ in range(self.max_iterations):
            adjusted_prices = prices + dual_signal
            fleet_result = self.fleet_scheduler.optimize_fleet(
                batteries,
                {
                    "prices": {
                        "buy": adjusted_prices,
                        "sell": adjusted_prices,
                        "regulation": system_constraints.get("regulation_prices", np.zeros_like(prices)),
                    },
                    "load": load,
                    "renewable": renewable,
                    "horizon_hours": int(system_constraints.get("horizon_hours", 24)),
                },
            )
            dr_result = self.dr_scheduler.schedule_dr(
                dr_resources,
                adjusted_prices,
                {
                    "max_net_load": net_limit,
                    "battery_schedule": fleet_result,
                },
            )
            net_grid = dr_result.adjusted_load_kw - renewable[: dr_result.adjusted_load_kw.size]
            residual = net_grid - net_limit[: net_grid.size]
            dual_signal[: residual.size] = np.maximum(0.0, dual_signal[: residual.size] + self.rho * residual)
            residual_history.append(float(np.linalg.norm(np.maximum(residual, 0.0))))
            if residual_history[-1] <= 1e-3:
                break

        if dr_result is None or fleet_result is None:
            raise RuntimeError("Coordination failed to produce a schedule.")
        status = "optimal" if dr_result.status == "optimal" and fleet_result.status == "optimal" else "degraded"
        return CoordinatedDispatchResult(
            dr_result=dr_result,
            fleet_result=fleet_result,
            residual_history=residual_history,
            dual_signal=dual_signal,
            status=status,
            message="Joint optimization completed.",
        )
