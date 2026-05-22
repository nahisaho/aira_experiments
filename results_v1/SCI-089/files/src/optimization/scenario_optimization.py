"""Stochastic and robust optimization models for power system scheduling."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog


ArrayLike = Sequence[float] | np.ndarray


@dataclass(slots=True)
class Scenario:
    """Single operating scenario for stochastic optimization."""

    demand: np.ndarray
    renewable: np.ndarray
    reserve: np.ndarray
    probability: float = 1.0
    name: str = "scenario"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.demand = np.asarray(self.demand, dtype=float)
        self.renewable = np.asarray(self.renewable, dtype=float)
        self.reserve = np.asarray(self.reserve, dtype=float)
        if not (len(self.demand) == len(self.renewable) == len(self.reserve)):
            raise ValueError("Scenario demand, renewable, and reserve arrays must share the same length.")
        if self.probability < 0.0:
            raise ValueError("Scenario probability must be non-negative.")

    @property
    def horizon(self) -> int:
        """Return the number of time periods in the scenario."""

        return int(self.demand.size)

    def feature_vector(self) -> np.ndarray:
        """Return a compact representation used for scenario reduction."""

        return np.concatenate((self.demand - self.renewable, self.reserve))


@dataclass(slots=True)
class ScenarioOptimizationResult:
    """Structured result for the stochastic optimization solve."""

    commitment: np.ndarray
    dispatch: np.ndarray
    reserve: np.ndarray
    load_shedding: np.ndarray
    renewable_spillage: np.ndarray
    reserve_shortfall: np.ndarray
    expected_cost: float
    cvar_value: float
    objective_value: float
    status: str
    message: str


@dataclass(slots=True)
class ScenarioReductionResult:
    """Reduced scenario set and redistributed probabilities."""

    scenarios: list[Scenario]
    retained_indices: list[int]
    probability_map: dict[int, float]
    distance_matrix: np.ndarray


@dataclass(slots=True)
class UncertaintySet:
    """Description of box or polyhedral uncertainty for robust optimization."""

    kind: str
    demand_nominal: np.ndarray
    renewable_nominal: np.ndarray
    demand_bounds: Optional[np.ndarray] = None
    renewable_bounds: Optional[np.ndarray] = None
    a_matrix: Optional[np.ndarray] = None
    b_vector: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        self.kind = self.kind.lower()
        self.demand_nominal = np.asarray(self.demand_nominal, dtype=float)
        self.renewable_nominal = np.asarray(self.renewable_nominal, dtype=float)
        if self.demand_nominal.shape != self.renewable_nominal.shape:
            raise ValueError("Demand and renewable nominal vectors must have the same shape.")
        if self.demand_bounds is not None:
            self.demand_bounds = np.asarray(self.demand_bounds, dtype=float)
        if self.renewable_bounds is not None:
            self.renewable_bounds = np.asarray(self.renewable_bounds, dtype=float)
        if self.a_matrix is not None:
            self.a_matrix = np.asarray(self.a_matrix, dtype=float)
        if self.b_vector is not None:
            self.b_vector = np.asarray(self.b_vector, dtype=float)


@dataclass(slots=True)
class RobustOptimizationResult:
    """Result container for robust optimization."""

    worst_case_demand: np.ndarray
    worst_case_renewable: np.ndarray
    scenario_result: ScenarioOptimizationResult
    status: str
    message: str


class ScenarioOptimizer:
    """Two-stage stochastic scenario optimizer with a CVaR risk term.

    The model uses an LP relaxation of unit commitment by relaxing binary
    commitment variables to the interval ``[0, 1]`` and solves the resulting
    problem with :func:`scipy.optimize.linprog`.
    """

    def __init__(
        self,
        risk_alpha: float = 0.95,
        cvar_weight: float = 1.0,
        load_shed_penalty: float = 10_000.0,
        reserve_penalty: float = 5_000.0,
        spill_penalty: float = 1.0,
    ) -> None:
        if not 0.0 < risk_alpha < 1.0:
            raise ValueError("risk_alpha must lie strictly between 0 and 1.")
        self.risk_alpha = risk_alpha
        self.cvar_weight = cvar_weight
        self.load_shed_penalty = load_shed_penalty
        self.reserve_penalty = reserve_penalty
        self.spill_penalty = spill_penalty
        self.network: dict[str, Any] | None = None
        self.scenarios: list[Scenario] = []
        self.weights: np.ndarray | None = None
        self.model: dict[str, Any] = {}
        self.solution: ScenarioOptimizationResult | None = None

    @staticmethod
    def generate_scenarios_from_forecast(
        probabilistic_forecast: Mapping[str, ArrayLike],
        n_scenarios: int,
        seed: Optional[int] = None,
    ) -> list[Scenario]:
        """Generate scenarios from Gaussian probabilistic forecasts."""

        rng = np.random.default_rng(seed)
        demand_mean = np.asarray(probabilistic_forecast["demand_mean"], dtype=float)
        demand_std = np.asarray(probabilistic_forecast.get("demand_std", np.zeros_like(demand_mean)), dtype=float)
        renewable_mean = np.asarray(
            probabilistic_forecast.get("renewable_mean", np.zeros_like(demand_mean)),
            dtype=float,
        )
        renewable_std = np.asarray(
            probabilistic_forecast.get("renewable_std", np.zeros_like(demand_mean)),
            dtype=float,
        )
        reserve_mean = np.asarray(
            probabilistic_forecast.get("reserve_mean", np.full_like(demand_mean, 0.05 * np.max(demand_mean))),
            dtype=float,
        )
        reserve_std = np.asarray(
            probabilistic_forecast.get("reserve_std", np.zeros_like(demand_mean)),
            dtype=float,
        )

        scenarios: list[Scenario] = []
        for idx in range(n_scenarios):
            demand = np.maximum(0.0, rng.normal(demand_mean, demand_std))
            renewable = np.maximum(0.0, rng.normal(renewable_mean, renewable_std))
            reserve = np.maximum(0.0, rng.normal(reserve_mean, reserve_std))
            scenarios.append(
                Scenario(
                    demand=demand,
                    renewable=renewable,
                    reserve=reserve,
                    probability=1.0 / n_scenarios,
                    name=f"scenario_{idx}",
                )
            )
        return scenarios

    def build_model(
        self,
        network: Mapping[str, Any],
        scenarios: Sequence[Scenario | Mapping[str, Any]],
        weights: Optional[Sequence[float]] = None,
    ) -> dict[str, Any]:
        """Build the LP relaxation for the stochastic scheduling model."""

        self.network = dict(network)
        self.scenarios = [self._coerce_scenario(scenario, idx) for idx, scenario in enumerate(scenarios)]
        if not self.scenarios:
            raise ValueError("At least one scenario is required.")

        horizon = self.scenarios[0].horizon
        if any(scenario.horizon != horizon for scenario in self.scenarios):
            raise ValueError("All scenarios must share the same horizon length.")

        generators = self.network.get("generators", [])
        if not generators:
            raise ValueError("Network must include a non-empty 'generators' list.")

        scenario_weights = np.asarray(
            weights if weights is not None else [scenario.probability for scenario in self.scenarios],
            dtype=float,
        )
        if scenario_weights.size != len(self.scenarios):
            raise ValueError("Scenario weights must match the number of scenarios.")
        if np.any(scenario_weights < 0.0):
            raise ValueError("Scenario weights must be non-negative.")
        total_weight = float(np.sum(scenario_weights))
        if total_weight <= 0.0:
            raise ValueError("Scenario weights must sum to a positive value.")
        self.weights = scenario_weights / total_weight

        n_gen = len(generators)
        n_scn = len(self.scenarios)
        n_t = horizon

        offsets: dict[str, Any] = {}
        cursor = 0
        offsets["u"] = (cursor, cursor + n_gen * n_t)
        cursor = offsets["u"][1]
        offsets["p"] = (cursor, cursor + n_scn * n_gen * n_t)
        cursor = offsets["p"][1]
        offsets["r"] = (cursor, cursor + n_scn * n_gen * n_t)
        cursor = offsets["r"][1]
        offsets["shed"] = (cursor, cursor + n_scn * n_t)
        cursor = offsets["shed"][1]
        offsets["spill"] = (cursor, cursor + n_scn * n_t)
        cursor = offsets["spill"][1]
        offsets["reserve_short"] = (cursor, cursor + n_scn * n_t)
        cursor = offsets["reserve_short"][1]
        offsets["eta"] = cursor
        cursor += 1
        offsets["z"] = (cursor, cursor + n_scn)
        cursor = offsets["z"][1]
        n_vars = cursor

        c = np.zeros(n_vars)
        for g_idx, generator in enumerate(generators):
            commitment_cost = float(generator.get("commitment_cost", generator.get("startup_cost", 0.0)))
            reserve_cost = float(generator.get("reserve_cost", 0.0))
            variable_cost = float(generator.get("variable_cost", generator.get("marginal_cost", generator.get("cost", 0.0))))
            for t_idx in range(n_t):
                c[self._u_index(offsets, g_idx, t_idx, n_t)] = commitment_cost
            for s_idx, weight in enumerate(self.weights):
                for t_idx in range(n_t):
                    c[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = weight * variable_cost
                    c[self._pg_index(offsets, "r", s_idx, g_idx, t_idx, n_gen, n_t)] = weight * reserve_cost
        for s_idx, weight in enumerate(self.weights):
            for t_idx in range(n_t):
                c[self._st_index(offsets, "shed", s_idx, t_idx, n_t)] = weight * self.load_shed_penalty
                c[self._st_index(offsets, "spill", s_idx, t_idx, n_t)] = weight * self.spill_penalty
                c[self._st_index(offsets, "reserve_short", s_idx, t_idx, n_t)] = weight * self.reserve_penalty
            c[offsets["z"][0] + s_idx] = self.cvar_weight * weight / (1.0 - self.risk_alpha)
        c[offsets["eta"]] = self.cvar_weight

        bounds: list[tuple[float | None, float | None]] = []
        bounds.extend([(0.0, 1.0)] * (n_gen * n_t))
        bounds.extend([(0.0, None)] * (2 * n_scn * n_gen * n_t))
        bounds.extend([(0.0, None)] * (3 * n_scn * n_t))
        bounds.append((None, None))
        bounds.extend([(0.0, None)] * n_scn)

        a_eq: list[np.ndarray] = []
        b_eq: list[float] = []
        a_ub: list[np.ndarray] = []
        b_ub: list[float] = []

        reserve_factor = float(self.network.get("reserve_factor", 1.0))
        for s_idx, scenario in enumerate(self.scenarios):
            net_load = scenario.demand - scenario.renewable
            for t_idx in range(n_t):
                row = np.zeros(n_vars)
                for g_idx in range(n_gen):
                    row[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = 1.0
                row[self._st_index(offsets, "shed", s_idx, t_idx, n_t)] = 1.0
                row[self._st_index(offsets, "spill", s_idx, t_idx, n_t)] = -1.0
                a_eq.append(row)
                b_eq.append(float(net_load[t_idx]))

                reserve_row = np.zeros(n_vars)
                for g_idx in range(n_gen):
                    reserve_row[self._pg_index(offsets, "r", s_idx, g_idx, t_idx, n_gen, n_t)] = -1.0
                reserve_row[self._st_index(offsets, "reserve_short", s_idx, t_idx, n_t)] = -1.0
                a_ub.append(reserve_row)
                b_ub.append(-float(reserve_factor * scenario.reserve[t_idx]))

            for g_idx, generator in enumerate(generators):
                min_output = float(generator.get("min_output", 0.0))
                max_output = float(generator.get("max_output", generator.get("capacity", 0.0)))
                ramp_rate = float(generator.get("ramp_rate", max_output))
                availability = np.asarray(generator.get("availability", np.ones(n_t)), dtype=float)
                if availability.size == 1:
                    availability = np.full(n_t, float(availability[0]))
                if availability.size != n_t:
                    raise ValueError("Generator availability vectors must match the scenario horizon.")

                for t_idx in range(n_t):
                    upper = np.zeros(n_vars)
                    upper[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = 1.0
                    upper[self._pg_index(offsets, "r", s_idx, g_idx, t_idx, n_gen, n_t)] = 1.0
                    upper[self._u_index(offsets, g_idx, t_idx, n_t)] = -max_output * float(availability[t_idx])
                    a_ub.append(upper)
                    b_ub.append(0.0)

                    lower = np.zeros(n_vars)
                    lower[self._u_index(offsets, g_idx, t_idx, n_t)] = min_output
                    lower[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = -1.0
                    a_ub.append(lower)
                    b_ub.append(0.0)

                    if t_idx == 0:
                        continue
                    ramp_up = np.zeros(n_vars)
                    ramp_up[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = 1.0
                    ramp_up[self._pg_index(offsets, "p", s_idx, g_idx, t_idx - 1, n_gen, n_t)] = -1.0
                    a_ub.append(ramp_up)
                    b_ub.append(ramp_rate)

                    ramp_down = np.zeros(n_vars)
                    ramp_down[self._pg_index(offsets, "p", s_idx, g_idx, t_idx - 1, n_gen, n_t)] = 1.0
                    ramp_down[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = -1.0
                    a_ub.append(ramp_down)
                    b_ub.append(ramp_rate)

        for s_idx, scenario in enumerate(self.scenarios):
            row = np.zeros(n_vars)
            for g_idx, generator in enumerate(generators):
                variable_cost = float(generator.get("variable_cost", generator.get("marginal_cost", generator.get("cost", 0.0))))
                reserve_cost = float(generator.get("reserve_cost", 0.0))
                commitment_cost = float(generator.get("commitment_cost", generator.get("startup_cost", 0.0)))
                for t_idx in range(n_t):
                    row[self._pg_index(offsets, "p", s_idx, g_idx, t_idx, n_gen, n_t)] = variable_cost
                    row[self._pg_index(offsets, "r", s_idx, g_idx, t_idx, n_gen, n_t)] = reserve_cost
                    row[self._u_index(offsets, g_idx, t_idx, n_t)] += commitment_cost
            for t_idx in range(n_t):
                row[self._st_index(offsets, "shed", s_idx, t_idx, n_t)] = self.load_shed_penalty
                row[self._st_index(offsets, "spill", s_idx, t_idx, n_t)] = self.spill_penalty
                row[self._st_index(offsets, "reserve_short", s_idx, t_idx, n_t)] = self.reserve_penalty
            row[offsets["eta"]] = -1.0
            row[offsets["z"][0] + s_idx] = -1.0
            a_ub.append(row)
            b_ub.append(0.0)

        self.model = {
            "c": c,
            "A_ub": np.vstack(a_ub) if a_ub else None,
            "b_ub": np.asarray(b_ub, dtype=float) if b_ub else None,
            "A_eq": np.vstack(a_eq) if a_eq else None,
            "b_eq": np.asarray(b_eq, dtype=float) if b_eq else None,
            "bounds": bounds,
            "offsets": offsets,
            "n_gen": n_gen,
            "n_scn": n_scn,
            "n_t": n_t,
            "generators": generators,
        }
        return self.model

    def solve(self) -> ScenarioOptimizationResult:
        """Solve the built stochastic optimization model."""

        if not self.model:
            raise RuntimeError("Call build_model before solve.")

        result = linprog(
            c=self.model["c"],
            A_ub=self.model["A_ub"],
            b_ub=self.model["b_ub"],
            A_eq=self.model["A_eq"],
            b_eq=self.model["b_eq"],
            bounds=self.model["bounds"],
            method="highs",
        )
        self.solution = self._decode_solution(result)
        return self.solution

    def get_solution(self) -> ScenarioOptimizationResult:
        """Return the latest optimization solution."""

        if self.solution is None:
            raise RuntimeError("No solution is available. Call solve first.")
        return self.solution

    def _decode_solution(self, result: Any) -> ScenarioOptimizationResult:
        offsets = self.model["offsets"]
        n_gen = self.model["n_gen"]
        n_scn = self.model["n_scn"]
        n_t = self.model["n_t"]
        weights = np.asarray(self.weights, dtype=float)

        if not result.success:
            empty_commitment = np.zeros((n_gen, n_t))
            empty_dispatch = np.zeros((n_scn, n_gen, n_t))
            empty_scalar = np.zeros((n_scn, n_t))
            return ScenarioOptimizationResult(
                commitment=empty_commitment,
                dispatch=empty_dispatch,
                reserve=empty_dispatch.copy(),
                load_shedding=empty_scalar.copy(),
                renewable_spillage=empty_scalar.copy(),
                reserve_shortfall=empty_scalar.copy(),
                expected_cost=float("inf"),
                cvar_value=float("inf"),
                objective_value=float("inf"),
                status="failed",
                message=result.message,
            )

        x = np.asarray(result.x, dtype=float)
        commitment = x[offsets["u"][0] : offsets["u"][1]].reshape(n_gen, n_t)
        dispatch = x[offsets["p"][0] : offsets["p"][1]].reshape(n_scn, n_gen, n_t)
        reserve = x[offsets["r"][0] : offsets["r"][1]].reshape(n_scn, n_gen, n_t)
        load_shedding = x[offsets["shed"][0] : offsets["shed"][1]].reshape(n_scn, n_t)
        renewable_spillage = x[offsets["spill"][0] : offsets["spill"][1]].reshape(n_scn, n_t)
        reserve_shortfall = x[offsets["reserve_short"][0] : offsets["reserve_short"][1]].reshape(n_scn, n_t)
        eta = float(x[offsets["eta"]])
        z = x[offsets["z"][0] : offsets["z"][1]]

        scenario_costs = np.zeros(n_scn)
        commitment_cost = 0.0
        for g_idx, generator in enumerate(self.model["generators"]):
            commitment_cost += float(generator.get("commitment_cost", generator.get("startup_cost", 0.0))) * float(
                np.sum(commitment[g_idx])
            )
            variable_cost = float(generator.get("variable_cost", generator.get("marginal_cost", generator.get("cost", 0.0))))
            reserve_cost = float(generator.get("reserve_cost", 0.0))
            scenario_costs += variable_cost * np.sum(dispatch[:, g_idx, :], axis=1)
            scenario_costs += reserve_cost * np.sum(reserve[:, g_idx, :], axis=1)
        scenario_costs += commitment_cost
        scenario_costs += self.load_shed_penalty * np.sum(load_shedding, axis=1)
        scenario_costs += self.spill_penalty * np.sum(renewable_spillage, axis=1)
        scenario_costs += self.reserve_penalty * np.sum(reserve_shortfall, axis=1)

        cvar_value = eta + float(np.dot(weights, z)) / (1.0 - self.risk_alpha)
        expected_cost = float(np.dot(weights, scenario_costs))
        return ScenarioOptimizationResult(
            commitment=commitment,
            dispatch=dispatch,
            reserve=reserve,
            load_shedding=load_shedding,
            renewable_spillage=renewable_spillage,
            reserve_shortfall=reserve_shortfall,
            expected_cost=expected_cost,
            cvar_value=cvar_value,
            objective_value=float(result.fun),
            status="optimal",
            message=result.message,
        )

    @staticmethod
    def _u_index(offsets: Mapping[str, Any], g_idx: int, t_idx: int, n_t: int) -> int:
        return offsets["u"][0] + g_idx * n_t + t_idx

    @staticmethod
    def _pg_index(
        offsets: Mapping[str, Any],
        key: str,
        s_idx: int,
        g_idx: int,
        t_idx: int,
        n_gen: int,
        n_t: int,
    ) -> int:
        return offsets[key][0] + s_idx * n_gen * n_t + g_idx * n_t + t_idx

    @staticmethod
    def _st_index(offsets: Mapping[str, Any], key: str, s_idx: int, t_idx: int, n_t: int) -> int:
        return offsets[key][0] + s_idx * n_t + t_idx

    @staticmethod
    def _coerce_scenario(scenario: Scenario | Mapping[str, Any], idx: int) -> Scenario:
        if isinstance(scenario, Scenario):
            return scenario
        demand = np.asarray(scenario["demand"], dtype=float)
        renewable = np.asarray(scenario.get("renewable", np.zeros_like(demand)), dtype=float)
        reserve = np.asarray(scenario.get("reserve", np.zeros_like(demand)), dtype=float)
        probability = float(scenario.get("probability", 1.0))
        return Scenario(
            demand=demand,
            renewable=renewable,
            reserve=reserve,
            probability=probability,
            name=str(scenario.get("name", f"scenario_{idx}")),
            metadata=dict(scenario.get("metadata", {})),
        )


class ScenarioReducer:
    """Scenario reduction using weighted distance-based forward/backward selection."""

    def compute_distances(self, scenarios: Sequence[Scenario | Mapping[str, Any]]) -> np.ndarray:
        """Compute a pairwise Kantorovich-style distance matrix."""

        scenario_objects = [ScenarioOptimizer._coerce_scenario(scenario, idx) for idx, scenario in enumerate(scenarios)]
        n_scenarios = len(scenario_objects)
        distance_matrix = np.zeros((n_scenarios, n_scenarios), dtype=float)
        for i in range(n_scenarios):
            vector_i = scenario_objects[i].feature_vector()
            for j in range(i + 1, n_scenarios):
                vector_j = scenario_objects[j].feature_vector()
                distance = float(np.mean(np.abs(vector_i - vector_j)))
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
        return distance_matrix

    def reduce(
        self,
        scenarios: Sequence[Scenario | Mapping[str, Any]],
        n_reduced: int,
        method: str = "forward",
    ) -> ScenarioReductionResult:
        """Reduce scenarios and redistribute probability mass to retained scenarios."""

        scenario_objects = [ScenarioOptimizer._coerce_scenario(scenario, idx) for idx, scenario in enumerate(scenarios)]
        if not scenario_objects:
            raise ValueError("At least one scenario is required.")
        if not 1 <= n_reduced <= len(scenario_objects):
            raise ValueError("n_reduced must be between 1 and the number of scenarios.")

        probabilities = np.asarray([scenario.probability for scenario in scenario_objects], dtype=float)
        probabilities = probabilities / np.sum(probabilities)
        distances = self.compute_distances(scenario_objects)
        method = method.lower()
        if method == "forward":
            retained = self._forward_select(distances, probabilities, n_reduced)
        elif method == "backward":
            retained = self._backward_select(distances, probabilities, n_reduced)
        else:
            raise ValueError("method must be either 'forward' or 'backward'.")

        redistributed = {index: 0.0 for index in retained}
        for idx, probability in enumerate(probabilities):
            if idx in redistributed:
                redistributed[idx] += float(probability)
                continue
            nearest = min(retained, key=lambda candidate: distances[idx, candidate])
            redistributed[nearest] += float(probability)

        reduced_scenarios: list[Scenario] = []
        for index in retained:
            scenario = scenario_objects[index]
            reduced_scenarios.append(
                Scenario(
                    demand=scenario.demand.copy(),
                    renewable=scenario.renewable.copy(),
                    reserve=scenario.reserve.copy(),
                    probability=redistributed[index],
                    name=scenario.name,
                    metadata=dict(scenario.metadata),
                )
            )

        return ScenarioReductionResult(
            scenarios=reduced_scenarios,
            retained_indices=retained,
            probability_map=redistributed,
            distance_matrix=distances,
        )

    @staticmethod
    def _forward_select(distances: np.ndarray, probabilities: np.ndarray, n_reduced: int) -> list[int]:
        n_scenarios = distances.shape[0]
        centroid_proxy = probabilities @ distances
        selected = [int(np.argmin(centroid_proxy))]
        while len(selected) < n_reduced:
            candidates = [idx for idx in range(n_scenarios) if idx not in selected]
            best = max(
                candidates,
                key=lambda idx: probabilities[idx] * min(distances[idx, ref] for ref in selected),
            )
            selected.append(int(best))
        return sorted(selected)

    @staticmethod
    def _backward_select(distances: np.ndarray, probabilities: np.ndarray, n_reduced: int) -> list[int]:
        retained = list(range(distances.shape[0]))
        while len(retained) > n_reduced:
            removable = min(
                retained,
                key=lambda idx: probabilities[idx]
                * min(
                    (distances[idx, other] for other in retained if other != idx),
                    default=0.0,
                ),
            )
            retained.remove(int(removable))
        return sorted(retained)


class RobustOptimizer:
    """Robust optimizer for box and polyhedral uncertainty sets."""

    def __init__(self, scenario_optimizer: Optional[ScenarioOptimizer] = None) -> None:
        self.scenario_optimizer = scenario_optimizer or ScenarioOptimizer(cvar_weight=0.0)

    def solve_robust(
        self,
        network: Mapping[str, Any],
        uncertainty_set: UncertaintySet | Mapping[str, Any],
    ) -> RobustOptimizationResult:
        """Solve a worst-case robust dispatch problem."""

        uncertainty = self._coerce_uncertainty_set(uncertainty_set)
        worst_case_demand = self._worst_case_series(
            uncertainty.demand_nominal,
            uncertainty.demand_bounds,
            1.0,
            uncertainty,
        )
        worst_case_renewable = self._worst_case_series(
            uncertainty.renewable_nominal,
            uncertainty.renewable_bounds,
            -1.0,
            uncertainty,
        )
        scenario = Scenario(
            demand=worst_case_demand,
            renewable=worst_case_renewable,
            reserve=np.asarray(network.get("reserve_requirement", np.zeros_like(worst_case_demand)), dtype=float),
            probability=1.0,
            name="worst_case",
        )
        self.scenario_optimizer.build_model(network, [scenario], [1.0])
        scenario_result = self.scenario_optimizer.solve()
        status = "optimal" if scenario_result.status == "optimal" else "failed"
        return RobustOptimizationResult(
            worst_case_demand=worst_case_demand,
            worst_case_renewable=worst_case_renewable,
            scenario_result=scenario_result,
            status=status,
            message=scenario_result.message,
        )

    def _worst_case_series(
        self,
        nominal: np.ndarray,
        bounds: Optional[np.ndarray],
        direction: float,
        uncertainty: UncertaintySet,
    ) -> np.ndarray:
        nominal = np.asarray(nominal, dtype=float)
        if bounds is None:
            bounds = np.zeros_like(nominal)
        bounds = np.asarray(bounds, dtype=float)
        if uncertainty.kind == "box":
            return nominal + direction * np.abs(bounds)
        if uncertainty.kind != "polyhedral":
            raise ValueError("Supported uncertainty kinds are 'box' and 'polyhedral'.")
        if uncertainty.a_matrix is None or uncertainty.b_vector is None:
            return nominal + direction * np.abs(bounds)

        objective = -direction * np.ones(nominal.size, dtype=float)
        lp_result = linprog(
            c=objective,
            A_ub=uncertainty.a_matrix,
            b_ub=uncertainty.b_vector,
            bounds=[(-float(abs(value)), float(abs(value))) for value in bounds],
            method="highs",
        )
        if not lp_result.success:
            return nominal + direction * np.abs(bounds)
        return nominal + lp_result.x

    @staticmethod
    def _coerce_uncertainty_set(uncertainty_set: UncertaintySet | Mapping[str, Any]) -> UncertaintySet:
        if isinstance(uncertainty_set, UncertaintySet):
            return uncertainty_set
        demand_nominal = np.asarray(uncertainty_set["demand_nominal"], dtype=float)
        renewable_nominal = np.asarray(
            uncertainty_set.get("renewable_nominal", np.zeros_like(demand_nominal)),
            dtype=float,
        )
        return UncertaintySet(
            kind=str(uncertainty_set["kind"]),
            demand_nominal=demand_nominal,
            renewable_nominal=renewable_nominal,
            demand_bounds=np.asarray(uncertainty_set.get("demand_bounds")) if uncertainty_set.get("demand_bounds") is not None else None,
            renewable_bounds=np.asarray(uncertainty_set.get("renewable_bounds")) if uncertainty_set.get("renewable_bounds") is not None else None,
            a_matrix=np.asarray(uncertainty_set.get("a_matrix")) if uncertainty_set.get("a_matrix") is not None else None,
            b_vector=np.asarray(uncertainty_set.get("b_vector")) if uncertainty_set.get("b_vector") is not None else None,
        )
