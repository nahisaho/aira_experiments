from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Union
import math
import random

import numpy as np
from scipy import optimize, stats


Number = Union[int, float, np.number]
ModelType = Any


def _coerce_rng(random_state: Optional[Union[int, np.random.Generator]] = None) -> np.random.Generator:
    if isinstance(random_state, np.random.Generator):
        return random_state
    return np.random.default_rng(random_state)


def _scalarize(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return float(array)
    return float(array.reshape(-1)[-1])


def _extract_output(raw_output: Any, output_species: Any) -> Any:
    if callable(output_species):
        return output_species(raw_output)
    if output_species is None:
        return raw_output
    if isinstance(raw_output, Mapping):
        return raw_output[output_species]
    if hasattr(raw_output, output_species):
        return getattr(raw_output, output_species)
    try:
        return raw_output[output_species]
    except (TypeError, KeyError, IndexError):
        pass
    raise KeyError(f"Could not extract output species {output_species!r} from model output.")


def _evaluate_model(model: ModelType, parameters: Mapping[str, float], output_species: Any) -> float:
    if callable(model):
        raw_output = model(dict(parameters))
    elif hasattr(model, "simulate"):
        raw_output = model.simulate(dict(parameters))
    elif hasattr(model, "evaluate"):
        raw_output = model.evaluate(dict(parameters))
    else:
        raise TypeError("Model must be callable or provide simulate()/evaluate() methods.")
    return _scalarize(_extract_output(raw_output, output_species))


def _passes_threshold(value: float, target: float, threshold: Any) -> bool:
    if callable(threshold):
        return bool(threshold(value, target))
    deviation = abs(value - target)
    if isinstance(threshold, Mapping):
        if "absolute" in threshold:
            return deviation <= float(threshold["absolute"])
        if "relative" in threshold:
            scale = max(abs(target), 1e-12)
            return deviation / scale <= float(threshold["relative"])
    if isinstance(threshold, (tuple, list)) and len(threshold) == 2:
        lower, upper = threshold
        return (target + float(lower)) <= value <= (target + float(upper))
    return deviation <= float(threshold)


def _design_signature(genotype: Mapping[str, int]) -> Tuple[Tuple[str, int], ...]:
    return tuple(sorted((key, int(value)) for key, value in genotype.items()))


_SOBOL_DIRECTION_DATA: Dict[int, Tuple[int, int, Sequence[int]]] = {
    1: (0, 0, []),
    2: (1, 0, [1]),
    3: (2, 1, [1, 3]),
    4: (3, 1, [1, 3, 1]),
    5: (3, 2, [1, 1, 1]),
    6: (4, 1, [1, 3, 5, 13]),
    7: (4, 4, [1, 1, 5, 5]),
    8: (5, 2, [1, 3, 3, 9, 7]),
    9: (5, 4, [1, 1, 5, 11, 27]),
    10: (5, 7, [1, 1, 7, 13, 3]),
    11: (5, 11, [1, 1, 5, 5, 17]),
    12: (5, 13, [1, 1, 7, 7, 1]),
    13: (5, 14, [1, 3, 1, 15, 11]),
    14: (6, 1, [1, 3, 5, 15, 17, 63]),
    15: (6, 13, [1, 1, 1, 15, 5, 49]),
    16: (6, 16, [1, 3, 3, 5, 19, 61]),
}


class _SobolSequence:
    """Minimal Sobol sequence generator for low-dimensional Saltelli sampling."""

    def __init__(self, dimensions: int, bits: int = 30) -> None:
        if dimensions < 1:
            raise ValueError("Sobol dimensions must be positive.")
        if dimensions > len(_SOBOL_DIRECTION_DATA):
            raise ValueError(
                f"Sobol sampler supports up to {len(_SOBOL_DIRECTION_DATA)} dimensions; got {dimensions}."
            )
        self.dimensions = dimensions
        self.bits = bits
        self._direction_numbers = self._build_direction_numbers()

    def _build_direction_numbers(self) -> np.ndarray:
        direction_numbers = np.zeros((self.dimensions, self.bits + 1), dtype=np.uint32)
        for bit in range(1, self.bits + 1):
            direction_numbers[0, bit] = 1 << (32 - bit)

        for dimension in range(2, self.dimensions + 1):
            degree, coefficient, initial_values = _SOBOL_DIRECTION_DATA[dimension]
            for bit in range(1, degree + 1):
                direction_numbers[dimension - 1, bit] = initial_values[bit - 1] << (32 - bit)
            for bit in range(degree + 1, self.bits + 1):
                value = direction_numbers[dimension - 1, bit - degree]
                value ^= value >> degree
                for offset in range(1, degree):
                    if (coefficient >> (degree - 1 - offset)) & 1:
                        value ^= direction_numbers[dimension - 1, bit - offset]
                direction_numbers[dimension - 1, bit] = value
        return direction_numbers

    def generate(self, n_samples: int) -> np.ndarray:
        if n_samples < 1:
            raise ValueError("n_samples must be positive.")
        samples = np.zeros((n_samples, self.dimensions), dtype=float)
        state = np.zeros(self.dimensions, dtype=np.uint32)
        scale = float(2**32)
        samples[0] = state / scale
        for index in range(1, n_samples):
            gray_index = index - 1
            bit = 1
            while gray_index & 1:
                gray_index >>= 1
                bit += 1
            state ^= self._direction_numbers[:, bit]
            samples[index] = state / scale
        return samples


@dataclass(frozen=True)
class ParameterUncertainty:
    parameter_name: str
    nominal_value: float
    distribution_type: str
    distribution_params: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        distribution_type = self.distribution_type.lower()
        object.__setattr__(self, "distribution_type", distribution_type)
        if distribution_type not in {"uniform", "lognormal", "normal"}:
            raise ValueError("distribution_type must be 'uniform', 'lognormal', or 'normal'.")

    @staticmethod
    def _lognormal_mu_sigma(mean: float, sigma: Optional[float] = None, cv: Optional[float] = None) -> Tuple[float, float]:
        if sigma is None:
            if cv is None:
                raise ValueError("Lognormal distribution requires either sigma or cv.")
            sigma = math.sqrt(math.log1p(float(cv) ** 2))
        sigma = float(sigma)
        mean = float(mean)
        if mean <= 0:
            raise ValueError("Lognormal mean must be positive.")
        mu = math.log(mean) - 0.5 * sigma**2
        return mu, sigma

    def _distribution(self) -> stats.rv_continuous:
        params = dict(self.distribution_params)
        if self.distribution_type == "uniform":
            if "low" in params and "high" in params:
                low = float(params["low"])
                high = float(params["high"])
            else:
                low = float(params.get("loc", self.nominal_value - params.get("scale", 0.0) / 2.0))
                high = low + float(params.get("scale", 0.0))
            if high <= low:
                raise ValueError(f"Uniform bounds invalid for {self.parameter_name}: ({low}, {high}).")
            return stats.uniform(loc=low, scale=high - low)

        if self.distribution_type == "normal":
            mean = float(params.get("mean", self.nominal_value))
            std = float(params.get("std", params.get("sigma", 0.0)))
            if std <= 0:
                raise ValueError(f"Normal std must be positive for {self.parameter_name}.")
            return stats.norm(loc=mean, scale=std)

        mean = float(params.get("mean", self.nominal_value))
        sigma = params.get("sigma")
        cv = params.get("cv")
        if "mu" in params and sigma is not None:
            mu = float(params["mu"])
            sigma = float(sigma)
        else:
            mu, sigma = self._lognormal_mu_sigma(mean=mean, sigma=None if sigma is None else float(sigma), cv=None if cv is None else float(cv))
        return stats.lognorm(s=sigma, scale=math.exp(mu))

    def sample(
        self,
        n: int,
        random_state: Optional[Union[int, np.random.Generator]] = None,
    ) -> np.ndarray:
        if n < 1:
            raise ValueError("n must be positive.")
        rng = _coerce_rng(random_state)
        distribution = self._distribution()
        return np.asarray(distribution.rvs(size=int(n), random_state=rng), dtype=float)

    def get_bounds(self, confidence: float = 0.95) -> Tuple[float, float]:
        if not 0 < confidence < 1:
            raise ValueError("confidence must lie in (0, 1).")
        alpha = (1.0 - confidence) / 2.0
        distribution = self._distribution()
        lower = float(distribution.ppf(alpha))
        upper = float(distribution.ppf(1.0 - alpha))
        return lower, upper

    def unit_to_sample(self, unit_samples: np.ndarray) -> np.ndarray:
        distribution = self._distribution()
        clipped = np.clip(unit_samples, np.finfo(float).eps, 1.0 - np.finfo(float).eps)
        return np.asarray(distribution.ppf(clipped), dtype=float)


class UncertaintySet:
    def __init__(self, uncertainties: Iterable[ParameterUncertainty]) -> None:
        self._uncertainties: List[ParameterUncertainty] = list(uncertainties)
        names = [unc.parameter_name for unc in self._uncertainties]
        if len(set(names)) != len(names):
            raise ValueError("Parameter names in an UncertaintySet must be unique.")

    def __iter__(self):
        return iter(self._uncertainties)

    def __len__(self) -> int:
        return len(self._uncertainties)

    def __getitem__(self, parameter_name: str) -> ParameterUncertainty:
        for uncertainty in self._uncertainties:
            if uncertainty.parameter_name == parameter_name:
                return uncertainty
        raise KeyError(parameter_name)

    @property
    def parameter_names(self) -> List[str]:
        return [uncertainty.parameter_name for uncertainty in self._uncertainties]

    def sample_all(
        self,
        n: int,
        random_state: Optional[Union[int, np.random.Generator]] = None,
    ) -> Dict[str, np.ndarray]:
        rng = _coerce_rng(random_state)
        return {
            uncertainty.parameter_name: uncertainty.sample(n, random_state=rng)
            for uncertainty in self._uncertainties
        }

    def get_nominal(self) -> Dict[str, float]:
        return {
            uncertainty.parameter_name: float(uncertainty.nominal_value)
            for uncertainty in self._uncertainties
        }

    @classmethod
    def default_gene_circuit(
        cls,
        nominal_parameters: Optional[Mapping[str, float]] = None,
        coefficient_of_variation: float = 0.30,
    ) -> "UncertaintySet":
        defaults: Dict[str, float] = {
            "transcription_rate": 1.0,
            "translation_rate": 5.0,
            "mrna_degradation_rate": 0.2,
            "protein_degradation_rate": 0.05,
            "binding_affinity": 10.0,
            "hill_coefficient": 2.0,
            "leakage_rate": 0.05,
        }
        if nominal_parameters is not None:
            defaults.update({key: float(value) for key, value in nominal_parameters.items()})
        uncertainties = [
            ParameterUncertainty(
                parameter_name=name,
                nominal_value=value,
                distribution_type="lognormal",
                distribution_params={"mean": value, "cv": coefficient_of_variation},
            )
            for name, value in defaults.items()
        ]
        return cls(uncertainties)


class RobustnessAnalyzer:
    def __init__(
        self,
        finite_difference_step: float = 1e-6,
        sobol_base_samples: int = 512,
        random_seed: Optional[int] = None,
    ) -> None:
        self.finite_difference_step = float(finite_difference_step)
        self.sobol_base_samples = int(sobol_base_samples)
        self.random_seed = random_seed

    def _resolve_uncertainties(self, param_uncertainties: Union[UncertaintySet, Iterable[ParameterUncertainty]]) -> UncertaintySet:
        if isinstance(param_uncertainties, UncertaintySet):
            return param_uncertainties
        return UncertaintySet(param_uncertainties)

    def _build_parameter_vector(self, names: Sequence[str], values: Sequence[float]) -> Dict[str, float]:
        return {name: float(value) for name, value in zip(names, values)}

    def _saltelli_samples(self, uncertainties: UncertaintySet, n_samples: int) -> Tuple[np.ndarray, np.ndarray]:
        dimensions = len(uncertainties)
        sobol = _SobolSequence(dimensions=2 * dimensions)
        unit_samples = sobol.generate(int(n_samples))
        a_unit = unit_samples[:, :dimensions]
        b_unit = unit_samples[:, dimensions:]
        a = np.zeros_like(a_unit)
        b = np.zeros_like(b_unit)
        for index, uncertainty in enumerate(uncertainties):
            a[:, index] = uncertainty.unit_to_sample(a_unit[:, index])
            b[:, index] = uncertainty.unit_to_sample(b_unit[:, index])
        return a, b

    def sensitivity_analysis(
        self,
        model: ModelType,
        param_uncertainties: Union[UncertaintySet, Iterable[ParameterUncertainty]],
        output_species: Any,
        method: str = "sobol",
    ) -> Dict[str, Any]:
        uncertainties = self._resolve_uncertainties(param_uncertainties)
        method = method.lower()
        if method == "local":
            nominal = uncertainties.get_nominal()
            baseline = _evaluate_model(model, nominal, output_species)
            sensitivities: Dict[str, Dict[str, float]] = {}
            for uncertainty in uncertainties:
                step = self.finite_difference_step * max(abs(uncertainty.nominal_value), 1.0)
                plus = dict(nominal)
                minus = dict(nominal)
                plus[uncertainty.parameter_name] += step
                minus[uncertainty.parameter_name] -= step
                value_plus = _evaluate_model(model, plus, output_species)
                value_minus = _evaluate_model(model, minus, output_species)
                derivative = (value_plus - value_minus) / (2.0 * step)
                normalized = derivative * uncertainty.nominal_value / max(abs(baseline), 1e-12)
                sensitivities[uncertainty.parameter_name] = {
                    "partial_derivative": float(derivative),
                    "normalized_sensitivity": float(normalized),
                }
            return {
                "method": "local",
                "baseline_output": float(baseline),
                "sensitivities": sensitivities,
            }

        if method != "sobol":
            raise ValueError("method must be 'local' or 'sobol'.")

        sample_size = max(self.sobol_base_samples, 64 * len(uncertainties))
        names = uncertainties.parameter_names
        a, b = self._saltelli_samples(uncertainties, sample_size)
        ya = np.array([
            _evaluate_model(model, self._build_parameter_vector(names, row), output_species)
            for row in a
        ])
        yb = np.array([
            _evaluate_model(model, self._build_parameter_vector(names, row), output_species)
            for row in b
        ])
        variance = float(np.var(np.concatenate([ya, yb]), ddof=1))
        if variance <= 1e-16:
            zero_indices = {name: 0.0 for name in names}
            return {
                "method": "sobol",
                "sample_size": sample_size,
                "first_order": zero_indices,
                "total": zero_indices.copy(),
            }

        first_order: Dict[str, float] = {}
        total_order: Dict[str, float] = {}
        for index, name in enumerate(names):
            mixed = np.array(a, copy=True)
            mixed[:, index] = b[:, index]
            ymixed = np.array([
                _evaluate_model(model, self._build_parameter_vector(names, row), output_species)
                for row in mixed
            ])
            s1 = np.mean(yb * (ymixed - ya)) / variance
            st = 0.5 * np.mean((ya - ymixed) ** 2) / variance
            first_order[name] = float(np.clip(s1, 0.0, 1.0))
            total_order[name] = float(np.clip(st, 0.0, 1.0))

        return {
            "method": "sobol",
            "sample_size": sample_size,
            "first_order": first_order,
            "total": total_order,
        }

    def robustness_score(
        self,
        model: ModelType,
        param_uncertainties: Union[UncertaintySet, Iterable[ParameterUncertainty]],
        output_species: Any,
        threshold: Any,
        n_samples: int = 500,
    ) -> float:
        uncertainties = self._resolve_uncertainties(param_uncertainties)
        nominal = uncertainties.get_nominal()
        target = _evaluate_model(model, nominal, output_species)
        samples = uncertainties.sample_all(n_samples, random_state=self.random_seed)
        successes = 0
        for sample_index in range(int(n_samples)):
            parameters = {
                name: float(values[sample_index])
                for name, values in samples.items()
            }
            value = _evaluate_model(model, parameters, output_species)
            if np.isfinite(value) and _passes_threshold(value, target, threshold):
                successes += 1
        return float(successes / max(int(n_samples), 1))

    def worst_case_analysis(
        self,
        model: ModelType,
        param_uncertainties: Union[UncertaintySet, Iterable[ParameterUncertainty]],
        output_species: Any,
    ) -> Dict[str, Any]:
        uncertainties = self._resolve_uncertainties(param_uncertainties)
        names = uncertainties.parameter_names
        nominal_parameters = uncertainties.get_nominal()
        nominal_output = _evaluate_model(model, nominal_parameters, output_species)
        bounds = [uncertainty.get_bounds() for uncertainty in uncertainties]

        def objective(x: np.ndarray) -> float:
            parameters = self._build_parameter_vector(names, x)
            value = _evaluate_model(model, parameters, output_species)
            if not np.isfinite(value):
                return 1e9
            return -abs(value - nominal_output)

        result = optimize.differential_evolution(
            objective,
            bounds=bounds,
            seed=self.random_seed,
            polish=True,
        )
        worst_parameters = self._build_parameter_vector(names, result.x)
        worst_output = _evaluate_model(model, worst_parameters, output_species)
        deviation = abs(worst_output - nominal_output)
        return {
            "worst_case_parameters": worst_parameters,
            "nominal_output": float(nominal_output),
            "worst_case_output": float(worst_output),
            "max_deviation": float(deviation),
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
        }


@dataclass
class OptimizationResult:
    best_design: Mapping[str, Any]
    robustness_score: float
    performance_metrics: Mapping[str, Any]
    pareto_front: List[Tuple[float, float]]
    convergence_history: List[Mapping[str, float]]


class RobustOptimizer:
    def __init__(
        self,
        n_generations: int = 50,
        population_size: int = 30,
        crossover_rate: float = 0.8,
        mutation_rate: float = 0.15,
        tournament_size: int = 3,
        analyzer: Optional[RobustnessAnalyzer] = None,
        random_seed: Optional[int] = None,
        robustness_samples: int = 500,
    ) -> None:
        self.n_generations = int(n_generations)
        self.population_size = int(population_size)
        self.crossover_rate = float(crossover_rate)
        self.mutation_rate = float(mutation_rate)
        self.tournament_size = int(tournament_size)
        self.random_seed = random_seed
        self.robustness_samples = int(robustness_samples)
        self.analyzer = analyzer or RobustnessAnalyzer(random_seed=random_seed)
        self._rng = random.Random(random_seed)

    def _roles(self, parts_catalog: Mapping[str, Sequence[Any]]) -> List[str]:
        roles = list(parts_catalog.keys())
        if not roles:
            raise ValueError("parts_catalog must contain at least one part family.")
        return roles

    def _random_genotype(self, parts_catalog: Mapping[str, Sequence[Any]]) -> Dict[str, int]:
        genotype: Dict[str, int] = {}
        for role in self._roles(parts_catalog):
            options = parts_catalog[role]
            if not options:
                raise ValueError(f"No parts available for role {role!r}.")
            genotype[role] = self._rng.randrange(len(options))
        return genotype

    def _decode_design(self, genotype: Mapping[str, int], parts_catalog: Mapping[str, Sequence[Any]]) -> Dict[str, Any]:
        return {role: parts_catalog[role][index] for role, index in genotype.items()}

    def _constraints_satisfied(
        self,
        design: Mapping[str, Any],
        constraints: Optional[Union[Callable[[Mapping[str, Any]], bool], Sequence[Callable[[Mapping[str, Any]], bool]]]],
    ) -> bool:
        if constraints is None:
            return True
        if callable(constraints):
            return bool(constraints(design))
        return all(bool(constraint(design)) for constraint in constraints)

    def _build_model(
        self,
        design: Mapping[str, Any],
        circuit_spec: Mapping[str, Any],
        objective_result: Any,
    ) -> ModelType:
        if isinstance(objective_result, Mapping) and objective_result.get("model") is not None:
            return objective_result["model"]
        if isinstance(objective_result, (tuple, list)) and len(objective_result) >= 2:
            return objective_result[1]
        model_builder = circuit_spec.get("model_builder")
        if model_builder is not None:
            return model_builder(design)
        if "model" in design:
            return design["model"]
        raise ValueError("A model is required from objective_func(...) or circuit_spec['model_builder'].")

    def _parse_performance(self, objective_result: Any) -> Tuple[float, Mapping[str, Any], bool]:
        if isinstance(objective_result, Mapping):
            performance = float(
                objective_result.get(
                    "performance_score",
                    objective_result.get("performance", objective_result.get("objective", 0.0)),
                )
            )
            metrics = dict(objective_result.get("performance_metrics", {}))
            if not metrics:
                metrics = {k: v for k, v in objective_result.items() if k not in {"model", "multi_objective"}}
            multi_objective = bool(objective_result.get("multi_objective", False))
            return performance, metrics, multi_objective
        if isinstance(objective_result, (tuple, list)):
            performance = float(objective_result[0])
            metrics = {"performance_score": performance}
            multi_objective = len(objective_result) >= 3 and bool(objective_result[2].get("multi_objective", False)) if isinstance(objective_result[2], Mapping) else False if len(objective_result) >= 3 else False
            return performance, metrics, multi_objective
        performance = float(objective_result)
        return performance, {"performance_score": performance}, False

    def _call_objective(
        self,
        objective_func: Callable[..., Any],
        design: Mapping[str, Any],
        circuit_spec: Mapping[str, Any],
        parts_catalog: Mapping[str, Sequence[Any]],
    ) -> Any:
        try:
            return objective_func(design, circuit_spec, parts_catalog)
        except TypeError:
            try:
                return objective_func(design, circuit_spec)
            except TypeError:
                return objective_func(design)

    def _evaluate_candidate(
        self,
        genotype: Mapping[str, int],
        circuit_spec: Mapping[str, Any],
        parts_catalog: Mapping[str, Sequence[Any]],
        uncertainties: Union[UncertaintySet, Iterable[ParameterUncertainty]],
        objective_func: Callable[..., Any],
        constraints: Optional[Union[Callable[[Mapping[str, Any]], bool], Sequence[Callable[[Mapping[str, Any]], bool]]]],
        cache: MutableMapping[Tuple[Tuple[str, int], ...], Dict[str, Any]],
    ) -> Dict[str, Any]:
        signature = _design_signature(genotype)
        if signature in cache:
            return cache[signature]

        design = self._decode_design(genotype, parts_catalog)
        if not self._constraints_satisfied(design, constraints):
            evaluation = {
                "design": design,
                "performance": 0.0,
                "robustness": 0.0,
                "fitness": 0.0,
                "metrics": {"constraint_satisfied": False},
            }
            cache[signature] = evaluation
            return evaluation

        objective_result = self._call_objective(objective_func, design, circuit_spec, parts_catalog)
        performance, metrics, _ = self._parse_performance(objective_result)
        model = self._build_model(design, circuit_spec, objective_result)
        threshold = circuit_spec.get("threshold", circuit_spec.get("robustness_threshold"))
        if threshold is None:
            raise ValueError("circuit_spec must define 'threshold' or 'robustness_threshold'.")
        output_species = circuit_spec.get("output_species")
        robustness = self.analyzer.robustness_score(
            model=model,
            param_uncertainties=uncertainties,
            output_species=output_species,
            threshold=threshold,
            n_samples=self.robustness_samples,
        )
        fitness = max(performance, 0.0) * robustness
        evaluation = {
            "design": design,
            "performance": float(performance),
            "robustness": float(robustness),
            "fitness": float(fitness),
            "metrics": metrics,
        }
        cache[signature] = evaluation
        return evaluation

    def _tournament_select(self, population: Sequence[Mapping[str, int]], evaluations: Sequence[Mapping[str, Any]]) -> Dict[str, int]:
        indices = self._rng.sample(range(len(population)), k=min(self.tournament_size, len(population)))
        best_index = max(indices, key=lambda idx: evaluations[idx]["fitness"])
        return dict(population[best_index])

    def _crossover(self, parent_a: Mapping[str, int], parent_b: Mapping[str, int]) -> Tuple[Dict[str, int], Dict[str, int]]:
        if self._rng.random() >= self.crossover_rate:
            return dict(parent_a), dict(parent_b)
        child_a = dict(parent_a)
        child_b = dict(parent_b)
        for role in parent_a:
            if self._rng.random() < 0.5:
                child_a[role], child_b[role] = child_b[role], child_a[role]
        return child_a, child_b

    def _mutate(self, genotype: Mapping[str, int], parts_catalog: Mapping[str, Sequence[Any]]) -> Dict[str, int]:
        mutated = dict(genotype)
        if self._rng.random() < self.mutation_rate:
            role = self._rng.choice(self._roles(parts_catalog))
            mutated[role] = self._rng.randrange(len(parts_catalog[role]))
        return mutated

    @staticmethod
    def _pareto_front(candidates: Sequence[Mapping[str, Any]]) -> List[Tuple[float, float]]:
        front: List[Tuple[float, float]] = []
        for candidate in candidates:
            point = (float(candidate["performance"]), float(candidate["robustness"]))
            dominated = False
            for other in candidates:
                other_point = (float(other["performance"]), float(other["robustness"]))
                if (
                    other_point[0] >= point[0]
                    and other_point[1] >= point[1]
                    and other_point != point
                ):
                    dominated = True
                    break
            if not dominated and point not in front:
                front.append(point)
        front.sort(key=lambda pair: (pair[0], pair[1]), reverse=True)
        return front

    def optimize(
        self,
        circuit_spec: Mapping[str, Any],
        parts_catalog: Mapping[str, Sequence[Any]],
        uncertainties: Union[UncertaintySet, Iterable[ParameterUncertainty]],
        objective_func: Callable[..., Any],
        constraints: Optional[Union[Callable[[Mapping[str, Any]], bool], Sequence[Callable[[Mapping[str, Any]], bool]]]] = None,
    ) -> OptimizationResult:
        population = [self._random_genotype(parts_catalog) for _ in range(self.population_size)]
        cache: Dict[Tuple[Tuple[str, int], ...], Dict[str, Any]] = {}
        evaluated_candidates: List[Dict[str, Any]] = []
        convergence_history: List[Mapping[str, float]] = []

        for generation in range(self.n_generations):
            evaluations = [
                self._evaluate_candidate(
                    genotype=genotype,
                    circuit_spec=circuit_spec,
                    parts_catalog=parts_catalog,
                    uncertainties=uncertainties,
                    objective_func=objective_func,
                    constraints=constraints,
                    cache=cache,
                )
                for genotype in population
            ]
            evaluated_candidates.extend(evaluations)
            best = max(evaluations, key=lambda item: item["fitness"])
            convergence_history.append(
                {
                    "generation": float(generation),
                    "best_fitness": float(best["fitness"]),
                    "best_performance": float(best["performance"]),
                    "best_robustness": float(best["robustness"]),
                    "mean_fitness": float(np.mean([item["fitness"] for item in evaluations])),
                }
            )

            next_population: List[Dict[str, int]] = [dict(population[evaluations.index(best)])]
            while len(next_population) < self.population_size:
                parent_a = self._tournament_select(population, evaluations)
                parent_b = self._tournament_select(population, evaluations)
                child_a, child_b = self._crossover(parent_a, parent_b)
                next_population.append(self._mutate(child_a, parts_catalog))
                if len(next_population) < self.population_size:
                    next_population.append(self._mutate(child_b, parts_catalog))
            population = next_population

        final_evaluations = [
            self._evaluate_candidate(
                genotype=genotype,
                circuit_spec=circuit_spec,
                parts_catalog=parts_catalog,
                uncertainties=uncertainties,
                objective_func=objective_func,
                constraints=constraints,
                cache=cache,
            )
            for genotype in population
        ]
        evaluated_candidates.extend(final_evaluations)
        best = max(final_evaluations, key=lambda item: item["fitness"])
        pareto_front = self._pareto_front(evaluated_candidates)
        return OptimizationResult(
            best_design=best["design"],
            robustness_score=float(best["robustness"]),
            performance_metrics=best["metrics"],
            pareto_front=pareto_front,
            convergence_history=convergence_history,
        )


__all__ = [
    "OptimizationResult",
    "ParameterUncertainty",
    "RobustOptimizer",
    "RobustnessAnalyzer",
    "UncertaintySet",
]
