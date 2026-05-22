from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.special import comb


PropensityCallable = Callable[..., float]


@dataclass(frozen=True)
class Reaction:
    """Chemical reaction with optional custom propensity."""

    reactants: Dict[str, int]
    products: Dict[str, int]
    rate_constant: float
    propensity_func: Optional[PropensityCallable] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "reactants", self._validate_stoichiometry(self.reactants, "reactants"))
        object.__setattr__(self, "products", self._validate_stoichiometry(self.products, "products"))
        if self.rate_constant < 0:
            raise ValueError("rate_constant must be non-negative.")

    @staticmethod
    def _validate_stoichiometry(stoich: Mapping[str, int], name: str) -> Dict[str, int]:
        cleaned: Dict[str, int] = {}
        for species, coefficient in stoich.items():
            coefficient = int(coefficient)
            if coefficient < 0:
                raise ValueError(f"{name} stoichiometry must be non-negative for {species!r}.")
            if coefficient > 0:
                cleaned[str(species)] = coefficient
        return cleaned

    def propensity(self, state: np.ndarray, species_index: Mapping[str, int]) -> float:
        if self.propensity_func is not None:
            try:
                return max(0.0, float(self.propensity_func(state, species_index)))
            except TypeError:
                return max(0.0, float(self.propensity_func(state)))

        propensity = float(self.rate_constant)
        for species, stoich in self.reactants.items():
            count = int(state[species_index[species]])
            if count < stoich:
                return 0.0
            propensity *= float(comb(count, stoich, exact=False))
            if propensity <= 0.0:
                return 0.0
        return propensity

    def stoichiometry_vector(self, species_index: Mapping[str, int]) -> np.ndarray:
        vector = np.zeros(len(species_index), dtype=np.int64)
        for species, stoich in self.products.items():
            vector[species_index[species]] += stoich
        for species, stoich in self.reactants.items():
            vector[species_index[species]] -= stoich
        return vector


class GeneCircuitModel:
    """Stochastic gene circuit model with SSA- and tau-leaping-ready state."""

    DEFAULT_PARAMETERS = {
        "transcription": 5.0,
        "translation": 1.0,
        "mrna_degradation": 0.2,
        "protein_degradation": 0.01,
        "binding_rate": 0.01,
        "unbinding_rate": 0.1,
        "hill_coefficient": 2.0,
        "dissociation_constant": 50.0,
        "activation_fold": 5.0,
        "leakiness": 0.0,
    }
    PROMOTER_STRENGTHS = {
        "very_weak": 0.1,
        "weak": 0.5,
        "medium": 5.0,
        "strong": 20.0,
        "very_strong": 50.0,
    }
    RBS_STRENGTHS = {
        "very_weak": 0.01,
        "weak": 0.05,
        "medium": 0.5,
        "strong": 2.0,
        "very_strong": 10.0,
    }

    def __init__(
        self,
        species: Mapping[str, int],
        reactions: Sequence[Reaction],
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> None:
        self.species = {str(name): int(count) for name, count in species.items()}
        self.reactions = list(reactions)
        self.metadata = dict(metadata or {})
        self.species_names = list(self.species.keys())
        self.species_index = {name: idx for idx, name in enumerate(self.species_names)}
        self.initial_state = np.array([self.species[name] for name in self.species_names], dtype=np.int64)
        if np.any(self.initial_state < 0):
            raise ValueError("Initial species counts must be non-negative.")
        if self.reactions:
            self.stoichiometry_matrix = np.column_stack(
                [reaction.stoichiometry_vector(self.species_index) for reaction in self.reactions]
            ).astype(np.int64, copy=False)
        else:
            self.stoichiometry_matrix = np.zeros((len(self.species_names), 0), dtype=np.int64)

    @classmethod
    def from_circuit_spec(
        cls,
        circuit_spec: Mapping[str, Any],
        part_assignments: Optional[Mapping[str, Any]] = None,
    ) -> "GeneCircuitModel":
        """
        Build a reaction network from a circuit specification.

        Expected structure:
        {
            "species": {"Regulator": 25},
            "defaults": {...},
            "genes": [
                {
                    "name": "GFP",
                    "promoter": "P1",
                    "promoter_strength": "strong",
                    "rbs": "B0034",
                    "rbs_strength": "medium",
                    "copy_number": 1,
                    "initial_mrna": 0,
                    "initial_protein": 0,
                    "regulators": [
                        {
                            "protein": "LacI",
                            "mode": "repression",
                            "Kd": 40,
                            "hill": 2,
                            "binding_rate": 0.02,
                            "unbinding_rate": 0.1,
                        }
                    ],
                }
            ],
        }
        """

        part_assignments = dict(part_assignments or {})
        defaults = dict(cls.DEFAULT_PARAMETERS)
        defaults.update(circuit_spec.get("defaults", {}))
        defaults.update(circuit_spec.get("parameters", {}))

        initial_species = dict(circuit_spec.get("species", {}))
        initial_species.update(circuit_spec.get("initial_species", {}))
        species: Dict[str, int] = {str(name): int(count) for name, count in initial_species.items()}
        reactions: List[Reaction] = []

        for gene in cls._normalize_genes(circuit_spec):
            name = str(gene.get("name") or gene.get("protein") or gene.get("id"))
            promoter_name = str(gene.get("promoter", f"{name}_promoter"))
            mrna_name = str(gene.get("mrna", f"{name}_mRNA"))
            protein_name = str(gene.get("protein", name))
            copy_number = int(gene.get("copy_number", gene.get("promoter_copy_number", 1)))
            if copy_number < 0:
                raise ValueError(f"copy_number must be non-negative for gene {name!r}.")

            species.setdefault(promoter_name, copy_number)
            species.setdefault(mrna_name, int(gene.get("initial_mrna", 0)))
            species.setdefault(protein_name, int(gene.get("initial_protein", 0)))

            promoter_rate = cls._resolve_rate(
                part_assignments=part_assignments,
                part_type="promoter",
                part_name=promoter_name,
                gene_value=gene.get("promoter_strength", gene.get("transcription_rate")),
                default_value=defaults["transcription"],
                strength_map=cls.PROMOTER_STRENGTHS,
                bounds=(0.1, 50.0),
            )
            translation_rate = cls._resolve_rate(
                part_assignments=part_assignments,
                part_type="rbs",
                part_name=str(gene.get("rbs", gene.get("utr", f"{name}_rbs"))),
                gene_value=gene.get("rbs_strength", gene.get("translation_rate")),
                default_value=defaults["translation"],
                strength_map=cls.RBS_STRENGTHS,
                bounds=(0.01, 10.0),
            )
            mrna_deg = float(gene.get("mrna_degradation", defaults["mrna_degradation"]))
            protein_deg = float(gene.get("protein_degradation", defaults["protein_degradation"]))
            if mrna_deg < 0 or protein_deg < 0:
                raise ValueError(f"Degradation rates must be non-negative for gene {name!r}.")

            regulators = cls._normalize_regulators(gene.get("regulators", []))
            for regulator in regulators:
                regulator_name = str(regulator.get("protein") or regulator.get("regulator"))
                if not regulator_name:
                    raise ValueError(f"Regulator protein missing for gene {name!r}.")
                species.setdefault(regulator_name, int(regulator.get("initial_count", 0)))
                bound_species = str(regulator.get("bound_species", f"{promoter_name}_bound_{regulator_name}"))
                species.setdefault(bound_species, 0)

            reactions.append(
                Reaction(
                    reactants={promoter_name: 1},
                    products={promoter_name: 1, mrna_name: 1},
                    rate_constant=promoter_rate,
                    propensity_func=cls._build_transcription_propensity(
                        promoter_name=promoter_name,
                        base_rate=promoter_rate,
                        regulators=regulators,
                        defaults=defaults,
                    ),
                )
            )
            reactions.append(
                Reaction(
                    reactants={mrna_name: 1},
                    products={mrna_name: 1, protein_name: 1},
                    rate_constant=translation_rate,
                )
            )
            reactions.append(Reaction(reactants={mrna_name: 1}, products={}, rate_constant=mrna_deg))
            reactions.append(Reaction(reactants={protein_name: 1}, products={}, rate_constant=protein_deg))

            for regulator in regulators:
                regulator_name = str(regulator.get("protein") or regulator.get("regulator"))
                mode = str(regulator.get("mode", regulator.get("type", "repression"))).lower()
                bound_species = str(regulator.get("bound_species", f"{promoter_name}_bound_{regulator_name}"))
                binding_rate = float(regulator.get("binding_rate", defaults["binding_rate"]))
                unbinding_rate = float(regulator.get("unbinding_rate", defaults["unbinding_rate"]))
                if binding_rate < 0 or unbinding_rate < 0:
                    raise ValueError(f"Binding rates must be non-negative for regulator {regulator_name!r}.")

                reactions.append(
                    Reaction(
                        reactants={promoter_name: 1, regulator_name: 1},
                        products={bound_species: 1},
                        rate_constant=binding_rate,
                    )
                )
                reactions.append(
                    Reaction(
                        reactants={bound_species: 1},
                        products={promoter_name: 1, regulator_name: 1},
                        rate_constant=unbinding_rate,
                    )
                )

                if mode.startswith("activ"):
                    activation_fold = float(regulator.get("activation_fold", defaults["activation_fold"]))
                    active_rate = promoter_rate * max(1.0, activation_fold)
                    reactions.append(
                        Reaction(
                            reactants={bound_species: 1},
                            products={bound_species: 1, mrna_name: 1},
                            rate_constant=active_rate,
                        )
                    )
                elif mode.startswith("repress"):
                    leakiness = float(regulator.get("leakiness", defaults["leakiness"]))
                    if leakiness > 0:
                        reactions.append(
                            Reaction(
                                reactants={bound_species: 1},
                                products={bound_species: 1, mrna_name: 1},
                                rate_constant=promoter_rate * leakiness,
                            )
                        )

        return cls(species=species, reactions=reactions, metadata={"circuit_spec": dict(circuit_spec)})

    @staticmethod
    def _normalize_genes(circuit_spec: Mapping[str, Any]) -> List[Dict[str, Any]]:
        genes = circuit_spec.get("genes")
        if genes is None:
            excluded = {"species", "initial_species", "defaults", "parameters", "parts"}
            candidates = {k: v for k, v in circuit_spec.items() if k not in excluded}
            if candidates and all(isinstance(v, Mapping) for v in candidates.values()):
                genes = candidates
            else:
                raise ValueError("circuit_spec must contain a 'genes' collection or gene mapping.")

        if isinstance(genes, Mapping):
            normalized: List[Dict[str, Any]] = []
            for gene_name, gene_spec in genes.items():
                if not isinstance(gene_spec, Mapping):
                    raise TypeError(f"Gene specification for {gene_name!r} must be a mapping.")
                entry = dict(gene_spec)
                entry.setdefault("name", gene_name)
                normalized.append(entry)
            return normalized

        normalized = []
        for gene_spec in genes:
            if not isinstance(gene_spec, Mapping):
                raise TypeError("Each gene specification must be a mapping.")
            normalized.append(dict(gene_spec))
        return normalized

    @staticmethod
    def _normalize_regulators(regulators: Any) -> List[Dict[str, Any]]:
        if regulators is None:
            return []
        if isinstance(regulators, Mapping) and "protein" not in regulators and "regulator" not in regulators:
            normalized: List[Dict[str, Any]] = []
            for protein_name, spec in regulators.items():
                if isinstance(spec, Mapping):
                    entry = dict(spec)
                    entry.setdefault("protein", protein_name)
                else:
                    entry = {"protein": protein_name, "mode": spec}
                normalized.append(entry)
            return normalized
        if isinstance(regulators, Mapping):
            return [dict(regulators)]
        normalized = []
        for regulator in regulators:
            if isinstance(regulator, Mapping):
                normalized.append(dict(regulator))
            else:
                normalized.append({"protein": str(regulator), "mode": "repression"})
        return normalized

    @classmethod
    def _resolve_rate(
        cls,
        part_assignments: Mapping[str, Any],
        part_type: str,
        part_name: str,
        gene_value: Any,
        default_value: float,
        strength_map: Mapping[str, float],
        bounds: Tuple[float, float],
    ) -> float:
        if gene_value is not None:
            raw_value = gene_value
        else:
            raw_value = cls._lookup_part_assignment(part_assignments, part_type, part_name)
        if isinstance(raw_value, Mapping):
            for key in ("strength", "rate", "rate_constant", f"{part_type}_strength"):
                if key in raw_value:
                    raw_value = raw_value[key]
                    break
        rate = cls._coerce_strength(raw_value, strength_map, default_value)
        return float(np.clip(rate, bounds[0], bounds[1]))

    @staticmethod
    def _lookup_part_assignment(part_assignments: Mapping[str, Any], part_type: str, part_name: str) -> Any:
        if part_name in part_assignments:
            return part_assignments[part_name]
        for container_name in (part_type, f"{part_type}s", "parts"):
            container = part_assignments.get(container_name)
            if isinstance(container, Mapping) and part_name in container:
                return container[part_name]
        return None

    @staticmethod
    def _coerce_strength(value: Any, strength_map: Mapping[str, float], default_value: float) -> float:
        if value is None:
            return float(default_value)
        if isinstance(value, (int, float, np.integer, np.floating)):
            return float(value)
        key = str(value).strip().lower().replace(" ", "_")
        if key in strength_map:
            return float(strength_map[key])
        try:
            return float(value)
        except (TypeError, ValueError):
            return float(default_value)

    @staticmethod
    def _build_transcription_propensity(
        promoter_name: str,
        base_rate: float,
        regulators: Sequence[Mapping[str, Any]],
        defaults: Mapping[str, float],
    ) -> PropensityCallable:
        def propensity(state: np.ndarray, species_index: Mapping[str, int]) -> float:
            promoter_count = float(state[species_index[promoter_name]])
            if promoter_count <= 0:
                return 0.0

            modulation = 1.0
            for regulator in regulators:
                regulator_name = str(regulator.get("protein") or regulator.get("regulator"))
                if regulator_name not in species_index:
                    continue
                regulator_count = max(0.0, float(state[species_index[regulator_name]]))
                kd = max(
                    1e-12,
                    float(
                        regulator.get(
                            "Kd",
                            regulator.get(
                                "dissociation_constant",
                                defaults["dissociation_constant"],
                            ),
                        )
                    ),
                )
                hill = max(
                    1.0,
                    float(
                        regulator.get(
                            "hill",
                            regulator.get("hill_coefficient", defaults["hill_coefficient"]),
                        )
                    ),
                )
                mode = str(regulator.get("mode", regulator.get("type", "repression"))).lower()
                scaled = (regulator_count / kd) ** hill if regulator_count > 0 else 0.0
                if mode.startswith("repress") or mode.startswith("inhib"):
                    modulation *= 1.0 / (1.0 + scaled)
                elif mode.startswith("activ"):
                    activation_fold = float(regulator.get("activation_fold", defaults["activation_fold"]))
                    hill_activation = scaled / (1.0 + scaled)
                    modulation *= 1.0 + max(0.0, activation_fold - 1.0) * hill_activation
            return max(0.0, base_rate * promoter_count * modulation)

        return propensity

    def propensities(self, state: np.ndarray) -> np.ndarray:
        if not self.reactions:
            return np.zeros(0, dtype=float)
        return np.fromiter(
            (reaction.propensity(state, self.species_index) for reaction in self.reactions),
            dtype=float,
            count=len(self.reactions),
        )


class GillespieSimulator:
    """Standard Gillespie direct-method SSA simulator."""

    def __init__(self, max_steps: int = int(1e6)) -> None:
        self.max_steps = int(max_steps)

    def simulate(
        self,
        model: GeneCircuitModel,
        t_end: float,
        seed: Optional[int] = None,
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        if t_end < 0:
            raise ValueError("t_end must be non-negative.")

        rng = np.random.default_rng(seed)
        state = model.initial_state.copy()
        time_points = [0.0]
        trajectory = [state.copy()]
        t = 0.0

        for _ in range(self.max_steps):
            if t >= t_end:
                break
            propensities = model.propensities(state)
            total_propensity = float(propensities.sum())
            if total_propensity <= 0.0:
                if t < t_end:
                    time_points.append(float(t_end))
                    trajectory.append(state.copy())
                break

            dt = float(rng.exponential(1.0 / total_propensity))
            if t + dt > t_end:
                time_points.append(float(t_end))
                trajectory.append(state.copy())
                break

            reaction_index = int(
                np.searchsorted(np.cumsum(propensities), rng.random() * total_propensity, side="right")
            )
            state = state + model.stoichiometry_matrix[:, reaction_index]
            state = np.maximum(state, 0)
            t += dt
            time_points.append(t)
            trajectory.append(state.copy())
        else:
            raise RuntimeError("Gillespie simulation exceeded max_steps.")

        if len(time_points) == 1 and t_end > 0:
            time_points.append(float(t_end))
            trajectory.append(state.copy())

        trajectory_array = np.vstack(trajectory)
        species_trajectories = {
            species: trajectory_array[:, idx].copy() for idx, species in enumerate(model.species_names)
        }
        return np.asarray(time_points, dtype=float), species_trajectories


class TauLeapingSimulator:
    """Tau-leaping simulator with adaptive tau and exact fallback at low counts."""

    def __init__(
        self,
        max_steps: int = int(1e6),
        epsilon: float = 0.03,
        exact_threshold: int = 10,
    ) -> None:
        self.max_steps = int(max_steps)
        self.epsilon = float(epsilon)
        self.exact_threshold = int(exact_threshold)

    def _exact_step(
        self,
        model: GeneCircuitModel,
        state: np.ndarray,
        t: float,
        t_end: float,
        rng: np.random.Generator,
    ) -> Tuple[float, np.ndarray, bool]:
        propensities = model.propensities(state)
        total_propensity = float(propensities.sum())
        if total_propensity <= 0.0:
            return t_end, state.copy(), False

        dt = float(rng.exponential(1.0 / total_propensity))
        if t + dt > t_end:
            return t_end, state.copy(), False

        reaction_index = int(
            np.searchsorted(np.cumsum(propensities), rng.random() * total_propensity, side="right")
        )
        new_state = state + model.stoichiometry_matrix[:, reaction_index]
        return t + dt, np.maximum(new_state, 0), True

    def simulate(
        self,
        model: GeneCircuitModel,
        t_end: float,
        tau: float = 0.1,
        seed: Optional[int] = None,
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        if t_end < 0:
            raise ValueError("t_end must be non-negative.")
        if tau <= 0:
            raise ValueError("tau must be positive.")

        rng = np.random.default_rng(seed)
        state = model.initial_state.copy()
        time_points = [0.0]
        trajectory = [state.copy()]
        t = 0.0
        steps = 0

        while t < t_end:
            if steps >= self.max_steps:
                raise RuntimeError("Tau-leaping simulation exceeded max_steps.")

            if np.any(state < self.exact_threshold):
                t, state, fired = self._exact_step(model, state, t, t_end, rng)
                if time_points[-1] != t:
                    time_points.append(float(t))
                    trajectory.append(state.copy())
                if not fired:
                    break
                steps += 1
                continue

            propensities = model.propensities(state)
            total_propensity = float(propensities.sum())
            if total_propensity <= 0.0:
                if t < t_end:
                    time_points.append(float(t_end))
                    trajectory.append(state.copy())
                break

            drift = model.stoichiometry_matrix @ propensities
            mask = (state > 0) & (np.abs(drift) > 1e-12)
            adaptive_tau = float(tau)
            if np.any(mask):
                adaptive_tau = float(np.min(self.epsilon * state[mask] / np.abs(drift[mask])))
                if not np.isfinite(adaptive_tau) or adaptive_tau <= 0.0:
                    adaptive_tau = float(tau)
            tau_step = min(float(tau), adaptive_tau, float(t_end - t))

            if tau_step <= 1e-12:
                t, state, fired = self._exact_step(model, state, t, t_end, rng)
                if time_points[-1] != t:
                    time_points.append(float(t))
                    trajectory.append(state.copy())
                if not fired:
                    break
                steps += 1
                continue

            firings = rng.poisson(propensities * tau_step)
            state = state + model.stoichiometry_matrix @ firings
            state = np.maximum(state, 0)
            t += tau_step
            time_points.append(float(t))
            trajectory.append(state.copy())
            steps += 1

        if len(time_points) == 1 and t_end > 0:
            time_points.append(float(t_end))
            trajectory.append(state.copy())
        elif time_points[-1] < t_end:
            time_points.append(float(t_end))
            trajectory.append(state.copy())

        trajectory_array = np.vstack(trajectory)
        species_trajectories = {
            species: trajectory_array[:, idx].copy() for idx, species in enumerate(model.species_names)
        }
        return np.asarray(time_points, dtype=float), species_trajectories


@dataclass
class SimulationResult:
    """Simulation outputs plus steady-state and switching analyses."""

    time_points: np.ndarray
    trajectories: Dict[str, np.ndarray]
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.time_points = np.asarray(self.time_points, dtype=float)
        self.trajectories = {
            species: np.asarray(values, dtype=float) for species, values in self.trajectories.items()
        }
        expected_length = len(self.time_points)
        for species, values in self.trajectories.items():
            if len(values) != expected_length:
                raise ValueError(
                    f"Trajectory length mismatch for {species!r}: expected {expected_length}, got {len(values)}."
                )

    def _steady_slice(self, window: float) -> slice:
        if not 0 < window <= 1:
            raise ValueError("window must be in the interval (0, 1].")
        start = int(np.floor((1.0 - window) * len(self.time_points)))
        return slice(min(start, len(self.time_points) - 1), len(self.time_points))

    def get_steady_state(self, window: float = 0.2) -> Dict[str, float]:
        steady_slice = self._steady_slice(window)
        return {
            species: float(np.mean(values[steady_slice])) for species, values in self.trajectories.items()
        }

    def get_statistics(self, window: float = 0.2) -> Dict[str, Dict[str, float]]:
        steady_slice = self._steady_slice(window)
        statistics: Dict[str, Dict[str, float]] = {}
        for species, values in self.trajectories.items():
            tail = values[steady_slice]
            mean = float(np.mean(tail))
            std = float(np.std(tail, ddof=0))
            cv = float(std / mean) if abs(mean) > 1e-12 else np.nan
            statistics[species] = {"mean": mean, "std": std, "cv": cv}
        return statistics

    def get_switching_times(self) -> Dict[str, Dict[str, Any]]:
        switching: Dict[str, Dict[str, Any]] = {}
        for species, values in self.trajectories.items():
            if len(values) < 3 or np.allclose(values, values[0]):
                continue

            low_threshold = float(np.quantile(values, 0.25))
            high_threshold = float(np.quantile(values, 0.75))
            if high_threshold <= low_threshold:
                continue

            midpoint = 0.5 * (low_threshold + high_threshold)
            current_state = 1 if values[0] >= midpoint else 0
            transitions: List[str] = []
            switching_times: List[float] = []

            for time, value in zip(self.time_points[1:], values[1:]):
                if current_state == 0 and value >= high_threshold:
                    current_state = 1
                    switching_times.append(float(time))
                    transitions.append("low_to_high")
                elif current_state == 1 and value <= low_threshold:
                    current_state = 0
                    switching_times.append(float(time))
                    transitions.append("high_to_low")

            if switching_times:
                switching[species] = {
                    "low_threshold": low_threshold,
                    "high_threshold": high_threshold,
                    "switching_times": switching_times,
                    "transitions": transitions,
                }
        return switching


def run_ensemble(
    simulator: Any,
    model: GeneCircuitModel,
    t_end: float,
    n_runs: int = 100,
    seed: Optional[int] = 42,
) -> Tuple[List[SimulationResult], Dict[str, Dict[str, float]]]:
    if n_runs <= 0:
        raise ValueError("n_runs must be positive.")

    seed_sequence = np.random.SeedSequence(seed)
    child_sequences = seed_sequence.spawn(n_runs)
    results: List[SimulationResult] = []
    steady_states: List[Dict[str, float]] = []
    final_values = {species: [] for species in model.species_names}

    for run_index, child_sequence in enumerate(child_sequences):
        child_seed = int(child_sequence.generate_state(1, dtype=np.uint64)[0])
        time_points, trajectories = simulator.simulate(model, t_end, seed=child_seed)
        result = SimulationResult(
            time_points=time_points,
            trajectories=trajectories,
            metadata={
                "run_index": run_index,
                "seed": child_seed,
                "simulator": simulator.__class__.__name__,
            },
        )
        results.append(result)
        steady_state = result.get_steady_state()
        steady_states.append(steady_state)
        for species in model.species_names:
            final_values[species].append(float(result.trajectories[species][-1]))

    ensemble_statistics: Dict[str, Dict[str, float]] = {
        "steady_state_mean": {},
        "steady_state_std": {},
        "steady_state_cv": {},
        "final_mean": {},
        "final_std": {},
    }
    for species in model.species_names:
        steady_array = np.array([state[species] for state in steady_states], dtype=float)
        final_array = np.array(final_values[species], dtype=float)
        steady_mean = float(np.mean(steady_array))
        steady_std = float(np.std(steady_array, ddof=0))
        ensemble_statistics["steady_state_mean"][species] = steady_mean
        ensemble_statistics["steady_state_std"][species] = steady_std
        ensemble_statistics["steady_state_cv"][species] = (
            float(steady_std / steady_mean) if abs(steady_mean) > 1e-12 else np.nan
        )
        ensemble_statistics["final_mean"][species] = float(np.mean(final_array))
        ensemble_statistics["final_std"][species] = float(np.std(final_array, ddof=0))

    return results, ensemble_statistics


__all__ = [
    "Reaction",
    "GeneCircuitModel",
    "GillespieSimulator",
    "TauLeapingSimulator",
    "SimulationResult",
    "run_ensemble",
]
