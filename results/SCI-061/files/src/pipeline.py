from __future__ import annotations

from dataclasses import dataclass, field
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    from .circuit_spec import CircuitSpec, Gate
    from .context_effects import ContextCorrector, ContextModel, InsulatorDesigner
    from .parts_catalog import (
        BioPart,
        CDS,
        CircuitAssembler,
        DesignCandidate,
        Insulator,
        PartsCatalog,
        Promoter,
        RBS,
        Terminator,
    )
    from .robust_design import ParameterUncertainty, RobustnessAnalyzer, UncertaintySet
    from .stochastic_sim import (
        GeneCircuitModel as StochasticGeneCircuitModel,
        GillespieSimulator,
        TauLeapingSimulator,
        run_ensemble,
    )
except ImportError:  # pragma: no cover - direct execution support
    from circuit_spec import CircuitSpec, Gate  # type: ignore
    from context_effects import ContextCorrector, ContextModel, InsulatorDesigner  # type: ignore
    from parts_catalog import (  # type: ignore
        BioPart,
        CDS,
        CircuitAssembler,
        DesignCandidate,
        Insulator,
        PartsCatalog,
        Promoter,
        RBS,
        Terminator,
    )
    from robust_design import ParameterUncertainty, RobustnessAnalyzer, UncertaintySet  # type: ignore
    from stochastic_sim import (  # type: ignore
        GeneCircuitModel as StochasticGeneCircuitModel,
        GillespieSimulator,
        TauLeapingSimulator,
        run_ensemble,
    )


@dataclass
class PipelineResult:
    ranked_designs: List[Tuple[DesignCandidate, Dict[str, Any], float, float]]
    best_design: Optional[DesignCandidate]
    summary_statistics: Dict[str, Any]
    execution_time: float
    best_entry: Optional[Tuple[DesignCandidate, Dict[str, Any], float, float]] = None


class AutoDesignPipeline:
    def __init__(self, catalog: Optional[PartsCatalog] = None):
        self.catalog = catalog or PartsCatalog.default_catalog()
        self.context_model = ContextModel()
        self.insulator_designer = InsulatorDesigner(self.context_model)
        self.context_corrector = ContextCorrector(
            context_model=self.context_model,
            insulator_designer=self.insulator_designer,
            max_iterations=6,
            convergence_threshold=0.08,
        )
        self.robustness_analyzer = RobustnessAnalyzer(sobol_base_samples=64, random_seed=17)

    def run(
        self,
        circuit_spec: Any,
        target_behavior: Optional[Mapping[str, Any]] = None,
        n_candidates: int = 10,
        n_sim_runs: int = 50,
        robust_optimization: bool = True,
    ) -> PipelineResult:
        """
        Full pipeline:
        1. Parse/validate circuit specification
        2. Enumerate design candidates from parts catalog
        3. Apply context effect predictions and corrections
        4. Run stochastic simulations for each candidate
        5. Evaluate robustness under parameter uncertainty
        6. Rank designs by combined score (performance * robustness * context_quality)
        7. Return PipelineResult
        """
        started = time.perf_counter()
        spec = self._coerce_circuit_spec(circuit_spec)
        validation = spec.validate_topology(raise_on_error=False)
        if not validation.get("is_valid", False):
            raise ValueError("Invalid circuit specification: " + "; ".join(validation.get("errors", [])))

        gate_requirements = self._build_gate_requirements(spec, target_behavior)
        assembler = CircuitAssembler(gate_requirements, self.catalog)
        candidates = assembler.enumerate_designs(max_designs=max(1, int(n_candidates)))
        if not candidates:
            raise ValueError("No feasible design candidates were produced by the parts catalog.")

        ranked_entries: List[Tuple[DesignCandidate, Dict[str, Any], float, float]] = []
        combined_scores: List[float] = []
        performance_scores: List[float] = []
        robustness_scores: List[float] = []
        context_scores: List[float] = []

        for candidate in candidates:
            context_summary = self.evaluate_context(candidate, target_behavior=target_behavior)
            sim_results = self.simulate_design(
                spec,
                candidate,
                n_runs=max(1, int(n_sim_runs)),
                t_end=self._screening_t_end(spec, target_behavior),
                seed=101,
            )
            performance_score = float(
                np.clip(0.5 * float(candidate.score) + 0.5 * float(sim_results.get("performance_score", 0.0)), 0.0, 1.0)
            )
            if robust_optimization:
                robustness_summary = self.evaluate_robustness(
                    spec,
                    candidate,
                    sim_results=sim_results,
                    target_behavior=target_behavior,
                )
                robustness_score = float(robustness_summary["robustness_score"])
            else:
                robustness_summary = {"robustness_score": 1.0, "uncertainty_parameters": []}
                robustness_score = 1.0
            context_score = float(context_summary["context_score"])
            combined_score = float(performance_score * robustness_score * context_score)

            sim_results["performance_score"] = performance_score
            sim_results["context_summary"] = context_summary
            sim_results["robustness_summary"] = robustness_summary
            sim_results["combined_score"] = combined_score

            ranked_entries.append((candidate, sim_results, robustness_score, context_score))
            combined_scores.append(combined_score)
            performance_scores.append(performance_score)
            robustness_scores.append(robustness_score)
            context_scores.append(context_score)

        ranked_entries.sort(key=lambda entry: float(entry[1].get("combined_score", 0.0)), reverse=True)
        best_entry = ranked_entries[0] if ranked_entries else None
        execution_time = time.perf_counter() - started
        summary_statistics = {
            "circuit_name": spec.name,
            "n_candidates_requested": int(n_candidates),
            "n_candidates_evaluated": len(ranked_entries),
            "validation": validation,
            "combined_score": self._summary_from_values(combined_scores),
            "performance_score": self._summary_from_values(performance_scores),
            "robustness_score": self._summary_from_values(robustness_scores),
            "context_score": self._summary_from_values(context_scores),
            "best_combined_score": float(best_entry[1].get("combined_score", 0.0)) if best_entry else 0.0,
            "best_design_length": len(best_entry[0].to_dna_sequence()) if best_entry else 0,
            "sbol_summary": spec.to_sbol_like_dict(),
        }
        return PipelineResult(
            ranked_designs=ranked_entries,
            best_design=best_entry[0] if best_entry else None,
            summary_statistics=summary_statistics,
            execution_time=float(execution_time),
            best_entry=best_entry,
        )

    def run_from_verilog(self, verilog_text: str, **kwargs: Any) -> PipelineResult:
        """Parse Verilog-like spec and run pipeline"""
        return self.run(CircuitSpec.from_verilog_like(verilog_text), **kwargs)

    def evaluate_context(
        self,
        design_candidate: DesignCandidate,
        target_behavior: Optional[Mapping[str, Any]] = None,
    ) -> Dict[str, Any]:
        design = self._candidate_to_design_dict(design_candidate)
        initial_report = self.context_corrector.generate_context_report(design)
        corrected_design = self.context_corrector.correct_design(
            design,
            self._target_expression_behavior(target_behavior),
        )
        corrected_report = self.context_corrector.generate_context_report(corrected_design)
        insulation_report = self.insulator_designer.evaluate_insulation(corrected_design)
        context_score = float(
            np.clip(
                0.5 * float(corrected_report.get("context_score", 0.0))
                + 0.5 * float(insulation_report.get("insulation_quality_score", 0.0)),
                0.0,
                1.0,
            )
        )
        return {
            "design": design,
            "initial_report": initial_report,
            "corrected_design": corrected_design,
            "corrected_report": corrected_report,
            "insulation_report": insulation_report,
            "context_score": context_score,
        }

    def build_simulation_model(self, circuit_spec: Any, design_candidate: DesignCandidate) -> StochasticGeneCircuitModel:
        spec = self._coerce_circuit_spec(circuit_spec)
        gates = self._spec_gate_sequence(spec)
        bundles = self._candidate_gate_bundles(design_candidate, len(gates))
        if not bundles:
            raise ValueError("Design candidate does not contain any recognizable transcription-unit bundles.")

        species: Dict[str, int] = {
            signal.name: int(max(0.0, float(signal.initial_value)))
            for signal in spec.signals.values()
            if signal.signal_type == "input"
        }
        genes: List[Dict[str, Any]] = []
        internal_outputs = {gate.output for gate in gates}

        for index, gate in enumerate(gates):
            bundle = bundles[min(index, len(bundles) - 1)]
            context_factors = self._bundle_context_factors(bundle, bundles[index + 1] if index + 1 < len(bundles) else None)
            promoter = bundle[1]
            rbs = bundle[2]
            cds = bundle[3]
            terminator = bundle[4]

            regulators: List[Dict[str, Any]] = []
            seen_regulators: set[str] = set()
            for regulator_name in spec.effective_inputs(gate.name):
                if regulator_name in seen_regulators or regulator_name not in internal_outputs:
                    continue
                seen_regulators.add(regulator_name)
                regulators.append(
                    {
                        "protein": regulator_name,
                        "mode": self._gate_regulation_mode(gate),
                        "Kd": max(5.0, 35.0 / max(context_factors["transcription_factor"], 0.2)),
                        "hill": float(gate.resolved_parameters["n"]),
                        "binding_rate": 0.01,
                        "unbinding_rate": 0.12,
                    }
                )

            genes.append(
                {
                    "name": gate.output,
                    "protein": gate.output,
                    "promoter": f"{gate.name}_{promoter.name}",
                    "transcription_rate": max(0.1, 4.0 * float(promoter.parameters.get("rpu", 1.0)) * context_factors["transcription_factor"]),
                    "rbs": f"{gate.name}_{rbs.name}",
                    "translation_rate": max(0.05, 2.0 * float(rbs.parameters.get("translation_rate", 1.0)) * context_factors["translation_factor"]),
                    "protein_degradation": max(
                        0.005,
                        float(cds.parameters.get("degradation_rate", 0.04)) + 0.02 * (1.0 - context_factors["termination_factor"]),
                    ),
                    "copy_number": 1,
                    "initial_protein": int(spec.signals.get(gate.output, spec.add_signal(gate.output)).initial_value),
                    "regulators": regulators,
                    "metadata": {
                        "selected_cds": cds.name,
                        "selected_terminator": terminator.name,
                        "context_factors": context_factors,
                    },
                }
            )

        return StochasticGeneCircuitModel.from_circuit_spec(
            {
                "species": species,
                "defaults": {
                    "transcription": 6.0,
                    "translation": 2.0,
                    "mrna_degradation": 0.25,
                    "protein_degradation": 0.04,
                    "binding_rate": 0.01,
                    "unbinding_rate": 0.12,
                    "hill_coefficient": 2.0,
                    "dissociation_constant": 35.0,
                    "activation_fold": 3.0,
                    "leakiness": 0.02,
                },
                "genes": genes,
            }
        )

    def simulate_design(
        self,
        circuit_spec: Any,
        design_candidate: DesignCandidate,
        n_runs: int = 10,
        t_end: Optional[float] = None,
        simulator: Optional[Any] = None,
        seed: Optional[int] = 42,
    ) -> Dict[str, Any]:
        spec = self._coerce_circuit_spec(circuit_spec)
        behavior_mode = self._behavior_mode(spec, None)
        model = self.build_simulation_model(spec, design_candidate)
        if t_end is None:
            t_end = self._screening_t_end(spec, None)
        if simulator is None:
            simulator = GillespieSimulator(max_steps=50000) if behavior_mode == "toggle" else TauLeapingSimulator(max_steps=50000)

        try:
            results, ensemble_stats = run_ensemble(simulator, model, float(t_end), n_runs=max(1, int(n_runs)), seed=seed)
        except RuntimeError:
            fallback = TauLeapingSimulator(max_steps=50000)
            results, ensemble_stats = run_ensemble(fallback, model, float(t_end), n_runs=max(1, int(n_runs)), seed=seed)
            simulator = fallback

        dynamic_metrics = self._dynamic_metrics(spec, results)
        return {
            "results": results,
            "ensemble_statistics": ensemble_stats,
            "dynamic_metrics": dynamic_metrics,
            "performance_score": float(dynamic_metrics.get("performance_score", 0.0)),
            "simulator": simulator.__class__.__name__,
            "t_end": float(t_end),
            "n_runs": int(n_runs),
            "output_species": self._primary_output_species(spec),
        }

    def evaluate_robustness(
        self,
        circuit_spec: Any,
        design_candidate: DesignCandidate,
        sim_results: Optional[Mapping[str, Any]] = None,
        target_behavior: Optional[Mapping[str, Any]] = None,
    ) -> Dict[str, Any]:
        spec = self._coerce_circuit_spec(circuit_spec)
        model = self.build_simulation_model(spec, design_candidate)
        nominal = self._nominal_uncertainty_parameters(model)
        uncertainties = UncertaintySet(
            [
                ParameterUncertainty(name, value, "lognormal", {"mean": value, "cv": 0.25})
                for name, value in nominal.items()
            ]
        )
        threshold = target_behavior.get("robustness_threshold", {"relative": 0.35}) if target_behavior else {"relative": 0.35}

        def proxy_model(params: Mapping[str, float]) -> float:
            tx = float(params["transcription_rate"])
            tl = float(params["translation_rate"])
            mrna_deg = float(params["mrna_degradation_rate"])
            protein_deg = float(params["protein_degradation_rate"])
            binding = float(params["binding_affinity"])
            hill = float(params["hill_coefficient"])
            leakage = float(params["leakage_rate"])
            repression = binding / (binding + 10.0)
            effective_gain = (tx * tl) / ((1.0 + 8.0 * mrna_deg) * (1.0 + 15.0 * protein_deg))
            nonlinear_bonus = 1.0 + 0.08 * max(hill - 1.0, 0.0)
            leakage_penalty = max(0.2, 1.0 - 2.0 * leakage)
            species_scale = max(len(model.species_names), 1)
            return float((effective_gain * repression * nonlinear_bonus * leakage_penalty) / species_scale)

        robustness_score = self.robustness_analyzer.robustness_score(
            proxy_model,
            uncertainties,
            output_species=None,
            threshold=threshold,
            n_samples=80,
        )
        sensitivity = self.robustness_analyzer.sensitivity_analysis(
            proxy_model,
            uncertainties,
            output_species=None,
            method="local",
        )
        return {
            "robustness_score": float(robustness_score),
            "sensitivity": sensitivity,
            "uncertainty_parameters": nominal,
            "reference_output_species": self._primary_output_species(spec),
            "simulation_context": dict(sim_results or {}),
        }

    def _coerce_circuit_spec(self, circuit_spec: Any) -> CircuitSpec:
        if isinstance(circuit_spec, CircuitSpec):
            return circuit_spec
        if isinstance(circuit_spec, Mapping):
            if "name" in circuit_spec and "gates" in circuit_spec:
                return CircuitSpec.from_dict(circuit_spec)
        raise TypeError("circuit_spec must be a CircuitSpec instance or CircuitSpec-compatible mapping.")

    def _build_gate_requirements(
        self,
        spec: CircuitSpec,
        target_behavior: Optional[Mapping[str, Any]],
    ) -> List[Dict[str, Any]]:
        mode = self._behavior_mode(spec, target_behavior)
        requirements: List[Dict[str, Any]] = []
        for gate in self._spec_gate_sequence(spec):
            gate_type = gate.gate_type.upper()
            feedback_sources = spec.feedback_sources_for_gate(gate.name)
            effective_inputs = [name for name in spec.effective_inputs(gate.name) if name not in spec.signals or spec.signals[name].signal_type != "input"]
            regulator = feedback_sources[0] if feedback_sources else (effective_inputs[0] if effective_inputs else None)
            is_reporter = spec.signals.get(gate.output) is not None and spec.signals[gate.output].signal_type == "output"
            target_rpu = 2.8 if gate_type in {"NOT", "NOR", "NAND"} else 1.8
            if mode == "repressilator" and not is_reporter:
                target_rpu = 2.4
            target_rbs = 0.8 if is_reporter else 0.55
            requirements.append(
                {
                    "name": gate.name,
                    "output": gate.output,
                    "regulator": regulator,
                    "target_rpu": target_rpu,
                    "target_rbs_strength": target_rbs,
                    "min_termination_efficiency": 0.9,
                    "min_insulation_score": 0.85,
                    "context": self._compatibility_family(regulator or gate.output),
                    "compatibility_family": self._compatibility_family(regulator or gate.output),
                }
            )
        return requirements

    def _spec_gate_sequence(self, circuit_spec: CircuitSpec) -> List[Gate]:
        return list(circuit_spec.gates.values())

    def _candidate_gate_bundles(self, candidate: DesignCandidate, n_gates: int) -> List[Tuple[BioPart, BioPart, BioPart, BioPart, BioPart]]:
        parts = list(candidate.parts)
        bundle_size = 5
        if len(parts) < bundle_size:
            return []
        bundles: List[Tuple[BioPart, BioPart, BioPart, BioPart, BioPart]] = []
        for index in range(0, min(len(parts), n_gates * bundle_size), bundle_size):
            chunk = parts[index : index + bundle_size]
            if len(chunk) == bundle_size:
                bundles.append(tuple(chunk))  # type: ignore[arg-type]
        return bundles

    def _candidate_to_design_dict(self, candidate: DesignCandidate) -> Dict[str, Any]:
        return {
            "parts": [self._part_to_context_dict(part) for part in candidate.parts],
            "score": float(candidate.score),
            "predicted_performance": dict(candidate.predicted_performance),
        }

    def _part_to_context_dict(self, part: BioPart) -> Dict[str, Any]:
        payload = {
            "name": part.name,
            "type": part.part_type,
            "sequence": part.sequence,
            "parameters": dict(part.parameters),
        }
        payload.update(dict(part.parameters))
        if part.part_type == "promoter":
            payload.setdefault("activity", float(part.parameters.get("rpu", 0.55)))
        elif part.part_type == "rbs":
            payload.setdefault("efficiency", float(part.parameters.get("translation_rate", 0.7)))
        elif part.part_type == "terminator":
            efficiency = float(part.parameters.get("termination_efficiency", 0.94))
            payload.setdefault("efficiency", efficiency)
            payload.setdefault("leakiness", 1.0 - efficiency)
        elif part.part_type == "insulator":
            payload.setdefault("insulation_strength", float(part.parameters.get("insulation_score", 0.8)))
        return payload

    def _bundle_context_factors(
        self,
        bundle: Sequence[BioPart],
        next_bundle: Optional[Sequence[BioPart]] = None,
    ) -> Dict[str, float]:
        promoter = self._part_to_context_dict(bundle[1])
        rbs = self._part_to_context_dict(bundle[2])
        cds = self._part_to_context_dict(bundle[3])
        terminator = self._part_to_context_dict(bundle[4])
        upstream = self._part_to_context_dict(bundle[0])
        downstream = self._part_to_context_dict(next_bundle[0]) if next_bundle else None

        promoter_base = self.context_model._base_promoter_activity(promoter)
        promoter_pred = self.context_model.predict_promoter_activity(promoter, upstream, rbs)
        rbs_base = self.context_model._base_rbs_efficiency(rbs)
        rbs_pred = self.context_model.predict_rbs_efficiency(rbs, promoter.get("sequence", ""), cds)
        terminator_base = self.context_model._base_terminator_efficiency(terminator)
        terminator_pred = self.context_model.predict_terminator_efficiency(terminator, downstream)
        return {
            "transcription_factor": float(np.clip(promoter_pred / max(promoter_base, 1e-9), 0.3, 1.5)),
            "translation_factor": float(np.clip(rbs_pred / max(rbs_base, 1e-9), 0.3, 1.5)),
            "termination_factor": float(np.clip(terminator_pred / max(terminator_base, 1e-9), 0.7, 1.1)),
        }

    def _dynamic_metrics(self, spec: CircuitSpec, results: Sequence[Any]) -> Dict[str, Any]:
        mode = self._behavior_mode(spec, None)
        if mode == "toggle":
            return self._toggle_metrics(spec, results)
        if mode == "repressilator":
            return self._repressilator_metrics(spec, results)
        return self._generic_metrics(spec, results)

    def _toggle_metrics(self, spec: CircuitSpec, results: Sequence[Any]) -> Dict[str, Any]:
        regulators = list(dict.fromkeys(loop.source_signal for loop in spec.feedback_loops))
        if len(regulators) < 2:
            return self._generic_metrics(spec, results)
        species_a, species_b = regulators[:2]
        final_pairs = np.array(
            [[float(run.trajectories[species_a][-1]), float(run.trajectories[species_b][-1])] for run in results],
            dtype=float,
        )
        diff = final_pairs[:, 0] - final_pairs[:, 1]
        total = np.maximum(final_pairs.sum(axis=1), 1.0)
        threshold = max(5.0, 0.2 * float(np.mean(total)))
        state_a = diff > threshold
        state_b = diff < -threshold
        separation = abs(float(np.mean(diff[state_a])) - float(np.mean(diff[state_b]))) / max(float(np.mean(total)), 1.0) if np.any(state_a) and np.any(state_b) else abs(float(np.mean(diff))) / max(float(np.mean(total)), 1.0)
        occupancy = 2.0 * min(float(np.mean(state_a)), float(np.mean(state_b))) if (np.any(state_a) and np.any(state_b)) else 0.0
        bistability_index = float(np.clip(0.6 * separation + 0.4 * occupancy, 0.0, 1.0))
        return {
            "mode": "toggle",
            "species": [species_a, species_b],
            "bistability_index": bistability_index,
            "state_a_fraction": float(np.mean(state_a)),
            "state_b_fraction": float(np.mean(state_b)),
            "performance_score": bistability_index,
        }

    def _repressilator_metrics(self, spec: CircuitSpec, results: Sequence[Any]) -> Dict[str, Any]:
        internal_species = [signal.name for signal in spec.signals.values() if signal.signal_type == "internal"]
        if not internal_species:
            return self._generic_metrics(spec, results)
        amplitudes: List[float] = []
        periods: List[float] = []
        zero_crossing_counts: List[int] = []
        probe = internal_species[0]
        for run in results:
            grid = np.linspace(float(run.time_points[0]), float(run.time_points[-1]), 400)
            signal = np.interp(grid, run.time_points, run.trajectories[probe])
            centered = signal - float(np.mean(signal))
            amplitudes.append(0.5 * float(np.quantile(signal, 0.95) - np.quantile(signal, 0.05)))
            crossings = np.where(np.diff(np.signbit(centered)))[0]
            zero_crossing_counts.append(int(len(crossings)))
            if len(crossings) >= 3:
                half_periods = np.diff(grid[crossings])
                periods.append(2.0 * float(np.mean(half_periods)))
        mean_amplitude = float(np.mean(amplitudes)) if amplitudes else 0.0
        mean_count = float(np.mean(zero_crossing_counts)) if zero_crossing_counts else 0.0
        normalized_amplitude = mean_amplitude / max(mean_amplitude + 5.0, 1.0)
        oscillation_score = float(np.clip(0.55 * normalized_amplitude + 0.45 * min(mean_count / 6.0, 1.0), 0.0, 1.0))
        return {
            "mode": "repressilator",
            "species": internal_species,
            "mean_amplitude": mean_amplitude,
            "mean_period": float(np.mean(periods)) if periods else float("nan"),
            "mean_zero_crossings": mean_count,
            "performance_score": oscillation_score,
        }

    def _generic_metrics(self, spec: CircuitSpec, results: Sequence[Any]) -> Dict[str, Any]:
        probe = self._primary_output_species(spec)
        if probe not in results[0].trajectories:
            probe = next(iter(results[0].trajectories))
        steady = np.array([run.get_steady_state().get(probe, 0.0) for run in results], dtype=float)
        mean_value = float(np.mean(steady))
        std_value = float(np.std(steady, ddof=0))
        cv = float(std_value / mean_value) if abs(mean_value) > 1e-12 else 1.0
        performance = float(np.clip(mean_value / (mean_value + 10.0), 0.0, 1.0) * np.clip(1.0 - min(cv, 1.0), 0.0, 1.0))
        return {
            "mode": "generic",
            "species": [probe],
            "steady_state_mean": mean_value,
            "steady_state_cv": cv,
            "performance_score": performance,
        }

    def _nominal_uncertainty_parameters(self, model: StochasticGeneCircuitModel) -> Dict[str, float]:
        genes = list(model.metadata.get("circuit_spec", {}).get("genes", []))
        transcription = [float(gene.get("transcription_rate", 1.0)) for gene in genes] or [1.0]
        translation = [float(gene.get("translation_rate", 1.0)) for gene in genes] or [1.0]
        protein_deg = [float(gene.get("protein_degradation", 0.04)) for gene in genes] or [0.04]
        mrna_deg = [float(model.metadata.get("circuit_spec", {}).get("defaults", {}).get("mrna_degradation", 0.25))]
        binding_terms = [float(reg.get("Kd", 35.0)) for gene in genes for reg in gene.get("regulators", [])] or [35.0]
        hill_terms = [float(reg.get("hill", 2.0)) for gene in genes for reg in gene.get("regulators", [])] or [2.0]
        leakage_terms = [float(model.metadata.get("circuit_spec", {}).get("defaults", {}).get("leakiness", 0.02))]
        return {
            "transcription_rate": float(np.mean(transcription)),
            "translation_rate": float(np.mean(translation)),
            "mrna_degradation_rate": float(np.mean(mrna_deg)),
            "protein_degradation_rate": float(np.mean(protein_deg)),
            "binding_affinity": float(np.mean(binding_terms)),
            "hill_coefficient": float(np.mean(hill_terms)),
            "leakage_rate": float(np.mean(leakage_terms)),
        }

    def _summary_from_values(self, values: Sequence[float]) -> Dict[str, float]:
        array = np.asarray(list(values), dtype=float)
        if array.size == 0:
            return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0}
        return {
            "mean": float(np.mean(array)),
            "std": float(np.std(array, ddof=0)),
            "min": float(np.min(array)),
            "max": float(np.max(array)),
        }

    def _screening_t_end(self, spec: CircuitSpec, target_behavior: Optional[Mapping[str, Any]]) -> float:
        mode = self._behavior_mode(spec, target_behavior)
        if mode == "toggle":
            return 60.0
        if mode == "repressilator":
            return 120.0
        return 80.0

    def _behavior_mode(self, spec: CircuitSpec, target_behavior: Optional[Mapping[str, Any]]) -> str:
        if target_behavior:
            for key in ("mode", "behavior", "target"):
                if key in target_behavior:
                    return str(target_behavior[key]).lower()
        name = spec.name.lower()
        if "toggle" in name:
            return "toggle"
        if "repress" in name or len(spec.feedback_loops) >= 3:
            return "repressilator"
        return "generic"

    def _compatibility_family(self, name: str) -> str:
        lowered = str(name).lower()
        if "lux" in lowered:
            return "lux"
        if lowered == "ci":
            return "lambda"
        return "sigma70"

    def _gate_regulation_mode(self, gate: Gate) -> str:
        if gate.gate_type.upper() in {"NOT", "NAND", "NOR"}:
            return "repression"
        return "activation"

    def _primary_output_species(self, spec: CircuitSpec) -> str:
        outputs = [signal.name for signal in spec.signals.values() if signal.signal_type == "output"]
        if outputs:
            return outputs[0]
        return next(iter(spec.signals))

    def _target_expression_behavior(self, target_behavior: Optional[Mapping[str, Any]]) -> Dict[str, float]:
        if not target_behavior:
            return {"expression": 1.0}
        if "expression" in target_behavior or "target_expression" in target_behavior:
            return {
                "expression": float(target_behavior.get("expression", target_behavior.get("target_expression", 1.0)))
            }
        return {"expression": 1.0}


__all__ = ["AutoDesignPipeline", "PipelineResult"]
