"""Utilities for predicting and correcting genetic context effects.

The heuristics in this module are lightweight design-stage proxies inspired by
characterization studies such as Lou et al. (2012) and Mutalik et al. (2013).
They are intended for ranking designs and proposing corrective edits rather than
replacing quantitative experimental calibration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np

Part = Mapping[str, Any]
Design = Dict[str, Any]

COMMON_PROMOTERS: tuple[str, ...] = ("J23100", "J23106", "J23119", "pTet", "pLacUV5")
PROMOTER_ACTIVITY_LOOKUP: Dict[str, float] = {
    "J23100": 1.00,
    "J23106": 0.72,
    "J23119": 0.65,
    "pTet": 0.84,
    "pLacUV5": 0.58,
}
PROMOTER_PROMOTER_INTERFERENCE_MATRIX = np.array(
    [
        [0.93, 0.89, 0.90, 0.88, 0.87],
        [0.90, 0.94, 0.91, 0.89, 0.88],
        [0.91, 0.90, 0.95, 0.90, 0.89],
        [0.88, 0.87, 0.88, 0.94, 0.86],
        [0.87, 0.86, 0.87, 0.85, 0.93],
    ],
    dtype=float,
)

RBS_EFFICIENCY_LOOKUP: Dict[str, float] = {
    "B0030": 0.35,
    "B0031": 0.52,
    "B0032": 0.62,
    "B0034": 1.00,
    "B0035": 0.78,
    "BCD2": 1.15,
}
RBS_CDS_COMPATIBILITY_SCORES: Dict[str, Dict[str, float]] = {
    "B0030": {"gfp": 0.90, "egfp": 0.89, "mcherry": 0.87, "lacz": 0.74, "tetr": 0.83},
    "B0031": {"gfp": 0.93, "egfp": 0.92, "mcherry": 0.90, "lacz": 0.79, "tetr": 0.86},
    "B0032": {"gfp": 0.96, "egfp": 0.95, "mcherry": 0.92, "lacz": 0.82, "tetr": 0.90},
    "B0034": {"gfp": 1.00, "egfp": 0.99, "mcherry": 0.96, "lacz": 0.78, "tetr": 0.94},
    "B0035": {"gfp": 0.98, "egfp": 0.97, "mcherry": 0.94, "lacz": 0.81, "tetr": 0.92},
    "BCD2": {"gfp": 1.04, "egfp": 1.03, "mcherry": 1.01, "lacz": 0.86, "tetr": 0.98},
}

TERMINATOR_READTHROUGH_RATES: Dict[str, float] = {
    "B0015": 0.015,
    "L3S2P21": 0.008,
    "ECK120033736": 0.030,
    "ECK120029600": 0.045,
    "T500": 0.090,
}
TERMINATOR_EFFICIENCY_LOOKUP: Dict[str, float] = {
    name: 1.0 - value for name, value in TERMINATOR_READTHROUGH_RATES.items()
}

INSULATOR_LIBRARY: Dict[str, Dict[str, Any]] = {
    "RiboJ": {
        "name": "RiboJ",
        "type": "insulator",
        "class": "ribozyme",
        "insulation_strength": 0.85,
        "target": "rbs",
    },
    "BydvJ": {
        "name": "BydvJ",
        "type": "insulator",
        "class": "junction_buffer",
        "insulation_strength": 0.72,
        "target": "promoter",
    },
    "RiboJ10": {
        "name": "RiboJ10",
        "type": "insulator",
        "class": "ribozyme",
        "insulation_strength": 0.80,
        "target": "general",
    },
}


def _clamp(value: float, low: float, high: float) -> float:
    return float(np.clip(value, low, high))


def _safe_name(part: Optional[Part]) -> str:
    if not part:
        return ""
    return str(part.get("name", "")).strip()


def _part_type(part: Optional[Part]) -> str:
    if not part:
        return ""
    return str(part.get("type", "")).strip().lower()


def _sequence(part_or_sequence: Optional[Any]) -> str:
    if part_or_sequence is None:
        return ""
    if isinstance(part_or_sequence, str):
        return part_or_sequence.upper()
    if isinstance(part_or_sequence, Mapping):
        return str(part_or_sequence.get("sequence", "")).upper()
    return ""


def _numeric(part: Optional[Part], key: str, default: float) -> float:
    if not part:
        return default
    value = part.get(key, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _gc_content(sequence: str) -> float:
    if not sequence:
        return 0.5
    encoded = np.fromiter(
        (1.0 if base in {"G", "C"} else 0.0 for base in sequence.upper()),
        dtype=float,
    )
    return float(encoded.mean()) if encoded.size else 0.5


def _window_sequence(upstream_sequence: str, rbs_sequence: str, cds_sequence: str, size: int = 30) -> str:
    combined = f"{upstream_sequence[-10:]}{rbs_sequence}{cds_sequence[:10]}"
    if len(combined) >= size:
        start = max((len(combined) - size) // 2, 0)
        return combined[start : start + size]
    return combined.ljust(size, "A")


def _copy_design(design: Mapping[str, Any]) -> Design:
    copied: Design = dict(design)
    copied["parts"] = [dict(part) for part in design.get("parts", [])]
    return copied


def _unique_strings(values: Sequence[str]) -> List[str]:
    seen: set[str] = set()
    unique: List[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            unique.append(value)
    return unique


@dataclass
class GeneCircuitModel:
    """Minimal gene circuit model with named reaction rates."""

    reaction_rates: Dict[str, float]
    species: Dict[str, float] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def copy(self) -> "GeneCircuitModel":
        return GeneCircuitModel(
            reaction_rates=dict(self.reaction_rates),
            species=dict(self.species),
            metadata=dict(self.metadata),
        )


class ContextModel:
    """Predicts context-dependent changes in promoter, RBS, and terminator performance."""

    def __init__(
        self,
        promoter_activity_lookup: Optional[Mapping[str, float]] = None,
        promoter_interference_matrix: Optional[np.ndarray] = None,
        rbs_efficiency_lookup: Optional[Mapping[str, float]] = None,
        rbs_cds_compatibility_scores: Optional[Mapping[str, Mapping[str, float]]] = None,
        terminator_readthrough_rates: Optional[Mapping[str, float]] = None,
    ) -> None:
        self.promoter_activity_lookup = dict(promoter_activity_lookup or PROMOTER_ACTIVITY_LOOKUP)
        self.promoter_interference_matrix = np.array(
            promoter_interference_matrix if promoter_interference_matrix is not None else PROMOTER_PROMOTER_INTERFERENCE_MATRIX,
            dtype=float,
        )
        self.common_promoters = tuple(self.promoter_activity_lookup.keys())
        self.promoter_index = {name: idx for idx, name in enumerate(self.common_promoters)}
        self.rbs_efficiency_lookup = dict(rbs_efficiency_lookup or RBS_EFFICIENCY_LOOKUP)
        self.rbs_cds_compatibility_scores = {
            key: dict(value) for key, value in (rbs_cds_compatibility_scores or RBS_CDS_COMPATIBILITY_SCORES).items()
        }
        self.terminator_readthrough_rates = dict(terminator_readthrough_rates or TERMINATOR_READTHROUGH_RATES)
        self.terminator_efficiency_lookup = {
            name: 1.0 - value for name, value in self.terminator_readthrough_rates.items()
        }

    def _base_promoter_activity(self, promoter: Part) -> float:
        name = _safe_name(promoter)
        return _numeric(promoter, "activity", self.promoter_activity_lookup.get(name, 0.55))

    def _base_rbs_efficiency(self, rbs: Part) -> float:
        name = _safe_name(rbs)
        return _numeric(rbs, "efficiency", self.rbs_efficiency_lookup.get(name, 0.70))

    def _base_terminator_efficiency(self, terminator: Part) -> float:
        name = _safe_name(terminator)
        return _numeric(terminator, "efficiency", self.terminator_efficiency_lookup.get(name, 0.94))

    def _promoter_interference(self, promoter_a: Optional[Part], promoter_b: Optional[Part]) -> float:
        name_a = _safe_name(promoter_a)
        name_b = _safe_name(promoter_b)
        if name_a in self.promoter_index and name_b in self.promoter_index:
            return float(
                self.promoter_interference_matrix[
                    self.promoter_index[name_a], self.promoter_index[name_b]
                ]
            )
        if name_a and name_b:
            return 0.94
        return 1.0

    def predict_promoter_activity(
        self,
        promoter: Part,
        upstream_part: Optional[Part],
        downstream_part: Optional[Part],
    ) -> float:
        """Predict promoter activity after upstream read-through and promoter interference."""

        base_activity = self._base_promoter_activity(promoter)
        upstream_type = _part_type(upstream_part)

        if upstream_type == "terminator":
            upstream_terminator_leakiness = _numeric(
                upstream_part,
                "leakiness",
                self.terminator_readthrough_rates.get(_safe_name(upstream_part), 0.03),
            )
            readthrough_factor = 1.0
        else:
            upstream_terminator_leakiness = 1.0
            readthrough_factor = 0.3 if upstream_part else 0.0

        effective_activity = base_activity * (1.0 - readthrough_factor * upstream_terminator_leakiness)
        if upstream_type == "promoter":
            effective_activity *= self._promoter_interference(upstream_part, promoter)
        if _part_type(downstream_part) == "promoter":
            effective_activity *= self._promoter_interference(promoter, downstream_part)

        return max(0.05 * base_activity, float(effective_activity))

    def predict_rbs_efficiency(self, rbs: Part, upstream_sequence: str, cds: Optional[Part]) -> float:
        """Predict RBS efficiency from local GC content and RBS-CDS compatibility."""

        base_efficiency = self._base_rbs_efficiency(rbs)
        local_window = _window_sequence(upstream_sequence, _sequence(rbs), _sequence(cds))
        gc_content = _gc_content(local_window)
        efficiency_factor = _clamp(1.0 - 0.5 * (gc_content - 0.4), 0.3, 1.0)

        rbs_name = _safe_name(rbs)
        cds_name = _safe_name(cds).lower()
        compatibility_factor = self.rbs_cds_compatibility_scores.get(rbs_name, {}).get(cds_name, 0.90)

        return float(base_efficiency * efficiency_factor * compatibility_factor)

    def predict_terminator_efficiency(self, terminator: Part, downstream_part: Optional[Part]) -> float:
        """Predict context-sensitive terminator efficiency from downstream sequence composition."""

        base_efficiency = self._base_terminator_efficiency(terminator)
        downstream_sequence = _sequence(downstream_part)[:30]
        downstream_gc = _gc_content(downstream_sequence) if downstream_sequence else 0.50
        downstream_factor = _clamp(1.0 - 0.15 * (downstream_gc - 0.45), 0.85, 1.05)
        if _part_type(downstream_part) == "promoter":
            downstream_factor *= 0.96
        return _clamp(base_efficiency * downstream_factor, 0.75, 0.999)


class InsulatorDesigner:
    """Detects strong context interactions and inserts genetic insulators."""

    def __init__(self, context_model: Optional[ContextModel] = None, significant_threshold: float = 0.20) -> None:
        self.context_model = context_model or ContextModel()
        self.significant_threshold = significant_threshold

    def _baseline_for_part(self, part: Part) -> float:
        part_type = _part_type(part)
        if part_type == "promoter":
            return self.context_model._base_promoter_activity(part)
        if part_type == "rbs":
            return self.context_model._base_rbs_efficiency(part)
        if part_type == "terminator":
            return self.context_model._base_terminator_efficiency(part)
        return 1.0

    def _analyze_part(self, parts: Sequence[Part], index: int) -> Optional[Dict[str, Any]]:
        part = parts[index]
        part_type = _part_type(part)
        left = parts[index - 1] if index > 0 else None
        right = parts[index + 1] if index + 1 < len(parts) else None

        if part_type == "promoter":
            predicted = self.context_model.predict_promoter_activity(part, left, right)
            recommendation = "Insert BydvJ-like insulation upstream of promoter."
        elif part_type == "rbs":
            upstream_sequence = _sequence(left)
            cds = right if _part_type(right) == "cds" else None
            predicted = self.context_model.predict_rbs_efficiency(part, upstream_sequence, cds)
            recommendation = "Insert RiboJ upstream of the RBS to standardize the 5' UTR."
        elif part_type == "terminator":
            predicted = self.context_model.predict_terminator_efficiency(part, right)
            recommendation = "Add a downstream buffer or stronger tandem terminator."
        else:
            return None

        baseline = self._baseline_for_part(part)
        deviation = abs(predicted - baseline) / max(baseline, 1e-9)
        return {
            "part_index": index,
            "part_name": _safe_name(part),
            "part_type": part_type,
            "predicted": predicted,
            "baseline": baseline,
            "deviation": float(deviation),
            "recommended_action": recommendation,
        }

    def _choose_insulator(self, affected_part: Part) -> Dict[str, Any]:
        part_type = _part_type(affected_part)
        if part_type == "rbs":
            return dict(INSULATOR_LIBRARY["RiboJ"])
        if part_type == "promoter":
            return dict(INSULATOR_LIBRARY["BydvJ"])
        return dict(INSULATOR_LIBRARY["RiboJ10"])

    def design_insulation(self, circuit_design: Mapping[str, Any]) -> Design:
        """Insert insulators around junctions with >20% predicted deviation."""

        design = _copy_design(circuit_design)
        original_parts = design.get("parts", [])
        modified_parts: List[Dict[str, Any]] = []
        additions: List[Dict[str, Any]] = []
        recommendations: List[str] = []

        for index, part in enumerate(original_parts):
            analysis = self._analyze_part(original_parts, index)
            needs_insulation = analysis is not None and analysis["deviation"] > self.significant_threshold
            previous_is_insulator = bool(modified_parts) and _part_type(modified_parts[-1]) == "insulator"

            if needs_insulation and not previous_is_insulator:
                insulator = self._choose_insulator(part)
                modified_parts.append(insulator)
                additions.append(
                    {
                        "insert_before": _safe_name(part),
                        "insulator": insulator["name"],
                        "predicted_deviation": analysis["deviation"],
                    }
                )
                recommendations.append(
                    f"{insulator['name']} before {_safe_name(part)} ({analysis['deviation']:.1%} deviation)."
                )

            modified_parts.append(dict(part))

        design["parts"] = modified_parts
        design["added_insulators"] = additions
        design["insulation_recommendations"] = recommendations
        return design

    def evaluate_insulation(self, design_with_insulators: Mapping[str, Any]) -> Dict[str, Any]:
        """Predict residual context burden after adding insulators."""

        parts = [dict(part) for part in design_with_insulators.get("parts", [])]
        residuals: List[Dict[str, Any]] = []

        for index, part in enumerate(parts):
            analysis = self._analyze_part(parts, index)
            if analysis is None:
                continue
            attenuation = 0.0
            if index > 0 and _part_type(parts[index - 1]) == "insulator":
                attenuation = max(attenuation, _numeric(parts[index - 1], "insulation_strength", 0.70))
            if index + 1 < len(parts) and _part_type(parts[index + 1]) == "insulator":
                attenuation = max(attenuation, _numeric(parts[index + 1], "insulation_strength", 0.70))
            residual = analysis["deviation"] * (1.0 - attenuation)
            residuals.append(
                {
                    "part_name": analysis["part_name"],
                    "part_type": analysis["part_type"],
                    "residual_deviation": float(residual),
                }
            )

        mean_residual = float(np.mean([entry["residual_deviation"] for entry in residuals])) if residuals else 0.0
        return {
            "insulation_quality_score": _clamp(1.0 - mean_residual, 0.0, 1.0),
            "mean_residual_deviation": mean_residual,
            "residual_context_effects": residuals,
        }


class ContextCorrector:
    """Iteratively retunes parts to compensate for predicted context effects."""

    def __init__(
        self,
        context_model: Optional[ContextModel] = None,
        insulator_designer: Optional[InsulatorDesigner] = None,
        max_iterations: int = 10,
        convergence_threshold: float = 0.05,
    ) -> None:
        self.context_model = context_model or ContextModel()
        self.insulator_designer = insulator_designer or InsulatorDesigner(self.context_model)
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold

    def _estimate_expression(self, design: Mapping[str, Any]) -> float:
        parts = design.get("parts", [])
        expression = 1.0

        for index, part in enumerate(parts):
            left = parts[index - 1] if index > 0 else None
            right = parts[index + 1] if index + 1 < len(parts) else None
            part_type = _part_type(part)

            if part_type == "promoter":
                expression *= self.context_model.predict_promoter_activity(part, left, right)
            elif part_type == "rbs":
                expression *= self.context_model.predict_rbs_efficiency(part, _sequence(left), right if _part_type(right) == "cds" else None)
            elif part_type == "terminator":
                expression *= self.context_model.predict_terminator_efficiency(part, right)

        return float(expression)

    def _select_replacement(self, part: Part, scale: float) -> Dict[str, Any]:
        part_type = _part_type(part)
        updated = dict(part)

        if part_type == "promoter":
            library = self.context_model.promoter_activity_lookup
            current = self.context_model._base_promoter_activity(part)
            key = "activity"
        elif part_type == "rbs":
            library = self.context_model.rbs_efficiency_lookup
            current = self.context_model._base_rbs_efficiency(part)
            key = "efficiency"
        elif part_type == "terminator":
            library = self.context_model.terminator_efficiency_lookup
            current = self.context_model._base_terminator_efficiency(part)
            key = "efficiency"
        else:
            return updated

        target = current * scale
        chosen_name, chosen_value = min(library.items(), key=lambda item: abs(item[1] - target))
        updated["name"] = chosen_name
        updated[key] = float(chosen_value)
        return updated

    def correct_design(self, design: Mapping[str, Any], target_behavior: Mapping[str, float]) -> Design:
        """Adjust parts until predicted behavior converges to the target or max iterations is reached."""

        corrected = _copy_design(design)
        target_expression = float(target_behavior.get("expression", target_behavior.get("target_expression", 1.0)))
        history: List[Dict[str, Any]] = []

        for iteration in range(1, self.max_iterations + 1):
            predicted_expression = self._estimate_expression(corrected)
            deviation = abs(predicted_expression - target_expression) / max(target_expression, 1e-9)
            report = self.generate_context_report(corrected)
            history.append(
                {
                    "iteration": iteration,
                    "predicted_expression": predicted_expression,
                    "target_expression": target_expression,
                    "deviation": deviation,
                }
            )

            if deviation <= self.convergence_threshold:
                break

            predicted_deviations = report["predicted_deviations"]
            if not predicted_deviations:
                break

            dominant = max(predicted_deviations, key=lambda entry: entry["deviation"])
            dominant_index = int(dominant["part_index"])
            scale = _clamp(target_expression / max(predicted_expression, 1e-9), 0.5, 1.5)
            corrected["parts"][dominant_index] = self._select_replacement(corrected["parts"][dominant_index], scale)

            if dominant["deviation"] > 0.20:
                corrected = self.insulator_designer.design_insulation(corrected)

        final_expression = self._estimate_expression(corrected)
        final_deviation = abs(final_expression - target_expression) / max(target_expression, 1e-9)
        corrected["correction_history"] = history
        corrected["predicted_expression"] = final_expression
        corrected["target_expression"] = target_expression
        corrected["final_deviation"] = final_deviation
        corrected["converged"] = final_deviation <= self.convergence_threshold
        return corrected

    def generate_context_report(self, design: Mapping[str, Any]) -> Dict[str, Any]:
        """Return per-junction context diagnostics and actionable recommendations."""

        parts = design.get("parts", [])
        junction_analysis: List[Dict[str, Any]] = []
        deviations: List[Dict[str, Any]] = []
        recommendations: List[str] = []

        for index in range(len(parts) - 1):
            left = parts[index]
            right = parts[index + 1]
            left_type = _part_type(left)
            right_type = _part_type(right)

            if right_type == "promoter":
                baseline = self.context_model._base_promoter_activity(right)
                predicted = self.context_model.predict_promoter_activity(
                    right,
                    left,
                    parts[index + 2] if index + 2 < len(parts) else None,
                )
                part_index = index + 1
                effect_type = "promoter_activity"
                recommendation = "Add BydvJ or strengthen upstream termination."
                part_name = _safe_name(right)
            elif right_type == "rbs":
                baseline = self.context_model._base_rbs_efficiency(right)
                predicted = self.context_model.predict_rbs_efficiency(
                    right,
                    _sequence(left),
                    parts[index + 2] if index + 2 < len(parts) else None,
                )
                part_index = index + 1
                effect_type = "rbs_accessibility"
                recommendation = "Use RiboJ or swap to a more compatible RBS/CDS pair."
                part_name = _safe_name(right)
            elif left_type == "terminator":
                baseline = self.context_model._base_terminator_efficiency(left)
                predicted = self.context_model.predict_terminator_efficiency(left, right)
                part_index = index
                effect_type = "termination"
                recommendation = "Upgrade terminator or insert spacer/insulator downstream."
                part_name = _safe_name(left)
            else:
                continue

            deviation = abs(predicted - baseline) / max(baseline, 1e-9)
            junction_analysis.append(
                {
                    "left_index": index,
                    "right_index": index + 1,
                    "junction": f"{_safe_name(left)}->{_safe_name(right)}",
                    "effect_type": effect_type,
                    "predicted": predicted,
                    "baseline": baseline,
                    "deviation": float(deviation),
                }
            )
            deviations.append(
                {
                    "part_index": part_index,
                    "part_name": part_name,
                    "effect_type": effect_type,
                    "predicted": predicted,
                    "baseline": baseline,
                    "deviation": float(deviation),
                }
            )
            if deviation > self.convergence_threshold:
                recommendations.append(f"{part_name}: {recommendation}")

        mean_deviation = float(np.mean([entry["deviation"] for entry in deviations])) if deviations else 0.0
        context_score = _clamp(1.0 - mean_deviation, 0.0, 1.0)
        return {
            "junction_analysis": junction_analysis,
            "predicted_deviations": sorted(deviations, key=lambda entry: entry["deviation"], reverse=True),
            "recommended_corrections": _unique_strings(recommendations),
            "overall_context_score": context_score,
            "context_score": context_score,
        }


class ContextSimulator:
    """Apply context-dependent rate scaling to a circuit model."""

    def simulate_with_context(
        self,
        model: GeneCircuitModel,
        context_corrections: Mapping[str, Any],
    ) -> GeneCircuitModel:
        adjusted_model = model.copy()
        rate_adjustments = dict(context_corrections.get("rate_adjustments", {}))
        transcription_factor = float(context_corrections.get("transcription_factor", 1.0))
        translation_factor = float(context_corrections.get("translation_factor", 1.0))
        termination_factor = float(context_corrections.get("termination_factor", 1.0))

        for rate_name, base_rate in adjusted_model.reaction_rates.items():
            lower_name = rate_name.lower()
            factor = float(rate_adjustments.get(rate_name, 1.0))

            if "transcription" in lower_name or lower_name.startswith("tx"):
                factor *= transcription_factor
            if "translation" in lower_name or lower_name.startswith("tl"):
                factor *= translation_factor
            if "termination" in lower_name or "readthrough" in lower_name:
                factor *= termination_factor

            adjusted_model.reaction_rates[rate_name] = float(base_rate * factor)

        adjusted_model.metadata["context_corrections"] = dict(context_corrections)
        adjusted_model.metadata["context_applied"] = True
        return adjusted_model


__all__ = [
    "COMMON_PROMOTERS",
    "PROMOTER_ACTIVITY_LOOKUP",
    "PROMOTER_PROMOTER_INTERFERENCE_MATRIX",
    "RBS_CDS_COMPATIBILITY_SCORES",
    "TERMINATOR_READTHROUGH_RATES",
    "ContextCorrector",
    "ContextModel",
    "ContextSimulator",
    "GeneCircuitModel",
    "InsulatorDesigner",
]
