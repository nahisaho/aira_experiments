from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import product
from typing import Any, Dict, List, Optional

import numpy as np

try:
    from .circuit_spec import CircuitSpec
except ImportError:  # pragma: no cover - supports direct import
    try:
        from circuit_spec import CircuitSpec  # type: ignore
    except ImportError:  # pragma: no cover - fallback for loosely coupled usage
        CircuitSpec = Any  # type: ignore


DNA_ALPHABET = frozenset({"A", "C", "G", "T"})
PART_TYPES = ("promoter", "rbs", "cds", "terminator", "insulator")
SBOL_ROLE_MAP: Dict[str, str] = {
    "promoter": "http://identifiers.org/so/SO:0000167",
    "rbs": "http://identifiers.org/so/SO:0000139",
    "cds": "http://identifiers.org/so/SO:0000316",
    "terminator": "http://identifiers.org/so/SO:0000141",
    "insulator": "http://identifiers.org/so/SO:0000627",
}


@dataclass(frozen=True)
class BioPart:
    name: str
    sequence: str
    parameters: Dict[str, float] = field(default_factory=dict)
    part_type: str = "part"
    sbol_role: str = "http://identifiers.org/so/SO:0000110"

    def __post_init__(self) -> None:
        cleaned = self.sequence.upper().replace(" ", "")
        if not cleaned:
            raise ValueError(f"Part {self.name} must have a non-empty DNA sequence")
        if not set(cleaned).issubset(DNA_ALPHABET):
            raise ValueError(f"Part {self.name} contains non-DNA characters")
        object.__setattr__(self, "sequence", cleaned)
        object.__setattr__(self, "parameters", dict(self.parameters))


@dataclass(frozen=True)
class Promoter(BioPart):
    part_type: str = "promoter"
    sbol_role: str = SBOL_ROLE_MAP["promoter"]


@dataclass(frozen=True)
class RBS(BioPart):
    part_type: str = "rbs"
    sbol_role: str = SBOL_ROLE_MAP["rbs"]


@dataclass(frozen=True)
class CDS(BioPart):
    part_type: str = "cds"
    sbol_role: str = SBOL_ROLE_MAP["cds"]


@dataclass(frozen=True)
class Terminator(BioPart):
    part_type: str = "terminator"
    sbol_role: str = SBOL_ROLE_MAP["terminator"]


@dataclass(frozen=True)
class Insulator(BioPart):
    part_type: str = "insulator"
    sbol_role: str = SBOL_ROLE_MAP["insulator"]


@dataclass
class DesignCandidate:
    parts: List[BioPart]
    score: float
    predicted_performance: Dict[str, Any]

    def to_dna_sequence(self) -> str:
        return "".join(part.sequence for part in self.parts)

    def to_sbol_dict(self) -> Dict[str, Any]:
        components: List[Dict[str, Any]] = []
        start = 1
        for index, part in enumerate(self.parts, start=1):
            end = start + len(part.sequence) - 1
            components.append(
                {
                    "displayId": f"{part.part_type}_{part.name}_{index}",
                    "name": part.name,
                    "type": part.part_type,
                    "role": part.sbol_role,
                    "sequence": part.sequence,
                    "start": start,
                    "end": end,
                    "parameters": dict(part.parameters),
                }
            )
            start = end + 1
        return {
            "type": "ComponentDefinition",
            "displayId": "assembled_design",
            "roles": ["http://identifiers.org/so/SO:0000804"],
            "sequence": self.to_dna_sequence(),
            "components": components,
            "score": self.score,
            "predicted_performance": dict(self.predicted_performance),
        }


class PartsCatalog:
    def __init__(self, parts: Optional[Iterable[BioPart]] = None, load_default: bool = True) -> None:
        self.library: Dict[str, List[BioPart]] = {part_type: [] for part_type in PART_TYPES}
        if load_default:
            for part in self._default_parts():
                self.add_part(part)
        if parts is not None:
            for part in parts:
                self.add_part(part)

    def add_part(self, part: BioPart) -> None:
        if part.part_type not in self.library:
            self.library[part.part_type] = []
        self.library[part.part_type].append(part)

    def get_parts(self, part_type: str) -> List[BioPart]:
        return list(self.library.get(part_type, []))

    def query(self, part_type: str, constraints: Optional[Mapping[str, Any]] = None) -> List[BioPart]:
        parts = self.get_parts(part_type)
        if not constraints:
            return parts
        return [part for part in parts if self._matches_constraints(part, constraints)]

    def get_compatible_parts(self, gate: Any) -> Dict[str, List[BioPart]]:
        gate_dict = self._normalize_gate(gate)
        compatible: Dict[str, List[BioPart]] = {}

        promoter_constraints: Dict[str, Any] = {}
        if "target_rpu" in gate_dict:
            promoter_constraints["rpu"] = ("approx", gate_dict["target_rpu"])
        if "regulator" in gate_dict:
            promoter_constraints["regulated_by"] = gate_dict["regulator"]
        compatible["promoter"] = self.query("promoter", promoter_constraints)

        rbs_constraints: Dict[str, Any] = {}
        if "target_rbs_strength" in gate_dict:
            rbs_constraints["translation_rate"] = ("approx", gate_dict["target_rbs_strength"])
        compatible["rbs"] = self.query("rbs", rbs_constraints)

        cds_constraints: Dict[str, Any] = {}
        output_name = gate_dict.get("output") or gate_dict.get("cds")
        if output_name is not None:
            cds_constraints["name"] = output_name
        compatible["cds"] = self.query("cds", cds_constraints)

        terminator_constraints: Dict[str, Any] = {}
        if "min_termination_efficiency" in gate_dict:
            terminator_constraints["termination_efficiency"] = (
                ">=",
                gate_dict["min_termination_efficiency"],
            )
        compatible["terminator"] = self.query("terminator", terminator_constraints)

        insulator_constraints: Dict[str, Any] = {}
        if "min_insulation_score" in gate_dict:
            insulator_constraints["insulation_score"] = (">=", gate_dict["min_insulation_score"])
        compatible["insulator"] = self.query("insulator", insulator_constraints)

        family = gate_dict.get("context") or gate_dict.get("compatibility_family")
        if family:
            compatible["promoter"] = [
                part
                for part in compatible["promoter"]
                if part.parameters.get("compatibility_family", family) == family
            ]
            compatible["rbs"] = [
                part
                for part in compatible["rbs"]
                if part.parameters.get("compatibility_family", family) == family
            ]

        return compatible

    @classmethod
    def default_catalog(cls) -> "PartsCatalog":
        return cls(load_default=True)

    @staticmethod
    def _matches_constraints(part: BioPart, constraints: Mapping[str, Any]) -> bool:
        for key, expected in constraints.items():
            actual = part.name if key == "name" else part.parameters.get(key)
            if isinstance(expected, tuple) and len(expected) == 2:
                op, value = expected
                if actual is None:
                    return False
                if op == ">=" and not actual >= value:
                    return False
                if op == "<=" and not actual <= value:
                    return False
                if op == ">" and not actual > value:
                    return False
                if op == "<" and not actual < value:
                    return False
                if op == "==" and not actual == value:
                    return False
                if op == "approx":
                    tolerance = max(abs(float(value)) * 0.35, 1e-9)
                    if abs(float(actual) - float(value)) > tolerance:
                        return False
                continue
            if callable(expected):
                if not expected(actual):
                    return False
            elif actual != expected:
                return False
        return True

    @staticmethod
    def _normalize_gate(gate: Any) -> Dict[str, Any]:
        if gate is None:
            return {}
        if isinstance(gate, Mapping):
            return dict(gate)
        gate_dict: Dict[str, Any] = {}
        for attr in dir(gate):
            if attr.startswith("_"):
                continue
            value = getattr(gate, attr)
            if callable(value):
                continue
            gate_dict[attr] = value
        return gate_dict

    @staticmethod
    def _default_parts() -> List[BioPart]:
        return [
            Promoter("pTac", "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA", {"rpu": 4.2, "regulated_by": "LacI", "compatibility_family": "sigma70", "basal_rpu": 0.3}),
            Promoter("pBAD", "TTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGG", {"rpu": 3.1, "regulated_by": "AraC", "compatibility_family": "sigma70", "basal_rpu": 0.15}),
            Promoter("pTet", "TCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGA", {"rpu": 2.7, "regulated_by": "TetR", "compatibility_family": "sigma70", "basal_rpu": 0.1}),
            Promoter("pLac", "AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATGAC", {"rpu": 1.8, "regulated_by": "LacI", "compatibility_family": "sigma70", "basal_rpu": 0.2}),
            Promoter("pLuxR", "ATAAATTCCTGTGTGAAATTGTTATCCGCTCACAATTCCAC", {"rpu": 2.2, "regulated_by": "LuxR", "compatibility_family": "lux", "basal_rpu": 0.12}),
            Promoter("pCI", "TTGACACTATCGTATAATGTGTGGATTATATCACCGCCAGAG", {"rpu": 1.4, "regulated_by": "CI", "compatibility_family": "lambda", "basal_rpu": 0.05}),
            RBS("B0030", "AAAGAGGAGAAA", {"translation_rate": 0.3, "compatibility_family": "sigma70"}),
            RBS("B0031", "AAAGAAGGAGAT", {"translation_rate": 0.45, "compatibility_family": "sigma70"}),
            RBS("B0032", "AAAGAGGTAGAC", {"translation_rate": 0.7, "compatibility_family": "lux"}),
            RBS("B0034", "AAAGAGGAGAAA", {"translation_rate": 1.0, "compatibility_family": "sigma70"}),
            CDS("GFP", "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAA", {"degradation_rate": 0.02, "output_role": "reporter", "preferred_context": "sigma70"}),
            CDS("RFP", "ATGGCCTCCTCCGAGGACGTCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAG", {"degradation_rate": 0.025, "output_role": "reporter", "preferred_context": "sigma70"}),
            CDS("BFP", "ATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGAACAACCACCACTTCAAGTGCACATCCGAGGACGGCAACATCCTGGGGCACAAGCTG", {"degradation_rate": 0.03, "output_role": "reporter", "preferred_context": "sigma70"}),
            CDS("LacI", "ATGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGAAACTCTTCAAGCGTTTCTCGCACGAGATGGTTTCGACGATGCCCTTG", {"degradation_rate": 0.05, "output_role": "repressor", "preferred_context": "sigma70"}),
            CDS("TetR", "ATGTCTAGATTAGATAAAAGTAAAGTGATTAACAGCGCATTAGAGCTGCTTAATGAGGTCGGAATCGAAATCAGGCTGATC", {"degradation_rate": 0.06, "output_role": "repressor", "preferred_context": "sigma70"}),
            CDS("AraC", "ATGGCTGAAGCGCAAATGATCCCGGCGATGAACATCAGCAGCGTTAACGCCAGCCCGGTTTTACCGTTGATGCCGATGTGG", {"degradation_rate": 0.04, "output_role": "activator", "preferred_context": "sigma70"}),
            CDS("LuxR", "ATGAAAAACATAAATGCCGACGACACATACAGAATAATTACCGCGACTGCCTTGGCGGATGCGGTGAAATCCGCCGAC", {"degradation_rate": 0.045, "output_role": "activator", "preferred_context": "lux"}),
            CDS("CI", "ATGAGCACAAAAAAGAAACCATTACCCGCGCCGTTGTTGCGGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCAC", {"degradation_rate": 0.035, "output_role": "repressor", "preferred_context": "lambda"}),
            Terminator("B0010", "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTC", {"termination_efficiency": 0.88}),
            Terminator("B0012", "TCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCT", {"termination_efficiency": 0.91}),
            Terminator("B1001", "GCGTTTTGCCCTGATAGTGACCTGTGCTCAGGAAAGGCCGATAAAGG", {"termination_efficiency": 0.95}),
            Terminator("B1002", "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG", {"termination_efficiency": 0.97}),
            Insulator("RiboJ", "AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAAC", {"insulation_score": 0.95, "context_family": "sigma70"}),
            Insulator("BydvJ", "GCTCGGATCCACTAGTCCAGTGTACAAGAAAGCTGGGTCTAGATCC", {"insulation_score": 0.89, "context_family": "lux"}),
            Insulator("SarJ", "GGCTCGAGCTGATCACTAGTGGTACCAGCGGCCGCATGCCTGCAGG", {"insulation_score": 0.86, "context_family": "lambda"}),
        ]


class CircuitAssembler:
    def __init__(self, circuit_spec: CircuitSpec, parts_catalog: PartsCatalog) -> None:
        self.circuit_spec = circuit_spec
        self.parts_catalog = parts_catalog

    def assemble(self) -> DesignCandidate:
        designs = self.enumerate_designs(max_designs=1)
        if not designs:
            raise ValueError("No feasible designs could be assembled from the provided circuit spec")
        return designs[0]

    def enumerate_designs(self, max_designs: int = 100) -> List[DesignCandidate]:
        gates = self._extract_gates(self.circuit_spec)
        if not gates:
            return []

        per_gate_choices = [self._enumerate_gate_assignments(gate) for gate in gates]
        if any(len(options) == 0 for options in per_gate_choices):
            return []

        limits = [min(len(options), max(1, int(np.ceil(max_designs ** (1.0 / max(len(gates), 1)))))) for options in per_gate_choices]
        trimmed_options = [options[:limit] for options, limit in zip(per_gate_choices, limits)]

        design_candidates: List[DesignCandidate] = []
        for assignment_combo in product(*trimmed_options):
            parts: List[BioPart] = []
            gate_scores: List[float] = []
            gate_predictions: List[Dict[str, Any]] = []
            for assignment in assignment_combo:
                parts.extend(assignment["parts"])
                gate_scores.append(float(assignment["score"]))
                gate_predictions.append(dict(assignment["predicted_performance"]))

            score = float(np.mean(gate_scores))
            predicted_performance = {
                "gate_predictions": gate_predictions,
                "mean_gate_score": score,
                "total_length": int(sum(len(part.sequence) for part in parts)),
            }
            design_candidates.append(DesignCandidate(parts=parts, score=score, predicted_performance=predicted_performance))

        design_candidates.sort(key=lambda candidate: candidate.score, reverse=True)
        return design_candidates[:max_designs]

    def _enumerate_gate_assignments(self, gate: Any) -> List[Dict[str, Any]]:
        gate_dict = self.parts_catalog._normalize_gate(gate)
        compatible = self.parts_catalog.get_compatible_parts(gate)

        promoters = compatible["promoter"] or self.parts_catalog.get_parts("promoter")
        rbss = compatible["rbs"] or self.parts_catalog.get_parts("rbs")
        cdss = compatible["cds"] or self.parts_catalog.get_parts("cds")
        terminators = compatible["terminator"] or self.parts_catalog.get_parts("terminator")
        insulators = compatible["insulator"] or self.parts_catalog.get_parts("insulator")

        assignments: List[Dict[str, Any]] = []
        for promoter, rbs, cds, terminator, insulator in product(promoters, rbss, cdss, terminators, insulators):
            compatibility_bonus = self._compatibility_score(promoter, rbs, cds, insulator, gate_dict)
            if compatibility_bonus <= 0.0:
                continue
            score = self._score_assignment(promoter, rbs, terminator, insulator, gate_dict, compatibility_bonus)
            predicted_performance = {
                "target_rpu": gate_dict.get("target_rpu", promoter.parameters.get("rpu")),
                "promoter_rpu": promoter.parameters.get("rpu"),
                "rbs_translation_rate": rbs.parameters.get("translation_rate"),
                "terminator_efficiency": terminator.parameters.get("termination_efficiency"),
                "insulation_score": insulator.parameters.get("insulation_score"),
                "compatibility_bonus": compatibility_bonus,
                "selected_output": cds.name,
            }
            assignments.append(
                {
                    "parts": [insulator, promoter, rbs, cds, terminator],
                    "score": score,
                    "predicted_performance": predicted_performance,
                }
            )

        assignments.sort(key=lambda item: float(item["score"]), reverse=True)
        return assignments

    @staticmethod
    def _extract_gates(circuit_spec: Any) -> List[Any]:
        if circuit_spec is None:
            return []
        if isinstance(circuit_spec, Mapping):
            gates = circuit_spec.get("gates")
            if gates is not None:
                return list(gates)
        for attr in ("gates", "logic_gates", "nodes", "stages"):
            if hasattr(circuit_spec, attr):
                value = getattr(circuit_spec, attr)
                if value is not None:
                    return list(value)
        if isinstance(circuit_spec, Sequence) and not isinstance(circuit_spec, (str, bytes, bytearray)):
            return list(circuit_spec)
        return []

    @staticmethod
    def _compatibility_score(
        promoter: Promoter,
        rbs: RBS,
        cds: CDS,
        insulator: Insulator,
        gate_dict: Mapping[str, Any],
    ) -> float:
        score = 1.0
        promoter_family = promoter.parameters.get("compatibility_family")
        rbs_family = rbs.parameters.get("compatibility_family")
        cds_context = cds.parameters.get("preferred_context")
        insulator_context = insulator.parameters.get("context_family")
        desired_family = gate_dict.get("compatibility_family") or gate_dict.get("context")

        if promoter_family is not None and rbs_family is not None and promoter_family != rbs_family:
            return 0.0
        if cds_context is not None and promoter_family is not None and cds_context != promoter_family:
            score *= 0.85
        if insulator_context is not None and promoter_family is not None and insulator_context != promoter_family:
            score *= 0.8
        if desired_family is not None and promoter_family is not None and desired_family != promoter_family:
            score *= 0.75
        regulator = gate_dict.get("regulator")
        if regulator is not None and promoter.parameters.get("regulated_by") not in (None, regulator):
            score *= 0.6
        return score

    @staticmethod
    def _score_assignment(
        promoter: Promoter,
        rbs: RBS,
        terminator: Terminator,
        insulator: Insulator,
        gate_dict: Mapping[str, Any],
        compatibility_bonus: float,
    ) -> float:
        promoter_target = float(gate_dict.get("target_rpu", promoter.parameters.get("rpu", 1.0)))
        rbs_target = float(gate_dict.get("target_rbs_strength", rbs.parameters.get("translation_rate", 1.0)))
        promoter_rpu = float(promoter.parameters.get("rpu", 0.0))
        rbs_strength = float(rbs.parameters.get("translation_rate", 0.0))
        insulation_score = float(insulator.parameters.get("insulation_score", 0.0))
        terminator_efficiency = float(terminator.parameters.get("termination_efficiency", 0.0))

        promoter_match = 1.0 / (1.0 + abs(promoter_rpu - promoter_target) / max(promoter_target, 1e-6))
        rbs_match = 1.0 / (1.0 + abs(rbs_strength - rbs_target) / max(rbs_target, 1e-6))

        weight_vector = np.array([0.4, 0.3, 0.15, 0.15], dtype=float)
        feature_vector = np.array([promoter_match, rbs_match, insulation_score, terminator_efficiency], dtype=float)
        base_score = float(np.dot(weight_vector, feature_vector))
        return base_score * float(compatibility_bonus)


__all__ = [
    "BioPart",
    "Promoter",
    "RBS",
    "CDS",
    "Terminator",
    "Insulator",
    "PartsCatalog",
    "DesignCandidate",
    "CircuitAssembler",
]
