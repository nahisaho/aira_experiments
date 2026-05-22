"""SARS-CoV-2 5'UTR structure prediction case study.

This module integrates RNA folding components from the :mod:`rna_structure`
package into a single reproducible case study focused on the SARS-CoV-2 5'UTR.
It is designed to work both with the package's native predictors (when present)
and with lightweight internal fallbacks for standalone execution.
"""

from __future__ import annotations

import json
import math
import random
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    from .deep_covariation import CovariationIntegrator as _CovariationIntegrator
except ImportError:  # pragma: no cover - fallback for direct script execution.
    try:
        from rna_structure.deep_covariation import CovariationIntegrator as _CovariationIntegrator
    except ImportError:  # pragma: no cover - optional dependency.
        _CovariationIntegrator = None

try:
    from .turner_model import ZukerMFE as _ZukerMFE
except ImportError:  # pragma: no cover - fallback for direct script execution.
    try:
        from rna_structure.turner_model import ZukerMFE as _ZukerMFE
    except ImportError:  # pragma: no cover - optional dependency.
        _ZukerMFE = None

try:
    from .chemical_probing import ProbeConstrainedFolding as _ProbeConstrainedFolding
except ImportError:  # pragma: no cover - optional dependency.
    try:
        from rna_structure.chemical_probing import ProbeConstrainedFolding as _ProbeConstrainedFolding
    except ImportError:  # pragma: no cover - module may not exist in lightweight installs.
        _ProbeConstrainedFolding = None

try:
    from .pseudoknot import HeuristicPseudoknot as _HeuristicPseudoknot
except ImportError:  # pragma: no cover - optional dependency.
    try:
        from rna_structure.pseudoknot import HeuristicPseudoknot as _HeuristicPseudoknot
    except ImportError:  # pragma: no cover - module may not exist in lightweight installs.
        _HeuristicPseudoknot = None

Pair = Tuple[int, int]
BRACKET_TYPES: List[Tuple[str, str]] = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
CANONICAL_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}

ZukerMFE = _ZukerMFE
ProbeConstrainedFolding = _ProbeConstrainedFolding
HeuristicPseudoknot = _HeuristicPseudoknot
CovariationIntegrator = _CovariationIntegrator


def _is_canonical(base_a: str, base_b: str) -> bool:
    return (base_a, base_b) in CANONICAL_PAIRS


def _structure_to_pairs(structure: str) -> List[Pair]:
    """Convert dot-bracket notation into a sorted list of 0-based base pairs."""

    opening_to_closing = {opening: closing for opening, closing in BRACKET_TYPES}
    closing_to_opening = {closing: opening for opening, closing in BRACKET_TYPES}
    stacks: Dict[str, List[int]] = {opening: [] for opening, _ in BRACKET_TYPES}
    pairs: List[Pair] = []

    for index, char in enumerate(structure):
        if char in opening_to_closing:
            stacks[char].append(index)
        elif char in closing_to_opening:
            opening = closing_to_opening[char]
            if stacks[opening]:
                pairs.append((stacks[opening].pop(), index))
    return sorted(pairs)


def _pairs_cross(first: Pair, second: Pair) -> bool:
    """Return True when two pairs form a pseudoknotted crossing."""

    i, j = first
    k, l = second
    return (i < k < j < l) or (k < i < l < j)


def _pairs_to_structure(length: int, pairs: Iterable[Pair]) -> str:
    """Convert a pair list into layered dot-bracket notation."""

    chars = ["."] * length
    layers: List[List[Pair]] = [[] for _ in BRACKET_TYPES]

    for i, j in sorted(set(pairs)):
        if not (0 <= i < j < length):
            continue
        for layer_index, brackets in enumerate(BRACKET_TYPES):
            if all(not _pairs_cross((i, j), existing) for existing in layers[layer_index]):
                layers[layer_index].append((i, j))
                chars[i] = brackets[0]
                chars[j] = brackets[1]
                break
    return "".join(chars)


def _score_structure_energy(sequence: str, pairs: Iterable[Pair]) -> float:
    """Estimate a simple pseudo-free-energy for ranking structures."""

    total = 0.0
    for i, j in pairs:
        total -= 1.6 if {sequence[i], sequence[j]} == {"G", "C"} else 1.0
        if (sequence[i], sequence[j]) in {("G", "U"), ("U", "G")}:
            total += 0.25
    return round(total, 3)


def _nussinov_fold(
    sequence: str,
    min_loop_length: int = 3,
    pair_bonus: Optional[Callable[[int, int], float]] = None,
) -> List[Pair]:
    """Fold an RNA with a lightweight Nussinov-style dynamic program."""

    n = len(sequence)
    dp = np.zeros((n, n), dtype=float)

    for span in range(min_loop_length + 1, n):
        for i in range(0, n - span):
            j = i + span
            best = max(dp[i + 1, j], dp[i, j - 1])
            if _is_canonical(sequence[i], sequence[j]):
                bonus = 1.0 if pair_bonus is None else max(0.0, pair_bonus(i, j))
                best = max(best, dp[i + 1, j - 1] + bonus)
            for k in range(i + 1, j):
                best = max(best, dp[i, k] + dp[k + 1, j])
            dp[i, j] = best

    pairs: List[Pair] = []

    def traceback(i: int, j: int) -> None:
        if i >= j:
            return
        if math.isclose(dp[i, j], dp[i + 1, j], abs_tol=1e-8):
            traceback(i + 1, j)
            return
        if math.isclose(dp[i, j], dp[i, j - 1], abs_tol=1e-8):
            traceback(i, j - 1)
            return
        if _is_canonical(sequence[i], sequence[j]):
            bonus = 1.0 if pair_bonus is None else max(0.0, pair_bonus(i, j))
            if math.isclose(dp[i, j], dp[i + 1, j - 1] + bonus, abs_tol=1e-8):
                pairs.append((i, j))
                traceback(i + 1, j - 1)
                return
        for k in range(i + 1, j):
            if math.isclose(dp[i, j], dp[i, k] + dp[k + 1, j], abs_tol=1e-8):
                traceback(i, k)
                traceback(k + 1, j)
                return

    if n:
        traceback(0, n - 1)
    return sorted(pairs)


def _extract_result_payload(result: Any, sequence: str) -> Dict[str, Any]:
    """Normalize predictor outputs from diverse APIs into a common payload."""

    if isinstance(result, tuple):
        structure = str(result[0])
        energy = float(result[1]) if len(result) > 1 else _score_structure_energy(sequence, _structure_to_pairs(structure))
        extra = result[2] if len(result) > 2 and isinstance(result[2], dict) else {}
        return {"structure": structure, "energy": energy, **extra}
    if isinstance(result, dict):
        structure = str(result.get("structure") or result.get("dot_bracket") or "")
        energy = float(result.get("energy", _score_structure_energy(sequence, _structure_to_pairs(structure))))
        payload = dict(result)
        payload["structure"] = structure
        payload["energy"] = energy
        return payload
    if isinstance(result, str):
        structure = result
        return {"structure": structure, "energy": _score_structure_energy(sequence, _structure_to_pairs(structure))}
    raise TypeError(f"Unsupported predictor output type: {type(result)!r}")


def _instantiate_external_model(model_class: Optional[type], *args: Any, **kwargs: Any) -> Optional[Any]:
    """Instantiate a predictor class while tolerating signature differences."""

    if model_class is None:
        return None
    constructor_candidates = [
        lambda: model_class(*args, **kwargs),
        lambda: model_class(sequence=kwargs.get("sequence")),
        lambda: model_class(kwargs.get("sequence")),
        model_class,
    ]
    for constructor in constructor_candidates:
        try:
            return constructor()
        except Exception:
            continue
    return None


def _call_external_predictor(model: Any, *args: Any, **kwargs: Any) -> Optional[Dict[str, Any]]:
    """Call an optional external predictor and normalize its output."""

    if model is None:
        return None

    methods = []
    for name in ("predict", "fold", "run", "fit_predict", "__call__"):
        method = getattr(model, name, None)
        if callable(method):
            methods.append(method)

    call_patterns = [
        lambda method: method(*args, **kwargs),
        lambda method: method(kwargs.get("sequence"), kwargs.get("shape_data"), kwargs.get("msa")),
        lambda method: method(kwargs.get("sequence"), kwargs.get("shape_data")),
        lambda method: method(kwargs.get("sequence"), kwargs.get("msa")),
        lambda method: method(kwargs.get("sequence")),
        lambda method: method(),
    ]

    for method in methods:
        for pattern in call_patterns:
            try:
                result = pattern(method)
                if result is not None:
                    return _extract_result_payload(result, kwargs.get("sequence", args[0] if args else ""))
            except Exception:
                continue
    return None


def _build_pairs_from_ranges(*range_specs: Tuple[range, range]) -> List[Pair]:
    """Create base pairs by zipping left and reversed right ranges."""

    pairs: List[Pair] = []
    for left_range, right_range in range_specs:
        right_positions = list(right_range)
        for left, right in zip(left_range, reversed(right_positions)):
            pairs.append((left - 1, right - 1))
    return pairs


@dataclass
class SARSCoV2Data:
    """Reference data and synthetic probing signals for the SARS-CoV-2 5'UTR."""

    random_seed: int = 7
    full_sequence: str = field(init=False)
    sequence: str = field(init=False)
    stem_loop_regions: Dict[str, Dict[str, Any]] = field(init=False)
    known_structure: str = field(init=False)
    shape_data: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        raw_sequence = (
            "AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUCGGCUGCAUGCUUAGUGCACUCACGCAGUAUAAUUAAUAACUAAUUACUGUCGUUGACAGGACACGAGUAACUCGUCUAUCUUCUGCAGGCUGCUUACGGUUUCGUCCGUGUUGCAGCCGAUCAUCAGCACAUCUAGGUUUCGUCCGGGUGUGACCGAAAGGUAAGAUGGAGAGCCUUGUCCCUGGUUUCAACGAGAAAACACACGUCCAACUCAGUUUGCCUGUUUUACAGGUUCGCGACGUGCUCGUACGUGGCUUUGGAGACUCCGUGGAGGAGGUCUUAUCAGAGGCACGUCAACAUCUUAAAGAUGGCACUUGUGGCUUAGUAGAAGUUGAAAAAGGCGUUUUGCCUCAACUUGAACAGCCCUAUGUGUUCAUCAAACGUUCGGAUGCUCGAACUGCACCUCAUGGUCAUGUUAUGGUUGAGCUGGUAGCAGAACUCGAAGGCAUUCAGUACGGUCGUAGUGGUGAGACACUUGGUGUCCUUGUCCCUCAUGUGGGCGAAAUACCAGUGGCUUACCGCAAGGUUCUUCUUCGUAAGAACGGUAAUAAAGGAGCUGGUGGCCAUAGUUACGGCGCCGAUCUAAAGUCAUUUGACUUAGGCGACGAGCUUGGCACUGAUCCUUAUGAAGAUUUUCAAGAAAACUGGAACACUAAACAUAGCAGUGGUGUUACCCGUGAACUCAUGCGUGAGCUUAACGGAGGG"
        )
        self.full_sequence = raw_sequence
        self.sequence = raw_sequence[:265]
        self.stem_loop_regions = {
            "SL1": {
                "start": 7,
                "end": 33,
                "description": "Stem-loop 1 with a CUCC-rich apical loop.",
            },
            "SL2": {
                "start": 45,
                "end": 59,
                "description": "Compact stem-loop 2.",
            },
            "SL3": {
                "start": 62,
                "end": 75,
                "description": "TRS-L-containing stem-loop 3.",
            },
            "SL4": {
                "start": 82,
                "end": 120,
                "description": "Extended stem-loop 4 with an internal loop.",
            },
            "SL5": {
                "start": 150,
                "end": 265,
                "description": "Large branched stem-loop 5 near the ORF1a start codon.",
            },
        }
        self.known_structure = self._build_known_structure()
        self.shape_data = self._simulate_shape_data()

    def _build_known_structure(self) -> str:
        """Construct a plausible dot-bracket model spanning the annotated stem-loops."""

        pairs: List[Pair] = []
        pairs.extend(_build_pairs_from_ranges((range(7, 18), range(23, 34))))
        pairs.extend(_build_pairs_from_ranges((range(45, 50), range(55, 60))))
        pairs.extend(_build_pairs_from_ranges((range(62, 67), range(71, 76))))
        pairs.extend(_build_pairs_from_ranges((range(82, 88), range(115, 121))))
        pairs.extend(_build_pairs_from_ranges((range(92, 97), range(104, 109))))
        pairs.extend(_build_pairs_from_ranges((range(150, 161), range(255, 266))))
        pairs.extend(_build_pairs_from_ranges((range(170, 177), range(194, 201))))
        pairs.extend(_build_pairs_from_ranges((range(206, 213), range(230, 237))))
        pairs.extend(_build_pairs_from_ranges((range(216, 220), range(224, 228))))
        return _pairs_to_structure(len(self.sequence), pairs)

    def _simulate_shape_data(self) -> np.ndarray:
        """Generate reproducible synthetic SHAPE reactivities from the consensus structure."""

        rng = np.random.default_rng(self.random_seed)
        paired_positions = {index for pair in _structure_to_pairs(self.known_structure) for index in pair}
        shape = np.empty(len(self.sequence), dtype=float)
        for idx in range(len(self.sequence)):
            base_mean = 0.18 if idx in paired_positions else 0.82
            jitter = rng.normal(loc=0.0, scale=0.08 if idx in paired_positions else 0.12)
            shape[idx] = np.clip(base_mean + jitter, 0.0, 1.4)
        return shape

    def get_sequence(self) -> str:
        """Return the 5'UTR-focused sequence used for structure prediction."""

        return self.sequence

    def get_known_structure(self) -> str:
        """Return the consensus-like reference structure in dot-bracket notation."""

        return self.known_structure

    def get_shape_data(self) -> np.ndarray:
        """Return simulated SHAPE reactivity data for the 5'UTR."""

        return self.shape_data.copy()

    def get_stem_loop_regions(self) -> Dict[str, Dict[str, Any]]:
        """Return annotated stem-loop regions (1-based inclusive coordinates)."""

        return {name: dict(values) for name, values in self.stem_loop_regions.items()}


class SARSCoV2Predictor:
    """Run and compare a panel of RNA secondary-structure prediction strategies."""

    def __init__(self, data: SARSCoV2Data) -> None:
        self.data = data
        self.sequence = data.get_sequence()
        self.shape_data = data.get_shape_data()
        self.known_structure = data.get_known_structure()
        self.known_pairs = _structure_to_pairs(self.known_structure)
        self.synthetic_msa = self._generate_synthetic_msa()

    def _generate_synthetic_msa(self, n_sequences: int = 24) -> List[str]:
        """Generate a small synthetic alignment that preserves reference pairing signals."""

        rng = random.Random(self.data.random_seed)
        paired_positions = {idx for pair in self.known_pairs for idx in pair}
        canonical_templates = [("G", "C"), ("C", "G"), ("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")]
        msa = [self.sequence]

        for _ in range(n_sequences - 1):
            chars = list(self.sequence)
            for i, j in self.known_pairs:
                if rng.random() < 0.22:
                    new_left, new_right = rng.choice(canonical_templates)
                    chars[i] = new_left
                    chars[j] = new_right
            for idx in range(len(chars)):
                if idx in paired_positions:
                    continue
                if rng.random() < 0.04:
                    chars[idx] = rng.choice(["A", "U", "G", "C"])
            msa.append("".join(chars))
        return msa

    def _predict_basic_mfe(self) -> Dict[str, Any]:
        if ZukerMFE is not None:
            try:
                model = ZukerMFE(self.sequence)
                structure = model.predict()
                energy = float(model.mfe()) if hasattr(model, "mfe") else _score_structure_energy(self.sequence, _structure_to_pairs(structure))
                return {"structure": structure, "energy": energy, "source": "external"}
            except Exception:
                pass

        model = _instantiate_external_model(ZukerMFE, sequence=self.sequence)
        external = _call_external_predictor(model, self.sequence, sequence=self.sequence)
        if external is not None:
            return {"structure": external["structure"], "energy": external["energy"], "source": "external"}

        pairs = _nussinov_fold(self.sequence)
        structure = _pairs_to_structure(len(self.sequence), pairs)
        return {"structure": structure, "energy": _score_structure_energy(self.sequence, pairs), "source": "fallback"}

    def _predict_shape_constrained(self) -> Dict[str, Any]:
        if ProbeConstrainedFolding is not None:
            try:
                model = ProbeConstrainedFolding()
                structure, energy = model.fold_with_shape(self.sequence, self.shape_data)
                return {"structure": structure, "energy": float(energy), "source": "external"}
            except Exception:
                pass

        model = _instantiate_external_model(ProbeConstrainedFolding, sequence=self.sequence, shape_data=self.shape_data)
        external = _call_external_predictor(
            model,
            self.sequence,
            self.shape_data,
            sequence=self.sequence,
            shape_data=self.shape_data,
        )
        if external is not None:
            return {"structure": external["structure"], "energy": external["energy"], "source": "external"}

        def bonus(i: int, j: int) -> float:
            return max(0.0, 1.25 - 0.75 * float(self.shape_data[i] + self.shape_data[j]))

        pairs = _nussinov_fold(self.sequence, pair_bonus=bonus)
        structure = _pairs_to_structure(len(self.sequence), pairs)
        return {"structure": structure, "energy": _score_structure_energy(self.sequence, pairs), "source": "fallback"}

    def _predict_pseudoknot(self) -> Dict[str, Any]:
        if HeuristicPseudoknot is not None:
            try:
                model = HeuristicPseudoknot(self.sequence)
                energy, pairs = model.predict_structure()
                structure = _pairs_to_structure(len(self.sequence), pairs)
                return {"structure": structure, "energy": float(energy), "source": "external"}
            except Exception:
                pass

        model = _instantiate_external_model(HeuristicPseudoknot, sequence=self.sequence)
        external = _call_external_predictor(model, self.sequence, sequence=self.sequence)
        if external is not None:
            return {"structure": external["structure"], "energy": external["energy"], "source": "external"}

        base_pairs = _structure_to_pairs(self._predict_shape_constrained()["structure"])
        used = {idx for pair in base_pairs for idx in pair}
        extra_pairs: List[Pair] = []
        for i in range(0, len(self.sequence) - 20):
            if i in used or self.shape_data[i] > 0.45:
                continue
            for j in range(len(self.sequence) - 1, i + 8, -1):
                if j in used or self.shape_data[j] > 0.45:
                    continue
                candidate = (i, j)
                if not _is_canonical(self.sequence[i], self.sequence[j]):
                    continue
                if any(_pairs_cross(candidate, pair) for pair in base_pairs) and all(idx not in used for idx in candidate):
                    extra_pairs.append(candidate)
                    used.update(candidate)
                    break
            if len(extra_pairs) >= 3:
                break
        combined_pairs = sorted(base_pairs + extra_pairs)
        structure = _pairs_to_structure(len(self.sequence), combined_pairs)
        return {"structure": structure, "energy": _score_structure_energy(self.sequence, combined_pairs), "source": "fallback"}

    def _covariation_matrix(self, msa: Sequence[str]) -> np.ndarray:
        """Build a simple covariation support matrix from a synthetic alignment."""

        n = len(msa[0])
        matrix = np.zeros((n, n), dtype=float)
        if not msa:
            return matrix

        for i in range(n):
            for j in range(i + 4, n):
                pair_types = Counter((sequence[i], sequence[j]) for sequence in msa)
                canonical_support = sum(count for pair, count in pair_types.items() if pair in CANONICAL_PAIRS) / len(msa)
                diversity = len([pair for pair in pair_types if pair in CANONICAL_PAIRS]) / 6.0
                matrix[i, j] = canonical_support * (0.7 + 0.3 * diversity)
                matrix[j, i] = matrix[i, j]
        return matrix

    def _predict_covariation_enhanced(self) -> Dict[str, Any]:
        cov = self._covariation_matrix(self.synthetic_msa)

        if CovariationIntegrator is not None and ZukerMFE is not None:
            try:
                integrator = CovariationIntegrator(lambda_weight=1.0)
                structure, energy = integrator.constrained_fold(self.sequence, -cov, ZukerMFE(self.sequence))
                return {
                    "structure": structure,
                    "energy": float(energy),
                    "source": "external",
                    "covariation_mean": float(np.mean(cov[np.triu_indices_from(cov, 1)])),
                }
            except Exception:
                pass

        model = _instantiate_external_model(CovariationIntegrator, sequence=self.sequence, msa=self.synthetic_msa)
        external = _call_external_predictor(
            model,
            self.sequence,
            self.synthetic_msa,
            sequence=self.sequence,
            msa=self.synthetic_msa,
        )
        if external is not None:
            external["covariation_mean"] = float(np.mean(cov[np.triu_indices_from(cov, 1)]))
            return {"structure": external["structure"], "energy": external["energy"], "source": "external", "covariation_mean": external["covariation_mean"]}

        def bonus(i: int, j: int) -> float:
            return 0.25 + 1.4 * cov[i, j]

        pairs = _nussinov_fold(self.sequence, pair_bonus=bonus)
        structure = _pairs_to_structure(len(self.sequence), pairs)
        return {
            "structure": structure,
            "energy": _score_structure_energy(self.sequence, pairs),
            "source": "fallback",
            "covariation_mean": float(np.mean(cov[np.triu_indices_from(cov, 1)])),
        }

    def _predict_combined(self, predictions: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        pair_counter: Counter[Pair] = Counter()
        pair_sources: Dict[Pair, List[str]] = {}
        for method, payload in predictions.items():
            for pair in _structure_to_pairs(payload["structure"]):
                pair_counter[pair] += 1
                pair_sources.setdefault(pair, []).append(method)

        selected: List[Pair] = []
        occupied: set[int] = set()
        for pair, votes in sorted(pair_counter.items(), key=lambda item: (-item[1], item[0][0], item[0][1])):
            if pair[0] in occupied or pair[1] in occupied:
                continue
            if votes < 2 and "covariation_enhanced" not in pair_sources[pair]:
                continue
            selected.append(pair)
            occupied.update(pair)

        structure = _pairs_to_structure(len(self.sequence), selected)
        mean_energy = float(np.mean([payload["energy"] for payload in predictions.values()]))
        return {
            "structure": structure,
            "energy": round(mean_energy, 3),
            "source": "ensemble",
            "consensus_pairs": len(selected),
        }

    def run_all_predictions(self) -> Dict[str, Dict[str, Any]]:
        """Run all available prediction strategies and return normalized results."""

        predictions = {
            "basic_mfe": self._predict_basic_mfe(),
            "shape_constrained": self._predict_shape_constrained(),
            "pseudoknot_aware": self._predict_pseudoknot(),
            "covariation_enhanced": self._predict_covariation_enhanced(),
        }
        predictions["combined"] = self._predict_combined(predictions)
        return predictions

    def compare_predictions(self, predictions: Dict[str, Dict[str, Any]], known_structure: str) -> Dict[str, Dict[str, float]]:
        """Compare predicted structures against the reference-like structure."""

        known_pairs = set(_structure_to_pairs(known_structure))
        paired_known = {idx for pair in known_pairs for idx in pair}
        metrics: Dict[str, Dict[str, float]] = {}

        for method, payload in predictions.items():
            predicted_pairs = set(_structure_to_pairs(payload["structure"]))
            paired_pred = {idx for pair in predicted_pairs for idx in pair}
            tp = len(predicted_pairs & known_pairs)
            fp = len(predicted_pairs - known_pairs)
            fn = len(known_pairs - predicted_pairs)
            sensitivity = tp / len(known_pairs) if known_pairs else 0.0
            ppv = tp / len(predicted_pairs) if predicted_pairs else 0.0
            f1 = 2 * sensitivity * ppv / (sensitivity + ppv) if (sensitivity + ppv) else 0.0
            paired_state_accuracy = sum(
                ((index in paired_pred) == (index in paired_known)) for index in range(len(self.sequence))
            ) / len(self.sequence)
            metrics[method] = {
                "tp": float(tp),
                "fp": float(fp),
                "fn": float(fn),
                "sensitivity": round(sensitivity, 3),
                "ppv": round(ppv, 3),
                "f1": round(f1, 3),
                "paired_state_accuracy": round(paired_state_accuracy, 3),
                "predicted_pairs": float(len(predicted_pairs)),
            }
            payload["metrics"] = metrics[method]
        return metrics


class SARSCoV2Analyzer:
    """Interpret predicted structures in terms of known SARS-CoV-2 leader features."""

    def __init__(self, data: SARSCoV2Data) -> None:
        self.data = data
        self.known_sls = data.get_stem_loop_regions()

    def analyze_stem_loops(self, structure: str, sequence: str) -> List[Dict[str, Any]]:
        """Summarize structure content across the expected SL1-SL5 regions."""

        pairs = _structure_to_pairs(structure)
        analyses: List[Dict[str, Any]] = []

        for name, region in self.known_sls.items():
            start_idx = region["start"] - 1
            end_idx = region["end"] - 1
            local_pairs = [(i, j) for i, j in pairs if start_idx <= i < j <= end_idx]
            local_structure = structure[start_idx : end_idx + 1]
            segment = sequence[start_idx : end_idx + 1]
            stems = sum(
                1
                for idx, char in enumerate(local_structure)
                if char in "([{<" and (idx == 0 or local_structure[idx - 1] == ".")
            )
            analyses.append(
                {
                    "name": name,
                    "start": region["start"],
                    "end": region["end"],
                    "length": region["end"] - region["start"] + 1,
                    "base_pairs": len(local_pairs),
                    "paired_fraction": round((2 * len(local_pairs)) / max(1, len(local_structure)), 3),
                    "unpaired_bases": local_structure.count("."),
                    "branching": "branched" if stems >= 2 else "hairpin",
                    "contains_trs_leader": "ACGAAC" in segment,
                    "contains_start_codon": "AUG" in segment,
                    "description": region["description"],
                }
            )
        return analyses

    def compare_with_known(self, predicted_sls: List[Dict[str, Any]], known_sls: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Match predicted stem-loops to known SL1-SL5 regions by interval overlap."""

        matches: Dict[str, Dict[str, Any]] = {}
        overlap_scores: List[float] = []
        for known_name, known_region in known_sls.items():
            known_start, known_end = known_region["start"], known_region["end"]
            best_match: Optional[Dict[str, Any]] = None
            best_score = 0.0
            for predicted in predicted_sls:
                overlap = max(0, min(known_end, predicted["end"]) - max(known_start, predicted["start"]) + 1)
                union = max(known_end, predicted["end"]) - min(known_start, predicted["start"]) + 1
                score = overlap / union if union else 0.0
                if score > best_score:
                    best_score = score
                    best_match = predicted
            matches[known_name] = {
                "best_match": best_match["name"] if best_match else None,
                "jaccard_overlap": round(best_score, 3),
                "description": known_region["description"],
            }
            overlap_scores.append(best_score)

        matched = sum(1 for info in matches.values() if info["jaccard_overlap"] >= 0.2)
        return {
            "matched_regions": matched,
            "total_known_regions": len(known_sls),
            "mean_jaccard_overlap": round(float(np.mean(overlap_scores)) if overlap_scores else 0.0, 3),
            "matches": matches,
        }

    def detect_pseudoknots(self, structure: str) -> List[Dict[str, int]]:
        """Identify crossing base-pair interactions in a predicted structure."""

        pairs = _structure_to_pairs(structure)
        pseudoknots: List[Dict[str, int]] = []
        for index, first in enumerate(pairs):
            for second in pairs[index + 1 :]:
                if _pairs_cross(first, second):
                    pseudoknots.append(
                        {
                            "left_start": first[0] + 1,
                            "left_end": first[1] + 1,
                            "right_start": second[0] + 1,
                            "right_end": second[1] + 1,
                        }
                    )
        return pseudoknots

    def analyze_frameshifting_signal(self, sequence: str) -> Dict[str, Any]:
        """Report whether the canonical frameshifting region is present in the analyzed window."""

        slippery_site = "UUUAAAC"
        index = sequence.find(slippery_site)
        return {
            "slippery_site_detected": index != -1,
            "position": index + 1 if index != -1 else None,
            "note": "The canonical ORF1a/1b frameshifting signal lies downstream of the 5'UTR and is usually absent here.",
        }


class SARSCoV2Visualization:
    """Text-only summaries for comparing structure predictions."""

    def format_comparison_table(self, results: Dict[str, Dict[str, Any]]) -> str:
        """Return a plain-text comparison table for prediction accuracy."""

        headers = ["Method", "Energy", "Sens", "PPV", "F1", "PairAcc"]
        rows = [headers]
        for method, payload in results.items():
            metrics = payload.get("metrics", {})
            rows.append(
                [
                    method,
                    f"{payload.get('energy', 0.0):.2f}",
                    f"{metrics.get('sensitivity', 0.0):.3f}",
                    f"{metrics.get('ppv', 0.0):.3f}",
                    f"{metrics.get('f1', 0.0):.3f}",
                    f"{metrics.get('paired_state_accuracy', 0.0):.3f}",
                ]
            )

        widths = [max(len(row[col]) for row in rows) for col in range(len(headers))]
        formatted_rows = []
        for row_index, row in enumerate(rows):
            formatted = " | ".join(value.ljust(widths[col]) for col, value in enumerate(row))
            formatted_rows.append(formatted)
            if row_index == 0:
                formatted_rows.append("-+-".join("-" * width for width in widths))
        return "\n".join(formatted_rows)

    def generate_arc_plot_data(self, structure: str) -> List[Tuple[int, int, int]]:
        """Return arc tuples as (1-based start, 1-based end, span)."""

        return [(i + 1, j + 1, j - i) for i, j in _structure_to_pairs(structure)]

    def format_stem_loop_summary(self, analysis: Sequence[Dict[str, Any]]) -> str:
        """Return a readable summary of predicted stem-loop features."""

        lines = []
        for item in analysis:
            lines.append(
                f"- {item['name']}: nt {item['start']}-{item['end']}, {item['branching']}, "
                f"{item['base_pairs']} bp, paired_fraction={item['paired_fraction']:.3f}, "
                f"TRS-L={'yes' if item['contains_trs_leader'] else 'no'}, "
                f"AUG={'yes' if item['contains_start_codon'] else 'no'}"
            )
        return "\n".join(lines)


def _make_json_safe(value: Any) -> Any:
    """Recursively convert numpy values into JSON-serializable Python objects."""

    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {key: _make_json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_make_json_safe(item) for item in value]
    return value


def run_case_study() -> Dict[str, Any]:
    """Execute the full SARS-CoV-2 5'UTR structure prediction case study."""

    data = SARSCoV2Data()
    predictor = SARSCoV2Predictor(data)
    analyzer = SARSCoV2Analyzer(data)
    visualization = SARSCoV2Visualization()

    predictions = predictor.run_all_predictions()
    comparisons = predictor.compare_predictions(predictions, data.get_known_structure())

    analyses: Dict[str, Any] = {}
    for method, payload in predictions.items():
        stem_loops = analyzer.analyze_stem_loops(payload["structure"], data.get_sequence())
        analyses[method] = {
            "stem_loops": stem_loops,
            "known_region_overlap": analyzer.compare_with_known(stem_loops, data.get_stem_loop_regions()),
            "pseudoknots": analyzer.detect_pseudoknots(payload["structure"]),
            "frameshifting_signal": analyzer.analyze_frameshifting_signal(data.get_sequence()),
            "arc_plot_data": visualization.generate_arc_plot_data(payload["structure"]),
        }

    summary_table = visualization.format_comparison_table(predictions)
    combined_summary = visualization.format_stem_loop_summary(analyses["combined"]["stem_loops"])

    results = {
        "sequence_length": len(data.get_sequence()),
        "sequence": data.get_sequence(),
        "known_structure": data.get_known_structure(),
        "shape_data": data.get_shape_data(),
        "stem_loop_regions": data.get_stem_loop_regions(),
        "predictions": predictions,
        "comparisons": comparisons,
        "analyses": analyses,
        "summary_table": summary_table,
        "combined_stem_loop_summary": combined_summary,
        "synthetic_msa_size": len(predictor.synthetic_msa),
    }

    workspace_root = Path(__file__).resolve().parent.parent
    results_dir = workspace_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    json_path = results_dir / "sars_cov2_case_study.json"
    text_path = results_dir / "sars_cov2_case_study_summary.txt"
    shape_path = results_dir / "sars_cov2_shape_data.csv"

    json_path.write_text(json.dumps(_make_json_safe(results), indent=2), encoding="utf-8")
    text_path.write_text(
        "SARS-CoV-2 5'UTR structure prediction case study\n"
        f"Sequence length: {len(data.get_sequence())} nt\n\n"
        "Comparison table\n"
        f"{summary_table}\n\n"
        "Combined stem-loop summary\n"
        f"{combined_summary}\n",
        encoding="utf-8",
    )
    np.savetxt(shape_path, data.get_shape_data(), delimiter=",", header="shape_reactivity", comments="")

    print(summary_table)
    print()
    print(combined_summary)

    return results


__all__ = [
    "SARSCoV2Analyzer",
    "SARSCoV2Data",
    "SARSCoV2Predictor",
    "SARSCoV2Visualization",
    "run_case_study",
]


if __name__ == "__main__":
    results = run_case_study()
