"""Heuristic riboswitch and functional RNA structure-function prediction tools.

The implementations in this module provide lightweight, self-contained predictors
for riboswitch-like RNAs without requiring an external RNA folding engine.
Predictions are based on canonical/Wobble pairing, simple dynamic programming,
and structure-aware heuristics that are suitable for prototyping and testing.

All reported energies are in kcal/mol.
"""

from __future__ import annotations

from difflib import SequenceMatcher
from typing import Any, Callable, Iterable, Mapping, Optional

import numpy as np
from scipy.spatial.distance import hamming
from scipy.special import expit, softmax

R_GAS_CONSTANT = 0.0019872041
DEFAULT_TEMPERATURE = 310.15
BASE_PAIR_ENERGIES: dict[tuple[str, str], float] = {
    ("G", "C"): -3.0,
    ("C", "G"): -3.0,
    ("A", "U"): -2.0,
    ("U", "A"): -2.0,
    ("G", "U"): -1.0,
    ("U", "G"): -1.0,
}
COMPLEMENTS: dict[str, tuple[str, ...]] = {
    "A": ("U",),
    "U": ("A", "G"),
    "G": ("C", "U"),
    "C": ("G",),
}


def _normalize_sequence(sequence: str) -> str:
    return "".join(base for base in sequence.upper().replace("T", "U") if base in "AUGC")


def _base_pair_energy(left: str, right: str) -> float:
    return BASE_PAIR_ENERGIES.get((left, right), 0.0)


def _pair_blocks_to_pairs(blocks: Iterable[tuple[int, int, int]]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for left_start, right_end, stem_len in blocks:
        for offset in range(stem_len):
            pairs.append((left_start + offset, right_end - offset))
    return sorted(pairs)


def _pairs_to_dot_bracket(length: int, pairs: Iterable[tuple[int, int]]) -> str:
    structure = ["."] * length
    for left, right in pairs:
        structure[left] = "("
        structure[right] = ")"
    return "".join(structure)


def _dot_bracket_to_pair_map(structure: str) -> dict[int, int]:
    stack: list[int] = []
    pair_map: dict[int, int] = {}
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                continue
            partner = stack.pop()
            pair_map[index] = partner
            pair_map[partner] = index
    return pair_map


def _pair_map_to_pairs(pair_map: Mapping[int, int]) -> list[tuple[int, int]]:
    return sorted((left, right) for left, right in pair_map.items() if left < right)


def _preferred_complement(base: str) -> str:
    return COMPLEMENTS.get(base, ("N",))[0]


def _sequence_from_pairs(
    length: int,
    pairs: list[tuple[int, int]],
    unpaired_pattern: str = "AUGCGAAU",
    fixed_bases: Optional[dict[int, str]] = None,
) -> str:
    seq = [""] * length
    pair_choices = [("G", "C"), ("C", "G"), ("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")]
    pair_map: dict[int, int] = {}
    for idx, (left, right) in enumerate(sorted(pairs)):
        pair_map[left] = right
        pair_map[right] = left
        seq[left], seq[right] = pair_choices[idx % len(pair_choices)]

    fixed = {pos: base.upper().replace("T", "U") for pos, base in (fixed_bases or {}).items()}
    for pos, base in fixed.items():
        if 0 <= pos < length:
            seq[pos] = base
            partner = pair_map.get(pos)
            if partner is not None and partner not in fixed:
                seq[partner] = _preferred_complement(base)

    pattern = _normalize_sequence(unpaired_pattern) or "AUGC"
    for idx in range(length):
        if not seq[idx]:
            seq[idx] = pattern[idx % len(pattern)]
    return "".join(seq)


def _replace_subsequence(sequence: str, start: int, subsequence: str) -> str:
    end = min(len(sequence), start + len(subsequence))
    repl = subsequence[: end - start]
    return f"{sequence[:start]}{repl}{sequence[end:]}"


def _find_unpaired_segments(structure: str) -> list[dict[str, Any]]:
    pair_map = _dot_bracket_to_pair_map(structure)
    segments: list[dict[str, Any]] = []
    start: Optional[int] = None
    for idx, char in enumerate(structure + "X"):
        if idx < len(structure) and char == ".":
            if start is None:
                start = idx
            continue
        if start is None:
            continue
        end = idx - 1
        left_flank = start - 1
        right_flank = end + 1
        loop_type = "external"
        if 0 <= left_flank < len(structure) and 0 <= right_flank < len(structure):
            if left_flank in pair_map and right_flank in pair_map and pair_map[left_flank] == right_flank:
                loop_type = "hairpin"
            elif left_flank in pair_map and right_flank in pair_map:
                loop_type = "junction"
            elif left_flank in pair_map or right_flank in pair_map:
                loop_type = "bulge"
        segments.append(
            {
                "start": start,
                "end": end,
                "positions": list(range(start, end + 1)),
                "size": end - start + 1,
                "loop_type": loop_type,
            }
        )
        start = None
    return segments


def _paired_fraction(structure: str, start: int, end: int) -> float:
    if end <= start:
        return 0.0
    region = structure[start:end]
    return 1.0 - (region.count(".") / max(1, len(region)))


def _contiguous_stems(structure: str, start: int = 0, end: Optional[int] = None) -> list[dict[str, Any]]:
    pair_map = _dot_bracket_to_pair_map(structure)
    end = len(structure) if end is None else end
    pairs = [(left, right) for left, right in _pair_map_to_pairs(pair_map) if start <= left < end and start < right <= end]
    stems: list[dict[str, Any]] = []
    idx = 0
    while idx < len(pairs):
        left, right = pairs[idx]
        stem_pairs = [(left, right)]
        idx += 1
        while idx < len(pairs) and pairs[idx][0] == stem_pairs[-1][0] + 1 and pairs[idx][1] == stem_pairs[-1][1] - 1:
            stem_pairs.append(pairs[idx])
            idx += 1
        stems.append(
            {
                "left_start": stem_pairs[0][0],
                "left_end": stem_pairs[-1][0],
                "right_start": stem_pairs[-1][1],
                "right_end": stem_pairs[0][1],
                "stem_length": len(stem_pairs),
                "pairs": stem_pairs,
            }
        )
    return stems


def _structure_energy(sequence: str, structure: str) -> float:
    pair_map = _dot_bracket_to_pair_map(structure)
    energy = 0.0
    for left, right in _pair_map_to_pairs(pair_map):
        energy += _base_pair_energy(sequence[left], sequence[right])
        if left + 1 not in pair_map and right - 1 not in pair_map:
            energy += 0.3
    for segment in _find_unpaired_segments(structure):
        if segment["loop_type"] == "hairpin":
            energy += 0.15 * abs(segment["size"] - 5)
        elif segment["loop_type"] == "junction":
            energy += 0.1 * max(0, segment["size"] - 6)
        else:
            energy += 0.03 * segment["size"]
    return float(energy)


def _base_pair_distance(structure_a: str, structure_b: str) -> int:
    pairs_a = set(_pair_map_to_pairs(_dot_bracket_to_pair_map(structure_a)))
    pairs_b = set(_pair_map_to_pairs(_dot_bracket_to_pair_map(structure_b)))
    return len(pairs_a.symmetric_difference(pairs_b))


def _structure_accuracy(predicted: str, reference: str) -> float:
    pred = np.array(list(predicted), dtype="U1")
    ref = np.array(list(reference), dtype="U1")
    return 1.0 - float(hamming(pred, ref))


def _pair_metrics(predicted: str, reference: str) -> dict[str, float]:
    pred_pairs = set(_pair_map_to_pairs(_dot_bracket_to_pair_map(predicted)))
    ref_pairs = set(_pair_map_to_pairs(_dot_bracket_to_pair_map(reference)))
    overlap = len(pred_pairs & ref_pairs)
    sensitivity = overlap / max(1, len(ref_pairs))
    precision = overlap / max(1, len(pred_pairs))
    f1 = 0.0 if sensitivity + precision == 0 else 2.0 * sensitivity * precision / (sensitivity + precision)
    return {
        "sensitivity": float(sensitivity),
        "precision": float(precision),
        "f1": float(f1),
    }


class RiboswitchDatabase:
    """Store curated riboswitch family prototypes and annotations.

    The stored consensus entries are compact, synthetic family representatives that
    preserve characteristic aptamer sizes, ligand annotations, and regulatory
    layout for downstream benchmarking and heuristic scoring.
    """

    def __init__(self) -> None:
        self._families = self._build_database()

    def _build_database(self) -> dict[str, dict[str, Any]]:
        families: dict[str, dict[str, Any]] = {}
        specs = [
            {
                "name": "TPP riboswitch",
                "key": "tpp",
                "ligand": "Thiamine pyrophosphate",
                "length": 104,
                "aptamer_end": 80,
                "blocks": [(0, 79, 6), (11, 28, 5), (34, 51, 5), (57, 71, 4)],
                "binding": [12, 13, 14, 36, 37, 38, 58, 59, 60],
                "pattern": "UGGAAUAC",
                "tail": "GGGCUUUUUUAGGAGGUAUGCCAA",
            },
            {
                "name": "SAM-I riboswitch",
                "key": "sam-i",
                "ligand": "S-adenosylmethionine",
                "length": 138,
                "aptamer_end": 110,
                "blocks": [(0, 109, 8), (14, 39, 6), (46, 71, 6), (78, 99, 5)],
                "binding": [18, 19, 20, 49, 50, 51, 80, 81, 82, 90],
                "pattern": "AGGAUUGC",
                "tail": "GCGCCCUUUUUUGGAGGAGUAUGCCAA",
            },
            {
                "name": "Adenine riboswitch",
                "key": "adenine",
                "ligand": "Adenine",
                "length": 94,
                "aptamer_end": 70,
                "blocks": [(0, 69, 5), (10, 27, 5), (33, 50, 5), (55, 66, 4)],
                "binding": [11, 12, 13, 35, 36, 37, 56, 57],
                "pattern": "GAAAGUCU",
                "tail": "AGGAGGUAAUGGCCUUUUUUGA",
            },
            {
                "name": "FMN riboswitch",
                "key": "fmn",
                "ligand": "Flavin mononucleotide",
                "length": 150,
                "aptamer_end": 120,
                "blocks": [(0, 119, 9), (16, 45, 7), (52, 81, 7), (88, 111, 6)],
                "binding": [17, 18, 19, 53, 54, 55, 89, 90, 91, 101],
                "pattern": "UGACGGAA",
                "tail": "GGCGCCUUUUUAGGAGGUAUGCGCGAAAUU",
            },
            {
                "name": "Glycine riboswitch",
                "key": "glycine",
                "ligand": "Glycine",
                "length": 118,
                "aptamer_end": 90,
                "blocks": [(0, 89, 7), (13, 36, 6), (43, 66, 6), (71, 84, 5)],
                "binding": [14, 15, 16, 45, 46, 47, 72, 73, 74],
                "pattern": "GGCAAAUU",
                "tail": "GGGCUUUUUUAGGAGGUAUGCCGAAU",
            },
        ]

        for spec in specs:
            pairs = _pair_blocks_to_pairs(spec["blocks"])
            fixed_bases = {pos: "A" if idx % 2 else "G" for idx, pos in enumerate(spec["binding"])}
            sequence = _sequence_from_pairs(spec["length"], pairs, spec["pattern"], fixed_bases)
            sequence = _replace_subsequence(sequence, spec["aptamer_end"], spec["tail"])
            structure = _pairs_to_dot_bracket(spec["length"], pairs)
            families[spec["key"]] = {
                "name": spec["name"],
                "consensus_sequence": sequence,
                "consensus_structure": structure,
                "aptamer_structure": structure[: spec["aptamer_end"]],
                "aptamer_range": (0, spec["aptamer_end"] - 1),
                "aptamer_length": spec["aptamer_end"],
                "ligand": spec["ligand"],
                "ligand_binding_site_positions": spec["binding"],
                "expression_platform_region": (spec["aptamer_end"], spec["length"] - 1),
            }
        return families

    def get_family(self, name: str) -> dict[str, Any]:
        """Return a riboswitch family entry by case-insensitive family name."""
        query = name.strip().lower()
        aliases = {
            "tpp riboswitch": "tpp",
            "thiamine pyrophosphate": "tpp",
            "sam-i riboswitch": "sam-i",
            "sam i riboswitch": "sam-i",
            "s-adenosylmethionine": "sam-i",
            "adenine riboswitch": "adenine",
            "purine riboswitch": "adenine",
            "fmn riboswitch": "fmn",
            "flavin mononucleotide": "fmn",
            "glycine riboswitch": "glycine",
        }
        key = aliases.get(query, query)
        if key not in self._families:
            raise KeyError(f"Unknown riboswitch family: {name}")
        return dict(self._families[key])

    def list_families(self) -> list[str]:
        """List the canonical family names available in the database."""
        return [self._families[key]["name"] for key in sorted(self._families)]


class StructuralSwitchPredictor:
    """Predict bistable riboswitch conformations using heuristic folding."""

    def __init__(self, database: Optional[RiboswitchDatabase] = None, temperature: float = DEFAULT_TEMPERATURE) -> None:
        self.database = database or RiboswitchDatabase()
        self.temperature = temperature

    def _best_family_match(self, sequence: str) -> tuple[dict[str, Any], float]:
        best_family = self.database.get_family(self.database.list_families()[0])
        best_score = -1.0
        for family_name in self.database.list_families():
            family = self.database.get_family(family_name)
            candidate = family["consensus_sequence"]
            score = SequenceMatcher(None, sequence[: min(len(sequence), len(candidate))], candidate[: min(len(sequence), len(candidate))]).ratio()
            if score > best_score:
                best_family = family
                best_score = score
        return best_family, float(best_score)

    def _score_matrix(
        self,
        sequence: str,
        bias: Optional[Callable[[int, int, int], float]] = None,
        suppress_pairs: Optional[set[tuple[int, int]]] = None,
    ) -> np.ndarray:
        n = len(sequence)
        matrix = np.zeros((n, n), dtype=float)
        suppressed = suppress_pairs or set()
        for left in range(n):
            for right in range(left + 1, n):
                energy = _base_pair_energy(sequence[left], sequence[right])
                if energy == 0.0:
                    continue
                score = -energy
                if bias is not None:
                    score *= bias(left, right, n)
                if (left, right) in suppressed or (right, left) in suppressed:
                    score *= 0.05
                matrix[left, right] = score
        return matrix

    def _fold_with_nussinov(
        self,
        sequence: str,
        min_loop: int = 3,
        bias: Optional[Callable[[int, int, int], float]] = None,
        suppress_pairs: Optional[set[tuple[int, int]]] = None,
    ) -> tuple[str, float, list[tuple[int, int]]]:
        sequence = _normalize_sequence(sequence)
        n = len(sequence)
        scores = self._score_matrix(sequence, bias=bias, suppress_pairs=suppress_pairs)
        dp = np.zeros((n, n), dtype=float)
        for span in range(min_loop + 1, n):
            for left in range(0, n - span):
                right = left + span
                best = max(dp[left + 1, right], dp[left, right - 1])
                if scores[left, right] > 0:
                    best = max(best, dp[left + 1, right - 1] + scores[left, right])
                if right - left > 1:
                    splits = dp[left, left:right] + dp[left + 1 : right + 1, right]
                    if len(splits):
                        best = max(best, float(np.max(splits)))
                dp[left, right] = best

        pairs: list[tuple[int, int]] = []
        epsilon = 1e-8

        def traceback(left: int, right: int) -> None:
            if left >= right:
                return
            if abs(dp[left, right] - dp[left + 1, right]) < epsilon:
                traceback(left + 1, right)
                return
            if abs(dp[left, right] - dp[left, right - 1]) < epsilon:
                traceback(left, right - 1)
                return
            if right - left > min_loop and scores[left, right] > 0 and abs(dp[left, right] - (dp[left + 1, right - 1] + scores[left, right])) < epsilon:
                pairs.append((left, right))
                traceback(left + 1, right - 1)
                return
            for split in range(left + 1, right):
                if abs(dp[left, right] - (dp[left, split] + dp[split + 1, right])) < epsilon:
                    traceback(left, split)
                    traceback(split + 1, right)
                    return

        if n:
            traceback(0, n - 1)
        structure = _pairs_to_dot_bracket(n, sorted(pairs))
        return structure, _structure_energy(sequence, structure), sorted(pairs)

    def _conformation_type(self, sequence: str, structure: str) -> tuple[str, dict[str, float]]:
        family, similarity = self._best_family_match(sequence)
        aptamer_end = min(len(sequence), family["aptamer_length"])
        aptamer_paired = _paired_fraction(structure, 0, aptamer_end)
        expression_paired = _paired_fraction(structure, aptamer_end, len(sequence))
        binding_positions = [pos for pos in family["ligand_binding_site_positions"] if pos < len(structure)]
        binding_open = np.mean([1.0 if structure[pos] == "." else 0.0 for pos in binding_positions]) if binding_positions else 0.0
        terminator_like = ExpressionPlatformAnalyzer(self).predict_regulatory_type(sequence, structure) == "transcriptional"
        aptamer_score = 0.55 * aptamer_paired + 0.30 * binding_open + 0.15 * similarity
        expression_score = 0.65 * expression_paired + 0.20 * float(terminator_like) + 0.15 * (1.0 - binding_open)
        label = "aptamer" if aptamer_score >= expression_score else "expression_platform"
        return label, {
            "aptamer_score": float(aptamer_score),
            "expression_score": float(expression_score),
        }

    def predict_conformations(self, sequence: str, n_suboptimal: int = 10, energy_range: float = 5.0) -> list[dict[str, Any]]:
        """Predict MFE and suboptimal conformations within an energy window."""
        sequence = _normalize_sequence(sequence)
        if not sequence:
            return []

        base_structure, _, base_pairs = self._fold_with_nussinov(sequence, min_loop=3)
        base_pair_set = set(base_pairs)
        pair_groups = [
            set(base_pairs[::2]),
            set(base_pairs[1::2]),
            {pair for pair in base_pairs if pair[0] < len(sequence) * 0.6},
            {pair for pair in base_pairs if pair[0] >= len(sequence) * 0.6},
        ]
        bias_functions: list[tuple[str, Callable[[int, int, int], float], int, set[tuple[int, int]]]] = [
            ("mfe", lambda i, j, n: 1.0, 3, set()),
            ("aptamer-biased", lambda i, j, n: 1.25 if j < int(0.72 * n) else 0.85, 3, set()),
            ("expression-biased", lambda i, j, n: 1.25 if i > int(0.40 * n) else 0.85, 3, set()),
            ("long-range", lambda i, j, n: 1.2 if (j - i) > int(0.30 * n) else 0.95, 3, set()),
            ("short-range", lambda i, j, n: 1.15 if (j - i) < int(0.20 * n) else 0.9, 2, set()),
            ("suppressed-even", lambda i, j, n: 1.0, 3, pair_groups[0]),
            ("suppressed-odd", lambda i, j, n: 1.0, 3, pair_groups[1]),
            ("suppressed-aptamer", lambda i, j, n: 1.0, 4, pair_groups[2]),
            ("suppressed-expression", lambda i, j, n: 1.0, 4, pair_groups[3]),
            ("relaxed-loop", lambda i, j, n: 1.05, 5, base_pair_set.intersection(pair_groups[0])),
        ]

        candidates: dict[str, dict[str, Any]] = {}
        for label, bias, min_loop, suppressed in bias_functions:
            structure, energy, pairs = self._fold_with_nussinov(sequence, min_loop=min_loop, bias=bias, suppress_pairs=suppressed)
            if structure in candidates:
                continue
            conformation_type, scores = self._conformation_type(sequence, structure)
            candidates[structure] = {
                "label": label,
                "structure": structure,
                "energy": float(energy),
                "pair_count": len(pairs),
                "pairs": pairs,
                "conformation_type": conformation_type,
                **scores,
            }

        ranked = sorted(candidates.values(), key=lambda item: item["energy"])
        best_energy = ranked[0]["energy"]
        windowed = [item for item in ranked if item["energy"] <= best_energy + energy_range]
        selected = windowed[: max(2, n_suboptimal)] if windowed else ranked[: max(2, n_suboptimal)]
        probabilities = softmax(-np.array([item["energy"] for item in selected]) / (R_GAS_CONSTANT * self.temperature))
        mfe_structure = selected[0]["structure"]
        for item, probability in zip(selected, probabilities, strict=False):
            item["ensemble_probability"] = float(probability)
            item["base_pair_distance_to_mfe"] = _base_pair_distance(item["structure"], mfe_structure)
        return selected

    def identify_switch_region(self, conformation1: Mapping[str, Any] | str, conformation2: Mapping[str, Any] | str) -> list[int]:
        """Return positions whose pairing partner changes between conformations."""
        structure1 = conformation1["structure"] if isinstance(conformation1, Mapping) else conformation1
        structure2 = conformation2["structure"] if isinstance(conformation2, Mapping) else conformation2
        map1 = _dot_bracket_to_pair_map(structure1)
        map2 = _dot_bracket_to_pair_map(structure2)
        return [idx for idx in range(min(len(structure1), len(structure2))) if map1.get(idx, -1) != map2.get(idx, -1)]

    def compute_switching_energy(self, sequence: str, struct1: str, struct2: str) -> float:
        """Compute the absolute energy gap required to switch between two structures."""
        sequence = _normalize_sequence(sequence)
        return float(abs(_structure_energy(sequence, struct2) - _structure_energy(sequence, struct1)))


class LigandBindingPredictor:
    """Predict riboswitch ligand-binding pockets from sequence and structure features."""

    def __init__(self, database: Optional[RiboswitchDatabase] = None) -> None:
        self.database = database or RiboswitchDatabase()

    def _best_family_match(self, sequence: str) -> tuple[dict[str, Any], float]:
        predictor = StructuralSwitchPredictor(self.database)
        return predictor._best_family_match(_normalize_sequence(sequence))

    def score_binding_site(self, sequence: str, structure: str, site_positions: list[int]) -> float:
        """Score a candidate binding site using geometry, conservation, and accessibility."""
        sequence = _normalize_sequence(sequence)
        if not site_positions:
            return 0.0
        family, family_similarity = self._best_family_match(sequence)
        site_positions = sorted(pos for pos in site_positions if 0 <= pos < len(sequence))
        if not site_positions:
            return 0.0
        site_sequence = "".join(sequence[pos] for pos in site_positions)
        loop_bonus = 0.0
        for segment in _find_unpaired_segments(structure):
            overlap = len(set(site_positions) & set(segment["positions"]))
            if overlap:
                ideal = 1.0 - min(1.0, abs(segment["size"] - 5) / 6.0)
                kind_bonus = {"hairpin": 1.0, "junction": 1.1, "bulge": 0.8, "external": 0.4}[segment["loop_type"]]
                loop_bonus = max(loop_bonus, ideal * kind_bonus)
        conserved = set(family["ligand_binding_site_positions"])
        conservation = len(conserved & set(site_positions)) / max(1, len(site_positions))
        sequence_bias = (site_sequence.count("A") + site_sequence.count("G")) / max(1, len(site_sequence))
        accessibility = np.mean([1.0 if structure[pos] == "." else 0.0 for pos in site_positions])
        local_context = np.mean([
            1.0
            if all(0 <= nbr < len(structure) and structure[nbr] == "." for nbr in (pos - 1, pos, pos + 1) if 0 <= nbr < len(structure))
            else 0.0
            for pos in site_positions
        ])
        raw_score = 2.0 * loop_bonus + 1.8 * conservation + 1.0 * accessibility + 0.8 * sequence_bias + 0.6 * local_context + 0.5 * family_similarity
        return float(expit(raw_score - 2.4))

    def predict_binding_sites(self, sequence: str, structure: str) -> list[dict[str, Any]]:
        """Identify loop and junction regions likely to host ligand binding."""
        sequence = _normalize_sequence(sequence)
        family, _ = self._best_family_match(sequence)
        candidates: list[dict[str, Any]] = []
        for segment in _find_unpaired_segments(structure):
            if segment["size"] < 3:
                continue
            score = self.score_binding_site(sequence, structure, segment["positions"])
            candidates.append(
                {
                    "site_positions": segment["positions"],
                    "loop_type": segment["loop_type"],
                    "score": score,
                    "sequence": sequence[segment["start"] : segment["end"] + 1],
                    "overlaps_family_binding_site": bool(set(segment["positions"]) & set(family["ligand_binding_site_positions"])),
                }
            )
        binding_positions = [pos for pos in family["ligand_binding_site_positions"] if pos < len(sequence)]
        if binding_positions:
            candidates.append(
                {
                    "site_positions": binding_positions,
                    "loop_type": "annotated_consensus",
                    "score": self.score_binding_site(sequence, structure, binding_positions),
                    "sequence": "".join(sequence[pos] for pos in binding_positions),
                    "overlaps_family_binding_site": True,
                }
            )
        unique: dict[tuple[int, ...], dict[str, Any]] = {}
        for candidate in sorted(candidates, key=lambda item: item["score"], reverse=True):
            key = tuple(candidate["site_positions"])
            unique.setdefault(key, candidate)
        return list(unique.values())[:10]


class ExpressionPlatformAnalyzer:
    """Analyze riboswitch expression platforms and regulatory outcomes."""

    def __init__(self, predictor: Optional[StructuralSwitchPredictor] = None) -> None:
        self.predictor = predictor

    @staticmethod
    def _find_shine_dalgarno(sequence: str) -> Optional[list[int]]:
        motifs = ["AGGAGG", "GGAGG", "AGGA"]
        for motif in motifs:
            position = sequence.find(motif)
            if position >= 0:
                return list(range(position, position + len(motif)))
        return None

    @staticmethod
    def _terminator_score(sequence: str, structure: str, aptamer_end: int) -> dict[str, Any]:
        best = {"score": 0.0, "positions": [], "stem_length": 0, "u_run": 0, "structure_fragment": ""}
        for stem in _contiguous_stems(structure, start=aptamer_end):
            tail_start = stem["right_end"] + 1
            u_run = 0
            while tail_start + u_run < len(sequence) and sequence[tail_start + u_run] == "U":
                u_run += 1
            gc_fraction = np.mean([
                1.0 if {sequence[left], sequence[right]} == {"G", "C"} else 0.0 for left, right in stem["pairs"]
            ])
            score = float(expit(0.8 * (stem["stem_length"] - 4) + 1.4 * gc_fraction + 0.5 * (u_run - 3)))
            if score > best["score"]:
                best = {
                    "score": score,
                    "positions": list(range(stem["left_start"], stem["right_end"] + 1)),
                    "stem_length": stem["stem_length"],
                    "u_run": u_run,
                    "structure_fragment": structure[stem["left_start"] : stem["right_end"] + 1],
                }
        return best

    @staticmethod
    def _anti_terminator_score(sequence: str, structure: str, aptamer_end: int, terminator_positions: set[int]) -> dict[str, Any]:
        downstream_segments = _find_unpaired_segments(structure[aptamer_end:])
        open_fraction = np.mean([1.0 if structure[pos] == "." else 0.0 for pos in terminator_positions]) if terminator_positions else 1.0
        best_loop = max((segment["size"] for segment in downstream_segments), default=0)
        score = float(expit(1.5 * open_fraction + 0.15 * best_loop - 1.0))
        return {"score": score, "open_terminator_fraction": float(open_fraction), "largest_loop": int(best_loop)}

    @staticmethod
    def _sd_state(structure: str, sd_positions: Optional[list[int]]) -> dict[str, Any]:
        if not sd_positions:
            return {"positions": [], "paired_fraction": 0.0, "sequestered": False}
        paired_fraction = np.mean([1.0 if structure[pos] != "." else 0.0 for pos in sd_positions if pos < len(structure)])
        return {
            "positions": sd_positions,
            "paired_fraction": float(paired_fraction),
            "sequestered": bool(paired_fraction >= 0.5),
        }

    def predict_regulatory_type(self, sequence: str, structure: str) -> str:
        """Classify the structure as transcriptional or translational regulation."""
        sequence = _normalize_sequence(sequence)
        aptamer_end = int(0.65 * len(sequence))
        terminator = self._terminator_score(sequence, structure, aptamer_end)
        sd_positions = self._find_shine_dalgarno(sequence)
        sd_state = self._sd_state(structure, sd_positions)
        return "transcriptional" if terminator["score"] >= max(0.55, sd_state["paired_fraction"]) else "translational"

    def analyze_regulation(self, sequence: str, aptamer_end: int, sd_positions: Optional[list[int]] = None) -> dict[str, Any]:
        """Analyze transcriptional and translational regulation from the ensemble."""
        sequence = _normalize_sequence(sequence)
        predictor = self.predictor or StructuralSwitchPredictor()
        sd_positions = sd_positions or self._find_shine_dalgarno(sequence)
        conformations = predictor.predict_conformations(sequence, n_suboptimal=8, energy_range=6.0)
        terminator_prob = 0.0
        anti_prob = 0.0
        sd_off_prob = 0.0
        best_terminator = {"score": 0.0, "positions": []}
        best_antiterminator = {"score": 0.0}
        for conformation in conformations:
            probability = conformation["ensemble_probability"]
            structure = conformation["structure"]
            terminator = self._terminator_score(sequence, structure, aptamer_end)
            anti = self._anti_terminator_score(sequence, structure, aptamer_end, set(terminator["positions"]))
            sd_state = self._sd_state(structure, sd_positions)
            terminator_prob += probability * terminator["score"]
            anti_prob += probability * anti["score"]
            sd_off_prob += probability * sd_state["paired_fraction"]
            if terminator["score"] > best_terminator["score"]:
                best_terminator = terminator
            if anti["score"] > best_antiterminator["score"]:
                best_antiterminator = anti
        regulatory_type = "transcriptional" if terminator_prob >= sd_off_prob else "translational"
        return {
            "regulatory_type": regulatory_type,
            "terminator": best_terminator,
            "anti_terminator": best_antiterminator,
            "shine_dalgarno": {
                "positions": sd_positions or [],
                "sequestration_probability": float(sd_off_prob),
            },
            "outcome_probabilities": {
                "transcription_off": float(terminator_prob),
                "transcription_on": float(max(0.0, 1.0 - terminator_prob)),
                "translation_off": float(sd_off_prob),
                "translation_on": float(max(0.0, 1.0 - sd_off_prob)),
                "anti_terminator": float(anti_prob),
            },
            "ensemble_size": len(conformations),
        }


class FunctionalMotifScanner:
    """Scan structures for common functional RNA motifs."""

    def score_motif(self, sequence: str, structure: str, motif_type: str, position: int) -> float:
        """Score a motif instance using sequence and structural context."""
        sequence = _normalize_sequence(sequence)
        motif_type = motif_type.lower()
        if position < 0 or position >= len(sequence):
            return 0.0
        loop_segments = _find_unpaired_segments(structure)
        segment = next((seg for seg in loop_segments if seg["start"] <= position <= seg["end"]), None)
        if motif_type == "gnra tetraloop" and segment and segment["size"] == 4:
            motif = sequence[segment["start"] : segment["end"] + 1]
            score = float(motif[0] == "G" and motif[-1] == "A" and motif[2] in {"A", "G"})
            return 0.6 + 0.4 * score
        if motif_type == "kink-turn":
            window = sequence[max(0, position - 2) : min(len(sequence), position + 5)]
            asymmetric = any(seg["loop_type"] in {"bulge", "junction"} and seg["start"] <= position <= seg["end"] for seg in loop_segments)
            ga_rich = "GA" in window or "AG" in window
            return float(expit(1.2 * float(asymmetric) + 1.0 * float(ga_rich) - 0.8))
        if motif_type == "sarcin-ricin loop":
            window = sequence[position : min(len(sequence), position + 7)]
            return float(expit(1.8 * float("GAGA" in window) + 0.8 * float(structure[position] == ".") - 1.0))
        if motif_type == "t-loop":
            motif = sequence[position : min(len(sequence), position + 5)]
            good = len(motif) == 5 and motif[0] == "U" and motif[2] in {"A", "G"} and motif[3] == "A"
            return 0.35 + 0.65 * float(good)
        if motif_type == "a-minor interaction site":
            window = sequence[position : min(len(sequence), position + 4)]
            a_fraction = window.count("A") / max(1, len(window))
            paired_neighbors = [idx for idx in range(max(0, position - 1), min(len(structure), position + len(window) + 1)) if structure[idx] in "()"]
            return float(expit(1.5 * a_fraction + 0.2 * len(paired_neighbors) - 1.0))
        return 0.0

    def scan_motifs(self, sequence: str, structure: str) -> list[dict[str, Any]]:
        """Scan for GNRA, kink-turn, sarcin-ricin loop, T-loop, and A-minor motifs."""
        sequence = _normalize_sequence(sequence)
        motifs: list[dict[str, Any]] = []
        for segment in _find_unpaired_segments(structure):
            start = segment["start"]
            if segment["loop_type"] == "hairpin" and segment["size"] == 4:
                score = self.score_motif(sequence, structure, "GNRA tetraloop", start)
                if score >= 0.6:
                    motifs.append({"motif_type": "GNRA tetraloop", "position": start, "score": score, "sequence": sequence[start : start + 4]})
            if segment["loop_type"] in {"hairpin", "junction"} and segment["size"] >= 4:
                score = self.score_motif(sequence, structure, "sarcin-ricin loop", start)
                if score >= 0.55:
                    motifs.append({"motif_type": "sarcin-ricin loop", "position": start, "score": score, "sequence": sequence[start : segment["end"] + 1]})
            if segment["loop_type"] == "hairpin" and segment["size"] >= 5:
                score = self.score_motif(sequence, structure, "T-loop", start)
                if score >= 0.55:
                    motifs.append({"motif_type": "T-loop", "position": start, "score": score, "sequence": sequence[start : start + 5]})
            score = self.score_motif(sequence, structure, "kink-turn", start)
            if score >= 0.55:
                motifs.append({"motif_type": "kink-turn", "position": start, "score": score, "sequence": sequence[start : segment["end"] + 1]})
            score = self.score_motif(sequence, structure, "A-minor interaction site", start)
            if score >= 0.55:
                motifs.append({"motif_type": "A-minor interaction site", "position": start, "score": score, "sequence": sequence[start : min(len(sequence), start + 4)]})
        motifs.sort(key=lambda item: item["score"], reverse=True)
        return motifs


class RiboswitchBenchmark:
    """Benchmark heuristic structure prediction across known riboswitch families."""

    def __init__(self, database: Optional[RiboswitchDatabase] = None) -> None:
        self.database = database or RiboswitchDatabase()

    def run_benchmark(self, predictor: StructuralSwitchPredictor) -> dict[str, Any]:
        """Benchmark predicted structures and switching behavior on database families."""
        per_family: dict[str, dict[str, Any]] = {}
        accuracies: list[float] = []
        switch_hits: list[float] = []
        pair_f1_values: list[float] = []
        for family_name in self.database.list_families()[:5]:
            family = self.database.get_family(family_name)
            conformations = predictor.predict_conformations(family["consensus_sequence"], n_suboptimal=10, energy_range=8.0)
            scored = []
            for conformation in conformations:
                accuracy = _structure_accuracy(conformation["structure"], family["consensus_structure"])
                metrics = _pair_metrics(conformation["structure"], family["consensus_structure"])
                scored.append((accuracy, metrics, conformation))
            best_accuracy, metrics, best = max(scored, key=lambda item: item[0])
            expression_candidates = [item[2] for item in scored if item[2]["conformation_type"] == "expression_platform"]
            expression = max(expression_candidates, key=lambda item: item["ensemble_probability"], default=None)
            switch_region = predictor.identify_switch_region(best, expression) if expression else []
            switch_detected = bool(expression and len(switch_region) >= max(6, len(family["consensus_sequence"]) // 12))
            accuracies.append(best_accuracy)
            pair_f1_values.append(metrics["f1"])
            switch_hits.append(float(switch_detected))
            per_family[family_name] = {
                "best_structure": best["structure"],
                "best_energy": best["energy"],
                "best_conformation_type": best["conformation_type"],
                "structure_prediction_accuracy": float(best_accuracy),
                "pair_metrics": metrics,
                "switch_detected": switch_detected,
                "switch_region_size": len(switch_region),
                "n_conformations": len(conformations),
            }
        return {
            "summary": {
                "n_cases": len(per_family),
                "structure_prediction_accuracy": float(np.mean(accuracies)) if accuracies else 0.0,
                "pair_f1": float(np.mean(pair_f1_values)) if pair_f1_values else 0.0,
                "switch_detection_rate": float(np.mean(switch_hits)) if switch_hits else 0.0,
            },
            "per_family": per_family,
        }


__all__ = [
    "RiboswitchDatabase",
    "StructuralSwitchPredictor",
    "LigandBindingPredictor",
    "ExpressionPlatformAnalyzer",
    "FunctionalMotifScanner",
    "RiboswitchBenchmark",
]
