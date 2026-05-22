"""RNA pseudoknot structure prediction algorithms.

This module implements a compact toolkit for pseudoknot-aware RNA secondary
structure prediction.  It includes:

* :class:`PseudoknotDetector` for crossing-pair detection and coarse taxonomy.
* :class:`AkutsuDP` for an exact dynamic-programming search over simple H-type
  pseudoknots using five DP matrices ``W``, ``V``, ``WK``, ``VK1``, and ``VK2``.
* :class:`HeuristicPseudoknot` for a faster two-pass predictor.
* :class:`PseudoknotScorer` for pseudoknot-specific free-energy terms.
* :class:`IterativeRelaxation` for Lagrangian relaxation over crossing pairs.
* :class:`BenchmarkPseudoknot` for lightweight runtime and accuracy comparison.

All energies are expressed in kcal/mol.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import time
from typing import Any, Dict, Iterable, List, Optional, Protocol, Sequence, Tuple

import numpy as np

Pair = Tuple[int, int]
INF = float("inf")
_CANONICAL_PAIRS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}


class TurnerParameters(Protocol):
    """Protocol for external Turner-like energy models."""

    def get_pair_energy(self, left: str, right: str) -> float:
        ...

    def get_stack_energy(self, left1: str, right1: str, left2: str, right2: str) -> float:
        ...


@dataclass(frozen=True)
class StemCandidate:
    """Contiguous stem candidate used by exact and heuristic predictors."""

    pairs: Tuple[Pair, ...]
    energy: float

    @property
    def length(self) -> int:
        return len(self.pairs)

    @property
    def left_start(self) -> int:
        return self.pairs[0][0]

    @property
    def left_end(self) -> int:
        return self.pairs[-1][0]

    @property
    def right_start(self) -> int:
        return self.pairs[-1][1]

    @property
    def right_end(self) -> int:
        return self.pairs[0][1]

    @property
    def span(self) -> Tuple[int, int]:
        return (self.left_start, self.right_end)


@dataclass(frozen=True)
class BenchmarkEntry:
    """Reference benchmark entry for pseudoknot-aware evaluation."""

    name: str
    sequence: str
    reference_pairs: Tuple[Pair, ...]
    pseudoknot_type: str
    note: str


def _normalize_pairs(pairs: Iterable[Pair]) -> List[Pair]:
    normalized = {tuple(sorted((int(i), int(j)))) for i, j in pairs if i != j}
    return sorted(normalized)


def _crosses(first: Pair, second: Pair) -> bool:
    i, j = first
    k, l = second
    return (i < k < j < l) or (k < i < l < j)


def _shares_index(first: Pair, second: Pair) -> bool:
    return len({first[0], first[1], second[0], second[1]}) < 4


def _conflicts(first: Pair, second: Pair) -> bool:
    return _shares_index(first, second) or _crosses(first, second)


def _pairs_to_index_set(pairs: Iterable[Pair]) -> set[int]:
    indices: set[int] = set()
    for i, j in pairs:
        indices.add(i)
        indices.add(j)
    return indices


def _can_pair(left: str, right: str) -> bool:
    return (left.upper().replace("T", "U"), right.upper().replace("T", "U")) in _CANONICAL_PAIRS


def _interval_value(matrix: np.ndarray, i: int, j: int) -> float:
    if i > j:
        return 0.0
    return float(matrix[i, j])


def _allowed(allowed_positions: Optional[set[int]], *positions: int) -> bool:
    return allowed_positions is None or all(position in allowed_positions for position in positions)


def _group_consecutive_stems(pairs: Sequence[Pair]) -> List[List[Pair]]:
    normalized = _normalize_pairs(pairs)
    if not normalized:
        return []
    stems: List[List[Pair]] = [[normalized[0]]]
    for pair in normalized[1:]:
        prev = stems[-1][-1]
        if pair[0] == prev[0] + 1 and pair[1] == prev[1] - 1:
            stems[-1].append(pair)
        else:
            stems.append([pair])
    return stems


def _valid_htype_configuration(stem1: StemCandidate, stem2: StemCandidate) -> bool:
    if _pairs_to_index_set(stem1.pairs) & _pairs_to_index_set(stem2.pairs):
        return False
    left, right = (stem1, stem2) if stem1.left_start <= stem2.left_start else (stem2, stem1)
    return (
        left.left_start < right.left_start < left.right_end < right.right_end
        and left.left_end < right.left_start
        and right.left_end < left.right_start
        and left.right_end < right.right_start
    )


class PseudoknotDetector:
    """Detect and classify crossing base-pair topologies."""

    def is_pseudoknotted(self, pairs: Sequence[Pair]) -> bool:
        normalized = _normalize_pairs(pairs)
        return any(_crosses(left, right) for index, left in enumerate(normalized) for right in normalized[index + 1 :])

    def classify_pseudoknot(self, pairs: Sequence[Pair]) -> str:
        normalized = _normalize_pairs(pairs)
        if not self.is_pseudoknotted(normalized):
            return "none"
        stems = self._group_stems(normalized)
        if len(stems) <= 2:
            return "H-type"

        edges = self._stem_crossing_edges(stems)
        degrees = [0 for _ in stems]
        for left, right in edges:
            degrees[left] += 1
            degrees[right] += 1
        if any(degree >= 2 for degree in degrees):
            return "kissing hairpin"
        if len(self.decompose_pseudoknot(normalized)) > 2 or self._connected_components(stems, edges) > 1:
            return "recursive"
        return "recursive"

    def decompose_pseudoknot(self, pairs: Sequence[Pair]) -> List[List[Pair]]:
        layers: List[List[Pair]] = []
        for pair in _normalize_pairs(pairs):
            placed = False
            for layer in layers:
                if all(not _crosses(pair, existing) for existing in layer):
                    layer.append(pair)
                    placed = True
                    break
            if not placed:
                layers.append([pair])
        return [sorted(layer) for layer in layers]

    def pseudoknotted_pairs(self, pairs: Sequence[Pair]) -> List[Pair]:
        normalized = _normalize_pairs(pairs)
        return [pair for pair in normalized if any(_crosses(pair, other) for other in normalized if pair != other)]

    def _group_stems(self, pairs: Sequence[Pair]) -> List[List[Pair]]:
        return _group_consecutive_stems(pairs)

    def _stem_crosses(self, first: Sequence[Pair], second: Sequence[Pair]) -> bool:
        return any(_crosses(left, right) for left in first for right in second)

    def _stem_crossing_edges(self, stems: Sequence[Sequence[Pair]]) -> List[Tuple[int, int]]:
        edges: List[Tuple[int, int]] = []
        for i, left in enumerate(stems):
            for j in range(i + 1, len(stems)):
                if self._stem_crosses(left, stems[j]):
                    edges.append((i, j))
        return edges

    def _connected_components(self, stems: Sequence[Sequence[Pair]], edges: Sequence[Tuple[int, int]]) -> int:
        if not stems:
            return 0
        adjacency: Dict[int, set[int]] = {index: set() for index in range(len(stems))}
        for left, right in edges:
            adjacency[left].add(right)
            adjacency[right].add(left)
        seen: set[int] = set()
        components = 0
        for start in range(len(stems)):
            if start in seen:
                continue
            components += 1
            stack = [start]
            seen.add(start)
            while stack:
                node = stack.pop()
                for neighbor in adjacency[node]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        stack.append(neighbor)
        return components


class PseudoknotScorer:
    """Energy model for pseudoknot stems, loop entropy, and coaxial stacking."""

    DEFAULT_PAIR_ENERGIES: Dict[Tuple[str, str], float] = {
        ("G", "C"): -3.3,
        ("C", "G"): -3.3,
        ("A", "U"): -2.1,
        ("U", "A"): -2.1,
        ("G", "U"): -1.1,
        ("U", "G"): -1.1,
    }

    def __init__(
        self,
        turner_params: Optional[TurnerParameters] = None,
        coaxial_bonus: float = -0.6,
        entropy_base: float = 3.2,
        loop_penalty: float = 0.32,
        pseudoknot_stack_bonus: float = -0.45,
    ) -> None:
        self.turner_params = turner_params
        self.coaxial_bonus = coaxial_bonus
        self.entropy_base = entropy_base
        self.loop_penalty = loop_penalty
        self.pseudoknot_stack_bonus = pseudoknot_stack_bonus
        self.detector = PseudoknotDetector()

    def pair_energy(self, left: str, right: str) -> float:
        left = left.upper().replace("T", "U")
        right = right.upper().replace("T", "U")
        if self.turner_params is not None:
            for attr in ("get_pair_energy", "pair_energy"):
                func = getattr(self.turner_params, attr, None)
                if callable(func):
                    try:
                        return float(func(left, right))
                    except TypeError:
                        return float(func((left, right)))
            energy_map = getattr(self.turner_params, "pair_energies", None)
            if isinstance(energy_map, dict) and (left, right) in energy_map:
                return float(energy_map[(left, right)])
        return self.DEFAULT_PAIR_ENERGIES.get((left, right), 0.6)

    def stack_energy(self, pair1: Pair, pair2: Pair, sequence: str) -> float:
        i, j = pair1
        k, l = pair2
        sequence = sequence.upper().replace("T", "U")
        if self.turner_params is not None:
            for attr in ("get_stack_energy", "stack_energy"):
                func = getattr(self.turner_params, attr, None)
                if callable(func):
                    return float(func(sequence[i], sequence[j], sequence[k], sequence[l]))
        if abs(i - k) == 1 and abs(j - l) == 1:
            return 0.18 * (
                self.pair_energy(sequence[i], sequence[j]) + self.pair_energy(sequence[k], sequence[l])
            ) - 0.15
        return 0.0

    def score_stem(self, stem_pairs: Sequence[Pair], sequence: str) -> float:
        ordered = list(sorted(stem_pairs))
        if not ordered:
            return 0.0
        total = sum(self.pair_energy(sequence[i], sequence[j]) for i, j in ordered)
        total += sum(
            self.stack_energy(ordered[index], ordered[index + 1], sequence) for index in range(len(ordered) - 1)
        )
        return float(total)

    def entropy_penalty(self, loop_lengths: Sequence[int]) -> float:
        positive = [max(0, loop) for loop in loop_lengths]
        if not positive:
            return 0.0
        linear = self.loop_penalty * float(sum(positive))
        logarithmic = 0.2 * float(sum(math.log1p(loop) for loop in positive))
        return self.entropy_base + linear + logarithmic

    def score_htype(self, stem1_pairs: Sequence[Pair], stem2_pairs: Sequence[Pair], sequence: str) -> float:
        stem1 = StemCandidate(tuple(sorted(stem1_pairs)), self.score_stem(stem1_pairs, sequence))
        stem2 = StemCandidate(tuple(sorted(stem2_pairs)), self.score_stem(stem2_pairs, sequence))
        if not _valid_htype_configuration(stem1, stem2):
            return stem1.energy + stem2.energy + 8.0
        if stem1.left_start > stem2.left_start:
            stem1, stem2 = stem2, stem1
        loop_lengths = [
            max(0, stem2.left_start - stem1.left_end - 1),
            max(0, stem1.right_start - stem2.left_end - 1),
            max(0, stem2.right_start - stem1.right_end - 1),
        ]
        coaxial = self._coaxial_bonus(stem1, stem2)
        return float(
            stem1.energy
            + stem2.energy
            + 2.0 * self.pseudoknot_stack_bonus
            + coaxial
            + self.entropy_penalty(loop_lengths)
        )

    def score_kissing(
        self,
        stem1: Sequence[Pair],
        stem2: Sequence[Pair],
        stem3: Sequence[Pair],
        sequence: str,
    ) -> float:
        stems = [
            StemCandidate(tuple(sorted(stem1)), self.score_stem(stem1, sequence)),
            StemCandidate(tuple(sorted(stem2)), self.score_stem(stem2, sequence)),
            StemCandidate(tuple(sorted(stem3)), self.score_stem(stem3, sequence)),
        ]
        stems.sort(key=lambda stem: stem.left_start)
        loop_lengths = [
            max(0, stems[1].left_start - stems[0].left_end - 1),
            max(0, stems[0].right_start - stems[1].left_end - 1),
            max(0, stems[2].left_start - stems[1].left_end - 1),
            max(0, stems[1].right_end - stems[0].right_end - 1),
            max(0, stems[2].right_start - stems[1].right_end - 1),
        ]
        coaxial = self._coaxial_bonus(stems[0], stems[1]) + self._coaxial_bonus(stems[1], stems[2])
        return float(
            sum(stem.energy for stem in stems)
            + 3.0 * self.pseudoknot_stack_bonus
            + coaxial
            + self.entropy_penalty(loop_lengths)
            + 1.2
        )

    def score_structure(self, pairs: Sequence[Pair], sequence: str) -> float:
        normalized = _normalize_pairs(pairs)
        if not normalized:
            return 0.0
        total = sum(self.pair_energy(sequence[i], sequence[j]) for i, j in normalized)
        detector = self.detector
        if detector.is_pseudoknotted(normalized):
            stems = detector._group_stems(detector.pseudoknotted_pairs(normalized))
            kind = detector.classify_pseudoknot(normalized)
            if kind == "H-type" and len(stems) >= 2:
                total += self.score_htype(stems[0], stems[1], sequence) - sum(
                    self.score_stem(stem, sequence) for stem in stems[:2]
                )
            elif kind == "kissing hairpin" and len(stems) >= 3:
                total += self.score_kissing(stems[0], stems[1], stems[2], sequence) - sum(
                    self.score_stem(stem, sequence) for stem in stems[:3]
                )
            else:
                total += 1.8
        return float(total)

    def _coaxial_bonus(self, stem1: StemCandidate, stem2: StemCandidate) -> float:
        close_left = abs(stem2.left_start - stem1.left_end - 1) <= 1
        close_right = abs(stem2.right_start - stem1.right_end - 1) <= 1
        if close_left or close_right:
            return self.coaxial_bonus
        return 0.5 * self.coaxial_bonus


class AkutsuDP:
    """Exact O(n^4) DP for simple H-type pseudoknots.

    The implementation explicitly maintains the five matrices commonly used in a
    compact Akutsu-style decomposition:

    * ``W``: best energy over an interval.
    * ``V``: best nested stem anchored at interval boundaries.
    * ``WK``: best H-type pseudoknot spanning an interval.
    * ``VK1``: first crossing-stem candidate matrix.
    * ``VK2``: second crossing-stem candidate matrix.
    """

    def __init__(
        self,
        sequence: str,
        turner_params: Optional[TurnerParameters] = None,
        min_stem: int = 2,
        max_stem: Optional[int] = 6,
        scorer: Optional[PseudoknotScorer] = None,
    ) -> None:
        self.sequence = sequence.upper().replace("T", "U")
        self.turner_params = turner_params
        self.min_stem = min_stem
        self.max_stem = max_stem
        self.scorer = scorer or PseudoknotScorer(turner_params=turner_params)
        self.detector = PseudoknotDetector()
        n = len(self.sequence)
        self.W = np.full((n, n), INF, dtype=float)
        self.V = np.full((n, n), INF, dtype=float)
        self.WK = np.full((n, n), INF, dtype=float)
        self.VK1 = np.full((n, n), INF, dtype=float)
        self.VK2 = np.full((n, n), INF, dtype=float)
        self.w_trace: Dict[Tuple[int, int], Tuple[str, Any]] = {}
        self.v_trace: Dict[Tuple[int, int], StemCandidate] = {}
        self.wk_trace: Dict[Tuple[int, int], Dict[str, Any]] = {}
        self.vk1_choice: Dict[Tuple[int, int], StemCandidate] = {}
        self.vk2_choice: Dict[Tuple[int, int], StemCandidate] = {}
        self.endpoint_candidates: Dict[Tuple[int, int], List[StemCandidate]] = {}
        self._computed = False

    def compute(self) -> float:
        if self._computed:
            return float(self.W[0, len(self.sequence) - 1]) if self.sequence else 0.0
        n = len(self.sequence)
        if n == 0:
            self._computed = True
            return 0.0
        self._prepare_candidates()
        for i in range(n):
            self.W[i, i] = 0.0
            self.w_trace[(i, i)] = ("empty", None)
        for span in range(2, n + 1):
            for i in range(0, n - span + 1):
                j = i + span - 1
                self._compute_v(i, j)
                self._compute_wk(i, j)
                best = 0.0
                action: Tuple[str, Any] = ("empty", None)
                for split in range(i, j):
                    value = _interval_value(self.W, i, split) + _interval_value(self.W, split + 1, j)
                    if value < best:
                        best = value
                        action = ("split", split)
                if self.V[i, j] < best:
                    best = float(self.V[i, j])
                    action = ("V", self.v_trace[(i, j)])
                if self.WK[i, j] < best:
                    best = float(self.WK[i, j])
                    action = ("WK", self.wk_trace[(i, j)])
                self.W[i, j] = best
                self.w_trace[(i, j)] = action
        self._computed = True
        return float(self.W[0, n - 1])

    def predict_structure(self) -> Tuple[float, List[Pair]]:
        energy = self.compute()
        if not self.sequence:
            return energy, []
        pairs = self._traceback_interval(0, len(self.sequence) - 1)
        return energy, _normalize_pairs(pairs)

    def _prepare_candidates(self) -> None:
        candidates = _enumerate_stem_candidates(
            self.sequence,
            self.scorer,
            self.min_stem,
            self.max_stem,
        )
        for candidate in candidates:
            self.endpoint_candidates.setdefault(candidate.span, []).append(candidate)
        for (i, j), endpoint_candidates in self.endpoint_candidates.items():
            endpoint_candidates.sort(key=lambda candidate: (candidate.energy, -candidate.length))
            best = endpoint_candidates[0]
            self.VK1[i, j] = best.energy
            self.VK2[i, j] = best.energy
            self.vk1_choice[(i, j)] = best
            self.vk2_choice[(i, j)] = best

    def _compute_v(self, i: int, j: int) -> None:
        best = INF
        best_candidate: Optional[StemCandidate] = None
        for candidate in self.endpoint_candidates.get((i, j), []):
            inside = _interval_value(self.W, candidate.left_end + 1, candidate.right_start - 1)
            value = candidate.energy + inside
            if value < best:
                best = value
                best_candidate = candidate
        self.V[i, j] = best
        if best_candidate is not None:
            self.v_trace[(i, j)] = best_candidate

    def _compute_wk(self, i: int, j: int) -> None:
        best = INF
        best_choice: Optional[Dict[str, Any]] = None
        for right1 in range(i + 2 * self.min_stem - 1, j):
            stem1 = self.vk1_choice.get((i, right1))
            if stem1 is None:
                continue
            for left2 in range(i + 1, right1):
                stem2 = self.vk2_choice.get((left2, j))
                if stem2 is None:
                    continue
                if not _valid_htype_configuration(stem1, stem2):
                    continue
                if stem1.left_start > stem2.left_start:
                    stem1, stem2 = stem2, stem1
                left_loop = (stem1.left_end + 1, stem2.left_start - 1)
                middle_loop = (stem2.left_end + 1, stem1.right_start - 1)
                right_loop = (stem1.right_end + 1, stem2.right_start - 1)
                value = self.scorer.score_htype(stem1.pairs, stem2.pairs, self.sequence)
                value += _interval_value(self.W, *left_loop)
                value += _interval_value(self.W, *middle_loop)
                value += _interval_value(self.W, *right_loop)
                if value < best:
                    best = value
                    best_choice = {
                        "stem1": stem1,
                        "stem2": stem2,
                        "left_loop": left_loop,
                        "middle_loop": middle_loop,
                        "right_loop": right_loop,
                    }
        self.WK[i, j] = best
        if best_choice is not None:
            self.wk_trace[(i, j)] = best_choice

    def _traceback_interval(self, i: int, j: int) -> List[Pair]:
        if i >= j:
            return []
        action, payload = self.w_trace.get((i, j), ("empty", None))
        if action == "split":
            split = int(payload)
            return self._traceback_interval(i, split) + self._traceback_interval(split + 1, j)
        if action == "V":
            candidate = payload
            inside = self._traceback_interval(candidate.left_end + 1, candidate.right_start - 1)
            return list(candidate.pairs) + inside
        if action == "WK":
            choice = payload
            stem1 = list(choice["stem1"].pairs)
            stem2 = list(choice["stem2"].pairs)
            left = self._traceback_interval(*choice["left_loop"])
            middle = self._traceback_interval(*choice["middle_loop"])
            right = self._traceback_interval(*choice["right_loop"])
            return stem1 + stem2 + left + middle + right
        return []


class HeuristicPseudoknot:
    """Efficient two-pass pseudoknot prediction.

    Pass 1 computes a standard nested MFE structure in O(n^3).  Pass 2 searches
    only the unpaired regions for crossing H-type stem pairs, scores the free
    energy gain, and greedily adds non-conflicting pseudoknots in O(n^2) over the
    reduced candidate pool.
    """

    def __init__(
        self,
        sequence: str,
        turner_params: Optional[TurnerParameters] = None,
        min_stem: int = 2,
        max_stem: Optional[int] = 5,
        min_gain: float = 0.25,
        scorer: Optional[PseudoknotScorer] = None,
    ) -> None:
        self.sequence = sequence.upper().replace("T", "U")
        self.turner_params = turner_params
        self.min_stem = min_stem
        self.max_stem = max_stem
        self.min_gain = min_gain
        self.scorer = scorer or PseudoknotScorer(turner_params=turner_params)
        self.detector = PseudoknotDetector()

    def predict_structure(self, regions: Optional[Sequence[Tuple[int, int]]] = None) -> Tuple[float, List[Pair]]:
        nested_energy, nested_pairs = self._nested_mfe()
        candidate_pseudoknots = self._candidate_pseudoknots(nested_pairs, regions)
        chosen_pairs: List[Pair] = []
        chosen_energy = 0.0
        occupied = _pairs_to_index_set(chosen_pairs)
        for candidate in sorted(candidate_pseudoknots, key=lambda item: (-item["gain"], item["energy"])):
            candidate_pairs = list(candidate["stem1"].pairs) + list(candidate["stem2"].pairs)
            candidate_indices = _pairs_to_index_set(candidate_pairs)
            if occupied & candidate_indices:
                continue
            if any(_conflicts(left, right) for left in chosen_pairs for right in candidate_pairs):
                continue
            chosen_pairs.extend(candidate_pairs)
            chosen_energy += float(candidate["energy"])
            occupied |= candidate_indices
        final_pairs = _normalize_pairs(nested_pairs + chosen_pairs)
        return nested_energy + chosen_energy, final_pairs

    def _nested_mfe(self) -> Tuple[float, List[Pair]]:
        n = len(self.sequence)
        if n == 0:
            return 0.0, []
        matrix = np.zeros((n, n), dtype=float)
        trace: Dict[Tuple[int, int], Tuple[str, Any]] = {}
        endpoint_candidates: Dict[Tuple[int, int], List[StemCandidate]] = {}
        for candidate in _enumerate_stem_candidates(self.sequence, self.scorer, self.min_stem, self.max_stem):
            endpoint_candidates.setdefault(candidate.span, []).append(candidate)
        for span in range(2, n + 1):
            for i in range(0, n - span + 1):
                j = i + span - 1
                best = 0.0
                action: Tuple[str, Any] = ("empty", None)
                for split in range(i, j):
                    value = _interval_value(matrix, i, split) + _interval_value(matrix, split + 1, j)
                    if value < best:
                        best = value
                        action = ("split", split)
                for candidate in endpoint_candidates.get((i, j), []):
                    inside = _interval_value(matrix, candidate.left_end + 1, candidate.right_start - 1)
                    value = candidate.energy + inside
                    if value < best:
                        best = value
                        action = ("stem", candidate)
                matrix[i, j] = best
                trace[(i, j)] = action
        return float(matrix[0, n - 1]), self._trace_nested(0, n - 1, trace)

    def _trace_nested(self, i: int, j: int, trace: Dict[Tuple[int, int], Tuple[str, Any]]) -> List[Pair]:
        if i >= j:
            return []
        action, payload = trace.get((i, j), ("empty", None))
        if action == "split":
            split = int(payload)
            return self._trace_nested(i, split, trace) + self._trace_nested(split + 1, j, trace)
        if action == "stem":
            candidate = payload
            return list(candidate.pairs) + self._trace_nested(candidate.left_end + 1, candidate.right_start - 1, trace)
        return []

    def _candidate_pseudoknots(
        self,
        nested_pairs: Sequence[Pair],
        regions: Optional[Sequence[Tuple[int, int]]] = None,
    ) -> List[Dict[str, Any]]:
        occupied = _pairs_to_index_set(nested_pairs)
        allowed = {index for index in range(len(self.sequence)) if index not in occupied}
        if regions is not None:
            region_positions: set[int] = set()
            for start, end in regions:
                region_positions.update(range(max(0, start), min(len(self.sequence), end + 1)))
            allowed &= region_positions
        stems = _enumerate_stem_candidates(
            self.sequence,
            self.scorer,
            self.min_stem,
            self.max_stem,
            allowed_positions=allowed,
        )
        candidates: List[Dict[str, Any]] = []
        for index, stem1 in enumerate(stems):
            for stem2 in stems[index + 1 :]:
                if not _valid_htype_configuration(stem1, stem2):
                    continue
                energy = self.scorer.score_htype(stem1.pairs, stem2.pairs, self.sequence)
                gain = -energy
                if gain < self.min_gain:
                    continue
                candidates.append({"stem1": stem1, "stem2": stem2, "energy": energy, "gain": gain})
        return candidates


class IterativeRelaxation:
    """Iterative relaxation for pseudoknot prediction with Lagrangian penalties."""

    def __init__(
        self,
        sequence: str,
        turner_params: Optional[TurnerParameters] = None,
        max_iterations: int = 25,
        crossing_penalty: float = 1.0,
        tolerance: float = 1e-3,
        scorer: Optional[PseudoknotScorer] = None,
    ) -> None:
        self.sequence = sequence.upper().replace("T", "U")
        self.turner_params = turner_params
        self.max_iterations = max_iterations
        self.crossing_penalty = crossing_penalty
        self.tolerance = tolerance
        self.scorer = scorer or PseudoknotScorer(turner_params=turner_params)
        self.detector = PseudoknotDetector()

    def predict_structure(self) -> Tuple[float, List[Pair]]:
        baseline_energy, baseline_pairs = HeuristicPseudoknot(
            self.sequence,
            turner_params=self.turner_params,
            scorer=self.scorer,
        )._nested_mfe()
        candidates = self._candidate_pairs()
        crossings = [(left, right) for index, left in enumerate(candidates) for right in candidates[index + 1 :] if _crosses(left, right)]
        multipliers: Dict[Tuple[Pair, Pair], float] = {tuple(sorted((left, right))): 0.0 for left, right in crossings}
        best_pairs = baseline_pairs
        best_energy = self.scorer.score_structure(best_pairs, self.sequence)
        previous_energy = baseline_energy
        for iteration in range(self.max_iterations):
            selected = self._solve_relaxed_problem(candidates, multipliers)
            feasible = self._repair_crossings(selected)
            energy = self.scorer.score_structure(feasible, self.sequence)
            if energy < best_energy - self.tolerance:
                best_energy = energy
                best_pairs = feasible
            violated = [edge for edge in multipliers if edge[0] in selected and edge[1] in selected]
            if not violated and abs(previous_energy - energy) < self.tolerance:
                break
            step = self.crossing_penalty / math.sqrt(iteration + 1.0)
            for edge in violated:
                multipliers[edge] += step
            previous_energy = energy
        return best_energy, _normalize_pairs(best_pairs)

    def _candidate_pairs(self) -> List[Pair]:
        candidates: List[Pair] = []
        for i in range(len(self.sequence)):
            for j in range(i + 4, len(self.sequence)):
                if _can_pair(self.sequence[i], self.sequence[j]):
                    candidates.append((i, j))
        return candidates

    def _solve_relaxed_problem(
        self,
        candidates: Sequence[Pair],
        multipliers: Dict[Tuple[Pair, Pair], float],
    ) -> List[Pair]:
        scored_candidates: List[Tuple[float, Pair]] = []
        for pair in candidates:
            reduced_cost = self.scorer.pair_energy(self.sequence[pair[0]], self.sequence[pair[1]])
            for edge, multiplier in multipliers.items():
                if pair in edge:
                    reduced_cost += multiplier
            scored_candidates.append((reduced_cost, pair))
        scored_candidates.sort(key=lambda item: item[0])
        chosen: List[Pair] = []
        occupied: set[int] = set()
        for score, pair in scored_candidates:
            if score >= 0.0:
                continue
            if pair[0] in occupied or pair[1] in occupied:
                continue
            chosen.append(pair)
            occupied.add(pair[0])
            occupied.add(pair[1])
        return chosen

    def _repair_crossings(self, pairs: Sequence[Pair]) -> List[Pair]:
        selected = list(_normalize_pairs(pairs))
        while True:
            crossing_pairs = [
                (left, right)
                for index, left in enumerate(selected)
                for right in selected[index + 1 :]
                if _crosses(left, right)
            ]
            if not crossing_pairs:
                return selected
            crossing_degree: Dict[Pair, int] = {pair: 0 for pair in selected}
            for left, right in crossing_pairs:
                crossing_degree[left] += 1
                crossing_degree[right] += 1
            worst = max(
                selected,
                key=lambda pair: (
                    crossing_degree[pair],
                    self.scorer.pair_energy(self.sequence[pair[0]], self.sequence[pair[1]]),
                ),
            )
            selected.remove(worst)


class BenchmarkPseudoknot:
    """Benchmark exact and heuristic pseudoknot predictors on motif panels."""

    def __init__(
        self,
        turner_params: Optional[TurnerParameters] = None,
        min_stem: int = 2,
        max_stem: Optional[int] = 5,
    ) -> None:
        self.turner_params = turner_params
        self.min_stem = min_stem
        self.max_stem = max_stem
        self.detector = PseudoknotDetector()
        self.scorer = PseudoknotScorer(turner_params=turner_params)
        self.benchmarks = self._build_benchmarks()

    def compare_algorithms(self) -> Dict[str, Any]:
        exact_rows: List[Dict[str, Any]] = []
        heuristic_rows: List[Dict[str, Any]] = []
        per_benchmark: List[Dict[str, Any]] = []
        for benchmark in self.benchmarks:
            exact_model = AkutsuDP(
                benchmark.sequence,
                turner_params=self.turner_params,
                min_stem=self.min_stem,
                max_stem=self.max_stem,
                scorer=self.scorer,
            )
            heuristic_model = HeuristicPseudoknot(
                benchmark.sequence,
                turner_params=self.turner_params,
                min_stem=self.min_stem,
                max_stem=self.max_stem,
                scorer=self.scorer,
            )

            start = time.perf_counter()
            exact_energy, exact_pairs = exact_model.predict_structure()
            exact_runtime = time.perf_counter() - start

            start = time.perf_counter()
            heuristic_energy, heuristic_pairs = heuristic_model.predict_structure()
            heuristic_runtime = time.perf_counter() - start

            exact_metrics = self._pseudoknot_metrics(exact_pairs, benchmark.reference_pairs)
            heuristic_metrics = self._pseudoknot_metrics(heuristic_pairs, benchmark.reference_pairs)
            exact_rows.append({"runtime_s": exact_runtime, **exact_metrics})
            heuristic_rows.append({"runtime_s": heuristic_runtime, **heuristic_metrics})
            per_benchmark.append(
                {
                    "name": benchmark.name,
                    "reference_type": benchmark.pseudoknot_type,
                    "sequence_length": len(benchmark.sequence),
                    "exact": {
                        "energy_kcal_mol": exact_energy,
                        "runtime_s": exact_runtime,
                        "pairs": exact_pairs,
                        "classification": self.detector.classify_pseudoknot(exact_pairs),
                        **exact_metrics,
                    },
                    "heuristic": {
                        "energy_kcal_mol": heuristic_energy,
                        "runtime_s": heuristic_runtime,
                        "pairs": heuristic_pairs,
                        "classification": self.detector.classify_pseudoknot(heuristic_pairs),
                        **heuristic_metrics,
                    },
                }
            )
        return {
            "benchmark_count": len(self.benchmarks),
            "benchmarks": per_benchmark,
            "exact_summary": self._summarize_rows(exact_rows),
            "heuristic_summary": self._summarize_rows(heuristic_rows),
        }

    def _pseudoknot_metrics(self, predicted: Sequence[Pair], reference: Sequence[Pair]) -> Dict[str, float]:
        predicted_pk = set(self.detector.pseudoknotted_pairs(predicted))
        reference_pk = set(self.detector.pseudoknotted_pairs(reference))
        true_positive = len(predicted_pk & reference_pk)
        false_positive = len(predicted_pk - reference_pk)
        false_negative = len(reference_pk - predicted_pk)
        sensitivity = true_positive / len(reference_pk) if reference_pk else 0.0
        ppv = true_positive / len(predicted_pk) if predicted_pk else 0.0
        f1 = (2.0 * sensitivity * ppv / (sensitivity + ppv)) if (sensitivity + ppv) else 0.0
        return {
            "tp": float(true_positive),
            "fp": float(false_positive),
            "fn": float(false_negative),
            "sensitivity": sensitivity,
            "ppv": ppv,
            "f1": f1,
        }

    def _summarize_rows(self, rows: Sequence[Dict[str, float]]) -> Dict[str, float]:
        if not rows:
            return {"runtime_s": 0.0, "sensitivity": 0.0, "ppv": 0.0, "f1": 0.0}
        return {
            "runtime_s": float(np.mean([row["runtime_s"] for row in rows])),
            "sensitivity": float(np.mean([row["sensitivity"] for row in rows])),
            "ppv": float(np.mean([row["ppv"] for row in rows])),
            "f1": float(np.mean([row["f1"] for row in rows])),
        }

    def _build_benchmarks(self) -> List[BenchmarkEntry]:
        return [
            BenchmarkEntry(
                name="tmRNA_PK1",
                sequence=self._build_sequence(26, [(0, 14), (1, 13), (2, 12), (6, 21), (7, 20), (8, 19)]),
                reference_pairs=((0, 14), (1, 13), (2, 12), (6, 21), (7, 20), (8, 19)),
                pseudoknot_type="H-type",
                note="Reduced tmRNA PK1-like H-type motif.",
            ),
            BenchmarkEntry(
                name="HDV_ribozyme",
                sequence=self._build_sequence(28, [(1, 15), (2, 14), (3, 13), (9, 24), (10, 23), (11, 22)]),
                reference_pairs=((1, 15), (2, 14), (3, 13), (9, 24), (10, 23), (11, 22)),
                pseudoknot_type="H-type",
                note="Reduced HDV ribozyme pseudoknot core.",
            ),
            BenchmarkEntry(
                name="BWYV_kissing",
                sequence=self._build_sequence(30, [(0, 11), (1, 10), (5, 18), (6, 17), (12, 27), (13, 26)]),
                reference_pairs=((0, 11), (1, 10), (5, 18), (6, 17), (12, 27), (13, 26)),
                pseudoknot_type="kissing hairpin",
                note="Beet western yellows virus kissing-loop inspired motif.",
            ),
            BenchmarkEntry(
                name="Human_telomerase_core",
                sequence=self._build_sequence(38, [(0, 13), (1, 12), (5, 18), (6, 17), (21, 32), (22, 31), (26, 37), (27, 36)]),
                reference_pairs=((0, 13), (1, 12), (5, 18), (6, 17), (21, 32), (22, 31), (26, 37), (27, 36)),
                pseudoknot_type="recursive",
                note="Two-tier telomerase-like recursive pseudoknot arrangement.",
            ),
            BenchmarkEntry(
                name="MMTV_frameshift",
                sequence=self._build_sequence(27, [(0, 15), (1, 14), (2, 13), (7, 23), (8, 22), (9, 21)]),
                reference_pairs=((0, 15), (1, 14), (2, 13), (7, 23), (8, 22), (9, 21)),
                pseudoknot_type="H-type",
                note="Frameshifting H-type pseudoknot inspired by MMTV/SARS class motifs.",
            ),
        ]

    def _build_sequence(self, length: int, pairs: Sequence[Pair]) -> str:
        sequence = ["" for _ in range(length)]
        motif_cycle = [("G", "C"), ("C", "G"), ("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")]
        for index, (left, right) in enumerate(_normalize_pairs(pairs)):
            nt_left, nt_right = motif_cycle[index % len(motif_cycle)]
            sequence[left] = sequence[left] or nt_left
            sequence[right] = sequence[right] or nt_right
        fillers = "ACGU"
        for index in range(length):
            if not sequence[index]:
                sequence[index] = fillers[index % len(fillers)]
        return "".join(sequence)


def _enumerate_stem_candidates(
    sequence: str,
    scorer: PseudoknotScorer,
    min_stem: int,
    max_stem: Optional[int] = None,
    allowed_positions: Optional[Iterable[int]] = None,
) -> List[StemCandidate]:
    allowed = set(allowed_positions) if allowed_positions is not None else None
    sequence = sequence.upper().replace("T", "U")
    candidates: List[StemCandidate] = []
    n = len(sequence)
    for i in range(n):
        for j in range(i + 4, n):
            if not _allowed(allowed, i, j):
                continue
            run_length = 0
            while i + run_length < j - run_length and _can_pair(sequence[i + run_length], sequence[j - run_length]):
                if not _allowed(allowed, i + run_length, j - run_length):
                    break
                run_length += 1
                if max_stem is not None and run_length >= max_stem:
                    break
            if run_length < min_stem:
                continue
            for length in range(min_stem, run_length + 1):
                pairs = tuple((i + offset, j - offset) for offset in range(length))
                candidates.append(StemCandidate(pairs=pairs, energy=scorer.score_stem(pairs, sequence)))
    candidates.sort(key=lambda candidate: (candidate.left_start, candidate.right_end, candidate.energy, -candidate.length))
    return candidates


__all__ = [
    "AkutsuDP",
    "BenchmarkEntry",
    "BenchmarkPseudoknot",
    "HeuristicPseudoknot",
    "IterativeRelaxation",
    "Pair",
    "PseudoknotDetector",
    "PseudoknotScorer",
    "StemCandidate",
    "TurnerParameters",
]
