"""RNA secondary structure prediction with Turner-style nearest-neighbor energetics.

This module implements two dynamic-programming predictors:

* :class:`NussinovDP` – a baseline O(n^3) maximum-base-pair model.
* :class:`ZukerMFE` – a Turner-inspired minimum-free-energy (MFE) model with
  Zuker-style dynamic programming matrices.

The implementation uses simplified but realistic thermodynamic parameters at
37°C based on the Turner 2004 nearest-neighbor model. Energies are reported in
kcal/mol throughout.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import exp, sqrt
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:  # Optional dependency used only for parameter fitting.
    from scipy.optimize import differential_evolution
except Exception:  # pragma: no cover - SciPy may be unavailable at runtime.
    differential_evolution = None

BASE_PAIRS: set[tuple[str, str]] = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}

INF = float("inf")
R_GAS_CONSTANT = 0.0019872041  # kcal mol^-1 K^-1
DEFAULT_TEMPERATURE_K = 310.15


Pair = Tuple[int, int]
BasePair = Tuple[str, str]
StackKey = Tuple[BasePair, BasePair]


def can_pair(base_i: str, base_j: str) -> bool:
    """Return ``True`` if two nucleotides form a canonical RNA pair."""

    return (base_i.upper(), base_j.upper()) in BASE_PAIRS


def dot_bracket_to_pairs(structure: str) -> List[Pair]:
    """Convert dot-bracket notation into a sorted list of base pairs."""

    stack: list[int] = []
    pairs: list[Pair] = []
    for idx, char in enumerate(structure):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket string: unexpected ')'.")
            start = stack.pop()
            pairs.append((start, idx))
        elif char != ".":
            raise ValueError(f"Unsupported dot-bracket character: {char!r}")
    if stack:
        raise ValueError("Unbalanced dot-bracket string: missing ')'.")
    return sorted(pairs)


def pairs_to_dot_bracket(pairs: Iterable[Pair], length: int) -> str:
    """Convert a base-pair list into dot-bracket notation."""

    chars = ["."] * length
    seen: set[int] = set()
    for i, j in pairs:
        if not (0 <= i < j < length):
            raise ValueError(f"Invalid base pair {(i, j)} for length {length}.")
        if i in seen or j in seen:
            raise ValueError("Each position may appear in at most one base pair.")
        seen.add(i)
        seen.add(j)
        chars[i] = "("
        chars[j] = ")"
    return "".join(chars)


def calculate_f1(predicted_pairs: Iterable[Pair], true_pairs: Iterable[Pair]) -> float:
    """Calculate the base-pair F1 score between two structures."""

    pred = set(predicted_pairs)
    true = set(true_pairs)
    if not pred and not true:
        return 1.0
    tp = len(pred & true)
    fp = len(pred - true)
    fn = len(true - pred)
    denom = 2 * tp + fp + fn
    return 0.0 if denom == 0 else (2 * tp) / denom


def calculate_mcc(
    predicted_pairs: Iterable[Pair], true_pairs: Iterable[Pair], length: int
) -> float:
    """Calculate the Matthews correlation coefficient for base-pair prediction."""

    pred = set(predicted_pairs)
    true = set(true_pairs)
    total_possible = length * (length - 1) // 2
    tp = len(pred & true)
    fp = len(pred - true)
    fn = len(true - pred)
    tn = total_possible - tp - fp - fn
    numerator = tp * tn - fp * fn
    denominator = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return 0.0
    return numerator / denominator


def calculate_sensitivity_ppv(
    predicted_pairs: Iterable[Pair], true_pairs: Iterable[Pair]
) -> Tuple[float, float]:
    """Return sensitivity and positive predictive value (PPV)."""

    pred = set(predicted_pairs)
    true = set(true_pairs)
    tp = len(pred & true)
    fn = len(true - pred)
    fp = len(pred - true)
    sensitivity = 0.0 if (tp + fn) == 0 else tp / (tp + fn)
    ppv = 0.0 if (tp + fp) == 0 else tp / (tp + fp)
    return sensitivity, ppv


def _ensure_rna(sequence: str) -> str:
    seq = sequence.upper().replace("T", "U")
    if any(base not in {"A", "C", "G", "U"} for base in seq):
        raise ValueError("RNA sequence may only contain A, C, G, U (or T).")
    return seq


def _reverse_stack(key: StackKey) -> StackKey:
    (i, j), (k, l) = key
    return ((l, k), (j, i))


def _build_stacking_energies() -> Dict[StackKey, float]:
    """Build a symmetric Turner-style stacking table in kcal/mol."""

    seeds: Dict[StackKey, float] = {
        (("A", "U"), ("A", "U")): -0.93,
        (("A", "U"), ("U", "A")): -1.10,
        (("A", "U"), ("G", "C")): -2.08,
        (("A", "U"), ("C", "G")): -2.24,
        (("A", "U"), ("G", "U")): -0.55,
        (("A", "U"), ("U", "G")): -1.36,
        (("U", "A"), ("A", "U")): -1.33,
        (("U", "A"), ("U", "A")): -0.93,
        (("U", "A"), ("G", "C")): -2.35,
        (("U", "A"), ("C", "G")): -2.11,
        (("U", "A"), ("G", "U")): -1.27,
        (("U", "A"), ("U", "G")): -0.50,
        (("G", "C"), ("A", "U")): -2.35,
        (("G", "C"), ("U", "A")): -3.26,
        (("G", "C"), ("G", "C")): -3.42,
        (("G", "C"), ("C", "G")): -2.36,
        (("G", "C"), ("G", "U")): -2.51,
        (("G", "C"), ("U", "G")): -1.41,
        (("C", "G"), ("A", "U")): -2.11,
        (("C", "G"), ("U", "A")): -2.08,
        (("C", "G"), ("G", "C")): -2.36,
        (("C", "G"), ("C", "G")): -3.26,
        (("C", "G"), ("G", "U")): -1.53,
        (("C", "G"), ("U", "G")): -2.51,
        (("G", "U"), ("A", "U")): -1.27,
        (("G", "U"), ("U", "A")): -1.36,
        (("G", "U"), ("G", "C")): -2.51,
        (("G", "U"), ("C", "G")): -1.53,
        (("G", "U"), ("G", "U")): -0.50,
        (("G", "U"), ("U", "G")): 0.30,
        (("U", "G"), ("A", "U")): -0.50,
        (("U", "G"), ("U", "A")): -0.55,
        (("U", "G"), ("G", "C")): -1.41,
        (("U", "G"), ("C", "G")): -2.51,
        (("U", "G"), ("G", "U")): 0.30,
        (("U", "G"), ("U", "G")): -0.50,
    }
    energies = dict(seeds)
    for key, value in list(seeds.items()):
        energies.setdefault(_reverse_stack(key), value)
    return energies


def _default_hairpin_initiation() -> Dict[int, float]:
    return {
        1: 100.0,
        2: 100.0,
        3: 5.40,
        4: 5.60,
        5: 5.70,
        6: 5.40,
        7: 5.90,
        8: 6.40,
        9: 6.50,
        10: 6.60,
        11: 6.70,
        12: 6.80,
        13: 6.90,
        14: 7.00,
        15: 7.10,
        16: 7.20,
        17: 7.30,
        18: 7.40,
        19: 7.50,
        20: 7.60,
        21: 7.70,
        22: 7.80,
        23: 7.90,
        24: 8.00,
        25: 8.10,
        26: 8.20,
        27: 8.30,
        28: 8.40,
        29: 8.50,
        30: 8.60,
    }


def _default_bulge_initiation() -> Dict[int, float]:
    return {
        1: 3.80,
        2: 2.80,
        3: 3.20,
        4: 3.60,
        5: 4.00,
        6: 4.40,
        7: 4.60,
        8: 4.80,
        9: 5.00,
        10: 5.20,
        11: 5.40,
        12: 5.60,
        13: 5.80,
        14: 6.00,
        15: 6.20,
        16: 6.40,
        17: 6.60,
        18: 6.80,
        19: 7.00,
        20: 7.20,
        21: 7.30,
        22: 7.40,
        23: 7.50,
        24: 7.60,
        25: 7.70,
        26: 7.80,
        27: 7.90,
        28: 8.00,
        29: 8.10,
        30: 8.20,
    }


def _default_internal_initiation() -> Dict[int, float]:
    return {
        1: 0.0,
        2: 0.8,
        3: 1.3,
        4: 1.7,
        5: 2.1,
        6: 2.5,
        7: 2.8,
        8: 3.0,
        9: 3.2,
        10: 3.4,
        11: 3.6,
        12: 3.8,
        13: 4.0,
        14: 4.2,
        15: 4.4,
        16: 4.6,
        17: 4.8,
        18: 5.0,
        19: 5.1,
        20: 5.2,
        21: 5.3,
        22: 5.4,
        23: 5.5,
        24: 5.6,
        25: 5.7,
        26: 5.8,
        27: 5.9,
        28: 6.0,
        29: 6.1,
        30: 6.2,
    }


def _default_symmetric_internal_bonus() -> Dict[int, float]:
    return {1: 0.0, 2: -0.4, 3: -0.3, 4: -0.2, 5: -0.1, 6: -0.1}


def _default_asymmetry_penalty() -> Dict[int, float]:
    return {0: 0.0, 1: 0.4, 2: 0.7, 3: 1.0, 4: 1.2, 5: 1.4, 6: 1.6}


def _default_tetraloop_bonuses() -> Dict[str, float]:
    return {
        "GGAA": -3.0,
        "GGAC": -3.0,
        "GGAG": -3.0,
        "GAAA": -2.5,
        "GUGA": -2.5,
        "CUUG": -2.2,
        "UUCG": -2.5,
        "UGAA": -2.0,
        "CGAA": -2.0,
        "AGAA": -2.0,
        "UACG": -2.0,
        "UUUG": -1.8,
    }


def _default_dangling_5p() -> Dict[Tuple[BasePair, str], float]:
    entries = {
        (("A", "U"), "A"): -0.3,
        (("A", "U"), "C"): -0.2,
        (("A", "U"), "G"): -0.4,
        (("A", "U"), "U"): -0.2,
        (("U", "A"), "A"): -0.2,
        (("U", "A"), "C"): -0.1,
        (("U", "A"), "G"): -0.3,
        (("U", "A"), "U"): -0.2,
        (("G", "C"), "A"): -0.8,
        (("G", "C"), "C"): -0.5,
        (("G", "C"), "G"): -0.8,
        (("G", "C"), "U"): -0.6,
        (("C", "G"), "A"): -1.1,
        (("C", "G"), "C"): -0.4,
        (("C", "G"), "G"): -0.7,
        (("C", "G"), "U"): -0.7,
        (("G", "U"), "A"): -0.4,
        (("G", "U"), "C"): -0.2,
        (("G", "U"), "G"): -0.5,
        (("G", "U"), "U"): -0.2,
        (("U", "G"), "A"): -0.2,
        (("U", "G"), "C"): -0.2,
        (("U", "G"), "G"): -0.4,
        (("U", "G"), "U"): -0.2,
    }
    return entries


def _default_dangling_3p() -> Dict[Tuple[BasePair, str], float]:
    entries = {
        (("A", "U"), "A"): -0.5,
        (("A", "U"), "C"): -0.3,
        (("A", "U"), "G"): -0.6,
        (("A", "U"), "U"): -0.4,
        (("U", "A"), "A"): -0.4,
        (("U", "A"), "C"): -0.2,
        (("U", "A"), "G"): -0.5,
        (("U", "A"), "U"): -0.3,
        (("G", "C"), "A"): -0.8,
        (("G", "C"), "C"): -0.6,
        (("G", "C"), "G"): -1.0,
        (("G", "C"), "U"): -0.7,
        (("C", "G"), "A"): -0.9,
        (("C", "G"), "C"): -0.6,
        (("C", "G"), "G"): -0.9,
        (("C", "G"), "U"): -0.8,
        (("G", "U"), "A"): -0.5,
        (("G", "U"), "C"): -0.3,
        (("G", "U"), "G"): -0.6,
        (("G", "U"), "U"): -0.4,
        (("U", "G"), "A"): -0.4,
        (("U", "G"), "C"): -0.2,
        (("U", "G"), "G"): -0.5,
        (("U", "G"), "U"): -0.3,
    }
    return entries


def _default_terminal_mismatches() -> Dict[Tuple[BasePair, str, str], float]:
    mismatches: Dict[Tuple[BasePair, str, str], float] = {}
    pair_bias = {
        ("A", "U"): 0.4,
        ("U", "A"): 0.4,
        ("G", "C"): 0.1,
        ("C", "G"): 0.1,
        ("G", "U"): 0.5,
        ("U", "G"): 0.5,
    }
    flank_bonus = {"A": 0.10, "C": 0.00, "G": -0.10, "U": 0.05}
    for pair, bias in pair_bias.items():
        for left in "ACGU":
            for right in "ACGU":
                mismatches[(pair, left, right)] = bias + flank_bonus[left] + flank_bonus[right]
    return mismatches


@dataclass
class TurnerParameters:
    """Container for Turner-style nearest-neighbor thermodynamic parameters."""

    stacking_energies: Dict[StackKey, float] = field(default_factory=_build_stacking_energies)
    hairpin_initiation: Dict[int, float] = field(default_factory=_default_hairpin_initiation)
    hairpin_tetraloop_bonuses: Dict[str, float] = field(default_factory=_default_tetraloop_bonuses)
    internal_loop_initiation: Dict[int, float] = field(default_factory=_default_internal_initiation)
    symmetric_internal_bonus: Dict[int, float] = field(default_factory=_default_symmetric_internal_bonus)
    asymmetry_penalty: Dict[int, float] = field(default_factory=_default_asymmetry_penalty)
    bulge_initiation: Dict[int, float] = field(default_factory=_default_bulge_initiation)
    multiloop_a: float = 3.4
    multiloop_b: float = 0.4
    multiloop_c: float = 0.1
    dangling_5p: Dict[Tuple[BasePair, str], float] = field(default_factory=_default_dangling_5p)
    dangling_3p: Dict[Tuple[BasePair, str], float] = field(default_factory=_default_dangling_3p)
    terminal_mismatch: Dict[Tuple[BasePair, str, str], float] = field(
        default_factory=_default_terminal_mismatches
    )
    temperature_k: float = DEFAULT_TEMPERATURE_K

    def get_stacking_energy(self, outer_pair: BasePair, inner_pair: BasePair) -> float:
        """Return the stacking free energy for adjacent base pairs."""

        key = ((outer_pair[0].upper(), outer_pair[1].upper()), (inner_pair[0].upper(), inner_pair[1].upper()))
        return self.stacking_energies.get(key, 0.0)

    def get_loop_initiation_energy(self, loop_type: str, size: int) -> float:
        """Return size-dependent initiation free energy for a loop."""

        if size < 0:
            raise ValueError("Loop size must be non-negative.")
        lookup = {
            "hairpin": self.hairpin_initiation,
            "bulge": self.bulge_initiation,
            "internal": self.internal_loop_initiation,
        }.get(loop_type.lower())
        if lookup is None:
            raise ValueError(f"Unsupported loop type: {loop_type!r}")
        if size in lookup:
            return lookup[size]
        max_size = max(lookup)
        if size <= max_size:
            return lookup[max_size]
        return lookup[max_size] + 1.75 * np.log(size / max_size)

    def get_dangling_energy(self, pair: BasePair, base: str, five_prime: bool = True) -> float:
        """Return the dangling-end stabilization for a nucleotide adjacent to a helix."""

        lookup = self.dangling_5p if five_prime else self.dangling_3p
        return lookup.get(((pair[0].upper(), pair[1].upper()), base.upper()), 0.0)

    def get_terminal_mismatch_energy(self, pair: BasePair, left: str, right: str) -> float:
        """Return the terminal mismatch penalty for nucleotides closing a loop."""

        return self.terminal_mismatch.get(
            ((pair[0].upper(), pair[1].upper()), left.upper(), right.upper()), 0.0
        )

    def hairpin_energy(self, sequence: str, i: int, j: int) -> float:
        """Calculate Turner-style hairpin loop energy closed by ``(i, j)``."""

        size = j - i - 1
        pair = (sequence[i], sequence[j])
        energy = self.get_loop_initiation_energy("hairpin", size)
        if size == 3:
            if sequence[i] in {"G", "C"} and sequence[j] in {"G", "C"}:
                energy += -0.4
        elif size >= 4:
            loop = sequence[i + 1 : j]
            if size == 4:
                energy += self.hairpin_tetraloop_bonuses.get(loop, 0.0)
            energy += self.get_terminal_mismatch_energy(pair, sequence[i + 1], sequence[j - 1])
        return energy

    def bulge_energy(self, size: int, closing_pair: Optional[BasePair] = None) -> float:
        """Calculate bulge loop energy."""

        energy = self.get_loop_initiation_energy("bulge", size)
        if size == 1 and closing_pair is not None:
            if closing_pair in {("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")}:
                energy += 0.5
        return energy

    def internal_loop_energy(
        self,
        sequence: str,
        i: int,
        j: int,
        p: int,
        q: int,
    ) -> float:
        """Calculate internal-loop energy for outer pair ``(i, j)`` and inner pair ``(p, q)``."""

        left = p - i - 1
        right = j - q - 1
        if left < 0 or right < 0 or (left == 0 and right == 0):
            raise ValueError("Invalid internal loop configuration.")
        if left == 0 or right == 0:
            size = left + right
            return self.bulge_energy(size, (sequence[i], sequence[j]))
        total = left + right
        energy = self.get_loop_initiation_energy("internal", total)
        if left == right:
            energy += self.symmetric_internal_bonus.get(left, 0.0)
        asym = abs(left - right)
        capped = min(asym, max(self.asymmetry_penalty))
        energy += self.asymmetry_penalty.get(capped, 0.0)
        outer_pair = (sequence[i], sequence[j])
        inner_pair = (sequence[p], sequence[q])
        energy += self.get_terminal_mismatch_energy(outer_pair, sequence[i + 1], sequence[j - 1])
        energy += self.get_terminal_mismatch_energy(inner_pair, sequence[p - 1], sequence[q + 1])
        if total == 2 and left == 1 and right == 1:
            energy -= 0.2
        return energy

    def scaled(
        self,
        stacking_scale: float = 1.0,
        hairpin_scale: float = 1.0,
        internal_scale: float = 1.0,
        bulge_scale: float = 1.0,
        multiloop_a: Optional[float] = None,
        multiloop_b: Optional[float] = None,
        multiloop_c: Optional[float] = None,
    ) -> "TurnerParameters":
        """Return a copy with globally scaled parameter subsets."""

        return TurnerParameters(
            stacking_energies={k: v * stacking_scale for k, v in self.stacking_energies.items()},
            hairpin_initiation={k: v * hairpin_scale for k, v in self.hairpin_initiation.items()},
            hairpin_tetraloop_bonuses=dict(self.hairpin_tetraloop_bonuses),
            internal_loop_initiation={k: v * internal_scale for k, v in self.internal_loop_initiation.items()},
            symmetric_internal_bonus=dict(self.symmetric_internal_bonus),
            asymmetry_penalty=dict(self.asymmetry_penalty),
            bulge_initiation={k: v * bulge_scale for k, v in self.bulge_initiation.items()},
            multiloop_a=self.multiloop_a if multiloop_a is None else multiloop_a,
            multiloop_b=self.multiloop_b if multiloop_b is None else multiloop_b,
            multiloop_c=self.multiloop_c if multiloop_c is None else multiloop_c,
            dangling_5p=dict(self.dangling_5p),
            dangling_3p=dict(self.dangling_3p),
            terminal_mismatch=dict(self.terminal_mismatch),
            temperature_k=self.temperature_k,
        )

    def optimize_parameters(
        self,
        training_set: Sequence[Tuple[str, str]],
        maxiter: int = 25,
        popsize: int = 10,
        min_loop_length: int = 3,
    ):
        """Optimize global Turner parameter scales with differential evolution.

        The optimization adjusts a compact set of scale factors and multi-loop
        coefficients to maximize mean base-pair F1 score across a training set.
        The function returns ``(optimized_parameters, scipy_result)``.
        """

        if differential_evolution is None:
            raise ImportError("scipy.optimize.differential_evolution is required for optimization.")
        dataset = [(_ensure_rna(seq), struct) for seq, struct in training_set]

        def objective(vector: np.ndarray) -> float:
            candidate = self.scaled(
                stacking_scale=float(vector[0]),
                hairpin_scale=float(vector[1]),
                internal_scale=float(vector[2]),
                bulge_scale=float(vector[3]),
                multiloop_a=float(vector[4]),
                multiloop_b=float(vector[5]),
                multiloop_c=float(vector[6]),
            )
            scores: list[float] = []
            for seq, known_structure in dataset:
                predictor = ZukerMFE(seq, params=candidate, min_loop_length=min_loop_length)
                predicted = predictor.predict()
                scores.append(
                    calculate_f1(dot_bracket_to_pairs(predicted), dot_bracket_to_pairs(known_structure))
                )
            return 1.0 - float(np.mean(scores))

        bounds = [
            (0.6, 1.4),
            (0.6, 1.4),
            (0.6, 1.5),
            (0.6, 1.5),
            (1.0, 6.0),
            (0.0, 2.0),
            (0.0, 1.0),
        ]
        result = differential_evolution(objective, bounds=bounds, maxiter=maxiter, popsize=popsize)
        optimized = self.scaled(
            stacking_scale=float(result.x[0]),
            hairpin_scale=float(result.x[1]),
            internal_scale=float(result.x[2]),
            bulge_scale=float(result.x[3]),
            multiloop_a=float(result.x[4]),
            multiloop_b=float(result.x[5]),
            multiloop_c=float(result.x[6]),
        )
        return optimized, result


class NussinovDP:
    """Baseline O(n^3) Nussinov dynamic program maximizing base-pair count."""

    def __init__(self, sequence: str, min_loop_length: int = 3) -> None:
        self.sequence = _ensure_rna(sequence)
        self.min_loop_length = min_loop_length
        self.n = len(self.sequence)
        self.dp = np.zeros((self.n, self.n), dtype=np.int32)

    def fill(self) -> np.ndarray:
        """Fill the Nussinov DP table."""

        for span in range(1, self.n):
            for i in range(self.n - span):
                j = i + span
                best = self.dp[i + 1, j] if i + 1 <= j else 0
                if i <= j - 1:
                    best = max(best, self.dp[i, j - 1])
                if j - i > self.min_loop_length and can_pair(self.sequence[i], self.sequence[j]):
                    best = max(best, (self.dp[i + 1, j - 1] if i + 1 <= j - 1 else 0) + 1)
                for k in range(i, j):
                    best = max(best, self.dp[i, k] + self.dp[k + 1, j])
                self.dp[i, j] = best
        return self.dp

    def traceback(self) -> List[Pair]:
        """Recover an optimal maximum-pairing structure from the DP table."""

        if self.n == 0:
            return []
        if not np.any(self.dp):
            self.fill()
        pairs: list[Pair] = []

        def _trace(i: int, j: int) -> None:
            if i >= j:
                return
            current = self.dp[i, j]
            if current == (self.dp[i + 1, j] if i + 1 <= j else 0):
                _trace(i + 1, j)
                return
            if current == (self.dp[i, j - 1] if i <= j - 1 else 0):
                _trace(i, j - 1)
                return
            if (
                j - i > self.min_loop_length
                and can_pair(self.sequence[i], self.sequence[j])
                and current == (self.dp[i + 1, j - 1] if i + 1 <= j - 1 else 0) + 1
            ):
                pairs.append((i, j))
                _trace(i + 1, j - 1)
                return
            for k in range(i, j):
                if current == self.dp[i, k] + self.dp[k + 1, j]:
                    _trace(i, k)
                    _trace(k + 1, j)
                    return

        _trace(0, self.n - 1)
        return sorted(pairs)

    def predict(self) -> str:
        """Return the Nussinov structure in dot-bracket notation."""

        self.fill()
        return pairs_to_dot_bracket(self.traceback(), self.n)


class ZukerMFE:
    """Turner-style Zuker dynamic program for RNA minimum free energy prediction."""

    def __init__(
        self,
        sequence: str,
        params: Optional[TurnerParameters] = None,
        min_loop_length: int = 3,
        max_internal_loop_size: int = 30,
    ) -> None:
        self.sequence = _ensure_rna(sequence)
        self.params = params or TurnerParameters()
        self.min_loop_length = min_loop_length
        self.max_internal_loop_size = max_internal_loop_size
        self.n = len(self.sequence)
        self.W = np.full((self.n, self.n), INF, dtype=float)
        self.V = np.full((self.n, self.n), INF, dtype=float)
        self.WM = np.full((self.n, self.n), INF, dtype=float)
        self.WM2 = np.full((self.n, self.n), INF, dtype=float)
        self._choice_W: Dict[Tuple[int, int], Tuple[str, Optional[int]]] = {}
        self._choice_V: Dict[Tuple[int, int], Tuple[str, Optional[Tuple[int, int]]]] = {}
        self._choice_WM: Dict[Tuple[int, int], Tuple[str, Optional[int]]] = {}
        self._choice_WM2: Dict[Tuple[int, int], Optional[int]] = {}
        self._filled = False

    def _w(self, i: int, j: int) -> float:
        return 0.0 if i > j else self.W[i, j]

    def _v(self, i: int, j: int) -> float:
        return INF if i >= self.n or j < 0 or i >= j else self.V[i, j]

    def _wm(self, i: int, j: int) -> float:
        return INF if i > j else self.WM[i, j]

    def _wm2(self, i: int, j: int) -> float:
        return INF if i > j else self.WM2[i, j]

    def fill_matrices(self) -> None:
        """Fill Zuker DP matrices ``W``, ``V``, ``WM``, and ``WM2``."""

        if self.n == 0:
            self._filled = True
            return
        for i in range(self.n):
            self.W[i, i] = 0.0
            self.WM[i, i] = self.params.multiloop_c
        for span in range(1, self.n):
            for i in range(self.n - span):
                j = i + span
                self._fill_v(i, j)
                self._fill_wm(i, j)
                self._fill_wm2(i, j)
                self._fill_w(i, j)
        self._filled = True

    def _fill_v(self, i: int, j: int) -> None:
        if j - i <= self.min_loop_length or not can_pair(self.sequence[i], self.sequence[j]):
            self.V[i, j] = INF
            return

        pair = (self.sequence[i], self.sequence[j])
        best = self.params.hairpin_energy(self.sequence, i, j)
        choice: Tuple[str, Optional[Tuple[int, int]]] = ("hairpin", None)

        if i + 1 < j - 1 and can_pair(self.sequence[i + 1], self.sequence[j - 1]):
            stacked = self.params.get_stacking_energy(pair, (self.sequence[i + 1], self.sequence[j - 1]))
            candidate = self.V[i + 1, j - 1] + stacked
            if candidate < best:
                best = candidate
                choice = ("stack", (i + 1, j - 1))

        p_max = min(j - self.min_loop_length - 1, i + self.max_internal_loop_size + 1)
        for p in range(i + 1, p_max + 1):
            q_min = max(p + self.min_loop_length + 1, j - self.max_internal_loop_size - 1)
            for q in range(q_min, j):
                if p >= q or not can_pair(self.sequence[p], self.sequence[q]):
                    continue
                left = p - i - 1
                right = j - q - 1
                total = left + right
                if total == 0 or total > self.max_internal_loop_size:
                    continue
                inner = self.V[p, q]
                if np.isinf(inner):
                    continue
                loop_energy = self.params.internal_loop_energy(self.sequence, i, j, p, q)
                candidate = inner + loop_energy
                if candidate < best:
                    best = candidate
                    choice = ("internal", (p, q))

        if i + 1 <= j - 1 and np.isfinite(self.WM2[i + 1, j - 1]):
            candidate = (
                self.params.multiloop_a
                + self.WM2[i + 1, j - 1]
                + self.params.get_terminal_mismatch_energy(pair, self.sequence[i + 1], self.sequence[j - 1])
            )
            if candidate < best:
                best = candidate
                choice = ("multiloop", (i + 1, j - 1))

        self.V[i, j] = best
        self._choice_V[(i, j)] = choice

    def _fill_wm(self, i: int, j: int) -> None:
        best = INF
        choice: Tuple[str, Optional[int]] = ("none", None)

        if i + 1 <= j and np.isfinite(self.WM[i + 1, j]):
            candidate = self.params.multiloop_c + self.WM[i + 1, j]
            if candidate < best:
                best = candidate
                choice = ("skip_i", None)

        if i <= j - 1 and np.isfinite(self.WM[i, j - 1]):
            candidate = self.params.multiloop_c + self.WM[i, j - 1]
            if candidate < best:
                best = candidate
                choice = ("skip_j", None)

        if np.isfinite(self.V[i, j]):
            candidate = self.params.multiloop_b + self.V[i, j]
            if candidate < best:
                best = candidate
                choice = ("branch", None)

        for k in range(i, j):
            left = self.WM[i, k]
            right = self.WM[k + 1, j]
            if np.isfinite(left) and np.isfinite(right):
                candidate = left + right
                if candidate < best:
                    best = candidate
                    choice = ("split", k)

        self.WM[i, j] = best
        self._choice_WM[(i, j)] = choice

    def _fill_wm2(self, i: int, j: int) -> None:
        best = INF
        best_k: Optional[int] = None
        for k in range(i, j):
            left = self.WM[i, k]
            right = self.WM[k + 1, j]
            if np.isfinite(left) and np.isfinite(right):
                candidate = left + right
                if candidate < best:
                    best = candidate
                    best_k = k
        self.WM2[i, j] = best
        self._choice_WM2[(i, j)] = best_k

    def _fill_w(self, i: int, j: int) -> None:
        best = self.W[i + 1, j] if i + 1 <= j else 0.0
        choice: Tuple[str, Optional[int]] = ("skip_i", None)

        if i <= j - 1 and self.W[i, j - 1] < best:
            best = self.W[i, j - 1]
            choice = ("skip_j", None)

        if np.isfinite(self.V[i, j]) and self.V[i, j] < best:
            best = self.V[i, j]
            choice = ("pair", None)

        for k in range(i, j):
            candidate = self.W[i, k] + self.W[k + 1, j]
            if candidate < best:
                best = candidate
                choice = ("split", k)

        self.W[i, j] = best
        self._choice_W[(i, j)] = choice

    def mfe(self) -> float:
        """Return the minimum free energy of the full sequence."""

        if not self._filled:
            self.fill_matrices()
        if self.n == 0:
            return 0.0
        return float(self.W[0, self.n - 1])

    def traceback(self) -> List[Pair]:
        """Recover the optimal MFE base pairs."""

        if not self._filled:
            self.fill_matrices()
        pairs: list[Pair] = []

        def trace_w(i: int, j: int) -> None:
            if i >= j or i < 0 or j >= self.n:
                return
            action, payload = self._choice_W.get((i, j), ("skip_i", None))
            if action == "skip_i":
                trace_w(i + 1, j)
            elif action == "skip_j":
                trace_w(i, j - 1)
            elif action == "pair":
                trace_v(i, j)
            elif action == "split" and payload is not None:
                trace_w(i, payload)
                trace_w(payload + 1, j)

        def trace_v(i: int, j: int) -> None:
            pairs.append((i, j))
            action, payload = self._choice_V.get((i, j), ("hairpin", None))
            if action in {"hairpin", "none"}:
                return
            if action in {"stack", "internal"} and payload is not None:
                p, q = payload
                trace_v(p, q)
            elif action == "multiloop" and payload is not None:
                p, q = payload
                trace_wm2(p, q)

        def trace_wm(i: int, j: int) -> None:
            if i > j or i < 0 or j >= self.n:
                return
            action, payload = self._choice_WM.get((i, j), ("none", None))
            if action == "skip_i":
                trace_wm(i + 1, j)
            elif action == "skip_j":
                trace_wm(i, j - 1)
            elif action == "branch":
                trace_v(i, j)
            elif action == "split" and payload is not None:
                trace_wm(i, payload)
                trace_wm(payload + 1, j)

        def trace_wm2(i: int, j: int) -> None:
            if i >= j or i < 0 or j >= self.n:
                return
            k = self._choice_WM2.get((i, j))
            if k is None:
                return
            trace_wm(i, k)
            trace_wm(k + 1, j)

        if self.n:
            trace_w(0, self.n - 1)
        return sorted(set(pairs))

    def predict(self) -> str:
        """Return the MFE structure in dot-bracket notation."""

        return pairs_to_dot_bracket(self.traceback(), self.n)

    def _boltzmann(self, energy: float) -> float:
        beta = 1.0 / (R_GAS_CONSTANT * self.params.temperature_k)
        return float(np.exp(-beta * energy))

    def compute_partition_function(self) -> Tuple[float, np.ndarray]:
        """Compute a simplified McCaskill-style partition function and pair probabilities.

        Returns
        -------
        tuple
            ``(Z, P)`` where ``Z`` is the partition function and ``P[i, j]`` is
            the estimated probability that ``i`` and ``j`` form a base pair.
        """

        n = self.n
        if n == 0:
            return 1.0, np.zeros((0, 0), dtype=float)

        QB = np.zeros((n, n), dtype=float)
        QM = np.zeros((n, n), dtype=float)
        QM2 = np.zeros((n, n), dtype=float)
        Q = np.zeros((n, n), dtype=float)

        def q(i: int, j: int) -> float:
            return 1.0 if i > j else Q[i, j]

        def qm(i: int, j: int) -> float:
            return 0.0 if i > j else QM[i, j]

        for span in range(1, n):
            for i in range(n - span):
                j = i + span
                if j - i > self.min_loop_length and can_pair(self.sequence[i], self.sequence[j]):
                    hairpin = self._boltzmann(self.params.hairpin_energy(self.sequence, i, j))
                    qb = hairpin
                    if i + 1 < j - 1 and can_pair(self.sequence[i + 1], self.sequence[j - 1]):
                        stack_energy = self.params.get_stacking_energy(
                            (self.sequence[i], self.sequence[j]),
                            (self.sequence[i + 1], self.sequence[j - 1]),
                        )
                        qb += self._boltzmann(stack_energy) * QB[i + 1, j - 1]
                    p_max = min(j - self.min_loop_length - 1, i + self.max_internal_loop_size + 1)
                    for p in range(i + 1, p_max + 1):
                        q_min = max(p + self.min_loop_length + 1, j - self.max_internal_loop_size - 1)
                        for q_idx in range(q_min, j):
                            if p >= q_idx or not can_pair(self.sequence[p], self.sequence[q_idx]):
                                continue
                            left = p - i - 1
                            right = j - q_idx - 1
                            total = left + right
                            if total == 0 or total > self.max_internal_loop_size:
                                continue
                            qb += self._boltzmann(
                                self.params.internal_loop_energy(self.sequence, i, j, p, q_idx)
                            ) * QB[p, q_idx]
                    if i + 1 <= j - 1:
                        loop_weight = self._boltzmann(self.params.multiloop_a)
                        qb += loop_weight * QM2[i + 1, j - 1]
                    QB[i, j] = qb

                qm_best = 0.0
                if i + 1 <= j:
                    qm_best += self._boltzmann(self.params.multiloop_c) * QM[i + 1, j]
                if i <= j - 1:
                    qm_best += self._boltzmann(self.params.multiloop_c) * QM[i, j - 1]
                if QB[i, j] > 0.0:
                    qm_best += self._boltzmann(self.params.multiloop_b) * QB[i, j]
                split_sum = 0.0
                for k in range(i, j):
                    split_sum += QM[i, k] * QM[k + 1, j]
                QM[i, j] = qm_best + split_sum

                split2 = 0.0
                for k in range(i, j):
                    split2 += QM[i, k] * QM[k + 1, j]
                QM2[i, j] = split2

                total = q(i + 1, j)
                for k in range(i + self.min_loop_length + 1, j + 1):
                    if QB[i, k] > 0.0:
                        total += QB[i, k] * q(k + 1, j)
                Q[i, j] = total if total > 0.0 else 1.0

        Z = Q[0, n - 1]
        prefix = np.ones(n + 1, dtype=float)
        for end in range(n):
            value = prefix[end]
            for start in range(end - self.min_loop_length):
                value += prefix[start] * QB[start, end]
            prefix[end + 1] = value

        suffix = np.ones(n + 1, dtype=float)
        for start in range(n - 1, -1, -1):
            value = suffix[start + 1]
            for end in range(start + self.min_loop_length + 1, n):
                value += QB[start, end] * suffix[end + 1]
            suffix[start] = value

        probabilities = np.zeros((n, n), dtype=float)
        normalizer = prefix[n] if prefix[n] > 0 else Z
        for i in range(n):
            for j in range(i + self.min_loop_length + 1, n):
                if QB[i, j] > 0.0 and normalizer > 0.0:
                    probabilities[i, j] = min(1.0, (prefix[i] * QB[i, j] * suffix[j + 1]) / normalizer)
                    probabilities[j, i] = probabilities[i, j]
        return float(normalizer), probabilities

    def base_pair_probabilities(self) -> np.ndarray:
        """Return the matrix of base-pair probabilities."""

        _, probabilities = self.compute_partition_function()
        return probabilities


class ParameterOptimizer:
    """Differential-evolution optimizer for Turner parameter calibration."""

    def __init__(
        self,
        training_set: Sequence[Tuple[str, str]],
        base_params: Optional[TurnerParameters] = None,
        min_loop_length: int = 3,
    ) -> None:
        self.training_set = [(_ensure_rna(seq), struct) for seq, struct in training_set]
        self.base_params = base_params or TurnerParameters()
        self.min_loop_length = min_loop_length

    def _score_dataset(self, params: TurnerParameters, dataset: Sequence[Tuple[str, str]]) -> float:
        scores: list[float] = []
        for sequence, known in dataset:
            predictor = ZukerMFE(sequence, params=params, min_loop_length=self.min_loop_length)
            predicted = predictor.predict()
            scores.append(calculate_f1(dot_bracket_to_pairs(predicted), dot_bracket_to_pairs(known)))
        return float(np.mean(scores)) if scores else 0.0

    def optimize(self, maxiter: int = 25, popsize: int = 10):
        """Optimize parameters to maximize mean base-pair F1 score."""

        return self.base_params.optimize_parameters(
            self.training_set,
            maxiter=maxiter,
            popsize=popsize,
            min_loop_length=self.min_loop_length,
        )

    def cross_validate(
        self,
        n_folds: int = 5,
        maxiter: int = 10,
        popsize: int = 8,
        random_seed: int = 0,
    ) -> Dict[str, object]:
        """Perform simple shuffled cross-validation for parameter tuning."""

        if n_folds < 2:
            raise ValueError("n_folds must be at least 2.")
        if len(self.training_set) < n_folds:
            raise ValueError("Training set must contain at least n_folds examples.")

        rng = np.random.default_rng(random_seed)
        indices = np.arange(len(self.training_set))
        rng.shuffle(indices)
        folds = np.array_split(indices, n_folds)
        fold_results: list[Dict[str, float]] = []

        for fold_idx in range(n_folds):
            test_idx = set(int(i) for i in folds[fold_idx])
            train = [self.training_set[i] for i in indices if int(i) not in test_idx]
            test = [self.training_set[i] for i in indices if int(i) in test_idx]
            optimizer = ParameterOptimizer(train, self.base_params, self.min_loop_length)
            params, result = optimizer.optimize(maxiter=maxiter, popsize=popsize)
            train_f1 = optimizer._score_dataset(params, train)
            test_f1 = optimizer._score_dataset(params, test)
            fold_results.append(
                {
                    "fold": float(fold_idx),
                    "train_f1": train_f1,
                    "test_f1": test_f1,
                    "objective": float(result.fun),
                }
            )

        return {
            "folds": fold_results,
            "mean_train_f1": float(np.mean([fold["train_f1"] for fold in fold_results])),
            "mean_test_f1": float(np.mean([fold["test_f1"] for fold in fold_results])),
        }


__all__ = [
    "TurnerParameters",
    "NussinovDP",
    "ZukerMFE",
    "ParameterOptimizer",
    "can_pair",
    "dot_bracket_to_pairs",
    "pairs_to_dot_bracket",
    "calculate_f1",
    "calculate_mcc",
    "calculate_sensitivity_ppv",
]
