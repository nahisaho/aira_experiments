"""Chemical probing utilities for RNA secondary-structure prediction.

This module provides lightweight processors for SHAPE and DMS probing data,
simplified probe-constrained folding routines, iterative constraint enforcement,
data simulation, and agreement metrics.
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
from scipy import stats

PROBE_DTYPE = np.dtype(
    [
        ("position", np.int64),
        ("nucleotide", "U1"),
        ("reactivity", np.float64),
    ]
)

PAIR_ENERGIES = {
    "AU": -2.0,
    "UA": -2.0,
    "GC": -3.0,
    "CG": -3.0,
    "GU": -1.0,
    "UG": -1.0,
}


def _parse_reactivity(value: Any) -> float:
    """Convert a value to a non-negative float or NaN."""
    if value is None:
        return float("nan")
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "" or stripped.lower() in {"na", "nan", "none", "missing"}:
            return float("nan")
        value = stripped
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid reactivity value: {value!r}") from exc
    if parsed < 0:
        return float("nan")
    return parsed


def _finalize_probe_records(records: Iterable[tuple[int, str, float]]) -> np.ndarray:
    """Aggregate records by position and return a sorted structured array."""
    grouped: dict[int, dict[str, Any]] = {}
    for position, nucleotide, reactivity in records:
        if position < 1:
            raise ValueError("Probe positions must be 1-based positive integers.")
        nucleotide = (nucleotide or "N").upper()
        if len(nucleotide) != 1:
            raise ValueError(f"Invalid nucleotide: {nucleotide!r}")
        entry = grouped.setdefault(position, {"nucleotide": nucleotide, "values": []})
        if entry["nucleotide"] == "N" and nucleotide != "N":
            entry["nucleotide"] = nucleotide
        if not np.isnan(reactivity):
            entry["values"].append(float(reactivity))
    data = np.empty(len(grouped), dtype=PROBE_DTYPE)
    for index, position in enumerate(sorted(grouped)):
        values = grouped[position]["values"]
        data[index] = (
            position,
            grouped[position]["nucleotide"],
            float(np.mean(values)) if values else float("nan"),
        )
    return data


def _load_probe_data(filepath_or_data: Any) -> np.ndarray:
    """Load probe data from a path or in-memory representation."""
    if isinstance(filepath_or_data, (str, Path)):
        path = Path(filepath_or_data)
        if not path.exists():
            raise FileNotFoundError(f"Probe data file not found: {path}")
        records: list[tuple[int, str, float]] = []
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.replace(",", " ").split()
                if len(fields) < 3:
                    raise ValueError(
                        f"Expected three columns at line {line_number}, got {len(fields)}."
                    )
                position = int(fields[0])
                nucleotide = fields[1]
                reactivity = _parse_reactivity(fields[2])
                records.append((position, nucleotide, reactivity))
        return _finalize_probe_records(records)

    if isinstance(filepath_or_data, np.ndarray) and filepath_or_data.dtype.names:
        names = {name.lower(): name for name in filepath_or_data.dtype.names}
        required = {"position", "nucleotide", "reactivity"}
        if not required.issubset(names):
            raise ValueError("Structured arrays must contain position, nucleotide, and reactivity fields.")
        records = []
        for row in filepath_or_data:
            records.append(
                (
                    int(row[names["position"]]),
                    str(row[names["nucleotide"]]),
                    _parse_reactivity(row[names["reactivity"]]),
                )
            )
        return _finalize_probe_records(records)

    if isinstance(filepath_or_data, Mapping):
        records = []
        for key, value in filepath_or_data.items():
            if isinstance(value, Mapping):
                nucleotide = str(value.get("nucleotide", "N"))
                reactivity = _parse_reactivity(value.get("reactivity"))
            elif isinstance(value, (tuple, list)) and len(value) >= 2:
                nucleotide = str(value[0])
                reactivity = _parse_reactivity(value[1])
            else:
                nucleotide = "N"
                reactivity = _parse_reactivity(value)
            records.append((int(key), nucleotide, reactivity))
        return _finalize_probe_records(records)

    records = []
    for row in filepath_or_data:
        if isinstance(row, Mapping):
            position = int(row["position"])
            nucleotide = str(row.get("nucleotide", "N"))
            reactivity = _parse_reactivity(row.get("reactivity"))
        else:
            if len(row) < 3:
                raise ValueError("Each probe record must contain position, nucleotide, and reactivity.")
            position = int(row[0])
            nucleotide = str(row[1])
            reactivity = _parse_reactivity(row[2])
        records.append((position, nucleotide, reactivity))
    return _finalize_probe_records(records)


def _winsorize_nonnegative(values: np.ndarray) -> np.ndarray:
    """Clip extreme values using an IQR rule while preserving NaNs."""
    result = np.asarray(values, dtype=float).copy()
    valid_mask = np.isfinite(result)
    valid = result[valid_mask]
    if valid.size < 4:
        return result
    q1, q3 = np.percentile(valid, [25.0, 75.0])
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    lower = max(0.0, q1 - 1.5 * iqr)
    result[valid_mask] = np.clip(valid, lower, upper)
    return result


def _two_to_eight_factor(values: np.ndarray) -> float:
    """Compute the 2-8% normalization factor."""
    valid = np.asarray(values, dtype=float)
    valid = valid[np.isfinite(valid) & (valid > 0)]
    if valid.size == 0:
        raise ValueError("At least one positive reactivity is required for normalization.")
    sorted_values = np.sort(valid)
    n = sorted_values.size
    start = max(int(math.floor(0.92 * n)), 0)
    end = max(int(math.ceil(0.98 * n)), start + 1)
    subset = sorted_values[start:end]
    if subset.size == 0:
        subset = sorted_values[max(n - max(1, n // 10), 0) :]
    factor = float(np.mean(subset))
    return factor if factor > 0 else float(np.max(sorted_values))


def _boxplot_factor(values: np.ndarray) -> float:
    """Compute an alternative normalization factor using a box-plot heuristic."""
    valid = np.asarray(values, dtype=float)
    valid = valid[np.isfinite(valid) & (valid > 0)]
    if valid.size == 0:
        raise ValueError("At least one positive reactivity is required for normalization.")
    q1, q3 = np.percentile(valid, [25.0, 75.0])
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    non_outliers = valid[valid <= upper]
    if non_outliers.size == 0:
        non_outliers = valid
    high_tail = non_outliers[non_outliers >= np.percentile(non_outliers, 90.0)]
    if high_tail.size == 0:
        high_tail = non_outliers
    factor = float(np.mean(high_tail))
    return factor if factor > 0 else float(np.max(non_outliers))


def _dot_bracket_pairs(structure: str) -> dict[int, int]:
    """Parse a dot-bracket structure into a pairing map."""
    stack: list[int] = []
    pairs: dict[int, int] = {}
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket structure.")
            partner = stack.pop()
            pairs[index] = partner
            pairs[partner] = index
        elif char != ".":
            raise ValueError("Only '.', '(' and ')' are supported in dot-bracket structures.")
    if stack:
        raise ValueError("Unbalanced dot-bracket structure.")
    return pairs


def _unpaired_indicator(structure: str) -> np.ndarray:
    """Return 1 for unpaired sites and 0 for paired sites."""
    _dot_bracket_pairs(structure)
    return np.fromiter((1.0 if char == "." else 0.0 for char in structure), dtype=float)


def _safe_correlation(fn: Any, x: np.ndarray, y: np.ndarray) -> float:
    """Compute a correlation coefficient while suppressing constant-input warnings."""
    if x.size < 2 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(fn(x, y).statistic)


class SHAPEProcessor:
    """Read, normalize, and classify SHAPE reactivity data."""

    def __init__(self) -> None:
        self.raw_data: np.ndarray | None = None
        self.normalized_data: np.ndarray | None = None

    def load_data(self, filepath_or_data: Any) -> np.ndarray:
        """Load SHAPE records in position, nucleotide, reactivity format."""
        self.raw_data = _load_probe_data(filepath_or_data)
        self.normalized_data = None
        return self.raw_data.copy()

    def normalize(self, method: str = "28") -> np.ndarray:
        """Normalize SHAPE data using the 2-8% or box-plot method."""
        if self.raw_data is None:
            raise ValueError("No SHAPE data loaded. Call load_data() first.")
        values = _winsorize_nonnegative(self.raw_data["reactivity"])
        factor = _two_to_eight_factor(values) if method == "28" else None
        if method == "boxplot":
            factor = _boxplot_factor(values)
        if factor is None:
            raise ValueError("Unsupported normalization method. Use '28' or 'boxplot'.")
        normalized = self.raw_data.copy()
        normalized_values = values / factor
        normalized_values[~np.isfinite(normalized_values)] = np.nan
        normalized["reactivity"] = normalized_values
        self.normalized_data = normalized
        return normalized_values.copy()

    def classify_positions(self, thresholds: tuple[float, float] = (0.3, 0.7)) -> dict[str, np.ndarray]:
        """Classify positions into low, medium, high, and missing reactivity groups."""
        if thresholds[0] >= thresholds[1]:
            raise ValueError("The lower threshold must be smaller than the upper threshold.")
        data = self.normalized_data if self.normalized_data is not None else self.raw_data
        if data is None:
            raise ValueError("No SHAPE data loaded. Call load_data() first.")
        reactivities = data["reactivity"]
        positions = data["position"]
        finite_mask = np.isfinite(reactivities)
        return {
            "low": positions[finite_mask & (reactivities < thresholds[0])].copy(),
            "medium": positions[
                finite_mask & (reactivities >= thresholds[0]) & (reactivities <= thresholds[1])
            ].copy(),
            "high": positions[finite_mask & (reactivities > thresholds[1])].copy(),
            "missing": positions[~finite_mask].copy(),
        }


class DMSProcessor:
    """Read, normalize, and interpret DMS probing data."""

    def __init__(self) -> None:
        self.raw_data: np.ndarray | None = None
        self.normalized_data: np.ndarray | None = None

    def load_data(self, filepath_or_data: Any) -> np.ndarray:
        """Load DMS records in position, nucleotide, reactivity/count format."""
        self.raw_data = _load_probe_data(filepath_or_data)
        self.normalized_data = None
        return self.raw_data.copy()

    def normalize(self) -> np.ndarray:
        """Normalize DMS counts to reactivities with separate A/C scaling."""
        if self.raw_data is None:
            raise ValueError("No DMS data loaded. Call load_data() first.")
        normalized = self.raw_data.copy()
        clipped = _winsorize_nonnegative(self.raw_data["reactivity"])
        nucleotides = np.char.upper(self.raw_data["nucleotide"].astype("U1"))
        result = np.full(clipped.shape, np.nan, dtype=float)

        scales: dict[str, float] = {}
        for nucleotide in ("A", "C"):
            mask = (nucleotides == nucleotide) & np.isfinite(clipped)
            if np.any(mask):
                scales[nucleotide] = _two_to_eight_factor(clipped[mask])
                result[mask] = clipped[mask] / scales[nucleotide]

        informative = np.isfinite(result)
        if np.any(informative):
            target = float(np.nanmedian(result[informative]))
            for nucleotide in ("A", "C"):
                mask = (nucleotides == nucleotide) & np.isfinite(result)
                if np.any(mask):
                    median_value = float(np.nanmedian(result[mask]))
                    if median_value > 0:
                        result[mask] *= target / median_value

        normalized["reactivity"] = result
        self.normalized_data = normalized
        return result.copy()

    def get_accessibility(self, threshold: float = 0.5) -> np.ndarray:
        """Return a boolean accessibility mask from normalized DMS reactivities."""
        if self.normalized_data is None:
            self.normalize()
        assert self.normalized_data is not None
        reactivities = self.normalized_data["reactivity"]
        return (np.isfinite(reactivities) & (reactivities >= threshold)).astype(bool)


class ProbeConstrainedFolding:
    """Simplified probe-constrained RNA secondary-structure prediction."""

    def __init__(
        self,
        min_loop_length: int = 3,
        shape_slope: float = 1.8,
        shape_intercept: float = -0.6,
        temperature_celsius: float = 37.0,
    ) -> None:
        self.min_loop_length = min_loop_length
        self.shape_slope = shape_slope
        self.shape_intercept = shape_intercept
        self.temperature_kelvin = temperature_celsius + 273.15
        self.gas_constant = 0.0019872041
        self.last_partition_function: np.ndarray | None = None

    def fold_with_shape(
        self,
        sequence: str,
        shape_data: Any,
        method: str = "soft",
    ) -> tuple[str, float]:
        """Fold an RNA using SHAPE-derived hard, soft, or probabilistic constraints."""
        reactivities = self._prepare_probe_vector(sequence, shape_data, probe_type="shape")
        return self._fold(sequence, reactivities, method)

    def fold_with_dms(
        self,
        sequence: str,
        dms_data: Any,
        method: str = "soft",
    ) -> tuple[str, float]:
        """Fold an RNA using DMS-derived hard, soft, or probabilistic constraints."""
        reactivities = self._prepare_probe_vector(sequence, dms_data, probe_type="dms")
        return self._fold(sequence, reactivities, method)

    def _prepare_probe_vector(self, sequence: str, probe_data: Any, probe_type: str) -> np.ndarray:
        sequence = self._validate_sequence(sequence)
        n = len(sequence)

        if probe_type == "shape" and isinstance(probe_data, SHAPEProcessor):
            if probe_data.raw_data is None:
                raise ValueError("The SHAPEProcessor has no loaded data.")
            if probe_data.normalized_data is None:
                probe_data.normalize()
            assert probe_data.normalized_data is not None
            return self._vector_from_structured(probe_data.normalized_data, n)

        if probe_type == "dms" and isinstance(probe_data, DMSProcessor):
            if probe_data.raw_data is None:
                raise ValueError("The DMSProcessor has no loaded data.")
            if probe_data.normalized_data is None:
                probe_data.normalize()
            assert probe_data.normalized_data is not None
            return self._vector_from_structured(probe_data.normalized_data, n, sequence, informative={"A", "C"})

        if isinstance(probe_data, (str, Path)):
            processor = SHAPEProcessor() if probe_type == "shape" else DMSProcessor()
            processor.load_data(probe_data)
            processor.normalize()  # type: ignore[call-arg]
            prepared = processor.normalized_data
            assert prepared is not None
            informative = {"A", "C"} if probe_type == "dms" else None
            return self._vector_from_structured(prepared, n, sequence, informative)

        if isinstance(probe_data, np.ndarray) and probe_data.dtype.names:
            informative = {"A", "C"} if probe_type == "dms" else None
            return self._vector_from_structured(probe_data, n, sequence, informative)

        if isinstance(probe_data, Mapping):
            vector = np.full(n, np.nan, dtype=float)
            for key, value in probe_data.items():
                index = int(key) - 1
                if not 0 <= index < n:
                    raise ValueError("Probe position is outside the sequence length.")
                vector[index] = _parse_reactivity(value)
            if probe_type == "dms":
                informative_mask = np.fromiter((nt in {"A", "C"} for nt in sequence), dtype=bool)
                vector[~informative_mask] = np.nan
            return vector

        array = np.asarray(probe_data, dtype=float)
        if array.ndim != 1 or array.size != n:
            raise ValueError("Probe reactivities must be a 1D array matching the sequence length.")
        array = array.astype(float, copy=True)
        array[array < 0] = np.nan
        if probe_type == "dms":
            informative_mask = np.fromiter((nt in {"A", "C"} for nt in sequence), dtype=bool)
            array[~informative_mask] = np.nan
        return array

    def _vector_from_structured(
        self,
        data: np.ndarray,
        sequence_length: int,
        sequence: str | None = None,
        informative: set[str] | None = None,
    ) -> np.ndarray:
        vector = np.full(sequence_length, np.nan, dtype=float)
        for row in data:
            index = int(row["position"]) - 1
            if not 0 <= index < sequence_length:
                raise ValueError("Probe position is outside the sequence length.")
            vector[index] = _parse_reactivity(row["reactivity"])
        if sequence is not None and informative is not None:
            mask = np.fromiter((nt in informative for nt in sequence), dtype=bool)
            vector[~mask] = np.nan
        return vector

    def _validate_sequence(self, sequence: str) -> str:
        sequence = sequence.upper().replace("T", "U")
        invalid = set(sequence) - {"A", "C", "G", "U"}
        if invalid:
            raise ValueError(f"Sequence contains invalid nucleotides: {sorted(invalid)}")
        return sequence

    def _can_pair(self, left: str, right: str) -> bool:
        return left + right in PAIR_ENERGIES

    def _base_pair_cost(self, left: str, right: str, span: int) -> float:
        return PAIR_ENERGIES[left + right] + 0.2 * math.log(span + 1.0)

    def _pseudo_energy(self, reactivity: float) -> float:
        if not np.isfinite(reactivity):
            return 0.0
        return self.shape_slope * math.log(reactivity + 1.0) + self.shape_intercept

    def _probabilistic_terms(self, reactivities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        unpaired = np.zeros_like(reactivities, dtype=float)
        paired = np.zeros_like(reactivities, dtype=float)
        valid = np.isfinite(reactivities)
        if not np.any(valid):
            return unpaired, paired
        probabilities = np.clip(reactivities[valid] / (reactivities[valid] + 1.0), 0.05, 0.95)
        rt = self.gas_constant * self.temperature_kelvin
        unpaired[valid] = -rt * np.log(probabilities)
        paired[valid] = -rt * np.log(1.0 - probabilities)
        return unpaired, paired

    def _fold(self, sequence: str, reactivities: np.ndarray, method: str) -> tuple[str, float]:
        method = method.lower()
        if method not in {"hard", "soft", "probabilistic"}:
            raise ValueError("method must be one of 'hard', 'soft', or 'probabilistic'.")

        sequence = self._validate_sequence(sequence)
        n = len(sequence)
        if reactivities.shape != (n,):
            raise ValueError("Probe reactivities must match the sequence length.")
        if n == 0:
            self.last_partition_function = np.empty((0, 0), dtype=float)
            return "", 0.0

        large_penalty = 1e6
        forced_unpaired = np.isfinite(reactivities) & (reactivities > 0.7)
        forced_paired = np.isfinite(reactivities) & (reactivities < 0.3)
        unpaired_costs = np.zeros(n, dtype=float)
        pair_prior_costs = np.zeros(n, dtype=float)
        if method == "hard":
            unpaired_costs[forced_paired] = large_penalty
        elif method == "probabilistic":
            unpaired_costs, pair_prior_costs = self._probabilistic_terms(reactivities)

        dp = np.zeros((n, n), dtype=float)
        choice: list[list[tuple[str, int | None]]] = [[("unpaired", None) for _ in range(n)] for _ in range(n)]

        for i in range(n - 1, -1, -1):
            for j in range(i, n):
                best = unpaired_costs[i] + (dp[i + 1, j] if i + 1 <= j else 0.0)
                best_choice: tuple[str, int | None] = ("unpaired", None)

                for k in range(i + self.min_loop_length + 1, j + 1):
                    if not self._can_pair(sequence[i], sequence[k]):
                        continue
                    if method == "hard" and (forced_unpaired[i] or forced_unpaired[k]):
                        continue

                    base_cost = self._base_pair_cost(sequence[i], sequence[k], k - i)
                    if method == "soft":
                        base_cost += self._pseudo_energy(reactivities[i]) + self._pseudo_energy(reactivities[k])
                    elif method == "probabilistic":
                        base_cost += pair_prior_costs[i] + pair_prior_costs[k]

                    inside = dp[i + 1, k - 1] if i + 1 <= k - 1 else 0.0
                    outside = dp[k + 1, j] if k + 1 <= j else 0.0
                    candidate = base_cost + inside + outside
                    if candidate < best:
                        best = candidate
                        best_choice = ("pair", k)

                dp[i, j] = best
                choice[i][j] = best_choice

        structure = ["."] * n
        self._traceback(0, n - 1, choice, structure)

        if method == "probabilistic":
            self.last_partition_function = self._compute_partition_function(
                sequence,
                reactivities,
                unpaired_costs,
                pair_prior_costs,
                forced_unpaired,
            )
        else:
            self.last_partition_function = None

        return "".join(structure), float(dp[0, n - 1])

    def _traceback(
        self,
        i: int,
        j: int,
        choice: list[list[tuple[str, int | None]]],
        structure: list[str],
    ) -> None:
        if i > j:
            return
        if i == j:
            structure[i] = "."
            return
        action, partner = choice[i][j]
        if action == "unpaired" or partner is None:
            structure[i] = "."
            self._traceback(i + 1, j, choice, structure)
            return
        structure[i] = "("
        structure[partner] = ")"
        self._traceback(i + 1, partner - 1, choice, structure)
        self._traceback(partner + 1, j, choice, structure)

    def _compute_partition_function(
        self,
        sequence: str,
        reactivities: np.ndarray,
        unpaired_costs: np.ndarray,
        pair_prior_costs: np.ndarray,
        forced_unpaired: np.ndarray,
    ) -> np.ndarray:
        n = len(sequence)
        z = np.zeros((n, n), dtype=float)
        beta = 1.0 / (self.gas_constant * self.temperature_kelvin)

        def boltzmann(energy: float) -> float:
            capped = float(np.clip(energy, -50.0, 50.0))
            return math.exp(-beta * capped)

        for i in range(n - 1, -1, -1):
            for j in range(i, n):
                total = boltzmann(unpaired_costs[i]) * (z[i + 1, j] if i + 1 <= j else 1.0)
                for k in range(i + self.min_loop_length + 1, j + 1):
                    if forced_unpaired[i] or forced_unpaired[k]:
                        continue
                    if not self._can_pair(sequence[i], sequence[k]):
                        continue
                    energy = self._base_pair_cost(sequence[i], sequence[k], k - i)
                    energy += pair_prior_costs[i] + pair_prior_costs[k]
                    left = z[i + 1, k - 1] if i + 1 <= k - 1 else 1.0
                    right = z[k + 1, j] if k + 1 <= j else 1.0
                    total += boltzmann(energy) * left * right
                z[i, j] = total
        return z


class ICEFold:
    """Iterative constraint enforcement for probe-constrained RNA folding."""

    def __init__(self, folder: ProbeConstrainedFolding | None = None) -> None:
        self.folder = folder or ProbeConstrainedFolding()

    def fold_iterative(
        self,
        sequence: str,
        probe_data: Any,
        max_iter: int = 10,
    ) -> tuple[str, float, list[dict[str, float | int | str]]]:
        """Alternately fold and update effective constraints until the structure stabilizes."""
        if max_iter < 1:
            raise ValueError("max_iter must be at least 1.")
        sequence = self.folder._validate_sequence(sequence)
        is_dms = isinstance(probe_data, DMSProcessor)
        working = self.folder._prepare_probe_vector(sequence, probe_data, "dms" if is_dms else "shape")
        original = working.copy()
        history: list[dict[str, float | int | str]] = []
        previous_structure: str | None = None
        structure = "." * len(sequence)
        energy = 0.0

        for iteration in range(1, max_iter + 1):
            if is_dms:
                structure, energy = self.folder.fold_with_dms(sequence, working, method="soft")
            else:
                structure, energy = self.folder.fold_with_shape(sequence, working, method="soft")
            agreement = self._agreement_score(structure, original)
            history.append(
                {
                    "iteration": iteration,
                    "energy": float(energy),
                    "agreement": float(agreement),
                    "structure": structure,
                }
            )
            if structure == previous_structure:
                break
            predicted = _unpaired_indicator(structure)
            valid = np.isfinite(original)
            updated = working.copy()
            updated[valid] = 0.75 * original[valid] + 0.25 * predicted[valid]
            disagreement = np.abs(predicted - np.clip(np.nan_to_num(original, nan=0.5), 0.0, 1.0))
            updated[valid & (disagreement > 0.5)] = 0.5 * (
                updated[valid & (disagreement > 0.5)] + predicted[valid & (disagreement > 0.5)]
            )
            working = updated
            previous_structure = structure

        return structure, float(energy), history

    def _agreement_score(self, structure: str, reactivities: np.ndarray) -> float:
        predicted = _unpaired_indicator(structure)
        valid = np.isfinite(reactivities)
        if not np.any(valid):
            return float("nan")
        target = np.clip(reactivities[valid], 0.0, 1.0)
        return float(1.0 - np.mean(np.abs(predicted[valid] - target)))


class ProbeDataSimulator:
    """Simulate SHAPE and DMS probing data from a known RNA structure."""

    def __init__(self, seed: int | None = None) -> None:
        self.rng = np.random.default_rng(seed)

    def simulate_shape(
        self,
        sequence: str,
        structure: str,
        noise_level: float = 0.2,
    ) -> np.ndarray:
        """Generate synthetic SHAPE data with log-normal noise."""
        sequence = self._validate_inputs(sequence, structure)
        pairs = _dot_bracket_pairs(structure)
        data = np.empty(len(sequence), dtype=PROBE_DTYPE)
        for index, nucleotide in enumerate(sequence, start=1):
            paired = (index - 1) in pairs
            base = self.rng.lognormal(mean=-1.5, sigma=0.35) if paired else self.rng.lognormal(mean=-0.1, sigma=0.45)
            noise = self.rng.lognormal(mean=0.0, sigma=max(noise_level, 1e-6))
            data[index - 1] = (index, nucleotide, float(base * noise))
        return data

    def simulate_dms(
        self,
        sequence: str,
        structure: str,
        noise_level: float = 0.2,
    ) -> np.ndarray:
        """Generate synthetic DMS data with A/C-specific accessibility patterns."""
        sequence = self._validate_inputs(sequence, structure)
        pairs = _dot_bracket_pairs(structure)
        data = np.empty(len(sequence), dtype=PROBE_DTYPE)
        for index, nucleotide in enumerate(sequence, start=1):
            paired = (index - 1) in pairs
            if nucleotide in {"A", "C"}:
                mean = -2.0 if paired else -0.2
                sigma = 0.35 if paired else 0.5
            else:
                mean = -3.0
                sigma = 0.25
            base = self.rng.lognormal(mean=mean, sigma=sigma)
            noise = self.rng.lognormal(mean=0.0, sigma=max(noise_level, 1e-6))
            data[index - 1] = (index, nucleotide, float(base * noise))
        return data

    def _validate_inputs(self, sequence: str, structure: str) -> str:
        sequence = sequence.upper().replace("T", "U")
        if len(sequence) != len(structure):
            raise ValueError("sequence and structure must have the same length.")
        invalid = set(sequence) - {"A", "C", "G", "U"}
        if invalid:
            raise ValueError(f"Sequence contains invalid nucleotides: {sorted(invalid)}")
        _dot_bracket_pairs(structure)
        return sequence


class ProbeEvaluator:
    """Evaluate agreement between probe data and a predicted RNA structure."""

    def compute_auc(self, structure: str, probe_data: Any) -> float:
        """Compute ROC AUC for classifying unpaired positions from probe reactivities."""
        indicator = _unpaired_indicator(structure)
        reactivities = self._prepare_probe_vector(probe_data, len(structure))
        valid = np.isfinite(reactivities)
        labels = indicator[valid]
        scores = reactivities[valid]
        if labels.size == 0 or np.unique(labels).size < 2:
            return float("nan")
        ranks = stats.rankdata(scores)
        positives = labels == 1.0
        n_pos = int(np.sum(positives))
        n_neg = int(labels.size - n_pos)
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        auc = (np.sum(ranks[positives]) - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
        return float(auc)

    def compute_correlation(self, structure: str, probe_data: Any) -> tuple[float, float]:
        """Compute Pearson and Spearman correlation with unpaired status."""
        indicator = _unpaired_indicator(structure)
        reactivities = self._prepare_probe_vector(probe_data, len(structure))
        valid = np.isfinite(reactivities)
        if np.sum(valid) < 2:
            return float("nan"), float("nan")
        x = reactivities[valid]
        y = indicator[valid]
        pearson_r = _safe_correlation(stats.pearsonr, x, y)
        spearman_r = _safe_correlation(stats.spearmanr, x, y)
        return pearson_r, spearman_r

    def _prepare_probe_vector(self, probe_data: Any, expected_length: int) -> np.ndarray:
        if isinstance(probe_data, (SHAPEProcessor, DMSProcessor)):
            data = probe_data.normalized_data if probe_data.normalized_data is not None else probe_data.raw_data
            if data is None:
                raise ValueError("The processor has no loaded data.")
            probe_data = data
        if isinstance(probe_data, np.ndarray) and probe_data.dtype.names:
            vector = np.full(expected_length, np.nan, dtype=float)
            for row in probe_data:
                index = int(row["position"]) - 1
                if not 0 <= index < expected_length:
                    raise ValueError("Probe position is outside the structure length.")
                vector[index] = _parse_reactivity(row["reactivity"])
            return vector
        array = np.asarray(probe_data, dtype=float)
        if array.ndim != 1 or array.size != expected_length:
            raise ValueError("Probe data must be a 1D array matching the structure length.")
        return array


__all__ = [
    "DMSProcessor",
    "ICEFold",
    "ProbeConstrainedFolding",
    "ProbeDataSimulator",
    "ProbeEvaluator",
    "SHAPEProcessor",
]
