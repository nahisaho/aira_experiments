"""Repeat-aware processing utilities for repetitive long-read SV loci."""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np


class RepeatRegionHandler:
    """Handle telomeric, centromeric and tandem-repeat-rich regions."""

    def __init__(self, telomere_motifs: Sequence[str] = ("TTAGGG", "CCCTAA")) -> None:
        self.telomere_motifs = tuple(motif.upper() for motif in telomere_motifs)

    def detect_telomere_repeats(self, sequence: str, min_copies: int = 4) -> List[Dict[str, Any]]:
        """Detect canonical telomere motifs (TTAGGG/CCCTAA) and report runs."""

        seq = sequence.upper()
        hits: List[Dict[str, Any]] = []
        for motif in self.telomere_motifs:
            motif_length = len(motif)
            index = 0
            while index <= len(seq) - motif_length * min_copies:
                copies = 0
                while seq[index + copies * motif_length : index + (copies + 1) * motif_length] == motif:
                    copies += 1
                if copies >= min_copies:
                    start = index
                    end = index + copies * motif_length
                    hits.append(
                        {
                            "motif": motif,
                            "start": start,
                            "end": end,
                            "copies": copies,
                            "fraction": (copies * motif_length) / max(len(seq), 1),
                        }
                    )
                    index = end
                else:
                    index += 1
        return hits

    def analyze_centromere(self, sequence: str, hor_unit_range: Tuple[int, int] = (150, 210)) -> Dict[str, Any]:
        """Estimate alpha-satellite enrichment and higher-order repeat structure."""

        seq = sequence.upper()
        if not seq:
            return {"alpha_satellite_score": 0.0, "estimated_hor_unit": None, "monomer_similarity": 0.0}

        monomer_scores = []
        for unit_size in range(hor_unit_range[0], hor_unit_range[1] + 1):
            blocks = [seq[index : index + unit_size] for index in range(0, len(seq) - unit_size + 1, unit_size)]
            if len(blocks) < 2:
                continue
            similarity = np.mean([self._sequence_identity(blocks[0], block) for block in blocks[1:]])
            monomer_scores.append((unit_size, float(similarity)))

        if not monomer_scores:
            return {"alpha_satellite_score": 0.0, "estimated_hor_unit": None, "monomer_similarity": 0.0}

        monomer_length, monomer_similarity = max(monomer_scores, key=lambda item: item[1])
        periodicities = self._autocorrelation_periods(seq, max_period=monomer_length * 20)
        hor_length = periodicities[0][0] if periodicities else monomer_length
        alpha_satellite_score = min(1.0, monomer_similarity * (hor_length / max(monomer_length, 1)) / 8.0)
        return {
            "alpha_satellite_score": alpha_satellite_score,
            "estimated_monomer_length": monomer_length,
            "estimated_hor_unit": hor_length,
            "monomer_similarity": monomer_similarity,
            "top_periodicities": periodicities[:5],
        }

    def detect_tandem_expansion(
        self,
        sequence: str,
        motif_lengths: Sequence[int] = (1, 2, 3, 4, 5, 6),
        min_repeats: int = 5,
    ) -> List[Dict[str, Any]]:
        """Detect candidate STR/VNTR expansions by scanning repeated motifs."""

        seq = sequence.upper()
        expansions: List[Dict[str, Any]] = []
        for motif_length in motif_lengths:
            for start in range(0, len(seq) - motif_length * min_repeats + 1):
                motif = seq[start : start + motif_length]
                if len(set(motif)) == 1 and motif_length > 3:
                    continue
                repeats = 1
                while seq[start + repeats * motif_length : start + (repeats + 1) * motif_length] == motif:
                    repeats += 1
                if repeats >= min_repeats:
                    end = start + repeats * motif_length
                    expansions.append(
                        {
                            "motif": motif,
                            "motif_length": motif_length,
                            "start": start,
                            "end": end,
                            "repeat_count": repeats,
                            "expansion_size": repeats * motif_length,
                        }
                    )
        return self._deduplicate_regions(expansions)

    def kmer_frequency_filter(
        self,
        sequence: str,
        candidate_score: float,
        k: int = 5,
        entropy_threshold: float = 1.5,
    ) -> Dict[str, Any]:
        """Filter repetitive false positives using k-mer complexity statistics."""

        seq = sequence.upper()
        if len(seq) < k:
            return {"pass_filter": True, "adjusted_score": candidate_score, "entropy": 0.0}
        kmers = [seq[index : index + k] for index in range(0, len(seq) - k + 1)]
        counts = Counter(kmers)
        total = sum(counts.values())
        probabilities = np.asarray([count / total for count in counts.values()], dtype=np.float64)
        entropy = float(-np.sum(probabilities * np.log2(probabilities + 1e-12)))
        dominant_fraction = max(probabilities) if probabilities.size else 0.0
        adjusted_score = candidate_score * max(entropy / max(entropy_threshold, 1e-6), 0.25) * (1.0 - 0.5 * dominant_fraction)
        return {
            "pass_filter": entropy >= entropy_threshold or dominant_fraction < 0.7,
            "adjusted_score": adjusted_score,
            "entropy": entropy,
            "dominant_kmer_fraction": float(dominant_fraction),
        }

    def repeat_aware_alignment_score(
        self,
        query: str,
        target: str,
        match_score: float = 2.0,
        mismatch_penalty: float = -2.5,
        gap_penalty: float = -3.0,
    ) -> float:
        """Compute a simple repeat-aware alignment score.

        Gaps traversing tandem-repeat tracts are penalized less harshly to avoid
        overcalling rearrangements in low-complexity sequence.
        """

        query_upper = query.upper()
        target_upper = target.upper()
        rows = len(query_upper) + 1
        cols = len(target_upper) + 1
        dp = np.zeros((rows, cols), dtype=np.float64)

        for index in range(1, rows):
            dp[index, 0] = dp[index - 1, 0] + gap_penalty
        for index in range(1, cols):
            dp[0, index] = dp[0, index - 1] + gap_penalty

        for row in range(1, rows):
            for col in range(1, cols):
                repeat_bonus = 0.6 if self._is_repetitive_context(query_upper, row - 1) or self._is_repetitive_context(target_upper, col - 1) else 1.0
                match = dp[row - 1, col - 1] + (match_score if query_upper[row - 1] == target_upper[col - 1] else mismatch_penalty)
                delete = dp[row - 1, col] + gap_penalty * repeat_bonus
                insert = dp[row, col - 1] + gap_penalty * repeat_bonus
                dp[row, col] = max(match, delete, insert)
        return float(dp[-1, -1])

    @staticmethod
    def _sequence_identity(first: str, second: str) -> float:
        length = min(len(first), len(second))
        if length == 0:
            return 0.0
        matches = sum(1 for left, right in zip(first[:length], second[:length]) if left == right)
        return matches / length

    @staticmethod
    def _autocorrelation_periods(sequence: str, max_period: int) -> List[Tuple[int, float]]:
        encoded = np.fromiter((ord(base) for base in sequence), dtype=np.float64)
        periods: List[Tuple[int, float]] = []
        for period in range(1, min(max_period, len(sequence) // 2) + 1):
            corr = np.corrcoef(encoded[:-period], encoded[period:])[0, 1]
            if np.isnan(corr):
                continue
            periods.append((period, float(corr)))
        return sorted(periods, key=lambda item: item[1], reverse=True)

    @staticmethod
    def _deduplicate_regions(expansions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        deduplicated: List[Dict[str, Any]] = []
        for expansion in sorted(expansions, key=lambda item: (item["start"], -item["expansion_size"])):
            if deduplicated and expansion["start"] < deduplicated[-1]["end"] and expansion["motif"] == deduplicated[-1]["motif"]:
                continue
            deduplicated.append(expansion)
        return deduplicated

    @staticmethod
    def _is_repetitive_context(sequence: str, index: int, window: int = 6) -> bool:
        start = max(0, index - window)
        end = min(len(sequence), index + window + 1)
        segment = sequence[start:end]
        return len(set(segment)) <= max(2, len(segment) // 4)


__all__ = ["RepeatRegionHandler"]
