#!/usr/bin/env python3
"""
Hybrid integration module for combining short-read and long-read SV calls.
Improves precision through cross-validation and breakpoint refinement.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class DataSource(Enum):
    LONG_READ = "long_read"
    SHORT_READ = "short_read"
    HYBRID = "hybrid"


@dataclass
class HybridSVCall:
    """SV call with evidence from multiple sequencing platforms."""
    sv_type: str
    chrom: str
    start: int
    end: int
    size: int
    quality: float
    genotype: str
    long_read_support: int
    short_read_support: int
    source: DataSource
    concordance_score: float = 0.0
    refined_start: Optional[int] = None
    refined_end: Optional[int] = None
    info: Dict = field(default_factory=dict)


@dataclass
class ShortReadEvidence:
    """Evidence from short-read sequencing for SV validation."""
    chrom: str
    start: int
    end: int
    sv_type: str
    split_reads: int
    discordant_pairs: int
    read_depth_ratio: float
    mapq_mean: float


class BreakpointRefiner:
    """Refine SV breakpoints using short-read precision."""

    def __init__(self, search_window: int = 200):
        self.search_window = search_window

    def refine(
        self, long_read_sv: Dict, short_read_evidence: List[ShortReadEvidence]
    ) -> Tuple[int, int]:
        best_start = long_read_sv["start"]
        best_end = long_read_sv["end"]
        best_score = 0

        for evidence in short_read_evidence:
            if evidence.sv_type != long_read_sv.get("sv_type", ""):
                continue
            if abs(evidence.start - long_read_sv["start"]) > self.search_window:
                continue

            score = (
                evidence.split_reads * 2 +
                evidence.discordant_pairs +
                evidence.mapq_mean / 60.0 * 10
            )
            if score > best_score:
                best_score = score
                best_start = evidence.start
                best_end = evidence.end

        return best_start, best_end


class ConcordanceCalculator:
    """Calculate concordance between long-read and short-read SV calls."""

    def __init__(
        self, position_tolerance: int = 500,
        size_tolerance: float = 0.25,
        reciprocal_overlap: float = 0.5
    ):
        self.position_tolerance = position_tolerance
        self.size_tolerance = size_tolerance
        self.reciprocal_overlap = reciprocal_overlap

    def compute_concordance(
        self, lr_call: Dict, sr_call: ShortReadEvidence
    ) -> float:
        if lr_call.get("sv_type") != sr_call.sv_type:
            return 0.0
        if lr_call.get("chrom") != sr_call.chrom:
            return 0.0

        # Position concordance
        start_diff = abs(lr_call["start"] - sr_call.start)
        end_diff = abs(lr_call["end"] - sr_call.end)
        pos_score = max(0, 1.0 - (start_diff + end_diff) / (2 * self.position_tolerance))

        # Size concordance
        lr_size = lr_call["end"] - lr_call["start"]
        sr_size = sr_call.end - sr_call.start
        if max(lr_size, sr_size) > 0:
            size_ratio = min(lr_size, sr_size) / max(lr_size, sr_size)
        else:
            size_ratio = 1.0
        size_score = size_ratio

        # Reciprocal overlap
        overlap_start = max(lr_call["start"], sr_call.start)
        overlap_end = min(lr_call["end"], sr_call.end)
        overlap = max(0, overlap_end - overlap_start)
        lr_overlap = overlap / max(lr_size, 1)
        sr_overlap = overlap / max(sr_size, 1)
        overlap_score = min(lr_overlap, sr_overlap) / self.reciprocal_overlap
        overlap_score = min(overlap_score, 1.0)

        return (pos_score * 0.3 + size_score * 0.3 + overlap_score * 0.4)

    def find_concordant_pairs(
        self, lr_calls: List[Dict], sr_calls: List[ShortReadEvidence]
    ) -> List[Tuple[int, int, float]]:
        pairs = []
        for i, lr in enumerate(lr_calls):
            for j, sr in enumerate(sr_calls):
                score = self.compute_concordance(lr, sr)
                if score > 0.3:
                    pairs.append((i, j, score))

        pairs.sort(key=lambda x: -x[2])

        used_lr: set = set()
        used_sr: set = set()
        matched = []
        for i, j, score in pairs:
            if i not in used_lr and j not in used_sr:
                matched.append((i, j, score))
                used_lr.add(i)
                used_sr.add(j)

        return matched


class HybridSVIntegrator:
    """Integrate long-read and short-read SV calls."""

    def __init__(
        self, concordance_threshold: float = 0.5,
        lr_only_min_support: int = 5,
        sr_validation_boost: float = 10.0
    ):
        self.refiner = BreakpointRefiner()
        self.concordance_calc = ConcordanceCalculator()
        self.concordance_threshold = concordance_threshold
        self.lr_only_min_support = lr_only_min_support
        self.sr_validation_boost = sr_validation_boost

    def integrate(
        self, lr_calls: List[Dict], sr_evidence: List[ShortReadEvidence]
    ) -> List[HybridSVCall]:
        matched_pairs = self.concordance_calc.find_concordant_pairs(
            lr_calls, sr_evidence
        )
        matched_lr = {i for i, _, _ in matched_pairs}
        matched_sr = {j for _, j, _ in matched_pairs}

        results = []

        # Process matched pairs (concordant calls)
        for lr_idx, sr_idx, concordance in matched_pairs:
            lr = lr_calls[lr_idx]
            sr = sr_evidence[sr_idx]
            refined_start, refined_end = self.refiner.refine(lr, [sr])

            results.append(HybridSVCall(
                sv_type=lr.get("sv_type", "UNK"),
                chrom=lr.get("chrom", ""),
                start=lr["start"], end=lr["end"],
                size=lr["end"] - lr["start"],
                quality=lr.get("quality", 0) + self.sr_validation_boost,
                genotype=lr.get("genotype", "./."),
                long_read_support=lr.get("support", 0),
                short_read_support=sr.split_reads + sr.discordant_pairs,
                source=DataSource.HYBRID,
                concordance_score=concordance,
                refined_start=refined_start,
                refined_end=refined_end,
            ))

        # Long-read-only calls (unique to long reads)
        for i, lr in enumerate(lr_calls):
            if i in matched_lr:
                continue
            if lr.get("support", 0) >= self.lr_only_min_support:
                results.append(HybridSVCall(
                    sv_type=lr.get("sv_type", "UNK"),
                    chrom=lr.get("chrom", ""),
                    start=lr["start"], end=lr["end"],
                    size=lr["end"] - lr["start"],
                    quality=lr.get("quality", 0),
                    genotype=lr.get("genotype", "./."),
                    long_read_support=lr.get("support", 0),
                    short_read_support=0,
                    source=DataSource.LONG_READ,
                ))

        # Short-read-only calls (high-confidence)
        for j, sr in enumerate(sr_evidence):
            if j in matched_sr:
                continue
            total_support = sr.split_reads + sr.discordant_pairs
            if total_support >= 10 and sr.mapq_mean >= 40:
                results.append(HybridSVCall(
                    sv_type=sr.sv_type,
                    chrom=sr.chrom,
                    start=sr.start, end=sr.end,
                    size=sr.end - sr.start,
                    quality=sr.mapq_mean,
                    genotype="0/1",
                    long_read_support=0,
                    short_read_support=total_support,
                    source=DataSource.SHORT_READ,
                ))

        results.sort(key=lambda x: (x.chrom, x.start))
        return results
