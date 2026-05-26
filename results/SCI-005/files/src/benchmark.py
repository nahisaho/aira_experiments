#!/usr/bin/env python3
"""
Benchmarking module for evaluating SV calls against GIAB Tier1 truth set.
Implements precision, recall, F1, and stratified evaluation metrics.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
import logging

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkResult:
    """Results from benchmarking against a truth set."""
    total_truth: int
    total_calls: int
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1_score: float
    genotype_concordance: float
    stratified_results: Dict[str, Dict] = field(default_factory=dict)


@dataclass
class TruthVariant:
    """A variant from the GIAB truth set."""
    chrom: str
    start: int
    end: int
    sv_type: str
    size: int
    genotype: str
    tier: str = "Tier1"
    region_class: str = "all"


class GIABEvaluator:
    """Evaluate SV calls against GIAB HG002 Tier1 truth set."""

    def __init__(
        self, position_tolerance: int = 1000,
        size_tolerance: float = 0.25,
        min_reciprocal_overlap: float = 0.5
    ):
        self.position_tolerance = position_tolerance
        self.size_tolerance = size_tolerance
        self.min_reciprocal_overlap = min_reciprocal_overlap

    def evaluate(
        self, calls: List[Dict], truth: List[TruthVariant]
    ) -> BenchmarkResult:
        matched_truth = set()
        matched_calls = set()
        gt_concordant = 0

        for i, call in enumerate(calls):
            for j, tv in enumerate(truth):
                if j in matched_truth:
                    continue
                if self._is_match(call, tv):
                    matched_truth.add(j)
                    matched_calls.add(i)
                    if call.get("genotype") == tv.genotype:
                        gt_concordant += 1
                    break

        tp = len(matched_truth)
        fp = len(calls) - len(matched_calls)
        fn = len(truth) - len(matched_truth)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        gt_conc = gt_concordant / tp if tp > 0 else 0.0

        # Stratified evaluation
        stratified = self._stratified_evaluation(calls, truth, matched_truth, matched_calls)

        return BenchmarkResult(
            total_truth=len(truth),
            total_calls=len(calls),
            true_positives=tp,
            false_positives=fp,
            false_negatives=fn,
            precision=precision,
            recall=recall,
            f1_score=f1,
            genotype_concordance=gt_conc,
            stratified_results=stratified
        )

    def _is_match(self, call: Dict, truth: TruthVariant) -> bool:
        if call.get("chrom") != truth.chrom:
            return False
        if call.get("sv_type") != truth.sv_type:
            return False

        # Position check
        start_diff = abs(call.get("start", 0) - truth.start)
        end_diff = abs(call.get("end", 0) - truth.end)
        if start_diff > self.position_tolerance or end_diff > self.position_tolerance:
            return False

        # Size check
        call_size = abs(call.get("size", call.get("end", 0) - call.get("start", 0)))
        if max(call_size, truth.size) > 0:
            size_ratio = min(call_size, truth.size) / max(call_size, truth.size)
            if size_ratio < (1 - self.size_tolerance):
                return False

        # Reciprocal overlap
        overlap_start = max(call.get("start", 0), truth.start)
        overlap_end = min(call.get("end", 0), truth.end)
        overlap = max(0, overlap_end - overlap_start)
        call_len = max(call.get("end", 0) - call.get("start", 0), 1)
        truth_len = max(truth.end - truth.start, 1)

        if min(overlap / call_len, overlap / truth_len) < self.min_reciprocal_overlap:
            return False

        return True

    def _stratified_evaluation(
        self, calls: List[Dict], truth: List[TruthVariant],
        matched_truth: set, matched_calls: set
    ) -> Dict[str, Dict]:
        strata = {}

        # By SV type
        for sv_type in ["DEL", "INS", "DUP", "INV"]:
            type_truth = [i for i, t in enumerate(truth) if t.sv_type == sv_type]
            type_calls = [i for i, c in enumerate(calls) if c.get("sv_type") == sv_type]
            type_tp = len(set(type_truth) & matched_truth)
            type_fp = len(set(type_calls) - matched_calls)
            type_fn = len(set(type_truth) - matched_truth)
            p = type_tp / (type_tp + type_fp) if (type_tp + type_fp) > 0 else 0
            r = type_tp / (type_tp + type_fn) if (type_tp + type_fn) > 0 else 0
            f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0
            strata[f"type_{sv_type}"] = {
                "precision": p, "recall": r, "f1": f1,
                "tp": type_tp, "fp": type_fp, "fn": type_fn
            }

        # By size range
        size_ranges = [
            ("50-300bp", 50, 300),
            ("300-1000bp", 300, 1000),
            ("1-10kb", 1000, 10000),
            ("10-100kb", 10000, 100000),
            (">100kb", 100000, float('inf'))
        ]
        for label, min_sz, max_sz in size_ranges:
            sz_truth = [i for i, t in enumerate(truth) if min_sz <= t.size < max_sz]
            sz_calls = [i for i, c in enumerate(calls)
                       if min_sz <= abs(c.get("size", 0)) < max_sz]
            sz_tp = len(set(sz_truth) & matched_truth)
            sz_fp = len(set(sz_calls) - matched_calls)
            sz_fn = len(set(sz_truth) - matched_truth)
            p = sz_tp / (sz_tp + sz_fp) if (sz_tp + sz_fp) > 0 else 0
            r = sz_tp / (sz_tp + sz_fn) if (sz_tp + sz_fn) > 0 else 0
            f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0
            strata[f"size_{label}"] = {
                "precision": p, "recall": r, "f1": f1,
                "tp": sz_tp, "fp": sz_fp, "fn": sz_fn
            }

        return strata


class BenchmarkSimulator:
    """Simulate benchmark data for pipeline validation."""

    def __init__(self, seed: int = 42):
        self.rng = np.random.RandomState(seed)

    def generate_truth_set(self, n_variants: int = 500) -> List[TruthVariant]:
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
        sv_types = ["DEL", "INS", "DUP", "INV"]
        type_weights = [0.4, 0.35, 0.15, 0.1]

        truth = []
        for _ in range(n_variants):
            chrom = self.rng.choice(chroms)
            sv_type = self.rng.choice(sv_types, p=type_weights)
            size = int(self.rng.lognormal(mean=7, sigma=1.5))
            size = max(50, min(size, 500000))
            start = self.rng.randint(1000000, 200000000)
            gt = self.rng.choice(["0/1", "1/1"], p=[0.6, 0.4])

            truth.append(TruthVariant(
                chrom=chrom, start=start, end=start + size,
                sv_type=sv_type, size=size, genotype=gt
            ))

        return truth

    def simulate_calls(
        self, truth: List[TruthVariant],
        sensitivity: float = 0.85,
        precision_rate: float = 0.90,
        position_noise: int = 100
    ) -> List[Dict]:
        calls = []

        # True positives (with noise)
        for tv in truth:
            if self.rng.random() < sensitivity:
                noise_start = self.rng.randint(-position_noise, position_noise)
                noise_end = self.rng.randint(-position_noise, position_noise)
                calls.append({
                    "chrom": tv.chrom,
                    "start": tv.start + noise_start,
                    "end": tv.end + noise_end,
                    "sv_type": tv.sv_type,
                    "size": tv.size + noise_start - noise_end,
                    "quality": self.rng.uniform(20, 60),
                    "genotype": tv.genotype if self.rng.random() < 0.9 else "0/1",
                    "support": self.rng.randint(3, 20),
                })

        # False positives
        n_fp = int(len(calls) * (1 - precision_rate) / precision_rate)
        chroms = [f"chr{i}" for i in range(1, 23)]
        for _ in range(n_fp):
            sv_type = self.rng.choice(["DEL", "INS", "DUP", "INV"])
            size = int(self.rng.lognormal(mean=6, sigma=1.5))
            start = self.rng.randint(1000000, 200000000)
            calls.append({
                "chrom": self.rng.choice(chroms),
                "start": start, "end": start + size,
                "sv_type": sv_type, "size": size,
                "quality": self.rng.uniform(10, 40),
                "genotype": "0/1",
                "support": self.rng.randint(2, 8),
            })

        return calls
