"""
GIAB Tier1 SV Benchmark Evaluation
=====================================
Evaluates SV caller performance against the GIAB (Genome in a Bottle)
HG002/NA24385 Tier1 SV truth set using truvari-style benchmarking.

Metrics:
  - Precision = TP / (TP + FP)
  - Recall    = TP / (TP + FN)
  - F1 score  = 2 × Precision × Recall / (Precision + Recall)
  - Genotype concordance

Reference:
  Zook et al. 2020, Nature Biotechnology: "A robust benchmark for germline
  structural variant detection"
"""

import json
import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Data Models
# ---------------------------------------------------------------------------

@dataclass
class TruthSVCall:
    """A ground-truth SV from GIAB Tier1 truth set."""
    sv_id: str
    sv_type: str
    chrom: str
    start: int
    end: int
    length: int
    gt: str = "0/1"
    tier: int = 1        # GIAB Tier1 = highest confidence


@dataclass
class BenchmarkResult:
    sv_type: str
    tp: int = 0
    fp: int = 0
    fn: int = 0
    gt_concordant: int = 0
    precision: float = 0.0
    recall: float = 0.0
    f1: float = 0.0
    gt_concordance: float = 0.0


# ---------------------------------------------------------------------------
# Truvari-style Matching
# ---------------------------------------------------------------------------

def match_sv_to_truth(
    query_sv: Dict,
    truth_calls: List[TruthSVCall],
    max_distance: int = 500,
    min_reciprocal_overlap: float = 0.5,
    min_size_similarity: float = 0.7,
    require_same_type: bool = True,
) -> Optional[TruthSVCall]:
    """
    Match a query SV call to a truth call (Truvari matching logic).

    Matching criteria:
      1. Same chromosome
      2. Same SV type (optional)
      3. |start_query - start_truth| ≤ max_distance
      4. Reciprocal overlap ≥ min_reciprocal_overlap (for DEL/DUP/INV)
      5. Size similarity = min(len_q, len_t) / max(len_q, len_t) ≥ min_size_similarity

    Returns the best matching truth call, or None if no match.
    """
    best_match: Optional[TruthSVCall] = None
    best_score = -1.0

    for truth in truth_calls:
        if truth.chrom != query_sv.get("chrom"):
            continue
        if require_same_type and truth.sv_type != query_sv.get("sv_type"):
            continue

        # Distance check
        dist = abs(query_sv["start"] - truth.start)
        if dist > max_distance:
            continue

        # Reciprocal overlap
        q_start, q_end = query_sv["start"], query_sv["end"]
        t_start, t_end = truth.start, truth.end
        overlap_start = max(q_start, t_start)
        overlap_end = min(q_end, t_end)
        if overlap_end > overlap_start:
            overlap_len = overlap_end - overlap_start
            q_len = max(q_end - q_start, 1)
            t_len = max(t_end - t_start, 1)
            recip_overlap = overlap_len / max(q_len, t_len)
        else:
            recip_overlap = 0.0

        if recip_overlap < min_reciprocal_overlap and query_sv.get("sv_type") != "INS":
            continue

        # Size similarity
        q_len = max(query_sv["end"] - query_sv["start"], 1)
        t_len = max(truth.length, 1)
        size_sim = min(q_len, t_len) / max(q_len, t_len)
        if size_sim < min_size_similarity:
            continue

        # Combined match score
        match_score = (
            0.4 * (1.0 - dist / max_distance) +
            0.4 * recip_overlap +
            0.2 * size_sim
        )

        if match_score > best_score:
            best_score = match_score
            best_match = truth

    return best_match


# ---------------------------------------------------------------------------
# Benchmark Evaluation Engine
# ---------------------------------------------------------------------------

def evaluate_against_giab(
    query_calls: List[Dict],
    truth_calls: List[TruthSVCall],
    min_size: int = 50,
    max_size: int = 1_000_000,
    sv_types: Optional[List[str]] = None,
    truvari_params: Optional[Dict] = None,
) -> Dict[str, BenchmarkResult]:
    """
    Main benchmark function. Returns per-SV-type BenchmarkResult.

    Parameters
    ----------
    query_calls  : list of dicts with keys: sv_id, sv_type, chrom, start, end, genotype
    truth_calls  : GIAB truth SVs
    min_size     : minimum SV size to evaluate
    max_size     : maximum SV size to evaluate
    sv_types     : list of SV types to evaluate; None = all
    truvari_params : override Truvari matching parameters

    Returns
    -------
    dict mapping sv_type → BenchmarkResult
    """
    params = {
        "max_distance": 500,
        "min_reciprocal_overlap": 0.5,
        "min_size_similarity": 0.7,
        "require_same_type": True,
    }
    if truvari_params:
        params.update(truvari_params)

    # Filter by size
    query_filtered = [
        c for c in query_calls
        if min_size <= (c["end"] - c["start"]) <= max_size
    ]
    truth_filtered = [
        t for t in truth_calls
        if min_size <= t.length <= max_size
    ]

    if sv_types is None:
        sv_types = list(set(
            [c["sv_type"] for c in query_filtered] +
            [t.sv_type for t in truth_filtered]
        ))

    results: Dict[str, BenchmarkResult] = {}

    for sv_type in sv_types:
        q_type = [c for c in query_filtered if c["sv_type"] == sv_type]
        t_type = [t for t in truth_filtered if t.sv_type == sv_type]

        tp_calls = []
        matched_truth = set()
        fp_calls = []

        for qc in q_type:
            truth_unmatched = [t for t in t_type if t.sv_id not in matched_truth]
            match = match_sv_to_truth(qc, truth_unmatched, **params)
            if match:
                tp_calls.append((qc, match))
                matched_truth.add(match.sv_id)
            else:
                fp_calls.append(qc)

        fn_calls = [t for t in t_type if t.sv_id not in matched_truth]

        tp = len(tp_calls)
        fp = len(fp_calls)
        fn = len(fn_calls)

        precision = tp / max(tp + fp, 1)
        recall = tp / max(tp + fn, 1)
        f1 = (2 * precision * recall) / max(precision + recall, 1e-9)

        # Genotype concordance
        gt_ok = sum(
            1 for qc, tc in tp_calls
            if qc.get("genotype", "./.")[:3] == tc.gt[:3]
        )
        gt_concordance = gt_ok / max(tp, 1)

        results[sv_type] = BenchmarkResult(
            sv_type=sv_type,
            tp=tp, fp=fp, fn=fn,
            gt_concordant=gt_ok,
            precision=precision,
            recall=recall,
            f1=f1,
            gt_concordance=gt_concordance,
        )

    return results


# ---------------------------------------------------------------------------
# Stratified Evaluation
# ---------------------------------------------------------------------------

def stratified_benchmark(
    query_calls: List[Dict],
    truth_calls: List[TruthSVCall],
    size_bins: Optional[List[Tuple[int, int]]] = None,
) -> Dict[str, Dict]:
    """
    Evaluate performance stratified by SV size bins.
    Default bins: 50–500, 500–5k, 5k–50k, 50k–1M bp.
    """
    if size_bins is None:
        size_bins = [
            (50, 500),
            (500, 5_000),
            (5_000, 50_000),
            (50_000, 1_000_000),
        ]

    stratified: Dict[str, Dict] = {}
    for min_s, max_s in size_bins:
        label = f"{min_s//1000 if min_s >= 1000 else min_s}bp-{max_s//1000}kb"
        results = evaluate_against_giab(
            query_calls, truth_calls, min_size=min_s, max_size=max_s
        )
        stratified[label] = {
            sv_type: {
                "precision": round(r.precision, 4),
                "recall": round(r.recall, 4),
                "f1": round(r.f1, 4),
                "tp": r.tp, "fp": r.fp, "fn": r.fn,
            }
            for sv_type, r in results.items()
        }
    return stratified


# ---------------------------------------------------------------------------
# GIAB HG002 Simulated Truth Set Generator
# ---------------------------------------------------------------------------

def generate_mock_giab_truth(
    n_del: int = 150,
    n_ins: int = 200,
    n_inv: int = 30,
    n_dup: int = 40,
    seed: int = 42,
) -> List[TruthSVCall]:
    """
    Generate a realistic mock GIAB Tier1 truth set for demonstration.
    Real GIAB: ~12,745 SVs in HG002 (Zook et al. 2020)
    """
    rng = np.random.default_rng(seed)
    calls: List[TruthSVCall] = []
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

    for i in range(n_del):
        chrom = rng.choice(chroms)
        pos = int(rng.integers(1_000_000, 100_000_000))
        length = int(rng.choice([
            rng.integers(50, 500),
            rng.integers(500, 5_000),
            rng.integers(5_000, 50_000),
        ]))
        calls.append(TruthSVCall(
            sv_id=f"GIAB_DEL_{i:06d}",
            sv_type="DEL",
            chrom=chrom,
            start=pos,
            end=pos + length,
            length=length,
            gt=rng.choice(["0/1", "1/1"], p=[0.7, 0.3]),
        ))

    for i in range(n_ins):
        chrom = rng.choice(chroms)
        pos = int(rng.integers(1_000_000, 100_000_000))
        length = int(rng.integers(50, 10_000))
        calls.append(TruthSVCall(
            sv_id=f"GIAB_INS_{i:06d}",
            sv_type="INS",
            chrom=chrom,
            start=pos,
            end=pos + 1,
            length=length,
            gt=rng.choice(["0/1", "1/1"], p=[0.75, 0.25]),
        ))

    for i in range(n_inv):
        chrom = rng.choice(chroms)
        pos = int(rng.integers(1_000_000, 100_000_000))
        length = int(rng.integers(1_000, 200_000))
        calls.append(TruthSVCall(
            sv_id=f"GIAB_INV_{i:06d}",
            sv_type="INV",
            chrom=chrom,
            start=pos,
            end=pos + length,
            length=length,
        ))

    for i in range(n_dup):
        chrom = rng.choice(chroms)
        pos = int(rng.integers(1_000_000, 100_000_000))
        length = int(rng.integers(1_000, 100_000))
        calls.append(TruthSVCall(
            sv_id=f"GIAB_DUP_{i:06d}",
            sv_type="DUP",
            chrom=chrom,
            start=pos,
            end=pos + length,
            length=length,
        ))

    return calls


def simulate_caller_output(
    truth_calls: List[TruthSVCall],
    sensitivity: float = 0.85,
    false_positive_rate: float = 0.08,
    position_noise_bp: int = 200,
    seed: int = 0,
) -> List[Dict]:
    """
    Simulate a realistic SV caller output based on a truth set.
    Models: sensitivity (recall), FP rate, and position noise.
    """
    rng = np.random.default_rng(seed)
    calls = []

    # True positives (subset of truth with noise)
    for truth in truth_calls:
        if rng.random() < sensitivity:
            noise = int(rng.integers(-position_noise_bp, position_noise_bp))
            length_noise = float(rng.uniform(0.9, 1.1))
            length = max(50, int(truth.length * length_noise))
            calls.append({
                "sv_id": f"CALL_{len(calls):06d}",
                "sv_type": truth.sv_type,
                "chrom": truth.chrom,
                "start": max(0, truth.start + noise),
                "end": max(0, truth.end + noise),
                "length": length,
                "combined_score": float(rng.uniform(0.6, 1.0)),
                "genotype": truth.gt if rng.random() < 0.9 else "./.",
                "read_support": int(rng.integers(5, 40)),
            })

    # False positives
    n_fp = int(len(truth_calls) * false_positive_rate)
    chroms = [f"chr{i}" for i in range(1, 23)]
    sv_types = ["DEL", "INS", "INV", "DUP"]
    for i in range(n_fp):
        chrom = rng.choice(chroms)
        pos = int(rng.integers(1_000_000, 100_000_000))
        sv_type = rng.choice(sv_types)
        length = int(rng.integers(50, 50_000))
        calls.append({
            "sv_id": f"FP_{i:06d}",
            "sv_type": sv_type,
            "chrom": chrom,
            "start": pos,
            "end": pos + length,
            "length": length,
            "combined_score": float(rng.uniform(0.3, 0.7)),
            "genotype": rng.choice(["0/1", "1/1"]),
            "read_support": int(rng.integers(2, 15)),
        })

    return calls


# ---------------------------------------------------------------------------
# Report Generation
# ---------------------------------------------------------------------------

def format_benchmark_report(
    results: Dict[str, BenchmarkResult],
    stratified: Dict[str, Dict],
    caller_name: str = "LongSV",
    truth_name: str = "GIAB_HG002_Tier1",
) -> str:
    lines = [
        f"# Benchmark Report: {caller_name} vs {truth_name}",
        "",
        "## Overall Performance by SV Type",
        "",
        f"| SV Type | Precision | Recall | F1 | TP | FP | FN | GT Concordance |",
        f"|---------|-----------|--------|----|----|----|-------|----------------|",
    ]
    for sv_type, r in sorted(results.items()):
        lines.append(
            f"| {sv_type} | {r.precision:.4f} | {r.recall:.4f} | {r.f1:.4f} | "
            f"{r.tp} | {r.fp} | {r.fn} | {r.gt_concordance:.4f} |"
        )

    lines += ["", "## Stratified by SV Size", ""]
    for size_bin, type_results in stratified.items():
        lines.append(f"### {size_bin}")
        lines.append("| SV Type | Precision | Recall | F1 |")
        lines.append("|---------|-----------|--------|----|")
        for sv_type, m in sorted(type_results.items()):
            lines.append(
                f"| {sv_type} | {m['precision']:.4f} | {m['recall']:.4f} | {m['f1']:.4f} |"
            )
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo():
    np.random.seed(42)
    print("=== GIAB Benchmark Demo ===\n")

    truth_calls = generate_mock_giab_truth(
        n_del=150, n_ins=200, n_inv=30, n_dup=40
    )
    print(f"Truth calls generated: {len(truth_calls)}")

    # Simulate three caller scenarios
    scenarios = {
        "LongSV_v1 (baseline)": simulate_caller_output(truth_calls, sensitivity=0.75, false_positive_rate=0.12),
        "LongSV_v2 (hybrid)": simulate_caller_output(truth_calls, sensitivity=0.88, false_positive_rate=0.06, seed=1),
        "LongSV_v3 (complex)": simulate_caller_output(truth_calls, sensitivity=0.91, false_positive_rate=0.05, seed=2),
    }

    all_benchmark = {}
    for scenario_name, calls in scenarios.items():
        results = evaluate_against_giab(calls, truth_calls)
        strat = stratified_benchmark(calls, truth_calls)
        all_benchmark[scenario_name] = {
            "n_calls": len(calls),
            "results_by_type": {
                sv_type: {
                    "precision": round(r.precision, 4),
                    "recall": round(r.recall, 4),
                    "f1": round(r.f1, 4),
                    "tp": r.tp, "fp": r.fp, "fn": r.fn,
                }
                for sv_type, r in results.items()
            },
        }
        # Print summary
        overall_f1 = np.mean([r.f1 for r in results.values()])
        overall_prec = np.mean([r.precision for r in results.values()])
        overall_rec = np.mean([r.recall for r in results.values()])
        print(f"\n{scenario_name}:")
        print(f"  Calls: {len(calls)}, Overall Precision: {overall_prec:.4f}, "
              f"Recall: {overall_rec:.4f}, F1: {overall_f1:.4f}")
        for sv_type, r in sorted(results.items()):
            print(f"  {sv_type}: P={r.precision:.4f} R={r.recall:.4f} F1={r.f1:.4f} "
                  f"TP={r.tp} FP={r.fp} FN={r.fn}")

        # Write per-scenario benchmark markdown
        report_text = format_benchmark_report(results, strat, caller_name=scenario_name)
        safe_name = scenario_name.split("(")[0].strip().replace(" ", "_").lower()
        report_path = f"/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/benchmark_{safe_name}.md"
        with open(report_path, "w") as f:
            f.write(report_text)

    with open("/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/benchmark_summary.json", "w") as f:
        json.dump(all_benchmark, f, indent=2)
    print("\nSaved: results/benchmark_summary.json")
    return all_benchmark


if __name__ == "__main__":
    demo()
