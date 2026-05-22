"""
Hybrid Short-Read + Long-Read SV Analysis
==========================================
Combines:
  - Long-read SV calls (ONT/PacBio) — high sensitivity for all SV types
  - Short-read SV calls (Illumina) — high precision for breakpoint resolution

Integration strategies:
  1. SURVIVOR-style merge (reciprocal overlap + distance)
  2. Short-read breakpoint refinement (±10 bp precision)
  3. Illumina genotyping at long-read-identified loci
  4. Discordant read-pair validation
"""

import json
import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


@dataclass
class HybridSVCall:
    sv_id: str
    sv_type: str
    chrom: str
    start: int
    end: int
    length: int
    # Evidence sources
    lr_calls: List[str] = field(default_factory=list)   # long-read SV IDs
    sr_calls: List[str] = field(default_factory=list)   # short-read SV IDs
    n_lr_support: int = 0
    n_sr_support: int = 0
    # Refined breakpoints from short reads
    start_refined: Optional[int] = None
    end_refined: Optional[int] = None
    breakpoint_precision_bp: int = 0
    # Quality metrics
    lr_score: float = 0.0
    sr_score: float = 0.0
    hybrid_score: float = 0.0
    genotype: str = "./."
    filter_flags: List[str] = field(default_factory=list)
    # Illumina validation
    illumina_dp: int = 0
    illumina_alt_reads: int = 0
    illumina_vaf: float = 0.0


# ---------------------------------------------------------------------------
# SURVIVOR-style SV Merge
# ---------------------------------------------------------------------------

def survivor_merge(
    lr_calls: List[Dict],     # long-read calls as dicts (sv_id, sv_type, chrom, start, end)
    sr_calls: List[Dict],     # short-read calls
    max_distance: int = 1000,
    min_overlap: float = 0.0,
    same_type_required: bool = True,
    same_strand_required: bool = False,
) -> List[HybridSVCall]:
    """
    Merge SV calls from multiple callers using SURVIVOR algorithm logic.

    For each pair of calls (lr, sr):
      - Same chromosome
      - Same SV type (optional)
      - |start_lr - start_sr| ≤ max_distance OR reciprocal overlap ≥ min_overlap
    → Merge into HybridSVCall with combined evidence.
    """
    merged: List[HybridSVCall] = []
    matched_sr = set()
    hybrid_id = 0

    for lr in lr_calls:
        best_sr = None
        best_dist = float("inf")

        for idx, sr in enumerate(sr_calls):
            if idx in matched_sr:
                continue
            if lr["chrom"] != sr["chrom"]:
                continue
            if same_type_required and lr["sv_type"] != sr["sv_type"]:
                continue

            dist = abs(lr["start"] - sr["start"])
            overlap = _reciprocal_overlap_dict(lr, sr)

            if dist <= max_distance or overlap >= min_overlap:
                if dist < best_dist:
                    best_dist = dist
                    best_sr = (idx, sr)

        if best_sr is not None:
            idx, sr = best_sr
            matched_sr.add(idx)

            # Breakpoint refinement: Illumina provides better resolution
            start_ref = sr["start"]  # Illumina-refined
            end_ref = sr["end"]
            bp_precision = min(
                abs(lr["start"] - sr["start"]),
                abs(lr["end"] - sr["end"]),
            )

            # Hybrid score: geometric mean weighted by caller quality
            lr_sc = lr.get("combined_score", 0.7)
            sr_sc = sr.get("score", 0.7)
            hybrid_sc = _hybrid_score(lr_sc, sr_sc, n_callers=2)

            merged.append(HybridSVCall(
                sv_id=f"HYB_{hybrid_id:06d}",
                sv_type=lr["sv_type"],
                chrom=lr["chrom"],
                start=lr["start"],
                end=lr["end"],
                length=lr["end"] - lr["start"],
                lr_calls=[lr["sv_id"]],
                sr_calls=[sr["sv_id"]],
                n_lr_support=lr.get("read_support", 1),
                n_sr_support=sr.get("read_support", 1),
                start_refined=start_ref,
                end_refined=end_ref,
                breakpoint_precision_bp=bp_precision,
                lr_score=lr_sc,
                sr_score=sr_sc,
                hybrid_score=hybrid_sc,
                genotype=lr.get("genotype", "0/1"),
                illumina_dp=sr.get("dp", 30),
                illumina_alt_reads=sr.get("alt_reads", 5),
                illumina_vaf=sr.get("vaf", 0.3),
            ))
            hybrid_id += 1
        else:
            # Long-read-only call
            lr_sc = lr.get("combined_score", 0.7)
            merged.append(HybridSVCall(
                sv_id=f"HYB_{hybrid_id:06d}",
                sv_type=lr["sv_type"],
                chrom=lr["chrom"],
                start=lr["start"],
                end=lr["end"],
                length=lr["end"] - lr["start"],
                lr_calls=[lr["sv_id"]],
                n_lr_support=lr.get("read_support", 1),
                lr_score=lr_sc,
                hybrid_score=lr_sc * 0.85,  # penalty for no SR support
                genotype=lr.get("genotype", "0/1"),
                filter_flags=["LR_ONLY"],
            ))
            hybrid_id += 1

    # Add SR-only calls not matched to any LR call
    for idx, sr in enumerate(sr_calls):
        if idx not in matched_sr:
            sr_sc = sr.get("score", 0.7)
            merged.append(HybridSVCall(
                sv_id=f"HYB_{hybrid_id:06d}",
                sv_type=sr["sv_type"],
                chrom=sr["chrom"],
                start=sr["start"],
                end=sr["end"],
                length=sr["end"] - sr["start"],
                sr_calls=[sr["sv_id"]],
                n_sr_support=sr.get("read_support", 1),
                sr_score=sr_sc,
                hybrid_score=sr_sc * 0.80,  # lower confidence for SR-only
                genotype=sr.get("genotype", "0/1"),
                filter_flags=["SR_ONLY"],
            ))
            hybrid_id += 1

    return merged


def _reciprocal_overlap_dict(a: Dict, b: Dict) -> float:
    if a.get("sv_type") != b.get("sv_type"):
        return 0.0
    overlap_start = max(a["start"], b["start"])
    overlap_end = min(a["end"], b["end"])
    if overlap_end <= overlap_start:
        return 0.0
    overlap_len = overlap_end - overlap_start
    return overlap_len / max(a["end"] - a["start"], b["end"] - b["start"], 1)


def _hybrid_score(lr_score: float, sr_score: float, n_callers: int) -> float:
    """
    Combine caller scores:
      - Geometric mean of individual caller scores
      - Multi-caller bonus: score * (1 + 0.1 * (n_callers - 1))
    """
    geo_mean = (lr_score * sr_score) ** 0.5
    bonus = 1 + 0.1 * (n_callers - 1)
    return min(1.0, geo_mean * bonus)


# ---------------------------------------------------------------------------
# Discordant Read-Pair Validation (Illumina)
# ---------------------------------------------------------------------------

@dataclass
class DiscordantReadPair:
    read1_chrom: str
    read1_pos: int
    read2_chrom: str
    read2_pos: int
    insert_size: int
    orientation: str  # "FF" | "RR" | "FR" | "RF"


def validate_with_discordant_pairs(
    sv_calls: List[HybridSVCall],
    discordant_pairs: List[DiscordantReadPair],
    flank: int = 500,
    min_pairs: int = 3,
) -> List[HybridSVCall]:
    """
    Validate SV calls using Illumina discordant read pairs.

    For each SV:
      - Count discordant pairs that span both breakpoints
      - Update sr_score and illumina_alt_reads accordingly
    """
    for sv in sv_calls:
        supporting_pairs = 0
        for dp in discordant_pairs:
            if dp.read1_chrom != sv.chrom:
                continue
            span1 = abs(dp.read1_pos - sv.start) <= flank
            span2 = abs(dp.read2_pos - sv.end) <= flank
            if span1 and span2:
                supporting_pairs += 1

        if supporting_pairs >= min_pairs:
            sv.illumina_alt_reads = max(sv.illumina_alt_reads, supporting_pairs)
            # Boost hybrid score with discordant pair evidence
            dp_evidence = min(1.0, supporting_pairs / 10.0)
            sv.hybrid_score = min(1.0, sv.hybrid_score + dp_evidence * 0.1)
        elif supporting_pairs == 0 and "LR_ONLY" in sv.filter_flags:
            sv.filter_flags.append("NO_SR_DISCORD")

    return sv_calls


# ---------------------------------------------------------------------------
# Illumina Genotyping (Re-genotyping at LR-called loci)
# ---------------------------------------------------------------------------

def regenotype_at_lr_loci(
    sv_calls: List[HybridSVCall],
    illumina_coverage: Dict[str, np.ndarray],  # chrom → per-base depth array
    bin_size: int = 1000,
    ploidy: int = 2,
) -> List[HybridSVCall]:
    """
    Re-genotype SV calls using Illumina coverage data.
    Uses a simple copy-number genotyper:
      CN = 2 × (within-SV depth) / (flanking depth)
      → CN ≈ 1 → HET deletion (0/1)
      → CN ≈ 0 → HOM deletion (1/1)
      → CN ≈ 3 → HET duplication (0/1 DUP)
    """
    for sv in sv_calls:
        cov = illumina_coverage.get(sv.chrom, np.array([]))
        if len(cov) == 0:
            continue

        sv_start_bin = sv.start // bin_size
        sv_end_bin = sv.end // bin_size
        flank_bins = max(5, (sv_end_bin - sv_start_bin) // 2)

        sv_bins = cov[sv_start_bin : sv_end_bin + 1]
        flank_left = cov[max(0, sv_start_bin - flank_bins) : sv_start_bin]
        flank_right = cov[sv_end_bin + 1 : sv_end_bin + 1 + flank_bins]

        if len(sv_bins) == 0:
            continue

        sv_depth = float(np.median(sv_bins)) if len(sv_bins) > 0 else 0
        flank_depth = float(np.median(np.concatenate([flank_left, flank_right]))) if (
            len(flank_left) + len(flank_right) > 0
        ) else sv_depth

        if flank_depth <= 0:
            continue

        relative_cn = (sv_depth / flank_depth) * ploidy

        if sv.sv_type == "DEL":
            if relative_cn < 0.3:
                sv.genotype = "1/1"
            elif relative_cn < 0.75:
                sv.genotype = "0/1"
            else:
                sv.filter_flags.append("GT_AMBIGUOUS")
        elif sv.sv_type == "DUP":
            if relative_cn > 2.5:
                sv.genotype = "1/1"
            elif relative_cn > 1.3:
                sv.genotype = "0/1"
            else:
                sv.filter_flags.append("GT_AMBIGUOUS")

        sv.illumina_dp = int(flank_depth)
        sv.illumina_vaf = sv_depth / max(flank_depth, 1)

    return sv_calls


# ---------------------------------------------------------------------------
# Quality Filters
# ---------------------------------------------------------------------------

def apply_hybrid_filters(
    sv_calls: List[HybridSVCall],
    min_hybrid_score: float = 0.5,
    min_lr_reads: int = 3,
    max_lr_only_score: float = 0.9,
) -> Tuple[List[HybridSVCall], List[HybridSVCall]]:
    """Return (PASS, FILTERED) call lists."""
    passing, filtered = [], []
    for sv in sv_calls:
        fail_reasons = []
        if sv.hybrid_score < min_hybrid_score:
            fail_reasons.append(f"LowScore:{sv.hybrid_score:.3f}")
        if sv.n_lr_support < min_lr_reads and "SR_ONLY" not in sv.filter_flags:
            fail_reasons.append(f"LowLRReads:{sv.n_lr_support}")
        if "LR_ONLY" in sv.filter_flags and sv.hybrid_score > max_lr_only_score:
            pass  # high-confidence LR-only → keep
        if fail_reasons:
            sv.filter_flags.extend(fail_reasons)
            filtered.append(sv)
        else:
            passing.append(sv)
    return passing, filtered


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo():
    np.random.seed(42)
    print("=== Hybrid SV Caller Demo ===\n")

    # Generate mock long-read SV calls
    lr_calls = []
    for i in range(40):
        chrom = f"chr{np.random.choice([1,2,7,8,17])}"
        start = int(np.random.randint(1_000_000, 100_000_000))
        sv_type = np.random.choice(["DEL", "INS", "INV", "DUP"])
        length = int(np.random.choice([500, 2000, 10000, 50000]))
        lr_calls.append({
            "sv_id": f"LR_{i:04d}",
            "sv_type": sv_type,
            "chrom": chrom,
            "start": start,
            "end": start + length,
            "combined_score": float(np.random.uniform(0.5, 1.0)),
            "read_support": int(np.random.randint(3, 50)),
            "genotype": np.random.choice(["0/1", "1/1"]),
        })

    # Generate Illumina calls — 70% overlap with LR calls
    sr_calls = []
    for i, lr in enumerate(lr_calls[:28]):
        noise = int(np.random.randint(-300, 300))
        sr_calls.append({
            "sv_id": f"SR_{i:04d}",
            "sv_type": lr["sv_type"],
            "chrom": lr["chrom"],
            "start": lr["start"] + noise,
            "end": lr["end"] + noise,
            "score": float(np.random.uniform(0.6, 1.0)),
            "read_support": int(np.random.randint(5, 30)),
            "genotype": lr["genotype"],
            "dp": int(np.random.randint(25, 60)),
            "alt_reads": int(np.random.randint(5, 20)),
            "vaf": float(np.random.uniform(0.3, 0.7)),
        })

    # Merge
    hybrid_calls = survivor_merge(lr_calls, sr_calls)

    # Mock Illumina coverage
    illumina_cov = {}
    for chrom in set(c["chrom"] for c in lr_calls):
        illumina_cov[chrom] = np.random.poisson(35, 200_000).astype(np.float32)

    # Re-genotype
    hybrid_calls = regenotype_at_lr_loci(hybrid_calls, illumina_cov)

    # Filter
    passing, filtered = apply_hybrid_filters(hybrid_calls)

    print(f"LR-only calls       : {len(lr_calls)}")
    print(f"SR-only calls       : {len(sr_calls)}")
    print(f"Merged hybrid calls : {len(hybrid_calls)}")
    print(f"  Supported by both : {sum(1 for c in hybrid_calls if c.lr_calls and c.sr_calls)}")
    print(f"  LR-only           : {sum(1 for c in hybrid_calls if 'LR_ONLY' in c.filter_flags)}")
    print(f"  SR-only           : {sum(1 for c in hybrid_calls if 'SR_ONLY' in c.filter_flags)}")
    print(f"PASS after filtering: {len(passing)}")
    print(f"Filtered out        : {len(filtered)}")

    mean_hybrid_score = float(np.mean([c.hybrid_score for c in passing]))
    print(f"Mean hybrid score (PASS): {mean_hybrid_score:.4f}")

    summary = {
        "n_lr_calls": len(lr_calls),
        "n_sr_calls": len(sr_calls),
        "n_hybrid_merged": len(hybrid_calls),
        "n_both_support": sum(1 for c in hybrid_calls if c.lr_calls and c.sr_calls),
        "n_pass": len(passing),
        "n_filtered": len(filtered),
        "mean_hybrid_score_pass": mean_hybrid_score,
    }
    with open("/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/hybrid_sv_demo.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved: results/hybrid_sv_demo.json")
    return summary


if __name__ == "__main__":
    demo()
