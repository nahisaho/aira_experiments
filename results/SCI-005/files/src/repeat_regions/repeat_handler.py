"""
Repeat Region Special Processing
==================================
Handles SVs in telomeres, centromeres, segmental duplications, and STRs.

Key challenges:
  - Multi-mapping reads → false-positive SV calls
  - Collapsed repeats in reference genome
  - Telomeric TTAGGG arrays require specialised detection
"""

import re
import math
import json
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Repeat Region Annotations
# ---------------------------------------------------------------------------

@dataclass
class RepeatRegion:
    chrom: str
    start: int
    end: int
    region_type: str   # "telomere" | "centromere" | "segdup" | "str" | "satellite"
    repeat_unit: str   # e.g. "TTAGGG" for telomere
    copy_number: float = 0.0
    gc_content: float = 0.0
    annotation: str = ""


# Known UCSC-style centromere coordinates (hg38, approximate)
HG38_CENTROMERES: Dict[str, Tuple[int, int]] = {
    "chr1": (121_500_000, 128_000_000),
    "chr2": (91_800_000, 96_000_000),
    "chr3": (87_500_000, 93_900_000),
    "chr4": (48_200_000, 52_700_000),
    "chr5": (45_800_000, 50_700_000),
    "chr6": (58_500_000, 62_600_000),
    "chr7": (58_900_000, 62_100_000),
    "chr8": (43_100_000, 48_100_000),
    "chr9": (42_200_000, 65_500_000),
    "chr10": (38_000_000, 42_300_000),
    "chr11": (51_000_000, 55_800_000),
    "chr12": (33_800_000, 38_200_000),
    "chrX": (60_600_000, 65_000_000),
    "chrY": (10_200_000, 12_500_000),
}

TELOMERE_MOTIF = "TTAGGG"
TELOMERE_MIN_COPIES = 5

# Segmental duplication identity threshold
SEGDUP_IDENTITY_THRESHOLD = 0.90


# ---------------------------------------------------------------------------
# Telomere Detection
# ---------------------------------------------------------------------------

def detect_telomeric_content(sequence: str, motif: str = TELOMERE_MOTIF) -> Dict:
    """
    Detect telomeric repeats in a read sequence.
    Returns fraction of sequence covered by TTAGGG / CCCTAA repeats.
    """
    fwd_motif = motif                           # TTAGGG
    rev_motif = _reverse_complement(motif)      # CCCTAA

    fwd_count = len(re.findall(fwd_motif, sequence))
    rev_count = len(re.findall(rev_motif, sequence))
    total_copies = fwd_count + rev_count
    covered_bases = total_copies * len(motif)
    telomere_fraction = covered_bases / max(len(sequence), 1)

    is_telomeric = (
        total_copies >= TELOMERE_MIN_COPIES or telomere_fraction > 0.5
    )

    return {
        "is_telomeric": is_telomeric,
        "fwd_copies": fwd_count,
        "rev_copies": rev_count,
        "total_copies": total_copies,
        "telomere_fraction": telomere_fraction,
    }


def estimate_telomere_length(
    reads: List[str],
    motif: str = TELOMERE_MOTIF,
) -> Dict:
    """
    Estimate telomere length distribution from a set of reads that
    map to telomeric regions.

    Method: TelomereHunter / Computel approach
      - Count canonical repeat copies per read
      - Distribution statistics → T-band estimate
    """
    lengths = []
    for read in reads:
        result = detect_telomeric_content(read, motif)
        if result["is_telomeric"]:
            estimated_len = result["total_copies"] * len(motif)
            lengths.append(estimated_len)

    if not lengths:
        return {"mean_bp": 0, "median_bp": 0, "std_bp": 0, "n_reads": 0}

    arr = np.array(lengths)
    return {
        "mean_bp": float(np.mean(arr)),
        "median_bp": float(np.median(arr)),
        "std_bp": float(np.std(arr)),
        "n_reads": len(lengths),
        "min_bp": float(np.min(arr)),
        "max_bp": float(np.max(arr)),
    }


# ---------------------------------------------------------------------------
# Centromere Handling
# ---------------------------------------------------------------------------

def is_centromeric(chrom: str, pos: int, flank: int = 500_000) -> bool:
    """Check whether a position is in or near a centromere."""
    if chrom not in HG38_CENTROMERES:
        return False
    cen_start, cen_end = HG38_CENTROMERES[chrom]
    return (cen_start - flank) <= pos <= (cen_end + flank)


def filter_centromeric_svs(
    sv_calls: List,
    min_evidence_in_centromere: float = 0.8,
) -> Tuple[List, List]:
    """
    Apply stricter evidence thresholds to SVs in centromeric regions.
    Returns (high_confidence_calls, flagged_centromeric_calls).
    """
    high_conf = []
    centromeric = []
    for sv in sv_calls:
        in_cen = is_centromeric(sv.chrom, sv.start) or is_centromeric(sv.chrom, sv.end)
        if in_cen:
            sv.filter_flags.append("CENTROMERE")
            sv.annotations["in_centromere"] = True
            if sv.combined_score >= min_evidence_in_centromere:
                high_conf.append(sv)
            else:
                centromeric.append(sv)
        else:
            high_conf.append(sv)
    return high_conf, centromeric


# ---------------------------------------------------------------------------
# Segmental Duplication Processing
# ---------------------------------------------------------------------------

@dataclass
class SegDupEntry:
    chrom: str
    start: int
    end: int
    chrom2: str
    start2: int
    end2: int
    identity: float


def resolve_segdup_multimap(
    sv_calls: List,
    segdup_db: List[SegDupEntry],
    identity_threshold: float = SEGDUP_IDENTITY_THRESHOLD,
) -> List:
    """
    Flag SV calls in segmental duplication regions that may be false positives
    due to multi-mapping reads.

    Strategy:
      - If SV breakpoint overlaps a segdup entry with identity > threshold,
        require assembly-based evidence (evidence_ab > 0) to pass.
    """
    for sv in sv_calls:
        for segdup in segdup_db:
            if segdup.identity < identity_threshold:
                continue
            if (segdup.chrom == sv.chrom and
                    segdup.start <= sv.start <= segdup.end):
                sv.filter_flags.append("SEGDUP")
                sv.annotations["segdup_identity"] = segdup.identity
                if sv.evidence_ab == 0.0:
                    sv.filter_flags.append("SEGDUP_NO_ASM")
    return sv_calls


# ---------------------------------------------------------------------------
# Short Tandem Repeat (STR) Expansion Detection
# ---------------------------------------------------------------------------

def detect_str_expansion(
    reads: List[str],
    motif: str,
    reference_copies: int,
    min_expansion_copies: int = 5,
) -> Dict:
    """
    Detect STR expansions (e.g. repeat expansions causing disease).
    Counts motif copies in each read and compares to reference.

    Applications: HTT (CAG), FMR1 (CGG), C9orf72 (GGGGCC).
    """
    copy_counts = []
    for read in reads:
        count = len(re.findall(motif, read))
        copy_counts.append(count)

    if not copy_counts:
        return {"expanded": False, "max_copies": 0, "reference_copies": reference_copies}

    max_copies = max(copy_counts)
    median_copies = float(np.median(copy_counts))
    expanded = max_copies >= reference_copies + min_expansion_copies

    return {
        "expanded": expanded,
        "motif": motif,
        "reference_copies": reference_copies,
        "max_copies": max_copies,
        "median_copies": median_copies,
        "n_reads_expanded": sum(1 for c in copy_counts if c >= reference_copies + min_expansion_copies),
        "total_reads": len(copy_counts),
    }


# ---------------------------------------------------------------------------
# Alpha Satellite / Satellite DNA
# ---------------------------------------------------------------------------

def compute_repeat_complexity(sequence: str, k: int = 4) -> float:
    """
    Linguistic complexity (Lempel-Ziv-like) of a sequence.
    Low values indicate repetitive / low-complexity regions.
    Range: 0 (maximally repetitive) → 1 (maximally complex).
    """
    n = len(sequence)
    if n <= k:
        return 1.0
    kmer_set = set()
    for i in range(n - k + 1):
        kmer_set.add(sequence[i : i + k])
    max_possible = min(4**k, n - k + 1)
    return len(kmer_set) / max_possible


def classify_repeat_region(sequence: str) -> str:
    """
    Classify a sequence into repeat category.
    Returns one of: "telomere", "alpha_satellite", "satellite", "str", "segdup", "normal"
    """
    tel = detect_telomeric_content(sequence)
    if tel["is_telomeric"]:
        return "telomere"

    complexity = compute_repeat_complexity(sequence)
    if complexity < 0.05:
        return "alpha_satellite"
    elif complexity < 0.15:
        return "satellite"
    elif complexity < 0.3:
        return "str"
    elif complexity < 0.5:
        return "segdup"
    return "normal"


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo():
    np.random.seed(42)
    print("=== Repeat Region Handler Demo ===\n")

    # 1. Telomere detection
    tel_read = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG" + "ACGTACGT" * 10
    tel_result = detect_telomeric_content(tel_read)
    print(f"Telomere detection: {tel_result}")

    # 2. Telomere length estimation
    reads = [
        (TELOMERE_MOTIF * (10 + i)) + "ACGT" * 20
        for i in range(20)
    ]
    tl_stats = estimate_telomere_length(reads)
    print(f"\nTelomere length estimate: mean={tl_stats['mean_bp']:.0f} bp, "
          f"median={tl_stats['median_bp']:.0f} bp")

    # 3. Centromere check
    print(f"\nIs chr7:60000000 centromeric? {is_centromeric('chr7', 60_000_000)}")
    print(f"Is chr7:10000000 centromeric? {is_centromeric('chr7', 10_000_000)}")

    # 4. STR expansion detection
    str_reads = [("CAG" * (20 + i)) + "ACGT" * 5 for i in range(15)]
    str_result = detect_str_expansion(str_reads, "CAG", reference_copies=20)
    print(f"\nSTR expansion (CAG/HTT): {str_result}")

    # 5. Complexity classification
    sequences = {
        "telomere": TELOMERE_MOTIF * 20,
        "alpha_sat": "TTCGTTGGAAACGGGA" * 20,
        "normal": "ACGTGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGCGATCGATCGATCGATCG",
        "str": "CACACACACACACACACACACACACACACACACACACACACA",
    }
    print("\nSequence complexity classification:")
    for name, seq in sequences.items():
        cplx = compute_repeat_complexity(seq)
        cls = classify_repeat_region(seq)
        print(f"  {name:12s}: complexity={cplx:.3f}, class={cls}")

    result = {
        "telomere_detection": tel_result,
        "telomere_length_stats": tl_stats,
        "str_expansion": str_result,
    }
    with open("/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/repeat_handler_demo.json", "w") as f:
        json.dump(result, f, indent=2)
    print("\nSaved: results/repeat_handler_demo.json")
    return result


if __name__ == "__main__":
    demo()
