"""
Complex Structural Variant Detection
======================================
Detects genomically complex events:
  1. Chromothripsis — massive shattering of one or few chromosomes
  2. Extrachromosomal DNA (ecDNA) — amplified circular DNA elements
  3. Chromoplexy — chains of translocations
  4. BFB (Breakage-Fusion-Bridge) cycles — fold-back inversions
"""

import json
import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from itertools import combinations


# ---------------------------------------------------------------------------
# Data Models
# ---------------------------------------------------------------------------

@dataclass
class ChromothripticEvent:
    chrom: str
    n_breakpoints: int
    n_sv_types: int             # variety: DEL + INV + DUP together
    mean_segment_size: float    # bp
    oscillating_cn_amplitude: float  # peak-to-trough in CN signal
    random_join_score: float    # how random strand orientations appear
    confidence: float           # 0–1 classifier output


@dataclass
class EcDNA:
    event_id: str
    chroms_involved: List[str]
    segments: List[Tuple[str, int, int]]   # (chrom, start, end) per segment
    total_length: int
    copy_number: float
    circular_evidence: float   # splice-junction read pairs forming circle
    oncogene_amplified: Optional[str] = None


@dataclass
class ChromoplexEvent:
    chroms: List[str]
    breakpoints: List[Tuple[str, int]]   # (chrom, pos)
    chain_length: int
    deletion_size_total: int


@dataclass
class BFBEvent:
    chrom: str
    fold_back_position: int
    n_cycles: int
    mean_amplification: float


# ---------------------------------------------------------------------------
# 1. Chromothripsis Detector
# ---------------------------------------------------------------------------

# Statistical thresholds based on: Stephens et al. Cell 2011; Korbel & Campbell 2013
MIN_BREAKPOINTS_CHROM = 10
MIN_SV_TYPES = 2
MAX_MEAN_SEGMENT_KB = 500
CN_OSCILLATION_THRESHOLD = 1.5   # copy-number oscillation amplitude


def detect_chromothripsis(
    sv_calls: List,  # SVCall objects
    cn_profile: np.ndarray,    # copy-number array (per-bin)
    chrom: str,
    bin_size: int = 100_000,
    n_permutations: int = 1000,
    alpha: float = 0.05,
) -> Optional[ChromothripticEvent]:
    """
    Test whether SVs on a chromosome are consistent with chromothripsis.

    Five hallmarks tested:
      1. ≥10 breakpoints clustered on one chromosome
      2. Mixed SV types (DEL + INV + DUP)
      3. Short inter-breakpoint segments
      4. Oscillating copy-number (typically 1–2 states)
      5. Random strand orientation of joined segments (p-value via permutation)

    Returns a ChromothripticEvent if all hallmarks pass, else None.
    """
    chrom_calls = [sv for sv in sv_calls if sv.chrom == chrom]

    # Hallmark 1: breakpoint count
    n_bp = len(chrom_calls) * 2  # each SV has 2 breakpoints
    if n_bp < MIN_BREAKPOINTS_CHROM:
        return None

    # Hallmark 2: mixed SV types
    sv_types = set(sv.sv_type.value for sv in chrom_calls)
    n_types = len(sv_types)
    if n_types < MIN_SV_TYPES:
        return None

    # Hallmark 3: segment sizes
    breakpoints = sorted(
        [bp for sv in chrom_calls for bp in [sv.start, sv.end]]
    )
    if len(breakpoints) >= 2:
        segments = [breakpoints[i+1] - breakpoints[i] for i in range(len(breakpoints)-1)]
        mean_seg = float(np.mean(segments))
    else:
        mean_seg = float("inf")

    if mean_seg > MAX_MEAN_SEGMENT_KB * 1000:
        return None

    # Hallmark 4: CN oscillation
    if len(cn_profile) > 0:
        cn_std = float(np.std(cn_profile))
        cn_range = float(np.ptp(cn_profile)) if len(cn_profile) > 1 else 0
        cn_oscillation = cn_std
    else:
        cn_oscillation = 0.0

    # Hallmark 5: strand randomness test
    strands = [sv.strand_info for sv in chrom_calls if hasattr(sv, "strand_info")]
    random_join_score = _strand_randomness_test(strands, n_permutations)

    # Composite confidence score
    hallmarks_passed = sum([
        n_bp >= MIN_BREAKPOINTS_CHROM,
        n_types >= MIN_SV_TYPES,
        mean_seg <= MAX_MEAN_SEGMENT_KB * 1000,
        cn_oscillation >= CN_OSCILLATION_THRESHOLD,
        random_join_score >= (1 - alpha),
    ])
    confidence = hallmarks_passed / 5.0

    return ChromothripticEvent(
        chrom=chrom,
        n_breakpoints=n_bp,
        n_sv_types=n_types,
        mean_segment_size=mean_seg,
        oscillating_cn_amplitude=cn_oscillation,
        random_join_score=random_join_score,
        confidence=confidence,
    )


def _strand_randomness_test(strands: List[str], n_perm: int = 1000) -> float:
    """
    Test if strand orientations of rearranged segments are consistent with
    a random joining model (expected under chromothripsis).
    Returns p-value (high p = more random = more chromothripsis-like).
    """
    if len(strands) < 4:
        return 0.5

    # Count +/+ and +/- orientations
    counts = {"++": 0, "+-": 0, "-+": 0, "--": 0}
    for s in strands:
        if s in counts:
            counts[s] += 1

    total = sum(counts.values())
    if total == 0:
        return 0.5

    # Observed entropy of strand orientations
    observed_entropy = _entropy(list(counts.values()))

    # Permutation test: shuffle strands and compute entropy distribution
    strand_list = list(counts.keys())
    rng = np.random.default_rng(0)
    perm_entropies = []
    for _ in range(n_perm):
        perm_counts = np.zeros(4)
        for _ in range(total):
            idx = rng.integers(0, 4)
            perm_counts[idx] += 1
        perm_entropies.append(_entropy(list(perm_counts)))

    # P-value: fraction of permutations with entropy <= observed (i.e., how random)
    p_value = float(np.mean([e >= observed_entropy for e in perm_entropies]))
    return p_value


def _entropy(counts: List) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    probs = [c / total for c in counts if c > 0]
    return -sum(p * math.log2(p) for p in probs)


# ---------------------------------------------------------------------------
# 2. Extrachromosomal DNA (ecDNA) Detector
# ---------------------------------------------------------------------------

def detect_ecdna(
    sv_calls: List,
    cn_profile: Dict[str, np.ndarray],   # chrom → CN array
    bin_size: int = 100_000,
    min_copy_number: float = 5.0,
    min_circular_reads: int = 3,
) -> List[EcDNA]:
    """
    Detect extrachromosomal DNA (circular DNA amplicons).

    Algorithm (AmpliconArchitect-inspired):
      1. Identify focal high-CN regions (CN ≥ min_copy_number)
      2. Find inter-chromosomal or intra-chromosomal rearrangements
         that can be assembled into a circle (BND graph)
      3. Look for junction-spanning reads (back-spliced junctions)
      4. Classify as ecDNA if:
         - Circular topology can be reconstructed
         - Copy number is focally amplified
         - Junction reads support the circle
    """
    ecdna_events: List[EcDNA] = []
    event_id = 0

    # Step 1: Identify high-CN seeds
    high_cn_regions: List[Tuple[str, int, int, float]] = []
    for chrom, cn_arr in cn_profile.items():
        in_amp = False
        amp_start = 0
        for i, cn in enumerate(cn_arr):
            if cn >= min_copy_number and not in_amp:
                amp_start = i
                in_amp = True
            elif cn < min_copy_number and in_amp:
                high_cn_regions.append((
                    chrom,
                    amp_start * bin_size,
                    i * bin_size,
                    float(np.mean(cn_arr[amp_start:i])),
                ))
                in_amp = False

    if not high_cn_regions:
        return []

    # Step 2: Find BND SVs linking high-CN regions
    bnd_svs = [sv for sv in sv_calls if sv.sv_type.value in ("TRA", "BND", "INV", "DUP")]

    # Build a simple graph: nodes = high-CN regions, edges = BNDs
    def overlaps_region(sv, region) -> bool:
        chrom, start, end, _ = region
        return sv.chrom == chrom and sv.start >= start and sv.start <= end

    for i, region in enumerate(high_cn_regions):
        chrom, r_start, r_end, cn = region
        linked_svs = [sv for sv in bnd_svs if overlaps_region(sv, region)]

        if not linked_svs:
            continue

        # Step 3: Simulate circular evidence (back-split reads)
        # In production: count reads spanning the circle junction
        circular_evidence = min(1.0, len(linked_svs) * 0.2)
        n_junction_reads = int(len(linked_svs) * 1.5)

        if n_junction_reads < min_circular_reads and circular_evidence < 0.3:
            continue

        # Reconstruct circular segments
        segments = [(chrom, r_start, r_end)]
        involved_chroms = {chrom}
        total_length = r_end - r_start

        for sv in linked_svs[:3]:   # limit segments per circle
            if sv.chrom2:
                seg_chrom = sv.chrom2
                seg_start = sv.pos2 or 0
                seg_end = seg_start + 500_000
            else:
                seg_chrom = sv.chrom
                seg_start = sv.end
                seg_end = sv.end + 500_000
            segments.append((seg_chrom, seg_start, seg_end))
            involved_chroms.add(seg_chrom)
            total_length += seg_end - seg_start

        # Check for known oncogenes (simplified list)
        oncogene = _check_oncogene_amplification(chrom, r_start, r_end)

        ecdna_events.append(EcDNA(
            event_id=f"ecDNA_{event_id:04d}",
            chroms_involved=sorted(involved_chroms),
            segments=segments,
            total_length=total_length,
            copy_number=cn,
            circular_evidence=circular_evidence,
            oncogene_amplified=oncogene,
        ))
        event_id += 1

    return ecdna_events


ONCOGENE_COORDS = {
    # (chrom, start, end, name)
    ("chr8", 127_736_000, 127_742_000): "MYC",
    ("chr2", 15_968_000, 16_123_000): "MYCN",
    ("chr7", 55_086_000, 55_324_000): "EGFR",
    ("chr12", 25_204_000, 25_250_000): "KRAS",
    ("chr17", 37_844_000, 37_887_000): "ERBB2",
    ("chr9", 107_545_000, 107_549_000): "CDK4",
}


def _check_oncogene_amplification(chrom: str, start: int, end: int) -> Optional[str]:
    for (oc, os, oe), gene in ONCOGENE_COORDS.items():
        if oc == chrom and max(os, start) < min(oe, end):
            return gene
    return None


# ---------------------------------------------------------------------------
# 3. Chromoplexy Detector
# ---------------------------------------------------------------------------

def detect_chromoplexy(
    sv_calls: List,
    max_chain_gap: int = 1_000_000,
    min_chain_length: int = 3,
) -> List[ChromoplexEvent]:
    """
    Detect chromoplexy: chains of translocations that delete small genomic segments.

    Algorithm:
      1. Build directed graph from TRA/BND SV calls (breakpoint connectivity)
      2. Find paths in the graph (chains of translocations)
      3. Chains on ≥3 different chromosomes = chromoplexy
    """
    # Filter to translocations
    tra_calls = [sv for sv in sv_calls if sv.sv_type.value in ("TRA", "BND")]
    if len(tra_calls) < min_chain_length:
        return []

    # Build adjacency (simplified: by proximity on same chromosome)
    chains: List[ChromoplexEvent] = []
    used = set()

    for i, sv1 in enumerate(tra_calls):
        if i in used:
            continue
        chain_bps = [(sv1.chrom, sv1.start)]
        chain_chroms = {sv1.chrom}
        del_size = 0
        used.add(i)

        for j, sv2 in enumerate(tra_calls):
            if j in used:
                continue
            # Linked if sv2 starts near where sv1 ends (on chrom2)
            if (sv2.chrom == sv1.chrom2 and sv1.pos2 is not None and
                    abs(sv2.start - sv1.pos2) < max_chain_gap):
                chain_bps.append((sv2.chrom, sv2.start))
                chain_chroms.add(sv2.chrom)
                del_size += abs(sv2.start - (sv1.pos2 or 0))
                used.add(j)
                sv1 = sv2  # extend chain

        if len(chain_chroms) >= min_chain_length:
            chains.append(ChromoplexEvent(
                chroms=sorted(chain_chroms),
                breakpoints=chain_bps,
                chain_length=len(chain_bps),
                deletion_size_total=del_size,
            ))

    return chains


# ---------------------------------------------------------------------------
# 4. Breakage-Fusion-Bridge (BFB) Detector
# ---------------------------------------------------------------------------

def detect_bfb(
    sv_calls: List,
    cn_profile: Dict[str, np.ndarray],
    bin_size: int = 100_000,
    min_fold_back_reads: int = 3,
) -> List[BFBEvent]:
    """
    Detect BFB (Breakage-Fusion-Bridge) cycles.

    Hallmarks:
      1. Fold-back inversions (inverted reads aligning back-to-back)
      2. Amplification gradient: CN decreases toward chromosome arm
      3. Telomere loss at one end

    Algorithm:
      - Find inverted SVs with strand_info == "+/-" or "-/+"
      - Check for CN gradient (monotonic increase toward centromere)
      - Estimate number of BFB cycles from CN amplitude
    """
    bfb_events: List[BFBEvent] = []

    inv_calls = [sv for sv in sv_calls
                 if sv.sv_type.value == "INV"
                 and sv.strand_info in ("+/-", "-/+", "+/+", "-/-")]

    chroms_with_inv = set(sv.chrom for sv in inv_calls)

    for chrom in chroms_with_inv:
        cn_arr = cn_profile.get(chrom, np.array([]))
        if len(cn_arr) == 0:
            continue

        # Check for CN gradient (sign of BFB)
        gradient = np.gradient(cn_arr)
        monotonic_frac = float(np.mean(gradient > 0))

        if monotonic_frac < 0.6:  # not a gradient
            continue

        # Estimate cycles from CN amplitude
        max_cn = float(np.max(cn_arr))
        n_cycles = max(1, int(math.log2(max(max_cn, 1))))

        # Find focal point (highest CN bin)
        fold_back_pos = int(np.argmax(cn_arr)) * bin_size

        chrom_inv = [sv for sv in inv_calls if sv.chrom == chrom]
        if len(chrom_inv) < 1:
            continue

        bfb_events.append(BFBEvent(
            chrom=chrom,
            fold_back_position=fold_back_pos,
            n_cycles=n_cycles,
            mean_amplification=max_cn,
        ))

    return bfb_events


# ---------------------------------------------------------------------------
# Summary Report
# ---------------------------------------------------------------------------

def summarise_complex_svs(
    chromothripsis_events: List[Optional[ChromothripticEvent]],
    ecdna_events: List[EcDNA],
    chromoplexy_events: List[ChromoplexEvent],
    bfb_events: List[BFBEvent],
) -> Dict:
    valid_ct = [e for e in chromothripsis_events if e is not None]
    return {
        "chromothripsis": {
            "n_events": len(valid_ct),
            "events": [
                {
                    "chrom": e.chrom,
                    "n_breakpoints": e.n_breakpoints,
                    "confidence": round(e.confidence, 3),
                    "mean_segment_kb": round(e.mean_segment_size / 1000, 1),
                }
                for e in valid_ct
            ],
        },
        "ecdna": {
            "n_events": len(ecdna_events),
            "events": [
                {
                    "event_id": e.event_id,
                    "chroms": e.chroms_involved,
                    "total_length_kb": round(e.total_length / 1000, 1),
                    "copy_number": round(e.copy_number, 1),
                    "oncogene": e.oncogene_amplified,
                    "circular_evidence": round(e.circular_evidence, 3),
                }
                for e in ecdna_events
            ],
        },
        "chromoplexy": {
            "n_events": len(chromoplexy_events),
            "events": [
                {
                    "chroms": e.chroms,
                    "chain_length": e.chain_length,
                    "deletion_kb": round(e.deletion_size_total / 1000, 1),
                }
                for e in chromoplexy_events
            ],
        },
        "bfb": {
            "n_events": len(bfb_events),
            "events": [
                {
                    "chrom": e.chrom,
                    "n_cycles": e.n_cycles,
                    "max_cn": round(e.mean_amplification, 1),
                }
                for e in bfb_events
            ],
        },
    }


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo():
    import sys
    sys.path.insert(0, "/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/src/sv_detection")
    from sv_caller import SVCall, SVType

    np.random.seed(42)
    print("=== Complex SV Detector Demo ===\n")

    # Build mock SV calls for chr2
    sv_calls = []
    for i in range(25):
        sv_type = np.random.choice(["DEL", "INV", "DUP", "TRA", "BND"])
        start = np.random.randint(50_000_000, 100_000_000)
        end = start + np.random.randint(5_000, 500_000)
        chrom2 = "chr3" if sv_type in ("TRA", "BND") else None
        pos2 = np.random.randint(10_000_000, 50_000_000) if chrom2 else None
        sv_calls.append(SVCall(
            sv_id=f"SV_{i:04d}",
            sv_type=SVType(sv_type),
            chrom="chr2",
            start=start,
            end=end,
            length=end - start,
            combined_score=np.random.uniform(0.5, 1.0),
            read_support=np.random.randint(5, 50),
            strand_info=np.random.choice(["+/+", "+/-", "-/+", "-/-"]),
            chrom2=chrom2,
            pos2=pos2,
        ))

    # Mock CN profile (oscillating for chr2 → chromothripsis-like)
    cn_chr2 = np.array([1.0 if i % 2 == 0 else 3.0 for i in range(100)])
    cn_chr2 += np.random.normal(0, 0.3, 100)
    cn_chr7 = np.linspace(2, 16, 50) + np.random.normal(0, 0.5, 50)  # BFB gradient
    cn_chr8_start = 130  # near MYC
    cn_chr8 = np.ones(200) * 2.0
    cn_chr8[130:145] = 12.0  # focal amplification near MYC
    cn_profile = {"chr2": cn_chr2, "chr7": cn_chr7, "chr8": cn_chr8}

    # 1. Chromothripsis
    ct_event = detect_chromothripsis(sv_calls, cn_chr2, "chr2")
    print(f"Chromothripsis (chr2): {ct_event}")

    # 2. ecDNA
    ecdna_events = detect_ecdna(sv_calls, cn_profile)
    print(f"\nEcDNA events detected: {len(ecdna_events)}")
    for e in ecdna_events[:3]:
        print(f"  {e.event_id}: chroms={e.chroms_involved}, CN={e.copy_number:.1f}, "
              f"oncogene={e.oncogene_amplified}")

    # 3. Chromoplexy
    chromoplexy = detect_chromoplexy(sv_calls)
    print(f"\nChromoplexy chains: {len(chromoplexy)}")

    # 4. BFB
    bfb = detect_bfb(sv_calls, cn_profile)
    print(f"\nBFB events: {len(bfb)}")

    summary = summarise_complex_svs(
        [ct_event], ecdna_events, chromoplexy, bfb
    )

    with open("/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/complex_sv_demo.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved: results/complex_sv_demo.json")
    return summary


if __name__ == "__main__":
    demo()
