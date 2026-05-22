"""
Integrated Structural Variant Detection
========================================
Three complementary strategies combined via a Bayesian evidence framework:

  1. Split-Read (SR)   — detect chimeric alignments indicating SV breakpoints
  2. Read-Depth (RD)   — copy-number changes from coverage depth
  3. Assembly-Based (AB)— local de-novo assembly + comparison to reference

SV types handled: DEL, INS, INV, DUP, TRA, BND, CNV
"""

import json
import math
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple
from enum import Enum


# ---------------------------------------------------------------------------
# Data Models
# ---------------------------------------------------------------------------

class SVType(str, Enum):
    DEL = "DEL"   # Deletion
    INS = "INS"   # Insertion
    INV = "INV"   # Inversion
    DUP = "DUP"   # Tandem duplication
    TRA = "TRA"   # Translocation
    BND = "BND"   # Breakend
    CNV = "CNV"   # Copy-number variant


@dataclass
class Breakpoint:
    chrom: str
    pos: int
    strand: str   # "+" | "-" | "."
    confidence: float = 1.0


@dataclass
class SVCall:
    sv_id: str
    sv_type: SVType
    chrom: str
    start: int
    end: int
    length: int
    chrom2: Optional[str] = None   # for TRA/BND
    pos2: Optional[int] = None
    evidence_sr: float = 0.0       # split-read support score
    evidence_rd: float = 0.0       # read-depth support score
    evidence_ab: float = 0.0       # assembly support score
    combined_score: float = 0.0
    genotype: str = "./."          # VCF GT field
    read_support: int = 0
    strand_info: str = "+/+"
    filter_flags: List[str] = field(default_factory=list)
    annotations: Dict = field(default_factory=dict)

    def passes_quality(self, min_score: float = 0.5, min_reads: int = 5) -> bool:
        return self.combined_score >= min_score and self.read_support >= min_reads


# ---------------------------------------------------------------------------
# 1. Split-Read Detector
# ---------------------------------------------------------------------------

@dataclass
class AlignmentRecord:
    """Simplified PAF/SAM alignment record."""
    read_id: str
    read_len: int
    chrom: str
    ref_start: int
    ref_end: int
    strand: str
    mapq: int
    cigar: str  # simplified: e.g. "100M 50S"
    sa_tag: Optional[str] = None  # supplementary alignments (SA:Z: field)


def parse_sa_tag(sa_tag: str) -> List[AlignmentRecord]:
    """Parse SA:Z supplementary alignment tag from SAM format."""
    records = []
    for part in sa_tag.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = part.split(",")
        if len(fields) >= 5:
            records.append(AlignmentRecord(
                read_id="supp",
                read_len=0,
                chrom=fields[0],
                ref_start=int(fields[1]),
                ref_end=int(fields[1]) + 100,  # approximate
                strand=fields[2],
                mapq=int(fields[4]),
                cigar=fields[3],
            ))
    return records


def detect_split_read_svs(
    alignments: List[AlignmentRecord],
    min_mapq: int = 20,
    min_sv_size: int = 50,
    max_sv_size: int = 100_000_000,
) -> List[SVCall]:
    """
    Detect SVs from split/chimeric reads.

    For each read with a supplementary alignment (SA tag), infer the SV type
    from the relationship between primary and supplementary alignments:
      - Same chrom, same strand, gap → DEL or INS
      - Same chrom, opposite strand → INV
      - Different chrom → TRA
    """
    sv_calls: List[SVCall] = []
    sv_counter = 0

    for primary in alignments:
        if primary.mapq < min_mapq or not primary.sa_tag:
            continue

        supps = parse_sa_tag(primary.sa_tag)
        for supp in supps:
            if supp.mapq < min_mapq:
                continue

            same_chrom = primary.chrom == supp.chrom
            same_strand = primary.strand == supp.strand

            if same_chrom and same_strand:
                gap = supp.ref_start - primary.ref_end
                if gap > min_sv_size:
                    sv_type = SVType.DEL if gap > 0 else SVType.INS
                    sv_len = abs(gap)
                    if sv_len > max_sv_size:
                        continue
                    sv_counter += 1
                    sv_calls.append(SVCall(
                        sv_id=f"SR_{sv_counter:06d}",
                        sv_type=sv_type,
                        chrom=primary.chrom,
                        start=min(primary.ref_end, supp.ref_start),
                        end=max(primary.ref_end, supp.ref_start),
                        length=sv_len,
                        evidence_sr=1.0,
                        read_support=1,
                        strand_info=f"{primary.strand}/{supp.strand}",
                    ))
                elif gap < -min_sv_size:
                    # Read sequence inserted relative to reference
                    sv_counter += 1
                    sv_calls.append(SVCall(
                        sv_id=f"SR_{sv_counter:06d}",
                        sv_type=SVType.INS,
                        chrom=primary.chrom,
                        start=primary.ref_end,
                        end=primary.ref_end + 1,
                        length=abs(gap),
                        evidence_sr=1.0,
                        read_support=1,
                    ))
            elif same_chrom and not same_strand:
                # Inversion
                start = min(primary.ref_start, supp.ref_start)
                end = max(primary.ref_end, supp.ref_end)
                sv_len = end - start
                if min_sv_size <= sv_len <= max_sv_size:
                    sv_counter += 1
                    sv_calls.append(SVCall(
                        sv_id=f"SR_{sv_counter:06d}",
                        sv_type=SVType.INV,
                        chrom=primary.chrom,
                        start=start,
                        end=end,
                        length=sv_len,
                        evidence_sr=1.0,
                        read_support=1,
                        strand_info=f"{primary.strand}/{supp.strand}",
                    ))
            else:
                # Translocation
                sv_counter += 1
                sv_calls.append(SVCall(
                    sv_id=f"SR_{sv_counter:06d}",
                    sv_type=SVType.TRA,
                    chrom=primary.chrom,
                    start=primary.ref_end,
                    end=primary.ref_end + 1,
                    length=0,
                    chrom2=supp.chrom,
                    pos2=supp.ref_start,
                    evidence_sr=1.0,
                    read_support=1,
                    strand_info=f"{primary.strand}/{supp.strand}",
                ))

    return sv_calls


def cluster_split_read_svs(
    calls: List[SVCall], max_distance: int = 200
) -> List[SVCall]:
    """
    Cluster nearby split-read SV calls of the same type to produce merged events.
    Uses a simple sweep-line algorithm.
    """
    if not calls:
        return []

    # Group by (chrom, sv_type)
    groups: Dict[Tuple, List[SVCall]] = {}
    for c in calls:
        key = (c.chrom, c.sv_type.value)
        groups.setdefault(key, []).append(c)

    merged: List[SVCall] = []
    cluster_id = 0

    for (chrom, svtype), group in groups.items():
        group.sort(key=lambda x: x.start)
        current_cluster: List[SVCall] = [group[0]]

        for call in group[1:]:
            if call.start - current_cluster[-1].end <= max_distance:
                current_cluster.append(call)
            else:
                merged.append(_merge_cluster(current_cluster, cluster_id, chrom, SVType(svtype)))
                cluster_id += 1
                current_cluster = [call]

        merged.append(_merge_cluster(current_cluster, cluster_id, chrom, SVType(svtype)))
        cluster_id += 1

    return merged


def _merge_cluster(cluster: List[SVCall], cid: int, chrom: str, sv_type: SVType) -> SVCall:
    start = min(c.start for c in cluster)
    end = max(c.end for c in cluster)
    total_support = sum(c.read_support for c in cluster)
    mean_score = float(np.mean([c.evidence_sr for c in cluster]))
    return SVCall(
        sv_id=f"SRC_{cid:06d}",
        sv_type=sv_type,
        chrom=chrom,
        start=start,
        end=end,
        length=end - start,
        evidence_sr=mean_score,
        read_support=total_support,
    )


# ---------------------------------------------------------------------------
# 2. Read-Depth (Coverage) Detector
# ---------------------------------------------------------------------------

def compute_coverage_depth(
    alignments: List[AlignmentRecord],
    chrom: str,
    chrom_len: int,
    bin_size: int = 1000,
) -> np.ndarray:
    """Compute binned read-depth coverage vector for a chromosome."""
    n_bins = math.ceil(chrom_len / bin_size)
    depth = np.zeros(n_bins, dtype=np.float32)
    for aln in alignments:
        if aln.chrom != chrom:
            continue
        b_start = aln.ref_start // bin_size
        b_end = min(aln.ref_end // bin_size, n_bins - 1)
        depth[b_start : b_end + 1] += 1
    return depth


def detect_cnv_from_depth(
    depth: np.ndarray,
    bin_size: int = 1000,
    chrom: str = "chr1",
    min_sv_bins: int = 5,
    z_threshold: float = 3.0,
) -> List[SVCall]:
    """
    Detect copy-number variants using a Z-score approach on GC-corrected depth.

    Algorithm:
      1. Compute median and MAD of depth
      2. Z-score normalisation
      3. Identify contiguous bins with |Z| > threshold
      4. Classify: Z < -threshold → DEL, Z > threshold → DUP
    """
    if len(depth) == 0:
        return []

    median = float(np.median(depth))
    mad = float(np.median(np.abs(depth - median))) + 1e-6
    z_scores = (depth - median) / (1.4826 * mad)

    calls: List[SVCall] = []
    sv_id = 0

    i = 0
    while i < len(z_scores):
        if abs(z_scores[i]) > z_threshold:
            direction = 1 if z_scores[i] > 0 else -1
            j = i
            while j < len(z_scores) and z_scores[j] * direction > z_threshold * 0.5:
                j += 1
            n_bins = j - i
            if n_bins >= min_sv_bins:
                sv_type = SVType.DUP if direction > 0 else SVType.DEL
                start = i * bin_size
                end = j * bin_size
                mean_z = float(np.mean(np.abs(z_scores[i:j])))
                # Normalised evidence score: clamp z to [0,1]
                evidence = min(1.0, (mean_z - z_threshold) / z_threshold)
                calls.append(SVCall(
                    sv_id=f"RD_{sv_id:06d}",
                    sv_type=sv_type,
                    chrom=chrom,
                    start=start,
                    end=end,
                    length=end - start,
                    evidence_rd=evidence,
                    read_support=int(np.sum(depth[i:j])),
                ))
                sv_id += 1
            i = j
        else:
            i += 1

    return calls


# ---------------------------------------------------------------------------
# 3. Assembly-Based Detector
# ---------------------------------------------------------------------------

@dataclass
class LocalAssembly:
    """Simulated local de-novo assembly result around a candidate breakpoint."""
    region: Tuple[str, int, int]
    contig_sequence: str
    alignment_score: float
    breakpoint_positions: List[int]


def local_assemble_region(
    reads: List[str],
    reference: str,
    region: Tuple[str, int, int],
) -> LocalAssembly:
    """
    Reference implementation of local de-novo assembly.
    Production: use wtdbg2 / hifiasm / miniasm for actual assembly.

    This simplified version:
      1. Build a de-Bruijn graph with k=15
      2. Find the Eulerian path (greedy)
      3. Align contig back to reference to find breakpoints
    """
    k = 15
    # Build k-mer graph
    kmer_counts: Dict[str, int] = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i : i + k]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    # Greedy contig assembly: extend from most frequent k-mer
    if not kmer_counts:
        return LocalAssembly(region, "", 0.0, [])

    start_kmer = max(kmer_counts, key=kmer_counts.get)
    contig = start_kmer
    for _ in range(500):  # max extension steps
        last_kmer = contig[-(k - 1):]
        best_ext, best_count = None, 0
        for base in "ACGT":
            candidate = last_kmer + base
            count = kmer_counts.get(candidate, 0)
            if count > best_count:
                best_count = count
                best_ext = base
        if best_ext is None or best_count < 2:
            break
        contig += best_ext

    # Detect breakpoints by comparing contig to reference (edit-distance heuristic)
    breakpoints = _find_breakpoints_in_contig(contig, reference)
    aln_score = len(contig) / max(len(reference), 1)

    return LocalAssembly(region, contig, aln_score, breakpoints)


def _find_breakpoints_in_contig(contig: str, reference: str, window: int = 50) -> List[int]:
    """
    Identify positions in contig where alignment to reference breaks.
    Uses a sliding-window mismatch rate heuristic.
    """
    breakpoints = []
    min_len = min(len(contig), len(reference))
    if min_len < window:
        return []
    for i in range(0, min_len - window, window // 2):
        c_chunk = contig[i : i + window]
        r_chunk = reference[i : i + window]
        mismatches = sum(a != b for a, b in zip(c_chunk, r_chunk))
        if mismatches / window > 0.3:  # >30% mismatch → breakpoint zone
            breakpoints.append(i)
    return breakpoints


def detect_assembly_based_svs(
    assemblies: List[LocalAssembly],
    min_sv_size: int = 50,
) -> List[SVCall]:
    """Convert local assembly breakpoints to SV calls."""
    calls: List[SVCall] = []
    for i, asm in enumerate(assemblies):
        if not asm.breakpoint_positions:
            continue
        chrom, reg_start, _ = asm.region
        for bp in asm.breakpoint_positions:
            abs_pos = reg_start + bp
            calls.append(SVCall(
                sv_id=f"AB_{i:06d}_{bp}",
                sv_type=SVType.INS,  # placeholder; refine with BLAST/minimap2
                chrom=chrom,
                start=abs_pos,
                end=abs_pos + min_sv_size,
                length=min_sv_size,
                evidence_ab=min(1.0, asm.alignment_score),
                read_support=len(asm.contig_sequence) // 100,
            ))
    return calls


# ---------------------------------------------------------------------------
# Evidence Integration (Bayesian Fusion)
# ---------------------------------------------------------------------------

EVIDENCE_WEIGHTS = {
    "sr": 0.4,   # split-read
    "rd": 0.3,   # read-depth
    "ab": 0.3,   # assembly-based
}


def bayesian_evidence_integration(
    sr_calls: List[SVCall],
    rd_calls: List[SVCall],
    ab_calls: List[SVCall],
    reciprocal_overlap: float = 0.5,
    distance_tolerance: int = 300,
) -> List[SVCall]:
    """
    Merge evidence from three SV calling strategies using Bayesian combination.

    For each cluster of overlapping calls:
      combined_score = sigmoid(
          w_sr * log_odds(sr) + w_rd * log_odds(rd) + w_ab * log_odds(ab)
      )

    Calls supported by ≥2 strategies receive a quality boost.
    """
    all_calls = sr_calls + rd_calls + ab_calls
    if not all_calls:
        return []

    # Index by chrom
    by_chrom: Dict[str, List[SVCall]] = {}
    for c in all_calls:
        by_chrom.setdefault(c.chrom, []).append(c)

    integrated: List[SVCall] = []
    global_id = 0

    for chrom, calls in by_chrom.items():
        calls.sort(key=lambda x: x.start)
        clusters: List[List[SVCall]] = []
        current: List[SVCall] = [calls[0]]

        for call in calls[1:]:
            # Check reciprocal overlap or proximity
            last = current[-1]
            overlap = _reciprocal_overlap(last, call)
            dist = call.start - last.end
            if overlap >= reciprocal_overlap or dist <= distance_tolerance:
                current.append(call)
            else:
                clusters.append(current)
                current = [call]
        clusters.append(current)

        for cluster in clusters:
            merged = _integrate_cluster(cluster, global_id, chrom)
            integrated.append(merged)
            global_id += 1

    return integrated


def _reciprocal_overlap(a: SVCall, b: SVCall) -> float:
    if a.sv_type != b.sv_type:
        return 0.0
    overlap_start = max(a.start, b.start)
    overlap_end = min(a.end, b.end)
    if overlap_end <= overlap_start:
        return 0.0
    overlap_len = overlap_end - overlap_start
    return overlap_len / max(a.end - a.start, b.end - b.start, 1)


def _integrate_cluster(cluster: List[SVCall], gid: int, chrom: str) -> SVCall:
    start = int(np.median([c.start for c in cluster]))
    end = int(np.median([c.end for c in cluster]))
    sv_type = max(set(c.sv_type for c in cluster), key=lambda t: sum(1 for c in cluster if c.sv_type == t))

    max_sr = max(c.evidence_sr for c in cluster)
    max_rd = max(c.evidence_rd for c in cluster)
    max_ab = max(c.evidence_ab for c in cluster)

    # Bayesian combination via log-odds fusion
    def log_odds(p: float) -> float:
        p = max(min(p, 0.9999), 0.0001)
        return math.log(p / (1 - p))

    fused_log_odds = (
        EVIDENCE_WEIGHTS["sr"] * log_odds(max_sr if max_sr > 0 else 0.01) +
        EVIDENCE_WEIGHTS["rd"] * log_odds(max_rd if max_rd > 0 else 0.01) +
        EVIDENCE_WEIGHTS["ab"] * log_odds(max_ab if max_ab > 0 else 0.01)
    )

    def sigmoid(x: float) -> float:
        return 1.0 / (1.0 + math.exp(-x))

    combined = sigmoid(fused_log_odds)

    # Multi-method bonus
    n_strategies = sum([max_sr > 0, max_rd > 0, max_ab > 0])
    if n_strategies >= 2:
        combined = min(1.0, combined * 1.15)

    total_reads = sum(c.read_support for c in cluster)

    return SVCall(
        sv_id=f"SV_{gid:06d}",
        sv_type=sv_type,
        chrom=chrom,
        start=start,
        end=end,
        length=end - start,
        evidence_sr=max_sr,
        evidence_rd=max_rd,
        evidence_ab=max_ab,
        combined_score=combined,
        read_support=total_reads,
        genotype="0/1" if combined < 0.85 else "1/1",
    )


# ---------------------------------------------------------------------------
# VCF Output
# ---------------------------------------------------------------------------

def write_vcf(
    sv_calls: List[SVCall],
    output_path: str,
    sample_name: str = "SAMPLE",
    ref_genome: str = "hg38",
) -> None:
    """Write SV calls to VCF 4.2 format."""
    header = [
        "##fileformat=VCFv4.2",
        f"##reference={ref_genome}",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=INS,Description=\"Insertion\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=TRA,Description=\"Translocation\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
        "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"Combined evidence score\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Read support count\">",
        "##INFO=<ID=SR,Number=1,Type=Float,Description=\"Split-read evidence\">",
        "##INFO=<ID=RD,Number=1,Type=Float,Description=\"Read-depth evidence\">",
        "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Assembly evidence\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}",
    ]
    with open(output_path, "w") as fh:
        for line in header:
            fh.write(line + "\n")
        for sv in sorted(sv_calls, key=lambda x: (x.chrom, x.start)):
            qual = int(sv.combined_score * 60)
            filt = "PASS" if sv.passes_quality() else "LowQual"
            info = (
                f"SVTYPE={sv.sv_type.value};"
                f"SVLEN={sv.length};"
                f"END={sv.end};"
                f"SCORE={sv.combined_score:.4f};"
                f"SUPPORT={sv.read_support};"
                f"SR={sv.evidence_sr:.3f};"
                f"RD={sv.evidence_rd:.3f};"
                f"AB={sv.evidence_ab:.3f}"
            )
            gt = sv.genotype
            fh.write(
                f"{sv.chrom}\t{sv.start}\t{sv.sv_id}\tN\t<{sv.sv_type.value}>\t"
                f"{qual}\t{filt}\t{info}\tGT\t{gt}\n"
            )


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo():
    np.random.seed(42)

    # Simulate split-read alignments
    alignments = []
    for i in range(200):
        chrom = f"chr{np.random.choice([1,2,3,7,22])}"
        pos = np.random.randint(1_000_000, 50_000_000)
        strand = np.random.choice(["+", "-"])
        # 30% of reads have supplementary alignment (SV-supporting)
        sa = None
        if np.random.random() < 0.3:
            gap = int(np.random.choice([500, 1000, 5000, 50000]))
            sa = f"{chrom},{pos + gap},{strand},100M,60,0;"
        alignments.append(AlignmentRecord(
            read_id=f"read_{i}",
            read_len=10000,
            chrom=chrom,
            ref_start=pos,
            ref_end=pos + 8000,
            strand=strand,
            mapq=np.random.randint(20, 60),
            cigar="8000M 2000S",
            sa_tag=sa,
        ))

    sr_calls = detect_split_read_svs(alignments)
    sr_calls = cluster_split_read_svs(sr_calls)

    # Simulate read-depth data
    depth = np.random.poisson(30, 3000).astype(np.float32)
    # Inject deletions and duplications
    depth[500:550] *= 0.1   # deletion
    depth[1200:1260] *= 2.5  # duplication
    rd_calls = detect_cnv_from_depth(depth, chrom="chr7")

    # Simulate assembly-based calls
    mock_assemblies = [
        LocalAssembly(
            ("chr1", 5_000_000, 5_020_000),
            "A" * 500 + "TTTTTTTTTTTTTTTT" + "C" * 300,
            0.75,
            [400, 516],
        )
    ]
    ab_calls = detect_assembly_based_svs(mock_assemblies)

    # Integrate evidence
    integrated = bayesian_evidence_integration(sr_calls, rd_calls, ab_calls)

    print("=== Integrated SV Caller Demo ===")
    print(f"Split-read calls (after clustering): {len(sr_calls)}")
    print(f"Read-depth CNV calls               : {len(rd_calls)}")
    print(f"Assembly-based calls               : {len(ab_calls)}")
    print(f"Integrated SV calls                : {len(integrated)}")
    passing = [c for c in integrated if c.passes_quality()]
    print(f"PASS-filter calls                  : {len(passing)}")

    vcf_path = "/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/sv_calls_demo.vcf"
    write_vcf(integrated, vcf_path)
    print(f"\nVCF written: {vcf_path}")

    # Summary JSON
    type_counts = {}
    for c in integrated:
        type_counts[c.sv_type.value] = type_counts.get(c.sv_type.value, 0) + 1

    summary = {
        "total_sv_calls": len(integrated),
        "passing_calls": len(passing),
        "by_type": type_counts,
        "mean_combined_score": float(np.mean([c.combined_score for c in integrated])),
    }
    with open("/app/projects/bf9f3f3c-3ec6-4692-a347-6ef4a8b2cc12/workspace/results/sv_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("Saved: results/sv_summary.json")
    return summary


if __name__ == "__main__":
    demo()
