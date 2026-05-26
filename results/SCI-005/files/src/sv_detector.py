#!/usr/bin/env python3
"""
Integrated structural variant detection module.
Combines Split-read, Read-depth, and Assembly-based strategies.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Set
from enum import Enum
import logging
from collections import defaultdict

logger = logging.getLogger(__name__)


class SVType(Enum):
    DELETION = "DEL"
    INSERTION = "INS"
    DUPLICATION = "DUP"
    INVERSION = "INV"
    TRANSLOCATION = "TRA"
    COMPLEX = "CPX"
    CHROMOTHRIPSIS = "CTH"
    ECDNA = "ecDNA"


@dataclass
class AlignmentSegment:
    """Represents a segment of a split-read alignment."""
    chrom: str
    start: int
    end: int
    strand: str
    mapq: int
    cigar: str
    query_start: int
    query_end: int


@dataclass
class SVCall:
    """Structural variant call with evidence and confidence."""
    sv_type: SVType
    chrom: str
    start: int
    end: int
    size: int
    quality: float
    genotype: str
    support_reads: int
    evidence_types: List[str] = field(default_factory=list)
    chrom2: Optional[str] = None
    info: Dict = field(default_factory=dict)

    def to_vcf_record(self) -> str:
        alt = f"<{self.sv_type.value}>"
        info_str = (
            f"SVTYPE={self.sv_type.value};SVLEN={self.size};"
            f"END={self.end};SUPPORT={self.support_reads};"
            f"EVIDENCE={','.join(self.evidence_types)}"
        )
        return (
            f"{self.chrom}\t{self.start}\t.\tN\t{alt}\t"
            f"{self.quality:.1f}\tPASS\t{info_str}\tGT\t{self.genotype}"
        )


@dataclass
class ReadDepthSignal:
    """Read depth signal for a genomic region."""
    chrom: str
    positions: np.ndarray
    depths: np.ndarray
    gc_content: np.ndarray


class SplitReadAnalyzer:
    """Detect SVs from split-read alignments."""

    def __init__(self, min_mapq: int = 20, min_sv_size: int = 50):
        self.min_mapq = min_mapq
        self.min_sv_size = min_sv_size

    def detect_from_splits(
        self, segments: List[AlignmentSegment]
    ) -> List[SVCall]:
        calls = []
        segments = [s for s in segments if s.mapq >= self.min_mapq]
        segments.sort(key=lambda s: (s.chrom, s.start))

        for i in range(len(segments) - 1):
            s1, s2 = segments[i], segments[i + 1]
            sv = self._classify_split_pair(s1, s2)
            if sv:
                calls.append(sv)

        return calls

    def _classify_split_pair(
        self, s1: AlignmentSegment, s2: AlignmentSegment
    ) -> Optional[SVCall]:
        if s1.chrom != s2.chrom:
            return SVCall(
                sv_type=SVType.TRANSLOCATION,
                chrom=s1.chrom, start=s1.end, end=s2.start,
                size=0, quality=min(s1.mapq, s2.mapq),
                genotype="0/1", support_reads=1,
                evidence_types=["split_read"], chrom2=s2.chrom
            )

        gap = s2.start - s1.end
        query_gap = s2.query_start - s1.query_end

        if gap > self.min_sv_size and abs(query_gap) < gap * 0.5:
            return SVCall(
                sv_type=SVType.DELETION, chrom=s1.chrom,
                start=s1.end, end=s2.start, size=gap,
                quality=min(s1.mapq, s2.mapq),
                genotype="0/1", support_reads=1,
                evidence_types=["split_read"]
            )

        if query_gap > self.min_sv_size and gap < query_gap * 0.5:
            return SVCall(
                sv_type=SVType.INSERTION, chrom=s1.chrom,
                start=s1.end, end=s1.end + 1, size=query_gap,
                quality=min(s1.mapq, s2.mapq),
                genotype="0/1", support_reads=1,
                evidence_types=["split_read"]
            )

        if s1.strand != s2.strand:
            return SVCall(
                sv_type=SVType.INVERSION, chrom=s1.chrom,
                start=min(s1.start, s2.start),
                end=max(s1.end, s2.end),
                size=max(s1.end, s2.end) - min(s1.start, s2.start),
                quality=min(s1.mapq, s2.mapq),
                genotype="0/1", support_reads=1,
                evidence_types=["split_read"]
            )

        return None


class ReadDepthAnalyzer:
    """Detect SVs from read depth changes."""

    def __init__(self, window_size: int = 1000, min_ratio: float = 0.5):
        self.window_size = window_size
        self.min_ratio = min_ratio

    def compute_depth_profile(
        self, alignments: List[Dict], chrom: str, chrom_length: int
    ) -> ReadDepthSignal:
        n_bins = chrom_length // self.window_size + 1
        depths = np.zeros(n_bins)
        positions = np.arange(n_bins) * self.window_size

        for aln in alignments:
            if aln.get("chrom") == chrom:
                start_bin = aln["start"] // self.window_size
                end_bin = min(aln["end"] // self.window_size, n_bins - 1)
                depths[start_bin:end_bin + 1] += 1

        gc_content = np.random.uniform(0.3, 0.7, n_bins)
        return ReadDepthSignal(chrom, positions, depths, gc_content)

    def gc_correct(self, signal: ReadDepthSignal) -> np.ndarray:
        """GC-content correction using LOESS-like approach."""
        corrected = signal.depths.copy()
        gc_bins = np.digitize(signal.gc_content, np.linspace(0, 1, 20))
        for b in range(1, 21):
            mask = gc_bins == b
            if np.sum(mask) > 0:
                median_depth = np.median(signal.depths[mask])
                if median_depth > 0:
                    corrected[mask] = signal.depths[mask] / median_depth
        return corrected

    def detect_cnv(self, signal: ReadDepthSignal) -> List[SVCall]:
        corrected = self.gc_correct(signal)
        median_depth = np.median(corrected)
        if median_depth == 0:
            return []

        calls = []
        ratio = corrected / median_depth
        segments = self._segment_cbs(ratio)

        for seg_start, seg_end, seg_mean in segments:
            if seg_mean < self.min_ratio:
                calls.append(SVCall(
                    sv_type=SVType.DELETION, chrom=signal.chrom,
                    start=int(signal.positions[seg_start]),
                    end=int(signal.positions[min(seg_end, len(signal.positions) - 1)]),
                    size=int((seg_end - seg_start) * self.window_size),
                    quality=30.0, genotype="0/1",
                    support_reads=0, evidence_types=["read_depth"]
                ))
            elif seg_mean > 1.5:
                calls.append(SVCall(
                    sv_type=SVType.DUPLICATION, chrom=signal.chrom,
                    start=int(signal.positions[seg_start]),
                    end=int(signal.positions[min(seg_end, len(signal.positions) - 1)]),
                    size=int((seg_end - seg_start) * self.window_size),
                    quality=30.0, genotype="0/1",
                    support_reads=0, evidence_types=["read_depth"]
                ))

        return calls

    def _segment_cbs(
        self, data: np.ndarray, min_segment: int = 5
    ) -> List[Tuple[int, int, float]]:
        """Circular Binary Segmentation (simplified)."""
        segments = []
        self._cbs_recursive(data, 0, len(data), segments, min_segment)
        return segments

    def _cbs_recursive(
        self, data: np.ndarray, start: int, end: int,
        segments: List, min_segment: int
    ):
        if end - start < min_segment:
            segments.append((start, end, float(np.mean(data[start:end]))))
            return

        best_t = -1
        best_stat = 0
        overall_mean = np.mean(data[start:end])

        for t in range(start + min_segment, end - min_segment):
            left_mean = np.mean(data[start:t])
            right_mean = np.mean(data[t:end])
            n_left = t - start
            n_right = end - t
            stat = abs(left_mean - right_mean) * math.sqrt(
                n_left * n_right / (n_left + n_right)
            )
            if stat > best_stat:
                best_stat = stat
                best_t = t

        if best_stat > 2.0 and best_t > 0:
            self._cbs_recursive(data, start, best_t, segments, min_segment)
            self._cbs_recursive(data, best_t, end, segments, min_segment)
        else:
            segments.append((start, end, float(overall_mean)))


class LocalAssembler:
    """Local assembly for SV breakpoint refinement."""

    def __init__(self, kmer_size: int = 21, min_coverage: int = 3):
        self.kmer_size = kmer_size
        self.min_coverage = min_coverage

    def build_dbg(self, reads: List[str]) -> Dict[str, List[str]]:
        """Build a de Bruijn graph from reads."""
        graph = defaultdict(list)
        kmer_counts = defaultdict(int)

        for read in reads:
            for i in range(len(read) - self.kmer_size):
                kmer = read[i:i + self.kmer_size]
                kmer_counts[kmer] += 1

        for read in reads:
            for i in range(len(read) - self.kmer_size - 1):
                kmer1 = read[i:i + self.kmer_size]
                kmer2 = read[i + 1:i + 1 + self.kmer_size]
                if (kmer_counts[kmer1] >= self.min_coverage and
                        kmer_counts[kmer2] >= self.min_coverage):
                    if kmer2 not in graph[kmer1]:
                        graph[kmer1].append(kmer2)

        return dict(graph)

    def assemble_contigs(self, reads: List[str]) -> List[str]:
        graph = self.build_dbg(reads)
        if not graph:
            return []

        contigs = []
        visited: Set[str] = set()

        for start_kmer in graph:
            if start_kmer in visited:
                continue
            contig = start_kmer
            current = start_kmer
            visited.add(current)

            while current in graph:
                neighbors = [n for n in graph[current] if n not in visited]
                if not neighbors:
                    break
                current = neighbors[0]
                visited.add(current)
                contig += current[-1]

            if len(contig) > self.kmer_size * 2:
                contigs.append(contig)

        return contigs

    def refine_breakpoints(
        self, sv: SVCall, contigs: List[str], reference: str
    ) -> SVCall:
        """Refine SV breakpoints using assembled contigs."""
        if not contigs:
            return sv

        best_contig = max(contigs, key=len)
        ref_region = reference[max(0, sv.start - 100):sv.end + 100]

        # Simple alignment to refine breakpoints
        best_score = 0
        best_offset = 0
        for offset in range(-50, 51):
            score = self._quick_align_score(best_contig, ref_region, offset)
            if score > best_score:
                best_score = score
                best_offset = offset

        refined = SVCall(
            sv_type=sv.sv_type, chrom=sv.chrom,
            start=sv.start + best_offset,
            end=sv.end + best_offset,
            size=sv.size, quality=sv.quality + 5,
            genotype=sv.genotype,
            support_reads=sv.support_reads,
            evidence_types=sv.evidence_types + ["assembly"]
        )
        return refined

    @staticmethod
    def _quick_align_score(query: str, ref: str, offset: int) -> int:
        score = 0
        start = max(0, offset)
        for i in range(min(len(query), len(ref) - start)):
            if i + start < len(ref) and i < len(query):
                if query[i] == ref[i + start]:
                    score += 1
        return score


class IntegratedSVDetector:
    """Integrates split-read, read-depth, and assembly-based SV detection."""

    def __init__(
        self, min_support: int = 3, min_quality: float = 20.0,
        merge_distance: int = 500
    ):
        self.split_analyzer = SplitReadAnalyzer()
        self.depth_analyzer = ReadDepthAnalyzer()
        self.assembler = LocalAssembler()
        self.min_support = min_support
        self.min_quality = min_quality
        self.merge_distance = merge_distance

    def detect(
        self, split_segments: List[List[AlignmentSegment]],
        alignments: List[Dict], chrom: str, chrom_length: int,
        reads: Optional[List[str]] = None,
        reference: Optional[str] = None
    ) -> List[SVCall]:
        # Phase 1: Split-read detection
        split_calls = []
        for segments in split_segments:
            split_calls.extend(self.split_analyzer.detect_from_splits(segments))

        # Phase 2: Read-depth detection
        depth_signal = self.depth_analyzer.compute_depth_profile(
            alignments, chrom, chrom_length
        )
        depth_calls = self.depth_analyzer.detect_cnv(depth_signal)

        # Phase 3: Local assembly refinement
        all_calls = split_calls + depth_calls
        if reads and reference:
            refined_calls = []
            for sv in all_calls:
                region_reads = self._extract_region_reads(
                    reads, sv.start, sv.end
                )
                if region_reads:
                    contigs = self.assembler.assemble_contigs(region_reads)
                    refined = self.assembler.refine_breakpoints(
                        sv, contigs, reference
                    )
                    refined_calls.append(refined)
                else:
                    refined_calls.append(sv)
            all_calls = refined_calls

        # Phase 4: Merge and filter
        merged = self._merge_calls(all_calls)
        filtered = [c for c in merged
                    if c.support_reads >= self.min_support
                    or len(c.evidence_types) >= 2]
        return filtered

    def _merge_calls(self, calls: List[SVCall]) -> List[SVCall]:
        if not calls:
            return []

        calls.sort(key=lambda c: (c.chrom, c.start))
        merged = [calls[0]]

        for call in calls[1:]:
            prev = merged[-1]
            if (call.chrom == prev.chrom and
                    call.sv_type == prev.sv_type and
                    abs(call.start - prev.start) < self.merge_distance):
                # Merge
                prev.support_reads += call.support_reads
                prev.quality = max(prev.quality, call.quality)
                prev.evidence_types = list(
                    set(prev.evidence_types + call.evidence_types)
                )
                prev.end = max(prev.end, call.end)
            else:
                merged.append(call)

        return merged

    @staticmethod
    def _extract_region_reads(
        reads: List[str], start: int, end: int
    ) -> List[str]:
        return reads[:min(100, len(reads))]
