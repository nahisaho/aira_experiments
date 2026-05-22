"""Structural variant calling engines used by DeepSV-LR."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, DefaultDict, Dict, Iterable, List, Mapping, MutableSequence, Optional, Sequence, Tuple

import numpy as np


class SVType(str, Enum):
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    BND = "BND"
    CNV = "CNV"


class EvidenceType(str, Enum):
    SPLIT_READ = "split_read"
    READ_DEPTH = "read_depth"
    ASSEMBLY = "assembly"
    SHORT_READ = "short_read"
    REPEAT = "repeat"


@dataclass(frozen=True)
class Breakpoint:
    """Representation of a genomic breakpoint."""

    chrom: str
    position: int
    orientation: str = "+"
    confidence_interval: Tuple[int, int] = (0, 0)


@dataclass
class SVEvidence:
    """Evidence item supporting a structural-variant candidate."""

    source: EvidenceType
    weight: float
    support_reads: int = 0
    score: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SVCandidate:
    """Canonical structural-variant representation used throughout the pipeline."""

    chrom: str
    start: int
    end: int
    svtype: SVType
    size: int
    left_breakpoint: Breakpoint
    right_breakpoint: Breakpoint
    evidence: List[SVEvidence] = field(default_factory=list)
    quality: float = 0.0
    sequence: Optional[str] = None
    copy_number: Optional[float] = None
    genotype: Optional[str] = None
    sample: str = "sample"
    info: Dict[str, Any] = field(default_factory=dict)

    def add_evidence(self, evidence_item: SVEvidence) -> None:
        self.evidence.append(evidence_item)

    @property
    def evidence_weight(self) -> float:
        return float(sum(item.weight for item in self.evidence))


class SplitReadCaller:
    """Detect structural variants from split and supplementary alignments."""

    def __init__(self, min_event_size: int = 30, min_mapq: int = 20) -> None:
        self.min_event_size = min_event_size
        self.min_mapq = min_mapq

    def call(self, alignments: Sequence[Mapping[str, Any]]) -> List[SVCandidate]:
        candidates: List[SVCandidate] = []
        for alignment in alignments:
            segments = sorted(alignment.get("segments", []), key=lambda segment: segment.get("query_start", 0))
            for left, right in zip(segments, segments[1:]):
                if min(left.get("mapq", 0), right.get("mapq", 0)) < self.min_mapq:
                    continue
                candidate = self._segment_pair_to_candidate(alignment.get("read_name", "read"), left, right)
                if candidate is not None:
                    candidates.append(candidate)
        return merge_sv_candidates(candidates)

    def _segment_pair_to_candidate(
        self,
        read_name: str,
        left: Mapping[str, Any],
        right: Mapping[str, Any],
    ) -> Optional[SVCandidate]:
        left_chrom = str(left.get("chrom", ""))
        right_chrom = str(right.get("chrom", ""))
        left_end = int(left.get("ref_end", left.get("ref_start", 0)))
        right_start = int(right.get("ref_start", right.get("ref_end", 0)))
        query_gap = int(right.get("query_start", 0)) - int(left.get("query_end", 0))
        ref_gap = right_start - left_end
        same_strand = left.get("strand", "+") == right.get("strand", "+")

        if left_chrom != right_chrom:
            svtype = SVType.BND
            start, end = sorted((left_end, right_start))
        elif not same_strand:
            svtype = SVType.INV
            start, end = sorted((left_end, right_start))
        elif ref_gap > self.min_event_size and query_gap <= self.min_event_size:
            svtype = SVType.DEL
            start, end = left_end, right_start
        elif query_gap > self.min_event_size and ref_gap <= self.min_event_size:
            svtype = SVType.INS
            start = left_end
            end = left_end + max(query_gap, 1)
        elif ref_gap < -self.min_event_size:
            svtype = SVType.DUP
            start, end = right_start, left_end
        else:
            return None

        start, end = int(min(start, end)), int(max(start, end))
        size = max(end - start, abs(query_gap), abs(ref_gap), 1)
        evidence = SVEvidence(
            source=EvidenceType.SPLIT_READ,
            weight=1.5,
            support_reads=1,
            score=float(min(left.get("mapq", 60), right.get("mapq", 60))) / 60.0,
            metadata={"read_name": read_name},
        )
        return SVCandidate(
            chrom=left_chrom,
            start=start,
            end=end,
            svtype=svtype,
            size=size,
            left_breakpoint=Breakpoint(left_chrom, start, str(left.get("strand", "+"))),
            right_breakpoint=Breakpoint(right_chrom, end, str(right.get("strand", "+"))),
            evidence=[evidence],
            quality=20.0 + 20.0 * evidence.score,
        )


class ReadDepthCaller:
    """Copy-number caller based on circular binary segmentation (CBS)."""

    def __init__(self, bin_size: int = 1000, min_bins: int = 4, z_threshold: float = 2.5) -> None:
        self.bin_size = bin_size
        self.min_bins = min_bins
        self.z_threshold = z_threshold

    def call(self, depth_profile: Mapping[str, Sequence[float]]) -> List[SVCandidate]:
        candidates: List[SVCandidate] = []
        for chrom, bins in depth_profile.items():
            values = np.asarray(bins, dtype=np.float64)
            if values.size < self.min_bins * 2:
                continue
            baseline = float(np.median(values)) if np.any(values) else 1.0
            for segment_start, segment_end, mean_depth in self._segment(values):
                copy_number = 2.0 * mean_depth / max(baseline, 1e-6)
                if 1.6 <= copy_number <= 2.4:
                    continue
                svtype = SVType.DEL if copy_number < 1.6 else SVType.DUP
                start = segment_start * self.bin_size
                end = segment_end * self.bin_size
                evidence = SVEvidence(
                    source=EvidenceType.READ_DEPTH,
                    weight=1.0,
                    support_reads=max(segment_end - segment_start, 1),
                    score=abs(copy_number - 2.0) / 2.0,
                    metadata={"mean_depth": mean_depth},
                )
                candidates.append(
                    SVCandidate(
                        chrom=chrom,
                        start=start,
                        end=end,
                        svtype=svtype,
                        size=max(end - start, self.bin_size),
                        left_breakpoint=Breakpoint(chrom, start),
                        right_breakpoint=Breakpoint(chrom, end),
                        evidence=[evidence],
                        quality=15.0 + 10.0 * evidence.score,
                        copy_number=copy_number,
                    )
                )
        return merge_sv_candidates(candidates)

    def _segment(self, values: np.ndarray) -> List[Tuple[int, int, float]]:
        segments: List[Tuple[int, int, float]] = []

        def recurse(start: int, end: int) -> None:
            segment = values[start:end]
            if end - start < self.min_bins * 2:
                segments.append((start, end, float(np.mean(segment))))
                return

            best_score = 0.0
            best_index: Optional[int] = None
            for split in range(start + self.min_bins, end - self.min_bins + 1):
                left = values[start:split]
                right = values[split:end]
                denominator = np.sqrt(np.var(left) / max(left.size, 1) + np.var(right) / max(right.size, 1) + 1e-6)
                score = abs(np.mean(left) - np.mean(right)) / denominator
                if score > best_score:
                    best_score = float(score)
                    best_index = split

            if best_index is not None and best_score >= self.z_threshold:
                recurse(start, best_index)
                recurse(best_index, end)
            else:
                segments.append((start, end, float(np.mean(segment))))

        recurse(0, len(values))
        return sorted(segments, key=lambda item: item[0])


class AssemblyCaller:
    """Local assembly around candidate breakpoints using greedy overlap assembly."""

    def __init__(self, min_overlap: int = 20, min_event_size: int = 30) -> None:
        self.min_overlap = min_overlap
        self.min_event_size = min_event_size

    def call(self, regions: Sequence[Mapping[str, Any]]) -> List[SVCandidate]:
        candidates: List[SVCandidate] = []
        for region in regions:
            reads = [str(read) for read in region.get("reads", []) if read]
            if len(reads) < 2:
                continue
            contig = self._assemble_reads(reads)
            reference_span = int(region.get("end", 0)) - int(region.get("start", 0))
            delta = len(contig) - max(reference_span, 1)
            if abs(delta) < self.min_event_size:
                continue
            svtype = SVType.INS if delta > 0 else SVType.DEL
            chrom = str(region.get("chrom", ""))
            start = int(region.get("start", 0))
            end = int(region.get("end", start + 1))
            evidence = SVEvidence(
                source=EvidenceType.ASSEMBLY,
                weight=1.3,
                support_reads=len(reads),
                score=min(abs(delta) / max(reference_span, 1), 5.0),
                metadata={"contig_length": len(contig)},
            )
            candidates.append(
                SVCandidate(
                    chrom=chrom,
                    start=start,
                    end=end if svtype != SVType.INS else start + max(abs(delta), 1),
                    svtype=svtype,
                    size=max(abs(delta), end - start, 1),
                    left_breakpoint=Breakpoint(chrom, start),
                    right_breakpoint=Breakpoint(chrom, end),
                    evidence=[evidence],
                    quality=18.0 + 8.0 * evidence.score,
                    sequence=contig if svtype == SVType.INS else None,
                )
            )
        return merge_sv_candidates(candidates)

    def _assemble_reads(self, reads: MutableSequence[str]) -> str:
        reads = list(reads)
        while len(reads) > 1:
            best_pair: Optional[Tuple[int, int, int, str]] = None
            for left_index, left in enumerate(reads):
                for right_index, right in enumerate(reads):
                    if left_index == right_index:
                        continue
                    overlap = self._suffix_prefix_overlap(left, right)
                    if overlap >= self.min_overlap:
                        merged = left + right[overlap:]
                        if best_pair is None or overlap > best_pair[2]:
                            best_pair = (left_index, right_index, overlap, merged)
            if best_pair is None:
                break
            left_index, right_index, _, merged = best_pair
            first, second = sorted((left_index, right_index), reverse=True)
            reads.pop(first)
            reads.pop(second)
            reads.append(merged)
        return max(reads, key=len)

    @staticmethod
    def _suffix_prefix_overlap(left: str, right: str) -> int:
        max_possible = min(len(left), len(right))
        for overlap in range(max_possible, 0, -1):
            if left[-overlap:] == right[:overlap]:
                return overlap
        return 0


class EnsembleSVCaller:
    """Merge calls from split-read, read-depth and assembly callers."""

    def __init__(self, weights: Optional[Mapping[EvidenceType, float]] = None, overlap_threshold: float = 0.5) -> None:
        self.weights = {
            EvidenceType.SPLIT_READ: 1.2,
            EvidenceType.READ_DEPTH: 1.0,
            EvidenceType.ASSEMBLY: 1.4,
            EvidenceType.SHORT_READ: 1.1,
            EvidenceType.REPEAT: 0.7,
        }
        if weights is not None:
            self.weights.update(weights)
        self.overlap_threshold = overlap_threshold

    def call(
        self,
        split_read_calls: Sequence[SVCandidate],
        read_depth_calls: Sequence[SVCandidate],
        assembly_calls: Sequence[SVCandidate],
    ) -> List[SVCandidate]:
        flattened = list(split_read_calls) + list(read_depth_calls) + list(assembly_calls)
        merged = merge_sv_candidates(flattened, overlap_threshold=self.overlap_threshold)
        for candidate in merged:
            weighted_support = 0.0
            total_weight = 0.0
            for evidence in candidate.evidence:
                model_weight = self.weights.get(evidence.source, 1.0)
                weighted_support += model_weight * max(evidence.score, 0.1)
                total_weight += model_weight
            candidate.quality = 10.0 + 20.0 * (weighted_support / max(total_weight, 1e-6))
        return merged


def reciprocal_overlap(first: SVCandidate, second: SVCandidate, insertion_tolerance: int = 50) -> float:
    """Compute reciprocal overlap used for call merging.

    Insertions are handled with a positional tolerance because they often have
    imprecise spans.
    """

    if first.chrom != second.chrom or first.svtype != second.svtype:
        return 0.0
    if first.svtype == SVType.INS:
        distance = abs(first.start - second.start)
        return 1.0 if distance <= insertion_tolerance else 0.0

    left = max(first.start, second.start)
    right = min(first.end, second.end)
    intersection = max(0, right - left)
    if intersection == 0:
        return 0.0
    len_first = max(first.end - first.start, 1)
    len_second = max(second.end - second.start, 1)
    return min(intersection / len_first, intersection / len_second)


def merge_sv_candidates(
    candidates: Iterable[SVCandidate],
    overlap_threshold: float = 0.5,
) -> List[SVCandidate]:
    """Merge highly overlapping SV calls using reciprocal overlap clustering."""

    merged_clusters: List[List[SVCandidate]] = []
    for candidate in sorted(candidates, key=lambda item: (item.chrom, item.start, item.end, item.svtype.value)):
        placed = False
        for cluster in merged_clusters:
            if reciprocal_overlap(candidate, cluster[0]) >= overlap_threshold:
                cluster.append(candidate)
                placed = True
                break
        if not placed:
            merged_clusters.append([candidate])

    merged: List[SVCandidate] = []
    for cluster in merged_clusters:
        merged.append(_collapse_cluster(cluster))
    return merged


def _collapse_cluster(cluster: Sequence[SVCandidate]) -> SVCandidate:
    weights = np.asarray([max(candidate.evidence_weight, 1.0) for candidate in cluster], dtype=np.float64)
    starts = np.asarray([candidate.start for candidate in cluster], dtype=np.float64)
    ends = np.asarray([candidate.end for candidate in cluster], dtype=np.float64)
    representative = max(cluster, key=lambda candidate: candidate.evidence_weight + candidate.quality)

    merged = SVCandidate(
        chrom=representative.chrom,
        start=int(np.average(starts, weights=weights)),
        end=int(np.average(ends, weights=weights)),
        svtype=representative.svtype,
        size=int(max(np.average(np.asarray([candidate.size for candidate in cluster], dtype=np.float64), weights=weights), 1)),
        left_breakpoint=Breakpoint(representative.chrom, int(np.average(starts, weights=weights))),
        right_breakpoint=Breakpoint(representative.chrom, int(np.average(ends, weights=weights))),
        evidence=[evidence for candidate in cluster for evidence in candidate.evidence],
        quality=float(np.average(np.asarray([candidate.quality for candidate in cluster], dtype=np.float64), weights=weights)),
        sequence=representative.sequence,
        copy_number=representative.copy_number,
        genotype=representative.genotype,
        sample=representative.sample,
        info=dict(representative.info),
    )
    merged.info["cluster_size"] = len(cluster)
    return merged


__all__ = [
    "AssemblyCaller",
    "Breakpoint",
    "EnsembleSVCaller",
    "EvidenceType",
    "ReadDepthCaller",
    "SVCandidate",
    "SVEvidence",
    "SVType",
    "SplitReadCaller",
    "merge_sv_candidates",
    "reciprocal_overlap",
]
