#!/usr/bin/env python3
"""
Repeat region handler for telomeric and centromeric regions.
Provides specialized SV detection in repetitive genomic regions.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
from collections import Counter
import re
import logging

logger = logging.getLogger(__name__)


@dataclass
class RepeatAnnotation:
    """Annotation for a repeat region."""
    chrom: str
    start: int
    end: int
    repeat_type: str  # telomere, centromere, satellite, LINE, SINE, etc.
    repeat_unit: str
    copy_number: float
    gc_content: float


@dataclass
class RepeatMaskedRegion:
    """Region with repeat masking information."""
    chrom: str
    start: int
    end: int
    masked_fraction: float
    repeat_classes: List[str]
    complexity_score: float


# Human telomere canonical repeat
TELOMERE_REPEAT = "TTAGGG"
TELOMERE_REVERSE = "CCCTAA"

# Human centromere alpha-satellite monomer length
ALPHA_SATELLITE_LENGTH = 171

# Known centromere positions (hg38, approximate)
CENTROMERE_REGIONS = {
    "chr1": (122026460, 125184587),
    "chr2": (92188146, 94090557),
    "chr3": (90772459, 93655574),
    "chr4": (49708101, 51743951),
    "chr5": (46485901, 50059807),
    "chr6": (58553889, 59829934),
    "chr7": (58169654, 60828234),
    "chr8": (44033745, 45877265),
    "chr9": (43236168, 45518558),
    "chr10": (39686683, 41593521),
    "chr11": (51078349, 54425074),
    "chr12": (34769408, 37185252),
    "chr13": (16000001, 18051248),
    "chr14": (16000001, 18173523),
    "chr15": (17083674, 19725254),
    "chr16": (36311159, 38280682),
    "chr17": (22813680, 26885980),
    "chr18": (15460900, 20861206),
    "chr19": (24498981, 27190874),
    "chr20": (26436233, 30038348),
    "chr21": (10864561, 12915808),
    "chr22": (12954789, 15054318),
    "chrX": (58605580, 62412542),
    "chrY": (10316945, 10544039),
}


class TelomereAnalyzer:
    """Analyze telomeric regions and detect telomere-associated SVs."""

    def __init__(self, min_repeat_count: int = 3):
        self.min_repeat_count = min_repeat_count
        self.telomere_pattern = re.compile(
            f"({TELOMERE_REPEAT}){{" + str(min_repeat_count) + ",}"
        )
        self.telomere_rev_pattern = re.compile(
            f"({TELOMERE_REVERSE}){{" + str(min_repeat_count) + ",}"
        )

    def detect_telomere_reads(self, sequences: List[str]) -> List[Dict]:
        """Identify reads containing telomeric repeats."""
        telomere_reads = []
        for i, seq in enumerate(sequences):
            fwd_matches = list(self.telomere_pattern.finditer(seq))
            rev_matches = list(self.telomere_rev_pattern.finditer(seq))

            if fwd_matches or rev_matches:
                all_matches = fwd_matches + rev_matches
                total_telomere_bp = sum(m.end() - m.start() for m in all_matches)
                telomere_reads.append({
                    "read_index": i,
                    "telomere_fraction": total_telomere_bp / len(seq),
                    "n_telomere_blocks": len(all_matches),
                    "is_terminal": (
                        any(m.start() < 50 for m in all_matches) or
                        any(m.end() > len(seq) - 50 for m in all_matches)
                    ),
                    "orientations": (
                        ["forward"] * len(fwd_matches) +
                        ["reverse"] * len(rev_matches)
                    ),
                })

        return telomere_reads

    def estimate_telomere_length(self, sequences: List[str]) -> Dict:
        """Estimate telomere length from long reads spanning telomeric regions."""
        lengths = []
        for seq in sequences:
            for match in self.telomere_pattern.finditer(seq):
                lengths.append(match.end() - match.start())
            for match in self.telomere_rev_pattern.finditer(seq):
                lengths.append(match.end() - match.start())

        if not lengths:
            return {"mean_length": 0, "median_length": 0, "n_estimates": 0}

        return {
            "mean_length": float(np.mean(lengths)),
            "median_length": float(np.median(lengths)),
            "std_length": float(np.std(lengths)),
            "n_estimates": len(lengths),
            "min_length": int(min(lengths)),
            "max_length": int(max(lengths)),
        }


class CentromereAnalyzer:
    """Analyze centromeric regions with alpha-satellite repeats."""

    def __init__(self):
        self.alpha_sat_length = ALPHA_SATELLITE_LENGTH
        self.regions = CENTROMERE_REGIONS

    def is_centromeric(self, chrom: str, start: int, end: int) -> bool:
        if chrom in self.regions:
            cen_start, cen_end = self.regions[chrom]
            return start < cen_end and end > cen_start
        return False

    def compute_repeat_complexity(self, sequence: str) -> float:
        """Compute linguistic complexity of a sequence (0=simple repeat, 1=complex)."""
        if len(sequence) < 4:
            return 0.0

        # Trigram complexity
        possible_trigrams = min(len(sequence) - 2, 64)
        observed = len(set(
            sequence[i:i + 3] for i in range(len(sequence) - 2)
        ))
        return observed / possible_trigrams if possible_trigrams > 0 else 0.0

    def detect_hor_variants(
        self, sequence: str, monomer_length: int = 171
    ) -> List[Dict]:
        """Detect Higher-Order Repeat (HOR) structural variants."""
        if len(sequence) < monomer_length * 2:
            return []

        monomers = []
        for i in range(0, len(sequence) - monomer_length, monomer_length):
            monomers.append(sequence[i:i + monomer_length])

        variants = []
        for i in range(len(monomers) - 1):
            divergence = self._sequence_divergence(monomers[i], monomers[i + 1])
            if divergence > 0.15:
                variants.append({
                    "position": i * monomer_length,
                    "divergence": divergence,
                    "type": "hor_variant",
                })

        return variants

    @staticmethod
    def _sequence_divergence(seq1: str, seq2: str) -> float:
        min_len = min(len(seq1), len(seq2))
        if min_len == 0:
            return 1.0
        mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
        return mismatches / min_len


class RepeatHandler:
    """Unified handler for repeat regions in SV detection."""

    def __init__(self):
        self.telomere_analyzer = TelomereAnalyzer()
        self.centromere_analyzer = CentromereAnalyzer()

    def classify_region(
        self, chrom: str, start: int, end: int, sequence: str
    ) -> RepeatAnnotation:
        """Classify a genomic region by its repeat content."""
        # Check centromere
        if self.centromere_analyzer.is_centromeric(chrom, start, end):
            complexity = self.centromere_analyzer.compute_repeat_complexity(sequence)
            gc = self._compute_gc(sequence)
            return RepeatAnnotation(
                chrom=chrom, start=start, end=end,
                repeat_type="centromere",
                repeat_unit="alpha-satellite",
                copy_number=(end - start) / ALPHA_SATELLITE_LENGTH,
                gc_content=gc
            )

        # Check telomere
        tel_reads = self.telomere_analyzer.detect_telomere_reads([sequence])
        if tel_reads and tel_reads[0]["telomere_fraction"] > 0.3:
            gc = self._compute_gc(sequence)
            return RepeatAnnotation(
                chrom=chrom, start=start, end=end,
                repeat_type="telomere",
                repeat_unit=TELOMERE_REPEAT,
                copy_number=tel_reads[0]["telomere_fraction"] * len(sequence) / 6,
                gc_content=gc
            )

        # Default: non-repeat
        gc = self._compute_gc(sequence)
        return RepeatAnnotation(
            chrom=chrom, start=start, end=end,
            repeat_type="unique",
            repeat_unit="",
            copy_number=1.0,
            gc_content=gc
        )

    def adjust_sv_confidence(
        self, sv_call, repeat_annotation: RepeatAnnotation
    ) -> float:
        """Adjust SV confidence based on repeat context."""
        base_quality = sv_call.quality

        if repeat_annotation.repeat_type == "centromere":
            return base_quality * 0.6  # Lower confidence in centromeric regions
        elif repeat_annotation.repeat_type == "telomere":
            return base_quality * 0.7
        elif repeat_annotation.repeat_type in ("LINE", "SINE", "satellite"):
            return base_quality * 0.8

        return base_quality

    @staticmethod
    def _compute_gc(sequence: str) -> float:
        if not sequence:
            return 0.0
        gc_count = sum(1 for c in sequence.upper() if c in ('G', 'C'))
        return gc_count / len(sequence)
