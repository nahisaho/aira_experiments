"""Short-read integration utilities for DeepSV-LR."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    from .sv_detector import Breakpoint, EvidenceType, SVCandidate, SVEvidence, SVType, merge_sv_candidates, reciprocal_overlap
except ImportError:  # pragma: no cover - fallback for flat module execution
    from sv_detector import Breakpoint, EvidenceType, SVCandidate, SVEvidence, SVType, merge_sv_candidates, reciprocal_overlap


class Genotype(str, Enum):
    HOM_REF = "0/0"
    HET = "0/1"
    HOM_ALT = "1/1"


@dataclass(frozen=True)
class ShortReadEvidence:
    chrom: str
    start: int
    end: int
    split_reads: int
    discordant_pairs: int
    reference_reads: int
    precision_positions: Tuple[int, int] = (0, 0)
    population_frequency: Optional[float] = None


class HybridIntegrator:
    """Combine short-read and long-read evidence for refined SV genotyping."""

    def __init__(self, overlap_threshold: float = 0.5) -> None:
        self.overlap_threshold = overlap_threshold

    def overlay_short_read_evidence(
        self,
        sv_calls: Sequence[SVCandidate],
        short_read_evidence: Sequence[ShortReadEvidence],
    ) -> List[SVCandidate]:
        """Attach short-read split/paired-end support to long-read calls."""

        augmented: List[SVCandidate] = []
        for candidate in sv_calls:
            updated = self._copy_candidate(candidate)
            for evidence in short_read_evidence:
                if evidence.chrom != candidate.chrom:
                    continue
                proxy = self._proxy_candidate(evidence)
                if reciprocal_overlap(candidate, proxy) < self.overlap_threshold:
                    continue
                support = evidence.split_reads + 0.5 * evidence.discordant_pairs
                updated.add_evidence(
                    SVEvidence(
                        source=EvidenceType.SHORT_READ,
                        weight=1.1,
                        support_reads=int(support),
                        score=float(support) / max(evidence.reference_reads + support, 1.0),
                        metadata={"precision_positions": evidence.precision_positions},
                    )
                )
                updated = self.enhance_breakpoint_precision(updated, evidence.precision_positions)
            augmented.append(updated)
        return merge_sv_candidates(augmented, overlap_threshold=self.overlap_threshold)

    def refine_genotype(
        self,
        candidate: SVCandidate,
        prior_alt_fraction: float = 0.01,
    ) -> Tuple[str, float]:
        """Infer genotype via a Bayesian binomial model over alt-support fractions."""

        alt_reads = float(sum(item.support_reads for item in candidate.evidence if item.source != EvidenceType.READ_DEPTH))
        ref_reads = float(candidate.info.get("reference_reads", max(2.0, alt_reads)))
        total_reads = max(alt_reads + ref_reads, 1.0)
        priors = {
            Genotype.HOM_REF: max(1.0 - prior_alt_fraction, 1e-3),
            Genotype.HET: 0.5,
            Genotype.HOM_ALT: max(prior_alt_fraction, 1e-3),
        }
        expectations = {
            Genotype.HOM_REF: 0.02,
            Genotype.HET: 0.5,
            Genotype.HOM_ALT: 0.98,
        }
        log_posteriors: Dict[Genotype, float] = {}
        for genotype, expected_alt_fraction in expectations.items():
            p = min(max(expected_alt_fraction, 1e-6), 1.0 - 1e-6)
            log_likelihood = alt_reads * np.log(p) + ref_reads * np.log(1.0 - p)
            log_posteriors[genotype] = np.log(priors[genotype]) + log_likelihood
        normalization = _logsumexp(list(log_posteriors.values()))
        posterior_probs = {genotype: float(np.exp(score - normalization)) for genotype, score in log_posteriors.items()}
        genotype = max(posterior_probs, key=posterior_probs.get)
        candidate.genotype = genotype.value
        candidate.info["genotype_posterior"] = posterior_probs[genotype]
        candidate.info["total_reads"] = total_reads
        return genotype.value, posterior_probs[genotype]

    def enhance_breakpoint_precision(
        self,
        candidate: SVCandidate,
        split_read_positions: Tuple[int, int],
    ) -> SVCandidate:
        """Refine long-read breakpoints using high-accuracy short-read split reads."""

        left, right = split_read_positions
        if left > 0:
            candidate.start = int(round((candidate.start + left) / 2))
            candidate.left_breakpoint = type(candidate.left_breakpoint)(
                candidate.left_breakpoint.chrom,
                candidate.start,
                candidate.left_breakpoint.orientation,
                (-3, 3),
            )
        if right > 0:
            candidate.end = int(round((candidate.end + right) / 2))
            candidate.right_breakpoint = type(candidate.right_breakpoint)(
                candidate.right_breakpoint.chrom,
                candidate.end,
                candidate.right_breakpoint.orientation,
                (-3, 3),
            )
        candidate.size = max(candidate.end - candidate.start, 1)
        return candidate

    def annotate_population_frequency(
        self,
        sv_calls: Sequence[SVCandidate],
        frequency_panel: Mapping[str, float],
    ) -> List[SVCandidate]:
        """Annotate population frequency using a keyed external SV panel."""

        annotated: List[SVCandidate] = []
        for candidate in sv_calls:
            updated = self._copy_candidate(candidate)
            key = self._frequency_key(candidate)
            if key in frequency_panel:
                updated.info["population_frequency"] = float(frequency_panel[key])
            annotated.append(updated)
        return annotated

    def call(
        self,
        sv_calls: Sequence[SVCandidate],
        short_read_evidence: Sequence[ShortReadEvidence],
        frequency_panel: Optional[Mapping[str, float]] = None,
    ) -> List[SVCandidate]:
        calls = self.overlay_short_read_evidence(sv_calls, short_read_evidence)
        if frequency_panel is not None:
            calls = self.annotate_population_frequency(calls, frequency_panel)
        for candidate in calls:
            self.refine_genotype(candidate, prior_alt_fraction=float(candidate.info.get("population_frequency", 0.01)))
        return calls

    @staticmethod
    def _frequency_key(candidate: SVCandidate) -> str:
        return f"{candidate.chrom}:{candidate.start}-{candidate.end}:{candidate.svtype.value}"

    @staticmethod
    def _proxy_candidate(evidence: ShortReadEvidence) -> SVCandidate:
        return SVCandidate(
            chrom=evidence.chrom,
            start=evidence.start,
            end=evidence.end,
            svtype=SVType.DEL if evidence.end > evidence.start else SVType.INS,
            size=max(evidence.end - evidence.start, 1),
            left_breakpoint=Breakpoint(evidence.chrom, evidence.start),
            right_breakpoint=Breakpoint(evidence.chrom, evidence.end),
        )

    @staticmethod
    def _copy_candidate(candidate: SVCandidate) -> SVCandidate:
        return SVCandidate(
            chrom=candidate.chrom,
            start=candidate.start,
            end=candidate.end,
            svtype=candidate.svtype,
            size=candidate.size,
            left_breakpoint=candidate.left_breakpoint,
            right_breakpoint=candidate.right_breakpoint,
            evidence=list(candidate.evidence),
            quality=candidate.quality,
            sequence=candidate.sequence,
            copy_number=candidate.copy_number,
            genotype=candidate.genotype,
            sample=candidate.sample,
            info=dict(candidate.info),
        )


def _logsumexp(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=np.float64)
    maximum = np.max(array)
    return float(maximum + np.log(np.sum(np.exp(array - maximum))))


__all__ = ["Genotype", "HybridIntegrator", "ShortReadEvidence"]
