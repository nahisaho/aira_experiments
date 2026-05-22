"""GIAB-style structural-variant benchmarking for DeepSV-LR."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from .sv_detector import Breakpoint, SVCandidate, SVType, reciprocal_overlap
except ImportError:  # pragma: no cover - fallback for flat module execution
    from sv_detector import Breakpoint, SVCandidate, SVType, reciprocal_overlap


class GIABBenchmark:
    """Evaluate predicted SVs against Genome in a Bottle style truth sets."""

    def __init__(self, overlap_threshold: float = 0.5, size_similarity: float = 0.7) -> None:
        self.overlap_threshold = overlap_threshold
        self.size_similarity = size_similarity
        self.truth_set: List[SVCandidate] = []

    def load_truth_set(self, path: str | Path) -> List[SVCandidate]:
        """Load truth calls from JSON or TSV/BED-like tabular input."""

        truth_path = Path(path)
        if truth_path.suffix.lower() == ".json":
            records = json.loads(truth_path.read_text(encoding="utf-8"))
        else:
            records = []
            for line in truth_path.read_text(encoding="utf-8").splitlines():
                if not line or line.startswith("#"):
                    continue
                chrom, start, end, svtype, *rest = line.split("\t")
                records.append({"chrom": chrom, "start": int(start), "end": int(end), "svtype": svtype, "id": rest[0] if rest else None})

        self.truth_set = [self._record_to_candidate(record) for record in records]
        return self.truth_set

    def evaluate(
        self,
        predictions: Sequence[SVCandidate],
        truth_set: Optional[Sequence[SVCandidate]] = None,
    ) -> Dict[str, Any]:
        """Run a Truvari-like evaluation with one-to-one truth matching."""

        truth = list(truth_set if truth_set is not None else self.truth_set)
        matches, false_positives, false_negatives = self._match_calls(predictions, truth)
        metrics = self.calculate_metrics(len(matches), len(false_positives), len(false_negatives))
        metrics["matches"] = [
            {
                "truth": f"{truth_call.chrom}:{truth_call.start}-{truth_call.end}:{truth_call.svtype.value}",
                "query": f"{query_call.chrom}:{query_call.start}-{query_call.end}:{query_call.svtype.value}",
            }
            for truth_call, query_call in matches
        ]
        metrics["stratified_by_svtype"] = self.stratify_by_svtype(predictions, truth)
        metrics["stratified_by_size"] = self.stratify_by_size(predictions, truth)
        return metrics

    @staticmethod
    def calculate_metrics(true_positives: int, false_positives: int, false_negatives: int) -> Dict[str, float]:
        precision = true_positives / max(true_positives + false_positives, 1)
        recall = true_positives / max(true_positives + false_negatives, 1)
        f1 = 2.0 * precision * recall / max(precision + recall, 1e-9)
        return {
            "true_positives": float(true_positives),
            "false_positives": float(false_positives),
            "false_negatives": float(false_negatives),
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }

    def stratify_by_region(
        self,
        predictions: Sequence[SVCandidate],
        regions: Mapping[str, Sequence[Tuple[str, int, int]]],
        truth_set: Optional[Sequence[SVCandidate]] = None,
    ) -> Dict[str, Dict[str, float]]:
        truth = list(truth_set if truth_set is not None else self.truth_set)
        results: Dict[str, Dict[str, float]] = {}
        for label, intervals in regions.items():
            pred_subset = [call for call in predictions if self._in_regions(call, intervals)]
            truth_subset = [call for call in truth if self._in_regions(call, intervals)]
            matches, fps, fns = self._match_calls(pred_subset, truth_subset)
            results[label] = self.calculate_metrics(len(matches), len(fps), len(fns))
        return results

    def stratify_by_svtype(
        self,
        predictions: Sequence[SVCandidate],
        truth_set: Optional[Sequence[SVCandidate]] = None,
    ) -> Dict[str, Dict[str, float]]:
        truth = list(truth_set if truth_set is not None else self.truth_set)
        results: Dict[str, Dict[str, float]] = {}
        for svtype in SVType:
            pred_subset = [call for call in predictions if call.svtype == svtype]
            truth_subset = [call for call in truth if call.svtype == svtype]
            matches, fps, fns = self._match_calls(pred_subset, truth_subset)
            results[svtype.value] = self.calculate_metrics(len(matches), len(fps), len(fns))
        return results

    def stratify_by_size(
        self,
        predictions: Sequence[SVCandidate],
        truth_set: Optional[Sequence[SVCandidate]] = None,
    ) -> Dict[str, Dict[str, float]]:
        truth = list(truth_set if truth_set is not None else self.truth_set)
        bins = {
            "50-100bp": (50, 100),
            "100-1000bp": (100, 1000),
            "1-10kb": (1000, 10_000),
            ">10kb": (10_000, float("inf")),
        }
        results: Dict[str, Dict[str, float]] = {}
        for label, (lower, upper) in bins.items():
            pred_subset = [call for call in predictions if lower <= call.size < upper]
            truth_subset = [call for call in truth if lower <= call.size < upper]
            matches, fps, fns = self._match_calls(pred_subset, truth_subset)
            results[label] = self.calculate_metrics(len(matches), len(fps), len(fns))
        return results

    def generate_report(self, evaluation: Mapping[str, Any]) -> str:
        """Generate a compact text report summarizing benchmark results."""

        lines = [
            "DeepSV-LR Benchmark Report",
            f"Precision: {evaluation['precision']:.3f}",
            f"Recall:    {evaluation['recall']:.3f}",
            f"F1 score:  {evaluation['f1']:.3f}",
            "",
            "Stratified by SV type:",
        ]
        for svtype, metrics in evaluation.get("stratified_by_svtype", {}).items():
            lines.append(f"  - {svtype}: P={metrics['precision']:.3f} R={metrics['recall']:.3f} F1={metrics['f1']:.3f}")
        lines.append("")
        lines.append("Stratified by size:")
        for size_bin, metrics in evaluation.get("stratified_by_size", {}).items():
            lines.append(f"  - {size_bin}: P={metrics['precision']:.3f} R={metrics['recall']:.3f} F1={metrics['f1']:.3f}")
        return "\n".join(lines)

    def _match_calls(
        self,
        predictions: Sequence[SVCandidate],
        truth: Sequence[SVCandidate],
    ) -> Tuple[List[Tuple[SVCandidate, SVCandidate]], List[SVCandidate], List[SVCandidate]]:
        used_truth: set[int] = set()
        matches: List[Tuple[SVCandidate, SVCandidate]] = []
        false_positives: List[SVCandidate] = []
        for prediction in predictions:
            best_index: Optional[int] = None
            best_score = 0.0
            for index, truth_call in enumerate(truth):
                if index in used_truth:
                    continue
                if truth_call.svtype != prediction.svtype or truth_call.chrom != prediction.chrom:
                    continue
                overlap = reciprocal_overlap(prediction, truth_call)
                size_ratio = min(prediction.size, truth_call.size) / max(prediction.size, truth_call.size, 1)
                if overlap >= self.overlap_threshold and size_ratio >= self.size_similarity and overlap > best_score:
                    best_score = overlap
                    best_index = index
            if best_index is None:
                false_positives.append(prediction)
            else:
                used_truth.add(best_index)
                matches.append((truth[best_index], prediction))
        false_negatives = [call for index, call in enumerate(truth) if index not in used_truth]
        return matches, false_positives, false_negatives

    @staticmethod
    def _record_to_candidate(record: Mapping[str, Any]) -> SVCandidate:
        chrom = str(record["chrom"])
        start = int(record["start"])
        end = int(record["end"])
        svtype = SVType(str(record["svtype"]))
        return SVCandidate(
            chrom=chrom,
            start=start,
            end=end,
            svtype=svtype,
            size=max(end - start, 1),
            left_breakpoint=Breakpoint(chrom, start),
            right_breakpoint=Breakpoint(chrom, end),
            info={"truth_id": record.get("id")},
        )

    @staticmethod
    def _in_regions(call: SVCandidate, regions: Sequence[Tuple[str, int, int]]) -> bool:
        return any(call.chrom == chrom and call.start < end and call.end > start for chrom, start, end in regions)


__all__ = ["GIABBenchmark"]
