"""Main orchestration layer for the DeepSV-LR long-read SV pipeline."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

try:
    from .benchmark import GIABBenchmark
    from .complex_sv import BreakpointGraph, ChromothripsisDector, EcDNADetector
    from .hybrid_integrator import HybridIntegrator, ShortReadEvidence
    from .repeat_handler import RepeatRegionHandler
    from .signal_basecaller import SignalBasecaller
    from .sv_detector import AssemblyCaller, EnsembleSVCaller, ReadDepthCaller, SVCandidate, SplitReadCaller
except ImportError:  # pragma: no cover - fallback for flat module execution
    from benchmark import GIABBenchmark
    from complex_sv import BreakpointGraph, ChromothripsisDector, EcDNADetector
    from hybrid_integrator import HybridIntegrator, ShortReadEvidence
    from repeat_handler import RepeatRegionHandler
    from signal_basecaller import SignalBasecaller
    from sv_detector import AssemblyCaller, EnsembleSVCaller, ReadDepthCaller, SVCandidate, SplitReadCaller


@dataclass(frozen=True)
class PipelineConfig:
    output_dir: str = "results"
    beam_width: int = 5
    min_sv_size: int = 30
    enable_repeat_analysis: bool = True
    enable_hybrid_integration: bool = True
    enable_benchmark: bool = False
    log_level: str = "INFO"
    progress_total: int = 7


@dataclass
class PipelineProgress:
    completed_steps: int = 0
    total_steps: int = 7
    current_stage: str = "initializing"

    @property
    def fraction(self) -> float:
        return self.completed_steps / max(self.total_steps, 1)


class DeepSVLRPipeline:
    """Coordinate DeepSV-LR basecalling, SV discovery and evaluation modules."""

    def __init__(self, config: Optional[PipelineConfig] = None) -> None:
        self.config = config or PipelineConfig()
        self.logger = logging.getLogger(self.__class__.__name__)
        self.progress = PipelineProgress(total_steps=self.config.progress_total)
        self.basecaller = SignalBasecaller()
        self.split_read_caller = SplitReadCaller(min_event_size=self.config.min_sv_size)
        self.read_depth_caller = ReadDepthCaller()
        self.assembly_caller = AssemblyCaller(min_event_size=self.config.min_sv_size)
        self.ensemble_caller = EnsembleSVCaller()
        self.repeat_handler = RepeatRegionHandler()
        self.hybrid_integrator = HybridIntegrator()
        self.chromothripsis_detector = ChromothripsisDector()
        self.ecdna_detector = EcDNADetector()
        self.benchmark = GIABBenchmark()
        self.configure()

    def configure(self, overrides: Optional[Mapping[str, Any]] = None) -> PipelineConfig:
        """Apply configuration overrides and initialize logging."""

        if overrides:
            self.config = replace(self.config, **dict(overrides))
            self.progress.total_steps = self.config.progress_total
        logging.basicConfig(level=getattr(logging, self.config.log_level.upper(), logging.INFO))
        return self.config

    def validate_inputs(self, inputs: Mapping[str, Any]) -> None:
        """Validate the minimum required input payload."""

        required = ["alignments", "depth_profile", "assembly_regions"]
        missing = [key for key in required if key not in inputs]
        if missing:
            raise ValueError(f"missing required inputs: {', '.join(missing)}")
        truth_path = inputs.get("truth_set_path")
        if truth_path is not None and not Path(truth_path).exists():
            raise FileNotFoundError(f"truth set not found: {truth_path}")

    def run(self, inputs: Mapping[str, Any]) -> Dict[str, Any]:
        """Execute the end-to-end DeepSV-LR pipeline."""

        self.validate_inputs(inputs)
        results: Dict[str, Any] = {"progress": self.progress}

        raw_signal = inputs.get("raw_signal")
        if raw_signal is not None:
            self._advance("basecalling")
            features = self.basecaller.preprocess_signal(raw_signal)
            logits = self.basecaller.forward_pass(features)
            results["basecall"] = self.basecaller.ctc_decode(logits, beam_width=self.config.beam_width)

        self._advance("split-read calling")
        split_calls = self.split_read_caller.call(inputs["alignments"])
        results["split_read_calls"] = split_calls

        self._advance("read-depth calling")
        depth_calls = self.read_depth_caller.call(inputs["depth_profile"])
        results["read_depth_calls"] = depth_calls

        self._advance("local assembly")
        assembly_calls = self.assembly_caller.call(inputs["assembly_regions"])
        results["assembly_calls"] = assembly_calls

        self._advance("ensembling")
        sv_calls = self.ensemble_caller.call(split_calls, depth_calls, assembly_calls)

        if self.config.enable_repeat_analysis and inputs.get("reference_sequences"):
            self._advance("repeat analysis")
            self._annotate_repeat_context(sv_calls, inputs["reference_sequences"])

        if self.config.enable_hybrid_integration and inputs.get("short_read_evidence"):
            self._advance("hybrid integration")
            evidence = [record if isinstance(record, ShortReadEvidence) else ShortReadEvidence(**record) for record in inputs["short_read_evidence"]]
            sv_calls = self.hybrid_integrator.call(sv_calls, evidence, inputs.get("population_frequency_panel"))

        graph = BreakpointGraph().build_from_calls(sv_calls)
        results["breakpoint_graph_components"] = graph.find_connected_components()
        results["complex_sv"] = {
            "chromothripsis": self.chromothripsis_detector.call(sv_calls, inputs.get("copy_number_segments")),
            "ecDNA": self.ecdna_detector.call(sv_calls, inputs.get("depth_profile")),
            "complex_components": graph.detect_complex_rearrangements(),
        }
        results["sv_calls"] = sv_calls

        if self.config.enable_benchmark and inputs.get("truth_set_path"):
            self._advance("benchmarking")
            truth = self.benchmark.load_truth_set(inputs["truth_set_path"])
            evaluation = self.benchmark.evaluate(sv_calls, truth)
            results["benchmark"] = evaluation
            results["benchmark_report"] = self.benchmark.generate_report(evaluation)

        Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)
        self.logger.info("DeepSV-LR pipeline completed with %d SV calls", len(results.get("sv_calls", [])))
        return results

    def _advance(self, stage: str) -> None:
        self.progress.completed_steps += 1
        self.progress.current_stage = stage
        self.logger.info("[%d/%d] %s", self.progress.completed_steps, self.progress.total_steps, stage)

    def _annotate_repeat_context(self, sv_calls: Sequence[SVCandidate], reference_sequences: Mapping[str, str]) -> None:
        for candidate in sv_calls:
            sequence = reference_sequences.get(candidate.chrom, "")
            if not sequence:
                continue
            start = max(candidate.start - 250, 0)
            end = min(candidate.end + 250, len(sequence))
            context = sequence[start:end]
            telomeres = self.repeat_handler.detect_telomere_repeats(context)
            tandem = self.repeat_handler.detect_tandem_expansion(context)
            filter_result = self.repeat_handler.kmer_frequency_filter(context, candidate.quality)
            candidate.info["telomere_hits"] = telomeres
            candidate.info["tandem_expansions"] = tandem
            candidate.info["repeat_entropy"] = filter_result["entropy"]
            candidate.quality = float(filter_result["adjusted_score"])


__all__ = ["DeepSVLRPipeline", "PipelineConfig", "PipelineProgress"]
