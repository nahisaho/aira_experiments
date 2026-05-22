"""Complex-structural-variant reconstruction for DeepSV-LR."""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any, DefaultDict, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

try:
    from .sv_detector import Breakpoint, SVCandidate, SVType
except ImportError:  # pragma: no cover - fallback for flat module execution
    from sv_detector import Breakpoint, SVCandidate, SVType


@dataclass
class GraphEdge:
    source: str
    target: str
    sv: SVCandidate


class BreakpointGraph:
    """Graph abstraction used to reconstruct complex SV patterns."""

    def __init__(self) -> None:
        self.adjacency: DefaultDict[str, List[GraphEdge]] = defaultdict(list)

    def build_from_calls(self, calls: Sequence[SVCandidate]) -> "BreakpointGraph":
        self.adjacency.clear()
        for call in calls:
            source = self._node_id(call.left_breakpoint.chrom, call.left_breakpoint.position, call.left_breakpoint.orientation)
            target = self._node_id(call.right_breakpoint.chrom, call.right_breakpoint.position, call.right_breakpoint.orientation)
            edge = GraphEdge(source=source, target=target, sv=call)
            reverse = GraphEdge(source=target, target=source, sv=call)
            self.adjacency[source].append(edge)
            self.adjacency[target].append(reverse)
        return self

    def find_connected_components(self) -> List[Set[str]]:
        components: List[Set[str]] = []
        seen: Set[str] = set()
        for node in self.adjacency:
            if node in seen:
                continue
            component: Set[str] = set()
            queue: deque[str] = deque([node])
            seen.add(node)
            while queue:
                current = queue.popleft()
                component.add(current)
                for edge in self.adjacency[current]:
                    if edge.target not in seen:
                        seen.add(edge.target)
                        queue.append(edge.target)
            components.append(component)
        return components

    def detect_cycles(self) -> List[List[str]]:
        cycles: List[List[str]] = []
        path: List[str] = []
        visited: Set[str] = set()

        def dfs(node: str, parent: Optional[str]) -> None:
            visited.add(node)
            path.append(node)
            for edge in self.adjacency[node]:
                if edge.target == parent:
                    continue
                if edge.target in path:
                    cycle_start = path.index(edge.target)
                    cycles.append(path[cycle_start:] + [edge.target])
                elif edge.target not in visited:
                    dfs(edge.target, node)
            path.pop()

        for node in self.adjacency:
            if node not in visited:
                dfs(node, None)
        return self._deduplicate_cycles(cycles)

    def detect_complex_rearrangements(self) -> List[Dict[str, Any]]:
        patterns: List[Dict[str, Any]] = []
        for component in self.find_connected_components():
            edges = [edge for node in component for edge in self.adjacency[node] if edge.source in component]
            svtypes = {edge.sv.svtype for edge in edges}
            if len(component) >= 4 or len(svtypes) >= 3:
                patterns.append(
                    {
                        "nodes": sorted(component),
                        "edge_count": len(edges) // 2,
                        "svtypes": sorted(svtype.value for svtype in svtypes),
                    }
                )
        return patterns

    @staticmethod
    def _node_id(chrom: str, position: int, orientation: str) -> str:
        return f"{chrom}:{position}:{orientation}"

    @staticmethod
    def _deduplicate_cycles(cycles: Sequence[Sequence[str]]) -> List[List[str]]:
        unique: Set[Tuple[str, ...]] = set()
        deduplicated: List[List[str]] = []
        for cycle in cycles:
            if len(cycle) < 3:
                continue
            canonical = tuple(sorted(set(cycle)))
            if canonical in unique:
                continue
            unique.add(canonical)
            deduplicated.append(list(cycle))
        return deduplicated


class ChromothripsisDector:
    """Detect chromothripsis-like signatures from clustered SV calls."""

    def __init__(self, breakpoint_window: int = 10_000_000) -> None:
        self.breakpoint_window = breakpoint_window

    def call(
        self,
        calls: Sequence[SVCandidate],
        copy_number_segments: Optional[Mapping[str, Sequence[Mapping[str, Any]]]] = None,
    ) -> List[Dict[str, Any]]:
        grouped: DefaultDict[str, List[SVCandidate]] = defaultdict(list)
        for call in calls:
            grouped[call.chrom].append(call)

        results: List[Dict[str, Any]] = []
        for chrom, chrom_calls in grouped.items():
            chrom_calls.sort(key=lambda call: call.start)
            cn_segments = list(copy_number_segments.get(chrom, [])) if copy_number_segments else []
            oscillation = self._oscillation_score(cn_segments)
            clustering = self._breakpoint_clustering(chrom_calls)
            orientation_entropy = self._orientation_entropy(chrom_calls)
            if oscillation >= 0.5 and clustering >= 0.5 and orientation_entropy >= 0.8:
                results.append(
                    {
                        "chrom": chrom,
                        "signature": "chromothripsis",
                        "oscillation_score": oscillation,
                        "breakpoint_clustering": clustering,
                        "orientation_entropy": orientation_entropy,
                        "event_count": len(chrom_calls),
                    }
                )
        return results

    @staticmethod
    def _oscillation_score(segments: Sequence[Mapping[str, Any]]) -> float:
        states = [int(round(segment.get("copy_number", 2))) for segment in segments]
        if len(states) < 4:
            return 0.0
        transitions = sum(1 for left, right in zip(states, states[1:]) if left != right)
        unique_states = len(set(states))
        return min(transitions / max(len(states) - 1, 1), 1.0) * (1.0 if unique_states <= 3 else 0.5)

    def _breakpoint_clustering(self, calls: Sequence[SVCandidate]) -> float:
        if len(calls) < 4:
            return 0.0
        positions = np.asarray([call.start for call in calls], dtype=np.float64)
        distances = np.diff(np.sort(positions))
        if distances.size == 0:
            return 0.0
        clustered = np.mean(distances < self.breakpoint_window / max(len(calls), 1))
        return float(clustered)

    @staticmethod
    def _orientation_entropy(calls: Sequence[SVCandidate]) -> float:
        orientations = [f"{call.left_breakpoint.orientation}{call.right_breakpoint.orientation}" for call in calls]
        counts = np.asarray(list({orientation: orientations.count(orientation) for orientation in set(orientations)}.values()), dtype=np.float64)
        if counts.size == 0:
            return 0.0
        probabilities = counts / counts.sum()
        entropy = -np.sum(probabilities * np.log2(probabilities + 1e-12))
        max_entropy = np.log2(max(len(probabilities), 1))
        return float(entropy / max(max_entropy, 1e-6))


ChromothripsisDetector = ChromothripsisDector


class EcDNADetector:
    """Detect extrachromosomal DNA signatures from SV cycles and amplification."""

    def __init__(self, amplification_threshold: float = 4.0) -> None:
        self.amplification_threshold = amplification_threshold

    def call(
        self,
        calls: Sequence[SVCandidate],
        coverage: Optional[Mapping[str, Sequence[float]]] = None,
    ) -> List[Dict[str, Any]]:
        graph = BreakpointGraph().build_from_calls(calls)
        cycles = graph.detect_cycles()
        results: List[Dict[str, Any]] = []
        for cycle in cycles:
            cycle_calls = self._calls_for_cycle(calls, cycle)
            amplification = self._amplification_score(cycle_calls, coverage)
            circular_read_support = self._circular_read_support(cycle_calls)
            if amplification >= self.amplification_threshold and circular_read_support >= 1:
                results.append(
                    {
                        "signature": "ecDNA",
                        "cycle_nodes": cycle,
                        "cycle_size": len(cycle_calls),
                        "amplification": amplification,
                        "circular_read_support": circular_read_support,
                    }
                )
        return results

    @staticmethod
    def _calls_for_cycle(calls: Sequence[SVCandidate], cycle: Sequence[str]) -> List[SVCandidate]:
        cycle_nodes = set(cycle)
        return [
            call
            for call in calls
            if f"{call.left_breakpoint.chrom}:{call.left_breakpoint.position}:{call.left_breakpoint.orientation}" in cycle_nodes
            or f"{call.right_breakpoint.chrom}:{call.right_breakpoint.position}:{call.right_breakpoint.orientation}" in cycle_nodes
        ]

    def _amplification_score(
        self,
        calls: Sequence[SVCandidate],
        coverage: Optional[Mapping[str, Sequence[float]]],
    ) -> float:
        cn_values = [call.copy_number for call in calls if call.copy_number is not None]
        if cn_values:
            return float(np.mean(cn_values))
        if coverage is None or not calls:
            return 0.0
        chrom = calls[0].chrom
        bins = np.asarray(coverage.get(chrom, []), dtype=np.float64)
        if bins.size == 0:
            return 0.0
        baseline = np.median(bins)
        local = bins[max(calls[0].start // 1000, 0) : max(calls[-1].end // 1000, 1)]
        return float(np.mean(local) / max(baseline, 1e-6) * 2.0) if local.size else 0.0

    @staticmethod
    def _circular_read_support(calls: Sequence[SVCandidate]) -> int:
        support = 0
        for call in calls:
            for evidence in call.evidence:
                support += int(evidence.metadata.get("circular_reads", 0))
        return support


__all__ = [
    "BreakpointGraph",
    "ChromothripsisDector",
    "ChromothripsisDetector",
    "EcDNADetector",
    "GraphEdge",
]
