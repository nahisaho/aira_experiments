#!/usr/bin/env python3
"""
Complex structural variant detection module.
Handles chromothripsis, extrachromosomal DNA (ecDNA), and other complex rearrangements.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Set
from collections import defaultdict
import math
import logging

logger = logging.getLogger(__name__)


@dataclass
class BreakpointEdge:
    """Edge in a breakpoint graph representing a rearrangement junction."""
    chrom1: str
    pos1: int
    strand1: str
    chrom2: str
    pos2: int
    strand2: str
    support: int
    edge_type: str  # reference, variant, or circular

    @property
    def is_interchromosomal(self) -> bool:
        return self.chrom1 != self.chrom2


@dataclass
class ChromothripsisEvent:
    """Detected chromothripsis event."""
    chromosomes: List[str]
    breakpoints: List[Tuple[str, int]]
    n_segments: int
    oscillating_cn: bool
    random_joins: bool
    score: float
    fragments: List[Dict]


@dataclass
class EcDNACandidate:
    """Candidate extrachromosomal DNA element."""
    segments: List[Dict]
    total_size: int
    circularity_score: float
    amplification_level: float
    breakpoint_edges: List[BreakpointEdge]
    is_circular: bool


class BreakpointGraph:
    """Graph-based representation of structural rearrangements."""

    def __init__(self):
        self.nodes: Dict[str, Dict] = {}
        self.edges: List[BreakpointEdge] = []
        self.adjacency: Dict[str, List[int]] = defaultdict(list)

    def add_breakpoint(self, chrom: str, pos: int, strand: str, info: Dict = None):
        node_id = f"{chrom}:{pos}:{strand}"
        self.nodes[node_id] = {
            "chrom": chrom, "pos": pos, "strand": strand,
            "info": info or {}
        }

    def add_edge(self, edge: BreakpointEdge):
        idx = len(self.edges)
        self.edges.append(edge)
        node1 = f"{edge.chrom1}:{edge.pos1}:{edge.strand1}"
        node2 = f"{edge.chrom2}:{edge.pos2}:{edge.strand2}"
        self.adjacency[node1].append(idx)
        self.adjacency[node2].append(idx)

    def find_cycles(self, max_cycle_length: int = 20) -> List[List[int]]:
        """Find cycles in the breakpoint graph (potential ecDNA)."""
        cycles = []
        visited_edges: Set[int] = set()

        for start_node in self.adjacency:
            cycle = self._dfs_cycle(
                start_node, start_node, [], visited_edges,
                max_cycle_length, set()
            )
            if cycle and len(cycle) >= 2:
                cycles.append(cycle)

        return cycles

    def _dfs_cycle(
        self, current: str, target: str, path: List[int],
        visited: Set[int], max_depth: int, visited_nodes: Set[str]
    ) -> Optional[List[int]]:
        if len(path) > 0 and current == target:
            return path.copy()
        if len(path) >= max_depth:
            return None

        visited_nodes.add(current)
        for edge_idx in self.adjacency.get(current, []):
            if edge_idx in visited:
                continue

            edge = self.edges[edge_idx]
            node1 = f"{edge.chrom1}:{edge.pos1}:{edge.strand1}"
            node2 = f"{edge.chrom2}:{edge.pos2}:{edge.strand2}"
            next_node = node2 if current == node1 else node1

            if next_node in visited_nodes and next_node != target:
                continue

            visited.add(edge_idx)
            result = self._dfs_cycle(
                next_node, target, path + [edge_idx],
                visited, max_depth, visited_nodes
            )
            if result:
                return result
            visited.discard(edge_idx)

        visited_nodes.discard(current)
        return None

    def get_connected_components(self) -> List[Set[str]]:
        visited: Set[str] = set()
        components = []

        for node in self.nodes:
            if node in visited:
                continue
            component: Set[str] = set()
            self._bfs_component(node, component, visited)
            components.append(component)

        return components

    def _bfs_component(
        self, start: str, component: Set[str], visited: Set[str]
    ):
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            for edge_idx in self.adjacency.get(node, []):
                edge = self.edges[edge_idx]
                node1 = f"{edge.chrom1}:{edge.pos1}:{edge.strand1}"
                node2 = f"{edge.chrom2}:{edge.pos2}:{edge.strand2}"
                next_node = node2 if node == node1 else node1
                if next_node not in visited:
                    queue.append(next_node)


class ChromothripsisDetector:
    """Detect chromothripsis events from breakpoint patterns."""

    CRITERIA = {
        "min_breakpoints": 10,
        "max_cn_states": 3,
        "min_random_join_ratio": 0.4,
        "clustering_distance": 50000,
    }

    def detect(
        self, breakpoints: List[Tuple[str, int, str]],
        copy_numbers: Dict[str, np.ndarray],
        window_size: int = 10000
    ) -> List[ChromothripsisEvent]:
        events = []
        chrom_breakpoints = defaultdict(list)
        for chrom, pos, strand in breakpoints:
            chrom_breakpoints[chrom].append((pos, strand))

        for chrom, bps in chrom_breakpoints.items():
            if len(bps) < self.CRITERIA["min_breakpoints"]:
                continue

            # Criterion 1: Clustering of breakpoints
            is_clustered = self._check_clustering(bps)
            if not is_clustered:
                continue

            # Criterion 2: Oscillating copy number states
            cn_data = copy_numbers.get(chrom, np.array([]))
            oscillating, n_states = self._check_cn_oscillation(cn_data)

            # Criterion 3: Random orientation of joins
            random_joins = self._check_random_joins(bps)

            score = self._compute_chromothripsis_score(
                len(bps), oscillating, n_states, random_joins
            )

            if score > 0.5:
                fragments = self._identify_fragments(bps, cn_data, window_size)
                events.append(ChromothripsisEvent(
                    chromosomes=[chrom],
                    breakpoints=[(chrom, pos) for pos, _ in bps],
                    n_segments=len(fragments),
                    oscillating_cn=oscillating,
                    random_joins=random_joins,
                    score=score,
                    fragments=fragments
                ))

        return events

    def _check_clustering(self, breakpoints: List[Tuple[int, str]]) -> bool:
        positions = sorted([p for p, _ in breakpoints])
        if len(positions) < 2:
            return False
        span = positions[-1] - positions[0]
        expected_span = len(positions) * self.CRITERIA["clustering_distance"]
        return span < expected_span

    def _check_cn_oscillation(
        self, cn_data: np.ndarray
    ) -> Tuple[bool, int]:
        if len(cn_data) == 0:
            return False, 0
        cn_rounded = np.round(cn_data)
        unique_states = len(np.unique(cn_rounded))
        changes = np.sum(np.abs(np.diff(cn_rounded)) > 0)
        oscillating = (
            unique_states <= self.CRITERIA["max_cn_states"] and
            changes > len(cn_data) * 0.3
        )
        return oscillating, unique_states

    def _check_random_joins(self, breakpoints: List[Tuple[int, str]]) -> bool:
        if len(breakpoints) < 4:
            return False
        strands = [s for _, s in breakpoints]
        plus_ratio = sum(1 for s in strands if s == '+') / len(strands)
        return 0.3 < plus_ratio < 0.7

    def _compute_chromothripsis_score(
        self, n_breakpoints: int, oscillating: bool,
        n_states: int, random_joins: bool
    ) -> float:
        score = 0.0
        score += min(n_breakpoints / 20.0, 0.3)
        if oscillating:
            score += 0.3
        if n_states <= 3:
            score += 0.1
        if random_joins:
            score += 0.3
        return min(score, 1.0)

    @staticmethod
    def _identify_fragments(
        breakpoints: List[Tuple[int, str]],
        cn_data: np.ndarray, window_size: int
    ) -> List[Dict]:
        positions = sorted([p for p, _ in breakpoints])
        fragments = []
        for i in range(len(positions) - 1):
            cn_val = 2.0
            if len(cn_data) > 0:
                bin_idx = positions[i] // window_size
                if bin_idx < len(cn_data):
                    cn_val = float(cn_data[bin_idx])
            fragments.append({
                "start": positions[i],
                "end": positions[i + 1],
                "size": positions[i + 1] - positions[i],
                "copy_number": cn_val,
            })
        return fragments


class EcDNADetector:
    """Detect extrachromosomal DNA (ecDNA) from long-read data."""

    def __init__(self, min_amplification: float = 4.0, min_size: int = 1000):
        self.min_amplification = min_amplification
        self.min_size = min_size

    def detect(
        self, graph: BreakpointGraph,
        copy_numbers: Dict[str, np.ndarray],
        window_size: int = 10000
    ) -> List[EcDNACandidate]:
        candidates = []
        cycles = graph.find_cycles()

        for cycle_edges in cycles:
            segments = self._extract_segments(graph, cycle_edges)
            total_size = sum(s.get("size", 0) for s in segments)

            if total_size < self.min_size:
                continue

            amp_level = self._estimate_amplification(
                segments, copy_numbers, window_size
            )
            circularity = self._compute_circularity_score(
                graph, cycle_edges
            )

            if amp_level >= self.min_amplification:
                candidates.append(EcDNACandidate(
                    segments=segments,
                    total_size=total_size,
                    circularity_score=circularity,
                    amplification_level=amp_level,
                    breakpoint_edges=[graph.edges[i] for i in cycle_edges],
                    is_circular=circularity > 0.7
                ))

        return candidates

    def _extract_segments(
        self, graph: BreakpointGraph, cycle_edges: List[int]
    ) -> List[Dict]:
        segments = []
        for edge_idx in cycle_edges:
            edge = graph.edges[edge_idx]
            if not edge.is_interchromosomal:
                segments.append({
                    "chrom": edge.chrom1,
                    "start": min(edge.pos1, edge.pos2),
                    "end": max(edge.pos1, edge.pos2),
                    "size": abs(edge.pos2 - edge.pos1),
                })
        return segments

    def _estimate_amplification(
        self, segments: List[Dict],
        copy_numbers: Dict[str, np.ndarray],
        window_size: int
    ) -> float:
        amp_values = []
        for seg in segments:
            chrom = seg.get("chrom", "")
            cn_data = copy_numbers.get(chrom, np.array([]))
            if len(cn_data) == 0:
                continue
            start_bin = seg["start"] // window_size
            end_bin = min(seg["end"] // window_size + 1, len(cn_data))
            if start_bin < end_bin:
                amp_values.append(float(np.mean(cn_data[start_bin:end_bin])))

        return float(np.mean(amp_values)) if amp_values else 1.0

    def _compute_circularity_score(
        self, graph: BreakpointGraph, cycle_edges: List[int]
    ) -> float:
        if not cycle_edges:
            return 0.0
        total_support = sum(graph.edges[i].support for i in cycle_edges)
        avg_support = total_support / len(cycle_edges)
        support_variance = np.var([
            graph.edges[i].support for i in cycle_edges
        ])
        uniformity = 1.0 / (1.0 + support_variance / (avg_support + 1))
        return min(uniformity * min(avg_support / 5.0, 1.0), 1.0)


class ComplexSVDetector:
    """Unified complex SV detection combining chromothripsis and ecDNA."""

    def __init__(self):
        self.chromothripsis_detector = ChromothripsisDetector()
        self.ecdna_detector = EcDNADetector()

    def analyze(
        self, breakpoints: List[Tuple[str, int, str]],
        copy_numbers: Dict[str, np.ndarray],
        graph: BreakpointGraph
    ) -> Dict:
        cth_events = self.chromothripsis_detector.detect(
            breakpoints, copy_numbers
        )
        ecdna_candidates = self.ecdna_detector.detect(
            graph, copy_numbers
        )

        return {
            "chromothripsis_events": cth_events,
            "ecdna_candidates": ecdna_candidates,
            "summary": {
                "n_chromothripsis": len(cth_events),
                "n_ecdna": len(ecdna_candidates),
                "affected_chromosomes": list(set(
                    chrom for e in cth_events for chrom in e.chromosomes
                )),
            }
        }
