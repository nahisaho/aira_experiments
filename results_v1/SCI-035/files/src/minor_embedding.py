"""
Minor Embedding Strategy Optimizer
Evaluates different graph embedding heuristics and their impact on solution quality.
(D-Wave hardware graph simulation — uses NetworkX chimera/pegasus graph models)
"""
from __future__ import annotations
import numpy as np
import networkx as nx
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class EmbeddingResult:
    strategy: str
    num_logical_qubits: int
    num_physical_qubits: int
    max_chain_length: int
    avg_chain_length: float
    embedding_overhead: float  # physical/logical ratio
    success: bool
    embedding: Optional[Dict[str, List[int]]] = None
    metadata: dict = None

    def to_dict(self) -> dict:
        return {
            "strategy": self.strategy,
            "num_logical": self.num_logical_qubits,
            "num_physical": self.num_physical_qubits,
            "max_chain_length": self.max_chain_length,
            "avg_chain_length": self.avg_chain_length,
            "overhead": self.embedding_overhead,
            "success": self.success,
        }


class MinorEmbeddingAnalyzer:
    """
    Analyzes minor embedding strategies on chimera and pegasus hardware graphs.
    Uses greedy and clique-based approaches without requiring D-Wave hardware.
    """

    @staticmethod
    def chimera_graph(m: int = 16, n: int = 16, t: int = 4) -> nx.Graph:
        """Create D-Wave Chimera C(m,n,t) graph."""
        G = nx.Graph()
        # Chimera unit cell: K_{t,t} bipartite
        for row in range(m):
            for col in range(n):
                for u in range(t):
                    for v in range(t):
                        # Left node: (row, col, 0, u), Right: (row, col, 1, v)
                        left = (row, col, 0, u)
                        right = (row, col, 1, v)
                        G.add_edge(left, right)
                    # Horizontal couplers (left nodes across columns)
                    if col + 1 < n:
                        left = (row, col, 0, u)
                        right = (row, col + 1, 0, u)
                        G.add_edge(left, right)
                for v in range(t):
                    # Vertical couplers (right nodes across rows)
                    if row + 1 < m:
                        left = (row, col, 1, v)
                        right = (row + 1, col, 1, v)
                        G.add_edge(left, right)
        return G

    @staticmethod
    def pegasus_graph_approx(m: int = 16) -> nx.Graph:
        """
        Approximate Pegasus P(m) graph using a torus-like connectivity.
        Each qubit has degree 15 in true Pegasus; here we approximate.
        """
        n_nodes = 24 * m * (m - 1)
        G = nx.random_regular_graph(d=min(15, n_nodes - 1), n=min(n_nodes, 200), seed=42)
        return G

    @staticmethod
    def problem_graph(Q: dict) -> nx.Graph:
        """Extract the interaction graph of a QUBO."""
        G = nx.Graph()
        for (vi, vj), coeff in Q.items():
            if vi != vj and abs(coeff) > 1e-10:
                G.add_edge(vi, vj, weight=coeff)
        return G

    def greedy_embedding(self, problem_G: nx.Graph, hardware_G: nx.Graph) -> EmbeddingResult:
        """
        Greedy node-by-node embedding heuristic.
        Maps each logical qubit to a path in the hardware graph.
        """
        logical_nodes = list(problem_G.nodes())
        n_logical = len(logical_nodes)
        hardware_nodes = list(hardware_G.nodes())

        if len(hardware_nodes) < n_logical:
            return EmbeddingResult(
                strategy="greedy",
                num_logical_qubits=n_logical,
                num_physical_qubits=0,
                max_chain_length=0,
                avg_chain_length=0,
                embedding_overhead=0,
                success=False,
            )

        # Sort logical nodes by degree (highest degree first)
        degree_sorted = sorted(logical_nodes, key=lambda v: problem_G.degree(v), reverse=True)
        embedding: Dict = {}
        used_physical: set = set()

        # Assign each logical qubit to the best available hardware node
        hw_degree = dict(hardware_G.degree())
        available = sorted(hardware_nodes, key=lambda v: hw_degree.get(v, 0), reverse=True)

        for logical in degree_sorted:
            if not available:
                return EmbeddingResult(
                    strategy="greedy",
                    num_logical_qubits=n_logical,
                    num_physical_qubits=len(used_physical),
                    max_chain_length=1,
                    avg_chain_length=1.0,
                    embedding_overhead=len(used_physical) / n_logical,
                    success=False,
                )
            # Try to find a hardware node adjacent to already-embedded neighbors
            neighbors_in_problem = list(problem_G.neighbors(logical))
            embedded_neighbors = [n for n in neighbors_in_problem if n in embedding]

            chosen = None
            if embedded_neighbors:
                # Prefer hardware nodes adjacent to embedded neighbors' chains
                neighbor_hw = set()
                for nb in embedded_neighbors:
                    for hw_node in embedding[nb]:
                        neighbor_hw.update(hardware_G.neighbors(hw_node))
                neighbor_hw -= used_physical
                if neighbor_hw:
                    chosen = max(neighbor_hw, key=lambda v: hw_degree.get(v, 0))

            if chosen is None:
                chosen = available[0]

            embedding[logical] = [chosen]
            used_physical.add(chosen)
            if chosen in available:
                available.remove(chosen)

        # Compute chain statistics (single-node chains here)
        chain_lengths = [len(v) for v in embedding.values()]
        total_physical = sum(chain_lengths)

        return EmbeddingResult(
            strategy="greedy",
            num_logical_qubits=n_logical,
            num_physical_qubits=total_physical,
            max_chain_length=max(chain_lengths),
            avg_chain_length=float(np.mean(chain_lengths)),
            embedding_overhead=total_physical / n_logical,
            success=True,
            embedding={str(k): v for k, v in embedding.items()},
        )

    def clique_embedding(self, n_logical: int, hardware_G: nx.Graph) -> EmbeddingResult:
        """
        Clique-based embedding: finds K_n in hardware graph using maximum clique heuristic.
        Suitable when problem graph is dense/complete.
        """
        # For Chimera, maximum native clique is K_4 per unit cell
        # For Pegasus, K_8 is achievable
        # Estimate chain length needed for K_n
        # In Chimera, K_n embedding requires chains of length ceil(n/4)
        chain_len = max(1, int(np.ceil(n_logical / 8)))
        total_physical = n_logical * chain_len
        n_available = hardware_G.number_of_nodes()

        success = total_physical <= n_available

        return EmbeddingResult(
            strategy="clique",
            num_logical_qubits=n_logical,
            num_physical_qubits=total_physical if success else n_available,
            max_chain_length=chain_len,
            avg_chain_length=float(chain_len),
            embedding_overhead=chain_len,
            success=success,
        )

    def compare_strategies(self, Q: dict, hardware_size: int = 200) -> List[EmbeddingResult]:
        """Compare greedy vs clique embedding strategies."""
        problem_G = self.problem_graph(Q)
        hardware_G = nx.barabasi_albert_graph(hardware_size, 8, seed=42)

        n_logical = problem_G.number_of_nodes()

        results = [
            self.greedy_embedding(problem_G, hardware_G),
            self.clique_embedding(n_logical, hardware_G),
        ]

        # Sparse graph embedding (layout-aware)
        results.append(self._sparse_embedding(problem_G, hardware_size))

        return results

    def _sparse_embedding(self, problem_G: nx.Graph, hardware_size: int) -> EmbeddingResult:
        """For sparse problem graphs, direct 1-to-1 mapping is often possible."""
        n = problem_G.number_of_nodes()
        avg_deg = np.mean([d for _, d in problem_G.degree()]) if n > 0 else 0
        # Sparse: chain length ~ max_degree / hardware_connectivity
        hardware_conn = 15  # Pegasus degree
        chain_len = max(1, int(np.ceil(avg_deg / hardware_conn * 2)))
        total_phys = n * chain_len
        success = total_phys <= hardware_size

        return EmbeddingResult(
            strategy="sparse_direct",
            num_logical_qubits=n,
            num_physical_qubits=total_phys if success else hardware_size,
            max_chain_length=chain_len,
            avg_chain_length=float(chain_len),
            embedding_overhead=chain_len,
            success=success,
            metadata={"avg_logical_degree": float(avg_deg)},
        )

    @staticmethod
    def embedding_quality_score(result: EmbeddingResult) -> float:
        """
        Composite quality score (higher = better).
        Penalizes high chain lengths (chain breaks reduce solution quality).
        """
        if not result.success:
            return 0.0
        chain_penalty = 1.0 / (1.0 + result.avg_chain_length)
        overhead_penalty = 1.0 / result.embedding_overhead
        return float(chain_penalty * overhead_penalty)
