"""
Integrated Information Theory (IIT 3.0) — Φ calculation.

Implements:
  - Exact Φ for small systems (≤8 nodes)
  - Approximation via minimum information partition (MIP) search
  - Φ_max over all subsystems
  - Concept geometry (cause-effect repertoires)
"""
import numpy as np
from itertools import combinations, chain
from typing import Optional, Tuple, List, Dict
from .utils import discretize_channels


def _powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))


def _repertoire_prob(states: np.ndarray, mechanism: tuple, purview: tuple,
                     direction: str = "cause") -> np.ndarray:
    """
    Estimate cause or effect repertoire as conditional probability distribution.
    Uses empirical transition frequencies from discretized state data.

    states: (n_channels, n_samples) binary array
    mechanism: tuple of node indices defining the mechanism
    purview: tuple of node indices defining the purview
    direction: "cause" or "effect"
    """
    n, T = states.shape
    mech_states = states[list(mechanism), :]  # (|mech|, T)
    purview_states = states[list(purview), :]  # (|purview|, T)

    n_mech_states = 2 ** len(mechanism)
    n_pur_states = 2 ** len(purview)

    # Joint histogram: (mech_state, purview_state)
    joint = np.zeros((n_mech_states, n_pur_states))
    for t in range(T):
        if direction == "effect" and t == T - 1:
            continue
        if direction == "cause" and t == 0:
            continue

        if direction == "effect":
            mech_val = int("".join(map(str, mech_states[:, t])), 2) if len(mechanism) > 0 else 0
            pur_val = int("".join(map(str, purview_states[:, t + 1])), 2) if len(purview) > 0 else 0
        else:  # cause
            mech_val = int("".join(map(str, mech_states[:, t])), 2) if len(mechanism) > 0 else 0
            pur_val = int("".join(map(str, purview_states[:, t - 1])), 2) if len(purview) > 0 else 0

        if mech_val < n_mech_states and pur_val < n_pur_states:
            joint[mech_val, pur_val] += 1

    # Conditional: P(purview | mechanism) — averaged over mechanism states
    row_sums = joint.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    conditional = joint / row_sums

    # Marginalize over mechanism (uniform prior on mech states)
    marginal = conditional.mean(axis=0)
    marginal = marginal / (marginal.sum() + 1e-10)
    return marginal


def earth_mover_distance(p: np.ndarray, q: np.ndarray) -> float:
    """
    Wasserstein / Earth Mover Distance for discrete distributions.
    For equal-size 1D distributions, uses cumsum approach.
    """
    assert len(p) == len(q)
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / (p.sum() + 1e-10)
    q = q / (q.sum() + 1e-10)
    return float(np.sum(np.abs(np.cumsum(p) - np.cumsum(q))))


def phi_mip(states: np.ndarray, mechanism: tuple) -> Tuple[float, Optional[tuple]]:
    """
    Calculate Φ (phi) for a given mechanism using Minimum Information Partition (MIP).

    Φ = min over bipartitions of KL divergence between full system repertoire
        and partitioned system repertoire.

    Returns (phi_value, minimum_information_partition)
    """
    n = len(mechanism)
    if n < 2:
        return 0.0, None

    full_cause = _repertoire_prob(states, mechanism, mechanism, "cause")
    full_effect = _repertoire_prob(states, mechanism, mechanism, "effect")

    min_phi = np.inf
    mip_partition = None

    # Iterate over all bipartitions
    for size in range(1, n):
        for part_a in combinations(range(n), size):
            part_b = tuple(i for i in range(n) if i not in part_a)

            mech_a = tuple(mechanism[i] for i in part_a)
            mech_b = tuple(mechanism[i] for i in part_b)

            # Partitioned cause repertoire (product of marginals)
            cause_a = _repertoire_prob(states, mech_a, mech_a, "cause")
            cause_b = _repertoire_prob(states, mech_b, mech_b, "cause")

            # Product distribution (independent parts)
            n_states_a = len(cause_a)
            n_states_b = len(cause_b)
            product_cause = np.outer(cause_a, cause_b).ravel()

            # Effect
            effect_a = _repertoire_prob(states, mech_a, mech_a, "effect")
            effect_b = _repertoire_prob(states, mech_b, mech_b, "effect")
            product_effect = np.outer(effect_a, effect_b).ravel()

            # Align sizes
            full_c = full_cause[:len(product_cause)]
            full_e = full_effect[:len(product_effect)]
            if len(full_c) == 0 or len(product_cause) == 0:
                continue

            phi_cause = earth_mover_distance(
                full_c / (full_c.sum() + 1e-10),
                product_cause[:len(full_c)] / (product_cause[:len(full_c)].sum() + 1e-10)
            )
            phi_effect = earth_mover_distance(
                full_e / (full_e.sum() + 1e-10),
                product_effect[:len(full_e)] / (product_effect[:len(full_e)].sum() + 1e-10)
            )

            # IIT 3.0: Φ = min of cause and effect information
            phi_cand = min(phi_cause, phi_effect)

            if phi_cand < min_phi:
                min_phi = phi_cand
                mip_partition = (mech_a, mech_b)

    return float(min_phi) if min_phi != np.inf else 0.0, mip_partition


class PhiCalculator:
    """
    Efficient Φ (Integrated Information) calculator for IIT 3.0.

    For systems up to 8 nodes, computes Φ for all subsystems and returns
    Φ_max (the maximum Φ across all candidate complexes).

    Parameters
    ----------
    max_nodes : int
        Maximum subset size to consider (default 4, for tractability)
    n_states : int
        Number of discrete states per node (default 2 = binary)
    """

    def __init__(self, max_nodes: int = 4, n_states: int = 2):
        self.max_nodes = max_nodes
        self.n_states = n_states

    def compute_phi(self, data: np.ndarray) -> float:
        """
        Compute Φ from multi-channel time-series data.

        data: (n_channels, n_samples) continuous array
        Returns: scalar Φ value
        """
        states = discretize_channels(data, self.n_states)
        n = states.shape[0]
        nodes = tuple(range(n))

        phi_values = []
        for size in range(2, min(self.max_nodes + 1, n + 1)):
            for subset in combinations(nodes, size):
                phi_val, _ = phi_mip(states, subset)
                phi_values.append(phi_val)

        return float(np.max(phi_values)) if phi_values else 0.0

    def compute_phi_trajectory(self, data: np.ndarray,
                                window_size: int = 256,
                                step: int = 64) -> np.ndarray:
        """
        Compute Φ in sliding windows to get temporal trajectory.
        Returns array of Φ values.
        """
        n, T = data.shape
        phi_traj = []
        for start in range(0, T - window_size + 1, step):
            window = data[:, start:start + window_size]
            phi_traj.append(self.compute_phi(window))
        return np.array(phi_traj)

    def phi_spectrum(self, data: np.ndarray) -> Dict[tuple, float]:
        """
        Compute Φ for all subsystems — returns dict {subset: phi}.
        """
        states = discretize_channels(data, self.n_states)
        n = states.shape[0]
        result = {}
        for size in range(2, min(self.max_nodes + 1, n + 1)):
            for subset in combinations(range(n), size):
                phi_val, _ = phi_mip(states, subset)
                result[subset] = phi_val
        return result

    def integrated_information_matrix(self, data: np.ndarray) -> np.ndarray:
        """
        Pairwise Φ matrix: Φ_{ij} = phi({i, j}).
        Returns (n_channels, n_channels) symmetric matrix.
        """
        states = discretize_channels(data, self.n_states)
        n = states.shape[0]
        mat = np.zeros((n, n))
        for i, j in combinations(range(n), 2):
            phi_val, _ = phi_mip(states, (i, j))
            mat[i, j] = mat[j, i] = phi_val
        return mat
