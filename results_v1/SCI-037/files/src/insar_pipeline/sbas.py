"""Small Baseline Subset processing routines."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy import linalg



def _to_days(acquisition_times: np.ndarray) -> np.ndarray:
    arr = np.asarray(acquisition_times)
    if np.issubdtype(arr.dtype, np.datetime64):
        return (arr - arr[0]).astype("timedelta64[D]").astype(float)
    return arr.astype(float)



def select_pairs(
    acquisition_times: np.ndarray,
    perpendicular_baselines: np.ndarray,
    max_temporal_days: float = 365.0,
    max_perpendicular_baseline: float = 200.0,
) -> list[tuple[int, int]]:
    """Select SBAS interferometric pairs using temporal and spatial baseline thresholds."""
    times = _to_days(acquisition_times)
    baselines = np.asarray(perpendicular_baselines, dtype=float)
    if times.ndim != 1 or baselines.ndim != 1 or times.size != baselines.size:
        raise ValueError("acquisition_times and perpendicular_baselines must be matching 1-D arrays")
    pairs: list[tuple[int, int]] = []
    for start in range(times.size - 1):
        delta_t = times[start + 1 :] - times[start]
        delta_b = np.abs(baselines[start + 1 :] - baselines[start])
        valid = np.where((delta_t <= max_temporal_days) & (delta_b <= max_perpendicular_baseline))[0]
        for item in valid:
            pairs.append((start, start + 1 + int(item)))
    return pairs



def construct_network(num_acquisitions: int, pairs: list[tuple[int, int]]) -> dict[str, np.ndarray | list[list[int]]]:
    """Build the SBAS design matrix and determine connected acquisition components."""
    if num_acquisitions < 2:
        raise ValueError("At least two acquisitions are required")
    incidence = np.zeros((len(pairs), num_acquisitions - 1), dtype=float)
    adjacency = np.zeros((num_acquisitions, num_acquisitions), dtype=int)
    for row, (start, stop) in enumerate(pairs):
        if not (0 <= start < stop < num_acquisitions):
            raise ValueError("Invalid acquisition pair encountered")
        incidence[row, start:stop] = 1.0
        adjacency[start, stop] = 1
        adjacency[stop, start] = 1
    graph = csr_matrix(adjacency)
    count, labels = connected_components(graph, directed=False, connection="weak")
    components = [np.where(labels == label)[0].tolist() for label in range(count)]
    return {"incidence": incidence, "components": components, "labels": labels}



def svd_inversion(incidence: np.ndarray, interferograms: np.ndarray, damping: float = 1.0e-3) -> dict[str, np.ndarray | int]:
    """Invert SBAS interferograms using damped SVD and minimum-norm regularization."""
    incidence = np.asarray(incidence, dtype=float)
    interferograms = np.asarray(interferograms, dtype=float)
    if incidence.ndim != 2:
        raise ValueError("incidence must be a 2-D design matrix")
    if interferograms.ndim == 1:
        interferograms = interferograms[:, None]
    if interferograms.shape[0] != incidence.shape[0]:
        raise ValueError("interferograms rows must match incidence rows")
    u, singular, vt = linalg.svd(incidence, full_matrices=False)
    filt = singular / (singular**2 + damping**2)
    solution = vt.T @ ((filt[:, None]) * (u.T @ interferograms))
    rank = int(np.sum(singular > 1.0e-8))
    return {"increments": solution, "rank": rank, "singular_values": singular}



def estimate_deformation(
    pair_displacements: np.ndarray,
    pairs: list[tuple[int, int]],
    num_acquisitions: int,
    damping: float = 1.0e-3,
) -> dict[str, np.ndarray | dict[str, np.ndarray | list[list[int]]]]:
    """Estimate displacement time series from SBAS interferograms."""
    network = construct_network(num_acquisitions, pairs)
    inversion = svd_inversion(network["incidence"], pair_displacements, damping=damping)
    increments = inversion["increments"]
    deformation = np.vstack([np.zeros((1, increments.shape[1])), np.cumsum(increments, axis=0)])
    return {"deformation": deformation, "network": network, "inversion": inversion}
