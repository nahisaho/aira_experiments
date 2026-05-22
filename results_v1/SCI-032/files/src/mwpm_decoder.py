"""
MWPM (Minimum Weight Perfect Matching) decoder using PyMatching.
Implements threshold analysis and provides decoded logical error rates.
"""

import numpy as np
import stim
import pymatching
from typing import Tuple


def build_detector_error_model(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    """Extract the detector error model from a Stim circuit."""
    return circuit.detector_error_model(decompose_errors=True)


def build_matching_graph(circuit: stim.Circuit) -> pymatching.Matching:
    """
    Build a PyMatching graph from a Stim circuit's detector error model.
    """
    dem = build_detector_error_model(circuit)
    matching = pymatching.Matching.from_detector_error_model(dem)
    return matching


def sample_and_decode_mwpm(
    circuit: stim.Circuit,
    num_shots: int,
    seed: int = 42,
) -> Tuple[int, int]:
    """
    Sample circuit, decode with MWPM, return (num_errors, num_shots).
    
    Returns:
        (logical_errors, total_shots)
    """
    sampler = circuit.compile_detector_sampler(seed=seed)
    detection_events, observable_flips = sampler.sample(
        num_shots, separate_observables=True
    )

    matching = build_matching_graph(circuit)
    predictions = matching.decode_batch(detection_events)

    num_errors = int(np.sum(predictions != observable_flips))
    return num_errors, num_shots


def logical_error_rate_mwpm(
    circuit: stim.Circuit,
    num_shots: int,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Compute logical error rate per round with Wilson confidence interval.
    
    Returns:
        (logical_error_rate_per_round, std_error)
    """
    errors, shots = sample_and_decode_mwpm(circuit, num_shots, seed)
    p_logical = errors / shots

    # Wilson confidence interval (95%)
    z = 1.96
    n = shots
    p_hat = p_logical
    denominator = 1 + z**2 / n
    center = (p_hat + z**2 / (2 * n)) / denominator
    half_width = (z * np.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2))) / denominator
    ci_low = max(0, center - half_width)
    ci_high = min(1, center + half_width)

    # Estimate rounds from circuit
    rounds = _estimate_rounds(circuit)
    p_per_round = 1 - (1 - p_logical) ** (1.0 / max(rounds, 1))

    return p_per_round, half_width / max(rounds, 1)


def _estimate_rounds(circuit: stim.Circuit) -> int:
    """Estimate number of rounds from circuit by counting REPEAT blocks."""
    circuit_str = str(circuit)
    rounds = 1
    for line in circuit_str.split("\n"):
        if line.strip().startswith("REPEAT"):
            try:
                r = int(line.strip().split()[1])
                rounds = max(rounds, r)
            except (IndexError, ValueError):
                pass
    return rounds


def sweep_error_rates_mwpm(
    distance: int,
    rounds: int,
    error_rates: np.ndarray,
    num_shots: int = 10_000,
    seed: int = 42,
) -> np.ndarray:
    """
    Sweep physical error rates and return logical error rates for MWPM.
    
    Returns:
        Array of logical error rates per round.
    """
    import sys
    sys.path.insert(0, '/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/src')
    from noise_models import build_depolarizing_circuit

    logical_rates = []
    for p in error_rates:
        if p < 1e-6:
            logical_rates.append(0.0)
            continue
        circuit = build_depolarizing_circuit(distance, rounds, p)
        p_l, _ = logical_error_rate_mwpm(circuit, num_shots, seed=seed)
        logical_rates.append(p_l)
    return np.array(logical_rates)
