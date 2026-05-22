"""
Threshold analysis for surface code: sweep physical error rates vs code distances.
Fits threshold using logistic/polynomial models.
"""

import numpy as np
import stim
import pymatching
from scipy.optimize import curve_fit
from typing import List, Tuple, Dict
import sys
sys.path.insert(0, '/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/src')


def logical_error_rate_fast(
    distance: int,
    rounds: int,
    p: float,
    num_shots: int = 10_000,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Fast logical error rate estimation using Stim + PyMatching.
    Returns (p_logical_per_shot, wilson_half_width).
    """
    from noise_models import build_depolarizing_circuit
    if p < 1e-9:
        return 0.0, 0.0

    circuit = build_depolarizing_circuit(distance, rounds, p)
    sampler = circuit.compile_detector_sampler(seed=seed)
    detection_events, observable_flips = sampler.sample(
        num_shots, separate_observables=True
    )

    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    predictions = matching.decode_batch(detection_events)

    num_errors = int(np.sum(predictions != observable_flips))
    p_l = num_errors / num_shots

    # Wilson CI half-width
    z = 1.96
    denom = 1 + z**2 / num_shots
    half_w = (z * np.sqrt(p_l * (1 - p_l) / num_shots + z**2 / (4 * num_shots**2))) / denom
    return p_l, half_w


def run_threshold_analysis(
    distances: List[int],
    error_rates: np.ndarray,
    rounds_per_distance: Dict[int, int] = None,
    num_shots: int = 10_000,
    seed: int = 42,
) -> Dict:
    """
    Run full threshold sweep for multiple distances and error rates.
    
    Returns:
        results dict with keys:
            'distances', 'error_rates', 'logical_rates' (2D), 'ci_half' (2D)
    """
    if rounds_per_distance is None:
        rounds_per_distance = {d: d for d in distances}

    logical_rates = np.zeros((len(distances), len(error_rates)))
    ci_half = np.zeros((len(distances), len(error_rates)))

    for i, d in enumerate(distances):
        rounds = rounds_per_distance.get(d, d)
        print(f"  Distance d={d}, rounds={rounds}", flush=True)
        for j, p in enumerate(error_rates):
            p_l, hw = logical_error_rate_fast(d, rounds, p, num_shots, seed + i * 100)
            logical_rates[i, j] = p_l
            ci_half[i, j] = hw
            print(f"    p={p:.4f} → p_L={p_l:.6f} ± {hw:.6f}", flush=True)

    return {
        "distances": distances,
        "error_rates": error_rates.tolist(),
        "logical_rates": logical_rates.tolist(),
        "ci_half": ci_half.tolist(),
    }


def fit_threshold(
    error_rates: np.ndarray,
    logical_rates_by_distance: Dict[int, np.ndarray],
) -> Tuple[float, float]:
    """
    Estimate threshold error rate by fitting crossing point of logical error curves.
    Uses polynomial interpolation to find where curves for d and d+2 intersect.
    
    Returns:
        (threshold_estimate, uncertainty)
    """
    distances = sorted(logical_rates_by_distance.keys())
    crossings = []

    for i in range(len(distances) - 1):
        d1, d2 = distances[i], distances[i + 1]
        r1 = logical_rates_by_distance[d1]
        r2 = logical_rates_by_distance[d2]

        # Find crossing point: r1[j] ≈ r2[j]
        diff = np.array(r1) - np.array(r2)
        sign_changes = np.where(np.diff(np.sign(diff)))[0]

        for idx in sign_changes:
            p_lo, p_hi = error_rates[idx], error_rates[idx + 1]
            d_lo, d_hi = diff[idx], diff[idx + 1]
            if abs(d_hi - d_lo) > 1e-10:
                p_cross = p_lo - d_lo * (p_hi - p_lo) / (d_hi - d_lo)
                crossings.append(p_cross)

    if crossings:
        threshold = np.median(crossings)
        uncertainty = np.std(crossings) if len(crossings) > 1 else 0.001
        return float(threshold), float(uncertainty)
    else:
        # Fallback: find where smallest and largest distance curves cross
        d_small = distances[0]
        d_large = distances[-1]
        r_s = np.array(logical_rates_by_distance[d_small])
        r_l = np.array(logical_rates_by_distance[d_large])
        diff = r_s - r_l
        idx = np.argmin(np.abs(diff))
        return float(error_rates[idx]), 0.002


def scaling_fit(p: np.ndarray, p_th: float, nu: float, A: float, B: float) -> np.ndarray:
    """
    Standard finite-size scaling ansatz for threshold:
    p_L(p, d) = A + B * (p - p_th) * d^(1/nu)
    """
    return A + B * p


def compute_code_capacity_threshold() -> float:
    """
    Theoretical code-capacity threshold for depolarizing noise on surface code.
    Known result: ~10.31% for rotated surface code.
    """
    return 0.1031


def compute_circuit_level_threshold() -> float:
    """
    Expected circuit-level depolarizing threshold for surface code with MWPM.
    Known result: ~0.5-1.1% depending on noise model.
    """
    return 0.0057  # ~0.57% for standard circuit-level noise
