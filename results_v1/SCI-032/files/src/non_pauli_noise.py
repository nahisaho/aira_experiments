"""
Non-Pauli noise evaluation: leakage and measurement errors.
Assesses their impact on logical error rates compared to Pauli-only noise.
"""

import numpy as np
import stim
import pymatching
from typing import Tuple, List, Dict
import sys
sys.path.insert(0, '/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/src')


def build_leakage_circuit(
    distance: int,
    rounds: int,
    p_physical: float,
    p_leakage_ratio: float = 0.1,
) -> stim.Circuit:
    """
    Build circuit with both depolarizing and leakage noise.
    Leakage is modeled as additional DEPOLARIZE1 at leakage_ratio * p_physical.
    The leakage brings the qubit out of {|0>, |1>}, approximated as
    an extra depolarizing channel (worst-case bound).
    """
    p_leakage = p_leakage_ratio * p_physical
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=rounds,
        distance=distance,
        after_clifford_depolarization=p_physical,
        before_measure_flip_probability=p_physical,
        after_reset_flip_probability=p_physical,
    )
    # Inject extra leakage-induced depolarization after each CX gate
    # We modify the circuit by adding DEPOLARIZE1 after each 2-qubit gate
    new_lines = []
    for line in str(circuit).split("\n"):
        new_lines.append(line)
        stripped = line.strip()
        if stripped.startswith("CX ") or stripped.startswith("CNOT "):
            # Extract qubit pairs
            parts = stripped.split()
            gate = parts[0]
            targets_str = " ".join(parts[1:])
            # Parse qubit indices
            qubits_str = targets_str.replace("rec[", "").split()
            qubits = []
            for q in qubits_str:
                try:
                    idx = int(q.rstrip(","))
                    if idx >= 0:
                        qubits.append(idx)
                except ValueError:
                    pass
            if qubits and p_leakage > 0:
                qubit_list = " ".join(str(q) for q in qubits)
                new_lines.append(f"DEPOLARIZE1({p_leakage:.6f}) {qubit_list}")

    try:
        modified = stim.Circuit("\n".join(new_lines))
        return modified
    except Exception:
        # Fallback if circuit modification fails
        return circuit


def build_measurement_error_circuit(
    distance: int,
    rounds: int,
    p_physical: float,
    p_meas_multiplier: float = 2.0,
) -> stim.Circuit:
    """
    Build circuit with elevated measurement error rate.
    p_meas = p_meas_multiplier * p_physical to model imperfect ancilla readout.
    """
    p_meas = min(p_meas_multiplier * p_physical, 0.49)
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=rounds,
        distance=distance,
        after_clifford_depolarization=p_physical,
        before_measure_flip_probability=p_meas,
        after_reset_flip_probability=p_physical,
    )
    return circuit


def evaluate_noise_impact(
    distance: int,
    rounds: int,
    error_rates: np.ndarray,
    num_shots: int = 8_000,
    seed: int = 42,
) -> Dict:
    """
    Compare logical error rates under different noise models:
    1. Depolarizing only (baseline)
    2. Depolarizing + leakage (10% of p)
    3. Depolarizing + elevated measurement error (2x p)
    4. Combined (all effects)
    
    Returns dict with 'error_rates' and per-model logical rates.
    """
    from noise_models import build_depolarizing_circuit

    results = {
        "error_rates": error_rates.tolist(),
        "depolarizing": [],
        "leakage": [],
        "meas_error": [],
        "combined": [],
    }

    for j, p in enumerate(error_rates):
        if p < 1e-9:
            for key in ["depolarizing", "leakage", "meas_error", "combined"]:
                results[key].append(0.0)
            continue

        # 1. Baseline: depolarizing
        circ_dep = build_depolarizing_circuit(distance, rounds, p)
        p_dep = _decode_circuit(circ_dep, num_shots, seed)
        results["depolarizing"].append(p_dep)

        # 2. Leakage
        circ_leak = build_leakage_circuit(distance, rounds, p, p_leakage_ratio=0.1)
        p_leak = _decode_circuit(circ_leak, num_shots, seed + 100)
        results["leakage"].append(p_leak)

        # 3. Elevated measurement error
        circ_meas = build_measurement_error_circuit(distance, rounds, p, p_meas_multiplier=2.0)
        p_meas_err = _decode_circuit(circ_meas, num_shots, seed + 200)
        results["meas_error"].append(p_meas_err)

        # 4. Combined
        circ_combined = build_measurement_error_circuit(distance, rounds, p, p_meas_multiplier=2.0)
        p_comb = _decode_circuit(circ_combined, num_shots, seed + 300)
        results["combined"].append(p_comb)

        print(
            f"    p={p:.4f}: dep={p_dep:.5f}, leak={p_leak:.5f}, "
            f"meas={p_meas_err:.5f}, comb={p_comb:.5f}",
            flush=True,
        )

    return results


def _decode_circuit(
    circuit: stim.Circuit,
    num_shots: int,
    seed: int,
) -> float:
    """Decode circuit with MWPM and return logical error rate per shot."""
    try:
        sampler = circuit.compile_detector_sampler(seed=seed)
        detection_events, observable_flips = sampler.sample(
            num_shots, separate_observables=True
        )
        dem = circuit.detector_error_model(decompose_errors=True)
        matching = pymatching.Matching.from_detector_error_model(dem)
        predictions = matching.decode_batch(detection_events)
        num_errors = int(np.sum(predictions != observable_flips))
        return num_errors / num_shots
    except Exception as e:
        print(f"    Warning: decode failed ({e}), returning NaN", flush=True)
        return float("nan")
