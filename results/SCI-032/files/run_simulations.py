"""
Main simulation runner for surface code logical error rate estimation.
Executes all simulation tasks and saves results to JSON files.
"""

import numpy as np
import json
import time
import sys
import os
sys.path.insert(0, '/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/src')

from noise_models import NoiseParams, build_depolarizing_circuit, build_independent_noise_circuit
from mwpm_decoder import sample_and_decode_mwpm, logical_error_rate_mwpm
from union_find_decoder import UnionFindDecoder, sample_and_decode_uf
from threshold_analysis import (
    run_threshold_analysis, fit_threshold,
    compute_circuit_level_threshold, compute_code_capacity_threshold
)
from non_pauli_noise import evaluate_noise_impact
from lattice_surgery import sweep_lattice_surgery, logical_cnot_error_rate

import stim
import pymatching

RESULTS_DIR = "/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/results"
LOGS_DIR = "/app/projects/9a7958af-1965-498d-ba8a-315793461ff6/workspace/logs"

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

log_entries = []

def log(phase, event_type, skill_or_tool, handoff_in=None, handoff_out=None,
        files_written=None, status="ok"):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files_written or [],
        "status": status,
    }
    log_entries.append(entry)
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


def save_results(filename, data):
    path = f"{RESULTS_DIR}/{filename}"
    with open(path, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"  Saved: {path}", flush=True)
    return path


def run_task1_noise_models():
    """Task 1: Noise model verification - compare depolarizing vs biased noise."""
    print("\n=== Task 1: Noise Models ===", flush=True)
    log("task1", "run_started", "noise_models")

    distance = 5
    rounds = 5
    num_shots = 8_000
    error_rates = np.array([0.001, 0.003, 0.005, 0.008, 0.01, 0.015, 0.02])

    results = {
        "distance": distance,
        "rounds": rounds,
        "error_rates": error_rates.tolist(),
        "depolarizing": [],
        "amplitude_damping": [],
        "phase_damping": [],
        "combined": [],
    }

    for p in error_rates:
        print(f"  p={p:.4f}", flush=True)

        # Pure depolarizing
        circ = build_depolarizing_circuit(distance, rounds, p)
        errors, shots = sample_and_decode_mwpm(circ, num_shots, seed=42)
        p_dep = errors / shots
        results["depolarizing"].append(p_dep)

        # Amplitude damping (T1-dominated: more X errors)
        params_amp = NoiseParams(p_amplitude_damp=p * 0.8, p_depolarize=p * 0.2)
        circ_amp = build_independent_noise_circuit(distance, rounds, params_amp)
        errors_amp, _ = sample_and_decode_mwpm(circ_amp, num_shots, seed=43)
        results["amplitude_damping"].append(errors_amp / shots)

        # Phase damping (T2*-dominated: more Z errors)
        params_phase = NoiseParams(p_phase_damp=p * 0.8, p_depolarize=p * 0.2)
        circ_phase = build_independent_noise_circuit(distance, rounds, params_phase)
        errors_phase, _ = sample_and_decode_mwpm(circ_phase, num_shots, seed=44)
        results["phase_damping"].append(errors_phase / shots)

        # Combined
        params_comb = NoiseParams(
            p_depolarize=p * 0.4,
            p_amplitude_damp=p * 0.3,
            p_phase_damp=p * 0.3,
        )
        circ_comb = build_independent_noise_circuit(distance, rounds, params_comb)
        errors_comb, _ = sample_and_decode_mwpm(circ_comb, num_shots, seed=45)
        results["combined"].append(errors_comb / shots)

    path = save_results("noise_model_comparison.json", results)
    log("task1", "run_completed", "noise_models", files_written=[path])
    return results


def run_task2_mwpm_threshold():
    """Task 2 & 3: MWPM decoder + threshold analysis."""
    print("\n=== Task 2&3: MWPM Threshold Analysis ===", flush=True)
    log("task2", "run_started", "mwpm_decoder+threshold_analysis")

    distances = [3, 5, 7, 9]
    error_rates = np.array([0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
                            0.008, 0.009, 0.01, 0.012, 0.015, 0.02])
    num_shots = 10_000

    results = run_threshold_analysis(
        distances, error_rates,
        rounds_per_distance={d: d for d in distances},
        num_shots=num_shots,
        seed=42,
    )

    # Fit threshold
    logical_rates_by_d = {
        distances[i]: np.array(results["logical_rates"][i])
        for i in range(len(distances))
    }
    threshold, threshold_unc = fit_threshold(error_rates, logical_rates_by_d)
    results["threshold_estimate"] = threshold
    results["threshold_uncertainty"] = threshold_unc
    results["theoretical_circuit_level_threshold"] = compute_circuit_level_threshold()
    results["theoretical_code_capacity_threshold"] = compute_code_capacity_threshold()

    print(f"\n  Fitted threshold: p_th = {threshold:.4f} ± {threshold_unc:.4f}", flush=True)
    print(f"  Theoretical circuit-level: ~{compute_circuit_level_threshold():.4f}", flush=True)

    path = save_results("mwpm_threshold.json", results)
    log("task2", "run_completed", "mwpm_decoder", files_written=[path])
    return results


def run_task4_decoder_comparison():
    """Task 4: MWPM vs Union-Find decoder comparison."""
    print("\n=== Task 4: Decoder Comparison (MWPM vs Union-Find) ===", flush=True)
    log("task4", "run_started", "decoder_comparison")

    distance = 5
    rounds = 5
    num_shots = 5_000
    error_rates = np.array([0.001, 0.003, 0.005, 0.007, 0.01, 0.015])

    results = {
        "distance": distance,
        "rounds": rounds,
        "error_rates": error_rates.tolist(),
        "mwpm_logical_rates": [],
        "mwpm_ci": [],
        "uf_logical_rates": [],
        "uf_ci": [],
        "mwpm_decode_time_ms": [],
        "uf_decode_time_ms": [],
    }

    for p in error_rates:
        print(f"  p={p:.4f}", flush=True)
        circ = build_depolarizing_circuit(distance, rounds, p)

        # MWPM timing
        sampler = circ.compile_detector_sampler(seed=42)
        detection_events, observable_flips = sampler.sample(
            num_shots, separate_observables=True
        )

        dem = circ.detector_error_model(decompose_errors=True)
        matching = pymatching.Matching.from_detector_error_model(dem)

        t0 = time.perf_counter()
        predictions_mwpm = matching.decode_batch(detection_events)
        t_mwpm = (time.perf_counter() - t0) * 1000 / num_shots  # ms per shot

        errors_mwpm = int(np.sum(predictions_mwpm != observable_flips))
        p_mwpm = errors_mwpm / num_shots
        ci_mwpm = 1.96 * np.sqrt(p_mwpm * (1 - p_mwpm) / num_shots)

        results["mwpm_logical_rates"].append(p_mwpm)
        results["mwpm_ci"].append(float(ci_mwpm))
        results["mwpm_decode_time_ms"].append(float(t_mwpm))

        # Union-Find timing
        uf_decoder = UnionFindDecoder(circ)
        t0 = time.perf_counter()
        predictions_uf = uf_decoder.decode_batch(detection_events[:500])  # smaller batch for timing
        t_uf = (time.perf_counter() - t0) * 1000 / 500

        errors_uf = int(np.sum(predictions_uf != observable_flips[:500]))
        p_uf = errors_uf / 500
        ci_uf = 1.96 * np.sqrt(p_uf * (1 - p_uf) / 500) if p_uf > 0 else 0.0

        results["uf_logical_rates"].append(p_uf)
        results["uf_ci"].append(float(ci_uf))
        results["uf_decode_time_ms"].append(float(t_uf))

        print(f"    MWPM: p_L={p_mwpm:.5f} ({t_mwpm*1000:.2f} μs/shot), "
              f"UF: p_L={p_uf:.5f} ({t_uf*1000:.2f} μs/shot)", flush=True)

    path = save_results("decoder_comparison.json", results)
    log("task4", "run_completed", "decoder_comparison", files_written=[path])
    return results


def run_task5_non_pauli():
    """Task 5: Non-Pauli noise (leakage + measurement errors)."""
    print("\n=== Task 5: Non-Pauli Noise Impact ===", flush=True)
    log("task5", "run_started", "non_pauli_noise")

    distance = 5
    rounds = 5
    error_rates = np.array([0.001, 0.003, 0.005, 0.008, 0.01, 0.015])
    num_shots = 8_000

    results = evaluate_noise_impact(distance, rounds, error_rates, num_shots, seed=42)

    path = save_results("non_pauli_noise.json", results)
    log("task5", "run_completed", "non_pauli_noise", files_written=[path])
    return results


def run_task6_lattice_surgery():
    """Task 6: Lattice surgery (logical CNOT) simulation."""
    print("\n=== Task 6: Lattice Surgery ===", flush=True)
    log("task6", "run_started", "lattice_surgery")

    distances = [3, 5, 7]
    error_rates = np.array([0.001, 0.003, 0.005, 0.008, 0.01])
    num_shots = 5_000

    results = sweep_lattice_surgery(distances, error_rates, num_shots, seed=42)

    # Add T-gate estimates
    t_gate_rates = {}
    for d in distances:
        t_rates = []
        for j, p in enumerate(error_rates):
            p_cnot = results["cnot_rates"][distances.index(d)][j]
            p_t = 35 * p_cnot ** 3
            t_rates.append(float(p_t))
        t_gate_rates[str(d)] = t_rates
    results["t_gate_rates"] = t_gate_rates

    path = save_results("lattice_surgery.json", results)
    log("task6", "run_completed", "lattice_surgery", files_written=[path])
    return results


def compute_decoder_scaling():
    """Compute decoder performance vs code distance."""
    print("\n=== Decoder Scaling Analysis ===", flush=True)

    distances = [3, 5, 7, 9]
    p_fixed = 0.005
    num_shots = 3_000
    results = {
        "distances": distances,
        "p_physical": p_fixed,
        "mwpm_logical_rates": [],
        "uf_logical_rates": [],
        "mwpm_time_ms": [],
        "uf_time_ms": [],
        "num_qubits": [],
    }

    for d in distances:
        print(f"  d={d}", flush=True)
        circ = build_depolarizing_circuit(d, d, p_fixed)
        n_qubits = 2 * d**2 - 2 * d + 1  # rotated surface code data qubits

        sampler = circ.compile_detector_sampler(seed=42)
        det_events, obs_flips = sampler.sample(num_shots, separate_observables=True)

        dem = circ.detector_error_model(decompose_errors=True)
        matching = pymatching.Matching.from_detector_error_model(dem)

        t0 = time.perf_counter()
        preds_mwpm = matching.decode_batch(det_events)
        t_mwpm = (time.perf_counter() - t0) * 1000 / num_shots

        errors_mwpm = int(np.sum(preds_mwpm != obs_flips))
        p_mwpm = errors_mwpm / num_shots

        uf = UnionFindDecoder(circ)
        batch = det_events[:500]
        t0 = time.perf_counter()
        preds_uf = uf.decode_batch(batch)
        t_uf = (time.perf_counter() - t0) * 1000 / 500

        errors_uf = int(np.sum(preds_uf != obs_flips[:500]))
        p_uf = errors_uf / 500

        results["mwpm_logical_rates"].append(float(p_mwpm))
        results["uf_logical_rates"].append(float(p_uf))
        results["mwpm_time_ms"].append(float(t_mwpm))
        results["uf_time_ms"].append(float(t_uf))
        results["num_qubits"].append(n_qubits)

        print(f"    n_qubits={n_qubits}, MWPM={p_mwpm:.5f} ({t_mwpm*1000:.1f}μs), "
              f"UF={p_uf:.5f} ({t_uf*1000:.1f}μs)", flush=True)

    path = save_results("decoder_scaling.json", results)
    return results


if __name__ == "__main__":
    np.random.seed(42)
    print("Surface Code Logical Error Rate Simulation Framework", flush=True)
    print("=" * 60, flush=True)

    log("main", "run_started", "main_simulation_runner")

    t_start = time.time()

    r1 = run_task1_noise_models()
    r2 = run_task2_mwpm_threshold()
    r4 = run_task4_decoder_comparison()
    r5 = run_task5_non_pauli()
    r6 = run_task6_lattice_surgery()
    r_scale = compute_decoder_scaling()

    total_time = time.time() - t_start
    print(f"\n=== All simulations complete in {total_time:.1f}s ===", flush=True)

    summary = {
        "total_runtime_seconds": total_time,
        "threshold_estimate": r2["threshold_estimate"],
        "threshold_uncertainty": r2["threshold_uncertainty"],
        "theoretical_threshold": r2["theoretical_circuit_level_threshold"],
        "tasks_completed": 6,
    }
    save_results("simulation_summary.json", summary)

    log("main", "run_completed", "main_simulation_runner",
        handoff_out=summary,
        files_written=[f"{RESULTS_DIR}/simulation_summary.json"],
        status="ok")

    print(f"\nThreshold: {r2['threshold_estimate']:.4f} ± {r2['threshold_uncertainty']:.4f}")
    print(f"Theoretical circuit-level: ~{r2['theoretical_circuit_level_threshold']:.4f}")
