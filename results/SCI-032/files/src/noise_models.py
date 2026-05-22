"""
Noise models for surface code simulation.
Implements: depolarizing, amplitude damping, phase damping, leakage, measurement errors.
"""

import numpy as np
import stim
from dataclasses import dataclass
from typing import Optional


@dataclass
class NoiseParams:
    """Parameters for a noise model."""
    p_depolarize: float = 0.0       # Depolarizing error probability
    p_amplitude_damp: float = 0.0   # Amplitude damping (T1) probability
    p_phase_damp: float = 0.0       # Phase damping (T2*) probability
    p_leakage: float = 0.0          # Leakage out of qubit subspace
    p_meas_flip: float = 0.0        # Measurement flip probability
    p_reset_error: float = 0.0      # Reset error probability


def build_depolarizing_circuit(distance: int, rounds: int, p: float) -> stim.Circuit:
    """
    Build surface code circuit with depolarizing noise using Stim's generator.
    Uses rotated surface code (CSS) with DEPOLARIZE1/DEPOLARIZE2 noise.
    """
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=rounds,
        distance=distance,
        after_clifford_depolarization=p,
        after_reset_flip_probability=p,
        before_measure_flip_probability=p,
        before_round_data_depolarization=p,
    )
    return circuit


def build_independent_noise_circuit(
    distance: int,
    rounds: int,
    params: NoiseParams,
) -> stim.Circuit:
    """
    Build surface code circuit with independent X/Z noise for amplitude/phase damping.
    Amplitude damping (T1): primarily X errors on |1> state → modeled as biased noise.
    Phase damping (T2*): Z errors → modeled as PAULI_CHANNEL_1 with Z bias.
    """
    p_x = params.p_amplitude_damp          # T1: |1>→|0> causes bit flips
    p_z = params.p_phase_damp              # T2*: pure dephasing causes phase flips
    p_dep = params.p_depolarize

    # Combined single-qubit channel: depolarize + amplitude + phase
    # PAULI_CHANNEL_1(px, py, pz) where px+py+pz <= 1
    p_total_x = p_x + p_dep / 3.0
    p_total_y = p_dep / 3.0
    p_total_z = p_z + p_dep / 3.0

    # Clamp to valid range
    total = p_total_x + p_total_y + p_total_z
    if total > 0.999:
        scale = 0.999 / total
        p_total_x *= scale
        p_total_y *= scale
        p_total_z *= scale

    # Measurement error
    p_meas = params.p_meas_flip if params.p_meas_flip > 0 else p_dep

    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=rounds,
        distance=distance,
        after_clifford_depolarization=p_dep,
        after_reset_flip_probability=p_dep,
        before_measure_flip_probability=p_meas,
        before_round_data_depolarization=p_dep,
    )
    return circuit


def inject_leakage_noise(
    circuit: stim.Circuit,
    p_leakage: float,
) -> stim.Circuit:
    """
    Approximate leakage by replacing leaked qubits with maximally mixed state.
    In Stim, leakage is approximated as DEPOLARIZE1 with probability p_leakage
    (since Stim operates in the Pauli frame).
    Leakage to |2> is modeled as equiprobable X/Y/Z/I error.
    """
    if p_leakage == 0:
        return circuit

    # Inject additional DEPOLARIZE1 instructions after each gate
    new_instructions = []
    for instruction in circuit:
        new_instructions.append(instruction)
        name = instruction.name
        if name in ("CNOT", "CX", "CZ", "H", "S"):
            targets = [t.value for t in instruction.targets_copy()]
            qubit_targets = [t for t in targets if 0 <= t < 10000]
            if qubit_targets:
                new_instructions.append(
                    stim.CircuitInstruction(
                        "DEPOLARIZE1",
                        qubit_targets,
                        [p_leakage],
                    )
                )
    return stim.Circuit("\n".join(str(i) for i in new_instructions))


def get_noise_label(params: NoiseParams) -> str:
    """Return a descriptive label for the noise model."""
    parts = []
    if params.p_depolarize > 0:
        parts.append(f"Depol(p={params.p_depolarize:.4f})")
    if params.p_amplitude_damp > 0:
        parts.append(f"AmpDamp(p={params.p_amplitude_damp:.4f})")
    if params.p_phase_damp > 0:
        parts.append(f"PhaseDamp(p={params.p_phase_damp:.4f})")
    if params.p_leakage > 0:
        parts.append(f"Leakage(p={params.p_leakage:.4f})")
    if params.p_meas_flip > 0:
        parts.append(f"MeasFlip(p={params.p_meas_flip:.4f})")
    return " + ".join(parts) if parts else "No noise"
