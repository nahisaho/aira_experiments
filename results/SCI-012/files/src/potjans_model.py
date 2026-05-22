"""
Potjans-Diesmann (2014) cortical microcircuit model.

Layers: L2/3, L4, L5, L6  x  cell types: E (excitatory), I (inhibitory)
→ 8 populations.

Reference: Potjans TC & Diesmann M (2014) Cerebral Cortex 24(3):785-806.

This implementation uses the Izhikevich neuron model as a computationally
efficient substitute for the LIF neurons in the original paper.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple
import time

from .gpu_architecture import IzhikevichPopulation, CSRConnectivity


# ---------------------------------------------------------------------------
# Model parameters (Potjans & Diesmann 2014, Table 1 & 2)
# ---------------------------------------------------------------------------

POPULATIONS = ["L23E", "L23I", "L4E", "L4I", "L5E", "L5I", "L6E", "L6I"]

# Neuron numbers (original 1mm² column)
N_FULL = {
    "L23E": 20683, "L23I": 5834,
    "L4E":  21915, "L4I":  5479,
    "L5E":  4850,  "L5I":  1065,
    "L6E":  14395, "L6I":  2948,
}

# Scale factor for faster simulation (use 0.1 for quick run, 1.0 for full)
DEFAULT_SCALE = 0.1

# Connection probabilities (8x8 matrix, rows=target, cols=source)
# Indices follow POPULATIONS order.
CONN_PROB = np.array([
    #  L23E    L23I    L4E     L4I     L5E     L5I     L6E     L6I
    [0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0.0,    0.0076, 0.0   ],  # L23E
    [0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0.0,    0.0042, 0.0   ],  # L23I
    [0.0077, 0.0059, 0.0497, 0.135,  0.0067, 0.0003, 0.0453, 0.0   ],  # L4E
    [0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0.0,    0.1057, 0.0   ],  # L4I
    [0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.0   ],  # L5E
    [0.0548, 0.0269, 0.0257, 0.0022, 0.06,   0.3158, 0.0086, 0.0   ],  # L5I
    [0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252],  # L6E
    [0.0364, 0.001,  0.0034, 0.0005, 0.0277, 0.008,  0.0658, 0.1443],  # L6I
], dtype=np.float32)

# Mean synaptic weight [mV] (excitatory=0.35, inhibitory=-0.35×g)
W_EXC = 0.35    # [a.u. compatible with Izhikevich current scaling]
W_INH = -0.35 * 4.0   # inhibitory ratio g=4

# Mean synaptic delay [ms]
DELAY_EXC = 1.5
DELAY_INH = 0.75

# Background Poisson input rate [spikes/s per neuron]
BG_RATE = 8.0

# Izhikevich parameters for excitatory (RS) and inhibitory (FS)
IZH_E = {"neuron_type": "RS"}
IZH_I = {"neuron_type": "FS"}


# ---------------------------------------------------------------------------
# Scaled model builder
# ---------------------------------------------------------------------------

def build_potjans_model(scale: float = DEFAULT_SCALE,
                        rng: Optional[np.random.Generator] = None,
                        backend: str = "auto") -> Dict:
    """
    Build the Potjans-Diesmann microcircuit.

    Args:
        scale:   fraction of full neuron numbers to simulate
        rng:     random number generator
        backend: computation backend

    Returns:
        dict with "populations", "connections", "N_scaled", "meta"
    """
    if rng is None:
        rng = np.random.default_rng(42)

    N_scaled = {k: max(1, int(v * scale)) for k, v in N_FULL.items()}

    # Create populations
    populations: Dict[str, IzhikevichPopulation] = {}
    for pop_name in POPULATIONS:
        is_exc = pop_name.endswith("E")
        ntype  = "RS" if is_exc else "FS"
        pop = IzhikevichPopulation(N_scaled[pop_name], neuron_type=ntype,
                                   backend=backend, rng=rng)
        populations[pop_name] = pop

    # Build CSR connectivity for each (target, source) pair
    connections: Dict[Tuple[str, str], CSRConnectivity] = {}
    for ti, tgt in enumerate(POPULATIONS):
        for si, src in enumerate(POPULATIONS):
            p_conn = float(CONN_PROB[ti, si])
            if p_conn == 0.0:
                continue
            n_pre  = N_scaled[src]
            n_post = N_scaled[tgt]
            is_exc_src = src.endswith("E")
            w = W_EXC if is_exc_src else W_INH
            # Small weight scatter
            w_std = abs(w) * 0.1
            conn = CSRConnectivity.random(
                n_pre, n_post, p_conn=p_conn,
                w_mean=abs(w), w_std=w_std, rng=rng
            )
            # Flip sign for inhibitory
            if not is_exc_src:
                conn.weights[:] *= -1.0
            connections[(tgt, src)] = conn

    return {
        "populations": populations,
        "connections":  connections,
        "N_scaled":     N_scaled,
        "meta": {
            "scale": scale,
            "total_neurons": sum(N_scaled.values()),
            "total_synapses": sum(c.n_syn for c in connections.values()),
        }
    }


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------

def run_potjans_simulation(model: Dict, T_ms: float = 500.0,
                           dt: float = 0.1,
                           rng: Optional[np.random.Generator] = None) -> Dict:
    """
    Simulate the Potjans-Diesmann microcircuit.

    Returns spike trains per population and summary statistics.
    """
    if rng is None:
        rng = np.random.default_rng(1)

    populations = model["populations"]
    connections  = model["connections"]
    N_scaled     = model["N_scaled"]
    T = int(T_ms / dt)
    p_bg = BG_RATE * dt * 1e-3   # Poisson background spike probability

    # Spike recorder: list of (time_ms, pop_index, neuron_index)
    spike_records: Dict[str, list] = {k: [] for k in POPULATIONS}

    from scipy.sparse import csr_matrix

    # Pre-build scipy CSR matrices for fast matrix-vector multiply
    csr_mats: Dict[Tuple[str,str], object] = {}
    for (tgt, src), conn in connections.items():
        mat = csr_matrix(
            (conn.weights, conn.indices,
             conn.indptr.astype(np.int32)),
            shape=(conn.n_post, conn.n_pre)
        )
        csr_mats[(tgt, src)] = mat

    # Calibrated background drive (tuned to produce ~2-5 Hz spontaneous rates
    # in the Izhikevich neuron model, accounting for its current-to-rate transfer)
    bg_drive = {
        "L23E": 6.5, "L23I": 5.8,
        "L4E":  7.2, "L4I":  6.0,
        "L5E":  6.0, "L5I":  5.5,
        "L6E":  5.5, "L6I":  5.0,
    }
    bg_noise_std = 2.5   # fluctuation SD

    t_start = time.perf_counter()
    for step in range(T):
        t_ms = step * dt
        # Compute I_syn for each population from all its sources
        I_dict: Dict[str, np.ndarray] = {}
        for tgt in POPULATIONS:
            N_t     = N_scaled[tgt]
            I_total = np.zeros(N_t, dtype=np.float32)
            for src in POPULATIONS:
                key = (tgt, src)
                if key not in csr_mats:
                    continue
                src_spk = populations[src].spikes.astype(np.float32)
                I_total += csr_mats[key].dot(src_spk)
            # Background drive + stochastic fluctuations
            noise = rng.standard_normal(N_t).astype(np.float32) * bg_noise_std
            I_total += bg_drive[tgt] + noise
            I_dict[tgt] = I_total

        # Step all populations
        for tgt in POPULATIONS:
            populations[tgt].step(dt, I_dict[tgt])

        # Record spikes
        for tgt in POPULATIONS:
            fired_idx = np.where(populations[tgt].spikes)[0]
            for idx in fired_idx:
                spike_records[tgt].append((t_ms, int(idx)))

    elapsed = time.perf_counter() - t_start

    # Compute per-layer firing rates
    rates = {}
    for pop_name in POPULATIONS:
        n_spikes = len(spike_records[pop_name])
        N = N_scaled[pop_name]
        rates[pop_name] = n_spikes / (N * T_ms * 1e-3)  # [Hz]

    return {
        "spike_records": spike_records,
        "rates":         rates,
        "elapsed_s":     elapsed,
        "T_ms":          T_ms,
        "N_scaled":      N_scaled,
        "meta":          model["meta"],
    }
