#!/usr/bin/env python3
"""
Quantum Internet Network Protocol Simulation
=============================================
Simulates QKD protocols (BB84/E91), quantum repeaters, entanglement distillation,
quantum routing, decoherence/loss effects, and a Tokyo QKD network case study.

This is a NetSquid/SimulaQron-inspired discrete-event simulation framework.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
import networkx as nx
from scipy.optimize import minimize_scalar
import os
import json

np.random.seed(42)
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# ==============================================================================
# 1. BB84/E91 Finite-Key Analysis
# ==============================================================================

@dataclass
class QKDParameters:
    """Parameters for QKD protocol simulation."""
    total_signals: int = 10**9
    dark_count_rate: float = 1e-6
    detector_efficiency: float = 0.1
    misalignment_error: float = 0.01
    fiber_loss_db_per_km: float = 0.2
    security_parameter: float = 1e-10
    error_correction_efficiency: float = 1.16

def bb84_finite_key_rate(N: int, distance_km: float, params: QKDParameters) -> dict:
    """
    Compute BB84 secure key rate with finite-key corrections.
    Uses Shor-Preskill + finite-key analysis (Lim et al. 2014).
    
    R = (s_{Z,0} + s_{Z,1}(1 - h(phi_Z))) / N - leak_EC / N - correction_terms / N
    """
    eta_channel = 10 ** (-params.fiber_loss_db_per_km * distance_km / 10)
    eta_total = eta_channel * params.detector_efficiency
    
    # Gains and error rates
    Q_mu = 1 - (1 - params.dark_count_rate) * np.exp(-eta_total * 0.5)  # decoy mean photon 0.5
    e_mu = (params.dark_count_rate / 2 + params.misalignment_error * eta_total * 0.5) / Q_mu if Q_mu > 0 else 0.5
    e_mu = min(e_mu, 0.5)
    
    # Single-photon contributions (decoy-state estimation)
    Y1 = eta_total + params.dark_count_rate
    e1 = (params.dark_count_rate / 2 + params.misalignment_error * eta_total) / Y1 if Y1 > 0 else 0.5
    e1 = min(e1, 0.5)
    
    # Binary entropy
    def h(x):
        if x <= 0 or x >= 1:
            return 0
        return -x * np.log2(x) - (1 - x) * np.log2(1 - x)
    
    # Finite-key correction terms
    n_Z = N * Q_mu * 0.5  # signals in Z basis
    n_X = N * Q_mu * 0.5  # signals in X basis
    
    if n_Z < 10 or n_X < 10:
        return {"rate": 0, "Q_mu": Q_mu, "e_mu": e_mu, "n_sifted": 0}
    
    # Statistical fluctuation term
    delta = 7 * np.sqrt(np.log2(2 / params.security_parameter) / n_X) if n_X > 0 else 1
    phi_X = min(e1 + delta, 0.5)
    
    # Key rate formula
    s1 = max(0, Y1 * N * 0.5 * np.exp(-0.5) * 0.5)  # single photon events in Z
    
    leak_EC = n_Z * params.error_correction_efficiency * h(e_mu)
    
    # Finite-key penalty
    finite_correction = 2 * np.log2(1 / (2 * params.security_parameter))
    
    key_length = max(0, s1 * (1 - h(phi_X)) - leak_EC - finite_correction)
    rate = key_length / N if N > 0 else 0
    
    return {
        "rate": rate,
        "Q_mu": Q_mu,
        "e_mu": e_mu,
        "e1": e1,
        "n_sifted": n_Z,
        "key_length": key_length,
        "phi_X": phi_X
    }


def e91_finite_key_rate(N: int, distance_km: float, params: QKDParameters) -> dict:
    """
    Compute E91 key rate with finite-key analysis.
    Based on CHSH violation and entanglement-based QKD security.
    
    Key rate: R = 1 - h(e_Z) - h(e_X) with finite-key corrections.
    """
    eta_channel = 10 ** (-params.fiber_loss_db_per_km * distance_km / 10)
    eta_total = eta_channel * params.detector_efficiency
    eta_pair = eta_total ** 2  # both sides must detect
    
    # Detection rate
    Q = eta_pair + 2 * params.dark_count_rate * (1 - eta_pair) + params.dark_count_rate**2
    
    # QBER
    e_det = (params.dark_count_rate * (1 - eta_pair) + params.misalignment_error * eta_pair) / Q if Q > 0 else 0.5
    e_det = min(e_det, 0.5)
    
    def h(x):
        if x <= 0 or x >= 1:
            return 0
        return -x * np.log2(x) - (1 - x) * np.log2(1 - x)
    
    # CHSH parameter S = 2*sqrt(2) * (1 - 2*e)
    S = 2 * np.sqrt(2) * (1 - 2 * e_det)
    
    n_sifted = N * Q * (3/8)  # 3/8 fraction for key generation in E91
    
    if n_sifted < 10:
        return {"rate": 0, "Q": Q, "QBER": e_det, "S": S, "n_sifted": 0}
    
    # Finite-key correction
    delta = np.sqrt(np.log2(2 / params.security_parameter) / n_sifted) if n_sifted > 0 else 1
    e_upper = min(e_det + delta, 0.5)
    
    finite_correction = 2 * np.log2(1 / (2 * params.security_parameter))
    
    key_length = max(0, n_sifted * (1 - h(e_upper) - params.error_correction_efficiency * h(e_det)) - finite_correction)
    rate = key_length / N if N > 0 else 0
    
    return {
        "rate": rate,
        "Q": Q,
        "QBER": e_det,
        "S": S,
        "n_sifted": n_sifted,
        "key_length": key_length
    }


def run_qkd_analysis():
    """Run BB84 and E91 finite-key analysis across distances and key lengths."""
    params = QKDParameters()
    
    # Distance sweep
    distances = np.linspace(1, 200, 100)
    N_values = [10**6, 10**8, 10**10, 10**12]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # BB84 key rate vs distance
    for N in N_values:
        rates = [bb84_finite_key_rate(N, d, params)["rate"] for d in distances]
        axes[0].semilogy(distances, [max(r, 1e-20) for r in rates], 
                        label=f'N = $10^{{{int(np.log10(N))}}}$')
    
    axes[0].set_xlabel('Distance (km)')
    axes[0].set_ylabel('Secure Key Rate (bits/pulse)')
    axes[0].set_title('BB84 Finite-Key Rate vs Distance')
    axes[0].legend()
    axes[0].set_ylim(1e-12, 1e-1)
    axes[0].grid(True, alpha=0.3)
    
    # E91 key rate vs distance  
    for N in N_values:
        rates = [e91_finite_key_rate(N, d, params)["rate"] for d in distances]
        axes[1].semilogy(distances, [max(r, 1e-20) for r in rates],
                        label=f'N = $10^{{{int(np.log10(N))}}}$')
    
    axes[1].set_xlabel('Distance (km)')
    axes[1].set_ylabel('Secure Key Rate (bits/pulse)')
    axes[1].set_title('E91 Finite-Key Rate vs Distance')
    axes[1].legend()
    axes[1].set_ylim(1e-12, 1e-1)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/qkd_finite_key_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Finite-key length dependence at fixed distance
    fig, ax = plt.subplots(figsize=(8, 5))
    N_range = np.logspace(4, 14, 50).astype(int)
    fixed_distances = [10, 50, 100, 150]
    
    for d in fixed_distances:
        rates_bb84 = [bb84_finite_key_rate(N, d, params)["rate"] for N in N_range]
        ax.loglog(N_range, [max(r, 1e-20) for r in rates_bb84], label=f'd = {d} km')
    
    ax.set_xlabel('Total Signals N')
    ax.set_ylabel('Secure Key Rate (bits/pulse)')
    ax.set_title('BB84 Key Rate Convergence: Finite to Asymptotic')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/finite_key_convergence.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Collect key results
    results = {}
    for d in [10, 50, 100]:
        bb84_r = bb84_finite_key_rate(10**10, d, params)
        e91_r = e91_finite_key_rate(10**10, d, params)
        results[d] = {"bb84": bb84_r, "e91": e91_r}
    
    return results


# ==============================================================================
# 2. Quantum Repeater Memory Requirements
# ==============================================================================

@dataclass
class RepeaterNode:
    """Quantum repeater node parameters."""
    memory_coherence_time_ms: float = 100.0
    gate_fidelity: float = 0.99
    swap_success_prob: float = 0.5
    memory_efficiency: float = 0.9
    n_memory_qubits: int = 10

def quantum_repeater_performance(
    total_distance_km: float,
    n_segments: int,
    repeater: RepeaterNode,
    fiber_loss_db_km: float = 0.2,
    source_rate_hz: float = 10e9
) -> dict:
    """
    Estimate quantum repeater chain performance.
    
    Uses nested purification + swapping protocol analysis.
    Entanglement generation rate: R = p_success * source_rate / attempts_needed
    """
    segment_length = total_distance_km / n_segments
    eta_segment = 10 ** (-fiber_loss_db_km * segment_length / 10)
    
    # Elementary link generation
    p_link = eta_segment * repeater.memory_efficiency
    t_link = segment_length * 1e3 / (2e8)  # round-trip time in seconds (fiber c ~ 2e8 m/s)
    
    # Average attempts for link generation
    n_attempts = 1 / p_link if p_link > 0 else float('inf')
    t_gen = n_attempts * t_link  # time to generate one link
    
    # Entanglement swapping through n_segments-1 swaps (tree structure)
    n_levels = int(np.ceil(np.log2(n_segments)))
    
    # Total time including swapping
    t_total = t_gen  # bottleneck is the link generation
    for level in range(n_levels):
        t_total *= (1 / repeater.swap_success_prob)
    
    # Fidelity after swapping
    F_link = repeater.gate_fidelity ** 2 * (1 - (1 - eta_segment) * 0.01)
    F_final = F_link ** n_segments * repeater.gate_fidelity ** (2 * (n_segments - 1))
    
    # Memory requirement: must hold state for t_total
    memory_required_ms = t_total * 1000
    memory_sufficient = memory_required_ms < repeater.memory_coherence_time_ms
    
    # Decoherence during waiting
    T2 = repeater.memory_coherence_time_ms / 1000  # convert to seconds
    if T2 > 0:
        decoherence_factor = np.exp(-t_total / T2)
        F_final *= decoherence_factor
    
    rate = 1 / t_total if t_total > 0 and t_total < float('inf') else 0
    
    # Memory qubits needed per node
    qubits_per_node = max(2, int(np.ceil(2 * n_attempts)))
    
    return {
        "segment_length_km": segment_length,
        "p_link": p_link,
        "t_gen_s": t_gen,
        "t_total_s": t_total,
        "F_final": max(0, min(1, F_final)),
        "rate_hz": rate,
        "memory_required_ms": memory_required_ms,
        "memory_sufficient": memory_sufficient,
        "qubits_per_node": qubits_per_node,
        "n_levels": n_levels
    }


def run_repeater_analysis():
    """Run quantum repeater performance analysis."""
    repeater = RepeaterNode()
    
    # Vary number of segments for fixed total distance
    total_distances = [100, 200, 500, 1000]
    segment_counts = range(2, 33)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Rate vs segments
    for d in total_distances:
        rates = []
        for n in segment_counts:
            result = quantum_repeater_performance(d, n, repeater)
            rates.append(result["rate_hz"])
        axes[0, 0].semilogy(list(segment_counts), rates, 'o-', label=f'{d} km', markersize=3)
    
    axes[0, 0].set_xlabel('Number of Segments')
    axes[0, 0].set_ylabel('Entanglement Rate (Hz)')
    axes[0, 0].set_title('Repeater Chain Rate vs Segments')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Fidelity vs segments
    for d in total_distances:
        fidelities = []
        for n in segment_counts:
            result = quantum_repeater_performance(d, n, repeater)
            fidelities.append(result["F_final"])
        axes[0, 1].plot(list(segment_counts), fidelities, 'o-', label=f'{d} km', markersize=3)
    
    axes[0, 1].set_xlabel('Number of Segments')
    axes[0, 1].set_ylabel('End-to-End Fidelity')
    axes[0, 1].set_title('Fidelity vs Number of Segments')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Memory requirements
    for d in total_distances:
        mem_req = []
        for n in segment_counts:
            result = quantum_repeater_performance(d, n, repeater)
            mem_req.append(result["memory_required_ms"])
        axes[1, 0].semilogy(list(segment_counts), mem_req, 'o-', label=f'{d} km', markersize=3)
    
    axes[1, 0].axhline(y=repeater.memory_coherence_time_ms, color='r', linestyle='--', label='Memory T₂')
    axes[1, 0].set_xlabel('Number of Segments')
    axes[1, 0].set_ylabel('Required Memory Time (ms)')
    axes[1, 0].set_title('Memory Coherence Requirements')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Qubits per node
    for d in total_distances:
        qubits = []
        for n in segment_counts:
            result = quantum_repeater_performance(d, n, repeater)
            qubits.append(result["qubits_per_node"])
        axes[1, 1].plot(list(segment_counts), qubits, 'o-', label=f'{d} km', markersize=3)
    
    axes[1, 1].set_xlabel('Number of Segments')
    axes[1, 1].set_ylabel('Memory Qubits per Node')
    axes[1, 1].set_title('Memory Qubit Requirements per Node')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/repeater_performance.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Coherence time sweep
    fig, ax = plt.subplots(figsize=(8, 5))
    T2_values = np.logspace(-1, 4, 50)  # 0.1 ms to 10 s
    for d in [100, 500, 1000]:
        rates = []
        for T2 in T2_values:
            rep = RepeaterNode(memory_coherence_time_ms=T2)
            result = quantum_repeater_performance(d, 10, rep)
            rates.append(result["rate_hz"])
        ax.loglog(T2_values, rates, label=f'{d} km')
    
    ax.set_xlabel('Memory Coherence Time T₂ (ms)')
    ax.set_ylabel('Entanglement Rate (Hz)')
    ax.set_title('Impact of Memory Coherence on Repeater Performance')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/memory_coherence_impact.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Collect results
    results = {}
    for d in [100, 500, 1000]:
        opt_n = 10
        results[d] = quantum_repeater_performance(d, opt_n, repeater)
    return results


# ==============================================================================
# 3. Entanglement Distillation Protocol Efficiency
# ==============================================================================

def bennett_distillation(F_in: float, n_rounds: int = 1) -> Tuple[float, float]:
    """
    BBPSSW (Bennett et al. 1996) entanglement distillation protocol.
    
    Takes two copies of a Werner state with fidelity F_in,
    returns one copy with higher fidelity F_out.
    
    F_out = (F_in^2 + (1-F_in)^2/9) / (F_in^2 + 2*F_in*(1-F_in)/3 + 5*(1-F_in)^2/9)
    p_success = F_in^2 + 2*F_in*(1-F_in)/3 + 5*(1-F_in)^2/9
    """
    F = F_in
    total_p = 1.0
    for _ in range(n_rounds):
        numerator = F**2 + ((1 - F) / 3)**2
        denominator = F**2 + 2 * F * (1 - F) / 3 + 5 * ((1 - F) / 3)**2
        if denominator <= 0:
            return F, 0
        p_success = denominator
        F = numerator / denominator
        total_p *= p_success
    
    return F, total_p


def dejmps_distillation(F_in: float, n_rounds: int = 1) -> Tuple[float, float]:
    """
    DEJMPS (Deutsch et al. 1996) entanglement distillation.
    Bilateral CNOT-based protocol.
    
    F_out = (F^2 + (1-F)^2/9) / (F^2 + 2F(1-F)/3 + 5(1-F)^2/9)
    with phase-flip correction applied.
    """
    F = F_in
    total_p = 1.0
    for _ in range(n_rounds):
        a = F
        b = (1 - F) / 3
        
        numerator = a**2 + b**2
        denominator = (a + b)**2 + 2 * b**2
        if denominator <= 0:
            return F, 0
        p_success = denominator
        F = numerator / denominator
        total_p *= p_success
    
    return F, total_p


def run_distillation_analysis():
    """Run entanglement distillation efficiency analysis."""
    
    # Fidelity improvement vs initial fidelity
    F_init = np.linspace(0.5, 0.99, 100)
    max_rounds = 10
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Single round comparison
    F_out_bbpssw = [bennett_distillation(f, 1)[0] for f in F_init]
    F_out_dejmps = [dejmps_distillation(f, 1)[0] for f in F_init]
    p_bbpssw = [bennett_distillation(f, 1)[1] for f in F_init]
    p_dejmps = [dejmps_distillation(f, 1)[1] for f in F_init]
    
    axes[0, 0].plot(F_init, F_out_bbpssw, label='BBPSSW')
    axes[0, 0].plot(F_init, F_out_dejmps, '--', label='DEJMPS')
    axes[0, 0].plot(F_init, F_init, ':', color='gray', label='No distillation')
    axes[0, 0].set_xlabel('Input Fidelity')
    axes[0, 0].set_ylabel('Output Fidelity')
    axes[0, 0].set_title('Single-Round Distillation')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Success probability
    axes[0, 1].plot(F_init, p_bbpssw, label='BBPSSW')
    axes[0, 1].plot(F_init, p_dejmps, '--', label='DEJMPS')
    axes[0, 1].set_xlabel('Input Fidelity')
    axes[0, 1].set_ylabel('Success Probability')
    axes[0, 1].set_title('Distillation Success Probability')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Multi-round distillation
    initial_fidelities = [0.6, 0.7, 0.8, 0.9]
    rounds = range(1, max_rounds + 1)
    
    for F0 in initial_fidelities:
        F_vals = [bennett_distillation(F0, r)[0] for r in rounds]
        axes[1, 0].plot(list(rounds), F_vals, 'o-', label=f'F₀ = {F0}', markersize=4)
    
    axes[1, 0].set_xlabel('Number of Rounds')
    axes[1, 0].set_ylabel('Output Fidelity')
    axes[1, 0].set_title('BBPSSW Multi-Round Distillation')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Resource cost: pairs consumed per output pair
    for F0 in initial_fidelities:
        costs = []
        for r in rounds:
            _, p = bennett_distillation(F0, r)
            pairs_consumed = 2**r
            effective_cost = pairs_consumed / p if p > 0 else float('inf')
            costs.append(effective_cost)
        axes[1, 1].semilogy(list(rounds), costs, 'o-', label=f'F₀ = {F0}', markersize=4)
    
    axes[1, 1].set_xlabel('Number of Rounds')
    axes[1, 1].set_ylabel('Input Pairs per Output Pair')
    axes[1, 1].set_title('Distillation Resource Cost')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/entanglement_distillation.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Yield vs fidelity tradeoff
    fig, ax = plt.subplots(figsize=(8, 5))
    for r in [1, 2, 3, 5]:
        fids = []
        yields = []
        for F0 in np.linspace(0.55, 0.99, 80):
            F_out, p = bennett_distillation(F0, r)
            fids.append(F_out)
            yields.append(p / (2**r))  # yield = p_success / pairs_consumed
        ax.plot(fids, yields, label=f'{r} round(s)')
    
    ax.set_xlabel('Output Fidelity')
    ax.set_ylabel('Yield (output pairs / input pairs)')
    ax.set_title('Distillation: Yield vs Output Fidelity Tradeoff')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/distillation_yield_tradeoff.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return {
        "F0_0.7_1round": bennett_distillation(0.7, 1),
        "F0_0.7_3round": bennett_distillation(0.7, 3),
        "F0_0.8_1round": bennett_distillation(0.8, 1),
        "F0_0.8_3round": bennett_distillation(0.8, 3),
        "F0_0.9_1round": bennett_distillation(0.9, 1),
        "F0_0.9_3round": bennett_distillation(0.9, 3),
    }


# ==============================================================================
# 4. Quantum Network Routing Algorithm
# ==============================================================================

@dataclass
class QuantumLink:
    """Properties of a quantum network link."""
    distance_km: float
    fiber_loss_db_km: float = 0.2
    fidelity: float = 0.95
    success_prob: float = 0.8
    has_repeater: bool = False

def create_quantum_network(topology: str = "tokyo") -> nx.Graph:
    """Create a quantum network topology."""
    G = nx.Graph()
    
    if topology == "tokyo":
        # Simplified Tokyo QKD Network topology (based on Sasaki et al. 2011)
        nodes = {
            "Otemachi": (0, 0),
            "Koganei": (25, 10),
            "Hakusan": (5, 8),
            "Hongo": (8, 6),
            "Nezu": (9, 7),
            "Shin-Ochanomizu": (6, 4),
            "Oshiage": (12, 3),
            "Tokiwabashi": (1, 1),
        }
        
        edges = [
            ("Otemachi", "Hakusan", 12),
            ("Otemachi", "Tokiwabashi", 3),
            ("Otemachi", "Shin-Ochanomizu", 8),
            ("Hakusan", "Hongo", 5),
            ("Hakusan", "Nezu", 7),
            ("Hongo", "Nezu", 3),
            ("Hongo", "Koganei", 24),
            ("Shin-Ochanomizu", "Hongo", 6),
            ("Shin-Ochanomizu", "Oshiage", 10),
            ("Tokiwabashi", "Shin-Ochanomizu", 7),
            ("Nezu", "Oshiage", 8),
            ("Oshiage", "Koganei", 20),
        ]
    elif topology == "grid":
        nodes = {f"N{i}{j}": (i * 20, j * 20) for i in range(4) for j in range(4)}
        edges = []
        for i in range(4):
            for j in range(4):
                if i < 3:
                    edges.append((f"N{i}{j}", f"N{i+1}{j}", 20))
                if j < 3:
                    edges.append((f"N{i}{j}", f"N{i}{j+1}", 20))
    else:
        raise ValueError(f"Unknown topology: {topology}")
    
    for name, pos in nodes.items():
        G.add_node(name, pos=pos)
    
    for n1, n2, d in edges:
        eta = 10 ** (-0.2 * d / 10)
        fidelity = 0.99 * eta**0.1
        G.add_edge(n1, n2, weight=d, 
                   link=QuantumLink(distance_km=d, fidelity=fidelity, success_prob=eta))
    
    return G


def quantum_dijkstra(G: nx.Graph, source: str, target: str, metric: str = "fidelity") -> dict:
    """
    Quantum-aware shortest path using modified Dijkstra.
    
    Metrics:
    - "fidelity": maximize end-to-end fidelity (multiplicative)
    - "success": maximize end-to-end success probability
    - "distance": minimize physical distance
    """
    if metric == "fidelity":
        # Convert to additive: -log(fidelity)
        weight_func = lambda u, v, d: -np.log(d['link'].fidelity) if d['link'].fidelity > 0 else float('inf')
    elif metric == "success":
        weight_func = lambda u, v, d: -np.log(d['link'].success_prob) if d['link'].success_prob > 0 else float('inf')
    elif metric == "distance":
        weight_func = lambda u, v, d: d['weight']
    else:
        raise ValueError(f"Unknown metric: {metric}")
    
    try:
        path = nx.dijkstra_path(G, source, target, weight=weight_func)
        path_length = nx.dijkstra_path_length(G, source, target, weight=weight_func)
    except nx.NetworkXNoPath:
        return {"path": [], "metric_value": 0, "total_distance": 0}
    
    # Compute path properties
    total_distance = 0
    total_fidelity = 1.0
    total_success = 1.0
    
    for i in range(len(path) - 1):
        edge_data = G[path[i]][path[i+1]]
        total_distance += edge_data['weight']
        total_fidelity *= edge_data['link'].fidelity
        total_success *= edge_data['link'].success_prob
    
    if metric == "fidelity":
        metric_value = total_fidelity
    elif metric == "success":
        metric_value = total_success
    else:
        metric_value = total_distance
    
    return {
        "path": path,
        "metric_value": metric_value,
        "total_distance": total_distance,
        "total_fidelity": total_fidelity,
        "total_success": total_success
    }


def k_shortest_paths(G: nx.Graph, source: str, target: str, k: int = 3, metric: str = "fidelity") -> List[dict]:
    """Find k-shortest quantum paths using Yen's algorithm."""
    if metric == "fidelity":
        weight_func = lambda u, v, d: -np.log(d['link'].fidelity) if d['link'].fidelity > 0 else float('inf')
    elif metric == "success":
        weight_func = lambda u, v, d: -np.log(d['link'].success_prob) if d['link'].success_prob > 0 else float('inf')
    else:
        weight_func = lambda u, v, d: d['weight']
    
    results = []
    try:
        paths = list(nx.shortest_simple_paths(G, source, target, weight=weight_func))
        for path in paths[:k]:
            total_distance = sum(G[path[i]][path[i+1]]['weight'] for i in range(len(path) - 1))
            total_fidelity = np.prod([G[path[i]][path[i+1]]['link'].fidelity for i in range(len(path) - 1)])
            total_success = np.prod([G[path[i]][path[i+1]]['link'].success_prob for i in range(len(path) - 1)])
            
            results.append({
                "path": path,
                "total_distance": total_distance,
                "total_fidelity": total_fidelity,
                "total_success": total_success,
                "hops": len(path) - 1
            })
    except nx.NetworkXNoPath:
        pass
    
    return results


def run_routing_analysis():
    """Run quantum routing algorithm analysis."""
    G = create_quantum_network("tokyo")
    
    # Visualize network
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    pos = nx.get_node_attributes(G, 'pos')
    
    # Draw network topology
    nx.draw_networkx(G, pos, ax=axes[0], node_color='lightblue', node_size=500,
                     font_size=7, font_weight='bold', edge_color='gray', width=2)
    
    # Highlight best path
    best_path = quantum_dijkstra(G, "Otemachi", "Koganei", "fidelity")
    if best_path["path"]:
        path_edges = list(zip(best_path["path"][:-1], best_path["path"][1:]))
        nx.draw_networkx_edges(G, pos, edgelist=path_edges, ax=axes[0],
                              edge_color='red', width=3)
    
    edge_labels = {(u, v): f"{d['weight']}km" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=axes[0], font_size=6)
    axes[0].set_title('Tokyo QKD Network\n(Red: Optimal Fidelity Path)')
    
    # Compare routing metrics
    metrics = ["fidelity", "success", "distance"]
    source, target = "Otemachi", "Koganei"
    
    paths_data = {}
    for metric in metrics:
        result = quantum_dijkstra(G, source, target, metric)
        paths_data[metric] = result
    
    # Bar chart comparing metrics
    metric_names = list(paths_data.keys())
    fidelities = [paths_data[m]["total_fidelity"] for m in metric_names]
    successes = [paths_data[m]["total_success"] for m in metric_names]
    distances = [paths_data[m]["total_distance"] for m in metric_names]
    
    x = np.arange(len(metric_names))
    width = 0.25
    
    ax2 = axes[1]
    bars1 = ax2.bar(x - width, fidelities, width, label='Fidelity', color='steelblue')
    bars2 = ax2.bar(x, successes, width, label='Success Prob', color='coral')
    bars3 = ax2.bar(x + width, [d/100 for d in distances], width, label='Distance/100 km', color='seagreen')
    
    ax2.set_xlabel('Optimization Metric')
    ax2.set_ylabel('Value')
    ax2.set_title(f'Path Properties: {source} → {target}')
    ax2.set_xticks(x)
    ax2.set_xticklabels(metric_names)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/quantum_routing.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # K-shortest paths analysis
    k_paths = k_shortest_paths(G, source, target, k=5, metric="fidelity")
    
    return {"best_paths": paths_data, "k_paths": k_paths}


# ==============================================================================
# 5. Decoherence and Channel Loss Simulation
# ==============================================================================

def simulate_decoherence_channel(
    distance_km: float,
    fiber_loss_db_km: float = 0.2,
    T1_ms: float = 1000,
    T2_ms: float = 100,
    gate_time_us: float = 1.0,
    n_gates: int = 10
) -> dict:
    """
    Simulate combined effects of fiber loss and qubit decoherence.
    
    Models:
    - Fiber attenuation: exponential loss
    - T1 relaxation (amplitude damping)
    - T2 dephasing (phase damping)
    - Gate errors accumulated over operations
    """
    # Fiber transmission
    eta_fiber = 10 ** (-fiber_loss_db_km * distance_km / 10)
    
    # Time for light to traverse fiber
    t_propagation = distance_km * 1e3 / (2e8)  # seconds
    
    # Decoherence during propagation (at memory nodes)
    t_total = t_propagation + n_gates * gate_time_us * 1e-6
    
    # Amplitude damping (T1)
    gamma1 = 1 - np.exp(-t_total / (T1_ms * 1e-3)) if T1_ms > 0 else 0
    
    # Phase damping (T2)
    gamma2 = 1 - np.exp(-t_total / (T2_ms * 1e-3)) if T2_ms > 0 else 0
    
    # Combined fidelity (for Bell state)
    F_amplitude = 1 - gamma1 / 2
    F_phase = (1 + np.exp(-t_total / (T2_ms * 1e-3))) / 2 if T2_ms > 0 else 0.5
    F_combined = F_amplitude * F_phase
    
    # Gate error contribution
    gate_fidelity = 0.999  # per gate
    F_gates = gate_fidelity ** n_gates
    
    # Total fidelity
    F_total = eta_fiber * F_combined * F_gates
    
    # Secret key rate bound
    def h(x):
        if x <= 0 or x >= 1:
            return 0
        return -x * np.log2(x) - (1 - x) * np.log2(1 - x)
    
    e_eff = max(0, min(0.5, (1 - F_total) / 2))
    skr = max(0, eta_fiber * (1 - 2 * h(e_eff)))
    
    return {
        "eta_fiber": eta_fiber,
        "F_amplitude": F_amplitude,
        "F_phase": F_phase,
        "F_gates": F_gates,
        "F_total": F_total,
        "QBER": e_eff,
        "secret_key_rate": skr,
        "t_propagation_ms": t_propagation * 1000,
        "gamma1": gamma1,
        "gamma2": gamma2
    }


def run_decoherence_analysis():
    """Run decoherence and channel loss analysis."""
    distances = np.linspace(1, 300, 200)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Channel transmission vs distance
    loss_values = [0.15, 0.2, 0.25, 0.35]
    for loss in loss_values:
        eta = [10 ** (-loss * d / 10) for d in distances]
        axes[0, 0].semilogy(distances, eta, label=f'{loss} dB/km')
    
    axes[0, 0].set_xlabel('Distance (km)')
    axes[0, 0].set_ylabel('Channel Transmission η')
    axes[0, 0].set_title('Fiber Channel Transmission')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Fidelity components vs distance
    F_amp = []
    F_pha = []
    F_gate = []
    F_tot = []
    for d in distances:
        result = simulate_decoherence_channel(d)
        F_amp.append(result["F_amplitude"])
        F_pha.append(result["F_phase"])
        F_gate.append(result["F_gates"])
        F_tot.append(result["F_total"])
    
    axes[0, 1].plot(distances, F_amp, label='Amplitude Damping')
    axes[0, 1].plot(distances, F_pha, label='Phase Damping')
    axes[0, 1].plot(distances, F_gate, label='Gate Errors')
    axes[0, 1].plot(distances, F_tot, 'k--', linewidth=2, label='Total')
    axes[0, 1].set_xlabel('Distance (km)')
    axes[0, 1].set_ylabel('Fidelity')
    axes[0, 1].set_title('Fidelity Degradation Components')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Secret key rate with different T2 values
    T2_values = [10, 50, 100, 500, 1000]
    for T2 in T2_values:
        skr = [simulate_decoherence_channel(d, T2_ms=T2)["secret_key_rate"] for d in distances]
        axes[1, 0].semilogy(distances, [max(s, 1e-20) for s in skr], label=f'T₂ = {T2} ms')
    
    axes[1, 0].set_xlabel('Distance (km)')
    axes[1, 0].set_ylabel('Secret Key Rate (bits/pulse)')
    axes[1, 0].set_title('SKR with Different Memory T₂')
    axes[1, 0].set_ylim(1e-15, 1)
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # QBER vs distance
    for T2 in T2_values:
        qber = [simulate_decoherence_channel(d, T2_ms=T2)["QBER"] for d in distances]
        axes[1, 1].plot(distances, qber, label=f'T₂ = {T2} ms')
    
    axes[1, 1].axhline(y=0.11, color='r', linestyle='--', alpha=0.7, label='QBER Threshold (11%)')
    axes[1, 1].set_xlabel('Distance (km)')
    axes[1, 1].set_ylabel('QBER')
    axes[1, 1].set_title('Quantum Bit Error Rate vs Distance')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/decoherence_channel_loss.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Repeater vs direct transmission comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Direct transmission (PLOB bound: -log2(1-eta))
    direct_rates = []
    repeater_rates = []
    rep = RepeaterNode(memory_coherence_time_ms=100)
    
    for d in distances:
        eta = 10 ** (-0.2 * d / 10)
        plob = -np.log2(1 - eta) if eta < 1 else 10
        direct_rates.append(eta)  # approximate SKR proportional to eta
        
        # With repeater (10 segments)
        result = quantum_repeater_performance(d, 10, rep)
        # Normalize to per-mode rate
        repeater_rates.append(result["rate_hz"] / 1e9 if result["rate_hz"] > 0 else 1e-20)
    
    ax.semilogy(distances, direct_rates, 'b-', linewidth=2, label='Direct Transmission')
    ax.semilogy(distances, [max(r, 1e-20) for r in repeater_rates], 'r-', linewidth=2, label='10-Segment Repeater')
    
    # PLOB bound
    plob_rates = [-np.log2(1 - 10**(-0.2*d/10)) for d in distances]
    ax.semilogy(distances, plob_rates, 'k--', linewidth=1, label='PLOB Bound')
    
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Rate (per mode)')
    ax.set_title('Direct vs Repeater-Assisted Quantum Communication')
    ax.legend()
    ax.set_ylim(1e-15, 10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/direct_vs_repeater.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return {d: simulate_decoherence_channel(d) for d in [10, 50, 100, 200]}


# ==============================================================================
# 6. Tokyo QKD Network Case Study
# ==============================================================================

def run_tokyo_case_study():
    """
    Simulate Tokyo QKD Network at metropolitan scale.
    Based on Sasaki et al. (2011) topology with modern parameters.
    """
    G = create_quantum_network("tokyo")
    params = QKDParameters()
    
    # Compute key rates for all links
    link_rates = {}
    for u, v, data in G.edges(data=True):
        d = data['weight']
        bb84_result = bb84_finite_key_rate(10**10, d, params)
        link_rates[(u, v)] = {
            "distance": d,
            "rate": bb84_result["rate"],
            "QBER": bb84_result.get("e_mu", 0)
        }
    
    # Network-wide analysis
    nodes = list(G.nodes())
    n_nodes = len(nodes)
    
    # All-pairs key rates
    rate_matrix = np.zeros((n_nodes, n_nodes))
    fidelity_matrix = np.zeros((n_nodes, n_nodes))
    
    for i, src in enumerate(nodes):
        for j, dst in enumerate(nodes):
            if i == j:
                rate_matrix[i, j] = 1.0
                fidelity_matrix[i, j] = 1.0
                continue
            result = quantum_dijkstra(G, src, dst, "fidelity")
            fidelity_matrix[i, j] = result["total_fidelity"]
            
            total_dist = result["total_distance"]
            bb84_r = bb84_finite_key_rate(10**10, total_dist, params)
            rate_matrix[i, j] = bb84_r["rate"]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Key rate heatmap
    im1 = axes[0, 0].imshow(np.log10(rate_matrix + 1e-20), cmap='viridis', aspect='auto')
    axes[0, 0].set_xticks(range(n_nodes))
    axes[0, 0].set_xticklabels(nodes, rotation=45, ha='right', fontsize=7)
    axes[0, 0].set_yticks(range(n_nodes))
    axes[0, 0].set_yticklabels(nodes, fontsize=7)
    axes[0, 0].set_title('log₁₀(Key Rate) Between All Pairs')
    plt.colorbar(im1, ax=axes[0, 0])
    
    # Fidelity heatmap
    im2 = axes[0, 1].imshow(fidelity_matrix, cmap='RdYlGn', aspect='auto', vmin=0.5, vmax=1.0)
    axes[0, 1].set_xticks(range(n_nodes))
    axes[0, 1].set_xticklabels(nodes, rotation=45, ha='right', fontsize=7)
    axes[0, 1].set_yticks(range(n_nodes))
    axes[0, 1].set_yticklabels(nodes, fontsize=7)
    axes[0, 1].set_title('End-to-End Path Fidelity')
    plt.colorbar(im2, ax=axes[0, 1])
    
    # Network throughput under load
    np.random.seed(42)
    n_requests = 100
    throughputs = []
    latencies = []
    
    for _ in range(n_requests):
        src, dst = np.random.choice(nodes, 2, replace=False)
        result = quantum_dijkstra(G, src, dst, "fidelity")
        d = result["total_distance"]
        bb84_r = bb84_finite_key_rate(10**10, d, params)
        throughputs.append(bb84_r["rate"] * 10**10)  # bits per block
        latencies.append(d * 1e3 / (2e8) * 1000)  # ms
    
    axes[1, 0].hist(throughputs, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('Key Bits per Block')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Network Key Generation Distribution')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Network topology with key rates
    pos = nx.get_node_attributes(G, 'pos')
    edge_rates = []
    for u, v, data in G.edges(data=True):
        d = data['weight']
        r = bb84_finite_key_rate(10**10, d, params)["rate"]
        edge_rates.append(r)
    
    edge_rates = np.array(edge_rates)
    edge_colors = plt.cm.viridis(edge_rates / max(edge_rates) if max(edge_rates) > 0 else edge_rates)
    
    nx.draw_networkx_nodes(G, pos, ax=axes[1, 1], node_color='lightcoral', node_size=400)
    nx.draw_networkx_labels(G, pos, ax=axes[1, 1], font_size=6, font_weight='bold')
    
    edges = list(G.edges())
    for idx, (u, v) in enumerate(edges):
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], ax=axes[1, 1],
                              edge_color=[edge_colors[idx]], width=3)
    
    axes[1, 1].set_title('Tokyo QKD Network\n(Color: Key Rate)')
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/tokyo_case_study.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Scalability analysis
    fig, ax = plt.subplots(figsize=(8, 5))
    
    node_counts = [4, 6, 8, 10, 12, 16, 20]
    avg_rates = []
    avg_fidelities = []
    
    for n in node_counts:
        # Create random network with n nodes
        G_rand = nx.random_geometric_graph(n, 0.6, seed=42)
        for u, v in G_rand.edges():
            d = np.random.uniform(5, 30)
            eta = 10 ** (-0.2 * d / 10)
            G_rand[u][v]['weight'] = d
            G_rand[u][v]['link'] = QuantumLink(distance_km=d, fidelity=0.99*eta**0.1, success_prob=eta)
        
        rates = []
        fids = []
        node_list = list(G_rand.nodes())
        for i in range(min(20, len(node_list))):
            for j in range(i+1, min(20, len(node_list))):
                try:
                    result = quantum_dijkstra(G_rand, node_list[i], node_list[j], "fidelity")
                    if result["total_fidelity"] > 0:
                        fids.append(result["total_fidelity"])
                        r = bb84_finite_key_rate(10**10, result["total_distance"], params)
                        rates.append(r["rate"])
                except:
                    pass
        
        avg_rates.append(np.mean(rates) if rates else 0)
        avg_fidelities.append(np.mean(fids) if fids else 0)
    
    ax.plot(node_counts, avg_rates, 'bo-', markersize=6, label='Avg Key Rate')
    ax2 = ax.twinx()
    ax2.plot(node_counts, avg_fidelities, 'rs-', markersize=6, label='Avg Fidelity')
    
    ax.set_xlabel('Number of Network Nodes')
    ax.set_ylabel('Average Key Rate (bits/pulse)', color='b')
    ax2.set_ylabel('Average Path Fidelity', color='r')
    ax.set_title('Network Scalability Analysis')
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='lower left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/network_scalability.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return {
        "link_rates": link_rates,
        "avg_key_rate": np.mean([r["rate"] for r in link_rates.values()]),
        "avg_fidelity": np.mean(fidelity_matrix[fidelity_matrix < 1]),
        "n_nodes": n_nodes,
        "n_links": G.number_of_edges()
    }


# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    print("=" * 70)
    print("Quantum Internet Network Protocol Simulation")
    print("=" * 70)
    
    print("\n[1/6] Running BB84/E91 Finite-Key Analysis...")
    qkd_results = run_qkd_analysis()
    for d, res in qkd_results.items():
        print(f"  d={d}km: BB84 rate={res['bb84']['rate']:.2e}, E91 rate={res['e91']['rate']:.2e}")
    
    print("\n[2/6] Running Quantum Repeater Analysis...")
    repeater_results = run_repeater_analysis()
    for d, res in repeater_results.items():
        print(f"  d={d}km: rate={res['rate_hz']:.2e} Hz, F={res['F_final']:.4f}, mem_req={res['memory_required_ms']:.1f}ms")
    
    print("\n[3/6] Running Entanglement Distillation Analysis...")
    distill_results = run_distillation_analysis()
    for key, (F, p) in distill_results.items():
        print(f"  {key}: F_out={F:.4f}, p_success={p:.4f}")
    
    print("\n[4/6] Running Quantum Routing Analysis...")
    routing_results = run_routing_analysis()
    for metric, res in routing_results["best_paths"].items():
        print(f"  {metric}: path={' → '.join(res['path'])}, F={res['total_fidelity']:.4f}, d={res['total_distance']}km")
    
    print("\n[5/6] Running Decoherence & Channel Loss Analysis...")
    decoherence_results = run_decoherence_analysis()
    for d, res in decoherence_results.items():
        print(f"  d={d}km: η={res['eta_fiber']:.2e}, F={res['F_total']:.4f}, SKR={res['secret_key_rate']:.2e}")
    
    print("\n[6/6] Running Tokyo QKD Network Case Study...")
    tokyo_results = run_tokyo_case_study()
    print(f"  Nodes: {tokyo_results['n_nodes']}, Links: {tokyo_results['n_links']}")
    print(f"  Avg key rate: {tokyo_results['avg_key_rate']:.2e} bits/pulse")
    print(f"  Avg path fidelity: {tokyo_results['avg_fidelity']:.4f}")
    
    # Save numerical results
    print("\n" + "=" * 70)
    print("All simulations complete. Figures saved to figures/")
    print("=" * 70)
    
    # Save summary data for report
    summary = {
        "qkd": {str(k): {
            "bb84_rate": v["bb84"]["rate"],
            "e91_rate": v["e91"]["rate"],
            "bb84_qber": v["bb84"].get("e_mu", 0),
            "e91_qber": v["e91"].get("QBER", 0)
        } for k, v in qkd_results.items()},
        "repeater": {str(k): {
            "rate_hz": v["rate_hz"],
            "fidelity": v["F_final"],
            "memory_ms": v["memory_required_ms"],
            "qubits_per_node": v["qubits_per_node"]
        } for k, v in repeater_results.items()},
        "distillation": {k: {"F_out": v[0], "p_success": v[1]} for k, v in distill_results.items()},
        "tokyo": {
            "avg_key_rate": tokyo_results["avg_key_rate"],
            "avg_fidelity": tokyo_results["avg_fidelity"],
            "n_nodes": tokyo_results["n_nodes"],
            "n_links": tokyo_results["n_links"]
        }
    }
    
    with open("simulation_results.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    
    return summary


if __name__ == "__main__":
    results = main()
