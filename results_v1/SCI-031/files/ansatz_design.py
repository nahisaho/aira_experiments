"""
Ansatz Design: Hardware-Efficient vs Chemically-Inspired (UCCSD)
Compares expressibility, entanglement capability, and convergence behavior.
"""

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
import json
import os
import time
from datetime import datetime

np.random.seed(42)

# ── Hardware-Efficient Ansatz ──────────────────────────────────────────────

def hardware_efficient_ansatz(params, n_qubits, n_layers, wires):
    """Ry-CNOT hardware-efficient ansatz."""
    for layer in range(n_layers):
        for i in range(n_qubits):
            qml.RY(params[layer, i, 0], wires=wires[i])
            qml.RZ(params[layer, i, 1], wires=wires[i])
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[wires[i], wires[i + 1]])

def he_param_shape(n_qubits, n_layers):
    return (n_layers, n_qubits, 2)

# ── UCCSD-inspired Ansatz ──────────────────────────────────────────────────

def uccsd_ansatz(params, n_qubits, singles, doubles, wires):
    """
    UCCSD-inspired ansatz using Givens rotations for singles
    and double excitations for doubles.
    """
    qml.BasisState(np.array([1, 1] + [0] * (n_qubits - 2)), wires=wires)
    idx = 0
    for (i, a) in singles:
        qml.SingleExcitation(params[idx], wires=[wires[i], wires[a]])
        idx += 1
    for (i, j, a, b) in doubles:
        qml.DoubleExcitation(params[idx], wires=[wires[i], wires[j], wires[a], wires[b]])
        idx += 1

# ── Expressibility Metric (Haar-distance approximation) ───────────────────

def compute_expressibility(circuit_fn, param_shape, n_samples=200, n_qubits=4):
    """
    Approximates expressibility as the average fidelity variance
    over random parameter pairs (lower Haar distance = higher expressibility).
    """
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def state_circuit(params):
        circuit_fn(params)
        return qml.state()

    fidelities = []
    for _ in range(n_samples):
        p1 = np.random.uniform(0, 2 * np.pi, param_shape)
        p2 = np.random.uniform(0, 2 * np.pi, param_shape)
        s1 = state_circuit(p1)
        s2 = state_circuit(p2)
        fid = abs(np.dot(s1.conj(), s2)) ** 2
        fidelities.append(float(fid))
    return float(np.mean(fidelities)), float(np.std(fidelities))

# ── Entanglement Capability ────────────────────────────────────────────────

def compute_entanglement_capability(circuit_fn, param_shape, n_samples=100, n_qubits=4):
    """
    Meyer-Wallach entanglement measure averaged over random parameters.
    """
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def state_circuit(params):
        circuit_fn(params)
        return qml.state()

    entanglements = []
    for _ in range(n_samples):
        params = np.random.uniform(0, 2 * np.pi, param_shape)
        state = state_circuit(params)
        # Simplified 2-qubit entanglement via reduced density matrix purity
        state_mat = state.reshape([2] * n_qubits)
        purities = []
        for q in range(n_qubits):
            # Trace out all other qubits
            axes = list(range(n_qubits))
            axes.remove(q)
            rho = np.tensordot(state_mat, state_mat.conj(), axes=(axes, axes))
            purity = float(np.real(np.trace(rho @ rho)))
            purities.append(purity)
        # Entanglement = 1 - avg purity (normalized)
        ent = 1.0 - float(np.mean(purities))
        entanglements.append(ent)
    return float(np.mean(entanglements)), float(np.std(entanglements))

# ── VQE Convergence Comparison ─────────────────────────────────────────────

def h2_hamiltonian():
    """Minimal H2 Hamiltonian (STO-3G, 2 qubits, JW mapping)."""
    coeffs = [
        -0.4804,  # identity
        +0.3435,  # Z0
        -0.4347,  # Z1
        +0.5716,  # Z0 Z1
        +0.0910,  # X0 X1
        +0.0910,  # Y0 Y1
    ]
    obs = [
        qml.Identity(0),
        qml.PauliZ(0),
        qml.PauliZ(1),
        qml.PauliZ(0) @ qml.PauliZ(1),
        qml.PauliX(0) @ qml.PauliX(1),
        qml.PauliY(0) @ qml.PauliY(1),
    ]
    return qml.Hamiltonian(coeffs, obs)

def run_vqe_convergence(ansatz_name, n_steps=80):
    """Run VQE and return energy history."""
    n_qubits = 2
    H = h2_hamiltonian()
    dev = qml.device("default.qubit", wires=n_qubits)

    if ansatz_name == "hardware_efficient":
        n_layers = 2
        shape = he_param_shape(n_qubits, n_layers)
        params = pnp.array(np.random.uniform(-np.pi, np.pi, shape), requires_grad=True)

        @qml.qnode(dev)
        def circuit(params):
            hardware_efficient_ansatz(params, n_qubits, n_layers, list(range(n_qubits)))
            return qml.expval(H)

    elif ansatz_name == "uccsd":
        singles = [(0, 1)]
        doubles = []
        n_params = len(singles) + len(doubles)
        params = pnp.array(np.zeros(n_params), requires_grad=True)

        @qml.qnode(dev)
        def circuit(params):
            # H2: HF reference |10⟩ (one electron in bonding orbital)
            qml.BasisState(np.array([1, 0]), wires=[0, 1])
            uccsd_ansatz(params, n_qubits, singles, doubles, list(range(n_qubits)))
            return qml.expval(H)

    opt = qml.AdamOptimizer(stepsize=0.05)
    energies = []
    t0 = time.time()
    for step in range(n_steps):
        params, cost = opt.step_and_cost(circuit, params)
        energies.append(float(cost))
    elapsed = time.time() - t0
    return energies, elapsed

# ── Main ───────────────────────────────────────────────────────────────────

def main():
    results = {}
    n_qubits = 4

    # Expressibility analysis
    print("Computing expressibility metrics...")
    he_shape = he_param_shape(n_qubits, 2)

    def he_fn(params):
        hardware_efficient_ansatz(params, n_qubits, 2, list(range(n_qubits)))

    he_expr_mean, he_expr_std = compute_expressibility(he_fn, he_shape, n_samples=150)
    he_ent_mean, he_ent_std = compute_entanglement_capability(he_fn, he_shape, n_samples=80)

    results["hardware_efficient"] = {
        "expressibility_fidelity_mean": he_expr_mean,
        "expressibility_fidelity_std": he_expr_std,
        "entanglement_mean": he_ent_mean,
        "entanglement_std": he_ent_std,
    }
    print(f"  HE ansatz: expr={he_expr_mean:.4f}±{he_expr_std:.4f}, ent={he_ent_mean:.4f}±{he_ent_std:.4f}")

    # UCCSD-like expressibility with 2-qubit circuit
    n_qubits2 = 2

    def uccsd_fn_expr(params):
        qml.BasisState(np.array([1, 0]), wires=[0, 1])
        qml.SingleExcitation(params[0], wires=[0, 1])

    uccsd_expr_mean, uccsd_expr_std = compute_expressibility(uccsd_fn_expr, (1,), n_samples=150, n_qubits=2)

    results["uccsd"] = {
        "expressibility_fidelity_mean": uccsd_expr_mean,
        "expressibility_fidelity_std": uccsd_expr_std,
    }
    print(f"  UCCSD ansatz: expr={uccsd_expr_mean:.4f}±{uccsd_expr_std:.4f}")

    # VQE convergence
    print("Running VQE convergence comparison...")
    he_energies, he_time = run_vqe_convergence("hardware_efficient", n_steps=80)
    uccsd_energies, uccsd_time = run_vqe_convergence("uccsd", n_steps=80)

    results["vqe_convergence"] = {
        "hardware_efficient": {
            "energies": he_energies,
            "final_energy": he_energies[-1],
            "time_sec": he_time,
            "min_energy": min(he_energies),
        },
        "uccsd": {
            "energies": uccsd_energies,
            "final_energy": uccsd_energies[-1],
            "time_sec": uccsd_time,
            "min_energy": min(uccsd_energies),
        },
        "exact_energy_h2": -1.8511,   # electronic eigenvalue (FCI); +0.7154 nuclear = -1.1357 Ha total
        "fci_total_energy": -1.1357,
    }

    print(f"  HE final energy: {he_energies[-1]:.6f} Ha (min: {min(he_energies):.6f})")
    print(f"  UCCSD final energy: {uccsd_energies[-1]:.6f} Ha (min: {min(uccsd_energies):.6f})")

    # Gate counts
    results["gate_counts"] = {
        "hardware_efficient_2layer_4qubit": {
            "single_qubit": n_qubits * 2 * 2,  # 2 layers * n_qubits * 2 rotations
            "two_qubit": (n_qubits - 1) * 2,
            "n_params": he_shape[0] * he_shape[1] * he_shape[2],
        },
        "uccsd_h2_2qubit": {
            "single_qubit": 4,
            "two_qubit": 1,
            "n_params": 1,
        }
    }

    os.makedirs("results", exist_ok=True)
    with open("results/ansatz_comparison.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Saved: results/ansatz_comparison.json")
    return results

if __name__ == "__main__":
    main()
