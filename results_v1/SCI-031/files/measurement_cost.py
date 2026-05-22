"""
Measurement Cost Reduction: Qubit Grouping and Classical Shadow Tomography
"""

import numpy as np
import pennylane as qml
import json
import os
from itertools import combinations

np.random.seed(42)

# ── Pauli Grouping (Qubit-Wise Commutativity) ──────────────────────────────

def pauli_strings_commute_qwc(p1_ops, p2_ops):
    """
    Check qubit-wise commutativity between two Pauli strings
    (represented as dicts {qubit: 'X'/'Y'/'Z'}).
    """
    for qubit in set(p1_ops) & set(p2_ops):
        if p1_ops[qubit] != p2_ops[qubit]:
            return False
    return True

def group_paulis_qwc(pauli_list):
    """
    Greedy qubit-wise commutativity grouping.
    Returns list of groups (each group is a list of Pauli string dicts).
    """
    groups = []
    for pauli in pauli_list:
        placed = False
        for group in groups:
            if all(pauli_strings_commute_qwc(pauli, member) for member in group):
                group.append(pauli)
                placed = True
                break
        if not placed:
            groups.append([pauli])
    return groups

def extract_pauli_ops(hamiltonian_terms):
    """
    Extract Pauli operators from Hamiltonian term list.
    Each term: (coeff, {qubit: 'X'/'Y'/'Z'}).
    """
    paulis = []
    for coeff, ops in hamiltonian_terms:
        if ops:  # skip identity
            paulis.append(ops)
    return paulis

def lih_hamiltonian_terms():
    """Simplified LiH Hamiltonian terms (selected Pauli strings for grouping demo)."""
    # Pauli strings represented as {qubit: operator}
    terms = [
        (-0.2, {0: 'Z'}),
        (+0.3, {1: 'Z'}),
        (-0.1, {2: 'Z'}),
        (+0.2, {3: 'Z'}),
        (+0.05, {0: 'X', 1: 'X'}),
        (+0.05, {0: 'Y', 1: 'Y'}),
        (+0.03, {1: 'X', 2: 'X'}),
        (+0.03, {1: 'Y', 2: 'Y'}),
        (+0.02, {2: 'X', 3: 'X'}),
        (+0.02, {2: 'Y', 3: 'Y'}),
        (+0.04, {0: 'Z', 1: 'Z'}),
        (+0.04, {1: 'Z', 2: 'Z'}),
        (+0.03, {0: 'Z', 2: 'Z'}),
        (+0.02, {0: 'X', 1: 'Z', 2: 'X'}),
        (+0.02, {0: 'Y', 1: 'Z', 2: 'Y'}),
        (+0.01, {0: 'X', 1: 'X', 2: 'Z', 3: 'Z'}),
        (+0.01, {0: 'Y', 1: 'Y', 2: 'Z', 3: 'Z'}),
    ]
    return terms

def analyze_grouping(hamiltonian_terms):
    """Analyze grouping efficiency."""
    paulis = [ops for (coeff, ops) in hamiltonian_terms if ops]
    n_total = len(paulis)
    groups = group_paulis_qwc(paulis)
    n_groups = len(groups)
    group_sizes = [len(g) for g in groups]
    reduction = 1.0 - n_groups / n_total
    return {
        "n_total_terms": n_total,
        "n_groups_qwc": n_groups,
        "group_sizes": group_sizes,
        "measurement_reduction_ratio": reduction,
        "speedup_factor": n_total / n_groups if n_groups > 0 else 1.0,
    }

# ── Classical Shadow Tomography ────────────────────────────────────────────

class ClassicalShadow:
    """
    Classical Shadow Tomography for efficient expectation value estimation.
    Uses random Pauli basis measurements.
    """
    def __init__(self, n_qubits, n_shadows=1000):
        self.n_qubits = n_qubits
        self.n_shadows = n_shadows
        self.dev = qml.device("default.qubit", wires=n_qubits)

    def take_shadows(self, state_prep_fn, params=None):
        """
        Take n_shadows snapshots with random Pauli basis measurements.
        Returns list of (basis_choice, measurement_outcomes).
        Uses probs-based simulation instead of sampling.
        """
        # Use shot-based device for sampling
        dev_shots = qml.device("default.qubit", wires=self.n_qubits, shots=1)
        shadows = []
        for _ in range(self.n_shadows):
            basis = np.random.choice(['X', 'Y', 'Z'], size=self.n_qubits)

            @qml.qnode(dev_shots)
            def circuit(b=basis):
                if params is not None:
                    state_prep_fn(params)
                else:
                    state_prep_fn()
                for q in range(self.n_qubits):
                    if b[q] == 'X':
                        qml.Hadamard(wires=q)
                    elif b[q] == 'Y':
                        qml.adjoint(qml.S)(wires=q)
                        qml.Hadamard(wires=q)
                return [qml.sample(qml.PauliZ(q)) for q in range(self.n_qubits)]

            outcomes = circuit()
            outcomes = np.array(outcomes).flatten()
            shadows.append((basis.copy(), outcomes))
        return shadows

    def estimate_pauli_expectation(self, shadows, pauli_ops):
        """
        Estimate ⟨P⟩ from classical shadows.
        pauli_ops: dict {qubit: 'X'/'Y'/'Z'}
        """
        estimates = []
        for basis, outcomes in shadows:
            # Check if the shadow's basis matches the Pauli string
            if all(basis[q] == pauli_ops.get(q, basis[q]) for q in range(self.n_qubits)):
                val = 1.0
                for q, op in pauli_ops.items():
                    val *= float(-outcomes[q])  # Z eigenvalue: +1 for 0, -1 for 1
                estimates.append(val)
        if len(estimates) == 0:
            return 0.0, float('inf')
        return float(np.mean(estimates)), float(np.std(estimates) / np.sqrt(len(estimates)) if len(estimates) > 1 else 0.0)

    def estimate_hamiltonian_energy(self, shadows, hamiltonian_terms):
        """Estimate full Hamiltonian expectation value from shadows."""
        total_energy = 0.0
        total_var = 0.0
        identity_coeff = sum(c for c, ops in hamiltonian_terms if not ops)
        total_energy += identity_coeff

        for coeff, ops in hamiltonian_terms:
            if not ops:
                continue
            mean_val, se = self.estimate_pauli_expectation(shadows, ops)
            total_energy += coeff * mean_val
            total_var += (coeff * se) ** 2
        return total_energy, float(np.sqrt(total_var))

# ── Demo: Shadow vs Direct Measurement Comparison ─────────────────────────

def demo_shadow_vs_direct(n_qubits=4, n_shadows_list=None):
    """Compare Classical Shadow vs direct measurement for energy estimation."""
    if n_shadows_list is None:
        n_shadows_list = [100, 300, 500, 1000, 2000]

    hamiltonian_terms = lih_hamiltonian_terms()
    dev_exact = qml.device("default.qubit", wires=n_qubits)

    # Build PennyLane Hamiltonian for exact reference
    coeffs = []
    obs_list = []
    pauli_map = {'X': qml.PauliX, 'Y': qml.PauliY, 'Z': qml.PauliZ}

    for coeff, ops in hamiltonian_terms:
        coeffs.append(coeff)
        if not ops:
            obs_list.append(qml.Identity(0))
        else:
            obs = None
            for q, p in ops.items():
                term = pauli_map[p](q)
                obs = term if obs is None else obs @ term
            obs_list.append(obs)

    H = qml.Hamiltonian(coeffs, obs_list)

    # Reference state (|0101⟩-like for 4 qubits)
    ref_params = np.array([0.3, 0.5, -0.2, 0.4])

    def state_prep(params):
        for q in range(n_qubits):
            qml.RY(params[q], wires=q)
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q + 1])

    @qml.qnode(dev_exact)
    def exact_energy(params):
        state_prep(params)
        return qml.expval(H)

    exact_val = float(exact_energy(ref_params))
    print(f"  Exact energy (direct): {exact_val:.6f}")

    shadow_results = []
    for n_s in n_shadows_list:
        cs = ClassicalShadow(n_qubits, n_shadows=n_s)
        shadows = cs.take_shadows(state_prep, ref_params)
        est_energy, est_error = cs.estimate_hamiltonian_energy(shadows, hamiltonian_terms)
        shadow_results.append({
            "n_shadows": n_s,
            "estimated_energy": est_energy,
            "std_error": est_error,
            "abs_error": abs(est_energy - exact_val),
        })
        print(f"  n_shadows={n_s:5d}: E={est_energy:.4f} ± {est_error:.4f}, "
              f"|err|={abs(est_energy - exact_val):.4f}")

    return exact_val, shadow_results

# ── Measurement Shot Analysis ──────────────────────────────────────────────

def measurement_shot_analysis(n_terms_list=None):
    """
    Compare required shots: naive vs grouped vs shadow.
    """
    if n_terms_list is None:
        n_terms_list = [10, 20, 50, 100, 200, 500]
    epsilon = 0.01  # target precision
    delta = 0.05    # failure probability

    results = []
    for n_terms in n_terms_list:
        # Naive: O(n_terms * 1/eps^2)
        naive_shots = int(n_terms * (1 / epsilon ** 2))
        # Grouped (assume ~sqrt(n_terms) groups)
        n_groups = max(1, int(np.sqrt(n_terms)))
        grouped_shots = int(n_groups * (1 / epsilon ** 2))
        # Classical shadow: O(log(n_terms) / eps^2)
        shadow_shots = int(np.log(n_terms + 1) * (1 / epsilon ** 2) * 3)  # 3x overhead
        results.append({
            "n_terms": n_terms,
            "naive_shots": naive_shots,
            "grouped_shots": grouped_shots,
            "shadow_shots": shadow_shots,
            "grouping_speedup": naive_shots / grouped_shots,
            "shadow_speedup": naive_shots / shadow_shots,
        })
    return results

# ── Main ───────────────────────────────────────────────────────────────────

def main():
    results = {}

    print("Analyzing qubit-wise commutativity grouping...")
    terms = lih_hamiltonian_terms()
    grouping_stats = analyze_grouping(terms)
    results["grouping"] = grouping_stats
    print(f"  {grouping_stats['n_total_terms']} terms → {grouping_stats['n_groups_qwc']} groups "
          f"({grouping_stats['measurement_reduction_ratio']*100:.1f}% reduction, "
          f"{grouping_stats['speedup_factor']:.1f}x speedup)")

    print("Running Classical Shadow comparison...")
    exact_energy, shadow_results = demo_shadow_vs_direct(n_shadows_list=[100, 300, 500, 1000])
    results["classical_shadow"] = {
        "exact_energy": exact_energy,
        "shadow_estimates": shadow_results,
    }

    print("Analyzing measurement shot scaling...")
    shot_analysis = measurement_shot_analysis()
    results["shot_scaling"] = shot_analysis

    os.makedirs("results", exist_ok=True)
    with open("results/measurement_cost.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Saved: results/measurement_cost.json")
    return results

if __name__ == "__main__":
    main()
