"""
Barren Plateau Analysis and Avoidance Strategies
"""

import numpy as np
import pennylane as qml
import json
import os

np.random.seed(42)

# ── Gradient Variance Analysis ─────────────────────────────────────────────

def compute_gradient_variance(circuit_fn, param_shape, n_samples=200, n_qubits=4):
    """
    Compute variance of gradients of the first parameter across random initializations.
    Barren plateau: variance exponentially vanishes with system size.
    Uses parameter-shift rule via finite differences.
    """
    dev = qml.device("default.qubit", wires=n_qubits)

    def cost_fn(p):
        @qml.qnode(dev)
        def circuit():
            circuit_fn(p)
            return qml.expval(qml.PauliZ(0))
        return float(circuit())

    gradients = []
    eps = 1e-4
    for _ in range(n_samples):
        params = np.random.uniform(0, 2 * np.pi, param_shape)
        # Finite-difference gradient for the first parameter element
        p_plus = params.copy()
        p_minus = params.copy()
        flat = p_plus.flatten()
        flat[0] += eps
        p_plus = flat.reshape(param_shape)
        flat2 = p_minus.flatten()
        flat2[0] -= eps
        p_minus = flat2.reshape(param_shape)
        grad = (cost_fn(p_plus) - cost_fn(p_minus)) / (2 * eps)
        gradients.append(float(grad))

    return float(np.var(gradients)), float(np.mean(gradients))

def random_circuit(params, n_qubits, n_layers):
    """Random deep circuit prone to barren plateaus."""
    for layer in range(n_layers):
        for q in range(n_qubits):
            qml.RY(params[layer, q, 0], wires=q)
            qml.RZ(params[layer, q, 1], wires=q)
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q + 1])
        qml.CNOT(wires=[n_qubits - 1, 0])  # periodic

def shallow_local_circuit(params, n_qubits, n_layers):
    """Shallow circuit with local operations - avoids barren plateaus."""
    for layer in range(min(n_layers, 2)):  # keep shallow
        for q in range(n_qubits):
            qml.RY(params[layer, q, 0], wires=q)
        for q in range(0, n_qubits - 1, 2):  # local 2-qubit gates only
            qml.CNOT(wires=[q, q + 1])

def analyze_barren_plateau_vs_depth(n_qubits=4):
    """
    Analyze gradient variance as a function of circuit depth.
    """
    depths = [1, 2, 3, 4, 6, 8, 10]
    results_deep = []
    results_shallow = []

    for depth in depths:
        shape = (depth, n_qubits, 2)
        # Deep random circuit
        def deep_fn(params, d=depth):
            random_circuit(params, n_qubits, d)
        var_deep, mean_deep = compute_gradient_variance(deep_fn, shape, n_samples=150)

        # Shallow local circuit
        def shallow_fn(params, d=depth):
            shallow_local_circuit(params, n_qubits, d)
        var_shallow, mean_shallow = compute_gradient_variance(shallow_fn, shape, n_samples=150)

        results_deep.append({"depth": depth, "grad_variance": var_deep, "grad_mean": mean_deep})
        results_shallow.append({"depth": depth, "grad_variance": var_shallow, "grad_mean": mean_shallow})
        print(f"  depth={depth:2d}: deep_var={var_deep:.2e}, shallow_var={var_shallow:.2e}")

    return results_deep, results_shallow

def analyze_barren_plateau_vs_qubits(depth=4):
    """
    Analyze gradient variance as a function of system size.
    Barren plateau: exponential decay in n_qubits.
    """
    qubit_sizes = [2, 3, 4, 5, 6]
    results = []

    for n_qubits in qubit_sizes:
        shape = (depth, n_qubits, 2)
        dev = qml.device("default.qubit", wires=n_qubits)
        gradients = []
        eps = 1e-4
        for _ in range(150):
            params = np.random.uniform(0, 2 * np.pi, shape)

            def cost_fn(p, nq=n_qubits, d=depth):
                @qml.qnode(dev)
                def circuit():
                    random_circuit(p, nq, d)
                    return qml.expval(qml.PauliZ(0))
                return float(circuit())

            p_plus = params.copy(); flat = p_plus.flatten(); flat[0] += eps
            p_plus = flat.reshape(shape)
            p_minus = params.copy(); flat2 = p_minus.flatten(); flat2[0] -= eps
            p_minus = flat2.reshape(shape)
            grad = (cost_fn(p_plus) - cost_fn(p_minus)) / (2 * eps)
            gradients.append(float(grad))

        var = float(np.var(gradients))
        results.append({
            "n_qubits": n_qubits,
            "depth": depth,
            "grad_variance": var,
            "log_variance": float(np.log(var + 1e-20)),
        })
        print(f"  n_qubits={n_qubits}: var={var:.2e}")

    return results

# ── Avoidance Strategies ──────────────────────────────────────────────────

def identity_blocks_init(params_shape, n_qubits):
    """
    Initialize parameters such that blocks form identity gates
    (avoids barren plateaus by starting near identity).
    """
    params = np.zeros(params_shape)
    # Small perturbation from identity
    params += np.random.normal(0, 0.01, params_shape)
    return params

def layerwise_training_demo(n_qubits=4, n_layers_total=6, n_steps_per_layer=30):
    """
    Layerwise (greedy) training: train one layer at a time,
    keeping previous layers fixed.
    """
    dev = qml.device("default.qubit", wires=n_qubits)
    shape_per_layer = (n_qubits, 2)
    all_params = [np.zeros(shape_per_layer) for _ in range(n_layers_total)]

    energy_history = []
    for layer_idx in range(n_layers_total):
        # Initialize new layer with small perturbation
        all_params[layer_idx] = np.random.normal(0, 0.05, shape_per_layer)
        trainable_params = all_params[layer_idx]

        frozen_params = [p.copy() for i, p in enumerate(all_params) if i != layer_idx]

        def make_cost(layer_i, frozen):
            @qml.qnode(dev)
            def cost(p):
                # Run all layers (frozen + trainable)
                for l in range(layer_i + 1):
                    if l == layer_i:
                        current = p
                    else:
                        offset = l if l < layer_i else l - 1
                        current = frozen[l]
                    for q in range(n_qubits):
                        qml.RY(current[q, 0], wires=q)
                        qml.RZ(current[q, 1], wires=q)
                    for q in range(n_qubits - 1):
                        qml.CNOT(wires=[q, q + 1])
                return qml.expval(qml.PauliZ(0))
            return cost

        cost_fn = make_cost(layer_idx, frozen_params)
        opt = qml.AdamOptimizer(stepsize=0.05)
        for step in range(n_steps_per_layer):
            trainable_params, cost = opt.step_and_cost(cost_fn, trainable_params)
            energy_history.append(float(cost))
        all_params[layer_idx] = trainable_params
        print(f"  Layer {layer_idx+1}/{n_layers_total} done, energy={cost:.4f}")

    return energy_history

def correlated_init_demo(n_qubits=4, n_layers=4, n_steps=60):
    """
    Correlated parameter initialization with identity-preserving blocks.
    Compare with random initialization.
    """
    dev = qml.device("default.qubit", wires=n_qubits)
    shape = (n_layers, n_qubits, 2)

    def circuit(params):
        random_circuit(params, n_qubits, n_layers)
        return qml.expval(qml.PauliZ(0))

    qnode = qml.QNode(circuit, dev)

    # Random init
    energies_random = []
    params_rand = np.random.uniform(0, 2 * np.pi, shape)
    opt = qml.AdamOptimizer(stepsize=0.05)
    for _ in range(n_steps):
        params_rand, cost = opt.step_and_cost(qnode, params_rand)
        energies_random.append(float(cost))

    # Identity-block init
    energies_ident = []
    params_ident = identity_blocks_init(shape, n_qubits)
    opt2 = qml.AdamOptimizer(stepsize=0.05)
    for _ in range(n_steps):
        params_ident, cost = opt2.step_and_cost(qnode, params_ident)
        energies_ident.append(float(cost))

    return energies_random, energies_ident

# ── Main ───────────────────────────────────────────────────────────────────

def main():
    results = {}

    print("Analyzing barren plateaus vs circuit depth...")
    depth_results_deep, depth_results_shallow = analyze_barren_plateau_vs_depth(n_qubits=4)
    results["depth_analysis"] = {
        "deep_circuit": depth_results_deep,
        "shallow_local": depth_results_shallow,
    }

    print("Analyzing barren plateaus vs system size...")
    qubit_results = analyze_barren_plateau_vs_qubits(depth=4)
    results["qubit_size_analysis"] = qubit_results

    print("Running identity-block init comparison...")
    e_rand, e_ident = correlated_init_demo(n_qubits=4, n_layers=4, n_steps=60)
    results["init_comparison"] = {
        "random_init": e_rand,
        "identity_block_init": e_ident,
        "random_final": e_rand[-1],
        "identity_final": e_ident[-1],
        "improvement": e_rand[-1] - e_ident[-1],
    }
    print(f"  Random init final: {e_rand[-1]:.4f}, Identity init final: {e_ident[-1]:.4f}")

    print("Running layerwise training...")
    lw_history = layerwise_training_demo(n_qubits=4, n_layers_total=4, n_steps_per_layer=25)
    results["layerwise_training"] = {
        "energy_history": lw_history,
        "final_energy": lw_history[-1] if lw_history else None,
    }

    # Theoretical scaling analysis
    n_q_vals = [2, 4, 6, 8, 10, 12]
    theory_variance = [2 ** (-n) for n in n_q_vals]
    results["theoretical_scaling"] = {
        "n_qubits": n_q_vals,
        "variance_scaling": theory_variance,
        "description": "Variance ~ 2^(-n) for global Pauli observable in random circuits",
    }

    os.makedirs("results", exist_ok=True)
    with open("results/barren_plateau.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Saved: results/barren_plateau.json")
    return results

if __name__ == "__main__":
    main()
