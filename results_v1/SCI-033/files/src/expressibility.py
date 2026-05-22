import json
from pathlib import Path

import numpy as np
import pennylane as qml
from scipy.stats import entropy

np.random.seed(42)


def meyer_wallach_entanglement(statevector, n_qubits):
    """Compute the Meyer-Wallach entanglement measure Q(psi)."""
    psi = np.asarray(statevector, dtype=np.complex128)
    norm = np.linalg.norm(psi)
    if norm == 0:
        return 0.0
    psi = psi / norm
    psi_tensor = psi.reshape([2] * n_qubits)
    total = 0.0
    for qubit in range(n_qubits):
        reduced = np.moveaxis(psi_tensor, qubit, 0).reshape(2, -1)
        rho = reduced @ reduced.conj().T
        purity = float(np.real(np.trace(rho @ rho)))
        total += 1.0 - purity
    return float((2.0 / n_qubits) * total)


def _ring_entanglers(n_qubits):
    if n_qubits < 2:
        return
    for wire in range(n_qubits - 1):
        qml.CNOT(wires=[wire, wire + 1])
    if n_qubits > 2:
        qml.CNOT(wires=[n_qubits - 1, 0])


def _hea(params, wires):
    depth = params.shape[0]
    for layer in range(depth):
        for wire in wires:
            qml.RY(params[layer, wire, 0], wires=wire)
            qml.RZ(params[layer, wire, 1], wires=wire)
            qml.RX(params[layer, wire, 2], wires=wire)
        _ring_entanglers(len(wires))


def _strongly_entangling(params, wires):
    qml.StronglyEntanglingLayers(params, wires=wires)


def _random_ansatz(params, wires):
    depth = params.shape[0]
    for layer in range(depth):
        for wire in wires:
            if (layer + wire) % 3 == 0:
                qml.RX(params[layer, wire, 0], wires=wire)
                qml.RY(params[layer, wire, 1], wires=wire)
                qml.RZ(params[layer, wire, 2], wires=wire)
            elif (layer + wire) % 3 == 1:
                qml.RZ(params[layer, wire, 0], wires=wire)
                qml.RX(params[layer, wire, 1], wires=wire)
                qml.RY(params[layer, wire, 2], wires=wire)
            else:
                qml.RY(params[layer, wire, 0], wires=wire)
                qml.RX(params[layer, wire, 1], wires=wire)
                qml.RZ(params[layer, wire, 2], wires=wire)
        for wire in range(0, len(wires) - 1, 2):
            qml.CZ(wires=[wire, wire + 1])
        for wire in range(1, len(wires) - 1, 2):
            qml.CZ(wires=[wire, wire + 1])


def _qaoa_style(params, wires):
    gammas = params[:, 0]
    betas = params[:, 1]
    for wire in wires:
        qml.Hadamard(wires=wire)
    for gamma, beta in zip(gammas, betas):
        for wire in range(len(wires) - 1):
            qml.IsingZZ(gamma, wires=[wire, wire + 1])
        if len(wires) > 2:
            qml.IsingZZ(gamma, wires=[len(wires) - 1, 0])
        for wire in wires:
            qml.RX(2.0 * beta, wires=wire)


def _build_state_circuit(ansatz_name, n_qubits, depth):
    dev = qml.device("default.qubit", wires=n_qubits)
    wires = list(range(n_qubits))
    if ansatz_name.startswith("HEA"):
        param_shape = (depth, n_qubits, 3)
        ansatz = _hea
    elif ansatz_name == "StronglyEntangling":
        param_shape = qml.StronglyEntanglingLayers.shape(n_layers=depth, n_wires=n_qubits)
        ansatz = _strongly_entangling
    elif ansatz_name == "Random":
        param_shape = (depth, n_qubits, 3)
        ansatz = _random_ansatz
    elif ansatz_name == "QAOA-style":
        param_shape = (depth, 2)
        ansatz = _qaoa_style
    else:
        raise ValueError(f"Unknown ansatz: {ansatz_name}")

    @qml.qnode(dev)
    def circuit(params):
        ansatz(params, wires)
        return qml.state()

    return circuit, param_shape


def compute_expressibility(circuit_fn, n_qubits, param_shape, n_samples=2000, n_bins=60, seed=42):
    """Compute expressibility via KL divergence from the Haar fidelity distribution."""
    rng = np.random.default_rng(seed)
    fidelities = np.empty(n_samples, dtype=float)
    entanglements = np.empty(n_samples, dtype=float)

    for idx in range(n_samples):
        params_1 = rng.uniform(0.0, 2.0 * np.pi, size=param_shape)
        params_2 = rng.uniform(0.0, 2.0 * np.pi, size=param_shape)
        state_1 = np.asarray(circuit_fn(params_1))
        state_2 = np.asarray(circuit_fn(params_2))
        fidelities[idx] = float(np.clip(np.abs(np.vdot(state_1, state_2)) ** 2, 0.0, 1.0))
        entanglements[idx] = meyer_wallach_entanglement(state_1, n_qubits)

    bins = np.linspace(0.0, 1.0, n_bins + 1)
    pqc_hist, edges = np.histogram(fidelities, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dim = 2 ** n_qubits
    haar_pdf = (dim - 1) * np.power(np.clip(1.0 - centers, 1e-12, 1.0), dim - 2)

    pqc_prob = pqc_hist + 1e-12
    pqc_prob /= pqc_prob.sum()
    haar_prob = haar_pdf + 1e-12
    haar_prob /= haar_prob.sum()

    return {
        "kl_divergence": float(entropy(pqc_prob, haar_prob)),
        "entanglement_capability": float(np.mean(entanglements)),
        "fidelity_mean": float(np.mean(fidelities)),
        "fidelity_std": float(np.std(fidelities)),
        "haar_l1_distance": float(np.mean(np.abs(pqc_prob - haar_prob))),
        "samples": int(n_samples),
    }


def run_expressibility_benchmark(output_path=None, n_qubits=4, n_samples=2000):
    output_path = Path(output_path or Path(__file__).resolve().parents[1] / "results" / "expressibility_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    circuit_specs = {
        "HEA-1": {"depth": 1},
        "HEA-2": {"depth": 2},
        "StronglyEntangling": {"depth": 2},
        "Random": {"depth": 2},
        "QAOA-style": {"depth": 2},
    }

    results = {
        "n_qubits": n_qubits,
        "n_samples": n_samples,
        "circuits": {},
    }

    for name, spec in circuit_specs.items():
        circuit_fn, param_shape = _build_state_circuit(name, n_qubits, spec["depth"])
        metrics = compute_expressibility(circuit_fn, n_qubits, param_shape, n_samples=n_samples)
        metrics["param_shape"] = list(param_shape)
        metrics["depth"] = spec["depth"]
        results["circuits"][name] = metrics

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    return results


if __name__ == "__main__":
    run_expressibility_benchmark()
