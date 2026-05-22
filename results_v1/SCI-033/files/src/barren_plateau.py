import json
from pathlib import Path

import numpy as np
import pennylane as qml

np.random.seed(42)


def _ansatz(theta, rest_params, n_qubits, depth, noise_strength=0.0):
    index = 0
    for layer in range(depth):
        for wire in range(n_qubits):
            angles = []
            for pos in range(3):
                if layer == 0 and wire == 0 and pos == 0:
                    angles.append(theta)
                else:
                    angles.append(rest_params[index])
                    index += 1
            qml.RY(angles[0], wires=wire)
            qml.RZ(angles[1], wires=wire)
            qml.RX(angles[2], wires=wire)
        for wire in range(n_qubits - 1):
            qml.CNOT(wires=[wire, wire + 1])
        if n_qubits > 2:
            qml.CNOT(wires=[n_qubits - 1, 0])
        if noise_strength > 0:
            for wire in range(n_qubits):
                qml.DepolarizingChannel(noise_strength, wires=wire)


def _observable(n_qubits, cost_type):
    if cost_type == "local":
        return qml.PauliZ(0)
    obs = qml.PauliZ(0)
    for wire in range(1, n_qubits):
        obs = obs @ qml.PauliZ(wire)
    return obs


def _gradient_samples(n_qubits, depth=3, cost_type="global", n_samples=24, noise_strength=0.0, seed=42):
    rng = np.random.default_rng(seed + n_qubits)
    dev_name = "default.mixed" if noise_strength > 0 else "default.qubit"
    dev = qml.device(dev_name, wires=n_qubits)
    obs = _observable(n_qubits, cost_type)
    rest_size = depth * n_qubits * 3 - 1

    @qml.qnode(dev)
    def cost_fn(theta, rest_params):
        _ansatz(theta, rest_params, n_qubits=n_qubits, depth=depth, noise_strength=noise_strength)
        return qml.expval(obs)

    eps = 1e-3
    grads = []
    for _ in range(n_samples):
        theta = rng.uniform(0.0, 2.0 * np.pi)
        rest = rng.uniform(0.0, 2.0 * np.pi, size=rest_size)
        forward = cost_fn(theta + eps, rest)
        backward = cost_fn(theta - eps, rest)
        grads.append(float((forward - backward) / (2.0 * eps)))
    return np.asarray(grads)


def run_barren_plateau_benchmark(output_path=None, qubit_list=None, depth=3):
    output_path = Path(output_path or Path(__file__).resolve().parents[1] / "results" / "barren_plateau_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    qubit_list = qubit_list or [2, 4, 6, 8, 10]

    global_variances = []
    local_variances = []
    by_qubits = {}

    for n_qubits in qubit_list:
        global_grads = _gradient_samples(n_qubits, depth=depth, cost_type="global", n_samples=20)
        local_grads = _gradient_samples(n_qubits, depth=depth, cost_type="local", n_samples=20)
        global_var = float(np.var(global_grads))
        local_var = float(np.var(local_grads))
        by_qubits[str(n_qubits)] = {
            "global_gradient_variance": global_var,
            "local_gradient_variance": local_var,
            "mean_abs_global_gradient": float(np.mean(np.abs(global_grads))),
            "mean_abs_local_gradient": float(np.mean(np.abs(local_grads))),
        }
        global_variances.append(global_var)
        local_variances.append(local_var)

    noise_levels = [0.0, 0.002, 0.005, 0.01, 0.02]
    noise_variances = {}
    for level in noise_levels:
        grads = _gradient_samples(6, depth=depth, cost_type="global", n_samples=16, noise_strength=level, seed=84)
        noise_variances[str(level)] = float(np.var(grads))

    fit_qubits = np.asarray(qubit_list, dtype=float)
    global_fit = np.polyfit(fit_qubits, np.log10(np.asarray(global_variances) + 1e-16), deg=1)
    local_fit = np.polyfit(np.log10(fit_qubits), np.log10(np.asarray(local_variances) + 1e-16), deg=1)

    results = {
        "depth": depth,
        "qubit_list": qubit_list,
        "by_qubits": by_qubits,
        "scaling": {
            "global_log10_slope_per_qubit": float(global_fit[0]),
            "local_loglog_slope": float(local_fit[0]),
            "global_expected_exponential": "O(1/2^n)",
            "local_expected_scaling": "poly(1/n)",
        },
        "noise_induced_plateau": {
            "noise_levels": noise_levels,
            "global_gradient_variance": noise_variances,
        },
    }

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    return results


if __name__ == "__main__":
    run_barren_plateau_benchmark()
