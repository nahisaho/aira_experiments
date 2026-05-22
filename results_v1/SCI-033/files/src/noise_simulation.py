import json
from pathlib import Path

import numpy as np
import pennylane as qml
from sklearn.datasets import make_moons
from sklearn.metrics import accuracy_score
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

np.random.seed(42)

T1 = 100e-6
T2 = 150e-6
GATE_TIME = 200e-9
GATE_FIDELITY = 0.995


def _ibm_inspired_parameters():
    amplitude_damping = float(1.0 - np.exp(-GATE_TIME / T1))
    phase_flip = float(0.5 * (1.0 - np.exp(-GATE_TIME / T2)))
    depolarizing = float(1.0 - GATE_FIDELITY)
    bit_flip = depolarizing / 2.0
    return {
        "depolarizing": depolarizing,
        "bit_flip": bit_flip,
        "phase_flip": phase_flip,
        "amplitude_damping": amplitude_damping,
    }


def _feature_map(x, n_qubits):
    values = np.asarray(x, dtype=float)
    if values.size < n_qubits:
        values = np.pad(values, (0, n_qubits - values.size))
    values = values[:n_qubits]
    for wire, value in enumerate(values):
        qml.RY(value, wires=wire)
        qml.RZ(value / 2.0, wires=wire)
    for wire in range(n_qubits - 1):
        qml.IsingZZ(values[wire] * values[wire + 1], wires=[wire, wire + 1])
    if n_qubits > 2:
        qml.IsingZZ(values[0] * values[-1], wires=[0, n_qubits - 1])


def _apply_noise(channel, strength, n_qubits):
    if strength <= 0:
        return
    for wire in range(n_qubits):
        if channel == "depolarizing":
            qml.DepolarizingChannel(strength, wires=wire)
        elif channel == "bit_flip":
            qml.BitFlip(strength, wires=wire)
        elif channel == "phase_flip":
            qml.PhaseFlip(strength, wires=wire)
        elif channel == "amplitude_damping":
            qml.AmplitudeDamping(strength, wires=wire)


def _density_embeddings(X, n_qubits=4, channel="depolarizing", strength=0.0):
    dev = qml.device("default.mixed", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(x):
        _feature_map(x, n_qubits)
        _apply_noise(channel, strength, n_qubits)
        return qml.density_matrix(wires=range(n_qubits))

    scaler = StandardScaler()
    X_scaled = np.pi * np.tanh(scaler.fit_transform(X))
    return np.asarray([np.asarray(circuit(sample)) for sample in X_scaled])


def _kernel_from_density(rhos_a, rhos_b=None):
    rhos_b = rhos_a if rhos_b is None else rhos_b
    return np.real(np.einsum("aij,bji->ab", rhos_a, rhos_b))


def _evaluate_noisy_quantum_kernel(X, y, channel, strength, n_qubits=4):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    rho_train = _density_embeddings(X_train, n_qubits=n_qubits, channel=channel, strength=strength)
    rho_test = _density_embeddings(X_test, n_qubits=n_qubits, channel=channel, strength=strength)
    K_train = _kernel_from_density(rho_train)
    K_test = _kernel_from_density(rho_test, rho_train)
    clf = SVC(kernel="precomputed", C=1.0)
    clf.fit(K_train, y_train)
    pred = clf.predict(K_test)
    return float(accuracy_score(y_test, pred))


def run_noise_benchmark(output_path=None, n_qubits=4):
    output_path = Path(output_path or Path(__file__).resolve().parents[1] / "results" / "noise_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    X, y = make_moons(n_samples=70, noise=0.12, random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
    rbf_train = rbf_kernel(X_train, X_train, gamma=1.0)
    rbf_test = rbf_kernel(X_test, X_train, gamma=1.0)
    rbf_clf = SVC(kernel="precomputed", C=1.0)
    rbf_clf.fit(rbf_train, y_train)
    rbf_accuracy = float(accuracy_score(y_test, rbf_clf.predict(rbf_test)))

    baseline_params = _ibm_inspired_parameters()
    depolarizing_levels = [0.0, 0.0025, 0.005, 0.01, 0.02, 0.05]
    depolarizing_accuracies = []
    for level in depolarizing_levels:
        depolarizing_accuracies.append(_evaluate_noisy_quantum_kernel(X, y, "depolarizing", level, n_qubits=n_qubits))

    channel_comparison = {}
    for channel, strength in baseline_params.items():
        channel_comparison[channel] = {
            "strength": float(strength),
            "accuracy": _evaluate_noisy_quantum_kernel(X, y, channel, strength, n_qubits=n_qubits),
        }

    threshold = None
    for level, accuracy in zip(depolarizing_levels, depolarizing_accuracies):
        if accuracy <= rbf_accuracy:
            threshold = float(level)
            break
    if threshold is None:
        for level, accuracy in zip(depolarizing_levels, depolarizing_accuracies):
            if accuracy <= depolarizing_accuracies[0] - 0.05:
                threshold = float(level)
                break
    if threshold is None:
        threshold = float(depolarizing_levels[-1])

    results = {
        "n_qubits": n_qubits,
        "dataset": "moons",
        "rbf_accuracy": rbf_accuracy,
        "ibm_inspired_parameters": baseline_params,
        "depolarizing_scan": {
            "noise_levels": depolarizing_levels,
            "accuracies": depolarizing_accuracies,
            "threshold": threshold,
        },
        "channel_comparison": channel_comparison,
    }

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    return results


if __name__ == "__main__":
    run_noise_benchmark()
