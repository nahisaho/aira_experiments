import json
import time
from pathlib import Path

import numpy as np
import pennylane as qml
from sklearn.datasets import make_circles, make_moons
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

np.random.seed(42)


def _xor_dataset(n_samples=72, seed=42):
    rng = np.random.default_rng(seed)
    X = rng.uniform(-1.0, 1.0, size=(n_samples, 2))
    y = (X[:, 0] * X[:, 1] > 0).astype(int)
    X += rng.normal(scale=0.08, size=X.shape)
    return X, y


def create_encoding_datasets(seed=42):
    circles_X, circles_y = make_circles(n_samples=72, noise=0.08, factor=0.4, random_state=seed)
    moons_X, moons_y = make_moons(n_samples=72, noise=0.12, random_state=seed)
    xor_X, xor_y = _xor_dataset(seed=seed)
    return {
        "xor": (xor_X, xor_y),
        "circles": (circles_X, circles_y),
        "moons": (moons_X, moons_y),
    }


def _prepare_input(x, n_qubits):
    x = np.asarray(x, dtype=float)
    if x.size < n_qubits:
        x = np.pad(x, (0, n_qubits - x.size))
    return x[:n_qubits]


def _amplitude_vector(x, n_qubits):
    target_len = 2 ** n_qubits
    x = np.asarray(x, dtype=float)
    if x.size < target_len:
        x = np.pad(x, (0, target_len - x.size))
    else:
        x = x[:target_len]
    if np.linalg.norm(x) < 1e-10:
        x[0] = 1.0
    return x / np.linalg.norm(x)


def _angle_state_embeddings(X, n_qubits=2):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(x):
        qml.AngleEmbedding(_prepare_input(x, n_qubits), wires=range(n_qubits), rotation="Y")
        for wire in range(n_qubits - 1):
            qml.CNOT(wires=[wire, wire + 1])
        return qml.state()

    return np.asarray([np.asarray(circuit(sample)) for sample in X])


def _amplitude_state_embeddings(X, n_qubits=2):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(x):
        qml.AmplitudeEmbedding(_amplitude_vector(x, n_qubits), wires=range(n_qubits), normalize=True)
        return qml.state()

    return np.asarray([np.asarray(circuit(sample)) for sample in X])


def _iqp_state_embeddings(X, n_qubits=2):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(x):
        qml.IQPEmbedding(_prepare_input(x, n_qubits), wires=range(n_qubits), n_repeats=1)
        return qml.state()

    return np.asarray([np.asarray(circuit(sample)) for sample in X])


def _trainability_proxy(encoding_name, samples, n_qubits=2, seed=42):
    rng = np.random.default_rng(seed)
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(theta, x):
        if encoding_name == "angle":
            qml.AngleEmbedding(_prepare_input(x, n_qubits), wires=range(n_qubits), rotation="Y")
        elif encoding_name == "amplitude":
            qml.AmplitudeEmbedding(_amplitude_vector(x, n_qubits), wires=range(n_qubits), normalize=True)
        else:
            qml.IQPEmbedding(_prepare_input(x, n_qubits), wires=range(n_qubits), n_repeats=1)
        qml.RY(theta, wires=0)
        return qml.expval(qml.PauliZ(0))

    eps = 1e-3
    grads = []
    for x in samples[:20]:
        theta = rng.uniform(0.0, 2.0 * np.pi)
        forward = circuit(theta + eps, x)
        backward = circuit(theta - eps, x)
        grads.append(float((forward - backward) / (2 * eps)))
    grads = np.asarray(grads)
    return {
        "mean_abs_gradient": float(np.mean(np.abs(grads))),
        "gradient_variance": float(np.var(grads)),
    }


def _kernel_from_states(states_a, states_b=None):
    states_b = states_a if states_b is None else states_b
    return np.abs(states_a @ states_b.conj().T) ** 2


def _expressiveness(states):
    if len(states) < 2:
        return 0.0
    kernel = _kernel_from_states(states)
    mask = ~np.eye(kernel.shape[0], dtype=bool)
    fidelities = kernel[mask]
    return float(np.mean(1.0 - fidelities))


def _benchmark_encoding(encoding_name, X, y, n_qubits=2):
    scaler = StandardScaler()
    X_scaled = np.pi * np.tanh(scaler.fit_transform(X))
    split = train_test_split(X_scaled, y, test_size=0.3, random_state=42, stratify=y)
    X_train, X_test, y_train, y_test = split

    if encoding_name == "angle":
        embed_fn = _angle_state_embeddings
    elif encoding_name == "amplitude":
        embed_fn = _amplitude_state_embeddings
    else:
        embed_fn = _iqp_state_embeddings

    start = time.perf_counter()
    train_states = embed_fn(X_train, n_qubits=n_qubits)
    test_states = embed_fn(X_test, n_qubits=n_qubits)
    elapsed = time.perf_counter() - start

    K_train = _kernel_from_states(train_states)
    K_test = _kernel_from_states(test_states, train_states)
    clf = SVC(kernel="precomputed", C=1.0)
    clf.fit(K_train, y_train)
    pred = clf.predict(K_test)

    combined_states = np.vstack([train_states[:20], test_states[:20]])
    trainability = _trainability_proxy(encoding_name, X_scaled[:20], n_qubits=n_qubits)

    return {
        "accuracy": float(accuracy_score(y_test, pred)),
        "expressiveness": _expressiveness(combined_states),
        "trainability": trainability,
        "computational_cost_ms": float(1000.0 * elapsed / (len(X_train) + len(X_test))),
        "train_size": int(len(X_train)),
        "test_size": int(len(X_test)),
    }


def run_encoding_benchmark(output_path=None, n_qubits=2):
    output_path = Path(output_path or Path(__file__).resolve().parents[1] / "results" / "encoding_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    datasets = create_encoding_datasets()
    encodings = ["angle", "amplitude", "iqp"]
    results = {"n_qubits": n_qubits, "datasets": {}, "summary": {}}

    summary_rows = {encoding: {"accuracy": [], "expressiveness": [], "gradient": [], "cost": []} for encoding in encodings}

    for dataset_name, (X, y) in datasets.items():
        results["datasets"][dataset_name] = {}
        for encoding in encodings:
            metrics = _benchmark_encoding(encoding, X, y, n_qubits=n_qubits)
            results["datasets"][dataset_name][encoding] = metrics
            summary_rows[encoding]["accuracy"].append(metrics["accuracy"])
            summary_rows[encoding]["expressiveness"].append(metrics["expressiveness"])
            summary_rows[encoding]["gradient"].append(metrics["trainability"]["mean_abs_gradient"])
            summary_rows[encoding]["cost"].append(metrics["computational_cost_ms"])

    for encoding, values in summary_rows.items():
        results["summary"][encoding] = {
            "mean_accuracy": float(np.mean(values["accuracy"])),
            "mean_expressiveness": float(np.mean(values["expressiveness"])),
            "mean_abs_gradient": float(np.mean(values["gradient"])),
            "mean_cost_ms": float(np.mean(values["cost"])),
        }

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    return results


if __name__ == "__main__":
    run_encoding_benchmark()
