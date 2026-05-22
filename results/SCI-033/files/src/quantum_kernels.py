import json
from pathlib import Path

import numpy as np
import pennylane as qml
from sklearn.datasets import make_circles, make_moons
from sklearn.metrics import accuracy_score
from sklearn.metrics.pairwise import rbf_kernel, polynomial_kernel
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

np.random.seed(42)


def _xor_dataset(n_samples=80, seed=42):
    rng = np.random.default_rng(seed)
    X = rng.uniform(-1.0, 1.0, size=(n_samples, 2))
    y = (X[:, 0] * X[:, 1] > 0).astype(int)
    X += rng.normal(scale=0.08, size=X.shape)
    return X, y


def create_datasets(seed=42):
    xor_X, xor_y = _xor_dataset(seed=seed)
    circles_X, circles_y = make_circles(n_samples=80, noise=0.08, factor=0.35, random_state=seed)
    moons_X, moons_y = make_moons(n_samples=80, noise=0.12, random_state=seed)
    return {
        "xor": (xor_X, xor_y),
        "circles": (circles_X, circles_y),
        "moons": (moons_X, moons_y),
    }


def _scale_features(X, n_qubits=4):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    if X_scaled.shape[1] < n_qubits:
        repeats = int(np.ceil(n_qubits / X_scaled.shape[1]))
        X_scaled = np.tile(X_scaled, (1, repeats))[:, :n_qubits]
    elif X_scaled.shape[1] > n_qubits:
        X_scaled = X_scaled[:, :n_qubits]
    X_scaled = np.pi * np.tanh(X_scaled)
    return X_scaled


def _feature_map(x, wires):
    for wire, value in enumerate(x):
        qml.RY(value, wires=wire)
        qml.RZ(0.5 * value, wires=wire)
    for wire in range(len(wires) - 1):
        qml.IsingZZ(x[wire] * x[wire + 1], wires=[wire, wire + 1])
    if len(wires) > 2:
        qml.IsingZZ(x[0] * x[-1], wires=[0, len(wires) - 1])


def compute_state_embeddings(X, n_qubits=4):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(x):
        _feature_map(x, list(range(n_qubits)))
        return qml.state()

    X_proc = _scale_features(X, n_qubits=n_qubits)
    return np.asarray([np.asarray(circuit(sample)) for sample in X_proc])


def compute_projected_features(X, n_qubits=4):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(x):
        _feature_map(x, list(range(n_qubits)))
        return [qml.expval(qml.PauliZ(wire)) for wire in range(n_qubits)]

    X_proc = _scale_features(X, n_qubits=n_qubits)
    return np.asarray([np.asarray(circuit(sample), dtype=float) for sample in X_proc])


def quantum_kernel(X1, X2=None, n_qubits=4):
    states_1 = compute_state_embeddings(X1, n_qubits=n_qubits)
    states_2 = states_1 if X2 is None else compute_state_embeddings(X2, n_qubits=n_qubits)
    return np.abs(states_1 @ states_2.conj().T) ** 2


def projected_quantum_kernel(X1, X2=None, n_qubits=4):
    features_1 = compute_projected_features(X1, n_qubits=n_qubits)
    features_2 = features_1 if X2 is None else compute_projected_features(X2, n_qubits=n_qubits)
    return (features_1 @ features_2.T) / n_qubits


def _center_kernel(K):
    n = K.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    return H @ K @ H


def kernel_alignment(K1, K2):
    A = _center_kernel(K1)
    B = _center_kernel(K2)
    numerator = float(np.sum(A * B))
    denominator = float(np.sqrt(np.sum(A * A) * np.sum(B * B)) + 1e-12)
    return numerator / denominator


def kernel_target_alignment(K, y):
    y_vec = np.where(np.asarray(y) > 0, 1.0, -1.0)
    yy = np.outer(y_vec, y_vec)
    numerator = float(np.sum(K * yy))
    denominator = float(np.sqrt(np.sum(K * K) * np.sum(yy * yy)) + 1e-12)
    return numerator / denominator


def analyze_rkhs_capacity(K):
    centered = _center_kernel(K)
    eigvals = np.linalg.eigvalsh(centered)
    eigvals = np.clip(eigvals, 0.0, None)
    total = float(np.sum(eigvals))
    if total <= 1e-12:
        return {
            "effective_dimension": 0.0,
            "spectral_entropy": 0.0,
            "normalized_rank": 0.0,
            "high_capacity": False,
        }
    probs = eigvals / total
    spectral_entropy = float(-np.sum(np.where(probs > 0, probs * np.log(probs), 0.0)))
    effective_dimension = float(total ** 2 / (np.sum(eigvals ** 2) + 1e-12))
    normalized_rank = float(np.exp(spectral_entropy) / len(eigvals))
    return {
        "effective_dimension": effective_dimension,
        "spectral_entropy": spectral_entropy,
        "normalized_rank": normalized_rank,
        "high_capacity": bool(normalized_rank > 0.35 and effective_dimension > 4.0),
    }


def _svm_scores(X, y, n_qubits=4):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    Kq_train = quantum_kernel(X_train, n_qubits=n_qubits)
    Kq_test = quantum_kernel(X_test, X_train, n_qubits=n_qubits)
    Kpq_train = projected_quantum_kernel(X_train, n_qubits=n_qubits)
    Kpq_test = projected_quantum_kernel(X_test, X_train, n_qubits=n_qubits)

    rbf_train = rbf_kernel(X_train, X_train, gamma=1.0)
    rbf_test = rbf_kernel(X_test, X_train, gamma=1.0)
    poly_train = polynomial_kernel(X_train, X_train, degree=3, gamma=1.0)
    poly_test = polynomial_kernel(X_test, X_train, degree=3, gamma=1.0)

    models = {
        "quantum_kernel": (Kq_train, Kq_test),
        "projected_quantum_kernel": (Kpq_train, Kpq_test),
        "rbf": (rbf_train, rbf_test),
        "polynomial": (poly_train, poly_test),
    }

    accuracies = {}
    for name, (train_kernel, test_kernel) in models.items():
        clf = SVC(kernel="precomputed", C=1.0)
        clf.fit(train_kernel, y_train)
        pred = clf.predict(test_kernel)
        accuracies[name] = float(accuracy_score(y_test, pred))

    return {
        "accuracies": accuracies,
        "kernel_alignment_with_rbf": float(kernel_alignment(Kq_train, rbf_train)),
        "projected_alignment_with_rbf": float(kernel_alignment(Kpq_train, rbf_train)),
        "kernel_target_alignment": float(kernel_target_alignment(Kq_train, y_train)),
        "projected_kernel_target_alignment": float(kernel_target_alignment(Kpq_train, y_train)),
        "rkhs_capacity": analyze_rkhs_capacity(Kq_train),
        "train_size": int(len(X_train)),
        "test_size": int(len(X_test)),
    }


def run_kernel_benchmark(output_path=None, n_qubits=4):
    output_path = Path(output_path or Path(__file__).resolve().parents[1] / "results" / "kernel_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    datasets = create_datasets()
    results = {
        "n_qubits": n_qubits,
        "datasets": {},
        "summary": {},
    }

    quantum_acc = []
    rbf_acc = []
    kta_values = []
    alignments = []
    high_capacity = []

    for name, (X, y) in datasets.items():
        metrics = _svm_scores(X, y, n_qubits=n_qubits)
        results["datasets"][name] = metrics
        quantum_acc.append(metrics["accuracies"]["quantum_kernel"])
        rbf_acc.append(metrics["accuracies"]["rbf"])
        kta_values.append(metrics["kernel_target_alignment"])
        alignments.append(metrics["kernel_alignment_with_rbf"])
        high_capacity.append(metrics["rkhs_capacity"]["high_capacity"])

    results["summary"] = {
        "mean_quantum_accuracy": float(np.mean(quantum_acc)),
        "mean_rbf_accuracy": float(np.mean(rbf_acc)),
        "mean_quantum_advantage": float(np.mean(np.asarray(quantum_acc) - np.asarray(rbf_acc))),
        "mean_kernel_target_alignment": float(np.mean(kta_values)),
        "mean_kernel_alignment_with_rbf": float(np.mean(alignments)),
        "high_capacity_fraction": float(np.mean(high_capacity)),
    }

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    return results


if __name__ == "__main__":
    run_kernel_benchmark()
