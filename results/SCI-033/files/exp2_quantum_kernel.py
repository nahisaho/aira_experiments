"""
Experiment 2: Quantum Kernel Methods vs Classical Kernels
Compares quantum kernel (angle, IQP encoding) with classical kernels (RBF, polynomial)
on classification tasks.
"""
import numpy as np
import pennylane as qml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.datasets import make_moons, make_circles
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score
import json

np.random.seed(42)

N_QUBITS = 4
N_SAMPLES = 200
N_LAYERS = 2

def generate_datasets():
    """Generate benchmark datasets."""
    datasets = {}
    # Moons
    X, y = make_moons(n_samples=N_SAMPLES, noise=0.15, random_state=42)
    scaler = MinMaxScaler(feature_range=(0, np.pi))
    X = scaler.fit_transform(X)
    datasets['Moons'] = (X, y)
    
    # Circles
    X, y = make_circles(n_samples=N_SAMPLES, noise=0.1, factor=0.5, random_state=42)
    X = scaler.fit_transform(X)
    datasets['Circles'] = (X, y)
    
    # XOR-like
    X = np.random.uniform(0, np.pi, (N_SAMPLES, 2))
    y = ((X[:, 0] > np.pi/2) ^ (X[:, 1] > np.pi/2)).astype(int)
    datasets['XOR'] = (X, y)
    
    # Checkerboard (harder for linear methods)
    X = np.random.uniform(0, np.pi, (N_SAMPLES, 2))
    grid = 3
    y = (((X[:, 0] * grid / np.pi).astype(int) + (X[:, 1] * grid / np.pi).astype(int)) % 2).astype(int)
    datasets['Checkerboard'] = (X, y)
    
    return datasets

def quantum_kernel_angle(x1, x2):
    """Angle encoding quantum kernel."""
    n_features = len(x1)
    n_qubits = n_features
    dev = qml.device('default.qubit', wires=n_qubits)
    
    @qml.qnode(dev)
    def kernel_circuit(x1, x2):
        for i in range(n_qubits):
            qml.RY(x1[i % n_features], wires=i)
            qml.RZ(x1[i % n_features], wires=i)
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i+1])
        # Adjoint
        for i in range(n_qubits - 1, 0, -1):
            qml.CNOT(wires=[i-1, i])
        for i in range(n_qubits - 1, -1, -1):
            qml.adjoint(qml.RZ)(x2[i % n_features], wires=i)
            qml.adjoint(qml.RY)(x2[i % n_features], wires=i)
        return qml.probs(wires=range(n_qubits))
    
    return kernel_circuit(x1, x2)[0]

def quantum_kernel_iqp(x1, x2):
    """IQP encoding quantum kernel."""
    n_features = len(x1)
    n_qubits = n_features
    dev = qml.device('default.qubit', wires=n_qubits)
    
    @qml.qnode(dev)
    def kernel_circuit(x1, x2):
        # IQP encoding for x1
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        for i in range(n_qubits):
            qml.RZ(x1[i % n_features], wires=i)
        for i in range(n_qubits - 1):
            qml.IsingZZ(x1[i % n_features] * x1[(i+1) % n_features], wires=[i, i+1])
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        # Adjoint IQP for x2
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        for i in range(n_qubits - 2, -1, -1):
            qml.IsingZZ(-x2[i % n_features] * x2[(i+1) % n_features], wires=[i, i+1])
        for i in range(n_qubits - 1, -1, -1):
            qml.RZ(-x2[i % n_features], wires=i)
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        return qml.probs(wires=range(n_qubits))
    
    return kernel_circuit(x1, x2)[0]

def compute_kernel_matrix(X, kernel_fn):
    """Compute full kernel matrix."""
    n = len(X)
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            k = kernel_fn(X[i], X[j])
            K[i, j] = k
            K[j, i] = k
    return K

print("Generating datasets...")
datasets = generate_datasets()

print("Computing quantum kernel matrices (this takes a while)...")
results = {}

for ds_name, (X, y) in datasets.items():
    print(f"\n  Dataset: {ds_name}")
    ds_results = {}
    
    # Classical kernels
    for kernel_name in ['rbf', 'poly']:
        svc = SVC(kernel=kernel_name, random_state=42)
        scores = cross_val_score(svc, X, y, cv=5, scoring='accuracy')
        ds_results[f'Classical-{kernel_name.upper()}'] = {
            'accuracy_mean': float(scores.mean()),
            'accuracy_std': float(scores.std())
        }
        print(f"    Classical {kernel_name.upper()}: {scores.mean():.4f}±{scores.std():.4f}")
    
    # Quantum kernels (use subset for speed)
    subset_size = min(80, len(X))
    indices = np.random.choice(len(X), subset_size, replace=False)
    X_sub, y_sub = X[indices], y[indices]
    
    for qk_name, qk_fn in [('Angle', quantum_kernel_angle), ('IQP', quantum_kernel_iqp)]:
        print(f"    Computing {qk_name} quantum kernel...")
        K = compute_kernel_matrix(X_sub, qk_fn)
        # Ensure PSD
        K = (K + K.T) / 2
        eigvals = np.linalg.eigvalsh(K)
        if eigvals.min() < 0:
            K += (-eigvals.min() + 1e-6) * np.eye(len(K))
        
        svc = SVC(kernel='precomputed', random_state=42)
        # Simple train/test split
        split = int(0.7 * len(X_sub))
        K_train = K[:split, :split]
        K_test = K[split:, :split]
        
        svc.fit(K_train, y_sub[:split])
        y_pred = svc.predict(K_test)
        acc = accuracy_score(y_sub[split:], y_pred)
        
        ds_results[f'Quantum-{qk_name}'] = {
            'accuracy_mean': float(acc),
            'accuracy_std': 0.0
        }
        print(f"    Quantum {qk_name}: {acc:.4f}")
    
    results[ds_name] = ds_results

# Plot: Accuracy comparison
fig, ax = plt.subplots(figsize=(12, 6))
ds_names = list(results.keys())
methods = ['Classical-RBF', 'Classical-POLY', 'Quantum-Angle', 'Quantum-IQP']
x = np.arange(len(ds_names))
width = 0.2

colors = ['#2196F3', '#4CAF50', '#FF9800', '#E91E63']
for i, method in enumerate(methods):
    accs = [results[ds].get(method, {}).get('accuracy_mean', 0) for ds in ds_names]
    stds = [results[ds].get(method, {}).get('accuracy_std', 0) for ds in ds_names]
    ax.bar(x + i*width, accs, width, yerr=stds, label=method, color=colors[i], alpha=0.85, capsize=3)

ax.set_xlabel('Dataset', fontsize=12)
ax.set_ylabel('Accuracy', fontsize=12)
ax.set_title('Quantum vs Classical Kernel Classification Accuracy', fontsize=14)
ax.set_xticks(x + width*1.5)
ax.set_xticklabels(ds_names)
ax.legend(fontsize=10)
ax.set_ylim(0.3, 1.05)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('figures/kernel_comparison.png', dpi=150)
plt.close()

with open('results_kernel.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nExperiment 2 complete.")
print(json.dumps(results, indent=2))
