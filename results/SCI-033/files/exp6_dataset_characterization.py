"""
Experiment 6: Dataset Characterization for Quantum Advantage
Identifies features of datasets where quantum models outperform classical ones.
"""
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
import json

np.random.seed(42)

N_SAMPLES = 150

def generate_varied_datasets():
    """Generate datasets with different geometric/algebraic structures."""
    datasets = {}
    
    # 1. Linear separable
    X = np.random.randn(N_SAMPLES, 2)
    y = (X[:, 0] + X[:, 1] > 0).astype(int)
    datasets['Linear'] = (X, y, 'Simple linear boundary')
    
    # 2. XOR (parity-like)
    X = np.random.uniform(-1, 1, (N_SAMPLES, 2))
    y = ((X[:, 0] > 0) ^ (X[:, 1] > 0)).astype(int)
    datasets['XOR'] = (X, y, 'Parity structure')
    
    # 3. Concentric circles
    r1 = np.random.uniform(0, 0.5, N_SAMPLES//2)
    r2 = np.random.uniform(0.7, 1.2, N_SAMPLES - N_SAMPLES//2)
    theta = np.random.uniform(0, 2*np.pi, N_SAMPLES)
    r = np.concatenate([r1, r2])
    X = np.column_stack([r*np.cos(theta), r*np.sin(theta)])
    y = np.concatenate([np.zeros(N_SAMPLES//2), np.ones(N_SAMPLES - N_SAMPLES//2)]).astype(int)
    datasets['Circles'] = (X, y, 'Radial symmetry')
    
    # 4. Spiral (hard for local methods)
    t1 = np.linspace(0, 4*np.pi, N_SAMPLES//2) + np.random.randn(N_SAMPLES//2)*0.2
    t2 = np.linspace(0, 4*np.pi, N_SAMPLES - N_SAMPLES//2) + np.random.randn(N_SAMPLES - N_SAMPLES//2)*0.2
    X1 = np.column_stack([t1*np.cos(t1), t1*np.sin(t1)])
    X2 = np.column_stack([t2*np.cos(t2+np.pi), t2*np.sin(t2+np.pi)])
    X = np.vstack([X1, X2])
    y = np.concatenate([np.zeros(N_SAMPLES//2), np.ones(N_SAMPLES - N_SAMPLES//2)]).astype(int)
    datasets['Spiral'] = (X, y, 'Complex topology')
    
    # 5. Checkerboard
    X = np.random.uniform(-1, 1, (N_SAMPLES, 2))
    y = (((X[:, 0]*3).astype(int) + (X[:, 1]*3).astype(int)) % 2).astype(int)
    datasets['Checkerboard'] = (X, y, 'High-frequency pattern')
    
    # 6. Gaussian clusters (easy)
    n_half = N_SAMPLES // 2
    X1 = np.random.randn(n_half, 2) * 0.5 + np.array([1, 1])
    X2 = np.random.randn(N_SAMPLES - n_half, 2) * 0.5 + np.array([-1, -1])
    X = np.vstack([X1, X2])
    y = np.concatenate([np.zeros(n_half), np.ones(N_SAMPLES - n_half)]).astype(int)
    datasets['Gaussian'] = (X, y, 'Well-separated clusters')
    
    return datasets

def quantum_classifier_accuracy(X, y, n_qubits=2, n_layers=2, n_epochs=30):
    """Train a simple VQC and return accuracy."""
    scaler = MinMaxScaler(feature_range=(0, np.pi))
    X_scaled = scaler.fit_transform(X)
    
    dev = qml.device('default.qubit', wires=n_qubits)
    
    @qml.qnode(dev, diff_method='parameter-shift')
    def circuit(weights, x):
        for i in range(n_qubits):
            qml.RY(x[i % X.shape[1]], wires=i)
        for l in range(n_layers):
            for i in range(n_qubits):
                qml.RY(weights[l * n_qubits * 2 + i * 2], wires=i)
                qml.RZ(weights[l * n_qubits * 2 + i * 2 + 1], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i+1])
        return qml.expval(qml.PauliZ(0))
    
    n_params = n_layers * n_qubits * 2
    weights = pnp.array(np.random.uniform(-np.pi, np.pi, n_params), requires_grad=True)
    
    split = int(0.7 * len(X_scaled))
    X_train, X_test = X_scaled[:split], X_scaled[split:]
    y_train_labels, y_test_labels = y[:split], y[split:]
    y_train_scaled = pnp.array(2*y_train_labels - 1, requires_grad=False)
    
    opt = qml.GradientDescentOptimizer(0.15)
    
    def cost(weights):
        total = pnp.array(0.0)
        for i in range(len(X_train)):
            pred = circuit(weights, X_train[i])
            total = total + (pred - y_train_scaled[i])**2
        return total / len(X_train)
    
    for _ in range(n_epochs):
        weights = opt.step(cost, weights)
    
    correct = 0
    for i in range(len(X_test)):
        pred = circuit(weights, X_test[i])
        if (pred > 0) == y_test_labels[i]:
            correct += 1
    return float(correct / len(X_test))

datasets = generate_varied_datasets()
results = {}

for ds_name, (X, y, desc) in datasets.items():
    print(f"\nDataset: {ds_name} ({desc})")
    
    # Classical models
    scaler = MinMaxScaler()
    X_s = scaler.fit_transform(X)
    
    svm_rbf = cross_val_score(SVC(kernel='rbf'), X_s, y, cv=5).mean()
    svm_lin = cross_val_score(SVC(kernel='linear'), X_s, y, cv=5).mean()
    mlp = cross_val_score(MLPClassifier(hidden_layer_sizes=(16, 8), max_iter=500, random_state=42), X_s, y, cv=5).mean()
    
    # Quantum
    qml_acc = quantum_classifier_accuracy(X, y)
    
    # Dataset complexity metrics
    # Geometric complexity: average distance ratio (inter/intra class)
    X0, X1 = X[y==0], X[y==1]
    intra_0 = np.mean(np.linalg.norm(X0 - X0.mean(axis=0), axis=1))
    intra_1 = np.mean(np.linalg.norm(X1 - X1.mean(axis=0), axis=1))
    inter = np.linalg.norm(X0.mean(axis=0) - X1.mean(axis=0))
    separability = inter / (intra_0 + intra_1 + 1e-10)
    
    results[ds_name] = {
        'description': desc,
        'SVM-Linear': float(svm_lin),
        'SVM-RBF': float(svm_rbf),
        'MLP': float(mlp),
        'Quantum-VQC': float(qml_acc),
        'separability_ratio': float(separability),
        'quantum_advantage': float(qml_acc - max(svm_rbf, mlp)),
    }
    print(f"  SVM-Lin={svm_lin:.3f}, SVM-RBF={svm_rbf:.3f}, MLP={mlp:.3f}, QML={qml_acc:.3f}")
    print(f"  Separability: {separability:.3f}, Q-advantage: {qml_acc - max(svm_rbf, mlp):.3f}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

ds_names = list(results.keys())
methods = ['SVM-Linear', 'SVM-RBF', 'MLP', 'Quantum-VQC']
colors = ['#9E9E9E', '#2196F3', '#4CAF50', '#FF9800']
x = np.arange(len(ds_names))
width = 0.2

for i, method in enumerate(methods):
    accs = [results[ds][method] for ds in ds_names]
    axes[0].bar(x + i*width, accs, width, label=method, color=colors[i], alpha=0.85)

axes[0].set_xlabel('Dataset', fontsize=12)
axes[0].set_ylabel('Accuracy', fontsize=12)
axes[0].set_title('Classification Accuracy by Method & Dataset', fontsize=13)
axes[0].set_xticks(x + width*1.5)
axes[0].set_xticklabels(ds_names, rotation=30, ha='right')
axes[0].legend(fontsize=9)
axes[0].set_ylim(0.3, 1.05)
axes[0].grid(True, alpha=0.3, axis='y')

# Quantum advantage vs dataset complexity
q_adv = [results[ds]['quantum_advantage'] for ds in ds_names]
sep = [results[ds]['separability_ratio'] for ds in ds_names]
axes[1].scatter(sep, q_adv, s=120, c='#FF9800', edgecolors='black', zorder=5)
for i, ds in enumerate(ds_names):
    axes[1].annotate(ds, (sep[i], q_adv[i]), textcoords="offset points", xytext=(5, 5), fontsize=9)
axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
axes[1].set_xlabel('Separability Ratio', fontsize=12)
axes[1].set_ylabel('Quantum Advantage (QML - best classical)', fontsize=12)
axes[1].set_title('Quantum Advantage vs Dataset Complexity', fontsize=13)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/dataset_characterization.png', dpi=150)
plt.close()

with open('results_dataset.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nExperiment 6 complete.")
print(json.dumps(results, indent=2))
