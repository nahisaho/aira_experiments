"""
Experiment 3: Data Encoding Strategy Comparison
Compares angle, amplitude, and IQP encoding strategies for quantum classification.
Uses parameter-shift gradient computation with PennyLane.
"""
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.datasets import make_moons, make_circles
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import json

np.random.seed(42)

N_QUBITS = 4
N_LAYERS = 3
N_EPOCHS = 50
LR = 0.15
N_SAMPLES = 120

def prepare_data():
    datasets = {}
    scaler = MinMaxScaler(feature_range=(0, np.pi))
    
    X, y = make_moons(n_samples=N_SAMPLES, noise=0.15, random_state=42)
    X = scaler.fit_transform(X)
    X = np.hstack([X, np.zeros((len(X), 2))])
    datasets['Moons'] = train_test_split(X, y, test_size=0.3, random_state=42)
    
    X, y = make_circles(n_samples=N_SAMPLES, noise=0.1, factor=0.5, random_state=42)
    X = scaler.fit_transform(X)
    X = np.hstack([X, np.zeros((len(X), 2))])
    datasets['Circles'] = train_test_split(X, y, test_size=0.3, random_state=42)
    
    return datasets

def train_and_evaluate(encoding_type, X_train, X_test, y_train, y_test):
    dev = qml.device('default.qubit', wires=N_QUBITS)
    
    @qml.qnode(dev, diff_method='parameter-shift')
    def circuit(weights, x):
        # Encoding
        if encoding_type == 'Angle':
            for i in range(N_QUBITS):
                qml.RY(x[i], wires=i)
        elif encoding_type == 'Amplitude':
            qml.AmplitudeEmbedding(x, wires=range(N_QUBITS), normalize=True, pad_with=0.1)
        elif encoding_type == 'IQP':
            for i in range(N_QUBITS):
                qml.Hadamard(wires=i)
            for i in range(N_QUBITS):
                qml.RZ(x[i], wires=i)
            for i in range(N_QUBITS - 1):
                qml.IsingZZ(x[i] * x[i+1], wires=[i, i+1])
            for i in range(N_QUBITS):
                qml.Hadamard(wires=i)
        
        # Variational layers
        for l in range(N_LAYERS):
            for i in range(N_QUBITS):
                qml.RY(weights[l * N_QUBITS * 2 + i * 2], wires=i)
                qml.RZ(weights[l * N_QUBITS * 2 + i * 2 + 1], wires=i)
            for i in range(N_QUBITS - 1):
                qml.CNOT(wires=[i, i+1])
        return qml.expval(qml.PauliZ(0))
    
    n_params = N_LAYERS * N_QUBITS * 2
    weights = pnp.array(np.random.uniform(-np.pi, np.pi, n_params), requires_grad=True)
    y_train_scaled = pnp.array(2*y_train - 1, requires_grad=False)
    
    opt = qml.GradientDescentOptimizer(LR)
    losses = []
    
    def cost(weights):
        total = pnp.array(0.0)
        for i in range(len(X_train)):
            pred = circuit(weights, X_train[i])
            total = total + (pred - y_train_scaled[i])**2
        return total / len(X_train)
    
    for epoch in range(N_EPOCHS):
        weights = opt.step(cost, weights)
        if (epoch+1) % 10 == 0:
            loss_val = float(cost(weights))
            losses.append(loss_val)
            print(f"      Epoch {epoch+1}: loss={loss_val:.4f}")
    
    # Evaluate
    correct = 0
    for i in range(len(X_test)):
        pred = circuit(weights, X_test[i])
        if (pred > 0) == y_test[i]:
            correct += 1
    acc = correct / len(X_test)
    
    return acc, losses

print("Preparing data...")
datasets = prepare_data()

encodings = ['Angle', 'Amplitude', 'IQP']

results = {}
all_losses = {}

for ds_name, (X_train, X_test, y_train, y_test) in datasets.items():
    print(f"\nDataset: {ds_name}")
    results[ds_name] = {}
    all_losses[ds_name] = {}
    
    for enc_name in encodings:
        print(f"  Encoding: {enc_name}")
        acc, losses = train_and_evaluate(enc_name, X_train, X_test, y_train, y_test)
        results[ds_name][enc_name] = float(acc)
        all_losses[ds_name][enc_name] = losses
        print(f"    Test Accuracy: {acc:.4f}")

# Plot: Loss curves
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
colors = {'Angle': '#2196F3', 'Amplitude': '#4CAF50', 'IQP': '#FF9800'}

for i, ds_name in enumerate(datasets.keys()):
    for enc_name in encodings:
        if all_losses[ds_name][enc_name]:
            axes[i].plot(all_losses[ds_name][enc_name], 'o-', label=enc_name, color=colors[enc_name], linewidth=2)
    axes[i].set_title(f'{ds_name} - Training Loss', fontsize=13)
    axes[i].set_xlabel('Epoch (x10)')
    axes[i].set_ylabel('Loss')
    axes[i].legend()
    axes[i].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/encoding_loss_curves.png', dpi=150)
plt.close()

# Plot: Accuracy comparison
fig, ax = plt.subplots(figsize=(8, 5))
ds_names = list(results.keys())
x = np.arange(len(ds_names))
width = 0.25

for i, enc in enumerate(encodings):
    accs = [results[ds][enc] for ds in ds_names]
    ax.bar(x + i*width, accs, width, label=enc, color=list(colors.values())[i], alpha=0.85)

ax.set_xlabel('Dataset', fontsize=12)
ax.set_ylabel('Test Accuracy', fontsize=12)
ax.set_title('Data Encoding Strategy Comparison', fontsize=14)
ax.set_xticks(x + width)
ax.set_xticklabels(ds_names)
ax.legend()
ax.set_ylim(0.3, 1.05)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('figures/encoding_accuracy.png', dpi=150)
plt.close()

with open('results_encoding.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nExperiment 3 complete.")
print(json.dumps(results, indent=2))
