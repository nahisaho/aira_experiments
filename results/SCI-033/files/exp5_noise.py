"""
Experiment 5: Noise Impact Analysis (Simulated IBM Quantum Noise)
Evaluates quantum model performance under depolarizing noise at various levels.
"""
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.datasets import make_moons
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import json

np.random.seed(42)

N_QUBITS = 4
N_LAYERS = 2
N_EPOCHS = 30
LR = 0.15
N_SAMPLES = 100

# Prepare data
X, y = make_moons(n_samples=N_SAMPLES, noise=0.15, random_state=42)
scaler = MinMaxScaler(feature_range=(0, np.pi))
X = scaler.fit_transform(X)
X = np.hstack([X, np.zeros((len(X), 2))])
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

noise_levels = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]

def create_noisy_circuit(noise_prob):
    dev = qml.device('default.mixed', wires=N_QUBITS)
    
    @qml.qnode(dev, diff_method='parameter-shift')
    def circuit(weights, x):
        for i in range(N_QUBITS):
            qml.RY(x[i], wires=i)
            if noise_prob > 0:
                qml.DepolarizingChannel(noise_prob, wires=i)
        for l in range(N_LAYERS):
            for i in range(N_QUBITS):
                qml.RY(weights[l * N_QUBITS * 2 + i * 2], wires=i)
                qml.RZ(weights[l * N_QUBITS * 2 + i * 2 + 1], wires=i)
                if noise_prob > 0:
                    qml.DepolarizingChannel(noise_prob, wires=i)
            for i in range(N_QUBITS - 1):
                qml.CNOT(wires=[i, i+1])
                if noise_prob > 0:
                    qml.DepolarizingChannel(noise_prob, wires=i)
                    qml.DepolarizingChannel(noise_prob, wires=i+1)
        return qml.expval(qml.PauliZ(0))
    return circuit

results = {}
n_params = N_LAYERS * N_QUBITS * 2

for noise in noise_levels:
    print(f"\nNoise level: {noise}")
    circuit = create_noisy_circuit(noise)
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
            print(f"  Epoch {epoch+1}: loss={loss_val:.4f}")
    
    # Evaluate
    correct = 0
    for i in range(len(X_test)):
        pred = circuit(weights, X_test[i])
        if (pred > 0) == y_test[i]:
            correct += 1
    acc = correct / len(X_test)
    
    # Gradient magnitude under noise
    grad_mags = []
    for _ in range(20):
        params = pnp.array(np.random.uniform(-np.pi, np.pi, n_params), requires_grad=True)
        x_sample = X_train[np.random.randint(len(X_train))]
        grad = qml.grad(circuit)(params, x_sample)
        grad_mags.append(float(np.mean(np.abs(grad))))
    
    results[str(noise)] = {
        'accuracy': float(acc),
        'final_loss': float(losses[-1]) if losses else 0.0,
        'losses': losses,
        'mean_grad_magnitude': float(np.mean(grad_mags)),
        'std_grad_magnitude': float(np.std(grad_mags)),
    }
    print(f"  Accuracy: {acc:.4f}, Mean |grad|: {np.mean(grad_mags):.6f}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

noise_vals = [float(n) for n in results.keys()]
accs = [results[str(n)]['accuracy'] for n in noise_vals]
final_losses = [results[str(n)]['final_loss'] for n in noise_vals]
grad_mags_plot = [results[str(n)]['mean_grad_magnitude'] for n in noise_vals]
grad_stds = [results[str(n)]['std_grad_magnitude'] for n in noise_vals]

axes[0].plot(range(len(noise_vals)), accs, 'o-', color='#2196F3', linewidth=2, markersize=8)
axes[0].set_xlabel('Noise Level (Depolarizing)', fontsize=12)
axes[0].set_ylabel('Test Accuracy', fontsize=12)
axes[0].set_title('Accuracy vs Noise', fontsize=13)
axes[0].set_xticks(range(len(noise_vals)))
axes[0].set_xticklabels([str(n) for n in noise_vals], rotation=45)
axes[0].grid(True, alpha=0.3)

axes[1].plot(range(len(noise_vals)), final_losses, 's-', color='#E91E63', linewidth=2, markersize=8)
axes[1].set_xlabel('Noise Level', fontsize=12)
axes[1].set_ylabel('Final Training Loss', fontsize=12)
axes[1].set_title('Loss vs Noise', fontsize=13)
axes[1].set_xticks(range(len(noise_vals)))
axes[1].set_xticklabels([str(n) for n in noise_vals], rotation=45)
axes[1].grid(True, alpha=0.3)

axes[2].errorbar(range(len(noise_vals)), grad_mags_plot, yerr=grad_stds, fmt='D-', color='#4CAF50', linewidth=2, markersize=8, capsize=5)
axes[2].set_xlabel('Noise Level', fontsize=12)
axes[2].set_ylabel('Mean |Gradient|', fontsize=12)
axes[2].set_title('Gradient Magnitude vs Noise', fontsize=13)
axes[2].set_xticks(range(len(noise_vals)))
axes[2].set_xticklabels([str(n) for n in noise_vals], rotation=45)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/noise_analysis.png', dpi=150)
plt.close()

# Loss curves comparison
fig, ax = plt.subplots(figsize=(10, 6))
cmap = plt.cm.viridis
for i, noise in enumerate(noise_vals):
    if results[str(noise)]['losses']:
        color = cmap(i / len(noise_vals))
        ax.plot(results[str(noise)]['losses'], 'o-', label=f'p={noise}', color=color, linewidth=1.5)
ax.set_xlabel('Checkpoint (every 10 epochs)', fontsize=12)
ax.set_ylabel('Training Loss', fontsize=12)
ax.set_title('Training Loss Curves Under Different Noise Levels', fontsize=14)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/noise_loss_curves.png', dpi=150)
plt.close()

save_results = {k: {kk: vv for kk, vv in v.items() if kk != 'losses'} for k, v in results.items()}
with open('results_noise.json', 'w') as f:
    json.dump(save_results, f, indent=2)

print("\nExperiment 5 complete.")
print(json.dumps(save_results, indent=2))
