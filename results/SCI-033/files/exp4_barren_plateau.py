"""
Experiment 4: Barren Plateau Analysis
Studies gradient variance decay with increasing qubits and circuit depth.
"""
import numpy as np
import pennylane as qml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

np.random.seed(42)

N_SAMPLES = 100

def hardware_efficient_circuit(params, n_qubits, n_layers):
    """Hardware-efficient ansatz."""
    idx = 0
    for l in range(n_layers):
        for i in range(n_qubits):
            qml.RY(params[idx], wires=i)
            qml.RZ(params[idx+1], wires=i)
            idx += 2
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i+1])

def compute_gradient_variance(n_qubits, n_layers, cost_type='global', n_samples=N_SAMPLES):
    """Compute variance of gradients for given circuit configuration."""
    n_params = 2 * n_qubits * n_layers
    dev = qml.device('default.qubit', wires=n_qubits)
    
    if cost_type == 'global':
        obs = qml.PauliZ(0)
        for i in range(1, n_qubits):
            obs = obs @ qml.PauliZ(i)
    else:
        obs = qml.PauliZ(0)
    
    @qml.qnode(dev, diff_method='parameter-shift')
    def circuit(params):
        hardware_efficient_circuit(params, n_qubits, n_layers)
        return qml.expval(obs)
    
    gradients = []
    for _ in range(n_samples):
        params = qml.numpy.array(np.random.uniform(0, 2*np.pi, n_params), requires_grad=True)
        grad = qml.grad(circuit)(params)
        gradients.append(float(grad[0]))
    
    return float(np.var(gradients)), float(np.mean(np.abs(gradients)))

# Study 1: Gradient variance vs number of qubits (fixed depth)
print("Study 1: Gradient variance vs qubits...")
qubit_range = [2, 3, 4, 5, 6, 7, 8]
n_layers_fixed = 2

var_global = []
var_local = []
mean_global = []
mean_local = []

for n_q in qubit_range:
    print(f"  n_qubits={n_q}...")
    vg, mg = compute_gradient_variance(n_q, n_layers_fixed, 'global', N_SAMPLES)
    vl, ml = compute_gradient_variance(n_q, n_layers_fixed, 'local', N_SAMPLES)
    var_global.append(vg)
    var_local.append(vl)
    mean_global.append(mg)
    mean_local.append(ml)
    print(f"    Global var: {vg:.6f}, Local var: {vl:.6f}")

# Study 2: Gradient variance vs depth (fixed qubits)
print("\nStudy 2: Gradient variance vs depth...")
n_qubits_fixed = 4
layer_range = [1, 2, 3, 4, 5, 6, 8, 10]

var_depth = []
mean_depth = []

for n_l in layer_range:
    print(f"  n_layers={n_l}...")
    v, m = compute_gradient_variance(n_qubits_fixed, n_l, 'global', N_SAMPLES)
    var_depth.append(v)
    mean_depth.append(m)
    print(f"    Variance: {v:.6f}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: Gradient variance vs qubits
axes[0].semilogy(qubit_range, var_global, 'o-', label='Global cost', color='#E91E63', linewidth=2, markersize=8)
axes[0].semilogy(qubit_range, var_local, 's-', label='Local cost', color='#2196F3', linewidth=2, markersize=8)
axes[0].set_xlabel('Number of Qubits', fontsize=12)
axes[0].set_ylabel('Gradient Variance (log scale)', fontsize=12)
axes[0].set_title('Barren Plateau: Qubits', fontsize=13)
axes[0].legend(fontsize=11)
axes[0].grid(True, alpha=0.3)

# Panel 2: Gradient variance vs depth
axes[1].semilogy(layer_range, var_depth, 'D-', color='#4CAF50', linewidth=2, markersize=8)
axes[1].set_xlabel('Number of Layers', fontsize=12)
axes[1].set_ylabel('Gradient Variance (log scale)', fontsize=12)
axes[1].set_title(f'Barren Plateau: Depth (n={n_qubits_fixed})', fontsize=13)
axes[1].grid(True, alpha=0.3)

# Panel 3: Mean absolute gradient vs qubits
axes[2].semilogy(qubit_range, mean_global, 'o-', label='Global cost', color='#E91E63', linewidth=2, markersize=8)
axes[2].semilogy(qubit_range, mean_local, 's-', label='Local cost', color='#2196F3', linewidth=2, markersize=8)
axes[2].set_xlabel('Number of Qubits', fontsize=12)
axes[2].set_ylabel('Mean |Gradient| (log scale)', fontsize=12)
axes[2].set_title('Trainability: Mean Gradient', fontsize=13)
axes[2].legend(fontsize=11)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/barren_plateau.png', dpi=150)
plt.close()

results = {
    'qubit_study': {
        'qubit_range': qubit_range,
        'global_variance': var_global,
        'local_variance': var_local,
        'global_mean_abs': mean_global,
        'local_mean_abs': mean_local,
    },
    'depth_study': {
        'n_qubits': n_qubits_fixed,
        'layer_range': layer_range,
        'variance': var_depth,
        'mean_abs': mean_depth,
    }
}

with open('results_barren.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nExperiment 4 complete.")
print(json.dumps(results, indent=2))
