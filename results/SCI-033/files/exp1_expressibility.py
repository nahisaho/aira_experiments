"""
Experiment 1: Expressibility and Entanglement Capability of Parameterized Quantum Circuits
Quantifies expressibility via KL divergence from Haar distribution and entanglement via Meyer-Wallach measure.
"""
import numpy as np
import pennylane as qml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import entropy
from itertools import combinations
import json

np.random.seed(42)

N_QUBITS = 4
N_SAMPLES = 500
N_BINS = 75

# Define circuit ansatze
def circuit_1(params, wires):
    """Hardware-efficient: RY-CNOT linear"""
    for i in wires:
        qml.RY(params[i], wires=i)
    for i in range(len(wires)-1):
        qml.CNOT(wires=[wires[i], wires[i+1]])

def circuit_2(params, wires):
    """Strongly entangling layers (ring CNOT)"""
    idx = 0
    for i in wires:
        qml.RY(params[idx], wires=i)
        qml.RZ(params[idx+1], wires=i)
        idx += 2
    for i in range(len(wires)):
        qml.CNOT(wires=[wires[i], wires[(i+1)%len(wires)]])

def circuit_3(params, wires):
    """IQP-style: H-ZZ-H"""
    for i in wires:
        qml.Hadamard(wires=i)
    idx = 0
    for i in range(len(wires)):
        for j in range(i+1, len(wires)):
            qml.IsingZZ(params[idx], wires=[wires[i], wires[j]])
            idx += 1
    for i in wires:
        qml.Hadamard(wires=i)

def circuit_4(params, wires):
    """Deep alternating: RX-RY-CNOT x2 layers"""
    n = len(wires)
    idx = 0
    for _ in range(2):
        for i in wires:
            qml.RX(params[idx], wires=i)
            qml.RY(params[idx+1], wires=i)
            idx += 2
        for i in range(n-1):
            qml.CNOT(wires=[wires[i], wires[i+1]])

def circuit_5(params, wires):
    """All-to-all CZ with RY"""
    idx = 0
    for i in wires:
        qml.RY(params[idx], wires=i)
        idx += 1
    for i in range(len(wires)):
        for j in range(i+1, len(wires)):
            qml.CZ(wires=[wires[i], wires[j]])
    for i in wires:
        qml.RY(params[idx], wires=i)
        idx += 1

circuits = {
    'C1: RY-CNOT linear': (circuit_1, N_QUBITS),
    'C2: RY-RZ ring CNOT': (circuit_2, 2*N_QUBITS),
    'C3: IQP (H-ZZ-H)': (circuit_3, N_QUBITS*(N_QUBITS-1)//2),
    'C4: Deep RX-RY-CNOT': (circuit_4, 4*N_QUBITS),
    'C5: All-to-all CZ': (circuit_5, 2*N_QUBITS),
}

def compute_fidelities(circuit_fn, n_params, n_qubits, n_samples):
    """Compute pairwise fidelities between random parameter samples."""
    dev = qml.device('default.qubit', wires=n_qubits)
    
    @qml.qnode(dev)
    def get_state(params):
        circuit_fn(params, wires=range(n_qubits))
        return qml.state()
    
    fidelities = []
    for _ in range(n_samples):
        p1 = np.random.uniform(0, 2*np.pi, n_params)
        p2 = np.random.uniform(0, 2*np.pi, n_params)
        s1 = get_state(p1)
        s2 = get_state(p2)
        fid = np.abs(np.dot(np.conj(s1), s2))**2
        fidelities.append(fid)
    return np.array(fidelities)

def haar_fidelity_pdf(f, n_qubits):
    """Haar-random fidelity distribution for n qubits."""
    d = 2**n_qubits
    return (d - 1) * (1 - f)**(d - 2)

def compute_expressibility(fidelities, n_qubits, n_bins):
    """KL divergence between sampled fidelity and Haar distribution."""
    hist, bin_edges = np.histogram(fidelities, bins=n_bins, range=(0, 1), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    haar = haar_fidelity_pdf(bin_centers, n_qubits)
    haar = haar / haar.sum()
    hist = hist / hist.sum()
    hist = np.clip(hist, 1e-10, None)
    haar = np.clip(haar, 1e-10, None)
    return entropy(hist, haar)

def compute_meyer_wallach(circuit_fn, n_params, n_qubits, n_samples=200):
    """Compute Meyer-Wallach entanglement measure."""
    dev = qml.device('default.qubit', wires=n_qubits)
    
    @qml.qnode(dev)
    def get_state(params):
        circuit_fn(params, wires=range(n_qubits))
        return qml.state()
    
    q_values = []
    for _ in range(n_samples):
        params = np.random.uniform(0, 2*np.pi, n_params)
        state = get_state(params)
        state = state.reshape([2]*n_qubits)
        
        entanglement = 0
        for k in range(n_qubits):
            # Partial trace over qubit k
            axes_to_trace = list(range(n_qubits))
            axes_to_trace.remove(k)
            # Compute reduced density matrix for qubit k
            rho_k = np.zeros((2, 2), dtype=complex)
            for idx_vals in np.ndindex(*([2]*(n_qubits-1))):
                for a in range(2):
                    for b in range(2):
                        idx_a = list(idx_vals)
                        idx_a.insert(k, a)
                        idx_b = list(idx_vals)
                        idx_b.insert(k, b)
                        rho_k[a, b] += state[tuple(idx_a)] * np.conj(state[tuple(idx_b)])
            purity = np.real(np.trace(rho_k @ rho_k))
            entanglement += 1 - purity
        
        q_values.append(2 * entanglement / n_qubits)
    
    return np.mean(q_values), np.std(q_values)

print("Computing expressibility and entanglement for 5 circuit ansatze...")
results = {}
for name, (circ, n_params) in circuits.items():
    print(f"  {name}...")
    fids = compute_fidelities(circ, n_params, N_QUBITS, N_SAMPLES)
    expr = compute_expressibility(fids, N_QUBITS, N_BINS)
    mw_mean, mw_std = compute_meyer_wallach(circ, n_params, N_QUBITS, 200)
    results[name] = {
        'expressibility_kl': float(expr),
        'entanglement_mw_mean': float(mw_mean),
        'entanglement_mw_std': float(mw_std),
        'fidelities': fids.tolist()
    }
    print(f"    Expr(KL): {expr:.4f}, Ent(MW): {mw_mean:.4f}±{mw_std:.4f}")

# Plot 1: Expressibility vs Entanglement scatter
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
names = list(results.keys())
expr_vals = [results[n]['expressibility_kl'] for n in names]
ent_vals = [results[n]['entanglement_mw_mean'] for n in names]
ent_stds = [results[n]['entanglement_mw_std'] for n in names]

colors = plt.cm.Set2(np.linspace(0, 1, len(names)))
for i, n in enumerate(names):
    ax.errorbar(expr_vals[i], ent_vals[i], yerr=ent_stds[i],
                fmt='o', markersize=12, color=colors[i], label=n, capsize=5)

ax.set_xlabel('Expressibility (KL Divergence, lower = more expressive)', fontsize=12)
ax.set_ylabel('Meyer-Wallach Entanglement', fontsize=12)
ax.set_title('Expressibility vs Entanglement Capability', fontsize=14)
ax.legend(fontsize=9, loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/expressibility_entanglement.png', dpi=150)
plt.close()

# Plot 2: Fidelity distributions
fig, axes = plt.subplots(1, len(names), figsize=(20, 4), sharey=True)
bins = np.linspace(0, 1, N_BINS)
bin_centers = (bins[:-1] + bins[1:]) / 2
haar = haar_fidelity_pdf(bin_centers, N_QUBITS)
haar = haar / (haar.sum() * (bins[1]-bins[0]))

for i, n in enumerate(names):
    axes[i].hist(results[n]['fidelities'], bins=bins, density=True, alpha=0.7, color=colors[i], label='Sampled')
    axes[i].plot(bin_centers, haar_fidelity_pdf(bin_centers, N_QUBITS), 'k--', label='Haar', linewidth=2)
    axes[i].set_title(n.split(':')[0], fontsize=10)
    axes[i].set_xlabel('Fidelity')
    if i == 0:
        axes[i].set_ylabel('Density')
        axes[i].legend(fontsize=8)

plt.suptitle('Fidelity Distributions vs Haar Random', fontsize=14)
plt.tight_layout()
plt.savefig('figures/fidelity_distributions.png', dpi=150)
plt.close()

# Save numerical results
save_results = {k: {kk: vv for kk, vv in v.items() if kk != 'fidelities'} for k, v in results.items()}
with open('results_expressibility.json', 'w') as f:
    json.dump(save_results, f, indent=2)

print("\nExperiment 1 complete. Figures saved.")
print(json.dumps(save_results, indent=2))
