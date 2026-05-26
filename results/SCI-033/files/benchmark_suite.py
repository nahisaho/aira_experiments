#!/usr/bin/env python3
"""
Quantum Machine Learning Benchmark Suite
=========================================
Systematic comparison of quantum and classical ML models.
Covers: expressibility, entanglement capability, quantum kernels,
data encoding strategies, barren plateaus, and noise analysis.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pennylane as qml
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from scipy.stats import entropy as kl_divergence
from scipy.linalg import expm
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

RESULTS = {}

# ============================================================
# 1. EXPRESSIBILITY & ENTANGLEMENT CAPABILITY
# ============================================================

def compute_expressibility(circuit_fn, n_qubits, n_layers, n_samples=500, n_bins=75):
    """Compute expressibility via KL divergence from Haar distribution."""
    dev = qml.device("default.qubit", wires=n_qubits)
    
    @qml.qnode(dev)
    def get_state(params):
        circuit_fn(params, n_qubits, n_layers)
        return qml.state()
    
    n_params = n_qubits * n_layers * 2
    fidelities = []
    for _ in range(n_samples):
        p1 = np.random.uniform(0, 2*np.pi, n_params)
        p2 = np.random.uniform(0, 2*np.pi, n_params)
        s1 = get_state(p1)
        s2 = get_state(p2)
        fid = np.abs(np.dot(np.conj(s1), s2))**2
        fidelities.append(fid)
    
    dim = 2**n_qubits
    haar_pdf = lambda f: (dim - 1) * (1 - f)**(dim - 2)
    
    hist, bin_edges = np.histogram(fidelities, bins=n_bins, range=(0, 1), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    haar_hist = np.array([haar_pdf(x) for x in bin_centers])
    
    hist = hist + 1e-10
    haar_hist = haar_hist + 1e-10
    hist = hist / hist.sum()
    haar_hist = haar_hist / haar_hist.sum()
    
    expr = kl_divergence(hist, haar_hist)
    return expr, fidelities

def compute_entanglement_capability(circuit_fn, n_qubits, n_layers, n_samples=200):
    """Compute average Meyer-Wallach entanglement measure."""
    dev = qml.device("default.qubit", wires=n_qubits)
    
    @qml.qnode(dev)
    def get_state(params):
        circuit_fn(params, n_qubits, n_layers)
        return qml.state()
    
    n_params = n_qubits * n_layers * 2
    ent_values = []
    
    for _ in range(n_samples):
        params = np.random.uniform(0, 2*np.pi, n_params)
        state = get_state(params)
        state = state.reshape([2]*n_qubits)
        
        total_ent = 0
        for k in range(n_qubits):
            axes_to_trace = list(range(n_qubits))
            axes_to_trace.remove(k)
            rho_k = np.tensordot(state, np.conj(state), axes=(axes_to_trace, axes_to_trace))
            purity = np.real(np.trace(rho_k @ rho_k))
            total_ent += 1 - purity
        
        mw = (2.0 / n_qubits) * total_ent
        ent_values.append(mw)
    
    return np.mean(ent_values), np.std(ent_values)

# Circuit ansatze
def circuit_1_hardware_efficient(params, n_qubits, n_layers):
    """Hardware-efficient ansatz with RY-RZ and CNOT."""
    idx = 0
    for l in range(n_layers):
        for q in range(n_qubits):
            qml.RY(params[idx], wires=q)
            idx += 1
            qml.RZ(params[idx], wires=q)
            idx += 1
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q+1])

def circuit_2_strongly_entangling(params, n_qubits, n_layers):
    """Strongly entangling layers with all-to-all CNOT."""
    idx = 0
    for l in range(n_layers):
        for q in range(n_qubits):
            qml.RX(params[idx], wires=q)
            idx += 1
            qml.RZ(params[idx], wires=q)
            idx += 1
        for q in range(n_qubits):
            qml.CNOT(wires=[q, (q+1) % n_qubits])

def circuit_3_simplified(params, n_qubits, n_layers):
    """Simplified two-design circuit with RY and CZ."""
    idx = 0
    for l in range(n_layers):
        for q in range(n_qubits):
            qml.RY(params[idx], wires=q)
            idx += 1
        for q in range(0, n_qubits - 1, 2):
            qml.CZ(wires=[q, q+1])
        for q in range(n_qubits):
            qml.RY(params[idx], wires=q)
            idx += 1
        for q in range(1, n_qubits - 1, 2):
            qml.CZ(wires=[q, q+1])

def circuit_4_iqp_inspired(params, n_qubits, n_layers):
    """IQP-inspired circuit."""
    idx = 0
    for l in range(n_layers):
        for q in range(n_qubits):
            qml.Hadamard(wires=q)
        for q in range(n_qubits):
            qml.RZ(params[idx], wires=q)
            idx += 1
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q+1])
            qml.RZ(params[idx], wires=q+1)
            idx += 1
            qml.CNOT(wires=[q, q+1])

def run_expressibility_benchmark():
    print("=" * 60)
    print("EXPERIMENT 1: Expressibility & Entanglement Capability")
    print("=" * 60)
    
    circuits = {
        'HW-Efficient': circuit_1_hardware_efficient,
        'Strongly-Ent': circuit_2_strongly_entangling,
        'Simplified': circuit_3_simplified,
        'IQP-Inspired': circuit_4_iqp_inspired,
    }
    
    n_qubits = 4
    layers_range = [1, 2, 3, 4]
    
    expr_results = {name: [] for name in circuits}
    ent_results = {name: [] for name in circuits}
    
    for name, cfn in circuits.items():
        print(f"  Circuit: {name}")
        for nl in layers_range:
            e, _ = compute_expressibility(cfn, n_qubits, nl, n_samples=300)
            ent_mean, ent_std = compute_entanglement_capability(cfn, n_qubits, nl, n_samples=150)
            expr_results[name].append(e)
            ent_results[name].append(ent_mean)
            print(f"    Layers={nl}: Expr={e:.4f}, Ent={ent_mean:.4f}±{ent_std:.4f}")
    
    # Plot expressibility
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    markers = ['o', 's', '^', 'D']
    for i, (name, vals) in enumerate(expr_results.items()):
        axes[0].plot(layers_range, vals, marker=markers[i], label=name, linewidth=2)
    axes[0].set_xlabel('Number of Layers', fontsize=12)
    axes[0].set_ylabel('Expressibility (KL Divergence)', fontsize=12)
    axes[0].set_title('Expressibility vs Circuit Depth', fontsize=14)
    axes[0].legend()
    axes[0].set_yscale('log')
    axes[0].grid(True, alpha=0.3)
    
    for i, (name, vals) in enumerate(ent_results.items()):
        axes[1].plot(layers_range, vals, marker=markers[i], label=name, linewidth=2)
    axes[1].set_xlabel('Number of Layers', fontsize=12)
    axes[1].set_ylabel('Entanglement Capability (MW)', fontsize=12)
    axes[1].set_title('Entanglement Capability vs Circuit Depth', fontsize=14)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/expressibility_entanglement.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    RESULTS['expressibility'] = expr_results
    RESULTS['entanglement'] = ent_results
    print("  → Saved figures/expressibility_entanglement.png")

# ============================================================
# 2. QUANTUM KERNEL METHODS
# ============================================================

def generate_dataset(name, n_samples=200, n_features=2):
    """Generate synthetic datasets for benchmarking."""
    if name == 'linear':
        X = np.random.randn(n_samples, n_features)
        y = (X[:, 0] + X[:, 1] > 0).astype(int)
    elif name == 'xor':
        X = np.random.randn(n_samples, n_features)
        y = ((X[:, 0] * X[:, 1]) > 0).astype(int)
    elif name == 'circle':
        r = np.random.randn(n_samples) * 0.3
        theta = np.random.uniform(0, 2*np.pi, n_samples)
        X = np.zeros((n_samples, n_features))
        labels = np.random.choice([0, 1], n_samples)
        for i in range(n_samples):
            radius = 1.0 + labels[i] * 1.5 + r[i]
            X[i, 0] = radius * np.cos(theta[i])
            X[i, 1] = radius * np.sin(theta[i])
        y = labels
    elif name == 'quantum_friendly':
        X = np.random.uniform(0, 2*np.pi, (n_samples, n_features))
        y = (np.sin(X[:, 0]) * np.cos(X[:, 1]) + np.sin(X[:, 0]*X[:, 1]) > 0).astype(int)
    elif name == 'checkerboard':
        X = np.random.uniform(-2, 2, (n_samples, n_features))
        y = ((np.floor(X[:, 0]) + np.floor(X[:, 1])) % 2 == 0).astype(int)
    else:
        raise ValueError(f"Unknown dataset: {name}")
    
    scaler = MinMaxScaler(feature_range=(0, np.pi))
    X = scaler.fit_transform(X)
    return X, y

def quantum_kernel_matrix(X1, X2, n_qubits=2, encoding='angle'):
    """Compute quantum kernel matrix."""
    dev = qml.device("default.qubit", wires=n_qubits)
    
    @qml.qnode(dev)
    def kernel_circuit(x1, x2):
        if encoding == 'angle':
            for i in range(n_qubits):
                qml.RY(x1[i % len(x1)], wires=i)
                qml.RZ(x1[(i+1) % len(x1)], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i+1])
            for i in range(n_qubits):
                qml.adjoint(qml.RZ)(x2[(i+1) % len(x2)], wires=i)
                qml.adjoint(qml.RY)(x2[i % len(x2)], wires=i)
        elif encoding == 'iqp':
            for i in range(n_qubits):
                qml.Hadamard(wires=i)
                qml.RZ(x1[i % len(x1)], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i+1])
                qml.RZ(x1[0] * x1[1], wires=i+1)
                qml.CNOT(wires=[i, i+1])
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i+1])
                qml.RZ(-x2[0] * x2[1], wires=i+1)
                qml.CNOT(wires=[i, i+1])
            for i in range(n_qubits):
                qml.adjoint(qml.RZ)(x2[i % len(x2)], wires=i)
                qml.Hadamard(wires=i)
        return qml.probs(wires=range(n_qubits))
    
    K = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            probs = kernel_circuit(X1[i], X2[j])
            K[i, j] = probs[0]
    return K

def run_kernel_benchmark():
    print("\n" + "=" * 60)
    print("EXPERIMENT 2: Quantum Kernel Methods Comparison")
    print("=" * 60)
    
    datasets = ['linear', 'xor', 'circle', 'quantum_friendly', 'checkerboard']
    results_kernel = {}
    
    for ds_name in datasets:
        print(f"  Dataset: {ds_name}")
        X, y = generate_dataset(ds_name, n_samples=120)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        
        # Classical SVM (RBF)
        svm_rbf = SVC(kernel='rbf', gamma='auto')
        svm_rbf.fit(X_train, y_train)
        acc_rbf = accuracy_score(y_test, svm_rbf.predict(X_test))
        
        # Classical SVM (Polynomial)
        svm_poly = SVC(kernel='poly', degree=3)
        svm_poly.fit(X_train, y_train)
        acc_poly = accuracy_score(y_test, svm_poly.predict(X_test))
        
        # Classical MLP
        mlp = MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=500, random_state=42)
        mlp.fit(X_train, y_train)
        acc_mlp = accuracy_score(y_test, mlp.predict(X_test))
        
        # Quantum Kernel (angle)
        K_train_q = quantum_kernel_matrix(X_train, X_train, encoding='angle')
        K_test_q = quantum_kernel_matrix(X_test, X_train, encoding='angle')
        svm_q = SVC(kernel='precomputed')
        svm_q.fit(K_train_q, y_train)
        acc_q_angle = accuracy_score(y_test, svm_q.predict(K_test_q))
        
        # Quantum Kernel (IQP)
        K_train_iqp = quantum_kernel_matrix(X_train, X_train, encoding='iqp')
        K_test_iqp = quantum_kernel_matrix(X_test, X_train, encoding='iqp')
        svm_iqp = SVC(kernel='precomputed')
        svm_iqp.fit(K_train_iqp, y_train)
        acc_q_iqp = accuracy_score(y_test, svm_iqp.predict(K_test_iqp))
        
        results_kernel[ds_name] = {
            'SVM-RBF': acc_rbf,
            'SVM-Poly': acc_poly,
            'MLP': acc_mlp,
            'QK-Angle': acc_q_angle,
            'QK-IQP': acc_q_iqp,
        }
        print(f"    RBF={acc_rbf:.3f} Poly={acc_poly:.3f} MLP={acc_mlp:.3f} "
              f"QK-Angle={acc_q_angle:.3f} QK-IQP={acc_q_iqp:.3f}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    x_pos = np.arange(len(datasets))
    width = 0.15
    methods = ['SVM-RBF', 'SVM-Poly', 'MLP', 'QK-Angle', 'QK-IQP']
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#E91E63', '#9C27B0']
    
    for i, method in enumerate(methods):
        vals = [results_kernel[ds][method] for ds in datasets]
        ax.bar(x_pos + i*width, vals, width, label=method, color=colors[i], alpha=0.85)
    
    ax.set_xlabel('Dataset', fontsize=12)
    ax.set_ylabel('Accuracy', fontsize=12)
    ax.set_title('Quantum vs Classical Kernel Methods', fontsize=14)
    ax.set_xticks(x_pos + width*2)
    ax.set_xticklabels(datasets, rotation=15)
    ax.legend(loc='lower right')
    ax.set_ylim(0.3, 1.05)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig('figures/kernel_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    RESULTS['kernel'] = results_kernel
    print("  → Saved figures/kernel_comparison.png")

# ============================================================
# 3. DATA ENCODING STRATEGIES
# ============================================================

def run_encoding_benchmark():
    print("\n" + "=" * 60)
    print("EXPERIMENT 3: Data Encoding Strategies")
    print("=" * 60)
    
    n_qubits = 4
    n_layers = 2
    n_samples_train = 80
    n_samples_test = 40
    
    X, y = generate_dataset('quantum_friendly', n_samples=n_samples_train + n_samples_test, n_features=4)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=n_samples_test, random_state=42)
    
    encoding_results = {}
    
    for enc_name in ['angle', 'amplitude', 'iqp']:
        print(f"  Encoding: {enc_name}")
        dev = qml.device("default.qubit", wires=n_qubits)
        
        @qml.qnode(dev, interface="autograd")
        def circuit(params, x):
            # Encoding layer
            if enc_name == 'angle':
                for i in range(n_qubits):
                    qml.RY(x[i % len(x)], wires=i)
            elif enc_name == 'amplitude':
                x_norm = x / (np.linalg.norm(x) + 1e-10)
                padded = np.zeros(2**n_qubits)
                padded[:len(x_norm)] = x_norm
                padded = padded / (np.linalg.norm(padded) + 1e-10)
                qml.StatePrep(padded, wires=range(n_qubits))
            elif enc_name == 'iqp':
                for i in range(n_qubits):
                    qml.Hadamard(wires=i)
                    qml.RZ(x[i % len(x)], wires=i)
                for i in range(n_qubits - 1):
                    qml.CNOT(wires=[i, i+1])
                    qml.RZ(x[i % len(x)] * x[(i+1) % len(x)], wires=i+1)
                    qml.CNOT(wires=[i, i+1])
            
            # Variational layers
            idx = 0
            for l in range(n_layers):
                for q in range(n_qubits):
                    qml.RY(params[idx], wires=q)
                    idx += 1
                    qml.RZ(params[idx], wires=q)
                    idx += 1
                for q in range(n_qubits - 1):
                    qml.CNOT(wires=[q, q+1])
            
            return qml.expval(qml.PauliZ(0))
        
        n_params = n_qubits * n_layers * 2
        params = np.random.uniform(0, 2*np.pi, n_params)
        
        opt = qml.AdamOptimizer(stepsize=0.15)
        costs = []
        
        batch_size = min(20, len(X_train))
        
        def cost_fn(params):
            indices = np.random.choice(len(X_train), batch_size, replace=False)
            total = 0
            for i in indices:
                pred = circuit(params, X_train[i])
                label = 2 * y_train[i] - 1
                total += (pred - label) ** 2
            return total / batch_size
        
        for epoch in range(40):
            params = opt.step(cost_fn, params)
            c = cost_fn(params)
            costs.append(float(c))
            if epoch % 10 == 0:
                print(f"    Epoch {epoch}: Cost={c:.4f}")
        
        # Evaluate
        preds = []
        for i in range(len(X_test)):
            p = circuit(params, X_test[i])
            preds.append(1 if p > 0 else 0)
        acc = accuracy_score(y_test, preds)
        
        encoding_results[enc_name] = {
            'accuracy': acc,
            'final_cost': costs[-1],
            'costs': costs,
        }
        print(f"    Final accuracy: {acc:.3f}")
    
    # Plot training curves and accuracy
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors_enc = {'angle': '#E91E63', 'amplitude': '#2196F3', 'iqp': '#4CAF50'}
    
    for enc, data in encoding_results.items():
        axes[0].plot(data['costs'], label=enc, color=colors_enc[enc], linewidth=2)
    axes[0].set_xlabel('Epoch', fontsize=12)
    axes[0].set_ylabel('Cost', fontsize=12)
    axes[0].set_title('Training Convergence by Encoding', fontsize=14)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    accs = [encoding_results[e]['accuracy'] for e in ['angle', 'amplitude', 'iqp']]
    bars = axes[1].bar(['Angle', 'Amplitude', 'IQP'], accs, 
                        color=[colors_enc[e] for e in ['angle', 'amplitude', 'iqp']], alpha=0.85)
    axes[1].set_ylabel('Test Accuracy', fontsize=12)
    axes[1].set_title('Classification Accuracy by Encoding', fontsize=14)
    axes[1].set_ylim(0, 1.1)
    for bar, acc in zip(bars, accs):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                    f'{acc:.3f}', ha='center', fontsize=11)
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/encoding_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    RESULTS['encoding'] = {k: {kk: vv for kk, vv in v.items() if kk != 'costs'} 
                           for k, v in encoding_results.items()}
    RESULTS['encoding_costs'] = {k: v['costs'] for k, v in encoding_results.items()}
    print("  → Saved figures/encoding_comparison.png")

# ============================================================
# 4. QUANTUM ADVANTAGE DATASET CHARACTERIZATION
# ============================================================

def run_dataset_characterization():
    print("\n" + "=" * 60)
    print("EXPERIMENT 4: Dataset Characterization for Quantum Advantage")
    print("=" * 60)
    
    datasets = ['linear', 'xor', 'circle', 'quantum_friendly', 'checkerboard']
    n_qubits = 2
    
    characteristics = {}
    
    for ds_name in datasets:
        X, y = generate_dataset(ds_name, n_samples=150)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        
        # Classical baseline
        svm = SVC(kernel='rbf', gamma='auto')
        svm.fit(X_train, y_train)
        acc_classical = accuracy_score(y_test, svm.predict(X_test))
        
        # Quantum kernel
        K_train = quantum_kernel_matrix(X_train, X_train, encoding='iqp')
        K_test = quantum_kernel_matrix(X_test, X_train, encoding='iqp')
        svm_q = SVC(kernel='precomputed')
        svm_q.fit(K_train, y_train)
        acc_quantum = accuracy_score(y_test, svm_q.predict(K_test))
        
        # Dataset complexity measures
        from sklearn.neighbors import KNeighborsClassifier
        knn = KNeighborsClassifier(n_neighbors=1)
        knn.fit(X_train, y_train)
        acc_1nn = accuracy_score(y_test, knn.predict(X_test))
        
        # Kernel alignment (quantum advantage proxy)
        q_advantage = acc_quantum - acc_classical
        
        characteristics[ds_name] = {
            'classical_acc': acc_classical,
            'quantum_acc': acc_quantum,
            'advantage': q_advantage,
            '1nn_acc': acc_1nn,
            'nonlinearity': 1 - acc_1nn,
        }
        print(f"  {ds_name}: Classical={acc_classical:.3f} Quantum={acc_quantum:.3f} "
              f"Δ={q_advantage:+.3f}")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    ds_names = list(characteristics.keys())
    classical_accs = [characteristics[d]['classical_acc'] for d in ds_names]
    quantum_accs = [characteristics[d]['quantum_acc'] for d in ds_names]
    advantages = [characteristics[d]['advantage'] for d in ds_names]
    
    x = np.arange(len(ds_names))
    axes[0].bar(x - 0.2, classical_accs, 0.4, label='Classical (RBF)', color='#2196F3', alpha=0.85)
    axes[0].bar(x + 0.2, quantum_accs, 0.4, label='Quantum (IQP)', color='#E91E63', alpha=0.85)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ds_names, rotation=15)
    axes[0].set_ylabel('Accuracy', fontsize=12)
    axes[0].set_title('Classical vs Quantum Accuracy', fontsize=14)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, axis='y')
    
    colors_adv = ['#4CAF50' if a > 0 else '#F44336' for a in advantages]
    axes[1].bar(ds_names, advantages, color=colors_adv, alpha=0.85)
    axes[1].set_ylabel('Quantum Advantage (Δ Accuracy)', fontsize=12)
    axes[1].set_title('Quantum Advantage by Dataset', fontsize=14)
    axes[1].axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    axes[1].set_xticklabels(ds_names, rotation=15)
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/dataset_characterization.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    RESULTS['dataset_char'] = characteristics
    print("  → Saved figures/dataset_characterization.png")

# ============================================================
# 5. BARREN PLATEAU ANALYSIS
# ============================================================

def run_barren_plateau_analysis():
    print("\n" + "=" * 60)
    print("EXPERIMENT 5: Barren Plateau & Trainability Analysis")
    print("=" * 60)
    
    qubit_range = [2, 3, 4, 5, 6]
    n_samples = 100
    n_layers = 3
    
    grad_results = {'global': {}, 'local': {}}
    
    for cost_type in ['global', 'local']:
        print(f"  Cost type: {cost_type}")
        grad_variances = []
        grad_means = []
        
        for n_q in qubit_range:
            dev = qml.device("default.qubit", wires=n_q)
            n_params = n_q * n_layers * 2
            
            if cost_type == 'global':
                @qml.qnode(dev, diff_method='parameter-shift')
                def cost_circuit(params):
                    idx = 0
                    for l in range(n_layers):
                        for q in range(n_q):
                            qml.RY(params[idx], wires=q)
                            idx += 1
                            qml.RZ(params[idx], wires=q)
                            idx += 1
                        for q in range(n_q - 1):
                            qml.CNOT(wires=[q, q+1])
                    # Global cost: sum of PauliZ on all qubits
                    return qml.expval(
                        sum(qml.PauliZ(i) for i in range(n_q))
                    )
            else:
                @qml.qnode(dev, diff_method='parameter-shift')
                def cost_circuit(params):
                    idx = 0
                    for l in range(n_layers):
                        for q in range(n_q):
                            qml.RY(params[idx], wires=q)
                            idx += 1
                            qml.RZ(params[idx], wires=q)
                            idx += 1
                        for q in range(n_q - 1):
                            qml.CNOT(wires=[q, q+1])
                    return qml.expval(qml.PauliZ(0))
            
            grads_first_param = []
            for _ in range(n_samples):
                params = np.random.uniform(0, 2*np.pi, n_params)
                # Use parameter-shift manually for first parameter
                shift = np.zeros(n_params)
                shift[0] = np.pi / 2
                p_plus = params + shift
                p_minus = params - shift
                grad_val = (cost_circuit(p_plus) - cost_circuit(p_minus)) / 2.0
                grads_first_param.append(float(grad_val))
            
            var_grad = float(np.var(grads_first_param))
            mean_grad = float(np.mean(np.abs(grads_first_param)))
            grad_variances.append(var_grad)
            grad_means.append(mean_grad)
            print(f"    n_qubits={n_q}: Var(∂C/∂θ₁)={var_grad:.6f}, <|∂C/∂θ₁|>={mean_grad:.6f}")
        
        grad_results[cost_type] = {
            'variances': grad_variances,
            'means': grad_means,
        }
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].semilogy(qubit_range, grad_results['global']['variances'], 'o-', 
                      label='Global Cost', color='#E91E63', linewidth=2, markersize=8)
    axes[0].semilogy(qubit_range, grad_results['local']['variances'], 's-', 
                      label='Local Cost', color='#2196F3', linewidth=2, markersize=8)
    axes[0].set_xlabel('Number of Qubits', fontsize=12)
    axes[0].set_ylabel('Var(∂C/∂θ₁)', fontsize=12)
    axes[0].set_title('Gradient Variance vs System Size', fontsize=14)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].semilogy(qubit_range, grad_results['global']['means'], 'o-',
                      label='Global Cost', color='#E91E63', linewidth=2, markersize=8)
    axes[1].semilogy(qubit_range, grad_results['local']['means'], 's-',
                      label='Local Cost', color='#2196F3', linewidth=2, markersize=8)
    axes[1].set_xlabel('Number of Qubits', fontsize=12)
    axes[1].set_ylabel('Mean |∂C/∂θ₁|', fontsize=12)
    axes[1].set_title('Mean Gradient Magnitude vs System Size', fontsize=14)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/barren_plateau.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    RESULTS['barren_plateau'] = {
        'qubit_range': qubit_range,
        'global_variances': grad_results['global']['variances'],
        'local_variances': grad_results['local']['variances'],
        'global_means': grad_results['global']['means'],
        'local_means': grad_results['local']['means'],
    }
    print("  → Saved figures/barren_plateau.png")

# ============================================================
# 6. NOISE ANALYSIS (Simulated IBM Quantum Noise)
# ============================================================

def run_noise_analysis():
    print("\n" + "=" * 60)
    print("EXPERIMENT 6: Noise Impact Analysis (Simulated)")
    print("=" * 60)
    
    n_qubits = 4
    n_layers = 2
    noise_levels = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    
    X, y = generate_dataset('quantum_friendly', n_samples=100, n_features=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    noise_results = {}
    
    for noise_p in noise_levels:
        print(f"  Noise level: {noise_p}")
        
        if noise_p == 0:
            dev = qml.device("default.qubit", wires=n_qubits)
        else:
            dev = qml.device("default.mixed", wires=n_qubits)
        
        @qml.qnode(dev)
        def noisy_circuit(params, x):
            for i in range(min(len(x), n_qubits)):
                qml.RY(x[i], wires=i)
            
            idx = 0
            for l in range(n_layers):
                for q in range(n_qubits):
                    qml.RY(params[idx], wires=q)
                    idx += 1
                    qml.RZ(params[idx], wires=q)
                    idx += 1
                    if noise_p > 0:
                        qml.DepolarizingChannel(noise_p, wires=q)
                for q in range(n_qubits - 1):
                    qml.CNOT(wires=[q, q+1])
                    if noise_p > 0:
                        qml.DepolarizingChannel(noise_p * 2, wires=q)
                        qml.DepolarizingChannel(noise_p * 2, wires=q+1)
            
            return qml.expval(qml.PauliZ(0))
        
        n_params = n_qubits * n_layers * 2
        best_acc = 0
        
        for trial in range(3):
            params = np.random.uniform(0, 2*np.pi, n_params)
            opt = qml.GradientDescentOptimizer(stepsize=0.1)
            
            def cost_fn(params):
                total = 0
                for i in range(len(X_train)):
                    pred = noisy_circuit(params, X_train[i])
                    label = 2 * y_train[i] - 1
                    total += (pred - label) ** 2
                return total / len(X_train)
            
            for epoch in range(20):
                params = opt.step(cost_fn, params)
            
            preds = []
            for i in range(len(X_test)):
                p = noisy_circuit(params, X_test[i])
                preds.append(1 if p > 0 else 0)
            acc = accuracy_score(y_test, preds)
            best_acc = max(best_acc, acc)
        
        # Fidelity degradation
        dev_ideal = qml.device("default.qubit", wires=n_qubits)
        
        @qml.qnode(dev_ideal)
        def ideal_state(params):
            idx = 0
            for l in range(n_layers):
                for q in range(n_qubits):
                    qml.RY(params[idx], wires=q)
                    idx += 1
                    qml.RZ(params[idx], wires=q)
                    idx += 1
                for q in range(n_qubits - 1):
                    qml.CNOT(wires=[q, q+1])
            return qml.state()
        
        fidelities = []
        for _ in range(20):
            rp = np.random.uniform(0, 2*np.pi, n_params)
            ideal = ideal_state(rp)
            
            if noise_p > 0:
                dev_noisy = qml.device("default.mixed", wires=n_qubits)
                @qml.qnode(dev_noisy)
                def noisy_state(params):
                    idx = 0
                    for l in range(n_layers):
                        for q in range(n_qubits):
                            qml.RY(params[idx], wires=q)
                            idx += 1
                            qml.RZ(params[idx], wires=q)
                            idx += 1
                            qml.DepolarizingChannel(noise_p, wires=q)
                        for q in range(n_qubits - 1):
                            qml.CNOT(wires=[q, q+1])
                            qml.DepolarizingChannel(noise_p*2, wires=q)
                    return qml.density_matrix(wires=range(n_qubits))
                
                rho = noisy_state(rp)
                fid = np.real(np.conj(ideal) @ rho @ ideal)
                fidelities.append(float(fid))
            else:
                fidelities.append(1.0)
        
        avg_fid = np.mean(fidelities)
        noise_results[noise_p] = {
            'accuracy': best_acc,
            'fidelity': avg_fid,
        }
        print(f"    Accuracy={best_acc:.3f}, Avg Fidelity={avg_fid:.4f}")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    nl = list(noise_results.keys())
    accs = [noise_results[n]['accuracy'] for n in nl]
    fids = [noise_results[n]['fidelity'] for n in nl]
    
    axes[0].plot(nl, accs, 'o-', color='#E91E63', linewidth=2, markersize=8)
    axes[0].set_xlabel('Depolarizing Noise Rate', fontsize=12)
    axes[0].set_ylabel('Classification Accuracy', fontsize=12)
    axes[0].set_title('Accuracy Degradation under Noise', fontsize=14)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xscale('symlog', linthresh=0.001)
    
    axes[1].plot(nl, fids, 's-', color='#2196F3', linewidth=2, markersize=8)
    axes[1].set_xlabel('Depolarizing Noise Rate', fontsize=12)
    axes[1].set_ylabel('Average State Fidelity', fontsize=12)
    axes[1].set_title('Fidelity Degradation under Noise', fontsize=14)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xscale('symlog', linthresh=0.001)
    
    plt.tight_layout()
    plt.savefig('figures/noise_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    RESULTS['noise'] = {str(k): v for k, v in noise_results.items()}
    print("  → Saved figures/noise_analysis.png")

# ============================================================
# COMPREHENSIVE SUMMARY FIGURE
# ============================================================

def create_summary_figure():
    """Create a comprehensive summary figure."""
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # Panel 1: Expressibility
    ax1 = fig.add_subplot(gs[0, 0])
    if 'expressibility' in RESULTS:
        for name, vals in RESULTS['expressibility'].items():
            ax1.plot([1,2,3,4], vals, 'o-', label=name, linewidth=1.5)
        ax1.set_yscale('log')
        ax1.set_xlabel('Layers')
        ax1.set_ylabel('Expressibility (KL)')
        ax1.set_title('(a) Expressibility')
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)
    
    # Panel 2: Kernel comparison
    ax2 = fig.add_subplot(gs[0, 1])
    if 'kernel' in RESULTS:
        datasets = list(RESULTS['kernel'].keys())
        methods = ['SVM-RBF', 'QK-Angle', 'QK-IQP']
        x = np.arange(len(datasets))
        w = 0.25
        for i, m in enumerate(methods):
            vals = [RESULTS['kernel'][d][m] for d in datasets]
            ax2.bar(x + i*w, vals, w, label=m, alpha=0.85)
        ax2.set_xticks(x + w)
        ax2.set_xticklabels(datasets, rotation=30, fontsize=7)
        ax2.set_ylabel('Accuracy')
        ax2.set_title('(b) Kernel Methods')
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3, axis='y')
    
    # Panel 3: Encoding
    ax3 = fig.add_subplot(gs[0, 2])
    if 'encoding' in RESULTS:
        encs = list(RESULTS['encoding'].keys())
        accs = [RESULTS['encoding'][e]['accuracy'] for e in encs]
        ax3.bar(encs, accs, color=['#E91E63', '#2196F3', '#4CAF50'], alpha=0.85)
        ax3.set_ylabel('Accuracy')
        ax3.set_title('(c) Encoding Strategies')
        ax3.grid(True, alpha=0.3, axis='y')
    
    # Panel 4: Dataset characterization
    ax4 = fig.add_subplot(gs[1, 0])
    if 'dataset_char' in RESULTS:
        ds = list(RESULTS['dataset_char'].keys())
        advs = [RESULTS['dataset_char'][d]['advantage'] for d in ds]
        colors = ['#4CAF50' if a > 0 else '#F44336' for a in advs]
        ax4.bar(ds, advs, color=colors, alpha=0.85)
        ax4.axhline(y=0, color='k', linewidth=0.5)
        ax4.set_ylabel('Δ Accuracy')
        ax4.set_title('(d) Quantum Advantage')
        ax4.set_xticklabels(ds, rotation=30, fontsize=7)
        ax4.grid(True, alpha=0.3, axis='y')
    
    # Panel 5: Barren plateau
    ax5 = fig.add_subplot(gs[1, 1])
    if 'barren_plateau' in RESULTS:
        bp = RESULTS['barren_plateau']
        ax5.semilogy(bp['qubit_range'], bp['global_variances'], 'o-', label='Global', color='#E91E63')
        ax5.semilogy(bp['qubit_range'], bp['local_variances'], 's-', label='Local', color='#2196F3')
        ax5.set_xlabel('Qubits')
        ax5.set_ylabel('Var(∂C/∂θ)')
        ax5.set_title('(e) Barren Plateaus')
        ax5.legend(fontsize=8)
        ax5.grid(True, alpha=0.3)
    
    # Panel 6: Noise
    ax6 = fig.add_subplot(gs[1, 2])
    if 'noise' in RESULTS:
        nl = sorted([float(k) for k in RESULTS['noise'].keys()])
        accs = [RESULTS['noise'][str(n)]['accuracy'] for n in nl]
        fids = [RESULTS['noise'][str(n)]['fidelity'] for n in nl]
        ax6.plot(nl, accs, 'o-', label='Accuracy', color='#E91E63')
        ax6.plot(nl, fids, 's-', label='Fidelity', color='#2196F3')
        ax6.set_xlabel('Noise Rate')
        ax6.set_title('(f) Noise Impact')
        ax6.set_xscale('symlog', linthresh=0.001)
        ax6.legend(fontsize=8)
        ax6.grid(True, alpha=0.3)
    
    plt.savefig('figures/summary.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\n  → Saved figures/summary.png")

# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Quantum ML Benchmark Suite")
    print("=" * 60)
    
    run_expressibility_benchmark()
    run_kernel_benchmark()
    run_encoding_benchmark()
    run_dataset_characterization()
    run_barren_plateau_analysis()
    run_noise_analysis()
    create_summary_figure()
    
    # Save results
    def convert_for_json(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, (np.int64, np.int32)):
            return int(obj)
        return obj
    
    def deep_convert(d):
        if isinstance(d, dict):
            return {k: deep_convert(v) for k, v in d.items()}
        if isinstance(d, list):
            return [deep_convert(i) for i in d]
        return convert_for_json(d)
    
    with open('results.json', 'w') as f:
        json.dump(deep_convert(RESULTS), f, indent=2)
    
    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETE")
    print("Results saved to results.json")
    print("Figures saved to figures/")
    print("=" * 60)
