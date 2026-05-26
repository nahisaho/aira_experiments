#!/usr/bin/env python3
"""
Quantum Annealing Performance Evaluation Framework
===================================================
Evaluates QUBO formulation strategies, annealing schedules,
minor embedding proxies, and classical solver comparisons
using OpenJij (simulated quantum annealing).

Case study: Capacitated Vehicle Routing Problem (CVRP)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openjij as oj
from scipy.optimize import minimize
import networkx as nx
import time
import json
import os
import warnings
warnings.filterwarnings('ignore')

FIGURES_DIR = "figures"
RESULTS = {}

np.random.seed(42)

# ============================================================
# 1. VRP Instance Generation
# ============================================================

def generate_vrp_instance(n_customers, n_vehicles, grid_size=100):
    """Generate a random CVRP instance."""
    depot = np.array([grid_size / 2, grid_size / 2])
    customers = np.random.rand(n_customers, 2) * grid_size
    locations = np.vstack([depot, customers])
    n = len(locations)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dist_matrix[i, j] = np.linalg.norm(locations[i] - locations[j])
    demands = np.random.randint(1, 10, size=n_customers)
    capacity = int(np.sum(demands) / n_vehicles * 1.5)
    return {
        'locations': locations,
        'dist_matrix': dist_matrix,
        'demands': demands,
        'capacity': capacity,
        'n_customers': n_customers,
        'n_vehicles': n_vehicles,
        'depot': depot
    }


# ============================================================
# 2. QUBO Formulation Strategies for VRP
# ============================================================

def vrp_to_qubo_standard(instance, penalty_weight=None):
    """Standard QUBO formulation for VRP with quadratic penalties."""
    n_c = instance['n_customers']
    n_v = instance['n_vehicles']
    n_steps = n_c  # max steps per vehicle
    D = instance['dist_matrix']

    if penalty_weight is None:
        penalty_weight = np.max(D) * 2.0

    n_vars = n_c * n_v * n_steps
    Q = {}

    def idx(c, v, s):
        return c * n_v * n_steps + v * n_steps + s

    # Objective: minimize travel distance (simplified)
    for v in range(n_v):
        for s in range(n_steps - 1):
            for c1 in range(n_c):
                for c2 in range(n_c):
                    i, j = idx(c1, v, s), idx(c2, v, s + 1)
                    key = (min(i, j), max(i, j)) if i != j else (i, i)
                    Q[key] = Q.get(key, 0) + D[c1 + 1, c2 + 1]

    # Constraint: each customer visited exactly once
    for c in range(n_c):
        indices = [idx(c, v, s) for v in range(n_v) for s in range(n_steps)]
        for i in indices:
            Q[(i, i)] = Q.get((i, i), 0) - penalty_weight
        for a in range(len(indices)):
            for b in range(a + 1, len(indices)):
                key = (min(indices[a], indices[b]), max(indices[a], indices[b]))
                Q[key] = Q.get(key, 0) + 2 * penalty_weight

    # Constraint: at most one customer per time step per vehicle
    for v in range(n_v):
        for s in range(n_steps):
            indices = [idx(c, v, s) for c in range(n_c)]
            for a in range(len(indices)):
                for b in range(a + 1, len(indices)):
                    key = (min(indices[a], indices[b]), max(indices[a], indices[b]))
                    Q[key] = Q.get(key, 0) + penalty_weight

    return Q, n_vars


def vrp_to_qubo_compact(instance, penalty_weight=None):
    """Compact QUBO: reduced variable count with aggregated constraints."""
    n_c = instance['n_customers']
    n_v = instance['n_vehicles']
    D = instance['dist_matrix']

    if penalty_weight is None:
        penalty_weight = np.max(D) * 1.5

    # Simplified: only assignment variables x_{c,v}
    n_vars = n_c * n_v
    Q = {}

    def idx(c, v):
        return c * n_v + v

    # Objective: approximate cost as assignment-based
    for v in range(n_v):
        for c1 in range(n_c):
            for c2 in range(c1 + 1, n_c):
                i, j = idx(c1, v), idx(c2, v)
                Q[(min(i, j), max(i, j))] = Q.get((min(i, j), max(i, j)), 0) + D[c1 + 1, c2 + 1] * 0.5

    # Constraint: each customer assigned to exactly one vehicle
    for c in range(n_c):
        indices = [idx(c, v) for v in range(n_v)]
        for i in indices:
            Q[(i, i)] = Q.get((i, i), 0) - penalty_weight
        for a in range(len(indices)):
            for b in range(a + 1, len(indices)):
                key = (min(indices[a], indices[b]), max(indices[a], indices[b]))
                Q[key] = Q.get(key, 0) + 2 * penalty_weight

    return Q, n_vars


def vrp_to_qubo_adaptive(instance):
    """Adaptive penalty QUBO: penalty weights tuned per constraint type."""
    n_c = instance['n_customers']
    n_v = instance['n_vehicles']
    D = instance['dist_matrix']

    obj_scale = np.mean(D[1:, 1:])
    penalty_assignment = obj_scale * 3.0
    penalty_capacity = obj_scale * 2.0

    n_vars = n_c * n_v
    Q = {}

    def idx(c, v):
        return c * n_v + v

    for v in range(n_v):
        for c1 in range(n_c):
            for c2 in range(c1 + 1, n_c):
                i, j = idx(c1, v), idx(c2, v)
                Q[(min(i, j), max(i, j))] = Q.get((min(i, j), max(i, j)), 0) + D[c1 + 1, c2 + 1] * 0.5

    for c in range(n_c):
        indices = [idx(c, v) for v in range(n_v)]
        for i in indices:
            Q[(i, i)] = Q.get((i, i), 0) - penalty_assignment
        for a in range(len(indices)):
            for b in range(a + 1, len(indices)):
                key = (min(indices[a], indices[b]), max(indices[a], indices[b]))
                Q[key] = Q.get(key, 0) + 2 * penalty_assignment

    return Q, n_vars


# ============================================================
# 3. Solver Wrappers
# ============================================================

def solve_sa(Q, n_vars, num_reads=100, beta_range=None, num_sweeps=1000):
    """Simulated Annealing via OpenJij."""
    sampler = oj.SASampler()
    if beta_range:
        sampler.beta_min = beta_range[0]
        sampler.beta_max = beta_range[1]
    sampler.num_sweeps = num_sweeps

    h = {}
    J = {}
    for (i, j), v in Q.items():
        if i == j:
            h[i] = h.get(i, 0) + v
        else:
            J[(i, j)] = J.get((i, j), 0) + v

    start = time.time()
    response = sampler.sample_qubo(Q, num_reads=num_reads)
    elapsed = time.time() - start

    energies = [s.energy for s in response.record]
    best_energy = min(energies)
    best_sample = response.first.sample

    return {
        'best_energy': best_energy,
        'best_sample': best_sample,
        'energies': energies,
        'time': elapsed,
        'method': 'SA'
    }


def solve_sqa(Q, n_vars, num_reads=100, trotter=4, beta=5.0, num_sweeps=1000):
    """Simulated Quantum Annealing via OpenJij."""
    sampler = oj.SQASampler()
    sampler.trotter = trotter
    sampler.beta = beta
    sampler.num_sweeps = num_sweeps

    start = time.time()
    response = sampler.sample_qubo(Q, num_reads=num_reads)
    elapsed = time.time() - start

    energies = [s.energy for s in response.record]
    best_energy = min(energies)
    best_sample = response.first.sample

    return {
        'best_energy': best_energy,
        'best_sample': best_sample,
        'energies': energies,
        'time': elapsed,
        'method': 'SQA'
    }


def solve_qaoa_inspired(Q, n_vars, p_layers=3, num_samples=50):
    """QAOA-inspired classical optimizer (variational approach)."""
    start = time.time()

    def qaoa_cost(params):
        gamma = params[:p_layers]
        beta_p = params[p_layers:]
        total_cost = 0
        for _ in range(num_samples):
            state = np.random.choice([0, 1], size=n_vars)
            for layer in range(p_layers):
                # Phase separation (probabilistic flip based on cost gradient)
                for (i, j), v in Q.items():
                    if i == j:
                        if np.random.rand() < abs(np.sin(gamma[layer] * v)):
                            state[i] = 1 - state[i]
                    else:
                        if state[i] == state[j] == 1:
                            if np.random.rand() < abs(np.sin(gamma[layer] * v * 0.1)):
                                idx_flip = np.random.choice([i, j])
                                state[idx_flip] = 1 - state[idx_flip]
                # Mixing
                for k in range(n_vars):
                    if np.random.rand() < abs(np.sin(beta_p[layer])):
                        state[k] = 1 - state[k]
            energy = sum(v * state[i] * (state[j] if i != j else 1)
                         for (i, j), v in Q.items())
            total_cost += energy
        return total_cost / num_samples

    init_params = np.random.rand(2 * p_layers) * np.pi
    result = minimize(qaoa_cost, init_params, method='COBYLA',
                      options={'maxiter': 100})

    # Final sampling
    best_energy = float('inf')
    best_sample = None
    gamma = result.x[:p_layers]
    beta_p = result.x[p_layers:]

    for _ in range(num_samples * 2):
        state = np.random.choice([0, 1], size=n_vars)
        for layer in range(p_layers):
            for (i, j), v in Q.items():
                if i == j:
                    if np.random.rand() < abs(np.sin(gamma[layer] * v)):
                        state[i] = 1 - state[i]
                else:
                    if state[i] == state[j] == 1:
                        if np.random.rand() < abs(np.sin(gamma[layer] * v * 0.1)):
                            idx_flip = np.random.choice([i, j])
                            state[idx_flip] = 1 - state[idx_flip]
            for k in range(n_vars):
                if np.random.rand() < abs(np.sin(beta_p[layer])):
                    state[k] = 1 - state[k]

        energy = sum(v * state[i] * (state[j] if i != j else 1)
                     for (i, j), v in Q.items())
        if energy < best_energy:
            best_energy = energy
            best_sample = dict(enumerate(state))

    elapsed = time.time() - start
    return {
        'best_energy': best_energy,
        'best_sample': best_sample,
        'energies': [best_energy],
        'time': elapsed,
        'method': 'QAOA-inspired'
    }


# ============================================================
# 4. Annealing Schedule Experiments
# ============================================================

def annealing_schedule_experiment(Q, n_vars):
    """Compare different annealing schedules."""
    schedules = {
        'fast (100 sweeps)': {'num_sweeps': 100},
        'standard (1000 sweeps)': {'num_sweeps': 1000},
        'slow (5000 sweeps)': {'num_sweeps': 5000},
        'cold-start (high beta)': {'num_sweeps': 1000, 'beta_range': (1.0, 100.0)},
        'warm-start (low beta)': {'num_sweeps': 1000, 'beta_range': (0.01, 10.0)},
    }

    results = {}
    for name, params in schedules.items():
        beta_range = params.pop('beta_range', None)
        r = solve_sa(Q, n_vars, num_reads=50, beta_range=beta_range, **params)
        r['schedule'] = name
        results[name] = r
        if 'beta_range' not in params and beta_range:
            params['beta_range'] = beta_range

    return results


def reverse_annealing_experiment(Q, n_vars):
    """Simulate reverse annealing: start from a solution, add noise, re-anneal."""
    # Forward anneal to get initial solution
    initial = solve_sa(Q, n_vars, num_reads=20, num_sweeps=500)
    initial_state = initial['best_sample']
    initial_energy = initial['best_energy']

    results = {'forward_only': initial}

    # Reverse annealing simulation: perturb + re-anneal
    noise_levels = [0.05, 0.1, 0.2, 0.3, 0.5]
    for noise in noise_levels:
        best_of_reverse = float('inf')
        energies = []
        start = time.time()

        for trial in range(30):
            # Perturb initial state
            state = dict(initial_state)
            for k in state:
                if np.random.rand() < noise:
                    state[k] = 1 - state[k]

            # Re-anneal from perturbed state using SA with warm start
            sampler = oj.SASampler()
            sampler.num_sweeps = 500
            response = sampler.sample_qubo(Q, num_reads=1, initial_state=state)
            e = response.first.energy
            energies.append(e)
            if e < best_of_reverse:
                best_of_reverse = e

        elapsed = time.time() - start
        results[f'reverse_noise_{noise}'] = {
            'best_energy': best_of_reverse,
            'energies': energies,
            'time': elapsed,
            'method': f'Reverse(noise={noise})',
            'improvement': initial_energy - best_of_reverse
        }

    return results


# ============================================================
# 5. Problem Scaling Experiment
# ============================================================

def scaling_experiment():
    """Measure performance scaling with problem size."""
    sizes = [4, 6, 8, 10, 12]
    results = {'sizes': sizes, 'sa': [], 'sqa': [], 'qaoa': []}

    for n in sizes:
        print(f"  Scaling: n_customers={n}")
        inst = generate_vrp_instance(n, 2)
        Q, nv = vrp_to_qubo_compact(inst)

        sa_r = solve_sa(Q, nv, num_reads=30, num_sweeps=1000)
        results['sa'].append({'energy': sa_r['best_energy'], 'time': sa_r['time']})

        sqa_r = solve_sqa(Q, nv, num_reads=30, num_sweeps=1000)
        results['sqa'].append({'energy': sqa_r['best_energy'], 'time': sqa_r['time']})

        if n <= 8:
            qaoa_r = solve_qaoa_inspired(Q, nv, p_layers=2, num_samples=20)
            results['qaoa'].append({'energy': qaoa_r['best_energy'], 'time': qaoa_r['time']})
        else:
            results['qaoa'].append({'energy': None, 'time': None})

    return results


# ============================================================
# 6. Minor Embedding Proxy Analysis
# ============================================================

def embedding_analysis(Q, n_vars):
    """Analyze QUBO graph structure as proxy for minor embedding difficulty."""
    G = nx.Graph()
    for (i, j), v in Q.items():
        if i != j and v != 0:
            G.add_edge(i, j, weight=abs(v))

    if len(G.nodes()) == 0:
        return {'density': 0, 'avg_degree': 0, 'max_degree': 0,
                'clustering': 0, 'components': 0, 'bandwidth': 0}

    degrees = [d for _, d in G.degree()]

    # Chimera-like embedding estimation
    chimera_unit = 8  # K_{4,4} unit cell
    n_logical = len(G.nodes())
    max_deg = max(degrees) if degrees else 0
    estimated_chain_length = max(1, max_deg / chimera_unit)
    estimated_physical_qubits = n_logical * estimated_chain_length

    return {
        'n_logical_qubits': n_logical,
        'n_edges': len(G.edges()),
        'density': nx.density(G),
        'avg_degree': np.mean(degrees) if degrees else 0,
        'max_degree': max_deg,
        'clustering_coeff': nx.average_clustering(G),
        'n_components': nx.number_connected_components(G),
        'est_chain_length': estimated_chain_length,
        'est_physical_qubits': estimated_physical_qubits
    }


# ============================================================
# 7. Visualization
# ============================================================

def plot_vrp_instance(instance, filename="vrp_instance.png"):
    """Plot VRP instance layout."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    locs = instance['locations']
    ax.plot(locs[0, 0], locs[0, 1], 'rs', markersize=15, label='Depot')
    ax.scatter(locs[1:, 0], locs[1:, 1], c='blue', s=80, zorder=5, label='Customers')
    for i in range(1, len(locs)):
        ax.annotate(f'C{i}', (locs[i, 0] + 1, locs[i, 1] + 1), fontsize=8)
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_title(f'VRP Instance: {instance["n_customers"]} customers, {instance["n_vehicles"]} vehicles')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_qubo_comparison(results_dict, filename="qubo_comparison.png"):
    """Compare QUBO formulation strategies."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    names = list(results_dict.keys())
    energies = [results_dict[n]['best_energy'] for n in names]
    times = [results_dict[n]['time'] for n in names]
    n_vars_list = [results_dict[n]['n_vars'] for n in names]

    axes[0].bar(names, energies, color=['#2196F3', '#4CAF50', '#FF9800'])
    axes[0].set_ylabel('Best Energy')
    axes[0].set_title('Solution Quality by QUBO Formulation')
    axes[0].tick_params(axis='x', rotation=15)

    axes[1].bar(names, times, color=['#2196F3', '#4CAF50', '#FF9800'])
    axes[1].set_ylabel('Time (s)')
    axes[1].set_title('Computation Time')
    axes[1].tick_params(axis='x', rotation=15)

    axes[2].bar(names, n_vars_list, color=['#2196F3', '#4CAF50', '#FF9800'])
    axes[2].set_ylabel('Number of Variables')
    axes[2].set_title('QUBO Size')
    axes[2].tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_solver_comparison(results_dict, filename="solver_comparison.png"):
    """Compare solvers."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    methods = list(results_dict.keys())
    energies = [results_dict[m]['best_energy'] for m in methods]
    times = [results_dict[m]['time'] for m in methods]

    colors = ['#E91E63', '#9C27B0', '#3F51B5'][:len(methods)]
    axes[0].bar(methods, energies, color=colors)
    axes[0].set_ylabel('Best Energy')
    axes[0].set_title('Solution Quality: SA vs SQA vs QAOA')
    axes[0].tick_params(axis='x', rotation=15)

    axes[1].bar(methods, times, color=colors)
    axes[1].set_ylabel('Time (s)')
    axes[1].set_title('Computation Time')
    axes[1].tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_annealing_schedules(results, filename="annealing_schedules.png"):
    """Plot annealing schedule comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    names = list(results.keys())
    best_energies = [results[n]['best_energy'] for n in names]
    times_list = [results[n]['time'] for n in names]

    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(names)))
    axes[0].barh(names, best_energies, color=colors)
    axes[0].set_xlabel('Best Energy')
    axes[0].set_title('Schedule vs Solution Quality')

    axes[1].barh(names, times_list, color=colors)
    axes[1].set_xlabel('Time (s)')
    axes[1].set_title('Schedule vs Computation Time')

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_reverse_annealing(results, filename="reverse_annealing.png"):
    """Plot reverse annealing results."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    forward_energy = results['forward_only']['best_energy']
    noise_levels = []
    improvements = []
    best_energies = []

    for key, val in results.items():
        if key.startswith('reverse_noise'):
            noise = float(key.split('_')[-1])
            noise_levels.append(noise)
            improvements.append(val['improvement'])
            best_energies.append(val['best_energy'])

    axes[0].plot(noise_levels, best_energies, 'bo-', linewidth=2, markersize=8, label='Reverse Anneal')
    axes[0].axhline(y=forward_energy, color='r', linestyle='--', linewidth=2, label='Forward Only')
    axes[0].set_xlabel('Perturbation Noise Level')
    axes[0].set_ylabel('Best Energy')
    axes[0].set_title('Reverse Annealing: Energy vs Noise')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].bar(noise_levels, improvements, width=0.03, color='green', alpha=0.7)
    axes[1].set_xlabel('Perturbation Noise Level')
    axes[1].set_ylabel('Energy Improvement')
    axes[1].set_title('Reverse Annealing: Improvement over Forward')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_scaling(results, filename="scaling_analysis.png"):
    """Plot scaling analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sizes = results['sizes']
    sa_times = [r['time'] for r in results['sa']]
    sqa_times = [r['time'] for r in results['sqa']]
    qaoa_times = [r['time'] if r['time'] else None for r in results['qaoa']]

    sa_energies = [r['energy'] for r in results['sa']]
    sqa_energies = [r['energy'] for r in results['sqa']]

    axes[0].plot(sizes, sa_times, 'ro-', linewidth=2, label='SA', markersize=8)
    axes[0].plot(sizes, sqa_times, 'bs-', linewidth=2, label='SQA', markersize=8)
    valid_qaoa = [(s, t) for s, t in zip(sizes, qaoa_times) if t is not None]
    if valid_qaoa:
        axes[0].plot([x[0] for x in valid_qaoa], [x[1] for x in valid_qaoa],
                     'g^-', linewidth=2, label='QAOA', markersize=8)
    axes[0].set_xlabel('Number of Customers')
    axes[0].set_ylabel('Time (s)')
    axes[0].set_title('Solver Time Scaling')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')

    axes[1].plot(sizes, sa_energies, 'ro-', linewidth=2, label='SA', markersize=8)
    axes[1].plot(sizes, sqa_energies, 'bs-', linewidth=2, label='SQA', markersize=8)
    axes[1].set_xlabel('Number of Customers')
    axes[1].set_ylabel('Best Energy')
    axes[1].set_title('Solution Quality Scaling')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_embedding_analysis(embed_results, filename="embedding_analysis.png"):
    """Plot embedding analysis results."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = list(embed_results.keys())
    logical = [embed_results[n]['n_logical_qubits'] for n in names]
    physical = [embed_results[n]['est_physical_qubits'] for n in names]
    chain = [embed_results[n]['est_chain_length'] for n in names]

    x = np.arange(len(names))
    w = 0.35
    axes[0].bar(x - w/2, logical, w, label='Logical Qubits', color='#2196F3')
    axes[0].bar(x + w/2, physical, w, label='Est. Physical Qubits', color='#FF5722')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(names, rotation=15)
    axes[0].set_ylabel('Qubit Count')
    axes[0].set_title('Embedding Overhead by Formulation')
    axes[0].legend()

    axes[1].bar(names, chain, color='#4CAF50')
    axes[1].set_ylabel('Estimated Chain Length')
    axes[1].set_title('Chain Length by Formulation')
    axes[1].tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_energy_distribution(results_dict, filename="energy_distribution.png"):
    """Plot energy distribution across reads."""
    fig, ax = plt.subplots(figsize=(10, 6))
    for method, result in results_dict.items():
        if len(result['energies']) > 1:
            ax.hist(result['energies'], bins=20, alpha=0.5, label=method)
    ax.set_xlabel('Energy')
    ax.set_ylabel('Frequency')
    ax.set_title('Energy Distribution Across Reads')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# Main Experiment Pipeline
# ============================================================

def main():
    print("=" * 60)
    print("Quantum Annealing Performance Evaluation Framework")
    print("=" * 60)

    # --- Generate VRP instance ---
    print("\n[1/7] Generating VRP instance...")
    instance = generate_vrp_instance(n_customers=8, n_vehicles=2)
    plot_vrp_instance(instance)
    RESULTS['instance'] = {
        'n_customers': instance['n_customers'],
        'n_vehicles': instance['n_vehicles'],
        'capacity': instance['capacity']
    }
    print(f"  Created: {instance['n_customers']} customers, {instance['n_vehicles']} vehicles")

    # --- QUBO Formulation Comparison ---
    print("\n[2/7] Comparing QUBO formulations...")
    qubo_results = {}

    Q_std, nv_std = vrp_to_qubo_standard(instance)
    r_std = solve_sa(Q_std, nv_std, num_reads=50)
    qubo_results['Standard'] = {**r_std, 'n_vars': nv_std}
    print(f"  Standard: {nv_std} vars, energy={r_std['best_energy']:.2f}")

    Q_comp, nv_comp = vrp_to_qubo_compact(instance)
    r_comp = solve_sa(Q_comp, nv_comp, num_reads=50)
    qubo_results['Compact'] = {**r_comp, 'n_vars': nv_comp}
    print(f"  Compact:  {nv_comp} vars, energy={r_comp['best_energy']:.2f}")

    Q_adap, nv_adap = vrp_to_qubo_adaptive(instance)
    r_adap = solve_sa(Q_adap, nv_adap, num_reads=50)
    qubo_results['Adaptive'] = {**r_adap, 'n_vars': nv_adap}
    print(f"  Adaptive: {nv_adap} vars, energy={r_adap['best_energy']:.2f}")

    plot_qubo_comparison(qubo_results)
    RESULTS['qubo_comparison'] = {
        k: {'energy': v['best_energy'], 'time': v['time'], 'n_vars': v['n_vars']}
        for k, v in qubo_results.items()
    }

    # --- Embedding Analysis ---
    print("\n[3/7] Analyzing embedding characteristics...")
    embed_results = {}
    embed_results['Standard'] = embedding_analysis(Q_std, nv_std)
    embed_results['Compact'] = embedding_analysis(Q_comp, nv_comp)
    embed_results['Adaptive'] = embedding_analysis(Q_adap, nv_adap)
    plot_embedding_analysis(embed_results)
    RESULTS['embedding'] = embed_results
    for name, er in embed_results.items():
        print(f"  {name}: logical={er['n_logical_qubits']}, est_physical={er['est_physical_qubits']:.0f}, chain={er['est_chain_length']:.1f}")

    # --- Solver Comparison ---
    print("\n[4/7] Comparing solvers (SA vs SQA vs QAOA)...")
    Q_test, nv_test = vrp_to_qubo_compact(instance)
    solver_results = {}

    r_sa = solve_sa(Q_test, nv_test, num_reads=100)
    solver_results['SA'] = r_sa
    print(f"  SA:   energy={r_sa['best_energy']:.2f}, time={r_sa['time']:.3f}s")

    r_sqa = solve_sqa(Q_test, nv_test, num_reads=100)
    solver_results['SQA'] = r_sqa
    print(f"  SQA:  energy={r_sqa['best_energy']:.2f}, time={r_sqa['time']:.3f}s")

    r_qaoa = solve_qaoa_inspired(Q_test, nv_test, p_layers=3, num_samples=30)
    solver_results['QAOA'] = r_qaoa
    print(f"  QAOA: energy={r_qaoa['best_energy']:.2f}, time={r_qaoa['time']:.3f}s")

    plot_solver_comparison(solver_results)
    plot_energy_distribution(solver_results)
    RESULTS['solver_comparison'] = {
        k: {'energy': v['best_energy'], 'time': v['time']}
        for k, v in solver_results.items()
    }

    # --- Annealing Schedule ---
    print("\n[5/7] Testing annealing schedules...")
    schedule_results = annealing_schedule_experiment(Q_test, nv_test)
    plot_annealing_schedules(schedule_results)
    RESULTS['schedules'] = {
        k: {'energy': v['best_energy'], 'time': v['time']}
        for k, v in schedule_results.items()
    }
    for name, sr in schedule_results.items():
        print(f"  {name}: energy={sr['best_energy']:.2f}, time={sr['time']:.3f}s")

    # --- Reverse Annealing ---
    print("\n[6/7] Reverse annealing experiment...")
    reverse_results = reverse_annealing_experiment(Q_test, nv_test)
    plot_reverse_annealing(reverse_results)
    RESULTS['reverse_annealing'] = {}
    for key, val in reverse_results.items():
        RESULTS['reverse_annealing'][key] = {
            'energy': val['best_energy'],
            'time': val.get('time', 0)
        }
        if key != 'forward_only':
            print(f"  {key}: energy={val['best_energy']:.2f}, improvement={val.get('improvement', 0):.2f}")

    # --- Scaling ---
    print("\n[7/7] Problem scaling experiment...")
    scaling_results = scaling_experiment()
    plot_scaling(scaling_results)
    RESULTS['scaling'] = scaling_results

    # --- Save results ---
    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            r = convert(obj)
            if r is not obj:
                return r
            return super().default(obj)

    with open('experiment_results.json', 'w') as f:
        json.dump(RESULTS, f, indent=2, cls=NumpyEncoder)

    print("\n" + "=" * 60)
    print("All experiments completed successfully!")
    print(f"Results saved to experiment_results.json")
    print(f"Figures saved to {FIGURES_DIR}/")
    print("=" * 60)


if __name__ == '__main__':
    main()
