"""
Module 2: RNA World Hypothesis — Self-Replicator Emergence
Stochastic model for template-directed RNA replication with error-prone copying.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

np.random.seed(123)

def simulate_rna_world(
    pool_size=1000,        # initial random RNA pool
    max_length=50,         # max RNA length
    init_length=10,        # initial strand length
    replication_rate=0.02, # base replication probability per step
    error_rate=0.01,       # per-nucleotide mutation rate
    degradation_rate=0.005,# degradation probability per step
    ligation_rate=0.005,   # ligation probability
    catalytic_threshold=0.6, # min "fitness" to be a ribozyme
    n_steps=5000,
    capacity=5000
):
    """
    Each RNA strand is represented by a fitness score (proxy for secondary structure).
    Strands above catalytic_threshold can catalyse replication of neighbours.
    """
    # Initialize random pool
    fitness = np.random.beta(2, 5, size=pool_size)
    lengths = np.full(pool_size, init_length, dtype=float)

    history = {
        'step': [], 'pop_size': [], 'mean_fitness': [], 'max_fitness': [],
        'n_ribozymes': [], 'mean_length': []
    }

    for step in range(n_steps):
        n = len(fitness)
        if n == 0:
            break

        # Record
        if step % 10 == 0:
            n_ribo = int(np.sum(fitness >= catalytic_threshold))
            history['step'].append(step)
            history['pop_size'].append(n)
            history['mean_fitness'].append(float(np.mean(fitness)))
            history['max_fitness'].append(float(np.max(fitness)))
            history['n_ribozymes'].append(n_ribo)
            history['mean_length'].append(float(np.mean(lengths)))

        # --- Degradation ---
        survive = np.random.random(n) > degradation_rate
        fitness = fitness[survive]
        lengths = lengths[survive]
        n = len(fitness)
        if n == 0:
            break

        # --- Replication (template-directed, catalysed by ribozymes) ---
        ribozyme_mask = fitness >= catalytic_threshold
        n_ribozymes = np.sum(ribozyme_mask)
        if n_ribozymes > 0 and n < capacity:
            # Ribozymes can replicate themselves or nearby strands
            repl_candidates = np.where(ribozyme_mask)[0]
            n_repl = min(len(repl_candidates), capacity - n)
            chosen = np.random.choice(repl_candidates, size=n_repl, replace=True)
            offspring_fitness = fitness[chosen].copy()
            offspring_lengths = lengths[chosen].copy()

            # Introduce mutations
            mutations = np.random.normal(0, error_rate * 5, size=n_repl)
            offspring_fitness = np.clip(offspring_fitness + mutations, 0, 1)
            length_changes = np.random.choice([-1, 0, 1], size=n_repl, p=[0.1, 0.8, 0.1])
            offspring_lengths = np.clip(offspring_lengths + length_changes, 3, max_length)

            fitness = np.concatenate([fitness, offspring_fitness])
            lengths = np.concatenate([lengths, offspring_lengths])

        # --- Non-enzymatic replication (low rate) ---
        if n < capacity:
            non_enz = int(n * replication_rate * 0.1)
            if non_enz > 0:
                idx = np.random.choice(len(fitness), size=min(non_enz, capacity - len(fitness)), replace=True)
                new_f = fitness[idx] + np.random.normal(0, error_rate * 10, size=len(idx))
                new_f = np.clip(new_f, 0, 1)
                new_l = lengths[idx] + np.random.choice([-1, 0, 1], size=len(idx))
                new_l = np.clip(new_l, 3, max_length)
                fitness = np.concatenate([fitness, new_f])
                lengths = np.concatenate([lengths, new_l])

        # --- Ligation (random joining increases length & may boost fitness) ---
        if len(fitness) >= 2:
            n_lig = int(len(fitness) * ligation_rate)
            if n_lig > 0:
                pairs = np.random.choice(len(fitness), size=(n_lig, 2), replace=False)
                for i, j in pairs:
                    if lengths[i] + lengths[j] <= max_length:
                        lengths[i] += lengths[j]
                        fitness[i] = min(1.0, fitness[i] + 0.05 * np.random.random())

        # Capacity control
        if len(fitness) > capacity:
            keep = np.random.choice(len(fitness), size=capacity, replace=False)
            fitness = fitness[keep]
            lengths = lengths[keep]

    return history

def phase_diagram(error_rates, catalytic_thresholds, n_trials=3):
    """Scan parameter space: does a self-replicating population emerge?"""
    results = np.zeros((len(error_rates), len(catalytic_thresholds)))
    for i, er in enumerate(error_rates):
        for j, ct in enumerate(catalytic_thresholds):
            successes = 0
            for _ in range(n_trials):
                h = simulate_rna_world(error_rate=er, catalytic_threshold=ct,
                                        n_steps=2000, pool_size=500)
                if len(h['n_ribozymes']) > 0 and h['n_ribozymes'][-1] > 10:
                    successes += 1
            results[i, j] = successes / n_trials
    return results

def plot_rna_dynamics(history):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('RNA World: Self-Replicator Emergence Dynamics', fontsize=14)

    ax = axes[0, 0]
    ax.plot(history['step'], history['pop_size'], color='#2196F3')
    ax.set_xlabel('Step')
    ax.set_ylabel('Population size')
    ax.set_title('A) RNA Pool Size')

    ax = axes[0, 1]
    ax.plot(history['step'], history['mean_fitness'], label='Mean', color='#4CAF50')
    ax.plot(history['step'], history['max_fitness'], label='Max', color='#F44336')
    ax.set_xlabel('Step')
    ax.set_ylabel('Fitness')
    ax.set_title('B) Fitness Evolution')
    ax.legend()

    ax = axes[1, 0]
    ax.plot(history['step'], history['n_ribozymes'], color='#9C27B0')
    ax.set_xlabel('Step')
    ax.set_ylabel('Count')
    ax.set_title('C) Ribozyme Count (self-replicators)')

    ax = axes[1, 1]
    ax.plot(history['step'], history['mean_length'], color='#FF9800')
    ax.set_xlabel('Step')
    ax.set_ylabel('Nucleotides')
    ax.set_title('D) Mean RNA Length')

    plt.tight_layout()
    plt.savefig('figures/fig3_rna_world.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig3_rna_world.svg', bbox_inches='tight')
    plt.close()

def plot_phase_diagram(error_rates, cat_thresholds, results):
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(results, origin='lower', aspect='auto', cmap='viridis',
                   extent=[cat_thresholds[0], cat_thresholds[-1],
                           error_rates[0], error_rates[-1]])
    ax.set_xlabel('Catalytic Threshold')
    ax.set_ylabel('Error Rate')
    ax.set_title('RNA World Phase Diagram: Self-Replicator Emergence Probability')
    plt.colorbar(im, ax=ax, label='Emergence Probability')
    plt.tight_layout()
    plt.savefig('figures/fig4_rna_phase_diagram.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("Running RNA World simulation...")
    history = simulate_rna_world(n_steps=5000)
    plot_rna_dynamics(history)

    print("Computing phase diagram...")
    error_rates = np.linspace(0.001, 0.05, 10)
    cat_thresholds = np.linspace(0.3, 0.8, 10)
    pd_results = phase_diagram(error_rates, cat_thresholds, n_trials=3)
    plot_phase_diagram(error_rates, cat_thresholds, pd_results)

    results = {
        'final_population': history['pop_size'][-1] if history['pop_size'] else 0,
        'final_ribozymes': history['n_ribozymes'][-1] if history['n_ribozymes'] else 0,
        'final_mean_fitness': history['mean_fitness'][-1] if history['mean_fitness'] else 0,
        'final_max_fitness': history['max_fitness'][-1] if history['max_fitness'] else 0,
        'final_mean_length': history['mean_length'][-1] if history['mean_length'] else 0,
        'error_threshold_estimate': float(error_rates[np.argmax(pd_results.mean(axis=1) < 0.3)])
            if np.any(pd_results.mean(axis=1) < 0.3) else None,
    }
    with open('results/rna_world_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("Module 2 complete.")
    print(f"Final pop={results['final_population']}, ribozymes={results['final_ribozymes']}")
    print(f"Mean fitness={results['final_mean_fitness']:.4f}, Max={results['final_max_fitness']:.4f}")
