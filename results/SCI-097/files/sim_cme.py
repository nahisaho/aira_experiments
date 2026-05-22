"""
Module 4: Chemical Master Equation (CME) — Biopolymer Emergence Probability
Stochastic chemical kinetics for estimating the probability of functional polymer emergence.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import comb
import json

np.random.seed(789)

def monomer_to_polymer_gillespie(
    n_monomers=10000,
    polymerisation_rate=0.001,
    hydrolysis_rate=0.0005,
    max_length=30,
    functional_length=20,
    n_steps=200000
):
    """
    Gillespie SSA for monomer → polymer chain extension.
    Tracks length distribution over time.
    """
    # State: count of chains of each length (1 = monomer)
    chains = np.zeros(max_length + 1, dtype=int)
    chains[1] = n_monomers

    history_t = [0.0]
    history_max_len = [1]
    history_total_polymers = [0]  # chains with length >= functional_length
    history_mean_len = [1.0]

    t = 0.0
    sample_interval = max(1, n_steps // 2000)

    for step in range(n_steps):
        # Compute propensities
        # Polymerisation: chain_i + monomer -> chain_{i+1}
        poly_props = np.zeros(max_length)
        for i in range(1, max_length):
            poly_props[i] = polymerisation_rate * chains[i] * chains[1]

        # Hydrolysis: chain_i -> chain_{i-1} + monomer
        hydro_props = np.zeros(max_length + 1)
        for i in range(2, max_length + 1):
            hydro_props[i] = hydrolysis_rate * chains[i] * (i - 1)

        a0 = poly_props.sum() + hydro_props.sum()
        if a0 == 0:
            break

        # Time step
        tau = -np.log(np.random.random()) / a0

        # Choose reaction
        r = np.random.random() * a0
        cumsum = 0.0
        fired = False

        # Check polymerisation
        for i in range(1, max_length):
            cumsum += poly_props[i]
            if cumsum > r:
                # chain_i + monomer -> chain_{i+1}
                chains[i] -= 1
                chains[1] -= 1
                chains[i + 1] += 1
                fired = True
                break

        if not fired:
            for i in range(2, max_length + 1):
                cumsum += hydro_props[i]
                if cumsum > r:
                    chains[i] -= 1
                    chains[i - 1] += 1
                    chains[1] += 1
                    break

        t += tau

        if step % sample_interval == 0:
            max_len = max(1, max(j for j in range(max_length + 1) if chains[j] > 0))
            n_functional = sum(chains[functional_length:])
            total_chains = sum(chains[i] * i for i in range(max_length + 1))
            wt_mean = sum(i * chains[i] for i in range(1, max_length + 1)) / max(1, sum(chains[1:]))
            history_t.append(t)
            history_max_len.append(max_len)
            history_total_polymers.append(int(n_functional))
            history_mean_len.append(float(wt_mean))

    return {
        'times': history_t,
        'max_length': history_max_len,
        'functional_count': history_total_polymers,
        'mean_length': history_mean_len,
        'final_distribution': chains.tolist(),
    }

def analytical_probability(n_monomers, chain_length, alphabet_size=20):
    """
    Estimate probability of a specific functional sequence appearing
    given random polymerisation.
    P(specific sequence of length L) = (1/alphabet)^L
    P(at least one in N trials) = 1 - (1 - (1/a)^L)^N
    """
    lengths = np.arange(5, chain_length + 1)
    p_specific = (1.0 / alphabet_size) ** lengths
    # Number of possible chains (rough estimate)
    n_trials = n_monomers  # simplification
    p_at_least_one = 1 - (1 - p_specific) ** n_trials
    return lengths, p_specific, p_at_least_one

def plot_cme_results(res, lengths, p_specific, p_at_least_one):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('CME: Biopolymer Emergence Probability', fontsize=14)

    ax = axes[0, 0]
    ax.plot(res['times'], res['max_length'], color='#E91E63')
    ax.set_xlabel('Time')
    ax.set_ylabel('Maximum chain length')
    ax.set_title('A) Longest Polymer Chain (Gillespie)')

    ax = axes[0, 1]
    ax.plot(res['times'], res['mean_length'], color='#009688')
    ax.set_xlabel('Time')
    ax.set_ylabel('Mean chain length')
    ax.set_title('B) Mean Polymer Length')

    ax = axes[1, 0]
    dist = res['final_distribution']
    ax.bar(range(len(dist)), dist, color='#3F51B5', alpha=0.7)
    ax.set_xlabel('Chain length')
    ax.set_ylabel('Count')
    ax.set_title('C) Final Length Distribution')
    ax.set_xlim(0, 35)

    ax = axes[1, 1]
    ax.semilogy(lengths, p_specific, label='P(specific seq)', color='#F44336')
    ax.semilogy(lengths, p_at_least_one, label=f'P(≥1 in {10000} trials)', color='#4CAF50')
    ax.set_xlabel('Sequence length')
    ax.set_ylabel('Probability')
    ax.set_title('D) Functional Sequence Probability')
    ax.legend()
    ax.axhline(y=1e-10, color='gray', linestyle='--', alpha=0.5, label='Detectability limit')

    plt.tight_layout()
    plt.savefig('figures/fig7_cme_biopolymer.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig7_cme_biopolymer.svg', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("Running CME biopolymer simulation...")
    res = monomer_to_polymer_gillespie()

    print("Computing analytical probabilities...")
    lengths, p_specific, p_at_least_one = analytical_probability(10000, 30)

    plot_cme_results(res, lengths, p_specific, p_at_least_one)

    results = {
        'final_max_length': int(max(i for i in range(len(res['final_distribution']))
                                     if res['final_distribution'][i] > 0)),
        'final_functional_count': res['functional_count'][-1],
        'final_mean_length': res['mean_length'][-1],
        'p_specific_20mer': float((1/20)**20),
        'p_at_least_one_20mer_in_10k': float(1 - (1 - (1/20)**20)**10000),
        'total_monomers_remaining': int(res['final_distribution'][1]),
    }
    with open('results/cme_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("Module 4 complete.")
    print(f"Max chain length: {results['final_max_length']}")
    print(f"P(specific 20-mer) = {results['p_specific_20mer']:.2e}")
