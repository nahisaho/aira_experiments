"""
Module 1: Primordial Soup Hypothesis — Extended Miller-Urey Reaction Network
Stochastic simulation of prebiotic chemistry in a reducing atmosphere.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import json, os
from datetime import datetime

np.random.seed(42)

# --- Species ---
SPECIES = [
    'CH4', 'NH3', 'H2O', 'H2', 'CO', 'CO2', 'N2',
    'HCN', 'HCHO', 'NH2CHO',          # small organics
    'Glycine', 'Alanine', 'Aspartate', # amino acids
    'Adenine', 'Uracil', 'Cytosine',   # nucleobases
    'Ribose', 'Glycerol',              # sugars/lipid precursors
    'AMP', 'Peptide2', 'Peptide5',     # polymers
    'FattyAcid'
]

INIT_COUNTS = {s: 0 for s in SPECIES}
INIT_COUNTS.update({'CH4': 50000, 'NH3': 40000, 'H2O': 100000, 'H2': 60000,
                     'CO': 10000, 'CO2': 5000, 'N2': 30000})

# --- Reactions (reactants, products, rate constant) ---
REACTIONS = [
    # Spark-discharge reactions
    (['CH4', 'NH3'],        ['HCN', 'H2', 'H2'],           0.005),
    (['CO', 'H2'],          ['HCHO'],                       0.008),
    (['HCN', 'H2O'],        ['NH2CHO'],                     0.003),
    # Strecker synthesis of amino acids
    (['HCN', 'HCHO', 'NH3'],['Glycine', 'H2O'],            0.002),
    (['HCN', 'HCHO', 'NH3'],['Alanine', 'H2O'],            0.001),
    (['HCN', 'HCHO', 'NH3'],['Aspartate', 'H2O'],          0.0005),
    # Nucleobase synthesis (HCN polymerisation)
    (['HCN','HCN','HCN','HCN','HCN'], ['Adenine'],         1e-6),
    (['HCN','CO','NH3'],     ['Uracil','H2O'],              5e-5),
    (['HCN','CO','NH3'],     ['Cytosine','H2O'],            3e-5),
    # Formose reaction (sugar)
    (['HCHO','HCHO','HCHO'],['Ribose'],                     1e-5),
    # Nucleotide formation
    (['Adenine','Ribose'],   ['AMP'],                        1e-4),
    # Peptide bond
    (['Glycine','Glycine'],  ['Peptide2','H2O'],             5e-5),
    (['Peptide2','Glycine','Alanine','Aspartate'], ['Peptide5','H2O','H2O','H2O'], 1e-7),
    # Lipid precursor
    (['CO','H2','H2','H2'],  ['FattyAcid','H2O'],           1e-5),
    (['HCHO','HCHO'],        ['Glycerol'],                   1e-5),
]

def propensity(counts, rxn_reactants, k):
    """Compute propensity for a reaction given species counts."""
    a = k
    from collections import Counter
    rc = Counter(rxn_reactants)
    for sp, n in rc.items():
        c = counts.get(sp, 0)
        for i in range(n):
            a *= max(c - i, 0)
    return a

def gillespie_step(counts, reactions):
    """One step of the Gillespie SSA."""
    props = []
    for reactants, products, k in reactions:
        props.append(propensity(counts, reactants, k))
    a0 = sum(props)
    if a0 == 0:
        return counts, float('inf'), -1
    r1, r2 = np.random.random(), np.random.random()
    tau = -np.log(r1) / a0
    cumsum = 0
    j = 0
    for i, a in enumerate(props):
        cumsum += a
        if cumsum > r2 * a0:
            j = i
            break
    reactants, products, _ = reactions[j]
    from collections import Counter
    for sp, n in Counter(reactants).items():
        counts[sp] -= n
    for sp, n in Counter(products).items():
        counts[sp] = counts.get(sp, 0) + n
    return counts, tau, j

def run_simulation(max_time=5000, max_steps=500000):
    counts = dict(INIT_COUNTS)
    t = 0.0
    record_times = []
    record_counts = {s: [] for s in SPECIES}
    rxn_freq = [0] * len(REACTIONS)

    step = 0
    sample_interval = max(1, max_steps // 2000)
    while t < max_time and step < max_steps:
        counts, tau, j = gillespie_step(counts, REACTIONS)
        if tau == float('inf'):
            break
        t += tau
        if j >= 0:
            rxn_freq[j] += 1
        if step % sample_interval == 0:
            record_times.append(t)
            for s in SPECIES:
                record_counts[s].append(counts.get(s, 0))
        step += 1

    return record_times, record_counts, rxn_freq, counts

def build_reaction_network():
    G = nx.DiGraph()
    for i, (reactants, products, k) in enumerate(REACTIONS):
        rxn_node = f"R{i}"
        G.add_node(rxn_node, node_type='reaction', rate=k)
        for r in set(reactants):
            G.add_node(r, node_type='species')
            G.add_edge(r, rxn_node)
        for p in set(products):
            G.add_node(p, node_type='species')
            G.add_edge(rxn_node, p)
    return G

def plot_results(times, counts):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Extended Miller-Urey: Primordial Soup Simulation (Gillespie SSA)', fontsize=14)

    # Panel A: Small organics
    ax = axes[0, 0]
    for s in ['HCN', 'HCHO', 'NH2CHO']:
        ax.plot(times, counts[s], label=s)
    ax.set_xlabel('Time (arbitrary units)')
    ax.set_ylabel('Molecule count')
    ax.set_title('A) Small Organic Intermediates')
    ax.legend()

    # Panel B: Amino acids
    ax = axes[0, 1]
    for s in ['Glycine', 'Alanine', 'Aspartate']:
        ax.plot(times, counts[s], label=s)
    ax.set_xlabel('Time')
    ax.set_ylabel('Molecule count')
    ax.set_title('B) Amino Acid Synthesis')
    ax.legend()

    # Panel C: Nucleobases & nucleotides
    ax = axes[1, 0]
    for s in ['Adenine', 'Uracil', 'Cytosine', 'AMP']:
        ax.plot(times, counts[s], label=s)
    ax.set_xlabel('Time')
    ax.set_ylabel('Molecule count')
    ax.set_title('C) Nucleobases & Nucleotides')
    ax.legend()

    # Panel D: Polymers & lipids
    ax = axes[1, 1]
    for s in ['Peptide2', 'Peptide5', 'FattyAcid', 'Glycerol']:
        ax.plot(times, counts[s], label=s)
    ax.set_xlabel('Time')
    ax.set_ylabel('Molecule count')
    ax.set_title('D) Polymers & Lipid Precursors')
    ax.legend()

    plt.tight_layout()
    plt.savefig('figures/fig1_primordial_soup.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig1_primordial_soup.svg', bbox_inches='tight')
    plt.close()

def plot_network(G):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    species_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'species']
    rxn_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'reaction']
    pos = nx.spring_layout(G, seed=42, k=1.5)
    nx.draw_networkx_nodes(G, pos, nodelist=species_nodes, node_color='#2196F3',
                           node_size=600, ax=ax, alpha=0.9)
    nx.draw_networkx_nodes(G, pos, nodelist=rxn_nodes, node_color='#FF9800',
                           node_shape='s', node_size=300, ax=ax, alpha=0.8)
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.4, arrows=True,
                           arrowsize=15, edge_color='gray')
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)
    ax.set_title('Prebiotic Reaction Network (Blue=Species, Orange=Reactions)')
    plt.tight_layout()
    plt.savefig('figures/fig2_reaction_network.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("Running Primordial Soup simulation...")
    times, counts, rxn_freq, final = run_simulation()
    print("Plotting results...")
    plot_results(times, counts)
    G = build_reaction_network()
    plot_network(G)

    # Save results
    results = {
        'final_counts': {k: int(v) for k, v in final.items()},
        'reaction_frequencies': rxn_freq,
        'network_stats': {
            'num_species': len([n for n, d in G.nodes(data=True) if d.get('node_type') == 'species']),
            'num_reactions': len([n for n, d in G.nodes(data=True) if d.get('node_type') == 'reaction']),
            'num_edges': G.number_of_edges(),
        }
    }
    with open('results/primordial_soup_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("Module 1 complete.")
    print(f"Final amino acids: Gly={final['Glycine']}, Ala={final['Alanine']}, Asp={final['Aspartate']}")
    print(f"Nucleobases: Ade={final['Adenine']}, Ura={final['Uracil']}, Cyt={final['Cytosine']}")
    print(f"AMP={final['AMP']}, Peptide2={final['Peptide2']}, Peptide5={final['Peptide5']}")
