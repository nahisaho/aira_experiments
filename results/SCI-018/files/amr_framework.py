#!/usr/bin/env python3
"""
Integrated Computational Framework for Predicting Antimicrobial Resistance Evolution
=====================================================================================
Combines population genetics simulation with epidemiological modeling to predict
AMR evolution trajectories, fitness landscapes, HGT networks, and optimal
antibiotic treatment strategies.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import integrate, optimize, stats
from scipy.spatial.distance import hamming
import networkx as nx
from itertools import product, combinations
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# ==============================================================================
# Module 1: Resistance Gene Detection Pipeline Simulation
# ==============================================================================

class ARGDetectionPipeline:
    """Simulates whole-genome sequence analysis for ARG detection."""

    GENE_FAMILIES = {
        'beta_lactamase': {'blaTEM': 0.35, 'blaSHV': 0.20, 'blaCTX-M': 0.30, 'blaKPC': 0.10, 'blaNDM': 0.05},
        'aminoglycoside': {'aac(6)': 0.40, 'aph(3)': 0.30, 'ant(2)': 0.20, 'armA': 0.10},
        'fluoroquinolone': {'qnrA': 0.25, 'qnrB': 0.30, 'qnrS': 0.25, 'aac(6)-Ib-cr': 0.20},
        'tetracycline': {'tetA': 0.30, 'tetB': 0.25, 'tetM': 0.25, 'tetX': 0.20},
        'colistin': {'mcr-1': 0.50, 'mcr-2': 0.20, 'mcr-3': 0.15, 'mcr-4': 0.10, 'mcr-5': 0.05},
    }

    ANTIBIOTICS = ['beta_lactam', 'aminoglycoside', 'fluoroquinolone', 'tetracycline', 'colistin']

    def __init__(self, n_genomes=200):
        self.n_genomes = n_genomes
        self.genomes = []
        self.arg_profiles = []

    def simulate_genomes(self):
        """Generate synthetic bacterial genomes with ARG profiles."""
        for i in range(self.n_genomes):
            genome = {
                'id': f'genome_{i:04d}',
                'species': np.random.choice(
                    ['E.coli', 'K.pneumoniae', 'A.baumannii', 'P.aeruginosa', 'S.aureus'],
                    p=[0.30, 0.25, 0.15, 0.15, 0.15]
                ),
                'genome_size_mb': np.random.normal(5.0, 0.8),
                'gc_content': np.random.normal(0.50, 0.05),
                'n_plasmids': np.random.poisson(2),
            }

            # Simulate ARG detection
            args_found = {}
            for family, genes in self.GENE_FAMILIES.items():
                family_present = np.random.random() < 0.4  # 40% chance per family
                if family_present:
                    for gene, freq in genes.items():
                        if np.random.random() < freq:
                            args_found[gene] = {
                                'identity': np.random.uniform(85, 100),
                                'coverage': np.random.uniform(80, 100),
                                'location': np.random.choice(['chromosome', 'plasmid'],
                                                              p=[0.3, 0.7]),
                                'family': family
                            }

            genome['args'] = args_found
            self.genomes.append(genome)

        return self.genomes

    def compute_resistance_profiles(self):
        """Compute resistance profiles (binary phenotype predictions)."""
        profiles = np.zeros((self.n_genomes, len(self.ANTIBIOTICS)))
        family_to_idx = {
            'beta_lactamase': 0, 'aminoglycoside': 1,
            'fluoroquinolone': 2, 'tetracycline': 3, 'colistin': 4
        }

        for i, g in enumerate(self.genomes):
            for gene_name, info in g['args'].items():
                fam = info['family']
                if fam in family_to_idx:
                    if info['identity'] > 90 and info['coverage'] > 85:
                        profiles[i, family_to_idx[fam]] = 1

        self.arg_profiles = profiles
        return profiles

    def compute_statistics(self):
        """Compute detection statistics."""
        n_args_per_genome = [len(g['args']) for g in self.genomes]
        species_counts = {}
        for g in self.genomes:
            sp = g['species']
            species_counts[sp] = species_counts.get(sp, 0) + 1

        mdr_count = np.sum(np.sum(self.arg_profiles, axis=1) >= 3)
        xdr_count = np.sum(np.sum(self.arg_profiles, axis=1) >= 4)

        return {
            'total_genomes': self.n_genomes,
            'mean_args': np.mean(n_args_per_genome),
            'max_args': np.max(n_args_per_genome),
            'species_distribution': species_counts,
            'mdr_count': int(mdr_count),
            'xdr_count': int(xdr_count),
            'resistance_rates': {
                ab: float(np.mean(self.arg_profiles[:, i]))
                for i, ab in enumerate(self.ANTIBIOTICS)
            }
        }


# ==============================================================================
# Module 2: Fitness Landscape Construction
# ==============================================================================

class FitnessLandscape:
    """Construct and analyze fitness landscapes for resistance mutations."""

    def __init__(self, n_loci=5, drug_concentrations=None):
        self.n_loci = n_loci
        self.n_genotypes = 2 ** n_loci
        self.drug_concs = drug_concentrations or [0, 0.5, 1.0, 2.0, 4.0, 8.0]
        self.genotypes = list(product([0, 1], repeat=n_loci))
        self.fitness_values = {}
        self.locus_names = [f'mut_{chr(65+i)}' for i in range(n_loci)]

    def construct_landscape(self):
        """Build NK-like fitness landscape with epistasis."""
        for conc in self.drug_concs:
            fitness = np.zeros(self.n_genotypes)
            for i, gt in enumerate(self.genotypes):
                n_mutations = sum(gt)
                # Base fitness cost of resistance mutations
                cost = 0.05 * n_mutations
                # Benefit under drug pressure
                benefit = 0
                if conc > 0:
                    for j, allele in enumerate(gt):
                        if allele == 1:
                            benefit += 0.15 * np.log2(1 + conc) * (1 + 0.1 * np.random.randn())

                # Epistatic interactions (pairwise)
                epistasis = 0
                for j, k in combinations(range(self.n_loci), 2):
                    if gt[j] == 1 and gt[k] == 1:
                        # Some pairs are synergistic, others antagonistic
                        if (j + k) % 3 == 0:
                            epistasis += 0.08 * conc / (conc + 1)  # synergistic
                        else:
                            epistasis -= 0.04  # antagonistic (fitness cost)

                fitness[i] = max(0.01, 1.0 - cost + benefit + epistasis)

            self.fitness_values[conc] = fitness

        return self.fitness_values

    def find_peaks(self, conc):
        """Find local fitness peaks at given drug concentration."""
        fitness = self.fitness_values[conc]
        peaks = []
        for i, gt in enumerate(self.genotypes):
            is_peak = True
            for j in range(self.n_loci):
                neighbor = list(gt)
                neighbor[j] = 1 - neighbor[j]
                neighbor_idx = self.genotypes.index(tuple(neighbor))
                if fitness[neighbor_idx] > fitness[i]:
                    is_peak = False
                    break
            if is_peak:
                peaks.append((gt, fitness[i]))
        return peaks

    def compute_ruggedness(self, conc):
        """Compute landscape ruggedness (correlation length)."""
        fitness = self.fitness_values[conc]
        correlations = []
        for d in range(1, self.n_loci + 1):
            pairs = []
            for i, gt1 in enumerate(self.genotypes):
                for j, gt2 in enumerate(self.genotypes):
                    if sum(a != b for a, b in zip(gt1, gt2)) == d:
                        pairs.append((fitness[i], fitness[j]))
            if len(pairs) > 1:
                f1, f2 = zip(*pairs)
                r, _ = stats.pearsonr(f1, f2)
                correlations.append(r)
        return correlations


# ==============================================================================
# Module 3: Evolutionary Path Prediction
# ==============================================================================

class EvolutionaryPathPredictor:
    """Enumerate and rank accessible evolutionary paths."""

    def __init__(self, landscape: FitnessLandscape):
        self.landscape = landscape

    def find_accessible_paths(self, conc, start=None, end=None):
        """Find all monotonically increasing fitness paths."""
        if start is None:
            start = tuple([0] * self.landscape.n_loci)
        if end is None:
            # Find global maximum
            fitness = self.landscape.fitness_values[conc]
            end_idx = np.argmax(fitness)
            end = self.landscape.genotypes[end_idx]

        paths = []
        self._dfs(conc, start, end, [start], paths)
        return paths

    def _dfs(self, conc, current, target, path, all_paths):
        if current == target:
            all_paths.append(list(path))
            return

        fitness = self.landscape.fitness_values[conc]
        current_fitness = fitness[self.landscape.genotypes.index(current)]

        for i in range(self.landscape.n_loci):
            if current[i] == 0 and target[i] == 1:
                neighbor = list(current)
                neighbor[i] = 1
                neighbor = tuple(neighbor)
                neighbor_fitness = fitness[self.landscape.genotypes.index(neighbor)]
                if neighbor_fitness > current_fitness:
                    path.append(neighbor)
                    self._dfs(conc, neighbor, target, path, all_paths)
                    path.pop()

    def rank_paths_by_probability(self, conc, paths):
        """Rank paths by transition probability (proportional to fitness gain)."""
        fitness = self.landscape.fitness_values[conc]
        path_probs = []
        for path in paths:
            prob = 1.0
            for k in range(len(path) - 1):
                current = path[k]
                current_idx = self.landscape.genotypes.index(current)
                current_fit = fitness[current_idx]

                # Find all beneficial neighbors
                neighbors = []
                for i in range(self.landscape.n_loci):
                    if current[i] == 0:
                        nb = list(current)
                        nb[i] = 1
                        nb = tuple(nb)
                        nb_idx = self.landscape.genotypes.index(nb)
                        if fitness[nb_idx] > current_fit:
                            neighbors.append((nb, fitness[nb_idx] - current_fit))

                total_gain = sum(g for _, g in neighbors)
                next_step = path[k + 1]
                next_idx = self.landscape.genotypes.index(next_step)
                step_gain = fitness[next_idx] - current_fit
                if total_gain > 0:
                    prob *= step_gain / total_gain

            path_probs.append(prob)
        return path_probs


# ==============================================================================
# Module 4: Horizontal Gene Transfer Network
# ==============================================================================

class HGTNetwork:
    """Model horizontal gene transfer dynamics."""

    def __init__(self, n_species=8, n_arg_types=5):
        self.n_species = n_species
        self.n_arg_types = n_arg_types
        self.species_names = ['E.coli', 'K.pneumoniae', 'E.faecium',
                              'S.aureus', 'A.baumannii', 'P.aeruginosa',
                              'E.cloacae', 'S.maltophilia'][:n_species]
        self.network = None
        self.transfer_history = []

    def build_transfer_network(self):
        """Build HGT network with species-specific transfer rates."""
        self.network = nx.DiGraph()

        for i, sp in enumerate(self.species_names):
            self.network.add_node(sp, args=set(), plasmid_count=np.random.poisson(3))

        # Transfer rates depend on phylogenetic distance and ecology
        for i, sp1 in enumerate(self.species_names):
            for j, sp2 in enumerate(self.species_names):
                if i != j:
                    # Closer species = higher transfer rate
                    base_rate = 0.01 * np.exp(-0.3 * abs(i - j))
                    # Gram-negative to gram-negative higher
                    if sp1 != 'S.aureus' and sp2 != 'S.aureus' and sp1 != 'E.faecium' and sp2 != 'E.faecium':
                        base_rate *= 2
                    self.network.add_edge(sp1, sp2, transfer_rate=base_rate)

        return self.network

    def simulate_transfer_dynamics(self, n_steps=500, initial_carriers=None):
        """Simulate HGT events over time."""
        if initial_carriers is None:
            # Seed some ARGs
            self.network.nodes['E.coli']['args'] = {0, 2}
            self.network.nodes['K.pneumoniae']['args'] = {1, 3}
            self.network.nodes['A.baumannii']['args'] = {0}

        arg_counts = {sp: [] for sp in self.species_names}
        transfer_events = []

        for t in range(n_steps):
            for sp in self.species_names:
                arg_counts[sp].append(len(self.network.nodes[sp]['args']))

            for u, v, data in self.network.edges(data=True):
                sender_args = self.network.nodes[u]['args']
                if len(sender_args) > 0:
                    for arg in list(sender_args):
                        if np.random.random() < data['transfer_rate']:
                            self.network.nodes[v]['args'].add(arg)
                            transfer_events.append({
                                'time': t,
                                'from': u,
                                'to': v,
                                'arg': arg
                            })

            # Random loss (plasmid curing)
            for sp in self.species_names:
                args = self.network.nodes[sp]['args']
                for arg in list(args):
                    if np.random.random() < 0.002:
                        args.discard(arg)

        self.transfer_history = transfer_events
        return arg_counts, transfer_events


# ==============================================================================
# Module 5: Spatiotemporal Dynamics Model (SIR + AMR)
# ==============================================================================

class SpatiotemporalAMRModel:
    """Coupled SIR-AMR model with antibiotic usage feedback."""

    def __init__(self, n_regions=4):
        self.n_regions = n_regions
        self.region_names = [f'Region_{i+1}' for i in range(n_regions)]

    def sir_amr_ode(self, y, t, params):
        """ODE system: S, I_s (sensitive), I_r (resistant), R for each region."""
        n = self.n_regions
        beta_s, beta_r, gamma, mu, tau, sigma, c = params[:7]
        migration = params[7]

        dydt = np.zeros(4 * n)

        for i in range(n):
            S = y[4*i]
            Is = y[4*i + 1]
            Ir = y[4*i + 2]
            R = y[4*i + 3]
            N = S + Is + Ir + R

            # Antibiotic usage drives resistance conversion
            usage_i = 0.3 + 0.2 * np.sin(2 * np.pi * t / 365)  # seasonal

            dS = mu * N - beta_s * S * Is / N - beta_r * S * Ir / N - mu * S
            dIs = beta_s * S * Is / N - gamma * Is - mu * Is - tau * usage_i * Is
            dIr = beta_r * S * Ir / N - gamma * Ir - mu * Ir + tau * usage_i * Is + sigma * Is
            dR_dt = gamma * (Is + Ir) - mu * R

            # Migration between regions
            for j in range(n):
                if j != i:
                    dS += migration * (y[4*j] - S) / n
                    dIs += migration * (y[4*j+1] - Is) / n
                    dIr += migration * (y[4*j+2] - Ir) / n

            dydt[4*i] = dS
            dydt[4*i + 1] = dIs
            dydt[4*i + 2] = dIr
            dydt[4*i + 3] = dR_dt

        return dydt

    def simulate(self, T=730, dt=1):
        """Run spatiotemporal simulation."""
        t = np.linspace(0, T, int(T/dt))

        # Initial conditions per region
        y0 = []
        for i in range(self.n_regions):
            pop = 100000 * (1 + 0.2 * i)
            y0.extend([pop * 0.95, pop * 0.03, pop * 0.01, pop * 0.01])

        params = [
            0.3,    # beta_s (transmission rate, sensitive)
            0.25,   # beta_r (transmission rate, resistant - fitness cost)
            0.1,    # gamma (recovery rate)
            0.0001, # mu (birth/death rate)
            0.05,   # tau (treatment-driven resistance conversion)
            0.001,  # sigma (spontaneous resistance)
            0.5,    # c (cost parameter)
            0.001,  # migration rate
        ]

        sol = integrate.odeint(self.sir_amr_ode, y0, t, args=(params,), mxstep=10000)
        return t, sol


# ==============================================================================
# Module 6: Treatment Strategy Optimization
# ==============================================================================

class TreatmentOptimizer:
    """Optimize antibiotic treatment strategies."""

    def __init__(self, n_antibiotics=3):
        self.n_antibiotics = n_antibiotics
        self.ab_names = ['Drug_A', 'Drug_B', 'Drug_C'][:n_antibiotics]

    def simulate_population(self, strategy, T=200, N0=1e6, n_steps=200):
        """Simulate bacterial population under treatment strategy."""
        dt = T / n_steps
        t = np.linspace(0, T, n_steps)

        # Define resistance genotypes (2^n_antibiotics)
        n_genotypes = 2 ** self.n_antibiotics
        genotypes = list(product([0, 1], repeat=self.n_antibiotics))

        # Population of each genotype
        pop = np.zeros((n_steps, n_genotypes))
        pop[0, 0] = N0 * 0.97  # mostly wild-type
        for i in range(1, n_genotypes):
            pop[0, i] = N0 * 0.03 / (n_genotypes - 1)

        # Fitness parameters
        growth_rate = 0.5
        carrying_capacity = 1e7
        mutation_rate = 1e-6

        for step in range(1, n_steps):
            current_drug = strategy(t[step], T)

            for i, gt in enumerate(genotypes):
                # Kill rate depends on drug activity vs resistance
                kill_rate = 0
                for d in range(self.n_antibiotics):
                    if current_drug[d] > 0 and gt[d] == 0:
                        kill_rate += current_drug[d] * 0.8
                    elif current_drug[d] > 0 and gt[d] == 1:
                        kill_rate += current_drug[d] * 0.1  # reduced efficacy

                # Fitness cost of resistance
                cost = 0.05 * sum(gt)
                total_pop = np.sum(pop[step-1])
                growth = growth_rate * (1 - cost) * (1 - total_pop / carrying_capacity)
                net_rate = growth - kill_rate

                pop[step, i] = pop[step-1, i] * np.exp(net_rate * dt)
                pop[step, i] = max(pop[step, i], 0)

            # Mutation events
            for i, gt in enumerate(genotypes):
                for j in range(self.n_antibiotics):
                    if gt[j] == 0:
                        mutant = list(gt)
                        mutant[j] = 1
                        mut_idx = genotypes.index(tuple(mutant))
                        transfer = pop[step, i] * mutation_rate * dt
                        pop[step, mut_idx] += transfer

        return t, pop, genotypes

    def monotherapy_strategy(self, drug_idx):
        def strategy(t, T):
            drugs = [0.0] * self.n_antibiotics
            drugs[drug_idx] = 1.0
            return drugs
        return strategy

    def cycling_strategy(self, cycle_length=30):
        def strategy(t, T):
            drugs = [0.0] * self.n_antibiotics
            phase = int(t / cycle_length) % self.n_antibiotics
            drugs[phase] = 1.0
            return drugs
        return strategy

    def combination_strategy(self):
        def strategy(t, T):
            return [0.5] * self.n_antibiotics
        return strategy

    def adaptive_strategy(self):
        """Simplified adaptive strategy - switch when resistance fraction is high."""
        state = {'current_drug': 0, 'last_switch': 0}
        def strategy(t, T):
            drugs = [0.0] * self.n_antibiotics
            if t - state['last_switch'] > 40:
                state['current_drug'] = (state['current_drug'] + 1) % self.n_antibiotics
                state['last_switch'] = t
            drugs[state['current_drug']] = 1.0
            return drugs
        return strategy


# ==============================================================================
# Module 7: Population Genetics Simulation (Wright-Fisher)
# ==============================================================================

class PopulationGeneticsSim:
    """Wright-Fisher model for resistance allele dynamics."""

    def __init__(self, N=1000, n_loci=4):
        self.N = N
        self.n_loci = n_loci

    def simulate(self, generations=500, s_vals=None, drug_on_off=None):
        """Simulate allele frequency trajectories."""
        if s_vals is None:
            s_vals = [0.02, 0.01, 0.03, 0.015][:self.n_loci]

        freqs = np.zeros((generations, self.n_loci))
        freqs[0, :] = 0.01  # initial frequency

        for gen in range(1, generations):
            for i in range(self.n_loci):
                p = freqs[gen-1, i]
                # Selection coefficient depends on drug presence
                if drug_on_off is not None and drug_on_off(gen):
                    s = s_vals[i]
                else:
                    s = -0.01 * (i + 1) * 0.5  # fitness cost without drug

                # Selection
                p_prime = p * (1 + s) / (1 + s * p)
                # Drift
                p_prime = np.random.binomial(2 * self.N, p_prime) / (2 * self.N)
                freqs[gen, i] = np.clip(p_prime, 0, 1)

        return freqs


# ==============================================================================
# Visualization Functions
# ==============================================================================

def plot_arg_detection(pipeline, stats):
    """Figure 1: ARG detection results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1a: Resistance rates by antibiotic
    ax = axes[0, 0]
    rates = stats['resistance_rates']
    bars = ax.bar(range(len(rates)), list(rates.values()),
                  color=['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6'])
    ax.set_xticks(range(len(rates)))
    ax.set_xticklabels([k.replace('_', '\n') for k in rates.keys()], fontsize=9)
    ax.set_ylabel('Resistance Rate')
    ax.set_title('(A) Predicted Resistance Rates by Antibiotic Class')
    ax.set_ylim(0, 0.5)
    for bar, val in zip(bars, rates.values()):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                f'{val:.2f}', ha='center', fontsize=9)

    # 1b: Species distribution
    ax = axes[0, 1]
    sp = stats['species_distribution']
    ax.pie(sp.values(), labels=sp.keys(), autopct='%1.1f%%',
           colors=['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6'])
    ax.set_title('(B) Species Distribution')

    # 1c: ARG count distribution
    ax = axes[1, 0]
    n_args = [len(g['args']) for g in pipeline.genomes]
    ax.hist(n_args, bins=range(0, max(n_args)+2), color='#3498db',
            edgecolor='black', alpha=0.7)
    ax.set_xlabel('Number of ARGs per Genome')
    ax.set_ylabel('Count')
    ax.set_title('(C) Distribution of ARGs per Genome')
    ax.axvline(np.mean(n_args), color='red', linestyle='--',
               label=f'Mean={np.mean(n_args):.1f}')
    ax.legend()

    # 1d: Co-resistance heatmap
    ax = axes[1, 1]
    profiles = pipeline.arg_profiles
    n_ab = profiles.shape[1]
    co_res = np.zeros((n_ab, n_ab))
    for i in range(n_ab):
        for j in range(n_ab):
            both = np.sum((profiles[:, i] == 1) & (profiles[:, j] == 1))
            co_res[i, j] = both / len(profiles)
    im = ax.imshow(co_res, cmap='YlOrRd', vmin=0)
    ax.set_xticks(range(n_ab))
    ax.set_yticks(range(n_ab))
    labels = [a.replace('_', '\n') for a in pipeline.ANTIBIOTICS]
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title('(D) Co-resistance Matrix')
    plt.colorbar(im, ax=ax, label='Co-resistance Rate')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig1_arg_detection.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig1_arg_detection.png")


def plot_fitness_landscape(landscape):
    """Figure 2: Fitness landscape visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 2a: Fitness vs number of mutations at different concentrations
    ax = axes[0, 0]
    for conc in [0, 1.0, 4.0, 8.0]:
        fitness = landscape.fitness_values[conc]
        n_muts = [sum(gt) for gt in landscape.genotypes]
        means = [np.mean([f for f, n in zip(fitness, n_muts) if n == k])
                 for k in range(landscape.n_loci + 1)]
        ax.plot(range(landscape.n_loci + 1), means, 'o-',
                label=f'[Drug]={conc}', markersize=6)
    ax.set_xlabel('Number of Resistance Mutations')
    ax.set_ylabel('Mean Fitness')
    ax.set_title('(A) Fitness vs Mutation Number')
    ax.legend()

    # 2b: Landscape ruggedness
    ax = axes[0, 1]
    for conc in [0, 1.0, 4.0, 8.0]:
        corr = landscape.compute_ruggedness(conc)
        ax.plot(range(1, len(corr)+1), corr, 'o-', label=f'[Drug]={conc}')
    ax.set_xlabel('Hamming Distance')
    ax.set_ylabel('Fitness Correlation')
    ax.set_title('(B) Landscape Ruggedness')
    ax.legend()
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

    # 2c: Fitness distribution
    ax = axes[1, 0]
    for conc in [0, 2.0, 8.0]:
        fitness = landscape.fitness_values[conc]
        ax.hist(fitness, bins=15, alpha=0.5, label=f'[Drug]={conc}')
    ax.set_xlabel('Fitness')
    ax.set_ylabel('Count')
    ax.set_title('(C) Fitness Distribution')
    ax.legend()

    # 2d: Number of peaks vs drug concentration
    ax = axes[1, 1]
    concs = landscape.drug_concs
    n_peaks = [len(landscape.find_peaks(c)) for c in concs]
    ax.plot(concs, n_peaks, 'rs-', markersize=8, linewidth=2)
    ax.set_xlabel('Drug Concentration')
    ax.set_ylabel('Number of Local Peaks')
    ax.set_title('(D) Landscape Complexity')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig2_fitness_landscape.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig2_fitness_landscape.png")


def plot_evolutionary_paths(landscape, predictor):
    """Figure 3: Evolutionary path analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 3a: Number of accessible paths at different concentrations
    ax = axes[0]
    concs = landscape.drug_concs[1:]  # skip 0
    n_paths_list = []
    for conc in concs:
        paths = predictor.find_accessible_paths(conc)
        n_paths_list.append(len(paths))
    ax.bar(range(len(concs)), n_paths_list, color='#2ecc71',
           edgecolor='black', alpha=0.8)
    ax.set_xticks(range(len(concs)))
    ax.set_xticklabels([f'{c}' for c in concs])
    ax.set_xlabel('Drug Concentration')
    ax.set_ylabel('Number of Accessible Paths')
    ax.set_title('(A) Accessible Evolutionary Paths')

    # 3b: Path probability distribution at concentration 4.0
    ax = axes[1]
    conc = 4.0
    paths = predictor.find_accessible_paths(conc)
    if len(paths) > 0:
        probs = predictor.rank_paths_by_probability(conc, paths)
        sorted_probs = sorted(probs, reverse=True)
        top_n = min(15, len(sorted_probs))
        ax.bar(range(top_n), sorted_probs[:top_n], color='#e74c3c',
               edgecolor='black', alpha=0.8)
        ax.set_xlabel('Path Rank')
        ax.set_ylabel('Path Probability')
        ax.set_title(f'(B) Path Probability Distribution ([Drug]={conc})')
    else:
        ax.text(0.5, 0.5, 'No accessible paths found', transform=ax.transAxes,
                ha='center')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig3_evolutionary_paths.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig3_evolutionary_paths.png")


def plot_hgt_network(hgt):
    """Figure 4: HGT network visualization."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 4a: Network graph
    ax = axes[0]
    G = hgt.network
    pos = nx.spring_layout(G, seed=42, k=2)
    edge_weights = [G[u][v]['transfer_rate'] * 100 for u, v in G.edges()]
    node_sizes = [300 + len(G.nodes[n]['args']) * 200 for n in G.nodes()]
    node_colors = [len(G.nodes[n]['args']) for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color=node_colors, cmap='YlOrRd',
                           edgecolors='black', linewidths=1)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_weights,
                           alpha=0.3, edge_color='gray',
                           arrows=True, arrowsize=10)
    ax.set_title('(A) HGT Network (node size ∝ ARGs)')
    ax.axis('off')

    # 4b: ARG spread over time
    ax = axes[1]
    arg_counts, _ = hgt.simulate_transfer_dynamics()
    for sp, counts in arg_counts.items():
        ax.plot(counts, label=sp, linewidth=1.5)
    ax.set_xlabel('Time Steps')
    ax.set_ylabel('Number of ARG Types')
    ax.set_title('(B) ARG Spread via HGT')
    ax.legend(fontsize=7, ncol=2)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig4_hgt_network.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig4_hgt_network.png")


def plot_spatiotemporal(model):
    """Figure 5: Spatiotemporal dynamics."""
    t, sol = model.simulate()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']

    for i in range(model.n_regions):
        ax = axes[i // 2, i % 2]
        S = sol[:, 4*i]
        Is = sol[:, 4*i + 1]
        Ir = sol[:, 4*i + 2]
        R = sol[:, 4*i + 3]

        ax.plot(t, S, label='Susceptible', color='#3498db', linewidth=1.5)
        ax.plot(t, Is, label='Infected (sensitive)', color='#2ecc71', linewidth=1.5)
        ax.plot(t, Ir, label='Infected (resistant)', color='#e74c3c', linewidth=1.5)
        ax.plot(t, R, label='Recovered', color='#95a5a6', linewidth=1.5)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Population')
        ax.set_title(f'Region {i+1}')
        ax.legend(fontsize=7)
        ax.set_xlim(0, 730)

    plt.suptitle('Spatiotemporal SIR-AMR Dynamics', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig5_spatiotemporal.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig5_spatiotemporal.png")


def plot_resistance_fraction(model):
    """Figure 5b: Resistance fraction across regions."""
    t, sol = model.simulate()
    fig, ax = plt.subplots(figsize=(10, 6))

    for i in range(model.n_regions):
        Is = sol[:, 4*i + 1]
        Ir = sol[:, 4*i + 2]
        total_infected = Is + Ir
        frac_r = np.where(total_infected > 1, Ir / total_infected, 0)
        ax.plot(t, frac_r, linewidth=2, label=f'Region {i+1}')

    ax.set_xlabel('Time (days)', fontsize=12)
    ax.set_ylabel('Resistance Fraction', fontsize=12)
    ax.set_title('Temporal Dynamics of Resistance Fraction Across Regions', fontsize=13)
    ax.legend()
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig5b_resistance_fraction.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig5b_resistance_fraction.png")


def plot_treatment_strategies(optimizer):
    """Figure 6: Treatment strategy comparison."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    strategies = {
        'Monotherapy (Drug A)': optimizer.monotherapy_strategy(0),
        'Cycling (30-day)': optimizer.cycling_strategy(30),
        'Combination': optimizer.combination_strategy(),
        'Adaptive Switching': optimizer.adaptive_strategy(),
    }

    colors_gt = plt.cm.Set2(np.linspace(0, 1, 8))
    results_summary = {}

    for idx, (name, strat) in enumerate(strategies.items()):
        ax = axes[idx // 2, idx % 2]
        t, pop, genotypes = optimizer.simulate_population(strat)

        total_pop = np.sum(pop, axis=1)
        resistant_pop = np.sum(pop[:, 1:], axis=1)  # all non-wildtype

        ax.semilogy(t, pop[:, 0], linewidth=2, label='Wild-type', color='#2ecc71')
        ax.semilogy(t, resistant_pop, linewidth=2, label='Total Resistant', color='#e74c3c')
        ax.semilogy(t, total_pop, linewidth=2, label='Total', color='#3498db', linestyle='--')

        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Population (log scale)')
        ax.set_title(name)
        ax.legend(fontsize=8)
        ax.set_ylim(1, 1e8)

        # Store final resistance fraction
        final_frac = resistant_pop[-1] / (total_pop[-1] + 1e-10) if total_pop[-1] > 0 else 0
        results_summary[name] = {
            'final_total': float(total_pop[-1]),
            'final_resistant': float(resistant_pop[-1]),
            'resistance_fraction': float(final_frac),
        }

    plt.suptitle('Comparison of Antibiotic Treatment Strategies', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig6_treatment_strategies.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig6_treatment_strategies.png")
    return results_summary


def plot_population_genetics(sim):
    """Figure 7: Population genetics simulation."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 7a: Constant drug pressure
    ax = axes[0]
    freqs = sim.simulate(generations=500, drug_on_off=lambda g: True)
    for i in range(sim.n_loci):
        ax.plot(freqs[:, i], linewidth=1.5, label=f'Locus {chr(65+i)}')
    ax.set_xlabel('Generation')
    ax.set_ylabel('Allele Frequency')
    ax.set_title('(A) Constant Drug Pressure')
    ax.legend()
    ax.set_ylim(0, 1)

    # 7b: Pulsed drug pressure
    ax = axes[1]
    freqs = sim.simulate(generations=500,
                         drug_on_off=lambda g: (g // 50) % 2 == 0)
    for i in range(sim.n_loci):
        ax.plot(freqs[:, i], linewidth=1.5, label=f'Locus {chr(65+i)}')
    ax.set_xlabel('Generation')
    ax.set_ylabel('Allele Frequency')
    ax.set_title('(B) Pulsed Drug Pressure (50-gen cycles)')
    ax.legend()
    ax.set_ylim(0, 1)

    # Add shading for drug-off periods
    for g_start in range(50, 500, 100):
        ax.axvspan(g_start, min(g_start+50, 500), alpha=0.1, color='gray')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig7_population_genetics.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig7_population_genetics.png")


def plot_strategy_comparison_summary(results):
    """Figure 8: Strategy comparison bar chart."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = list(results.keys())
    short_names = ['Mono', 'Cycling', 'Combo', 'Adaptive']
    fracs = [results[n]['resistance_fraction'] for n in names]
    totals = [results[n]['final_total'] for n in names]

    colors = ['#e74c3c', '#f39c12', '#3498db', '#2ecc71']

    ax = axes[0]
    bars = ax.bar(range(len(names)), fracs, color=colors, edgecolor='black')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(short_names)
    ax.set_ylabel('Final Resistance Fraction')
    ax.set_title('(A) Resistance Fraction at End of Treatment')
    for bar, val in zip(bars, fracs):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', fontsize=9)

    ax = axes[1]
    bars = ax.bar(range(len(names)), [np.log10(t+1) for t in totals],
                  color=colors, edgecolor='black')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(short_names)
    ax.set_ylabel('log₁₀(Final Population)')
    ax.set_title('(B) Total Bacterial Population at End')
    for bar, val in zip(bars, totals):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
                f'{val:.1e}', ha='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig8_strategy_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved fig8_strategy_comparison.png")


# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    print("=" * 70)
    print("AMR Evolution Prediction Framework")
    print("=" * 70)

    # --- Module 1: ARG Detection ---
    print("\n[1/7] Running ARG Detection Pipeline...")
    pipeline = ARGDetectionPipeline(n_genomes=200)
    pipeline.simulate_genomes()
    pipeline.compute_resistance_profiles()
    stats = pipeline.compute_statistics()
    print(f"  Genomes analyzed: {stats['total_genomes']}")
    print(f"  Mean ARGs/genome: {stats['mean_args']:.2f}")
    print(f"  MDR isolates: {stats['mdr_count']}")
    print(f"  XDR isolates: {stats['xdr_count']}")
    plot_arg_detection(pipeline, stats)

    # --- Module 2: Fitness Landscape ---
    print("\n[2/7] Constructing Fitness Landscape...")
    landscape = FitnessLandscape(n_loci=5)
    landscape.construct_landscape()
    for conc in [0, 4.0, 8.0]:
        peaks = landscape.find_peaks(conc)
        print(f"  [Drug]={conc}: {len(peaks)} local peaks")
    plot_fitness_landscape(landscape)

    # --- Module 3: Evolutionary Paths ---
    print("\n[3/7] Predicting Evolutionary Paths...")
    predictor = EvolutionaryPathPredictor(landscape)
    for conc in [1.0, 4.0, 8.0]:
        paths = predictor.find_accessible_paths(conc)
        print(f"  [Drug]={conc}: {len(paths)} accessible paths")
    plot_evolutionary_paths(landscape, predictor)

    # --- Module 4: HGT Network ---
    print("\n[4/7] Modeling HGT Network...")
    hgt = HGTNetwork(n_species=8, n_arg_types=5)
    hgt.build_transfer_network()
    print(f"  Network nodes: {hgt.network.number_of_nodes()}")
    print(f"  Network edges: {hgt.network.number_of_edges()}")
    plot_hgt_network(hgt)

    # --- Module 5: Spatiotemporal Dynamics ---
    print("\n[5/7] Running Spatiotemporal SIR-AMR Model...")
    st_model = SpatiotemporalAMRModel(n_regions=4)
    t, sol = st_model.simulate()
    for i in range(st_model.n_regions):
        Is_final = sol[-1, 4*i+1]
        Ir_final = sol[-1, 4*i+2]
        total = Is_final + Ir_final
        frac = Ir_final / total if total > 0 else 0
        print(f"  Region {i+1} final resistance fraction: {frac:.4f}")
    plot_spatiotemporal(st_model)
    plot_resistance_fraction(st_model)

    # --- Module 6: Treatment Optimization ---
    print("\n[6/7] Optimizing Treatment Strategies...")
    optimizer = TreatmentOptimizer(n_antibiotics=3)
    results = plot_treatment_strategies(optimizer)
    for name, res in results.items():
        print(f"  {name}: resistance fraction = {res['resistance_fraction']:.4f}")
    plot_strategy_comparison_summary(results)

    # --- Module 7: Population Genetics ---
    print("\n[7/7] Running Population Genetics Simulation...")
    pg_sim = PopulationGeneticsSim(N=1000, n_loci=4)
    plot_population_genetics(pg_sim)

    # Save results
    all_results = {
        'arg_detection': stats,
        'fitness_landscape': {
            'n_loci': landscape.n_loci,
            'n_genotypes': landscape.n_genotypes,
            'peaks_by_conc': {str(c): len(landscape.find_peaks(c)) for c in landscape.drug_concs},
        },
        'evolutionary_paths': {
            str(c): len(predictor.find_accessible_paths(c))
            for c in [1.0, 4.0, 8.0]
        },
        'hgt_network': {
            'n_nodes': hgt.network.number_of_nodes(),
            'n_edges': hgt.network.number_of_edges(),
            'n_transfer_events': len(hgt.transfer_history),
        },
        'treatment_strategies': results,
    }

    with open('results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print("\nResults saved to results.json")
    print(f"Figures saved to {FIGURES_DIR}/")

    return all_results


if __name__ == '__main__':
    results = main()
