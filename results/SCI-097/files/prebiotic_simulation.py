#!/usr/bin/env python3
"""
Integrated Prebiotic Chemistry Simulation Framework
====================================================
Combines stochastic chemical kinetics (Gillespie SSA / Chemical Master Equation),
reaction network analysis, and protocell formation modeling to simulate
chemical evolution scenarios relevant to the origin of life.

Modules:
1. Primordial Soup (Extended Miller-Urey) reaction network
2. RNA World self-replication dynamics
3. Metabolism-First (hydrothermal vent) autocatalytic cycles
4. Stochastic Chemical Kinetics (CME / Gillespie SSA)
5. Membrane self-assembly and protocell formation
6. Enceladus/Titan environmental adaptation
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import odeint
from scipy.stats import poisson
import networkx as nx
from collections import defaultdict
import json
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})

FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)
RESULTS = {}

# =============================================================================
# 1. PRIMORDIAL SOUP — Extended Miller-Urey Reaction Network
# =============================================================================

class MillerUreyExtended:
    """
    Extended Miller-Urey reaction network modeling the synthesis of
    amino acids, nucleobases, and sugars from simple precursors
    (CH4, NH3, H2O, H2, HCN) under energy input (lightning/UV).
    """
    
    SPECIES = [
        'CH4', 'NH3', 'H2O', 'H2', 'HCN', 'HCHO',  # precursors
        'Glycine', 'Alanine', 'Aspartate', 'Valine',  # amino acids
        'Adenine', 'Guanine', 'Cytosine', 'Uracil',   # nucleobases
        'Ribose', 'Glycerol',                          # sugars/lipid precursors
        'Formate', 'Acetate',                          # organic acids
        'Cyanamide', 'Urea',                           # intermediates
    ]
    
    REACTIONS = [
        # (reactants, products, rate_constant)
        (['CH4', 'NH3'], ['HCN', 'H2', 'H2'], 0.005),
        (['HCN', 'HCN', 'HCN', 'HCN', 'HCN'], ['Adenine'], 0.0001),
        (['HCN', 'H2O'], ['HCHO', 'NH3'], 0.01),
        (['HCHO', 'HCHO', 'HCHO', 'HCHO', 'HCHO'], ['Ribose'], 0.00005),
        (['HCN', 'H2O', 'CH4'], ['Glycine', 'H2'], 0.008),
        (['HCN', 'CH4', 'H2O', 'H2'], ['Alanine'], 0.004),
        (['HCN', 'HCHO', 'NH3'], ['Glycine'], 0.006),
        (['HCN', 'H2O', 'HCHO'], ['Aspartate'], 0.002),
        (['CH4', 'NH3', 'H2O'], ['Valine'], 0.001),
        (['HCN', 'HCN', 'HCN'], ['Cytosine'], 0.0003),
        (['HCN', 'HCHO', 'H2O'], ['Uracil'], 0.0002),
        (['HCN', 'HCN', 'HCN', 'HCN'], ['Guanine'], 0.00015),
        (['HCHO', 'HCHO', 'HCHO'], ['Glycerol'], 0.001),
        (['HCN', 'H2O'], ['Formate', 'NH3'], 0.015),
        (['CH4', 'H2O'], ['Acetate', 'H2'], 0.003),
        (['HCN', 'NH3'], ['Cyanamide', 'H2'], 0.007),
        (['Cyanamide', 'H2O'], ['Urea'], 0.01),
    ]
    
    def __init__(self, initial_conc=None):
        self.species_idx = {s: i for i, s in enumerate(self.SPECIES)}
        n = len(self.SPECIES)
        if initial_conc is None:
            self.initial = np.zeros(n)
            self.initial[self.species_idx['CH4']] = 100.0
            self.initial[self.species_idx['NH3']] = 80.0
            self.initial[self.species_idx['H2O']] = 500.0
            self.initial[self.species_idx['H2']] = 50.0
            self.initial[self.species_idx['HCN']] = 20.0
            self.initial[self.species_idx['HCHO']] = 10.0
        else:
            self.initial = initial_conc
    
    def build_network_graph(self):
        G = nx.DiGraph()
        for s in self.SPECIES:
            G.add_node(s)
        for reactants, products, k in self.REACTIONS:
            for r in set(reactants):
                for p in set(products):
                    G.add_edge(r, p, weight=k)
        return G
    
    def ode_system(self, y, t):
        dydt = np.zeros(len(self.SPECIES))
        for reactants, products, k in self.REACTIONS:
            rate = k
            for r in reactants:
                idx = self.species_idx[r]
                rate *= max(y[idx], 0)
            for r in reactants:
                idx = self.species_idx[r]
                dydt[idx] -= rate
            for p in products:
                idx = self.species_idx[p]
                dydt[idx] += rate
        return dydt
    
    def simulate(self, t_span=(0, 500), n_points=2000):
        t = np.linspace(t_span[0], t_span[1], n_points)
        sol = odeint(self.ode_system, self.initial, t, mxstep=10000)
        sol = np.maximum(sol, 0)
        return t, sol
    
    def run_and_plot(self):
        t, sol = self.simulate()
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Extended Miller-Urey Reaction Network Simulation', fontsize=14, fontweight='bold')
        
        # Precursors
        ax = axes[0, 0]
        for s in ['CH4', 'NH3', 'H2O', 'HCN', 'HCHO']:
            ax.plot(t, sol[:, self.species_idx[s]], label=s, linewidth=1.5)
        ax.set_xlabel('Time (arbitrary units)')
        ax.set_ylabel('Concentration')
        ax.set_title('(a) Precursor Depletion')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # Amino acids
        ax = axes[0, 1]
        for s in ['Glycine', 'Alanine', 'Aspartate', 'Valine']:
            ax.plot(t, sol[:, self.species_idx[s]], label=s, linewidth=1.5)
        ax.set_xlabel('Time (arbitrary units)')
        ax.set_ylabel('Concentration')
        ax.set_title('(b) Amino Acid Synthesis')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # Nucleobases
        ax = axes[1, 0]
        for s in ['Adenine', 'Guanine', 'Cytosine', 'Uracil']:
            ax.plot(t, sol[:, self.species_idx[s]], label=s, linewidth=1.5)
        ax.set_xlabel('Time (arbitrary units)')
        ax.set_ylabel('Concentration')
        ax.set_title('(c) Nucleobase Formation')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # Sugars & intermediates
        ax = axes[1, 1]
        for s in ['Ribose', 'Glycerol', 'Formate', 'Acetate', 'Urea']:
            ax.plot(t, sol[:, self.species_idx[s]], label=s, linewidth=1.5)
        ax.set_xlabel('Time (arbitrary units)')
        ax.set_ylabel('Concentration')
        ax.set_title('(d) Sugars & Intermediates')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'miller_urey_extended.png')
        fig.savefig(path)
        plt.close(fig)
        
        # Final concentrations
        final = {s: float(sol[-1, self.species_idx[s]]) for s in self.SPECIES}
        
        # Network analysis
        G = self.build_network_graph()
        fig2, ax2 = plt.subplots(1, 1, figsize=(12, 10))
        pos = nx.spring_layout(G, k=2.0, seed=42)
        
        node_colors = []
        for n in G.nodes():
            if n in ['CH4', 'NH3', 'H2O', 'H2', 'HCN', 'HCHO']:
                node_colors.append('#3498db')
            elif n in ['Glycine', 'Alanine', 'Aspartate', 'Valine']:
                node_colors.append('#e74c3c')
            elif n in ['Adenine', 'Guanine', 'Cytosine', 'Uracil']:
                node_colors.append('#2ecc71')
            else:
                node_colors.append('#f39c12')
        
        nx.draw(G, pos, ax=ax2, with_labels=True, node_color=node_colors,
                node_size=1200, font_size=8, font_weight='bold',
                edge_color='gray', arrows=True, arrowsize=15,
                connectionstyle='arc3,rad=0.1')
        ax2.set_title('Miller-Urey Extended Reaction Network', fontsize=14, fontweight='bold')
        
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#3498db', label='Precursors'),
            Patch(facecolor='#e74c3c', label='Amino Acids'),
            Patch(facecolor='#2ecc71', label='Nucleobases'),
            Patch(facecolor='#f39c12', label='Other Products'),
        ]
        ax2.legend(handles=legend_elements, loc='upper left', fontsize=9)
        
        path2 = os.path.join(FIGURES_DIR, 'reaction_network.png')
        fig2.savefig(path2)
        plt.close(fig2)
        
        return {
            'final_concentrations': final,
            'network_nodes': G.number_of_nodes(),
            'network_edges': G.number_of_edges(),
            'network_density': nx.density(G),
        }


# =============================================================================
# 2. RNA WORLD — Self-Replicating RNA Emergence
# =============================================================================

class RNAWorldSimulation:
    """
    Stochastic simulation of RNA World dynamics:
    - Template-directed polymerization
    - Ribozyme catalysis emergence
    - Error-prone replication with mutation
    - Selection via differential degradation
    """
    
    def __init__(self, n_sequences=200, seq_length=20, n_steps=5000):
        self.n_sequences = n_sequences
        self.seq_length = seq_length
        self.n_steps = n_steps
        self.mutation_rate = 0.02
        self.replication_base_rate = 0.1
        self.degradation_rate = 0.05
        self.catalytic_threshold = 0.6  # fitness threshold for catalytic activity
        self.alphabet = ['A', 'U', 'G', 'C']
    
    def generate_random_sequence(self):
        return np.random.choice(4, size=self.seq_length)
    
    def compute_fitness(self, seq):
        """Fitness based on GC content, secondary structure potential, and catalytic motifs."""
        gc_content = np.mean((seq == 2) | (seq == 3))
        # Check for stem-loop potential (palindromic subsequences)
        stem_score = 0
        for i in range(len(seq) - 5):
            complement = 3 - seq[i:i+3]
            if np.array_equal(complement[::-1], seq[i+3:i+6]):
                stem_score += 0.1
        # Catalytic motif (consecutive GUG or similar)
        motif_score = 0
        for i in range(len(seq) - 2):
            if seq[i] == 2 and seq[i+1] == 1 and seq[i+2] == 2:  # GUG
                motif_score += 0.15
            if seq[i] == 2 and seq[i+1] == 0 and seq[i+2] == 0:  # GAA
                motif_score += 0.1
        
        fitness = 0.3 * gc_content + 0.4 * min(stem_score, 0.5) + 0.3 * min(motif_score, 0.5)
        return fitness
    
    def replicate_with_error(self, seq):
        new_seq = seq.copy()
        for i in range(len(new_seq)):
            if np.random.random() < self.mutation_rate:
                new_seq[i] = np.random.randint(0, 4)
        return new_seq
    
    def simulate(self):
        population = [self.generate_random_sequence() for _ in range(self.n_sequences)]
        
        history = {
            'mean_fitness': [],
            'max_fitness': [],
            'population_size': [],
            'catalytic_fraction': [],
            'diversity': [],
            'gc_content': [],
        }
        
        for step in range(self.n_steps):
            fitnesses = np.array([self.compute_fitness(s) for s in population])
            
            catalytic = np.sum(fitnesses > self.catalytic_threshold) / max(len(population), 1)
            
            # Diversity (mean pairwise Hamming distance of a sample)
            if len(population) > 1:
                sample_size = min(20, len(population))
                sample_idx = np.random.choice(len(population), sample_size, replace=False)
                dists = []
                for i in range(sample_size):
                    for j in range(i+1, sample_size):
                        dists.append(np.mean(population[sample_idx[i]] != population[sample_idx[j]]))
                diversity = np.mean(dists) if dists else 0
            else:
                diversity = 0
            
            gc = np.mean([(np.mean((s == 2) | (s == 3))) for s in population]) if population else 0
            
            history['mean_fitness'].append(float(np.mean(fitnesses)) if len(fitnesses) > 0 else 0)
            history['max_fitness'].append(float(np.max(fitnesses)) if len(fitnesses) > 0 else 0)
            history['population_size'].append(len(population))
            history['catalytic_fraction'].append(float(catalytic))
            history['diversity'].append(float(diversity))
            history['gc_content'].append(float(gc))
            
            # Replication (fitness-proportional)
            new_population = []
            for i, seq in enumerate(population):
                rep_rate = self.replication_base_rate * (1 + 2 * fitnesses[i])
                if np.random.random() < rep_rate:
                    new_population.append(self.replicate_with_error(seq))
            
            # Degradation
            surviving = []
            for i, seq in enumerate(population):
                deg_rate = self.degradation_rate * (1.5 - fitnesses[i])
                if np.random.random() > deg_rate:
                    surviving.append(seq)
            
            population = surviving + new_population
            
            # Carrying capacity
            if len(population) > 500:
                fitnesses_all = np.array([self.compute_fitness(s) for s in population])
                top_idx = np.argsort(fitnesses_all)[-500:]
                population = [population[i] for i in top_idx]
            
            if len(population) == 0:
                population = [self.generate_random_sequence() for _ in range(10)]
        
        return history
    
    def run_and_plot(self):
        history = self.simulate()
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        fig.suptitle('RNA World: Self-Replicator Emergence Simulation', fontsize=14, fontweight='bold')
        
        steps = range(len(history['mean_fitness']))
        
        axes[0, 0].plot(steps, history['mean_fitness'], 'b-', alpha=0.7, label='Mean')
        axes[0, 0].plot(steps, history['max_fitness'], 'r-', alpha=0.7, label='Max')
        axes[0, 0].set_title('(a) Fitness Evolution')
        axes[0, 0].set_xlabel('Generation')
        axes[0, 0].set_ylabel('Fitness')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        axes[0, 1].plot(steps, history['population_size'], 'g-', linewidth=1.5)
        axes[0, 1].set_title('(b) Population Size')
        axes[0, 1].set_xlabel('Generation')
        axes[0, 1].set_ylabel('Count')
        axes[0, 1].grid(True, alpha=0.3)
        
        axes[0, 2].plot(steps, history['catalytic_fraction'], 'm-', linewidth=1.5)
        axes[0, 2].set_title('(c) Catalytic RNA Fraction')
        axes[0, 2].set_xlabel('Generation')
        axes[0, 2].set_ylabel('Fraction')
        axes[0, 2].grid(True, alpha=0.3)
        
        axes[1, 0].plot(steps, history['diversity'], 'c-', linewidth=1.5)
        axes[1, 0].set_title('(d) Sequence Diversity')
        axes[1, 0].set_xlabel('Generation')
        axes[1, 0].set_ylabel('Mean Hamming Distance')
        axes[1, 0].grid(True, alpha=0.3)
        
        axes[1, 1].plot(steps, history['gc_content'], color='orange', linewidth=1.5)
        axes[1, 1].set_title('(e) GC Content')
        axes[1, 1].set_xlabel('Generation')
        axes[1, 1].set_ylabel('GC Fraction')
        axes[1, 1].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
        axes[1, 1].grid(True, alpha=0.3)
        
        # Phase portrait: fitness vs diversity
        axes[1, 2].scatter(history['diversity'], history['mean_fitness'],
                          c=range(len(history['diversity'])), cmap='viridis', s=2, alpha=0.5)
        axes[1, 2].set_title('(f) Fitness–Diversity Phase Space')
        axes[1, 2].set_xlabel('Diversity')
        axes[1, 2].set_ylabel('Mean Fitness')
        axes[1, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'rna_world.png')
        fig.savefig(path)
        plt.close(fig)
        
        return {
            'final_mean_fitness': history['mean_fitness'][-1],
            'final_max_fitness': history['max_fitness'][-1],
            'final_population_size': history['population_size'][-1],
            'final_catalytic_fraction': history['catalytic_fraction'][-1],
            'final_diversity': history['diversity'][-1],
            'emergence_generation': next((i for i, f in enumerate(history['catalytic_fraction']) if f > 0.1), -1),
        }


# =============================================================================
# 3. METABOLISM-FIRST — Hydrothermal Vent Autocatalytic Cycles
# =============================================================================

class HydrothermalVentModel:
    """
    Simulates autocatalytic reaction cycles at hydrothermal vents:
    - Reverse citric acid cycle (rTCA) analogue
    - FeS-catalyzed reactions
    - Temperature and pH gradients
    - Energy coupling from proton/thermal gradients
    """
    
    SPECIES = [
        'CO2', 'H2', 'H2S', 'FeS', 'Fe2S2',
        'Acetate', 'Pyruvate', 'Oxaloacetate', 'Succinate',
        'Fumarate', 'Malate', 'Citrate',
        'Acetyl_CoA_analog', 'ATP_analog',
    ]
    
    def __init__(self, temp_vent=350, temp_ocean=4, pH_vent=9, pH_ocean=6):
        self.temp_vent = temp_vent
        self.temp_ocean = temp_ocean
        self.pH_vent = pH_vent
        self.pH_ocean = pH_ocean
        self.species_idx = {s: i for i, s in enumerate(self.SPECIES)}
        self.n_species = len(self.SPECIES)
        
    def temp_profile(self, x):
        """Temperature as function of distance from vent (0=vent, 1=ocean)."""
        return self.temp_vent * np.exp(-3 * x) + self.temp_ocean * (1 - np.exp(-3 * x))
    
    def pH_profile(self, x):
        return self.pH_vent * np.exp(-2 * x) + self.pH_ocean * (1 - np.exp(-2 * x))
    
    def arrhenius_factor(self, T, Ea=50.0):
        R = 8.314e-3  # kJ/(mol·K)
        T_ref = 300
        return np.exp(-Ea / R * (1/T - 1/T_ref))
    
    def ode_system(self, y, t, distance=0.3):
        dydt = np.zeros(self.n_species)
        T = self.temp_profile(distance) + 273.15
        pH = self.pH_profile(distance)
        af = self.arrhenius_factor(T)
        
        idx = self.species_idx
        
        # CO2 + H2 -> Formate -> Acetate (Wood-Ljungdahl analog)
        r1 = 0.01 * af * y[idx['CO2']] * y[idx['H2']]
        dydt[idx['CO2']] -= r1
        dydt[idx['H2']] -= r1
        dydt[idx['Acetate']] += r1
        
        # FeS catalysis: Acetate -> Pyruvate
        r2 = 0.005 * af * y[idx['Acetate']] * y[idx['FeS']] * y[idx['CO2']]
        dydt[idx['Acetate']] -= r2
        dydt[idx['CO2']] -= r2
        dydt[idx['Pyruvate']] += r2
        
        # rTCA cycle steps (simplified)
        r3 = 0.003 * af * y[idx['Pyruvate']] * y[idx['CO2']]
        dydt[idx['Pyruvate']] -= r3
        dydt[idx['Oxaloacetate']] += r3
        
        r4 = 0.004 * af * y[idx['Oxaloacetate']] * y[idx['H2']]
        dydt[idx['Oxaloacetate']] -= r4
        dydt[idx['Malate']] += r4
        
        r5 = 0.004 * af * y[idx['Malate']]
        dydt[idx['Malate']] -= r5
        dydt[idx['Fumarate']] += r5
        
        r6 = 0.003 * af * y[idx['Fumarate']] * y[idx['H2']]
        dydt[idx['Fumarate']] -= r6
        dydt[idx['Succinate']] += r6
        
        r7 = 0.002 * af * y[idx['Succinate']] * y[idx['CO2']]
        dydt[idx['Succinate']] -= r7
        dydt[idx['Citrate']] += r7
        
        # Citrate cleavage -> Acetyl-CoA analog + Oxaloacetate (autocatalytic closure)
        r8 = 0.003 * af * y[idx['Citrate']] * y[idx['FeS']]
        dydt[idx['Citrate']] -= r8
        dydt[idx['Acetyl_CoA_analog']] += r8
        dydt[idx['Oxaloacetate']] += r8  # autocatalytic feedback
        
        # FeS cluster formation
        r9 = 0.002 * y[idx['H2S']] * y[idx['FeS']]
        dydt[idx['FeS']] -= r9
        dydt[idx['Fe2S2']] += r9
        
        # Energy currency: proton gradient -> ATP analog
        proton_gradient = abs(self.pH_vent - self.pH_ocean)
        r10 = 0.001 * proton_gradient * y[idx['Fe2S2']]
        dydt[idx['ATP_analog']] += r10
        
        # Continuous supply from vent
        dydt[idx['CO2']] += 0.5
        dydt[idx['H2']] += 0.3
        dydt[idx['H2S']] += 0.1
        dydt[idx['FeS']] += 0.05
        
        return dydt
    
    def simulate(self, t_span=(0, 300), n_points=1500):
        t = np.linspace(t_span[0], t_span[1], n_points)
        y0 = np.zeros(self.n_species)
        y0[self.species_idx['CO2']] = 50.0
        y0[self.species_idx['H2']] = 30.0
        y0[self.species_idx['H2S']] = 10.0
        y0[self.species_idx['FeS']] = 5.0
        
        sol = odeint(self.ode_system, y0, t, mxstep=10000)
        sol = np.maximum(sol, 0)
        return t, sol
    
    def simulate_distance_sweep(self, distances=None, t_end=200):
        if distances is None:
            distances = np.linspace(0.05, 1.0, 10)
        results = {}
        for d in distances:
            t = np.linspace(0, t_end, 500)
            y0 = np.zeros(self.n_species)
            y0[self.species_idx['CO2']] = 50.0
            y0[self.species_idx['H2']] = 30.0
            y0[self.species_idx['H2S']] = 10.0
            y0[self.species_idx['FeS']] = 5.0
            sol = odeint(lambda y, t: self.ode_system(y, t, distance=d), y0, t, mxstep=10000)
            sol = np.maximum(sol, 0)
            results[d] = sol[-1, :]
        return distances, results
    
    def run_and_plot(self):
        t, sol = self.simulate()
        distances, dist_results = self.simulate_distance_sweep()
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Metabolism-First: Hydrothermal Vent Autocatalytic Network', fontsize=14, fontweight='bold')
        
        # rTCA cycle intermediates
        ax = axes[0, 0]
        for s in ['Pyruvate', 'Oxaloacetate', 'Malate', 'Fumarate', 'Succinate', 'Citrate']:
            ax.plot(t, sol[:, self.species_idx[s]], label=s, linewidth=1.5)
        ax.set_title('(a) rTCA Cycle Intermediates')
        ax.set_xlabel('Time')
        ax.set_ylabel('Concentration')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        
        # Energy carriers
        ax = axes[0, 1]
        ax.plot(t, sol[:, self.species_idx['Acetyl_CoA_analog']], label='Acetyl-CoA analog', linewidth=2)
        ax.plot(t, sol[:, self.species_idx['ATP_analog']], label='ATP analog', linewidth=2)
        ax.plot(t, sol[:, self.species_idx['Acetate']], label='Acetate', linewidth=1.5)
        ax.set_title('(b) Energy Carriers & Key Products')
        ax.set_xlabel('Time')
        ax.set_ylabel('Concentration')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Temperature & pH profiles
        ax = axes[1, 0]
        x = np.linspace(0, 1, 100)
        ax2 = ax.twinx()
        l1, = ax.plot(x, [self.temp_profile(xi) for xi in x], 'r-', linewidth=2, label='Temperature (°C)')
        l2, = ax2.plot(x, [self.pH_profile(xi) for xi in x], 'b-', linewidth=2, label='pH')
        ax.set_title('(c) Vent Environment Gradients')
        ax.set_xlabel('Distance from vent (normalized)')
        ax.set_ylabel('Temperature (°C)', color='r')
        ax2.set_ylabel('pH', color='b')
        ax.legend(handles=[l1, l2], loc='center right')
        ax.grid(True, alpha=0.3)
        
        # Distance sweep: key product yields
        ax = axes[1, 1]
        key_species = ['Acetate', 'Pyruvate', 'Citrate', 'ATP_analog']
        for s in key_species:
            vals = [dist_results[d][self.species_idx[s]] for d in distances]
            ax.plot(distances, vals, 'o-', label=s, linewidth=1.5)
        ax.set_title('(d) Product Yield vs. Distance from Vent')
        ax.set_xlabel('Distance (normalized)')
        ax.set_ylabel('Final Concentration')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'hydrothermal_vent.png')
        fig.savefig(path)
        plt.close(fig)
        
        final = {s: float(sol[-1, self.species_idx[s]]) for s in self.SPECIES}
        autocatalytic_ratio = final.get('Oxaloacetate', 0) / max(final.get('Pyruvate', 1), 0.01)
        
        return {
            'final_concentrations': final,
            'autocatalytic_ratio': autocatalytic_ratio,
            'optimal_distance': float(distances[np.argmax([dist_results[d][self.species_idx['Citrate']] for d in distances])]),
        }


# =============================================================================
# 4. STOCHASTIC CHEMICAL KINETICS — Gillespie SSA / CME
# =============================================================================

class GillespieSSA:
    """
    Gillespie Stochastic Simulation Algorithm for prebiotic reaction systems.
    Computes exact stochastic trajectories and estimates biopolymer emergence probability.
    """
    
    def __init__(self):
        self.species = ['Monomer_A', 'Monomer_B', 'Dimer_AB', 'Trimer', 'Tetramer',
                       'Pentamer', 'Hexamer', 'Polymer_7', 'Polymer_8']
        self.n_species = len(self.species)
        self.species_idx = {s: i for i, s in enumerate(self.species)}
    
    def propensities(self, state):
        """Calculate reaction propensities for polymerization system."""
        a = []
        s = self.species_idx
        
        # Monomer_A + Monomer_B -> Dimer_AB
        a.append(0.01 * state[s['Monomer_A']] * state[s['Monomer_B']])
        # Dimer_AB + Monomer_A -> Trimer
        a.append(0.008 * state[s['Dimer_AB']] * state[s['Monomer_A']])
        # Dimer_AB + Monomer_B -> Trimer
        a.append(0.008 * state[s['Dimer_AB']] * state[s['Monomer_B']])
        # Trimer + Monomer -> Tetramer
        a.append(0.005 * state[s['Trimer']] * (state[s['Monomer_A']] + state[s['Monomer_B']]))
        # Tetramer + Monomer -> Pentamer
        a.append(0.003 * state[s['Tetramer']] * (state[s['Monomer_A']] + state[s['Monomer_B']]))
        # Pentamer + Monomer -> Hexamer
        a.append(0.002 * state[s['Pentamer']] * (state[s['Monomer_A']] + state[s['Monomer_B']]))
        # Hexamer + Monomer -> Polymer_7
        a.append(0.001 * state[s['Hexamer']] * (state[s['Monomer_A']] + state[s['Monomer_B']]))
        # Polymer_7 + Monomer -> Polymer_8
        a.append(0.0005 * state[s['Polymer_7']] * (state[s['Monomer_A']] + state[s['Monomer_B']]))
        
        # Hydrolysis reactions (reverse)
        a.append(0.002 * state[s['Dimer_AB']])
        a.append(0.003 * state[s['Trimer']])
        a.append(0.004 * state[s['Tetramer']])
        a.append(0.005 * state[s['Pentamer']])
        
        return np.array(a)
    
    def stoichiometry(self):
        """Stoichiometry matrix for each reaction."""
        n_rxn = 12
        S = np.zeros((n_rxn, self.n_species), dtype=int)
        s = self.species_idx
        
        # Forward polymerization
        S[0, s['Monomer_A']] = -1; S[0, s['Monomer_B']] = -1; S[0, s['Dimer_AB']] = 1
        S[1, s['Dimer_AB']] = -1; S[1, s['Monomer_A']] = -1; S[1, s['Trimer']] = 1
        S[2, s['Dimer_AB']] = -1; S[2, s['Monomer_B']] = -1; S[2, s['Trimer']] = 1
        S[3, s['Trimer']] = -1; S[3, s['Monomer_A']] = -1; S[3, s['Tetramer']] = 1
        S[4, s['Tetramer']] = -1; S[4, s['Monomer_A']] = -1; S[4, s['Pentamer']] = 1
        S[5, s['Pentamer']] = -1; S[5, s['Monomer_A']] = -1; S[5, s['Hexamer']] = 1
        S[6, s['Hexamer']] = -1; S[6, s['Monomer_A']] = -1; S[6, s['Polymer_7']] = 1
        S[7, s['Polymer_7']] = -1; S[7, s['Monomer_A']] = -1; S[7, s['Polymer_8']] = 1
        
        # Hydrolysis
        S[8, s['Dimer_AB']] = -1; S[8, s['Monomer_A']] = 1; S[8, s['Monomer_B']] = 1
        S[9, s['Trimer']] = -1; S[9, s['Dimer_AB']] = 1; S[9, s['Monomer_A']] = 1
        S[10, s['Tetramer']] = -1; S[10, s['Trimer']] = 1; S[10, s['Monomer_A']] = 1
        S[11, s['Pentamer']] = -1; S[11, s['Tetramer']] = 1; S[11, s['Monomer_A']] = 1
        
        return S
    
    def run_single(self, initial_state=None, t_max=500):
        if initial_state is None:
            initial_state = np.zeros(self.n_species, dtype=int)
            initial_state[self.species_idx['Monomer_A']] = 500
            initial_state[self.species_idx['Monomer_B']] = 500
        
        S = self.stoichiometry()
        state = initial_state.copy()
        t = 0
        times = [0]
        states = [state.copy()]
        
        while t < t_max:
            props = self.propensities(state)
            a0 = np.sum(props)
            if a0 == 0:
                break
            
            tau = np.random.exponential(1.0 / a0)
            t += tau
            
            j = np.searchsorted(np.cumsum(props), np.random.random() * a0)
            j = min(j, len(S) - 1)
            
            new_state = state + S[j]
            if np.all(new_state >= 0):
                state = new_state
            
            times.append(t)
            states.append(state.copy())
        
        return np.array(times), np.array(states)
    
    def run_ensemble(self, n_runs=100, t_max=500):
        max_polymer_lengths = []
        final_dimer_counts = []
        final_trimer_counts = []
        emergence_times = []
        
        for _ in range(n_runs):
            times, states = self.run_single(t_max=t_max)
            
            max_len = 2
            for sp_name in ['Trimer', 'Tetramer', 'Pentamer', 'Hexamer', 'Polymer_7', 'Polymer_8']:
                if states[-1, self.species_idx[sp_name]] > 0:
                    max_len = max(max_len, self.species.index(sp_name) + 1)
            
            max_polymer_lengths.append(max_len)
            final_dimer_counts.append(states[-1, self.species_idx['Dimer_AB']])
            final_trimer_counts.append(states[-1, self.species_idx['Trimer']])
            
            # First time a tetramer+ appears
            tetramer_col = states[:, self.species_idx['Tetramer']]
            tetramer_times = np.where(tetramer_col > 0)[0]
            if len(tetramer_times) > 0:
                emergence_times.append(times[tetramer_times[0]])
            else:
                emergence_times.append(np.nan)
        
        return {
            'max_polymer_lengths': max_polymer_lengths,
            'final_dimer_counts': final_dimer_counts,
            'final_trimer_counts': final_trimer_counts,
            'emergence_times': emergence_times,
        }
    
    def run_and_plot(self):
        # Single trajectory
        times, states = self.run_single(t_max=500)
        
        # Ensemble
        ensemble = self.run_ensemble(n_runs=200, t_max=500)
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Stochastic Chemical Kinetics: Biopolymer Emergence (Gillespie SSA)',
                     fontsize=14, fontweight='bold')
        
        # Single trajectory
        ax = axes[0, 0]
        for sp in ['Monomer_A', 'Dimer_AB', 'Trimer', 'Tetramer', 'Pentamer']:
            ax.plot(times, states[:, self.species_idx[sp]], label=sp, linewidth=1)
        ax.set_title('(a) Single Stochastic Trajectory')
        ax.set_xlabel('Time')
        ax.set_ylabel('Copy Number')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        
        # Polymer length distribution
        ax = axes[0, 1]
        lengths = ensemble['max_polymer_lengths']
        bins = np.arange(1.5, 10.5, 1)
        ax.hist(lengths, bins=bins, edgecolor='black', color='steelblue', alpha=0.7, density=True)
        ax.set_title('(b) Max Polymer Length Distribution (N=200)')
        ax.set_xlabel('Maximum Polymer Length')
        ax.set_ylabel('Probability')
        ax.grid(True, alpha=0.3)
        
        # Emergence time distribution
        ax = axes[1, 0]
        valid_times = [t for t in ensemble['emergence_times'] if not np.isnan(t)]
        if valid_times:
            ax.hist(valid_times, bins=30, edgecolor='black', color='coral', alpha=0.7)
        ax.set_title(f'(c) Tetramer+ Emergence Time (appeared in {len(valid_times)}/200 runs)')
        ax.set_xlabel('First Appearance Time')
        ax.set_ylabel('Count')
        ax.grid(True, alpha=0.3)
        
        # Probability of polymer emergence vs chain length
        ax = axes[1, 1]
        chain_lengths = list(range(2, 9))
        probs = []
        for cl in chain_lengths:
            sp_name = self.species[cl]
            count = sum(1 for l in lengths if l >= cl + 1)
            probs.append(count / len(lengths))
        ax.bar(chain_lengths, probs, color='teal', alpha=0.7, edgecolor='black')
        ax.set_title('(d) Emergence Probability by Chain Length')
        ax.set_xlabel('Minimum Chain Length')
        ax.set_ylabel('P(emergence)')
        ax.set_xticks(chain_lengths)
        ax.set_xticklabels([self.species[i] for i in chain_lengths], rotation=45, fontsize=7)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'gillespie_ssa.png')
        fig.savefig(path)
        plt.close(fig)
        
        prob_tetramer = len(valid_times) / 200
        mean_emergence = np.mean(valid_times) if valid_times else float('nan')
        
        return {
            'prob_tetramer_emergence': prob_tetramer,
            'mean_emergence_time': float(mean_emergence) if not np.isnan(mean_emergence) else None,
            'mean_max_polymer_length': float(np.mean(lengths)),
            'polymer_length_distribution': {str(k): int(v) for k, v in
                                            zip(*np.unique(lengths, return_counts=True))},
        }


# =============================================================================
# 5. MEMBRANE SELF-ASSEMBLY — Protocell Formation
# =============================================================================

class ProtocellFormation:
    """
    Models the self-assembly of amphiphilic molecules into vesicles/protocells.
    Uses a coarse-grained lattice-based approach with:
    - Amphiphile aggregation dynamics
    - Critical micelle concentration (CMC)
    - Vesicle growth and division
    - Encapsulation of catalytic polymers
    """
    
    def __init__(self, grid_size=50, n_amphiphiles=800, n_polymers=50):
        self.grid_size = grid_size
        self.n_amphiphiles = n_amphiphiles
        self.n_polymers = n_polymers
        self.cmc = 0.05  # Critical Micelle Concentration (fraction)
    
    def initialize(self):
        """Place amphiphiles and polymers randomly on grid."""
        positions_a = np.random.randint(0, self.grid_size, (self.n_amphiphiles, 2))
        positions_p = np.random.randint(0, self.grid_size, (self.n_polymers, 2))
        return positions_a, positions_p
    
    def compute_local_density(self, positions, grid_size):
        """Compute local density on grid."""
        density = np.zeros((grid_size, grid_size))
        for p in positions:
            density[p[0], p[1]] += 1
        return density / (grid_size * grid_size)
    
    def identify_clusters(self, positions, threshold=2.5):
        """Identify clusters of amphiphiles using distance-based grouping."""
        from scipy.spatial.distance import pdist, squareform
        if len(positions) < 2:
            return [list(range(len(positions)))]
        
        dist_matrix = squareform(pdist(positions))
        visited = set()
        clusters = []
        
        for i in range(len(positions)):
            if i in visited:
                continue
            cluster = []
            queue = [i]
            while queue:
                node = queue.pop(0)
                if node in visited:
                    continue
                visited.add(node)
                cluster.append(node)
                neighbors = np.where(dist_matrix[node] < threshold)[0]
                for n in neighbors:
                    if n not in visited:
                        queue.append(n)
            if len(cluster) > 0:
                clusters.append(cluster)
        
        return clusters
    
    def simulate(self, n_steps=300):
        positions_a, positions_p = self.initialize()
        
        history = {
            'n_clusters': [],
            'largest_cluster': [],
            'mean_cluster_size': [],
            'vesicle_count': [],
            'encapsulated_polymers': [],
        }
        
        density_snapshots = []
        
        for step in range(n_steps):
            density = self.compute_local_density(positions_a, self.grid_size)
            
            # Move amphiphiles: biased toward high-density regions (aggregation)
            new_positions = []
            for i, pos in enumerate(positions_a):
                local_d = density[pos[0], pos[1]]
                if local_d > self.cmc:
                    # Stay near cluster (small random walk)
                    dx = np.random.randint(-1, 2, 2)
                else:
                    # Random walk with drift toward high density
                    dx = np.random.randint(-2, 3, 2)
                    # Check neighbors for gradient
                    for direction in [np.array([1,0]), np.array([-1,0]),
                                     np.array([0,1]), np.array([0,-1])]:
                        np_ = (pos + direction) % self.grid_size
                        if density[np_[0], np_[1]] > local_d:
                            dx += direction
                            break
                
                new_pos = (pos + dx) % self.grid_size
                new_positions.append(new_pos)
            positions_a = np.array(new_positions)
            
            # Move polymers toward nearby clusters
            new_positions_p = []
            for pos in positions_p:
                local_d = density[pos[0], pos[1]]
                if local_d > self.cmc * 2:
                    dx = np.random.randint(-1, 2, 2)
                else:
                    dx = np.random.randint(-2, 3, 2)
                new_pos = (pos + dx) % self.grid_size
                new_positions_p.append(new_pos)
            positions_p = np.array(new_positions_p)
            
            # Identify clusters
            clusters = self.identify_clusters(positions_a, threshold=3.0)
            cluster_sizes = [len(c) for c in clusters]
            
            vesicles = [c for c in clusters if len(c) >= 15]  # Vesicle threshold
            
            # Count encapsulated polymers
            encapsulated = 0
            for v_cluster in vesicles:
                v_positions = positions_a[v_cluster]
                center = np.mean(v_positions, axis=0)
                radius = np.max(np.sqrt(np.sum((v_positions - center)**2, axis=1))) + 1
                for pp in positions_p:
                    if np.sqrt(np.sum((pp - center)**2)) < radius:
                        encapsulated += 1
            
            history['n_clusters'].append(len(clusters))
            history['largest_cluster'].append(max(cluster_sizes) if cluster_sizes else 0)
            history['mean_cluster_size'].append(np.mean(cluster_sizes) if cluster_sizes else 0)
            history['vesicle_count'].append(len(vesicles))
            history['encapsulated_polymers'].append(encapsulated)
            
            if step in [0, n_steps//4, n_steps//2, n_steps-1]:
                density_snapshots.append((step, density.copy(), positions_a.copy(), positions_p.copy()))
        
        return history, density_snapshots
    
    def run_and_plot(self):
        history, snapshots = self.simulate()
        
        fig = plt.figure(figsize=(16, 12))
        gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)
        fig.suptitle('Protocell Formation: Membrane Self-Assembly Simulation',
                     fontsize=14, fontweight='bold')
        
        # Density snapshots
        for i, (step, density, pos_a, pos_p) in enumerate(snapshots):
            ax = fig.add_subplot(gs[0, i])
            ax.imshow(density.T, origin='lower', cmap='YlOrRd', vmin=0, vmax=0.02)
            ax.scatter(pos_p[:, 0], pos_p[:, 1], c='blue', s=5, alpha=0.5, label='Polymers')
            ax.set_title(f'Step {step}', fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_ylabel('Amphiphile Density')
        
        # Time series
        steps = range(len(history['n_clusters']))
        
        ax1 = fig.add_subplot(gs[1, 0:2])
        ax1.plot(steps, history['n_clusters'], 'b-', linewidth=1.5, label='Total clusters')
        ax1.plot(steps, history['vesicle_count'], 'r-', linewidth=2, label='Vesicles (≥15)')
        ax1.set_title('(e) Cluster & Vesicle Count')
        ax1.set_xlabel('Time Step')
        ax1.set_ylabel('Count')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        ax2 = fig.add_subplot(gs[1, 2:4])
        ax2.plot(steps, history['largest_cluster'], 'g-', linewidth=1.5, label='Largest')
        ax2.plot(steps, history['mean_cluster_size'], 'm-', linewidth=1.5, label='Mean')
        ax2.set_title('(f) Cluster Size Dynamics')
        ax2.set_xlabel('Time Step')
        ax2.set_ylabel('Size (# amphiphiles)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        ax3 = fig.add_subplot(gs[2, 0:2])
        ax3.plot(steps, history['encapsulated_polymers'], 'k-', linewidth=2)
        ax3.set_title('(g) Encapsulated Polymers in Vesicles')
        ax3.set_xlabel('Time Step')
        ax3.set_ylabel('Count')
        ax3.grid(True, alpha=0.3)
        
        # Final cluster size distribution
        ax4 = fig.add_subplot(gs[2, 2:4])
        final_density = snapshots[-1][1]
        ax4.hist(final_density.flatten(), bins=50, edgecolor='black', color='steelblue', alpha=0.7)
        ax4.axvline(x=self.cmc, color='red', linestyle='--', label=f'CMC={self.cmc}')
        ax4.set_title('(h) Final Density Distribution')
        ax4.set_xlabel('Local Amphiphile Density')
        ax4.set_ylabel('Count')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        path = os.path.join(FIGURES_DIR, 'protocell_formation.png')
        fig.savefig(path)
        plt.close(fig)
        
        return {
            'final_vesicle_count': history['vesicle_count'][-1],
            'final_largest_cluster': history['largest_cluster'][-1],
            'final_encapsulated': history['encapsulated_polymers'][-1],
            'mean_cluster_size_final': history['mean_cluster_size'][-1],
        }


# =============================================================================
# 6. ENCELADUS / TITAN ENVIRONMENTAL CONDITIONS
# =============================================================================

class ExoplanetaryChemistry:
    """
    Adapts prebiotic chemistry simulations to Enceladus and Titan conditions.
    - Enceladus: subsurface ocean, hydrothermal activity, ~0°C water
    - Titan: cryogenic (-179°C), liquid methane/ethane, tholin chemistry
    """
    
    def simulate_enceladus(self, t_max=300, n_points=1000):
        """
        Enceladus subsurface ocean chemistry:
        H2 + CO2 -> organics in alkaline hydrothermal conditions
        """
        species = ['H2', 'CO2', 'CH4_aq', 'Formate', 'Acetate', 'Amino_acids',
                   'NH3', 'H2S', 'Methanol', 'Formaldehyde']
        n = len(species)
        idx = {s: i for i, s in enumerate(species)}
        
        def ode(y, t):
            dydt = np.zeros(n)
            T = 273 + 50  # ~50°C at vent
            af = np.exp(-40 / (8.314e-3 * T))
            
            r1 = 0.005 * af * y[idx['CO2']] * y[idx['H2']]
            dydt[idx['CO2']] -= r1; dydt[idx['H2']] -= 4*r1; dydt[idx['CH4_aq']] += r1
            
            r2 = 0.008 * af * y[idx['CO2']] * y[idx['H2']]
            dydt[idx['CO2']] -= r2; dydt[idx['H2']] -= r2; dydt[idx['Formate']] += r2
            
            r3 = 0.003 * af * y[idx['Formate']] * y[idx['H2']]
            dydt[idx['Formate']] -= r3; dydt[idx['Methanol']] += r3
            
            r4 = 0.002 * af * y[idx['Formate']] * y[idx['Formate']]
            dydt[idx['Formate']] -= 2*r4; dydt[idx['Acetate']] += r4
            
            r5 = 0.001 * af * y[idx['Formaldehyde']] * y[idx['NH3']] * y[idx['Formate']]
            dydt[idx['Amino_acids']] += r5
            dydt[idx['Formaldehyde']] -= r5; dydt[idx['NH3']] -= r5; dydt[idx['Formate']] -= r5
            
            r6 = 0.004 * af * y[idx['CO2']] * y[idx['H2']]
            dydt[idx['Formaldehyde']] += r6; dydt[idx['CO2']] -= r6; dydt[idx['H2']] -= r6
            
            # Continuous supply
            dydt[idx['H2']] += 0.3
            dydt[idx['CO2']] += 0.2
            dydt[idx['NH3']] += 0.05
            dydt[idx['H2S']] += 0.02
            
            return dydt
        
        t = np.linspace(0, t_max, n_points)
        y0 = np.zeros(n)
        y0[idx['H2']] = 30; y0[idx['CO2']] = 40; y0[idx['NH3']] = 5; y0[idx['H2S']] = 3
        
        sol = odeint(ode, y0, t, mxstep=10000)
        sol = np.maximum(sol, 0)
        
        return t, sol, species, idx
    
    def simulate_titan(self, t_max=500, n_points=1000):
        """
        Titan surface/atmosphere chemistry:
        N2 + CH4 + UV -> HCN + tholins (complex organics)
        Liquid methane/ethane solvent chemistry
        """
        species = ['N2', 'CH4', 'C2H6', 'HCN', 'C2H2', 'C2H4',
                   'Tholins', 'HC3N', 'CH2NH', 'Adenine_precursor']
        n = len(species)
        idx = {s: i for i, s in enumerate(species)}
        
        def ode(y, t):
            dydt = np.zeros(n)
            T = 94  # Titan surface temperature in K
            uv = 0.01  # UV flux factor (attenuated by atmosphere)
            
            # N2 + CH4 -> HCN + H2 (photolysis)
            r1 = 0.01 * uv * y[idx['N2']] * y[idx['CH4']]
            dydt[idx['N2']] -= r1; dydt[idx['CH4']] -= r1; dydt[idx['HCN']] += r1
            
            # CH4 -> C2H6 (ethane production)
            r2 = 0.005 * uv * y[idx['CH4']] ** 2
            dydt[idx['CH4']] -= 2*r2; dydt[idx['C2H6']] += r2
            
            # CH4 -> C2H2 + C2H4
            r3 = 0.003 * uv * y[idx['CH4']]
            dydt[idx['CH4']] -= r3; dydt[idx['C2H2']] += 0.5*r3; dydt[idx['C2H4']] += 0.5*r3
            
            # HCN + C2H2 -> HC3N
            r4 = 0.002 * y[idx['HCN']] * y[idx['C2H2']]
            dydt[idx['HCN']] -= r4; dydt[idx['C2H2']] -= r4; dydt[idx['HC3N']] += r4
            
            # HCN polymerization -> Tholins
            r5 = 0.0005 * y[idx['HCN']] ** 2 * uv
            dydt[idx['HCN']] -= 2*r5; dydt[idx['Tholins']] += r5
            
            # CH4 + N2 -> CH2NH (methylenimine)
            r6 = 0.001 * uv * y[idx['CH4']] * y[idx['N2']]
            dydt[idx['CH4']] -= r6; dydt[idx['N2']] -= r6; dydt[idx['CH2NH']] += r6
            
            # HCN pentamerization -> Adenine precursor (very slow at 94K)
            cryo_factor = np.exp(-20 / (8.314e-3 * T))  # very slow at 94K
            r7 = 0.0001 * cryo_factor * y[idx['HCN']] ** 3
            dydt[idx['HCN']] -= 3*r7; dydt[idx['Adenine_precursor']] += r7
            
            # Atmospheric replenishment
            dydt[idx['N2']] += 0.1
            dydt[idx['CH4']] += 0.2
            
            return dydt
        
        t = np.linspace(0, t_max, n_points)
        y0 = np.zeros(n)
        y0[idx['N2']] = 100; y0[idx['CH4']] = 50
        
        sol = odeint(ode, y0, t, mxstep=10000)
        sol = np.maximum(sol, 0)
        
        return t, sol, species, idx
    
    def compare_environments(self):
        """Compare chemical evolution potential across Earth, Enceladus, and Titan."""
        environments = {
            'Early Earth': {
                'temperature': 80, 'pH': 7.0, 'energy_sources': 4,
                'liquid_water': True, 'organic_precursors': 3,
                'redox_gradient': True, 'mineral_catalysts': True,
                'score': 0
            },
            'Enceladus': {
                'temperature': 50, 'pH': 9.0, 'energy_sources': 2,
                'liquid_water': True, 'organic_precursors': 2,
                'redox_gradient': True, 'mineral_catalysts': True,
                'score': 0
            },
            'Titan': {
                'temperature': -179, 'pH': None, 'energy_sources': 1,
                'liquid_water': False, 'organic_precursors': 5,
                'redox_gradient': False, 'mineral_catalysts': False,
                'score': 0
            },
        }
        
        for env in environments.values():
            s = 0
            s += 3 if env['liquid_water'] else 0
            s += min(env['energy_sources'], 4)
            s += min(env['organic_precursors'], 5)
            s += 2 if env['redox_gradient'] else 0
            s += 2 if env['mineral_catalysts'] else 0
            if env['temperature'] is not None:
                if 0 <= env['temperature'] <= 120:
                    s += 3
                elif -50 <= env['temperature'] < 0 or 120 < env['temperature'] <= 200:
                    s += 1
            env['score'] = s
        
        return environments
    
    def run_and_plot(self):
        t_e, sol_e, sp_e, idx_e = self.simulate_enceladus()
        t_t, sol_t, sp_t, idx_t = self.simulate_titan()
        env_comparison = self.compare_environments()
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        fig.suptitle('Extraterrestrial Chemical Evolution: Enceladus & Titan',
                     fontsize=14, fontweight='bold')
        
        # Enceladus: reactants
        ax = axes[0, 0]
        for s in ['H2', 'CO2', 'NH3']:
            ax.plot(t_e, sol_e[:, idx_e[s]], label=s, linewidth=1.5)
        ax.set_title('(a) Enceladus: Reactants')
        ax.set_xlabel('Time'); ax.set_ylabel('Concentration')
        ax.legend(); ax.grid(True, alpha=0.3)
        
        # Enceladus: products
        ax = axes[0, 1]
        for s in ['Formate', 'Acetate', 'Methanol', 'Amino_acids']:
            ax.plot(t_e, sol_e[:, idx_e[s]], label=s, linewidth=1.5)
        ax.set_title('(b) Enceladus: Organic Products')
        ax.set_xlabel('Time'); ax.set_ylabel('Concentration')
        ax.legend(); ax.grid(True, alpha=0.3)
        
        # Titan: atmosphere
        ax = axes[0, 2]
        for s in ['N2', 'CH4', 'C2H6', 'C2H2']:
            ax.plot(t_t, sol_t[:, idx_t[s]], label=s, linewidth=1.5)
        ax.set_title('(c) Titan: Atmospheric Species')
        ax.set_xlabel('Time'); ax.set_ylabel('Concentration')
        ax.legend(); ax.grid(True, alpha=0.3)
        
        # Titan: complex organics
        ax = axes[1, 0]
        for s in ['HCN', 'Tholins', 'HC3N', 'CH2NH', 'Adenine_precursor']:
            ax.plot(t_t, sol_t[:, idx_t[s]], label=s, linewidth=1.5)
        ax.set_title('(d) Titan: Complex Organics')
        ax.set_xlabel('Time'); ax.set_ylabel('Concentration')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
        
        # Environment comparison
        ax = axes[1, 1]
        envs = list(env_comparison.keys())
        scores = [env_comparison[e]['score'] for e in envs]
        colors = ['#2ecc71', '#3498db', '#e67e22']
        bars = ax.bar(envs, scores, color=colors, edgecolor='black', alpha=0.8)
        ax.set_title('(e) Chemical Evolution Potential Score')
        ax.set_ylabel('Score')
        for bar, score in zip(bars, scores):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                   str(score), ha='center', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Comparative organic yield (normalized)
        ax = axes[1, 2]
        categories = ['Amino Acids\n/ Tholins', 'Simple Organics\n(Formate/HCN)', 'Energy Carriers\n(Acetate/HC3N)']
        earth_vals = [0.8, 0.9, 0.7]
        enceladus_vals = [
            float(sol_e[-1, idx_e['Amino_acids']]) / 50,
            float(sol_e[-1, idx_e['Formate']]) / 50,
            float(sol_e[-1, idx_e['Acetate']]) / 50,
        ]
        titan_vals = [
            float(sol_t[-1, idx_t['Tholins']]) / 50,
            float(sol_t[-1, idx_t['HCN']]) / 50,
            float(sol_t[-1, idx_t['HC3N']]) / 50,
        ]
        
        x = np.arange(len(categories))
        w = 0.25
        ax.bar(x - w, earth_vals, w, label='Early Earth', color='#2ecc71', alpha=0.8)
        ax.bar(x, enceladus_vals, w, label='Enceladus', color='#3498db', alpha=0.8)
        ax.bar(x + w, titan_vals, w, label='Titan', color='#e67e22', alpha=0.8)
        ax.set_title('(f) Comparative Organic Yield')
        ax.set_ylabel('Normalized Yield')
        ax.set_xticks(x)
        ax.set_xticklabels(categories, fontsize=8)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'exoplanetary_chemistry.png')
        fig.savefig(path)
        plt.close(fig)
        
        return {
            'enceladus_final': {s: float(sol_e[-1, idx_e[s]]) for s in sp_e},
            'titan_final': {s: float(sol_t[-1, idx_t[s]]) for s in sp_t},
            'environment_scores': {e: env_comparison[e]['score'] for e in envs},
        }


# =============================================================================
# INTEGRATED ANALYSIS
# =============================================================================

def run_integrated_analysis(results):
    """Cross-module comparative analysis and summary figure."""
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle('Integrated Analysis: Chemical Evolution Pathways', fontsize=14, fontweight='bold')
    
    # 1. Comparison of organic yields across models
    ax = axes[0]
    models = ['Miller-Urey', 'Hydrothermal\nVent', 'Enceladus', 'Titan']
    
    mu_aa = sum(results['miller_urey']['final_concentrations'].get(aa, 0)
                for aa in ['Glycine', 'Alanine', 'Aspartate', 'Valine'])
    hv_org = sum(results['hydrothermal']['final_concentrations'].get(s, 0)
                 for s in ['Pyruvate', 'Citrate', 'Acetate'])
    en_org = sum(results['exoplanetary']['enceladus_final'].get(s, 0)
                 for s in ['Amino_acids', 'Formate', 'Acetate'])
    ti_org = sum(results['exoplanetary']['titan_final'].get(s, 0)
                 for s in ['Tholins', 'HCN', 'HC3N'])
    
    yields = [mu_aa, hv_org, en_org, ti_org]
    colors = ['#e74c3c', '#2ecc71', '#3498db', '#e67e22']
    ax.bar(models, yields, color=colors, edgecolor='black', alpha=0.8)
    ax.set_title('(a) Total Organic Yield by Model')
    ax.set_ylabel('Total Concentration')
    ax.grid(True, alpha=0.3, axis='y')
    
    # 2. Stochastic vs deterministic comparison
    ax = axes[1]
    polymer_dist = results['gillespie']['polymer_length_distribution']
    lengths = sorted([int(k) for k in polymer_dist.keys()])
    counts = [polymer_dist[str(l)] for l in lengths]
    ax.bar(lengths, counts, color='teal', edgecolor='black', alpha=0.7)
    ax.set_title('(b) Polymer Length Distribution\n(Gillespie SSA, N=200)')
    ax.set_xlabel('Max Polymer Length')
    ax.set_ylabel('Count')
    ax.grid(True, alpha=0.3)
    
    # 3. Summary radar chart (as bar chart alternative)
    ax = axes[2]
    metrics = ['Amino Acid\nYield', 'Nucleobase\nYield', 'Autocatalytic\nStrength',
               'Polymer\nEmergence', 'Protocell\nFormation']
    
    mu_nb = sum(results['miller_urey']['final_concentrations'].get(nb, 0)
                for nb in ['Adenine', 'Guanine', 'Cytosine', 'Uracil'])
    
    values = [
        min(mu_aa / 100, 1.0),
        min(mu_nb / 10, 1.0),
        min(results['hydrothermal']['autocatalytic_ratio'] / 5, 1.0),
        results['gillespie']['prob_tetramer_emergence'],
        min(results['protocell']['final_vesicle_count'] / 10, 1.0),
    ]
    
    ax.barh(metrics, values, color=['#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#3498db'],
            edgecolor='black', alpha=0.8)
    ax.set_title('(c) Chemical Evolution Metrics\n(Normalized)')
    ax.set_xlabel('Score (0-1)')
    ax.set_xlim(0, 1.1)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, 'integrated_analysis.png')
    fig.savefig(path)
    plt.close(fig)
    
    return {
        'total_organic_yields': dict(zip(models, [float(v) for v in yields])),
        'evolution_metrics': dict(zip([m.replace('\n', ' ') for m in metrics],
                                      [float(v) for v in values])),
    }


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    print("=" * 70)
    print("PREBIOTIC CHEMISTRY SIMULATION FRAMEWORK")
    print("Integrated Stochastic-Deterministic Chemical Evolution Models")
    print("=" * 70)
    
    np.random.seed(42)
    
    # 1. Miller-Urey Extended
    print("\n[1/6] Running Extended Miller-Urey simulation...")
    mu = MillerUreyExtended()
    RESULTS['miller_urey'] = mu.run_and_plot()
    print(f"  → Generated {RESULTS['miller_urey']['network_nodes']} species, "
          f"{RESULTS['miller_urey']['network_edges']} reactions")
    print(f"  → Glycine yield: {RESULTS['miller_urey']['final_concentrations']['Glycine']:.4f}")
    print(f"  → Adenine yield: {RESULTS['miller_urey']['final_concentrations']['Adenine']:.4f}")
    
    # 2. RNA World
    print("\n[2/6] Running RNA World self-replication simulation...")
    rna = RNAWorldSimulation(n_sequences=200, seq_length=20, n_steps=5000)
    RESULTS['rna_world'] = rna.run_and_plot()
    print(f"  → Final mean fitness: {RESULTS['rna_world']['final_mean_fitness']:.4f}")
    print(f"  → Final catalytic fraction: {RESULTS['rna_world']['final_catalytic_fraction']:.4f}")
    print(f"  → Emergence generation: {RESULTS['rna_world']['emergence_generation']}")
    
    # 3. Hydrothermal Vent
    print("\n[3/6] Running Hydrothermal Vent metabolism-first simulation...")
    hv = HydrothermalVentModel()
    RESULTS['hydrothermal'] = hv.run_and_plot()
    print(f"  → Autocatalytic ratio: {RESULTS['hydrothermal']['autocatalytic_ratio']:.4f}")
    print(f"  → Optimal distance: {RESULTS['hydrothermal']['optimal_distance']:.2f}")
    
    # 4. Gillespie SSA
    print("\n[4/6] Running Gillespie SSA stochastic polymerization...")
    ssa = GillespieSSA()
    RESULTS['gillespie'] = ssa.run_and_plot()
    print(f"  → Tetramer emergence probability: {RESULTS['gillespie']['prob_tetramer_emergence']:.4f}")
    print(f"  → Mean max polymer length: {RESULTS['gillespie']['mean_max_polymer_length']:.2f}")
    
    # 5. Protocell Formation
    print("\n[5/6] Running Protocell formation simulation...")
    pc = ProtocellFormation(grid_size=50, n_amphiphiles=800, n_polymers=50)
    RESULTS['protocell'] = pc.run_and_plot()
    print(f"  → Final vesicle count: {RESULTS['protocell']['final_vesicle_count']}")
    print(f"  → Encapsulated polymers: {RESULTS['protocell']['final_encapsulated']}")
    
    # 6. Exoplanetary Chemistry
    print("\n[6/6] Running Enceladus/Titan chemistry simulation...")
    exo = ExoplanetaryChemistry()
    RESULTS['exoplanetary'] = exo.run_and_plot()
    for env, score in RESULTS['exoplanetary']['environment_scores'].items():
        print(f"  → {env} score: {score}")
    
    # Integrated Analysis
    print("\n[+] Running integrated cross-model analysis...")
    RESULTS['integrated'] = run_integrated_analysis(RESULTS)
    
    # Save results
    results_serializable = {}
    for k, v in RESULTS.items():
        results_serializable[k] = {}
        for k2, v2 in v.items():
            if isinstance(v2, (dict, list, int, float, str, bool, type(None))):
                results_serializable[k][k2] = v2
            elif isinstance(v2, np.floating):
                results_serializable[k][k2] = float(v2)
            elif isinstance(v2, np.integer):
                results_serializable[k][k2] = int(v2)
            else:
                results_serializable[k][k2] = str(v2)
    
    with open('simulation_results.json', 'w') as f:
        json.dump(results_serializable, f, indent=2, default=str)
    
    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE")
    print(f"Results saved to simulation_results.json")
    print(f"Figures saved to {FIGURES_DIR}/")
    print("=" * 70)
    
    return RESULTS


if __name__ == '__main__':
    results = main()
