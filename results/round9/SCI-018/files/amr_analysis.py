"""
AMR Evolution Prediction Computational Framework
=================================================
A comprehensive computational framework integrating:
1. ARG detection from whole-genome sequencing data (simulated)
2. Fitness landscape construction
3. Evolutionary path prediction (accessible mutational paths)
4. HGT network modeling
5. Spatiotemporal antibiotic resistance dynamics (SIR-based)
6. Antibiotic treatment strategy optimization

Reproducibility: random_state=42 fixed throughout.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.integrate import odeint
from scipy.spatial.distance import hamming
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler
from itertools import combinations, product
import networkx as nx
import warnings
import os
import random

warnings.filterwarnings('ignore')

# ─── Reproducibility ────────────────────────────────────────────────────────
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

# ─── Directories ─────────────────────────────────────────────────────────────
os.makedirs('figures', exist_ok=True)
os.makedirs('data/raw', exist_ok=True)

print("=" * 70)
print("AMR Evolution Prediction Computational Framework")
print("=" * 70)

# =============================================================================
# MODULE 1: ARG Detection from Whole-Genome Sequencing (Simulated)
# =============================================================================
print("\n[Module 1] ARG Detection from WGS Data")

# Simulate WGS data: 500 bacterial isolates, 50 potential ARG loci
N_ISOLATES = 500
N_LOCI = 50
N_RESISTANCE_GENES = 12  # known ARG families

np.random.seed(SEED)
# Simulate presence/absence matrix of genes (1=present, 0=absent)
# Higher prevalence for true ARGs (first 12 loci)
true_arg_prevalence = np.array([0.15, 0.22, 0.08, 0.31, 0.18, 0.25,
                                  0.12, 0.20, 0.09, 0.27, 0.14, 0.19])
noise_gene_prevalence = np.random.uniform(0.02, 0.06, N_LOCI - N_RESISTANCE_GENES)
all_prevalences = np.concatenate([true_arg_prevalence, noise_gene_prevalence])

genome_matrix = np.zeros((N_ISOLATES, N_LOCI))
for j in range(N_LOCI):
    genome_matrix[:, j] = np.random.binomial(1, all_prevalences[j], N_ISOLATES)

# Drug resistance phenotypes: 8 antibiotics
# Based on combinations of ARGs (with noise)
ANTIBIOTICS = ['Ampicillin', 'Ciprofloxacin', 'Tetracycline', 'Gentamicin',
               'Cefotaxime', 'Meropenem', 'Azithromycin', 'Trimethoprim']
ARG_DRUG_LINKS = {
    'Ampicillin':     [0, 4],
    'Ciprofloxacin':  [1, 5],
    'Tetracycline':   [2, 6],
    'Gentamicin':     [3, 7],
    'Cefotaxime':     [0, 4, 8],
    'Meropenem':      [9, 10],
    'Azithromycin':   [5, 11],
    'Trimethoprim':   [2, 3],
}

resistance_labels = {}
for drug, loci in ARG_DRUG_LINKS.items():
    # Resistant if any linked ARG present; add 10% noise
    base = (genome_matrix[:, loci].sum(axis=1) > 0).astype(float)
    noise = np.random.binomial(1, 0.10, N_ISOLATES)
    resistance_labels[drug] = ((base + noise) > 0).astype(int)

df_genomes = pd.DataFrame(genome_matrix,
                           columns=[f'locus_{i:02d}' for i in range(N_LOCI)])
df_phenotype = pd.DataFrame(resistance_labels)

# Save raw data
df_genomes.to_csv('data/raw/genome_matrix.csv', index=False)
df_phenotype.to_csv('data/raw/resistance_phenotypes.csv', index=False)

print(f"  Simulated {N_ISOLATES} isolates × {N_LOCI} genomic loci")
print(f"  Resistance prevalence per drug:")
for drug in ANTIBIOTICS:
    prev = df_phenotype[drug].mean()
    print(f"    {drug:20s}: {prev:.1%}")

# ARG Detection Model: Random Forest (5-fold CV)
print("\n  Training ARG-based resistance classifiers (5-fold CV):")
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
X = df_genomes.values

auroc_results = {}
for drug in ANTIBIOTICS:
    y = df_phenotype[drug].values
    clf = RandomForestClassifier(n_estimators=100, random_state=SEED, n_jobs=-1)
    scores = cross_val_score(clf, X, y, cv=skf, scoring='roc_auc')
    auroc_results[drug] = {'mean': scores.mean(), 'std': scores.std()}
    print(f"    {drug:20s}: AUROC = {scores.mean():.3f} ± {scores.std():.3f}")

# Save AUROC table
df_auroc = pd.DataFrame(auroc_results).T.reset_index()
df_auroc.columns = ['Drug', 'AUROC_mean', 'AUROC_std']
df_auroc.to_csv('data/raw/auroc_results.csv', index=False)

# Plot ARG heatmap
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
# Presence/absence heatmap (first 80 isolates)
sns.heatmap(genome_matrix[:80, :N_RESISTANCE_GENES],
            ax=axes[0], cmap='YlOrRd', cbar_kws={'label': 'Gene Present'},
            xticklabels=[f'ARG-{i+1}' for i in range(N_RESISTANCE_GENES)],
            yticklabels=False)
axes[0].set_title('ARG Presence/Absence Matrix\n(80 representative isolates)')
axes[0].set_xlabel('Antibiotic Resistance Genes')
axes[0].set_ylabel('Isolates')

# AUROC bar chart
drugs = [d['Drug'] for _, d in df_auroc.iterrows()]
means = [d['AUROC_mean'] for _, d in df_auroc.iterrows()]
stds = [d['AUROC_std'] for _, d in df_auroc.iterrows()]
colors = ['#2ecc71' if m >= 0.85 else '#e74c3c' for m in means]
bars = axes[1].bar(range(len(ANTIBIOTICS)), means, yerr=stds,
                    color=colors, alpha=0.8, capsize=4, edgecolor='black')
axes[1].axhline(0.8, color='navy', linestyle='--', label='AUROC=0.80 threshold')
axes[1].set_xticks(range(len(ANTIBIOTICS)))
axes[1].set_xticklabels(ANTIBIOTICS, rotation=45, ha='right', fontsize=9)
axes[1].set_ylim(0.5, 1.05)
axes[1].set_ylabel('AUROC (5-fold CV)')
axes[1].set_title('ARG Detection Performance\n(Random Forest, 5-fold CV)')
axes[1].legend()

plt.tight_layout()
plt.savefig('figures/fig1_arg_detection.png', dpi=150, bbox_inches='tight')
plt.close()
print("  → figures/fig1_arg_detection.png saved")

# =============================================================================
# MODULE 2: Fitness Landscape Construction
# =============================================================================
print("\n[Module 2] Fitness Landscape Construction")

# TEM β-lactamase: 4-site combinatorial fitness landscape
# Inspired by Weinreich et al. (2006) / Ogbunugafor et al.
# Each site: wildtype (0) or resistant mutation (1)
N_SITES = 4
N_GENOTYPES = 2 ** N_SITES  # 16 genotypes

np.random.seed(SEED)
# Epistatic fitness model: f(g) = baseline + Σ α_i*g_i + Σ β_ij*g_i*g_j + noise
# Parameters inspired by empirical β-lactamase landscapes
alpha = np.array([0.10, 0.08, 0.15, 0.06])  # additive effects
# Pairwise epistatic interactions (some positive, some negative)
beta = np.array([
    [0,    0.12, -0.05, 0.08],
    [0.12, 0,    0.07, -0.03],
    [-0.05, 0.07, 0,    0.11],
    [0.08, -0.03, 0.11, 0],
])

genotypes = list(product([0, 1], repeat=N_SITES))
fitness_values = {}
for g in genotypes:
    g_arr = np.array(g)
    f = 1.0  # wildtype fitness
    f += g_arr @ alpha
    f += 0.5 * g_arr @ beta @ g_arr
    f += np.random.normal(0, 0.02)  # measurement noise
    fitness_values[g] = max(0.01, f)  # floor at 0.01

df_fitness = pd.DataFrame([
    {'genotype': ''.join(map(str, g)),
     'n_mutations': sum(g),
     'fitness': fitness_values[g]}
    for g in genotypes
])
df_fitness.to_csv('data/raw/fitness_landscape.csv', index=False)

print(f"  Fitness landscape: {N_GENOTYPES} genotypes, {N_SITES} mutational sites")
print(f"  Fitness range: {df_fitness['fitness'].min():.3f} – {df_fitness['fitness'].max():.3f}")
print(f"  Most fit genotype: {df_fitness.loc[df_fitness['fitness'].idxmax(), 'genotype']} "
      f"(f={df_fitness['fitness'].max():.3f})")
print(f"  Wildtype fitness: {fitness_values[(0,0,0,0)]:.3f}")

# =============================================================================
# MODULE 3: Evolutionary Path Prediction
# =============================================================================
print("\n[Module 3] Evolutionary Path Prediction (Accessible Mutational Paths)")

# Enumerate all monotonically increasing fitness paths
# from wildtype (0000) to most-resistant (1111)
wildtype = (0, 0, 0, 0)
target = (1, 1, 1, 1)

def hamming_dist(g1, g2):
    return sum(a != b for a, b in zip(g1, g2))

def get_neighbors_towards_target(g, tgt):
    """Single-step mutations moving closer to target."""
    neighbors = []
    for i in range(len(g)):
        if g[i] != tgt[i]:
            new_g = list(g)
            new_g[i] = tgt[i]
            neighbors.append(tuple(new_g))
    return neighbors

def enumerate_accessible_paths(start, end, fitness_dict, min_steps=None):
    """Enumerate paths where each step increases fitness."""
    if min_steps is None:
        min_steps = hamming_dist(start, end)
    
    accessible = []
    inaccessible = []
    
    def dfs(path):
        current = path[-1]
        if current == end:
            accessible.append(list(path))
            return
        for neighbor in get_neighbors_towards_target(current, end):
            if fitness_dict[neighbor] >= fitness_dict[current]:
                dfs(path + [neighbor])
            else:
                # Dead end – record as inaccessible branch
                pass
    
    # All N_SITES! permutations of mutational order
    from itertools import permutations
    all_paths = []
    for perm in permutations(range(N_SITES)):
        path = [start]
        for site in perm:
            prev = path[-1]
            new = list(prev)
            new[site] = 1
            path.append(tuple(new))
        all_paths.append(path)
    
    for path in all_paths:
        accessible_flag = all(
            fitness_values[path[k+1]] >= fitness_values[path[k]]
            for k in range(len(path)-1)
        )
        if accessible_flag:
            accessible.append(path)
        else:
            inaccessible.append(path)
    
    return accessible, inaccessible

accessible_paths, inaccessible_paths = enumerate_accessible_paths(
    wildtype, target, fitness_values)
total_paths = len(accessible_paths) + len(inaccessible_paths)
accessibility_fraction = len(accessible_paths) / total_paths

print(f"  Total paths (wildtype→full resistance): {total_paths}")
print(f"  Accessible paths (monotone fitness): {len(accessible_paths)}")
print(f"  Inaccessible paths: {len(inaccessible_paths)}")
print(f"  Accessibility fraction: {accessibility_fraction:.1%}")

# Compute path probabilities (proportional to fitness differences)
def path_probability(path, fitness_dict, temperature=1.0):
    """Sella-Hirsh fixation probability model."""
    log_prob = 0.0
    for k in range(len(path) - 1):
        neighbors = get_neighbors_towards_target(path[k], target)
        delta_f = [fitness_dict[n] - fitness_dict[path[k]] for n in neighbors]
        # Fixation probability ∝ (1 - e^{-2Δf}) / (1 - e^{-2NΔf}) ~ Δf for Δf>0
        # Simplified: prob ∝ max(Δf, 0) / Σ max(Δf_j, 0)
        pos_delta = [max(d, 1e-10) for d in delta_f]
        total = sum(pos_delta)
        next_node_idx = neighbors.index(path[k+1]) if path[k+1] in neighbors else None
        if next_node_idx is not None and total > 0:
            log_prob += np.log(pos_delta[next_node_idx] / total)
        else:
            log_prob += -np.inf
    return np.exp(log_prob)

path_probs = []
for path in accessible_paths:
    p = path_probability(path, fitness_values)
    path_probs.append(p)

total_prob = sum(path_probs)
normalized_probs = [p / total_prob for p in path_probs]
most_likely_idx = np.argmax(normalized_probs)
most_likely_path = accessible_paths[most_likely_idx]
most_likely_prob = normalized_probs[most_likely_idx]

print(f"  Most probable evolutionary path: "
      f"{'→'.join([''.join(map(str,g)) for g in most_likely_path])}")
print(f"  Path probability: {most_likely_prob:.3f}")

# =============================================================================
# MODULE 4: HGT Network Modeling
# =============================================================================
print("\n[Module 4] HGT Network Modeling")

np.random.seed(SEED)
N_STRAINS = 30
N_PLASMID_TYPES = 5

# Generate strain metadata
strain_types = np.random.choice(['E.coli', 'K.pneumoniae', 'P.aeruginosa',
                                   'S.aureus', 'E.faecalis'], N_STRAINS)
strain_resistance = np.random.binomial(1, 0.45, N_STRAINS)

# Create HGT network: preferential attachment + strain-type bias
G_hgt = nx.DiGraph()
for i in range(N_STRAINS):
    G_hgt.add_node(i, species=strain_types[i], resistant=strain_resistance[i])

# Add edges: HGT events (plasmid transfer)
# Higher probability within same species, lower between species
hgt_events = []
for i in range(N_STRAINS):
    for j in range(N_STRAINS):
        if i == j:
            continue
        same_species = strain_types[i] == strain_types[j]
        p_transfer = 0.18 if same_species else 0.05
        # Donor must be resistant
        if strain_resistance[i] == 1 and np.random.random() < p_transfer:
            G_hgt.add_edge(i, j, plasmid_type=np.random.randint(N_PLASMID_TYPES))
            hgt_events.append({'donor': i, 'recipient': j,
                                'same_species': same_species})

df_hgt = pd.DataFrame(hgt_events)
df_hgt.to_csv('data/raw/hgt_events.csv', index=False)

# Network statistics
in_degrees = dict(G_hgt.in_degree())
out_degrees = dict(G_hgt.out_degree())
betweenness = nx.betweenness_centrality(G_hgt)
top_donors = sorted(out_degrees.items(), key=lambda x: x[1], reverse=True)[:5]

print(f"  HGT network: {G_hgt.number_of_nodes()} strains, "
      f"{G_hgt.number_of_edges()} transfer events")
print(f"  Within-species events: {df_hgt['same_species'].sum()} "
      f"({df_hgt['same_species'].mean():.1%})")
print(f"  Top 5 donor strains (out-degree): "
      f"{[(n, d) for n,d in top_donors]}")
print(f"  Network density: {nx.density(G_hgt):.4f}")

# Communities (Louvain-like using greedy modularity on undirected)
G_undirected = G_hgt.to_undirected()
communities = list(nx.community.greedy_modularity_communities(G_undirected))
print(f"  Detected {len(communities)} HGT transmission communities")

# =============================================================================
# MODULE 5: Spatiotemporal Resistance Dynamics (SIR + Resistance)
# =============================================================================
print("\n[Module 5] Spatiotemporal Resistance Dynamics")

# Extended SIR model with two bacterial strains:
# S: susceptible host, I_S: infected with susceptible bacteria,
# I_R: infected with resistant bacteria, R: recovered
# Parameters based on published AMR epidemiology literature

def amr_sir_model(y, t, beta_s, beta_r, gamma, mu, phi, kappa):
    """
    Extended SIR with resistance dynamics.
    
    y = [S, I_S, I_R, R]
    beta_s: transmission rate, susceptible strain
    beta_r: transmission rate, resistant strain  
    gamma: recovery rate
    mu: antibiotic usage (promotes resistance selection)
    phi: rate of resistance acquisition (de novo + HGT)
    kappa: fitness cost of resistance (0=no cost, 1=full cost)
    """
    S, I_S, I_R, R = y
    N = S + I_S + I_R + R
    
    # Force of infection
    lambda_s = beta_s * I_S / N
    lambda_r = beta_r * (1 - kappa) * I_R / N  # fitness cost reduces spread
    
    dS  = -lambda_s * S - lambda_r * S
    dI_S = lambda_s * S - gamma * I_S - mu * phi * I_S  # antibiotics clear susceptible
    dI_R = lambda_r * S + mu * phi * I_S - gamma * I_R  # resistance emerges
    dR  = gamma * (I_S + I_R)
    
    return [dS, dI_S, dI_R, dR]

# Parameter values (literature-informed)
params_base = {
    'beta_s': 0.35,   # transmission rate
    'beta_r': 0.30,   # resistant strain slightly less fit
    'gamma': 0.10,    # 10-day infection duration
    'mu': 0.40,       # antibiotic usage rate (moderate)
    'phi': 0.02,      # resistance emergence rate
    'kappa': 0.15,    # 15% fitness cost of resistance
}

# Initial conditions: 10% infected, 1% resistant
N_POP = 10000
y0 = [0.88 * N_POP, 0.10 * N_POP, 0.01 * N_POP, 0.01 * N_POP]
t = np.linspace(0, 365, 730)  # 1 year, 0.5-day steps

sol_base = odeint(amr_sir_model, y0, t,
                  args=(params_base['beta_s'], params_base['beta_r'],
                        params_base['gamma'], params_base['mu'],
                        params_base['phi'], params_base['kappa']))

# Compute resistance fraction over time
S_base, IS_base, IR_base, R_base = sol_base.T
# Clip to non-negative (ODE can drift slightly negative at very small values)
S_base  = np.clip(S_base,  0, None)
IS_base = np.clip(IS_base, 0, None)
IR_base = np.clip(IR_base, 0, None)
R_base  = np.clip(R_base,  0, None)
resistance_frac_base = np.clip(IR_base / (IS_base + IR_base + 1e-10), 0, 1)

# Scenario analysis: different antibiotic usage levels
mu_values = [0.1, 0.2, 0.4, 0.6, 0.8]
resistance_outcomes = {}
for mu_val in mu_values:
    sol = odeint(amr_sir_model, y0, t,
                 args=(params_base['beta_s'], params_base['beta_r'],
                       params_base['gamma'], mu_val,
                       params_base['phi'], params_base['kappa']))
    IR_t = np.clip(sol[:, 2], 0, None)
    IS_t = np.clip(sol[:, 1], 0, None)
    rf = np.clip(IR_t / (IS_t + IR_t + 1e-10), 0, 1)
    resistance_outcomes[mu_val] = rf

# Compute R_eff over time
R_eff_s = params_base['beta_s'] / params_base['gamma'] * S_base / N_POP
R_eff_r = params_base['beta_r'] * (1 - params_base['kappa']) / params_base['gamma'] * S_base / N_POP

print(f"  Baseline scenario (μ={params_base['mu']}):")
print(f"    Peak infected (susceptible): {IS_base.max():.0f} on day {t[IS_base.argmax()]:.0f}")
print(f"    Peak infected (resistant): {IR_base.max():.0f} on day {t[IR_base.argmax()]:.0f}")
print(f"    Final resistance fraction: {resistance_frac_base[-1]:.3f}")
print(f"    R_eff(susceptible) at peak: {R_eff_s[IS_base.argmax()]:.2f}")

df_dynamics = pd.DataFrame({
    'time': t,
    'S': S_base, 'I_S': IS_base, 'I_R': IR_base, 'R': R_base,
    'resistance_fraction': resistance_frac_base,
    'R_eff_s': R_eff_s, 'R_eff_r': R_eff_r
})
df_dynamics.to_csv('data/raw/sir_dynamics.csv', index=False)

# =============================================================================
# MODULE 6: Antibiotic Treatment Strategy Optimization
# =============================================================================
print("\n[Module 6] Treatment Strategy Optimization")

# Simulate: cycling vs combination vs monotherapy
# Metric: resistance fraction at day 365, total bacterial burden

def simulate_strategy(strategy, t_span, y0, base_params, seed=42):
    """
    Simulate treatment strategies over time.
    strategy: 'monotherapy', 'cycling', 'combination'
    """
    np.random.seed(seed)
    n_steps = len(t_span)
    dt = t_span[1] - t_span[0]
    
    results = {'resistance_frac': [], 'total_infected': []}
    y = list(y0)
    
    for k in range(n_steps - 1):
        t_k = t_span[k]
        
        # Determine effective mu based on strategy
        if strategy == 'monotherapy':
            mu_eff = base_params['mu']
            phi_eff = base_params['phi']
        elif strategy == 'cycling':
            # Alternate antibiotics every 30 days
            cycle_phase = int(t_k / 30) % 2
            mu_eff = base_params['mu']
            phi_eff = base_params['phi'] * (0.5 if cycle_phase == 0 else 1.5)
        elif strategy == 'combination':
            # Two antibiotics together: lower emergence rate but higher clearance
            mu_eff = base_params['mu'] * 1.2
            phi_eff = base_params['phi'] * 0.3  # synergistic effect
        
        dy = amr_sir_model(y, t_k, base_params['beta_s'], base_params['beta_r'],
                           base_params['gamma'], mu_eff, phi_eff, base_params['kappa'])
        y = [max(0, y[i] + dy[i] * dt) for i in range(4)]
        
        S_, IS_, IR_, R_ = y
        rf = IR_ / (IS_ + IR_ + 1e-10)
        results['resistance_frac'].append(rf)
        results['total_infected'].append(IS_ + IR_)
    
    return results

t_opt = np.linspace(0, 365, 365)
strategies = ['monotherapy', 'cycling', 'combination']
strategy_results = {}
for strat in strategies:
    strategy_results[strat] = simulate_strategy(strat, t_opt, y0, params_base)

print("  Strategy comparison at day 365:")
for strat in strategies:
    final_rf = strategy_results[strat]['resistance_frac'][-1]
    final_burden = strategy_results[strat]['total_infected'][-1]
    print(f"    {strat:15s}: resistance={final_rf:.3f}, burden={final_burden:.0f}")

# Optimization: find optimal cycling period
print("\n  Optimizing cycling period (7–90 days):")
cycle_periods = [7, 14, 21, 30, 45, 60, 90]
cycle_outcomes = {}
for period in cycle_periods:
    total_rf = 0
    n_pts = 0
    for k, t_k in enumerate(t_opt[1:]):
        pass  # simplified – use final resistance fraction
    
    # Simulate with variable cycling period
    y = list(y0)
    dt = t_opt[1] - t_opt[0]
    rf_series = []
    for k in range(len(t_opt) - 1):
        t_k = t_opt[k]
        cycle_phase = int(t_k / period) % 2
        phi_eff = params_base['phi'] * (0.5 if cycle_phase == 0 else 1.5)
        dy = amr_sir_model(y, t_k, params_base['beta_s'], params_base['beta_r'],
                           params_base['gamma'], params_base['mu'], phi_eff,
                           params_base['kappa'])
        y = [max(0, y[i] + dy[i] * dt) for i in range(4)]
        S_, IS_, IR_, R_ = y
        rf_series.append(IR_ / (IS_ + IR_ + 1e-10))
    cycle_outcomes[period] = {'final_rf': rf_series[-1],
                               'mean_rf': np.mean(rf_series)}
    print(f"    Period={period:3d}d: final_resistance={rf_series[-1]:.4f}, "
          f"mean_resistance={np.mean(rf_series):.4f}")

best_period = min(cycle_outcomes.keys(), key=lambda p: cycle_outcomes[p]['final_rf'])
print(f"  → Optimal cycling period: {best_period} days "
      f"(final resistance = {cycle_outcomes[best_period]['final_rf']:.4f})")

df_cycle = pd.DataFrame([
    {'period': p, **v} for p, v in cycle_outcomes.items()
])
df_cycle.to_csv('data/raw/cycling_optimization.csv', index=False)

# =============================================================================
# VISUALIZATION: Comprehensive Figure Panel
# =============================================================================
print("\n[Visualization] Creating figure panels...")

# Figure 2: Fitness Landscape
fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5))

# 2a: Fitness by number of mutations (violin)
ax = axes2[0]
data_by_nmut = [df_fitness[df_fitness['n_mutations'] == n]['fitness'].values
                for n in range(N_SITES + 1)]
parts = ax.violinplot(data_by_nmut, positions=range(N_SITES + 1),
                       showmeans=True, showmedians=False)
for pc in parts['bodies']:
    pc.set_facecolor('#3498db')
    pc.set_alpha(0.7)
ax.set_xlabel('Number of Mutations')
ax.set_ylabel('Relative Fitness')
ax.set_title('Fitness Distribution\nby Mutation Count')
ax.set_xticks(range(N_SITES + 1))

# 2b: Fitness landscape heatmap (2D slice: sites 0 & 1, sites 2 & 3 fixed)
ax = axes2[1]
landscape_2d = np.zeros((2, 2, 2, 2))
for g in genotypes:
    landscape_2d[g[0], g[1], g[2], g[3]] = fitness_values[g]
# Mean over sites 2,3
landscape_slice = landscape_2d.mean(axis=(2, 3))
im = ax.imshow(landscape_slice, cmap='RdYlGn', aspect='auto',
               interpolation='nearest', vmin=0.9, vmax=1.55)
plt.colorbar(im, ax=ax, label='Mean Fitness')
ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
ax.set_xticklabels(['WT', 'Mut']); ax.set_yticklabels(['WT', 'Mut'])
ax.set_xlabel('Site 1'); ax.set_ylabel('Site 0')
ax.set_title('2D Fitness Landscape Slice\n(Sites 0 & 1, averaged over Sites 2 & 3)')
for i in range(2):
    for j in range(2):
        ax.text(j, i, f'{landscape_slice[i, j]:.3f}', ha='center', va='center',
                fontsize=11, fontweight='bold')

# 2c: Accessible paths visualization
ax = axes2[2]
path_lengths = [len(p) for p in accessible_paths]
path_fitness_gains = [fitness_values[p[-1]] - fitness_values[p[0]]
                       for p in accessible_paths]
scatter = ax.scatter(path_lengths, path_fitness_gains,
                      c=normalized_probs, cmap='viridis', s=80, alpha=0.8)
plt.colorbar(scatter, ax=ax, label='Normalized Probability')
ax.set_xlabel('Path Length (steps)')
ax.set_ylabel('Total Fitness Gain')
ax.set_title(f'Accessible Evolutionary Paths\n({len(accessible_paths)}/{total_paths} paths, '
             f'{accessibility_fraction:.0%} accessible)')
ax.axvline(np.mean(path_lengths), color='red', linestyle='--', alpha=0.7,
           label=f'Mean length={np.mean(path_lengths):.1f}')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('figures/fig2_fitness_landscape.png', dpi=150, bbox_inches='tight')
plt.close()
print("  → figures/fig2_fitness_landscape.png saved")

# Figure 3: HGT Network
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 6))

ax = axes3[0]
species_colors = {
    'E.coli': '#e74c3c', 'K.pneumoniae': '#3498db',
    'P.aeruginosa': '#2ecc71', 'S.aureus': '#f39c12',
    'E.faecalis': '#9b59b6'
}
node_colors = [species_colors[G_hgt.nodes[n]['species']] for n in G_hgt.nodes()]
node_sizes = [200 + 50 * G_hgt.out_degree(n) for n in G_hgt.nodes()]
node_shapes = ['s' if G_hgt.nodes[n]['resistant'] else 'o' for n in G_hgt.nodes()]

pos = nx.spring_layout(G_hgt, seed=SEED, k=1.5)
# Draw by species
for species, color in species_colors.items():
    nodes_sp = [n for n in G_hgt.nodes() if G_hgt.nodes[n]['species'] == species]
    nx.draw_networkx_nodes(G_hgt, pos, nodelist=nodes_sp, node_color=color,
                            node_size=150, ax=ax, label=species, alpha=0.8)

nx.draw_networkx_edges(G_hgt, pos, ax=ax, alpha=0.3, arrows=True,
                        arrowsize=10, edge_color='gray', connectionstyle='arc3,rad=0.1')
ax.set_title(f'HGT Transmission Network\n({G_hgt.number_of_nodes()} strains, '
              f'{G_hgt.number_of_edges()} transfer events)')
ax.legend(fontsize=8, loc='upper left')
ax.axis('off')

# Degree distribution
ax = axes3[1]
in_deg_vals = [d for _, d in G_hgt.in_degree()]
out_deg_vals = [d for _, d in G_hgt.out_degree()]
ax.hist(in_deg_vals, bins=range(max(in_deg_vals)+2), alpha=0.6,
        label='In-degree (ARG acquisition)', color='#e74c3c', edgecolor='black')
ax.hist(out_deg_vals, bins=range(max(out_deg_vals)+2), alpha=0.6,
        label='Out-degree (ARG donation)', color='#3498db', edgecolor='black')
ax.set_xlabel('Node Degree')
ax.set_ylabel('Frequency')
ax.set_title('HGT Network Degree Distribution')
ax.legend()

plt.tight_layout()
plt.savefig('figures/fig3_hgt_network.png', dpi=150, bbox_inches='tight')
plt.close()
print("  → figures/fig3_hgt_network.png saved")

# Figure 4: Spatiotemporal Dynamics + Treatment Optimization
fig4, axes4 = plt.subplots(2, 2, figsize=(14, 10))

# 4a: SIR dynamics
ax = axes4[0, 0]
ax.plot(t, S_base / N_POP, 'b-', linewidth=2, label='Susceptible hosts (S)')
ax.plot(t, IS_base / N_POP, 'g-', linewidth=2, label='Infected (susceptible strain)')
ax.plot(t, IR_base / N_POP, 'r-', linewidth=2, label='Infected (resistant strain)')
ax.plot(t, R_base / N_POP, 'gray', linewidth=2, label='Recovered (R)')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Population Fraction')
ax.set_title(f'AMR-SIR Dynamics (μ={params_base["mu"]})')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# 4b: Resistance fraction vs antibiotic usage
ax = axes4[0, 1]
for mu_val, rf_series in resistance_outcomes.items():
    ax.plot(t, rf_series, linewidth=2, label=f'μ={mu_val}')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Resistance Fraction I_R / (I_S + I_R)')
ax.set_title('Resistance Dynamics vs\nAntibiotic Usage Intensity')
ax.legend(title='Antibiotic usage (μ)', fontsize=8)
ax.grid(alpha=0.3)

# 4c: Treatment strategy comparison
ax = axes4[1, 0]
colors_strat = {'monotherapy': '#e74c3c', 'cycling': '#f39c12', 'combination': '#2ecc71'}
for strat in strategies:
    ax.plot(t_opt[1:], strategy_results[strat]['resistance_frac'],
            linewidth=2.5, label=strat.capitalize(),
            color=colors_strat[strat])
ax.set_xlabel('Time (days)')
ax.set_ylabel('Resistance Fraction')
ax.set_title('Treatment Strategy Comparison')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# 4d: Cycling period optimization
ax = axes4[1, 1]
periods = list(cycle_outcomes.keys())
final_rfs = [cycle_outcomes[p]['final_rf'] for p in periods]
mean_rfs = [cycle_outcomes[p]['mean_rf'] for p in periods]
ax.bar([str(p) for p in periods], final_rfs, alpha=0.7, color='#3498db',
        label='Final resistance fraction', edgecolor='black')
ax.plot([str(p) for p in periods], mean_rfs, 'ro-', linewidth=2,
         label='Mean resistance fraction', markersize=8)
ax.axvline(str(best_period), color='green', linestyle='--', linewidth=2,
            label=f'Optimal period ({best_period}d)')
ax.set_xlabel('Cycling Period (days)')
ax.set_ylabel('Resistance Fraction at Day 365')
ax.set_title('Cycling Period Optimization')
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('figures/fig4_dynamics_optimization.png', dpi=150, bbox_inches='tight')
plt.close()
print("  → figures/fig4_dynamics_optimization.png saved")

# Figure 5: Population Genetics Integration
print("\n[Module 7] Population Genetics Integration")

np.random.seed(SEED)
N_POP_GEN = 1000   # effective population size
N_GENERATIONS = 200
MUTATION_RATE = 1e-4
SELECTION_COEFF = 0.05  # selective advantage of resistance under antibiotics
ANTIBIOTIC_PRES = 0.6   # fraction of time antibiotics are present

# Wright-Fisher simulation with selection
def wright_fisher_amr(N, n_gen, mu, s, p_antibiotic, seed=42):
    """
    Wright-Fisher model with selection for AMR.
    Tracks frequency of resistance allele over generations.
    """
    np.random.seed(seed)
    freq = np.zeros(n_gen)
    freq[0] = 1 / N  # start with single resistant mutant
    
    for t in range(1, n_gen):
        p = freq[t-1]
        # Selection: antibiotic present with probability p_antibiotic
        antibiotic_present = np.random.random() < p_antibiotic
        if antibiotic_present:
            # Resistant strain has advantage s
            p_selected = p * (1 + s) / (p * (1 + s) + (1 - p))
        else:
            # Slightly costly without antibiotics
            p_selected = p * (1 - 0.02) / (p * (1 - 0.02) + (1 - p))
        
        # Mutation
        p_mut = p_selected + mu * (1 - p_selected) - mu * p_selected
        
        # Genetic drift
        new_count = np.random.binomial(N, np.clip(p_mut, 0, 1))
        freq[t] = new_count / N
    
    return freq

# Run multiple simulations
N_SIMS = 50
trajectories = []
for i in range(N_SIMS):
    traj = wright_fisher_amr(N_POP_GEN, N_GENERATIONS, MUTATION_RATE,
                               SELECTION_COEFF, ANTIBIOTIC_PRES, seed=SEED+i)
    trajectories.append(traj)

trajectories = np.array(trajectories)
mean_traj = trajectories.mean(axis=0)
std_traj = trajectories.std(axis=0)

# Fixation probability
fixation_count = sum(t[-1] > 0.5 for t in trajectories)
fixation_prob = fixation_count / N_SIMS

print(f"  Wright-Fisher AMR simulation:")
print(f"    N_eff={N_POP_GEN}, s={SELECTION_COEFF}, "
      f"p_antibiotic={ANTIBIOTIC_PRES}")
print(f"    Fixation probability: {fixation_prob:.2f} ({fixation_count}/{N_SIMS})")
print(f"    Mean final frequency: {mean_traj[-1]:.3f} ± {std_traj[-1]:.3f}")

# Compute Tajima's D analog (simplified for simulated data)
# Using final frequency distribution across simulations
final_freqs = trajectories[:, -1]
tajima_stat, tajima_p = stats.ttest_1samp(final_freqs, 0.5)
print(f"  Frequency vs 0.5 (t-test): t={tajima_stat:.3f}, p={tajima_p:.4f}")

fig5, axes5 = plt.subplots(1, 3, figsize=(16, 5))

# 5a: WF trajectories
ax = axes5[0]
gen_idx = np.arange(N_GENERATIONS)
for i in range(min(20, N_SIMS)):
    ax.plot(gen_idx, trajectories[i], alpha=0.2, color='steelblue', linewidth=0.8)
ax.plot(gen_idx, mean_traj, 'r-', linewidth=2.5, label='Mean trajectory')
ax.fill_between(gen_idx, mean_traj - std_traj, mean_traj + std_traj,
                 alpha=0.3, color='red', label='±1 SD')
ax.set_xlabel('Generation')
ax.set_ylabel('Resistance Allele Frequency')
ax.set_title(f'Wright-Fisher AMR Simulation\n(N={N_POP_GEN}, s={SELECTION_COEFF}, '
              f'p_ab={ANTIBIOTIC_PRES})')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# 5b: Final frequency distribution
ax = axes5[1]
ax.hist(final_freqs, bins=20, color='#9b59b6', alpha=0.8, edgecolor='black')
ax.axvline(mean_traj[-1], color='red', linestyle='--', linewidth=2,
            label=f'Mean={mean_traj[-1]:.3f}')
ax.axvline(0.5, color='gray', linestyle=':', linewidth=2, label='0.5 threshold')
ax.set_xlabel('Final Resistance Frequency')
ax.set_ylabel('Count')
ax.set_title(f'Distribution of Final Frequencies\n(Fixation: {fixation_prob:.0%})')
ax.legend(fontsize=9)

# 5c: Selection coefficient sensitivity
s_values = np.arange(0.01, 0.16, 0.01)
fixation_probs_by_s = []
for s_val in s_values:
    trajs = [wright_fisher_amr(N_POP_GEN, N_GENERATIONS, MUTATION_RATE,
                                s_val, ANTIBIOTIC_PRES, seed=SEED+k)
             for k in range(20)]
    fp = sum(tr[-1] > 0.5 for tr in trajs) / 20
    fixation_probs_by_s.append(fp)

ax = axes5[2]
ax.plot(s_values, fixation_probs_by_s, 'bo-', linewidth=2, markersize=6)
# Kimura formula: P_fix ≈ (1-e^{-2s}) / (1-e^{-2Ns}) for large N
kimura_fix = [(1 - np.exp(-2 * s)) / (1 - np.exp(-2 * N_POP_GEN * s))
              for s in s_values]
ax.plot(s_values, kimura_fix, 'r--', linewidth=2, label='Kimura formula')
ax.set_xlabel('Selection Coefficient (s)')
ax.set_ylabel('Fixation Probability')
ax.set_title('Fixation Probability vs\nSelection Coefficient')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('figures/fig5_population_genetics.png', dpi=150, bbox_inches='tight')
plt.close()
print("  → figures/fig5_population_genetics.png saved")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY OF RESULTS")
print("=" * 70)

mean_auroc = np.mean([v['mean'] for v in auroc_results.values()])
std_auroc = np.std([v['mean'] for v in auroc_results.values()])
print(f"\nModule 1 (ARG Detection):")
print(f"  Mean AUROC across 8 drugs: {mean_auroc:.3f} ± {std_auroc:.3f}")

print(f"\nModule 2 (Fitness Landscape):")
print(f"  WT fitness: {fitness_values[(0,0,0,0)]:.3f}")
print(f"  Max fitness (MGRG): {df_fitness['fitness'].max():.3f}")
print(f"  Fitness gain (WT→full-res): "
      f"{fitness_values[target] - fitness_values[wildtype]:.3f}")

print(f"\nModule 3 (Evolutionary Paths):")
print(f"  Accessible paths: {len(accessible_paths)}/{total_paths} ({accessibility_fraction:.1%})")
print(f"  Most probable path probability: {most_likely_prob:.3f}")

print(f"\nModule 4 (HGT Network):")
print(f"  Network density: {nx.density(G_hgt):.4f}")
print(f"  Communities: {len(communities)}")

print(f"\nModule 5 (Spatiotemporal Dynamics):")
print(f"  Final resistance fraction (μ=0.4): {resistance_frac_base[-1]:.3f}")
print(f"  R_eff(susceptible) at epidemic peak: "
      f"{R_eff_s[IS_base.argmax()]:.2f}")

print(f"\nModule 6 (Treatment Optimization):")
for strat in strategies:
    rf = strategy_results[strat]['resistance_frac'][-1]
    print(f"  {strat:15s}: final resistance fraction = {rf:.4f}")
print(f"  Optimal cycling period: {best_period} days")

print(f"\nModule 7 (Population Genetics):")
print(f"  Fixation probability: {fixation_prob:.2f}")
print(f"  Mean final frequency: {mean_traj[-1]:.3f} ± {std_traj[-1]:.3f}")
print(f"  t-test vs 0.5: t={tajima_stat:.3f}, p={tajima_p:.4f}")

print("\n  All figures saved to figures/")
print("  All raw data saved to data/raw/")
print("\n[DONE]")
