"""
Module 3: Species Interaction Network Model
Lotka-Volterra based community dynamics with predation, competition, and symbiosis.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import networkx as nx


# === Reef Community Functional Groups ===
SPECIES = {
    0: {'name': 'Hard Coral (Acropora)',     'type': 'coral',     'K': 0.60, 'r': 0.10},
    1: {'name': 'Hard Coral (Porites)',      'type': 'coral',     'K': 0.40, 'r': 0.06},
    2: {'name': 'Soft Coral',                'type': 'coral',     'K': 0.20, 'r': 0.08},
    3: {'name': 'Crustose Coralline Algae',  'type': 'CCA',      'K': 0.30, 'r': 0.05},
    4: {'name': 'Macroalgae',                'type': 'algae',     'K': 0.80, 'r': 0.15},
    5: {'name': 'Turf Algae',                'type': 'algae',     'K': 0.70, 'r': 0.20},
    6: {'name': 'Zooxanthellae',             'type': 'symbiont',  'K': 1.00, 'r': 0.12},
    7: {'name': 'Herbivorous Fish',          'type': 'herbivore', 'K': 0.50, 'r': 0.08},
    8: {'name': 'Corallivorous Fish',        'type': 'predator',  'K': 0.30, 'r': 0.05},
    9: {'name': 'Crown-of-Thorns Starfish',  'type': 'predator',  'K': 0.15, 'r': 0.12},
    10:{'name': 'Sea Urchins',               'type': 'herbivore', 'K': 0.25, 'r': 0.07},
    11:{'name': 'Reef Fish (general)',       'type': 'fish',      'K': 0.40, 'r': 0.06},
}

N_SPECIES = len(SPECIES)

def build_interaction_matrix():
    """
    Build the species interaction matrix.
    Positive: benefit (mutualism/prey benefit to predator)
    Negative: cost (competition/predation on prey)
    """
    # alpha[i,j]: effect of species j on species i
    alpha = np.zeros((N_SPECIES, N_SPECIES))

    # Competition between corals
    alpha[0, 1] = -0.3;  alpha[1, 0] = -0.3   # Acropora-Porites
    alpha[0, 2] = -0.2;  alpha[2, 0] = -0.2   # Acropora-Soft coral
    alpha[1, 2] = -0.15; alpha[2, 1] = -0.15  # Porites-Soft coral

    # Coral-algae competition (space)
    for coral in [0, 1, 2]:
        for algae in [4, 5]:
            alpha[coral, algae] = -0.4  # Algae harm coral
            alpha[algae, coral] = -0.2  # Coral limits algae space

    # CCA-algae competition
    alpha[3, 4] = -0.3; alpha[4, 3] = -0.2
    alpha[3, 5] = -0.25; alpha[5, 3] = -0.15

    # CCA facilitates coral recruitment
    alpha[0, 3] = 0.15; alpha[1, 3] = 0.12

    # Coral-Zooxanthellae mutualism
    for coral in [0, 1, 2]:
        alpha[coral, 6] = 0.5   # Zoox benefit coral
        alpha[6, coral] = 0.4   # Coral benefit zoox

    # Herbivore-algae (predation)
    alpha[7, 4] = 0.3;  alpha[4, 7] = -0.5   # Herbivorous fish eat macroalgae
    alpha[7, 5] = 0.25; alpha[5, 7] = -0.4   # Herbivorous fish eat turf
    alpha[10, 4] = 0.2; alpha[4, 10] = -0.3  # Urchins eat macroalgae
    alpha[10, 5] = 0.2; alpha[5, 10] = -0.3  # Urchins eat turf

    # Corallivore-coral (predation)
    alpha[8, 0] = 0.3;  alpha[0, 8] = -0.3   # Corallivores eat Acropora
    alpha[8, 1] = 0.2;  alpha[1, 8] = -0.2   # Corallivores eat Porites
    alpha[9, 0] = 0.5;  alpha[0, 9] = -0.5   # COTS eat Acropora (strong)
    alpha[9, 1] = 0.3;  alpha[1, 9] = -0.3   # COTS eat Porites

    # Reef fish depend on coral habitat
    alpha[11, 0] = 0.3; alpha[11, 1] = 0.2; alpha[11, 2] = 0.15

    return alpha


def build_network_graph(alpha):
    """Build NetworkX directed graph from interaction matrix."""
    G = nx.DiGraph()
    for i in range(N_SPECIES):
        G.add_node(i, **SPECIES[i])
    
    for i in range(N_SPECIES):
        for j in range(N_SPECIES):
            if i != j and abs(alpha[i, j]) > 0.01:
                interaction_type = 'positive' if alpha[i, j] > 0 else 'negative'
                G.add_edge(j, i, weight=alpha[i, j], interaction=interaction_type)
    
    return G


def reef_community_ode(t, N, r, K, alpha, stress_factors=None):
    """
    Generalized Lotka-Volterra community dynamics.
    
    dN_i/dt = r_i * N_i * (1 - (N_i + Σ α_ij * N_j) / K_i) * stress_i
    """
    N = np.maximum(N, 0)
    dNdt = np.zeros(N_SPECIES)
    
    for i in range(N_SPECIES):
        if N[i] < 1e-6:
            dNdt[i] = 0
            continue
        
        interaction_sum = 0
        for j in range(N_SPECIES):
            if i != j:
                interaction_sum += alpha[i, j] * N[j]
        
        stress = 1.0
        if stress_factors is not None:
            stress = stress_factors.get(i, 1.0)
        
        dNdt[i] = r[i] * N[i] * (1 - (N[i] + interaction_sum) / K[i]) * stress
    
    return dNdt


def simulate_community(years=80, scenario='baseline', stress_profile=None):
    """
    Simulate reef community dynamics over time.
    
    Parameters
    ----------
    years : int - Simulation duration (years from 2020)
    scenario : str - Stress scenario name
    stress_profile : dict - Time-varying stress factors
    """
    alpha = build_interaction_matrix()
    r = np.array([SPECIES[i]['r'] for i in range(N_SPECIES)])
    K = np.array([SPECIES[i]['K'] for i in range(N_SPECIES)])
    
    # Initial conditions (proportional cover/abundance)
    N0 = np.array([0.30, 0.20, 0.08, 0.15, 0.05, 0.05, 0.50, 0.30, 0.10, 0.02, 0.08, 0.25])
    
    t_span = (0, years)
    t_eval = np.linspace(0, years, years * 12)  # monthly
    
    if stress_profile is None:
        # Default: no additional stress
        sol = solve_ivp(
            reef_community_ode, t_span, N0,
            args=(r, K, alpha, None),
            t_eval=t_eval, method='RK45',
            max_step=0.1
        )
    else:
        # Time-varying stress: use piecewise integration
        dt = 1.0  # yearly steps
        t_all = [0]
        N_all = [N0.copy()]
        N_current = N0.copy()
        
        for year in range(years):
            stress = {}
            for sp_id, sp_stress_func in stress_profile.items():
                stress[sp_id] = sp_stress_func(year)
            
            sol_step = solve_ivp(
                reef_community_ode, (0, dt), N_current,
                args=(r, K, alpha, stress),
                t_eval=np.linspace(0, dt, 12),
                method='RK45', max_step=0.1
            )
            N_current = np.maximum(sol_step.y[:, -1], 0)
            for month_idx in range(12):
                t_all.append(year + (month_idx + 1) / 12)
                N_all.append(sol_step.y[:, min(month_idx, sol_step.y.shape[1]-1)])
        
        class SolContainer:
            pass
        sol = SolContainer()
        sol.t = np.array(t_all)
        sol.y = np.array(N_all).T
    
    return sol, alpha


def run_network_analysis():
    """Run full network analysis and save results."""
    alpha = build_interaction_matrix()
    G = build_network_graph(alpha)
    
    # Network metrics
    metrics = {
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        'density': round(nx.density(G), 4),
        'species_centrality': {},
        'interaction_summary': {
            'positive': sum(1 for _, _, d in G.edges(data=True) if d['weight'] > 0),
            'negative': sum(1 for _, _, d in G.edges(data=True) if d['weight'] < 0),
        }
    }
    
    # Degree centrality
    in_deg = nx.in_degree_centrality(G)
    out_deg = nx.out_degree_centrality(G)
    betw = nx.betweenness_centrality(G)
    
    for i in range(N_SPECIES):
        metrics['species_centrality'][SPECIES[i]['name']] = {
            'in_degree': round(in_deg[i], 3),
            'out_degree': round(out_deg[i], 3),
            'betweenness': round(betw[i], 3)
        }
    
    # Run baseline simulation
    sol_base, _ = simulate_community(years=80, scenario='baseline')
    
    # Define OA stress profile (RCP8.5)
    def make_oa_stress(sp_type, severity):
        def stress_func(year):
            # Linearly increasing stress
            oa_factor = max(0.2, 1.0 - severity * year / 80)
            return oa_factor
        return stress_func
    
    stress_rcp85 = {
        0: make_oa_stress('coral', 0.6),    # Acropora most sensitive
        1: make_oa_stress('coral', 0.4),    # Porites more tolerant
        2: make_oa_stress('coral', 0.5),    # Soft coral
        3: make_oa_stress('CCA', 0.7),      # CCA very sensitive
        6: make_oa_stress('symbiont', 0.3), # Zooxanthellae
    }
    
    sol_rcp85, _ = simulate_community(years=80, scenario='RCP8.5', stress_profile=stress_rcp85)
    
    # Extract key timepoints
    community_results = {
        'baseline': {},
        'RCP8.5': {}
    }
    
    for i in range(N_SPECIES):
        name = SPECIES[i]['name']
        community_results['baseline'][name] = {
            'initial': round(float(sol_base.y[i, 0]), 4),
            'year_2060': round(float(sol_base.y[i, min(40*12, sol_base.y.shape[1]-1)]), 4),
            'year_2100': round(float(sol_base.y[i, -1]), 4)
        }
        community_results['RCP8.5'][name] = {
            'initial': round(float(sol_rcp85.y[i, 0]), 4),
            'year_2060': round(float(sol_rcp85.y[i, min(40*12, sol_rcp85.y.shape[1]-1)]), 4),
            'year_2100': round(float(sol_rcp85.y[i, -1]), 4)
        }
    
    # Interaction matrix for export
    alpha_list = alpha.tolist()
    species_names = [SPECIES[i]['name'] for i in range(N_SPECIES)]
    
    all_results = {
        'network_metrics': metrics,
        'community_dynamics': community_results,
        'interaction_matrix': alpha_list,
        'species_names': species_names
    }
    
    with open('results/network_model.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    return all_results, sol_base, sol_rcp85


if __name__ == '__main__':
    results, sol_base, sol_rcp85 = run_network_analysis()
    print("=== Network Analysis Results ===")
    print(f"Nodes: {results['network_metrics']['n_nodes']}")
    print(f"Edges: {results['network_metrics']['n_edges']}")
    print(f"Positive interactions: {results['network_metrics']['interaction_summary']['positive']}")
    print(f"Negative interactions: {results['network_metrics']['interaction_summary']['negative']}")
    print("\nResults saved to results/network_model.json")
