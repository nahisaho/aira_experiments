"""
Visualization script: Generate all figures for the integrated reef model.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json
import os

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

COLORS = {
    'RCP2.6': '#2166ac',
    'RCP4.5': '#fdae61',
    'RCP8.5': '#d73027',
}

SPECIES_COLORS = {
    'Acropora': '#e41a1c',
    'Porites': '#377eb8',
    'Montipora': '#4daf4a',
    'Stylophora': '#984ea3',
    'Pavona': '#ff7f00',
    'CCA': '#a65628',
}


def fig1_carbonate_chemistry():
    """Figure 1: Carbonate chemistry projections under RCP scenarios."""
    with open('results/carbonate_chemistry.json') as f:
        data = json.load(f)
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    
    for scenario in ['RCP2.6', 'RCP4.5', 'RCP8.5']:
        years = [d['year'] for d in data[scenario]]
        pH_vals = [d['pH'] for d in data[scenario]]
        omega_vals = [d['Omega_aragonite'] for d in data[scenario]]
        pCO2_vals = [d['pCO2_uatm'] for d in data[scenario]]
        CO3_vals = [d['CO3_umol_kg'] for d in data[scenario]]
        
        axes[0,0].plot(years, pH_vals, color=COLORS[scenario], label=scenario, linewidth=2)
        axes[0,1].plot(years, omega_vals, color=COLORS[scenario], label=scenario, linewidth=2)
        axes[1,0].plot(years, pCO2_vals, color=COLORS[scenario], label=scenario, linewidth=2)
        axes[1,1].plot(years, CO3_vals, color=COLORS[scenario], label=scenario, linewidth=2)
    
    axes[0,0].set_ylabel('pH (total scale)')
    axes[0,0].set_title('(a) Seawater pH')
    axes[0,0].axhline(y=7.8, color='gray', linestyle='--', alpha=0.5, label='Critical threshold')
    
    axes[0,1].set_ylabel('Ω aragonite')
    axes[0,1].set_title('(b) Aragonite Saturation State')
    axes[0,1].axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Dissolution threshold')
    axes[0,1].axhline(y=3.0, color='gray', linestyle=':', alpha=0.5, label='Optimal calcification')
    
    axes[1,0].set_ylabel('pCO₂ (µatm)')
    axes[1,0].set_title('(c) Seawater pCO₂')
    
    axes[1,1].set_ylabel('[CO₃²⁻] (µmol/kg)')
    axes[1,1].set_title('(d) Carbonate Ion Concentration')
    
    for ax in axes.flat:
        ax.set_xlabel('Year')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    fig.suptitle('Projected Carbonate Chemistry (GBR Conditions)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig('figures/fig1_carbonate_chemistry.png')
    fig.savefig('figures/fig1_carbonate_chemistry.svg')
    plt.close()
    print("  Fig 1 saved")


def fig2_calcification():
    """Figure 2: Calcification rate response to pH and Omega."""
    from model_calcification import calcification_rate_IpCC, calcification_pH_response
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    Omega = np.linspace(0.5, 5.0, 200)
    pH = np.linspace(7.2, 8.4, 200)
    
    for sp, color in SPECIES_COLORS.items():
        G = calcification_rate_IpCC(Omega, sp)
        axes[0].plot(Omega, G, color=color, label=sp, linewidth=2)
        
        G_pH = calcification_pH_response(pH, sp)
        axes[1].plot(pH, G_pH, color=color, label=sp, linewidth=2)
    
    axes[0].axvline(x=3.5, color='green', linestyle='--', alpha=0.5, label='Present Ω')
    axes[0].axvline(x=1.8, color='red', linestyle='--', alpha=0.5, label='RCP8.5 2100 Ω')
    axes[0].axvline(x=1.0, color='black', linestyle=':', alpha=0.5, label='Dissolution')
    axes[0].set_xlabel('Ω aragonite')
    axes[0].set_ylabel('Calcification Rate (µmol CaCO₃ cm⁻² h⁻¹)')
    axes[0].set_title('(a) Calcification vs. Aragonite Saturation')
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].axvline(x=8.1, color='green', linestyle='--', alpha=0.5, label='Present pH')
    axes[1].axvline(x=7.7, color='red', linestyle='--', alpha=0.5, label='RCP8.5 2100 pH')
    axes[1].set_xlabel('pH')
    axes[1].set_ylabel('Relative Calcification')
    axes[1].set_title('(b) Calcification Response to pH')
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)
    
    fig.suptitle('Coral Calcification Rate Models', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig('figures/fig2_calcification.png')
    fig.savefig('figures/fig2_calcification.svg')
    plt.close()
    print("  Fig 2 saved")


def fig3_network():
    """Figure 3: Species interaction network."""
    import networkx as nx
    from model_network import build_interaction_matrix, build_network_graph, SPECIES, N_SPECIES
    
    alpha = build_interaction_matrix()
    G = build_network_graph(alpha)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Network graph
    type_colors = {
        'coral': '#e41a1c', 'CCA': '#a65628', 'algae': '#4daf4a',
        'symbiont': '#ff7f00', 'herbivore': '#377eb8',
        'predator': '#984ea3', 'fish': '#999999'
    }
    node_colors = [type_colors.get(SPECIES[n]['type'], '#666') for n in G.nodes()]
    
    pos = nx.spring_layout(G, seed=42, k=2.5)
    
    pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] > 0]
    neg_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] < 0]
    
    nx.draw_networkx_edges(G, pos, edgelist=pos_edges, edge_color='green',
                           alpha=0.4, width=1.5, ax=axes[0], arrows=True,
                           connectionstyle='arc3,rad=0.1')
    nx.draw_networkx_edges(G, pos, edgelist=neg_edges, edge_color='red',
                           alpha=0.4, width=1.5, ax=axes[0], arrows=True,
                           connectionstyle='arc3,rad=0.1')
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=600,
                           alpha=0.9, ax=axes[0])
    labels = {i: SPECIES[i]['name'].split('(')[0].strip()[:12] for i in range(N_SPECIES)}
    nx.draw_networkx_labels(G, pos, labels, font_size=7, ax=axes[0])
    
    axes[0].set_title('(a) Reef Community Interaction Network')
    # Legend
    for type_name, color in type_colors.items():
        axes[0].scatter([], [], c=color, s=80, label=type_name)
    axes[0].legend(fontsize=7, loc='lower left')
    axes[0].axis('off')
    
    # Interaction matrix heatmap
    species_names = [SPECIES[i]['name'].split('(')[0].strip()[:10] for i in range(N_SPECIES)]
    im = axes[1].imshow(alpha, cmap='RdBu_r', vmin=-0.6, vmax=0.6, aspect='auto')
    axes[1].set_xticks(range(N_SPECIES))
    axes[1].set_xticklabels(species_names, rotation=45, ha='right', fontsize=7)
    axes[1].set_yticks(range(N_SPECIES))
    axes[1].set_yticklabels(species_names, fontsize=7)
    axes[1].set_title('(b) Interaction Matrix (α)')
    plt.colorbar(im, ax=axes[1], shrink=0.8, label='Interaction strength')
    
    fig.suptitle('Species Interaction Network Model', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig('figures/fig3_network.png')
    fig.savefig('figures/fig3_network.svg')
    plt.close()
    print("  Fig 3 saved")


def fig4_compound_stress():
    """Figure 4: Temperature-pH compound stress heatmaps."""
    from model_compound_stress import compound_stress_model
    
    T_range = np.linspace(22, 34, 60)
    pH_range = np.linspace(7.2, 8.4, 60)
    T_grid, pH_grid = np.meshgrid(T_range, pH_range)
    
    species_params = {
        'Acropora': {'T_opt': 27.0, 'synergy': 1.8},
        'Porites': {'T_opt': 27.5, 'synergy': 1.2},
        'Montipora': {'T_opt': 27.0, 'synergy': 1.5},
        'Stylophora': {'T_opt': 26.5, 'synergy': 1.6},
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    
    for idx, (sp_name, params) in enumerate(species_params.items()):
        ax = axes[idx // 2, idx % 2]
        perf, _ = compound_stress_model(
            T_grid, pH_grid, T_opt=params['T_opt'],
            synergy_factor=params['synergy'], model_type='synergistic'
        )
        
        im = ax.contourf(T_range, pH_range, perf, levels=20, cmap='viridis')
        ax.contour(T_range, pH_range, perf, levels=[0.3], colors='red', linewidths=2)
        ax.contour(T_range, pH_range, perf, levels=[0.5], colors='orange', linewidths=1.5)
        
        # Mark scenarios
        scenarios = {'2020': (25.5, 8.07), '2050': (27.0, 7.90), '2100': (29.2, 7.70)}
        for label, (t, ph) in scenarios.items():
            ax.plot(t, ph, 'w*', markersize=10, markeredgecolor='black')
            ax.annotate(label, (t, ph), fontsize=7, color='white',
                       fontweight='bold', ha='center', va='bottom',
                       xytext=(0, 5), textcoords='offset points')
        
        plt.colorbar(im, ax=ax, shrink=0.8, label='Performance')
        ax.set_xlabel('Temperature (°C)')
        ax.set_ylabel('pH')
        ax.set_title(f'({chr(97+idx)}) {sp_name}')
    
    fig.suptitle('Compound Temperature-pH Stress (Synergistic Model)\n'
                 'Red contour = 30% performance threshold', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig('figures/fig4_compound_stress.png')
    fig.savefig('figures/fig4_compound_stress.svg')
    plt.close()
    print("  Fig 4 saved")


def fig5_popgen():
    """Figure 5: Population genetics trajectories."""
    data = np.load('results/popgen_trajectories.npz')
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    species_list = ['Acropora', 'Porites', 'Stylophora']
    sp_colors = {'Acropora': '#e41a1c', 'Porites': '#377eb8', 'Stylophora': '#984ea3'}
    
    for scenario, ls in [('RCP26', '-'), ('RCP45', '--'), ('RCP85', ':')]:
        for sp in species_list:
            key_freq = f"{sp}_{scenario}_freq"
            key_fit = f"{sp}_{scenario}_fitness"
            key_va = f"{sp}_{scenario}_Va"
            
            if key_freq in data:
                freq = data[key_freq]
                gens = np.arange(len(freq))
                axes[0].plot(gens, freq, color=sp_colors[sp], linestyle=ls,
                           linewidth=1.5, alpha=0.8)
            
            if key_fit in data:
                fit = data[key_fit]
                gens = np.arange(len(fit))
                axes[1].plot(gens, fit, color=sp_colors[sp], linestyle=ls,
                           linewidth=1.5, alpha=0.8)
            
            if key_va in data:
                va = data[key_va]
                gens = np.arange(len(va))
                axes[2].plot(gens, va, color=sp_colors[sp], linestyle=ls,
                           linewidth=1.5, alpha=0.8)
    
    axes[0].set_xlabel('Generation')
    axes[0].set_ylabel('Mean Tolerance Allele Frequency')
    axes[0].set_title('(a) Allele Frequency Trajectories')
    
    axes[1].set_xlabel('Generation')
    axes[1].set_ylabel('Mean Population Fitness')
    axes[1].set_title('(b) Fitness Trajectories')
    
    axes[2].set_xlabel('Generation')
    axes[2].set_ylabel('Additive Genetic Variance (Va)')
    axes[2].set_title('(c) Genetic Variance')
    
    # Custom legend
    from matplotlib.lines import Line2D
    sp_handles = [Line2D([0], [0], color=c, linewidth=2) for c in sp_colors.values()]
    scen_handles = [Line2D([0], [0], color='gray', linestyle=ls, linewidth=2)
                    for ls in ['-', '--', ':']]
    
    axes[0].legend(sp_handles + scen_handles,
                   list(sp_colors.keys()) + ['RCP2.6', 'RCP4.5', 'RCP8.5'],
                   fontsize=7, ncol=2)
    
    for ax in axes:
        ax.grid(True, alpha=0.3)
    
    fig.suptitle('Evolutionary Response to Ocean Acidification', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig('figures/fig5_popgen.png')
    fig.savefig('figures/fig5_popgen.svg')
    plt.close()
    print("  Fig 5 saved")


def fig6_gbr_projection():
    """Figure 6: GBR 2100 projection."""
    with open('results/gbr_projection.json') as f:
        data = json.load(f)
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    
    regions = list(data['projection_timelines'].keys())
    
    for col, region in enumerate(regions):
        for scenario in ['RCP2.6', 'RCP4.5', 'RCP8.5']:
            timeline = data['projection_timelines'][region][scenario]
            years = [e['year'] for e in timeline]
            covers = [e['coral_cover'] * 100 for e in timeline]
            temps = [e['T_C'] for e in timeline]
            
            axes[0, col].plot(years, covers, color=COLORS[scenario],
                            label=scenario, linewidth=2, marker='o', markersize=3)
            axes[1, col].plot(years, temps, color=COLORS[scenario],
                            label=scenario, linewidth=2)
        
        axes[0, col].axhline(y=10, color='red', linestyle='--', alpha=0.5,
                            label='Functional collapse')
        axes[0, col].set_title(region)
        axes[0, col].set_ylabel('Coral Cover (%)')
        axes[0, col].set_ylim(0, 45)
        axes[0, col].legend(fontsize=7)
        axes[0, col].grid(True, alpha=0.3)
        
        axes[1, col].set_xlabel('Year')
        axes[1, col].set_ylabel('SST (°C)')
        axes[1, col].axhline(y=29, color='red', linestyle=':', alpha=0.5,
                            label='Bleaching threshold')
        axes[1, col].legend(fontsize=7)
        axes[1, col].grid(True, alpha=0.3)
    
    axes[0, 0].set_ylabel('Coral Cover (%)')
    axes[1, 0].set_ylabel('SST (°C)')
    
    fig.suptitle('Great Barrier Reef 2100 Projections', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig('figures/fig6_gbr_projection.png')
    fig.savefig('figures/fig6_gbr_projection.svg')
    plt.close()
    print("  Fig 6 saved")


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    print("Generating figures...")
    fig1_carbonate_chemistry()
    fig2_calcification()
    fig3_network()
    fig4_compound_stress()
    fig5_popgen()
    fig6_gbr_projection()
    print("All figures saved to figures/")
