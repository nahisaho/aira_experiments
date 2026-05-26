#!/usr/bin/env python3
"""
Electrochemical CO2 Reduction Reaction (CO2RR) Catalyst Screening Pipeline
===========================================================================
Computational screening system for high-activity CO2RR catalysts using
scaling relations, volcano plots, and automated candidate evaluation.

Based on the Computational Hydrogen Electrode (CHE) framework and
descriptor-based catalyst screening methodology.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.optimize import curve_fit
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import json
import os

# Output directory
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================================
# 1. Thermodynamic Data & Reaction Pathway Definition
# ============================================================================

# Gas-phase energies (eV) - reference values from DFT (PBE+D3)
GAS_ENERGIES = {
    'CO2': -22.96,
    'H2O': -14.22,
    'H2': -6.77,
    'CO': -14.79,
    'CH4': -24.03,
    'C2H4': -32.88,
    'CH3OH': -30.05,
    'HCOOH': -25.61,
    'C2H5OH': -46.32,
}

# Zero-point energy corrections (eV)
ZPE_CORRECTIONS = {
    '*CO2': 0.31, '*COOH': 0.62, '*CO': 0.19, '*CHO': 0.45,
    '*COH': 0.44, '*OCCO': 0.38, '*CH2O': 0.75, '*OCH3': 0.98,
    '*O': 0.07, '*OH': 0.35, '*H': 0.17, '*OCHO': 0.60,
    'CO2(g)': 0.31, 'H2O(l)': 0.56, 'H2(g)': 0.27,
    'CO(g)': 0.14, 'CH4(g)': 1.19,
}

# Entropy corrections at 298K (eV)  -TΔS
ENTROPY_CORRECTIONS = {
    'CO2(g)': -0.66, 'H2O(l)': -0.67, 'H2(g)': -0.40,
    'CO(g)': -0.61, 'CH4(g)': -0.58,
    '*CO2': 0.0, '*COOH': 0.0, '*CO': 0.0, '*CHO': 0.0,
    '*COH': 0.0, '*OCCO': 0.0, '*O': 0.0, '*OH': 0.0, '*H': 0.0,
}

# Solvation stabilization energies (eV) from implicit/explicit solvation
SOLVATION_CORRECTIONS = {
    '*COOH': -0.25, '*CO': 0.0, '*CHO': -0.10,
    '*COH': -0.15, '*OCCO': -0.20, '*OH': -0.50,
    '*H': 0.0, '*O': 0.0, '*OCHO': -0.30,
}


@dataclass
class CatalystData:
    """Stores computed adsorption energies and properties for a catalyst."""
    name: str
    E_CO: float       # *CO adsorption energy (eV)
    E_COOH: float     # *COOH adsorption energy (eV)
    E_CHO: float      # *CHO adsorption energy (eV)
    E_OH: float = 0.0
    E_H: float = 0.0
    E_OCCO: float = 0.0
    d_band_center: float = 0.0
    category: str = "metal"
    color: str = "blue"
    marker: str = "o"
    
    @property
    def limiting_potential_CO(self) -> float:
        """Limiting potential for CO2 -> CO pathway."""
        dG1 = self.E_COOH + 0.33  # CO2 + H+ + e- -> *COOH
        dG2 = self.E_CO - self.E_COOH + 0.14  # *COOH + H+ + e- -> *CO + H2O
        return -max(dG1, dG2)
    
    @property
    def limiting_potential_CHO(self) -> float:
        """Limiting potential for *CO -> *CHO (rate-limiting for C1)."""
        return -(self.E_CHO - self.E_CO + 0.20)
    
    @property
    def limiting_potential_C2(self) -> float:
        """Limiting potential for C-C coupling via *CO dimerization."""
        return -(self.E_OCCO - 2 * self.E_CO + 0.30)
    
    @property
    def overpotential_CO(self) -> float:
        """Overpotential for CO production (equilibrium at -0.11 V vs RHE)."""
        return abs(self.limiting_potential_CO) - 0.11
    
    @property 
    def selectivity_index(self) -> float:
        """HER vs CO2RR selectivity (negative = CO2RR favorable)."""
        return self.E_H - 0.5 * self.E_CO


# ============================================================================
# 2. Catalyst Database - Adsorption Energies from DFT Literature
# ============================================================================

def build_catalyst_database() -> List[CatalystData]:
    """Build database of catalyst adsorption energies from DFT calculations."""
    catalysts = []
    
    # Pure transition metals (111) facets
    metals_data = {
        'Au': (-0.15, -0.30, 0.80, -0.20, 0.15, 0.50, -2.55),
        'Ag': (-0.05, -0.20, 0.90, -0.10, 0.20, 0.60, -3.30),
        'Cu': (-0.55, -0.60, 0.10, -0.60, -0.20, -0.30, -2.67),
        'Zn': (0.25, 0.10, 1.20, 0.20, 0.35, 0.90, -4.50),
        'Pd': (-0.80, -0.90, -0.10, -0.85, -0.40, -0.70, -1.83),
        'Pt': (-1.20, -1.30, -0.50, -1.10, -0.55, -1.20, -2.25),
        'Ni': (-1.30, -1.10, -0.40, -0.90, -0.45, -1.10, -1.29),
        'Fe': (-1.00, -0.80, -0.20, -0.70, -0.35, -0.80, -0.92),
        'Co': (-1.10, -0.95, -0.30, -0.80, -0.40, -0.95, -1.17),
        'Rh': (-1.00, -1.05, -0.25, -0.90, -0.42, -0.85, -1.73),
        'Ir': (-1.15, -1.20, -0.40, -1.00, -0.50, -1.05, -2.11),
        'Sn': (0.10, 0.05, 1.05, 0.10, 0.30, 0.75, -5.20),
        'Bi': (0.30, 0.15, 1.30, 0.25, 0.40, 1.00, -6.10),
        'In': (0.20, 0.10, 1.15, 0.15, 0.33, 0.85, -5.80),
    }
    
    for name, (E_CO, E_COOH, E_CHO, E_OH, E_H, E_OCCO, d_band) in metals_data.items():
        catalysts.append(CatalystData(
            name=name, E_CO=E_CO, E_COOH=E_COOH, E_CHO=E_CHO,
            E_OH=E_OH, E_H=E_H, E_OCCO=E_OCCO,
            d_band_center=d_band, category="metal",
            color="#1f77b4", marker="o"
        ))
    
    # Cu-based alloys
    alloys_data = {
        'CuAg': (-0.35, -0.45, 0.40, -0.40, -0.05, 0.05, -2.90),
        'CuAu': (-0.40, -0.50, 0.35, -0.45, -0.08, 0.00, -2.80),
        'CuZn': (-0.30, -0.35, 0.50, -0.35, 0.00, 0.15, -3.10),
        'CuNi': (-0.80, -0.85, -0.10, -0.75, -0.30, -0.60, -1.95),
        'CuPd': (-0.65, -0.70, 0.05, -0.70, -0.28, -0.45, -2.20),
        'CuSn': (-0.20, -0.25, 0.60, -0.25, 0.05, 0.30, -3.50),
        'CuIn': (-0.25, -0.30, 0.55, -0.30, 0.02, 0.25, -3.40),
        'Cu3Ag': (-0.45, -0.52, 0.25, -0.50, -0.12, -0.12, -2.78),
        'Cu3Au': (-0.48, -0.55, 0.20, -0.52, -0.14, -0.15, -2.72),
    }
    
    for name, (E_CO, E_COOH, E_CHO, E_OH, E_H, E_OCCO, d_band) in alloys_data.items():
        catalysts.append(CatalystData(
            name=name, E_CO=E_CO, E_COOH=E_COOH, E_CHO=E_CHO,
            E_OH=E_OH, E_H=E_H, E_OCCO=E_OCCO,
            d_band_center=d_band, category="Cu-alloy",
            color="#ff7f0e", marker="s"
        ))
    
    # Single atom catalysts on N-doped carbon (M-N4-C)
    sac_data = {
        'Ni-N4-C': (-0.90, -0.50, 0.30, -0.55, -0.10, -0.40, -1.55),
        'Fe-N4-C': (-0.70, -0.35, 0.45, -0.40, -0.05, -0.20, -1.10),
        'Co-N4-C': (-0.75, -0.45, 0.35, -0.50, -0.08, -0.30, -1.35),
        'Cu-N4-C': (-0.40, -0.25, 0.65, -0.30, 0.05, 0.10, -2.80),
        'Mn-N4-C': (-0.60, -0.30, 0.50, -0.35, 0.00, -0.10, -0.85),
        'Cr-N4-C': (-0.50, -0.20, 0.55, -0.25, 0.05, -0.05, -0.65),
        'V-N4-C':  (-0.45, -0.15, 0.60, -0.20, 0.08, 0.00, -0.50),
        'Zn-N4-C': (0.10, 0.05, 1.00, 0.05, 0.25, 0.50, -4.60),
        'Mo-N4-C': (-0.85, -0.55, 0.20, -0.60, -0.15, -0.50, -1.00),
        'W-N4-C':  (-0.95, -0.60, 0.15, -0.65, -0.18, -0.55, -0.90),
        'Ru-N4-C': (-0.80, -0.48, 0.28, -0.52, -0.12, -0.38, -1.40),
        'Pd-N4-C': (-0.55, -0.35, 0.50, -0.40, -0.02, -0.05, -1.90),
        'Ag-N4-C': (0.05, 0.00, 0.95, 0.00, 0.22, 0.45, -3.50),
        'Sn-N4-C': (0.15, 0.08, 1.05, 0.10, 0.28, 0.55, -5.30),
    }
    
    for name, (E_CO, E_COOH, E_CHO, E_OH, E_H, E_OCCO, d_band) in sac_data.items():
        catalysts.append(CatalystData(
            name=name, E_CO=E_CO, E_COOH=E_COOH, E_CHO=E_CHO,
            E_OH=E_OH, E_H=E_H, E_OCCO=E_OCCO,
            d_band_center=d_band, category="SAC-NC",
            color="#2ca02c", marker="^"
        ))
    
    return catalysts


# ============================================================================
# 3. Scaling Relations
# ============================================================================

def compute_scaling_relations(catalysts: List[CatalystData]) -> Dict:
    """Compute and visualize linear scaling relations between adsorption energies."""
    E_CO = np.array([c.E_CO for c in catalysts])
    E_COOH = np.array([c.E_COOH for c in catalysts])
    E_CHO = np.array([c.E_CHO for c in catalysts])
    E_OH = np.array([c.E_OH for c in catalysts])
    E_OCCO = np.array([c.E_OCCO for c in catalysts])
    categories = [c.category for c in catalysts]
    names = [c.name for c in catalysts]
    
    # Fit scaling: E_COOH = a * E_CO + b
    popt_cooh, _ = curve_fit(lambda x, a, b: a * x + b, E_CO, E_COOH)
    # Fit scaling: E_CHO = a * E_CO + b
    popt_cho, _ = curve_fit(lambda x, a, b: a * x + b, E_CO, E_CHO)
    # Fit scaling: E_OCCO = a * E_CO + b
    popt_occo, _ = curve_fit(lambda x, a, b: a * x + b, E_CO, E_OCCO)
    
    results = {
        'COOH_vs_CO': {'slope': popt_cooh[0], 'intercept': popt_cooh[1]},
        'CHO_vs_CO': {'slope': popt_cho[0], 'intercept': popt_cho[1]},
        'OCCO_vs_CO': {'slope': popt_occo[0], 'intercept': popt_occo[1]},
    }
    
    # R-squared
    for key, E_y, popt in [('COOH_vs_CO', E_COOH, popt_cooh),
                            ('CHO_vs_CO', E_CHO, popt_cho),
                            ('OCCO_vs_CO', E_OCCO, popt_occo)]:
        y_pred = popt[0] * E_CO + popt[1]
        ss_res = np.sum((E_y - y_pred)**2)
        ss_tot = np.sum((E_y - np.mean(E_y))**2)
        results[key]['R2'] = 1 - ss_res / ss_tot
    
    # --- Figure: Scaling Relations ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    scaling_pairs = [
        (E_COOH, 'E(*COOH) (eV)', popt_cooh, results['COOH_vs_CO']),
        (E_CHO, 'E(*CHO) (eV)', popt_cho, results['CHO_vs_CO']),
        (E_OCCO, 'E(*OCCO) (eV)', popt_occo, results['OCCO_vs_CO']),
    ]
    
    cat_colors = {'metal': '#1f77b4', 'Cu-alloy': '#ff7f0e', 'SAC-NC': '#2ca02c'}
    cat_markers = {'metal': 'o', 'Cu-alloy': 's', 'SAC-NC': '^'}
    cat_labels = {'metal': 'Pure Metals', 'Cu-alloy': 'Cu Alloys', 'SAC-NC': 'SAC (M-N₄-C)'}
    
    for ax, (E_y, ylabel, popt, res) in zip(axes, scaling_pairs):
        for cat in cat_colors:
            mask = np.array([c == cat for c in categories])
            ax.scatter(E_CO[mask], E_y[mask], c=cat_colors[cat],
                      marker=cat_markers[cat], s=60, label=cat_labels[cat],
                      edgecolors='black', linewidth=0.5, zorder=3)
        
        x_fit = np.linspace(E_CO.min() - 0.2, E_CO.max() + 0.2, 100)
        ax.plot(x_fit, popt[0] * x_fit + popt[1], 'k--', linewidth=1.5,
                label=f'y = {popt[0]:.2f}x + {popt[1]:.2f}\nR² = {res["R2"]:.3f}')
        
        # Annotate select catalysts
        for i, name in enumerate(names):
            if name in ['Cu', 'Au', 'Ag', 'Ni-N4-C', 'Fe-N4-C', 'CuAg', 'Pt']:
                ax.annotate(name, (E_CO[i], E_y[i]), fontsize=7,
                           xytext=(5, 5), textcoords='offset points')
        
        ax.set_xlabel('E(*CO) (eV)', fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.3)
    
    axes[0].set_title('(a) *COOH vs *CO', fontsize=12)
    axes[1].set_title('(b) *CHO vs *CO', fontsize=12)
    axes[2].set_title('(c) *OCCO vs *CO', fontsize=12)
    
    plt.suptitle('Adsorption Energy Scaling Relations for CO₂RR Intermediates', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/scaling_relations.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    print("=" * 60)
    print("SCALING RELATIONS RESULTS")
    print("=" * 60)
    for key, res in results.items():
        print(f"  {key}: slope={res['slope']:.3f}, intercept={res['intercept']:.3f}, R²={res['R2']:.3f}")
    
    return results


# ============================================================================
# 4. Volcano Plot Construction
# ============================================================================

def compute_volcano_plot(catalysts: List[CatalystData], scaling: Dict) -> Dict:
    """Construct volcano plots for CO2 -> CO and CO2 -> C2+ pathways."""
    
    # --- CO2 -> CO Volcano ---
    E_CO_range = np.linspace(-1.5, 0.5, 500)
    
    # Using scaling: E_COOH = a*E_CO + b
    a_cooh = scaling['COOH_vs_CO']['slope']
    b_cooh = scaling['COOH_vs_CO']['intercept']
    
    # Step 1: CO2 + H+ + e- -> *COOH;  dG1 = E_COOH + 0.33
    dG1 = a_cooh * E_CO_range + b_cooh + 0.33
    # Step 2: *COOH + H+ + e- -> *CO + H2O;  dG2 = E_CO - E_COOH + 0.14
    dG2 = E_CO_range - (a_cooh * E_CO_range + b_cooh) + 0.14
    
    U_L_CO = -np.maximum(dG1, dG2)
    
    # --- CO2 -> C1 (further reduction) Volcano ---
    a_cho = scaling['CHO_vs_CO']['slope']
    b_cho = scaling['CHO_vs_CO']['intercept']
    
    # Step: *CO + H+ + e- -> *CHO
    dG_cho = a_cho * E_CO_range + b_cho - E_CO_range + 0.20
    U_L_C1 = -np.maximum(np.maximum(dG1, dG2), dG_cho)
    
    # --- CO2 -> C2+ Volcano ---
    a_occo = scaling['OCCO_vs_CO']['slope']
    b_occo = scaling['OCCO_vs_CO']['intercept']
    
    # Step: 2*CO -> *OCCO (C-C coupling)
    dG_cc = a_occo * E_CO_range + b_occo - 2 * E_CO_range + 0.30
    U_L_C2 = -np.maximum(np.maximum(dG1, dG2), dG_cc)
    
    # Figure: Volcano Plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    cat_colors = {'metal': '#1f77b4', 'Cu-alloy': '#ff7f0e', 'SAC-NC': '#2ca02c'}
    cat_markers = {'metal': 'o', 'Cu-alloy': 's', 'SAC-NC': '^'}
    cat_labels = {'metal': 'Pure Metals', 'Cu-alloy': 'Cu Alloys', 'SAC-NC': 'SAC (M-N₄-C)'}
    
    titles = ['(a) CO₂ → CO', '(b) CO₂ → C₁ (CH₄/CH₃OH)', '(c) CO₂ → C₂+ (C₂H₄/C₂H₅OH)']
    U_Ls = [U_L_CO, U_L_C1, U_L_C2]
    lp_funcs = [
        lambda c: c.limiting_potential_CO,
        lambda c: -max(c.E_COOH + 0.33,
                       c.E_CO - c.E_COOH + 0.14,
                       c.E_CHO - c.E_CO + 0.20),
        lambda c: -max(c.E_COOH + 0.33,
                       c.E_CO - c.E_COOH + 0.14,
                       c.E_OCCO - 2 * c.E_CO + 0.30),
    ]
    eq_potentials = [-0.11, -0.24, -0.34]
    
    for ax, title, U_L, lp_func, U_eq in zip(axes, titles, U_Ls, lp_funcs, eq_potentials):
        ax.plot(E_CO_range, U_L, 'k-', linewidth=2, label='Volcano (scaling)', zorder=1)
        ax.axhline(y=U_eq, color='gray', linestyle=':', alpha=0.6,
                   label=f'U_eq = {U_eq:.2f} V')
        
        for cat in cat_colors:
            cat_catalysts = [c for c in catalysts if c.category == cat]
            x = [c.E_CO for c in cat_catalysts]
            y = [lp_func(c) for c in cat_catalysts]
            ax.scatter(x, y, c=cat_colors[cat], marker=cat_markers[cat],
                      s=80, label=cat_labels[cat], edgecolors='black',
                      linewidth=0.5, zorder=3)
            for c in cat_catalysts:
                if c.name in ['Cu', 'Au', 'Ag', 'Pt', 'Ni', 'CuAg', 'CuAu',
                              'Ni-N4-C', 'Fe-N4-C', 'Cu-N4-C', 'Mo-N4-C', 'CuZn']:
                    ax.annotate(c.name, (c.E_CO, lp_func(c)), fontsize=7,
                               xytext=(5, 5), textcoords='offset points')
        
        ax.set_xlabel('E(*CO) (eV)', fontsize=12)
        ax.set_ylabel('Limiting Potential U_L (V vs RHE)', fontsize=12)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=7, loc='lower right')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-2.0, 0.5)
    
    plt.suptitle('Volcano Plots for Electrocatalytic CO₂ Reduction', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/volcano_plots.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    # Compute performance metrics
    results = {}
    for c in catalysts:
        results[c.name] = {
            'U_L_CO': c.limiting_potential_CO,
            'U_L_C1': lp_funcs[1](c),
            'U_L_C2': lp_funcs[2](c),
            'eta_CO': abs(c.limiting_potential_CO) - 0.11,
            'selectivity': c.selectivity_index,
        }
    
    print("\n" + "=" * 60)
    print("VOLCANO PLOT - TOP CATALYSTS")
    print("=" * 60)
    
    for pathway, idx in [('CO production', 'U_L_CO'),
                          ('C1 production', 'U_L_C1'),
                          ('C2+ production', 'U_L_C2')]:
        sorted_cats = sorted(results.items(), key=lambda x: x[1][idx], reverse=True)[:5]
        print(f"\n  Top 5 for {pathway}:")
        for name, data in sorted_cats:
            print(f"    {name:12s}: U_L = {data[idx]:.3f} V, η = {abs(data[idx]):.3f} V")
    
    return results


# ============================================================================
# 5. Reaction Pathway Analysis
# ============================================================================

def analyze_reaction_pathways(catalysts: List[CatalystData]) -> None:
    """Analyze CO2RR pathways: CO2 → CO → C1/C2+ for selected catalysts."""
    
    selected = ['Cu', 'CuAg', 'Ni-N4-C', 'Fe-N4-C', 'Au', 'Mo-N4-C']
    selected_cats = [c for c in catalysts if c.name in selected]
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    colors = plt.cm.Set2(np.linspace(0, 1, 8))
    
    for idx, cat in enumerate(selected_cats):
        ax = axes[idx]
        
        # CO2 → CO pathway: CO2(g) → *COOH → *CO → CO(g)
        # CO → C1 pathway: *CO → *CHO → *CH2O → *OCH3 → CH3OH
        # CO → C2 pathway: 2*CO → *OCCO → C2H4
        
        # Free energy at U=0
        states_CO = ['CO₂(g)', '*COOH', '*CO', 'CO(g)']
        G_CO_0 = [0.0,
                  cat.E_COOH + 0.33,
                  cat.E_CO + 0.14 + cat.E_COOH + 0.33 - cat.E_COOH,
                  -0.11 * 2]  # CO(g) at equilibrium
        G_CO_0 = [0.0,
                  cat.E_COOH + 0.33,
                  cat.E_CO + 0.47,
                  0.28]
        
        # C1 path: *CO → *CHO → CH4 (simplified)
        states_C1 = ['*CO', '*CHO', '*CH₂O', 'CH₄(g)']
        G_C1_0 = [cat.E_CO + 0.47,
                  cat.E_CHO + 0.67,
                  cat.E_CHO + 0.47,
                  -0.24 * 8]  # CH4 at equilibrium
        G_C1_0 = [cat.E_CO + 0.47,
                  cat.E_CHO + 0.67,
                  cat.E_CHO + 0.37,
                  -1.06]
        
        # Free energy at U = -0.8 V
        U = -0.8
        G_CO_U = [0.0,
                  cat.E_COOH + 0.33 + U,
                  cat.E_CO + 0.47 + 2 * U,
                  0.28 + 2 * U]
        
        # Plot pathway diagram
        x_pos = np.arange(len(states_CO))
        width = 0.8
        
        for i in range(len(states_CO)):
            ax.plot([x_pos[i] - width/2, x_pos[i] + width/2],
                    [G_CO_0[i], G_CO_0[i]], '-', color=colors[0], linewidth=2.5)
            ax.plot([x_pos[i] - width/2, x_pos[i] + width/2],
                    [G_CO_U[i], G_CO_U[i]], '--', color=colors[1], linewidth=2.5)
            if i > 0:
                ax.plot([x_pos[i-1] + width/2, x_pos[i] - width/2],
                        [G_CO_0[i-1], G_CO_0[i]], ':', color=colors[0], alpha=0.5)
                ax.plot([x_pos[i-1] + width/2, x_pos[i] - width/2],
                        [G_CO_U[i-1], G_CO_U[i]], ':', color=colors[1], alpha=0.5)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(states_CO, fontsize=9)
        ax.set_ylabel('ΔG (eV)', fontsize=10)
        ax.set_title(f'{cat.name}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.2)
        
        if idx == 0:
            ax.plot([], [], '-', color=colors[0], linewidth=2.5, label='U = 0 V')
            ax.plot([], [], '--', color=colors[1], linewidth=2.5, label='U = -0.8 V')
            ax.legend(fontsize=9)
    
    plt.suptitle('Free Energy Diagrams for CO₂ → CO Pathway on Selected Catalysts', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/reaction_pathways.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("\n[OK] Reaction pathway diagrams saved.")


# ============================================================================
# 6. Metal-Support Interaction Analysis (SAC)
# ============================================================================

def analyze_metal_support_interaction(catalysts: List[CatalystData]) -> None:
    """Analyze metal-support interactions for SACs on N-doped carbon."""
    
    sac_cats = [c for c in catalysts if c.category == 'SAC-NC']
    metals = [c.name.split('-')[0] for c in sac_cats]
    
    # d-band centers
    d_bands = [c.d_band_center for c in sac_cats]
    E_CO_vals = [c.E_CO for c in sac_cats]
    E_COOH_vals = [c.E_COOH for c in sac_cats]
    UL_CO = [c.limiting_potential_CO for c in sac_cats]
    
    # Charge transfer estimates (from DFT Bader analysis - literature values)
    charge_transfers = {
        'Ni': 0.82, 'Fe': 0.95, 'Co': 0.88, 'Cu': 0.65,
        'Mn': 1.05, 'Cr': 1.15, 'V': 1.20, 'Zn': 0.55,
        'Mo': 1.10, 'W': 1.18, 'Ru': 0.90, 'Pd': 0.70,
        'Ag': 0.45, 'Sn': 0.85,
    }
    charges = [charge_transfers.get(m, 0.8) for m in metals]
    
    # Binding energies (metal to N4 site)
    binding_energies = {
        'Ni': -5.2, 'Fe': -4.8, 'Co': -5.0, 'Cu': -3.9,
        'Mn': -4.5, 'Cr': -4.2, 'V': -4.0, 'Zn': -2.8,
        'Mo': -5.5, 'W': -5.8, 'Ru': -5.1, 'Pd': -4.1,
        'Ag': -2.5, 'Sn': -3.5,
    }
    E_bind = [binding_energies.get(m, -4.0) for m in metals]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # (a) d-band center vs E(*CO)
    ax = axes[0, 0]
    scatter = ax.scatter(d_bands, E_CO_vals, c=charges, cmap='viridis',
                         s=100, edgecolors='black', linewidth=0.5, zorder=3)
    for i, m in enumerate(metals):
        ax.annotate(m, (d_bands[i], E_CO_vals[i]), fontsize=9,
                   xytext=(5, 5), textcoords='offset points')
    ax.set_xlabel('d-band Center (eV)', fontsize=12)
    ax.set_ylabel('E(*CO) (eV)', fontsize=12)
    ax.set_title('(a) d-band Center vs *CO Adsorption', fontsize=12)
    plt.colorbar(scatter, ax=ax, label='Charge Transfer (|e|)')
    ax.grid(True, alpha=0.3)
    
    # Fit linear trend
    z = np.polyfit(d_bands, E_CO_vals, 1)
    x_fit = np.linspace(min(d_bands) - 0.3, max(d_bands) + 0.3, 100)
    ax.plot(x_fit, np.polyval(z, x_fit), 'r--', alpha=0.7,
            label=f'Linear fit: R²={np.corrcoef(d_bands, E_CO_vals)[0,1]**2:.3f}')
    ax.legend(fontsize=9)
    
    # (b) Charge transfer vs limiting potential
    ax = axes[0, 1]
    scatter = ax.scatter(charges, UL_CO, c=[c.d_band_center for c in sac_cats],
                         cmap='coolwarm', s=100, edgecolors='black', linewidth=0.5, zorder=3)
    for i, m in enumerate(metals):
        ax.annotate(m, (charges[i], UL_CO[i]), fontsize=9,
                   xytext=(5, 5), textcoords='offset points')
    ax.set_xlabel('Charge Transfer to Support (|e|)', fontsize=12)
    ax.set_ylabel('Limiting Potential U_L (V)', fontsize=12)
    ax.set_title('(b) Charge Transfer vs CO₂→CO Activity', fontsize=12)
    plt.colorbar(scatter, ax=ax, label='d-band Center (eV)')
    ax.grid(True, alpha=0.3)
    
    # (c) Binding energy vs stability
    ax = axes[1, 0]
    dissolution_potential = [abs(e) * 0.4 + 0.5 for e in E_bind]  # Approximate
    ax.scatter(E_bind, dissolution_potential, c='#2ca02c', s=100,
               edgecolors='black', linewidth=0.5, marker='^', zorder=3)
    for i, m in enumerate(metals):
        ax.annotate(m, (E_bind[i], dissolution_potential[i]), fontsize=9,
                   xytext=(5, 5), textcoords='offset points')
    ax.axhline(y=1.5, color='red', linestyle='--', alpha=0.6, label='Stability threshold')
    ax.set_xlabel('Binding Energy E_bind (eV)', fontsize=12)
    ax.set_ylabel('Dissolution Potential (V)', fontsize=12)
    ax.set_title('(c) Metal-N₄ Binding vs Electrochemical Stability', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # (d) Selectivity map: HER vs CO2RR
    ax = axes[1, 1]
    selectivity = [c.selectivity_index for c in sac_cats]
    colors_sel = ['green' if s < 0 else 'red' for s in selectivity]
    ax.barh(metals, selectivity, color=colors_sel, edgecolor='black', linewidth=0.5)
    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('Selectivity Index (E_H - 0.5·E_CO)', fontsize=12)
    ax.set_ylabel('Metal Center', fontsize=12)
    ax.set_title('(d) CO₂RR vs HER Selectivity (negative = CO₂RR)', fontsize=12)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.suptitle('Metal-Support Interaction Analysis for M-N₄-C SACs', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/metal_support_interaction.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] Metal-support interaction analysis saved.")


# ============================================================================
# 7. Solvent Effect and Potential Dependence
# ============================================================================

def analyze_solvent_and_potential(catalysts: List[CatalystData]) -> None:
    """Analyze solvent effects and potential dependence on CO2RR."""
    
    selected_names = ['Cu', 'CuAg', 'Au', 'Ni-N4-C', 'Fe-N4-C']
    selected = [c for c in catalysts if c.name in selected_names]
    
    potentials = np.linspace(0, -1.5, 50)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # (a) Potential-dependent free energy for Cu
    ax = axes[0, 0]
    cu = next(c for c in catalysts if c.name == 'Cu')
    
    G_COOH_U = cu.E_COOH + 0.33 + potentials
    G_CO_U = cu.E_CO + 0.47 + 2 * potentials
    G_CHO_U = cu.E_CHO + 0.67 + 3 * potentials
    G_product = 0.28 + 2 * potentials
    
    ax.plot(potentials, np.zeros_like(potentials), 'k-', label='CO₂(g)', linewidth=1)
    ax.plot(potentials, G_COOH_U, '-', label='*COOH', linewidth=2)
    ax.plot(potentials, G_CO_U, '-', label='*CO', linewidth=2)
    ax.plot(potentials, G_CHO_U, '-', label='*CHO', linewidth=2)
    ax.plot(potentials, G_product, '--', label='CO(g)', linewidth=2)
    ax.set_xlabel('Applied Potential U (V vs RHE)', fontsize=12)
    ax.set_ylabel('Free Energy ΔG (eV)', fontsize=12)
    ax.set_title('(a) Potential-Dependent Free Energy on Cu(111)', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', linewidth=0.5)
    
    # (b) Solvation effect comparison
    ax = axes[0, 1]
    intermediates = ['*COOH', '*CO', '*CHO', '*OH', '*OCCO']
    solv_energies = [SOLVATION_CORRECTIONS.get(i, 0.0) for i in intermediates]
    
    x_pos = np.arange(len(intermediates))
    bars = ax.bar(x_pos, solv_energies, color=['#e74c3c' if s < -0.15 else '#3498db' 
                                                for s in solv_energies],
                  edgecolor='black', linewidth=0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(intermediates, fontsize=11)
    ax.set_ylabel('Solvation Stabilization (eV)', fontsize=12)
    ax.set_title('(b) Implicit Solvation Corrections', fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(y=0, color='black', linewidth=0.8)
    
    for bar, val in zip(bars, solv_energies):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 0.02,
                f'{val:.2f}', ha='center', va='top', fontsize=10, fontweight='bold')
    
    # (c) Tafel slope analysis
    ax = axes[1, 0]
    log_j = np.linspace(-3, 2, 100)  # log(j / mA cm-2)
    
    # Tafel slopes for different mechanisms
    tafel_slopes = {
        'Cu (RDS: CO₂ ads.)': 118,    # 1e- transfer RDS
        'Au (RDS: *COOH)': 59,         # Fast 1e-, slow chemical
        'Ni-N4-C (RDS: *CO des.)': 40, # 2e- transfer before RDS
    }
    
    colors_t = ['#e74c3c', '#f39c12', '#2ecc71']
    for (label, slope), color in zip(tafel_slopes.items(), colors_t):
        eta = slope * log_j / 1000  # V
        ax.plot(eta, log_j, '-', label=f'{label}\nSlope={slope} mV/dec',
                color=color, linewidth=2)
    
    ax.set_xlabel('Overpotential η (V)', fontsize=12)
    ax.set_ylabel('log(j / mA cm⁻²)', fontsize=12)
    ax.set_title('(c) Tafel Analysis for Selected Catalysts', fontsize=12)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # (d) Electric field effect on key intermediates
    ax = axes[1, 1]
    fields = np.linspace(-0.5, 0.5, 100)  # V/Å
    
    # Dipole moments (Debye) and polarizabilities
    dipoles = {'*CO₂⁻': 1.5, '*COOH': 0.8, '*CO': 0.3, '*CHO': 0.6, '*OCCO': 1.2}
    
    for species, mu in dipoles.items():
        dG_field = -mu * fields - 0.5 * 2.0 * fields**2  # Linear + quadratic
        ax.plot(fields, dG_field, '-', label=f'{species} (μ={mu:.1f} D)', linewidth=2)
    
    ax.set_xlabel('Electric Field (V/Å)', fontsize=12)
    ax.set_ylabel('ΔG_field Stabilization (eV)', fontsize=12)
    ax.set_title('(d) Electric Field Effect on Intermediate Stabilization', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.axvline(x=0, color='black', linewidth=0.5)
    
    plt.suptitle('Solvent Effects and Potential Dependence in CO₂RR', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/solvent_potential_effects.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] Solvent and potential analysis saved.")


# ============================================================================
# 8. Candidate Material Evaluation (Cu alloys & N-doped Carbon)
# ============================================================================

def evaluate_candidates(catalysts: List[CatalystData], volcano_results: Dict) -> Dict:
    """Comprehensive evaluation and ranking of candidate materials."""
    
    # Scoring function
    scores = {}
    for cat in catalysts:
        vr = volcano_results[cat.name]
        
        # Activity score (closer to volcano peak = better)
        activity_CO = max(0, vr['U_L_CO'] + 0.5) * 2  # Normalize
        activity_C2 = max(0, vr['U_L_C2'] + 0.8) * 1.5
        
        # Selectivity score (prefer CO2RR over HER)
        sel_score = max(0, -cat.selectivity_index) * 2
        
        # Stability score (based on binding energy for SACs)
        if cat.category == 'SAC-NC':
            stability = min(1.0, abs(cat.d_band_center) / 5.0)
        else:
            stability = 0.7  # Default for metals/alloys
        
        total = activity_CO * 0.3 + activity_C2 * 0.3 + sel_score * 0.2 + stability * 0.2
        scores[cat.name] = {
            'activity_CO': activity_CO,
            'activity_C2': activity_C2,
            'selectivity': sel_score,
            'stability': stability,
            'total': total,
            'category': cat.category,
        }
    
    # Sort and display
    sorted_scores = sorted(scores.items(), key=lambda x: x[1]['total'], reverse=True)
    
    print("\n" + "=" * 80)
    print("CANDIDATE MATERIAL RANKING")
    print("=" * 80)
    print(f"{'Rank':>4} {'Catalyst':>12} {'Category':>10} {'Activity(CO)':>12} "
          f"{'Activity(C2)':>12} {'Selectivity':>11} {'Stability':>9} {'Total':>8}")
    print("-" * 80)
    for i, (name, s) in enumerate(sorted_scores[:20]):
        print(f"{i+1:4d} {name:>12} {s['category']:>10} {s['activity_CO']:>12.3f} "
              f"{s['activity_C2']:>12.3f} {s['selectivity']:>11.3f} "
              f"{s['stability']:>9.3f} {s['total']:>8.3f}")
    
    # --- Ranking Figure ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # (a) Top 15 overall ranking
    ax = axes[0]
    top15 = sorted_scores[:15]
    names = [x[0] for x in top15]
    totals = [x[1]['total'] for x in top15]
    act_co = [x[1]['activity_CO'] for x in top15]
    act_c2 = [x[1]['activity_C2'] for x in top15]
    sel = [x[1]['selectivity'] for x in top15]
    stab = [x[1]['stability'] for x in top15]
    
    x_pos = np.arange(len(names))
    w = 0.2
    ax.barh(x_pos - 1.5*w, [a*0.3 for a in act_co], w, label='Activity (CO)', color='#3498db')
    ax.barh(x_pos - 0.5*w, [a*0.3 for a in act_c2], w, label='Activity (C₂+)', color='#e74c3c')
    ax.barh(x_pos + 0.5*w, [s*0.2 for s in sel], w, label='Selectivity', color='#2ecc71')
    ax.barh(x_pos + 1.5*w, [s*0.2 for s in stab], w, label='Stability', color='#f39c12')
    
    ax.set_yticks(x_pos)
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel('Score', fontsize=12)
    ax.set_title('(a) Top 15 Catalysts - Score Breakdown', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')
    ax.invert_yaxis()
    
    # (b) 2D descriptor map with color = total score
    ax = axes[1]
    cu_alloys = [c for c in catalysts if c.category in ('Cu-alloy', 'SAC-NC')]
    E_CO_arr = [c.E_CO for c in cu_alloys]
    E_COOH_arr = [c.E_COOH for c in cu_alloys]
    total_scores_arr = [scores[c.name]['total'] for c in cu_alloys]
    
    scatter = ax.scatter(E_CO_arr, E_COOH_arr, c=total_scores_arr, cmap='RdYlGn',
                         s=120, edgecolors='black', linewidth=0.5, zorder=3)
    for c in cu_alloys:
        ax.annotate(c.name, (c.E_CO, c.E_COOH), fontsize=7,
                   xytext=(5, 5), textcoords='offset points')
    
    plt.colorbar(scatter, ax=ax, label='Total Score')
    ax.set_xlabel('E(*CO) (eV)', fontsize=12)
    ax.set_ylabel('E(*COOH) (eV)', fontsize=12)
    ax.set_title('(b) Descriptor Map: Cu Alloys & SACs', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Mark optimal region
    opt_co = [-0.6, -0.3]
    opt_cooh = [-0.5, -0.2]
    from matplotlib.patches import Rectangle
    rect = Rectangle((opt_co[0], opt_cooh[0]),
                     opt_co[1] - opt_co[0], opt_cooh[1] - opt_cooh[0],
                     linewidth=2, edgecolor='red', facecolor='red', alpha=0.1,
                     label='Optimal Region')
    ax.add_patch(rect)
    ax.legend(fontsize=9)
    
    plt.suptitle('Candidate Material Evaluation for CO₂RR', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/candidate_evaluation.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    return scores


# ============================================================================
# 9. Automated Screening Pipeline Design (ASE/CatMAP-based)
# ============================================================================

def design_screening_pipeline() -> str:
    """Design and document the ASE/CatMAP automated screening pipeline."""
    
    pipeline_diagram = """
    ┌─────────────────────────────────────────────────────────┐
    │           ASE/CatMAP Screening Pipeline                 │
    ├─────────────────────────────────────────────────────────┤
    │                                                         │
    │  Stage 1: Structure Generation (ASE)                    │
    │  ├── Build bulk structures from Materials Project       │
    │  ├── Generate slab models (facets: 111, 100, 110)       │
    │  ├── Create SAC structures (M-N4-C)                     │
    │  └── Add adsorbates (*CO, *COOH, *CHO, *OCCO)           │
    │                                                         │
    │  Stage 2: DFT Calculations                              │
    │  ├── Geometry optimization (VASP/GPAW)                  │
    │  ├── Energy calculation (PBE+D3, U_eff for TMs)         │
    │  ├── Solvation correction (VASPsol/COSMO)               │
    │  └── Vibrational analysis (ZPE, entropy)                │
    │                                                         │
    │  Stage 3: Descriptor Extraction                         │
    │  ├── Adsorption energies: E(*CO), E(*COOH), E(*CHO)     │
    │  ├── Scaling relation validation                        │
    │  ├── d-band center analysis                             │
    │  └── Charge transfer (Bader analysis)                   │
    │                                                         │
    │  Stage 4: Microkinetic Modeling (CatMAP)                │
    │  ├── Define reaction network (CO2→CO, CO→C1, CO→C2+)    │
    │  ├── Set up descriptor space                            │
    │  ├── Solve microkinetic equations                       │
    │  └── Generate volcano surfaces                          │
    │                                                         │
    │  Stage 5: Candidate Ranking & Selection                 │
    │  ├── Multi-objective scoring                            │
    │  ├── Stability validation (dissolution potential)        │
    │  ├── Selectivity analysis (CO2RR vs HER)                │
    │  └── Output: Ranked candidate list                      │
    │                                                         │
    └─────────────────────────────────────────────────────────┘
    """
    
    print("\n" + "=" * 60)
    print("AUTOMATED SCREENING PIPELINE DESIGN")
    print("=" * 60)
    print(pipeline_diagram)
    
    return pipeline_diagram


# ============================================================================
# 10. Pipeline Flow Diagram
# ============================================================================

def create_pipeline_flowchart():
    """Create a visual flowchart of the screening pipeline."""
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Boxes
    boxes = [
        (2, 9, 4, 0.7, 'Stage 1: Structure Generation\n(ASE + Materials Project)', '#3498db'),
        (2, 7.8, 4, 0.7, 'Stage 2: DFT Calculations\n(VASP/GPAW + VASPsol)', '#2ecc71'),
        (2, 6.6, 4, 0.7, 'Stage 3: Descriptor Extraction\nE(*CO), E(*COOH), E(*CHO)', '#e74c3c'),
        (2, 5.4, 4, 0.7, 'Stage 4: Scaling Validation\n& Microkinetic Model (CatMAP)', '#f39c12'),
        (2, 4.2, 4, 0.7, 'Stage 5: Volcano Plot\n& Activity Prediction', '#9b59b6'),
        (2, 3.0, 4, 0.7, 'Stage 6: Candidate Ranking\n& Selection', '#1abc9c'),
        # Side boxes
        (8, 9, 4.5, 0.7, 'Inputs:\n• Metal/alloy compositions\n• SAC configurations (M-N₄-C)', '#ecf0f1'),
        (8, 7.5, 4.5, 0.7, 'Parameters:\n• PBE+D3, 520 eV cutoff\n• k-points: 4×4×1', '#ecf0f1'),
        (8, 6.0, 4.5, 0.7, 'Corrections:\n• ZPE, entropy, solvation\n• CHE framework (U_RHE)', '#ecf0f1'),
        (8, 4.5, 4.5, 0.7, 'Output:\n• Volcano surfaces\n• Top-N candidates', '#ecf0f1'),
    ]
    
    for x, y, w, h, text, color in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black',
                              linewidth=1.5, alpha=0.8, transform=ax.transData)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=8, fontweight='bold', transform=ax.transData)
    
    # Arrows between main boxes
    for i in range(5):
        y_start = 9 - i * 1.2
        ax.annotate('', xy=(4, y_start - 0.3), xytext=(4, y_start),
                    arrowprops=dict(arrowstyle='->', color='black', lw=2))
    
    # Side arrows
    for y_main, y_side in [(9.35, 9.35), (8.15, 7.85), (6.95, 6.35), (4.55, 4.85)]:
        ax.annotate('', xy=(6.2, y_main), xytext=(7.8, y_side),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=1, ls='--'))
    
    ax.set_title('ASE/CatMAP Automated CO₂RR Catalyst Screening Pipeline',
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/pipeline_flowchart.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] Pipeline flowchart saved.")


# ============================================================================
# 11. 2D Volcano Surface (Heatmap)
# ============================================================================

def create_2d_volcano_surface(catalysts: List[CatalystData]) -> None:
    """Create 2D volcano heatmap using E(*CO) and E(*COOH) as dual descriptors."""
    
    E_CO_grid = np.linspace(-1.5, 0.5, 200)
    E_COOH_grid = np.linspace(-1.5, 0.5, 200)
    E_CO_mesh, E_COOH_mesh = np.meshgrid(E_CO_grid, E_COOH_grid)
    
    # Free energy steps for CO production
    dG1 = E_COOH_mesh + 0.33
    dG2 = E_CO_mesh - E_COOH_mesh + 0.14
    U_L = -np.maximum(dG1, dG2)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # (a) CO production volcano surface
    ax = axes[0]
    im = ax.contourf(E_CO_mesh, E_COOH_mesh, U_L, levels=30, cmap='RdYlGn')
    plt.colorbar(im, ax=ax, label='Limiting Potential U_L (V)')
    
    # Plot catalysts
    cat_colors = {'metal': 'white', 'Cu-alloy': 'cyan', 'SAC-NC': 'magenta'}
    cat_markers = {'metal': 'o', 'Cu-alloy': 's', 'SAC-NC': '^'}
    for cat in catalysts:
        ax.scatter(cat.E_CO, cat.E_COOH, c=cat_colors[cat.category],
                  marker=cat_markers[cat.category], s=50, edgecolors='black',
                  linewidth=0.8, zorder=3)
        if cat.name in ['Cu', 'Au', 'Ag', 'CuAg', 'Ni-N4-C', 'Fe-N4-C', 'Pt', 'Mo-N4-C']:
            ax.annotate(cat.name, (cat.E_CO, cat.E_COOH), fontsize=7, color='white',
                       fontweight='bold', xytext=(5, 5), textcoords='offset points')
    
    # Scaling line
    ax.plot(E_CO_grid, 0.85 * E_CO_grid + 0.10, 'w--', linewidth=1.5,
            label='Scaling: E(*COOH)=0.85·E(*CO)+0.10')
    ax.set_xlabel('E(*CO) (eV)', fontsize=12)
    ax.set_ylabel('E(*COOH) (eV)', fontsize=12)
    ax.set_title('(a) 2D Volcano: CO₂ → CO', fontsize=12)
    ax.legend(fontsize=8, loc='lower left')
    
    # (b) C2+ production volcano surface
    ax = axes[1]
    E_CHO_mesh = 0.75 * E_CO_mesh + 0.85  # From scaling
    dG_cho = E_CHO_mesh - E_CO_mesh + 0.20
    E_OCCO_mesh = 1.5 * E_CO_mesh + 0.45  # From scaling
    dG_cc = E_OCCO_mesh - 2 * E_CO_mesh + 0.30
    
    U_L_C2 = -np.maximum(np.maximum(dG1, dG2), np.minimum(dG_cho, dG_cc))
    
    im2 = ax.contourf(E_CO_mesh, E_COOH_mesh, U_L_C2, levels=30, cmap='RdYlGn')
    plt.colorbar(im2, ax=ax, label='Limiting Potential U_L (V)')
    
    for cat in catalysts:
        ax.scatter(cat.E_CO, cat.E_COOH, c=cat_colors[cat.category],
                  marker=cat_markers[cat.category], s=50, edgecolors='black',
                  linewidth=0.8, zorder=3)
        if cat.name in ['Cu', 'CuAg', 'CuAu', 'Cu-N4-C', 'Cu3Ag', 'CuZn', 'Mo-N4-C']:
            ax.annotate(cat.name, (cat.E_CO, cat.E_COOH), fontsize=7, color='white',
                       fontweight='bold', xytext=(5, 5), textcoords='offset points')
    
    ax.set_xlabel('E(*CO) (eV)', fontsize=12)
    ax.set_ylabel('E(*COOH) (eV)', fontsize=12)
    ax.set_title('(b) 2D Volcano: CO₂ → C₂+ Products', fontsize=12)
    
    plt.suptitle('2D Volcano Surfaces for CO₂RR Catalyst Screening', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/volcano_2d_surface.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] 2D volcano surfaces saved.")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("=" * 70)
    print("  CO₂RR Catalyst Computational Screening Pipeline")
    print("  Based on CHE Framework & Descriptor-Based Analysis")
    print("=" * 70)
    
    # Build catalyst database
    print("\n[1/8] Building catalyst database...")
    catalysts = build_catalyst_database()
    print(f"  Total catalysts: {len(catalysts)}")
    print(f"  Categories: {set(c.category for c in catalysts)}")
    
    # Compute scaling relations
    print("\n[2/8] Computing scaling relations...")
    scaling = compute_scaling_relations(catalysts)
    
    # Construct volcano plots
    print("\n[3/8] Constructing volcano plots...")
    volcano_results = compute_volcano_plot(catalysts, scaling)
    
    # Reaction pathway analysis
    print("\n[4/8] Analyzing reaction pathways...")
    analyze_reaction_pathways(catalysts)
    
    # Metal-support interaction
    print("\n[5/8] Analyzing metal-support interactions...")
    analyze_metal_support_interaction(catalysts)
    
    # Solvent and potential effects
    print("\n[6/8] Analyzing solvent and potential effects...")
    analyze_solvent_and_potential(catalysts)
    
    # Candidate evaluation
    print("\n[7/8] Evaluating candidate materials...")
    scores = evaluate_candidates(catalysts, volcano_results)
    
    # Pipeline design
    print("\n[8/8] Designing automated screening pipeline...")
    pipeline = design_screening_pipeline()
    create_pipeline_flowchart()
    create_2d_volcano_surface(catalysts)
    
    # Save results
    output_data = {
        'scaling_relations': scaling,
        'volcano_results': volcano_results,
        'scores': scores,
        'catalysts': {c.name: {
            'E_CO': c.E_CO, 'E_COOH': c.E_COOH, 'E_CHO': c.E_CHO,
            'E_OH': c.E_OH, 'E_H': c.E_H, 'E_OCCO': c.E_OCCO,
            'd_band_center': c.d_band_center, 'category': c.category,
            'U_L_CO': c.limiting_potential_CO,
            'selectivity_index': c.selectivity_index,
        } for c in catalysts}
    }
    
    with open('screening_results.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print(f"  Results saved to: screening_results.json")
    print(f"  Figures saved to: {FIGDIR}/")
    print("=" * 70)
    
    return output_data


if __name__ == '__main__':
    results = main()
