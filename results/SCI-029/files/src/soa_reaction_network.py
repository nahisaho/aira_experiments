"""
SOA Reaction Network Analysis System
Automated Chemical Reaction Network Generation and Analysis for 
Secondary Organic Aerosol Formation in Urban Atmospheres
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import r2_score, mean_squared_error
import json
import os

np.random.seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# =============================================================================
# Module 1: RMG-Based Automated Reaction Pathway Generation
# =============================================================================

class ReactionPathwayGenerator:
    """RMG-inspired automated reaction pathway generator for VOC oxidation."""
    
    # VOC species with their molecular properties
    VOC_SPECIES = {
        'isoprene': {'formula': 'C5H8', 'MW': 68.12, 'C_num': 5, 'H_num': 8, 'O_num': 0,
                     'double_bonds': 2, 'category': 'biogenic'},
        'alpha_pinene': {'formula': 'C10H16', 'MW': 136.23, 'C_num': 10, 'H_num': 16, 'O_num': 0,
                        'double_bonds': 1, 'category': 'biogenic'},
        'beta_caryophyllene': {'formula': 'C15H24', 'MW': 204.35, 'C_num': 15, 'H_num': 24, 'O_num': 0,
                              'double_bonds': 2, 'category': 'biogenic'},
        'toluene': {'formula': 'C7H8', 'MW': 92.14, 'C_num': 7, 'H_num': 8, 'O_num': 0,
                   'double_bonds': 4, 'category': 'anthropogenic'},
        'limonene': {'formula': 'C10H16', 'MW': 136.23, 'C_num': 10, 'H_num': 16, 'O_num': 0,
                    'double_bonds': 2, 'category': 'biogenic'},
    }
    
    OXIDANTS = ['OH', 'O3', 'NO3']
    
    REACTION_TYPES = [
        'H_abstraction', 'OH_addition', 'O3_ozonolysis', 'NO3_addition',
        'peroxy_radical_formation', 'alkoxy_radical_decomposition',
        'isomerization', 'autoxidation', 'accretion', 'hydrolysis'
    ]
    
    def __init__(self):
        self.reactions = []
        self.species = {}
        self.graph = nx.DiGraph()
        
    def generate_oxidation_products(self, voc_name):
        """Generate oxidation products for a given VOC."""
        if voc_name not in self.VOC_SPECIES:
            raise ValueError(f"Unknown VOC: {voc_name}")
        
        voc = self.VOC_SPECIES[voc_name]
        products = []
        
        # First-generation OH oxidation products
        for i in range(3):
            prod = {
                'name': f"{voc_name}_OH_prod{i+1}",
                'C_num': voc['C_num'],
                'H_num': voc['H_num'] - 1 + (1 if i == 0 else 0),
                'O_num': 1 + i,
                'MW': voc['MW'] + 16 * (1 + i),
                'generation': 1,
                'volatility': 'SVOC' if i < 2 else 'LVOC',
                'parent': voc_name,
                'oxidant': 'OH'
            }
            products.append(prod)
        
        # O3 ozonolysis products (for species with double bonds)
        if voc['double_bonds'] > 0:
            for i in range(2):
                c_split = max(1, voc['C_num'] // 2 + (i * (voc['C_num'] % 2)))
                prod = {
                    'name': f"{voc_name}_O3_prod{i+1}",
                    'C_num': c_split,
                    'H_num': max(2, c_split * 2),
                    'O_num': 2 + i,
                    'MW': c_split * 14 + 32 + 16 * i,
                    'generation': 1,
                    'volatility': 'IVOC' if c_split < 5 else 'SVOC',
                    'parent': voc_name,
                    'oxidant': 'O3'
                }
                products.append(prod)
        
        # NO3 oxidation products (nighttime chemistry)
        prod_no3 = {
            'name': f"{voc_name}_NO3_prod1",
            'C_num': voc['C_num'],
            'H_num': voc['H_num'],
            'O_num': 4,
            'MW': voc['MW'] + 62,
            'generation': 1,
            'volatility': 'SVOC',
            'parent': voc_name,
            'oxidant': 'NO3'
        }
        products.append(prod_no3)
        
        # Second-generation products (further oxidation)
        second_gen = []
        for prod in products[:3]:
            sg = {
                'name': f"{prod['name']}_2nd",
                'C_num': prod['C_num'],
                'H_num': max(2, prod['H_num'] - 2),
                'O_num': prod['O_num'] + 2,
                'MW': prod['MW'] + 32,
                'generation': 2,
                'volatility': 'LVOC' if prod['O_num'] >= 2 else 'SVOC',
                'parent': prod['name'],
                'oxidant': 'OH'
            }
            second_gen.append(sg)
        
        # ELVOC from autoxidation
        elvoc = {
            'name': f"{voc_name}_ELVOC",
            'C_num': voc['C_num'],
            'H_num': max(2, voc['H_num'] - 4),
            'O_num': 6 + voc['double_bonds'],
            'MW': voc['MW'] + 96 + 16 * voc['double_bonds'],
            'generation': 1,
            'volatility': 'ELVOC',
            'parent': voc_name,
            'oxidant': 'autoxidation'
        }
        products.append(elvoc)
        products.extend(second_gen)
        
        return products
    
    def build_reaction_network(self):
        """Build the full reaction network for all VOCs."""
        all_products = {}
        reaction_id = 0
        
        for voc_name in self.VOC_SPECIES:
            self.graph.add_node(voc_name, type='VOC', **self.VOC_SPECIES[voc_name])
            products = self.generate_oxidation_products(voc_name)
            
            for prod in products:
                self.graph.add_node(prod['name'], type='product', **prod)
                
                # Estimate rate constant (cm³ molecule⁻¹ s⁻¹)
                if prod['oxidant'] == 'OH':
                    k = np.random.uniform(1e-11, 3e-10)
                elif prod['oxidant'] == 'O3':
                    k = np.random.uniform(1e-18, 5e-16)
                elif prod['oxidant'] == 'NO3':
                    k = np.random.uniform(1e-14, 1e-11)
                else:
                    k = np.random.uniform(1e-3, 1e-1)  # autoxidation (s⁻¹)
                
                self.graph.add_edge(prod['parent'], prod['name'], 
                                   reaction_id=reaction_id,
                                   rate_constant=k,
                                   oxidant=prod['oxidant'],
                                   reaction_type=self._classify_reaction(prod))
                
                self.reactions.append({
                    'id': reaction_id,
                    'reactant': prod['parent'],
                    'product': prod['name'],
                    'oxidant': prod['oxidant'],
                    'rate_constant': k,
                    'type': self._classify_reaction(prod)
                })
                reaction_id += 1
                all_products[prod['name']] = prod
        
        self.species = all_products
        return self.graph
    
    def _classify_reaction(self, product):
        """Classify reaction type based on product properties."""
        if product['oxidant'] == 'OH':
            return 'OH_addition' if product['generation'] == 1 else 'multigenerational_oxidation'
        elif product['oxidant'] == 'O3':
            return 'ozonolysis'
        elif product['oxidant'] == 'NO3':
            return 'NO3_addition'
        else:
            return 'autoxidation'
    
    def get_network_statistics(self):
        """Calculate network statistics."""
        stats = {
            'num_species': self.graph.number_of_nodes(),
            'num_reactions': self.graph.number_of_edges(),
            'num_vocs': len(self.VOC_SPECIES),
            'avg_degree': np.mean([d for n, d in self.graph.degree()]),
            'max_path_length': 0,
            'connected_components': nx.number_weakly_connected_components(self.graph),
            'reaction_types': {},
        }
        
        for rxn in self.reactions:
            rt = rxn['type']
            stats['reaction_types'][rt] = stats['reaction_types'].get(rt, 0) + 1
        
        # Find longest pathway
        for voc in self.VOC_SPECIES:
            for node in self.graph.nodes():
                if node != voc and nx.has_path(self.graph, voc, node):
                    path_len = nx.shortest_path_length(self.graph, voc, node)
                    stats['max_path_length'] = max(stats['max_path_length'], path_len)
        
        return stats
    
    def visualize_network(self, filename='reaction_network.png'):
        """Visualize the reaction network."""
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        color_map = []
        size_map = []
        for node in self.graph.nodes():
            ndata = self.graph.nodes[node]
            if ndata.get('type') == 'VOC':
                color_map.append('#E74C3C')
                size_map.append(800)
            elif ndata.get('volatility') == 'ELVOC':
                color_map.append('#8E44AD')
                size_map.append(500)
            elif ndata.get('volatility') == 'LVOC':
                color_map.append('#2980B9')
                size_map.append(400)
            elif ndata.get('volatility') == 'SVOC':
                color_map.append('#27AE60')
                size_map.append(350)
            else:
                color_map.append('#F39C12')
                size_map.append(300)
        
        pos = nx.spring_layout(self.graph, k=2.5, iterations=100, seed=42)
        
        edge_colors = []
        for u, v in self.graph.edges():
            edata = self.graph.edges[u, v]
            if edata.get('oxidant') == 'OH':
                edge_colors.append('#E74C3C')
            elif edata.get('oxidant') == 'O3':
                edge_colors.append('#3498DB')
            elif edata.get('oxidant') == 'NO3':
                edge_colors.append('#2ECC71')
            else:
                edge_colors.append('#9B59B6')
        
        nx.draw_networkx_edges(self.graph, pos, edge_color=edge_colors, 
                              alpha=0.6, arrows=True, arrowsize=15, ax=ax)
        nx.draw_networkx_nodes(self.graph, pos, node_color=color_map, 
                              node_size=size_map, alpha=0.85, ax=ax)
        
        # Label only VOC nodes
        voc_labels = {n: n for n in self.graph.nodes() if self.graph.nodes[n].get('type') == 'VOC'}
        nx.draw_networkx_labels(self.graph, pos, voc_labels, font_size=9, 
                               font_weight='bold', ax=ax)
        
        # Legend
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        legend_elements = [
            Patch(facecolor='#E74C3C', label='VOC precursors'),
            Patch(facecolor='#8E44AD', label='ELVOC'),
            Patch(facecolor='#2980B9', label='LVOC'),
            Patch(facecolor='#27AE60', label='SVOC'),
            Patch(facecolor='#F39C12', label='IVOC'),
            Line2D([0], [0], color='#E74C3C', label='OH oxidation'),
            Line2D([0], [0], color='#3498DB', label='O₃ ozonolysis'),
            Line2D([0], [0], color='#2ECC71', label='NO₃ oxidation'),
            Line2D([0], [0], color='#9B59B6', label='Autoxidation'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=8)
        ax.set_title('VOC Oxidation Reaction Network for SOA Formation', fontsize=14, fontweight='bold')
        ax.axis('off')
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")


# =============================================================================
# Module 2: Gas-Particle Partitioning (UNIFAC/AIOMFAC-based)
# =============================================================================

class ThermodynamicPartitioning:
    """UNIFAC/AIOMFAC-inspired thermodynamic model for gas-particle partitioning."""
    
    # UNIFAC group interaction parameters (simplified)
    UNIFAC_GROUPS = {
        'CH3': {'R': 0.9011, 'Q': 0.848},
        'CH2': {'R': 0.6744, 'Q': 0.540},
        'CH': {'R': 0.4469, 'Q': 0.228},
        'C=C': {'R': 1.1167, 'Q': 0.867},
        'OH': {'R': 1.0000, 'Q': 1.200},
        'CHO': {'R': 0.9980, 'Q': 0.948},
        'COOH': {'R': 1.3013, 'Q': 1.224},
        'C=O': {'R': 1.6724, 'Q': 1.488},
        'ONO2': {'R': 1.6000, 'Q': 1.420},
        'OOH': {'R': 1.2000, 'Q': 1.100},
    }
    
    def __init__(self, T=298.15, RH=0.5, C_OA=10.0):
        """
        T: Temperature (K)
        RH: Relative humidity (0-1)
        C_OA: Organic aerosol concentration (µg/m³)
        """
        self.T = T
        self.RH = RH
        self.C_OA = C_OA
        self.R = 8.314  # J/(mol·K)
        
    def estimate_saturation_vapor_pressure(self, species):
        """Estimate p° using SIMPOL.1 group contribution method."""
        C_num = species.get('C_num', 5)
        O_num = species.get('O_num', 1)
        MW = species.get('MW', 100)
        
        # SIMPOL.1 parameterization (simplified)
        log10_p = 1.79 - 0.438 * C_num - 2.84 * (O_num / max(C_num, 1))
        
        # Temperature correction (Clausius-Clapeyron)
        dH_vap = (50 + 6.5 * C_num) * 1000  # J/mol
        T_ref = 298.15
        log10_p += (dH_vap / (2.303 * self.R)) * (1/T_ref - 1/self.T)
        
        return 10**log10_p  # in atm
    
    def calculate_activity_coefficient(self, species, composition=None):
        """Calculate activity coefficient using simplified UNIFAC model."""
        O_C_ratio = species.get('O_num', 1) / max(species.get('C_num', 5), 1)
        
        # Combinatorial contribution
        ln_gamma_c = 0.1 * (1 - O_C_ratio)
        
        # Residual contribution (simplified interaction)
        ln_gamma_r = -0.5 * O_C_ratio + 0.3 * O_C_ratio**2
        
        # Water interaction (AIOMFAC-like)
        if self.RH > 0.3:
            ln_gamma_w = 0.2 * self.RH * (1 - O_C_ratio)
        else:
            ln_gamma_w = 0.0
        
        gamma = np.exp(ln_gamma_c + ln_gamma_r + ln_gamma_w)
        return max(0.1, min(gamma, 10.0))
    
    def calculate_partitioning_coefficient(self, species):
        """Calculate Kp (m³/µg) for gas-particle partitioning."""
        MW = species.get('MW', 100)
        p_sat = self.estimate_saturation_vapor_pressure(species)
        gamma = self.calculate_activity_coefficient(species)
        
        # Kp = (RT) / (MW * gamma * p_sat * 10^6)
        # Converting to m³/µg
        p_sat_Pa = p_sat * 101325
        Kp = (self.R * self.T) / (MW * gamma * p_sat_Pa * 1e6)
        
        return Kp
    
    def calculate_particle_fraction(self, species):
        """Calculate fraction in particle phase (Fp)."""
        Kp = self.calculate_partitioning_coefficient(species)
        Fp = (Kp * self.C_OA) / (1 + Kp * self.C_OA)
        return Fp
    
    def calculate_effective_saturation_concentration(self, species):
        """Calculate C* (µg/m³) - effective saturation concentration."""
        Kp = self.calculate_partitioning_coefficient(species)
        if Kp > 0:
            C_star = 1.0 / Kp
        else:
            C_star = 1e10
        return C_star
    
    def run_partitioning_analysis(self, species_list):
        """Run partitioning analysis for all species."""
        results = []
        for sp in species_list:
            Kp = self.calculate_partitioning_coefficient(sp)
            Fp = self.calculate_particle_fraction(sp)
            C_star = self.calculate_effective_saturation_concentration(sp)
            gamma = self.calculate_activity_coefficient(sp)
            p_sat = self.estimate_saturation_vapor_pressure(sp)
            
            results.append({
                'name': sp['name'],
                'MW': sp.get('MW', 100),
                'O_C': sp.get('O_num', 1) / max(sp.get('C_num', 5), 1),
                'Kp': Kp,
                'Fp': Fp,
                'C_star': C_star,
                'gamma': gamma,
                'p_sat': p_sat,
                'volatility': sp.get('volatility', 'unknown'),
            })
        return results
    
    def plot_volatility_distribution(self, results, filename='volatility_distribution.png'):
        """Plot the volatility basis set distribution."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # VBS distribution
        log_cstar = [np.log10(max(r['C_star'], 1e-10)) for r in results]
        bins = np.arange(-6, 8, 1)
        
        colors_by_vol = {'ELVOC': '#8E44AD', 'LVOC': '#2980B9', 'SVOC': '#27AE60', 'IVOC': '#F39C12'}
        for vol_class in ['ELVOC', 'LVOC', 'SVOC', 'IVOC']:
            vals = [lc for lc, r in zip(log_cstar, results) if r['volatility'] == vol_class]
            if vals:
                axes[0].hist(vals, bins=bins, alpha=0.6, label=vol_class, 
                           color=colors_by_vol.get(vol_class, '#95A5A6'))
        
        axes[0].set_xlabel('log₁₀(C*) [µg/m³]', fontsize=11)
        axes[0].set_ylabel('Number of species', fontsize=11)
        axes[0].set_title('Volatility Basis Set Distribution', fontsize=12, fontweight='bold')
        axes[0].legend()
        axes[0].axvline(x=0, color='red', linestyle='--', alpha=0.5, label='C*=1 µg/m³')
        
        # Fp vs C*
        c_stars = [r['C_star'] for r in results]
        fps = [r['Fp'] for r in results]
        vol_colors = [colors_by_vol.get(r['volatility'], '#95A5A6') for r in results]
        
        axes[1].scatter(c_stars, fps, c=vol_colors, s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        # Theoretical curve
        c_star_range = np.logspace(-6, 6, 200)
        Fp_theory = (self.C_OA / c_star_range) / (1 + self.C_OA / c_star_range)
        axes[1].plot(c_star_range, Fp_theory, 'k--', alpha=0.5, label=f'Theory (COA={self.C_OA} µg/m³)')
        
        axes[1].set_xscale('log')
        axes[1].set_xlabel('C* [µg/m³]', fontsize=11)
        axes[1].set_ylabel('Particle fraction (Fp)', fontsize=11)
        axes[1].set_title('Gas-Particle Partitioning', fontsize=12, fontweight='bold')
        axes[1].legend()
        axes[1].set_ylim(-0.05, 1.05)
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")


# =============================================================================
# Module 3: ML Prediction of Photochemical Rate Constants
# =============================================================================

class RateConstantPredictor:
    """ML-based prediction of photochemical reaction rate constants."""
    
    def __init__(self):
        self.model_oh = None
        self.model_o3 = None
        self.model_no3 = None
        self.training_data = None
        
    def generate_training_data(self, n_samples=500):
        """Generate synthetic training data based on Evans-Polanyi relationships."""
        np.random.seed(42)
        
        data = []
        for _ in range(n_samples):
            C_num = np.random.randint(1, 16)
            H_num = max(2, int(2.1 * C_num + np.random.normal(0, 1)))
            O_num = np.random.randint(0, 8)
            double_bonds = np.random.randint(0, min(4, C_num // 2 + 1))
            MW = C_num * 12 + H_num + O_num * 16
            O_C = O_num / max(C_num, 1)
            H_C = H_num / max(C_num, 1)
            
            # Evans-Polanyi: log(k) = a * BDE + b
            # Features: C_num, O_C ratio, double bonds, MW
            BDE_eff = 400 - 10 * double_bonds - 5 * O_num + 2 * C_num + np.random.normal(0, 10)
            
            # OH rate constant
            log_k_OH = -10.5 + 0.8 * np.log(C_num + 1) + 0.3 * double_bonds - 0.001 * BDE_eff
            log_k_OH += np.random.normal(0, 0.15)
            
            # O3 rate constant
            log_k_O3 = -17.5 + 1.5 * double_bonds + 0.2 * np.log(C_num + 1)
            log_k_O3 += np.random.normal(0, 0.3)
            
            # NO3 rate constant
            log_k_NO3 = -13.0 + 0.7 * double_bonds + 0.15 * np.log(C_num + 1) - 0.1 * O_C
            log_k_NO3 += np.random.normal(0, 0.2)
            
            data.append({
                'C_num': C_num, 'H_num': H_num, 'O_num': O_num,
                'double_bonds': double_bonds, 'MW': MW,
                'O_C': O_C, 'H_C': H_C, 'BDE_eff': BDE_eff,
                'log_k_OH': log_k_OH, 'log_k_O3': log_k_O3, 'log_k_NO3': log_k_NO3,
            })
        
        self.training_data = data
        return data
    
    def train_models(self):
        """Train GBR models for each oxidant."""
        if self.training_data is None:
            self.generate_training_data()
        
        features = ['C_num', 'H_num', 'O_num', 'double_bonds', 'MW', 'O_C', 'H_C']
        X = np.array([[d[f] for f in features] for d in self.training_data])
        
        results = {}
        for target, name in [('log_k_OH', 'OH'), ('log_k_O3', 'O3'), ('log_k_NO3', 'NO3')]:
            y = np.array([d[target] for d in self.training_data])
            
            model = GradientBoostingRegressor(
                n_estimators=200, max_depth=5, learning_rate=0.1,
                min_samples_split=5, random_state=42
            )
            
            cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2')
            model.fit(X, y)
            y_pred = model.predict(X)
            
            r2 = r2_score(y, y_pred)
            rmse = np.sqrt(mean_squared_error(y, y_pred))
            
            results[name] = {
                'model': model,
                'r2_train': r2,
                'r2_cv': cv_scores.mean(),
                'r2_cv_std': cv_scores.std(),
                'rmse': rmse,
                'feature_importance': dict(zip(features, model.feature_importances_)),
            }
            
            if name == 'OH':
                self.model_oh = model
            elif name == 'O3':
                self.model_o3 = model
            else:
                self.model_no3 = model
        
        self.ml_results = results
        return results
    
    def predict(self, species, oxidant='OH'):
        """Predict rate constant for a species."""
        features = [
            species.get('C_num', 5),
            species.get('H_num', 8),
            species.get('O_num', 0),
            species.get('double_bonds', 0),
            species.get('MW', 100),
            species.get('O_num', 0) / max(species.get('C_num', 5), 1),
            species.get('H_num', 8) / max(species.get('C_num', 5), 1),
        ]
        X = np.array([features])
        
        if oxidant == 'OH' and self.model_oh:
            return 10**self.model_oh.predict(X)[0]
        elif oxidant == 'O3' and self.model_o3:
            return 10**self.model_o3.predict(X)[0]
        elif oxidant == 'NO3' and self.model_no3:
            return 10**self.model_no3.predict(X)[0]
        return None
    
    def plot_ml_results(self, filename='ml_rate_prediction.png'):
        """Plot ML model performance."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        features = ['C_num', 'H_num', 'O_num', 'double_bonds', 'MW', 'O_C', 'H_C']
        X = np.array([[d[f] for f in features] for d in self.training_data])
        
        # Parity plots for each oxidant
        for idx, (target, name, color) in enumerate([
            ('log_k_OH', 'OH', '#E74C3C'),
            ('log_k_O3', 'O₃', '#3498DB'),
            ('log_k_NO3', 'NO₃', '#2ECC71')
        ]):
            ax = axes[idx // 2][idx % 2]
            y = np.array([d[target] for d in self.training_data])
            model = self.ml_results[name.replace('₃', '3')]['model']
            y_pred = model.predict(X)
            
            ax.scatter(y, y_pred, alpha=0.3, s=20, c=color, edgecolors='none')
            lim = [min(y.min(), y_pred.min()) - 0.5, max(y.max(), y_pred.max()) + 0.5]
            ax.plot(lim, lim, 'k--', alpha=0.5)
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.set_xlabel(f'Actual log₁₀(k_{name})', fontsize=10)
            ax.set_ylabel(f'Predicted log₁₀(k_{name})', fontsize=10)
            
            r2 = self.ml_results[name.replace('₃', '3')]['r2_train']
            rmse = self.ml_results[name.replace('₃', '3')]['rmse']
            ax.set_title(f'{name} Rate Constants (R²={r2:.3f}, RMSE={rmse:.3f})', fontsize=11)
        
        # Feature importance
        ax = axes[1][1]
        importance_data = {}
        for name in ['OH', 'O3', 'NO3']:
            for feat, imp in self.ml_results[name]['feature_importance'].items():
                if feat not in importance_data:
                    importance_data[feat] = {}
                importance_data[feat][name] = imp
        
        x_pos = np.arange(len(features))
        width = 0.25
        for i, (name, color) in enumerate([('OH', '#E74C3C'), ('O3', '#3498DB'), ('NO3', '#2ECC71')]):
            vals = [importance_data[f][name] for f in features]
            ax.bar(x_pos + i * width, vals, width, label=name, color=color, alpha=0.7)
        
        ax.set_xticks(x_pos + width)
        ax.set_xticklabels(features, rotation=45, ha='right')
        ax.set_ylabel('Feature Importance', fontsize=10)
        ax.set_title('Feature Importance by Oxidant', fontsize=11)
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")


# =============================================================================
# Module 4: Atmospheric Box Model
# =============================================================================

class AtmosphericBoxModel:
    """Zero-dimensional atmospheric chemistry box model."""
    
    def __init__(self, T=298.15, P=101325, RH=0.5):
        self.T = T
        self.P = P
        self.RH = RH
        self.species_names = []
        self.initial_concentrations = {}
        self.reactions = []
        self.results = None
        
    def setup_urban_scenario(self, voc_concentrations=None):
        """Set up a typical urban atmospheric scenario."""
        # Default concentrations (molecules/cm³)
        self.initial_concentrations = {
            'OH': 1e6,
            'O3': 1e12,
            'NO': 5e10,
            'NO2': 2e11,
            'NO3': 1e8,
            'HO2': 1e8,
        }
        
        if voc_concentrations:
            self.initial_concentrations.update(voc_concentrations)
        else:
            self.initial_concentrations.update({
                'isoprene': 5e10,
                'alpha_pinene': 2e10,
                'toluene': 3e10,
                'limonene': 1e10,
                'beta_caryophyllene': 5e9,
            })
        
        self.species_names = list(self.initial_concentrations.keys())
        
    def add_reactions_from_network(self, reaction_list, rate_predictor=None):
        """Add reactions from the generated network."""
        for rxn in reaction_list:
            self.reactions.append(rxn)
            
            # Add product species if not present
            product = rxn['product']
            if product not in self.initial_concentrations:
                self.initial_concentrations[product] = 0.0
                self.species_names.append(product)
    
    def _build_ode_system(self):
        """Build the ODE system for the box model."""
        n_species = len(self.species_names)
        idx = {name: i for i, name in enumerate(self.species_names)}
        
        def dydt(t, y):
            rates = np.zeros(n_species)
            
            # Diurnal cycle for photolysis (simplified)
            hour = (t / 3600) % 24
            solar_factor = max(0, np.sin(np.pi * (hour - 6) / 12)) if 6 < hour < 18 else 0.0
            
            # OH production/loss (simplified HOx cycling)
            if 'OH' in idx:
                oh_prod = 1e6 * solar_factor  # photolysis source
                oh_loss = -1e-5 * y[idx['OH']]
                rates[idx['OH']] += oh_prod + oh_loss
            
            for rxn in self.reactions:
                reactant = rxn['reactant']
                product = rxn['product']
                oxidant = rxn['oxidant']
                k = rxn['rate_constant']
                
                if reactant not in idx or product not in idx:
                    continue
                
                # Calculate reaction rate
                if oxidant in idx and oxidant != 'autoxidation':
                    rate = k * y[idx[reactant]] * y[idx[oxidant]]
                    if oxidant == 'OH':
                        rate *= solar_factor
                    elif oxidant == 'NO3':
                        rate *= (1 - solar_factor + 0.01)  # nighttime chemistry
                else:
                    rate = k * y[idx[reactant]]
                
                rate = max(0, rate)
                
                rates[idx[reactant]] -= rate
                rates[idx[product]] += rate
            
            return rates
        
        return dydt
    
    def run_simulation(self, t_hours=48, dt_output=0.1):
        """Run the box model simulation."""
        t_span = (0, t_hours * 3600)
        t_eval = np.arange(0, t_hours * 3600, dt_output * 3600)
        
        y0 = np.array([self.initial_concentrations.get(sp, 0.0) for sp in self.species_names])
        
        dydt = self._build_ode_system()
        
        sol = solve_ivp(dydt, t_span, y0, t_eval=t_eval, method='BDF',
                       rtol=1e-6, atol=1e-8, max_step=60)
        
        self.results = {
            'time_hours': sol.t / 3600,
            'concentrations': {},
        }
        
        for i, name in enumerate(self.species_names):
            self.results['concentrations'][name] = sol.y[i]
        
        return self.results
    
    def plot_simulation_results(self, filename='box_model_results.png'):
        """Plot box model simulation results."""
        if self.results is None:
            raise ValueError("Run simulation first")
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        time = self.results['time_hours']
        
        # VOC decay
        ax = axes[0][0]
        vocs = ['isoprene', 'alpha_pinene', 'toluene', 'limonene', 'beta_caryophyllene']
        colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6']
        for voc, color in zip(vocs, colors):
            if voc in self.results['concentrations']:
                conc = self.results['concentrations'][voc]
                conc0 = conc[0] if conc[0] > 0 else 1
                ax.plot(time, conc / conc0, label=voc.replace('_', '-'), color=color, linewidth=1.5)
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Normalized concentration (C/C₀)')
        ax.set_title('VOC Decay over Time', fontweight='bold')
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1.1)
        ax.axvspan(0, 6, alpha=0.1, color='gray')
        ax.axvspan(18, 30, alpha=0.1, color='gray')
        ax.axvspan(42, 48, alpha=0.1, color='gray')
        
        # Oxidant concentrations
        ax = axes[0][1]
        for ox, color, label in [('OH', '#E74C3C', 'OH'), ('O3', '#3498DB', 'O₃'), 
                                   ('NO3', '#2ECC71', 'NO₃')]:
            if ox in self.results['concentrations']:
                ax.plot(time, self.results['concentrations'][ox], label=label, color=color, linewidth=1.5)
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (molec/cm³)')
        ax.set_title('Oxidant Concentrations', fontweight='bold')
        ax.legend()
        ax.set_yscale('log')
        
        # SOA product formation
        ax = axes[1][0]
        soa_products = [name for name in self.results['concentrations'] 
                       if any(x in name for x in ['ELVOC', 'LVOC', '2nd']) 
                       and name in self.results['concentrations']]
        
        # Group by parent VOC
        voc_soa = {}
        for prod in self.results['concentrations']:
            for voc in vocs:
                if voc in prod and prod != voc:
                    if voc not in voc_soa:
                        voc_soa[voc] = np.zeros_like(time)
                    voc_soa[voc] += np.maximum(0, self.results['concentrations'][prod])
        
        for voc, color in zip(vocs, colors):
            if voc in voc_soa:
                ax.plot(time, voc_soa[voc], label=f'{voc.replace("_", "-")} products', 
                       color=color, linewidth=1.5)
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Product concentration (molec/cm³)')
        ax.set_title('SOA Precursor Product Formation', fontweight='bold')
        ax.legend(fontsize=8)
        
        # Total SOA mass estimation
        ax = axes[1][1]
        total_soa = np.zeros_like(time)
        for voc in voc_soa:
            total_soa += voc_soa[voc]
        
        # Convert to µg/m³ (rough estimate)
        avg_MW = 200  # average molecular weight
        Na = 6.022e23
        soa_mass = total_soa * avg_MW / Na * 1e12  # µg/m³
        
        ax.plot(time, soa_mass, 'k-', linewidth=2, label='Total SOA')
        ax.fill_between(time, 0, soa_mass, alpha=0.2, color='blue')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('SOA mass concentration (µg/m³)')
        ax.set_title('Estimated Total SOA Mass', fontweight='bold')
        ax.legend()
        ax.axvspan(0, 6, alpha=0.1, color='gray', label='Night')
        ax.axvspan(18, 30, alpha=0.1, color='gray')
        ax.axvspan(42, 48, alpha=0.1, color='gray')
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")
        
        return soa_mass


# =============================================================================
# Module 5: Sensitivity Analysis
# =============================================================================

class SensitivityAnalyzer:
    """Sensitivity analysis for identifying key SOA formation pathways."""
    
    def __init__(self, box_model):
        self.box_model = box_model
        self.sensitivity_results = {}
    
    def one_at_a_time_analysis(self, perturbation_factor=0.1):
        """Perform OAT sensitivity analysis on reaction rate constants."""
        # Baseline run
        baseline_results = self.box_model.run_simulation()
        baseline_soa = self._calculate_total_soa(baseline_results)
        
        sensitivities = {}
        
        for i, rxn in enumerate(self.box_model.reactions):
            # Perturb rate constant
            original_k = rxn['rate_constant']
            rxn['rate_constant'] = original_k * (1 + perturbation_factor)
            
            perturbed_results = self.box_model.run_simulation()
            perturbed_soa = self._calculate_total_soa(perturbed_results)
            
            # Normalized sensitivity
            delta_soa = (perturbed_soa - baseline_soa) / baseline_soa if baseline_soa > 0 else 0
            sensitivity = delta_soa / perturbation_factor
            
            sensitivities[f"R{i}_{rxn['reactant']}->{rxn['product']}"] = {
                'sensitivity': sensitivity,
                'reaction_id': i,
                'reactant': rxn['reactant'],
                'product': rxn['product'],
                'oxidant': rxn['oxidant'],
                'rate_constant': original_k,
            }
            
            # Restore
            rxn['rate_constant'] = original_k
        
        self.sensitivity_results = sensitivities
        return sensitivities
    
    def _calculate_total_soa(self, results):
        """Calculate total SOA from simulation results."""
        total = 0
        vocs = ['isoprene', 'alpha_pinene', 'toluene', 'limonene', 'beta_caryophyllene']
        for name, conc in results['concentrations'].items():
            for voc in vocs:
                if voc in name and name != voc:
                    total += max(0, conc[-1])
        return total
    
    def morris_screening(self, n_trajectories=10):
        """Simplified Morris method for global sensitivity analysis."""
        n_params = len(self.box_model.reactions)
        
        # Generate Morris trajectories
        sensitivities_elementary = {i: [] for i in range(n_params)}
        
        for traj in range(n_trajectories):
            # Random base point
            perturbations = np.random.choice([-0.2, 0.2], size=n_params)
            
            # Baseline
            baseline = self._run_with_perturbations(np.zeros(n_params))
            
            for i in range(min(n_params, 30)):  # Limit for computational efficiency
                pert = np.zeros(n_params)
                pert[i] = perturbations[i]
                perturbed = self._run_with_perturbations(pert)
                
                if baseline > 0:
                    ee = (perturbed - baseline) / (baseline * perturbations[i])
                else:
                    ee = 0
                sensitivities_elementary[i].append(ee)
        
        morris_results = {}
        for i in range(min(n_params, 30)):
            ees = sensitivities_elementary[i]
            if ees:
                rxn = self.box_model.reactions[i]
                morris_results[f"R{i}"] = {
                    'mu_star': np.mean(np.abs(ees)),
                    'sigma': np.std(ees),
                    'reactant': rxn['reactant'],
                    'product': rxn['product'],
                    'oxidant': rxn['oxidant'],
                }
        
        return morris_results
    
    def _run_with_perturbations(self, perturbation_vector):
        """Run model with perturbed rate constants."""
        original_rates = []
        for i, rxn in enumerate(self.box_model.reactions):
            original_rates.append(rxn['rate_constant'])
            if i < len(perturbation_vector):
                rxn['rate_constant'] *= (1 + perturbation_vector[i])
        
        results = self.box_model.run_simulation()
        total_soa = self._calculate_total_soa(results)
        
        # Restore
        for i, rxn in enumerate(self.box_model.reactions):
            rxn['rate_constant'] = original_rates[i]
        
        return total_soa
    
    def plot_sensitivity_results(self, filename='sensitivity_analysis.png'):
        """Plot sensitivity analysis results."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Top sensitivities (OAT)
        ax = axes[0]
        sorted_sens = sorted(self.sensitivity_results.items(), 
                           key=lambda x: abs(x[1]['sensitivity']), reverse=True)[:15]
        
        names = [s[0].split('_', 1)[1] if '_' in s[0] else s[0] for s in sorted_sens]
        values = [s[1]['sensitivity'] for s in sorted_sens]
        colors = ['#E74C3C' if v > 0 else '#3498DB' for v in values]
        
        y_pos = np.arange(len(names))
        ax.barh(y_pos, values, color=colors, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names, fontsize=7)
        ax.set_xlabel('Normalized Sensitivity Coefficient')
        ax.set_title('Top 15 Sensitive Reactions (OAT)', fontweight='bold')
        ax.axvline(x=0, color='black', linewidth=0.5)
        
        # Sensitivity by oxidant
        ax = axes[1]
        oxidant_sens = {}
        for name, data in self.sensitivity_results.items():
            ox = data['oxidant']
            if ox not in oxidant_sens:
                oxidant_sens[ox] = []
            oxidant_sens[ox].append(abs(data['sensitivity']))
        
        ox_names = list(oxidant_sens.keys())
        ox_means = [np.mean(v) for v in oxidant_sens.values()]
        ox_stds = [np.std(v) for v in oxidant_sens.values()]
        
        colors = {'OH': '#E74C3C', 'O3': '#3498DB', 'NO3': '#2ECC71', 'autoxidation': '#9B59B6'}
        bar_colors = [colors.get(ox, '#95A5A6') for ox in ox_names]
        
        ax.bar(ox_names, ox_means, yerr=ox_stds, color=bar_colors, alpha=0.7, capsize=5)
        ax.set_ylabel('Mean |Sensitivity|')
        ax.set_title('Sensitivity by Oxidant Type', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")


# =============================================================================
# Module 6: SOA Yield Prediction
# =============================================================================

class SOAYieldPredictor:
    """SOA yield prediction for terpene/isoprene systems."""
    
    def __init__(self):
        self.yield_data = {}
    
    def two_product_model(self, C_OA, alpha1, Kp1, alpha2, Kp2):
        """Two-product Odum model for SOA yield."""
        Y = C_OA * (alpha1 * Kp1 / (1 + Kp1 * C_OA) + alpha2 * Kp2 / (1 + Kp2 * C_OA))
        return Y
    
    def vbs_yield(self, C_OA, alpha_bins, C_star_bins):
        """Volatility Basis Set SOA yield calculation."""
        Y = 0
        for alpha, C_star in zip(alpha_bins, C_star_bins):
            xi = 1 / (1 + C_star / max(C_OA, 1e-10))
            Y += alpha * xi
        return Y
    
    def calculate_yields_for_vocs(self, C_OA_range=None):
        """Calculate SOA yields for different VOC systems."""
        if C_OA_range is None:
            C_OA_range = np.logspace(-1, 3, 100)
        
        # VBS parameters from literature
        voc_params = {
            'α-pinene (high NOx)': {
                'alpha': [0.0, 0.038, 0.0, 0.140, 0.0, 0.050, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
            'α-pinene (low NOx)': {
                'alpha': [0.012, 0.122, 0.0, 0.210, 0.0, 0.030, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
            'isoprene (high NOx)': {
                'alpha': [0.0, 0.01, 0.0, 0.023, 0.0, 0.015, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
            'isoprene (low NOx)': {
                'alpha': [0.0, 0.009, 0.0, 0.030, 0.0, 0.015, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
            'β-caryophyllene': {
                'alpha': [0.05, 0.10, 0.0, 0.250, 0.0, 0.100, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
            'toluene': {
                'alpha': [0.0, 0.065, 0.0, 0.0, 0.168, 0.0, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
            'limonene': {
                'alpha': [0.03, 0.15, 0.0, 0.200, 0.0, 0.060, 0.0],
                'C_star': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
            },
        }
        
        yields = {}
        for voc_name, params in voc_params.items():
            Y = np.array([self.vbs_yield(c, params['alpha'], params['C_star']) for c in C_OA_range])
            yields[voc_name] = Y
        
        self.yield_data = {'C_OA': C_OA_range, 'yields': yields, 'params': voc_params}
        return yields
    
    def plot_soa_yields(self, filename='soa_yields.png'):
        """Plot SOA yield curves."""
        if not self.yield_data:
            self.calculate_yields_for_vocs()
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        C_OA = self.yield_data['C_OA']
        yields = self.yield_data['yields']
        
        colors = ['#E74C3C', '#C0392B', '#3498DB', '#2980B9', '#9B59B6', '#2ECC71', '#F39C12']
        linestyles = ['-', '--', '-', '--', '-', '-', '-']
        
        ax = axes[0]
        for (name, Y), color, ls in zip(yields.items(), colors, linestyles):
            ax.plot(C_OA, Y, color=color, linestyle=ls, linewidth=1.5, label=name)
        
        ax.set_xscale('log')
        ax.set_xlabel('Organic Aerosol Concentration C_OA (µg/m³)', fontsize=10)
        ax.set_ylabel('SOA Yield (Y)', fontsize=10)
        ax.set_title('SOA Yield Curves (VBS Model)', fontsize=12, fontweight='bold')
        ax.legend(fontsize=7, loc='upper left')
        ax.set_ylim(0, 0.6)
        ax.grid(True, alpha=0.3)
        ax.axvline(x=10, color='gray', linestyle=':', alpha=0.5, label='Typical urban COA')
        
        # VBS distribution at C_OA = 10 µg/m³
        ax = axes[1]
        C_star_bins = self.yield_data['params']['α-pinene (high NOx)']['C_star']
        bar_width = 0.12
        x = np.arange(len(C_star_bins))
        
        vocs_to_plot = ['α-pinene (high NOx)', 'isoprene (high NOx)', 'β-caryophyllene', 'toluene']
        plot_colors = ['#E74C3C', '#3498DB', '#9B59B6', '#2ECC71']
        
        for i, (voc, color) in enumerate(zip(vocs_to_plot, plot_colors)):
            alphas = self.yield_data['params'][voc]['alpha']
            ax.bar(x + i * bar_width, alphas, bar_width, label=voc, color=color, alpha=0.7)
        
        ax.set_xticks(x + 1.5 * bar_width)
        ax.set_xticklabels([f'{c}' for c in C_star_bins], fontsize=8)
        ax.set_xlabel('C* (µg/m³)', fontsize=10)
        ax.set_ylabel('Mass yield coefficient (α)', fontsize=10)
        ax.set_title('VBS Mass Yield Coefficients', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8)
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")

    def plot_temperature_rh_sensitivity(self, filename='temp_rh_sensitivity.png'):
        """Plot temperature and RH sensitivity of SOA yields."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Temperature sensitivity
        ax = axes[0]
        temps = np.arange(270, 320, 5)
        C_OA = 10.0  # µg/m³
        
        voc_params = {
            'α-pinene': {'alpha': [0.012, 0.122, 0.0, 0.210, 0.0, 0.030, 0.0],
                        'C_star_298': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
                        'dH': 30e3},
            'isoprene': {'alpha': [0.0, 0.009, 0.0, 0.030, 0.0, 0.015, 0.0],
                        'C_star_298': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
                        'dH': 40e3},
            'β-caryophyllene': {'alpha': [0.05, 0.10, 0.0, 0.250, 0.0, 0.100, 0.0],
                               'C_star_298': [0.01, 0.1, 1.0, 10, 100, 1000, 10000],
                               'dH': 35e3},
        }
        
        colors = ['#E74C3C', '#3498DB', '#9B59B6']
        for (voc, params), color in zip(voc_params.items(), colors):
            yields_T = []
            for T in temps:
                # Clausius-Clapeyron adjustment
                C_star_T = [c * np.exp(params['dH']/8.314 * (1/298.15 - 1/T)) 
                           for c in params['C_star_298']]
                Y = self.vbs_yield(C_OA, params['alpha'], C_star_T)
                yields_T.append(Y)
            ax.plot(temps - 273.15, yields_T, color=color, linewidth=2, label=voc, marker='o', markersize=4)
        
        ax.set_xlabel('Temperature (°C)', fontsize=10)
        ax.set_ylabel('SOA Yield', fontsize=10)
        ax.set_title('Temperature Dependence of SOA Yield', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # RH sensitivity (via COA effect)
        ax = axes[1]
        rh_values = np.arange(0.1, 1.0, 0.05)
        
        for (voc, params), color in zip(voc_params.items(), colors):
            yields_RH = []
            for rh in rh_values:
                # RH affects activity coefficients and water content
                gamma_correction = 1 + 0.3 * (rh - 0.5)
                C_star_RH = [c * gamma_correction for c in params['C_star_298']]
                Y = self.vbs_yield(C_OA, params['alpha'], C_star_RH)
                yields_RH.append(Y)
            ax.plot(rh_values * 100, yields_RH, color=color, linewidth=2, label=voc, marker='s', markersize=4)
        
        ax.set_xlabel('Relative Humidity (%)', fontsize=10)
        ax.set_ylabel('SOA Yield', fontsize=10)
        ax.set_title('RH Dependence of SOA Yield', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")


# =============================================================================
# Main Experiment Runner
# =============================================================================

def run_full_experiment():
    """Run the complete SOA reaction network analysis experiment."""
    results = {}
    
    print("=" * 70)
    print("SOA Reaction Network Analysis System - Full Experiment")
    print("=" * 70)
    
    # 1. Reaction Pathway Generation
    print("\n[1/6] Generating VOC Oxidation Reaction Network...")
    rpg = ReactionPathwayGenerator()
    graph = rpg.build_reaction_network()
    net_stats = rpg.get_network_statistics()
    rpg.visualize_network()
    print(f"  Species: {net_stats['num_species']}, Reactions: {net_stats['num_reactions']}")
    print(f"  Reaction types: {net_stats['reaction_types']}")
    results['network_stats'] = net_stats
    
    # 2. Thermodynamic Partitioning
    print("\n[2/6] Running Gas-Particle Partitioning Analysis...")
    thermo = ThermodynamicPartitioning(T=298.15, RH=0.5, C_OA=10.0)
    all_products = list(rpg.species.values())
    part_results = thermo.run_partitioning_analysis(all_products)
    thermo.plot_volatility_distribution(part_results)
    
    # Summary statistics
    fp_values = [r['Fp'] for r in part_results]
    print(f"  Mean particle fraction: {np.mean(fp_values):.3f}")
    print(f"  Species with Fp > 0.5: {sum(1 for fp in fp_values if fp > 0.5)}/{len(fp_values)}")
    results['partitioning'] = {
        'mean_Fp': float(np.mean(fp_values)),
        'n_particle_phase': int(sum(1 for fp in fp_values if fp > 0.5)),
        'n_total': len(fp_values),
    }
    
    # 3. ML Rate Constant Prediction
    print("\n[3/6] Training ML Rate Constant Predictors...")
    predictor = RateConstantPredictor()
    predictor.generate_training_data(n_samples=500)
    ml_results = predictor.train_models()
    predictor.plot_ml_results()
    
    for ox, res in ml_results.items():
        print(f"  {ox}: R²(CV) = {res['r2_cv']:.3f} ± {res['r2_cv_std']:.3f}, RMSE = {res['rmse']:.3f}")
    results['ml_performance'] = {ox: {'r2_cv': float(res['r2_cv']), 'rmse': float(res['rmse'])} 
                                  for ox, res in ml_results.items()}
    
    # 4. Box Model Simulation
    print("\n[4/6] Running Atmospheric Box Model Simulation...")
    box_model = AtmosphericBoxModel(T=298.15, P=101325, RH=0.5)
    box_model.setup_urban_scenario()
    box_model.add_reactions_from_network(rpg.reactions)
    sim_results = box_model.run_simulation(t_hours=48)
    soa_mass = box_model.plot_simulation_results()
    
    print(f"  Peak SOA mass: {np.max(soa_mass):.3f} µg/m³")
    print(f"  Final SOA mass: {soa_mass[-1]:.3f} µg/m³")
    results['box_model'] = {
        'peak_soa': float(np.max(soa_mass)),
        'final_soa': float(soa_mass[-1]),
        'simulation_hours': 48,
    }
    
    # 5. Sensitivity Analysis
    print("\n[5/6] Running Sensitivity Analysis...")
    sa = SensitivityAnalyzer(box_model)
    sens_results = sa.one_at_a_time_analysis(perturbation_factor=0.1)
    sa.plot_sensitivity_results()
    
    top_sensitive = sorted(sens_results.items(), 
                          key=lambda x: abs(x[1]['sensitivity']), reverse=True)[:5]
    print("  Top 5 most sensitive reactions:")
    for name, data in top_sensitive:
        print(f"    {name}: S = {data['sensitivity']:.4f} ({data['oxidant']})")
    results['sensitivity'] = {name: float(data['sensitivity']) for name, data in top_sensitive}
    
    # 6. SOA Yield Prediction
    print("\n[6/6] Predicting SOA Yields...")
    yield_pred = SOAYieldPredictor()
    yields = yield_pred.calculate_yields_for_vocs()
    yield_pred.plot_soa_yields()
    yield_pred.plot_temperature_rh_sensitivity()
    
    # Report yields at C_OA = 10 µg/m³
    C_OA_10_idx = np.argmin(np.abs(yield_pred.yield_data['C_OA'] - 10))
    print("  SOA Yields at C_OA = 10 µg/m³:")
    yield_at_10 = {}
    for name, Y in yields.items():
        print(f"    {name}: Y = {Y[C_OA_10_idx]:.4f}")
        yield_at_10[name] = float(Y[C_OA_10_idx])
    results['soa_yields'] = yield_at_10
    
    print("\n" + "=" * 70)
    print("Experiment Complete!")
    print("=" * 70)
    
    # Save results
    results_path = os.path.join(os.path.dirname(FIGURES_DIR), 'experiment_results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")
    
    return results


if __name__ == '__main__':
    results = run_full_experiment()
