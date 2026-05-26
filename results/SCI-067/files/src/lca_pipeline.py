"""
AutoLCA: AI-Driven Life Cycle Assessment Automation Pipeline
=============================================================
Modules:
1. NLP-based process tree construction
2. Ecoinvent database auto-matching (TF-IDF + cosine similarity)
3. Uncertainty propagation (Monte Carlo & Taylor expansion)
4. Hotspot analysis & scenario comparison
5. Scope 3 emissions estimation (ML-based)
6. EV battery manufacturing case study
"""

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
import json
import re
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# Module 1: NLP-based Process Tree Construction
# ============================================================

@dataclass
class ProcessNode:
    name: str
    category: str
    inputs: Dict[str, float] = field(default_factory=dict)
    outputs: Dict[str, float] = field(default_factory=dict)
    emissions: Dict[str, float] = field(default_factory=dict)
    unit: str = "kg"
    uncertainty: float = 0.1  # coefficient of variation

class NLPProcessTreeBuilder:
    """Builds process trees from textual descriptions using NLP techniques."""

    PROCESS_PATTERNS = [
        r'(?:production|manufacturing|processing|extraction|refining|assembly)\s+of\s+(\w[\w\s]+)',
        r'(\w[\w\s]+)\s+(?:production|manufacturing|processing|extraction|refining)',
        r'(?:convert|transform|synthesize)\s+(\w[\w\s]+)',
    ]

    FLOW_PATTERNS = [
        r'(\d+\.?\d*)\s*(kg|kWh|MJ|L|m3|t)\s+(?:of\s+)?(\w[\w\s]+)',
        r'(\w[\w\s]+):\s*(\d+\.?\d*)\s*(kg|kWh|MJ|L|m3|t)',
    ]

    EMISSION_KEYWORDS = {
        'CO2': ['carbon dioxide', 'co2', 'CO2'],
        'CH4': ['methane', 'ch4'],
        'N2O': ['nitrous oxide', 'n2o'],
        'SO2': ['sulfur dioxide', 'so2'],
        'NOx': ['nitrogen oxides', 'nox'],
        'PM': ['particulate matter', 'pm2.5', 'pm10'],
    }

    def __init__(self):
        self.graph = nx.DiGraph()
        self.processes = {}

    def extract_processes(self, text: str) -> List[str]:
        processes = []
        for pattern in self.PROCESS_PATTERNS:
            matches = re.findall(pattern, text, re.IGNORECASE)
            processes.extend([m.strip() for m in matches])
        return list(set(processes))

    def extract_flows(self, text: str) -> List[Tuple[str, float, str]]:
        flows = []
        for pattern in self.FLOW_PATTERNS:
            matches = re.findall(pattern, text, re.IGNORECASE)
            for m in matches:
                if len(m) == 3:
                    try:
                        flows.append((m[2].strip(), float(m[0]), m[1]))
                    except ValueError:
                        pass
        return flows

    def build_ev_battery_tree(self) -> nx.DiGraph:
        """Build a complete process tree for EV battery (NMC 811) manufacturing."""
        processes = {
            'lithium_mining': ProcessNode(
                name='Lithium Carbonate Mining & Refining',
                category='raw_material',
                inputs={'brine_water': 500.0, 'electricity': 15.0, 'chemicals': 2.5},
                outputs={'lithium_carbonate': 1.0},
                emissions={'CO2': 5.3, 'SO2': 0.02, 'water_consumption': 469.0},
                uncertainty=0.15
            ),
            'nickel_mining': ProcessNode(
                name='Nickel Sulfate Production',
                category='raw_material',
                inputs={'nickel_ore': 80.0, 'electricity': 25.0, 'sulfuric_acid': 3.0},
                outputs={'nickel_sulfate': 1.0},
                emissions={'CO2': 12.4, 'SO2': 0.15, 'NOx': 0.08},
                uncertainty=0.20
            ),
            'cobalt_mining': ProcessNode(
                name='Cobalt Sulfate Production',
                category='raw_material',
                inputs={'cobalt_ore': 120.0, 'electricity': 30.0},
                outputs={'cobalt_sulfate': 1.0},
                emissions={'CO2': 18.2, 'SO2': 0.25, 'PM': 0.05},
                uncertainty=0.25
            ),
            'manganese_mining': ProcessNode(
                name='Manganese Sulfate Production',
                category='raw_material',
                inputs={'manganese_ore': 15.0, 'electricity': 5.0},
                outputs={'manganese_sulfate': 1.0},
                emissions={'CO2': 2.1, 'SO2': 0.01},
                uncertainty=0.12
            ),
            'graphite_production': ProcessNode(
                name='Synthetic Graphite Production',
                category='raw_material',
                inputs={'petroleum_coke': 3.5, 'electricity': 35.0, 'natural_gas': 8.0},
                outputs={'graphite': 1.0},
                emissions={'CO2': 15.8, 'NOx': 0.04, 'PM': 0.02},
                uncertainty=0.18
            ),
            'cathode_production': ProcessNode(
                name='NMC 811 Cathode Production',
                category='component',
                inputs={'nickel_sulfate': 0.65, 'manganese_sulfate': 0.08,
                        'cobalt_sulfate': 0.08, 'lithium_carbonate': 0.22,
                        'electricity': 20.0, 'natural_gas': 5.0},
                outputs={'cathode_material': 1.0},
                emissions={'CO2': 8.5, 'NOx': 0.03},
                uncertainty=0.15
            ),
            'anode_production': ProcessNode(
                name='Graphite Anode Production',
                category='component',
                inputs={'graphite': 1.1, 'copper_foil': 0.15, 'binder': 0.05,
                        'electricity': 8.0},
                outputs={'anode': 1.0},
                emissions={'CO2': 4.2, 'PM': 0.01},
                uncertainty=0.12
            ),
            'electrolyte_production': ProcessNode(
                name='Electrolyte Production',
                category='component',
                inputs={'lithium_salt': 0.15, 'organic_solvent': 0.85,
                        'electricity': 3.0},
                outputs={'electrolyte': 1.0},
                emissions={'CO2': 2.8, 'VOC': 0.05},
                uncertainty=0.20
            ),
            'separator_production': ProcessNode(
                name='Separator Production',
                category='component',
                inputs={'polyethylene': 0.8, 'electricity': 5.0},
                outputs={'separator': 1.0},
                emissions={'CO2': 3.5},
                uncertainty=0.10
            ),
            'cell_assembly': ProcessNode(
                name='Cell Assembly',
                category='manufacturing',
                inputs={'cathode_material': 0.35, 'anode': 0.25,
                        'electrolyte': 0.12, 'separator': 0.03,
                        'aluminum_casing': 0.15, 'electricity': 45.0,
                        'dry_room_energy': 20.0},
                outputs={'battery_cell': 1.0},
                emissions={'CO2': 22.0, 'NOx': 0.05, 'PM': 0.02},
                uncertainty=0.15
            ),
            'module_assembly': ProcessNode(
                name='Battery Module Assembly',
                category='manufacturing',
                inputs={'battery_cell': 12.0, 'bms_electronics': 0.5,
                        'cooling_system': 2.0, 'connectors': 0.3,
                        'electricity': 10.0},
                outputs={'battery_module': 1.0},
                emissions={'CO2': 5.0, 'PM': 0.01},
                uncertainty=0.10
            ),
            'pack_assembly': ProcessNode(
                name='Battery Pack Assembly',
                category='manufacturing',
                inputs={'battery_module': 8.0, 'pack_housing': 15.0,
                        'thermal_management': 5.0, 'electricity': 15.0},
                outputs={'battery_pack': 1.0},
                emissions={'CO2': 8.0},
                uncertainty=0.10
            ),
        }

        G = nx.DiGraph()
        edges = [
            ('lithium_mining', 'cathode_production'),
            ('nickel_mining', 'cathode_production'),
            ('cobalt_mining', 'cathode_production'),
            ('manganese_mining', 'cathode_production'),
            ('graphite_production', 'anode_production'),
            ('cathode_production', 'cell_assembly'),
            ('anode_production', 'cell_assembly'),
            ('electrolyte_production', 'cell_assembly'),
            ('separator_production', 'cell_assembly'),
            ('cell_assembly', 'module_assembly'),
            ('module_assembly', 'pack_assembly'),
        ]

        for name, proc in processes.items():
            G.add_node(name, process=proc)

        for src, dst in edges:
            G.add_edge(src, dst)

        self.graph = G
        self.processes = processes
        return G


# ============================================================
# Module 2: Ecoinvent Database Auto-Matching
# ============================================================

class EcoinventMatcher:
    """TF-IDF based fuzzy matching to Ecoinvent database entries."""

    ECOINVENT_DB = [
        {'id': 'ei_001', 'name': 'lithium carbonate production', 'category': 'chemical',
         'location': 'GLO', 'unit': 'kg', 'gwp': 5.3},
        {'id': 'ei_002', 'name': 'nickel sulfate production', 'category': 'chemical',
         'location': 'GLO', 'unit': 'kg', 'gwp': 12.4},
        {'id': 'ei_003', 'name': 'cobalt sulfate production', 'category': 'chemical',
         'location': 'GLO', 'unit': 'kg', 'gwp': 18.2},
        {'id': 'ei_004', 'name': 'manganese sulfate production', 'category': 'chemical',
         'location': 'GLO', 'unit': 'kg', 'gwp': 2.1},
        {'id': 'ei_005', 'name': 'synthetic graphite production', 'category': 'chemical',
         'location': 'GLO', 'unit': 'kg', 'gwp': 15.8},
        {'id': 'ei_006', 'name': 'electricity production, hard coal', 'category': 'energy',
         'location': 'CN', 'unit': 'kWh', 'gwp': 1.1},
        {'id': 'ei_007', 'name': 'electricity production, natural gas', 'category': 'energy',
         'location': 'EU', 'unit': 'kWh', 'gwp': 0.5},
        {'id': 'ei_008', 'name': 'electricity production, wind', 'category': 'energy',
         'location': 'GLO', 'unit': 'kWh', 'gwp': 0.012},
        {'id': 'ei_009', 'name': 'aluminum production, primary', 'category': 'metal',
         'location': 'GLO', 'unit': 'kg', 'gwp': 8.0},
        {'id': 'ei_010', 'name': 'copper production, primary', 'category': 'metal',
         'location': 'GLO', 'unit': 'kg', 'gwp': 3.5},
        {'id': 'ei_011', 'name': 'polyethylene production, HDPE', 'category': 'plastic',
         'location': 'GLO', 'unit': 'kg', 'gwp': 2.0},
        {'id': 'ei_012', 'name': 'natural gas, burned in furnace', 'category': 'energy',
         'location': 'GLO', 'unit': 'MJ', 'gwp': 0.067},
        {'id': 'ei_013', 'name': 'sulfuric acid production', 'category': 'chemical',
         'location': 'GLO', 'unit': 'kg', 'gwp': 0.09},
        {'id': 'ei_014', 'name': 'transport, freight, lorry', 'category': 'transport',
         'location': 'GLO', 'unit': 'tkm', 'gwp': 0.11},
        {'id': 'ei_015', 'name': 'battery cell production, Li-ion', 'category': 'manufacturing',
         'location': 'CN', 'unit': 'kg', 'gwp': 12.5},
    ]

    def __init__(self):
        self.db_names = [e['name'] for e in self.ECOINVENT_DB]
        self.vectorizer = TfidfVectorizer(ngram_range=(1, 3), analyzer='char_wb')
        self.tfidf_matrix = self.vectorizer.fit_transform(self.db_names)

    def match(self, query: str, top_k: int = 3) -> List[Dict]:
        query_vec = self.vectorizer.transform([query.lower()])
        similarities = cosine_similarity(query_vec, self.tfidf_matrix)[0]
        top_indices = similarities.argsort()[-top_k:][::-1]
        results = []
        for idx in top_indices:
            results.append({
                'match': self.ECOINVENT_DB[idx],
                'similarity': float(similarities[idx]),
                'confidence': 'high' if similarities[idx] > 0.5 else
                             'medium' if similarities[idx] > 0.3 else 'low'
            })
        return results

    def auto_match_process_tree(self, processes: Dict[str, ProcessNode]) -> Dict:
        matching_results = {}
        for proc_id, proc in processes.items():
            matches = self.match(proc.name)
            matching_results[proc_id] = {
                'process_name': proc.name,
                'best_match': matches[0] if matches else None,
                'all_matches': matches
            }
        return matching_results


# ============================================================
# Module 3: Uncertainty Propagation
# ============================================================

class UncertaintyPropagation:
    """Monte Carlo and Taylor expansion methods for LCA uncertainty."""

    def __init__(self, n_simulations: int = 10000):
        self.n_simulations = n_simulations

    def monte_carlo(self, processes: Dict[str, ProcessNode],
                    graph: nx.DiGraph) -> Dict:
        """Run Monte Carlo simulation across the process tree."""
        results = {
            'gwp_samples': np.zeros(self.n_simulations),
            'process_contributions': {},
        }

        for proc_id, proc in processes.items():
            co2 = proc.emissions.get('CO2', 0)
            cv = proc.uncertainty
            samples = np.random.lognormal(
                mean=np.log(max(co2, 0.001)),
                sigma=cv,
                size=self.n_simulations
            )
            results['process_contributions'][proc_id] = samples
            results['gwp_samples'] += samples

        results['statistics'] = {
            'mean': float(np.mean(results['gwp_samples'])),
            'std': float(np.std(results['gwp_samples'])),
            'cv': float(np.std(results['gwp_samples']) / np.mean(results['gwp_samples'])),
            'p5': float(np.percentile(results['gwp_samples'], 5)),
            'p25': float(np.percentile(results['gwp_samples'], 25)),
            'median': float(np.percentile(results['gwp_samples'], 50)),
            'p75': float(np.percentile(results['gwp_samples'], 75)),
            'p95': float(np.percentile(results['gwp_samples'], 95)),
        }
        return results

    def taylor_expansion(self, processes: Dict[str, ProcessNode]) -> Dict:
        """First-order Taylor expansion for uncertainty propagation."""
        mean_total = sum(p.emissions.get('CO2', 0) for p in processes.values())
        var_total = sum(
            (p.emissions.get('CO2', 0) * p.uncertainty) ** 2
            for p in processes.values()
        )
        std_total = np.sqrt(var_total)
        return {
            'mean': mean_total,
            'std': std_total,
            'cv': std_total / mean_total if mean_total > 0 else 0,
            'ci_95_lower': mean_total - 1.96 * std_total,
            'ci_95_upper': mean_total + 1.96 * std_total,
        }

    def compare_methods(self, processes, graph) -> Dict:
        mc = self.monte_carlo(processes, graph)
        taylor = self.taylor_expansion(processes)
        return {
            'monte_carlo': mc['statistics'],
            'taylor': taylor,
            'relative_difference_mean': abs(mc['statistics']['mean'] - taylor['mean']) / taylor['mean'] * 100,
            'relative_difference_std': abs(mc['statistics']['std'] - taylor['std']) / taylor['std'] * 100,
        }


# ============================================================
# Module 4: Hotspot Analysis & Scenario Comparison
# ============================================================

class HotspotAnalyzer:
    """Identify environmental hotspots and compare scenarios."""

    def analyze(self, processes: Dict[str, ProcessNode]) -> pd.DataFrame:
        data = []
        total_co2 = sum(p.emissions.get('CO2', 0) for p in processes.values())
        for pid, proc in processes.items():
            co2 = proc.emissions.get('CO2', 0)
            data.append({
                'process_id': pid,
                'process_name': proc.name,
                'category': proc.category,
                'co2_kg': co2,
                'contribution_pct': co2 / total_co2 * 100 if total_co2 > 0 else 0,
                'electricity_kwh': proc.inputs.get('electricity', 0),
            })
        df = pd.DataFrame(data).sort_values('co2_kg', ascending=False)
        df['cumulative_pct'] = df['contribution_pct'].cumsum()
        return df

    def scenario_comparison(self, base_processes: Dict[str, ProcessNode]) -> Dict:
        """Compare baseline vs. improvement scenarios."""
        scenarios = {}

        # Baseline (China grid, coal-heavy)
        scenarios['baseline_china'] = self._calc_total_gwp(base_processes, grid_factor=1.1)

        # Scenario 1: EU grid (natural gas mix)
        scenarios['eu_grid'] = self._calc_total_gwp(base_processes, grid_factor=0.5)

        # Scenario 2: Renewable energy
        scenarios['renewable'] = self._calc_total_gwp(base_processes, grid_factor=0.05)

        # Scenario 3: LFP chemistry (lower material impacts)
        lfp_processes = self._modify_to_lfp(base_processes)
        scenarios['lfp_chemistry'] = self._calc_total_gwp(lfp_processes, grid_factor=1.1)

        # Scenario 4: Recycled materials (30% recycled content)
        scenarios['recycled_30pct'] = self._calc_total_gwp(
            base_processes, grid_factor=1.1, recycling_factor=0.7)

        # Scenario 5: Best case (renewable + LFP + recycling)
        scenarios['best_case'] = self._calc_total_gwp(
            lfp_processes, grid_factor=0.05, recycling_factor=0.7)

        return scenarios

    def _calc_total_gwp(self, processes, grid_factor=1.1, recycling_factor=1.0):
        total = 0
        breakdown = {}
        for pid, proc in processes.items():
            co2 = proc.emissions.get('CO2', 0)
            elec = proc.inputs.get('electricity', 0)
            elec_co2 = elec * grid_factor
            material_co2 = co2 * recycling_factor
            proc_total = material_co2 + elec_co2
            total += proc_total
            breakdown[pid] = proc_total
        return {'total_gwp': total, 'breakdown': breakdown}

    def _modify_to_lfp(self, processes):
        """Modify NMC processes to approximate LFP chemistry."""
        lfp = {}
        for pid, proc in processes.items():
            new_proc = ProcessNode(
                name=proc.name, category=proc.category,
                inputs=dict(proc.inputs), outputs=dict(proc.outputs),
                emissions=dict(proc.emissions), unit=proc.unit,
                uncertainty=proc.uncertainty
            )
            if pid == 'cobalt_mining':
                new_proc.emissions['CO2'] = 0  # No cobalt in LFP
            if pid == 'nickel_mining':
                new_proc.emissions['CO2'] *= 0.1  # Minimal nickel
            if pid == 'cathode_production':
                new_proc.emissions['CO2'] *= 0.6  # Simpler process
            lfp[pid] = new_proc
        return lfp


# ============================================================
# Module 5: Scope 3 Emissions Estimation
# ============================================================

class Scope3Estimator:
    """ML-based Scope 3 emissions estimation using proxy data."""

    SCOPE3_CATEGORIES = [
        'purchased_goods', 'capital_goods', 'fuel_energy',
        'upstream_transport', 'waste', 'business_travel',
        'employee_commuting', 'upstream_leased', 'downstream_transport',
        'processing', 'use_phase', 'end_of_life',
        'downstream_leased', 'franchises', 'investments'
    ]

    def __init__(self):
        self.models = {}

    def generate_training_data(self, n_samples=500):
        """Generate synthetic training data for Scope 3 estimation."""
        np.random.seed(42)
        data = {
            'revenue_musd': np.random.lognormal(5, 1.5, n_samples),
            'employees': np.random.lognormal(6, 1.2, n_samples).astype(int),
            'scope1_tco2': np.random.lognormal(8, 1.0, n_samples),
            'scope2_tco2': np.random.lognormal(7, 1.0, n_samples),
            'industry_code': np.random.choice([20, 25, 28, 29, 35], n_samples),
            'supplier_count': np.random.lognormal(4, 1.0, n_samples).astype(int),
        }
        df = pd.DataFrame(data)

        # Scope 3 as function of other variables (with noise)
        df['scope3_total'] = (
            df['scope1_tco2'] * 3.5 +
            df['scope2_tco2'] * 2.0 +
            df['revenue_musd'] * 0.8 +
            df['employees'] * 0.05 +
            np.random.lognormal(6, 0.5, n_samples)
        )

        # Category breakdown (approximate percentages)
        df['cat1_purchased'] = df['scope3_total'] * np.random.uniform(0.3, 0.5, n_samples)
        df['cat4_transport'] = df['scope3_total'] * np.random.uniform(0.05, 0.15, n_samples)
        df['cat11_use'] = df['scope3_total'] * np.random.uniform(0.15, 0.35, n_samples)
        df['cat12_eol'] = df['scope3_total'] * np.random.uniform(0.02, 0.08, n_samples)

        return df

    def train_models(self, df: pd.DataFrame) -> Dict:
        features = ['revenue_musd', 'employees', 'scope1_tco2',
                     'scope2_tco2', 'industry_code', 'supplier_count']
        X = df[features]

        targets = ['scope3_total', 'cat1_purchased', 'cat4_transport',
                    'cat11_use', 'cat12_eol']

        results = {}
        for target in targets:
            y = df[target]
            rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
            gb = GradientBoostingRegressor(n_estimators=100, random_state=42)

            rf_scores = cross_val_score(rf, X, y, cv=5, scoring='r2')
            gb_scores = cross_val_score(gb, X, y, cv=5, scoring='r2')

            best_model = rf if rf_scores.mean() > gb_scores.mean() else gb
            best_model.fit(X, y)
            self.models[target] = best_model

            results[target] = {
                'rf_r2_mean': float(rf_scores.mean()),
                'rf_r2_std': float(rf_scores.std()),
                'gb_r2_mean': float(gb_scores.mean()),
                'gb_r2_std': float(gb_scores.std()),
                'best_model': 'RandomForest' if rf_scores.mean() > gb_scores.mean() else 'GradientBoosting',
                'feature_importance': dict(zip(features,
                    best_model.feature_importances_.tolist()))
            }

        return results

    def estimate_ev_battery_scope3(self) -> Dict:
        """Estimate Scope 3 for a representative EV battery manufacturer."""
        company = {
            'revenue_musd': 5000,
            'employees': 15000,
            'scope1_tco2': 50000,
            'scope2_tco2': 120000,
            'industry_code': 29,
            'supplier_count': 500,
        }
        X = pd.DataFrame([company])
        estimates = {}
        for target, model in self.models.items():
            estimates[target] = float(model.predict(X)[0])
        return estimates


# ============================================================
# Module 6: Full Pipeline Runner
# ============================================================

class AutoLCAPipeline:
    """Orchestrates the full AutoLCA pipeline."""

    def __init__(self):
        self.tree_builder = NLPProcessTreeBuilder()
        self.matcher = EcoinventMatcher()
        self.uncertainty = UncertaintyPropagation(n_simulations=10000)
        self.hotspot = HotspotAnalyzer()
        self.scope3 = Scope3Estimator()
        self.results = {}

    def run(self) -> Dict:
        print("=" * 60)
        print("AutoLCA Pipeline - EV Battery Manufacturing Case Study")
        print("=" * 60)

        # Step 1: Build process tree
        print("\n[1/6] Building process tree...")
        graph = self.tree_builder.build_ev_battery_tree()
        processes = self.tree_builder.processes
        print(f"  -> {len(processes)} processes, {len(graph.edges())} edges")
        self.results['process_tree'] = {
            'n_processes': len(processes),
            'n_edges': len(graph.edges()),
            'processes': {k: v.name for k, v in processes.items()},
        }

        # Step 2: Ecoinvent matching
        print("\n[2/6] Auto-matching to Ecoinvent database...")
        matches = self.matcher.auto_match_process_tree(processes)
        match_stats = {
            'high_confidence': sum(1 for m in matches.values()
                                   if m['best_match'] and m['best_match']['confidence'] == 'high'),
            'medium_confidence': sum(1 for m in matches.values()
                                     if m['best_match'] and m['best_match']['confidence'] == 'medium'),
            'low_confidence': sum(1 for m in matches.values()
                                  if m['best_match'] and m['best_match']['confidence'] == 'low'),
        }
        print(f"  -> High: {match_stats['high_confidence']}, "
              f"Medium: {match_stats['medium_confidence']}, "
              f"Low: {match_stats['low_confidence']}")
        self.results['matching'] = {
            'statistics': match_stats,
            'details': {k: {
                'process': v['process_name'],
                'best_match': v['best_match']['match']['name'] if v['best_match'] else 'None',
                'similarity': v['best_match']['similarity'] if v['best_match'] else 0,
                'confidence': v['best_match']['confidence'] if v['best_match'] else 'none',
            } for k, v in matches.items()}
        }

        # Step 3: Uncertainty propagation
        print("\n[3/6] Running uncertainty analysis...")
        unc_results = self.uncertainty.compare_methods(processes, graph)
        print(f"  -> MC mean GWP: {unc_results['monte_carlo']['mean']:.2f} kg CO2-eq")
        print(f"  -> MC std: {unc_results['monte_carlo']['std']:.2f}")
        print(f"  -> Taylor mean: {unc_results['taylor']['mean']:.2f} kg CO2-eq")
        print(f"  -> Method diff (mean): {unc_results['relative_difference_mean']:.1f}%")
        self.results['uncertainty'] = unc_results

        # Step 4: Hotspot analysis
        print("\n[4/6] Performing hotspot analysis...")
        hotspot_df = self.hotspot.analyze(processes)
        print(f"  -> Top hotspot: {hotspot_df.iloc[0]['process_name']} "
              f"({hotspot_df.iloc[0]['contribution_pct']:.1f}%)")
        self.results['hotspot'] = hotspot_df.to_dict('records')

        # Step 5: Scenario comparison
        print("\n[5/6] Running scenario comparison...")
        scenarios = self.hotspot.scenario_comparison(processes)
        for name, data in scenarios.items():
            print(f"  -> {name}: {data['total_gwp']:.1f} kg CO2-eq")
        self.results['scenarios'] = {k: v['total_gwp'] for k, v in scenarios.items()}
        self.results['scenario_details'] = scenarios

        # Step 6: Scope 3 estimation
        print("\n[6/6] Training Scope 3 ML models...")
        training_data = self.scope3.generate_training_data()
        model_results = self.scope3.train_models(training_data)
        scope3_est = self.scope3.estimate_ev_battery_scope3()
        print(f"  -> Estimated total Scope 3: {scope3_est.get('scope3_total', 0):,.0f} tCO2-eq")
        for target, metrics in model_results.items():
            print(f"  -> {target}: R²={metrics[metrics['best_model'].lower().replace('forest','_forest').replace('boosting','_boosting').replace('random_forest','rf').replace('gradient_boosting','gb')+'_r2_mean'] if False else metrics['rf_r2_mean' if metrics['best_model']=='RandomForest' else 'gb_r2_mean']:.3f}")
        self.results['scope3'] = {
            'model_performance': model_results,
            'ev_battery_estimates': scope3_est,
        }

        # Store full MC samples for plotting
        mc_full = self.uncertainty.monte_carlo(processes, graph)
        self.results['mc_samples'] = mc_full['gwp_samples']
        self.results['mc_process_contributions'] = mc_full['process_contributions']

        print("\n" + "=" * 60)
        print("Pipeline complete!")
        print("=" * 60)

        return self.results


if __name__ == '__main__':
    pipeline = AutoLCAPipeline()
    results = pipeline.run()
    print(f"\nTotal processes analyzed: {results['process_tree']['n_processes']}")
    print(f"Total scenarios compared: {len(results['scenarios'])}")
