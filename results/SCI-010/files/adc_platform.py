#!/usr/bin/env python3
"""
ADC Payload-Linker Optimization Computational Platform
=====================================================
Implements ODE-based PK/PD models and Monte Carlo simulations for
Antibody-Drug Conjugate optimization, with HER2-targeted T-DXd analog case study.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import json
import os
from datetime import datetime

sns.set_theme(style="whitegrid", font_scale=1.1)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'DejaVu Sans'

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
DATA_DIR = "data"
LOGS_DIR = "logs"

# ============================================================================
# 1. DAR Distribution Model
# ============================================================================

@dataclass
class DARParameters:
    """Parameters for DAR distribution modeling."""
    mean_dar: float = 4.0          # Target average DAR
    max_dar: int = 8               # Maximum possible DAR (IgG1 has 8 interchain cysteines)
    conjugation_efficiency: float = 0.85
    batch_variability: float = 0.15


class DARDistributionModel:
    """Models DAR distribution using binomial/Poisson mixture and relates to therapeutic window."""

    def __init__(self, params: DARParameters = None):
        self.params = params or DARParameters()

    def binomial_dar_distribution(self, p_conjugation: float, n_sites: int = 8) -> np.ndarray:
        """Binomial DAR distribution based on conjugation probability per site."""
        from scipy.stats import binom
        dar_values = np.arange(0, n_sites + 1)
        probs = binom.pmf(dar_values, n_sites, p_conjugation)
        return dar_values, probs

    def poisson_dar_distribution(self, mean_dar: float, max_dar: int = 8) -> np.ndarray:
        """Poisson-approximated DAR distribution (for stochastic conjugation)."""
        from scipy.stats import poisson
        dar_values = np.arange(0, max_dar + 1)
        probs = poisson.pmf(dar_values, mean_dar)
        probs = probs / probs.sum()  # Normalize for truncation
        return dar_values, probs

    def monte_carlo_dar_sampling(self, n_molecules: int = 10000,
                                  n_batches: int = 50) -> Dict:
        """Monte Carlo simulation of DAR distribution across batches."""
        np.random.seed(42)
        all_dar = []
        batch_means = []
        batch_stds = []

        for batch in range(n_batches):
            # Vary conjugation efficiency per batch
            eff = np.clip(
                np.random.normal(self.params.conjugation_efficiency,
                                 self.params.batch_variability),
                0.1, 1.0
            )
            p_site = eff * self.params.mean_dar / self.params.max_dar
            p_site = np.clip(p_site, 0, 1)
            dar_batch = np.random.binomial(self.params.max_dar, p_site, n_molecules)
            all_dar.extend(dar_batch)
            batch_means.append(np.mean(dar_batch))
            batch_stds.append(np.std(dar_batch))

        return {
            'all_dar': np.array(all_dar),
            'batch_means': np.array(batch_means),
            'batch_stds': np.array(batch_stds),
            'overall_mean': np.mean(all_dar),
            'overall_std': np.std(all_dar),
        }

    def therapeutic_window_model(self, dar_values: np.ndarray) -> Dict:
        """Model relationship between DAR and efficacy/toxicity."""
        # Efficacy: sigmoidal increase with DAR
        ec50_dar = 3.5
        emax = 1.0
        hill_eff = 2.5
        efficacy = emax * dar_values**hill_eff / (ec50_dar**hill_eff + dar_values**hill_eff)

        # Toxicity: exponential increase at high DAR
        tox_threshold = 4.0
        tox_slope = 0.8
        toxicity = 1.0 / (1.0 + np.exp(-tox_slope * (dar_values - tox_threshold - 2)))

        # Aggregation propensity increases with DAR
        aggregation = 0.01 * np.exp(0.4 * dar_values)

        # Clearance rate increases with DAR (hydrophobicity)
        clearance_factor = 1.0 + 0.15 * (dar_values - 2)**2 * (dar_values > 2)

        # Therapeutic index
        ti = efficacy / (toxicity + 0.01)

        return {
            'efficacy': efficacy,
            'toxicity': toxicity,
            'aggregation': aggregation,
            'clearance_factor': clearance_factor,
            'therapeutic_index': ti,
        }

    def plot_dar_analysis(self, mc_results: Dict):
        """Generate comprehensive DAR analysis plots."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))

        # (a) DAR distribution histogram
        ax = axes[0, 0]
        dar_counts = np.bincount(mc_results['all_dar'], minlength=9)
        dar_fracs = dar_counts / dar_counts.sum()
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, 9))
        ax.bar(range(9), dar_fracs, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_xlabel('DAR Value')
        ax.set_ylabel('Fraction')
        ax.set_title(f'(a) DAR Distribution (Mean={mc_results["overall_mean"]:.2f}±{mc_results["overall_std"]:.2f})')
        ax.set_xticks(range(9))

        # (b) Batch-to-batch variability
        ax = axes[0, 1]
        ax.errorbar(range(len(mc_results['batch_means'])),
                     mc_results['batch_means'],
                     yerr=mc_results['batch_stds'],
                     fmt='o', markersize=4, capsize=3, color='#2196F3',
                     ecolor='#90CAF9')
        ax.axhline(y=self.params.mean_dar, color='red', linestyle='--',
                    label=f'Target DAR={self.params.mean_dar}')
        ax.fill_between(range(len(mc_results['batch_means'])),
                         self.params.mean_dar - 0.5, self.params.mean_dar + 0.5,
                         alpha=0.15, color='red', label='±0.5 spec')
        ax.set_xlabel('Batch Number')
        ax.set_ylabel('Mean DAR')
        ax.set_title('(b) Batch-to-Batch DAR Variability')
        ax.legend(fontsize=9)

        # (c) Therapeutic window
        ax = axes[1, 0]
        dar_cont = np.linspace(0, 8, 200)
        tw = self.therapeutic_window_model(dar_cont)
        ax.plot(dar_cont, tw['efficacy'], 'g-', linewidth=2, label='Efficacy')
        ax.plot(dar_cont, tw['toxicity'], 'r-', linewidth=2, label='Toxicity')
        ax.fill_between(dar_cont, 0, 1, where=(tw['efficacy'] > 0.5) & (tw['toxicity'] < 0.3),
                         alpha=0.2, color='green', label='Therapeutic Window')
        ax.set_xlabel('DAR')
        ax.set_ylabel('Normalized Response')
        ax.set_title('(c) DAR–Efficacy/Toxicity Relationship')
        ax.legend(fontsize=9)

        # (d) Therapeutic index
        ax = axes[1, 1]
        ax.plot(dar_cont, tw['therapeutic_index'], 'b-', linewidth=2)
        optimal_dar = dar_cont[np.argmax(tw['therapeutic_index'])]
        max_ti = np.max(tw['therapeutic_index'])
        ax.axvline(x=optimal_dar, color='red', linestyle='--',
                    label=f'Optimal DAR={optimal_dar:.1f}')
        ax.annotate(f'TI_max={max_ti:.1f}', xy=(optimal_dar, max_ti),
                     xytext=(optimal_dar + 1.5, max_ti * 0.9),
                     arrowprops=dict(arrowstyle='->', color='red'),
                     fontsize=10, color='red')
        ax.set_xlabel('DAR')
        ax.set_ylabel('Therapeutic Index')
        ax.set_title('(d) Therapeutic Index vs DAR')
        ax.legend(fontsize=9)

        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'fig1_dar_analysis.png')
        plt.savefig(path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {path}")
        return path


# ============================================================================
# 2. Linker Cleavage Mechanism Simulation
# ============================================================================

@dataclass
class LinkerParameters:
    """Parameters for different linker types."""
    # Acid-sensitive (hydrazone)
    acid_k_tumor: float = 0.15     # h^-1, tumor pH ~6.5
    acid_k_plasma: float = 0.005   # h^-1, plasma pH 7.4
    acid_k_lysosome: float = 0.8   # h^-1, lysosome pH ~5.0

    # Enzyme-cleavable (Val-Cit, GGFG)
    enzyme_vmax: float = 0.5       # h^-1
    enzyme_km: float = 5.0         # μM
    cathepsin_conc_tumor: float = 10.0   # μM
    cathepsin_conc_plasma: float = 0.1   # μM

    # Disulfide (reducible)
    disulfide_k_intracellular: float = 0.3  # h^-1 (GSH ~10 mM)
    disulfide_k_plasma: float = 0.01        # h^-1 (GSH ~2 μM)


class LinkerCleavageSimulator:
    """Simulates linker cleavage kinetics for different mechanisms."""

    def __init__(self, params: LinkerParameters = None):
        self.params = params or LinkerParameters()

    def acid_cleavage_ode(self, t, y, pH):
        """ODE for acid-sensitive linker cleavage."""
        intact, released = y
        # pH-dependent rate: k = k0 * 10^(7.4 - pH)
        k_base = self.params.acid_k_plasma
        k = k_base * 10**(7.4 - pH)
        d_intact = -k * intact
        d_released = k * intact
        return [d_intact, d_released]

    def enzyme_cleavage_ode(self, t, y, cathepsin_conc):
        """ODE for enzyme-cleavable linker (Michaelis-Menten)."""
        intact, released = y
        v = self.params.enzyme_vmax * cathepsin_conc / (self.params.enzyme_km + cathepsin_conc)
        d_intact = -v * intact
        d_released = v * intact
        return [d_intact, d_released]

    def disulfide_cleavage_ode(self, t, y, gsh_conc):
        """ODE for disulfide linker cleavage."""
        intact, released = y
        # Rate proportional to GSH concentration
        k = self.params.disulfide_k_intracellular * (gsh_conc / 10.0)
        d_intact = -k * intact
        d_released = k * intact
        return [d_intact, d_released]

    def simulate_all_linkers(self, t_span=(0, 72), n_points=500):
        """Simulate cleavage for all linker types in different environments."""
        t_eval = np.linspace(t_span[0], t_span[1], n_points)
        y0 = [1.0, 0.0]
        results = {}

        environments = {
            'Plasma (pH 7.4)': {'pH': 7.4, 'cathepsin': 0.1, 'gsh': 0.002},
            'Tumor ECM (pH 6.5)': {'pH': 6.5, 'cathepsin': 2.0, 'gsh': 0.5},
            'Endosome (pH 5.5)': {'pH': 5.5, 'cathepsin': 5.0, 'gsh': 1.0},
            'Lysosome (pH 4.5)': {'pH': 4.5, 'cathepsin': 10.0, 'gsh': 10.0},
        }

        for env_name, env_params in environments.items():
            results[env_name] = {}

            # Acid-sensitive
            sol = solve_ivp(self.acid_cleavage_ode, t_span, y0,
                           args=(env_params['pH'],), t_eval=t_eval, method='RK45')
            results[env_name]['acid'] = {'t': sol.t, 'released': sol.y[1]}

            # Enzyme-cleavable
            sol = solve_ivp(self.enzyme_cleavage_ode, t_span, y0,
                           args=(env_params['cathepsin'],), t_eval=t_eval, method='RK45')
            results[env_name]['enzyme'] = {'t': sol.t, 'released': sol.y[1]}

            # Disulfide
            sol = solve_ivp(self.disulfide_cleavage_ode, t_span, y0,
                           args=(env_params['gsh'],), t_eval=t_eval, method='RK45')
            results[env_name]['disulfide'] = {'t': sol.t, 'released': sol.y[1]}

        return results

    def plot_linker_cleavage(self, results: Dict):
        """Plot linker cleavage kinetics."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))
        linker_types = ['acid', 'enzyme', 'disulfide']
        linker_labels = ['Acid-Sensitive (Hydrazone)', 'Enzyme-Cleavable (Val-Cit)',
                         'Disulfide (Reducible)']
        colors = ['#E53935', '#1E88E5', '#43A047']
        linestyles = ['-', '--', '-.', ':']
        env_names = list(results.keys())

        # Panel per linker type
        for i, (lt, ll) in enumerate(zip(linker_types, linker_labels)):
            ax = axes[i // 2, i % 2]
            for j, env in enumerate(env_names):
                data = results[env][lt]
                ax.plot(data['t'], data['released'] * 100,
                        color=colors[i], linestyle=linestyles[j],
                        linewidth=2, label=env, alpha=0.7 + 0.1 * j)
            ax.set_xlabel('Time (h)')
            ax.set_ylabel('Payload Released (%)')
            ax.set_title(f'({chr(97+i)}) {ll}')
            ax.legend(fontsize=8)
            ax.set_ylim(0, 105)

        # Panel (d): Selectivity ratio
        ax = axes[1, 1]
        t_points = [6, 12, 24, 48, 72]
        bar_width = 0.25
        x = np.arange(len(t_points))

        for i, (lt, ll, c) in enumerate(zip(linker_types, ['Acid', 'Enzyme', 'Disulfide'], colors)):
            selectivity = []
            for tp in t_points:
                lyso_data = results['Lysosome (pH 4.5)'][lt]
                plasma_data = results['Plasma (pH 7.4)'][lt]
                idx = np.argmin(np.abs(lyso_data['t'] - tp))
                lyso_rel = lyso_data['released'][idx]
                plasma_rel = plasma_data['released'][idx]
                ratio = lyso_rel / (plasma_rel + 1e-6)
                selectivity.append(min(ratio, 200))
            ax.bar(x + i * bar_width, selectivity, bar_width, label=ll, color=c, alpha=0.8)

        ax.set_xlabel('Time (h)')
        ax.set_ylabel('Selectivity Ratio (Lysosome/Plasma)')
        ax.set_title('(d) Cleavage Selectivity')
        ax.set_xticks(x + bar_width)
        ax.set_xticklabels([f'{t}h' for t in t_points])
        ax.legend(fontsize=9)
        ax.set_yscale('log')

        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'fig2_linker_cleavage.png')
        plt.savefig(path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {path}")
        return path


# ============================================================================
# 3. Bystander Effect Model (Tumor Diffusion)
# ============================================================================

@dataclass
class BustanderParameters:
    """Parameters for bystander effect diffusion model."""
    D_payload: float = 1e-7        # cm²/s, diffusion coefficient
    k_uptake: float = 0.05         # s^-1, cellular uptake rate
    k_efflux: float = 0.005        # s^-1, efflux rate
    k_kill: float = 0.01           # s^-1, cell killing rate
    tumor_radius: float = 0.05     # cm (500 μm)
    cell_spacing: float = 0.001    # cm (10 μm)
    antigen_pos_fraction: float = 0.7  # fraction of Ag+ cells
    membrane_permeability: float = 1e-5  # cm/s


class BustanderEffectModel:
    """2D reaction-diffusion model for bystander killing effect."""

    def __init__(self, params: BustanderParameters = None):
        self.params = params or BustanderParameters()

    def solve_1d_diffusion(self, t_max: float = 3600, nx: int = 100, nt: int = 500):
        """Solve 1D reaction-diffusion PDE for payload spread."""
        L = self.params.tumor_radius * 2
        dx = L / nx
        dt = t_max / nt
        x = np.linspace(0, L, nx)

        # CFL condition check
        D = self.params.D_payload
        cfl = D * dt / dx**2
        if cfl > 0.5:
            dt = 0.4 * dx**2 / D
            nt = int(t_max / dt) + 1

        # Initialize: payload released at center (Ag+ cell that internalized ADC)
        c_extra = np.zeros(nx)   # Extracellular payload
        c_intra = np.zeros(nx)   # Intracellular payload
        viability = np.ones(nx)  # Cell viability

        # Ag+ cells: clustered near center with some random positions
        np.random.seed(42)
        ag_positive = np.zeros(nx, dtype=bool)
        center = nx // 2
        ag_range = int(nx * self.params.antigen_pos_fraction * 0.6)
        ag_positive[center - ag_range:center + ag_range] = True
        # Add some scattered Ag+ cells
        random_pos = np.random.choice(nx, size=int(nx * 0.1), replace=False)
        ag_positive[random_pos] = True

        # Source: Ag+ cells release payload after ADC internalization
        source_rate = np.zeros(nx)
        source_rate[ag_positive] = 0.001  # Continuous release from internalized ADC

        snapshots = {'time': [], 'x': x * 1e4, 'c_extra': [], 'c_intra': [],
                     'viability': [], 'ag_positive': ag_positive}
        save_times = [0, t_max * 0.1, t_max * 0.25, t_max * 0.5, t_max]

        t = 0
        for step in range(nt):
            # Diffusion (explicit finite difference)
            c_new = c_extra.copy()
            c_new[1:-1] += D * dt / dx**2 * (c_extra[2:] - 2 * c_extra[1:-1] + c_extra[:-2])

            # Source from Ag+ cells
            c_new += source_rate * dt

            # Cellular uptake and efflux
            uptake = self.params.k_uptake * c_new * dt
            efflux = self.params.k_efflux * c_intra * dt
            c_extra_new = c_new - uptake + efflux
            c_intra += uptake - efflux

            # Cell killing (depends on intracellular concentration)
            kill_prob = self.params.k_kill * c_intra * dt
            viability *= (1 - np.clip(kill_prob, 0, 0.99))

            c_extra = np.clip(c_extra_new, 0, None)
            t += dt

            if any(abs(t - st) < dt * 1.5 for st in save_times):
                snapshots['time'].append(t)
                snapshots['c_extra'].append(c_extra.copy())
                snapshots['c_intra'].append(c_intra.copy())
                snapshots['viability'].append(viability.copy())

        return snapshots

    def simulate_permeability_comparison(self):
        """Compare bystander effect for different payload membrane permeabilities."""
        permeabilities = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4]
        results = {}

        for perm in permeabilities:
            self.params.membrane_permeability = perm
            self.params.k_uptake = perm * 100
            snap = self.solve_1d_diffusion()
            final_viability = snap['viability'][-1] if snap['viability'] else np.ones(100)
            results[f'{perm:.0e}'] = {
                'viability': final_viability,
                'x': snap['x'],
                'mean_kill': 1 - np.mean(final_viability),
                'bystander_kill': 1 - np.mean(final_viability[~snap['ag_positive']]),
            }
        return results

    def plot_bystander_effect(self, snapshots: Dict, perm_results: Dict):
        """Plot bystander effect analysis."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))

        # (a) Extracellular payload diffusion over time
        ax = axes[0, 0]
        colors_time = plt.cm.plasma(np.linspace(0.1, 0.9, len(snapshots['time'])))
        for i, (t, c) in enumerate(zip(snapshots['time'], snapshots['c_extra'])):
            ax.plot(snapshots['x'], c, color=colors_time[i],
                    linewidth=2, label=f't={t/60:.0f} min')
        # Mark Ag+ region
        ag_mask = snapshots['ag_positive']
        ax.fill_between(snapshots['x'], 0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.01,
                         where=ag_mask, alpha=0.1, color='green', label='Ag+ cells')
        ax.set_xlabel('Position (μm)')
        ax.set_ylabel('Extracellular [Payload] (a.u.)')
        ax.set_title('(a) Payload Diffusion Profile')
        ax.legend(fontsize=8)

        # (b) Cell viability over time
        ax = axes[0, 1]
        for i, (t, v) in enumerate(zip(snapshots['time'], snapshots['viability'])):
            ax.plot(snapshots['x'], v * 100, color=colors_time[i],
                    linewidth=2, label=f't={t/60:.0f} min')
        ax.fill_between(snapshots['x'], 0, 100,
                         where=ag_mask, alpha=0.1, color='green')
        ax.set_xlabel('Position (μm)')
        ax.set_ylabel('Cell Viability (%)')
        ax.set_title('(b) Spatial Viability Profile')
        ax.legend(fontsize=8)
        ax.set_ylim(0, 105)

        # (c) Bystander vs target cell killing
        ax = axes[1, 0]
        if len(snapshots['viability']) > 0:
            final_v = snapshots['viability'][-1]
            target_kill = (1 - np.mean(final_v[ag_mask])) * 100
            bystander_kill = (1 - np.mean(final_v[~ag_mask])) * 100

            categories = ['Target (Ag+)\nCells', 'Bystander (Ag-)\nCells', 'Overall']
            kills = [target_kill, bystander_kill,
                     (1 - np.mean(final_v)) * 100]
            bars = ax.bar(categories, kills,
                          color=['#2196F3', '#FF9800', '#4CAF50'],
                          edgecolor='black', linewidth=0.5)
            for bar, val in zip(bars, kills):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f'{val:.1f}%', ha='center', fontsize=10)
        ax.set_ylabel('Cell Killing (%)')
        ax.set_title('(c) Target vs Bystander Killing')
        ax.set_ylim(0, 100)

        # (d) Permeability effect on bystander killing
        ax = axes[1, 1]
        perms = list(perm_results.keys())
        bystander_kills = [perm_results[p]['bystander_kill'] * 100 for p in perms]
        total_kills = [perm_results[p]['mean_kill'] * 100 for p in perms]
        x_pos = np.arange(len(perms))
        ax.bar(x_pos - 0.15, total_kills, 0.3, label='Total Kill', color='#2196F3', alpha=0.8)
        ax.bar(x_pos + 0.15, bystander_kills, 0.3, label='Bystander Kill', color='#FF9800', alpha=0.8)
        ax.set_xlabel('Membrane Permeability (cm/s)')
        ax.set_ylabel('Cell Killing (%)')
        ax.set_title('(d) Effect of Payload Permeability')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(perms, rotation=45, fontsize=8)
        ax.legend(fontsize=9)

        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'fig3_bystander_effect.png')
        plt.savefig(path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {path}")
        return path


# ============================================================================
# 4. Plasma Stability vs Tumor Release Optimization
# ============================================================================

class StabilityReleaseOptimizer:
    """Optimizes balance between plasma stability and intratumoral release."""

    def __init__(self):
        self.optimization_history = []

    def objective_function(self, params, verbose=False):
        """
        Multi-objective: maximize tumor release, minimize plasma release.
        params: [k_cleavage_base, pH_sensitivity, enzyme_sensitivity, hydrophobicity]
        """
        k_base, pH_sens, enz_sens, hydrophob = params

        # Plasma stability (want high = low release)
        k_plasma = k_base * 10**(pH_sens * (7.4 - 7.4))  # pH 7.4
        k_plasma += enz_sens * 0.1  # Low cathepsin in plasma
        plasma_release_24h = 1 - np.exp(-k_plasma * 24)

        # Tumor release (want high)
        k_tumor = k_base * 10**(pH_sens * (7.4 - 6.0))  # pH ~6.0 lysosome-like
        k_tumor += enz_sens * 10.0  # High cathepsin in tumor
        tumor_release_24h = 1 - np.exp(-k_tumor * 24)

        # Clearance penalty (hydrophobicity increases clearance)
        clearance_penalty = 0.3 * hydrophob**2

        # Aggregation penalty
        agg_penalty = 0.2 * max(0, hydrophob - 0.5)**2

        # Objective: maximize selectivity while maintaining efficacy
        selectivity = tumor_release_24h / (plasma_release_24h + 0.01)
        score = -(selectivity * tumor_release_24h - clearance_penalty - agg_penalty)

        if verbose:
            return {
                'score': -score,
                'plasma_release': plasma_release_24h * 100,
                'tumor_release': tumor_release_24h * 100,
                'selectivity': selectivity,
                'clearance_penalty': clearance_penalty,
            }
        return score

    def optimize(self, n_iterations: int = 200):
        """Run differential evolution optimization."""
        bounds = [(0.001, 0.5), (0.1, 3.0), (0.01, 1.0), (0.1, 1.0)]
        result = differential_evolution(
            self.objective_function, bounds,
            maxiter=n_iterations, seed=42, tol=1e-8,
            popsize=30
        )

        optimal_result = self.objective_function(result.x, verbose=True)
        return {
            'optimal_params': {
                'k_cleavage_base': result.x[0],
                'pH_sensitivity': result.x[1],
                'enzyme_sensitivity': result.x[2],
                'hydrophobicity': result.x[3],
            },
            'optimal_metrics': optimal_result,
            'convergence': result.fun,
        }

    def parameter_sensitivity_analysis(self, n_samples: int = 5000):
        """Monte Carlo sensitivity analysis of linker parameters."""
        np.random.seed(42)
        param_names = ['k_base', 'pH_sens', 'enz_sens', 'hydrophob']
        bounds = [(0.001, 0.5), (0.1, 3.0), (0.01, 1.0), (0.1, 1.0)]

        samples = np.column_stack([
            np.random.uniform(low, high, n_samples)
            for low, high in bounds
        ])

        scores = []
        plasma_releases = []
        tumor_releases = []

        for s in samples:
            result = self.objective_function(s, verbose=True)
            scores.append(result['score'])
            plasma_releases.append(result['plasma_release'])
            tumor_releases.append(result['tumor_release'])

        return {
            'samples': samples,
            'param_names': param_names,
            'scores': np.array(scores),
            'plasma_releases': np.array(plasma_releases),
            'tumor_releases': np.array(tumor_releases),
        }

    def plot_optimization(self, opt_result: Dict, sensitivity: Dict):
        """Plot optimization results."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))

        # (a) Pareto front: plasma stability vs tumor release
        ax = axes[0, 0]
        sc = ax.scatter(sensitivity['plasma_releases'],
                        sensitivity['tumor_releases'],
                        c=sensitivity['scores'], cmap='viridis',
                        alpha=0.4, s=10, edgecolors='none')
        plt.colorbar(sc, ax=ax, label='Objective Score')
        # Mark optimal
        opt_metrics = opt_result['optimal_metrics']
        ax.scatter([opt_metrics['plasma_release']], [opt_metrics['tumor_release']],
                   c='red', s=200, marker='*', zorder=5, label='Optimal')
        ax.set_xlabel('Plasma Release at 24h (%)')
        ax.set_ylabel('Tumor Release at 24h (%)')
        ax.set_title('(a) Stability-Release Pareto Space')
        ax.legend(fontsize=9)

        # (b) Sensitivity tornado chart
        ax = axes[0, 1]
        from scipy.stats import spearmanr
        correlations = []
        for i, name in enumerate(sensitivity['param_names']):
            r, _ = spearmanr(sensitivity['samples'][:, i], sensitivity['scores'])
            correlations.append(r)
        sorted_idx = np.argsort(np.abs(correlations))
        names_sorted = [sensitivity['param_names'][i] for i in sorted_idx]
        corr_sorted = [correlations[i] for i in sorted_idx]
        colors = ['#E53935' if c < 0 else '#1E88E5' for c in corr_sorted]
        ax.barh(names_sorted, corr_sorted, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_xlabel('Spearman Correlation with Score')
        ax.set_title('(b) Parameter Sensitivity (Tornado)')
        ax.axvline(x=0, color='black', linewidth=0.5)

        # (c) Release kinetics for optimal linker
        ax = axes[1, 0]
        opt_p = opt_result['optimal_params']
        t_array = np.linspace(0, 72, 200)

        k_plasma = opt_p['k_cleavage_base']
        k_plasma += opt_p['enzyme_sensitivity'] * 0.1
        plasma_curve = (1 - np.exp(-k_plasma * t_array)) * 100

        k_tumor = opt_p['k_cleavage_base'] * 10**(opt_p['pH_sensitivity'] * 1.4)
        k_tumor += opt_p['enzyme_sensitivity'] * 10.0
        tumor_curve = (1 - np.exp(-k_tumor * t_array)) * 100

        ax.plot(t_array, tumor_curve, 'g-', linewidth=2, label='Tumor (pH 6.0)')
        ax.plot(t_array, plasma_curve, 'r--', linewidth=2, label='Plasma (pH 7.4)')
        ax.fill_between(t_array, plasma_curve, tumor_curve, alpha=0.15, color='green')
        ax.set_xlabel('Time (h)')
        ax.set_ylabel('Payload Released (%)')
        ax.set_title('(c) Optimized Release Kinetics')
        ax.legend(fontsize=9)

        # (d) Optimal parameters radar chart
        ax = axes[1, 1]
        params = opt_result['optimal_params']
        param_labels = ['k_base\n(Cleavage)', 'pH\nSensitivity', 'Enzyme\nSensitivity', 'Hydro-\nphobicity']
        param_values = list(params.values())
        # Normalize to [0, 1]
        bounds_arr = np.array([(0.001, 0.5), (0.1, 3.0), (0.01, 1.0), (0.1, 1.0)])
        norm_values = [(v - b[0]) / (b[1] - b[0]) for v, b in zip(param_values, bounds_arr)]
        norm_values.append(norm_values[0])

        angles = np.linspace(0, 2 * np.pi, len(param_labels), endpoint=False).tolist()
        angles.append(angles[0])

        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect('equal')

        # Draw radar manually
        for i, (angle, label) in enumerate(zip(angles[:-1], param_labels)):
            ax.plot([0, np.cos(angle)], [0, np.sin(angle)], 'gray', linewidth=0.5)
            ax.text(1.2 * np.cos(angle), 1.2 * np.sin(angle), label,
                    ha='center', va='center', fontsize=9)

        # Plot values
        radar_x = [v * np.cos(a) for v, a in zip(norm_values, angles)]
        radar_y = [v * np.sin(a) for v, a in zip(norm_values, angles)]
        ax.fill(radar_x, radar_y, alpha=0.25, color='#2196F3')
        ax.plot(radar_x, radar_y, 'o-', color='#2196F3', linewidth=2, markersize=6)
        ax.set_title('(d) Optimal Parameter Profile')
        ax.axis('off')

        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'fig4_optimization.png')
        plt.savefig(path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {path}")
        return path


# ============================================================================
# 5. PK Model Integration
# ============================================================================

@dataclass
class PKParameters:
    """Two-compartment PK model parameters for ADC."""
    # Central compartment
    V1: float = 3.0            # L, central volume
    CL: float = 0.01           # L/h, clearance
    # Peripheral compartment
    V2: float = 4.0            # L, peripheral volume
    Q: float = 0.02            # L/h, intercompartmental clearance
    # ADC-specific
    k_deconj: float = 0.005    # h^-1, deconjugation rate
    k_internalization: float = 0.05  # h^-1
    # Payload
    CL_payload: float = 0.5    # L/h
    V_payload: float = 50.0    # L
    # Tumor
    k_tumor_uptake: float = 0.002   # h^-1
    k_tumor_release: float = 0.01   # h^-1 (payload release in tumor)
    tumor_volume: float = 0.01      # L
    # Target-mediated disposition
    k_on: float = 1.0          # nM^-1 h^-1
    k_off: float = 0.01        # h^-1
    R0: float = 10.0           # nM, receptor concentration
    k_int: float = 0.1         # h^-1, internalization of complex


class PKModel:
    """Two-compartment PK model with target-mediated drug disposition (TMDD)."""

    def __init__(self, params: PKParameters = None):
        self.params = params or PKParameters()

    def pk_ode_system(self, t, y, dose_times=None, dose_amount=0):
        """
        Full PK/PD ODE system.
        States:
        y[0] = ADC in central (nM)
        y[1] = ADC in peripheral (nM)
        y[2] = Free payload in plasma (nM)
        y[3] = ADC in tumor (nM)
        y[4] = Payload in tumor (nM)
        y[5] = Free receptor (nM)
        y[6] = ADC-receptor complex (nM)
        y[7] = Tumor cell fraction
        """
        p = self.params
        adc_c, adc_p, payload_p, adc_t, payload_t, R, AR, tumor = y

        # Target-mediated disposition
        binding = p.k_on * adc_c * R
        unbinding = p.k_off * AR
        internalization = p.k_int * AR

        # ADC in central compartment
        d_adc_c = (-p.CL / p.V1 * adc_c           # Clearance
                   - p.Q / p.V1 * (adc_c - adc_p)  # Distribution
                   - p.k_deconj * adc_c             # Deconjugation
                   - binding + unbinding             # TMDD
                   - p.k_tumor_uptake * adc_c)       # Tumor uptake

        # ADC in peripheral
        d_adc_p = p.Q / p.V2 * (adc_c - adc_p) - p.CL / p.V2 * 0.3 * adc_p

        # Free payload in plasma (from deconjugation)
        d_payload_p = (p.k_deconj * adc_c * 4.0     # DAR=4 average
                       - p.CL_payload / p.V_payload * payload_p)

        # ADC in tumor
        d_adc_t = (p.k_tumor_uptake * adc_c * (p.V1 / p.tumor_volume)
                   - p.k_internalization * adc_t
                   - p.k_tumor_release * adc_t)

        # Payload in tumor
        d_payload_t = (p.k_internalization * adc_t * 4.0  # DAR=4
                       + p.k_tumor_release * adc_t * 2.0
                       - 0.05 * payload_t)   # Payload clearance from tumor

        # Receptor dynamics
        k_syn = p.k_int * p.R0  # Synthesis rate (steady state)
        d_R = k_syn - p.k_int * R - binding + unbinding
        d_AR = binding - unbinding - internalization

        # Tumor dynamics (logistic growth + drug killing)
        k_growth = 0.003  # h^-1
        k_kill = 0.0005   # h^-1 nM^-1
        d_tumor = k_growth * tumor * (1 - tumor) - k_kill * payload_t * tumor

        return [d_adc_c, d_adc_p, d_payload_p, d_adc_t, d_payload_t, d_R, d_AR, d_tumor]

    def simulate_dosing_regimen(self, dose_mg_kg: float = 5.4,
                                 interval_h: float = 504,  # 3 weeks
                                 n_doses: int = 6,
                                 body_weight: float = 70):
        """Simulate multiple-dose PK with typical T-DXd regimen."""
        mw_adc = 150000  # Da
        dose_mg = dose_mg_kg * body_weight
        dose_nmol = dose_mg * 1e6 / mw_adc
        dose_nM = dose_nmol / self.params.V1  # nM in central compartment

        t_total = interval_h * n_doses + interval_h
        t_eval = np.linspace(0, t_total, 5000)

        # Initial conditions
        y0 = [0, 0, 0, 0, 0, self.params.R0, 0, 1.0]

        # Simulate with event-driven dosing
        all_t = []
        all_y = []
        current_y = np.array(y0)

        for dose_n in range(n_doses):
            t_start = dose_n * interval_h
            t_end = (dose_n + 1) * interval_h if dose_n < n_doses - 1 else t_total

            # Add dose
            current_y[0] += dose_nM

            t_span = (t_start, t_end)
            t_eval_segment = np.linspace(t_start, t_end, 1000)

            sol = solve_ivp(self.pk_ode_system, t_span, current_y,
                           t_eval=t_eval_segment, method='LSODA',
                           rtol=1e-8, atol=1e-10)

            all_t.extend(sol.t)
            all_y.append(sol.y)
            current_y = sol.y[:, -1].copy()

        t_array = np.array(all_t)
        y_array = np.hstack(all_y)

        return {
            't': t_array,
            'adc_central': y_array[0],
            'adc_peripheral': y_array[1],
            'payload_plasma': y_array[2],
            'adc_tumor': y_array[3],
            'payload_tumor': y_array[4],
            'free_receptor': y_array[5],
            'adc_receptor': y_array[6],
            'tumor_fraction': y_array[7],
            'dose_times': [i * interval_h for i in range(n_doses)],
            'dose_nM': dose_nM,
        }

    def dose_response_simulation(self, doses_mg_kg: List[float] = None):
        """Simulate dose-response for different dose levels."""
        if doses_mg_kg is None:
            doses_mg_kg = [1.6, 3.2, 5.4, 6.4, 8.0]

        results = {}
        for dose in doses_mg_kg:
            sim = self.simulate_dosing_regimen(dose_mg_kg=dose, n_doses=4)
            final_tumor = sim['tumor_fraction'][-1]
            max_payload_plasma = np.max(sim['payload_plasma'])
            auc_payload_tumor = np.trapz(sim['payload_tumor'], sim['t'])

            results[dose] = {
                'final_tumor': final_tumor,
                'max_payload_plasma': max_payload_plasma,
                'auc_payload_tumor': auc_payload_tumor,
                'tumor_response': (1 - final_tumor) * 100,
            }
        return results

    def plot_pk_simulation(self, sim_result: Dict, dose_response: Dict):
        """Plot PK simulation results."""
        fig, axes = plt.subplots(2, 3, figsize=(18, 11))

        t_days = sim_result['t'] / 24

        # (a) ADC concentration in plasma
        ax = axes[0, 0]
        ax.semilogy(t_days, sim_result['adc_central'], 'b-', linewidth=2, label='Central')
        ax.semilogy(t_days, sim_result['adc_peripheral'], 'b--', linewidth=1.5, label='Peripheral')
        for dt in sim_result['dose_times']:
            ax.axvline(x=dt/24, color='red', linestyle=':', alpha=0.5)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('ADC Concentration (nM)')
        ax.set_title('(a) ADC Pharmacokinetics')
        ax.legend(fontsize=9)

        # (b) Payload in plasma vs tumor
        ax = axes[0, 1]
        ax.plot(t_days, sim_result['payload_plasma'], 'r-', linewidth=2, label='Plasma')
        ax.plot(t_days, sim_result['payload_tumor'], 'g-', linewidth=2, label='Tumor')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Payload Concentration (nM)')
        ax.set_title('(b) Payload Distribution')
        ax.legend(fontsize=9)

        # (c) Receptor occupancy
        ax = axes[0, 2]
        total_receptor = sim_result['free_receptor'] + sim_result['adc_receptor']
        occupancy = sim_result['adc_receptor'] / (total_receptor + 1e-10) * 100
        ax.plot(t_days, occupancy, 'purple', linewidth=2)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('HER2 Receptor Occupancy (%)')
        ax.set_title('(c) Target Engagement')
        ax.set_ylim(0, 105)

        # (d) Tumor growth inhibition
        ax = axes[1, 0]
        ax.plot(t_days, sim_result['tumor_fraction'] * 100, 'k-', linewidth=2)
        ax.axhline(y=100, color='gray', linestyle='--', alpha=0.5, label='Baseline')
        ax.axhline(y=30, color='green', linestyle='--', alpha=0.5, label='PR threshold')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Tumor Size (% baseline)')
        ax.set_title('(d) Tumor Growth Inhibition')
        ax.legend(fontsize=9)

        # (e) Dose-response
        ax = axes[1, 1]
        doses = sorted(dose_response.keys())
        responses = [dose_response[d]['tumor_response'] for d in doses]
        max_payloads = [dose_response[d]['max_payload_plasma'] for d in doses]

        ax.bar(range(len(doses)), responses, color='#4CAF50', edgecolor='black',
               linewidth=0.5, alpha=0.8)
        ax.set_xlabel('Dose (mg/kg)')
        ax.set_ylabel('Tumor Response (%)')
        ax.set_title('(e) Dose-Response Relationship')
        ax.set_xticks(range(len(doses)))
        ax.set_xticklabels([f'{d}' for d in doses])

        ax2 = ax.twinx()
        ax2.plot(range(len(doses)), max_payloads, 'ro-', linewidth=2, label='Max Plasma Payload')
        ax2.set_ylabel('Max Plasma Payload (nM)', color='red')
        ax2.tick_params(axis='y', labelcolor='red')

        # (f) Therapeutic index by dose
        ax = axes[1, 2]
        ti_values = []
        for d in doses:
            efficacy = dose_response[d]['tumor_response'] / 100
            toxicity_proxy = dose_response[d]['max_payload_plasma'] / max(max_payloads)
            ti = efficacy / (toxicity_proxy + 0.01)
            ti_values.append(ti)

        colors = ['#4CAF50' if ti > np.median(ti_values) else '#FF9800' for ti in ti_values]
        ax.bar(range(len(doses)), ti_values, color=colors, edgecolor='black', linewidth=0.5)
        optimal_dose = doses[np.argmax(ti_values)]
        ax.set_xlabel('Dose (mg/kg)')
        ax.set_ylabel('Therapeutic Index (a.u.)')
        ax.set_title(f'(f) Therapeutic Index (Optimal: {optimal_dose} mg/kg)')
        ax.set_xticks(range(len(doses)))
        ax.set_xticklabels([f'{d}' for d in doses])

        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'fig5_pk_simulation.png')
        plt.savefig(path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {path}")
        return path


# ============================================================================
# 6. HER2-Targeted ADC Case Study (T-DXd Analog)
# ============================================================================

class TDXdCaseStudy:
    """Integrated case study for HER2-targeted ADC (T-DXd-like)."""

    def __init__(self):
        # T-DXd specific parameters
        self.dar_params = DARParameters(mean_dar=8.0, max_dar=8,
                                         conjugation_efficiency=0.95,
                                         batch_variability=0.05)
        self.linker_params = LinkerParameters(
            enzyme_vmax=0.6,       # GGFG tetrapeptide linker
            enzyme_km=3.0,
            cathepsin_conc_tumor=12.0,
        )
        self.pk_params = PKParameters(
            V1=2.83, CL=0.0088, V2=3.5, Q=0.015,
            k_deconj=0.003,        # Low deconjugation (stable linker)
            k_internalization=0.06,
            k_tumor_release=0.015,  # DXd release rate
        )

    def run_comprehensive_analysis(self):
        """Run full analysis pipeline."""
        results = {}

        # 1. DAR analysis (T-DXd has DAR ~8)
        print("=" * 60)
        print("1. DAR Distribution Analysis (T-DXd: DAR≈8)")
        print("=" * 60)
        dar_model = DARDistributionModel(self.dar_params)
        mc_dar = dar_model.monte_carlo_dar_sampling(n_molecules=20000, n_batches=100)
        results['dar'] = mc_dar
        print(f"  Mean DAR: {mc_dar['overall_mean']:.2f} ± {mc_dar['overall_std']:.2f}")

        # Compare with conventional ADC (DAR 2-4)
        conventional_params = DARParameters(mean_dar=3.5, max_dar=8,
                                             conjugation_efficiency=0.7,
                                             batch_variability=0.2)
        dar_conv = DARDistributionModel(conventional_params)
        mc_conv = dar_conv.monte_carlo_dar_sampling(n_molecules=20000, n_batches=100)
        results['dar_conventional'] = mc_conv
        print(f"  Conventional Mean DAR: {mc_conv['overall_mean']:.2f} ± {mc_conv['overall_std']:.2f}")

        # 2. Linker analysis (GGFG peptide)
        print("\n" + "=" * 60)
        print("2. Linker Cleavage Analysis (GGFG peptide)")
        print("=" * 60)
        linker_sim = LinkerCleavageSimulator(self.linker_params)
        linker_results = linker_sim.simulate_all_linkers()
        results['linker'] = linker_results

        # 3. PK simulation
        print("\n" + "=" * 60)
        print("3. PK/PD Simulation (T-DXd 5.4 mg/kg Q3W)")
        print("=" * 60)
        pk_model = PKModel(self.pk_params)
        pk_sim = pk_model.simulate_dosing_regimen(dose_mg_kg=5.4, n_doses=6)
        results['pk'] = pk_sim
        print(f"  Peak ADC conc: {np.max(pk_sim['adc_central']):.1f} nM")
        print(f"  Final tumor: {pk_sim['tumor_fraction'][-1]*100:.1f}% of baseline")

        # 4. Dose comparison
        dose_response = pk_model.dose_response_simulation([1.6, 3.2, 5.4, 6.4, 8.0])
        results['dose_response'] = dose_response

        return results

    def plot_case_study(self, results: Dict):
        """Generate comprehensive case study figure."""
        fig, axes = plt.subplots(2, 3, figsize=(18, 11))

        # (a) DAR comparison: T-DXd vs conventional
        ax = axes[0, 0]
        dar_tdxd = np.bincount(results['dar']['all_dar'], minlength=9)
        dar_conv = np.bincount(results['dar_conventional']['all_dar'], minlength=9)
        x = np.arange(9)
        ax.bar(x - 0.2, dar_tdxd / dar_tdxd.sum(), 0.4,
               label=f'T-DXd (DAR≈{results["dar"]["overall_mean"]:.1f})',
               color='#2196F3', alpha=0.8)
        ax.bar(x + 0.2, dar_conv / dar_conv.sum(), 0.4,
               label=f'Conventional (DAR≈{results["dar_conventional"]["overall_mean"]:.1f})',
               color='#FF9800', alpha=0.8)
        ax.set_xlabel('DAR')
        ax.set_ylabel('Fraction')
        ax.set_title('(a) DAR Distribution Comparison')
        ax.legend(fontsize=9)

        # (b) T-DXd therapeutic window
        ax = axes[0, 1]
        dar_model = DARDistributionModel(self.dar_params)
        dar_cont = np.linspace(0, 8, 200)
        tw = dar_model.therapeutic_window_model(dar_cont)
        ax.plot(dar_cont, tw['efficacy'], 'g-', linewidth=2, label='Efficacy')
        ax.plot(dar_cont, tw['toxicity'], 'r-', linewidth=2, label='Toxicity')
        ax.plot(dar_cont, tw['aggregation'], 'b--', linewidth=1.5, label='Aggregation')
        ax.axvline(x=results['dar']['overall_mean'], color='purple',
                   linestyle=':', linewidth=2, label=f'T-DXd DAR={results["dar"]["overall_mean"]:.1f}')
        ax.set_xlabel('DAR')
        ax.set_ylabel('Normalized Score')
        ax.set_title('(b) T-DXd Therapeutic Window')
        ax.legend(fontsize=8)

        # (c) DXd release kinetics
        ax = axes[0, 2]
        pk_t = results['pk']['t'] / 24
        ax.plot(pk_t, results['pk']['payload_tumor'], 'g-', linewidth=2, label='Tumor DXd')
        ax.plot(pk_t, results['pk']['payload_plasma'], 'r--', linewidth=1.5, label='Plasma DXd')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('DXd Concentration (nM)')
        ax.set_title('(c) DXd Release Profile')
        ax.legend(fontsize=9)

        # (d) ADC PK profile
        ax = axes[1, 0]
        ax.semilogy(pk_t, results['pk']['adc_central'], 'b-', linewidth=2, label='Intact ADC')
        for dt in results['pk']['dose_times']:
            ax.axvline(x=dt/24, color='red', linestyle=':', alpha=0.4)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('ADC Concentration (nM)')
        ax.set_title('(d) T-DXd Plasma PK (5.4 mg/kg Q3W)')
        ax.legend(fontsize=9)

        # (e) Tumor growth inhibition
        ax = axes[1, 1]
        tumor = np.clip(results['pk']['tumor_fraction'] * 100, 0, 200)
        ax.plot(pk_t, tumor, 'k-', linewidth=2)
        # Add RECIST criteria lines
        ax.axhline(y=100, color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=70, color='orange', linestyle='--', alpha=0.5, label='SD boundary')
        ax.axhline(y=30, color='green', linestyle='--', alpha=0.5, label='PR boundary')
        ax.fill_between(pk_t, 0, 30, alpha=0.1, color='green')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Tumor Size (% baseline)')
        ax.set_title('(e) Tumor Growth Inhibition')
        ax.legend(fontsize=9)

        # (f) Dose-response with therapeutic index
        ax = axes[1, 2]
        doses = sorted(results['dose_response'].keys())
        responses = [results['dose_response'][d]['tumor_response'] for d in doses]
        ax.bar(range(len(doses)), responses, color='#4CAF50', edgecolor='black',
               linewidth=0.5, alpha=0.8)
        ax.axhline(y=70, color='green', linestyle='--', alpha=0.5, label='PR threshold')
        ax.set_xlabel('Dose (mg/kg)')
        ax.set_ylabel('Tumor Response (%)')
        ax.set_title('(f) T-DXd Dose-Response')
        ax.set_xticks(range(len(doses)))
        ax.set_xticklabels([f'{d}' for d in doses])
        ax.legend(fontsize=9)

        plt.suptitle('HER2-Targeted ADC (T-DXd Analog) — Integrated Analysis',
                     fontsize=14, fontweight='bold', y=1.01)
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, 'fig6_tdxd_case_study.png')
        plt.savefig(path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {path}")
        return path

    def generate_summary_table(self, results: Dict) -> pd.DataFrame:
        """Generate summary metrics table."""
        metrics = {
            'Parameter': [
                'Mean DAR (T-DXd)', 'Mean DAR (Conventional)',
                'DAR CV% (T-DXd)', 'DAR CV% (Conventional)',
                'Peak ADC Conc (nM)', 'Final Tumor Size (%)',
                'Tumor Response 1.6 mg/kg (%)', 'Tumor Response 3.2 mg/kg (%)',
                'Tumor Response 5.4 mg/kg (%)', 'Tumor Response 6.4 mg/kg (%)',
                'Tumor Response 8.0 mg/kg (%)',
            ],
            'Value': [
                f"{results['dar']['overall_mean']:.2f}",
                f"{results['dar_conventional']['overall_mean']:.2f}",
                f"{results['dar']['overall_std']/results['dar']['overall_mean']*100:.1f}",
                f"{results['dar_conventional']['overall_std']/results['dar_conventional']['overall_mean']*100:.1f}",
                f"{np.max(results['pk']['adc_central']):.1f}",
                f"{results['pk']['tumor_fraction'][-1]*100:.1f}",
            ] + [f"{results['dose_response'][d]['tumor_response']:.1f}"
                 for d in sorted(results['dose_response'].keys())],
        }
        df = pd.DataFrame(metrics)
        return df


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Execute full ADC optimization platform analysis."""
    timestamp = datetime.now().isoformat()
    log_entries = []

    def log_event(phase, event_type, **kwargs):
        entry = {
            'timestamp': datetime.now().isoformat(),
            'phase': phase,
            'event_type': event_type,
            'actor': 'co-scientist',
            'skill_or_tool': 'adc-platform',
            **kwargs
        }
        log_entries.append(entry)

    log_event('init', 'run_started')
    print("=" * 70)
    print("ADC Payload-Linker Optimization Computational Platform")
    print("=" * 70)
    print(f"Timestamp: {timestamp}\n")

    # ---- Module 1: DAR Distribution ----
    print("\n>>> MODULE 1: DAR Distribution Analysis")
    dar_model = DARDistributionModel()
    mc_results = dar_model.monte_carlo_dar_sampling(n_molecules=20000, n_batches=100)
    fig1_path = dar_model.plot_dar_analysis(mc_results)
    log_event('module1', 'file_written', files_written=[fig1_path])

    # Save DAR data
    dar_df = pd.DataFrame({
        'batch': range(len(mc_results['batch_means'])),
        'mean_dar': mc_results['batch_means'],
        'std_dar': mc_results['batch_stds'],
    })
    dar_df.to_csv(os.path.join(DATA_DIR, 'dar_batch_statistics.csv'), index=False)

    # ---- Module 2: Linker Cleavage ----
    print("\n>>> MODULE 2: Linker Cleavage Simulation")
    linker_sim = LinkerCleavageSimulator()
    linker_results = linker_sim.simulate_all_linkers()
    fig2_path = linker_sim.plot_linker_cleavage(linker_results)
    log_event('module2', 'file_written', files_written=[fig2_path])

    # ---- Module 3: Bystander Effect ----
    print("\n>>> MODULE 3: Bystander Effect Model")
    bystander = BustanderEffectModel()
    snapshots = bystander.solve_1d_diffusion(t_max=3600, nx=100, nt=2000)
    perm_results = bystander.simulate_permeability_comparison()
    fig3_path = bystander.plot_bystander_effect(snapshots, perm_results)
    log_event('module3', 'file_written', files_written=[fig3_path])

    # Save bystander metrics
    perm_df = pd.DataFrame([
        {'permeability': k, 'mean_kill_pct': v['mean_kill']*100,
         'bystander_kill_pct': v['bystander_kill']*100}
        for k, v in perm_results.items()
    ])
    perm_df.to_csv(os.path.join(RESULTS_DIR, 'bystander_permeability.csv'), index=False)

    # ---- Module 4: Stability-Release Optimization ----
    print("\n>>> MODULE 4: Stability-Release Optimization")
    optimizer = StabilityReleaseOptimizer()
    opt_result = optimizer.optimize()
    sensitivity = optimizer.parameter_sensitivity_analysis(n_samples=5000)
    fig4_path = optimizer.plot_optimization(opt_result, sensitivity)
    log_event('module4', 'file_written', files_written=[fig4_path])

    # Save optimization results
    opt_df = pd.DataFrame([opt_result['optimal_params']])
    opt_df.to_csv(os.path.join(RESULTS_DIR, 'optimal_linker_params.csv'), index=False)

    print(f"\n  Optimal parameters: {opt_result['optimal_params']}")
    print(f"  Plasma release (24h): {opt_result['optimal_metrics']['plasma_release']:.1f}%")
    print(f"  Tumor release (24h): {opt_result['optimal_metrics']['tumor_release']:.1f}%")
    print(f"  Selectivity: {opt_result['optimal_metrics']['selectivity']:.1f}x")

    # ---- Module 5: PK Model ----
    print("\n>>> MODULE 5: PK/PD Simulation")
    pk_model = PKModel()
    pk_sim = pk_model.simulate_dosing_regimen(dose_mg_kg=5.4, n_doses=6)
    dose_response = pk_model.dose_response_simulation()
    fig5_path = pk_model.plot_pk_simulation(pk_sim, dose_response)
    log_event('module5', 'file_written', files_written=[fig5_path])

    # Save PK data
    pk_df = pd.DataFrame({
        'time_h': pk_sim['t'],
        'adc_central_nM': pk_sim['adc_central'],
        'payload_plasma_nM': pk_sim['payload_plasma'],
        'payload_tumor_nM': pk_sim['payload_tumor'],
        'tumor_fraction': pk_sim['tumor_fraction'],
    })
    pk_df.to_csv(os.path.join(DATA_DIR, 'pk_simulation_data.csv'), index=False)

    # ---- Module 6: T-DXd Case Study ----
    print("\n>>> MODULE 6: T-DXd Case Study")
    case_study = TDXdCaseStudy()
    cs_results = case_study.run_comprehensive_analysis()
    fig6_path = case_study.plot_case_study(cs_results)
    summary_table = case_study.generate_summary_table(cs_results)
    log_event('module6', 'file_written', files_written=[fig6_path])

    # Save summary
    summary_table.to_csv(os.path.join(RESULTS_DIR, 'tdxd_summary_metrics.csv'), index=False)
    print("\n  Summary Metrics:")
    print(summary_table.to_string(index=False))

    # ---- Save process log ----
    log_event('final', 'run_completed', status='ok')
    log_path = os.path.join(LOGS_DIR, 'process-log.jsonl')
    with open(log_path, 'w') as f:
        for entry in log_entries:
            f.write(json.dumps(entry, ensure_ascii=False) + '\n')
    print(f"\nProcess log saved: {log_path}")

    # ---- Collect all results for report generation ----
    all_results = {
        'mc_dar': mc_results,
        'opt_result': opt_result,
        'pk_sim_summary': {
            'peak_adc': float(np.max(pk_sim['adc_central'])),
            'final_tumor_pct': float(pk_sim['tumor_fraction'][-1] * 100),
        },
        'dose_response': {str(k): v for k, v in dose_response.items()},
        'bystander_metrics': {k: {'mean_kill': v['mean_kill'], 'bystander_kill': v['bystander_kill']}
                              for k, v in perm_results.items()},
        'summary_table': summary_table.to_dict(),
    }

    with open(os.path.join(RESULTS_DIR, 'all_results.json'), 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\n" + "=" * 70)
    print("Analysis Complete! All figures and data saved.")
    print("=" * 70)

    return all_results


if __name__ == '__main__':
    results = main()
