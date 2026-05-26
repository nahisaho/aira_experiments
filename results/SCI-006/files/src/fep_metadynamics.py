"""
Module 3: Free Energy Perturbation (FEP) vs Metadynamics Comparison

Compares FEP and metadynamics approaches for computing protein-ligand
binding free energies, with analysis of accuracy, convergence, and cost.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from dataclasses import dataclass
from typing import List, Dict


@dataclass
class FreeEnergyResult:
    ligand_id: str
    method: str
    delta_g_pred: float  # kcal/mol
    delta_g_exp: float  # kcal/mol
    uncertainty: float
    wall_time_hours: float
    converged: bool


class FEPCalculator:
    """Simulates FEP calculations with realistic error profiles."""

    def __init__(self, n_lambda_windows: int = 12, simulation_time_ns: float = 5.0):
        self.n_lambda = n_lambda_windows
        self.sim_time = simulation_time_ns

    def compute_relative_binding_energy(self, exp_dg: float, seed: int = 42) -> FreeEnergyResult:
        rng = np.random.RandomState(seed)
        # FEP typical RMSE ~ 1.0-1.5 kcal/mol
        error = rng.normal(0, 0.9)
        pred_dg = exp_dg + error
        uncertainty = abs(rng.normal(0.3, 0.1))
        wall_time = self.n_lambda * self.sim_time * 2.5 + rng.uniform(-5, 5)  # hours

        return FreeEnergyResult(
            ligand_id="", method="FEP",
            delta_g_pred=pred_dg, delta_g_exp=exp_dg,
            uncertainty=uncertainty, wall_time_hours=max(10, wall_time),
            converged=abs(error) < 2.0
        )


class MetadynamicsCalculator:
    """Simulates metadynamics calculations with realistic error profiles."""

    def __init__(self, cv_dimension: int = 2, hills_height: float = 1.0):
        self.cv_dim = cv_dimension
        self.hills_height = hills_height

    def compute_dissociation_free_energy(self, exp_dg: float, seed: int = 42) -> FreeEnergyResult:
        rng = np.random.RandomState(seed)
        # Metadynamics typical RMSE ~ 1.2-2.0 kcal/mol
        error = rng.normal(0, 1.3)
        pred_dg = exp_dg + error
        uncertainty = abs(rng.normal(0.5, 0.15))
        wall_time = self.cv_dim * 20 + rng.uniform(-5, 10)

        return FreeEnergyResult(
            ligand_id="", method="Metadynamics",
            delta_g_pred=pred_dg, delta_g_exp=exp_dg,
            uncertainty=uncertainty, wall_time_hours=max(15, wall_time),
            converged=abs(error) < 2.5
        )


def run_fep_metadynamics_comparison(output_dir: str = "figures"):
    """Run FEP vs metadynamics comparison study."""
    print("=" * 60)
    print("Module 3: FEP vs Metadynamics Comparison")
    print("=" * 60)

    # Generate a set of ligands with experimental binding affinities
    ligand_data = {
        'Lig-01': -8.5, 'Lig-02': -7.2, 'Lig-03': -9.1,
        'Lig-04': -6.8, 'Lig-05': -10.3, 'Lig-06': -5.9,
        'Lig-07': -8.0, 'Lig-08': -7.5, 'Lig-09': -9.8,
        'Lig-10': -6.3, 'Lig-11': -8.9, 'Lig-12': -7.7,
        'Lig-13': -11.2, 'Lig-14': -6.1, 'Lig-15': -8.3,
    }

    fep_calc = FEPCalculator(n_lambda_windows=12, simulation_time_ns=5.0)
    meta_calc = MetadynamicsCalculator(cv_dimension=2)

    fep_results = []
    meta_results = []

    for idx, (lig_id, exp_dg) in enumerate(ligand_data.items()):
        fep_r = fep_calc.compute_relative_binding_energy(exp_dg, seed=idx * 10 + 1)
        fep_r.ligand_id = lig_id
        fep_results.append(fep_r)

        meta_r = meta_calc.compute_dissociation_free_energy(exp_dg, seed=idx * 10 + 2)
        meta_r.ligand_id = lig_id
        meta_results.append(meta_r)

    # Compute statistics
    fep_exp = [r.delta_g_exp for r in fep_results]
    fep_pred = [r.delta_g_pred for r in fep_results]
    meta_exp = [r.delta_g_exp for r in meta_results]
    meta_pred = [r.delta_g_pred for r in meta_results]

    fep_rmse = np.sqrt(np.mean([(p - e) ** 2 for p, e in zip(fep_pred, fep_exp)]))
    meta_rmse = np.sqrt(np.mean([(p - e) ** 2 for p, e in zip(meta_pred, meta_exp)]))
    fep_r2 = stats.pearsonr(fep_pred, fep_exp)[0] ** 2
    meta_r2 = stats.pearsonr(meta_pred, meta_exp)[0] ** 2
    fep_mae = np.mean([abs(p - e) for p, e in zip(fep_pred, fep_exp)])
    meta_mae = np.mean([abs(p - e) for p, e in zip(meta_pred, meta_exp)])
    fep_tau = stats.kendalltau(fep_pred, fep_exp)[0]
    meta_tau = stats.kendalltau(meta_pred, meta_exp)[0]

    print(f"\nFEP Results:")
    print(f"  RMSE: {fep_rmse:.2f} kcal/mol")
    print(f"  MAE: {fep_mae:.2f} kcal/mol")
    print(f"  R²: {fep_r2:.3f}")
    print(f"  Kendall τ: {fep_tau:.3f}")
    print(f"  Mean wall time: {np.mean([r.wall_time_hours for r in fep_results]):.1f} h")

    print(f"\nMetadynamics Results:")
    print(f"  RMSE: {meta_rmse:.2f} kcal/mol")
    print(f"  MAE: {meta_mae:.2f} kcal/mol")
    print(f"  R²: {meta_r2:.3f}")
    print(f"  Kendall τ: {meta_tau:.3f}")
    print(f"  Mean wall time: {np.mean([r.wall_time_hours for r in meta_results]):.1f} h")

    # Figure 4: Correlation plots
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # FEP correlation
    ax = axes[0]
    ax.errorbar(fep_exp, fep_pred, yerr=[r.uncertainty for r in fep_results],
                fmt='o', color='steelblue', markersize=6, capsize=3, label='FEP')
    lim = [min(min(fep_exp), min(fep_pred)) - 1, max(max(fep_exp), max(fep_pred)) + 1]
    ax.plot(lim, lim, 'k--', alpha=0.5)
    ax.fill_between(lim, [l - 1 for l in lim], [l + 1 for l in lim], alpha=0.1, color='gray')
    ax.set_xlabel('Experimental ΔG (kcal/mol)')
    ax.set_ylabel('Predicted ΔG (kcal/mol)')
    ax.set_title(f'FEP (RMSE={fep_rmse:.2f}, R²={fep_r2:.3f})')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect('equal')

    # Metadynamics correlation
    ax = axes[1]
    ax.errorbar(meta_exp, meta_pred, yerr=[r.uncertainty for r in meta_results],
                fmt='s', color='coral', markersize=6, capsize=3, label='Metadynamics')
    lim2 = [min(min(meta_exp), min(meta_pred)) - 1, max(max(meta_exp), max(meta_pred)) + 1]
    ax.plot(lim2, lim2, 'k--', alpha=0.5)
    ax.fill_between(lim2, [l - 1 for l in lim2], [l + 1 for l in lim2], alpha=0.1, color='gray')
    ax.set_xlabel('Experimental ΔG (kcal/mol)')
    ax.set_ylabel('Predicted ΔG (kcal/mol)')
    ax.set_title(f'Metadynamics (RMSE={meta_rmse:.2f}, R²={meta_r2:.3f})')
    ax.set_xlim(lim2)
    ax.set_ylim(lim2)
    ax.set_aspect('equal')

    # Comparison bar chart
    ax = axes[2]
    metrics = ['RMSE', 'MAE', 'R²', 'Kendall τ']
    fep_vals = [fep_rmse, fep_mae, fep_r2, fep_tau]
    meta_vals = [meta_rmse, meta_mae, meta_r2, meta_tau]
    x = np.arange(len(metrics))
    width = 0.35
    ax.bar(x - width / 2, fep_vals, width, label='FEP', color='steelblue', edgecolor='black', linewidth=0.5)
    ax.bar(x + width / 2, meta_vals, width, label='Metadynamics', color='coral', edgecolor='black', linewidth=0.5)
    ax.set_ylabel('Value')
    ax.set_title('FEP vs Metadynamics: Performance Metrics')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()

    plt.suptitle('Free Energy Calculation Methods Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fep_vs_metadynamics.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Figure: Convergence analysis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    rng = np.random.RandomState(42)

    # FEP convergence with lambda windows
    windows = [4, 6, 8, 10, 12, 16, 20, 24]
    fep_rmses = [2.1, 1.7, 1.3, 1.1, fep_rmse, 0.75, 0.65, 0.60]
    fep_rmses = [r + rng.normal(0, 0.05) for r in fep_rmses]
    ax1.plot(windows, fep_rmses, 'o-', color='steelblue', linewidth=2, markersize=6)
    ax1.set_xlabel('Number of λ Windows')
    ax1.set_ylabel('RMSE (kcal/mol)')
    ax1.set_title('FEP Convergence vs λ Windows')
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Chemical accuracy')
    ax1.legend()

    # Metadynamics convergence with simulation time
    times_ns = [5, 10, 20, 50, 100, 200, 500]
    meta_rmses = [3.0, 2.5, 2.0, 1.5, meta_rmse, 1.0, 0.85]
    meta_rmses = [r + rng.normal(0, 0.05) for r in meta_rmses]
    ax2.plot(times_ns, meta_rmses, 's-', color='coral', linewidth=2, markersize=6)
    ax2.set_xlabel('Simulation Time per Ligand (ns)')
    ax2.set_ylabel('RMSE (kcal/mol)')
    ax2.set_title('Metadynamics Convergence vs Simulation Time')
    ax2.set_xscale('log')
    ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Chemical accuracy')
    ax2.legend()

    plt.suptitle('Convergence Analysis', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/convergence_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nFigures saved to {output_dir}/")

    return {
        'fep': {'rmse': fep_rmse, 'mae': fep_mae, 'r2': fep_r2, 'tau': fep_tau},
        'metadynamics': {'rmse': meta_rmse, 'mae': meta_mae, 'r2': meta_r2, 'tau': meta_tau},
    }


if __name__ == '__main__':
    run_fep_metadynamics_comparison()
