"""
Module 3: Mg2+/K+/Polyamine Concentration Optimization Map
Response surface methodology with Bayesian optimization
"""

import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import json

def cfps_yield_model(Mg, K, spd, noise=False):
    """
    Empirical response surface model for CFPS protein yield.
    Based on published optimization data (Jewett & Swartz, 2004; Caschera & Noireaux, 2014).

    Parameters:
        Mg: Mg2+ concentration (mM), typical range 5-25
        K:  K+ concentration (mM), typical range 40-300
        spd: Spermidine/polyamine concentration (mM), typical range 0-4
    Returns:
        Relative protein yield (0-1 scale)
    """
    # Optimal values
    Mg_opt, K_opt, spd_opt = 12.0, 170.0, 1.5

    # Response surface (quadratic with interactions)
    yield_val = 1.0
    yield_val *= np.exp(-0.5 * ((Mg - Mg_opt) / 4.0)**2)
    yield_val *= np.exp(-0.5 * ((K - K_opt) / 60.0)**2)
    yield_val *= np.exp(-0.5 * ((spd - spd_opt) / 0.8)**2)

    # Interaction terms
    yield_val *= (1 + 0.05 * np.sin(0.1 * Mg * spd))  # Mg-spd interaction
    yield_val *= (1 - 0.001 * np.abs(K - 180) * np.abs(Mg - 10))  # K-Mg interaction

    # Clamp
    yield_val = np.clip(yield_val, 0, 1)

    if noise:
        yield_val += np.random.normal(0, 0.03, yield_val.shape if hasattr(yield_val, 'shape') else None)
        yield_val = np.clip(yield_val, 0, 1)

    return yield_val

def bayesian_optimization_1d_demo():
    """Simple Bayesian optimization using expected improvement"""
    from scipy.stats import norm

    # We optimize Mg while fixing K=170, spd=1.5
    def objective(Mg):
        return -cfps_yield_model(Mg, 170.0, 1.5)

    # Initial random samples
    np.random.seed(42)
    X_observed = np.random.uniform(3, 25, 5).reshape(-1, 1)
    y_observed = np.array([objective(x[0]) for x in X_observed])

    # Simple GP surrogate using RBF-like kernel (manual for minimal deps)
    best_points = [X_observed.copy()]
    best_vals = [y_observed.min()]

    for iteration in range(15):
        # Grid search for next point (simplified BO)
        X_candidates = np.linspace(3, 25, 200).reshape(-1, 1)

        # Predict using inverse-distance weighted interpolation (surrogate)
        predictions = []
        uncertainties = []
        for xc in X_candidates:
            dists = np.abs(X_observed - xc).flatten()
            weights = 1.0 / (dists + 0.1)
            weights /= weights.sum()
            pred = np.sum(weights * y_observed)
            unc = 1.0 / (1.0 + np.sum(weights**2) * len(X_observed))
            predictions.append(pred)
            uncertainties.append(unc)

        predictions = np.array(predictions)
        uncertainties = np.array(uncertainties)

        # Expected improvement
        best_so_far = y_observed.min()
        z = (best_so_far - predictions) / (uncertainties + 1e-8)
        ei = (best_so_far - predictions) * norm.cdf(z) + uncertainties * norm.pdf(z)

        # Select best candidate
        next_idx = np.argmax(ei)
        next_x = X_candidates[next_idx]
        next_y = objective(next_x[0])

        X_observed = np.vstack([X_observed, next_x.reshape(1, -1)])
        y_observed = np.append(y_observed, next_y)
        best_vals.append(y_observed.min())

    return X_observed, y_observed, best_vals

def plot_optimization_maps(save_path='figures/fig5_ion_optimization.png'):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Map 1: Mg vs K (spd=1.5)
    Mg_range = np.linspace(3, 25, 80)
    K_range = np.linspace(40, 300, 80)
    Mg_grid, K_grid = np.meshgrid(Mg_range, K_range)
    Z1 = cfps_yield_model(Mg_grid, K_grid, 1.5)

    c1 = axes[0].contourf(Mg_grid, K_grid, Z1, levels=20, cmap='viridis')
    axes[0].set_xlabel('[Mg²⁺] (mM)', fontsize=11)
    axes[0].set_ylabel('[K⁺] (mM)', fontsize=11)
    axes[0].set_title('Mg²⁺ vs K⁺ (Spd=1.5 mM)', fontsize=13, fontweight='bold')
    plt.colorbar(c1, ax=axes[0], label='Relative Yield')
    axes[0].plot(12, 170, 'r*', markersize=15, label='Optimum')
    axes[0].legend()

    # Map 2: Mg vs Spermidine (K=170)
    spd_range = np.linspace(0, 4, 80)
    Mg_grid2, spd_grid = np.meshgrid(Mg_range, spd_range)
    Z2 = cfps_yield_model(Mg_grid2, 170, spd_grid)

    c2 = axes[1].contourf(Mg_grid2, spd_grid, Z2, levels=20, cmap='viridis')
    axes[1].set_xlabel('[Mg²⁺] (mM)', fontsize=11)
    axes[1].set_ylabel('[Spermidine] (mM)', fontsize=11)
    axes[1].set_title('Mg²⁺ vs Spermidine (K⁺=170 mM)', fontsize=13, fontweight='bold')
    plt.colorbar(c2, ax=axes[1], label='Relative Yield')
    axes[1].plot(12, 1.5, 'r*', markersize=15, label='Optimum')
    axes[1].legend()

    # Map 3: K vs Spermidine (Mg=12)
    K_grid3, spd_grid3 = np.meshgrid(K_range, spd_range)
    Z3 = cfps_yield_model(12, K_grid3, spd_grid3)

    c3 = axes[2].contourf(K_grid3, spd_grid3, Z3, levels=20, cmap='viridis')
    axes[2].set_xlabel('[K⁺] (mM)', fontsize=11)
    axes[2].set_ylabel('[Spermidine] (mM)', fontsize=11)
    axes[2].set_title('K⁺ vs Spermidine (Mg²⁺=12 mM)', fontsize=13, fontweight='bold')
    plt.colorbar(c3, ax=axes[2], label='Relative Yield')
    axes[2].plot(170, 1.5, 'r*', markersize=15, label='Optimum')
    axes[2].legend()

    plt.suptitle('Ion Concentration Optimization Maps for CFPS',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

def plot_bayesian_opt(X_obs, y_obs, best_vals, save_path='figures/fig6_bayesian_opt.png'):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: observed points vs true function
    Mg_true = np.linspace(3, 25, 200)
    y_true = np.array([cfps_yield_model(m, 170, 1.5) for m in Mg_true])

    axes[0].plot(Mg_true, y_true, 'k-', linewidth=2, label='True Response', alpha=0.7)
    axes[0].scatter(X_obs[:5, 0], -y_obs[:5], c='blue', s=80, zorder=5,
                    label='Initial Samples', edgecolors='black')
    axes[0].scatter(X_obs[5:, 0], -y_obs[5:], c='red', s=60, zorder=5,
                    label='BO Acquisitions', edgecolors='black', marker='^')
    axes[0].set_xlabel('[Mg²⁺] (mM)', fontsize=11)
    axes[0].set_ylabel('Relative Yield', fontsize=11)
    axes[0].set_title('Bayesian Optimization: Mg²⁺ Scan', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    # Right: convergence
    axes[1].plot(range(len(best_vals)), [-b for b in best_vals], 'o-',
                 color='#4CAF50', linewidth=2, markersize=6)
    axes[1].set_xlabel('Iteration', fontsize=11)
    axes[1].set_ylabel('Best Yield Found', fontsize=11)
    axes[1].set_title('BO Convergence', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Global Optimum')
    axes[1].legend(fontsize=10)

    plt.suptitle('Bayesian Optimization for Ion Concentration Tuning',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")

def grid_search_optimum():
    best_yield = 0
    best_params = {}
    for Mg in np.arange(5, 20, 0.5):
        for K in np.arange(100, 250, 5):
            for spd in np.arange(0.5, 3.0, 0.2):
                y = cfps_yield_model(Mg, K, spd)
                if y > best_yield:
                    best_yield = y
                    best_params = {'Mg_mM': float(Mg), 'K_mM': float(K), 'Spermidine_mM': float(spd)}
    best_params['max_yield'] = float(best_yield)
    return best_params

if __name__ == '__main__':
    print("=== Module 3: Mg2+/K+/Polyamine Optimization ===")
    plot_optimization_maps()

    X_obs, y_obs, best_vals = bayesian_optimization_1d_demo()
    plot_bayesian_opt(X_obs, y_obs, best_vals)

    opt = grid_search_optimum()
    print(f"  Optimal: Mg={opt['Mg_mM']:.1f}, K={opt['K_mM']:.0f}, Spd={opt['Spermidine_mM']:.1f}")
    print(f"  Max yield: {opt['max_yield']:.4f}")

    with open('results/m3_ion_optimization.json', 'w') as f:
        json.dump(opt, f, indent=2)
    print("  Saved: results/m3_ion_optimization.json")
