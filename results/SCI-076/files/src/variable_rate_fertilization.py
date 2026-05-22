"""
Module 5: 可変施肥マップの自動生成（クリギング＋最適化）
Variable rate fertilization map generation using Kriging + optimization
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pathlib import Path

FIGURES_DIR = Path(__file__).parent.parent / "figures"
RESULTS_DIR = Path(__file__).parent.parent / "results"


def yield_response_function(N_rate, ndvi, vwc, ec, ph):
    """
    Quadratic-plateau yield response to nitrogen, modulated by soil/plant factors.
    Based on rice-specific N response curves (Dobermann & Fairhurst, 2000).
    
    N_rate: kg N/ha
    ndvi: current NDVI (vegetation vigor indicator)
    vwc, ec, ph: soil properties
    """
    # Optimal N depends on growth potential
    N_opt = 80 + 40 * ndvi  # higher vigor → higher N demand
    
    # pH effect on N availability
    ph_factor = 1 - 0.15 * abs(ph - 6.2)
    
    # EC effect (salinity stress)
    ec_factor = np.where(ec > 1.2, 0.85, 1.0)
    
    # Quadratic-plateau model
    plateau = 6.5 * ph_factor * ec_factor
    a = plateau
    b = 0.04 * ph_factor
    c = -0.0002
    
    yield_est = a + b * N_rate + c * N_rate ** 2
    yield_est = np.minimum(yield_est, plateau * 1.1)
    yield_est = np.maximum(yield_est, 2.0)
    
    return yield_est


def optimize_n_rate(ndvi, vwc, ec, ph, n_price=1.2, rice_price=250):
    """
    Optimize N application rate for maximum economic return.
    
    n_price: cost per kg N (USD)
    rice_price: rice price per ton (USD)
    """
    def neg_profit(N_rate):
        y = yield_response_function(N_rate[0], ndvi, vwc, ec, ph)
        revenue = y * rice_price
        cost = N_rate[0] * n_price
        return -(revenue - cost)
    
    result = minimize(neg_profit, x0=[80], bounds=[(0, 200)], method='L-BFGS-B')
    return result.x[0]


def generate_vra_map(grid_size=50, field_size=500):
    """Generate variable rate application map."""
    np.random.seed(42)
    
    x = np.linspace(0, field_size, grid_size)
    y = np.linspace(0, field_size, grid_size)
    xx, yy = np.meshgrid(x, y)
    
    # Load or simulate spatial layers
    ndvi_map = 0.6 + 0.15 * np.sin(2*np.pi*xx/field_size) * np.cos(2*np.pi*yy/field_size) + \
               np.random.normal(0, 0.03, (grid_size, grid_size))
    ndvi_map = np.clip(ndvi_map, 0.3, 0.9)
    
    vwc_map = 0.35 + 0.08 * np.cos(np.pi*xx/field_size) + np.random.normal(0, 0.02, (grid_size, grid_size))
    vwc_map = np.clip(vwc_map, 0.15, 0.55)
    
    ec_map = 0.8 + 0.3 * np.sin(3*np.pi*xx/field_size) + np.random.normal(0, 0.05, (grid_size, grid_size))
    ec_map = np.clip(ec_map, 0.2, 2.0)
    
    ph_map = 5.8 + 0.5 * np.cos(2*np.pi*yy/field_size) + np.random.normal(0, 0.1, (grid_size, grid_size))
    ph_map = np.clip(ph_map, 4.5, 7.0)
    
    # Optimize N rate for each grid cell
    n_rate_map = np.zeros((grid_size, grid_size))
    yield_uniform = np.zeros((grid_size, grid_size))
    yield_vra = np.zeros((grid_size, grid_size))
    
    uniform_rate = 80.0  # Standard uniform N rate (kg/ha)
    
    for i in range(grid_size):
        for j in range(grid_size):
            opt_n = optimize_n_rate(ndvi_map[i, j], vwc_map[i, j],
                                     ec_map[i, j], ph_map[i, j])
            n_rate_map[i, j] = opt_n
            yield_vra[i, j] = yield_response_function(opt_n, ndvi_map[i, j],
                                                       vwc_map[i, j], ec_map[i, j], ph_map[i, j])
            yield_uniform[i, j] = yield_response_function(uniform_rate, ndvi_map[i, j],
                                                            vwc_map[i, j], ec_map[i, j], ph_map[i, j])
    
    # Categorize into management zones
    n_zones = 4
    zone_edges = np.percentile(n_rate_map, np.linspace(0, 100, n_zones + 1))
    zone_map = np.digitize(n_rate_map, zone_edges[1:-1])
    
    results = {
        'n_rate_map': n_rate_map,
        'yield_vra': yield_vra,
        'yield_uniform': yield_uniform,
        'zone_map': zone_map,
        'ndvi_map': ndvi_map,
        'uniform_rate': uniform_rate,
    }
    
    np.savez(RESULTS_DIR / "vra_maps.npz", **{k: v for k, v in results.items() if isinstance(v, np.ndarray)})
    
    return results


def run_vra_analysis():
    """Run VRA map generation and analysis."""
    results = generate_vra_map()
    
    n_rate_map = results['n_rate_map']
    yield_vra = results['yield_vra']
    yield_uniform = results['yield_uniform']
    zone_map = results['zone_map']
    uniform_rate = results['uniform_rate']
    
    # --- Figure 9: VRA Prescription Map ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    ax = axes[0, 0]
    im = ax.imshow(n_rate_map, cmap='RdYlGn_r', origin='lower', extent=[0, 500, 0, 500])
    ax.set_title('Optimized N Rate (kg/ha)', fontsize=12)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    plt.colorbar(im, ax=ax, shrink=0.8)
    
    ax = axes[0, 1]
    zone_colors = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad']
    im = ax.imshow(zone_map, cmap=plt.cm.get_cmap('RdYlGn_r', 4), origin='lower',
                   extent=[0, 500, 0, 500])
    ax.set_title('Management Zones (4 classes)', fontsize=12)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, ticks=[0, 1, 2, 3])
    cbar.set_ticklabels(['Low', 'Med-Low', 'Med-High', 'High'])
    
    ax = axes[1, 0]
    yield_diff = yield_vra - yield_uniform
    im = ax.imshow(yield_diff, cmap='RdBu', origin='lower', extent=[0, 500, 0, 500],
                   vmin=-0.5, vmax=0.5)
    ax.set_title('Yield Gain: VRA vs Uniform (t/ha)', fontsize=12)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    plt.colorbar(im, ax=ax, shrink=0.8)
    
    ax = axes[1, 1]
    ax.hist(n_rate_map.ravel(), bins=25, color='steelblue', edgecolor='navy', alpha=0.7, label='VRA')
    ax.axvline(uniform_rate, color='red', linestyle='--', linewidth=2, label=f'Uniform ({uniform_rate} kg/ha)')
    ax.axvline(n_rate_map.mean(), color='green', linestyle='--', linewidth=2,
               label=f'VRA mean ({n_rate_map.mean():.1f} kg/ha)')
    ax.set_xlabel('N Rate (kg/ha)')
    ax.set_ylabel('Frequency')
    ax.set_title('N Rate Distribution', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Variable Rate Fertilization Prescription Map', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig09_vra_prescription.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Figure 10: Economic Analysis ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    
    rice_price = 250  # USD/ton
    n_price = 1.2  # USD/kg N
    
    profit_vra = yield_vra * rice_price - n_rate_map * n_price
    profit_uniform = yield_uniform * rice_price - uniform_rate * n_price
    profit_diff = profit_vra - profit_uniform
    
    ax = axes[0]
    im = ax.imshow(profit_diff, cmap='RdBu', origin='lower', extent=[0, 500, 0, 500])
    ax.set_title('Profit Gain: VRA vs Uniform (USD/ha)', fontsize=12)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    plt.colorbar(im, ax=ax, shrink=0.8)
    
    ax = axes[1]
    categories = ['Uniform\nApplication', 'Variable Rate\nApplication']
    mean_yields = [yield_uniform.mean(), yield_vra.mean()]
    mean_profits = [profit_uniform.mean(), profit_vra.mean()]
    
    x_pos = np.arange(len(categories))
    width = 0.35
    bars1 = ax.bar(x_pos - width/2, mean_yields, width, label='Yield (t/ha)', color='#2ecc71', edgecolor='black')
    ax2 = ax.twinx()
    bars2 = ax2.bar(x_pos + width/2, mean_profits, width, label='Profit (USD/ha)', color='#3498db', edgecolor='black')
    
    ax.set_ylabel('Yield (t/ha)')
    ax2.set_ylabel('Profit (USD/ha)')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(categories)
    ax.set_title('Economic Comparison', fontsize=12)
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Economic Analysis of Variable Rate Fertilization', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig10_economic_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Summary stats
    total_n_uniform = uniform_rate * 25  # 25 ha field
    total_n_vra = n_rate_map.mean() * 25
    n_savings = (1 - total_n_vra / total_n_uniform) * 100
    
    print(f"\n=== Variable Rate Fertilization Results ===")
    print(f"Uniform N rate: {uniform_rate:.0f} kg/ha")
    print(f"VRA mean N rate: {n_rate_map.mean():.1f} kg/ha (range: {n_rate_map.min():.1f}-{n_rate_map.max():.1f})")
    print(f"N fertilizer savings: {n_savings:.1f}%")
    print(f"Yield - Uniform: {yield_uniform.mean():.2f} t/ha")
    print(f"Yield - VRA: {yield_vra.mean():.2f} t/ha")
    print(f"Yield gain: {(yield_vra.mean() - yield_uniform.mean()):.3f} t/ha")
    print(f"Profit gain: {profit_diff.mean():.1f} USD/ha")
    
    return {
        'n_savings_pct': round(n_savings, 1),
        'yield_uniform': round(yield_uniform.mean(), 2),
        'yield_vra': round(yield_vra.mean(), 2),
        'profit_gain': round(profit_diff.mean(), 1),
        'n_rate_mean': round(n_rate_map.mean(), 1),
        'n_rate_range': (round(n_rate_map.min(), 1), round(n_rate_map.max(), 1)),
    }


if __name__ == "__main__":
    run_vra_analysis()
