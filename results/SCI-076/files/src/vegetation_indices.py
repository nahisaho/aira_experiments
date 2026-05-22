"""
Module 1: 衛星/ドローンマルチスペクトル画像からの植生指数計算
Vegetation Index calculation from multispectral imagery
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path

FIGURES_DIR = Path(__file__).parent.parent / "figures"
RESULTS_DIR = Path(__file__).parent.parent / "results"
DATA_DIR = Path(__file__).parent.parent / "data"


def calculate_ndvi(nir: np.ndarray, red: np.ndarray) -> np.ndarray:
    """Normalized Difference Vegetation Index"""
    return np.where((nir + red) == 0, 0, (nir - red) / (nir + red))


def calculate_evi(nir: np.ndarray, red: np.ndarray, blue: np.ndarray,
                  G=2.5, C1=6.0, C2=7.5, L=1.0) -> np.ndarray:
    """Enhanced Vegetation Index"""
    denom = nir + C1 * red - C2 * blue + L
    return np.where(denom == 0, 0, G * (nir - red) / denom)


def calculate_savi(nir: np.ndarray, red: np.ndarray, L=0.5) -> np.ndarray:
    """Soil Adjusted Vegetation Index"""
    return np.where((nir + red + L) == 0, 0,
                    (1 + L) * (nir - red) / (nir + red + L))


def calculate_ndre(nir: np.ndarray, red_edge: np.ndarray) -> np.ndarray:
    """Normalized Difference Red Edge Index"""
    return np.where((nir + red_edge) == 0, 0,
                    (nir - red_edge) / (nir + red_edge))


def calculate_gndvi(nir: np.ndarray, green: np.ndarray) -> np.ndarray:
    """Green NDVI"""
    return np.where((nir + green) == 0, 0, (nir - green) / (nir + green))


def calculate_lai_from_ndvi(ndvi: np.ndarray) -> np.ndarray:
    """Empirical LAI estimation from NDVI (rice-specific)"""
    return np.clip(0.57 * np.exp(2.33 * ndvi), 0, 8.0)


def simulate_multispectral_timeseries(n_pixels=100, n_dates=18):
    """
    Simulate rice paddy multispectral imagery across growing season.
    Dates represent ~10-day intervals from transplanting (June) to harvest (Oct).
    """
    np.random.seed(42)
    doy = np.linspace(150, 290, n_dates)  # June to mid-Oct
    
    # Rice growth curve (double sigmoid)
    growth_phase = 1 / (1 + np.exp(-0.08 * (doy - 190)))
    senescence = 1 / (1 + np.exp(0.06 * (doy - 270)))
    green_fraction = growth_phase * senescence
    
    bands = {}
    for i in range(n_pixels):
        pixel_var = np.random.normal(1.0, 0.1)
        noise = np.random.normal(0, 0.02, n_dates)
        
        nir = 0.15 + 0.45 * green_fraction * pixel_var + noise
        red = 0.10 - 0.06 * green_fraction * pixel_var + np.random.normal(0, 0.01, n_dates)
        blue = 0.08 - 0.03 * green_fraction * pixel_var + np.random.normal(0, 0.01, n_dates)
        green = 0.12 - 0.04 * green_fraction * pixel_var + np.random.normal(0, 0.01, n_dates)
        red_edge = 0.12 + 0.20 * green_fraction * pixel_var + np.random.normal(0, 0.01, n_dates)
        
        nir = np.clip(nir, 0.01, 0.95)
        red = np.clip(red, 0.01, 0.95)
        blue = np.clip(blue, 0.01, 0.95)
        green = np.clip(green, 0.01, 0.95)
        red_edge = np.clip(red_edge, 0.01, 0.95)
        
        if i == 0:
            bands = {k: [] for k in ['nir', 'red', 'blue', 'green', 'red_edge']}
        bands['nir'].append(nir)
        bands['red'].append(red)
        bands['blue'].append(blue)
        bands['green'].append(green)
        bands['red_edge'].append(red_edge)
    
    for k in bands:
        bands[k] = np.array(bands[k])
    
    return bands, doy


def run_vegetation_analysis():
    """Run full vegetation index analysis and generate figures."""
    bands, doy = simulate_multispectral_timeseries(n_pixels=200, n_dates=18)
    
    ndvi = calculate_ndvi(bands['nir'], bands['red'])
    evi = calculate_evi(bands['nir'], bands['red'], bands['blue'])
    savi = calculate_savi(bands['nir'], bands['red'])
    ndre = calculate_ndre(bands['nir'], bands['red_edge'])
    gndvi = calculate_gndvi(bands['nir'], bands['green'])
    lai = calculate_lai_from_ndvi(ndvi)
    
    # Summary stats
    stats = {
        'NDVI': {'mean': ndvi.mean(axis=0), 'std': ndvi.std(axis=0)},
        'EVI': {'mean': evi.mean(axis=0), 'std': evi.std(axis=0)},
        'SAVI': {'mean': savi.mean(axis=0), 'std': savi.std(axis=0)},
        'NDRE': {'mean': ndre.mean(axis=0), 'std': ndre.std(axis=0)},
        'GNDVI': {'mean': gndvi.mean(axis=0), 'std': gndvi.std(axis=0)},
        'LAI': {'mean': lai.mean(axis=0), 'std': lai.std(axis=0)},
    }
    
    # Save numeric results
    np.savez(RESULTS_DIR / "vegetation_indices.npz",
             doy=doy, ndvi_mean=stats['NDVI']['mean'], ndvi_std=stats['NDVI']['std'],
             evi_mean=stats['EVI']['mean'], lai_mean=stats['LAI']['mean'])
    
    # --- Figure 1: VI Time Series ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    months = ['Jun', 'Jul', 'Aug', 'Sep', 'Oct']
    month_doy = [152, 182, 213, 244, 274]
    
    for ax, (name, s) in zip(axes.flat, stats.items()):
        ax.fill_between(doy, s['mean'] - s['std'], s['mean'] + s['std'], alpha=0.3)
        ax.plot(doy, s['mean'], linewidth=2)
        ax.set_title(f'{name} Time Series', fontsize=12)
        ax.set_xlabel('Day of Year')
        ax.set_ylabel(name)
        ax.set_xticks(month_doy)
        ax.set_xticklabels(months)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('Rice Paddy Vegetation Index Temporal Profiles', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig01_vegetation_indices_timeseries.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Figure 2: Spatial NDVI map (simulated) ---
    np.random.seed(123)
    grid_size = 50
    x, y = np.meshgrid(np.linspace(0, 1, grid_size), np.linspace(0, 1, grid_size))
    ndvi_spatial = 0.6 + 0.2 * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + np.random.normal(0, 0.05, (grid_size, grid_size))
    ndvi_spatial = np.clip(ndvi_spatial, 0, 1)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    titles = ['NDVI (Heading Stage)', 'EVI (Heading Stage)', 'LAI Estimate']
    evi_spatial = 0.4 + 0.15 * np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + np.random.normal(0, 0.03, (grid_size, grid_size))
    lai_spatial = calculate_lai_from_ndvi(ndvi_spatial)
    
    for ax, data, title, cmap in zip(axes,
                                      [ndvi_spatial, evi_spatial, lai_spatial],
                                      titles,
                                      ['YlGn', 'YlGn', 'viridis']):
        im = ax.imshow(data, cmap=cmap, origin='lower', extent=[0, 500, 0, 500])
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        plt.colorbar(im, ax=ax, shrink=0.8)
    
    plt.suptitle('Spatial Distribution of Vegetation Indices (Heading Stage)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig02_spatial_vegetation_map.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    peak_ndvi = stats['NDVI']['mean'].max()
    peak_doy = doy[np.argmax(stats['NDVI']['mean'])]
    peak_lai = stats['LAI']['mean'].max()
    
    print(f"=== Vegetation Index Analysis Results ===")
    print(f"Peak NDVI: {peak_ndvi:.3f} at DOY {peak_doy:.0f}")
    print(f"Peak LAI:  {peak_lai:.2f}")
    print(f"NDVI range: {stats['NDVI']['mean'].min():.3f} - {peak_ndvi:.3f}")
    print(f"EVI range:  {stats['EVI']['mean'].min():.3f} - {stats['EVI']['mean'].max():.3f}")
    print(f"Figures saved to {FIGURES_DIR}")
    
    return stats, doy, ndvi_spatial


if __name__ == "__main__":
    run_vegetation_analysis()
