"""
Module 3: 土壌センサーデータの空間補間
Soil sensor data spatial interpolation using Kriging
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RBFInterpolator
from scipy.spatial.distance import pdist, squareform
from pathlib import Path

FIGURES_DIR = Path(__file__).parent.parent / "figures"
RESULTS_DIR = Path(__file__).parent.parent / "results"
DATA_DIR = Path(__file__).parent.parent / "data"


def generate_soil_sensor_data(n_sensors=30, field_size=500):
    """
    Generate synthetic soil sensor network data for a rice paddy field.
    Parameters measured: volumetric water content (VWC), EC, pH
    """
    np.random.seed(42)
    
    # Sensor locations (stratified random within field)
    grid_n = int(np.sqrt(n_sensors))
    cell_size = field_size / grid_n
    x_sensors, y_sensors = [], []
    for i in range(grid_n):
        for j in range(grid_n):
            x_sensors.append(cell_size * (i + np.random.uniform(0.2, 0.8)))
            y_sensors.append(cell_size * (j + np.random.uniform(0.2, 0.8)))
    
    # Add a few extra random sensors
    extra = n_sensors - len(x_sensors)
    if extra > 0:
        x_sensors.extend(np.random.uniform(20, field_size - 20, extra))
        y_sensors.extend(np.random.uniform(20, field_size - 20, extra))
    
    x_sensors = np.array(x_sensors[:n_sensors])
    y_sensors = np.array(y_sensors[:n_sensors])
    
    # Spatially correlated soil properties using Gaussian fields
    def generate_spatial_field(x, y, mean, sill, range_param, nugget=0):
        dist_matrix = squareform(pdist(np.column_stack([x, y])))
        # Exponential variogram model
        cov_matrix = (sill - nugget) * np.exp(-3 * dist_matrix / range_param)
        np.fill_diagonal(cov_matrix, sill)
        L = np.linalg.cholesky(cov_matrix + 1e-6 * np.eye(len(x)))
        z = np.random.randn(len(x))
        return mean + L @ z
    
    vwc = generate_spatial_field(x_sensors, y_sensors, mean=0.35, sill=0.003, range_param=200, nugget=0.0005)
    ec = generate_spatial_field(x_sensors, y_sensors, mean=0.8, sill=0.04, range_param=150, nugget=0.005)
    ph = generate_spatial_field(x_sensors, y_sensors, mean=5.8, sill=0.15, range_param=250, nugget=0.02)
    
    vwc = np.clip(vwc, 0.15, 0.55)
    ec = np.clip(ec, 0.2, 2.0)
    ph = np.clip(ph, 4.5, 7.5)
    
    sensor_df = pd.DataFrame({
        'sensor_id': [f'S{i+1:02d}' for i in range(n_sensors)],
        'x_m': np.round(x_sensors, 1),
        'y_m': np.round(y_sensors, 1),
        'vwc': np.round(vwc, 3),
        'ec_dSm': np.round(ec, 2),
        'ph': np.round(ph, 2),
    })
    
    sensor_df.to_csv(DATA_DIR / "soil_sensor_data.csv", index=False)
    return sensor_df


def ordinary_kriging(x_obs, y_obs, z_obs, x_pred, y_pred, variogram_model='exponential',
                     sill=None, range_param=None, nugget=0):
    """
    Simple Ordinary Kriging implementation.
    """
    coords_obs = np.column_stack([x_obs, y_obs])
    coords_pred = np.column_stack([x_pred, y_pred])
    
    # Estimate variogram parameters if not provided
    dists = pdist(coords_obs)
    if sill is None:
        sill = np.var(z_obs)
    if range_param is None:
        range_param = np.percentile(dists, 60)
    
    def variogram(h):
        if variogram_model == 'exponential':
            return nugget + (sill - nugget) * (1 - np.exp(-3 * h / range_param))
        elif variogram_model == 'spherical':
            result = np.where(h == 0, 0,
                              np.where(h <= range_param,
                                       nugget + (sill - nugget) * (1.5 * h / range_param - 0.5 * (h / range_param) ** 3),
                                       sill))
            return result
        return nugget + (sill - nugget) * (1 - np.exp(-3 * h / range_param))
    
    n_obs = len(z_obs)
    # Build kriging matrix
    dist_obs = squareform(pdist(coords_obs))
    K = variogram(dist_obs)
    K_ext = np.zeros((n_obs + 1, n_obs + 1))
    K_ext[:n_obs, :n_obs] = K
    K_ext[n_obs, :n_obs] = 1
    K_ext[:n_obs, n_obs] = 1
    K_ext[n_obs, n_obs] = 0
    
    K_inv = np.linalg.solve(K_ext + 1e-8 * np.eye(n_obs + 1), np.eye(n_obs + 1))
    
    z_pred = np.zeros(len(coords_pred))
    z_var = np.zeros(len(coords_pred))
    
    for i, cp in enumerate(coords_pred):
        dist_to_obs = np.sqrt(np.sum((coords_obs - cp) ** 2, axis=1))
        k_vec = np.zeros(n_obs + 1)
        k_vec[:n_obs] = variogram(dist_to_obs)
        k_vec[n_obs] = 1
        
        weights = K_inv @ k_vec
        z_pred[i] = np.sum(weights[:n_obs] * z_obs)
        z_var[i] = np.sum(weights[:n_obs] * k_vec[:n_obs]) + weights[n_obs]
    
    return z_pred, z_var


def run_soil_interpolation():
    """Run soil data interpolation and generate figures."""
    sensor_df = generate_soil_sensor_data(n_sensors=30)
    
    # Create prediction grid
    grid_res = 5  # 5m resolution
    x_grid = np.arange(0, 500, grid_res)
    y_grid = np.arange(0, 500, grid_res)
    xx, yy = np.meshgrid(x_grid, y_grid)
    x_pred = xx.ravel()
    y_pred = yy.ravel()
    
    # Interpolate each soil property
    results = {}
    for prop, col in [('VWC', 'vwc'), ('EC', 'ec_dSm'), ('pH', 'ph')]:
        z_pred, z_var = ordinary_kriging(
            sensor_df['x_m'].values, sensor_df['y_m'].values,
            sensor_df[col].values, x_pred, y_pred
        )
        results[prop] = {
            'pred': z_pred.reshape(xx.shape),
            'var': z_var.reshape(xx.shape),
        }
    
    # Save interpolated data
    np.savez(RESULTS_DIR / "soil_interpolation.npz",
             x_grid=x_grid, y_grid=y_grid,
             vwc=results['VWC']['pred'], ec=results['EC']['pred'], ph=results['pH']['pred'])
    
    # --- Figure 5: Soil Property Maps ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    
    props = [('VWC', 'Volumetric Water Content (m³/m³)', 'Blues'),
             ('EC', 'Electrical Conductivity (dS/m)', 'YlOrRd'),
             ('pH', 'Soil pH', 'RdYlGn')]
    col_names = ['vwc', 'ec_dSm', 'ph']
    
    for j, (prop, title, cmap) in enumerate(props):
        # Interpolated map
        ax = axes[0, j]
        im = ax.imshow(results[prop]['pred'], cmap=cmap, origin='lower',
                       extent=[0, 500, 0, 500], aspect='equal')
        ax.scatter(sensor_df['x_m'], sensor_df['y_m'], c='black', s=20, marker='^',
                   label='Sensors', zorder=5)
        ax.set_title(f'{title}\n(Kriging Interpolation)', fontsize=11)
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.legend(fontsize=8)
        
        # Prediction variance
        ax = axes[1, j]
        im = ax.imshow(np.sqrt(np.abs(results[prop]['var'])), cmap='Reds', origin='lower',
                       extent=[0, 500, 0, 500], aspect='equal')
        ax.scatter(sensor_df['x_m'], sensor_df['y_m'], c='black', s=20, marker='^', zorder=5)
        ax.set_title(f'{title}\n(Kriging Std. Dev.)', fontsize=11)
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        plt.colorbar(im, ax=ax, shrink=0.8)
    
    plt.suptitle('Soil Property Spatial Interpolation — Rice Paddy Field', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig05_soil_interpolation.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Figure 6: Experimental Variogram ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    
    for ax, (prop, title, _), col in zip(axes, props, col_names):
        coords = sensor_df[['x_m', 'y_m']].values
        z = sensor_df[col].values
        dists = squareform(pdist(coords))
        
        # Compute experimental variogram
        max_dist = 350
        n_lags = 15
        lag_edges = np.linspace(0, max_dist, n_lags + 1)
        lag_centers = (lag_edges[:-1] + lag_edges[1:]) / 2
        gamma = np.zeros(n_lags)
        counts = np.zeros(n_lags)
        
        n = len(z)
        for i in range(n):
            for jj in range(i + 1, n):
                d = dists[i, jj]
                for k in range(n_lags):
                    if lag_edges[k] <= d < lag_edges[k + 1]:
                        gamma[k] += (z[i] - z[jj]) ** 2
                        counts[k] += 1
                        break
        
        valid = counts > 0
        gamma[valid] /= (2 * counts[valid])
        
        ax.scatter(lag_centers[valid], gamma[valid], c='steelblue', s=40, zorder=5)
        
        # Fit exponential model
        sill = np.var(z)
        range_est = 200
        h_fit = np.linspace(0, max_dist, 100)
        gamma_fit = sill * (1 - np.exp(-3 * h_fit / range_est))
        ax.plot(h_fit, gamma_fit, 'r-', linewidth=2, label=f'Exp. model (sill={sill:.4f})')
        
        ax.set_xlabel('Lag Distance (m)')
        ax.set_ylabel('Semivariance')
        ax.set_title(f'{prop} Experimental Variogram', fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('Experimental Variograms for Soil Properties', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig06_variograms.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"=== Soil Interpolation Results ===")
    for prop, col in [('VWC', 'vwc'), ('EC', 'ec_dSm'), ('pH', 'ph')]:
        print(f"{prop}: mean={sensor_df[col].mean():.3f}, std={sensor_df[col].std():.3f}")
    print(f"Grid resolution: {grid_res}m × {grid_res}m")
    print(f"Prediction grid: {xx.shape}")
    
    return sensor_df, results


if __name__ == "__main__":
    run_soil_interpolation()
