"""
Visualization utilities for volcanic deformation inversion results.

Generates publication-quality figures with colorblind-friendly palettes.
"""

import numpy as np
from typing import Dict, Optional, List
import os


def plot_displacement_comparison(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    models: Dict[str, np.ndarray],
    obs_disp: Optional[np.ndarray] = None,
    station_names: Optional[List[str]] = None,
    title: str = "Displacement Comparison",
    output_path: str = "figures/displacement_comparison.png"
):
    """
    Plot displacement vectors for multiple source models.

    Parameters
    ----------
    obs_x, obs_y : station coordinates
    models : dict of model_name -> (N, 3) displacement
    obs_disp : observed displacements (optional)
    station_names : station labels
    output_path : save path
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    components = ['East', 'North', 'Up']
    colors = {'mogi': '#0072B2', 'spheroid': '#D55E00', 'fem': '#009E73'}

    for k, (comp, ax) in enumerate(zip(components, axes)):
        if obs_disp is not None:
            ax.scatter(obs_x, obs_y, c=obs_disp[:, k] * 1000,
                      cmap='cividis', s=80, edgecolors='k',
                      linewidths=0.5, zorder=5, label='Observed')
            cb = plt.colorbar(ax.collections[0], ax=ax, shrink=0.8)
            cb.set_label(f'{comp} displacement [mm]')

        for name, disp in models.items():
            scale_val = max(np.max(np.abs(disp[:, k])), 1e-10)
            ax.scatter(obs_x, obs_y, c=disp[:, k] * 1000,
                      cmap='viridis', s=40, marker='x',
                      alpha=0.7, label=name)

        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')
        ax.set_title(f'{comp} Component')
        ax.set_aspect('equal')
        ax.legend(fontsize=8)

        if station_names is not None:
            for i, name in enumerate(station_names):
                ax.annotate(name, (obs_x[i], obs_y[i]),
                           fontsize=6, ha='left')

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_insar_map(
    x: np.ndarray,
    y: np.ndarray,
    los_obs: np.ndarray,
    los_pred: np.ndarray,
    title: str = "InSAR LOS Displacement",
    output_path: str = "figures/insar_map.png"
):
    """Plot InSAR observed vs predicted LOS displacement maps."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    vmin = min(los_obs.min(), los_pred.min()) * 1000
    vmax = max(los_obs.max(), los_pred.max()) * 1000
    resid = (los_obs - los_pred) * 1000

    n_side = int(np.sqrt(len(x)))
    if n_side**2 == len(x):
        reshape = lambda a: a.reshape(n_side, n_side)
        extent = [x.min(), x.max(), y.min(), y.max()]

        im0 = axes[0].imshow(reshape(los_obs * 1000), extent=extent,
                            origin='lower', cmap='cividis',
                            vmin=vmin, vmax=vmax)
        axes[0].set_title('Observed LOS [mm]')

        im1 = axes[1].imshow(reshape(los_pred * 1000), extent=extent,
                            origin='lower', cmap='cividis',
                            vmin=vmin, vmax=vmax)
        axes[1].set_title('Predicted LOS [mm]')

        im2 = axes[2].imshow(reshape(resid), extent=extent,
                            origin='lower', cmap='RdBu_r')
        axes[2].set_title('Residual [mm]')
    else:
        axes[0].scatter(x, y, c=los_obs*1000, cmap='cividis', s=5)
        axes[0].set_title('Observed LOS [mm]')
        axes[1].scatter(x, y, c=los_pred*1000, cmap='cividis', s=5)
        axes[1].set_title('Predicted LOS [mm]')
        axes[2].scatter(x, y, c=resid, cmap='RdBu_r', s=5)
        axes[2].set_title('Residual [mm]')

    for ax in axes:
        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')
        ax.set_aspect('equal')
        plt.colorbar(ax.collections[0] if hasattr(ax, 'images') and not ax.images
                     else (ax.images[0] if ax.images else ax.collections[0]),
                     ax=ax, shrink=0.8)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_posterior_distributions(
    trace_data: Dict[str, np.ndarray],
    true_values: Optional[Dict[str, float]] = None,
    output_path: str = "figures/posterior_distributions.png"
):
    """Plot marginal posterior distributions from MCMC trace."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    n_params = len(trace_data)
    n_cols = min(3, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_params == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for i, (name, samples) in enumerate(trace_data.items()):
        ax = axes[i]
        ax.hist(samples.flatten(), bins=50, density=True,
                color='#0072B2', alpha=0.7, edgecolor='white')
        ax.set_xlabel(name)
        ax.set_ylabel('Density')

        mean = np.mean(samples)
        hdi_lo = np.percentile(samples, 3)
        hdi_hi = np.percentile(samples, 97)
        ax.axvline(mean, color='#D55E00', linestyle='-', linewidth=2,
                   label=f'Mean={mean:.2g}')
        ax.axvspan(hdi_lo, hdi_hi, alpha=0.2, color='#D55E00',
                   label=f'94% HDI')

        if true_values and name in true_values:
            ax.axvline(true_values[name], color='#009E73', linestyle='--',
                       linewidth=2, label=f'True={true_values[name]:.2g}')

        ax.legend(fontsize=7)

    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Posterior Distributions', fontsize=14, fontweight='bold')
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_kalman_timeseries(
    times: np.ndarray,
    true_dV: np.ndarray,
    filtered_dV: np.ndarray,
    filtered_std: np.ndarray,
    smoothed_dV: Optional[np.ndarray] = None,
    smoothed_std: Optional[np.ndarray] = None,
    title: str = "Volume Change Time Series",
    output_path: str = "figures/kalman_timeseries.png"
):
    """Plot Kalman filter/smoother volume change estimates."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Volume change
    ax = axes[0]
    ax.plot(times, true_dV / 1e6, 'k-', linewidth=2, label='True', zorder=3)
    ax.plot(times[:len(filtered_dV)], filtered_dV / 1e6,
            color='#0072B2', linewidth=1.5, label='EKF filtered')
    ax.fill_between(
        times[:len(filtered_dV)],
        (filtered_dV - 2*filtered_std) / 1e6,
        (filtered_dV + 2*filtered_std) / 1e6,
        color='#0072B2', alpha=0.2, label='±2σ (filter)'
    )

    if smoothed_dV is not None:
        ax.plot(times[:len(smoothed_dV)], smoothed_dV / 1e6,
                color='#D55E00', linewidth=1.5, label='RTS smoothed')
        if smoothed_std is not None:
            ax.fill_between(
                times[:len(smoothed_dV)],
                (smoothed_dV - 2*smoothed_std) / 1e6,
                (smoothed_dV + 2*smoothed_std) / 1e6,
                color='#D55E00', alpha=0.15, label='±2σ (smoother)'
            )

    ax.set_ylabel('Volume Change [×10⁶ m³]')
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Rate
    ax2 = axes[1]
    true_rate = np.gradient(true_dV, times)
    ax2.plot(times, true_rate, 'k-', linewidth=2, label='True rate')
    filt_rate = np.gradient(filtered_dV, times[:len(filtered_dV)])
    ax2.plot(times[:len(filt_rate)], filt_rate,
             color='#0072B2', linewidth=1.5, label='Filtered rate')
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Volume Rate [m³/day]')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_viscoelastic_correction(
    times_days: np.ndarray,
    correction_factors: Dict[str, np.ndarray],
    title: str = "Viscoelastic Correction Factors",
    output_path: str = "figures/viscoelastic_correction.png"
):
    """Plot viscoelastic correction factors for different rheologies."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'maxwell': '#0072B2', 'sls': '#D55E00', 'burgers': '#009E73'}

    for model, C in correction_factors.items():
        ax.plot(times_days, C, color=colors.get(model, '#666'),
                linewidth=2, label=model.upper())

    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time since pressurization [days]')
    ax.set_ylabel('Correction Factor C(t)')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_model_comparison_residuals(
    residuals: Dict[str, Dict[str, float]],
    output_path: str = "figures/model_comparison.png"
):
    """Bar chart comparing model residuals (RMS, WRMS, chi²)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    models = list(residuals.keys())
    metrics = ['rms', 'wrms']

    fig, axes = plt.subplots(1, len(metrics), figsize=(5*len(metrics), 5))
    colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']

    for j, metric in enumerate(metrics):
        vals = [residuals[m][metric] * 1000 for m in models]  # to mm
        axes[j].bar(models, vals, color=colors[:len(models)], alpha=0.8)
        axes[j].set_ylabel(f'{metric.upper()} [mm]')
        axes[j].set_title(f'{metric.upper()} by Model')

    fig.suptitle('Source Model Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
