"""
Main experiment: Train and evaluate GraphWeatherNet.
Compares against persistence, climatology, and linear regression baselines.
Evaluates at 6h, 24h, and 120h lead times across multiple resolution scales.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import json
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from model import GraphWeatherNet, SphericalGraphBuilder
from data_generator import ERA5SyntheticGenerator
from baselines import PersistenceModel, ClimatologyModel, LinearRegressionModel
from metrics import (compute_rmse_per_variable, compute_acc_per_variable,
                     compute_physics_metrics, compute_skill_score,
                     VARIABLE_NAMES, SURFACE_NAMES, PRESSURE_LEVELS)


def setup_graph(resolution=2.5):
    """Build graph structure for the given resolution."""
    builder = SphericalGraphBuilder()
    graph = builder.build_graph(resolution)

    edge_index = torch.from_numpy(graph['edge_index']).long()
    pos = torch.from_numpy(graph['pos']).float()

    # Edge features: relative position
    src, dst = edge_index
    edge_attr = pos[dst] - pos[src]

    return edge_index, edge_attr, graph['n_nodes']


def train_model(model, train_data, edge_index, edge_attr, n_epochs=30,
                lr=1e-3, physics_weight=0.1, device='cpu'):
    """Train GraphWeatherNet."""
    model = model.to(device)
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    edge_index = edge_index.to(device)
    edge_attr = edge_attr.to(device)

    history = {'total': [], 'mse_pressure': [], 'mse_surface': [],
               'mass_conservation': [], 'energy_conservation': []}

    batch_size = 4
    n_batches = max(1, len(train_data) // batch_size)

    for epoch in range(n_epochs):
        model.train()
        epoch_losses = {k: [] for k in history}
        np.random.shuffle(train_data)

        for b in range(n_batches):
            batch = train_data[b * batch_size:(b + 1) * batch_size]
            if len(batch) == 0:
                continue

            in_p = torch.stack([d['input_pressure'] for d in batch]).to(device)
            in_s = torch.stack([d['input_surface'] for d in batch]).to(device)
            tgt_p = torch.stack([d['target_pressure'] for d in batch]).to(device)
            tgt_s = torch.stack([d['target_surface'] for d in batch]).to(device)

            optimizer.zero_grad()
            pred_p, pred_s = model(in_p, in_s, edge_index, edge_attr)
            loss, loss_dict = model.compute_loss(
                pred_p, pred_s, tgt_p, tgt_s, in_p, physics_weight)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            for k, v in loss_dict.items():
                epoch_losses[k].append(v)

        scheduler.step()

        for k in history:
            history[k].append(np.mean(epoch_losses[k]))

        if (epoch + 1) % 5 == 0 or epoch == 0:
            print(f"  Epoch {epoch+1}/{n_epochs}: "
                  f"total={history['total'][-1]:.6f}, "
                  f"mse_p={history['mse_pressure'][-1]:.6f}, "
                  f"mass={history['mass_conservation'][-1]:.6f}")

    return history


def evaluate_model(model, test_data, edge_index, edge_attr, clim_p, device='cpu'):
    """Evaluate model on test set."""
    model.eval()
    edge_index = edge_index.to(device)
    edge_attr = edge_attr.to(device)

    all_pred_p, all_pred_s = [], []
    all_tgt_p, all_tgt_s = [], []
    all_in_p = []

    with torch.no_grad():
        for d in test_data:
            in_p = d['input_pressure'].unsqueeze(0).to(device)
            in_s = d['input_surface'].unsqueeze(0).to(device)
            pred_p, pred_s = model(in_p, in_s, edge_index, edge_attr)
            all_pred_p.append(pred_p.cpu())
            all_pred_s.append(pred_s.cpu())
            all_tgt_p.append(d['target_pressure'].unsqueeze(0))
            all_tgt_s.append(d['target_surface'].unsqueeze(0))
            all_in_p.append(d['input_pressure'].unsqueeze(0))

    pred_p = torch.cat(all_pred_p)
    pred_s = torch.cat(all_pred_s)
    tgt_p = torch.cat(all_tgt_p)
    tgt_s = torch.cat(all_tgt_s)
    in_p = torch.cat(all_in_p)

    rmse = compute_rmse_per_variable(pred_p, tgt_p, pred_s, tgt_s)
    acc = compute_acc_per_variable(pred_p, tgt_p, clim_p.unsqueeze(0).expand_as(tgt_p))
    physics = compute_physics_metrics(pred_p, in_p)

    return rmse, acc, physics, pred_p, pred_s, tgt_p, tgt_s


def evaluate_baselines(test_data, clim_p):
    """Evaluate all baseline models."""
    results = {}

    # Persistence
    persist = PersistenceModel()
    in_p = torch.stack([d['input_pressure'] for d in test_data])
    in_s = torch.stack([d['input_surface'] for d in test_data])
    tgt_p = torch.stack([d['target_pressure'] for d in test_data])
    tgt_s = torch.stack([d['target_surface'] for d in test_data])

    pp, ps = persist.predict(in_p, in_s)
    rmse = compute_rmse_per_variable(pp, tgt_p, ps, tgt_s)
    acc = compute_acc_per_variable(pp, tgt_p, clim_p.unsqueeze(0).expand_as(tgt_p))
    physics = compute_physics_metrics(pp, in_p)
    results['Persistence'] = {'rmse': rmse, 'acc': acc, 'physics': physics}

    # Climatology
    clim_model = ClimatologyModel()
    clim_model.fit(test_data)
    cp, cs = clim_model.predict(in_p, in_s)
    rmse = compute_rmse_per_variable(cp, tgt_p, cs, tgt_s)
    acc = compute_acc_per_variable(cp, tgt_p, clim_p.unsqueeze(0).expand_as(tgt_p))
    physics_clim = compute_physics_metrics(cp, in_p)
    results['Climatology'] = {'rmse': rmse, 'acc': acc, 'physics': physics_clim}

    # Linear Regression
    lr_model = LinearRegressionModel()
    lr_model.fit(test_data[:len(test_data)//2])  # Use half for training
    test_subset = test_data[len(test_data)//2:]
    in_p_lr = torch.stack([d['input_pressure'] for d in test_subset])
    in_s_lr = torch.stack([d['input_surface'] for d in test_subset])
    tgt_p_lr = torch.stack([d['target_pressure'] for d in test_subset])
    tgt_s_lr = torch.stack([d['target_surface'] for d in test_subset])

    lp, ls = lr_model.predict(in_p_lr, in_s_lr)
    rmse = compute_rmse_per_variable(lp, tgt_p_lr, ls, tgt_s_lr)
    acc = compute_acc_per_variable(lp, tgt_p_lr, clim_p.unsqueeze(0).expand_as(tgt_p_lr))
    physics_lr = compute_physics_metrics(lp, in_p_lr)
    results['Linear Regression'] = {'rmse': rmse, 'acc': acc, 'physics': physics_lr}

    return results


# ============ Plotting Functions ============

def plot_training_curves(history, save_path):
    """Plot training loss curves."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].plot(history['total'], 'b-', linewidth=2)
    axes[0, 0].set_title('Total Loss', fontsize=14)
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Loss')
    axes[0, 0].set_yscale('log')
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(history['mse_pressure'], 'r-', linewidth=2, label='Pressure')
    axes[0, 1].plot(history['mse_surface'], 'g-', linewidth=2, label='Surface')
    axes[0, 1].set_title('MSE Loss Components', fontsize=14)
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('MSE')
    axes[0, 1].legend()
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(history['mass_conservation'], 'm-', linewidth=2)
    axes[1, 0].set_title('Mass Conservation Loss', fontsize=14)
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('Loss')
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(history['energy_conservation'], 'c-', linewidth=2)
    axes[1, 1].set_title('Energy Conservation Loss', fontsize=14)
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('Loss')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_rmse_comparison(all_results, lead_times, save_path):
    """Plot RMSE comparison across models and lead times."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    variables = ['Temperature', 'U-wind', 'V-wind', 'Specific Humidity']
    colors = {'GraphWeatherNet': '#2196F3', 'Persistence': '#FF5722',
              'Climatology': '#4CAF50', 'Linear Regression': '#FF9800'}

    for idx, var in enumerate(variables):
        ax = axes[idx // 2, idx % 2]
        for model_name, model_results in all_results.items():
            rmses = []
            for lt in lead_times:
                if var in model_results[lt]['rmse']:
                    r = model_results[lt]['rmse'][var]
                    rmses.append(r['mean'] if isinstance(r, dict) else r)
                else:
                    rmses.append(0)
            ax.plot(lead_times, rmses, 'o-', color=colors.get(model_name, 'gray'),
                    linewidth=2, markersize=8, label=model_name)
        ax.set_title(f'{var} RMSE', fontsize=14)
        ax.set_xlabel('Lead Time (hours)')
        ax.set_ylabel('RMSE')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(lead_times)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_acc_comparison(all_results, lead_times, save_path):
    """Plot ACC comparison."""
    fig, ax = plt.subplots(figsize=(10, 7))
    variables = ['Temperature', 'U-wind', 'V-wind', 'Specific Humidity']
    colors_model = {'GraphWeatherNet': '#2196F3', 'Persistence': '#FF5722',
                    'Climatology': '#4CAF50', 'Linear Regression': '#FF9800'}

    x = np.arange(len(variables))
    width = 0.2
    offsets = np.arange(len(all_results)) - (len(all_results) - 1) / 2

    for i, (model_name, model_results) in enumerate(all_results.items()):
        # Use 24h lead time for ACC comparison
        accs = [model_results[24]['acc'].get(v, 0) for v in variables]
        bars = ax.bar(x + offsets[i] * width, accs, width,
                      label=model_name, color=colors_model.get(model_name, 'gray'),
                      alpha=0.85)

    ax.set_ylabel('Anomaly Correlation Coefficient', fontsize=12)
    ax.set_title('ACC at 24h Lead Time', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(variables, fontsize=11)
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.axhline(y=0.6, color='red', linestyle='--', alpha=0.5, label='Useful skill threshold')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_vertical_profile(rmse_results, save_path):
    """Plot RMSE vertical profiles."""
    fig, axes = plt.subplots(1, 4, figsize=(18, 8), sharey=True)
    variables = ['Temperature', 'U-wind', 'V-wind', 'Specific Humidity']
    units = ['K', 'm/s', 'm/s', 'kg/kg']

    for idx, (var, unit) in enumerate(zip(variables, units)):
        ax = axes[idx]
        if var in rmse_results and isinstance(rmse_results[var], dict) and 'levels' in rmse_results[var]:
            levels = list(rmse_results[var]['levels'].keys())
            rmses = list(rmse_results[var]['levels'].values())
            ax.plot(rmses, levels, 'b-o', linewidth=2, markersize=5)
        ax.set_title(f'{var} ({unit})', fontsize=13)
        ax.set_xlabel('RMSE', fontsize=11)
        if idx == 0:
            ax.set_ylabel('Pressure Level (hPa)', fontsize=11)
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_yticks(PRESSURE_LEVELS)
        ax.set_yticklabels([str(p) for p in PRESSURE_LEVELS])
        ax.grid(True, alpha=0.3)

    plt.suptitle('GraphWeatherNet RMSE Vertical Profiles (24h Forecast)', fontsize=15)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_physics_constraints(all_results, lead_times, save_path):
    """Plot physics constraint metrics."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    model_names = list(all_results.keys())
    colors = {'GraphWeatherNet': '#2196F3', 'Persistence': '#FF5722',
              'Climatology': '#4CAF50', 'Linear Regression': '#FF9800'}

    # Mass conservation
    for model_name in model_names:
        mass_errors = [all_results[model_name][lt]['physics']['mass_error_kg_m2']
                       for lt in lead_times]
        axes[0].plot(lead_times, mass_errors, 'o-',
                     color=colors.get(model_name, 'gray'),
                     linewidth=2, markersize=8, label=model_name)
    axes[0].set_title('Column Moisture Conservation Error', fontsize=14)
    axes[0].set_xlabel('Lead Time (hours)')
    axes[0].set_ylabel('Error (kg/m²)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xticks(lead_times)

    # Energy conservation
    for model_name in model_names:
        energy_errors = [all_results[model_name][lt]['physics']['energy_relative_error']
                         for lt in lead_times]
        axes[1].plot(lead_times, energy_errors, 'o-',
                     color=colors.get(model_name, 'gray'),
                     linewidth=2, markersize=8, label=model_name)
    axes[1].set_title('Column Energy Relative Error', fontsize=14)
    axes[1].set_xlabel('Lead Time (hours)')
    axes[1].set_ylabel('Relative Error')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xticks(lead_times)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_spatial_prediction(pred_p, tgt_p, n_lat, n_lon, save_path):
    """Plot spatial maps of prediction vs truth."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # Temperature at 500 hPa (level index 7)
    level_idx = 7
    sample_idx = 0

    pred_field = pred_p[sample_idx, :, 0, level_idx].reshape(n_lat, n_lon).numpy()
    tgt_field = tgt_p[sample_idx, :, 0, level_idx].reshape(n_lat, n_lon).numpy()
    diff_field = pred_field - tgt_field

    im1 = axes[0, 0].imshow(tgt_field, aspect='auto', cmap='RdBu_r',
                             extent=[0, 360, -90, 90])
    axes[0, 0].set_title('Target: T500 (K)', fontsize=13)
    plt.colorbar(im1, ax=axes[0, 0])

    im2 = axes[0, 1].imshow(pred_field, aspect='auto', cmap='RdBu_r',
                             extent=[0, 360, -90, 90],
                             vmin=tgt_field.min(), vmax=tgt_field.max())
    axes[0, 1].set_title('Predicted: T500 (K)', fontsize=13)
    plt.colorbar(im2, ax=axes[0, 1])

    vmax_diff = np.abs(diff_field).max()
    im3 = axes[0, 2].imshow(diff_field, aspect='auto', cmap='bwr',
                             extent=[0, 360, -90, 90],
                             vmin=-vmax_diff, vmax=vmax_diff)
    axes[0, 2].set_title('Error: T500 (K)', fontsize=13)
    plt.colorbar(im3, ax=axes[0, 2])

    # U-wind at 250 hPa (jet stream level, index 4)
    level_idx = 4
    pred_u = pred_p[sample_idx, :, 1, level_idx].reshape(n_lat, n_lon).numpy()
    tgt_u = tgt_p[sample_idx, :, 1, level_idx].reshape(n_lat, n_lon).numpy()
    diff_u = pred_u - tgt_u

    im4 = axes[1, 0].imshow(tgt_u, aspect='auto', cmap='coolwarm',
                             extent=[0, 360, -90, 90])
    axes[1, 0].set_title('Target: U250 (m/s)', fontsize=13)
    plt.colorbar(im4, ax=axes[1, 0])

    im5 = axes[1, 1].imshow(pred_u, aspect='auto', cmap='coolwarm',
                             extent=[0, 360, -90, 90],
                             vmin=tgt_u.min(), vmax=tgt_u.max())
    axes[1, 1].set_title('Predicted: U250 (m/s)', fontsize=13)
    plt.colorbar(im5, ax=axes[1, 1])

    vmax_u = np.abs(diff_u).max()
    im6 = axes[1, 2].imshow(diff_u, aspect='auto', cmap='bwr',
                             extent=[0, 360, -90, 90],
                             vmin=-vmax_u, vmax=vmax_u)
    axes[1, 2].set_title('Error: U250 (m/s)', fontsize=13)
    plt.colorbar(im6, ax=axes[1, 2])

    for ax in axes.flat:
        ax.set_xlabel('Longitude (°)')
        ax.set_ylabel('Latitude (°)')

    plt.suptitle('GraphWeatherNet Spatial Prediction (24h Forecast)', fontsize=15)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_skill_scores(all_results, lead_times, save_path):
    """Plot skill scores relative to persistence."""
    fig, ax = plt.subplots(figsize=(10, 7))
    variables = ['Temperature', 'U-wind', 'V-wind', 'Specific Humidity']
    colors_var = ['#E53935', '#1E88E5', '#43A047', '#FB8C00']

    for i, var in enumerate(variables):
        skills = []
        for lt in lead_times:
            gw_rmse = all_results['GraphWeatherNet'][lt]['rmse'][var]
            gw_rmse = gw_rmse['mean'] if isinstance(gw_rmse, dict) else gw_rmse
            pers_rmse = all_results['Persistence'][lt]['rmse'][var]
            pers_rmse = pers_rmse['mean'] if isinstance(pers_rmse, dict) else pers_rmse
            skills.append(compute_skill_score(gw_rmse, pers_rmse))
        ax.plot(lead_times, skills, 'o-', color=colors_var[i],
                linewidth=2, markersize=8, label=var)

    ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax.set_title('GraphWeatherNet Skill Score vs Persistence', fontsize=14)
    ax.set_xlabel('Lead Time (hours)', fontsize=12)
    ax.set_ylabel('Skill Score', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(lead_times)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_multi_resolution(resolution_results, save_path):
    """Plot multi-resolution comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    resolutions = list(resolution_results.keys())
    res_labels = [f"{r}°" for r in resolutions]

    variables = ['Temperature', 'U-wind']
    colors_var = ['#E53935', '#1E88E5']

    for idx, var in enumerate(variables):
        rmses = []
        nodes = []
        for r in resolutions:
            rmse = resolution_results[r]['rmse'].get(var, {})
            rmses.append(rmse.get('mean', 0) if isinstance(rmse, dict) else rmse)
            nodes.append(resolution_results[r]['n_nodes'])

        axes[0].bar(np.arange(len(resolutions)) + idx * 0.3, rmses, 0.3,
                    label=var, color=colors_var[idx], alpha=0.85)

    axes[0].set_title('RMSE by Resolution (24h)', fontsize=14)
    axes[0].set_xticks(np.arange(len(resolutions)) + 0.15)
    axes[0].set_xticklabels(res_labels)
    axes[0].set_ylabel('RMSE')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, axis='y')

    # Nodes per resolution
    nodes_list = [resolution_results[r]['n_nodes'] for r in resolutions]
    axes[1].bar(res_labels, nodes_list, color='#7B1FA2', alpha=0.85)
    axes[1].set_title('Grid Points by Resolution', fontsize=14)
    axes[1].set_ylabel('Number of Grid Points')
    axes[1].set_yscale('log')
    axes[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_architecture_diagram(save_path):
    """Create architecture diagram."""
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Boxes
    boxes = [
        (1, 3, 2.5, 2, 'Input\nAtmospheric\nState\n(T, u, v, q)', '#BBDEFB'),
        (4.5, 3, 2.5, 2, 'Pressure Level\nEncoder\n(MLP + Embedding)', '#C8E6C9'),
        (8, 3, 2.5, 2, 'Multi-Scale\nGNN Processor\n(4 MP Blocks)', '#FFF9C4'),
        (11.5, 3, 2.5, 2, 'Pressure Level\nDecoder\n(MLP)', '#FFCCBC'),
        (8, 0.5, 2.5, 1.5, 'Physics\nConstraints\n(Mass + Energy)', '#E1BEE7'),
    ]

    for x, y, w, h, text, color in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color,
                              edgecolor='black', linewidth=2, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=10, fontweight='bold', zorder=3)

    # Arrows
    arrows = [(3.5, 4), (7, 4), (10.5, 4)]
    for x, y in arrows:
        ax.annotate('', xy=(x + 0.5, y), xytext=(x, y),
                    arrowprops=dict(arrowstyle='->', lw=2.5, color='#333'))

    # Physics constraint arrow
    ax.annotate('', xy=(9.25, 3), xytext=(9.25, 2),
                arrowprops=dict(arrowstyle='->', lw=2, color='#9C27B0', linestyle='dashed'))

    # Output label
    ax.annotate('Predicted\nState (t+Δt)', xy=(14.5, 4), fontsize=12,
                fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#E3F2FD', edgecolor='#1565C0'))
    ax.annotate('', xy=(14, 4), xytext=(14, 4),
                arrowprops=dict(arrowstyle='->', lw=2.5, color='#333'))

    # Title
    ax.text(8, 7.2, 'GraphWeatherNet Architecture', fontsize=18,
            fontweight='bold', ha='center', va='center')
    ax.text(8, 6.5, 'Encoder-Processor-Decoder with Physics-Informed Constraints',
            fontsize=12, ha='center', va='center', color='#555')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


# ============ Main Experiment ============

def main():
    print("=" * 70)
    print("GraphWeatherNet: Data-Driven Weather Prediction Experiment")
    print("=" * 70)

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")

    figures_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    results_all = {}
    lead_times = [6, 24, 120]

    # ---- Architecture Diagram ----
    print("\n[1/7] Generating architecture diagram...")
    plot_architecture_diagram(os.path.join(figures_dir, 'architecture.png'))

    # ---- Data Generation ----
    print("\n[2/7] Generating synthetic ERA5-like data...")
    resolution = 10.0  # Use coarser grid for tractable CPU training
    gen = ERA5SyntheticGenerator(resolution=resolution, seed=42)
    datasets = gen.generate_dataset(n_samples=40, lead_times=lead_times)

    for lt in lead_times:
        print(f"  Lead time {lt}h: {len(datasets[lt])} samples, "
              f"nodes={gen.n_nodes}, shape=({gen.n_nodes}, 4, 13)")

    # ---- Graph Setup ----
    print("\n[3/7] Building spherical graph...")
    edge_index, edge_attr, n_nodes = setup_graph(resolution)
    print(f"  Nodes: {n_nodes}, Edges: {edge_index.shape[1]}")

    # ---- Climatology ----
    clim_p_all = torch.stack([d['input_pressure'] for d in datasets[6]])
    clim_p = clim_p_all.mean(dim=0)

    # ---- Train & Evaluate GraphWeatherNet at each lead time ----
    print("\n[4/7] Training and evaluating GraphWeatherNet...")
    results_all['GraphWeatherNet'] = {}
    training_history = {}

    for lt in lead_times:
        print(f"\n  --- Lead Time: {lt}h ---")
        data = datasets[lt]
        train_data = data[:30]
        test_data = data[30:]

        model = GraphWeatherNet(embed_dim=64, n_pressure_levels=13,
                                n_variables=4, n_blocks=3)
        n_params = sum(p.numel() for p in model.parameters())
        print(f"  Model parameters: {n_params:,}")

        t0 = time.time()
        history = train_model(model, train_data, edge_index, edge_attr,
                              n_epochs=20, lr=5e-4, physics_weight=0.01,
                              device=device)
        train_time = time.time() - t0
        print(f"  Training time: {train_time:.1f}s")

        training_history[lt] = history

        rmse, acc, physics, pred_p, pred_s, tgt_p, tgt_s = evaluate_model(
            model, test_data, edge_index, edge_attr, clim_p, device)

        results_all['GraphWeatherNet'][lt] = {
            'rmse': rmse, 'acc': acc, 'physics': physics,
            'train_time': train_time, 'n_params': n_params
        }

        # Save spatial prediction for 24h
        if lt == 24:
            pred_p_24, tgt_p_24 = pred_p, tgt_p
            rmse_24 = rmse

    # ---- Evaluate Baselines ----
    print("\n[5/7] Evaluating baselines...")
    for baseline_name in ['Persistence', 'Climatology', 'Linear Regression']:
        results_all[baseline_name] = {}

    for lt in lead_times:
        print(f"  --- Lead Time: {lt}h ---")
        data = datasets[lt]
        test_data = data[30:]
        baseline_results = evaluate_baselines(test_data, clim_p)
        for bname, bresult in baseline_results.items():
            results_all[bname][lt] = bresult

    # ---- Multi-Resolution Experiment ----
    print("\n[6/7] Multi-resolution experiment...")
    resolution_results = {}
    for res in [10.0, 15.0, 30.0]:
        print(f"  Resolution: {res}°")
        gen_r = ERA5SyntheticGenerator(resolution=res, seed=42)
        data_r = gen_r.generate_dataset(n_samples=20, lead_times=[24])
        edge_r, attr_r, nodes_r = setup_graph(res)

        model_r = GraphWeatherNet(embed_dim=32, n_pressure_levels=13,
                                  n_variables=4, n_blocks=2)
        train_model(model_r, data_r[24][:15], edge_r, attr_r,
                    n_epochs=10, lr=1e-3, device=device)

        clim_p_r = torch.stack([d['input_pressure'] for d in data_r[24]]).mean(dim=0)
        rmse_r, acc_r, phys_r, _, _, _, _ = evaluate_model(
            model_r, data_r[24][15:], edge_r, attr_r, clim_p_r, device)
        resolution_results[res] = {
            'rmse': rmse_r, 'acc': acc_r, 'physics': phys_r, 'n_nodes': nodes_r
        }

    # ---- Generate Figures ----
    print("\n[7/7] Generating figures...")
    plot_training_curves(training_history[24],
                         os.path.join(figures_dir, 'training_curves.png'))
    plot_rmse_comparison(results_all, lead_times,
                         os.path.join(figures_dir, 'rmse_comparison.png'))
    plot_acc_comparison(results_all, lead_times,
                        os.path.join(figures_dir, 'acc_comparison.png'))
    plot_vertical_profile(rmse_24,
                          os.path.join(figures_dir, 'vertical_profile.png'))
    plot_physics_constraints(results_all, lead_times,
                             os.path.join(figures_dir, 'physics_constraints.png'))
    plot_spatial_prediction(pred_p_24, tgt_p_24, gen.n_lat, gen.n_lon,
                            os.path.join(figures_dir, 'spatial_prediction.png'))
    plot_skill_scores(results_all, lead_times,
                      os.path.join(figures_dir, 'skill_scores.png'))
    plot_multi_resolution(resolution_results,
                          os.path.join(figures_dir, 'multi_resolution.png'))

    # ---- Save results ----
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)

    summary = {}
    for model_name in results_all:
        summary[model_name] = {}
        for lt in lead_times:
            if lt not in results_all[model_name]:
                continue
            r = results_all[model_name][lt]
            lt_summary = {'acc': r['acc'], 'physics': r['physics']}
            rmse_summary = {}
            for var_name in VARIABLE_NAMES + SURFACE_NAMES:
                if var_name in r['rmse']:
                    val = r['rmse'][var_name]
                    rmse_summary[var_name] = val['mean'] if isinstance(val, dict) else val
            lt_summary['rmse'] = rmse_summary
            summary[model_name][str(lt)] = lt_summary

    # Print key results
    for lt in lead_times:
        print(f"\n  Lead Time: {lt}h")
        print(f"  {'Model':<20} {'T RMSE':>10} {'U RMSE':>10} {'T ACC':>10}")
        print(f"  {'-'*52}")
        for model_name in results_all:
            if lt in results_all[model_name]:
                r = results_all[model_name][lt]
                t_rmse = r['rmse']['Temperature']
                t_rmse = t_rmse['mean'] if isinstance(t_rmse, dict) else t_rmse
                u_rmse = r['rmse']['U-wind']
                u_rmse = u_rmse['mean'] if isinstance(u_rmse, dict) else u_rmse
                t_acc = r['acc'].get('Temperature', 0)
                print(f"  {model_name:<20} {t_rmse:>10.4f} {u_rmse:>10.4f} {t_acc:>10.4f}")

    # Save JSON
    results_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                'results.json')
    with open(results_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")

    print("\nExperiment complete!")
    return summary, resolution_results


if __name__ == '__main__':
    main()
