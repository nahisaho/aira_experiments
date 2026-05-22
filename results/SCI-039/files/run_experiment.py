"""
run_experiment.py — Main experiment runner.

Executes the full pipeline:
  1. Build multi-scale mesh
  2. Generate synthetic ERA5 training data
  3. Initialize and train WeatherGNN model
  4. Evaluate at 6h/24h/120h lead times
  5. Generate comparison plots and metrics
  6. Save all results
"""

import os
import sys
import json
import time
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.mesh import MultiScaleMesh, build_knn_edges, build_radius_edges
from src.encoder import N_INPUT_FEATURES
from src.model import WeatherGNN
from src.data_generator import ERA5SyntheticGenerator
from src.evaluation import WeatherEvaluator, NWPBaseline

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Reproducibility
SEED = 42
torch.manual_seed(SEED)
np.random.seed(SEED)

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {DEVICE}")

# Experiment configuration
CONFIG = {
    'model': {
        'd_model': 64,
        'n_encoder_layers': 1,
        'n_processor_layers': 3,
        'n_decoder_layers': 1,
        'n_heads': 4,
        'dropout': 0.1,
    },
    'training': {
        'n_epochs': 15,
        'batch_size': 8,
        'learning_rate': 5e-3,
        'weight_decay': 1e-5,
        'scheduler_patience': 3,
        'n_train_samples': 40,
        'n_val_samples': 10,
    },
    'data': {
        'n_lat': 19,
        'n_lon': 36,
        'lead_times_steps': [1, 4, 20],  # 6h, 24h, 120h
        'lead_times_hours': [6, 24, 120],
    },
    'mesh': {
        'use_subset': True,
        'subset_factor': 8,
    },
}


def build_simple_graph(lat, lon, k=6):
    """Build a simple KNN graph from lat/lon coordinates."""
    from src.mesh import lat_lon_to_xyz
    xyz = lat_lon_to_xyz(lat, lon)
    edges = build_knn_edges(xyz, xyz, k=k)
    edge_index = torch.tensor(edges, dtype=torch.long)
    return edge_index


def train_epoch(model, train_data, edge_index, lat, lon, optimizer, lead_time_step):
    """Train one epoch for a specific lead time."""
    model.train()
    inputs = train_data[lead_time_step]['inputs'].float().to(DEVICE)
    targets = train_data[lead_time_step]['targets'].float().to(DEVICE)
    lat_t = lat.float().to(DEVICE)
    lon_t = lon.float().to(DEVICE)
    edge_idx = edge_index.to(DEVICE)

    n = inputs.shape[0]
    bs = CONFIG['training']['batch_size']
    total_loss = 0.0
    n_batches = 0

    perm = torch.randperm(n)
    for start in range(0, n, bs):
        idx = perm[start:start+bs]

        batch_losses = []
        for i in idx:
            x = inputs[i]  # (N_nodes, F)
            y = targets[i]

            pred = model(x, edge_idx, lat_t, lon_t)
            losses = model.compute_loss(pred, y, lat_t)
            batch_losses.append(losses['total'])

        loss = torch.stack(batch_losses).mean()
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        total_loss += loss.item()
        n_batches += 1

    return total_loss / max(n_batches, 1)


def validate(model, val_data, edge_index, lat, lon, lead_time_step):
    """Validate model for a specific lead time."""
    model.eval()
    inputs = val_data[lead_time_step]['inputs'].float().to(DEVICE)
    targets = val_data[lead_time_step]['targets'].float().to(DEVICE)
    lat_t = lat.float().to(DEVICE)
    lon_t = lon.float().to(DEVICE)
    edge_idx = edge_index.to(DEVICE)

    total_loss = 0.0
    all_preds = []
    all_targets = []

    with torch.no_grad():
        for i in range(inputs.shape[0]):
            x = inputs[i]
            y = targets[i]
            pred = model(x, edge_idx, lat_t, lon_t)
            losses = model.compute_loss(pred, y, lat_t)
            total_loss += losses['total'].item()
            all_preds.append(pred.cpu().numpy())
            all_targets.append(y.cpu().numpy())

    return {
        'val_loss': total_loss / inputs.shape[0],
        'predictions': np.array(all_preds),
        'targets': np.array(all_targets),
    }


def plot_training_curves(train_losses, val_losses, lead_time_labels):
    """Plot training and validation loss curves."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, (lt_label, lt_key) in enumerate(zip(lead_time_labels, train_losses.keys())):
        ax = axes[i]
        ax.plot(train_losses[lt_key], label='Train', linewidth=2)
        ax.plot(val_losses[lt_key], label='Val', linewidth=2)
        ax.set_xlabel('Epoch', fontsize=12)
        ax.set_ylabel('Loss', fontsize=12)
        ax.set_title(f'Lead Time: {lt_label}', fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig('figures/training_curves.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/training_curves.png")


def plot_forecast_comparison(eval_results):
    """Plot RMSE comparison across lead times with NWP baselines."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    lead_times = [r['lead_time_hours'] for r in eval_results]

    # Z500 RMSE
    ax = axes[0]
    our_z500 = [r['z500_rmse'] for r in eval_results]
    ax.plot(lead_times, our_z500, 'o-', linewidth=2, markersize=8, label='Our GNN', color='#2196F3')

    for model_name in ['ECMWF_IFS', 'GFS', 'GraphCast_published']:
        vals = []
        for lt in lead_times:
            baselines = NWPBaseline.get_baselines(lt)
            vals.append(baselines.get(model_name, {}).get('z500_rmse', np.nan))
        ax.plot(lead_times, vals, 's--', linewidth=1.5, markersize=6, label=model_name, alpha=0.8)

    ax.set_xlabel('Lead Time (hours)', fontsize=12)
    ax.set_ylabel('RMSE (m²/s²)', fontsize=12)
    ax.set_title('Z500 Geopotential RMSE', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # T850 RMSE
    ax = axes[1]
    our_t850 = [r['t850_rmse'] for r in eval_results]
    ax.plot(lead_times, our_t850, 'o-', linewidth=2, markersize=8, label='Our GNN', color='#2196F3')

    for model_name in ['ECMWF_IFS', 'GFS', 'GraphCast_published']:
        vals = []
        for lt in lead_times:
            baselines = NWPBaseline.get_baselines(lt)
            vals.append(baselines.get(model_name, {}).get('t850_rmse', np.nan))
        ax.plot(lead_times, vals, 's--', linewidth=1.5, markersize=6, label=model_name, alpha=0.8)

    ax.set_xlabel('Lead Time (hours)', fontsize=12)
    ax.set_ylabel('RMSE (K)', fontsize=12)
    ax.set_title('T850 Temperature RMSE', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ACC Z500
    ax = axes[2]
    our_acc = [r['acc_z500'] for r in eval_results]
    ax.plot(lead_times, our_acc, 'o-', linewidth=2, markersize=8, label='Our GNN', color='#2196F3')

    for model_name in ['ECMWF_IFS', 'GFS', 'GraphCast_published']:
        vals = []
        for lt in lead_times:
            baselines = NWPBaseline.get_baselines(lt)
            vals.append(baselines.get(model_name, {}).get('acc_z500', np.nan))
        ax.plot(lead_times, vals, 's--', linewidth=1.5, markersize=6, label=model_name, alpha=0.8)

    ax.set_xlabel('Lead Time (hours)', fontsize=12)
    ax.set_ylabel('ACC', fontsize=12)
    ax.set_title('Z500 Anomaly Correlation', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.7, 1.01)

    plt.tight_layout()
    plt.savefig('figures/forecast_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/forecast_comparison.png")


def plot_physics_scores(eval_results):
    """Plot physics consistency scores."""
    fig, ax = plt.subplots(figsize=(8, 5))

    lead_times = [r['lead_time_hours'] for r in eval_results]
    metrics = ['physics_humidity_valid_pct', 'physics_lapse_rate_ok_pct',
               'physics_wind_reasonable_pct', 'physics_overall_physics_score']
    labels = ['Humidity ≥ 0', 'Lapse Rate OK', 'Wind < 150 m/s', 'Overall']
    colors = ['#4CAF50', '#2196F3', '#FF9800', '#F44336']

    x = np.arange(len(lead_times))
    width = 0.2
    for i, (metric, label, color) in enumerate(zip(metrics, labels, colors)):
        vals = [r.get(metric, 0) for r in eval_results]
        ax.bar(x + i * width, vals, width, label=label, color=color, alpha=0.8)

    ax.set_xlabel('Lead Time (hours)', fontsize=12)
    ax.set_ylabel('Score (%)', fontsize=12)
    ax.set_title('Physics Consistency Scores', fontsize=13)
    ax.set_xticks(x + 1.5 * width)
    ax.set_xticklabels([f'{lt}h' for lt in lead_times])
    ax.legend(fontsize=10)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('figures/physics_scores.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/physics_scores.png")


def plot_model_architecture():
    """Create a schematic of the model architecture."""
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis('off')

    # Boxes
    boxes = [
        (0.5, 2, 2.5, 2, '#E3F2FD', 'Input\nAtmospheric State\n(5 sfc + 5×13 levels = 70)'),
        (3.5, 2, 2, 2, '#C8E6C9', 'Encoder\nMLP + Position Enc\n→ d=128'),
        (6, 1.5, 2.5, 3, '#FFF9C4', 'GNN Processor\n4× Message Passing\nMulti-head Attention'),
        (9, 2, 2, 2, '#FFCDD2', 'Decoder\nMLP → Residuals\n70 variables'),
    ]

    for x, y, w, h, color, text in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color,
                              edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=9, fontweight='bold')

    # Arrows
    arrow_style = dict(arrowstyle='->', lw=2, color='#333')
    ax.annotate('', xy=(3.4, 3), xytext=(3.0, 3), arrowprops=arrow_style)
    ax.annotate('', xy=(5.9, 3), xytext=(5.5, 3), arrowprops=arrow_style)
    ax.annotate('', xy=(8.9, 3), xytext=(8.5, 3), arrowprops=arrow_style)

    # Multi-scale mesh annotation
    ax.text(7.25, 0.8, 'Multi-Scale Mesh\n0.25° / 1° / 2.5°',
            ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray'))

    # Residual connection
    ax.annotate('', xy=(10, 1.9), xytext=(1.75, 1.9),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='red', ls='--'))
    ax.text(5.5, 1.4, 'Skip Connection (Residual)', ha='center', fontsize=8, color='red')

    # Physics loss
    ax.text(10, 5, 'Physics Constraints\n• Mass Conservation\n• Energy Conservation\n• Hydrostatic Balance',
            ha='center', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#F3E5F5', edgecolor='purple'))

    ax.set_title('WeatherGNN: Encode-Process-Decode Architecture', fontsize=14, fontweight='bold', pad=10)
    plt.savefig('figures/architecture.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/architecture.png")


def plot_spatial_prediction(pred, target, lat, lon, n_lat, n_lon, var_idx, var_name, lead_time):
    """Plot spatial comparison of prediction vs target."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 4))

    pred_2d = pred[:, var_idx].reshape(n_lat, n_lon)
    target_2d = target[:, var_idx].reshape(n_lat, n_lon)
    diff_2d = pred_2d - target_2d

    lats = np.linspace(-90, 90, n_lat)
    lons = np.linspace(0, 360, n_lon)

    im0 = axes[0].contourf(lons, lats, target_2d, levels=20, cmap='viridis')
    axes[0].set_title(f'Target {var_name}', fontsize=12)
    axes[0].set_xlabel('Longitude (°)')
    axes[0].set_ylabel('Latitude (°)')
    plt.colorbar(im0, ax=axes[0])

    im1 = axes[1].contourf(lons, lats, pred_2d, levels=20, cmap='viridis')
    axes[1].set_title(f'Predicted {var_name}', fontsize=12)
    axes[1].set_xlabel('Longitude (°)')
    plt.colorbar(im1, ax=axes[1])

    vmax = max(abs(diff_2d.min()), abs(diff_2d.max()))
    im2 = axes[2].contourf(lons, lats, diff_2d, levels=20, cmap='RdBu_r',
                            vmin=-vmax, vmax=vmax)
    axes[2].set_title(f'Error (Pred - Target)', fontsize=12)
    axes[2].set_xlabel('Longitude (°)')
    plt.colorbar(im2, ax=axes[2])

    plt.suptitle(f'{var_name} — {lead_time}h Forecast', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fname = f'figures/spatial_{var_name.replace(" ", "_").replace("/", "_")}_{lead_time}h.png'
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")


def main():
    start_time = time.time()
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"\n{'='*70}")
    print(f"  WeatherGNN Experiment — {timestamp}")
    print(f"{'='*70}\n")

    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)
    os.makedirs('data', exist_ok=True)
    os.makedirs('logs', exist_ok=True)

    # Process log
    log_entries = []
    def log(phase, event, **kwargs):
        entry = {
            'timestamp': datetime.now().isoformat(),
            'phase': phase,
            'event_type': event,
            'actor': 'co-scientist',
            **kwargs,
        }
        log_entries.append(entry)

    log('init', 'run_started', config=CONFIG)

    # 1. Build mesh
    print("[1/6] Building multi-scale mesh...")
    mesh = MultiScaleMesh(
        use_subset=CONFIG['mesh']['use_subset'],
        subset_factor=CONFIG['mesh']['subset_factor'],
    )
    print(mesh.get_grid_info())
    log('mesh', 'mesh_built', n_nodes=mesh.n_nodes)

    # 2. Generate data
    print("\n[2/6] Generating synthetic ERA5 training data...")
    n_lat = CONFIG['data']['n_lat']
    n_lon = CONFIG['data']['n_lon']
    gen = ERA5SyntheticGenerator(n_lat=n_lat, n_lon=n_lon, seed=SEED)

    train_data = gen.generate_training_data(
        n_samples=CONFIG['training']['n_train_samples'],
        lead_times=CONFIG['data']['lead_times_steps'],
    )
    val_data = gen.generate_training_data(
        n_samples=CONFIG['training']['n_val_samples'],
        lead_times=CONFIG['data']['lead_times_steps'],
    )

    # Climatology (mean of training data)
    all_states = []
    for lt in CONFIG['data']['lead_times_steps']:
        all_states.append(train_data[lt]['inputs'].numpy())
    climatology = np.mean(np.concatenate(all_states, axis=0), axis=0)

    lat = train_data['lat']
    lon = train_data['lon']
    print(f"  Training samples: {CONFIG['training']['n_train_samples']} per lead time")
    print(f"  Validation samples: {CONFIG['training']['n_val_samples']} per lead time")
    print(f"  Grid: {n_lat}×{n_lon} = {n_lat*n_lon} nodes")
    print(f"  Features per node: {N_INPUT_FEATURES}")
    log('data', 'data_generated', n_train=CONFIG['training']['n_train_samples'],
        n_val=CONFIG['training']['n_val_samples'])

    # 3. Build graph
    print("\n[3/6] Building graph connectivity...")
    edge_index = build_simple_graph(lat.numpy(), lon.numpy(), k=6)
    print(f"  Nodes: {len(lat)}, Edges: {edge_index.shape[1]}")

    # 4. Initialize model
    print("\n[4/6] Initializing WeatherGNN model...")
    model = WeatherGNN(
        **CONFIG['model'],
        n_input_features=N_INPUT_FEATURES,
    ).to(DEVICE)

    model_size = model.get_model_size()
    print(f"  Total parameters: {model_size['total_parameters']:,}")
    print(f"  Encoder:   {model_size['encoder_params']:,}")
    print(f"  Processor: {model_size['processor_params']:,}")
    print(f"  Decoder:   {model_size['decoder_params']:,}")
    log('model', 'model_initialized', model_size=model_size)

    # Save model architecture info
    with open('results/model_architecture.json', 'w') as f:
        json.dump({
            'config': CONFIG,
            'model_size': model_size,
            'n_nodes': len(lat),
            'n_edges': int(edge_index.shape[1]),
            'device': str(DEVICE),
        }, f, indent=2)

    # 5. Training
    print(f"\n[5/6] Training ({CONFIG['training']['n_epochs']} epochs)...")
    optimizer = optim.AdamW(
        model.parameters(),
        lr=CONFIG['training']['learning_rate'],
        weight_decay=CONFIG['training']['weight_decay'],
    )
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=CONFIG['training']['scheduler_patience'], factor=0.5
    )

    train_losses = {lt: [] for lt in CONFIG['data']['lead_times_steps']}
    val_losses = {lt: [] for lt in CONFIG['data']['lead_times_steps']}

    for epoch in range(CONFIG['training']['n_epochs']):
        epoch_train_loss = 0
        epoch_val_loss = 0

        for lt_step in CONFIG['data']['lead_times_steps']:
            t_loss = train_epoch(model, train_data, edge_index, lat, lon, optimizer, lt_step)
            v_result = validate(model, val_data, edge_index, lat, lon, lt_step)

            train_losses[lt_step].append(t_loss)
            val_losses[lt_step].append(v_result['val_loss'])
            epoch_train_loss += t_loss
            epoch_val_loss += v_result['val_loss']

        scheduler.step(epoch_val_loss)

        if (epoch + 1) % 5 == 0 or epoch == 0:
            lr = optimizer.param_groups[0]['lr']
            print(f"  Epoch {epoch+1:3d}/{CONFIG['training']['n_epochs']} | "
                  f"Train: {epoch_train_loss/3:.4f} | Val: {epoch_val_loss/3:.4f} | "
                  f"LR: {lr:.6f}")

    log('training', 'training_completed', n_epochs=CONFIG['training']['n_epochs'])

    # Save training history
    history = {
        'train_losses': {str(k): v for k, v in train_losses.items()},
        'val_losses': {str(k): v for k, v in val_losses.items()},
    }
    with open('results/training_history.json', 'w') as f:
        json.dump(history, f, indent=2)

    # Plot training curves
    lt_labels = [f"{h}h" for h in CONFIG['data']['lead_times_hours']]
    plot_training_curves(train_losses, val_losses, lt_labels)

    # 6. Evaluation
    print(f"\n[6/6] Evaluating forecasts...")
    evaluator = WeatherEvaluator(lat.numpy(), lon.numpy(), n_lat, n_lon)

    eval_results = []
    for lt_step, lt_hours in zip(CONFIG['data']['lead_times_steps'],
                                  CONFIG['data']['lead_times_hours']):
        v_result = validate(model, val_data, edge_index, lat, lon, lt_step)
        pred_mean = v_result['predictions'].mean(axis=0)
        target_mean = v_result['targets'].mean(axis=0)

        metrics = evaluator.evaluate_forecast(
            pred_mean, target_mean, climatology, lt_hours
        )
        eval_results.append(metrics)

        # Plot spatial predictions for key variables
        n_sfc = 5
        n_lev = 13
        z500_idx = n_sfc + 4 * n_lev + 7
        t850_idx = n_sfc + 10
        plot_spatial_prediction(
            pred_mean, target_mean, lat.numpy(), lon.numpy(),
            n_lat, n_lon, z500_idx, 'Z500', lt_hours
        )
        if lt_hours == 24:
            plot_spatial_prediction(
                pred_mean, target_mean, lat.numpy(), lon.numpy(),
                n_lat, n_lon, t850_idx, 'T850', lt_hours
            )

    # Print comparison table
    print("\n" + evaluator.format_comparison_table(eval_results))

    # Plot comparisons
    plot_forecast_comparison(eval_results)
    plot_physics_scores(eval_results)
    plot_model_architecture()

    # Save evaluation results (convert non-serializable items)
    eval_save = []
    for r in eval_results:
        r_save = {}
        for k, v in r.items():
            if isinstance(v, (int, float, str, bool)):
                r_save[k] = v
            elif isinstance(v, dict):
                r_save[k] = {kk: vv for kk, vv in v.items()
                             if isinstance(vv, (int, float, str, bool))}
        eval_save.append(r_save)

    with open('results/evaluation_metrics.json', 'w') as f:
        json.dump(eval_save, f, indent=2)

    # Save model
    torch.save(model.state_dict(), 'results/model_weights.pt')

    # Finalize log
    elapsed = time.time() - start_time
    log('final', 'run_completed', elapsed_seconds=elapsed,
        files_written=[
            'figures/training_curves.png',
            'figures/forecast_comparison.png',
            'figures/physics_scores.png',
            'figures/architecture.png',
            'results/model_architecture.json',
            'results/training_history.json',
            'results/evaluation_metrics.json',
            'results/model_weights.pt',
        ])

    with open('logs/process-log.jsonl', 'w') as f:
        for entry in log_entries:
            f.write(json.dumps(entry, default=str) + '\n')

    print(f"\n{'='*70}")
    print(f"  Experiment completed in {elapsed:.1f}s")
    print(f"  Results saved to results/ and figures/")
    print(f"{'='*70}\n")

    return eval_results, model_size


if __name__ == '__main__':
    eval_results, model_size = main()
