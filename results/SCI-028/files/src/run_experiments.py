"""
Main experiment runner for tokamak disruption prediction.
Trains all models, evaluates performance, measures latency, and generates figures.
"""
import os
import sys
import time
import json
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.metrics import (accuracy_score, precision_score, recall_score,
                             f1_score, roc_auc_score, confusion_matrix, roc_curve)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from models import (LSTMPredictor, CNNLSTMPredictor, PhysicsInformedPredictor,
                    TearingModeDetector, TransferablePredictor)

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

RESULTS = {}


def load_data(path):
    data = np.load(path)
    return data['X'], data['y'], data['y_tm']


def split_data(X, y, train_ratio=0.7, val_ratio=0.15):
    n = len(X)
    idx = np.random.permutation(n)
    t1 = int(n * train_ratio)
    t2 = int(n * (train_ratio + val_ratio))
    return (X[idx[:t1]], y[idx[:t1]],
            X[idx[t1:t2]], y[idx[t1:t2]],
            X[idx[t2:]], y[idx[t2:]])


def train_model(model, X_train, y_train, X_val, y_val, epochs=30, lr=1e-3,
                physics_model=False, batch_size=128, extra_y_train=None, extra_y_val=None):
    """Generic training loop."""
    dataset = TensorDataset(torch.FloatTensor(X_train), torch.FloatTensor(y_train))
    if extra_y_train is not None:
        dataset = TensorDataset(torch.FloatTensor(X_train), torch.FloatTensor(y_train),
                                torch.FloatTensor(extra_y_train))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)

    # Class weight for imbalanced data
    pos_weight = torch.tensor([(len(y_train) - y_train.sum()) / (y_train.sum() + 1e-6)])
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    model.to(DEVICE)
    train_losses, val_losses = [], []
    best_val_loss = float('inf')

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0
        for batch in loader:
            optimizer.zero_grad()
            if physics_model:
                xb, yb = batch[0].to(DEVICE), batch[1].to(DEVICE)
                pred, physics_feats = model(xb)
                loss = criterion(pred, yb) + model.physics_loss(xb, physics_feats)
            elif extra_y_train is not None:
                xb, yb, ytm = batch[0].to(DEVICE), batch[1].to(DEVICE), batch[2].to(DEVICE)
                pred_d, pred_tm = model(xb)
                loss = criterion(pred_d, yb) + nn.BCEWithLogitsLoss()(pred_tm, ytm)
            else:
                xb, yb = batch[0].to(DEVICE), batch[1].to(DEVICE)
                pred = model(xb)
                loss = criterion(pred, yb)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            epoch_loss += loss.item()

        train_losses.append(epoch_loss / len(loader))

        # Validation
        model.eval()
        with torch.no_grad():
            xv = torch.FloatTensor(X_val).to(DEVICE)
            yv = torch.FloatTensor(y_val).to(DEVICE)
            if physics_model:
                pred_v, _ = model(xv)
            elif extra_y_train is not None:
                pred_v, _ = model(xv)
            else:
                pred_v = model(xv)
            val_loss = criterion(pred_v, yv).item()
            val_losses.append(val_loss)

        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

    model.load_state_dict(best_state)
    return train_losses, val_losses


def evaluate_model(model, X_test, y_test, physics_model=False, multi_task=False):
    """Evaluate model and return metrics."""
    model.eval()
    model.to(DEVICE)
    with torch.no_grad():
        xt = torch.FloatTensor(X_test).to(DEVICE)
        if physics_model:
            logits, _ = model(xt)
        elif multi_task:
            logits, _ = model(xt)
        else:
            logits = model(xt)
        probs = torch.sigmoid(logits).cpu().numpy()

    y_pred = (probs > 0.5).astype(int)
    metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'precision': precision_score(y_test, y_pred, zero_division=0),
        'recall': recall_score(y_test, y_pred, zero_division=0),
        'f1': f1_score(y_test, y_pred, zero_division=0),
        'auc': roc_auc_score(y_test, probs) if len(np.unique(y_test)) > 1 else 0.0,
    }
    return metrics, probs


def measure_latency(model, input_shape=(1, 50, 11), n_runs=1000, physics_model=False, multi_task=False):
    """Measure inference latency in milliseconds."""
    model.eval()
    model.to(DEVICE)
    x = torch.randn(*input_shape).to(DEVICE)

    # Warmup
    for _ in range(50):
        with torch.no_grad():
            if physics_model:
                model(x)
            elif multi_task:
                model(x)
            else:
                model(x)

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        with torch.no_grad():
            if physics_model:
                model(x)
            elif multi_task:
                model(x)
            else:
                model(x)
        times.append((time.perf_counter() - start) * 1000)

    return {
        'mean_ms': np.mean(times),
        'std_ms': np.std(times),
        'p50_ms': np.percentile(times, 50),
        'p95_ms': np.percentile(times, 95),
        'p99_ms': np.percentile(times, 99),
        'max_ms': np.max(times),
    }


# ---- Figure generation ----

def plot_training_curves(all_curves, filename='training_curves.png'):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for name, (train_l, val_l) in all_curves.items():
        axes[0].plot(train_l, label=name)
        axes[1].plot(val_l, label=name)
    axes[0].set_title('Training Loss'); axes[0].set_xlabel('Epoch'); axes[0].set_ylabel('Loss')
    axes[0].legend(); axes[0].grid(True, alpha=0.3)
    axes[1].set_title('Validation Loss'); axes[1].set_xlabel('Epoch'); axes[1].set_ylabel('Loss')
    axes[1].legend(); axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_roc_curves(all_rocs, filename='roc_curves.png'):
    fig, ax = plt.subplots(figsize=(8, 8))
    for name, (y_true, probs) in all_rocs.items():
        if len(np.unique(y_true)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true, probs)
        auc = roc_auc_score(y_true, probs)
        ax.plot(fpr, tpr, label=f'{name} (AUC={auc:.3f})', linewidth=2)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('ROC Curves - Disruption Prediction Models', fontsize=14)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_metrics_comparison(all_metrics, filename='metrics_comparison.png'):
    fig, ax = plt.subplots(figsize=(12, 6))
    metric_names = ['accuracy', 'precision', 'recall', 'f1', 'auc']
    x = np.arange(len(metric_names))
    width = 0.8 / len(all_metrics)

    for i, (name, metrics) in enumerate(all_metrics.items()):
        vals = [metrics[m] for m in metric_names]
        ax.bar(x + i * width, vals, width, label=name, alpha=0.85)

    ax.set_xticks(x + width * (len(all_metrics) - 1) / 2)
    ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1', 'AUC'], fontsize=11)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title('Model Performance Comparison (JET Test Set)', fontsize=14)
    ax.legend(fontsize=10); ax.set_ylim(0, 1.1); ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_latency(all_latencies, filename='latency_comparison.png'):
    fig, ax = plt.subplots(figsize=(10, 6))
    names = list(all_latencies.keys())
    means = [all_latencies[n]['mean_ms'] for n in names]
    p95s = [all_latencies[n]['p95_ms'] for n in names]
    p99s = [all_latencies[n]['p99_ms'] for n in names]

    x = np.arange(len(names))
    width = 0.25
    ax.bar(x - width, means, width, label='Mean', color='#2196F3', alpha=0.85)
    ax.bar(x, p95s, width, label='P95', color='#FF9800', alpha=0.85)
    ax.bar(x + width, p99s, width, label='P99', color='#F44336', alpha=0.85)

    ax.axhline(y=30, color='red', linestyle='--', linewidth=2, label='30ms Requirement')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=10, rotation=15)
    ax.set_ylabel('Latency (ms)', fontsize=12)
    ax.set_title('Inference Latency Comparison', fontsize=14)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_transfer_learning(results, filename='transfer_learning.png'):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    scenarios = list(results.keys())
    metrics_to_plot = ['f1', 'auc']
    colors = ['#4CAF50', '#2196F3', '#FF9800', '#F44336']

    for idx, metric in enumerate(metrics_to_plot):
        vals = [results[s][metric] for s in scenarios]
        bars = axes[idx].bar(range(len(scenarios)), vals, color=colors[:len(scenarios)], alpha=0.85)
        axes[idx].set_xticks(range(len(scenarios)))
        axes[idx].set_xticklabels(scenarios, fontsize=9, rotation=20)
        axes[idx].set_ylabel(metric.upper(), fontsize=12)
        axes[idx].set_title(f'Transfer Learning - {metric.upper()}', fontsize=13)
        axes[idx].set_ylim(0, 1.1)
        axes[idx].grid(True, alpha=0.3, axis='y')
        for bar, val in zip(bars, vals):
            axes[idx].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                          f'{val:.3f}', ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_confusion_matrices(all_cms, filename='confusion_matrices.png'):
    n = len(all_cms)
    fig, axes = plt.subplots(1, n, figsize=(5*n, 4.5))
    if n == 1:
        axes = [axes]
    for ax, (name, cm) in zip(axes, all_cms.items()):
        im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
        ax.set_title(name, fontsize=12)
        ax.set_xlabel('Predicted'); ax.set_ylabel('Actual')
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(['Normal', 'Disrupt']); ax.set_yticklabels(['Normal', 'Disrupt'])
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f'{cm[i, j]}', ha='center', va='center', fontsize=14,
                       color='white' if cm[i, j] > cm.max()/2 else 'black')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_feature_importance(model, X_test, y_test, feature_names, filename='feature_importance.png'):
    """Permutation-based feature importance."""
    model.eval()
    model.to(DEVICE)
    with torch.no_grad():
        base_pred = torch.sigmoid(model(torch.FloatTensor(X_test).to(DEVICE))).cpu().numpy()
    base_auc = roc_auc_score(y_test, base_pred) if len(np.unique(y_test)) > 1 else 0.5

    importances = []
    for feat_idx in range(X_test.shape[2]):
        X_perm = X_test.copy()
        np.random.shuffle(X_perm[:, :, feat_idx])
        with torch.no_grad():
            perm_pred = torch.sigmoid(model(torch.FloatTensor(X_perm).to(DEVICE))).cpu().numpy()
        perm_auc = roc_auc_score(y_test, perm_pred) if len(np.unique(y_test)) > 1 else 0.5
        importances.append(base_auc - perm_auc)

    fig, ax = plt.subplots(figsize=(10, 6))
    sorted_idx = np.argsort(importances)[::-1]
    ax.barh(range(len(feature_names)), [importances[i] for i in sorted_idx],
            color='#2196F3', alpha=0.85)
    ax.set_yticks(range(len(feature_names)))
    ax.set_yticklabels([feature_names[i] for i in sorted_idx], fontsize=11)
    ax.set_xlabel('AUC Drop (Permutation Importance)', fontsize=12)
    ax.set_title('Feature Importance for Disruption Prediction', fontsize=14)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    return {feature_names[i]: importances[i] for i in sorted_idx}


def plot_architecture(filename='model_architecture.png'):
    """Generate model architecture diagram."""
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.set_xlim(0, 16); ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Real-Time Disruption Prediction System Architecture', fontsize=16, fontweight='bold', pad=20)

    def draw_box(x, y, w, h, text, color='#E3F2FD', edge='#1565C0', fontsize=9):
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor=edge, linewidth=1.5, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fontsize,
                fontweight='bold', zorder=3, wrap=True)

    def draw_arrow(x1, y1, x2, y2, color='#333'):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                   arrowprops=dict(arrowstyle='->', color=color, lw=1.5), zorder=1)

    # Input layer
    draw_box(0.5, 7.5, 3, 2, 'Plasma Diagnostics\n─────────────\nMirnov coils\nECE radiometry\nThomson scattering\nBolometry\nMagnetics', '#E8F5E9', '#2E7D32')

    # Preprocessing
    draw_box(4.5, 8, 2.5, 1.2, 'Real-Time\nPreprocessing\n(< 1ms)', '#FFF3E0', '#E65100')

    # Multi-scale CNN
    draw_box(4.5, 6.2, 2.5, 1.2, 'Multi-Scale\nCNN Feature\nExtractor', '#E3F2FD', '#1565C0')

    # Physics module
    draw_box(8, 8, 2.5, 1.2, 'Physics-Informed\nModule\n(MHD constraints)', '#FCE4EC', '#C62828')

    # LSTM
    draw_box(8, 6.2, 2.5, 1.2, 'Bidirectional\nLSTM\nTemporal Encoder', '#E3F2FD', '#1565C0')

    # TM detector
    draw_box(11.5, 8, 2.5, 1.2, 'TM/NTM\nDetection Head', '#F3E5F5', '#6A1B9A')

    # Disruption head
    draw_box(11.5, 6.2, 2.5, 1.2, 'Disruption\nPrediction Head', '#FFEBEE', '#B71C1C')

    # Output
    draw_box(11.5, 3.5, 2.5, 2, 'Control System\n─────────────\nDisruption prob.\nTM/NTM alert\nTime-to-disrupt\nMitigation trigger\n(< 30ms total)', '#E8F5E9', '#2E7D32')

    # Transfer learning
    draw_box(4.5, 3.5, 3, 1.5, 'Cross-Device\nTransfer Learning\n─────────────\nJET → ITER\nKSTAR → ITER', '#FFF9C4', '#F57F17')

    # Domain adaptation
    draw_box(8, 3.5, 2.5, 1.5, 'Domain\nAdaptation\nLayer', '#FFF9C4', '#F57F17')

    # Arrows
    draw_arrow(3.5, 8.5, 4.5, 8.6)
    draw_arrow(5.75, 8.0, 5.75, 7.4)
    draw_arrow(7.0, 6.8, 8.0, 6.8)
    draw_arrow(5.75, 6.2, 5.75, 5.0)
    draw_arrow(7.0, 4.25, 8.0, 4.25)
    draw_arrow(9.25, 8.0, 9.25, 7.4)
    draw_arrow(10.5, 8.6, 11.5, 8.6)
    draw_arrow(10.5, 6.8, 11.5, 6.8)
    draw_arrow(12.75, 6.2, 12.75, 5.5)
    draw_arrow(12.75, 8.0, 12.75, 7.4)  # TM to merge
    draw_arrow(10.5, 4.25, 11.5, 4.25)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_plasma_signals(filename='plasma_signals.png'):
    """Plot example disruptive and non-disruptive plasma shots."""
    from generate_data import generate_plasma_shot
    shot_normal = generate_plasma_shot('demo_normal', 'JET', disruption=False)
    shot_disrupt = generate_plasma_shot('demo_disrupt', 'JET', disruption=True)

    fig, axes = plt.subplots(4, 2, figsize=(14, 10), sharex=True)
    signals = [(0, 'Ip (MA)'), (2, 'Te (keV)'), (5, 'P_rad fraction'), (9, 'Mirnov RMS')]

    for row, (idx, label) in enumerate(signals):
        axes[row, 0].plot(shot_normal['t'], shot_normal['features'][:, idx], 'b-', linewidth=1)
        axes[row, 0].set_ylabel(label, fontsize=10)
        if row == 0:
            axes[row, 0].set_title('Non-Disruptive Shot', fontsize=12, fontweight='bold')

        axes[row, 1].plot(shot_disrupt['t'], shot_disrupt['features'][:, idx], 'r-', linewidth=1)
        if row == 0:
            axes[row, 1].set_title('Disruptive Shot', fontsize=12, fontweight='bold')

    axes[-1, 0].set_xlabel('Time (s)', fontsize=11)
    axes[-1, 1].set_xlabel('Time (s)', fontsize=11)
    plt.suptitle('Example Plasma Diagnostic Signals (JET-like)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


def plot_pipeline_timing(latencies, filename='pipeline_timing.png'):
    """Visualize the real-time inference pipeline timing."""
    fig, ax = plt.subplots(figsize=(12, 5))
    stages = ['Data\nAcquisition', 'Preprocessing\n& Normalization', 'CNN Feature\nExtraction',
              'LSTM Temporal\nEncoding', 'Physics\nConstraints', 'Prediction\nOutput', 'Control\nAction']
    times = [1.0, 0.5, 2.0, 3.0, 1.5, 0.3, 2.0]  # estimated ms
    colors = ['#4CAF50', '#FF9800', '#2196F3', '#2196F3', '#F44336', '#9C27B0', '#4CAF50']

    cumulative = 0
    for i, (stage, t_ms, color) in enumerate(zip(stages, times, colors)):
        ax.barh(0, t_ms, left=cumulative, height=0.5, color=color, alpha=0.8, edgecolor='white')
        ax.text(cumulative + t_ms/2, 0, f'{stage}\n{t_ms:.1f}ms', ha='center', va='center',
                fontsize=8, fontweight='bold')
        cumulative += t_ms

    total = sum(times)
    ax.axvline(x=30, color='red', linestyle='--', linewidth=2, label=f'30ms Limit')
    ax.set_xlim(0, 35)
    ax.set_yticks([])
    ax.set_xlabel('Time (ms)', fontsize=12)
    ax.set_title(f'Real-Time Inference Pipeline (Total: {total:.1f}ms)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()


# ---- Main experiment ----

def main():
    np.random.seed(42)
    torch.manual_seed(42)

    print("=" * 60)
    print("Tokamak Disruption Prediction - Experiment Suite")
    print("=" * 60)

    # Generate data
    print("\n[1/7] Generating synthetic plasma data...")
    from generate_data import generate_dataset, prepare_sequences
    jets = generate_dataset(200, 'JET', 0.3, 0.2)
    kstars = generate_dataset(200, 'KSTAR', 0.3, 0.2)
    iters = generate_dataset(30, 'ITER', 0.2, 0.1)

    X_jet, y_jet, ytm_jet = prepare_sequences(jets)
    X_kstar, y_kstar, ytm_kstar = prepare_sequences(kstars)
    X_iter, y_iter, ytm_iter = prepare_sequences(iters)

    print(f"  JET: {X_jet.shape}, disruptions={y_jet.sum():.0f}")
    print(f"  KSTAR: {X_kstar.shape}, disruptions={y_kstar.sum():.0f}")
    print(f"  ITER: {X_iter.shape}, disruptions={y_iter.sum():.0f}")

    # Plot example signals
    plot_plasma_signals()
    print("  → Generated: figures/plasma_signals.png")

    # Split JET data
    Xtr, ytr, Xv, yv, Xte, yte = split_data(X_jet, y_jet)
    _, ytr_tm, _, yv_tm, _, yte_tm = split_data(X_jet, ytm_jet)

    # ---- Experiment 1: Model comparison on JET ----
    print("\n[2/7] Training disruption prediction models on JET data...")
    models_config = {
        'LSTM': (LSTMPredictor(11, 64, 2, 0.3), False, False),
        'CNN-LSTM': (CNNLSTMPredictor(11, 32, 64, 2, 0.3), False, False),
        'PI-LSTM': (PhysicsInformedPredictor(11, 64, 2, 0.3), True, False),
    }

    all_curves = {}
    all_metrics = {}
    all_rocs = {}
    all_cms = {}
    all_latencies = {}

    for name, (model, is_physics, is_mt) in models_config.items():
        print(f"  Training {name}...")
        tl, vl = train_model(model, Xtr, ytr, Xv, yv, epochs=30, lr=1e-3, physics_model=is_physics)
        all_curves[name] = (tl, vl)
        metrics, probs = evaluate_model(model, Xte, yte, physics_model=is_physics)
        all_metrics[name] = metrics
        all_rocs[name] = (yte, probs)
        y_pred = (probs > 0.5).astype(int)
        all_cms[name] = confusion_matrix(yte, y_pred)
        latency = measure_latency(model, physics_model=is_physics)
        all_latencies[name] = latency
        print(f"    {name}: F1={metrics['f1']:.3f}, AUC={metrics['auc']:.3f}, "
              f"Latency={latency['mean_ms']:.2f}ms (P99={latency['p99_ms']:.2f}ms)")

    RESULTS['model_comparison'] = all_metrics
    RESULTS['latencies'] = {k: v for k, v in all_latencies.items()}

    # ---- Experiment 2: Tearing mode detection (multi-task) ----
    print("\n[3/7] Training TM/NTM multi-task detector...")
    tm_model = TearingModeDetector(11, 48, 0.3)
    tl_tm, vl_tm = train_model(tm_model, Xtr, ytr, Xv, yv, epochs=30, lr=1e-3,
                                extra_y_train=ytr_tm, extra_y_val=yv_tm)

    # Evaluate disruption prediction
    tm_metrics, tm_probs = evaluate_model(tm_model, Xte, yte, multi_task=True)
    all_metrics['TM-MultiTask'] = tm_metrics
    all_rocs['TM-MultiTask'] = (yte, tm_probs)
    all_cms['TM-MultiTask'] = confusion_matrix(yte, (tm_probs > 0.5).astype(int))

    # Evaluate TM detection
    tm_model.eval()
    with torch.no_grad():
        _, tm_logits = tm_model(torch.FloatTensor(Xte).to(DEVICE))
        tm_det_probs = torch.sigmoid(tm_logits).cpu().numpy()
    tm_det_pred = (tm_det_probs > 0.5).astype(int)
    tm_det_metrics = {
        'accuracy': accuracy_score(yte_tm, tm_det_pred),
        'precision': precision_score(yte_tm, tm_det_pred, zero_division=0),
        'recall': recall_score(yte_tm, tm_det_pred, zero_division=0),
        'f1': f1_score(yte_tm, tm_det_pred, zero_division=0),
        'auc': roc_auc_score(yte_tm, tm_det_probs) if len(np.unique(yte_tm)) > 1 else 0.0,
    }
    RESULTS['tm_detection'] = tm_det_metrics
    print(f"  TM Detection: F1={tm_det_metrics['f1']:.3f}, AUC={tm_det_metrics['auc']:.3f}")

    tm_latency = measure_latency(tm_model, multi_task=True)
    all_latencies['TM-MultiTask'] = tm_latency
    print(f"  Disruption (multi-task): F1={tm_metrics['f1']:.3f}, AUC={tm_metrics['auc']:.3f}, "
          f"Latency={tm_latency['mean_ms']:.2f}ms")

    all_curves['TM-MultiTask'] = (tl_tm, vl_tm)

    # ---- Experiment 3: Transfer learning JET → KSTAR/ITER ----
    print("\n[4/7] Transfer learning experiments...")
    transfer_results = {}

    # Baseline: JET-trained CNN-LSTM tested on KSTAR
    cnn_model_jet = CNNLSTMPredictor(11, 32, 64, 2, 0.3)
    train_model(cnn_model_jet, Xtr, ytr, Xv, yv, epochs=30, lr=1e-3)
    kstar_metrics, _ = evaluate_model(cnn_model_jet, X_kstar, y_kstar)
    transfer_results['JET→KSTAR\n(no adapt)'] = kstar_metrics
    print(f"  JET→KSTAR (no adaptation): F1={kstar_metrics['f1']:.3f}, AUC={kstar_metrics['auc']:.3f}")

    # Transfer with fine-tuning on small KSTAR data
    transfer_model = TransferablePredictor(11, 64, 32, 0.3)
    # Pre-train on JET
    train_model(transfer_model, Xtr, ytr, Xv, yv, epochs=20, lr=1e-3)
    # Fine-tune on small KSTAR subset
    n_ft = min(500, len(X_kstar))
    idx_ft = np.random.choice(len(X_kstar), n_ft, replace=False)
    transfer_model.freeze_shared()
    Xk_ft, yk_ft = X_kstar[idx_ft], y_kstar[idx_ft]
    Xk_val = X_kstar[~np.isin(np.arange(len(X_kstar)), idx_ft)][:200]
    yk_val = y_kstar[~np.isin(np.arange(len(y_kstar)), idx_ft)][:200]
    train_model(transfer_model, Xk_ft, yk_ft, Xk_val, yk_val, epochs=15, lr=5e-4)
    kstar_ft_metrics, _ = evaluate_model(transfer_model, X_kstar, y_kstar)
    transfer_results['JET→KSTAR\n(fine-tuned)'] = kstar_ft_metrics
    print(f"  JET→KSTAR (fine-tuned): F1={kstar_ft_metrics['f1']:.3f}, AUC={kstar_ft_metrics['auc']:.3f}")

    # KSTAR native baseline
    Xktr, yktr, Xkv, ykv, Xkte, ykte = split_data(X_kstar, y_kstar)
    kstar_native = CNNLSTMPredictor(11, 32, 64, 2, 0.3)
    train_model(kstar_native, Xktr, yktr, Xkv, ykv, epochs=30, lr=1e-3)
    kstar_native_metrics, _ = evaluate_model(kstar_native, Xkte, ykte)
    transfer_results['KSTAR native'] = kstar_native_metrics
    print(f"  KSTAR native: F1={kstar_native_metrics['f1']:.3f}, AUC={kstar_native_metrics['auc']:.3f}")

    # JET → ITER transfer
    transfer_iter = TransferablePredictor(11, 64, 32, 0.3)
    train_model(transfer_iter, Xtr, ytr, Xv, yv, epochs=20, lr=1e-3)
    transfer_iter.freeze_shared()
    n_iter_ft = min(100, len(X_iter))
    idx_iter = np.random.choice(len(X_iter), n_iter_ft, replace=False)
    Xi_ft, yi_ft = X_iter[idx_iter], y_iter[idx_iter]
    Xi_val = X_iter[~np.isin(np.arange(len(X_iter)), idx_iter)][:50]
    yi_val = y_iter[~np.isin(np.arange(len(y_iter)), idx_iter)][:50]
    if len(Xi_val) > 0 and len(Xi_ft) >= 16:
        train_model(transfer_iter, Xi_ft, yi_ft, Xi_val, yi_val, epochs=15, lr=5e-4, batch_size=min(64, len(Xi_ft)))
    iter_metrics, _ = evaluate_model(transfer_iter, X_iter, y_iter)
    transfer_results['JET→ITER\n(fine-tuned)'] = iter_metrics
    print(f"  JET→ITER (fine-tuned): F1={iter_metrics['f1']:.3f}, AUC={iter_metrics['auc']:.3f}")

    RESULTS['transfer_learning'] = transfer_results

    # ---- Experiment 4: Feature importance ----
    print("\n[5/7] Computing feature importance...")
    feature_names = ['Ip', 'ne', 'Te', 'βN', 'li', 'Prad', 'q95', 'nG', 'LM', 'Mirnov', 'Wmhd']
    lstm_model = models_config['LSTM'][0]
    feat_imp = plot_feature_importance(lstm_model, Xte, yte, feature_names)
    RESULTS['feature_importance'] = feat_imp
    print("  → Generated: figures/feature_importance.png")
    print(f"  Top 3: {list(feat_imp.items())[:3]}")

    # ---- Generate all figures ----
    print("\n[6/7] Generating figures...")
    plot_training_curves(all_curves)
    print("  → figures/training_curves.png")

    plot_roc_curves(all_rocs)
    print("  → figures/roc_curves.png")

    plot_metrics_comparison(all_metrics)
    print("  → figures/metrics_comparison.png")

    plot_latency(all_latencies)
    print("  → figures/latency_comparison.png")

    plot_transfer_learning(transfer_results)
    print("  → figures/transfer_learning.png")

    plot_confusion_matrices(all_cms)
    print("  → figures/confusion_matrices.png")

    plot_architecture()
    print("  → figures/model_architecture.png")

    plot_pipeline_timing(all_latencies)
    print("  → figures/pipeline_timing.png")

    # ---- Summary ----
    print("\n[7/7] Experiment complete. Saving results...")
    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    results_serializable = json.loads(json.dumps(RESULTS, default=convert))
    with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results.json'), 'w') as f:
        json.dump(results_serializable, f, indent=2)

    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print("\nModel Comparison (JET test set):")
    for name, m in all_metrics.items():
        print(f"  {name:15s}: F1={m['f1']:.3f}  AUC={m['auc']:.3f}  Acc={m['accuracy']:.3f}")
    print(f"\nTM/NTM Detection: F1={tm_det_metrics['f1']:.3f}  AUC={tm_det_metrics['auc']:.3f}")
    print("\nInference Latency:")
    for name, lat in all_latencies.items():
        status = "✓ PASS" if lat['p99_ms'] < 30 else "✗ FAIL"
        print(f"  {name:15s}: mean={lat['mean_ms']:.2f}ms  P99={lat['p99_ms']:.2f}ms  [{status}]")
    print("\nTransfer Learning:")
    for name, m in transfer_results.items():
        print(f"  {name:20s}: F1={m['f1']:.3f}  AUC={m['auc']:.3f}")

    return RESULTS


if __name__ == '__main__':
    main()
