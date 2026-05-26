"""
Run all experiments and generate figures for the tactile sensing research.
"""

import torch
import torch.nn.functional as F
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from sklearn.metrics import confusion_matrix
import os
import sys
import warnings
warnings.filterwarnings('ignore')

from tactile_framework import (
    TactileSimulator, TactileDataset,
    ContactEstimationNet, TextureClassifier, MultimodalFusionNet,
    GraspStabilityNet, SlipDetectionNet, ExploratoryGraspPolicy,
    ForceController, TactileGraspEnv,
    train_contact_estimation, train_texture_classifier,
    train_multimodal, train_slip_detector, train_grasp_stability,
    run_force_control_simulation, run_exploratory_grasp
)
from torch.utils.data import DataLoader, random_split

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)
np.random.seed(42)
torch.manual_seed(42)

sns.set_theme(style='whitegrid', font_scale=1.1)
COLORS = sns.color_palette('Set2', 10)


def plot_tactile_samples():
    """Generate sample tactile images for different shapes and textures."""
    print("Generating tactile sample visualizations...")
    sim = TactileSimulator(resolution=64)
    fig, axes = plt.subplots(4, 5, figsize=(18, 14))
    shapes = ['sphere', 'cylinder', 'edge', 'flat']
    
    for i, shape in enumerate(shapes):
        params = {'center': (0, 0), 'radius': 0.4, 'force': 1.0, 'angle': 0.3}
        depth = sim.generate_contact_geometry(shape, params)
        tactile = sim.render_tactile_image(depth)
        force = sim.compute_force_distribution(depth)
        normal = sim.depth_to_normal(depth)
        
        axes[i, 0].imshow(tactile)
        axes[i, 0].set_title(f'{shape.capitalize()} - Tactile', fontsize=11)
        axes[i, 0].axis('off')
        
        im1 = axes[i, 1].imshow(depth, cmap='viridis')
        axes[i, 1].set_title('Depth Map', fontsize=11)
        axes[i, 1].axis('off')
        plt.colorbar(im1, ax=axes[i, 1], fraction=0.046)
        
        im2 = axes[i, 2].imshow(force, cmap='hot')
        axes[i, 2].set_title('Force Distribution', fontsize=11)
        axes[i, 2].axis('off')
        plt.colorbar(im2, ax=axes[i, 2], fraction=0.046)
        
        normal_vis = (normal + 1) / 2
        axes[i, 3].imshow(normal_vis)
        axes[i, 3].set_title('Normal Map', fontsize=11)
        axes[i, 3].axis('off')
        
        # Texture overlay
        tex = sim.generate_texture('striped' if i < 2 else 'dotted', frequency=8)
        depth_tex = depth + tex * (depth > 0.01).astype(float)
        tactile_tex = sim.render_tactile_image(depth_tex)
        axes[i, 4].imshow(tactile_tex)
        axes[i, 4].set_title('With Texture', fontsize=11)
        axes[i, 4].axis('off')

    fig.suptitle('GelSight/DIGIT Tactile Sensor Simulation: Contact Geometry & Force Estimation',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/tactile_samples.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/tactile_samples.png")


def plot_texture_gallery():
    """Generate a gallery of all texture types."""
    print("Generating texture gallery...")
    sim = TactileSimulator(resolution=64)
    textures = ['smooth', 'rough', 'striped', 'dotted', 'crosshatch', 'wavy', 'grid', 'random_bumps']
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    
    for i, tex_name in enumerate(textures):
        row, col = i // 4, i % 4
        params = {'center': (0, 0), 'radius': 0.5, 'force': 1.0, 'angle': 0}
        depth = sim.generate_contact_geometry('sphere', params)
        texture = sim.generate_texture(tex_name, frequency=8)
        depth_tex = depth + texture * (depth > 0.01).astype(float)
        tactile = sim.render_tactile_image(depth_tex)
        axes[row, col].imshow(tactile)
        axes[row, col].set_title(tex_name.capitalize(), fontsize=12, fontweight='bold')
        axes[row, col].axis('off')
    
    fig.suptitle('Texture Classification: 8 Material Categories', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/texture_gallery.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/texture_gallery.png")


def run_experiment_1():
    """Contact shape and force distribution estimation."""
    print("\n=== Experiment 1: Contact Shape & Force Estimation ===")
    dataset = TactileDataset(num_samples=1500, resolution=64, include_visual=False)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_ds, test_ds = random_split(dataset, [train_size, test_size])
    train_loader = DataLoader(train_ds, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=32, shuffle=False)

    model = ContactEstimationNet()
    history = train_contact_estimation(model, train_loader, epochs=25, lr=1e-3, device=DEVICE)

    # Evaluate
    model.eval()
    depth_errors, force_errors = [], []
    sample_results = []
    with torch.no_grad():
        for batch in test_loader:
            tactile = batch['tactile_image'].to(DEVICE)
            gt_depth = batch['depth_map'].to(DEVICE)
            gt_force = batch['force_distribution'].to(DEVICE)
            pred_depth, pred_force = model(tactile)
            depth_errors.append(F.mse_loss(pred_depth, gt_depth, reduction='none').mean(dim=[1,2,3]).cpu().numpy())
            force_errors.append(F.mse_loss(pred_force, gt_force, reduction='none').mean(dim=[1,2,3]).cpu().numpy())
            if len(sample_results) == 0:
                sample_results = {
                    'tactile': tactile[:4].cpu(), 'gt_depth': gt_depth[:4].cpu(),
                    'pred_depth': pred_depth[:4].cpu(), 'gt_force': gt_force[:4].cpu(),
                    'pred_force': pred_force[:4].cpu()
                }

    depth_errors = np.concatenate(depth_errors)
    force_errors = np.concatenate(force_errors)

    # Plot training curves
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].plot(history['depth_loss'], color=COLORS[0], linewidth=2, label='Depth Loss')
    axes[0].plot(history['force_loss'], color=COLORS[1], linewidth=2, label='Force Loss')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('MSE Loss')
    axes[0].set_title('Training Loss Curves')
    axes[0].legend()
    axes[0].set_yscale('log')

    axes[1].hist(depth_errors, bins=30, color=COLORS[0], alpha=0.7, edgecolor='black')
    axes[1].axvline(np.mean(depth_errors), color='red', linestyle='--',
                    label=f'Mean={np.mean(depth_errors):.6f}')
    axes[1].set_xlabel('Depth MSE')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Depth Estimation Error Distribution')
    axes[1].legend()

    axes[2].hist(force_errors, bins=30, color=COLORS[1], alpha=0.7, edgecolor='black')
    axes[2].axvline(np.mean(force_errors), color='red', linestyle='--',
                    label=f'Mean={np.mean(force_errors):.6f}')
    axes[2].set_xlabel('Force MSE')
    axes[2].set_ylabel('Count')
    axes[2].set_title('Force Estimation Error Distribution')
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/contact_estimation_training.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/contact_estimation_training.png")

    # Plot qualitative results
    fig, axes = plt.subplots(4, 5, figsize=(20, 16))
    col_titles = ['Tactile Input', 'GT Depth', 'Pred Depth', 'GT Force', 'Pred Force']
    for j, title in enumerate(col_titles):
        axes[0, j].set_title(title, fontsize=12, fontweight='bold')
    for i in range(4):
        axes[i, 0].imshow(sample_results['tactile'][i].permute(1, 2, 0).numpy())
        axes[i, 1].imshow(sample_results['gt_depth'][i, 0].numpy(), cmap='viridis')
        axes[i, 2].imshow(sample_results['pred_depth'][i, 0].numpy(), cmap='viridis')
        axes[i, 3].imshow(sample_results['gt_force'][i, 0].numpy(), cmap='hot')
        axes[i, 4].imshow(sample_results['pred_force'][i, 0].numpy(), cmap='hot')
        for j in range(5):
            axes[i, j].axis('off')

    fig.suptitle('Contact Shape & Force Distribution Estimation Results', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/contact_estimation_qualitative.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/contact_estimation_qualitative.png")

    results = {
        'depth_mse_mean': float(np.mean(depth_errors)),
        'depth_mse_std': float(np.std(depth_errors)),
        'force_mse_mean': float(np.mean(force_errors)),
        'force_mse_std': float(np.std(force_errors)),
        'final_train_loss': history['total_loss'][-1],
    }
    print(f"  Depth MSE: {results['depth_mse_mean']:.6f} ± {results['depth_mse_std']:.6f}")
    print(f"  Force MSE: {results['force_mse_mean']:.6f} ± {results['force_mse_std']:.6f}")
    return results


def run_experiment_2():
    """Texture classification."""
    print("\n=== Experiment 2: Texture Classification ===")
    dataset = TactileDataset(num_samples=2000, resolution=64, include_visual=False)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_ds, test_ds = random_split(dataset, [train_size, test_size])
    train_loader = DataLoader(train_ds, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=32, shuffle=False)

    model = TextureClassifier(num_classes=8)
    history = train_texture_classifier(model, train_loader, epochs=30, lr=1e-3, device=DEVICE)

    # Evaluate
    model.eval()
    all_preds, all_labels = [], []
    with torch.no_grad():
        for batch in test_loader:
            tactile = batch['tactile_image'].to(DEVICE)
            labels = batch['texture_label']
            logits = model(tactile)
            all_preds.extend(logits.argmax(1).cpu().numpy())
            all_labels.extend(labels.numpy())

    all_preds, all_labels = np.array(all_preds), np.array(all_labels)
    test_acc = (all_preds == all_labels).mean()

    # Per-class accuracy
    texture_names = ['Smooth', 'Rough', 'Striped', 'Dotted', 'CrossH', 'Wavy', 'Grid', 'RndBump']
    cm = confusion_matrix(all_labels, all_preds, labels=range(8))
    per_class_acc = cm.diagonal() / (cm.sum(axis=1) + 1e-8)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    axes[0].plot(history['loss'], color=COLORS[0], linewidth=2, label='Loss')
    ax_acc = axes[0].twinx()
    ax_acc.plot(history['accuracy'], color=COLORS[1], linewidth=2, label='Accuracy')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Loss')
    ax_acc.set_ylabel('Accuracy')
    axes[0].set_title('Training Loss & Accuracy')
    lines1, labels1 = axes[0].get_legend_handles_labels()
    lines2, labels2 = ax_acc.get_legend_handles_labels()
    axes[0].legend(lines1 + lines2, labels1 + labels2)

    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1],
                xticklabels=texture_names, yticklabels=texture_names)
    axes[1].set_xlabel('Predicted')
    axes[1].set_ylabel('True')
    axes[1].set_title(f'Confusion Matrix (Acc={test_acc:.3f})')

    bars = axes[2].bar(texture_names, per_class_acc, color=COLORS[:8], edgecolor='black')
    axes[2].set_ylabel('Accuracy')
    axes[2].set_title('Per-Class Accuracy')
    axes[2].set_ylim(0, 1.1)
    for bar, acc in zip(bars, per_class_acc):
        axes[2].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                     f'{acc:.2f}', ha='center', fontsize=9)
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/texture_classification.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/texture_classification.png")

    results = {
        'test_accuracy': float(test_acc),
        'per_class_accuracy': {n: float(a) for n, a in zip(texture_names, per_class_acc)},
        'final_train_loss': history['loss'][-1],
        'final_train_acc': history['accuracy'][-1],
    }
    print(f"  Test Accuracy: {test_acc:.4f}")
    return results


def run_experiment_3():
    """Multimodal fusion experiment."""
    print("\n=== Experiment 3: Multimodal Tactile-Visual Fusion ===")
    dataset = TactileDataset(num_samples=2000, resolution=64, include_visual=True)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_ds, test_ds = random_split(dataset, [train_size, test_size])
    train_loader = DataLoader(train_ds, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=32, shuffle=False)

    # Train multimodal
    mm_model = MultimodalFusionNet()
    mm_history = train_multimodal(mm_model, train_loader, epochs=30, lr=1e-3, device=DEVICE)

    # Train tactile-only baseline (TextureClassifier for shape)
    tactile_only = TextureClassifier(num_classes=4)
    to_history = {'loss': [], 'accuracy': []}
    tactile_only.to(DEVICE)
    optimizer = torch.optim.Adam(tactile_only.parameters(), lr=1e-3)
    criterion = torch.nn.CrossEntropyLoss()
    for epoch in range(30):
        tactile_only.train()
        total_loss, correct, total = 0, 0, 0
        for batch in train_loader:
            x = batch['tactile_image'].to(DEVICE)
            y = batch['shape_label'].to(DEVICE)
            logits = tactile_only(x)
            loss = criterion(logits, y)
            optimizer.zero_grad(); loss.backward(); optimizer.step()
            total_loss += loss.item()
            correct += (logits.argmax(1) == y).sum().item()
            total += y.size(0)
        to_history['loss'].append(total_loss / len(train_loader))
        to_history['accuracy'].append(correct / total)

    # Evaluate
    mm_model.eval()
    tactile_only.eval()
    mm_shape_preds, mm_tex_preds = [], []
    to_preds = []
    shape_labels_all, tex_labels_all = [], []

    with torch.no_grad():
        for batch in test_loader:
            t = batch['tactile_image'].to(DEVICE)
            v = batch['visual_image'].to(DEVICE)
            s_logits, t_logits = mm_model(t, v)
            mm_shape_preds.extend(s_logits.argmax(1).cpu().numpy())
            mm_tex_preds.extend(t_logits.argmax(1).cpu().numpy())
            to_logits = tactile_only(t)
            to_preds.extend(to_logits.argmax(1).cpu().numpy())
            shape_labels_all.extend(batch['shape_label'].numpy())
            tex_labels_all.extend(batch['texture_label'].numpy())

    mm_shape_acc = (np.array(mm_shape_preds) == np.array(shape_labels_all)).mean()
    mm_tex_acc = (np.array(mm_tex_preds) == np.array(tex_labels_all)).mean()
    to_shape_acc = (np.array(to_preds) == np.array(shape_labels_all)).mean()

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].plot(mm_history['shape_acc'], color=COLORS[0], linewidth=2, label='MM Shape')
    axes[0].plot(mm_history['texture_acc'], color=COLORS[1], linewidth=2, label='MM Texture')
    axes[0].plot(to_history['accuracy'], color=COLORS[2], linewidth=2, linestyle='--',
                 label='Tactile-Only Shape')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Accuracy')
    axes[0].set_title('Training Accuracy Comparison')
    axes[0].legend()

    methods = ['Multimodal\n(Shape)', 'Multimodal\n(Texture)', 'Tactile-Only\n(Shape)']
    accs = [mm_shape_acc, mm_tex_acc, to_shape_acc]
    bars = axes[1].bar(methods, accs, color=[COLORS[0], COLORS[1], COLORS[2]], edgecolor='black')
    axes[1].set_ylabel('Test Accuracy')
    axes[1].set_title('Modality Comparison')
    axes[1].set_ylim(0, 1.15)
    for bar, acc in zip(bars, accs):
        axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                     f'{acc:.3f}', ha='center', fontsize=11, fontweight='bold')

    axes[2].plot(mm_history['loss'], color=COLORS[0], linewidth=2, label='Multimodal')
    axes[2].plot(to_history['loss'], color=COLORS[2], linewidth=2, linestyle='--', label='Tactile-Only')
    axes[2].set_xlabel('Epoch')
    axes[2].set_ylabel('Loss')
    axes[2].set_title('Training Loss Comparison')
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/multimodal_fusion.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/multimodal_fusion.png")

    results = {
        'mm_shape_acc': float(mm_shape_acc),
        'mm_texture_acc': float(mm_tex_acc),
        'tactile_only_shape_acc': float(to_shape_acc),
        'improvement': float(mm_shape_acc - to_shape_acc),
    }
    print(f"  Multimodal Shape Acc: {mm_shape_acc:.4f}")
    print(f"  Multimodal Texture Acc: {mm_tex_acc:.4f}")
    print(f"  Tactile-Only Shape Acc: {to_shape_acc:.4f}")
    return results


def run_experiment_4():
    """Grasp stability evaluation."""
    print("\n=== Experiment 4: Grasp Stability Evaluation ===")
    dataset = TactileDataset(num_samples=1500, resolution=64, include_visual=False)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_ds, test_ds = random_split(dataset, [train_size, test_size])
    train_loader = DataLoader(train_ds, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=32, shuffle=False)

    model = GraspStabilityNet()
    history = train_grasp_stability(model, train_loader, epochs=25, lr=1e-3, device=DEVICE)

    # Evaluate
    model.eval()
    all_preds, all_gts = [], []
    with torch.no_grad():
        for batch in test_loader:
            tactile = batch['tactile_image'].to(DEVICE)
            gt = batch['stability_score']
            pred = model(tactile)
            all_preds.extend(pred.cpu().numpy())
            all_gts.extend(gt.numpy())

    all_preds, all_gts = np.array(all_preds), np.array(all_gts)
    mae = np.abs(all_preds - all_gts).mean()
    correlation = np.corrcoef(all_preds, all_gts)[0, 1]

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].plot(history['loss'], color=COLORS[0], linewidth=2, label='MSE Loss')
    axes[0].plot(history['mae'], color=COLORS[1], linewidth=2, label='MAE')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Loss')
    axes[0].set_title('Training Curves')
    axes[0].legend()

    axes[1].scatter(all_gts, all_preds, alpha=0.4, s=20, color=COLORS[0])
    axes[1].plot([0, 1], [0, 1], 'r--', linewidth=2, label='Ideal')
    axes[1].set_xlabel('Ground Truth Stability')
    axes[1].set_ylabel('Predicted Stability')
    axes[1].set_title(f'Prediction vs GT (r={correlation:.3f})')
    axes[1].legend()
    axes[1].set_xlim(-0.1, 1.1)
    axes[1].set_ylim(-0.1, 1.1)

    errors = all_preds - all_gts
    axes[2].hist(errors, bins=40, color=COLORS[2], alpha=0.7, edgecolor='black')
    axes[2].axvline(0, color='red', linestyle='--')
    axes[2].set_xlabel('Prediction Error')
    axes[2].set_ylabel('Count')
    axes[2].set_title(f'Error Distribution (MAE={mae:.4f})')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/grasp_stability.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/grasp_stability.png")

    results = {
        'mae': float(mae),
        'correlation': float(correlation),
        'final_train_loss': history['loss'][-1],
    }
    print(f"  MAE: {mae:.4f}")
    print(f"  Correlation: {correlation:.4f}")
    return results


def run_experiment_5():
    """Slip detection and force control."""
    print("\n=== Experiment 5: Slip Detection & Force Control ===")
    dataset = TactileDataset(num_samples=2000, resolution=64, include_visual=False)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_ds, test_ds = random_split(dataset, [train_size, test_size])
    train_loader = DataLoader(train_ds, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=32, shuffle=False)

    model = SlipDetectionNet()
    history = train_slip_detector(model, train_loader, epochs=30, lr=1e-3, device=DEVICE)

    # Force control simulation
    fc_results = run_force_control_simulation(n_steps=300)

    # Evaluate slip detection
    model.eval()
    all_preds, all_labels = [], []
    with torch.no_grad():
        for batch in test_loader:
            tactile = batch['tactile_image'].to(DEVICE)
            labels = batch['slip_label']
            logits = model(tactile)
            all_preds.extend(logits.argmax(1).cpu().numpy())
            all_labels.extend(labels.numpy())

    all_preds, all_labels = np.array(all_preds), np.array(all_labels)
    test_acc = (all_preds == all_labels).mean()
    tp = ((all_preds == 1) & (all_labels == 1)).sum()
    fp = ((all_preds == 1) & (all_labels == 0)).sum()
    fn = ((all_preds == 0) & (all_labels == 1)).sum()
    precision = tp / (tp + fp + 1e-8)
    recall = tp / (tp + fn + 1e-8)
    f1 = 2 * precision * recall / (precision + recall + 1e-8)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Slip detection training
    axes[0, 0].plot(history['accuracy'], color=COLORS[0], linewidth=2, label='Accuracy')
    axes[0, 0].plot(history['precision'], color=COLORS[1], linewidth=2, label='Precision')
    axes[0, 0].plot(history['recall'], color=COLORS[2], linewidth=2, label='Recall')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Score')
    axes[0, 0].set_title('Slip Detection Training Metrics')
    axes[0, 0].legend()

    # Confusion matrix
    cm = confusion_matrix(all_labels, all_preds)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Oranges', ax=axes[0, 1],
                xticklabels=['No Slip', 'Slip'], yticklabels=['No Slip', 'Slip'])
    axes[0, 1].set_xlabel('Predicted')
    axes[0, 1].set_ylabel('True')
    axes[0, 1].set_title(f'Slip Detection CM (F1={f1:.3f})')

    # Force control
    steps = range(len(fc_results['forces']))
    axes[1, 0].plot(steps, fc_results['forces'], color=COLORS[0], linewidth=1.5, label='Actual Force')
    axes[1, 0].plot(steps, fc_results['targets'], color=COLORS[1], linewidth=1.5,
                    linestyle='--', label='Target Force')
    slip_regions = np.array(fc_results['slips'])
    axes[1, 0].fill_between(steps, 0, max(fc_results['targets']) * 1.2,
                             where=slip_regions, alpha=0.2, color='red', label='Slip Region')
    axes[1, 0].set_xlabel('Time Step')
    axes[1, 0].set_ylabel('Force (N)')
    axes[1, 0].set_title('Force Control with Slip Compensation')
    axes[1, 0].legend()

    # Control signal
    axes[1, 1].plot(steps, fc_results['controls'], color=COLORS[3], linewidth=1.5)
    axes[1, 1].fill_between(steps, 0, max(fc_results['controls']) * 1.2,
                             where=slip_regions, alpha=0.2, color='red', label='Slip Region')
    axes[1, 1].set_xlabel('Time Step')
    axes[1, 1].set_ylabel('Control Signal')
    axes[1, 1].set_title('PID Control Output')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/slip_detection_force_control.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/slip_detection_force_control.png")

    results = {
        'test_accuracy': float(test_acc),
        'precision': float(precision),
        'recall': float(recall),
        'f1_score': float(f1),
    }
    print(f"  Accuracy: {test_acc:.4f}, Precision: {precision:.4f}, "
          f"Recall: {recall:.4f}, F1: {f1:.4f}")
    return results


def run_experiment_6():
    """Exploratory grasping strategy."""
    print("\n=== Experiment 6: Exploratory Grasping Strategy ===")
    env = TactileGraspEnv(resolution=64)
    policy = ExploratoryGraspPolicy(action_dim=6)
    results = run_exploratory_grasp(env, policy, n_episodes=80, device=DEVICE)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    episodes = range(len(results['rewards']))
    axes[0, 0].plot(episodes, results['rewards'], color=COLORS[0], alpha=0.5, linewidth=1)
    window = 10
    if len(results['rewards']) > window:
        moving_avg = np.convolve(results['rewards'], np.ones(window)/window, mode='valid')
        axes[0, 0].plot(range(window-1, len(results['rewards'])), moving_avg,
                        color=COLORS[0], linewidth=2.5, label=f'{window}-ep MA')
    axes[0, 0].set_xlabel('Episode')
    axes[0, 0].set_ylabel('Cumulative Reward')
    axes[0, 0].set_title('Exploration Reward per Episode')
    axes[0, 0].legend()

    axes[0, 1].plot(episodes, results['stabilities'], color=COLORS[1], alpha=0.5, linewidth=1)
    if len(results['stabilities']) > window:
        ma = np.convolve(results['stabilities'], np.ones(window)/window, mode='valid')
        axes[0, 1].plot(range(window-1, len(results['stabilities'])), ma,
                        color=COLORS[1], linewidth=2.5, label=f'{window}-ep MA')
    axes[0, 1].set_xlabel('Episode')
    axes[0, 1].set_ylabel('Max Stability')
    axes[0, 1].set_title('Grasp Stability Achievement')
    axes[0, 1].legend()

    axes[1, 0].scatter(results['forces'], results['stabilities'],
                        c=results['rewards'], cmap='viridis', alpha=0.6, s=40)
    axes[1, 0].set_xlabel('Applied Force')
    axes[1, 0].set_ylabel('Stability Score')
    axes[1, 0].set_title('Force vs Stability (color=reward)')
    plt.colorbar(axes[1, 0].collections[0], ax=axes[1, 0], label='Reward')

    axes[1, 1].hist(results['steps'], bins=20, color=COLORS[3], edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('Steps per Episode')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title('Episode Length Distribution')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/exploratory_grasping.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/exploratory_grasping.png")

    res = {
        'mean_reward': float(np.mean(results['rewards'])),
        'mean_stability': float(np.mean(results['stabilities'])),
        'mean_steps': float(np.mean(results['steps'])),
        'max_stability': float(np.max(results['stabilities'])),
    }
    print(f"  Mean Reward: {res['mean_reward']:.4f}")
    print(f"  Mean Stability: {res['mean_stability']:.4f}")
    return res


def plot_architecture_diagram():
    """Generate system architecture diagram."""
    print("\nGenerating architecture diagram...")
    fig, ax = plt.subplots(1, 1, figsize=(18, 10))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 10)
    ax.axis('off')

    def draw_box(x, y, w, h, text, color, fontsize=9):
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                              facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.85)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', wrap=True)

    def draw_arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    # Sensors
    draw_box(0.5, 7.5, 3, 1.8, 'GelSight/DIGIT\nTactile Sensor', '#FFB3BA', 10)
    draw_box(0.5, 5.2, 3, 1.8, 'RGB-D Camera\n(Visual)', '#BAFFC9', 10)

    # Processing
    draw_box(5, 8, 3.5, 1.2, 'Contact Shape &\nForce Estimation', '#BAE1FF', 9)
    draw_box(5, 6.3, 3.5, 1.2, 'Texture\nClassification', '#FFFFBA', 9)
    draw_box(5, 4.6, 3.5, 1.2, 'Multimodal\nFusion (Attention)', '#E8BAFF', 9)
    draw_box(5, 2.9, 3.5, 1.2, 'Slip Detection\nNetwork', '#FFD4BA', 9)

    # Control
    draw_box(10, 7.2, 3.5, 1.5, 'Grasp Stability\nEvaluator', '#C9FFE5', 10)
    draw_box(10, 5, 3.5, 1.5, 'Force Controller\n(PID + Slip Comp.)', '#FFC9E5', 10)
    draw_box(10, 2.8, 3.5, 1.5, 'Exploratory Grasp\nPolicy (PPO)', '#E5C9FF', 10)

    # Output
    draw_box(14.5, 5.5, 3, 2.5, 'Robot Gripper\nAction\nExecution', '#C9D4FF', 11)

    # Arrows
    draw_arrow(3.5, 8.4, 5, 8.6)
    draw_arrow(3.5, 8.0, 5, 6.9)
    draw_arrow(3.5, 6.5, 5, 5.2)
    draw_arrow(3.5, 6.0, 5, 3.5)
    draw_arrow(8.5, 8.4, 10, 8.0)
    draw_arrow(8.5, 6.9, 10, 7.5)
    draw_arrow(8.5, 5.2, 10, 5.5)
    draw_arrow(8.5, 3.5, 10, 3.5)
    draw_arrow(13.5, 7.9, 14.5, 7.0)
    draw_arrow(13.5, 5.7, 14.5, 6.5)
    draw_arrow(13.5, 3.5, 14.5, 5.8)

    ax.set_title('Tactile Sensing System Architecture for Object Recognition & Manipulation',
                 fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/system_architecture.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/system_architecture.png")


def plot_summary_comparison(all_results):
    """Generate summary comparison figure."""
    print("\nGenerating summary comparison...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Performance metrics
    metrics = ['Contact\nDepth MSE', 'Contact\nForce MSE', 'Texture\nAccuracy',
               'MM Shape\nAccuracy', 'MM Texture\nAccuracy', 'Slip\nF1-Score',
               'Stability\nCorrelation']
    values = [
        1 - min(1, all_results['exp1']['depth_mse_mean'] * 100),
        1 - min(1, all_results['exp1']['force_mse_mean'] * 100),
        all_results['exp2']['test_accuracy'],
        all_results['exp3']['mm_shape_acc'],
        all_results['exp3']['mm_texture_acc'],
        all_results['exp5']['f1_score'],
        all_results['exp4']['correlation'],
    ]
    colors_bar = [COLORS[i] for i in range(len(metrics))]
    bars = axes[0].barh(metrics, values, color=colors_bar, edgecolor='black', height=0.6)
    axes[0].set_xlim(0, 1.15)
    axes[0].set_xlabel('Performance Score')
    axes[0].set_title('System Performance Summary', fontweight='bold')
    for bar, val in zip(bars, values):
        axes[0].text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2,
                     f'{val:.3f}', va='center', fontsize=10)

    # Model complexity
    models = {
        'ContactEst': sum(p.numel() for p in ContactEstimationNet().parameters()),
        'TextureClf': sum(p.numel() for p in TextureClassifier().parameters()),
        'Multimodal': sum(p.numel() for p in MultimodalFusionNet().parameters()),
        'GraspStab': sum(p.numel() for p in GraspStabilityNet().parameters()),
        'SlipDet': sum(p.numel() for p in SlipDetectionNet().parameters()),
        'ExplPolicy': sum(p.numel() for p in ExploratoryGraspPolicy().parameters()),
    }
    bars2 = axes[1].bar(models.keys(), [v/1000 for v in models.values()],
                        color=COLORS[:6], edgecolor='black')
    axes[1].set_ylabel('Parameters (×1000)')
    axes[1].set_title('Model Complexity', fontweight='bold')
    for bar, val in zip(bars2, models.values()):
        axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 5,
                     f'{val/1000:.0f}K', ha='center', fontsize=9)
    plt.xticks(rotation=30, ha='right')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/summary_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> figures/summary_comparison.png")


if __name__ == '__main__':
    print("=" * 60)
    print("TACTILE SENSING RESEARCH: FULL EXPERIMENT SUITE")
    print("=" * 60)

    # Generate sample visualizations
    plot_tactile_samples()
    plot_texture_gallery()
    plot_architecture_diagram()

    # Run all experiments
    all_results = {}
    all_results['exp1'] = run_experiment_1()
    all_results['exp2'] = run_experiment_2()
    all_results['exp3'] = run_experiment_3()
    all_results['exp4'] = run_experiment_4()
    all_results['exp5'] = run_experiment_5()
    all_results['exp6'] = run_experiment_6()

    # Summary
    plot_summary_comparison(all_results)

    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 60)

    # Print summary
    for key, val in all_results.items():
        print(f"\n{key}:")
        for k, v in val.items():
            if isinstance(v, dict):
                continue
            print(f"  {k}: {v}")

    # Save results
    import json
    with open('experiment_results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print("\nResults saved to experiment_results.json")
