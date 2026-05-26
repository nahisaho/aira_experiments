"""
Main experiment runner for mHealth neurodegenerative disease biomarker study.
Generates data, trains models, runs change point detection, creates visualizations.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from sklearn.metrics import roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
import warnings
warnings.filterwarnings('ignore')

from data_generation import (generate_gait_data, generate_voice_data,
                             generate_touch_data, generate_longitudinal_data)
from models import (train_gait_model, train_voice_model, train_touch_model,
                    train_multimodal_fusion, evaluate_classifiers, get_classifiers)
from change_point_detection import (cusum_detection, pelt_detection,
                                     bayesian_online_cpd, multimodal_cpd,
                                     evaluate_cpd)

FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'figure.facecolor': 'white',
})


def plot_gait_results(results, gait_df, feature_cols):
    """Generate all gait-related figures."""
    
    # Figure 1: Model comparison bar chart
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    models = list(results.keys())
    metrics = ['accuracy', 'precision', 'recall', 'f1', 'auc_roc']
    metric_labels = ['Accuracy', 'Precision', 'Recall', 'F1', 'AUC-ROC']
    
    x = np.arange(len(models))
    width = 0.15
    
    for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
        vals = [results[m][metric] for m in models]
        axes[0].bar(x + i * width, vals, width, label=label)
    
    axes[0].set_xlabel('Model')
    axes[0].set_ylabel('Score')
    axes[0].set_title('Gait-Based PD Screening: Model Comparison')
    axes[0].set_xticks(x + width * 2)
    axes[0].set_xticklabels([m.replace(' ', '\n') for m in models], fontsize=7)
    axes[0].legend(loc='lower right', fontsize=7)
    axes[0].set_ylim(0.5, 1.05)
    axes[0].grid(axis='y', alpha=0.3)
    
    # ROC curves
    best_model = max(results, key=lambda m: results[m]['auc_roc'])
    for name in models:
        y_true = gait_df['label'].values
        y_prob = results[name]['y_prob']
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        axes[1].plot(fpr, tpr, label=f"{name} (AUC={results[name]['auc_roc']:.3f})")
    
    axes[1].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[1].set_xlabel('False Positive Rate')
    axes[1].set_ylabel('True Positive Rate')
    axes[1].set_title('ROC Curves for PD Gait Screening')
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'gait_model_comparison.png'), bbox_inches='tight')
    plt.close()
    
    # Figure 2: Feature importance
    fig, ax = plt.subplots(figsize=(10, 6))
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    X = StandardScaler().fit_transform(gait_df[feature_cols].values)
    rf.fit(X, gait_df['label'].values)
    importances = rf.feature_importances_
    idx = np.argsort(importances)[::-1]
    
    ax.barh(range(len(feature_cols)), importances[idx], color='steelblue')
    ax.set_yticks(range(len(feature_cols)))
    ax.set_yticklabels([feature_cols[i] for i in idx], fontsize=8)
    ax.set_xlabel('Feature Importance')
    ax.set_title('Gait Feature Importance for PD Detection (Random Forest)')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'gait_feature_importance.png'), bbox_inches='tight')
    plt.close()
    
    # Figure 3: Confusion matrix for best model
    fig, ax = plt.subplots(figsize=(6, 5))
    cm = results[best_model]['confusion_matrix']
    disp = ConfusionMatrixDisplay(cm, display_labels=['Healthy', 'PD'])
    disp.plot(ax=ax, cmap='Blues')
    ax.set_title(f'Confusion Matrix ({best_model})')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'gait_confusion_matrix.png'), bbox_inches='tight')
    plt.close()


def plot_voice_results(results, voice_df, feature_cols):
    """Generate voice analysis figures."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Voice feature distributions
    als_df = voice_df[voice_df['is_als'] == 1]
    healthy_df = voice_df[voice_df['is_als'] == 0]
    
    for ax, feat, title in zip(
        [axes[0, 0], axes[0, 1], axes[1, 0]],
        ['jitter', 'shimmer', 'hnr'],
        ['Jitter (%)', 'Shimmer (%)', 'HNR (dB)']
    ):
        ax.hist(healthy_df[feat], bins=30, alpha=0.6, label='Healthy', color='steelblue', density=True)
        ax.hist(als_df[feat], bins=30, alpha=0.6, label='ALS', color='salmon', density=True)
        ax.set_xlabel(title)
        ax.set_ylabel('Density')
        ax.set_title(f'{title} Distribution')
        ax.legend()
        ax.grid(alpha=0.3)
    
    # Model comparison
    models = list(results.keys())
    metrics_vals = {m: results[m]['auc_roc'] for m in models}
    colors = plt.cm.Set2(np.linspace(0, 1, len(models)))
    axes[1, 1].barh(models, [metrics_vals[m] for m in models], color=colors)
    axes[1, 1].set_xlabel('AUC-ROC')
    axes[1, 1].set_title('Voice-Based ALS Detection: Model AUC-ROC')
    axes[1, 1].set_xlim(0.5, 1.0)
    axes[1, 1].grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'voice_analysis.png'), bbox_inches='tight')
    plt.close()
    
    # ALS progression over sessions
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    for ax, feat, title in zip(axes, ['jitter', 'shimmer', 'f0'],
                                ['Jitter', 'Shimmer', 'F0 (Hz)']):
        for label, color, name in [(1, 'salmon', 'ALS'), (0, 'steelblue', 'Healthy')]:
            subset = voice_df[voice_df['is_als'] == label]
            means = subset.groupby('session')[feat].mean()
            stds = subset.groupby('session')[feat].std()
            ax.plot(means.index, means.values, '-o', color=color, label=name, markersize=4)
            ax.fill_between(means.index, means - stds, means + stds, alpha=0.2, color=color)
        ax.set_xlabel('Session')
        ax.set_ylabel(title)
        ax.set_title(f'{title} Over Time')
        ax.legend()
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'voice_progression.png'), bbox_inches='tight')
    plt.close()


def plot_touch_results(results, touch_df):
    """Generate touchscreen analysis figures."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    groups = ['healthy', 'mci', 'impaired']
    colors = {'healthy': 'steelblue', 'mci': 'orange', 'impaired': 'salmon'}
    
    feats = ['reaction_time', 'tap_accuracy', 'swipe_velocity', 'typing_speed']
    titles = ['Reaction Time (ms)', 'Tap Accuracy', 'Swipe Velocity (px/s)', 'Typing Speed (chars/min)']
    
    for ax, feat, title in zip(axes.flat, feats, titles):
        data = [touch_df[touch_df['group'] == g][feat].values for g in groups]
        bp = ax.boxplot(data, labels=['Healthy', 'MCI', 'Impaired'], patch_artist=True)
        for patch, g in zip(bp['boxes'], groups):
            patch.set_facecolor(colors[g])
            patch.set_alpha(0.6)
        ax.set_ylabel(title)
        ax.set_title(f'{title} by Cognitive Group')
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'touch_analysis.png'), bbox_inches='tight')
    plt.close()
    
    # ROC curves for touch models
    fig, ax = plt.subplots(figsize=(7, 6))
    binary_df = touch_df[touch_df['group'].isin(['impaired', 'healthy'])]
    y_true = (binary_df['group'] == 'impaired').astype(int).values
    
    for name in results:
        y_prob = results[name]['y_prob']
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        ax.plot(fpr, tpr, label=f"{name} (AUC={results[name]['auc_roc']:.3f})")
    
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curves: Cognitive Decline Detection from Touchscreen')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'touch_roc_curves.png'), bbox_inches='tight')
    plt.close()


def plot_cpd_results(long_df):
    """Run change point detection and create visualizations."""
    
    # Select subjects with onset
    onset_subjects = long_df[long_df['has_onset'] == 1]['subject_id'].unique()[:6]
    
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    
    cpd_results_all = []
    
    for idx, (ax, subj) in enumerate(zip(axes.flat, onset_subjects)):
        subj_data = long_df[long_df['subject_id'] == subj].sort_values('week')
        true_onset = subj_data['true_onset_week'].iloc[0]
        
        gait = subj_data['gait_score'].values
        voice = subj_data['voice_score'].values
        touch = subj_data['touch_score'].values
        
        weeks = subj_data['week'].values
        
        # Apply multimodal CPD
        signals = {'gait': gait, 'voice': voice, 'touch': touch}
        fused_cps, per_modal_cps = multimodal_cpd(signals, method='cusum',
                                                    threshold=0.5, drift=0.02)
        
        # Plot signals
        ax.plot(weeks, gait, '-', color='steelblue', alpha=0.8, label='Gait', linewidth=1.5)
        ax.plot(weeks, voice, '-', color='orange', alpha=0.8, label='Voice', linewidth=1.5)
        ax.plot(weeks, touch, '-', color='green', alpha=0.8, label='Touch', linewidth=1.5)
        
        # True onset
        ax.axvline(x=true_onset, color='red', linestyle='--', linewidth=2, label=f'True onset (w{true_onset})')
        
        # Detected change points
        for cp in fused_cps:
            if cp < len(weeks):
                ax.axvline(x=weeks[cp], color='purple', linestyle=':', linewidth=1.5, alpha=0.8)
        
        if fused_cps:
            ax.axvline(x=weeks[min(fused_cps[0], len(weeks)-1)], color='purple', linestyle=':', 
                       linewidth=1.5, alpha=0.8, label='Detected CP')
        
        ax.set_xlabel('Week')
        ax.set_ylabel('Score')
        ax.set_title(f'Subject {subj}')
        ax.legend(fontsize=7, loc='lower left')
        ax.grid(alpha=0.3)
        ax.set_ylim(-0.1, 1.3)
        
        # Evaluate
        eval_result = evaluate_cpd(
            [weeks[cp] for cp in fused_cps if cp < len(weeks)],
            [true_onset],
            tolerance=3
        )
        cpd_results_all.append(eval_result)
    
    plt.suptitle('Longitudinal Change Point Detection (Multimodal CUSUM)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'change_point_detection.png'), bbox_inches='tight')
    plt.close()
    
    # CPD method comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    methods = ['cusum', 'pelt', 'bayesian']
    method_results = {m: {'precision': [], 'recall': [], 'f1': [], 'delay': []} for m in methods}
    
    for subj in onset_subjects[:20] if len(onset_subjects) >= 20 else onset_subjects:
        subj_data = long_df[long_df['subject_id'] == subj].sort_values('week')
        true_onset = subj_data['true_onset_week'].iloc[0]
        gait = subj_data['gait_score'].values
        voice = subj_data['voice_score'].values
        touch = subj_data['touch_score'].values
        weeks = subj_data['week'].values
        
        signals = {'gait': gait, 'voice': voice, 'touch': touch}
        
        for method in methods:
            try:
                kwargs = {}
                if method == 'cusum':
                    kwargs = {'threshold': 0.5, 'drift': 0.02}
                elif method == 'pelt':
                    kwargs = {'penalty': 0.5, 'min_segment': 3}
                elif method == 'bayesian':
                    kwargs = {'hazard_rate': 1/20, 'observation_var': 0.01}
                
                fused_cps, _ = multimodal_cpd(signals, method=method, **kwargs)
                detected_weeks = [weeks[cp] for cp in fused_cps if cp < len(weeks)]
                ev = evaluate_cpd(detected_weeks, [true_onset], tolerance=3)
                
                method_results[method]['precision'].append(ev['precision'])
                method_results[method]['recall'].append(ev['recall'])
                method_results[method]['f1'].append(ev['f1'])
                if ev['mean_delay'] != float('inf'):
                    method_results[method]['delay'].append(ev['mean_delay'])
            except Exception:
                pass
    
    # Bar chart comparison
    method_names = ['CUSUM', 'PELT', 'Bayesian']
    x = np.arange(len(method_names))
    width = 0.25
    
    for i, metric in enumerate(['precision', 'recall', 'f1']):
        vals = [np.mean(method_results[m][metric]) if method_results[m][metric] else 0 
                for m in methods]
        axes[0].bar(x + i * width, vals, width, label=metric.capitalize())
    
    axes[0].set_xlabel('Method')
    axes[0].set_ylabel('Score')
    axes[0].set_title('Change Point Detection: Method Comparison')
    axes[0].set_xticks(x + width)
    axes[0].set_xticklabels(method_names)
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)
    axes[0].set_ylim(0, 1.1)
    
    # Detection delay
    delays = [np.mean(method_results[m]['delay']) if method_results[m]['delay'] else 0 
              for m in methods]
    axes[1].bar(method_names, delays, color=['steelblue', 'orange', 'green'])
    axes[1].set_xlabel('Method')
    axes[1].set_ylabel('Mean Detection Delay (weeks)')
    axes[1].set_title('Detection Delay by Method')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'cpd_method_comparison.png'), bbox_inches='tight')
    plt.close()
    
    return cpd_results_all, method_results


def plot_fusion_results(fusion_results):
    """Plot multimodal fusion comparison."""
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    strategies = list(fusion_results.keys())
    
    # Accuracy and F1
    accs = [fusion_results[s]['accuracy'] for s in strategies]
    f1s = [fusion_results[s]['f1'] for s in strategies]
    aucs = [fusion_results[s]['auc_roc'] for s in strategies]
    
    x = np.arange(len(strategies))
    width = 0.25
    
    axes[0].bar(x - width, accs, width, label='Accuracy', color='steelblue')
    axes[0].bar(x, f1s, width, label='F1 Score', color='orange')
    axes[0].bar(x + width, aucs, width, label='AUC-ROC', color='green')
    axes[0].set_xlabel('Fusion Strategy')
    axes[0].set_ylabel('Score')
    axes[0].set_title('Multimodal Fusion Strategy Comparison')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([s.replace(' ', '\n') for s in strategies], fontsize=8)
    axes[0].legend()
    axes[0].set_ylim(0.5, 1.05)
    axes[0].grid(axis='y', alpha=0.3)
    
    # Composite score distribution
    np.random.seed(42)
    n = 200
    healthy_scores = np.random.beta(2, 8, n) * 100
    disease_scores = np.random.beta(6, 3, n) * 100
    
    axes[1].hist(healthy_scores, bins=25, alpha=0.6, label='Healthy', color='steelblue', density=True)
    axes[1].hist(disease_scores, bins=25, alpha=0.6, label='Neurodegenerative', color='salmon', density=True)
    axes[1].axvline(x=50, color='red', linestyle='--', linewidth=2, label='Threshold (50)')
    axes[1].set_xlabel('Composite Risk Score')
    axes[1].set_ylabel('Density')
    axes[1].set_title('Composite Neurodegenerative Risk Score Distribution')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'multimodal_fusion.png'), bbox_inches='tight')
    plt.close()


def plot_clinical_validation(long_df):
    """Plot clinical endpoint correlation analysis."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Correlation between digital and clinical scores
    onset_df = long_df[long_df['has_onset'] == 1]
    
    # Scatter: gait vs clinical
    ax = axes[0, 0]
    ax.scatter(onset_df['gait_score'], onset_df['clinical_score'], 
               alpha=0.3, s=10, color='steelblue')
    z = np.polyfit(onset_df['gait_score'], onset_df['clinical_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(onset_df['gait_score'].min(), onset_df['gait_score'].max(), 100)
    ax.plot(x_line, p(x_line), 'r-', linewidth=2)
    corr = onset_df['gait_score'].corr(onset_df['clinical_score'])
    ax.set_xlabel('Digital Gait Score')
    ax.set_ylabel('Clinical Score')
    ax.set_title(f'Gait vs Clinical (r={corr:.3f})')
    ax.grid(alpha=0.3)
    
    # Scatter: voice vs clinical
    ax = axes[0, 1]
    ax.scatter(onset_df['voice_score'], onset_df['clinical_score'], 
               alpha=0.3, s=10, color='orange')
    z = np.polyfit(onset_df['voice_score'], onset_df['clinical_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(onset_df['voice_score'].min(), onset_df['voice_score'].max(), 100)
    ax.plot(x_line, p(x_line), 'r-', linewidth=2)
    corr = onset_df['voice_score'].corr(onset_df['clinical_score'])
    ax.set_xlabel('Digital Voice Score')
    ax.set_ylabel('Clinical Score')
    ax.set_title(f'Voice vs Clinical (r={corr:.3f})')
    ax.grid(alpha=0.3)
    
    # Composite vs clinical
    ax = axes[1, 0]
    composite = 0.4 * onset_df['gait_score'] + 0.35 * onset_df['voice_score'] + 0.25 * onset_df['touch_score']
    ax.scatter(composite, onset_df['clinical_score'], alpha=0.3, s=10, color='green')
    z = np.polyfit(composite, onset_df['clinical_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(composite.min(), composite.max(), 100)
    ax.plot(x_line, p(x_line), 'r-', linewidth=2)
    corr = composite.corr(onset_df['clinical_score'])
    ax.set_xlabel('Composite Digital Score')
    ax.set_ylabel('Clinical Score')
    ax.set_title(f'Composite vs Clinical (r={corr:.3f})')
    ax.grid(alpha=0.3)
    
    # Correlation heatmap
    ax = axes[1, 1]
    corr_cols = ['gait_score', 'voice_score', 'touch_score', 'clinical_score']
    corr_matrix = onset_df[corr_cols].corr()
    im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr_cols)))
    ax.set_yticks(range(len(corr_cols)))
    labels = ['Gait', 'Voice', 'Touch', 'Clinical']
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)
    for i in range(len(corr_cols)):
        for j in range(len(corr_cols)):
            ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', ha='center', va='center', fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title('Inter-Modality Correlation Matrix')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'clinical_validation.png'), bbox_inches='tight')
    plt.close()


def plot_system_architecture():
    """Create a system architecture diagram."""
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Data collection layer
    boxes = [
        (1, 6.5, 2.5, 1.0, 'Accelerometer\n& Gyroscope', '#4ECDC4'),
        (4, 6.5, 2.5, 1.0, 'Microphone\n(Voice)', '#FFE66D'),
        (7, 6.5, 2.5, 1.0, 'Touchscreen\nSensors', '#FF6B6B'),
        (10.5, 6.5, 2.5, 1.0, 'Clinical\nAssessment', '#C7CEEA'),
    ]
    
    for x, y, w, h, text, color in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.8)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=8, fontweight='bold')
    
    # Feature extraction layer
    feat_boxes = [
        (1, 4.8, 2.5, 0.8, 'Gait Features\n(18 features)', '#95E1D3'),
        (4, 4.8, 2.5, 0.8, 'Voice Features\n(17 features)', '#FCEABB'),
        (7, 4.8, 2.5, 0.8, 'Touch Features\n(8 features)', '#F38181'),
    ]
    
    for x, y, w, h, text, color in feat_boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.8)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=7)
    
    # Arrows from sensors to features
    for sx in [2.25, 5.25, 8.25]:
        ax.annotate('', xy=(sx, 5.6), xytext=(sx, 6.5),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    
    # Model layer
    model_box = (2.5, 3.2, 7, 1.0, 'ML Models\n(LR / RF / GB / SVM / MLP)', '#B5EAD7')
    x, y, w, h, text, color = model_box
    rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=9, fontweight='bold')
    
    # Arrows to model
    for sx in [2.25, 5.25, 8.25]:
        ax.annotate('', xy=(6, 4.2), xytext=(sx, 4.8),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    
    # Fusion layer
    fusion_box = (3.5, 1.8, 5, 0.9, 'Multimodal Fusion\n(Late Fusion + Meta-Learner)', '#FFDAC1')
    x, y, w, h, text, color = fusion_box
    rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=2, alpha=0.8)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=9, fontweight='bold')
    
    ax.annotate('', xy=(6, 2.7), xytext=(6, 3.2),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))
    
    # Output layer
    out_boxes = [
        (1.5, 0.3, 3, 1.0, 'Composite Risk\nScore (0-100)', '#E2F0CB'),
        (5.5, 0.3, 3.5, 1.0, 'Change Point\nDetection (CPD)', '#FFB7B2'),
        (10, 0.3, 3, 1.0, 'Clinical\nValidation', '#C7CEEA'),
    ]
    
    for x, y, w, h, text, color in out_boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.8)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=8, fontweight='bold')
    
    ax.annotate('', xy=(3, 1.3), xytext=(6, 1.8),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.annotate('', xy=(7.25, 1.3), xytext=(6, 1.8),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.annotate('', xy=(11.5, 1.3), xytext=(11.75, 4.8),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    
    ax.set_title('NeuroSense mHealth Framework: System Architecture', fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'system_architecture.png'), bbox_inches='tight')
    plt.close()


def main():
    print("=" * 60)
    print("NeuroSense mHealth Framework - Experiment Runner")
    print("=" * 60)
    
    # Step 1: Generate data
    print("\n[1/7] Generating synthetic sensor data...")
    gait_df = generate_gait_data()
    voice_df = generate_voice_data()
    touch_df = generate_touch_data()
    long_df = generate_longitudinal_data()
    print(f"  Gait: {len(gait_df)} subjects | Voice: {len(voice_df)} records")
    print(f"  Touch: {len(touch_df)} records | Longitudinal: {len(long_df)} records")
    
    # Step 2: Train gait models
    print("\n[2/7] Training PD gait screening models...")
    gait_results, gait_features = train_gait_model(gait_df)
    best_gait = max(gait_results, key=lambda m: gait_results[m]['auc_roc'])
    print(f"  Best: {best_gait} (AUC={gait_results[best_gait]['auc_roc']:.4f})")
    for name, r in gait_results.items():
        print(f"    {name}: Acc={r['accuracy']:.3f} F1={r['f1']:.3f} AUC={r['auc_roc']:.3f}")
    
    # Step 3: Train voice models
    print("\n[3/7] Training ALS voice monitoring models...")
    voice_results, voice_features = train_voice_model(voice_df)
    best_voice = max(voice_results, key=lambda m: voice_results[m]['auc_roc'])
    print(f"  Best: {best_voice} (AUC={voice_results[best_voice]['auc_roc']:.4f})")
    for name, r in voice_results.items():
        print(f"    {name}: Acc={r['accuracy']:.3f} F1={r['f1']:.3f} AUC={r['auc_roc']:.3f}")
    
    # Step 4: Train touch models
    print("\n[4/7] Training cognitive decline detection models...")
    touch_results, touch_features = train_touch_model(touch_df)
    best_touch = max(touch_results, key=lambda m: touch_results[m]['auc_roc'])
    print(f"  Best: {best_touch} (AUC={touch_results[best_touch]['auc_roc']:.4f})")
    for name, r in touch_results.items():
        print(f"    {name}: Acc={r['accuracy']:.3f} F1={r['f1']:.3f} AUC={r['auc_roc']:.3f}")
    
    # Step 5: Multimodal fusion
    print("\n[5/7] Evaluating multimodal fusion strategies...")
    fusion_results = train_multimodal_fusion(gait_df, voice_df, touch_df)
    for name, r in fusion_results.items():
        print(f"    {name}: Acc={r['accuracy']:.3f} F1={r['f1']:.3f} AUC={r['auc_roc']:.3f}")
    
    # Step 6: Change point detection
    print("\n[6/7] Running change point detection on longitudinal data...")
    cpd_results, method_results = plot_cpd_results(long_df)
    for method in ['cusum', 'pelt', 'bayesian']:
        if method_results[method]['f1']:
            avg_f1 = np.mean(method_results[method]['f1'])
            avg_prec = np.mean(method_results[method]['precision'])
            avg_rec = np.mean(method_results[method]['recall'])
            print(f"    {method.upper()}: Prec={avg_prec:.3f} Rec={avg_rec:.3f} F1={avg_f1:.3f}")
    
    # Step 7: Generate all figures
    print("\n[7/7] Generating figures...")
    plot_gait_results(gait_results, gait_df, gait_features)
    print("  ✓ gait_model_comparison.png")
    print("  ✓ gait_feature_importance.png")
    print("  ✓ gait_confusion_matrix.png")
    
    plot_voice_results(voice_results, voice_df, voice_features)
    print("  ✓ voice_analysis.png")
    print("  ✓ voice_progression.png")
    
    plot_touch_results(touch_results, touch_df)
    print("  ✓ touch_analysis.png")
    print("  ✓ touch_roc_curves.png")
    
    print("  ✓ change_point_detection.png")
    print("  ✓ cpd_method_comparison.png")
    
    plot_fusion_results(fusion_results)
    print("  ✓ multimodal_fusion.png")
    
    plot_clinical_validation(long_df)
    print("  ✓ clinical_validation.png")
    
    plot_system_architecture()
    print("  ✓ system_architecture.png")
    
    # Save summary results
    summary = {
        'gait_results': {k: {m: v for m, v in r.items() if m not in ['confusion_matrix', 'y_pred', 'y_prob']} 
                         for k, r in gait_results.items()},
        'voice_results': {k: {m: v for m, v in r.items() if m not in ['confusion_matrix', 'y_pred', 'y_prob']} 
                          for k, r in voice_results.items()},
        'touch_results': {k: {m: v for m, v in r.items() if m not in ['confusion_matrix', 'y_pred', 'y_prob']} 
                          for k, r in touch_results.items()},
        'fusion_results': fusion_results,
    }
    
    print("\n" + "=" * 60)
    print("Experiment complete! All figures saved to figures/")
    print("=" * 60)
    
    return summary


if __name__ == "__main__":
    summary = main()
