"""
Main experiment pipeline for epigenetic clock development and evaluation.
Runs all models, generates figures, and outputs results.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, KFold
from scipy import stats
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from data_generator import generate_dataset
from models import (
    HorvathBaseline, TissueAwareClock, DeepClockTrainer,
    GradientBoostClock, evaluate_model
)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

plt.rcParams.update({
    'font.size': 11,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.dpi': 150,
})


def run_experiment():
    print("=" * 60)
    print("Epigenetic Clock Improvement: Full Experiment Pipeline")
    print("=" * 60)
    
    # --- Step 1: Data Generation ---
    print("\n[1/7] Generating synthetic DNA methylation data...")
    df_meth, df_meta = generate_dataset()
    X = df_meth.values
    y_chrono = df_meta['chronological_age'].values
    y_bio = df_meta['biological_age'].values
    tissues = df_meta['tissue'].values
    interventions = df_meta['intervention'].values
    
    print(f"  Samples: {X.shape[0]}, CpG sites: {X.shape[1]}")
    print(f"  Tissues: {np.unique(tissues)}")
    print(f"  Age range: {y_chrono.min():.1f} - {y_chrono.max():.1f}")
    
    # Train/test split
    idx = np.arange(len(X))
    train_idx, test_idx = train_test_split(idx, test_size=0.2, random_state=42)
    
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y_bio[train_idx], y_bio[test_idx]
    y_chrono_train, y_chrono_test = y_chrono[train_idx], y_chrono[test_idx]
    tissues_train, tissues_test = tissues[train_idx], tissues[test_idx]
    interventions_test = interventions[test_idx]
    
    results = {}
    
    # --- Step 2: Horvath Baseline (ElasticNet) ---
    print("\n[2/7] Training Horvath-like baseline (ElasticNet)...")
    horvath = HorvathBaseline(alpha=0.1, l1_ratio=0.5)
    horvath.fit(X_train, y_train)
    pred_horvath = horvath.predict(X_test)
    results['ElasticNet (Horvath-like)'] = evaluate_model(y_test, pred_horvath)
    n_selected = len(horvath.get_selected_cpgs())
    print(f"  Selected CpGs: {n_selected}/{X.shape[1]}")
    print(f"  MAE: {results['ElasticNet (Horvath-like)']['MAE']:.2f} years")
    
    # --- Step 3: Gradient Boosting Clock ---
    print("\n[3/7] Training Gradient Boosting Clock...")
    gb_clock = GradientBoostClock()
    gb_clock.fit(X_train, y_train)
    pred_gb = gb_clock.predict(X_test)
    results['Gradient Boosting'] = evaluate_model(y_test, pred_gb)
    print(f"  MAE: {results['Gradient Boosting']['MAE']:.2f} years")
    
    # --- Step 4: Tissue-Aware Clock ---
    print("\n[4/7] Training Tissue-Aware Clock...")
    tissue_clock = TissueAwareClock()
    tissue_clock.fit(X_train, y_train, tissues_train)
    pred_tissue = tissue_clock.predict(X_test, tissues_test)
    results['Tissue-Aware ElasticNet'] = evaluate_model(y_test, pred_tissue)
    print(f"  MAE: {results['Tissue-Aware ElasticNet']['MAE']:.2f} years")
    
    # --- Step 5: Deep Learning Clock ---
    print("\n[5/7] Training DeepEpiClock (Neural Network)...")
    # Further split train into train/val
    train_sub_idx, val_idx = train_test_split(
        np.arange(len(X_train)), test_size=0.15, random_state=42
    )
    
    deep_clock = DeepClockTrainer(
        n_features=X.shape[1], n_tissues=5, lr=0.002, epochs=150, batch_size=32
    )
    deep_clock.fit(
        X_train[train_sub_idx], y_train[train_sub_idx], tissues_train[train_sub_idx],
        X_train[val_idx], y_train[val_idx], tissues_train[val_idx]
    )
    pred_deep = deep_clock.predict(X_test, tissues_test)
    results['DeepEpiClock'] = evaluate_model(y_test, pred_deep)
    print(f"  MAE: {results['DeepEpiClock']['MAE']:.2f} years")
    
    # --- Step 6: Cross-validated results ---
    print("\n[6/7] Running 5-fold cross-validation for DeepEpiClock...")
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_maes = []
    for fold, (cv_train, cv_test) in enumerate(kf.split(X)):
        cv_sub_train, cv_val = train_test_split(cv_train, test_size=0.15, random_state=fold)
        dc = DeepClockTrainer(n_features=X.shape[1], n_tissues=5, lr=0.002, epochs=120, batch_size=32)
        dc.fit(X[cv_sub_train], y_bio[cv_sub_train], tissues[cv_sub_train],
               X[cv_val], y_bio[cv_val], tissues[cv_val])
        cv_pred = dc.predict(X[cv_test], tissues[cv_test])
        fold_mae = np.mean(np.abs(y_bio[cv_test] - cv_pred))
        cv_maes.append(fold_mae)
        print(f"  Fold {fold+1}: MAE = {fold_mae:.2f}")
    print(f"  Mean CV MAE: {np.mean(cv_maes):.2f} ± {np.std(cv_maes):.2f}")
    
    # --- Step 7: Generate Figures ---
    print("\n[7/7] Generating figures...")
    
    predictions_dict = {
        'ElasticNet (Horvath-like)': pred_horvath,
        'Gradient Boosting': pred_gb,
        'Tissue-Aware ElasticNet': pred_tissue,
        'DeepEpiClock': pred_deep,
    }
    
    # Figure 1: Model comparison scatter plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (name, pred) in zip(axes.flatten(), predictions_dict.items()):
        ax.scatter(y_test, pred, alpha=0.4, s=15, c='steelblue')
        lims = [min(y_test.min(), pred.min()) - 2, max(y_test.max(), pred.max()) + 2]
        ax.plot(lims, lims, 'r--', lw=1.5, label='Perfect prediction')
        ax.set_xlabel('True Biological Age')
        ax.set_ylabel('Predicted Age')
        r = results[name]
        ax.set_title(f'{name}\nMAE={r["MAE"]:.2f}, R²={r["R2"]:.3f}')
        ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'model_comparison_scatter.png'))
    plt.close()
    print("  Saved: model_comparison_scatter.png")
    
    # Figure 2: Performance bar chart
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    model_names = list(results.keys())
    metrics_to_plot = ['MAE', 'RMSE', 'R2']
    colors = ['#2196F3', '#FF9800', '#4CAF50', '#E91E63']
    
    for ax, metric in zip(axes, metrics_to_plot):
        vals = [results[m][metric] for m in model_names]
        bars = ax.bar(range(len(model_names)), vals, color=colors)
        ax.set_xticks(range(len(model_names)))
        ax.set_xticklabels([n.replace(' ', '\n') for n in model_names], fontsize=8)
        ax.set_title(metric, fontweight='bold')
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                    f'{v:.3f}', ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'performance_metrics.png'))
    plt.close()
    print("  Saved: performance_metrics.png")
    
    # Figure 3: Training curves for DeepEpiClock
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(deep_clock.train_losses, label='Training Loss', color='steelblue')
    if deep_clock.val_losses:
        ax.plot(deep_clock.val_losses, label='Validation Loss', color='coral')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('MSE Loss')
    ax.set_title('DeepEpiClock Training Curves')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig(os.path.join(FIGURES_DIR, 'training_curves.png'))
    plt.close()
    print("  Saved: training_curves.png")
    
    # Figure 4: Tissue-specific performance
    tissue_results = {}
    for tissue in np.unique(tissues_test):
        mask = tissues_test == tissue
        tissue_results[tissue] = {
            name: np.mean(np.abs(y_test[mask] - pred[mask]))
            for name, pred in predictions_dict.items()
        }
    
    tissue_df = pd.DataFrame(tissue_results).T
    fig, ax = plt.subplots(figsize=(10, 6))
    tissue_df.plot(kind='bar', ax=ax, colormap='Set2')
    ax.set_xlabel('Tissue Type')
    ax.set_ylabel('MAE (years)')
    ax.set_title('Tissue-Specific Model Performance')
    ax.legend(fontsize=8, loc='upper right')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'tissue_specific_performance.png'))
    plt.close()
    print("  Saved: tissue_specific_performance.png")
    
    # Figure 5: Age acceleration distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Age acceleration by intervention
    age_accel_deep = pred_deep - y_chrono_test
    ax = axes[0]
    for intervention in np.unique(interventions_test):
        mask = interventions_test == intervention
        ax.hist(age_accel_deep[mask], bins=20, alpha=0.5, label=intervention, density=True)
    ax.set_xlabel('Age Acceleration (years)')
    ax.set_ylabel('Density')
    ax.set_title('Age Acceleration by Intervention (DeepEpiClock)')
    ax.legend()
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    
    # Age acceleration by tissue
    ax = axes[1]
    data_for_box = []
    labels_for_box = []
    for tissue in np.unique(tissues_test):
        mask = tissues_test == tissue
        data_for_box.append(age_accel_deep[mask])
        labels_for_box.append(tissue)
    ax.boxplot(data_for_box, labels=labels_for_box)
    ax.set_xlabel('Tissue Type')
    ax.set_ylabel('Age Acceleration (years)')
    ax.set_title('Age Acceleration by Tissue (DeepEpiClock)')
    ax.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'age_acceleration_analysis.png'))
    plt.close()
    print("  Saved: age_acceleration_analysis.png")
    
    # Figure 6: Intervention sensitivity analysis
    fig, ax = plt.subplots(figsize=(10, 6))
    sensitivity_data = {}
    for name, pred in predictions_dict.items():
        age_accel = pred - y_chrono_test
        sensitivities = {}
        baseline_accel = np.mean(age_accel[interventions_test == 'none'])
        for intervention in ['exercise', 'diet', 'drug']:
            mask = interventions_test == intervention
            if mask.sum() > 0:
                int_accel = np.mean(age_accel[mask])
                sensitivities[intervention] = baseline_accel - int_accel
        sensitivity_data[name] = sensitivities
    
    sens_df = pd.DataFrame(sensitivity_data).T
    sens_df.plot(kind='bar', ax=ax, colormap='coolwarm')
    ax.set_xlabel('Model')
    ax.set_ylabel('Age Deceleration Effect (years)')
    ax.set_title('Intervention Detection Sensitivity')
    ax.legend(title='Intervention')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'intervention_sensitivity.png'))
    plt.close()
    print("  Saved: intervention_sensitivity.png")
    
    # Figure 7: Residual analysis
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (name, pred) in zip(axes.flatten(), predictions_dict.items()):
        residuals = pred - y_test
        ax.scatter(y_test, residuals, alpha=0.3, s=15, c='steelblue')
        ax.axhline(y=0, color='r', linestyle='--', lw=1.5)
        ax.set_xlabel('True Biological Age')
        ax.set_ylabel('Residual (years)')
        ax.set_title(f'{name}')
        z = np.polyfit(y_test, residuals, 1)
        p = np.poly1d(z)
        x_line = np.linspace(y_test.min(), y_test.max(), 100)
        ax.plot(x_line, p(x_line), 'g-', alpha=0.7, label=f'Trend')
        ax.legend(fontsize=8)
    plt.suptitle('Residual Analysis', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'residual_analysis.png'))
    plt.close()
    print("  Saved: residual_analysis.png")
    
    # Figure 8: Cross-validation boxplot
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(range(len(cv_maes)), cv_maes, color='steelblue', alpha=0.8)
    ax.axhline(y=np.mean(cv_maes), color='r', linestyle='--', label=f'Mean: {np.mean(cv_maes):.2f}')
    ax.set_xlabel('Fold')
    ax.set_ylabel('MAE (years)')
    ax.set_title('DeepEpiClock 5-Fold Cross-Validation')
    ax.set_xticks(range(len(cv_maes)))
    ax.set_xticklabels([f'Fold {i+1}' for i in range(len(cv_maes))])
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'cross_validation.png'))
    plt.close()
    print("  Saved: cross_validation.png")
    
    # --- Summary ---
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    results_df = pd.DataFrame(results).T
    print(results_df.to_string())
    
    # Save results
    results_df.to_csv(os.path.join(os.path.dirname(FIGURES_DIR), 'results_summary.csv'))
    
    # Save detailed results as JSON
    detailed = {
        'model_results': {k: {mk: float(mv) for mk, mv in v.items()} for k, v in results.items()},
        'cv_maes': [float(m) for m in cv_maes],
        'cv_mean_mae': float(np.mean(cv_maes)),
        'cv_std_mae': float(np.std(cv_maes)),
        'tissue_results': {t: {m: float(v) for m, v in tv.items()} for t, tv in tissue_results.items()},
        'sensitivity_data': {k: {sk: float(sv) for sk, sv in v.items()} for k, v in sensitivity_data.items()},
        'n_selected_cpgs': int(n_selected),
        'dataset_info': {
            'n_samples': int(X.shape[0]),
            'n_cpg_sites': int(X.shape[1]),
            'n_train': int(len(train_idx)),
            'n_test': int(len(test_idx)),
        }
    }
    with open(os.path.join(os.path.dirname(FIGURES_DIR), 'results_detailed.json'), 'w') as f:
        json.dump(detailed, f, indent=2)
    
    print("\nAll experiments completed successfully!")
    return detailed


if __name__ == '__main__':
    results = run_experiment()
