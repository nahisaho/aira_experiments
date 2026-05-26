#!/usr/bin/env python3
"""
Module 5: RA drug response prediction model.
Uses multi-omics features to predict treatment response via ML.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (roc_curve, auc, precision_recall_curve,
                             classification_report, confusion_matrix, f1_score)
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
import os

np.random.seed(42)

def generate_drug_response_data(n_patients=200):
    """Generate multi-omics features and drug response labels."""
    n_features = 80  # Combined omics features
    
    X = np.random.normal(0, 1, (n_patients, n_features))
    
    # Create response-associated signal in first 20 features
    # Responders: ~60%, Non-responders: ~40%
    response_weights = np.zeros(n_features)
    response_weights[:20] = np.random.uniform(0.3, 0.8, 20)
    
    logits = X @ response_weights + np.random.normal(0, 1, n_patients)
    prob = 1 / (1 + np.exp(-logits))
    y = (prob > 0.5).astype(int)  # 1 = responder, 0 = non-responder
    
    feature_names = (
        [f'Gene_{i}' for i in range(30)] +
        [f'Protein_{i}' for i in range(20)] +
        [f'Metabolite_{i}' for i in range(15)] +
        [f'CellFrac_{i}' for i in range(10)] +
        [f'Clinical_{i}' for i in range(5)]
    )
    # Replace some with meaningful names
    meaningful = ['TNF_expr', 'IL6_expr', 'IL17A_expr', 'IL10_expr', 'FOXP3_expr',
                  'DAS28_baseline', 'CRP_level', 'RF_titer', 'ACPA_status', 'Th17_frac',
                  'Treg_frac', 'M1_frac', 'PD1_expr', 'CTLA4_expr', 'HAQ_score']
    for i, name in enumerate(meaningful):
        feature_names[i] = name
    
    X_df = pd.DataFrame(X, columns=feature_names)
    return X_df, y, feature_names

def train_and_evaluate(X, y):
    """Train multiple classifiers and evaluate via cross-validation."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Feature selection
    selector = SelectKBest(f_classif, k=30)
    X_selected = selector.fit_transform(X_scaled, y)
    selected_mask = selector.get_support()
    
    models = {
        'Random Forest': RandomForestClassifier(n_estimators=200, max_depth=6, random_state=42),
        'Gradient Boosting': GradientBoostingClassifier(n_estimators=150, max_depth=4, random_state=42),
        'Logistic Regression': LogisticRegression(C=1.0, max_iter=1000, random_state=42),
        'SVM (RBF)': SVC(kernel='rbf', probability=True, random_state=42)
    }
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results = {}
    
    for name, model in models.items():
        # AUC scores
        auc_scores = cross_val_score(model, X_selected, y, cv=cv, scoring='roc_auc')
        f1_scores = cross_val_score(model, X_selected, y, cv=cv, scoring='f1')
        acc_scores = cross_val_score(model, X_selected, y, cv=cv, scoring='accuracy')
        
        results[name] = {
            'AUC': (auc_scores.mean(), auc_scores.std()),
            'F1': (f1_scores.mean(), f1_scores.std()),
            'Accuracy': (acc_scores.mean(), acc_scores.std()),
            'model': model,
        }
        print(f"  {name}: AUC={auc_scores.mean():.3f}±{auc_scores.std():.3f}, "
              f"F1={f1_scores.mean():.3f}±{f1_scores.std():.3f}")
    
    return results, X_selected, selected_mask, scaler, selector

def plot_roc_curves(X, y, results):
    """Plot ROC curves for all models."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12']
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    for (name, res), color in zip(results.items(), colors):
        model = res['model']
        tprs, aucs = [], []
        mean_fpr = np.linspace(0, 1, 100)
        
        for train_idx, test_idx in cv.split(X, y):
            model.fit(X[train_idx], y[train_idx])
            y_prob = model.predict_proba(X[test_idx])[:, 1]
            fpr, tpr, _ = roc_curve(y[test_idx], y_prob)
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            aucs.append(auc(fpr, tpr))
        
        mean_tpr = np.mean(tprs, axis=0)
        mean_auc = np.mean(aucs)
        std_auc = np.std(aucs)
        
        axes[0].plot(mean_fpr, mean_tpr, color=color, linewidth=2,
                    label=f'{name} (AUC={mean_auc:.3f}±{std_auc:.3f})')
    
    axes[0].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[0].set_xlabel('False Positive Rate')
    axes[0].set_ylabel('True Positive Rate')
    axes[0].set_title('ROC Curves: Drug Response Prediction')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)
    
    # Bar chart comparison
    model_names = list(results.keys())
    metrics = ['AUC', 'F1', 'Accuracy']
    x = np.arange(len(model_names))
    width = 0.25
    
    for i, metric in enumerate(metrics):
        vals = [results[m][metric][0] for m in model_names]
        errs = [results[m][metric][1] for m in model_names]
        axes[1].bar(x + i*width, vals, width, yerr=errs, label=metric, alpha=0.8)
    
    axes[1].set_xticks(x + width)
    axes[1].set_xticklabels(model_names, rotation=15, fontsize=9)
    axes[1].set_ylabel('Score')
    axes[1].set_title('Model Performance Comparison')
    axes[1].legend()
    axes[1].set_ylim(0.5, 1.0)
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/drug_response_roc.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_feature_importance(X_df, y, selected_mask, feature_names):
    """Plot feature importance from Random Forest."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_df)
    selector = SelectKBest(f_classif, k=30)
    X_sel = selector.fit_transform(X_scaled, y)
    
    rf = RandomForestClassifier(n_estimators=200, max_depth=6, random_state=42)
    rf.fit(X_sel, y)
    
    selected_names = [feature_names[i] for i in range(len(feature_names)) if selected_mask[i]]
    importances = rf.feature_importances_
    sorted_idx = np.argsort(importances)[-15:]
    
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(range(len(sorted_idx)), importances[sorted_idx], color='#3498DB', edgecolor='black')
    ax.set_yticks(range(len(sorted_idx)))
    ax.set_yticklabels([selected_names[i] for i in sorted_idx])
    ax.set_xlabel('Feature Importance (Gini)')
    ax.set_title('Top 15 Predictive Features for Drug Response')
    plt.tight_layout()
    plt.savefig('figures/feature_importance.png', dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    
    print("Generating drug response dataset...")
    X_df, y, feature_names = generate_drug_response_data(200)
    print(f"  Responders: {y.sum()}, Non-responders: {(1-y).sum()}")
    
    print("Training models...")
    results, X_selected, selected_mask, scaler, selector = train_and_evaluate(X_df, y)
    
    print("Plotting ROC curves...")
    plot_roc_curves(X_selected, y, results)
    
    print("Plotting feature importance...")
    plot_feature_importance(X_df, y, selected_mask, feature_names)
    
    print("Module 5 complete.")
