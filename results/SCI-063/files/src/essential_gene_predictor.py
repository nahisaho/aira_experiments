#!/usr/bin/env python3
"""
Module 1: Essential Gene Prediction
Machine learning + transposon mutagenesis data integration
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import (roc_auc_score, roc_curve, confusion_matrix)
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import mutual_info_classif
import json
import os

class EssentialGenePredictor:
    """Predict essential genes using ML + transposon insertion data."""

    FEATURE_NAMES = [
        'tn_insertion_density', 'phyletic_retention', 'codon_adaptation_index',
        'protein_length', 'gc_content', 'network_degree', 'network_betweenness',
        'expression_level', 'operonic_context', 'functional_redundancy',
        'metabolic_flux', 'subcellular_location', 'domain_count',
        'evolutionary_rate', 'gene_strand_bias',
    ]

    def __init__(self, n_genes=525, seed=42):
        self.n_genes = n_genes
        self.seed = seed
        self.rng = np.random.RandomState(seed)
        self.scaler = StandardScaler()

    def generate_synthetic_dataset(self):
        """Generate realistic dataset based on M. genitalium properties. ~382/525 essential (Glass et al. 2006)."""
        n = self.n_genes
        rng = self.rng
        n_essential = 382
        labels = np.zeros(n, dtype=int)
        labels[rng.choice(n, n_essential, replace=False)] = 1
        features = np.zeros((n, len(self.FEATURE_NAMES)))

        for i in range(n):
            e = labels[i]
            features[i, 0] = rng.exponential(0.02 if e else 0.15)
            features[i, 1] = rng.beta(8, 2) if e else rng.beta(3, 5)
            features[i, 2] = np.clip(rng.normal(0.72 if e else 0.55, 0.08 if e else 0.12), 0, 1)
            features[i, 3] = rng.lognormal(5.5 if e else 5.2, 0.7 if e else 0.9)
            features[i, 4] = np.clip(rng.normal(0.317, 0.04), 0.15, 0.50)
            features[i, 5] = rng.poisson(8 if e else 3)
            features[i, 6] = rng.exponential(0.05 if e else 0.02)
            features[i, 7] = rng.lognormal(3.0 if e else 2.0, 1.0 if e else 1.2)
            features[i, 8] = rng.poisson(2.5 if e else 1.5)
            features[i, 9] = rng.poisson(0.3 if e else 1.5)
            features[i, 10] = rng.exponential(5.0 if e else 2.0)
            features[i, 11] = rng.choice([0, 1, 2], p=[0.6, 0.3, 0.1])
            features[i, 12] = rng.poisson(2.0 if e else 1.2)
            features[i, 13] = rng.exponential(0.08 if e else 0.25)
            features[i, 14] = rng.binomial(1, 0.75 if e else 0.55)

        gene_ids = [f"MG_{i+1:04d}" for i in range(n)]
        df = pd.DataFrame(features, columns=self.FEATURE_NAMES, index=gene_ids)
        df['essential'] = labels
        df.index.name = 'gene_id'
        self.dataset = df
        return df

    def train_ensemble(self):
        X = self.scaler.fit_transform(self.dataset[self.FEATURE_NAMES].values)
        y = self.dataset['essential'].values
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=self.seed)

        self.model_rf = RandomForestClassifier(n_estimators=500, max_depth=10, min_samples_leaf=5,
                                                class_weight='balanced', random_state=self.seed, n_jobs=-1)
        self.model_gb = GradientBoostingClassifier(n_estimators=300, max_depth=5, learning_rate=0.05,
                                                    subsample=0.8, random_state=self.seed)

        rf_scores = cross_val_score(self.model_rf, X, y, cv=cv, scoring='roc_auc')
        gb_scores = cross_val_score(self.model_gb, X, y, cv=cv, scoring='roc_auc')

        self.model_rf.fit(X, y)
        self.model_gb.fit(X, y)

        rf_proba = self.model_rf.predict_proba(X)[:, 1]
        gb_proba = self.model_gb.predict_proba(X)[:, 1]
        self.dataset['rf_prob'] = rf_proba
        self.dataset['gb_prob'] = gb_proba
        self.dataset['ensemble_prob'] = 0.5 * rf_proba + 0.5 * gb_proba
        self.dataset['predicted_essential'] = (self.dataset['ensemble_prob'] >= 0.5).astype(int)

        return {
            'rf_auc_cv_mean': float(rf_scores.mean()), 'rf_auc_cv_std': float(rf_scores.std()),
            'gb_auc_cv_mean': float(gb_scores.mean()), 'gb_auc_cv_std': float(gb_scores.std()),
        }

    def feature_importance_analysis(self):
        X = self.dataset[self.FEATURE_NAMES].values
        y = self.dataset['essential'].values
        mi = mutual_info_classif(X, y, random_state=self.seed)
        return pd.DataFrame({
            'feature': self.FEATURE_NAMES,
            'rf_importance': self.model_rf.feature_importances_,
            'mutual_information': mi
        }).sort_values('rf_importance', ascending=False)

    def identify_essential_set(self, threshold=0.5):
        d = self.dataset
        pred = d['predicted_essential']
        tp = int(((pred == 1) & (d['essential'] == 1)).sum())
        fp = int(((pred == 1) & (d['essential'] == 0)).sum())
        fn = int(((pred == 0) & (d['essential'] == 1)).sum())
        tn = int(((pred == 0) & (d['essential'] == 0)).sum())
        n_ess = int(d['essential'].sum())
        summary = {
            'total_genes': len(d), 'predicted_essential': int(pred.sum()),
            'true_essential': n_ess,
            'true_positive': tp, 'false_positive': fp, 'false_negative': fn, 'true_negative': tn,
            'sensitivity': tp / max(n_ess, 1), 'specificity': tn / max(tn + fp, 1),
            'ppv': tp / max(tp + fp, 1),
        }
        return d[d['ensemble_prob'] >= threshold].sort_values('ensemble_prob', ascending=False), summary

    def generate_plots(self, output_dir='figures'):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        os.makedirs(output_dir, exist_ok=True)

        # Feature importance
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        imp = self.feature_importance_analysis()
        sns.barplot(data=imp, y='feature', x='rf_importance', ax=axes[0], palette='viridis')
        axes[0].set_title('Random Forest Feature Importance')
        sns.barplot(data=imp, y='feature', x='mutual_information', ax=axes[1], palette='cividis')
        axes[1].set_title('Mutual Information Score')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig1_feature_importance.png', dpi=300, bbox_inches='tight')
        plt.close()

        # ROC
        fig, ax = plt.subplots(figsize=(8, 8))
        y = self.dataset['essential'].values
        for name, col in [('Random Forest','rf_prob'),('Gradient Boosting','gb_prob'),('Ensemble','ensemble_prob')]:
            fpr, tpr, _ = roc_curve(y, self.dataset[col])
            auc = roc_auc_score(y, self.dataset[col])
            ax.plot(fpr, tpr, label=f'{name} (AUC={auc:.3f})', linewidth=2)
        ax.plot([0,1],[0,1],'k--',alpha=0.5)
        ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curves — Essential Gene Prediction'); ax.legend(); ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig2_roc_curves.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Tn insertion density
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(self.dataset[self.dataset['essential']==1]['tn_insertion_density'], bins=40, alpha=0.6, label='Essential', color='#2196F3', density=True)
        ax.hist(self.dataset[self.dataset['essential']==0]['tn_insertion_density'], bins=40, alpha=0.6, label='Non-essential', color='#FF5722', density=True)
        ax.set_xlabel('Transposon Insertion Density'); ax.set_ylabel('Density')
        ax.set_title('Tn Insertion Density: Essential vs Non-essential'); ax.legend()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig3_tn_insertion_density.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Confusion matrix
        cm = confusion_matrix(y, self.dataset['predicted_essential'].values)
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['Non-ess','Essential'], yticklabels=['Non-ess','Essential'], ax=ax)
        ax.set_xlabel('Predicted'); ax.set_ylabel('Actual'); ax.set_title('Confusion Matrix — Ensemble')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig4_confusion_matrix.png', dpi=300, bbox_inches='tight')
        plt.close()

        imp.to_csv('results/feature_importance.csv', index=False)
        return [f'{output_dir}/fig{i}_{n}.png' for i, n in [(1,'feature_importance'),(2,'roc_curves'),(3,'tn_insertion_density'),(4,'confusion_matrix')]]

def run_module1():
    print("="*60); print("MODULE 1: Essential Gene Prediction"); print("="*60)
    p = EssentialGenePredictor(525, 42)
    df = p.generate_synthetic_dataset(); df.to_csv('data/essential_genes_dataset.csv')
    print(f"  Dataset: {len(df)} genes, {int(df['essential'].sum())} essential")
    cv = p.train_ensemble()
    print(f"  RF AUC: {cv['rf_auc_cv_mean']:.4f}±{cv['rf_auc_cv_std']:.4f}")
    print(f"  GB AUC: {cv['gb_auc_cv_mean']:.4f}±{cv['gb_auc_cv_std']:.4f}")
    ess, s = p.identify_essential_set()
    ess.to_csv('results/predicted_essential_genes.csv')
    with open('results/module1_summary.json','w') as f: json.dump({**cv,**s}, f, indent=2)
    print(f"  Predicted essential: {s['predicted_essential']}/{s['total_genes']}")
    print(f"  Sensitivity={s['sensitivity']:.4f} Specificity={s['specificity']:.4f} PPV={s['ppv']:.4f}")
    plots = p.generate_plots()
    for pl in plots: print(f"  Saved: {pl}")
    return p, s

if __name__ == '__main__': run_module1()
