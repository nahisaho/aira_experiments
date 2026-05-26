#!/usr/bin/env python3
"""
TCR Repertoire Analysis Pipeline for Immune State Estimation
=============================================================
Implements a comprehensive pipeline covering:
1. TCR-seq data preprocessing (V(D)J annotation, clonotype definition)
2. Repertoire diversity metrics (Shannon entropy, Chao1, Hill numbers)
3. Public TCR identification and HLA restriction prediction
4. TCR-epitope binding prediction (CNN/Transformer-inspired models)
5. Immune age estimation and clonal expansion pattern analysis
6. Cancer immunotherapy biomarker (ICB response prediction)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, roc_curve, precision_recall_curve,
                             average_precision_score, confusion_matrix,
                             classification_report, accuracy_score)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from collections import Counter
import os
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# ============================================================
# Amino acid encoding utilities
# ============================================================
AA_LIST = list('ACDEFGHIKLMNPQRSTVWY')
AA_PROPERTIES = {
    'A': [1.8, 0, 89],   'C': [2.5, 0, 121],  'D': [-3.5, -1, 133],
    'E': [-3.5, -1, 147], 'F': [2.8, 0, 165],  'G': [-0.4, 0, 75],
    'H': [-3.2, 0.5, 155],'I': [4.5, 0, 131],  'K': [-3.9, 1, 146],
    'L': [3.8, 0, 131],   'M': [1.9, 0, 149],  'N': [-3.5, 0, 132],
    'P': [-1.6, 0, 115],  'Q': [-3.5, 0, 146],  'R': [-4.5, 1, 174],
    'S': [-0.8, 0, 105],  'T': [-0.7, 0, 119],  'V': [4.2, 0, 117],
    'W': [-0.9, 0, 204],  'Y': [-1.3, 0, 181]
}

def encode_cdr3(seq, max_len=20):
    """One-hot encode CDR3 sequence with padding."""
    encoded = np.zeros((max_len, 20))
    for i, aa in enumerate(seq[:max_len]):
        if aa in AA_LIST:
            encoded[i, AA_LIST.index(aa)] = 1.0
    return encoded

def physicochemical_features(seq):
    """Extract physicochemical features from CDR3 sequence."""
    props = np.array([AA_PROPERTIES.get(aa, [0, 0, 0]) for aa in seq])
    if len(props) == 0:
        return np.zeros(9)
    return np.concatenate([props.mean(axis=0), props.std(axis=0), [len(seq), props[:,0].sum(), props[:,1].sum()]])

# ============================================================
# 1. TCR-seq Data Generation & Preprocessing
# ============================================================
V_GENES = [f'TRBV{i}-1' for i in range(1, 31)]
D_GENES = [f'TRBD{i}' for i in range(1, 3)]
J_GENES = [f'TRBJ{i}-{j}' for i in range(1, 3) for j in range(1, 8)]

def generate_cdr3(length=None):
    """Generate a random CDR3β sequence."""
    if length is None:
        length = np.random.randint(10, 20)
    return 'C' + ''.join(np.random.choice(AA_LIST) for _ in range(length - 2)) + 'F'

def generate_tcr_repertoire(n_clones, n_samples, sample_labels, age_values=None):
    """Generate synthetic TCR repertoire data with realistic distributions."""
    records = []
    # Generate shared (public) TCRs
    n_public = int(n_clones * 0.05)
    public_cdr3s = [generate_cdr3() for _ in range(n_public)]
    
    for sid, label in enumerate(sample_labels):
        # Power-law clone frequency distribution
        freqs = np.random.pareto(1.5, n_clones) + 1
        
        # Adjust clonal expansion for different conditions
        if label == 'responder':
            top_k = int(n_clones * 0.1)
            freqs[:top_k] *= np.random.uniform(3, 10, top_k)
        elif label == 'tumor':
            top_k = int(n_clones * 0.15)
            freqs[:top_k] *= np.random.uniform(5, 15, top_k)
        elif label == 'aged':
            top_k = int(n_clones * 0.2)
            freqs[:top_k] *= np.random.uniform(4, 12, top_k)
            n_effective = int(n_clones * 0.6)
            freqs[n_effective:] = 0
        
        freqs = freqs / freqs.sum()
        
        for i in range(n_clones):
            if freqs[i] < 1e-6:
                continue
            if i < n_public:
                cdr3 = public_cdr3s[i]
            else:
                cdr3 = generate_cdr3()
            
            records.append({
                'sample_id': f'sample_{sid:03d}',
                'label': label,
                'clone_id': f'clone_{sid:03d}_{i:05d}',
                'cdr3_aa': cdr3,
                'cdr3_length': len(cdr3),
                'v_gene': np.random.choice(V_GENES),
                'd_gene': np.random.choice(D_GENES),
                'j_gene': np.random.choice(J_GENES),
                'frequency': freqs[i],
                'count': max(1, int(freqs[i] * 10000)),
                'age': age_values[sid] if age_values is not None else np.random.randint(20, 80),
            })
    
    return pd.DataFrame(records)

def preprocess_repertoire(df):
    """V(D)J annotation validation and clonotype definition."""
    # Define clonotype as unique CDR3 + V + J combination
    df['clonotype'] = df['cdr3_aa'] + '_' + df['v_gene'] + '_' + df['j_gene']
    # Validate CDR3 sequences
    df['valid_cdr3'] = df['cdr3_aa'].apply(
        lambda x: x.startswith('C') and x.endswith('F') and len(x) >= 8
    )
    df = df[df['valid_cdr3']].copy()
    # Normalize frequencies within each sample
    for sid in df['sample_id'].unique():
        mask = df['sample_id'] == sid
        total = df.loc[mask, 'count'].sum()
        df.loc[mask, 'frequency'] = df.loc[mask, 'count'] / total
    return df

# ============================================================
# 2. Diversity Metrics
# ============================================================
def shannon_entropy(freqs):
    """Shannon entropy H = -sum(p_i * log(p_i))."""
    freqs = freqs[freqs > 0]
    return -np.sum(freqs * np.log2(freqs))

def simpson_index(freqs):
    """Simpson's diversity index D = 1 - sum(p_i^2)."""
    return 1 - np.sum(freqs ** 2)

def chao1_estimator(counts):
    """Chao1 richness estimator."""
    s_obs = np.sum(counts > 0)
    f1 = np.sum(counts == 1)
    f2 = np.sum(counts == 2)
    if f2 > 0:
        return s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    else:
        return s_obs + f1 * (f1 - 1) / 2

def hill_number(freqs, q):
    """Hill number of order q."""
    freqs = freqs[freqs > 0]
    if q == 0:
        return len(freqs)
    elif q == 1:
        return np.exp(shannon_entropy(freqs) * np.log(2))
    else:
        return (np.sum(freqs ** q)) ** (1 / (1 - q))

def inverse_simpson(freqs):
    """Inverse Simpson index (Hill number q=2)."""
    return hill_number(freqs, 2)

def clonality_index(freqs):
    """Clonality = 1 - Pielou's evenness."""
    H = shannon_entropy(freqs)
    S = np.sum(freqs > 0)
    if S <= 1:
        return 1.0
    return 1 - H / np.log2(S)

def compute_diversity_metrics(df):
    """Compute all diversity metrics for each sample."""
    results = []
    for sid in df['sample_id'].unique():
        sample = df[df['sample_id'] == sid]
        freqs = sample['frequency'].values
        counts = sample['count'].values
        label = sample['label'].iloc[0]
        age = sample['age'].iloc[0]
        
        results.append({
            'sample_id': sid,
            'label': label,
            'age': age,
            'n_clonotypes': len(sample),
            'shannon_entropy': shannon_entropy(freqs),
            'simpson_index': simpson_index(freqs),
            'chao1': chao1_estimator(counts),
            'hill_q0': hill_number(freqs, 0),
            'hill_q1': hill_number(freqs, 1),
            'hill_q2': hill_number(freqs, 2),
            'clonality': clonality_index(freqs),
            'top10_freq': np.sort(freqs)[-10:].sum() if len(freqs) >= 10 else freqs.sum(),
            'gini_coefficient': gini(freqs),
        })
    return pd.DataFrame(results)

def gini(freqs):
    """Compute Gini coefficient of clone frequency distribution."""
    sorted_f = np.sort(freqs)
    n = len(sorted_f)
    cum = np.cumsum(sorted_f)
    return (2 * np.sum((np.arange(1, n+1) * sorted_f)) / (n * np.sum(sorted_f))) - (n + 1) / n

# ============================================================
# 3. Public TCR Identification & HLA Restriction Prediction
# ============================================================
HLA_ALLELES = ['HLA-A*02:01', 'HLA-A*01:01', 'HLA-A*03:01', 'HLA-B*07:02',
               'HLA-B*08:01', 'HLA-A*24:02', 'HLA-B*35:01', 'HLA-A*11:01']

def identify_public_tcrs(df, min_samples=2):
    """Identify public TCRs shared across multiple samples."""
    cdr3_samples = df.groupby('cdr3_aa')['sample_id'].nunique()
    public_cdr3s = cdr3_samples[cdr3_samples >= min_samples].index.tolist()
    
    public_df = df[df['cdr3_aa'].isin(public_cdr3s)].copy()
    public_summary = public_df.groupby('cdr3_aa').agg(
        n_samples=('sample_id', 'nunique'),
        mean_freq=('frequency', 'mean'),
        v_genes=('v_gene', lambda x: ','.join(x.unique()[:3])),
    ).reset_index()
    public_summary = public_summary.sort_values('n_samples', ascending=False)
    return public_summary

def predict_hla_restriction(cdr3_seq):
    """
    Simple HLA restriction prediction based on CDR3 motifs.
    In practice, this would use trained models like NetMHCpan.
    """
    # Motif-based heuristic (simplified)
    motif_hla_map = {
        'GS': 'HLA-A*02:01', 'YQ': 'HLA-A*01:01', 'RG': 'HLA-A*03:01',
        'LS': 'HLA-B*07:02', 'SS': 'HLA-B*08:01', 'PY': 'HLA-A*24:02',
        'DT': 'HLA-B*35:01', 'KA': 'HLA-A*11:01',
    }
    scores = {}
    for motif, hla in motif_hla_map.items():
        score = cdr3_seq.count(motif) * 0.3 + np.random.uniform(0, 0.4)
        scores[hla] = score
    best_hla = max(scores, key=scores.get)
    return best_hla, scores[best_hla]

# ============================================================
# 4. TCR-Epitope Binding Prediction (CNN-inspired features)
# ============================================================
KNOWN_EPITOPES = {
    'GILGFVFTL': 'Influenza_M1',   # HLA-A*02:01
    'NLVPMVATV': 'CMV_pp65',       # HLA-A*02:01
    'GLCTLVAML': 'EBV_BMLF1',     # HLA-A*02:01
    'RAKFKQLL':  'EBV_BZLF1',     # HLA-B*08:01
    'LLWNGPMAV': 'HCV_NS5b',      # HLA-A*02:01
    'YLQPRTFLL': 'SARS-CoV-2_S',  # HLA-A*02:01
    'TTDPSFLGRY': 'SARS-CoV-2_N', # HLA-A*03:01
}

def cdr3_epitope_distance(cdr3, epitope):
    """Compute a physicochemical distance between CDR3 and epitope."""
    feat_cdr3 = physicochemical_features(cdr3)
    feat_epi = physicochemical_features(epitope)
    min_len = min(len(feat_cdr3), len(feat_epi))
    return np.sqrt(np.sum((feat_cdr3[:min_len] - feat_epi[:min_len])**2))

def predict_binding_scores(df, epitopes=None):
    """Predict TCR-epitope binding scores for all CDR3 sequences."""
    if epitopes is None:
        epitopes = KNOWN_EPITOPES
    
    results = []
    for _, row in df.iterrows():
        cdr3 = row['cdr3_aa']
        best_score = 0
        best_epitope = None
        for epi_seq, epi_name in epitopes.items():
            dist = cdr3_epitope_distance(cdr3, epi_seq)
            # Convert distance to binding probability
            score = np.exp(-dist / 10.0) + np.random.uniform(0, 0.15)
            score = min(score, 1.0)
            if score > best_score:
                best_score = score
                best_epitope = epi_name
        results.append({
            'clone_id': row['clone_id'],
            'cdr3_aa': cdr3,
            'predicted_epitope': best_epitope,
            'binding_score': best_score,
        })
    return pd.DataFrame(results)

def build_cnn_features(df, max_len=20):
    """Build CNN-like feature matrix from CDR3 sequences."""
    features = []
    for _, row in df.iterrows():
        onehot = encode_cdr3(row['cdr3_aa'], max_len).flatten()
        physico = physicochemical_features(row['cdr3_aa'])
        features.append(np.concatenate([onehot, physico]))
    return np.array(features)

# ============================================================
# 5. Immune Age Estimation & Clonal Expansion Analysis
# ============================================================
def estimate_immune_age(diversity_df):
    """Estimate immune age from repertoire diversity features."""
    features = diversity_df[['shannon_entropy', 'simpson_index', 'clonality',
                             'n_clonotypes', 'hill_q1', 'hill_q2', 'gini_coefficient',
                             'top10_freq']].values
    actual_ages = diversity_df['age'].values
    
    scaler = StandardScaler()
    X = scaler.fit_transform(features)
    
    # Linear model for immune age
    from sklearn.linear_model import Ridge
    model = Ridge(alpha=1.0)
    model.fit(X, actual_ages)
    predicted_ages = model.predict(X)
    
    # Immune age acceleration = predicted - actual
    age_acceleration = predicted_ages - actual_ages
    
    return predicted_ages, age_acceleration, model

def analyze_clonal_expansion(df):
    """Analyze clonal expansion patterns across conditions."""
    results = []
    for sid in df['sample_id'].unique():
        sample = df[df['sample_id'] == sid].sort_values('frequency', ascending=False)
        label = sample['label'].iloc[0]
        freqs = sample['frequency'].values
        
        n_expanded = np.sum(freqs > 0.01)
        n_hyperexpanded = np.sum(freqs > 0.05)
        top1_freq = freqs[0]
        top10_freq = freqs[:10].sum()
        
        # Expansion evenness
        if n_expanded > 1:
            exp_freqs = freqs[freqs > 0.01]
            expansion_evenness = shannon_entropy(exp_freqs / exp_freqs.sum()) / np.log2(len(exp_freqs))
        else:
            expansion_evenness = 0
        
        results.append({
            'sample_id': sid,
            'label': label,
            'n_expanded': n_expanded,
            'n_hyperexpanded': n_hyperexpanded,
            'top1_freq': top1_freq,
            'top10_freq': top10_freq,
            'expansion_evenness': expansion_evenness,
        })
    return pd.DataFrame(results)

# ============================================================
# 6. ICB Response Prediction
# ============================================================
def build_icb_features(diversity_df, expansion_df):
    """Build feature matrix for ICB response prediction."""
    merged = diversity_df.merge(expansion_df, on=['sample_id', 'label'])
    feature_cols = ['shannon_entropy', 'simpson_index', 'chao1', 'clonality',
                    'hill_q1', 'hill_q2', 'gini_coefficient', 'top10_freq_x',
                    'n_expanded', 'n_hyperexpanded', 'top1_freq',
                    'expansion_evenness', 'n_clonotypes']
    X = merged[feature_cols].values
    return X, merged

def train_icb_predictor(X, y):
    """Train and evaluate ICB response predictors."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models = {
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'Gradient Boosting': GradientBoostingClassifier(n_estimators=100, random_state=42),
        'Logistic Regression': LogisticRegression(random_state=42, max_iter=1000),
    }
    
    results = {}
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    for name, model in models.items():
        auc_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
        acc_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='accuracy')
        
        model.fit(X_scaled, y)
        y_prob = model.predict_proba(X_scaled)[:, 1]
        y_pred = model.predict(X_scaled)
        
        results[name] = {
            'model': model,
            'scaler': scaler,
            'cv_auc_mean': auc_scores.mean(),
            'cv_auc_std': auc_scores.std(),
            'cv_acc_mean': acc_scores.mean(),
            'cv_acc_std': acc_scores.std(),
            'train_auc': roc_auc_score(y, y_prob),
            'y_prob': y_prob,
            'y_pred': y_pred,
        }
    
    return results

# ============================================================
# Visualization Functions
# ============================================================
def plot_diversity_comparison(diversity_df):
    """Plot diversity metrics comparison across conditions."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    metrics = ['shannon_entropy', 'simpson_index', 'chao1', 'clonality', 'hill_q1', 'gini_coefficient']
    titles = ['Shannon Entropy', 'Simpson Index', 'Chao1 Richness', 'Clonality', 'Hill Number (q=1)', 'Gini Coefficient']
    
    palette = {'healthy': '#2ecc71', 'responder': '#3498db', 'non_responder': '#e74c3c',
               'tumor': '#9b59b6', 'aged': '#f39c12', 'young': '#1abc9c'}
    
    for ax, metric, title in zip(axes.flatten(), metrics, titles):
        sns.boxplot(data=diversity_df, x='label', y=metric, ax=ax, palette=palette)
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.set_xlabel('')
        ax.tick_params(axis='x', rotation=30)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'diversity_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_clonal_expansion(expansion_df):
    """Plot clonal expansion patterns."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    palette = {'healthy': '#2ecc71', 'responder': '#3498db', 'non_responder': '#e74c3c',
               'tumor': '#9b59b6', 'aged': '#f39c12', 'young': '#1abc9c'}
    
    sns.boxplot(data=expansion_df, x='label', y='n_expanded', ax=axes[0], palette=palette)
    axes[0].set_title('Number of Expanded Clones (>1%)', fontweight='bold')
    axes[0].tick_params(axis='x', rotation=30)
    
    sns.boxplot(data=expansion_df, x='label', y='top1_freq', ax=axes[1], palette=palette)
    axes[1].set_title('Top Clone Frequency', fontweight='bold')
    axes[1].tick_params(axis='x', rotation=30)
    
    sns.boxplot(data=expansion_df, x='label', y='expansion_evenness', ax=axes[2], palette=palette)
    axes[2].set_title('Expansion Evenness', fontweight='bold')
    axes[2].tick_params(axis='x', rotation=30)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'clonal_expansion.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_hill_diversity_profile(diversity_df):
    """Plot Hill diversity profiles for different conditions."""
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'healthy': '#2ecc71', 'responder': '#3498db', 'non_responder': '#e74c3c',
              'tumor': '#9b59b6', 'aged': '#f39c12', 'young': '#1abc9c'}
    
    for label in diversity_df['label'].unique():
        subset = diversity_df[diversity_df['label'] == label]
        q_values = [0, 1, 2]
        means = [subset['hill_q0'].mean(), subset['hill_q1'].mean(), subset['hill_q2'].mean()]
        stds = [subset['hill_q0'].std(), subset['hill_q1'].std(), subset['hill_q2'].std()]
        ax.errorbar(q_values, means, yerr=stds, label=label, marker='o',
                    capsize=5, linewidth=2, color=colors.get(label, '#333'))
    
    ax.set_xlabel('Order q', fontsize=12)
    ax.set_ylabel('Hill Number', fontsize=12)
    ax.set_title('Hill Diversity Profile Across Conditions', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['q=0 (Richness)', 'q=1 (Shannon)', 'q=2 (Simpson)'])
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'hill_diversity_profile.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_immune_age(actual_ages, predicted_ages, labels):
    """Plot immune age estimation results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    colors_map = {'healthy': '#2ecc71', 'responder': '#3498db', 'non_responder': '#e74c3c',
                  'tumor': '#9b59b6', 'aged': '#f39c12', 'young': '#1abc9c'}
    colors = [colors_map.get(l, '#333') for l in labels]
    
    axes[0].scatter(actual_ages, predicted_ages, c=colors, alpha=0.7, s=60, edgecolors='white')
    min_age, max_age = min(actual_ages.min(), predicted_ages.min()), max(actual_ages.max(), predicted_ages.max())
    axes[0].plot([min_age, max_age], [min_age, max_age], 'k--', alpha=0.5, linewidth=2)
    r = np.corrcoef(actual_ages, predicted_ages)[0, 1]
    axes[0].set_xlabel('Chronological Age', fontsize=12)
    axes[0].set_ylabel('Predicted Immune Age', fontsize=12)
    axes[0].set_title(f'Immune Age Estimation (r={r:.3f})', fontsize=14, fontweight='bold')
    
    # Add legend
    for label, color in colors_map.items():
        if label in labels.values:
            axes[0].scatter([], [], c=color, label=label, s=60)
    axes[0].legend(fontsize=9)
    
    # Age acceleration by group
    acceleration = predicted_ages - actual_ages
    accel_df = pd.DataFrame({'label': labels, 'acceleration': acceleration})
    sns.boxplot(data=accel_df, x='label', y='acceleration', ax=axes[1], palette=colors_map)
    axes[1].axhline(y=0, color='black', linestyle='--', alpha=0.5)
    axes[1].set_title('Immune Age Acceleration by Condition', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('')
    axes[1].set_ylabel('Age Acceleration (years)', fontsize=12)
    axes[1].tick_params(axis='x', rotation=30)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'immune_age.png'), dpi=150, bbox_inches='tight')
    plt.close()
    return r

def plot_icb_roc_curves(results, y_true):
    """Plot ROC curves for ICB response prediction models."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = {'Random Forest': '#3498db', 'Gradient Boosting': '#e74c3c', 'Logistic Regression': '#2ecc71'}
    
    for name, res in results.items():
        fpr, tpr, _ = roc_curve(y_true, res['y_prob'])
        axes[0].plot(fpr, tpr, label=f"{name} (AUC={res['train_auc']:.3f})",
                     color=colors[name], linewidth=2)
    
    axes[0].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[0].set_xlabel('False Positive Rate', fontsize=12)
    axes[0].set_ylabel('True Positive Rate', fontsize=12)
    axes[0].set_title('ROC Curves — ICB Response Prediction', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=10)
    
    # Cross-validation AUC comparison
    model_names = list(results.keys())
    cv_aucs = [results[n]['cv_auc_mean'] for n in model_names]
    cv_stds = [results[n]['cv_auc_std'] for n in model_names]
    bars = axes[1].bar(model_names, cv_aucs, yerr=cv_stds, capsize=5,
                       color=[colors[n] for n in model_names], alpha=0.8)
    axes[1].set_ylabel('Cross-Validated AUC', fontsize=12)
    axes[1].set_title('5-Fold CV AUC Comparison', fontsize=14, fontweight='bold')
    axes[1].set_ylim(0.4, 1.0)
    for bar, auc in zip(bars, cv_aucs):
        axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                     f'{auc:.3f}', ha='center', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'icb_prediction.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_feature_importance(model, feature_names):
    """Plot feature importance for the best ICB prediction model."""
    if hasattr(model, 'feature_importances_'):
        importances = model.feature_importances_
    else:
        importances = np.abs(model.coef_[0])
    
    idx = np.argsort(importances)[::-1]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(importances)))
    ax.barh(range(len(importances)), importances[idx][::-1], color=colors)
    ax.set_yticks(range(len(importances)))
    ax.set_yticklabels([feature_names[i] for i in idx][::-1], fontsize=10)
    ax.set_xlabel('Feature Importance', fontsize=12)
    ax.set_title('Feature Importance — ICB Response Prediction (Gradient Boosting)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'feature_importance.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_public_tcr_network(public_df, df):
    """Visualize public TCR sharing network."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Frequency distribution of public TCRs
    axes[0].hist(public_df['n_samples'], bins=range(2, public_df['n_samples'].max()+2),
                 color='#3498db', edgecolor='white', alpha=0.8)
    axes[0].set_xlabel('Number of Samples Sharing TCR', fontsize=12)
    axes[0].set_ylabel('Count of Public TCRs', fontsize=12)
    axes[0].set_title('Public TCR Sharing Distribution', fontsize=14, fontweight='bold')
    
    # CDR3 length distribution by condition
    for label in df['label'].unique():
        subset = df[df['label'] == label]
        axes[1].hist(subset['cdr3_length'], bins=range(8, 22), alpha=0.5,
                     label=label, density=True, edgecolor='white')
    axes[1].set_xlabel('CDR3 Length', fontsize=12)
    axes[1].set_ylabel('Density', fontsize=12)
    axes[1].set_title('CDR3 Length Distribution by Condition', fontsize=14, fontweight='bold')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'public_tcr_analysis.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_binding_prediction(binding_df):
    """Plot TCR-epitope binding prediction results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Binding score distribution by epitope
    top_epitopes = binding_df['predicted_epitope'].value_counts().head(7).index
    subset = binding_df[binding_df['predicted_epitope'].isin(top_epitopes)]
    sns.boxplot(data=subset, x='predicted_epitope', y='binding_score', ax=axes[0],
                palette='Set2')
    axes[0].set_title('Binding Score Distribution by Predicted Epitope', fontsize=13, fontweight='bold')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].set_xlabel('')
    
    # Binding score vs CDR3 length
    axes[1].scatter(binding_df['cdr3_aa'].str.len(), binding_df['binding_score'],
                    alpha=0.3, s=20, c='#3498db')
    axes[1].set_xlabel('CDR3 Length', fontsize=12)
    axes[1].set_ylabel('Binding Score', fontsize=12)
    axes[1].set_title('Binding Score vs CDR3 Length', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'binding_prediction.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_vgene_usage(df):
    """Plot V gene usage across conditions."""
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Get top V genes
    top_v = df['v_gene'].value_counts().head(15).index
    usage = df[df['v_gene'].isin(top_v)].groupby(['label', 'v_gene']).size().unstack(fill_value=0)
    usage = usage.div(usage.sum(axis=1), axis=0)
    
    usage.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.8)
    ax.set_title('V Gene Usage Distribution by Condition', fontsize=14, fontweight='bold')
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_xlabel('')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.tick_params(axis='x', rotation=30)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'vgene_usage.png'), dpi=150, bbox_inches='tight')
    plt.close()

def plot_clone_frequency_distribution(df):
    """Plot clone frequency rank distribution (Zipf plot)."""
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = {'healthy': '#2ecc71', 'responder': '#3498db', 'non_responder': '#e74c3c',
              'tumor': '#9b59b6', 'aged': '#f39c12', 'young': '#1abc9c'}
    
    for label in df['label'].unique():
        subset = df[df['label'] == label]
        # Average across samples
        sample_ids = subset['sample_id'].unique()[:3]
        for sid in sample_ids:
            sample = subset[subset['sample_id'] == sid].sort_values('frequency', ascending=False)
            freqs = sample['frequency'].values
            ranks = np.arange(1, len(freqs) + 1)
            ax.plot(ranks, freqs, alpha=0.4, color=colors.get(label, '#333'), linewidth=1)
        # Plot mean line
        ax.plot([], [], color=colors.get(label, '#333'), linewidth=2, label=label)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Clone Rank', fontsize=12)
    ax.set_ylabel('Clone Frequency', fontsize=12)
    ax.set_title('Clone Frequency Rank Distribution (Zipf Plot)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'clone_frequency_rank.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Main Pipeline Execution
# ============================================================
def main():
    print("=" * 70)
    print("TCR Repertoire Analysis Pipeline — Immune State Estimation")
    print("=" * 70)
    
    # --- Generate synthetic data ---
    print("\n[1/6] Generating synthetic TCR-seq data and preprocessing...")
    n_clones = 500
    n_per_group = 15
    
    labels = (['healthy'] * n_per_group +
              ['responder'] * n_per_group +
              ['non_responder'] * n_per_group +
              ['tumor'] * n_per_group +
              ['aged'] * n_per_group +
              ['young'] * n_per_group)
    
    ages = (list(np.random.randint(30, 60, n_per_group)) +  # healthy
            list(np.random.randint(40, 70, n_per_group)) +  # responder
            list(np.random.randint(40, 70, n_per_group)) +  # non_responder
            list(np.random.randint(35, 65, n_per_group)) +  # tumor
            list(np.random.randint(65, 85, n_per_group)) +  # aged
            list(np.random.randint(20, 35, n_per_group)))   # young
    
    df = generate_tcr_repertoire(n_clones, len(labels), labels, ages)
    df = preprocess_repertoire(df)
    
    print(f"  Total records: {len(df):,}")
    print(f"  Samples: {df['sample_id'].nunique()}")
    print(f"  Groups: {df['label'].unique().tolist()}")
    print(f"  Valid CDR3 rate: {df['valid_cdr3'].mean()*100:.1f}%")
    
    # Save preprocessed data
    df.to_csv(os.path.join(DATA_DIR, 'preprocessed_repertoire.csv'), index=False)
    
    # --- Diversity Metrics ---
    print("\n[2/6] Computing repertoire diversity metrics...")
    diversity_df = compute_diversity_metrics(df)
    diversity_df.to_csv(os.path.join(DATA_DIR, 'diversity_metrics.csv'), index=False)
    
    print("  Diversity metrics by group (mean ± std):")
    for label in diversity_df['label'].unique():
        sub = diversity_df[diversity_df['label'] == label]
        print(f"    {label:15s}: Shannon={sub['shannon_entropy'].mean():.2f}±{sub['shannon_entropy'].std():.2f}, "
              f"Clonality={sub['clonality'].mean():.3f}±{sub['clonality'].std():.3f}, "
              f"Chao1={sub['chao1'].mean():.0f}±{sub['chao1'].std():.0f}")
    
    plot_diversity_comparison(diversity_df)
    plot_hill_diversity_profile(diversity_df)
    plot_clone_frequency_distribution(df)
    print("  → Figures saved: diversity_comparison.png, hill_diversity_profile.png, clone_frequency_rank.png")
    
    # --- Public TCR & HLA ---
    print("\n[3/6] Identifying public TCRs and predicting HLA restriction...")
    public_df = identify_public_tcrs(df, min_samples=2)
    print(f"  Public TCRs found: {len(public_df)}")
    if len(public_df) > 0:
        print(f"  Max sharing: {public_df['n_samples'].max()} samples")
        public_df['predicted_hla'], public_df['hla_score'] = zip(
            *public_df['cdr3_aa'].apply(predict_hla_restriction)
        )
        public_df.to_csv(os.path.join(DATA_DIR, 'public_tcrs.csv'), index=False)
    
    plot_public_tcr_network(public_df, df)
    plot_vgene_usage(df)
    print("  → Figures saved: public_tcr_analysis.png, vgene_usage.png")
    
    # --- Binding Prediction ---
    print("\n[4/6] Predicting TCR-epitope binding...")
    # Sample for binding prediction
    sample_df = df.groupby('sample_id').head(50)
    binding_df = predict_binding_scores(sample_df)
    binding_df.to_csv(os.path.join(DATA_DIR, 'binding_predictions.csv'), index=False)
    
    print(f"  Predictions: {len(binding_df):,}")
    print(f"  Epitope distribution:")
    for epi, count in binding_df['predicted_epitope'].value_counts().head(5).items():
        print(f"    {epi}: {count} ({count/len(binding_df)*100:.1f}%)")
    
    plot_binding_prediction(binding_df)
    print("  → Figure saved: binding_prediction.png")
    
    # --- Immune Age ---
    print("\n[5/6] Estimating immune age and analyzing clonal expansion...")
    predicted_ages, age_accel, age_model = estimate_immune_age(diversity_df)
    expansion_df = analyze_clonal_expansion(df)
    expansion_df.to_csv(os.path.join(DATA_DIR, 'clonal_expansion.csv'), index=False)
    
    r = plot_immune_age(diversity_df['age'].values, predicted_ages, diversity_df['label'])
    print(f"  Immune age correlation (r): {r:.4f}")
    print(f"  Mean age acceleration by group:")
    for label in diversity_df['label'].unique():
        mask = diversity_df['label'] == label
        print(f"    {label:15s}: {age_accel[mask].mean():+.2f} ± {age_accel[mask].std():.2f} years")
    
    plot_clonal_expansion(expansion_df)
    print("  → Figures saved: immune_age.png, clonal_expansion.png")
    
    # --- ICB Response Prediction ---
    print("\n[6/6] Training ICB response prediction models...")
    # Create binary labels: responder vs non_responder
    icb_mask = diversity_df['label'].isin(['responder', 'non_responder'])
    icb_diversity = diversity_df[icb_mask].copy()
    icb_expansion = expansion_df[expansion_df['label'].isin(['responder', 'non_responder'])].copy()
    
    X_icb, merged_df = build_icb_features(icb_diversity, icb_expansion)
    y_icb = (merged_df['label'] == 'responder').astype(int).values
    
    feature_names = ['Shannon Entropy', 'Simpson Index', 'Chao1', 'Clonality',
                     'Hill q=1', 'Hill q=2', 'Gini', 'Top10 Freq',
                     'N Expanded', 'N Hyperexpanded', 'Top1 Freq',
                     'Expansion Evenness', 'N Clonotypes']
    
    icb_results = train_icb_predictor(X_icb, y_icb)
    
    print("  Model performance (5-fold CV):")
    for name, res in icb_results.items():
        print(f"    {name:25s}: AUC={res['cv_auc_mean']:.3f}±{res['cv_auc_std']:.3f}, "
              f"Acc={res['cv_acc_mean']:.3f}±{res['cv_acc_std']:.3f}")
    
    plot_icb_roc_curves(icb_results, y_icb)
    plot_feature_importance(icb_results['Gradient Boosting']['model'], feature_names)
    print("  → Figures saved: icb_prediction.png, feature_importance.png")
    
    # --- Summary ---
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    
    summary = {
        'total_records': len(df),
        'n_samples': df['sample_id'].nunique(),
        'n_clonotypes_total': df['clonotype'].nunique(),
        'n_public_tcrs': len(public_df),
        'binding_predictions': len(binding_df),
        'immune_age_r': float(r),
        'icb_models': {name: {'cv_auc': float(res['cv_auc_mean']), 'cv_acc': float(res['cv_acc_mean'])}
                       for name, res in icb_results.items()},
        'diversity_summary': {
            label: {
                'shannon_mean': float(diversity_df[diversity_df['label']==label]['shannon_entropy'].mean()),
                'clonality_mean': float(diversity_df[diversity_df['label']==label]['clonality'].mean()),
            }
            for label in diversity_df['label'].unique()
        }
    }
    
    with open(os.path.join(DATA_DIR, 'pipeline_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nFiles generated:")
    print(f"  Data: {os.listdir(DATA_DIR)}")
    print(f"  Figures: {os.listdir(FIGURES_DIR)}")
    
    return summary

if __name__ == '__main__':
    summary = main()
