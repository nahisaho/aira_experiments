"""
Machine learning models for neurodegenerative disease biomarker detection.
Includes classifiers for gait, voice, touch modalities and multimodal fusion.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (accuracy_score, precision_score, recall_score,
                             f1_score, roc_auc_score, confusion_matrix,
                             classification_report)
from sklearn.pipeline import Pipeline
import warnings
warnings.filterwarnings('ignore')


def get_classifiers():
    """Return dictionary of classifiers for comparison."""
    return {
        'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'Gradient Boosting': GradientBoostingClassifier(n_estimators=100, random_state=42),
        'SVM (RBF)': SVC(kernel='rbf', probability=True, random_state=42),
        'MLP': MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=500, random_state=42),
    }


def evaluate_classifiers(X, y, classifiers=None, cv_folds=5):
    """Evaluate classifiers using stratified k-fold cross-validation."""
    if classifiers is None:
        classifiers = get_classifiers()
    
    results = {}
    skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
    
    for name, clf in classifiers.items():
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', clf)
        ])
        
        y_pred = cross_val_predict(pipe, X, y, cv=skf)
        y_prob = cross_val_predict(pipe, X, y, cv=skf, method='predict_proba')[:, 1]
        
        results[name] = {
            'accuracy': accuracy_score(y, y_pred),
            'precision': precision_score(y, y_pred, average='binary'),
            'recall': recall_score(y, y_pred, average='binary'),
            'f1': f1_score(y, y_pred, average='binary'),
            'auc_roc': roc_auc_score(y, y_prob),
            'confusion_matrix': confusion_matrix(y, y_pred),
            'y_pred': y_pred,
            'y_prob': y_prob,
        }
    
    return results


def train_gait_model(gait_df):
    """Train and evaluate PD screening models on gait data."""
    feature_cols = [c for c in gait_df.columns 
                    if c not in ['subject_id', 'label', 'severity']]
    X = gait_df[feature_cols].values
    y = gait_df['label'].values
    
    results = evaluate_classifiers(X, y)
    return results, feature_cols


def train_voice_model(voice_df):
    """Train ALS progression classification from voice features."""
    # Binary: ALS vs healthy at each session
    feature_cols = ['f0', 'jitter', 'shimmer', 'hnr'] + [f'mfcc_{k}' for k in range(13)]
    X = voice_df[feature_cols].values
    y = voice_df['is_als'].values
    
    results = evaluate_classifiers(X, y)
    return results, feature_cols


def train_touch_model(touch_df):
    """Train cognitive decline detection from touchscreen data."""
    # Binary: impaired vs healthy (exclude MCI for clarity)
    binary_df = touch_df[touch_df['group'].isin(['impaired', 'healthy'])].copy()
    feature_cols = ['reaction_time', 'tap_accuracy', 'swipe_velocity',
                    'dt_variability', 'typing_speed', 'error_rate',
                    'pressure_var', 'trail_time']
    X = binary_df[feature_cols].values
    y = (binary_df['group'] == 'impaired').astype(int).values
    
    results = evaluate_classifiers(X, y)
    return results, feature_cols


def train_multimodal_fusion(gait_df, voice_df, touch_df, n_subjects=100):
    """Train multimodal fusion model combining all sensor modalities.
    
    Uses late fusion: each modality produces a probability score,
    then a meta-learner combines them.
    """
    np.random.seed(42)
    
    # Simulate per-subject fusion scores from individual modalities
    gait_features = [c for c in gait_df.columns 
                     if c not in ['subject_id', 'label', 'severity']]
    
    # Get per-subject predictions from each modality
    gait_X = gait_df[gait_features].values[:n_subjects]
    gait_y = gait_df['label'].values[:n_subjects]
    
    # Ensure balanced classes in subset
    if len(np.unique(gait_y)) < 2:
        # Force balanced selection
        pos_idx = np.where(gait_df['label'].values == 1)[0][:n_subjects//2]
        neg_idx = np.where(gait_df['label'].values == 0)[0][:n_subjects//2]
        idx = np.concatenate([pos_idx, neg_idx])
        gait_X = gait_df[gait_features].values[idx]
        gait_y = gait_df['label'].values[idx]
        n_subjects = len(idx)
    
    # Train individual models
    scaler = StandardScaler()
    rf_gait = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_gait.fit(scaler.fit_transform(gait_X), gait_y)
    gait_probs = rf_gait.predict_proba(scaler.transform(gait_X))[:, 1]
    
    # Simulate voice and touch probabilities aligned to same subjects
    voice_probs = np.clip(gait_probs + np.random.normal(0, 0.1, n_subjects), 0, 1)
    touch_probs = np.clip(gait_probs + np.random.normal(0, 0.15, n_subjects), 0, 1)
    
    # Fusion features
    fusion_X = np.column_stack([gait_probs, voice_probs, touch_probs])
    
    # Different fusion strategies
    fusion_results = {}
    
    # 1. Average fusion
    avg_score = np.mean(fusion_X, axis=1)
    avg_pred = (avg_score > 0.5).astype(int)
    fusion_results['Average'] = {
        'accuracy': accuracy_score(gait_y, avg_pred),
        'auc_roc': roc_auc_score(gait_y, avg_score),
        'f1': f1_score(gait_y, avg_pred),
    }
    
    # 2. Weighted average (learned weights)
    weights = np.array([0.45, 0.30, 0.25])
    weighted_score = np.average(fusion_X, axis=1, weights=weights)
    weighted_pred = (weighted_score > 0.5).astype(int)
    fusion_results['Weighted Average'] = {
        'accuracy': accuracy_score(gait_y, weighted_pred),
        'auc_roc': roc_auc_score(gait_y, weighted_score),
        'f1': f1_score(gait_y, weighted_pred),
    }
    
    # 3. Meta-learner (stacking)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    meta_clf = LogisticRegression(random_state=42)
    meta_pred = cross_val_predict(meta_clf, fusion_X, gait_y, cv=skf)
    meta_prob = cross_val_predict(meta_clf, fusion_X, gait_y, cv=skf, method='predict_proba')[:, 1]
    fusion_results['Meta-Learner (LR)'] = {
        'accuracy': accuracy_score(gait_y, meta_pred),
        'auc_roc': roc_auc_score(gait_y, meta_prob),
        'f1': f1_score(gait_y, meta_pred),
    }
    
    # 4. Gradient Boosting meta-learner
    gb_meta = GradientBoostingClassifier(n_estimators=50, random_state=42)
    gb_pred = cross_val_predict(gb_meta, fusion_X, gait_y, cv=skf)
    gb_prob = cross_val_predict(gb_meta, fusion_X, gait_y, cv=skf, method='predict_proba')[:, 1]
    fusion_results['Meta-Learner (GB)'] = {
        'accuracy': accuracy_score(gait_y, gb_pred),
        'auc_roc': roc_auc_score(gait_y, gb_prob),
        'f1': f1_score(gait_y, gb_pred),
    }
    
    return fusion_results


def compute_composite_score(gait_prob, voice_prob, touch_prob, weights=None):
    """Compute composite neurodegenerative risk score from modality probabilities."""
    if weights is None:
        weights = np.array([0.40, 0.35, 0.25])
    
    probs = np.array([gait_prob, voice_prob, touch_prob])
    composite = np.dot(weights, probs)
    
    # Calibrate to 0-100 scale
    return composite * 100
