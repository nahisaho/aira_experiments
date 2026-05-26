#!/usr/bin/env python3
"""
Automated Quality Control and Anomaly Detection for Large-Scale Scientific Data
================================================================================
This experiment implements a streaming anomaly detection pipeline with:
1. Changepoint detection (PELT / BOCPD)
2. Multivariate outlier detection (Isolation Forest / Deep SVDD-like)
3. Physics-constrained anomaly scoring
4. Concept drift detection and model retraining triggers
5. Explainable anomaly detection (feature attribution)
6. CERN/LIGO-style large-scale data application design
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats, signal
from scipy.special import logsumexp
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (precision_score, recall_score, f1_score,
                             roc_auc_score, confusion_matrix, roc_curve,
                             precision_recall_curve, average_precision_score)
from sklearn.decomposition import PCA
import ruptures
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
FIGDIR = 'figures'
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================================
# 1. Synthetic Scientific Data Generation
# ============================================================================
def generate_scientific_data(n_samples=5000, n_features=8):
    """Generate synthetic data mimicking CERN/LIGO-style sensor readings."""
    t = np.linspace(0, 100, n_samples)
    
    # Base physics signals
    temperature = 20 + 0.5 * np.sin(2 * np.pi * t / 24) + np.random.normal(0, 0.3, n_samples)
    pressure = 1013.25 + 2 * np.sin(2 * np.pi * t / 12) + np.random.normal(0, 1.0, n_samples)
    voltage = 220 + 5 * np.sin(2 * np.pi * t / 6) + np.random.normal(0, 1.5, n_samples)
    current = voltage / (50 + np.random.normal(0, 0.5, n_samples))  # Ohm's law
    magnetic_field = 1.5 + 0.1 * np.sin(2 * np.pi * t / 48) + np.random.normal(0, 0.02, n_samples)
    beam_intensity = 1e6 * (1 + 0.05 * np.sin(2 * np.pi * t / 8)) + np.random.normal(0, 1e4, n_samples)
    luminosity = beam_intensity * (1 + np.random.normal(0, 0.01, n_samples)) * 1e-4
    event_rate = luminosity * 0.1 + np.random.poisson(5, n_samples)
    
    data = np.column_stack([temperature, pressure, voltage, current,
                            magnetic_field, beam_intensity, luminosity, event_rate])
    feature_names = ['Temperature', 'Pressure', 'Voltage', 'Current',
                     'MagneticField', 'BeamIntensity', 'Luminosity', 'EventRate']
    
    # Inject anomalies
    labels = np.zeros(n_samples, dtype=int)
    
    # Point anomalies (sensor spikes)
    spike_idx = np.random.choice(n_samples, 30, replace=False)
    for idx in spike_idx:
        feat = np.random.randint(0, n_features)
        data[idx, feat] += 10 * np.std(data[:, feat]) * np.random.choice([-1, 1])
    labels[spike_idx] = 1
    
    # Contextual anomalies (physics violations: current not matching voltage/resistance)
    ctx_start = 2000
    ctx_end = 2050
    data[ctx_start:ctx_end, 3] = data[ctx_start:ctx_end, 2] / 10  # Wrong resistance
    labels[ctx_start:ctx_end] = 1
    
    # Collective anomalies (gradual drift in a segment)
    drift_start = 3500
    drift_end = 3600
    drift_ramp = np.linspace(0, 15, drift_end - drift_start)
    data[drift_start:drift_end, 0] += drift_ramp  # Temperature drift
    labels[drift_start:drift_end] = 1
    
    # Changepoints
    cp1, cp2 = 1500, 3000
    data[cp1:, 1] += 5  # Pressure shift
    data[cp2:, 4] -= 0.3  # Magnetic field shift
    
    return t, data, labels, feature_names, [cp1, cp2]


# ============================================================================
# 2. Changepoint Detection
# ============================================================================
def run_pelt_detection(signal_data, pen=10):
    """PELT changepoint detection."""
    algo = ruptures.Pelt(model="rbf", min_size=50).fit(signal_data)
    result = algo.predict(pen=pen)
    return result[:-1]  # Remove the last element (end of signal)

def run_bocpd(data, hazard_rate=1/200, mu0=0, kappa0=1, alpha0=1, beta0=1):
    """Bayesian Online Changepoint Detection."""
    n = len(data)
    R = np.zeros((n + 1, n + 1))
    R[0, 0] = 1.0
    
    maxes = np.zeros(n)
    muT = np.array([mu0])
    kappaT = np.array([kappa0])
    alphaT = np.array([alpha0])
    betaT = np.array([beta0])
    
    changepoint_prob = np.zeros(n)
    
    for t in range(n):
        x = data[t]
        predprobs = stats.t.pdf(x, 2 * alphaT, loc=muT,
                                scale=np.sqrt(betaT * (kappaT + 1) / (alphaT * kappaT)))
        
        H = 1 / hazard_rate
        R[1:t+2, t+1] = R[:t+1, t] * predprobs * (1 - 1/H)
        R[0, t+1] = np.sum(R[:t+1, t] * predprobs * (1/H))
        
        evidence = np.sum(R[:t+2, t+1])
        if evidence > 0:
            R[:t+2, t+1] /= evidence
        
        changepoint_prob[t] = R[0, t+1]
        
        # Update sufficient statistics
        new_mu = np.concatenate([[mu0], (kappaT * muT + x) / (kappaT + 1)])
        new_kappa = np.concatenate([[kappa0], kappaT + 1])
        new_alpha = np.concatenate([[alpha0], alphaT + 0.5])
        new_beta = np.concatenate([[beta0], betaT + kappaT * (x - muT)**2 / (2 * (kappaT + 1))])
        
        muT = new_mu
        kappaT = new_kappa
        alphaT = new_alpha
        betaT = new_beta
    
    return changepoint_prob

def evaluate_changepoints(detected, true_cps, tolerance=50):
    """Evaluate changepoint detection with tolerance window."""
    tp = 0
    for tcp in true_cps:
        if any(abs(d - tcp) <= tolerance for d in detected):
            tp += 1
    fp = len(detected) - tp
    fn = len(true_cps) - tp
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return {'precision': precision, 'recall': recall, 'f1': f1, 'tp': tp, 'fp': fp, 'fn': fn}


# ============================================================================
# 3. Multivariate Outlier Detection
# ============================================================================
def isolation_forest_detection(data, contamination=0.05):
    """Standard Isolation Forest."""
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    clf = IsolationForest(contamination=contamination, random_state=42, n_estimators=200)
    clf.fit(data_scaled)
    scores = -clf.decision_function(data_scaled)  # Higher = more anomalous
    predictions = clf.predict(data_scaled)
    predictions = (predictions == -1).astype(int)
    return scores, predictions, clf, scaler

def deep_svdd_like_detection(data, nu=0.05):
    """Deep SVDD-inspired detection using PCA + hypersphere approach."""
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    
    # Use PCA as a learned representation (simulating deep network)
    pca = PCA(n_components=min(4, data.shape[1]))
    data_embedded = pca.fit_transform(data_scaled)
    
    # Find center of normal data (use majority of data)
    center = np.median(data_embedded, axis=0)
    
    # Compute distances to center
    distances = np.sqrt(np.sum((data_embedded - center) ** 2, axis=1))
    
    # Determine threshold using quantile
    threshold = np.quantile(distances, 1 - nu)
    predictions = (distances > threshold).astype(int)
    scores = distances / threshold  # Normalized scores
    
    return scores, predictions, center, pca, scaler


# ============================================================================
# 4. Physics-Constrained Anomaly Scoring
# ============================================================================
def physics_constrained_scoring(data, feature_names):
    """Compute physics constraint violations."""
    n = data.shape[0]
    constraint_scores = np.zeros(n)
    violation_details = []
    
    # Constraint 1: Ohm's Law (V = I * R, expect R ~ 50 ohms)
    voltage_idx = feature_names.index('Voltage')
    current_idx = feature_names.index('Current')
    expected_current = data[:, voltage_idx] / 50.0
    ohm_violation = np.abs(data[:, current_idx] - expected_current) / np.std(expected_current)
    
    # Constraint 2: Luminosity ~ BeamIntensity * constant
    beam_idx = feature_names.index('BeamIntensity')
    lumi_idx = feature_names.index('Luminosity')
    expected_lumi = data[:, beam_idx] * 1e-4
    lumi_violation = np.abs(data[:, lumi_idx] - expected_lumi) / np.std(expected_lumi)
    
    # Constraint 3: EventRate ~ Luminosity * cross_section
    event_idx = feature_names.index('EventRate')
    expected_events = data[:, lumi_idx] * 0.1
    event_violation = np.abs(data[:, event_idx] - expected_events) / np.std(expected_events)
    
    # Constraint 4: Temperature physical bounds
    temp_idx = feature_names.index('Temperature')
    temp_violation = np.maximum(0, np.abs(data[:, temp_idx] - 20) - 5) / 5.0
    
    # Combine constraint violations
    constraint_scores = (ohm_violation + lumi_violation + event_violation + temp_violation) / 4.0
    
    return constraint_scores, {
        'ohm_law': ohm_violation,
        'luminosity': lumi_violation,
        'event_rate': event_violation,
        'temperature': temp_violation
    }


# ============================================================================
# 5. Concept Drift Detection
# ============================================================================
def adwin_drift_detection(data, delta=0.002):
    """Simplified ADWIN-like drift detection."""
    n = len(data)
    drift_points = []
    window = []
    
    for i in range(n):
        window.append(data[i])
        if len(window) < 30:
            continue
        
        # Check for drift by comparing sub-windows
        best_cut = -1
        best_stat = 0
        w = np.array(window)
        
        for cut in range(max(10, len(w)//4), len(w) - max(10, len(w)//4)):
            w1, w2 = w[:cut], w[cut:]
            mu1, mu2 = np.mean(w1), np.mean(w2)
            n1, n2 = len(w1), len(w2)
            
            # Hoeffding bound
            m = 1.0 / (1.0/n1 + 1.0/n2)
            eps = np.sqrt(np.log(4.0/delta) / (2.0 * m))
            
            if abs(mu1 - mu2) >= eps:
                if abs(mu1 - mu2) > best_stat:
                    best_stat = abs(mu1 - mu2)
                    best_cut = cut
        
        if best_cut >= 0:
            drift_points.append(i)
            window = window[best_cut:]  # Shrink window
    
    return drift_points

def page_hinkley_test(data, threshold=50, alpha=0.005):
    """Page-Hinkley test for drift detection."""
    n = len(data)
    drift_points = []
    m_t = 0  # Cumulative sum
    M_t = 0  # Minimum of cumulative sum
    sum_x = 0
    
    for t in range(1, n):
        sum_x += data[t]
        mean_t = sum_x / t
        m_t += (data[t] - mean_t - alpha)
        M_t = min(M_t, m_t)
        
        if m_t - M_t > threshold:
            drift_points.append(t)
            m_t = 0
            M_t = 0
            sum_x = 0
    
    return drift_points


# ============================================================================
# 6. Explainable Anomaly Detection
# ============================================================================
def compute_feature_attribution(data, scores, feature_names, top_k=3):
    """Compute feature-level attribution for anomalies (SHAP-inspired)."""
    n, p = data.shape
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    
    # Marginal contribution of each feature
    attributions = np.zeros((n, p))
    baseline = np.mean(data_scaled, axis=0)
    
    for j in range(p):
        # Measure the contribution of each feature to anomaly score
        perturbed = data_scaled.copy()
        perturbed[:, j] = baseline[j]
        
        clf = IsolationForest(contamination=0.05, random_state=42, n_estimators=100)
        clf.fit(data_scaled)
        original_scores = -clf.decision_function(data_scaled)
        perturbed_scores = -clf.decision_function(perturbed)
        
        attributions[:, j] = original_scores - perturbed_scores
    
    # For each anomaly, find top contributing features
    anomaly_mask = scores > np.quantile(scores, 0.95)
    explanations = []
    for i in np.where(anomaly_mask)[0]:
        top_features = np.argsort(attributions[i])[-top_k:][::-1]
        explanation = {
            'index': int(i),
            'score': float(scores[i]),
            'top_features': [(feature_names[f], float(attributions[i, f])) for f in top_features]
        }
        explanations.append(explanation)
    
    return attributions, explanations


# ============================================================================
# 7. Combined Anomaly Score
# ============================================================================
def combined_anomaly_pipeline(data, feature_names, labels):
    """Full streaming-compatible anomaly detection pipeline."""
    # Isolation Forest
    if_scores, if_preds, if_clf, if_scaler = isolation_forest_detection(data)
    
    # Deep SVDD-like
    svdd_scores, svdd_preds, center, pca, svdd_scaler = deep_svdd_like_detection(data)
    
    # Physics constraints
    phys_scores, phys_details = physics_constrained_scoring(data, feature_names)
    
    # Normalize all scores to [0, 1]
    if_norm = (if_scores - if_scores.min()) / (if_scores.max() - if_scores.min())
    svdd_norm = (svdd_scores - svdd_scores.min()) / (svdd_scores.max() - svdd_scores.min())
    phys_norm = (phys_scores - phys_scores.min()) / (phys_scores.max() - phys_scores.min() + 1e-10)
    
    # Weighted combination
    combined = 0.35 * if_norm + 0.30 * svdd_norm + 0.35 * phys_norm
    
    # Threshold selection using Otsu's method
    hist, bin_edges = np.histogram(combined, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    total = np.sum(hist)
    best_thresh = 0.5
    best_var = 0
    
    for i in range(1, len(hist)):
        w0 = np.sum(hist[:i]) / total
        w1 = np.sum(hist[i:]) / total
        if w0 == 0 or w1 == 0:
            continue
        mu0 = np.sum(bin_centers[:i] * hist[:i]) / np.sum(hist[:i])
        mu1 = np.sum(bin_centers[i:] * hist[i:]) / np.sum(hist[i:])
        var_between = w0 * w1 * (mu0 - mu1) ** 2
        if var_between > best_var:
            best_var = var_between
            best_thresh = bin_centers[i]
    
    combined_preds = (combined > best_thresh).astype(int)
    
    return {
        'if_scores': if_scores, 'if_preds': if_preds,
        'svdd_scores': svdd_scores, 'svdd_preds': svdd_preds,
        'phys_scores': phys_scores, 'phys_details': phys_details,
        'combined': combined, 'combined_preds': combined_preds,
        'threshold': best_thresh,
        'if_norm': if_norm, 'svdd_norm': svdd_norm, 'phys_norm': phys_norm
    }


# ============================================================================
# Plotting Functions
# ============================================================================
def plot_data_overview(t, data, labels, feature_names):
    fig, axes = plt.subplots(4, 2, figsize=(16, 14))
    fig.suptitle('Scientific Data Overview with Anomaly Labels', fontsize=14, fontweight='bold')
    for i, (ax, name) in enumerate(zip(axes.flat, feature_names)):
        ax.plot(t, data[:, i], 'b-', alpha=0.5, linewidth=0.5)
        anomaly_mask = labels == 1
        ax.scatter(t[anomaly_mask], data[anomaly_mask, i], c='red', s=8, zorder=5, label='Anomaly')
        ax.set_title(name, fontsize=10)
        ax.set_xlabel('Time')
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/data_overview.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_changepoint_results(t, data, feature_names, pelt_cps, bocpd_prob, true_cps):
    fig, axes = plt.subplots(3, 1, figsize=(16, 10))
    
    # Pressure with PELT
    ax = axes[0]
    ax.plot(t, data[:, 1], 'b-', alpha=0.7, label='Pressure')
    for cp in pelt_cps:
        ax.axvline(t[cp], color='red', linestyle='--', alpha=0.8, label='PELT CP' if cp == pelt_cps[0] else '')
    for cp in true_cps:
        ax.axvline(t[cp], color='green', linestyle=':', alpha=0.8, label='True CP' if cp == true_cps[0] else '')
    ax.set_title('PELT Changepoint Detection (Pressure)', fontweight='bold')
    ax.legend()
    
    # Magnetic field with PELT
    ax = axes[1]
    ax.plot(t, data[:, 4], 'b-', alpha=0.7, label='MagneticField')
    for cp in pelt_cps:
        ax.axvline(t[cp], color='red', linestyle='--', alpha=0.8)
    for cp in true_cps:
        ax.axvline(t[cp], color='green', linestyle=':', alpha=0.8)
    ax.set_title('PELT Changepoint Detection (Magnetic Field)', fontweight='bold')
    
    # BOCPD probability
    ax = axes[2]
    ax.plot(t, bocpd_prob, 'purple', alpha=0.8, label='CP Probability')
    for cp in true_cps:
        ax.axvline(t[cp], color='green', linestyle=':', alpha=0.8, label='True CP' if cp == true_cps[0] else '')
    ax.set_title('BOCPD Changepoint Probability', fontweight='bold')
    ax.set_xlabel('Time')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/changepoint_detection.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_anomaly_scores(t, results, labels):
    fig, axes = plt.subplots(4, 1, figsize=(16, 12))
    
    for ax, (name, scores) in zip(axes, [
        ('Isolation Forest', results['if_norm']),
        ('Deep SVDD-like', results['svdd_norm']),
        ('Physics Constraints', results['phys_norm']),
        ('Combined Score', results['combined'])
    ]):
        ax.plot(t, scores, 'b-', alpha=0.5, linewidth=0.5)
        anomaly_mask = labels == 1
        ax.scatter(t[anomaly_mask], scores[anomaly_mask], c='red', s=8, zorder=5)
        if name == 'Combined Score':
            ax.axhline(results['threshold'], color='orange', linestyle='--', label=f'Threshold={results["threshold"]:.3f}')
            ax.legend()
        ax.set_title(f'{name} Anomaly Scores', fontweight='bold')
        ax.set_xlabel('Time')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/anomaly_scores.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_roc_pr_curves(labels, results):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for name, scores in [('IsolationForest', results['if_norm']),
                         ('DeepSVDD-like', results['svdd_norm']),
                         ('Physics', results['phys_norm']),
                         ('Combined', results['combined'])]:
        fpr, tpr, _ = roc_curve(labels, scores)
        auc = roc_auc_score(labels, scores)
        axes[0].plot(fpr, tpr, label=f'{name} (AUC={auc:.3f})')
        
        prec, rec, _ = precision_recall_curve(labels, scores)
        ap = average_precision_score(labels, scores)
        axes[1].plot(rec, prec, label=f'{name} (AP={ap:.3f})')
    
    axes[0].plot([0,1], [0,1], 'k--', alpha=0.3)
    axes[0].set_xlabel('False Positive Rate')
    axes[0].set_ylabel('True Positive Rate')
    axes[0].set_title('ROC Curves', fontweight='bold')
    axes[0].legend(fontsize=9)
    
    axes[1].set_xlabel('Recall')
    axes[1].set_ylabel('Precision')
    axes[1].set_title('Precision-Recall Curves', fontweight='bold')
    axes[1].legend(fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/roc_pr_curves.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_physics_constraints(t, phys_details):
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    fig.suptitle('Physics Constraint Violations', fontsize=14, fontweight='bold')
    
    for ax, (name, violation) in zip(axes.flat, phys_details.items()):
        ax.plot(t, violation, 'b-', alpha=0.5, linewidth=0.5)
        high_mask = violation > np.quantile(violation, 0.95)
        ax.scatter(t[high_mask], violation[high_mask], c='red', s=8, zorder=5)
        ax.set_title(name.replace('_', ' ').title(), fontweight='bold')
        ax.set_xlabel('Time')
        ax.set_ylabel('Violation Score')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/physics_constraints.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_drift_detection(t, data, adwin_drifts, ph_drifts, true_cps):
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    
    ax = axes[0]
    ax.plot(t, data[:, 0], 'b-', alpha=0.7, label='Temperature')
    for d in adwin_drifts[:20]:
        ax.axvline(t[d], color='red', alpha=0.3, linewidth=0.5)
    ax.set_title('ADWIN Drift Detection (Temperature)', fontweight='bold')
    ax.legend()
    
    ax = axes[1]
    ax.plot(t, data[:, 1], 'b-', alpha=0.7, label='Pressure')
    for d in ph_drifts[:20]:
        ax.axvline(t[d], color='orange', alpha=0.3, linewidth=0.5)
    for cp in true_cps:
        ax.axvline(t[cp], color='green', linestyle=':', linewidth=2, label='True CP' if cp == true_cps[0] else '')
    ax.set_title('Page-Hinkley Drift Detection (Pressure)', fontweight='bold')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/drift_detection.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_feature_attribution(attributions, feature_names, labels):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    anomaly_mask = labels == 1
    normal_mask = labels == 0
    
    # Mean attribution for anomalies vs normal
    mean_attr_anomaly = np.mean(np.abs(attributions[anomaly_mask]), axis=0)
    mean_attr_normal = np.mean(np.abs(attributions[normal_mask]), axis=0)
    
    x = np.arange(len(feature_names))
    width = 0.35
    axes[0].barh(x - width/2, mean_attr_anomaly, width, label='Anomaly', color='red', alpha=0.7)
    axes[0].barh(x + width/2, mean_attr_normal, width, label='Normal', color='blue', alpha=0.7)
    axes[0].set_yticks(x)
    axes[0].set_yticklabels(feature_names, fontsize=9)
    axes[0].set_xlabel('Mean |Attribution|')
    axes[0].set_title('Feature Attribution: Anomaly vs Normal', fontweight='bold')
    axes[0].legend()
    
    # Heatmap of top anomalies
    top_anomalies = np.argsort(-np.max(np.abs(attributions), axis=1))[:20]
    im = axes[1].imshow(attributions[top_anomalies].T, aspect='auto', cmap='RdBu_r', interpolation='nearest')
    axes[1].set_yticks(range(len(feature_names)))
    axes[1].set_yticklabels(feature_names, fontsize=9)
    axes[1].set_xlabel('Top Anomaly Index')
    axes[1].set_title('Attribution Heatmap (Top 20 Anomalies)', fontweight='bold')
    plt.colorbar(im, ax=axes[1], shrink=0.8)
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/feature_attribution.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return mean_attr_anomaly, mean_attr_normal

def plot_streaming_pipeline():
    """Create a pipeline architecture diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(18, 8))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 8)
    ax.axis('off')
    ax.set_title('Streaming Anomaly Detection Pipeline Architecture', fontsize=14, fontweight='bold')
    
    boxes = [
        (1, 5.5, 'Data Ingestion\n(Kafka/Flink)', '#3498db'),
        (4.5, 5.5, 'Preprocessing\n& Validation', '#2ecc71'),
        (8, 5.5, 'Feature\nExtraction', '#e67e22'),
        (11.5, 5.5, 'Anomaly\nDetection', '#e74c3c'),
        (15, 5.5, 'Alert &\nExplanation', '#9b59b6'),
        (4.5, 2.5, 'Physics\nConstraints', '#1abc9c'),
        (8, 2.5, 'Drift\nDetection', '#f39c12'),
        (11.5, 2.5, 'Model\nRetraining', '#c0392b'),
        (15, 2.5, 'Dashboard\n& Storage', '#34495e'),
    ]
    
    for x, y, text, color in boxes:
        rect = plt.Rectangle((x-1.2, y-0.7), 2.4, 1.4, 
                             facecolor=color, alpha=0.8, edgecolor='black', linewidth=1.5,
                             zorder=2)
        ax.add_patch(rect)
        ax.text(x, y, text, ha='center', va='center', fontsize=9,
                fontweight='bold', color='white', zorder=3)
    
    # Arrows for main pipeline
    for x1, x2 in [(2.2, 3.3), (5.7, 6.8), (9.2, 10.3), (12.7, 13.8)]:
        ax.annotate('', xy=(x2, 5.5), xytext=(x1, 5.5),
                   arrowprops=dict(arrowstyle='->', color='black', lw=2))
    
    # Arrows for sub-pipeline
    ax.annotate('', xy=(4.5, 4.8), xytext=(4.5, 3.2),
               arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5, linestyle='--'))
    ax.annotate('', xy=(8, 4.8), xytext=(8, 3.2),
               arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5, linestyle='--'))
    ax.annotate('', xy=(11.5, 4.8), xytext=(11.5, 3.2),
               arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5, linestyle='--'))
    ax.annotate('', xy=(15, 4.8), xytext=(15, 3.2),
               arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5, linestyle='--'))
    
    # Bottom arrow
    for x1, x2 in [(5.7, 6.8), (9.2, 10.3), (12.7, 13.8)]:
        ax.annotate('', xy=(x2, 2.5), xytext=(x1, 2.5),
                   arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    
    # Labels
    ax.text(9, 7.5, 'Real-Time Streaming Layer (Apache Kafka + Apache Flink)', 
            ha='center', fontsize=12, fontstyle='italic', color='#2c3e50')
    ax.text(9, 1.2, 'Feedback & Adaptation Layer', 
            ha='center', fontsize=12, fontstyle='italic', color='#2c3e50')
    
    plt.savefig(f'{FIGDIR}/pipeline_architecture.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_confusion_matrices(labels, results):
    fig, axes = plt.subplots(1, 4, figsize=(18, 4))
    fig.suptitle('Confusion Matrices', fontsize=13, fontweight='bold')
    
    methods = [
        ('Isolation Forest', results['if_preds']),
        ('Deep SVDD-like', results['svdd_preds']),
        ('Physics', (results['phys_norm'] > np.quantile(results['phys_norm'], 0.95)).astype(int)),
        ('Combined', results['combined_preds'])
    ]
    
    for ax, (name, preds) in zip(axes, methods):
        cm = confusion_matrix(labels, preds)
        im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
        ax.set_title(name, fontsize=10, fontweight='bold')
        ax.set_xlabel('Predicted')
        ax.set_ylabel('Actual')
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(['Normal', 'Anomaly'], fontsize=8)
        ax.set_yticklabels(['Normal', 'Anomaly'], fontsize=8)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm[i, j]), ha='center', va='center',
                       fontsize=12, fontweight='bold', 
                       color='white' if cm[i, j] > cm.max()/2 else 'black')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/confusion_matrices.png', dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================
def main():
    print("=" * 70)
    print("EXPERIMENT: Automated Quality Control & Anomaly Detection")
    print("         for Large-Scale Scientific Data")
    print("=" * 70)
    
    # 1. Generate data
    print("\n[1/7] Generating synthetic scientific data...")
    t, data, labels, feature_names, true_cps = generate_scientific_data()
    print(f"  Data shape: {data.shape}")
    print(f"  Anomaly ratio: {labels.mean():.3f}")
    print(f"  True changepoints at: {true_cps}")
    plot_data_overview(t, data, labels, feature_names)
    print("  -> figures/data_overview.png saved")
    
    # 2. Changepoint detection
    print("\n[2/7] Running changepoint detection...")
    # PELT on pressure
    pelt_cps_pressure = run_pelt_detection(data[:, 1].reshape(-1, 1), pen=10)
    print(f"  PELT detected changepoints (Pressure): {pelt_cps_pressure}")
    pelt_eval = evaluate_changepoints(pelt_cps_pressure, true_cps)
    print(f"  PELT evaluation: {pelt_eval}")
    
    # BOCPD on pressure
    bocpd_prob = run_bocpd(data[:, 1])
    bocpd_cps = np.where(bocpd_prob > 0.3)[0].tolist()
    bocpd_eval = evaluate_changepoints(bocpd_cps, true_cps, tolerance=100)
    print(f"  BOCPD high-probability changepoints: {len(bocpd_cps)} detected")
    print(f"  BOCPD evaluation: {bocpd_eval}")
    
    plot_changepoint_results(t, data, feature_names, pelt_cps_pressure, bocpd_prob, true_cps)
    print("  -> figures/changepoint_detection.png saved")
    
    # 3. Multivariate anomaly detection
    print("\n[3/7] Running multivariate anomaly detection...")
    results = combined_anomaly_pipeline(data, feature_names, labels)
    
    for name, preds in [('IsolationForest', results['if_preds']),
                        ('DeepSVDD-like', results['svdd_preds']),
                        ('Combined', results['combined_preds'])]:
        p = precision_score(labels, preds, zero_division=0)
        r = recall_score(labels, preds, zero_division=0)
        f = f1_score(labels, preds, zero_division=0)
        auc = roc_auc_score(labels, results[f'{"if" if "Iso" in name else "svdd" if "SVDD" in name else "combined"}_{"norm" if name != "Combined" else ""}'.rstrip('_')])
        print(f"  {name}: P={p:.3f}, R={r:.3f}, F1={f:.3f}, AUC={auc:.3f}")
    
    plot_anomaly_scores(t, results, labels)
    print("  -> figures/anomaly_scores.png saved")
    
    plot_roc_pr_curves(labels, results)
    print("  -> figures/roc_pr_curves.png saved")
    
    plot_confusion_matrices(labels, results)
    print("  -> figures/confusion_matrices.png saved")
    
    # 4. Physics constraints
    print("\n[4/7] Evaluating physics constraint violations...")
    plot_physics_constraints(t, results['phys_details'])
    print("  -> figures/physics_constraints.png saved")
    
    phys_preds = (results['phys_norm'] > np.quantile(results['phys_norm'], 0.95)).astype(int)
    phys_p = precision_score(labels, phys_preds, zero_division=0)
    phys_r = recall_score(labels, phys_preds, zero_division=0)
    phys_f = f1_score(labels, phys_preds, zero_division=0)
    print(f"  Physics-only: P={phys_p:.3f}, R={phys_r:.3f}, F1={phys_f:.3f}")
    
    # 5. Drift detection
    print("\n[5/7] Running drift detection...")
    adwin_drifts = adwin_drift_detection(data[:, 0])  # Temperature
    ph_drifts = page_hinkley_test(data[:, 1])  # Pressure
    print(f"  ADWIN drifts detected (Temperature): {len(adwin_drifts)}")
    print(f"  Page-Hinkley drifts detected (Pressure): {len(ph_drifts)}")
    
    plot_drift_detection(t, data, adwin_drifts, ph_drifts, true_cps)
    print("  -> figures/drift_detection.png saved")
    
    # 6. Explainable anomaly detection
    print("\n[6/7] Computing feature attributions...")
    attributions, explanations = compute_feature_attribution(data, results['combined'], feature_names)
    mean_attr_a, mean_attr_n = plot_feature_attribution(attributions, feature_names, labels)
    print("  -> figures/feature_attribution.png saved")
    print(f"  Top 3 anomaly explanations:")
    for exp in explanations[:3]:
        feat_str = ", ".join([f"{f}:{v:.3f}" for f, v in exp['top_features']])
        print(f"    idx={exp['index']}, score={exp['score']:.3f} -> [{feat_str}]")
    
    # 7. Pipeline architecture
    print("\n[7/7] Generating streaming pipeline architecture...")
    plot_streaming_pipeline()
    print("  -> figures/pipeline_architecture.png saved")
    
    # ========================================================================
    # Summary Results
    # ========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY OF RESULTS")
    print("=" * 70)
    
    # Compute final metrics
    methods_results = {}
    for name, scores_key, preds_key in [
        ('Isolation Forest', 'if_norm', 'if_preds'),
        ('Deep SVDD-like', 'svdd_norm', 'svdd_preds'),
        ('Physics Constraints', 'phys_norm', None),
        ('Combined Pipeline', 'combined', 'combined_preds')
    ]:
        scores = results[scores_key]
        preds = results[preds_key] if preds_key else (scores > np.quantile(scores, 0.95)).astype(int)
        auc = roc_auc_score(labels, scores)
        ap = average_precision_score(labels, scores)
        p = precision_score(labels, preds, zero_division=0)
        r = recall_score(labels, preds, zero_division=0)
        f = f1_score(labels, preds, zero_division=0)
        methods_results[name] = {'AUC': auc, 'AP': ap, 'Precision': p, 'Recall': r, 'F1': f}
    
    print(f"\n{'Method':<22} {'AUC':>6} {'AP':>6} {'Prec':>6} {'Rec':>6} {'F1':>6}")
    print("-" * 56)
    for name, m in methods_results.items():
        print(f"{name:<22} {m['AUC']:>6.3f} {m['AP']:>6.3f} {m['Precision']:>6.3f} {m['Recall']:>6.3f} {m['F1']:>6.3f}")
    
    print(f"\nChangepoint Detection:")
    print(f"  PELT: {pelt_eval}")
    print(f"  BOCPD: {bocpd_eval}")
    
    print(f"\nDrift Detection:")
    print(f"  ADWIN points: {len(adwin_drifts)}")
    print(f"  Page-Hinkley points: {len(ph_drifts)}")
    
    # Save results as JSON
    output = {
        'data_info': {
            'n_samples': int(data.shape[0]),
            'n_features': int(data.shape[1]),
            'anomaly_ratio': float(labels.mean()),
            'true_changepoints': true_cps
        },
        'methods_results': methods_results,
        'changepoint_results': {
            'pelt': pelt_eval,
            'bocpd': bocpd_eval,
            'pelt_detected': pelt_cps_pressure
        },
        'drift_results': {
            'adwin_count': len(adwin_drifts),
            'page_hinkley_count': len(ph_drifts)
        },
        'threshold': float(results['threshold']),
        'feature_attribution': {
            'anomaly_mean': {fn: float(v) for fn, v in zip(feature_names, mean_attr_a)},
            'normal_mean': {fn: float(v) for fn, v in zip(feature_names, mean_attr_n)}
        }
    }
    
    with open('results.json', 'w') as f:
        json.dump(output, f, indent=2)
    print("\n  -> results.json saved")
    
    print("\n" + "=" * 70)
    print("Experiment completed successfully!")
    print(f"Generated figures: {os.listdir(FIGDIR)}")
    print("=" * 70)
    
    return output

if __name__ == '__main__':
    main()
