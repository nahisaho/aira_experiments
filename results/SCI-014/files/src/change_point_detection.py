"""
Change point detection algorithms for longitudinal mHealth data.
Implements CUSUM, PELT-inspired, and Bayesian approaches.
"""

import numpy as np
import pandas as pd
from scipy import stats


def cusum_detection(signal, threshold=1.0, drift=0.0):
    """Cumulative Sum (CUSUM) change point detection.
    
    Args:
        signal: 1D array of observations
        threshold: detection threshold (h)
        drift: allowance parameter (k)
    
    Returns:
        change_points: list of detected change point indices
        cusum_pos: positive CUSUM statistics
        cusum_neg: negative CUSUM statistics
    """
    n = len(signal)
    mean_val = np.mean(signal[:min(10, n)])
    
    cusum_pos = np.zeros(n)
    cusum_neg = np.zeros(n)
    change_points = []
    
    for i in range(1, n):
        cusum_pos[i] = max(0, cusum_pos[i-1] + (signal[i] - mean_val) - drift)
        cusum_neg[i] = max(0, cusum_neg[i-1] - (signal[i] - mean_val) - drift)
        
        if cusum_pos[i] > threshold or cusum_neg[i] > threshold:
            change_points.append(i)
            cusum_pos[i] = 0
            cusum_neg[i] = 0
            mean_val = signal[i]
    
    return change_points, cusum_pos, cusum_neg


def pelt_detection(signal, penalty=1.0, min_segment=5):
    """Pruned Exact Linear Time (PELT) inspired change point detection.
    
    Uses a simplified dynamic programming approach with L2 cost.
    """
    n = len(signal)
    
    def segment_cost(start, end):
        seg = signal[start:end]
        if len(seg) < 2:
            return 0
        return np.sum((seg - np.mean(seg))**2)
    
    # Dynamic programming
    F = np.full(n + 1, np.inf)
    F[0] = -penalty
    cp_sets = [[] for _ in range(n + 1)]
    
    for t in range(min_segment, n + 1):
        candidates = range(max(0, t - 50), t - min_segment + 1)
        for s in candidates:
            cost = F[s] + segment_cost(s, t) + penalty
            if cost < F[t]:
                F[t] = cost
                cp_sets[t] = cp_sets[s] + [s] if s > 0 else []
    
    return cp_sets[n]


def bayesian_online_cpd(signal, hazard_rate=1/50, observation_var=1.0):
    """Bayesian Online Change Point Detection (Adams & MacKay, 2007).
    
    Args:
        signal: 1D observation array
        hazard_rate: prior probability of change at each step
        observation_var: observation noise variance
    
    Returns:
        change_points: detected change points
        run_length_probs: run length probability matrix
    """
    n = len(signal)
    R = np.zeros((n + 1, n + 1))
    R[0, 0] = 1.0
    
    mu_prior = np.mean(signal[:min(10, n)])
    var_prior = np.var(signal[:min(10, n)]) + 1e-6
    
    change_points = []
    max_run_lengths = np.zeros(n)
    
    for t in range(n):
        # Predictive probability under each run length
        pred_probs = np.zeros(t + 1)
        for r in range(t + 1):
            if R[r, t] > 1e-10:
                start = max(0, t - r)
                segment = signal[start:t+1]
                mu = np.mean(segment) if len(segment) > 0 else mu_prior
                var = max(np.var(segment), observation_var) if len(segment) > 1 else var_prior
                pred_probs[r] = stats.norm.pdf(signal[t], mu, np.sqrt(var))
        
        # Growth probabilities
        growth = R[:t+1, t] * pred_probs * (1 - hazard_rate)
        
        # Change point probability
        cp_prob = np.sum(R[:t+1, t] * pred_probs * hazard_rate)
        
        # Update run length distribution
        R[1:t+2, t+1] = growth
        R[0, t+1] = cp_prob
        
        # Normalize
        total = np.sum(R[:t+2, t+1])
        if total > 0:
            R[:t+2, t+1] /= total
        
        max_run_lengths[t] = np.argmax(R[:t+2, t+1])
        
        # Detect change point when run length drops significantly
        if t > 0 and max_run_lengths[t] < max_run_lengths[t-1] * 0.3 and max_run_lengths[t-1] > 3:
            change_points.append(t)
    
    return change_points, R, max_run_lengths


def multimodal_cpd(signals_dict, method='cusum', **kwargs):
    """Apply change point detection across multiple modalities and fuse results.
    
    Args:
        signals_dict: dict of {modality_name: signal_array}
        method: 'cusum', 'pelt', or 'bayesian'
    
    Returns:
        fused_change_points: consensus change points
        per_modality_cps: dict of change points per modality
    """
    per_modality_cps = {}
    
    for name, signal in signals_dict.items():
        if method == 'cusum':
            cps, _, _ = cusum_detection(signal, **kwargs)
        elif method == 'pelt':
            cps = pelt_detection(signal, **kwargs)
        elif method == 'bayesian':
            cps, _, _ = bayesian_online_cpd(signal, **kwargs)
        else:
            raise ValueError(f"Unknown method: {method}")
        per_modality_cps[name] = cps
    
    # Fuse: find time points where at least 2 modalities agree (within ±2 window)
    all_cps = []
    for cps in per_modality_cps.values():
        all_cps.extend(cps)
    
    if not all_cps:
        return [], per_modality_cps
    
    all_cps = sorted(set(all_cps))
    
    fused = []
    for cp in all_cps:
        agreement = 0
        for name, cps in per_modality_cps.items():
            if any(abs(cp - c) <= 2 for c in cps):
                agreement += 1
        if agreement >= 2:
            fused.append(cp)
    
    # Merge nearby fused points
    if fused:
        merged = [fused[0]]
        for cp in fused[1:]:
            if cp - merged[-1] > 3:
                merged.append(cp)
        fused = merged
    
    return fused, per_modality_cps


def evaluate_cpd(detected_cps, true_cps, tolerance=3):
    """Evaluate change point detection performance.
    
    Args:
        detected_cps: list of detected change points
        true_cps: list of true change points
        tolerance: window for matching (±tolerance)
    
    Returns:
        dict with precision, recall, F1, mean_delay
    """
    if not true_cps:
        return {'precision': 1.0 if not detected_cps else 0.0,
                'recall': 1.0,
                'f1': 1.0 if not detected_cps else 0.0,
                'mean_delay': 0.0}
    
    if not detected_cps:
        return {'precision': 0.0, 'recall': 0.0, 'f1': 0.0, 'mean_delay': float('inf')}
    
    # Match detected to true
    tp = 0
    delays = []
    matched_true = set()
    
    for d in detected_cps:
        for i, t in enumerate(true_cps):
            if i not in matched_true and abs(d - t) <= tolerance:
                tp += 1
                delays.append(d - t)
                matched_true.add(i)
                break
    
    precision = tp / len(detected_cps) if detected_cps else 0
    recall = tp / len(true_cps) if true_cps else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    mean_delay = np.mean(delays) if delays else float('inf')
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'mean_delay': mean_delay,
    }
