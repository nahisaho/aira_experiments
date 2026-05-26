"""
Consciousness Level Estimation and Disorder of Consciousness (DoC) Classification

Implements information-theoretic metrics for consciousness level estimation
and differential diagnosis of VS/UWS vs MCS.
"""

import numpy as np
from scipy.stats import entropy as scipy_entropy
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix


def compute_shannon_entropy(signal, n_bins=50):
    """Compute Shannon entropy of a signal."""
    hist, _ = np.histogram(signal, bins=n_bins, density=True)
    hist = hist[hist > 0]
    return scipy_entropy(hist, base=2)


def compute_spectral_entropy(signal, fs=256):
    """Compute spectral entropy of a signal."""
    fft_vals = np.abs(np.fft.rfft(signal))
    psd = fft_vals ** 2
    psd_norm = psd / psd.sum()
    psd_norm = psd_norm[psd_norm > 0]
    return scipy_entropy(psd_norm, base=2)


def compute_sample_entropy(signal, m=2, r_factor=0.2, max_len=1000):
    """Compute sample entropy of a time series."""
    signal = signal[:max_len]
    N = len(signal)
    r = r_factor * np.std(signal)
    
    if r == 0 or N < m + 1:
        return 0
    
    def count_matches(template_len):
        count = 0
        templates = np.array([signal[i:i+template_len] for i in range(N - template_len)])
        for i in range(len(templates)):
            for j in range(i + 1, len(templates)):
                if np.max(np.abs(templates[i] - templates[j])) < r:
                    count += 1
        return count
    
    A = count_matches(m + 1)
    B = count_matches(m)
    
    if B == 0:
        return 0
    
    return -np.log(A / B) if A > 0 else 0


def compute_permutation_entropy(signal, order=3, delay=1):
    """Compute permutation entropy."""
    n = len(signal)
    permutations = {}
    count = 0
    
    for i in range(n - (order - 1) * delay):
        indices = [i + j * delay for j in range(order)]
        values = [signal[idx] for idx in indices]
        perm = tuple(np.argsort(values))
        permutations[perm] = permutations.get(perm, 0) + 1
        count += 1
    
    if count == 0:
        return 0
    
    probs = np.array(list(permutations.values())) / count
    return scipy_entropy(probs, base=2)


def compute_lempel_ziv_complexity(signal, threshold=None):
    """Compute normalized Lempel-Ziv complexity."""
    if threshold is None:
        threshold = np.median(signal)
    
    binary = (signal > threshold).astype(int)
    s = ''.join(map(str, binary))
    n = len(s)
    
    if n <= 1:
        return 0
    
    complexity = 1
    i = 0
    k = 1
    
    while i + k < n:
        substr = s[i+1:i+k+1]
        if substr in s[:i+k]:
            k += 1
        else:
            complexity += 1
            i += k
            k = 1
    
    norm = n / np.log2(n)
    return complexity / norm


def extract_consciousness_features(eeg_data, fs=256):
    """
    Extract information-theoretic features from multi-channel EEG.
    
    Parameters:
    -----------
    eeg_data: ndarray (n_channels, n_timepoints)
    fs: sampling frequency
    
    Returns:
    --------
    features: dict of feature names and values
    """
    n_channels = eeg_data.shape[0]
    
    features = {}
    
    # Per-channel metrics
    shannon_vals = [compute_shannon_entropy(eeg_data[ch]) for ch in range(n_channels)]
    spectral_vals = [compute_spectral_entropy(eeg_data[ch], fs) for ch in range(n_channels)]
    perm_vals = [compute_permutation_entropy(eeg_data[ch]) for ch in range(n_channels)]
    lzc_vals = [compute_lempel_ziv_complexity(eeg_data[ch]) for ch in range(n_channels)]
    
    features['shannon_entropy_mean'] = np.mean(shannon_vals)
    features['shannon_entropy_std'] = np.std(shannon_vals)
    features['spectral_entropy_mean'] = np.mean(spectral_vals)
    features['spectral_entropy_std'] = np.std(spectral_vals)
    features['permutation_entropy_mean'] = np.mean(perm_vals)
    features['permutation_entropy_std'] = np.std(perm_vals)
    features['lzc_mean'] = np.mean(lzc_vals)
    features['lzc_std'] = np.std(lzc_vals)
    
    # Cross-channel connectivity
    corr_matrix = np.corrcoef(eeg_data)
    upper_triangle = corr_matrix[np.triu_indices(n_channels, k=1)]
    features['mean_connectivity'] = np.mean(np.abs(upper_triangle))
    features['connectivity_variance'] = np.var(upper_triangle)
    
    # Power spectral features
    for ch in range(min(n_channels, 4)):
        fft_vals = np.abs(np.fft.rfft(eeg_data[ch]))
        freqs = np.fft.rfftfreq(eeg_data.shape[1], 1.0/fs)
        
        delta_power = np.sum(fft_vals[(freqs >= 0.5) & (freqs < 4)] ** 2)
        theta_power = np.sum(fft_vals[(freqs >= 4) & (freqs < 8)] ** 2)
        alpha_power = np.sum(fft_vals[(freqs >= 8) & (freqs < 13)] ** 2)
        beta_power = np.sum(fft_vals[(freqs >= 13) & (freqs < 30)] ** 2)
        total_power = delta_power + theta_power + alpha_power + beta_power
        
        if total_power > 0:
            features[f'delta_ratio_ch{ch}'] = delta_power / total_power
            features[f'alpha_ratio_ch{ch}'] = alpha_power / total_power
    
    return features


def simulate_doc_dataset(n_subjects_per_class=30, n_channels=16, n_timepoints=1000, fs=256):
    """
    Simulate EEG dataset for disorders of consciousness classification.
    
    Classes: 0=VS/UWS, 1=MCS, 2=Healthy
    """
    from src.pci_simulation import simulate_neural_mass_model
    
    X_list = []
    y_list = []
    
    class_map = {0: 'vegetative', 1: 'mcs', 2: 'awake'}
    
    for class_label, condition in class_map.items():
        for subj in range(n_subjects_per_class):
            np.random.seed(class_label * 1000 + subj)
            E, _ = simulate_neural_mass_model(n_channels, n_timepoints, consciousness_level=condition)
            noise = np.random.randn(n_channels, n_timepoints) * 0.05
            eeg = E + noise
            
            features = extract_consciousness_features(eeg, fs)
            X_list.append(list(features.values()))
            y_list.append(class_label)
    
    feature_names = list(features.keys())
    X = np.array(X_list)
    y = np.array(y_list)
    
    return X, y, feature_names


def classify_consciousness_states(X, y, feature_names=None):
    """
    Train and evaluate classifiers for consciousness state classification.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    results = {}
    
    # SVM classifier
    svm = SVC(kernel='rbf', C=10, gamma='scale', random_state=42)
    svm_scores = cross_val_score(svm, X_scaled, y, cv=5, scoring='accuracy')
    results['SVM'] = {'mean_accuracy': svm_scores.mean(), 'std_accuracy': svm_scores.std(), 'scores': svm_scores}
    
    # Random Forest
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_scores = cross_val_score(rf, X_scaled, y, cv=5, scoring='accuracy')
    results['RandomForest'] = {'mean_accuracy': rf_scores.mean(), 'std_accuracy': rf_scores.std(), 'scores': rf_scores}
    
    # Feature importance from RF
    rf.fit(X_scaled, y)
    if feature_names:
        importance = dict(zip(feature_names, rf.feature_importances_))
        results['feature_importance'] = importance
    
    # Full model evaluation
    svm.fit(X_scaled, y)
    y_pred = svm.predict(X_scaled)
    results['confusion_matrix'] = confusion_matrix(y, y_pred)
    results['classification_report'] = classification_report(
        y, y_pred, target_names=['VS/UWS', 'MCS', 'Healthy'], output_dict=True
    )
    
    return results
