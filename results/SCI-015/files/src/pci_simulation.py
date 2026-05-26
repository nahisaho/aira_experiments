"""
Perturbational Complexity Index (PCI) Simulation Module

Implements PCI computation based on Lempel-Ziv complexity of EEG responses
to simulated transcranial magnetic stimulation (TMS).
"""

import numpy as np
from scipy.signal import hilbert


def lempel_ziv_complexity(binary_sequence):
    """
    Compute Lempel-Ziv complexity of a binary sequence.
    Based on Kaspar & Schuster (1987) algorithm.
    """
    s = ''.join(map(str, binary_sequence.astype(int)))
    n = len(s)
    
    if n == 0:
        return 0
    
    c = 1
    l = 1
    i = 0
    k = 1
    k_max = 1
    
    while i + k <= n:
        substring = s[i + 1:i + k + 1] if i + k + 1 <= n else s[i + 1:]
        if k <= len(substring) and s[i + k - 1 + 1:i + k + 1] if False else True:
            pass
        
        sub = s[i:i + k]
        history = s[0:i + k - 1]
        
        if sub in history or k == 1:
            k += 1
            if k > k_max:
                k_max = k
        else:
            c += 1
            i += k_max
            k = 1
            k_max = 1
        
        if i + k > n:
            break
    
    # Normalize
    b = n / np.log2(n) if n > 1 else 1
    return c / b


def lempel_ziv_complexity_v2(binary_string):
    """Improved LZ complexity using the incremental parsing algorithm."""
    s = binary_string
    n = len(s)
    if n <= 1:
        return n
    
    complexity = 1
    prefix_len = 1
    
    i = 0
    while i + prefix_len < n:
        found = False
        for l in range(1, n - i - prefix_len + 1):
            substr = s[i + prefix_len:i + prefix_len + l]
            if substr in s[i:i + prefix_len + l - 1]:
                found = True
                continue
            else:
                complexity += 1
                i = i + prefix_len
                prefix_len = l
                found = False
                break
        if found:
            break
    
    norm = n / np.log2(n) if n > 1 else 1
    return complexity / norm


def simulate_neural_mass_model(n_channels, n_timepoints, coupling_strength=0.3, 
                                noise_level=0.1, consciousness_level='awake'):
    """
    Simulate a neural mass model with Wilson-Cowan dynamics.
    
    Parameters:
    -----------
    consciousness_level: str
        'awake', 'light_sedation', 'deep_anesthesia', 'vegetative', 'mcs'
    """
    params = {
        'awake': {'coupling': coupling_strength, 'noise': noise_level, 'inhibition': 0.3},
        'light_sedation': {'coupling': coupling_strength * 0.7, 'noise': noise_level * 0.7, 'inhibition': 0.5},
        'deep_anesthesia': {'coupling': coupling_strength * 0.3, 'noise': noise_level * 0.3, 'inhibition': 0.8},
        'vegetative': {'coupling': coupling_strength * 0.2, 'noise': noise_level * 0.2, 'inhibition': 0.9},
        'mcs': {'coupling': coupling_strength * 0.5, 'noise': noise_level * 0.5, 'inhibition': 0.6},
    }
    
    p = params[consciousness_level]
    dt = 0.001
    tau_e = 0.01
    tau_i = 0.02
    
    E = np.zeros((n_channels, n_timepoints))
    I = np.zeros((n_channels, n_timepoints))
    
    W_ee = np.random.randn(n_channels, n_channels) * p['coupling']
    np.fill_diagonal(W_ee, 0)
    W_ei = np.random.randn(n_channels, n_channels) * p['inhibition']
    
    E[:, 0] = np.random.rand(n_channels) * 0.5
    I[:, 0] = np.random.rand(n_channels) * 0.3
    
    def sigmoid(x, threshold=0.5, steepness=5):
        return 1 / (1 + np.exp(-steepness * (x - threshold)))
    
    for t in range(1, n_timepoints):
        noise_e = np.random.randn(n_channels) * p['noise']
        noise_i = np.random.randn(n_channels) * p['noise'] * 0.5
        
        dE = (-E[:, t-1] + sigmoid(W_ee @ E[:, t-1] - W_ei @ I[:, t-1] + noise_e)) / tau_e * dt
        dI = (-I[:, t-1] + sigmoid(E[:, t-1] + noise_i)) / tau_i * dt
        
        E[:, t] = np.clip(E[:, t-1] + dE, 0, 1)
        I[:, t] = np.clip(I[:, t-1] + dI, 0, 1)
    
    return E, I


def apply_tms_perturbation(E, I, stim_channel, stim_time, stim_strength=0.8):
    """Apply a simulated TMS pulse perturbation."""
    E_perturbed = E.copy()
    I_perturbed = I.copy()
    
    E_perturbed[stim_channel, stim_time:stim_time+5] += stim_strength
    E_perturbed = np.clip(E_perturbed, 0, 1)
    
    return E_perturbed


def compute_pci(eeg_response, pre_stim_samples=100):
    """
    Compute Perturbational Complexity Index from EEG response.
    
    PCI = LZ complexity of the spatiotemporal binary matrix of significant responses.
    """
    n_channels, n_timepoints = eeg_response.shape
    
    baseline = eeg_response[:, :pre_stim_samples]
    baseline_mean = baseline.mean(axis=1, keepdims=True)
    baseline_std = baseline.std(axis=1, keepdims=True)
    baseline_std[baseline_std < 1e-10] = 1e-10
    
    z_scores = (eeg_response - baseline_mean) / baseline_std
    binary_matrix = (np.abs(z_scores[:, pre_stim_samples:]) > 2.0).astype(int)
    
    binary_sequence = binary_matrix.flatten()
    
    lzc = lempel_ziv_complexity_v2(''.join(map(str, binary_sequence)))
    
    source_entropy = np.mean(binary_matrix)
    if source_entropy > 0 and source_entropy < 1:
        normalization = -source_entropy * np.log2(source_entropy) - (1-source_entropy) * np.log2(1-source_entropy)
    else:
        normalization = 1.0
    
    pci = lzc * normalization if normalization > 0 else 0
    
    return pci, binary_matrix


def compute_pci_across_conditions(n_channels=32, n_timepoints=500, n_trials=10):
    """Run PCI computation across multiple consciousness conditions."""
    conditions = ['awake', 'light_sedation', 'deep_anesthesia', 'vegetative', 'mcs']
    results = {}
    
    for cond in conditions:
        pci_values = []
        for trial in range(n_trials):
            np.random.seed(trial * 100 + hash(cond) % 1000)
            E, I = simulate_neural_mass_model(n_channels, n_timepoints, consciousness_level=cond)
            
            stim_ch = n_channels // 2
            stim_t = 100
            E_stim = apply_tms_perturbation(E, I, stim_ch, stim_t)
            
            response = E_stim - E
            pci, _ = compute_pci(E_stim, pre_stim_samples=stim_t)
            pci_values.append(pci)
        
        results[cond] = {
            'mean': np.mean(pci_values),
            'std': np.std(pci_values),
            'values': pci_values
        }
    
    return results
