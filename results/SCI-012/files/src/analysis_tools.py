"""
Analysis Tools: Firing rate, phase synchrony, and information transfer.
"""
import numpy as np
from scipy import signal


def compute_firing_rates(spike_trains, N_neurons, T, bin_size=50):
    """Compute population firing rates in time bins."""
    n_bins = int(T / bin_size)
    rates = np.zeros(n_bins)
    
    for t_spike, _ in spike_trains:
        bin_idx = min(int(t_spike / bin_size), n_bins - 1)
        rates[bin_idx] += 1
    
    rates = rates / N_neurons / (bin_size / 1000)  # Hz
    return rates, np.arange(n_bins) * bin_size


def compute_cv_isi(spike_trains, N_neurons, T):
    """Compute coefficient of variation of inter-spike intervals."""
    neuron_spikes = {}
    for t, idx in spike_trains:
        neuron_spikes.setdefault(idx, []).append(t)
    
    cvs = []
    for idx, times in neuron_spikes.items():
        if len(times) > 2:
            isis = np.diff(sorted(times))
            if isis.mean() > 0:
                cvs.append(isis.std() / isis.mean())
    
    return np.mean(cvs) if cvs else 0.0, np.std(cvs) if cvs else 0.0


def compute_phase_synchrony(pop_rates_1, pop_rates_2, dt):
    """Compute phase synchrony between two population rate signals."""
    analytic_1 = signal.hilbert(pop_rates_1 - pop_rates_1.mean())
    analytic_2 = signal.hilbert(pop_rates_2 - pop_rates_2.mean())
    
    phase_1 = np.angle(analytic_1)
    phase_2 = np.angle(analytic_2)
    
    phase_diff = phase_1 - phase_2
    pli = np.abs(np.mean(np.sign(np.sin(phase_diff))))  # Phase Lag Index
    plv = np.abs(np.mean(np.exp(1j * phase_diff)))  # Phase Locking Value
    
    return {'PLI': pli, 'PLV': plv, 'mean_phase_diff': np.mean(phase_diff)}


def compute_transfer_entropy(source, target, k=1, bins=10):
    """Compute transfer entropy from source to target time series."""
    source_d = np.digitize(source, np.linspace(source.min(), source.max(), bins))
    target_d = np.digitize(target, np.linspace(target.min(), target.max(), bins))
    
    n = len(source_d) - k
    if n <= 0:
        return 0.0
    
    # Estimate joint and conditional probabilities
    joint_counts = {}
    cond_counts = {}
    target_counts = {}
    
    for t in range(k, len(source_d)):
        target_future = target_d[t]
        target_past = tuple(target_d[t-k:t])
        source_past = tuple(source_d[t-k:t])
        
        key_joint = (target_future, target_past, source_past)
        key_cond = (target_past, source_past)
        key_target = (target_future, target_past)
        
        joint_counts[key_joint] = joint_counts.get(key_joint, 0) + 1
        cond_counts[key_cond] = cond_counts.get(key_cond, 0) + 1
        target_counts[key_target] = target_counts.get(key_target, 0) + 1
    
    te = 0.0
    for key, count in joint_counts.items():
        tf, tp, sp = key
        p_joint = count / n
        p_cond = cond_counts.get((tp, sp), 1) / n
        p_target = target_counts.get((tf, tp), 1) / n
        tp_only = sum(v for k, v in cond_counts.items() if k[0] == tp) / n
        
        if p_cond > 0 and tp_only > 0:
            te += p_joint * np.log2(
                (p_joint * tp_only) / (p_target * p_cond + 1e-12) + 1e-12
            )
    
    return max(te, 0.0)


def compute_power_spectrum(rate_signal, dt):
    """Compute power spectrum of firing rate signal."""
    fs = 1000.0 / dt
    freqs, psd = signal.welch(rate_signal, fs=fs, nperseg=min(256, len(rate_signal)//2))
    return freqs, psd


def compute_fano_factor(spike_trains, N_neurons, T, window=100):
    """Compute Fano factor (variance/mean of spike counts)."""
    n_windows = int(T / window)
    counts = np.zeros((N_neurons, n_windows))
    
    for t, idx in spike_trains:
        w = min(int(t / window), n_windows - 1)
        if idx < N_neurons:
            counts[idx, w] += 1
    
    means = counts.mean(axis=1)
    variances = counts.var(axis=1)
    
    active = means > 0
    if active.sum() == 0:
        return 1.0
    
    return np.mean(variances[active] / means[active])
