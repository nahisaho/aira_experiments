"""
Integrated Information Theory (IIT) - Φ (Phi) Computation Module

Efficient algorithms for computing integrated information in neural network models.
Supports both exact and approximate methods for scalability.
"""

import numpy as np
from itertools import combinations
from scipy.stats import entropy as scipy_entropy


def compute_tpm(data, n_nodes, tau=1):
    """Compute Transition Probability Matrix from time-series data."""
    n_states = 2 ** n_nodes
    tpm = np.zeros((n_states, n_states))
    
    for t in range(len(data) - tau):
        current = int(''.join(map(str, data[t].astype(int))), 2)
        future = int(''.join(map(str, data[t + tau].astype(int))), 2)
        tpm[current, future] += 1
    
    row_sums = tpm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    tpm = tpm / row_sums
    return tpm


def mutual_information(p_xy):
    """Compute mutual information from joint distribution."""
    p_x = p_xy.sum(axis=1)
    p_y = p_xy.sum(axis=0)
    
    mi = 0.0
    for i in range(p_xy.shape[0]):
        for j in range(p_xy.shape[1]):
            if p_xy[i, j] > 1e-12 and p_x[i] > 1e-12 and p_y[j] > 1e-12:
                mi += p_xy[i, j] * np.log2(p_xy[i, j] / (p_x[i] * p_y[j]))
    return mi


def effective_information(tpm):
    """Compute effective information (EI) of a system's TPM."""
    n_states = tpm.shape[0]
    uniform = np.ones(n_states) / n_states
    
    ei = 0.0
    for i in range(n_states):
        for j in range(n_states):
            if tpm[i, j] > 1e-12:
                marginal = uniform @ tpm[:, j]
                if marginal > 1e-12:
                    ei += uniform[i] * tpm[i, j] * np.log2(tpm[i, j] / marginal)
    return ei


def partition_system(n_nodes):
    """Generate all bipartitions of a system (for MIP search)."""
    partitions = []
    nodes = list(range(n_nodes))
    for r in range(1, n_nodes):
        for combo in combinations(nodes, r):
            part_a = list(combo)
            part_b = [n for n in nodes if n not in combo]
            if len(part_a) <= len(part_b):
                partitions.append((part_a, part_b))
    return partitions


def compute_stochastic_interaction(data, n_nodes):
    """
    Compute stochastic interaction (Φ-SI) from multivariate time series.
    Φ_SI = H(whole) - Σ H(parts) — measures how much the whole exceeds
    the sum of its parts in terms of entropy.
    """
    n_samples = data.shape[0]
    
    # Whole system entropy (using covariance-based Gaussian approximation)
    cov_whole = np.cov(data.T) + np.eye(n_nodes) * 1e-6
    sign, logdet_whole = np.linalg.slogdet(cov_whole)
    h_whole = 0.5 * logdet_whole + 0.5 * n_nodes * np.log(2 * np.pi * np.e)
    
    # Sum of individual entropies
    h_parts = 0.0
    for i in range(n_nodes):
        var_i = np.var(data[:, i]) + 1e-6
        h_parts += 0.5 * np.log(var_i) + 0.5 * np.log(2 * np.pi * np.e)
    
    # Stochastic interaction: redundancy measure
    phi_si = h_parts - h_whole
    return max(0, phi_si)


def compute_phi_geometric(data, n_nodes):
    """
    Compute geometric integrated information (Φ_G).
    Based on Oizumi et al. (2016) geometric approach using KL divergence.
    """
    cov = np.cov(data.T) + np.eye(n_nodes) * 1e-6
    
    partitions = partition_system(n_nodes)
    if not partitions:
        return 0.0, None
    
    min_phi = float('inf')
    mip = None
    
    for part_a, part_b in partitions:
        # Disconnected covariance (block diagonal)
        cov_disconnected = np.zeros_like(cov)
        idx_a = np.array(part_a)
        idx_b = np.array(part_b)
        cov_disconnected[np.ix_(idx_a, idx_a)] = cov[np.ix_(idx_a, idx_a)]
        cov_disconnected[np.ix_(idx_b, idx_b)] = cov[np.ix_(idx_b, idx_b)]
        
        # KL divergence between whole and disconnected (Gaussian)
        try:
            cov_disc_inv = np.linalg.inv(cov_disconnected + np.eye(n_nodes) * 1e-6)
            sign_d, logdet_disc = np.linalg.slogdet(cov_disconnected + np.eye(n_nodes) * 1e-6)
            sign_w, logdet_whole = np.linalg.slogdet(cov)
            
            kl = 0.5 * (np.trace(cov_disc_inv @ cov) - n_nodes + logdet_disc - logdet_whole)
            kl = max(0, kl)
        except np.linalg.LinAlgError:
            kl = float('inf')
        
        if kl < min_phi:
            min_phi = kl
            mip = (part_a, part_b)
    
    return min_phi, mip


def compute_phi_approximate(connectivity_matrix, n_samples=1000, noise_level=0.1):
    """
    Approximate Φ computation using stochastic sampling.
    Uses geometric Φ (Φ_G) for more robust measurement.
    """
    n_nodes = connectivity_matrix.shape[0]
    
    # Simulate dynamics
    data = np.zeros((n_samples, n_nodes))
    data[0] = np.random.randn(n_nodes) * 0.5
    
    for t in range(1, n_samples):
        data[t] = 0.5 * np.tanh(connectivity_matrix @ data[t-1]) + np.random.randn(n_nodes) * noise_level
    
    # Compute multiple Φ measures
    phi_g, mip = compute_phi_geometric(data, n_nodes)
    phi_si = compute_stochastic_interaction(data, n_nodes)
    
    # Use geometric Φ as primary measure
    phi = phi_g
    
    return phi, mip, data, phi_si


def generate_network(n_nodes, connectivity_type='integrated'):
    """Generate connectivity matrices for different network types."""
    if connectivity_type == 'integrated':
        W = np.random.randn(n_nodes, n_nodes) * 0.5
        W += np.eye(n_nodes) * 0.8
        W = (W + W.T) / 2
    elif connectivity_type == 'modular':
        W = np.zeros((n_nodes, n_nodes))
        half = n_nodes // 2
        W[:half, :half] = np.random.randn(half, half) * 0.4 + np.eye(half) * 0.7
        W[half:, half:] = np.random.randn(n_nodes - half, n_nodes - half) * 0.4 + np.eye(n_nodes - half) * 0.7
        W[:half, half:] = np.random.randn(half, n_nodes - half) * 0.02
        W[half:, :half] = np.random.randn(n_nodes - half, half) * 0.02
    elif connectivity_type == 'feedforward':
        W = np.zeros((n_nodes, n_nodes))
        for i in range(n_nodes - 1):
            W[i+1, i] = np.random.randn() * 0.4 + 0.6
    elif connectivity_type == 'disconnected':
        W = np.eye(n_nodes) * 0.5
    else:
        raise ValueError(f"Unknown connectivity type: {connectivity_type}")
    
    return W
