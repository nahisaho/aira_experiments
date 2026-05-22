"""
BB84 / E91 Finite-Key Analysis
================================
Finite-key security bounds following:
  - Scarani & Renner (2008) PRL 100, 200501
  - Tomamichel et al. (2012) Nature Comm 3, 634
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import rel_entr
import json, os, time

np.random.seed(42)

# ─── helpers ──────────────────────────────────────────────────────────────────
def h(p):
    """Binary entropy."""
    p = np.asarray(p, dtype=float)
    out = np.zeros_like(p)
    mask = (p > 0) & (p < 1)
    pm = p[mask]
    out[mask] = -pm * np.log2(pm) - (1 - pm) * np.log2(1 - pm)
    return out


def h2(p):
    """Binary entropy (scalar-safe)."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


# ─── BB84 finite-key ──────────────────────────────────────────────────────────
def bb84_key_length(n, q, epsilon=1e-8):
    """
    Finite-key secret key length for BB84 (simple Shor-Preskill + finite-sample correction).

    Parameters
    ----------
    n       : int   – number of raw sifted bits
    q       : float – quantum bit error rate (QBER)
    epsilon : float – security parameter

    Returns
    -------
    l : float – secret key bits (≥0)
    """
    # Privacy amplification correction (composite epsilon)
    leak_ec = n * h2(q) + (np.log2(1 / epsilon)) * np.sqrt(n)  # error-correction leakage
    # Statistical correction on phase error estimate (Chernoff bound)
    delta = np.sqrt((2 * np.log(1 / epsilon) + np.log(2 * n + 1)) / (2 * n))
    q_phase = min(q + delta, 0.5)
    # Key length
    l = n * (1 - h2(q_phase)) - leak_ec - 2 * np.log2(1 / epsilon)
    return max(l, 0.0)


def bb84_key_rate(n_values, qber=0.03, epsilon=1e-8):
    lengths = [bb84_key_length(n, qber, epsilon) / n for n in n_values]
    return np.array(lengths)


# ─── E91 finite-key ──────────────────────────────────────────────────────────
def e91_chsh_violation(noise_param):
    """
    CHSH value S for Werner state rho = (1-p)|Phi+><Phi+| + p/4 * I.
    S = 2*sqrt(2)*(1-p)
    """
    return 2 * np.sqrt(2) * (1 - noise_param)


def e91_key_rate_from_chsh(S, n, epsilon=1e-8):
    """
    Secret key rate from CHSH inequality violation S ∈ (2, 2√2].
    Using Acin et al. device-independent bound approximation.
    """
    # Tolerated QBER from S: q ≤ (1 - sqrt(S^2/4 - 1))/2 approximately
    # Using Pironio et al. bound: r ≥ 1 - h(e) where e derived from S
    if S <= 2.0:
        return 0.0
    # Critical noise threshold
    eta = (S - 2) / (2 * np.sqrt(2) - 2)  # 0→no violation, 1→pure Bell state
    q_eff = 0.5 * (1 - eta)
    # Finite-key correction
    delta = np.sqrt(np.log(1 / epsilon) / (2 * n))
    q_eff = min(q_eff + delta, 0.5)
    return max(1 - h2(q_eff), 0.0)


# ─── simulation sweep ─────────────────────────────────────────────────────────
n_values = np.logspace(3, 8, 60).astype(int)
qber_levels = [0.01, 0.03, 0.05, 0.08, 0.11]

# BB84 key-rate vs block size for different QBER
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax1, ax2 = axes

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(qber_levels)))
for q, c in zip(qber_levels, colors):
    rates = bb84_key_rate(n_values, qber=q)
    ax1.semilogx(n_values, rates, color=c, linewidth=2, label=f'QBER={q:.2f}')

asymptotic = [1 - h2(q) for q in qber_levels]
for q, a, c in zip(qber_levels, asymptotic, colors):
    ax1.axhline(a, color=c, linestyle='--', alpha=0.4)

ax1.set_xlabel('Block size n (raw sifted bits)', fontsize=12)
ax1.set_ylabel('Secret key rate (bits per sifted bit)', fontsize=12)
ax1.set_title('BB84 Finite-Key Rate vs Block Size', fontsize=13)
ax1.legend(fontsize=9)
ax1.set_ylim(0, 1)
ax1.grid(True, alpha=0.3)
ax1.fill_between([n_values.min(), n_values.max()], 0, 0, alpha=0)

# E91 CHSH vs noise, key rate
noise_params = np.linspace(0, 0.45, 200)
S_values = np.array([e91_chsh_violation(p) for p in noise_params])
rates_e91_large = np.array([e91_key_rate_from_chsh(S, n=int(1e6)) for S in S_values])
rates_e91_small = np.array([e91_key_rate_from_chsh(S, n=int(1e4)) for S in S_values])

ax2.plot(noise_params, S_values / (2 * np.sqrt(2)), 'b-', linewidth=2, label='CHSH S/(2√2)')
ax2.plot(noise_params, rates_e91_large, 'r-', linewidth=2, label='E91 key rate (n=10⁶)')
ax2.plot(noise_params, rates_e91_small, 'r--', linewidth=2, label='E91 key rate (n=10⁴)')
ax2.axvline(1 - 1 / np.sqrt(2), color='gray', linestyle=':', label='Classical CHSH threshold')
ax2.set_xlabel('Werner state noise parameter p', fontsize=12)
ax2.set_ylabel('Normalized value', fontsize=12)
ax2.set_title('E91 Protocol: CHSH Violation & Key Rate vs Noise', fontsize=13)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 0.45)

plt.tight_layout()
plt.savefig('figures/bb84_e91_analysis.png', dpi=200, bbox_inches='tight')
plt.savefig('figures/bb84_e91_analysis.svg', bbox_inches='tight')
plt.close()

# ─── key length vs epsilon (security parameter) ─────────────────────────────
epsilons = np.logspace(-12, -3, 50)
n_fixed = 100_000
q_fixed = 0.04

key_lengths = [bb84_key_length(n_fixed, q_fixed, eps) for eps in epsilons]

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(epsilons, key_lengths, 'navy', linewidth=2.5)
ax.fill_between(epsilons, 0, key_lengths, alpha=0.15, color='navy')
ax.set_xlabel('Security parameter ε', fontsize=12)
ax.set_ylabel('Secret key length (bits)', fontsize=12)
ax.set_title(f'BB84 Key Length vs Security Parameter\n(n={n_fixed:,}, QBER={q_fixed})', fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_ylim(0)
plt.tight_layout()
plt.savefig('figures/bb84_security_parameter.png', dpi=200, bbox_inches='tight')
plt.close()

# ─── collect results ─────────────────────────────────────────────────────────
results = {
    "bb84": {
        "qber_levels_tested": qber_levels,
        "asymptotic_key_rates": asymptotic,
        "finite_key_rate_at_n1e6_qber003": float(bb84_key_rate([int(1e6)], qber=0.03)[0]),
        "finite_key_rate_at_n1e4_qber003": float(bb84_key_rate([int(1e4)], qber=0.03)[0]),
        "minimum_block_size_qber003": int(n_values[np.where(bb84_key_rate(n_values, qber=0.03) > 0.01)[0][0]]) if len(np.where(bb84_key_rate(n_values, qber=0.03) > 0.01)[0]) > 0 else int(n_values[-1]),
        "key_length_n100k_q004_eps1e8": float(bb84_key_length(100_000, 0.04, 1e-8)),
    },
    "e91": {
        "S_max_pure": float(2 * np.sqrt(2)),
        "S_classical_threshold": 2.0,
        "key_rate_pure_state_n1e6": float(e91_key_rate_from_chsh(2 * np.sqrt(2), int(1e6))),
        "noise_threshold_for_positive_key": float(noise_params[np.where(rates_e91_large > 0)[0][0]]) if any(rates_e91_large > 0) else "N/A",
    }
}

os.makedirs('results', exist_ok=True)
with open('results/bb84_e91_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("BB84/E91 analysis complete.")
print(json.dumps(results, indent=2))
