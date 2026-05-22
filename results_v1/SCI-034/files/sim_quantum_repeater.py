"""
Quantum Repeater Memory Requirements & Performance Estimation
=============================================================
Model: Elementary link generation → entanglement swapping → purification
Reference: Sangouard et al. Rev. Mod. Phys. 83, 33 (2011)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import expon
import json, os

np.random.seed(42)

# ─── physical parameters ──────────────────────────────────────────────────────
C_FIBER = 2e8          # m/s  (speed of light in fiber)
ALPHA_DB_KM = 0.2      # dB/km fiber loss
ETA_DETECTOR = 0.85    # detector efficiency
ETA_COUPLING = 0.9     # memory coupling efficiency

def fiber_transmittance(L_km):
    """Fiber transmittance η(L) = 10^(-αL/10)."""
    return 10 ** (-ALPHA_DB_KM * L_km / 10)


def elementary_link_success_prob(L_km, p_det=ETA_DETECTOR, p_coup=ETA_COUPLING):
    """Probability of successful Bell-state measurement on an elementary link."""
    eta = fiber_transmittance(L_km / 2)  # half-link each side
    return (eta * p_det * p_coup) ** 2


def waiting_time_geometric(p):
    """Expected waiting time in units of L/(2c) for geometric trials."""
    return 1.0 / p  # in channel round-trip units


# ─── nested repeater chain ────────────────────────────────────────────────────
def repeater_chain_analysis(L_total_km, n_segments_list):
    """
    Analyze performance of nested quantum repeater chain.
    Returns: generation rate, fidelity, memory requirement.
    """
    results = []
    for n_seg in n_segments_list:
        L_seg = L_total_km / n_seg
        p_link = elementary_link_success_prob(L_seg)
        t_link = (L_seg / 2) / C_FIBER * 1e3  # ms, round-trip

        # Expected attempts until success (geometric)
        attempts_per_link = 1.0 / p_link

        # Entanglement swapping success probability
        p_swap = 0.5  # linear optics BSM

        # Total chain generation rate (nested, n_levels = log2(n_seg))
        n_levels = int(np.ceil(np.log2(n_seg)))
        rate_seg = p_link / t_link  # Hz per elementary link
        chain_rate = rate_seg / (n_seg * (1 / p_swap) ** n_levels)

        # Memory coherence time required: t_coh ≥ t_link * attempts_per_link
        t_coh_required_ms = t_link * attempts_per_link * 2  # safety factor 2

        # Dephasing fidelity loss (assuming T2 = t_coh_required * 10)
        # F ≈ exp(-t_wait/T2), t_wait = t_coh_required
        F_loss_per_level = np.exp(-1 / 10)  # 1/10 of T2
        F_chain = F_loss_per_level ** n_levels

        results.append({
            'n_segments': n_seg,
            'L_seg_km': L_seg,
            'link_success_prob': p_link,
            'attempts_per_link': attempts_per_link,
            't_link_ms': t_link,
            'chain_rate_Hz': chain_rate,
            't_coh_required_ms': t_coh_required_ms,
            'n_levels': n_levels,
            'fidelity_estimate': F_chain,
            'min_memory_per_node': 2 * int(np.ceil(attempts_per_link)),
        })
    return results


# ─── Tokyo-scale analysis: L_total = 250 km ──────────────────────────────────
L_TOTAL = 250  # km  (Tokyo QKD network approximate span)
n_segs_list = [1, 2, 4, 8, 16, 32, 64]

chain_data = repeater_chain_analysis(L_TOTAL, n_segs_list)

# Memory time vs number of segments
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

n_segs = [d['n_segments'] for d in chain_data]
t_cohs = [d['t_coh_required_ms'] for d in chain_data]
rates = [d['chain_rate_Hz'] for d in chain_data]
fids = [d['fidelity_estimate'] for d in chain_data]
mem_reqs = [d['min_memory_per_node'] for d in chain_data]

axes[0, 0].semilogy(n_segs, t_cohs, 'bo-', linewidth=2, markersize=8)
axes[0, 0].set_xlabel('Number of repeater segments', fontsize=11)
axes[0, 0].set_ylabel('Required coherence time (ms)', fontsize=11)
axes[0, 0].set_title(f'Memory Coherence Time Required\n(Total link: {L_TOTAL} km)', fontsize=12)
axes[0, 0].grid(True, alpha=0.3)
# Add horizontal lines for current memory technologies
axes[0, 0].axhline(1, color='red', linestyle='--', label='Solid-state (~1 ms)')
axes[0, 0].axhline(100, color='orange', linestyle='--', label='Rare-earth (~100 ms)')
axes[0, 0].axhline(1000, color='green', linestyle='--', label='Trapped ion (~1 s)')
axes[0, 0].legend(fontsize=8)

axes[0, 1].semilogy(n_segs, [max(r, 1e-10) for r in rates], 'rs-', linewidth=2, markersize=8)
axes[0, 1].set_xlabel('Number of repeater segments', fontsize=11)
axes[0, 1].set_ylabel('Entanglement generation rate (Hz)', fontsize=11)
axes[0, 1].set_title('E2E Entanglement Rate vs Segments', fontsize=12)
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(n_segs, fids, 'g^-', linewidth=2, markersize=8)
axes[1, 0].set_xlabel('Number of repeater segments', fontsize=11)
axes[1, 0].set_ylabel('Estimated end-to-end fidelity', fontsize=11)
axes[1, 0].set_title('Fidelity vs Number of Segments\n(Dephasing model)', fontsize=12)
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].axhline(0.5, color='red', linestyle=':', label='Classical threshold')
axes[1, 0].legend(fontsize=9)

axes[1, 1].bar(range(len(n_segs)), mem_reqs, color=plt.cm.plasma(np.linspace(0.2, 0.8, len(n_segs))))
axes[1, 1].set_xticks(range(len(n_segs)))
axes[1, 1].set_xticklabels([str(n) for n in n_segs])
axes[1, 1].set_xlabel('Number of repeater segments', fontsize=11)
axes[1, 1].set_ylabel('Min quantum memories per node', fontsize=11)
axes[1, 1].set_title('Memory Requirement per Repeater Node', fontsize=12)
axes[1, 1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('figures/quantum_repeater_analysis.png', dpi=200, bbox_inches='tight')
plt.savefig('figures/quantum_repeater_analysis.svg', bbox_inches='tight')
plt.close()

# ─── rate-fidelity tradeoff curve ────────────────────────────────────────────
T2_values_ms = [1, 10, 100, 1000]  # different memory coherence times
n_seg_range = np.array([2, 4, 8, 16, 32, 64])

fig, ax = plt.subplots(figsize=(8, 6))
colors = plt.cm.cool(np.linspace(0.1, 0.9, len(T2_values_ms)))

for T2, c in zip(T2_values_ms, colors):
    rates_T2 = []
    fids_T2 = []
    for n_seg in n_seg_range:
        L_seg = L_TOTAL / n_seg
        p_link = elementary_link_success_prob(L_seg)
        t_link = (L_seg / 2) / C_FIBER * 1e3
        n_levels = int(np.ceil(np.log2(n_seg)))
        rate = p_link / (t_link * n_seg * (2 ** n_levels))
        t_wait = t_link / p_link
        F = np.exp(-t_wait / T2) ** n_levels
        rates_T2.append(max(rate, 1e-12))
        fids_T2.append(F)
    ax.loglog(rates_T2, fids_T2, 'o-', color=c, linewidth=2, markersize=7,
              label=f'T₂ = {T2} ms')

ax.axhline(0.5, color='gray', linestyle=':', label='Classical fidelity bound')
ax.set_xlabel('Entanglement generation rate (Hz)', fontsize=12)
ax.set_ylabel('End-to-end fidelity', fontsize=12)
ax.set_title('Rate–Fidelity Tradeoff for Quantum Repeater Chain\n(250 km total distance)', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('figures/repeater_rate_fidelity_tradeoff.png', dpi=200, bbox_inches='tight')
plt.close()

# ─── save results ────────────────────────────────────────────────────────────
results_out = {
    "total_distance_km": L_TOTAL,
    "chain_analysis": chain_data,
    "optimal_segments_rate": n_segs[int(np.argmax(rates))],
    "optimal_segments_fidelity": n_segs[int(np.argmax(fids))],
    "memory_requirements_summary": {
        "2_segments": chain_data[1],
        "8_segments": chain_data[3],
        "32_segments": chain_data[5],
    }
}

with open('results/quantum_repeater_results.json', 'w') as f:
    json.dump(results_out, f, indent=2)

print("Quantum repeater analysis complete.")
for d in chain_data:
    print(f"  n={d['n_segments']:3d}: rate={d['chain_rate_Hz']:.3e} Hz, "
          f"T_coh≥{d['t_coh_required_ms']:.1f} ms, F≈{d['fidelity_estimate']:.3f}, "
          f"mem/node≥{d['min_memory_per_node']}")
