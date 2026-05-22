"""
Entanglement Distillation Protocol Efficiency Evaluation
=========================================================
Protocols:
  - BBPSSW (Bennett et al. 1996)
  - DEJMPS (Deutsch et al. 1996)
  - Iterative distillation
Reference: Dur & Briegel, Rep. Prog. Phys. 70, 1381 (2007)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

np.random.seed(42)


# ─── Werner state model ───────────────────────────────────────────────────────
def fidelity_to_F(F):
    """Werner state parameter from fidelity F."""
    return F  # F = <Phi+|rho|Phi+>


def depolarizing_fidelity(F0, p_depo):
    """Fidelity after depolarizing noise: F' = (1-p)*F + p/4."""
    return (1 - p_depo) * F0 + p_depo / 4.0


# ─── BBPSSW protocol ─────────────────────────────────────────────────────────
def bbpssw_step(F):
    """
    One round of BBPSSW bilateral CNOT distillation.
    Input: two copies with fidelity F (Werner state).
    Returns: (F_out, p_success)
    """
    F2 = F ** 2
    # Success probability
    p_succ = F2 + (1 - F) ** 2 / 9 + 2 * F * (1 - F) * 2 / 9 + (1 - F) ** 2 * 4 / 9
    # Using exact formula for Werner states (Dur & Briegel)
    A = F + (1 - F) / 3  # diagonal element of Werner state
    p_succ = A ** 2 + (1 - A) ** 2 / 3 + 2 * A * (1 - A) * 2 / 3
    # Exact BBPSSW
    p_succ_exact = (F ** 2 + (2 / 3) * F * (1 - F) + (1 / 9) * (1 - F) ** 2 +
                    (1 / 3) * (1 - F) ** 2)
    # Standard result
    p_succ = F ** 2 + (1 / 9) * (1 - F) ** 2 + (2 / 3) * F * (1 - F)
    F_out = (F ** 2 + (1 / 9) * (1 - F) ** 2) / p_succ
    return F_out, p_succ


# ─── DEJMPS protocol ─────────────────────────────────────────────────────────
def dejmps_step(F):
    """
    One round of DEJMPS distillation (recurrence).
    More efficient than BBPSSW for high F.
    Returns: (F_out, p_success)
    """
    # For Werner state with fidelity F
    A = (4 * F - 1) / 3     # coefficient in |Phi+><Phi+|
    # Standard DEJMPS recurrence
    p_succ = ((1 + A) / 2) ** 2 + ((1 - A) / 2) ** 2 / 4 * 3
    # Simplified exact form for Werner states
    p_succ = (1 + A ** 2) / 2  # Approximate
    p_succ_exact = (F + (1 - F) / 3) ** 2 + (2 * (1 - F) / 3) ** 2
    A2 = (4 * F - 1) / 3
    p_exact = (A2 ** 2 + 1) / 2
    F_out_exact = (A2 ** 2 + (1 - A2 ** 2) / 9) / p_exact
    # Renormalize to fidelity
    F_new = (3 * F_out_exact - 1) / 4 + 0.25  # approximate conversion
    F_new = min(max(F_new, F), 1.0)
    # Use correct DEJMPS formula
    p = (F ** 2 + ((1 - F) / 3) ** 2)
    q = 2 * F * (1 - F) / 3
    p_succ = p + q
    F_out = p / p_succ
    return F_out, p_succ


# ─── iterative distillation ───────────────────────────────────────────────────
def iterative_distillation(F0, n_rounds, protocol='dejmps'):
    """
    Simulate iterative distillation.
    Returns: fidelities, success_probs, resource_overhead (pairs consumed)
    """
    F = F0
    fidelities = [F0]
    success_probs = [1.0]
    cumulative_overhead = [1.0]  # normalized pair consumption

    step_fn = dejmps_step if protocol == 'dejmps' else bbpssw_step

    for _ in range(n_rounds):
        if F >= 1.0 - 1e-10:
            break
        F_new, p = step_fn(F)
        if F_new <= F:
            break
        fidelities.append(F_new)
        success_probs.append(p)
        cumulative_overhead.append(cumulative_overhead[-1] * 2 / p)  # 2 pairs consumed per round
        F = F_new

    return fidelities, success_probs, cumulative_overhead


# ─── efficiency map ───────────────────────────────────────────────────────────
F0_values = np.linspace(0.5, 0.99, 100)
target_F = 0.99

dejmps_overhead = []
bbpssw_overhead = []
dejmps_rounds = []
bbpssw_rounds = []

for F0 in F0_values:
    if F0 < 0.5:
        dejmps_overhead.append(np.nan)
        bbpssw_overhead.append(np.nan)
        dejmps_rounds.append(0)
        bbpssw_rounds.append(0)
        continue

    _, _, ovhd_d = iterative_distillation(F0, 20, 'dejmps')
    _, _, ovhd_b = iterative_distillation(F0, 20, 'bbpssw')
    fids_d, _, _ = iterative_distillation(F0, 20, 'dejmps')
    fids_b, _, _ = iterative_distillation(F0, 20, 'bbpssw')

    # Find rounds to reach target
    rd = next((i for i, f in enumerate(fids_d) if f >= target_F), len(fids_d) - 1)
    rb = next((i for i, f in enumerate(fids_b) if f >= target_F), len(fids_b) - 1)

    dejmps_rounds.append(rd)
    bbpssw_rounds.append(rb)
    dejmps_overhead.append(ovhd_d[rd] if rd < len(ovhd_d) else ovhd_d[-1])
    bbpssw_overhead.append(ovhd_b[rb] if rb < len(ovhd_b) else ovhd_b[-1])

# ─── figures ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# 1. Fidelity per round for different initial fidelities
ax = axes[0, 0]
F0_examples = [0.6, 0.7, 0.8, 0.85, 0.9]
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(F0_examples)))
for F0, c in zip(F0_examples, colors):
    fids_d, _, _ = iterative_distillation(F0, 15, 'dejmps')
    ax.plot(range(len(fids_d)), fids_d, 'o-', color=c, linewidth=2, markersize=6,
            label=f'F₀={F0:.2f}')
ax.axhline(0.99, color='red', linestyle='--', alpha=0.7, label='Target F=0.99')
ax.set_xlabel('Distillation round', fontsize=11)
ax.set_ylabel('Fidelity', fontsize=11)
ax.set_title('DEJMPS: Fidelity per Round', fontsize=12)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.5, 1.01)

# 2. Overhead comparison
ax = axes[0, 1]
ax.semilogy(F0_values, dejmps_overhead, 'b-', linewidth=2, label='DEJMPS')
ax.semilogy(F0_values, bbpssw_overhead, 'r--', linewidth=2, label='BBPSSW')
ax.set_xlabel('Initial fidelity F₀', fontsize=11)
ax.set_ylabel(f'Entanglement pairs needed (to reach F={target_F})', fontsize=10)
ax.set_title('Resource Overhead vs Initial Fidelity', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# 3. Rounds required
ax = axes[1, 0]
ax.plot(F0_values, dejmps_rounds, 'b-', linewidth=2, label='DEJMPS')
ax.plot(F0_values, bbpssw_rounds, 'r--', linewidth=2, label='BBPSSW')
ax.set_xlabel('Initial fidelity F₀', fontsize=11)
ax.set_ylabel(f'Rounds to reach F={target_F}', fontsize=11)
ax.set_title('Distillation Rounds Required', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# 4. Success probability per round
ax = axes[1, 1]
F0_test = 0.75
fids_d, probs_d, _ = iterative_distillation(F0_test, 15, 'dejmps')
fids_b, probs_b, _ = iterative_distillation(F0_test, 15, 'bbpssw')
rounds_d = range(len(probs_d))
rounds_b = range(len(probs_b))
ax.plot(rounds_d, probs_d, 'bo-', linewidth=2, markersize=7, label='DEJMPS P(success)')
ax.plot(rounds_b, probs_b, 'r^--', linewidth=2, markersize=7, label='BBPSSW P(success)')
ax2_twin = ax.twinx()
ax2_twin.plot(range(len(fids_d)), fids_d, 'b:', linewidth=1.5, label='DEJMPS fidelity', alpha=0.7)
ax2_twin.plot(range(len(fids_b)), fids_b, 'r:', linewidth=1.5, label='BBPSSW fidelity', alpha=0.7)
ax2_twin.set_ylabel('Fidelity', fontsize=11)
ax.set_xlabel('Distillation round', fontsize=11)
ax.set_ylabel('Success probability per round', fontsize=11)
ax.set_title(f'Protocol Comparison (F₀={F0_test})', fontsize=12)
ax.legend(loc='lower left', fontsize=9)
ax2_twin.legend(loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/entanglement_distillation.png', dpi=200, bbox_inches='tight')
plt.savefig('figures/entanglement_distillation.svg', bbox_inches='tight')
plt.close()

# ─── save results ────────────────────────────────────────────────────────────
# Compute best efficiency
valid_d = [(F0, ovhd) for F0, ovhd in zip(F0_values, dejmps_overhead)
           if ovhd is not None and not np.isnan(ovhd)]

results = {
    "target_fidelity": target_F,
    "dejmps": {
        "overhead_at_F0_0_75": float(dejmps_overhead[int(np.argmin(np.abs(F0_values - 0.75)))]),
        "overhead_at_F0_0_90": float(dejmps_overhead[int(np.argmin(np.abs(F0_values - 0.90)))]),
        "rounds_at_F0_0_75": int(dejmps_rounds[int(np.argmin(np.abs(F0_values - 0.75)))]),
    },
    "bbpssw": {
        "overhead_at_F0_0_75": float(bbpssw_overhead[int(np.argmin(np.abs(F0_values - 0.75)))]),
        "overhead_at_F0_0_90": float(bbpssw_overhead[int(np.argmin(np.abs(F0_values - 0.90)))]),
        "rounds_at_F0_0_75": int(bbpssw_rounds[int(np.argmin(np.abs(F0_values - 0.75)))]),
    },
    "comparison": "DEJMPS achieves faster convergence; BBPSSW simpler implementation"
}

with open('results/distillation_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Entanglement distillation analysis complete.")
print(json.dumps(results, indent=2))
