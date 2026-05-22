"""
Decoherence and Channel Loss Simulation
=========================================
Models:
  - T1/T2 decoherence (amplitude/phase damping)
  - Fiber channel loss
  - Combined QKD channel fidelity vs distance
  - Monte Carlo simulation of qubit transmission
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import expm
import json

np.random.seed(42)


# ─── quantum channel models ───────────────────────────────────────────────────
def amplitude_damping_kraus(gamma):
    """Amplitude damping Kraus operators. gamma = 1 - exp(-t/T1)."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    return [K0, K1]


def phase_damping_kraus(lambda_):
    """Phase damping (dephasing) Kraus operators. lambda = 1 - exp(-t/T2*)."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - lambda_)]])
    K1 = np.array([[0, 0], [0, np.sqrt(lambda_)]])
    return [K0, K1]


def depolarizing_kraus(p):
    """Depolarizing channel Kraus operators."""
    K0 = np.sqrt(1 - p) * np.eye(2)
    K1 = np.sqrt(p / 3) * np.array([[0, 1], [1, 0]])   # X
    K2 = np.sqrt(p / 3) * np.array([[0, -1j], [1j, 0]])  # Y
    K3 = np.sqrt(p / 3) * np.array([[1, 0], [0, -1]])   # Z
    return [K0, K1, K2, K3]


def apply_channel(rho, kraus_ops):
    """Apply quantum channel to density matrix rho."""
    out = np.zeros_like(rho, dtype=complex)
    for K in kraus_ops:
        out += K @ rho @ K.conj().T
    return out


def qubit_fidelity(rho1, rho2):
    """Fidelity F(rho1, rho2) = Tr(sqrt(sqrt(rho1)*rho2*sqrt(rho1)))."""
    # For pure state rho1 = |psi><psi|, F = <psi|rho2|psi>
    from scipy.linalg import sqrtm
    sqrt_r1 = sqrtm(rho1)
    M = sqrt_r1 @ rho2 @ sqrt_r1
    F = np.real(np.trace(sqrtm(M)))
    return min(max(F, 0.0), 1.0)


# ─── Bell state (ideal) ───────────────────────────────────────────────────────
phi_plus = np.array([1, 0, 0, 1]) / np.sqrt(2)
rho_bell = np.outer(phi_plus, phi_plus.conj())


def bell_state_fidelity_after_decoherence(t_us, T1_us, T2_us, system='fiber_memory'):
    """
    Fidelity of Bell pair after memory storage time t.
    Combined T1 (amplitude damping) and T2 (phase damping) on each qubit.
    """
    gamma = 1 - np.exp(-t_us / T1_us)
    lambda_ = 1 - np.exp(-t_us / T2_us)

    # Single qubit channel
    K_ad = amplitude_damping_kraus(gamma)
    K_pd = phase_damping_kraus(lambda_)

    # Single qubit density matrix for tracing (start with |+> state for phase test)
    rho_q = np.array([[0.5, 0.5], [0.5, 0.5]])  # |+> state
    rho_ideal = np.array([[0.5, 0.5], [0.5, 0.5]])

    rho_q_ad = apply_channel(rho_q, K_ad)
    rho_q_pd = apply_channel(rho_q_ad, K_pd)

    F = qubit_fidelity(rho_ideal, rho_q_pd)
    # Bell pair: both qubits decohere → approximate F_bell ≈ F_qubit^2 for product noise
    return F ** 2


# ─── fiber channel model ─────────────────────────────────────────────────────
def fiber_qkd_channel(L_km, alpha_dB_km=0.2, n_dark=100, detector_eff=0.85,
                       source_rate_MHz=1):
    """
    Comprehensive fiber QKD channel model.
    Returns: transmittance, QBER, key_rate_kHz
    """
    # Transmittance
    eta_fiber = 10 ** (-alpha_dB_km * L_km / 10)
    eta_total = eta_fiber * detector_eff

    # QBER from dark counts
    mu = 0.1  # mean photon number per pulse
    # Signal counts per second
    R_signal = source_rate_MHz * 1e6 * mu * eta_total
    R_dark = n_dark  # dark counts per second per detector
    # QBER ≈ dark_counts / (signal + dark)
    if R_signal + R_dark > 0:
        QBER = R_dark / (R_signal + R_dark)
    else:
        QBER = 0.5

    # Add misalignment error
    QBER_align = 0.01  # 1% alignment error
    QBER_total = min(QBER + QBER_align, 0.5)

    # Asymptotic BB84 key rate
    def h2(p):
        if p <= 0 or p >= 1:
            return 0.0
        return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

    key_rate = max(source_rate_MHz * 1e3 * mu * eta_total * (1 - h2(QBER_total) - h2(QBER_total)), 0)

    return eta_total, QBER_total, key_rate


# ─── simulation sweeps ────────────────────────────────────────────────────────
distances = np.linspace(1, 300, 200)
eta_arr, qber_arr, rate_arr = zip(*[fiber_qkd_channel(L) for L in distances])
eta_arr = np.array(eta_arr)
qber_arr = np.array(qber_arr)
rate_arr = np.array(rate_arr)

# Max distance where positive key rate
max_dist_idx = np.where(rate_arr > 0)[0]
max_dist = distances[max_dist_idx[-1]] if len(max_dist_idx) > 0 else 0

# T1/T2 decoherence sweep
times = np.logspace(-1, 4, 100)  # 0.1 μs to 10 ms
T1T2_configs = [
    (10, 5, 'Superconducting (T1=10μs, T2=5μs)', 'red'),
    (1000, 300, 'NV center (T1=1ms, T2=300μs)', 'blue'),
    (100000, 10000, 'Rare-earth (T1=100ms, T2=10ms)', 'green'),
    (1e9, 1e7, 'Trapped ion (T1=1000s, T2=10s)', 'purple'),
]

# ─── figures ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# 1. Transmittance vs distance
ax = axes[0, 0]
ax.semilogy(distances, eta_arr, 'b-', linewidth=2.5, label='Total η (fiber + detector)')
ax.semilogy(distances, 10 ** (-0.2 * distances / 10), 'r--', linewidth=1.5,
            label='Fiber transmittance (α=0.2 dB/km)')
ax.set_xlabel('Fiber length (km)', fontsize=11)
ax.set_ylabel('Transmittance η', fontsize=11)
ax.set_title('QKD Channel Transmittance vs Distance', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, which='both')
ax.axvline(100, color='gray', linestyle=':', alpha=0.6, label='100 km reference')
ax.set_xlim(0, 300)

# 2. QBER and key rate vs distance
ax = axes[0, 1]
ax_twin = ax.twinx()
l1, = ax.plot(distances, qber_arr * 100, 'r-', linewidth=2, label='QBER (%)')
l2, = ax_twin.semilogy(distances, np.maximum(rate_arr, 1e-6), 'b-', linewidth=2,
                        label='Key rate (kHz)')
ax_twin.axvline(max_dist, color='black', linestyle='--', alpha=0.7,
                label=f'Max dist = {max_dist:.0f} km')
ax.set_xlabel('Fiber length (km)', fontsize=11)
ax.set_ylabel('QBER (%)', fontsize=11, color='red')
ax_twin.set_ylabel('Key rate (kHz)', fontsize=11, color='blue')
ax.set_title('QBER & Secret Key Rate vs Distance\n(BB84, 1 GHz source, 100 dark counts/s)', fontsize=11)
ax.tick_params(axis='y', labelcolor='red')
ax_twin.tick_params(axis='y', labelcolor='blue')
lines = [l1, l2]
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, fontsize=9, loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 300)

# 3. Decoherence fidelity vs storage time
ax = axes[1, 0]
for T1, T2, label, color in T1T2_configs:
    fids = [bell_state_fidelity_after_decoherence(t, T1, T2) for t in times]
    ax.semilogx(times, fids, color=color, linewidth=2, label=label)

ax.axhline(0.5, color='gray', linestyle=':', label='Classical threshold F=0.5')
ax.axhline(2/3, color='black', linestyle='--', alpha=0.5, label='Distillation threshold F=2/3')
ax.set_xlabel('Storage time (μs)', fontsize=11)
ax.set_ylabel('Bell pair fidelity', fontsize=11)
ax.set_title('Bell Pair Fidelity vs Memory Storage Time\n(T1+T2 decoherence)', fontsize=12)
ax.legend(fontsize=8, loc='lower left')
ax.grid(True, alpha=0.3, which='both')
ax.set_ylim(0, 1.05)
ax.set_xlim(times[0], times[-1])

# 4. Combined distance-decoherence heatmap (fidelity)
ax = axes[1, 1]
distances_2d = np.linspace(10, 200, 50)
storage_times = np.logspace(1, 5, 50)  # 10 μs to 100 ms
T1_ref, T2_ref = 100000, 10000  # Rare-earth crystal reference

fidelity_map = np.zeros((len(storage_times), len(distances_2d)))
for i, t in enumerate(storage_times):
    for j, L in enumerate(distances_2d):
        # Memory decoherence
        F_mem = bell_state_fidelity_after_decoherence(t, T1_ref, T2_ref)
        # Channel fidelity (depolarizing approximation)
        eta = 10 ** (-0.2 * L / 10) * 0.85
        F_ch = eta + (1 - eta) * 0.25
        fidelity_map[i, j] = F_mem * F_ch

im = ax.contourf(distances_2d, np.log10(storage_times), fidelity_map,
                  levels=20, cmap='viridis')
ax.contour(distances_2d, np.log10(storage_times), fidelity_map,
           levels=[0.5, 0.6, 0.7, 0.8], colors='white', linewidths=1.5)
plt.colorbar(im, ax=ax, label='End-to-end fidelity')
ax.set_xlabel('Channel distance (km)', fontsize=11)
ax.set_ylabel('Storage time log₁₀(μs)', fontsize=11)
ax.set_title('Fidelity Map: Distance vs Storage Time\n(Rare-earth crystal memory, α=0.2 dB/km)', fontsize=11)

plt.tight_layout()
plt.savefig('figures/decoherence_channel_loss.png', dpi=200, bbox_inches='tight')
plt.savefig('figures/decoherence_channel_loss.svg', bbox_inches='tight')
plt.close()

# ─── Monte Carlo error simulation ────────────────────────────────────────────
N_SHOTS = 10000
n_qubits = 100
L_test = 50  # km

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Simulate BB84 raw key generation
def simulate_bb84_channel(n_bits, L_km, qber):
    """Monte Carlo simulation of BB84 key sifting."""
    alice_bits = np.random.randint(0, 2, n_bits)
    alice_bases = np.random.randint(0, 2, n_bits)
    bob_bases = np.random.randint(0, 2, n_bits)

    # Transmission losses: some photons don't arrive
    eta, _, _ = fiber_qkd_channel(L_km)
    arrived = np.random.random(n_bits) < eta

    # Bob measures
    bob_bits = alice_bits.copy()
    # Add QBER errors
    errors = np.random.random(n_bits) < qber
    bob_bits[errors] = 1 - bob_bits[errors]

    # Sifting: keep matching bases where photon arrived
    sifted = (alice_bases == bob_bases) & arrived
    alice_sifted = alice_bits[sifted]
    bob_sifted = bob_bits[sifted]

    # Measured QBER
    measured_qber = np.mean(alice_sifted != bob_sifted) if len(alice_sifted) > 0 else 0.5

    return len(alice_sifted), measured_qber

N_RUNS = 500
lengths = [10, 20, 50, 100, 150]
results_mc = {L: [] for L in lengths}

for L in lengths:
    _, true_qber, _ = fiber_qkd_channel(L)
    for _ in range(N_RUNS):
        n_sifted, meas_qber = simulate_bb84_channel(10000, L, true_qber)
        results_mc[L].append((n_sifted, meas_qber))

ax = axes[0]
for L in lengths:
    sifted_counts = [r[0] for r in results_mc[L]]
    ax.hist(sifted_counts, bins=30, alpha=0.6, label=f'{L} km', density=True)
ax.set_xlabel('Sifted key length (bits per 10k sent)', fontsize=11)
ax.set_ylabel('Probability density', fontsize=11)
ax.set_title('Monte Carlo: Sifted Key Length Distribution', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

ax = axes[1]
qbers_by_L = {L: [r[1] for r in results_mc[L]] for L in lengths}
ax.boxplot([qbers_by_L[L] for L in lengths], labels=[f'{L}km' for L in lengths],
           patch_artist=True, boxprops=dict(facecolor='lightblue', alpha=0.7))
ax.set_xlabel('Fiber distance', fontsize=11)
ax.set_ylabel('Measured QBER', fontsize=11)
ax.set_title('Monte Carlo: QBER Distribution per Distance\n(500 runs, 10k pulses each)', fontsize=12)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('figures/monte_carlo_qkd.png', dpi=200, bbox_inches='tight')
plt.close()

# ─── save results ────────────────────────────────────────────────────────────
results_out = {
    "fiber_channel": {
        "max_distance_km": float(max_dist),
        "qber_at_50km": float(fiber_qkd_channel(50)[1]),
        "key_rate_at_50km_kHz": float(fiber_qkd_channel(50)[2]),
        "key_rate_at_100km_kHz": float(fiber_qkd_channel(100)[2]),
        "transmittance_at_100km": float(fiber_qkd_channel(100)[0]),
    },
    "decoherence": {
        "bell_fidelity_SC_at_100us": float(bell_state_fidelity_after_decoherence(100, 10, 5)),
        "bell_fidelity_NV_at_1ms": float(bell_state_fidelity_after_decoherence(1000, 1000, 300)),
        "bell_fidelity_RE_at_1ms": float(bell_state_fidelity_after_decoherence(1000, 100000, 10000)),
        "bell_fidelity_Ion_at_1ms": float(bell_state_fidelity_after_decoherence(1000, 1e9, 1e7)),
    }
}

with open('results/decoherence_results.json', 'w') as f:
    json.dump(results_out, f, indent=2)

print("Decoherence & channel loss simulation complete.")
print(json.dumps(results_out, indent=2))
