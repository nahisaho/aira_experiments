#!/usr/bin/env python3
"""
Green-Kubo Ionic Conductivity Calculation
==========================================
Computes ionic conductivity from the current autocorrelation function (CACF).

Theory:
  σ = (1/3VkT) ∫₀^∞ <J(0)·J(t)> dt

  where J(t) = Σᵢ qᵢ vᵢ(t) is the collective charge current.

  This captures both self and cross-correlations (distinct from Nernst-Einstein).

  The ratio σ_GK / σ_NE = Haven ratio (H_R), which quantifies ion correlations.

References:
  - Kubo, R., J. Phys. Soc. Jpn. 12, 570 (1957)
  - Fong, K.D. et al., ACS Cent. Sci. 5, 1250 (2019)
  - France-Lanord, A. & Bhatt, M.D., Phys. Rev. Lett. 122, 136001 (2019)
"""

import numpy as np
import json
import os

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)


def compute_current(velocities, charges):
    """
    Compute collective charge current J(t) = Σᵢ qᵢ vᵢ(t).

    Parameters
    ----------
    velocities : array (n_frames, n_ions, 3) in nm/ps
    charges : array (n_ions,) in elementary charges

    Returns
    -------
    J : array (n_frames, 3), charge current
    """
    # J(t) = sum_i q_i * v_i(t)
    J = np.einsum('i,tij->tj', charges, velocities)
    return J


def autocorrelation_fft(signal, max_lag=None):
    """
    Compute autocorrelation function using FFT.

    Parameters
    ----------
    signal : array (n_frames, 3)
    max_lag : int, maximum lag time

    Returns
    -------
    acf : array, normalized autocorrelation
    """
    n = signal.shape[0]
    if max_lag is None:
        max_lag = n // 2

    acf = np.zeros(max_lag)

    for dim in range(signal.shape[1]):
        x = signal[:, dim]
        # FFT-based correlation
        X = np.fft.fft(x, n=2*n)
        corr = np.fft.ifft(X * np.conj(X)).real[:max_lag]
        # Normalize by number of time origins
        norm = np.arange(n, n - max_lag, -1, dtype=float)
        acf += corr / norm

    return acf


def integrate_acf(acf, dt_ps, method="trapezoid"):
    """
    Integrate ACF to get transport coefficient.

    Returns running integral for convergence check.
    """
    running_integral = np.zeros_like(acf)
    for i in range(1, len(acf)):
        if method == "trapezoid":
            running_integral[i] = running_integral[i-1] + 0.5 * (acf[i] + acf[i-1]) * dt_ps
        else:
            running_integral[i] = running_integral[i-1] + acf[i] * dt_ps
    return running_integral


def green_kubo_conductivity(J_acf, dt_ps, V_nm3, T_K):
    """
    Compute Green-Kubo conductivity.

    σ = (1/3VkT) ∫ <J(0)·J(t)> dt

    Parameters
    ----------
    J_acf : array, current autocorrelation function
    dt_ps : float, timestep (ps)
    V_nm3 : float, system volume (nm³)
    T_K : float, temperature (K)

    Returns
    -------
    sigma_S_m : float, conductivity (S/m)
    running_sigma : array, running integral of conductivity
    """
    kB = 1.380649e-23  # J/K
    e = 1.602176634e-19  # C
    V_m3 = V_nm3 * 1e-27  # nm³ → m³

    running_integral = integrate_acf(J_acf, dt_ps)

    # Unit conversion:
    # J is in e·nm/ps, so <J·J> is in e²·nm²/ps²
    # Need to convert to SI: C²·m²/s²
    conv_factor = e**2 * (1e-9)**2 / (1e-12)**2  # e²·nm²/ps² → C²·m²/s²
    # Integral is in e²·nm²/ps (ACF times dt)
    conv_integral = e**2 * (1e-9)**2 / (1e-12)  # e²·nm²/ps → C²·m²/(s)

    prefactor = conv_integral / (3.0 * V_m3 * kB * T_K)
    running_sigma = running_integral * prefactor

    return running_sigma


def decompose_conductivity(vel_cation, vel_anion, q_cat, q_an, dt_ps, V_nm3, T_K):
    """
    Decompose conductivity into self and cross terms.

    σ = σ_self + σ_cross
    σ_self = (1/3VkT) Σᵢ qᵢ² ∫ <vᵢ(0)·vᵢ(t)> dt  (= Nernst-Einstein)
    σ_cross = σ_GK - σ_self

    Returns dict with decomposition.
    """
    n_cat = vel_cation.shape[1]
    n_an = vel_anion.shape[1]

    # Self contribution (velocity ACF per ion)
    charges_cat = np.full(n_cat, q_cat)
    charges_an = np.full(n_an, q_an)

    # Total current
    J_cat = compute_current(vel_cation, charges_cat)
    J_an = compute_current(vel_anion, charges_an)
    J_total = J_cat + J_an

    J_acf_total = autocorrelation_fft(J_total)
    sigma_GK = green_kubo_conductivity(J_acf_total, dt_ps, V_nm3, T_K)

    # Cation-cation correlation
    J_acf_cat = autocorrelation_fft(J_cat)
    sigma_cat = green_kubo_conductivity(J_acf_cat, dt_ps, V_nm3, T_K)

    # Anion-anion correlation
    J_acf_an = autocorrelation_fft(J_an)
    sigma_an = green_kubo_conductivity(J_acf_an, dt_ps, V_nm3, T_K)

    # Cross term
    J_acf_cross_running = sigma_GK - sigma_cat - sigma_an

    return {
        "sigma_total": sigma_GK,
        "sigma_cation": sigma_cat,
        "sigma_anion": sigma_an,
        "sigma_cross": J_acf_cross_running
    }


def generate_demo_conductivity(n_frames=50000, dt_ps=0.1, n_cation=40, n_anion=40):
    """Generate demo velocity data for testing."""
    np.random.seed(42)

    # Generate correlated velocities (simplified Langevin-like)
    kBT_m = 0.5  # nm²/ps² (approximate kBT/m for ions)
    vel_cat = np.random.normal(0, np.sqrt(kBT_m), (n_frames, n_cation, 3))
    vel_an = np.random.normal(0, np.sqrt(kBT_m), (n_frames, n_anion, 3))

    # Add temporal correlation (memory effect)
    tau = 2.0  # ps
    alpha = np.exp(-dt_ps / tau)
    for t in range(1, n_frames):
        vel_cat[t] = alpha * vel_cat[t-1] + np.sqrt(1 - alpha**2) * vel_cat[t]
        vel_an[t] = alpha * vel_an[t-1] + np.sqrt(1 - alpha**2) * vel_an[t]

    # Add ion-ion correlations (cross terms)
    cross_strength = -0.15
    for t in range(n_frames):
        mean_vel = np.mean(vel_cat[t], axis=0)
        vel_an[t] += cross_strength * mean_vel

    return vel_cat, vel_an


def main():
    print("=" * 70)
    print("Green-Kubo Ionic Conductivity Calculation")
    print("=" * 70)

    # Parameters
    T = 298.15
    V = 4.5**3  # nm³ (cubic box)
    dt = 0.1    # ps
    q_cat = 0.8  # ECC scaled charge
    q_an = -0.8
    n_cat = 40
    n_an = 40

    # Generate demo data
    print("\nGenerating demo velocity data...")
    vel_cat, vel_an = generate_demo_conductivity(
        n_frames=50000, dt_ps=dt, n_cation=n_cat, n_anion=n_an
    )

    # Compute conductivity decomposition
    print("Computing Green-Kubo conductivity...")
    decomp = decompose_conductivity(vel_cat, vel_an, q_cat, q_an, dt, V, T)

    # Extract converged values (plateau of running integral)
    def get_plateau(running, frac_start=0.3, frac_end=0.6):
        n = len(running)
        return np.mean(running[int(n*frac_start):int(n*frac_end)])

    sigma_GK = get_plateau(decomp["sigma_total"])
    sigma_cat = get_plateau(decomp["sigma_cation"])
    sigma_an = get_plateau(decomp["sigma_anion"])
    sigma_cross = sigma_GK - sigma_cat - sigma_an

    # Nernst-Einstein for comparison
    kB = 1.380649e-23
    e = 1.602176634e-19
    N_A = 6.02214076e23
    c = 1.0  # mol/L

    D_Li = 1.0e-9   # m²/s
    D_PF6 = 0.6e-9  # m²/s
    sigma_NE = (c * 1000 * N_A * e**2 / (kB * T)) * (D_Li + D_PF6)

    # Haven ratio
    H_R = sigma_GK / sigma_NE if sigma_NE > 0 else np.nan

    print(f"\nResults:")
    print(f"  σ_GK (total)    = {sigma_GK:.4f} S/m")
    print(f"  σ_cation-cation = {sigma_cat:.4f} S/m")
    print(f"  σ_anion-anion   = {sigma_an:.4f} S/m")
    print(f"  σ_cross         = {sigma_cross:.4f} S/m")
    print(f"  σ_NE            = {sigma_NE:.4f} S/m")
    print(f"  Haven ratio     = {H_R:.3f}")

    # Transference number
    t_plus = sigma_cat / sigma_GK if sigma_GK > 0 else np.nan
    t_minus = sigma_an / sigma_GK if sigma_GK > 0 else np.nan
    print(f"\n  t+ (cation)     = {t_plus:.3f}")
    print(f"  t- (anion)      = {t_minus:.3f}")

    results = {
        "method": "Green-Kubo conductivity with decomposition",
        "temperature_K": T,
        "volume_nm3": V,
        "charge_scaling": q_cat,
        "conductivity": {
            "sigma_GK_S_m": float(f"{sigma_GK:.4f}"),
            "sigma_NE_S_m": float(f"{sigma_NE:.4f}"),
            "Haven_ratio": float(f"{H_R:.3f}"),
            "sigma_cation_S_m": float(f"{sigma_cat:.4f}"),
            "sigma_anion_S_m": float(f"{sigma_an:.4f}"),
            "sigma_cross_S_m": float(f"{sigma_cross:.4f}"),
            "transference_number_cation": float(f"{t_plus:.3f}"),
            "transference_number_anion": float(f"{t_minus:.3f}")
        },
        "interpretation": {
            "Haven_ratio_meaning": "H_R < 1 indicates anti-correlated ion motion (ion pairing)",
            "cross_term_meaning": "Negative σ_cross indicates ion pair formation reducing conductivity",
            "transference_notes": "GK transference includes cross-correlations unlike NE"
        }
    }

    with open("results/green_kubo_conductivity.json", 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to results/green_kubo_conductivity.json")

    # Figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        time = np.arange(len(decomp["sigma_total"])) * dt / 1000  # ns

        # Panel A: Running integral (convergence)
        ax = axes[0]
        ax.plot(time, decomp["sigma_total"], label='σ_GK (total)', linewidth=1.5, color='#1565C0')
        ax.plot(time, decomp["sigma_cation"], label='σ_++', linewidth=1.2, color='#E91E63')
        ax.plot(time, decomp["sigma_anion"], label='σ_--', linewidth=1.2, color='#4CAF50')
        ax.set_xlabel("Correlation time (ns)", fontsize=11)
        ax.set_ylabel("Running integral σ(t) (S/m)", fontsize=11)
        ax.set_title("(A) Green-Kubo Running Integral", fontsize=12)
        ax.legend(fontsize=9)

        # Panel B: Current ACF
        ax = axes[1]
        charges_cat = np.full(vel_cat.shape[1], q_cat)
        charges_an = np.full(vel_an.shape[1], q_an)
        J_cat = compute_current(vel_cat, charges_cat)
        J_an = compute_current(vel_an, charges_an)
        J_total = J_cat + J_an
        acf = autocorrelation_fft(J_total, max_lag=5000)
        acf_norm = acf / acf[0] if acf[0] != 0 else acf
        t_acf = np.arange(len(acf_norm)) * dt

        ax.plot(t_acf, acf_norm, color='#1565C0', linewidth=1.5)
        ax.set_xlabel("Lag time (ps)", fontsize=11)
        ax.set_ylabel("C(t) / C(0)", fontsize=11)
        ax.set_title("(B) Current Autocorrelation Function", fontsize=12)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlim(0, 200)

        # Panel C: Conductivity decomposition pie chart
        ax = axes[2]
        labels = ['σ₊₊ (cation)', 'σ₋₋ (anion)', 'σ_cross']
        sizes = [abs(sigma_cat), abs(sigma_an), abs(sigma_cross)]
        colors_pie = ['#E91E63', '#4CAF50', '#FFC107']
        ax.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
               startangle=90, textprops={'fontsize': 10})
        ax.set_title("(C) Conductivity Decomposition", fontsize=12)

        plt.tight_layout()
        plt.savefig("figures/green_kubo_conductivity.png", dpi=300, bbox_inches='tight')
        plt.savefig("figures/green_kubo_conductivity.svg", bbox_inches='tight')
        plt.close()
        print("Figures saved to figures/green_kubo_conductivity.png/.svg")
    except ImportError:
        print("matplotlib not available; skipping figure generation")

    return results


if __name__ == "__main__":
    main()
