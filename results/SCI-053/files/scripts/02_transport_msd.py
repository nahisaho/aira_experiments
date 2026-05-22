#!/usr/bin/env python3
"""
Transport Properties: MSD-based Diffusion Coefficient Analysis
===============================================================
Computes self-diffusion coefficients from Mean Square Displacement (MSD).

Theory:
  D = lim(t→∞) 1/(6t) * <|r(t) - r(0)|²>

  With finite-size correction (Yeh-Hummer):
    D_corr = D_PBC + (k_B T ξ) / (6π η L)
    where ξ = 2.837297 for cubic box

References:
  - Allen, M.P. & Tildesley, D.J., Computer Simulation of Liquids (2017)
  - Yeh, I.-C. & Hummer, G., J. Phys. Chem. B 108, 15873 (2004)
"""

import numpy as np
import json
import os

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)


def compute_msd(positions, dt_ps):
    """
    Compute MSD using multiple time origins for better statistics.

    Parameters
    ----------
    positions : array (n_frames, n_atoms, 3), unwrapped coordinates (nm)
    dt_ps : float, time between frames (ps)

    Returns
    -------
    time : array, lag times (ps)
    msd : array, mean square displacement (nm²)
    """
    n_frames = positions.shape[0]
    max_lag = n_frames // 2

    time = np.arange(max_lag) * dt_ps
    msd = np.zeros(max_lag)
    counts = np.zeros(max_lag)

    for lag in range(1, max_lag):
        displacements = positions[lag:] - positions[:-lag]
        sq_disp = np.sum(displacements**2, axis=-1)  # per atom
        msd[lag] = np.mean(sq_disp)
        counts[lag] = sq_disp.size

    return time, msd


def fit_diffusion(time_ps, msd_nm2, fit_start_frac=0.2, fit_end_frac=0.8):
    """
    Fit diffusion coefficient from linear region of MSD.

    D = slope / 6 (3D diffusion)

    Parameters
    ----------
    time_ps : array
    msd_nm2 : array
    fit_start_frac, fit_end_frac : fraction of data for fitting

    Returns
    -------
    D_nm2_ps : float, diffusion coefficient (nm²/ps)
    D_cm2_s : float, diffusion coefficient (cm²/s)
    r_squared : float, R² of linear fit
    """
    n = len(time_ps)
    i_start = int(n * fit_start_frac)
    i_end = int(n * fit_end_frac)

    t_fit = time_ps[i_start:i_end]
    msd_fit = msd_nm2[i_start:i_end]

    coeffs = np.polyfit(t_fit, msd_fit, 1)
    slope = coeffs[0]

    D_nm2_ps = slope / 6.0
    D_cm2_s = D_nm2_ps * 1e-7  # nm²/ps → cm²/s

    # R²
    msd_pred = np.polyval(coeffs, t_fit)
    ss_res = np.sum((msd_fit - msd_pred)**2)
    ss_tot = np.sum((msd_fit - np.mean(msd_fit))**2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return D_nm2_ps, D_cm2_s, r_squared


def yeh_hummer_correction(D_pbc, T_K, eta_Pa_s, L_nm):
    """
    Yeh-Hummer finite-size correction for diffusion coefficient.

    D_inf = D_PBC + k_B T ξ / (6π η L)

    Parameters
    ----------
    D_pbc : float, PBC diffusion coefficient (cm²/s)
    T_K : float, temperature (K)
    eta_Pa_s : float, viscosity (Pa·s)
    L_nm : float, box length (nm)

    Returns
    -------
    D_corrected : float, corrected diffusion coefficient (cm²/s)
    correction : float, correction term (cm²/s)
    """
    kB = 1.380649e-23  # J/K
    xi = 2.837297
    L_m = L_nm * 1e-9

    correction = (kB * T_K * xi) / (6.0 * np.pi * eta_Pa_s * L_m)
    correction_cm2_s = correction * 1e4  # m²/s → cm²/s
    D_corrected = D_pbc + correction_cm2_s

    return D_corrected, correction_cm2_s


def generate_demo_msd(species="Li+", n_frames=10000, dt_ps=0.1):
    """Generate realistic demo MSD data."""
    time = np.arange(n_frames) * dt_ps

    if species == "Li+":
        D_true = 1.0e-5  # cm²/s (typical for Li+ in EC/DMC)
        D_nm2_ps = D_true / 1e-7
    elif species == "PF6-":
        D_true = 0.6e-5
        D_nm2_ps = D_true / 1e-7
    elif species == "EC":
        D_true = 2.5e-5
        D_nm2_ps = D_true / 1e-7
    else:
        D_true = 3.0e-5
        D_nm2_ps = D_true / 1e-7

    # MSD = 6Dt + noise + subdiffusive early regime
    msd = 6.0 * D_nm2_ps * time
    # Add subdiffusive cage effect at short times
    cage = 0.02 * (1 - np.exp(-time / 5.0))
    msd = msd + cage
    # Add noise
    noise = np.random.normal(0, 0.005 * msd.max(), len(time)) * (time / time.max())
    msd += noise
    msd = np.maximum(msd, 0)
    msd[0] = 0

    return time, msd


def main():
    print("=" * 70)
    print("Transport Properties: MSD-based Diffusion Coefficients")
    print("=" * 70)

    species_list = ["Li+", "PF6-", "EC", "DMC"]
    T = 298.15
    L_box = 4.5  # nm
    eta_solvent = 0.0008  # Pa·s (typical for EC/DMC mixture)

    results = {}

    print(f"\n{'Species':<10} {'D_PBC (cm²/s)':<18} {'D_corr (cm²/s)':<18} {'R²':<8}")
    print("-" * 54)

    for species in species_list:
        time, msd = generate_demo_msd(species)
        D_nm2_ps, D_cm2_s, r2 = fit_diffusion(time, msd)
        D_corr, correction = yeh_hummer_correction(D_cm2_s, T, eta_solvent, L_box)

        results[species] = {
            "D_PBC_cm2_s": float(f"{D_cm2_s:.4e}"),
            "D_corrected_cm2_s": float(f"{D_corr:.4e}"),
            "YH_correction_cm2_s": float(f"{correction:.4e}"),
            "R_squared": float(f"{r2:.4f}"),
            "fit_range_frac": [0.2, 0.8]
        }
        print(f"{species:<10} {D_cm2_s:<18.4e} {D_corr:<18.4e} {r2:<8.4f}")

    # Nernst-Einstein conductivity estimate
    e = 1.602176634e-19  # C
    kB = 1.380649e-23
    N_A = 6.02214076e23
    c_salt = 1.0  # mol/L = 1000 mol/m³

    D_Li = results["Li+"]["D_corrected_cm2_s"] * 1e-4  # cm²/s → m²/s
    D_PF6 = results["PF6-"]["D_corrected_cm2_s"] * 1e-4

    sigma_NE = (c_salt * 1000 * N_A * e**2 / (kB * T)) * (D_Li + D_PF6)

    # Haven ratio (typically 0.3-0.7 for concentrated electrolytes)
    H_R = 0.45  # estimated
    sigma_actual = sigma_NE * H_R

    conductivity_results = {
        "Nernst_Einstein_conductivity_S_m": float(f"{sigma_NE:.4f}"),
        "Haven_ratio": H_R,
        "estimated_actual_conductivity_S_m": float(f"{sigma_actual:.4f}"),
        "note": "Haven ratio < 1 indicates correlated ion motion"
    }
    results["conductivity"] = conductivity_results

    print(f"\nNernst-Einstein conductivity: σ_NE = {sigma_NE:.4f} S/m")
    print(f"Estimated actual conductivity: σ = {sigma_actual:.4f} S/m (H_R = {H_R})")

    output = {
        "method": "MSD analysis with Yeh-Hummer correction",
        "temperature_K": T,
        "box_length_nm": L_box,
        "viscosity_Pa_s": eta_solvent,
        "diffusion_coefficients": results
    }

    with open("results/diffusion_results.json", 'w') as f:
        json.dump(output, f, indent=2)

    # Figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Panel A: MSD vs time
        ax = axes[0]
        colors = {'Li+': '#E91E63', 'PF6-': '#2196F3', 'EC': '#4CAF50', 'DMC': '#FF9800'}
        for species in species_list:
            time, msd = generate_demo_msd(species)
            ax.plot(time / 1000, msd, label=species, color=colors.get(species, 'gray'),
                    linewidth=1.5)
        ax.set_xlabel("Time (ns)", fontsize=11)
        ax.set_ylabel("MSD (nm²)", fontsize=11)
        ax.set_title("(A) Mean Square Displacement", fontsize=12)
        ax.legend(fontsize=9)

        # Panel B: Log-log MSD (check anomalous diffusion)
        ax = axes[1]
        for species in species_list:
            time, msd = generate_demo_msd(species)
            mask = (time > 0) & (msd > 0)
            ax.loglog(time[mask] / 1000, msd[mask], label=species,
                      color=colors.get(species, 'gray'), linewidth=1.5)
        # Reference slope = 1 line
        t_ref = np.logspace(-2, 1, 100)
        ax.loglog(t_ref, 0.1 * t_ref, '--', color='gray', alpha=0.5, label='slope = 1')
        ax.set_xlabel("Time (ns)", fontsize=11)
        ax.set_ylabel("MSD (nm²)", fontsize=11)
        ax.set_title("(B) Log-Log MSD (Anomalous Check)", fontsize=12)
        ax.legend(fontsize=9)

        # Panel C: Diffusion coefficients bar chart
        ax = axes[2]
        species_names = list(results.keys())
        species_names = [s for s in species_names if s != "conductivity"]
        D_pbc = [results[s]["D_PBC_cm2_s"] * 1e5 for s in species_names]
        D_corr_vals = [results[s]["D_corrected_cm2_s"] * 1e5 for s in species_names]

        x = np.arange(len(species_names))
        width = 0.35
        ax.bar(x - width/2, D_pbc, width, label='D (PBC)', color='#90CAF9')
        ax.bar(x + width/2, D_corr_vals, width, label='D (YH-corrected)', color='#1565C0')
        ax.set_xticks(x)
        ax.set_xticklabels(species_names)
        ax.set_ylabel("D (×10⁻⁵ cm²/s)", fontsize=11)
        ax.set_title("(C) Diffusion Coefficients", fontsize=12)
        ax.legend(fontsize=9)

        plt.tight_layout()
        plt.savefig("figures/diffusion_analysis.png", dpi=300, bbox_inches='tight')
        plt.savefig("figures/diffusion_analysis.svg", bbox_inches='tight')
        plt.close()
        print("\nFigures saved to figures/diffusion_analysis.png/.svg")
    except ImportError:
        print("matplotlib not available; skipping figure generation")

    return output


if __name__ == "__main__":
    main()
