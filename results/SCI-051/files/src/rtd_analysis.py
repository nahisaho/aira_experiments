"""
Residence Time Distribution (RTD) Analysis
============================================
Theoretical RTD models (PFR, CSTR, axial dispersion, tanks-in-series)
and experimental RTD pulse/step response analysis.
"""

import numpy as np
import json, os

def rtd_pfr(t, tau):
    """Ideal PFR: Dirac delta at t = tau."""
    E = np.zeros_like(t, dtype=float)
    idx = np.argmin(np.abs(t - tau))
    dt = t[1] - t[0] if len(t) > 1 else 1.0
    E[idx] = 1.0 / dt
    return E

def rtd_cstr(t, tau):
    """Single ideal CSTR: E(t) = (1/τ) exp(-t/τ)."""
    return (1.0 / tau) * np.exp(-t / tau)

def rtd_tanks_in_series(t, tau, N):
    """Tanks-in-series model: E(t) = N/τ * (Nt/τ)^(N-1) / (N-1)! * exp(-Nt/τ)."""
    from math import factorial
    tau_i = tau / N
    theta = t / tau_i
    E = (1.0 / tau_i) * (theta ** (N - 1)) * np.exp(-theta) / factorial(N - 1)
    return E

def rtd_axial_dispersion(t, tau, Pe):
    """Axial dispersion model (open-open boundary): 
    E(t) = 1/(2*sqrt(π*t*D/uL)) * exp(-(1-t/τ)²/(4t/(τ*Pe)))"""
    theta = t / tau
    mask = theta > 0
    E = np.zeros_like(t, dtype=float)
    E[mask] = (1.0 / tau) * np.sqrt(Pe / (4 * np.pi * theta[mask])) * \
              np.exp(-Pe * (1 - theta[mask])**2 / (4 * theta[mask]))
    return E

def compute_rtd_moments(t, E):
    """Compute mean residence time and variance from RTD."""
    dt = np.diff(t, prepend=t[0] - (t[1] - t[0]))
    norm = np.trapz(E, t)
    if norm > 0:
        E_norm = E / norm
    else:
        E_norm = E
    t_mean = np.trapz(t * E_norm, t)
    sigma2 = np.trapz((t - t_mean)**2 * E_norm, t)
    return t_mean, sigma2

def generate_synthetic_pulse_response(tau, Pe, t_max_factor=3, n_points=500, noise_level=0.02):
    """Generate synthetic experimental pulse response with noise."""
    t = np.linspace(0, tau * t_max_factor, n_points)
    E_true = rtd_axial_dispersion(t, tau, Pe)
    np.random.seed(42)
    noise = noise_level * np.max(E_true) * np.random.randn(n_points)
    E_noisy = np.maximum(E_true + noise, 0)
    return t, E_true, E_noisy

def fit_tanks_in_series(t_mean, sigma2):
    """Estimate N from moments: N = t_mean² / σ²."""
    if sigma2 > 0:
        N = t_mean**2 / sigma2
    else:
        N = float('inf')
    return N

def fit_peclet_from_variance(t_mean, sigma2):
    """Estimate Pe from σ²_θ = 2/Pe - 2/Pe²*(1-exp(-Pe))."""
    sigma2_theta = sigma2 / t_mean**2
    if sigma2_theta <= 0 or sigma2_theta >= 1:
        return 1.0
    Pe_est = 2.0 / sigma2_theta
    # Newton refinement
    for _ in range(20):
        f = 2.0 / Pe_est - 2.0 / Pe_est**2 * (1 - np.exp(-Pe_est)) - sigma2_theta
        df = -2.0 / Pe_est**2 + 4.0 / Pe_est**3 * (1 - np.exp(-Pe_est)) - 2.0 / Pe_est**2 * np.exp(-Pe_est)
        if abs(df) > 1e-12:
            Pe_est -= f / df
            Pe_est = max(Pe_est, 0.1)
    return Pe_est

def run_rtd_analysis():
    tau_design = 30.0  # s (design residence time at 1 mL/min)
    Pe_true = 50.0     # Peclet number

    t, E_true, E_noisy = generate_synthetic_pulse_response(tau_design, Pe_true)

    t_mean_true, sigma2_true = compute_rtd_moments(t, E_true)
    t_mean_exp, sigma2_exp = compute_rtd_moments(t, E_noisy)

    N_est = fit_tanks_in_series(t_mean_exp, sigma2_exp)
    Pe_est = fit_peclet_from_variance(t_mean_exp, sigma2_exp)

    # Compare models
    E_cstr = rtd_cstr(t, t_mean_exp)
    N_int = max(1, int(round(N_est)))
    E_tis = rtd_tanks_in_series(t, t_mean_exp, N_int)
    E_ad = rtd_axial_dispersion(t, t_mean_exp, Pe_est)

    results = {
        "design_parameters": {
            "design_tau_s": tau_design,
            "true_peclet": Pe_true,
        },
        "experimental_moments": {
            "mean_residence_time_s": round(t_mean_exp, 3),
            "variance_s2": round(sigma2_exp, 3),
            "dimensionless_variance": round(sigma2_exp / t_mean_exp**2, 5),
        },
        "model_fitting": {
            "tanks_in_series_N": round(N_est, 1),
            "tanks_in_series_N_int": N_int,
            "estimated_peclet": round(Pe_est, 1),
        },
        "model_comparison": {},
    }

    # Goodness of fit (R²)
    for name, E_model in [("CSTR", E_cstr), ("Tanks_in_Series", E_tis), ("Axial_Dispersion", E_ad)]:
        ss_res = np.sum((E_noisy - E_model)**2)
        ss_tot = np.sum((E_noisy - np.mean(E_noisy))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        results["model_comparison"][name] = {"R_squared": round(r2, 4)}

    rtd_data = {
        "time_s": t.tolist(),
        "E_true": E_true.tolist(),
        "E_experimental": E_noisy.tolist(),
        "E_cstr": E_cstr.tolist(),
        "E_tanks_in_series": E_tis.tolist(),
        "E_axial_dispersion": E_ad.tolist(),
    }

    return results, rtd_data

if __name__ == "__main__":
    results, rtd_data = run_rtd_analysis()
    os.makedirs("results", exist_ok=True)
    os.makedirs("data", exist_ok=True)

    with open("results/rtd_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    with open("data/rtd_curves.json", "w") as f:
        json.dump(rtd_data, f, indent=2)

    print("=== RTD Analysis Results ===")
    print(f"Mean Residence Time: {results['experimental_moments']['mean_residence_time_s']:.3f} s")
    print(f"Variance: {results['experimental_moments']['variance_s2']:.3f} s²")
    print(f"σ²_θ: {results['experimental_moments']['dimensionless_variance']:.5f}")
    print(f"\nEstimated N (tanks-in-series): {results['model_fitting']['tanks_in_series_N']:.1f}")
    print(f"Estimated Pe: {results['model_fitting']['estimated_peclet']:.1f}")
    print(f"\nModel Comparison (R²):")
    for name, val in results["model_comparison"].items():
        print(f"  {name}: {val['R_squared']:.4f}")
