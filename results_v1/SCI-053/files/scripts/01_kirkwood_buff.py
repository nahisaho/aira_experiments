#!/usr/bin/env python3
"""
Kirkwood-Buff Integral Analysis for Activity Coefficients and Osmotic Pressure
===============================================================================
Computes KB integrals G_ij from radial distribution functions g(r),
then derives thermodynamic properties of concentrated electrolyte solutions.

Theory:
  G_ij = 4π ∫₀^∞ [g_ij(r) - 1] r² dr

  Activity coefficient derivative:
    ∂ln(γ±)/∂c = -(G++ + G-- - 2G+-) / (2 + c(G++ + G-- - 2G+-))

  Osmotic coefficient:
    φ = 1 - c_s * (G_ss - G_si) / (1 + c_s * G_ss)

References:
  - Kirkwood, J.G. & Buff, F.P., J. Chem. Phys. 19, 774 (1951)
  - Ben-Naim, A., Molecular Theory of Solutions, Oxford (2006)
  - Ganguly, P. & van der Vegt, N.F.A., JCTC 9, 1347 (2013)
"""

import numpy as np
import json
import os

# Output directories
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)


def load_rdf(filename, skiprows=0):
    """Load RDF data from GROMACS xvg or plain text format."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(('#', '@', ';')) or not line:
                continue
            vals = line.split()
            if len(vals) >= 2:
                data.append([float(v) for v in vals[:2]])
    return np.array(data)


def compute_kb_integral(r, g_r, r_max=None):
    """
    Compute Kirkwood-Buff integral G = 4π ∫₀^R [g(r)-1] r² dr

    Parameters
    ----------
    r : array, radial distance (nm)
    g_r : array, radial distribution function
    r_max : float, upper integration limit (nm). If None, use full range.

    Returns
    -------
    R_values : array, running upper limit
    G_values : array, running KB integral (nm³ → L/mol after conversion)
    """
    if r_max is not None:
        mask = r <= r_max
        r = r[mask]
        g_r = g_r[mask]

    integrand = 4.0 * np.pi * (g_r - 1.0) * r**2
    G_running = np.zeros_like(r)
    for i in range(1, len(r)):
        dr = r[i] - r[i-1]
        G_running[i] = G_running[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * dr

    return r, G_running


def kb_to_activity_coefficient(G_pp, G_mm, G_pm, c_salt):
    """
    Compute activity coefficient derivative from KB integrals.

    Parameters
    ----------
    G_pp : float, cation-cation KB integral (nm³)
    G_mm : float, anion-anion KB integral (nm³)
    G_pm : float, cation-anion KB integral (nm³)
    c_salt : float, salt concentration (mol/L = M)

    Returns
    -------
    dln_gamma_dc : float, ∂ln(γ±)/∂c
    """
    Delta = G_pp + G_mm - 2.0 * G_pm
    # Convert nm³ to L/mol: 1 nm³ = 6.022e-4 L/mol (with Avogadro)
    N_A = 6.02214076e23
    Delta_Lmol = Delta * 1e-24 * N_A  # nm³ → cm³ → L/mol

    denominator = 2.0 + c_salt * Delta_Lmol
    if abs(denominator) < 1e-10:
        return np.nan
    return -Delta_Lmol / denominator


def kb_to_osmotic_coefficient(G_ss, G_si, c_solvent):
    """
    Compute osmotic coefficient from KB integrals.

    Parameters
    ----------
    G_ss : float, solvent-solvent KB integral (nm³)
    G_si : float, solvent-ion KB integral (nm³)
    c_solvent : float, solvent concentration (mol/L)

    Returns
    -------
    phi : float, osmotic coefficient
    """
    N_A = 6.02214076e23
    G_ss_Lmol = G_ss * 1e-24 * N_A
    G_si_Lmol = G_si * 1e-24 * N_A

    numerator = c_solvent * (G_ss_Lmol - G_si_Lmol)
    denominator = 1.0 + c_solvent * G_ss_Lmol
    if abs(denominator) < 1e-10:
        return np.nan
    return 1.0 - numerator / denominator


def analyze_convergence(r, G_running, window_nm=0.5):
    """Check KB integral convergence by plateau detection."""
    plateau_start = len(r) // 2
    G_plateau = G_running[plateau_start:]
    r_plateau = r[plateau_start:]

    mean_G = np.mean(G_plateau)
    std_G = np.std(G_plateau)
    relative_fluct = std_G / abs(mean_G) if abs(mean_G) > 1e-10 else np.inf

    return {
        "converged_value_nm3": float(mean_G),
        "std_nm3": float(std_G),
        "relative_fluctuation": float(relative_fluct),
        "plateau_start_nm": float(r_plateau[0]),
        "converged": bool(relative_fluct < 0.1)
    }


def generate_demo_rdf(pair_type="cation-anion", n_points=500, r_max=2.0):
    """
    Generate realistic demo RDF for testing.
    Based on typical MD results for concentrated NaCl.
    """
    r = np.linspace(0.01, r_max, n_points)

    if pair_type == "cation-anion":
        # Strong first peak (contact ion pair)
        g = 1.0 + 8.0 * np.exp(-((r - 0.28)**2) / (2 * 0.01**2))
        # Solvent-separated peak
        g += 2.5 * np.exp(-((r - 0.50)**2) / (2 * 0.02**2))
        # Depletion
        g -= 0.3 * np.exp(-((r - 0.38)**2) / (2 * 0.015**2))
        g = np.maximum(g, 0)
        # Decay to 1
        g = 1.0 + (g - 1.0) * np.exp(-r / 0.8)
    elif pair_type == "cation-cation":
        # Exclusion zone, then peak
        g = np.zeros_like(r)
        mask = r > 0.3
        g[mask] = 1.0 + 1.5 * np.exp(-((r[mask] - 0.45)**2) / (2 * 0.02**2))
        g = 1.0 + (g - 1.0) * np.exp(-r / 0.6)
    elif pair_type == "solvent-solvent":
        g = 1.0 + 2.0 * np.exp(-((r - 0.28)**2) / (2 * 0.012**2))
        g += 1.0 * np.exp(-((r - 0.45)**2) / (2 * 0.02**2))
        g = 1.0 + (g - 1.0) * np.exp(-r / 0.5)
    else:
        g = np.ones_like(r)

    return r, g


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("Kirkwood-Buff Integral Analysis for Concentrated Electrolytes")
    print("=" * 70)

    # Generate demo RDFs (replace with actual data loading)
    pairs = {
        "cation-anion": "cation-anion",
        "cation-cation": "cation-cation",
        "anion-anion": "cation-cation",  # Similar structure
        "solvent-solvent": "solvent-solvent",
        "solvent-cation": "cation-anion",
    }

    results = {}
    kb_integrals = {}

    for pair_name, pair_type in pairs.items():
        r, g = generate_demo_rdf(pair_type)
        r_kb, G_kb = compute_kb_integral(r, g, r_max=1.8)
        convergence = analyze_convergence(r_kb, G_kb)
        kb_integrals[pair_name] = convergence["converged_value_nm3"]
        results[pair_name] = {
            "KB_integral_nm3": convergence["converged_value_nm3"],
            "convergence": convergence
        }
        print(f"\n{pair_name}:")
        print(f"  G = {convergence['converged_value_nm3']:.4f} nm³")
        print(f"  Converged: {convergence['converged']}")

    # Compute thermodynamic properties at multiple concentrations
    concentrations = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
    thermo_results = []

    print("\n" + "=" * 70)
    print("Thermodynamic Properties vs. Concentration")
    print("=" * 70)
    print(f"{'c (mol/L)':<12} {'∂ln(γ±)/∂c':<15} {'φ':<10}")
    print("-" * 37)

    for c in concentrations:
        dln_gamma = kb_to_activity_coefficient(
            kb_integrals["cation-cation"],
            kb_integrals["anion-anion"],
            kb_integrals["cation-anion"],
            c
        )
        phi = kb_to_osmotic_coefficient(
            kb_integrals["solvent-solvent"],
            kb_integrals["solvent-cation"],
            55.5 - c  # approximate solvent conc
        )
        thermo_results.append({
            "concentration_mol_L": c,
            "dln_gamma_dc": float(dln_gamma) if not np.isnan(dln_gamma) else None,
            "osmotic_coefficient": float(phi) if not np.isnan(phi) else None
        })
        print(f"{c:<12.1f} {dln_gamma:<15.4f} {phi:<10.4f}")

    # Save results
    output = {
        "method": "Kirkwood-Buff integral analysis",
        "kb_integrals": results,
        "thermodynamic_properties": thermo_results,
        "notes": [
            "Demo RDFs used for illustration; replace with actual MD trajectory data",
            "KB integrals converged from plateau region (r > r_max/2)",
            "Activity coefficients derived via Ben-Naim KB theory",
            "Unit conversion: 1 nm³ = 6.022e-4 L/mol"
        ]
    }

    with open("results/kb_analysis_results.json", 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to results/kb_analysis_results.json")

    # Generate figure (matplotlib)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Panel A: RDFs
        ax = axes[0, 0]
        for pair_name, pair_type in [("Li⁺–PF₆⁻", "cation-anion"),
                                      ("Li⁺–Li⁺", "cation-cation"),
                                      ("Solvent–Solvent", "solvent-solvent")]:
            r, g = generate_demo_rdf(pair_type)
            ax.plot(r, g, label=pair_name, linewidth=1.5)
        ax.set_xlabel("r (nm)", fontsize=11)
        ax.set_ylabel("g(r)", fontsize=11)
        ax.set_title("(A) Radial Distribution Functions", fontsize=12)
        ax.legend(fontsize=9)
        ax.set_xlim(0, 1.5)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)

        # Panel B: Running KB integrals
        ax = axes[0, 1]
        for pair_name, pair_type in [("G₊₋", "cation-anion"),
                                      ("G₊₊", "cation-cation"),
                                      ("G_ss", "solvent-solvent")]:
            r, g = generate_demo_rdf(pair_type)
            r_kb, G_kb = compute_kb_integral(r, g)
            ax.plot(r_kb, G_kb, label=pair_name, linewidth=1.5)
        ax.set_xlabel("R (nm)", fontsize=11)
        ax.set_ylabel("G(R) (nm³)", fontsize=11)
        ax.set_title("(B) Running KB Integrals", fontsize=12)
        ax.legend(fontsize=9)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        # Panel C: Activity coefficient derivative
        ax = axes[1, 0]
        concs = [t["concentration_mol_L"] for t in thermo_results]
        dlng = [t["dln_gamma_dc"] for t in thermo_results]
        ax.plot(concs, dlng, 'o-', color='#2196F3', linewidth=2, markersize=8)
        ax.set_xlabel("Concentration (mol/L)", fontsize=11)
        ax.set_ylabel("∂ln(γ±)/∂c", fontsize=11)
        ax.set_title("(C) Activity Coefficient Derivative", fontsize=12)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        # Panel D: Osmotic coefficient
        ax = axes[1, 1]
        phi_vals = [t["osmotic_coefficient"] for t in thermo_results]
        ax.plot(concs, phi_vals, 's-', color='#FF5722', linewidth=2, markersize=8)
        ax.set_xlabel("Concentration (mol/L)", fontsize=11)
        ax.set_ylabel("Osmotic Coefficient φ", fontsize=11)
        ax.set_title("(D) Osmotic Coefficient", fontsize=12)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)

        plt.tight_layout()
        plt.savefig("figures/kirkwood_buff_analysis.png", dpi=300, bbox_inches='tight')
        plt.savefig("figures/kirkwood_buff_analysis.svg", bbox_inches='tight')
        plt.close()
        print("Figures saved to figures/kirkwood_buff_analysis.png/.svg")
    except ImportError:
        print("matplotlib not available; skipping figure generation")

    return output


if __name__ == "__main__":
    main()
