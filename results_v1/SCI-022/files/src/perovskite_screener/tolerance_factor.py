"""
Module 1: Extended Goldschmidt Tolerance Factor & Structural Stability
=======================================================================
Implements:
  - Classical Goldschmidt tolerance factor  t = (rA + rX) / (√2 (rB + rX))
  - Octahedral factor  μ = rB / rX
  - New tolerance factor  τ (Bartel et al. 2019)
  - Stability windows and phase prediction
  - Mixed-halide / double perovskite extensions
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
from .materials_database import get_ionic_radius, IONIC_RADII


@dataclass
class ToleranceResult:
    formula: str
    A: str
    B: str
    X: str
    rA: float
    rB: float
    rX: float
    goldschmidt_t: float
    octahedral_mu: float
    bartel_tau: float
    stability_class: str        # "perovskite", "hexagonal", "ilmenite", "unstable"
    distortion: str             # "cubic", "tetragonal", "orthorhombic", "rhombohedral"
    decomposition_risk: float   # 0–1 (higher = less stable)
    notes: list = field(default_factory=list)


def goldschmidt_tolerance(rA: float, rB: float, rX: float) -> float:
    """Classic Goldschmidt tolerance factor t."""
    return (rA + rX) / (np.sqrt(2) * (rB + rX))


def octahedral_factor(rB: float, rX: float) -> float:
    """Octahedral factor μ = rB/rX. Stable range: 0.414–0.732."""
    return rB / rX


def bartel_tau(rA: float, rB: float, rX: float, nA: float = 1.0) -> float:
    """
    New tolerance factor τ (Bartel et al. Science Adv. 2019).
    τ = rX/rB - nA*(nA - rA/rB / ln(rA/rB))
    nA = oxidation state of A-site (typically 1)
    Stable perovskite: τ < 4.18
    """
    ratio = rA / rB
    tau = rX / rB - nA * (nA - ratio / np.log(ratio))
    return tau


def classify_stability(t: float, mu: float, tau: float) -> tuple:
    """
    Classify structural stability based on tolerance factors.
    Returns (stability_class, distortion, decomposition_risk)
    """
    # Octahedral stability window
    if mu < 0.414:
        return "unstable", "amorphous", 0.95
    if mu > 0.732:
        return "hexagonal", "hexagonal", 0.75

    # Bartel τ-based classification (primary predictor)
    if tau < 3.5:
        return "unstable", "amorphous", 0.90
    elif tau < 4.18:
        # Stable perovskite – determine distortion from t
        if 0.95 <= t <= 1.02:
            return "perovskite", "cubic", 0.05
        elif 0.89 <= t < 0.95:
            return "perovskite", "tetragonal", 0.12
        elif 0.80 <= t < 0.89:
            return "perovskite", "orthorhombic", 0.22
        elif 0.71 <= t < 0.80:
            return "perovskite", "rhombohedral", 0.35
        else:
            return "ilmenite", "rhombohedral", 0.55
    elif tau < 6.0:
        if t > 1.02:
            return "hexagonal", "hexagonal", 0.65
        else:
            return "ilmenite", "trigonal", 0.60
    else:
        return "unstable", "amorphous", 0.90


def estimate_decomposition_enthalpy(A: str, B: str, X: str, B_ox: int) -> float:
    """
    Empirical estimate of decomposition enthalpy ΔH_decomp (eV/formula unit).
    Negative = thermodynamically stable. Uses electronegativity differences.
    Based on: Brgoch et al. parametrization, corrected for Sn/Ge/Bi.
    """
    from .materials_database import ELECTRONEGATIVITY
    try:
        chi_A = ELECTRONEGATIVITY.get(A, 1.5)
        chi_B = ELECTRONEGATIVITY.get(B, 2.0)
        chi_X = ELECTRONEGATIVITY.get(X, 2.7)

        # Ionicity of B-X bond
        delta_BX = abs(chi_B - chi_X)
        # Ionicity of A-X bond
        delta_AX = abs(chi_A - chi_X)

        # Empirical formula calibrated to DFT hull distances
        # Sn penalty (oxidation instability)
        sn_penalty = 0.15 if B == "Sn" else 0.0
        ge_penalty = 0.10 if B == "Ge" else 0.0
        bi_factor  = 0.08 if B == "Bi" else 0.0  # Bi more stable vs decomp

        dH = -0.45 * delta_BX - 0.20 * delta_AX + sn_penalty + ge_penalty - bi_factor
        return round(dH, 3)
    except Exception:
        return 0.0


def analyze_perovskite(A: str, B: str, X: str, B_ox: int = 2) -> ToleranceResult:
    """
    Full tolerance factor analysis for ABX3 perovskite.
    """
    # A-site: 12-coordinate; B-site: 6-coordinate; X-site: 6-coordinate
    try:
        rA = get_ionic_radius(A, 1, cn=12)
    except Exception:
        rA = get_ionic_radius(A, 1, cn=6) * 1.12  # fallback scaling

    rB = get_ionic_radius(B, B_ox, cn=6)
    rX = get_ionic_radius(X, -1, cn=6)

    t   = goldschmidt_tolerance(rA, rB, rX)
    mu  = octahedral_factor(rB, rX)
    tau = bartel_tau(rA, rB, rX, nA=1.0)

    stab_class, distortion, decomp_risk = classify_stability(t, mu, tau)

    notes = []
    if B == "Sn":
        notes.append("Sn²⁺→Sn⁴⁺ oxidation risk; requires reducing atmosphere or Sn-excess")
    if B == "Ge":
        notes.append("Ge²⁺ susceptible to oxidation and moisture; consider encapsulation")
    if B == "Bi":
        notes.append("Bi³⁺ forms layered A₃Bi₂X₉; indirect band gap may limit Jsc")
    if t > 1.05:
        notes.append("Over-tolerance: likely hexagonal or 2D phases")
    if tau > 5.5:
        notes.append("τ > 5.5: high decomposition probability on hull")

    formula = f"{A}{B}{X}3"
    return ToleranceResult(
        formula=formula, A=A, B=B, X=X,
        rA=round(rA, 4), rB=round(rB, 4), rX=round(rX, 4),
        goldschmidt_t=round(t, 4),
        octahedral_mu=round(mu, 4),
        bartel_tau=round(tau, 4),
        stability_class=stab_class,
        distortion=distortion,
        decomposition_risk=round(decomp_risk, 3),
        notes=notes,
    )


def mixed_halide_tolerance(A: str, B: str, halides: dict, B_ox: int = 2) -> dict:
    """
    Tolerance factors for mixed-halide ABX₃₋ₓYₓ perovskites.
    halides: {"I": 0.7, "Br": 0.3} – mole fractions summing to 1.
    """
    assert abs(sum(halides.values()) - 1.0) < 1e-6, "Halide fractions must sum to 1"
    rX_mix = sum(get_ionic_radius(h, -1, cn=6) * frac for h, frac in halides.items())
    try:
        rA = get_ionic_radius(A, 1, cn=12)
    except Exception:
        rA = get_ionic_radius(A, 1, cn=6) * 1.12
    rB = get_ionic_radius(B, B_ox, cn=6)

    t   = goldschmidt_tolerance(rA, rB, rX_mix)
    mu  = octahedral_factor(rB, rX_mix)
    tau = bartel_tau(rA, rB, rX_mix)
    stab_class, distortion, decomp_risk = classify_stability(t, mu, tau)

    halide_str = "".join(f"{h}{frac:.1f}" for h, frac in halides.items())
    return {
        "formula": f"{A}{B}({halide_str})3",
        "rX_eff": round(rX_mix, 4),
        "goldschmidt_t": round(t, 4),
        "octahedral_mu": round(mu, 4),
        "bartel_tau": round(tau, 4),
        "stability_class": stab_class,
        "distortion": distortion,
        "decomposition_risk": round(decomp_risk, 3),
    }
