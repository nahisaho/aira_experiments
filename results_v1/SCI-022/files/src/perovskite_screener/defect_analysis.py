"""
Module 3: Defect Formation Energy & Non-Radiative Recombination
===============================================================
Implements:
  - Defect formation energy ΔHf(q, EF) for vacancies, interstitials, antisites
  - Charge transition levels ε(q/q')
  - Shockley-Read-Hall (SRH) recombination rate estimate
  - Non-radiative recombination loss: ΔVoc,nr
  - Defect tolerance classification (defect-tolerant vs. fatal)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from .materials_database import ELECTRONEGATIVITY, get_ionic_radius


# ── Constants ─────────────────────────────────────────────────────────────────
kB = 8.617e-5   # eV/K
T  = 300        # K
kT = kB * T     # ≈ 0.02585 eV


# ── Empirical Defect Parameters ───────────────────────────────────────────────

# Empirical correction energies (eV) for Sn/Ge/Bi perovskites
# Compiled from DFT literature (PBE+SOC, HSE06 where available)
DEFECT_PARAMS = {
    "Sn": {
        "V_Sn":    {"formation_base": 0.15, "charge_states": [-2, -1, 0], "type": "vacancy",     "deep": False},
        "V_I":     {"formation_base": 0.30, "charge_states": [0, 1],      "type": "vacancy",     "deep": False},
        "Sn_i":    {"formation_base": 0.45, "charge_states": [0, 1, 2],   "type": "interstitial","deep": False},
        "I_i":     {"formation_base": 0.55, "charge_states": [-1, 0],     "type": "interstitial","deep": False},
        "Sn_Pb":   {"formation_base": 0.00, "charge_states": [0],         "type": "antisite",    "deep": False},
        "Sn4_Sn2": {"formation_base": 0.20, "charge_states": [2],         "type": "oxidation",   "deep": True},
        "V_MA":    {"formation_base": 0.55, "charge_states": [-1, 0],     "type": "vacancy",     "deep": False},
    },
    "Ge": {
        "V_Ge":    {"formation_base": 0.30, "charge_states": [-2, -1, 0], "type": "vacancy",     "deep": True},
        "V_I":     {"formation_base": 0.35, "charge_states": [0, 1],      "type": "vacancy",     "deep": False},
        "Ge_i":    {"formation_base": 0.60, "charge_states": [0, 1, 2],   "type": "interstitial","deep": False},
        "I_i":     {"formation_base": 0.50, "charge_states": [-1, 0],     "type": "interstitial","deep": False},
        "Ge4_Ge2": {"formation_base": 0.25, "charge_states": [2],         "type": "oxidation",   "deep": True},
    },
    "Bi": {
        "V_Bi":    {"formation_base": 0.55, "charge_states": [-3,-2,-1,0],"type": "vacancy",     "deep": False},
        "V_I":     {"formation_base": 0.40, "charge_states": [0, 1],      "type": "vacancy",     "deep": False},
        "Bi_i":    {"formation_base": 0.80, "charge_states": [0, 1, 2, 3],"type": "interstitial","deep": False},
        "I_Bi":    {"formation_base": 0.65, "charge_states": [-2,-1, 0],  "type": "antisite",    "deep": True},
        "Bi_vac":  {"formation_base": 0.45, "charge_states": [-1, 0],     "type": "vacancy",     "deep": False},
    },
}

# Defect formation energy modifiers for halide choice
HALIDE_MODIFIER = {"I": 0.0, "Br": +0.12, "Cl": +0.25, "F": +0.40}

# A-site modifier (cation size / organic vs inorganic)
ASITE_MODIFIER = {
    "Cs": +0.05, "Rb": +0.03, "MA": 0.0, "FA": -0.05, "EA": -0.08, "DMA": -0.10
}


@dataclass
class DefectResult:
    formula: str
    B_site: str
    defect_name: str
    defect_type: str
    formation_energy_eV: float
    dominant_charge_state: int
    is_deep_trap: bool
    concentration_cm3: float   # at T=300K, EF at midgap
    srh_rate_relative: float   # relative SRH rate (1 = MAPbI3 reference)
    notes: str = ""


@dataclass
class DefectSummary:
    formula: str
    B_site: str
    X_site: str
    defects: List[DefectResult] = field(default_factory=list)
    n_deep_traps: int = 0
    defect_tolerance_score: float = 0.0   # 0 (intolerant) – 1 (tolerant)
    Voc_nr_loss_mV: float = 0.0
    dominant_defect: str = ""


def formation_energy(base: float, charge: int, EF: float, Eg: float,
                     correction: float = 0.0) -> float:
    """
    ΔHf(q, EF) = ΔHf(0) + q·EF + correction
    EF measured from VBM. Typical range: 0 → Eg.
    """
    return base + charge * EF + correction


def equilibrium_concentration(dHf: float, g: int = 1) -> float:
    """
    Boltzmann concentration: N_sites * g * exp(-dHf / kT)
    N_sites ≈ 5e21 cm⁻³ (typical perovskite site density)
    """
    N_sites = 5e21
    c = N_sites * g * np.exp(-max(dHf, 0) / kT)
    return min(c, N_sites)


def srh_recombination_rate(trap_density: float, Eg: float,
                            sigma: float = 1e-14) -> float:
    """
    Relative SRH lifetime. τ_SRH ∝ 1 / (σ·v_th·N_trap)
    sigma: capture cross-section (cm²); v_th ≈ 1e7 cm/s
    Returns rate normalized to MAPbI3 reference (N_trap = 1e15 cm⁻³).
    """
    v_th = 1e7
    tau_inv = sigma * v_th * trap_density
    tau_inv_ref = 1e-14 * 1e7 * 1e15  # reference
    return tau_inv / tau_inv_ref


def voc_nr_loss(srh_rate_relative: float, Eg: float) -> float:
    """
    Non-radiative Voc loss (mV) from SRH recombination.
    ΔVoc,nr ≈ kT * ln(SRH_rate_relative) / q
    Capped at 300 mV.
    """
    if srh_rate_relative <= 0:
        return 0.0
    loss = kT * np.log(max(srh_rate_relative, 1.0)) * 1000  # mV
    return min(round(loss, 1), 300.0)


def analyze_defects(A: str, B: str, X: str, Eg: float, B_ox: int = 2) -> DefectSummary:
    """
    Full defect analysis for ABX3 perovskite.
    """
    params = DEFECT_PARAMS.get(B, DEFECT_PARAMS["Sn"])
    halide_mod = HALIDE_MODIFIER.get(X, 0.0)
    asite_mod  = ASITE_MODIFIER.get(A, 0.0)
    formula    = f"{A}{B}{X}3"

    defects = []
    total_deep_conc = 0.0
    total_srh = 0.0

    EF_midgap = Eg / 2   # Fermi level at midgap (intrinsic)

    for name, info in params.items():
        base = info["formation_base"] + halide_mod + asite_mod
        charges = info["charge_states"]

        # Find dominant (lowest formation energy) charge state at midgap
        min_dHf = 1e10
        dom_q   = 0
        for q in charges:
            dHf = formation_energy(base, q, EF_midgap, Eg)
            if dHf < min_dHf:
                min_dHf = dHf
                dom_q   = q

        min_dHf = max(min_dHf, 0.01)
        conc    = equilibrium_concentration(min_dHf)
        srh_r   = srh_recombination_rate(conc, Eg) if info["deep"] else 0.01
        note    = "deep trap – recombination center" if info["deep"] else "shallow defect"

        if info["deep"]:
            total_deep_conc += conc
        total_srh += srh_r

        defects.append(DefectResult(
            formula=formula, B_site=B, defect_name=name,
            defect_type=info["type"],
            formation_energy_eV=round(min_dHf, 3),
            dominant_charge_state=dom_q,
            is_deep_trap=info["deep"],
            concentration_cm3=round(conc, 3),
            srh_rate_relative=round(srh_r, 4),
            notes=note,
        ))

    # Defect tolerance score
    n_deep = sum(1 for d in defects if d.is_deep_trap and d.concentration_cm3 > 1e14)
    deep_conc_total = max(total_deep_conc, 1.0)
    # Higher score = more tolerant (low deep trap density)
    tol_score = np.exp(-deep_conc_total / 1e16)
    tol_score = round(float(np.clip(tol_score, 0, 1)), 3)

    avg_srh   = total_srh / max(len(defects), 1)
    dvoc_nr   = voc_nr_loss(avg_srh, Eg)

    dom_defect = min(defects, key=lambda d: d.formation_energy_eV).defect_name

    return DefectSummary(
        formula=formula, B_site=B, X_site=X,
        defects=defects,
        n_deep_traps=n_deep,
        defect_tolerance_score=tol_score,
        Voc_nr_loss_mV=dvoc_nr,
        dominant_defect=dom_defect,
    )


def defect_tolerance_classification(summary: DefectSummary) -> str:
    """Classify material as defect tolerant / moderate / intolerant."""
    if summary.defect_tolerance_score >= 0.7 and summary.Voc_nr_loss_mV < 80:
        return "defect-tolerant"
    elif summary.defect_tolerance_score >= 0.4 and summary.Voc_nr_loss_mV < 150:
        return "moderate"
    else:
        return "defect-intolerant"
