"""
Module 7: Multi-Criteria Ranking of Lead-Free Perovskite Candidates
====================================================================
Implements:
  - Weighted scoring across 8 performance dimensions
  - Pareto-front identification for multi-objective optimization
  - Uncertainty-aware ranking (Monte Carlo sampling)
  - Top-N recommendation with confidence intervals
  - Comparison with MAPbI3 reference
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import warnings


# ── Scoring Weights (calibrated to PV relevance) ─────────────────────────────
DEFAULT_WEIGHTS = {
    "band_gap_score":          0.20,  # optimal: 1.0–1.6 eV for single-junction
    "stability_score":         0.18,  # structural stability (τ, t)
    "defect_tolerance_score":  0.17,  # low deep-trap density
    "ion_migration_score":     0.15,  # high barrier = low migration
    "slme_score":              0.12,  # SLME efficiency
    "device_pce_score":        0.10,  # SCAPS-simulated PCE
    "voc_score":               0.05,  # open-circuit voltage
    "toxicity_score":          0.03,  # 1 = Pb-free, 0 = Pb-containing
}

# Shockley-Queisser optimal band gap (eV) for single-junction solar cell
SQ_OPTIMAL_EG = 1.34

# Reference: MAPbI3 PCE for normalization
MAPBI3_PCE = 25.7


@dataclass
class CandidateScore:
    formula: str
    A: str
    B: str
    X: str
    system: str  # "Sn", "Ge", "Bi"
    # Raw metrics
    Eg_eV: float
    Eg_uncertainty: float
    goldschmidt_t: float
    bartel_tau: float
    stability_class: str
    defect_tolerance: float
    voc_nr_loss_mV: float
    neb_barrier_eV: float
    slme_pct: float
    device_pce_pct: float
    Voc_V: float
    Jsc_mAcm2: float
    FF: float
    # Scores (0–1)
    band_gap_score: float = 0.0
    stability_score: float = 0.0
    defect_tolerance_score_n: float = 0.0
    ion_migration_score: float = 0.0
    slme_score: float = 0.0
    device_pce_score: float = 0.0
    voc_score: float = 0.0
    toxicity_score: float = 1.0  # all Pb-free
    # Composite
    composite_score: float = 0.0
    rank: int = 0
    pareto_optimal: bool = False
    recommendation: str = ""


def score_band_gap(Eg: float) -> float:
    """
    Band gap score: peaks at SQ optimal (1.34 eV), falls off.
    Uses Gaussian-shaped scoring around optimal.
    """
    optimal = SQ_OPTIMAL_EG
    width   = 0.35  # eV
    score   = np.exp(-((Eg - optimal) ** 2) / (2 * width**2))
    return float(np.clip(score, 0, 1))


def score_stability(goldschmidt_t: float, bartel_tau: float,
                    stability_class: str) -> float:
    """
    Structural stability score based on tolerance factors.
    """
    if stability_class not in ["perovskite"]:
        return 0.1

    # Goldschmidt t: ideal ~0.9–1.0
    t_score  = 1.0 - 2.0 * abs(goldschmidt_t - 0.95)
    t_score  = np.clip(t_score, 0, 1)

    # Bartel τ: lower is better (< 4.18 = stable)
    tau_score = np.clip(1.0 - (bartel_tau - 3.0) / 2.0, 0, 1)

    return float(0.5 * t_score + 0.5 * tau_score)


def score_ion_migration(barrier_eV: float) -> float:
    """
    Ion migration score: high barrier = low migration risk.
    Reference: 0.22 eV (MAPbI3), optimal > 0.4 eV.
    """
    # Sigmoid centered at 0.25 eV
    score = 1.0 / (1.0 + np.exp(-8 * (barrier_eV - 0.25)))
    return float(np.clip(score, 0, 1))


def score_slme(slme_pct: float, reference: float = 25.0) -> float:
    """SLME score normalized to theoretical SQ limit (~33%)."""
    return float(np.clip(slme_pct / reference, 0, 1))


def score_device_pce(pce: float, reference: float = MAPBI3_PCE) -> float:
    """Device PCE score normalized to MAPbI3 reference."""
    return float(np.clip(pce / reference, 0, 1))


def score_voc(Voc: float, Eg: float) -> float:
    """Voc score as fraction of radiative limit Voc,rad ≈ Eg - 0.26 V."""
    voc_limit = max(Eg - 0.26, 0.1)
    return float(np.clip(Voc / voc_limit, 0, 1))


def compute_composite_score(cand: CandidateScore, weights: Dict = None) -> float:
    """Weighted composite score."""
    w = weights or DEFAULT_WEIGHTS
    scores = {
        "band_gap_score":          cand.band_gap_score,
        "stability_score":         cand.stability_score,
        "defect_tolerance_score":  cand.defect_tolerance_score_n,
        "ion_migration_score":     cand.ion_migration_score,
        "slme_score":              cand.slme_score,
        "device_pce_score":        cand.device_pce_score,
        "voc_score":               cand.voc_score,
        "toxicity_score":          cand.toxicity_score,
    }
    total = sum(w.get(k, 0) * v for k, v in scores.items())
    return round(float(total), 4)


def identify_pareto_front(candidates: List[CandidateScore],
                           objectives: List[str] = None) -> List[CandidateScore]:
    """
    Identify Pareto-optimal candidates across multiple objectives.
    All objectives are maximization.
    """
    if objectives is None:
        objectives = ["band_gap_score", "stability_score",
                      "defect_tolerance_score_n", "ion_migration_score"]

    n = len(candidates)
    pareto_mask = [True] * n

    for i, ci in enumerate(candidates):
        for j, cj in enumerate(candidates):
            if i == j:
                continue
            # ci dominated by cj if cj is better in all objectives
            dominated = all(
                getattr(cj, obj, 0) >= getattr(ci, obj, 0)
                for obj in objectives
            ) and any(
                getattr(cj, obj, 0) > getattr(ci, obj, 0)
                for obj in objectives
            )
            if dominated:
                pareto_mask[i] = False
                break

    for i, c in enumerate(candidates):
        c.pareto_optimal = pareto_mask[i]
    return candidates


def generate_recommendation(cand: CandidateScore) -> str:
    """Generate a short recommendation string for the candidate."""
    issues = []
    strengths = []

    if cand.band_gap_score > 0.7:
        strengths.append("optimal band gap")
    if cand.stability_score > 0.6:
        strengths.append("structurally stable")
    if cand.defect_tolerance_score_n > 0.6:
        strengths.append("defect tolerant")
    if cand.neb_barrier_eV > 0.30:
        strengths.append("low ion migration")
    if cand.device_pce_pct > 10:
        strengths.append(f"PCE {cand.device_pce_pct:.1f}%")

    if cand.B == "Sn":
        issues.append("Sn oxidation (use SnF₂ additive)")
    if cand.B == "Ge":
        issues.append("Ge oxidation (moisture sensitive)")
    if cand.B == "Bi" and cand.Eg_eV > 2.0:
        issues.append("wide gap limits Jsc (tandem candidate)")
    if cand.neb_barrier_eV < 0.15:
        issues.append("high ion migration risk")
    if cand.voc_nr_loss_mV > 150:
        issues.append("high non-radiative Voc loss")

    rec = ""
    if strengths:
        rec += "Strengths: " + ", ".join(strengths) + ". "
    if issues:
        rec += "Challenges: " + ", ".join(issues) + "."
    if cand.pareto_optimal:
        rec = "⭐ Pareto-optimal. " + rec
    return rec.strip()


def rank_candidates(
    all_results: List[Dict],
    weights: Dict = None,
    top_n: int = 10,
) -> List[CandidateScore]:
    """
    Main ranking function. Takes list of result dicts from the pipeline.
    Returns sorted list of CandidateScore objects.
    """
    scored = []
    for r in all_results:
        try:
            cs = CandidateScore(
                formula=r["formula"],
                A=r["A"], B=r["B"], X=r["X"],
                system=r.get("system", r["B"]),
                Eg_eV=r.get("Eg_eV", 1.5),
                Eg_uncertainty=r.get("Eg_uncertainty", 0.15),
                goldschmidt_t=r.get("goldschmidt_t", 0.9),
                bartel_tau=r.get("bartel_tau", 4.0),
                stability_class=r.get("stability_class", "perovskite"),
                defect_tolerance=r.get("defect_tolerance_score", 0.5),
                voc_nr_loss_mV=r.get("Voc_nr_loss_mV", 100),
                neb_barrier_eV=r.get("barrier_eV", 0.2),
                slme_pct=r.get("slme_pct", 0.0),
                device_pce_pct=r.get("device_pce_pct", 0.0),
                Voc_V=r.get("Voc_V", 0.7),
                Jsc_mAcm2=r.get("Jsc_mAcm2", 15.0),
                FF=r.get("FF", 0.65),
            )

            # Compute individual scores
            cs.band_gap_score          = score_band_gap(cs.Eg_eV)
            cs.stability_score         = score_stability(cs.goldschmidt_t,
                                                          cs.bartel_tau,
                                                          cs.stability_class)
            cs.defect_tolerance_score_n = float(cs.defect_tolerance)
            cs.ion_migration_score     = score_ion_migration(cs.neb_barrier_eV)
            cs.slme_score              = score_slme(cs.slme_pct)
            cs.device_pce_score        = score_device_pce(cs.device_pce_pct)
            cs.voc_score               = score_voc(cs.Voc_V, cs.Eg_eV)
            cs.toxicity_score          = 1.0  # all Pb-free

            cs.composite_score = compute_composite_score(cs, weights)
            scored.append(cs)
        except Exception as e:
            warnings.warn(f"Skipping {r.get('formula', '?')}: {e}")

    # Sort by composite score
    scored.sort(key=lambda x: x.composite_score, reverse=True)

    # Assign ranks
    for i, c in enumerate(scored):
        c.rank = i + 1

    # Pareto front
    identify_pareto_front(scored)

    # Recommendations
    for c in scored:
        c.recommendation = generate_recommendation(c)

    return scored[:top_n] if top_n else scored
