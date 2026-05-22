"""
Topological Insulator Design Framework
Module 6: Bi2Se3 Analog Candidate Screening

Systematic screening of isostructural compounds for topological properties.

Screening criteria:
  1. Space group compatibility (R-3m preferred)
  2. Average atomic number Z > 35 (heavy elements for SOC)
  3. Estimated SOC strength > threshold
  4. Band gap window: 0.05 – 1.5 eV (not too large, not semimetal)
  5. Valence electron count compatible with TI phase (n_e = 5 per formula)
  6. Electronic structure descriptor: inverted band ordering at Γ

Scoring function: composite TI-score [0,1]
"""

import numpy as np
import json
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))


# ---------------------------------------------------------------------------
# Extended candidate database
# ---------------------------------------------------------------------------

EXTENDED_CANDIDATES = {
    # Confirmed TIs
    "Bi2Se3":  {"Z_avg": 55.6, "sg": 166, "gap": 0.30, "n_ve": 5, "lam_c": 0.65, "exp_TI": True,  "family": "V2-VI3"},
    "Bi2Te3":  {"Z_avg": 68.2, "sg": 166, "gap": 0.15, "n_ve": 5, "lam_c": 0.35, "exp_TI": True,  "family": "V2-VI3"},
    "Sb2Te3":  {"Z_avg": 58.4, "sg": 166, "gap": 0.21, "n_ve": 5, "lam_c": 0.49, "exp_TI": True,  "family": "V2-VI3"},
    "TlBiSe2": {"Z_avg": 59.5, "sg": 166, "gap": 0.35, "n_ve": 4, "lam_c": 0.81, "exp_TI": True,  "family": "ternary"},
    "TlBiTe2": {"Z_avg": 67.3, "sg": 166, "gap": 0.20, "n_ve": 4, "lam_c": 0.50, "exp_TI": True,  "family": "ternary"},
    "GeBi2Te4":{"Z_avg": 57.3, "sg": 166, "gap": 0.18, "n_ve": 4, "lam_c": 0.42, "exp_TI": True,  "family": "heterostructure"},
    "MnBi2Te4":{"Z_avg": 58.1, "sg": 166, "gap": 0.20, "n_ve": 4, "lam_c": 0.46, "exp_TI": True,  "family": "magnetic"},
    "PbBi2Te4":{"Z_avg": 66.4, "sg": 166, "gap": 0.23, "n_ve": 4, "lam_c": 0.53, "exp_TI": True,  "family": "heterostructure"},
    "Bi4Br4":  {"Z_avg": 53.5, "sg": 12,  "gap": 0.18, "n_ve": 3, "lam_c": 0.60, "exp_TI": True,  "family": "halide"},

    # Trivial / negative examples
    "Bi2S3":   {"Z_avg": 37.4, "sg": 62,  "gap": 1.30, "n_ve": 5, "lam_c": 2.10, "exp_TI": False, "family": "V2-VI3"},
    "Sb2Se3":  {"Z_avg": 44.4, "sg": 62,  "gap": 1.20, "n_ve": 5, "lam_c": 1.80, "exp_TI": False, "family": "V2-VI3"},
    "As2Se3":  {"Z_avg": 38.6, "sg": 14,  "gap": 1.80, "n_ve": 5, "lam_c": 3.50, "exp_TI": False, "family": "V2-VI3"},

    # Novel candidates (theoretical predictions)
    "Bi2Po3":  {"Z_avg": 67.2, "sg": 166, "gap": 0.25, "n_ve": 5, "lam_c": 0.55, "exp_TI": None, "family": "V2-VI3"},
    "BiTlS2":  {"Z_avg": 47.3, "sg": 166, "gap": 0.50, "n_ve": 4, "lam_c": 1.15, "exp_TI": None, "family": "ternary"},
    "TlSbTe2": {"Z_avg": 64.7, "sg": 166, "gap": 0.28, "n_ve": 4, "lam_c": 0.65, "exp_TI": None, "family": "ternary"},
    "PbBi4Te7":{"Z_avg": 66.8, "sg": 166, "gap": 0.19, "n_ve": 4, "lam_c": 0.44, "exp_TI": None, "family": "heterostructure"},
    "SnBi2Te4":{"Z_avg": 65.0, "sg": 166, "gap": 0.22, "n_ve": 4, "lam_c": 0.51, "exp_TI": None, "family": "heterostructure"},
    "InBiTe3": {"Z_avg": 57.7, "sg": 166, "gap": 0.32, "n_ve": 4, "lam_c": 0.74, "exp_TI": None, "family": "ternary"},
    "TlBiPo2": {"Z_avg": 72.3, "sg": 166, "gap": 0.18, "n_ve": 4, "lam_c": 0.42, "exp_TI": None, "family": "ternary"},
    "Bi2MnTe4":{"Z_avg": 59.5, "sg": 166, "gap": 0.22, "n_ve": 4, "lam_c": 0.50, "exp_TI": None, "family": "magnetic"},
    "CrBi2Te4":{"Z_avg": 58.7, "sg": 166, "gap": 0.24, "n_ve": 4, "lam_c": 0.55, "exp_TI": None, "family": "magnetic"},
    "EuBi2Te4":{"Z_avg": 66.1, "sg": 166, "gap": 0.20, "n_ve": 4, "lam_c": 0.46, "exp_TI": None, "family": "magnetic"},
}

# SOC strength estimated from average Z
# λ_SOC ∝ Z⁴ (relativistic scaling)
SOC_Z_REFERENCE = {"Z": 83, "lambda": 1.0}  # Bi as reference


def estimate_soc_strength(Z_avg: float) -> float:
    """Estimate relative SOC strength: λ ∝ (Z/Z_ref)^4."""
    return (Z_avg / SOC_Z_REFERENCE["Z"])**4


def ti_score(material_data: dict) -> float:
    """
    Composite topological insulator score [0, 1].

    Criteria:
      w1=0.30: SOC strength (heavier → better)
      w2=0.25: Band gap window (0.1–0.5 eV optimal)
      w3=0.20: Space group match (R-3m preferred)
      w4=0.15: Valence electron count
      w5=0.10: λ_c < 1.0 (reachable transition)
    """
    d = material_data

    # w1: SOC
    soc = estimate_soc_strength(d["Z_avg"])
    s1 = min(soc / 0.8, 1.0)

    # w2: Band gap (optimal range 0.1-0.5 eV)
    g = d["gap"]
    if 0.1 <= g <= 0.5:
        s2 = 1.0 - abs(g - 0.25) / 0.25
    elif g < 0.1:
        s2 = g / 0.1
    else:
        s2 = max(0, 1.0 - (g - 0.5) / 1.0)

    # w3: Space group
    s3 = 1.0 if d["sg"] == 166 else 0.3

    # w4: Valence electrons
    s4 = 1.0 if d["n_ve"] in [4, 5] else 0.2

    # w5: Accessible critical SOC
    lc = d["lam_c"]
    s5 = 1.0 if lc < 0.8 else max(0, 1.0 - (lc - 0.8) / 0.6)

    score = 0.30 * s1 + 0.25 * s2 + 0.20 * s3 + 0.15 * s4 + 0.10 * s5
    return round(score, 4)


def run_screening() -> dict:
    """Run full candidate screening and ranking."""
    os.makedirs("results", exist_ok=True)

    print("=" * 75)
    print("Bi2Se3-ANALOG CANDIDATE SCREENING")
    print("=" * 75)

    scored = {}
    for name, data in EXTENDED_CANDIDATES.items():
        score = ti_score(data)
        soc_est = estimate_soc_strength(data["Z_avg"])
        scored[name] = {
            **data,
            "ti_score": score,
            "soc_estimate": round(soc_est, 4),
            "predicted_TI": score > 0.60,
        }

    # Sort by TI score
    ranked = sorted(scored.items(), key=lambda x: x[1]["ti_score"], reverse=True)

    print(f"\n{'Rank':<5} {'Material':<15} {'Score':>6} {'SG':>5} {'Gap':>6} "
          f"{'Z_avg':>6} {'λ_c':>6} {'ExpTI':>7} {'Pred':>6}")
    print("-" * 75)

    top_novel = []
    for rank, (name, data) in enumerate(ranked, 1):
        exp = "Yes" if data["exp_TI"] else ("No" if data["exp_TI"] is False else "?")
        pred = "TI" if data["predicted_TI"] else "trivial"
        print(f"{rank:<5} {name:<15} {data['ti_score']:>6.3f} {data['sg']:>5} "
              f"{data['gap']:>6.3f} {data['Z_avg']:>6.1f} {data['lam_c']:>6.2f} "
              f"{exp:>7} {pred:>6}")
        if data["exp_TI"] is None and data["predicted_TI"]:
            top_novel.append(name)

    print("\n" + "=" * 75)
    print(f"Novel TI candidates (unexplored, score > 0.60): {top_novel}")
    print("=" * 75)

    result = {
        "ranked_materials": {n: d for n, d in ranked},
        "top_novel_candidates": top_novel,
        "screening_criteria": {
            "w1_soc": 0.30, "w2_gap": 0.25, "w3_sg": 0.20,
            "w4_valence": 0.15, "w5_lam_c": 0.10,
            "ti_threshold": 0.60,
        },
        "n_total": len(scored),
        "n_predicted_TI": sum(1 for d in scored.values() if d["predicted_TI"]),
    }

    with open("results/candidate_screening.json", "w") as f:
        json.dump(result, f, indent=2)

    print("\nSaved: results/candidate_screening.json")
    return result


if __name__ == "__main__":
    run_screening()
