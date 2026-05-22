"""
SOA Yield Prediction for Terpene/Isoprene Systems
Implements:
  - Literature-based yield parameterizations
  - Two-product model (Odum et al. 1996)
  - VBS-based yield prediction
  - NOx dependence modeling
"""
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple
import logging

logger = logging.getLogger(__name__)


@dataclass
class SOAYieldResult:
    voc: str
    oxidant: str
    NOx_regime: str       # "low" or "high"
    T: float              # K
    Coa: float            # μg/m3 OA loading
    Y_2prod: float        # Two-product model yield
    Y_VBS: float          # VBS-based yield
    Y_literature: float   # Literature reference yield
    Y_predicted: float    # Best estimate (weighted average)
    uncertainty: float    # ± fractional uncertainty


# ── Two-product model parameters (Odum 1996) ─────────────────────────────────
# alpha_i: stoichiometric mass yield
# K_om_i:  partitioning coefficient [m3/μg]
# Source: Griffin et al. (1999), Ng et al. (2007), Presto & Donahue (2006)
TWO_PROD_PARAMS = {
    # (VOC, oxidant, NOx_regime): [(alpha1, Kom1), (alpha2, Kom2)]
    ("alpha_pinene", "OH",  "low"):  [(0.232, 0.0357), (0.038, 1.621)],
    ("alpha_pinene", "OH",  "high"): [(0.091, 0.0609), (0.154, 0.0538)],
    ("alpha_pinene", "O3",  "any"):  [(0.300, 0.0420), (0.049, 2.140)],
    ("beta_pinene",  "OH",  "low"):  [(0.117, 0.0832), (0.013, 2.820)],
    ("beta_pinene",  "OH",  "high"): [(0.055, 0.0681), (0.097, 0.0464)],
    ("beta_pinene",  "O3",  "any"):  [(0.100, 0.0550), (0.018, 3.100)],
    ("limonene",     "OH",  "low"):  [(0.239, 0.0500), (0.363, 0.0180)],
    ("limonene",     "OH",  "high"): [(0.150, 0.0750), (0.250, 0.0300)],
    ("limonene",     "O3",  "any"):  [(0.400, 0.0300), (0.200, 0.0100)],
    ("isoprene",     "OH",  "low"):  [(0.050, 0.0600), (0.007, 5.000)],
    ("isoprene",     "OH",  "high"): [(0.015, 0.0500), (0.035, 0.0200)],
    ("isoprene",     "O3",  "any"):  [(0.018, 0.1200), (0.002, 8.000)],
    ("toluene",      "OH",  "low"):  [(0.360, 0.0320), (0.136, 0.0053)],
    ("toluene",      "OH",  "high"): [(0.074, 0.0280), (0.217, 0.0058)],
}

# ── VBS parameters: volatility distribution of products ──────────────────────
# bins: log10(C*) → mass fraction in each bin
VBS_YIELDS = {
    # (VOC, oxidant, NOx_regime): {log_Cstar_bin: yield_fraction}
    ("alpha_pinene", "OH",  "low"):  {-3: 0.05, -1: 0.12, 1: 0.15, 3: 0.08},
    ("alpha_pinene", "OH",  "high"): {-3: 0.02, -1: 0.08, 1: 0.10, 3: 0.05},
    ("alpha_pinene", "O3",  "any"):  {-3: 0.08, -1: 0.15, 1: 0.12, 3: 0.05},
    ("beta_pinene",  "OH",  "low"):  {-3: 0.03, -1: 0.08, 1: 0.10, 3: 0.06},
    ("beta_pinene",  "O3",  "any"):  {-3: 0.04, -1: 0.09, 1: 0.08, 3: 0.04},
    ("limonene",     "OH",  "low"):  {-3: 0.07, -1: 0.18, 1: 0.18, 3: 0.10},
    ("limonene",     "O3",  "any"):  {-3: 0.10, -1: 0.22, 1: 0.15, 3: 0.08},
    ("isoprene",     "OH",  "low"):  {-3: 0.00, -1: 0.01, 1: 0.04, 3: 0.02},
    ("isoprene",     "OH",  "high"): {-3: 0.00, -1: 0.01, 1: 0.02, 3: 0.01},
    ("toluene",      "OH",  "low"):  {-3: 0.03, -1: 0.12, 1: 0.20, 3: 0.08},
    ("toluene",      "OH",  "high"): {-3: 0.01, -1: 0.05, 1: 0.10, 3: 0.04},
}

# ── Literature reference yields at Coa = 10 μg/m3, 298 K ─────────────────────
LITERATURE_YIELDS = {
    ("alpha_pinene", "OH",  "low"):  0.30,  # Presto & Donahue 2006
    ("alpha_pinene", "OH",  "high"): 0.14,
    ("alpha_pinene", "O3",  "any"):  0.40,  # Pathak et al. 2007
    ("beta_pinene",  "OH",  "low"):  0.15,  # Griffin et al. 1999
    ("beta_pinene",  "OH",  "high"): 0.09,
    ("beta_pinene",  "O3",  "any"):  0.13,
    ("limonene",     "OH",  "low"):  0.39,  # Ng et al. 2006
    ("limonene",     "OH",  "high"): 0.25,
    ("limonene",     "O3",  "any"):  0.50,  # Saathoff et al. 2009
    ("isoprene",     "OH",  "low"):  0.03,  # Kroll et al. 2006
    ("isoprene",     "OH",  "high"): 0.02,
    ("isoprene",     "O3",  "any"):  0.01,
    ("toluene",      "OH",  "low"):  0.28,  # Ng et al. 2007
    ("toluene",      "OH",  "high"): 0.12,
}


def two_product_yield(
    voc: str,
    oxidant: str,
    NOx_regime: str,
    Coa: float,
) -> float:
    """
    Odum two-product absorptive partitioning yield:
    Y = Coa * sum_i [ alpha_i * Kom_i / (1 + Kom_i * Coa) ]
    """
    key = (voc, oxidant, NOx_regime)
    if key not in TWO_PROD_PARAMS:
        key_any = (voc, oxidant, "any")
        if key_any in TWO_PROD_PARAMS:
            key = key_any
        else:
            return 0.0

    params = TWO_PROD_PARAMS[key]
    yield_val = Coa * sum(
        alpha * Kom / (1.0 + Kom * Coa)
        for alpha, Kom in params
    )
    return float(np.clip(yield_val, 0, 1))


def vbs_yield(
    voc: str,
    oxidant: str,
    NOx_regime: str,
    Coa: float,
    T: float = 298.15,
) -> float:
    """
    VBS-based yield: sum over bins of (alpha_bin * Fpart_bin)
    with temperature correction for Cstar.
    """
    key = (voc, oxidant, NOx_regime)
    if key not in VBS_YIELDS:
        key_any = (voc, oxidant, "any")
        if key_any in VBS_YIELDS:
            key = key_any
        else:
            return 0.0

    bins = VBS_YIELDS[key]
    yield_val = 0.0
    for log_Cstar, alpha_bin in bins.items():
        # Temperature correction: dHvap ~ 80 kJ/mol
        Cstar = (10 ** log_Cstar) * np.exp(80e3 / 8.314 * (1 / 298.15 - 1 / T))
        Fpart = Coa / (Coa + Cstar)
        yield_val += alpha_bin * Fpart
    return float(np.clip(yield_val, 0, 1))


def predict_soa_yield(
    voc: str,
    oxidant: str,
    NOx_ppb: float = 5.0,
    Coa: float = 10.0,
    T: float = 298.15,
) -> SOAYieldResult:
    """Predict SOA yield using combined two-product + VBS models."""
    NOx_regime = "low" if NOx_ppb < 10.0 else "high"

    y_2p  = two_product_yield(voc, oxidant, NOx_regime, Coa)
    y_vbs = vbs_yield(voc, oxidant, NOx_regime, Coa, T)
    y_lit = LITERATURE_YIELDS.get(
        (voc, oxidant, NOx_regime),
        LITERATURE_YIELDS.get((voc, oxidant, "any"), np.nan)
    )

    # Temperature correction for two-product model
    T_correction = np.exp(-3000.0 * (1 / T - 1 / 298.15))
    y_2p_T = y_2p * T_correction

    # Best estimate: weighted average
    weights = [0.4, 0.6]  # 2-product, VBS
    y_pred  = weights[0] * y_2p_T + weights[1] * y_vbs

    # Uncertainty: coefficient of variation between methods
    y_vals = [v for v in [y_2p_T, y_vbs, y_lit] if not np.isnan(v)]
    unc    = float(np.std(y_vals) / max(np.mean(y_vals), 1e-6)) if len(y_vals) > 1 else 0.2

    return SOAYieldResult(
        voc=voc, oxidant=oxidant, NOx_regime=NOx_regime,
        T=T, Coa=Coa,
        Y_2prod=float(y_2p_T), Y_VBS=float(y_vbs),
        Y_literature=float(y_lit) if not np.isnan(y_lit) else -1.0,
        Y_predicted=float(y_pred), uncertainty=unc,
    )


def yield_vs_coa(
    voc: str,
    oxidant: str,
    NOx_ppb: float = 5.0,
    T: float = 298.15,
    Coa_range: np.ndarray = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (Coa, Y_2prod, Y_VBS) over a range of OA loadings."""
    if Coa_range is None:
        Coa_range = np.logspace(-1, 2, 50)

    NOx_regime = "low" if NOx_ppb < 10.0 else "high"
    Y_2p  = np.array([two_product_yield(voc, oxidant, NOx_regime, c) for c in Coa_range])
    Y_vbs = np.array([vbs_yield(voc, oxidant, NOx_regime, c, T) for c in Coa_range])
    return Coa_range, Y_2p, Y_vbs


def yield_temperature_sensitivity(
    voc: str, oxidant: str, NOx_ppb: float = 5.0, Coa: float = 10.0,
    T_range: np.ndarray = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """SOA yield vs temperature."""
    if T_range is None:
        T_range = np.linspace(270, 320, 30)
    NOx_regime = "low" if NOx_ppb < 10.0 else "high"
    Y = np.array([vbs_yield(voc, oxidant, NOx_regime, Coa, T) for T in T_range])
    return T_range, Y


def generate_yield_table(
    voc_list: List[str],
    oxidant_list: List[str] = ["OH", "O3"],
    NOx_ppb: float = 5.0,
    Coa: float = 10.0,
    T: float = 298.15,
) -> List[SOAYieldResult]:
    """Generate yield predictions for all VOC × oxidant combinations."""
    results = []
    for voc in voc_list:
        for oxidant in oxidant_list:
            res = predict_soa_yield(voc, oxidant, NOx_ppb, Coa, T)
            results.append(res)
    return results
