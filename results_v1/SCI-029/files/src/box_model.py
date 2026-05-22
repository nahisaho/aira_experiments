"""
Atmospheric Box Model for SOA Formation
Solves coupled ODE system for gas-phase chemistry + partitioning.
Uses explicit Euler + adaptive Runge-Kutta 4(5) integration.
"""
import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Dict, List, Callable, Tuple
import logging

logger = logging.getLogger(__name__)

# ── Physical constants ────────────────────────────────────────────────────────
AVOGADRO     = 6.022e23
R_GAS        = 8.314      # J/(mol·K)
M_AIR        = 2.46e19    # molecules/cm3 at 298K, 1 atm


@dataclass
class AtmosphericConditions:
    """Urban atmosphere initial/boundary conditions."""
    T: float       = 298.15   # K
    P: float       = 101325.0 # Pa
    RH: float      = 0.50     # relative humidity
    NOx_ppb: float = 5.0      # ppbv NOx
    O3_ppb:  float = 30.0     # ppbv O3
    SO2_ppb: float = 1.0      # ppbv SO2
    JNO2:    float = 8.0e-3   # s-1  photolysis rate (midday)
    JO3:     float = 3.0e-5   # s-1  O3 + hv -> O(1D)
    emission_time: float = 0.0  # s


@dataclass
class VOCEmission:
    name: str
    conc_ppb: float    # initial ppbv


# ── Helper: ppb to molec/cm3 ─────────────────────────────────────────────────
def ppb2molec(ppb: float, M: float = M_AIR) -> float:
    return ppb * 1e-9 * M


def molec2ppb(molec: float, M: float = M_AIR) -> float:
    return molec / M * 1e9


# ── Rate constants (298 K, cm3 molecule-1 s-1 unless noted) ──────────────────
# Source: MCM v3.3.1, Atkinson et al. (2006)
def k_arrhenius(A: float, Ea_R: float, T: float) -> float:
    """Arrhenius rate constant."""
    return A * np.exp(-Ea_R / T)


def k_troe(k0: float, kinf: float, Fc: float, M: float) -> float:
    """Troe fall-off expression."""
    kf = k0 * M / (1 + k0 * M / kinf)
    N  = 0.75 - 1.27 * np.log10(Fc)
    f  = Fc ** (1 / (1 + (np.log10(k0 * M / kinf) / N) ** 2))
    return kf * f


class SimplifiedSOABoxModel:
    """
    Simplified atmospheric box model tracking:
    - OH, O3, NO, NO2, HO2 (fast photochemistry)
    - Parent VOC (alpha-pinene as representative)
    - 4 lumped product classes: ELVOC, LVOC, SVOC, IVOC
    - SOA mass via VBS partitioning

    State vector: [OH, O3, NO, NO2, HO2, VOC, ELVOC_g, LVOC_g, SVOC_g, IVOC_g,
                   ELVOC_p, LVOC_p, SVOC_p, IVOC_p]
    Units: all in molecules cm-3 except SOA mass [μg m-3]
    """

    NAMES = [
        "OH", "O3", "NO", "NO2", "HO2",
        "VOC",
        "ELVOC_g", "LVOC_g", "SVOC_g", "IVOC_g",      # gas-phase products
        "ELVOC_p", "LVOC_p", "SVOC_p", "IVOC_p",      # particle-phase products
    ]

    # VBS bin properties: [log10(C*), MW, yield_OH, yield_O3]
    VBS_BINS = {
        "ELVOC": {"log_Cstar": -3.0, "MW": 250.0, "y_OH": 0.12, "y_O3": 0.15},
        "LVOC":  {"log_Cstar": -1.0, "MW": 220.0, "y_OH": 0.18, "y_O3": 0.20},
        "SVOC":  {"log_Cstar":  1.0, "MW": 190.0, "y_OH": 0.15, "y_O3": 0.10},
        "IVOC":  {"log_Cstar":  3.0, "MW": 160.0, "y_OH": 0.08, "y_O3": 0.03},
    }

    def __init__(self, cond: AtmosphericConditions, voc: VOCEmission):
        self.cond = cond
        self.voc  = voc
        self.M    = M_AIR * cond.P / 101325 * 298.15 / cond.T

    def _calc_rates(self, T: float) -> Dict[str, float]:
        """Temperature-dependent rate constants."""
        return {
            # OH + NO2 → HNO3
            "k_OH_NO2": k_troe(
                k0=2.5e-30 * (T / 300) ** -4.4 * self.M,
                kinf=1.6e-11,
                Fc=0.3, M=self.M
            ),
            # OH + NO → HONO
            "k_OH_NO":  k_troe(
                k0=7.0e-31 * (T / 300) ** -2.6 * self.M,
                kinf=3.6e-11,
                Fc=0.6, M=self.M
            ),
            # HO2 + NO → OH + NO2
            "k_HO2_NO":  8.1e-12 * np.exp(270 / T),
            # O3 + NO → NO2 + O2
            "k_O3_NO":   1.9e-14 * np.exp(-1310 / T),
            # OH + VOC (alpha-pinene)
            "k_OH_VOC":  1.2e-11 * np.exp(440 / T),
            # O3 + VOC
            "k_O3_VOC":  6.3e-16 * np.exp(-580 / T),
            # OH source from O(1D) + H2O
            "S_OH":      2 * self.cond.JO3 * ppb2molec(self.cond.O3_ppb, self.M)
                         * 2.2e-10 * ppb2molec(0.47 * self.cond.RH * 33.0, self.M),  # [H2O]
        }

    def _partitioning_coeff(self, log_Cstar: float, Coa: float, T: float) -> float:
        """Equilibrium particle-phase fraction from VBS."""
        Cstar = 10 ** log_Cstar * np.exp(10e3 / 8.314e-3 * (1 / 298.15 - 1 / T))
        return Coa / (Coa + Cstar) if Coa > 0 else 0.0

    def _molec_to_ugm3(self, molec_cm3: float, MW: float) -> float:
        """Convert molec/cm3 → μg/m3."""
        return molec_cm3 / AVOGADRO * MW * 1e6 * 1e6  # g/mol * 1e6 cm3/m3 * 1e6 μg/g

    def _ugm3_to_molec(self, ugm3: float, MW: float) -> float:
        return ugm3 / 1e12 * AVOGADRO / MW

    def odes(self, t: float, y: np.ndarray) -> np.ndarray:
        """RHS of the ODE system."""
        (OH, O3, NO, NO2, HO2, VOC,
         ELVOC_g, LVOC_g, SVOC_g, IVOC_g,
         ELVOC_p, LVOC_p, SVOC_p, IVOC_p) = np.maximum(y, 0.0)

        T   = self.cond.T
        JNO2 = self.cond.JNO2
        JO3  = self.cond.JO3
        kr   = self._calc_rates(T)

        # Total SOA mass for partitioning
        bin_names = ["ELVOC", "LVOC", "SVOC", "IVOC"]
        g_phases  = [ELVOC_g, LVOC_g, SVOC_g, IVOC_g]
        p_phases  = [ELVOC_p, LVOC_p, SVOC_p, IVOC_p]
        MWs       = [self.VBS_BINS[b]["MW"] for b in bin_names]

        Coa = sum(
            self._molec_to_ugm3(p, MW)
            for p, MW in zip(p_phases, MWs)
        ) + 1.0  # seed aerosol [μg m-3]

        # O3 photolysis → O(1D) + O2 (→ OH)
        dO3_phot = -JO3 * O3

        # NOx photochemistry
        dNO2_phot = -JNO2 * NO2           # NO2 + hv → NO + O
        dNO_phot  = +JNO2 * NO2
        dO3_prod  = +JNO2 * NO2 * 0.95    # O + O2 + M → O3

        # Bimolecular
        r_OH_NO2  = kr["k_OH_NO2"] * OH * NO2
        r_OH_NO   = kr["k_OH_NO"]  * OH * NO
        r_HO2_NO  = kr["k_HO2_NO"] * HO2 * NO
        r_O3_NO   = kr["k_O3_NO"]  * O3  * NO

        # VOC oxidation
        r_OH_VOC  = kr["k_OH_VOC"] * OH  * VOC
        r_O3_VOC  = kr["k_O3_VOC"] * O3  * VOC

        r_VOC_tot = r_OH_VOC + r_O3_VOC

        # Product formation and partitioning
        dg, dp = [], []
        for b in bin_names:
            props  = self.VBS_BINS[b]
            MW_b   = props["MW"]
            # Production rate [molec cm-3 s-1] from OH and O3 pathways
            prod_rate = (
                r_OH_VOC * props["y_OH"] +
                r_O3_VOC * props["y_O3"]
            )
            Cstar_b    = 10 ** props["log_Cstar"]
            Fpart_b    = self._partitioning_coeff(props["log_Cstar"], Coa, T)
            # Condensation rate constant [s-1]
            k_cond     = 0.1  # s-1 (mass transfer limited)
            # Gas-phase production minus condensation
            g_idx  = ["ELVOC", "LVOC", "SVOC", "IVOC"].index(b)
            g_conc = g_phases[g_idx]
            p_conc = p_phases[g_idx]
            p_conc_eq  = g_conc * Fpart_b / max(1 - Fpart_b, 1e-6)
            dg.append(prod_rate - k_cond * (g_conc - (1 - Fpart_b) * (g_conc + p_conc)))
            dp.append(             k_cond * (g_conc - (1 - Fpart_b) * (g_conc + p_conc)))

        # Differential equations
        dOH   = (+kr["S_OH"]
                 - r_OH_NO2 - r_OH_NO - r_OH_VOC
                 + r_HO2_NO * 0.9
                 + 2 * JO3 * O3 * 0.1)

        dO3   = (dO3_prod + dO3_phot - r_O3_NO - r_O3_VOC)

        dNO   = (dNO_phot + r_HO2_NO + r_O3_NO
                 - r_OH_NO - r_O3_NO * 1.0
                 + JNO2 * NO2)

        dNO2  = (-JNO2 * NO2
                 - r_OH_NO2
                 + r_O3_NO
                 - r_HO2_NO)

        dHO2  = (r_OH_VOC * 0.5 - r_HO2_NO + 1e3)   # simplified HO2 source

        dVOC  = -r_VOC_tot

        return np.array([
            dOH, dO3, dNO, dNO2, dHO2, dVOC,
            *dg, *dp,
        ])

    def run(self, t_end: float = 3600.0 * 8, n_points: int = 500) -> Dict:
        """Integrate the box model over t_end seconds."""
        t_span  = (0.0, t_end)
        t_eval  = np.linspace(0, t_end, n_points)

        # Initial conditions
        y0 = np.array([
            ppb2molec(2e-4, self.M),          # OH  ~ 0.0002 ppt (elevated urban)
            ppb2molec(self.cond.O3_ppb, self.M),
            ppb2molec(self.cond.NOx_ppb * 0.6, self.M),
            ppb2molec(self.cond.NOx_ppb * 0.4, self.M),
            ppb2molec(0.05, self.M),           # HO2
            ppb2molec(self.voc.conc_ppb, self.M),
            0.0, 0.0, 0.0, 0.0,               # gas-phase products
            0.0, 0.0, 0.0, 0.0,               # particle-phase products
        ])

        sol = solve_ivp(
            self.odes,
            t_span,
            y0,
            method="RK45",
            t_eval=t_eval,
            rtol=1e-6,
            atol=1e-12,
            max_step=60.0,
        )

        if not sol.success:
            logger.warning(f"ODE solver warning: {sol.message}")

        # Compute SOA mass [μg m-3]
        bin_names = ["ELVOC", "LVOC", "SVOC", "IVOC"]
        SOA = np.zeros_like(sol.t)
        for i, b in enumerate(bin_names):
            MW_b = self.VBS_BINS[b]["MW"]
            SOA += np.array([
                self._molec_to_ugm3(sol.y[10 + i, j], MW_b)
                for j in range(len(sol.t))
            ])

        t_h = sol.t / 3600.0

        return {
            "t":       sol.t,
            "t_h":     t_h,
            "y":       sol.y,
            "names":   self.NAMES,
            "SOA":     SOA,
            "OH":      molec2ppb(sol.y[0]),
            "O3":      molec2ppb(sol.y[1]),
            "NO":      molec2ppb(sol.y[2]),
            "NO2":     molec2ppb(sol.y[3]),
            "VOC":     molec2ppb(sol.y[5]),
            "VOC0":    self.voc.conc_ppb,
            "success": sol.success,
        }


def run_multiple_scenarios(
    voc_list: List[VOCEmission],
    cond: AtmosphericConditions,
    t_end: float = 3600 * 8,
) -> Dict[str, Dict]:
    """Run box model for multiple VOC scenarios."""
    results = {}
    for voc in voc_list:
        logger.info(f"Running box model: {voc.name} @ {voc.conc_ppb} ppb")
        model = SimplifiedSOABoxModel(cond, voc)
        res   = model.run(t_end=t_end)
        results[voc.name] = res
    return results
