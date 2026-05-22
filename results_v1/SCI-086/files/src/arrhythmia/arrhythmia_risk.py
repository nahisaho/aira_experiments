"""
arrhythmia_risk.py
==================
Module 5: Arrhythmia risk assessment via simulation.
Implements re-entry vulnerability analysis, APD dispersion mapping,
and arrhythmia inducibility protocols.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class ArrhythmiaType(Enum):
    VT_MONOMORPHIC = "monomorphic_vt"
    VT_POLYMORPHIC = "polymorphic_vt"
    VF = "ventricular_fibrillation"
    ATRIAL_FLUTTER = "atrial_flutter"
    AF = "atrial_fibrillation"
    AVNRT = "av_nodal_reentry"


@dataclass
class ArrhythmiaRiskScore:
    """Comprehensive arrhythmia risk assessment output."""
    overall_risk: float          # 0-1 normalized risk score
    risk_category: str           # "low", "moderate", "high", "very_high"
    sub_scores: Dict[str, float] = field(default_factory=dict)
    vulnerable_regions: List[int] = field(default_factory=list)
    critical_wavelength: float = 0.0     # mm
    reentry_inducible: bool = False
    dominant_frequency: float = 0.0      # Hz
    details: Dict = field(default_factory=dict)


class ReentryVulnerabilityAnalysis:
    """
    Analyze vulnerability to re-entrant arrhythmias.

    Methods:
    1. Vulnerability window (VW) assessment via S1-S2 protocol
    2. APD restitution curve analysis
    3. Conduction velocity restitution
    4. Wavelength analysis (λ = CV × ERP)
    5. Phase singularity detection
    """

    def __init__(self, cv_long: float = 0.6, cv_trans: float = 0.2,
                 apd_90: float = 280.0):
        self.cv_long = cv_long      # m/s
        self.cv_trans = cv_trans    # m/s
        self.apd_90 = apd_90       # ms

    def compute_apd_restitution(self, di_range: np.ndarray = None
                                  ) -> Dict[str, np.ndarray]:
        """
        Compute APD restitution curve: APD = f(DI).

        Steep restitution (slope > 1) promotes alternans and wave break.
        """
        if di_range is None:
            di_range = np.linspace(20, 600, 100)

        # Mono-exponential restitution model
        APD_max = self.apd_90
        APD_min = 150.0
        tau_restitution = 80.0  # ms

        apd = APD_min + (APD_max - APD_min) * (1 - np.exp(-di_range / tau_restitution))

        # Restitution slope
        slope = (APD_max - APD_min) / tau_restitution * np.exp(-di_range / tau_restitution)

        # Critical DI where slope = 1
        critical_di_idx = np.argmin(np.abs(slope - 1.0))
        critical_di = di_range[critical_di_idx]

        return {
            "di": di_range,
            "apd": apd,
            "slope": slope,
            "max_slope": np.max(slope),
            "critical_di": critical_di,
            "alternans_prone": np.max(slope) > 1.0,
        }

    def compute_vulnerability_window(self, erp: float = 230.0,
                                       substrate_size: float = 30.0
                                       ) -> Dict[str, float]:
        """
        Compute the vulnerability window for re-entry induction.

        VW = time interval during which a premature stimulus can
        initiate unidirectional block and re-entry.
        """
        wavelength = self.cv_long * 1000 * erp / 1000  # Convert to mm

        # Critical circuit length for sustained re-entry
        critical_length = wavelength  # λ = CV × ERP

        # Vulnerability window width
        # Wider VW → easier to induce arrhythmia
        vw_width = max(0, erp - 0.8 * self.apd_90)

        # Can substrate support re-entry?
        reentry_possible = substrate_size > critical_length

        return {
            "wavelength_mm": wavelength,
            "critical_circuit_length_mm": critical_length,
            "vw_width_ms": vw_width,
            "erp_ms": erp,
            "reentry_possible": reentry_possible,
            "substrate_size_mm": substrate_size,
        }

    def s1s2_protocol(self, bcl_s1: float = 600.0,
                       s2_range: Tuple[float, float] = (200.0, 500.0),
                       n_s1_beats: int = 8
                       ) -> Dict[str, any]:
        """
        Simulate S1-S2 programmed stimulation protocol.

        S1: Drive train at fixed BCL
        S2: Premature extra stimulus at varying coupling intervals
        """
        s2_intervals = np.linspace(s2_range[0], s2_range[1], 50)

        responses = []
        for s2_ci in s2_intervals:
            # Check if S2 falls within vulnerable window
            di = s2_ci - self.apd_90
            if di < 0:
                response = "no_capture"
            elif di < 20:
                response = "unidirectional_block"
            elif di < 50:
                response = "slow_conduction"
            else:
                response = "normal_capture"

            responses.append({
                "s2_ci": s2_ci,
                "di": max(0, di),
                "response": response,
            })

        # Find ERP (shortest S2 CI with capture)
        captured = [r for r in responses if r["response"] != "no_capture"]
        erp = captured[0]["s2_ci"] if captured else s2_range[1]

        # Find unidirectional block window
        ub_responses = [r for r in responses if r["response"] == "unidirectional_block"]
        vw_start = ub_responses[0]["s2_ci"] if ub_responses else 0
        vw_end = ub_responses[-1]["s2_ci"] if ub_responses else 0

        return {
            "responses": responses,
            "erp_ms": erp,
            "vw_start_ms": vw_start,
            "vw_end_ms": vw_end,
            "vw_width_ms": vw_end - vw_start,
            "n_s1_beats": n_s1_beats,
            "bcl_s1_ms": bcl_s1,
        }


class APDDispersionAnalyzer:
    """
    Analyze spatial dispersion of repolarization.

    Elevated APD dispersion creates substrate for functional re-entry.
    """

    def compute_dispersion_map(self, apd_field: np.ndarray,
                                 element_coords: np.ndarray
                                 ) -> Dict[str, any]:
        """
        Compute APD dispersion metrics from a spatial APD field.

        apd_field: APD value per mesh element
        element_coords: (n_elements, 3) coordinates
        """
        mean_apd = np.mean(apd_field)
        std_apd = np.std(apd_field)
        range_apd = np.max(apd_field) - np.min(apd_field)

        # Dispersion index
        dispersion_index = std_apd / mean_apd

        # Gradient-based analysis
        gradients = self._compute_apd_gradients(apd_field, element_coords)
        max_gradient = np.max(gradients) if len(gradients) > 0 else 0

        # Identify high-gradient regions (potential re-entry sites)
        gradient_threshold = np.percentile(gradients, 95) if len(gradients) > 0 else 0
        high_gradient_regions = np.where(gradients > gradient_threshold)[0]

        return {
            "mean_apd_ms": mean_apd,
            "std_apd_ms": std_apd,
            "range_apd_ms": range_apd,
            "dispersion_index": dispersion_index,
            "max_gradient_ms_per_mm": max_gradient,
            "high_gradient_regions": high_gradient_regions.tolist(),
            "n_vulnerable_regions": len(high_gradient_regions),
        }

    def _compute_apd_gradients(self, apd_field: np.ndarray,
                                 coords: np.ndarray) -> np.ndarray:
        """Compute spatial gradients of APD field."""
        n = len(apd_field)
        gradients = np.zeros(n)

        for i in range(n):
            dists = np.linalg.norm(coords - coords[i], axis=1)
            neighbors = np.where((dists > 0) & (dists < 5.0))[0]
            if len(neighbors) > 0:
                apd_diffs = np.abs(apd_field[neighbors] - apd_field[i])
                neighbor_dists = dists[neighbors]
                gradients[i] = np.max(apd_diffs / neighbor_dists)

        return gradients


class FibrosisMappingAnalyzer:
    """Analyze fibrosis patterns and their arrhythmogenic potential."""

    def analyze_fibrosis_pattern(self, fibrosis_map: np.ndarray,
                                   element_coords: np.ndarray
                                   ) -> Dict[str, any]:
        """
        Analyze fibrosis distribution for arrhythmia substrate.

        fibrosis_map: 0/1 per element (1 = fibrotic)
        """
        total_elements = len(fibrosis_map)
        n_fibrotic = np.sum(fibrosis_map)
        fibrosis_burden = n_fibrotic / total_elements * 100

        # Classify fibrosis pattern
        if fibrosis_burden < 5:
            pattern = "minimal"
            risk_modifier = 0.2
        elif fibrosis_burden < 15:
            pattern = "patchy"
            risk_modifier = 0.6
        elif fibrosis_burden < 30:
            pattern = "dense_patchy"
            risk_modifier = 1.0
        else:
            pattern = "dense_confluent"
            risk_modifier = 0.8  # Paradoxically lower risk (complete block)

        # Border zone analysis (most arrhythmogenic)
        border_zone = self._compute_border_zone(fibrosis_map, element_coords)

        return {
            "fibrosis_burden_pct": fibrosis_burden,
            "pattern": pattern,
            "risk_modifier": risk_modifier,
            "n_border_zone_elements": len(border_zone),
            "border_zone_fraction": len(border_zone) / total_elements * 100,
        }

    def _compute_border_zone(self, fibrosis_map: np.ndarray,
                               coords: np.ndarray) -> np.ndarray:
        """Identify border zone elements (fibrotic-healthy interface)."""
        border = []
        for i in range(len(fibrosis_map)):
            if fibrosis_map[i] == 0:
                dists = np.linalg.norm(coords - coords[i], axis=1)
                neighbors = np.where((dists > 0) & (dists < 3.0))[0]
                if np.any(fibrosis_map[neighbors] == 1):
                    border.append(i)
        return np.array(border)


class ArrhythmiaRiskAssessor:
    """
    Comprehensive arrhythmia risk assessment integrating multiple factors.

    Combines:
    - Restitution properties
    - APD dispersion
    - Fibrosis substrate
    - Conduction abnormalities
    - Chamber geometry
    """

    def __init__(self):
        self.weights = {
            "restitution": 0.2,
            "dispersion": 0.25,
            "fibrosis": 0.25,
            "conduction": 0.15,
            "geometry": 0.15,
        }

    def assess_risk(self, restitution_data: Dict,
                     dispersion_data: Dict,
                     fibrosis_data: Dict,
                     conduction_data: Dict,
                     geometry_data: Optional[Dict] = None
                     ) -> ArrhythmiaRiskScore:
        """Compute comprehensive arrhythmia risk score."""

        sub_scores = {}

        # 1. Restitution risk
        max_slope = restitution_data.get("max_slope", 0)
        sub_scores["restitution"] = min(1.0, max_slope / 1.5)

        # 2. Dispersion risk
        disp_idx = dispersion_data.get("dispersion_index", 0)
        sub_scores["dispersion"] = min(1.0, disp_idx / 0.15)

        # 3. Fibrosis risk
        sub_scores["fibrosis"] = fibrosis_data.get("risk_modifier", 0)

        # 4. Conduction risk
        vw = conduction_data.get("vw_width_ms", 0)
        sub_scores["conduction"] = min(1.0, vw / 50.0)

        # 5. Geometry risk
        if geometry_data:
            wall_thin = geometry_data.get("min_wall_thickness", 10)
            sub_scores["geometry"] = max(0, 1.0 - wall_thin / 15.0)
        else:
            sub_scores["geometry"] = 0.3

        # Weighted overall score
        overall = sum(self.weights[k] * sub_scores[k] for k in self.weights)

        # Risk category
        if overall < 0.25:
            category = "low"
        elif overall < 0.50:
            category = "moderate"
        elif overall < 0.75:
            category = "high"
        else:
            category = "very_high"

        # Detailed analysis
        reentry_inducible = (conduction_data.get("reentry_possible", False) and
                             sub_scores["restitution"] > 0.5)

        wavelength = conduction_data.get("wavelength_mm", 0)

        return ArrhythmiaRiskScore(
            overall_risk=overall,
            risk_category=category,
            sub_scores=sub_scores,
            vulnerable_regions=dispersion_data.get("high_gradient_regions", []),
            critical_wavelength=wavelength,
            reentry_inducible=reentry_inducible,
            details={
                "restitution_max_slope": max_slope,
                "apd_dispersion_ms": dispersion_data.get("std_apd_ms", 0),
                "fibrosis_burden_pct": fibrosis_data.get("fibrosis_burden_pct", 0),
                "vulnerability_window_ms": vw,
            },
        )
