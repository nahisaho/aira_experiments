"""
Module 3: Geometric descriptor – adsorption relationship analysis.

Analyzes correlations between pore geometry and gas adsorption capacity.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np

logger = logging.getLogger(__name__)


class GeometricAdsorptionAnalyzer:
    """Analyze relationships between geometric descriptors and adsorption."""

    OPTIMAL_RANGES = {
        "CO2_LCD_min": 3.5, "CO2_LCD_max": 12.0,
        "CO2_PLD_min": 3.3, "CO2_PLD_max": 8.0,
        "CO2_ASA_min": 500, "CO2_ASA_max": 3000,
        "CO2_porosity_min": 0.3, "CO2_porosity_max": 0.8,
        "H2_LCD_min": 2.9, "H2_LCD_max": 10.0,
        "H2_PLD_min": 2.89,
        "DAC_LCD_optimal": (5.0, 10.0),
    }

    def compute_correlations(self, descriptors: np.ndarray,
                              loadings: np.ndarray,
                              feature_names: List[str]) -> Dict[str, float]:
        correlations = {}
        for i, name in enumerate(feature_names):
            feat = descriptors[:, i]
            mask = np.isfinite(feat) & np.isfinite(loadings)
            if mask.sum() < 10:
                continue
            r = np.corrcoef(feat[mask], loadings[mask])[0, 1]
            correlations[name] = round(r, 4)
        return dict(sorted(correlations.items(), key=lambda x: abs(x[1]), reverse=True))

    def identify_optimal_windows(self, descriptors: np.ndarray,
                                  loadings: np.ndarray,
                                  feature_names: List[str],
                                  target_percentile: float = 90) -> Dict[str, Tuple]:
        top_mask = loadings >= np.percentile(loadings, target_percentile)
        windows = {}
        for i, name in enumerate(feature_names):
            feat = descriptors[:, i]
            if not np.isfinite(feat).all():
                continue
            top_vals = feat[top_mask]
            if len(top_vals) < 5:
                continue
            windows[name] = (
                round(float(np.percentile(top_vals, 5)), 3),
                round(float(np.percentile(top_vals, 95)), 3),
            )
        return windows

    def pore_size_adsorption_profile(self, lcd_values: np.ndarray,
                                      loadings: np.ndarray,
                                      n_bins: int = 30) -> Dict:
        mask = np.isfinite(lcd_values) & np.isfinite(loadings) & (lcd_values > 0)
        lcd, load = lcd_values[mask], loadings[mask]
        edges = np.linspace(lcd.min(), lcd.max(), n_bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2
        means = np.zeros(n_bins)
        stds = np.zeros(n_bins)
        counts = np.zeros(n_bins, dtype=int)
        for j in range(n_bins):
            in_bin = (lcd >= edges[j]) & (lcd < edges[j + 1])
            if in_bin.sum() > 0:
                means[j] = load[in_bin].mean()
                stds[j] = load[in_bin].std()
                counts[j] = in_bin.sum()
        return {"bin_centers_angstrom": centers.tolist(),
                "mean_loading_mmol_g": means.tolist(),
                "std_loading": stds.tolist(), "counts": counts.tolist()}

    def surface_area_analysis(self, asa_values: np.ndarray,
                               loadings: np.ndarray) -> Dict:
        mask = np.isfinite(asa_values) & np.isfinite(loadings) & (asa_values > 0)
        asa, load = asa_values[mask], loadings[mask]
        if len(asa) < 10:
            return {"slope": 0, "intercept": 0, "r_squared": 0}
        coeffs = np.polyfit(asa, load, 1)
        pred = np.polyval(coeffs, asa)
        ss_res = np.sum((load - pred) ** 2)
        ss_tot = np.sum((load - load.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return {"slope": round(coeffs[0], 6), "intercept": round(coeffs[1], 4),
                "r_squared": round(r2, 4), "n_samples": int(len(asa))}

    def void_fraction_analysis(self, porosity: np.ndarray,
                                loadings_low_p: np.ndarray,
                                loadings_high_p: np.ndarray) -> Dict:
        results = {}
        for label, load in [("low_pressure", loadings_low_p),
                            ("high_pressure", loadings_high_p)]:
            mask = np.isfinite(porosity) & np.isfinite(load) & (porosity > 0)
            p, l = porosity[mask], load[mask]
            if len(p) < 10:
                results[label] = {"r_squared": 0}
                continue
            coeffs = np.polyfit(p, l, 1)
            pred = np.polyval(coeffs, p)
            ss_res = np.sum((l - pred) ** 2)
            ss_tot = np.sum((l - l.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            results[label] = {"slope": round(coeffs[0], 4), "r_squared": round(r2, 4)}
        return results

    def generate_screening_criteria(self, gas: str = "CO2") -> Dict[str, Tuple]:
        if gas == "CO2":
            return {
                "LCD": (self.OPTIMAL_RANGES["CO2_LCD_min"],
                        self.OPTIMAL_RANGES["CO2_LCD_max"]),
                "PLD": (self.OPTIMAL_RANGES["CO2_PLD_min"],
                        self.OPTIMAL_RANGES["CO2_PLD_max"]),
                "ASA": (self.OPTIMAL_RANGES["CO2_ASA_min"],
                        self.OPTIMAL_RANGES["CO2_ASA_max"]),
                "porosity": (self.OPTIMAL_RANGES["CO2_porosity_min"],
                             self.OPTIMAL_RANGES["CO2_porosity_max"]),
            }
        elif gas == "H2":
            return {
                "LCD": (self.OPTIMAL_RANGES["H2_LCD_min"],
                        self.OPTIMAL_RANGES["H2_LCD_max"]),
                "PLD": (self.OPTIMAL_RANGES["H2_PLD_min"], 8.0),
                "ASA": (500, 5000), "porosity": (0.4, 0.9),
            }
        return {}

    def apply_geometric_filter(self, mof_data: List[Dict],
                                criteria: Dict[str, Tuple]) -> List[Dict]:
        passed = []
        for mof in mof_data:
            if all(criteria[p][0] <= mof.get(p, 0) <= criteria[p][1]
                   for p in criteria):
                passed.append(mof)
        logger.info(f"Geometric filter: {len(passed)}/{len(mof_data)} passed")
        return passed
