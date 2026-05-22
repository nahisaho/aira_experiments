"""
Heat Stress Risk Assessment — WBGT Prediction Module

Computes WBGT and associated heat stroke risk from UCM/WRF output.
WBGT_outdoor = 0.7 × T_w + 0.2 × T_g + 0.1 × T_a
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class PopulationExposure:
    total_population: int = 1_000_000
    outdoor_fraction: float = 0.15
    elderly_fraction: float = 0.28
    elderly_relative_risk: float = 3.0
    child_fraction: float = 0.12
    child_relative_risk: float = 1.5
    outdoor_worker_fraction: float = 0.05
    worker_relative_risk: float = 2.5


class WBGTCalculator:
    @staticmethod
    def wet_bulb_temperature(T_air, rh):
        """Stull (2011) approximation."""
        return (T_air * np.arctan(0.151977 * np.sqrt(rh + 8.313659))
                + np.arctan(T_air + rh) - np.arctan(rh - 1.676331)
                + 0.00391838 * rh**(3/2) * np.arctan(0.023101 * rh) - 4.686035)

    @staticmethod
    def globe_temperature(T_air, T_mrt, wind_speed, solar_rad=800.0):
        T_g = 1.01 * T_air + 0.006 * solar_rad - 0.35 * np.sqrt(max(wind_speed, 0.1)) + 2.8
        if T_mrt > 0:
            T_g = 0.5 * (T_g + 0.7 * T_mrt + 0.3 * T_air)
        return T_g

    def compute_wbgt_outdoor(self, T_air, rh, wind_speed, solar_rad, T_mrt=0.0):
        T_w = self.wet_bulb_temperature(T_air, rh)
        T_g = self.globe_temperature(T_air, T_mrt, wind_speed, solar_rad)
        return 0.7 * T_w + 0.2 * T_g + 0.1 * T_air

    @staticmethod
    def classify_risk(wbgt):
        if wbgt >= 31:
            return {"level": "危険", "level_en": "Danger", "color": "red",
                    "action": "運動は原則中止。外出を避ける"}
        elif wbgt >= 28:
            return {"level": "厳重警戒", "level_en": "Severe Warning", "color": "orange",
                    "action": "激しい運動は中止。こまめな水分補給"}
        elif wbgt >= 25:
            return {"level": "警戒", "level_en": "Warning", "color": "yellow",
                    "action": "積極的に休憩。水分・塩分補給"}
        elif wbgt >= 21:
            return {"level": "注意", "level_en": "Caution", "color": "lightblue",
                    "action": "適度な休憩と水分補給"}
        else:
            return {"level": "ほぼ安全", "level_en": "Almost Safe", "color": "blue",
                    "action": "通常の注意"}


class HeatStrokeRiskAssessor:
    def __init__(self, population=None):
        self.pop = population or PopulationExposure()
        self.wbgt_calc = WBGTCalculator()

    def estimate_patient_count(self, wbgt, duration_hours=1.0):
        if wbgt < 25:
            base_rate = 0.1
        else:
            base_rate = np.exp(0.12 * wbgt - 2.5)

        outdoor_pop = self.pop.total_population * self.pop.outdoor_fraction
        general_pop = outdoor_pop * (1 - self.pop.elderly_fraction
                                      - self.pop.child_fraction
                                      - self.pop.outdoor_worker_fraction)
        elderly_eq = outdoor_pop * self.pop.elderly_fraction * self.pop.elderly_relative_risk
        child_eq = outdoor_pop * self.pop.child_fraction * self.pop.child_relative_risk
        worker_eq = outdoor_pop * self.pop.outdoor_worker_fraction * self.pop.worker_relative_risk
        effective_100k = (general_pop + elderly_eq + child_eq + worker_eq) / 100_000
        patients = base_rate * effective_100k * duration_hours
        return {
            "wbgt": wbgt,
            "risk_level": self.wbgt_calc.classify_risk(wbgt),
            "estimated_patients": round(patients, 1),
            "patients_per_100k_hr": round(base_rate, 2),
        }

    def daily_risk_profile(self, hourly_T, hourly_rh, hourly_wind, hourly_solar,
                            hourly_Tmrt=None):
        if hourly_Tmrt is None:
            hourly_Tmrt = np.zeros(24)
        wbgt = np.zeros(24)
        risk_levels = []
        patients = np.zeros(24)
        for h in range(24):
            wbgt[h] = self.wbgt_calc.compute_wbgt_outdoor(
                hourly_T[h], hourly_rh[h], hourly_wind[h], hourly_solar[h], hourly_Tmrt[h])
            risk = self.estimate_patient_count(wbgt[h])
            risk_levels.append(risk["risk_level"]["level"])
            patients[h] = risk["estimated_patients"]
        return {
            "hourly_wbgt": wbgt, "hourly_risk_level": risk_levels,
            "hourly_patients": patients, "daily_total_patients": patients.sum(),
            "peak_wbgt": wbgt.max(), "peak_hour": int(np.argmax(wbgt)),
            "danger_hours": int(np.sum(wbgt >= 31)),
            "severe_warning_hours": int(np.sum((wbgt >= 28) & (wbgt < 31))),
        }


TOKYO_WARD_POPULATIONS = {
    "chiyoda": PopulationExposure(total_population=67_000, outdoor_fraction=0.25,
                                   elderly_fraction=0.22, outdoor_worker_fraction=0.08),
    "chuo": PopulationExposure(total_population=175_000, outdoor_fraction=0.20,
                                elderly_fraction=0.20, outdoor_worker_fraction=0.06),
    "minato": PopulationExposure(total_population=260_000, outdoor_fraction=0.20,
                                  elderly_fraction=0.18, outdoor_worker_fraction=0.05),
    "shinjuku": PopulationExposure(total_population=350_000, outdoor_fraction=0.20,
                                    elderly_fraction=0.22, outdoor_worker_fraction=0.04),
    "nerima": PopulationExposure(total_population=740_000, outdoor_fraction=0.12,
                                  elderly_fraction=0.25, outdoor_worker_fraction=0.03),
    "setagaya": PopulationExposure(total_population=920_000, outdoor_fraction=0.12,
                                    elderly_fraction=0.22, outdoor_worker_fraction=0.03),
}
