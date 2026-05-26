"""
WBGT (Wet Bulb Globe Temperature) Heat Stress Risk Model
Based on Liljegren et al. (2008) simplified approach.
"""
import numpy as np


class WBGTModel:
    """Compute WBGT and heat stress risk categories."""

    # Risk thresholds (°C)
    THRESHOLDS = {
        'low': 25.0,
        'moderate': 28.0,
        'high': 31.0,
        'very_high': 33.0,
        'extreme': 35.0
    }

    RISK_LABELS = ['Low', 'Moderate', 'High', 'Very High', 'Extreme']

    @staticmethod
    def compute_wbgt_outdoor(T_air, T_globe, T_wet, wind_speed=1.0):
        """Compute outdoor WBGT.
        WBGT_outdoor = 0.7 * T_wet + 0.2 * T_globe + 0.1 * T_air
        """
        return 0.7 * T_wet + 0.2 * T_globe + 0.1 * T_air

    @staticmethod
    def estimate_globe_temperature(T_air, solar_rad, wind_speed):
        """Estimate globe temperature from meteorological variables.
        Simplified Liljegren approximation.
        """
        T_globe = T_air + 0.01 * solar_rad - 0.5 * np.sqrt(wind_speed)
        return T_globe

    @staticmethod
    def estimate_wet_bulb(T_air, humidity):
        """Estimate natural wet bulb temperature.
        Stull (2011) approximation.
        """
        T_c = T_air - 273.15 if T_air > 200 else T_air  # handle K or C
        RH = humidity

        T_wet = T_c * np.arctan(0.151977 * np.sqrt(RH + 8.313659)) + \
                np.arctan(T_c + RH) - np.arctan(RH - 1.676331) + \
                0.00391838 * RH**1.5 * np.arctan(0.023101 * RH) - 4.686035

        if T_air > 200:  # was in Kelvin
            T_wet += 273.15

        return T_wet

    def compute_wbgt_from_meteo(self, T_air_K, humidity, solar_rad, wind_speed):
        """Compute WBGT from standard meteorological variables."""
        T_air_C = T_air_K - 273.15
        T_globe = self.estimate_globe_temperature(T_air_C, solar_rad, wind_speed)
        T_wet = self.estimate_wet_bulb(T_air_C, humidity)
        wbgt = self.compute_wbgt_outdoor(T_air_C, T_globe, T_wet, wind_speed)
        return wbgt

    def classify_risk(self, wbgt):
        """Classify heat stress risk level."""
        if isinstance(wbgt, np.ndarray):
            risk = np.zeros_like(wbgt, dtype=int)
            risk[wbgt >= self.THRESHOLDS['moderate']] = 1
            risk[wbgt >= self.THRESHOLDS['high']] = 2
            risk[wbgt >= self.THRESHOLDS['very_high']] = 3
            risk[wbgt >= self.THRESHOLDS['extreme']] = 4
            return risk
        else:
            if wbgt >= self.THRESHOLDS['extreme']:
                return 4
            elif wbgt >= self.THRESHOLDS['very_high']:
                return 3
            elif wbgt >= self.THRESHOLDS['high']:
                return 2
            elif wbgt >= self.THRESHOLDS['moderate']:
                return 1
            return 0

    def compute_heatstroke_risk_index(self, wbgt, population_density, age_fraction_elderly):
        """Compute composite heatstroke risk index.
        Combines WBGT with demographic vulnerability.
        """
        thermal_risk = np.clip((wbgt - 25) / 10, 0, 1)
        vulnerability = 1.0 + 2.0 * age_fraction_elderly
        exposure = np.log1p(population_density) / 10.0

        risk_index = thermal_risk * vulnerability * exposure
        return np.clip(risk_index, 0, 1)
