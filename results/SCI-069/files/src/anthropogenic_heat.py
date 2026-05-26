"""
Anthropogenic Heat Emission Model
Spatiotemporal distribution of traffic, HVAC, and industrial heat sources.
"""
import numpy as np


class AnthropogenicHeatModel:
    """Models spatiotemporal anthropogenic heat emissions for urban areas."""

    def __init__(self, grid_size=(50, 50)):
        self.nx, self.ny = grid_size
        # Total AH [W/m2] for each source category
        self.Q_traffic = np.zeros((self.nx, self.ny))
        self.Q_hvac = np.zeros((self.nx, self.ny))
        self.Q_industry = np.zeros((self.nx, self.ny))
        self.Q_metabolism = np.zeros((self.nx, self.ny))

    def initialize_tokyo(self, building_fraction, building_height):
        """Set spatial distribution based on urban morphology."""
        cx, cy = self.nx // 2, self.ny // 2

        for i in range(self.nx):
            for j in range(self.ny):
                dist = np.sqrt((i - cx)**2 + (j - cy)**2)

                # Traffic: peaks along major roads and CBD
                if dist < 10:
                    self.Q_traffic[i, j] = 45.0
                elif dist < 20:
                    self.Q_traffic[i, j] = 25.0
                else:
                    self.Q_traffic[i, j] = 10.0

                # HVAC: proportional to building volume
                vol = building_fraction[i, j] * building_height[i, j]
                self.Q_hvac[i, j] = vol * 0.08

                # Industry: concentrated in bay area (south)
                if j < 10 and dist > 15:
                    self.Q_industry[i, j] = 30.0
                else:
                    self.Q_industry[i, j] = 2.0

                # Metabolism: population density proxy
                self.Q_metabolism[i, j] = building_fraction[i, j] * 5.0

    def get_diurnal_profile(self, hour):
        """Return diurnal scaling factors for each source."""
        # Traffic: bimodal peaks at 8 and 18
        traffic_profile = 0.3 + 0.7 * (
            np.exp(-0.5 * ((hour - 8) / 2)**2) +
            np.exp(-0.5 * ((hour - 18) / 2)**2)
        )
        traffic_profile = min(traffic_profile, 1.0)

        # HVAC: peaks during afternoon (cooling demand)
        hvac_profile = 0.2 + 0.8 * np.exp(-0.5 * ((hour - 14) / 4)**2)

        # Industry: daytime operation
        if 6 <= hour <= 22:
            industry_profile = 0.8 + 0.2 * np.exp(-0.5 * ((hour - 12) / 4)**2)
        else:
            industry_profile = 0.2

        # Metabolism: follows activity patterns
        if 6 <= hour <= 23:
            metab_profile = 0.7 + 0.3 * np.exp(-0.5 * ((hour - 15) / 5)**2)
        else:
            metab_profile = 0.4

        return traffic_profile, hvac_profile, industry_profile, metab_profile

    def compute_total(self, hour, climate_factor=1.0):
        """Compute total anthropogenic heat at given hour.

        Args:
            hour: Hour of day (0-23)
            climate_factor: Scaling factor for future climate (e.g., 1.3 for 2050)
        """
        tp, hp, ip, mp = self.get_diurnal_profile(hour)

        Q_total = (self.Q_traffic * tp +
                   self.Q_hvac * hp * climate_factor +
                   self.Q_industry * ip +
                   self.Q_metabolism * mp)

        return Q_total

    def get_component_breakdown(self, hour):
        """Return individual components for analysis."""
        tp, hp, ip, mp = self.get_diurnal_profile(hour)
        return {
            'traffic': self.Q_traffic * tp,
            'hvac': self.Q_hvac * hp,
            'industry': self.Q_industry * ip,
            'metabolism': self.Q_metabolism * mp
        }
