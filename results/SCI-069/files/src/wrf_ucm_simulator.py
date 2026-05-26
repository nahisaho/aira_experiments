"""
WRF-UCM Coupling Simulation Framework
Mesoscale simulation driver for Tokyo UHI prediction.
"""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from ucm_model import UrbanCanopyModel
from anthropogenic_heat import AnthropogenicHeatModel
from mitigation import MitigationScenario
from wbgt_model import WBGTModel


class WRFUCMSimulator:
    """Simplified WRF-UCM coupling framework for UHI simulation."""

    def __init__(self, grid_size=(50, 50), dx=500.0):
        self.grid_size = grid_size
        self.dx = dx
        self.ucm = UrbanCanopyModel(grid_size, dx)
        self.ah_model = AnthropogenicHeatModel(grid_size)
        self.wbgt_model = WBGTModel()

        # Meteorological forcing (simplified)
        self.T_air_base = 303.0  # K (30°C summer baseline)
        self.humidity = 70.0  # %
        self.wind_speed = 3.0  # m/s
        self.solar_max = 900.0  # W/m2

    def initialize(self):
        """Initialize all model components."""
        self.ucm.initialize_tokyo_morphology()
        self.ah_model.initialize_tokyo(
            self.ucm.building_fraction,
            self.ucm.building_height
        )

    def get_forcing(self, hour):
        """Generate meteorological forcing for given hour."""
        # Solar radiation (diurnal cycle)
        if 5 <= hour <= 19:
            solar_angle = np.pi * (hour - 5) / 14
            S_down = self.solar_max * np.sin(solar_angle)
        else:
            S_down = 0.0

        # Longwave from atmosphere
        L_down = 350.0 + 20 * np.sin(2 * np.pi * (hour - 6) / 24)

        # Air temperature (diurnal)
        T_air = self.T_air_base + 5.0 * np.sin(2 * np.pi * (hour - 14) / 24)

        # Wind (weaker at night)
        wind = self.wind_speed * (0.5 + 0.5 * np.sin(2 * np.pi * (hour - 14) / 24))
        wind = max(wind, 0.5)

        return S_down, L_down, T_air, wind

    def run_24h(self, mitigation=None, climate_factor=1.0, year_label="2020"):
        """Run 24-hour simulation."""
        if mitigation is not None:
            mitigation.apply_to_ucm(self.ucm)

        results = {
            'hours': [],
            'T_canyon_mean': [],
            'T_canyon_max': [],
            'UHI_mean': [],
            'UHI_max': [],
            'Q_anthro_mean': [],
            'WBGT_mean': [],
            'WBGT_max': [],
            'spatial_T': [],
            'spatial_UHI': [],
            'spatial_WBGT': [],
            'risk_map': [],
        }

        # Reset temperatures
        self.ucm.T_roof = np.full(self.grid_size, 300.0)
        self.ucm.T_wall = np.full(self.grid_size, 300.0)
        self.ucm.T_road = np.full(self.grid_size, 300.0)

        # Spin-up (3 days)
        for day in range(3):
            for hour_step in range(24 * 60):
                hour = hour_step / 60.0
                S_down, L_down, T_air, wind = self.get_forcing(hour)
                Q_anthro = self.ah_model.compute_total(hour, climate_factor)
                self.ucm.step(S_down, L_down, T_air, wind, self.humidity, Q_anthro, hour)

        # Production run (24h)
        for hour in range(24):
            S_down, L_down, T_air, wind = self.get_forcing(hour)
            Q_anthro = self.ah_model.compute_total(hour, climate_factor)

            # Sub-hourly stepping (60 steps per hour)
            for substep in range(60):
                h = hour + substep / 60.0
                S, L, Ta, w = self.get_forcing(h)
                Qa = self.ah_model.compute_total(h, climate_factor)
                state = self.ucm.step(S, L, Ta, w, self.humidity, Qa, h)

            T_canyon = state['T_canyon']
            uhi = state['UHI_intensity']

            # WBGT computation
            S_rad = max(S_down, 0)
            wbgt_field = np.zeros(self.grid_size)
            for i in range(self.grid_size[0]):
                for j in range(self.grid_size[1]):
                    wbgt_field[i, j] = self.wbgt_model.compute_wbgt_from_meteo(
                        T_canyon[i, j], self.humidity, S_rad, wind
                    )

            risk = self.wbgt_model.classify_risk(wbgt_field)

            results['hours'].append(hour)
            results['T_canyon_mean'].append(np.mean(T_canyon) - 273.15)
            results['T_canyon_max'].append(np.max(T_canyon) - 273.15)
            results['UHI_mean'].append(np.mean(uhi))
            results['UHI_max'].append(np.max(uhi))
            results['Q_anthro_mean'].append(np.mean(Q_anthro))
            results['WBGT_mean'].append(np.mean(wbgt_field))
            results['WBGT_max'].append(np.max(wbgt_field))
            results['spatial_T'].append(T_canyon.copy())
            results['spatial_UHI'].append(uhi.copy())
            results['spatial_WBGT'].append(wbgt_field.copy())
            results['risk_map'].append(risk.copy())

        results['year'] = year_label
        return results


def run_all_scenarios():
    """Run complete experiment suite."""
    grid_size = (50, 50)
    scenarios = {}

    # 1. Baseline 2020
    print("Running baseline 2020...")
    sim = WRFUCMSimulator(grid_size)
    sim.initialize()
    scenarios['baseline_2020'] = sim.run_24h(year_label="2020")

    # 2. Green infrastructure
    print("Running green infrastructure scenario...")
    sim2 = WRFUCMSimulator(grid_size)
    sim2.initialize()
    green = MitigationScenario.create_green_scenario(grid_size)
    scenarios['green'] = sim2.run_24h(mitigation=green, year_label="2020-Green")

    # 3. Cool roofs
    print("Running cool roof scenario...")
    sim3 = WRFUCMSimulator(grid_size)
    sim3.initialize()
    cool = MitigationScenario.create_cool_roof_scenario(grid_size)
    scenarios['cool_roof'] = sim3.run_24h(mitigation=cool, year_label="2020-CoolRoof")

    # 4. Combined
    print("Running combined scenario...")
    sim4 = WRFUCMSimulator(grid_size)
    sim4.initialize()
    combined = MitigationScenario.create_combined_scenario(grid_size)
    scenarios['combined'] = sim4.run_24h(mitigation=combined, year_label="2020-Combined")

    # 5. Future 2050 (no mitigation)
    print("Running 2050 baseline...")
    sim5 = WRFUCMSimulator(grid_size)
    sim5.T_air_base = 305.0  # +2K from climate change
    sim5.initialize()
    scenarios['baseline_2050'] = sim5.run_24h(climate_factor=1.3, year_label="2050")

    # 6. Future 2050 + combined mitigation
    print("Running 2050 + combined mitigation...")
    sim6 = WRFUCMSimulator(grid_size)
    sim6.T_air_base = 305.0
    sim6.initialize()
    combined2 = MitigationScenario.create_combined_scenario(grid_size)
    scenarios['mitigated_2050'] = sim6.run_24h(
        mitigation=combined2, climate_factor=1.3, year_label="2050-Mitigated"
    )

    return scenarios


if __name__ == "__main__":
    scenarios = run_all_scenarios()

    # Print summary
    for name, result in scenarios.items():
        peak_hour = np.argmax(result['UHI_max'])
        print(f"\n=== {name} ({result['year']}) ===")
        print(f"  Peak UHI: {result['UHI_max'][peak_hour]:.2f} K at {peak_hour}:00")
        print(f"  Mean daily UHI: {np.mean(result['UHI_mean']):.2f} K")
        print(f"  Peak WBGT: {max(result['WBGT_max']):.1f} °C")
        print(f"  Peak canyon T: {max(result['T_canyon_max']):.1f} °C")
