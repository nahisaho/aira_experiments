"""
Daylighting Simulation Engine
Simplified Radiance/Honeybee-compatible daylighting analysis.
"""
import numpy as np
import math
from typing import Dict, List, Tuple


class DaylightSimulation:
    """Simplified daylighting simulation using radiosity-based method."""

    def __init__(self, params: Dict, zone_params: Dict):
        self.params = params
        self.zone = zone_params
        self.grid_spacing = params['analysis_grid']['spacing']
        self.work_plane_height = params['analysis_grid']['height']

        # Room dimensions (from zone)
        self.length = 25.0  # m
        self.width = 20.0   # m
        self.height = 3.5   # m

        # Grid
        self.nx = int(self.length / self.grid_spacing)
        self.ny = int(self.width / self.grid_spacing)
        self.illuminance_grid = np.zeros((self.ny, self.nx))

    def _sky_luminance(self, month: int, hour: int) -> float:
        """Calculate exterior horizontal illuminance (simplified Perez model)."""
        day_of_year = (month - 1) * 30 + 15
        solar_declination = 23.45 * math.sin(math.radians(360 / 365 * (day_of_year - 81)))
        latitude = self.params['location']['latitude']

        hour_angle = 15 * (hour - 12)
        sin_alt = (math.sin(math.radians(latitude)) * math.sin(math.radians(solar_declination)) +
                   math.cos(math.radians(latitude)) * math.cos(math.radians(solar_declination)) *
                   math.cos(math.radians(hour_angle)))
        solar_altitude = math.degrees(math.asin(max(-1, min(1, sin_alt))))

        if solar_altitude <= 0:
            return 0

        # Exterior illuminance (lux) - simplified
        luminous_efficacy = 110  # lm/W typical
        solar_constant = 1367  # W/m²
        clearness = 0.65
        ext_illuminance = solar_constant * luminous_efficacy * math.sin(math.radians(solar_altitude)) * clearness
        return min(ext_illuminance, 100000)

    def _window_daylight_factor(self, x: float, y: float, win_x: float, win_y: float,
                                  win_width: float, win_height: float) -> float:
        """Calculate daylight factor at a point from a window."""
        dx = abs(x - win_x)
        dy = abs(y - win_y)
        dist = math.sqrt(dx**2 + dy**2 + self.work_plane_height**2)

        if dist < 0.5:
            dist = 0.5

        vlt = self.params['materials']['glass_transmittance']
        win_area = win_width * win_height

        # Simplified daylight factor (BRE split-flux method approximation)
        sc = win_area / (2 * math.pi * dist**2)
        wall_ref = self.params['materials']['wall_reflectance']
        ceil_ref = self.params['materials']['ceiling_reflectance']
        floor_ref = self.params['materials']['floor_reflectance']
        avg_ref = (wall_ref + ceil_ref + floor_ref) / 3

        # Direct component
        df_direct = sc * vlt * 100

        # Inter-reflected component
        room_area = 2 * (self.length * self.width + self.length * self.height + self.width * self.height)
        total_win_area = self.zone.get('total_window_area', 30)
        df_irc = (total_win_area * vlt * avg_ref) / (room_area * (1 - avg_ref)) * 100

        return df_direct + df_irc * 0.5

    def calculate_daylight_grid(self, month: int = 9, hour: int = 12) -> np.ndarray:
        """Calculate illuminance distribution across the floor plan."""
        ext_illum = self._sky_luminance(month, hour)
        if ext_illum <= 0:
            return self.illuminance_grid

        # Define window positions
        windows = []
        # South windows (y=0)
        for i in range(5):
            windows.append({
                'x': 3 + i * 4.5, 'y': 0, 'width': 1.5, 'height': 1.8,
                'orientation': 'south'
            })
        # North windows (y=width)
        for i in range(4):
            windows.append({
                'x': 4 + i * 5, 'y': self.width, 'width': 1.5, 'height': 1.8,
                'orientation': 'north'
            })
        # East windows
        for i in range(3):
            windows.append({
                'x': self.length, 'y': 3 + i * 5, 'width': 1.5, 'height': 1.8,
                'orientation': 'east'
            })
        # West windows
        for i in range(3):
            windows.append({
                'x': 0, 'y': 3 + i * 5, 'width': 1.5, 'height': 1.8,
                'orientation': 'west'
            })

        orientation_factors = {
            'south': 1.0, 'north': 0.5, 'east': 0.7, 'west': 0.7
        }

        for j in range(self.ny):
            for i in range(self.nx):
                x = (i + 0.5) * self.grid_spacing
                y = (j + 0.5) * self.grid_spacing

                total_df = 0
                for win in windows:
                    df = self._window_daylight_factor(x, y, win['x'], win['y'],
                                                     win['width'], win['height'])
                    of = orientation_factors[win['orientation']]
                    total_df += df * of

                self.illuminance_grid[j, i] = ext_illum * total_df / 100

        return self.illuminance_grid

    def calculate_annual_metrics(self) -> Dict:
        """Calculate annual daylighting metrics: sDA, ASE, UDI."""
        occupied_hours = 0
        sda_count = np.zeros((self.ny, self.nx))
        ase_count = np.zeros((self.ny, self.nx))
        udi_count = np.zeros((self.ny, self.nx))
        da_count = np.zeros((self.ny, self.nx))
        total_annual_illum = np.zeros((self.ny, self.nx))

        for month in range(1, 13):
            for hour in range(self.params['schedules']['occupancy_start'],
                              self.params['schedules']['occupancy_end'] + 1):
                grid = self.calculate_daylight_grid(month, hour)
                occupied_hours += 1

                # sDA300/50: illuminance >= 300 lux
                sda_count += (grid >= 300).astype(float)
                # ASE1000/250: illuminance >= 1000 lux
                ase_count += (grid >= 1000).astype(float)
                # UDI100-3000: 100 <= illuminance <= 3000
                udi_count += ((grid >= 100) & (grid <= 3000)).astype(float)
                # DA300: illuminance >= 300
                da_count += (grid >= 300).astype(float)
                total_annual_illum += grid

        # Calculate spatial metrics
        sda_pct = sda_count / occupied_hours * 100
        ase_pct = ase_count / occupied_hours * 100
        udi_pct = udi_count / occupied_hours * 100
        da_pct = da_count / occupied_hours * 100

        # Spatial daylight autonomy (sDA300/50): % of floor where DA >= 50%
        sda_value = np.mean(sda_pct >= 50) * 100
        # Annual sunlight exposure (ASE1000/250): % of floor where ASE >= 250 hours
        ase_threshold_hours = 250 / occupied_hours * 100
        ase_value = np.mean(ase_pct >= ase_threshold_hours) * 100

        mean_da = np.mean(da_pct)
        mean_udi = np.mean(udi_pct)
        mean_annual_illum = np.mean(total_annual_illum) / occupied_hours

        # LEED v4.1 daylight credit assessment
        leed_daylight = "Option 2"
        if sda_value >= 75 and ase_value <= 10:
            leed_points = 3
        elif sda_value >= 55 and ase_value <= 10:
            leed_points = 2
        elif sda_value >= 40:
            leed_points = 1
        else:
            leed_points = 0

        return {
            'sDA300_50': round(float(sda_value), 1),
            'ASE1000_250': round(float(ase_value), 1),
            'mean_DA300': round(float(mean_da), 1),
            'mean_UDI_100_3000': round(float(mean_udi), 1),
            'mean_annual_illuminance_lux': round(float(mean_annual_illum), 0),
            'occupied_hours': occupied_hours,
            'sda_grid': sda_pct.tolist(),
            'ase_grid': ase_pct.tolist(),
            'udi_grid': udi_pct.tolist(),
            'leed_daylight_points': leed_points,
            'leed_sda_pass': sda_value >= 55,
            'leed_ase_pass': ase_value <= 10,
            'illuminance_grid_sample': self.illuminance_grid.tolist(),
        }


if __name__ == "__main__":
    from ifc_converter import IFCConverter
    import json

    converter = IFCConverter()
    model = converter.create_reference_building()
    rad_params = converter.generate_radiance_params()
    ep_params = converter.generate_energyplus_params()

    sim = DaylightSimulation(rad_params, ep_params['zones'][0])
    metrics = sim.calculate_annual_metrics()
    print("=== Daylighting Metrics ===")
    for k, v in metrics.items():
        if 'grid' not in k:
            print(f"  {k}: {v}")
