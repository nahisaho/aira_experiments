"""
Urban Canopy Model (UCM) for Tokyo Heat Island Simulation
Implements a simplified single-layer UCM with building morphology parameterization.
"""
import numpy as np

class UrbanCanopyModel:
    """Single-layer Urban Canopy Model with building morphology parameters."""

    def __init__(self, grid_size=(50, 50), dx=500.0):
        self.nx, self.ny = grid_size
        self.dx = dx  # grid spacing in meters
        self.dt = 60.0  # time step in seconds

        # Building morphology parameters (Tokyo central ward defaults)
        self.building_height = np.zeros((self.nx, self.ny))  # mean building height [m]
        self.building_fraction = np.zeros((self.nx, self.ny))  # plan area fraction
        self.wall_area_ratio = np.zeros((self.nx, self.ny))  # wall-to-plan area ratio
        self.canyon_aspect = np.zeros((self.nx, self.ny))  # H/W ratio
        self.sky_view_factor = np.zeros((self.nx, self.ny))

        # Surface properties
        self.albedo_roof = np.full((self.nx, self.ny), 0.20)
        self.albedo_wall = np.full((self.nx, self.ny), 0.25)
        self.albedo_road = np.full((self.nx, self.ny), 0.08)
        self.emissivity = np.full((self.nx, self.ny), 0.95)
        self.green_fraction = np.zeros((self.nx, self.ny))

        # State variables
        self.T_roof = np.full((self.nx, self.ny), 300.0)  # K
        self.T_wall = np.full((self.nx, self.ny), 300.0)
        self.T_road = np.full((self.nx, self.ny), 300.0)
        self.T_canyon = np.full((self.nx, self.ny), 300.0)

    def initialize_tokyo_morphology(self):
        """Initialize building morphology for Tokyo central wards."""
        np.random.seed(42)

        # Define urban zones: CBD, commercial, residential, suburban, parks
        cx, cy = self.nx // 2, self.ny // 2
        for i in range(self.nx):
            for j in range(self.ny):
                dist = np.sqrt((i - cx)**2 + (j - cy)**2)
                if dist < 8:  # CBD (Marunouchi/Otemachi)
                    self.building_height[i, j] = 80 + np.random.normal(0, 20)
                    self.building_fraction[i, j] = 0.65
                    self.green_fraction[i, j] = 0.05
                elif dist < 15:  # Commercial
                    self.building_height[i, j] = 40 + np.random.normal(0, 10)
                    self.building_fraction[i, j] = 0.50
                    self.green_fraction[i, j] = 0.10
                elif dist < 22:  # Dense residential
                    self.building_height[i, j] = 15 + np.random.normal(0, 5)
                    self.building_fraction[i, j] = 0.40
                    self.green_fraction[i, j] = 0.15
                else:  # Suburban
                    self.building_height[i, j] = 8 + np.random.normal(0, 3)
                    self.building_fraction[i, j] = 0.25
                    self.green_fraction[i, j] = 0.30

        self.building_height = np.clip(self.building_height, 3, 200)
        self.building_fraction = np.clip(self.building_fraction, 0.05, 0.80)

        # Derived parameters
        road_width = self.building_height / np.clip(self.building_fraction * 3, 0.3, 3.0)
        self.canyon_aspect = self.building_height / np.clip(road_width, 5, 100)
        self.wall_area_ratio = 2 * self.canyon_aspect * (1 - self.building_fraction)
        self.sky_view_factor = 1.0 / (1.0 + self.canyon_aspect)

    def compute_radiation_balance(self, S_down, L_down, hour):
        """Compute radiation balance for canyon surfaces."""
        sigma = 5.67e-8
        zenith = np.abs(hour - 12) / 12.0 * np.pi / 2

        # Direct + diffuse shortwave on roof
        S_roof = S_down * (1 - self.albedo_roof) * np.cos(zenith)
        S_roof = np.clip(S_roof, 0, None)

        # Canyon shortwave (multiple reflections)
        tau_canyon = self.sky_view_factor
        S_road = S_down * tau_canyon * (1 - self.albedo_road) * np.cos(zenith)
        S_wall = S_down * (1 - tau_canyon) * 0.5 * (1 - self.albedo_wall)
        S_road = np.clip(S_road, 0, None)
        S_wall = np.clip(S_wall, 0, None)

        # Longwave
        L_roof = self.emissivity * (L_down - sigma * self.T_roof**4)
        L_road = self.emissivity * (tau_canyon * L_down +
                 (1 - tau_canyon) * sigma * self.T_wall**4 - sigma * self.T_road**4)
        L_wall = self.emissivity * (0.5 * L_down +
                 0.5 * sigma * self.T_road**4 - sigma * self.T_wall**4)

        return {
            'Q_roof': S_roof + L_roof,
            'Q_road': S_road + L_road,
            'Q_wall': S_wall + L_wall
        }

    def compute_turbulent_fluxes(self, T_air, wind_speed, humidity):
        """Compute sensible and latent heat fluxes."""
        rho = 1.2  # air density
        cp = 1005  # specific heat
        Ch = 0.005  # bulk transfer coefficient

        # Sensible heat
        H_roof = rho * cp * Ch * wind_speed * (self.T_roof - T_air)
        H_road = rho * cp * Ch * wind_speed * 0.5 * (self.T_road - T_air)
        H_wall = rho * cp * Ch * wind_speed * 0.3 * (self.T_wall - T_air)

        # Latent heat (from green fraction)
        Lv = 2.5e6
        LE = self.green_fraction * Lv * rho * Ch * wind_speed * humidity * 0.001

        return {
            'H_roof': H_roof, 'H_road': H_road, 'H_wall': H_wall,
            'LE': LE,
            'H_total': (self.building_fraction * H_roof +
                       (1 - self.building_fraction) * (H_road + H_wall * self.wall_area_ratio * 0.5))
        }

    def step(self, S_down, L_down, T_air, wind_speed, humidity, Q_anthro, hour):
        """Advance model by one time step."""
        C_roof = 1.0e5  # heat capacity J/m2/K
        C_wall = 1.5e5
        C_road = 2.0e5

        rad = self.compute_radiation_balance(S_down, L_down, hour)
        turb = self.compute_turbulent_fluxes(T_air, wind_speed, humidity)

        # Anthropogenic heat contribution
        Q_a = Q_anthro * (1 - self.building_fraction)

        # Update surface temperatures
        self.T_roof += self.dt / C_roof * (rad['Q_roof'] - turb['H_roof'])
        self.T_wall += self.dt / C_wall * (rad['Q_wall'] - turb['H_wall'])
        self.T_road += self.dt / C_road * (rad['Q_road'] - turb['H_road'] + Q_a)

        # Canyon air temperature
        self.T_canyon = (self.building_fraction * self.T_roof +
                        (1 - self.building_fraction) * 0.5 * (self.T_road + self.T_wall))

        # Evaporative cooling
        cooling = self.green_fraction * 2.0  # K reduction
        self.T_canyon -= cooling

        return {
            'T_canyon': self.T_canyon.copy(),
            'T_roof': self.T_roof.copy(),
            'H_total': turb['H_total'].copy(),
            'UHI_intensity': self.T_canyon - T_air
        }
