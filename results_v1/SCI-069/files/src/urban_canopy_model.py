"""
Urban Canopy Model (UCM) — Building Morphology Parameterization

Implements a single-layer UCM following Kusaka et al. (2001) with extensions
for Tokyo's high-density built environment. Parameterizes:
  - Street canyon geometry (H/W ratio, sky view factor)
  - Building surface energy balance
  - In-canyon radiation trapping
  - Anthropogenic heat injection interface
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple


@dataclass
class BuildingMorphology:
    """Building morphology parameters for a grid cell."""
    building_height_mean: float = 20.0        # [m]
    building_height_std: float = 10.0         # [m]
    building_width_mean: float = 15.0         # [m]
    street_width_mean: float = 12.0           # [m]
    building_fraction: float = 0.45           # [-] plan area fraction
    wall_area_ratio: float = 2.5              # [-] wall/plan area
    roof_fraction: float = 0.45              # [-]
    impervious_fraction: float = 0.40        # [-]
    green_fraction: float = 0.15             # [-]

    @property
    def hw_ratio(self) -> float:
        """Canyon height-to-width ratio."""
        return self.building_height_mean / self.street_width_mean

    @property
    def sky_view_factor(self) -> float:
        """Sky view factor from canyon floor (Oke, 1981)."""
        hw = self.hw_ratio
        return np.sqrt(1 + hw**2) - hw

    @property
    def canyon_aspect_ratio(self) -> float:
        return self.building_height_mean / (self.street_width_mean + self.building_width_mean)


# Tokyo district morphology presets
TOKYO_MORPHOLOGY = {
    "marunouchi": BuildingMorphology(
        building_height_mean=120.0, building_height_std=50.0,
        building_width_mean=50.0, street_width_mean=25.0,
        building_fraction=0.55, wall_area_ratio=4.0,
        roof_fraction=0.55, impervious_fraction=0.35, green_fraction=0.10
    ),
    "shinjuku": BuildingMorphology(
        building_height_mean=80.0, building_height_std=40.0,
        building_width_mean=30.0, street_width_mean=20.0,
        building_fraction=0.50, wall_area_ratio=3.5,
        roof_fraction=0.50, impervious_fraction=0.35, green_fraction=0.15
    ),
    "residential_23ku": BuildingMorphology(
        building_height_mean=8.0, building_height_std=3.0,
        building_width_mean=8.0, street_width_mean=6.0,
        building_fraction=0.40, wall_area_ratio=2.0,
        roof_fraction=0.40, impervious_fraction=0.40, green_fraction=0.20
    ),
    "suburban": BuildingMorphology(
        building_height_mean=6.0, building_height_std=2.0,
        building_width_mean=10.0, street_width_mean=8.0,
        building_fraction=0.25, wall_area_ratio=1.2,
        roof_fraction=0.25, impervious_fraction=0.30, green_fraction=0.45
    ),
}


@dataclass
class SurfaceThermalProperties:
    """Thermal properties for urban surfaces."""
    albedo: float = 0.15                   # [-]
    emissivity: float = 0.95               # [-]
    thermal_conductivity: float = 1.0      # [W/m/K]
    heat_capacity: float = 1.5e6           # [J/m³/K]
    roughness_length: float = 1.0          # [m]
    thickness: float = 0.5                 # [m] representative layer

    @property
    def thermal_diffusivity(self) -> float:
        return self.thermal_conductivity / self.heat_capacity


SURFACE_MATERIALS = {
    "concrete":     SurfaceThermalProperties(0.20, 0.90, 1.4, 2.0e6, 0.5),
    "asphalt":      SurfaceThermalProperties(0.10, 0.95, 0.75, 1.9e6, 0.3),
    "glass_facade": SurfaceThermalProperties(0.25, 0.85, 1.0, 1.8e6, 0.05),
    "green_roof":   SurfaceThermalProperties(0.22, 0.95, 0.3, 1.2e6, 0.2),
    "cool_roof":    SurfaceThermalProperties(0.60, 0.90, 1.4, 2.0e6, 0.5),
    "cool_pavement": SurfaceThermalProperties(0.40, 0.90, 0.75, 1.9e6, 0.3),
    "soil":         SurfaceThermalProperties(0.20, 0.95, 0.5, 1.4e6, 0.5),
    "water":        SurfaceThermalProperties(0.06, 0.97, 0.6, 4.2e6, 1.0),
}


class UrbanCanopyModel:
    """
    Single-layer Urban Canopy Model.

    Solves the canyon energy balance:
      Q* + QF = QH + QE + ΔQS
    where:
      Q*  = net all-wave radiation
      QF  = anthropogenic heat flux
      QH  = sensible heat flux
      QE  = latent heat flux
      ΔQS = storage heat flux
    """

    STEFAN_BOLTZMANN = 5.67e-8
    CP_AIR = 1005.0
    RHO_AIR = 1.2
    LATENT_HEAT_VAP = 2.45e6

    def __init__(self, morphology: BuildingMorphology,
                 roof_props: Optional[SurfaceThermalProperties] = None,
                 wall_props: Optional[SurfaceThermalProperties] = None,
                 road_props: Optional[SurfaceThermalProperties] = None):
        self.morph = morphology
        self.roof = roof_props or SURFACE_MATERIALS["concrete"]
        self.wall = wall_props or SURFACE_MATERIALS["concrete"]
        self.road = road_props or SURFACE_MATERIALS["asphalt"]

        self.T_roof = 300.0
        self.T_wall = 300.0
        self.T_road = 300.0
        self.T_canyon = 300.0

    def compute_radiation_trapping(self, sw_down: float, lw_down: float) -> dict:
        """Multi-reflection radiation in street canyon."""
        svf = self.morph.sky_view_factor
        sw_roof = sw_down * (1 - self.roof.albedo)
        sw_road_direct = sw_down * svf * (1 - self.road.albedo)
        view_wall = (1 - svf) / 2
        alpha_w = self.wall.albedo
        alpha_r = self.road.albedo
        sw_wall = sw_down * view_wall * (1 - alpha_w) * (1 + alpha_r * (1 - svf))
        eps_w = self.wall.emissivity
        eps_r = self.road.emissivity
        lw_roof = lw_down * self.roof.emissivity
        lw_road = (lw_down * svf * eps_r
                   + eps_r * eps_w * (1 - svf) * self.STEFAN_BOLTZMANN * self.T_wall**4)
        lw_wall = (lw_down * view_wall * eps_w
                   + eps_w * eps_r * 0.5 * self.STEFAN_BOLTZMANN * self.T_road**4
                   + eps_w**2 * (1 - svf - view_wall) * self.STEFAN_BOLTZMANN * self.T_wall**4)

        return {
            "sw_roof": sw_roof, "sw_wall": sw_wall, "sw_road": sw_road_direct,
            "lw_roof": lw_roof, "lw_wall": lw_wall, "lw_road": lw_road,
            "net_roof": sw_roof + lw_roof - self.roof.emissivity * self.STEFAN_BOLTZMANN * self.T_roof**4,
            "net_wall": sw_wall + lw_wall - eps_w * self.STEFAN_BOLTZMANN * self.T_wall**4,
            "net_road": sw_road_direct + lw_road - eps_r * self.STEFAN_BOLTZMANN * self.T_road**4,
        }

    def compute_turbulent_fluxes(self, T_surface: float, T_air: float,
                                 u_star: float = 0.3, q_surface: float = 0.005,
                                 q_air: float = 0.008) -> Tuple[float, float]:
        """Sensible and latent heat fluxes via bulk aerodynamic method."""
        C_H = 0.003 * (1 + 0.5 * self.morph.hw_ratio)
        QH = self.RHO_AIR * self.CP_AIR * C_H * u_star * (T_surface - T_air)
        QE = self.RHO_AIR * self.LATENT_HEAT_VAP * C_H * u_star * (q_surface - q_air)
        return QH, max(QE, 0)

    def compute_storage_flux(self, Q_net: float, dT_dt: float = 0.0) -> float:
        """Objective Hysteresis Model (OHM) for storage heat flux."""
        bf = self.morph.building_fraction
        a1 = 0.20 + 0.30 * bf
        a2 = 0.10
        a3 = -15.0 * bf
        return a1 * Q_net + a2 * dT_dt + a3

    def energy_balance(self, sw_down: float, lw_down: float, T_air: float,
                       u_star: float, q_air: float,
                       Q_anthropogenic: float = 0.0,
                       dt: float = 3600.0) -> dict:
        """
        Solve the complete urban energy balance for one timestep.

        Q* + QF = QH + QE + ΔQS
        """
        rad = self.compute_radiation_trapping(sw_down, lw_down)
        bf = self.morph.building_fraction
        rf = self.morph.roof_fraction
        war = self.morph.wall_area_ratio
        Q_star = (rf * rad["net_roof"]
                  + (1 - rf) * rad["net_road"]
                  + war * (1 - rf) * rad["net_wall"])

        # Surface temperature: relax toward equilibrium with implicit damping
        C_eff = 2.0e6 * 0.5  # effective thermal mass [J/m²/K]
        # Damping coefficient: linearized LW + turbulent coupling
        lambda_damp = (4 * self.STEFAN_BOLTZMANN * 0.92 * T_air**3
                       + self.RHO_AIR * self.CP_AIR * 0.003 * u_star)

        for surf, net_key, qf_frac in [
            ("T_roof", "net_roof", 0.2),
            ("T_wall", "net_wall", 0.0),
            ("T_road", "net_road", 0.3),
        ]:
            T_s = getattr(self, surf)
            Q_net = rad[net_key] + Q_anthropogenic * qf_frac
            # Implicit update: C dT/dt = Q_net - λ(T_s - T_air)
            T_eq = T_air + Q_net / lambda_damp
            tau = C_eff / lambda_damp
            alpha = dt / (tau + dt)  # implicit blending
            T_new = T_s + alpha * (T_eq - T_s)
            setattr(self, surf, T_new)

        T_urban = rf * self.T_roof + (1 - rf) * 0.5 * (self.T_road + self.T_wall)
        QH, QE = self.compute_turbulent_fluxes(T_urban, T_air, u_star, q_air=q_air)
        gf = self.morph.green_fraction
        QH_adj = QH * (1 - 0.6 * gf)
        QE_adj = QE + QH * 0.4 * gf
        dQS = self.compute_storage_flux(Q_star)
        residual = Q_star + Q_anthropogenic - QH_adj - QE_adj - dQS

        # Canyon air: implicit relaxation toward T_air with heat addition
        rho_cp = self.RHO_AIR * self.CP_AIR
        effective_height = max(self.morph.building_height_mean, 10.0)
        C_canyon = rho_cp * effective_height  # [J/m²/K]
        ventilation_lambda = rho_cp * u_star * 0.05  # [W/m²/K]
        Q_canyon = QH_adj * 0.3 + Q_anthropogenic * 0.3
        T_eq_canyon = T_air + Q_canyon / (ventilation_lambda + 1e-6)
        tau_canyon = C_canyon / (ventilation_lambda + 1e-6)
        alpha_c = dt / (tau_canyon + dt)
        self.T_canyon = self.T_canyon + alpha_c * (T_eq_canyon - self.T_canyon)

        return {
            "Q_star": Q_star, "QF": Q_anthropogenic,
            "QH": QH_adj, "QE": QE_adj, "dQS": dQS,
            "residual": residual,
            "T_roof": self.T_roof, "T_wall": self.T_wall,
            "T_road": self.T_road, "T_canyon": self.T_canyon,
            "UHI_intensity": self.T_canyon - T_air,
        }

    def reset(self, T_init: float = 300.0):
        self.T_roof = T_init
        self.T_wall = T_init
        self.T_road = T_init
        self.T_canyon = T_init


def compute_local_climate_zone(morphology: BuildingMorphology) -> str:
    """Classify grid cell into Local Climate Zone (Stewart & Oke, 2012)."""
    h = morphology.building_height_mean
    bf = morphology.building_fraction
    if h > 50 and bf > 0.4:
        return "LCZ-1 (Compact High-rise)"
    elif h > 25 and bf > 0.4:
        return "LCZ-2 (Compact Mid-rise)"
    elif h > 25 and bf <= 0.4:
        return "LCZ-4 (Open High-rise)"
    elif 10 < h <= 25 and bf > 0.4:
        return "LCZ-2 (Compact Mid-rise)"
    elif 10 < h <= 25:
        return "LCZ-5 (Open Mid-rise)"
    elif 3 < h <= 10 and bf > 0.4:
        return "LCZ-3 (Compact Low-rise)"
    elif 3 < h <= 10:
        return "LCZ-6 (Open Low-rise)"
    else:
        return "LCZ-9 (Sparsely Built)"
