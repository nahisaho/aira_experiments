"""
WRF-UCM Coupling Framework for Mesoscale UHI Simulation

Implements:
  - WRF namelist generation for urban physics options
  - UCM parameter tables for WRF-Urban
  - Domain nesting configuration for Tokyo region
  - Offline coupling mode for scenario analysis
  - ENVI-met microscale interface
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple


@dataclass
class WRFDomain:
    name: str
    dx: float
    dy: float
    nx: int
    ny: int
    nz: int = 40
    dt: float = None
    parent_ratio: int = 1
    center_lat: float = 35.6762
    center_lon: float = 139.6503

    def __post_init__(self):
        if self.dt is None:
            self.dt = min(6 * self.dx / 1000, 180)


TOKYO_DOMAINS = [
    WRFDomain("d01_kanto", dx=9000, dy=9000, nx=100, ny=100, nz=45),
    WRFDomain("d02_tokyo", dx=3000, dy=3000, nx=121, ny=121, nz=45, parent_ratio=3),
    WRFDomain("d03_23ku", dx=1000, dy=1000, nx=151, ny=151, nz=45, parent_ratio=3),
    WRFDomain("d04_focus", dx=333, dy=333, nx=151, ny=151, nz=45, parent_ratio=3),
]


@dataclass
class WRFUrbanConfig:
    sf_urban_physics: int = 1
    num_urban_categories: int = 3
    use_wudapt_lcz: bool = True
    num_urban_hi: int = 15


class WRFNamelistGenerator:
    def generate_physics_namelist(self, config, domains):
        nd = len(domains)
        lines = [
            "&physics",
            f" sf_urban_physics = {', '.join([str(config.sf_urban_physics)] * nd)},",
            f" num_urban_categories = {config.num_urban_categories},",
            f" num_land_cat = 21,",
            f" sf_surface_physics = {', '.join(['2'] * nd)},",
            f" bl_pbl_physics = {', '.join(['1'] * nd)},",
            f" ra_lw_physics = {', '.join(['4'] * nd)},",
            f" ra_sw_physics = {', '.join(['4'] * nd)},",
            f" sf_sfclay_physics = {', '.join(['1'] * nd)},",
            f" cu_physics = {', '.join(['1'] + ['0'] * (nd - 1))},",
            f" mp_physics = {', '.join(['8'] * nd)},",
            "/",
        ]
        return "\n".join(lines)

    def generate_domains_namelist(self, domains):
        nd = len(domains)
        lines = ["&domains", f" max_dom = {nd},"]
        for attr in ["dx", "dy"]:
            vals = ", ".join([f"{getattr(d, attr):.0f}" for d in domains])
            lines.append(f" {attr} = {vals},")
        for attr, mapped in [("nx", "e_we"), ("ny", "e_sn"), ("nz", "e_vert")]:
            vals = ", ".join([str(getattr(d, attr)) for d in domains])
            lines.append(f" {mapped} = {vals},")
        ratios = ", ".join([str(d.parent_ratio) for d in domains])
        lines.append(f" parent_grid_ratio = {ratios},")
        lines.append(f" time_step = {int(domains[0].dt)},")
        lines.append("/")
        return "\n".join(lines)


LCZ_URBAN_PARAMS = {
    1: {"ZR": 50, "ZW": 30, "SW": 15, "BF": 0.55, "WAR": 4.0,
        "ALB_R": 0.15, "ALB_W": 0.20, "ALB_G": 0.10,
        "EM_R": 0.90, "EM_W": 0.90, "EM_G": 0.95,
        "LAM_R": 1.4, "LAM_W": 1.4, "LAM_G": 0.75,
        "C_R": 2.0e6, "C_W": 2.0e6, "C_G": 1.9e6},
    2: {"ZR": 18, "ZW": 12, "SW": 10, "BF": 0.50, "WAR": 2.5,
        "ALB_R": 0.15, "ALB_W": 0.20, "ALB_G": 0.10,
        "EM_R": 0.90, "EM_W": 0.90, "EM_G": 0.95,
        "LAM_R": 1.4, "LAM_W": 1.4, "LAM_G": 0.75,
        "C_R": 2.0e6, "C_W": 2.0e6, "C_G": 1.9e6},
    3: {"ZR": 7, "ZW": 6, "SW": 5, "BF": 0.55, "WAR": 1.5,
        "ALB_R": 0.15, "ALB_W": 0.25, "ALB_G": 0.12,
        "EM_R": 0.90, "EM_W": 0.90, "EM_G": 0.95,
        "LAM_R": 1.0, "LAM_W": 1.0, "LAM_G": 0.75,
        "C_R": 1.8e6, "C_W": 1.8e6, "C_G": 1.9e6},
}


class OfflineCouplingEngine:
    """Offline WRF-UCM coupling for scenario analysis."""

    def __init__(self, ucm, anthropogenic_model, cooling_model=None):
        self.ucm = ucm
        self.anthro = anthropogenic_model
        self.cooling = cooling_model

    def run_diurnal_cycle(self, forcing, dt=3600.0):
        nt = len(forcing['sw_down'])
        results = {key: np.zeros(nt) for key in
                   ['Q_star', 'QF', 'QH', 'QE', 'dQS', 'T_canyon',
                    'T_roof', 'T_wall', 'T_road', 'UHI_intensity',
                    'QF_traffic', 'QF_building', 'QF_industry', 'dT_cooling']}
        self.ucm.reset(forcing['T_air'][0])

        for t in range(nt):
            hour = t % 24
            qf = self.anthro.total_anthropogenic_heat(hour, forcing['T_air'][t], day_of_week=2, month=8)
            dT_cool = 0.0
            if self.cooling:
                cool = self.cooling.total_cooling_effect(forcing['sw_down'][t], forcing['T_air'][t])
                dT_cool = cool['dT_net']
            eb = self.ucm.energy_balance(
                sw_down=forcing['sw_down'][t], lw_down=forcing['lw_down'][t],
                T_air=forcing['T_air'][t] + dT_cool,
                u_star=forcing['u_star'][t], q_air=forcing['q_air'][t],
                Q_anthropogenic=qf['QF_total'], dt=dt)
            for key in ['Q_star', 'QF', 'QH', 'QE', 'dQS',
                        'T_canyon', 'T_roof', 'T_wall', 'T_road', 'UHI_intensity']:
                results[key][t] = eb[key]
            results['QF_traffic'][t] = qf['QF_traffic']
            results['QF_building'][t] = qf['QF_building']
            results['QF_industry'][t] = qf['QF_industry']
            results['dT_cooling'][t] = dT_cool
        return results

    def run_scenario_comparison(self, forcing, scenarios, dt=3600.0):
        all_results = {}
        for name, cooling_model in scenarios.items():
            self.cooling = cooling_model
            self.ucm.reset(forcing['T_air'][0])
            all_results[name] = self.run_diurnal_cycle(forcing, dt)
        return all_results


class ENVImetInterface:
    def __init__(self, domain_size=(100, 100, 30), resolution=2.0):
        self.nx, self.ny, self.nz = domain_size
        self.resolution = resolution

    def generate_area_input(self, building_array, vegetation_array, surface_type):
        return {
            "model_area": {"nx": self.nx, "ny": self.ny, "nz": self.nz,
                           "dx": self.resolution, "dy": self.resolution},
            "buildings": building_array.tolist(),
            "vegetation": vegetation_array.tolist(),
            "surfaces": surface_type.tolist(),
        }

    def generate_simulation_config(self, start_date="2024-08-01", duration_hours=48):
        return {
            "simulation": {
                "start_date": start_date, "duration": duration_hours,
                "output_interval": 3600, "turbulence_model": "prognostic",
                "radiation_model": "IVS", "plant_model": True,
            },
            "output_variables": ["T_air", "T_surface", "T_mrt", "wind_speed",
                                  "relative_humidity", "UTCI", "PET", "WBGT"],
        }
