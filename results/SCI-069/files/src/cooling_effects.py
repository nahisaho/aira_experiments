"""
Green Infrastructure and High-Albedo Material Cooling Effect Quantification

Models cooling contributions from:
  - Urban greening (street trees, parks, green roofs)
  - High-albedo (cool) roofs and pavements
  - Water features (misting, retention ponds)
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict


@dataclass
class GreenInfrastructure:
    tree_coverage: float = 0.15
    grass_coverage: float = 0.10
    green_roof_coverage: float = 0.0
    park_area_fraction: float = 0.05
    leaf_area_index: float = 3.0
    stomatal_conductance: float = 0.01


@dataclass
class CoolMaterials:
    cool_roof_fraction: float = 0.0
    cool_roof_albedo: float = 0.60
    original_roof_albedo: float = 0.15
    cool_pavement_fraction: float = 0.0
    cool_pavement_albedo: float = 0.40
    original_pavement_albedo: float = 0.10


@dataclass
class WaterFeatures:
    misting_coverage: float = 0.0
    misting_rate: float = 0.5
    retention_pond_fraction: float = 0.0
    pervious_pavement_fraction: float = 0.0


class CoolingEffectModel:
    CP_AIR = 1005.0
    RHO_AIR = 1.2
    LATENT_HEAT = 2.45e6

    def __init__(self, green=None, cool=None, water=None):
        self.green = green or GreenInfrastructure()
        self.cool = cool or CoolMaterials()
        self.water = water or WaterFeatures()

    def tree_shade_cooling(self, sw_down, T_air):
        coverage = self.green.tree_coverage
        lai = self.green.leaf_area_index
        transmissivity = np.exp(-0.5 * lai)
        shade_fraction = (1 - transmissivity) * coverage
        sw_reduction = sw_down * shade_fraction
        vpd = 1.5
        ga = 0.05
        gs = self.green.stomatal_conductance
        et_rate = coverage * gs * ga / (gs + ga) * vpd * 0.622 / 101.3
        latent_flux = et_rate * self.LATENT_HEAT
        dT_shade = -shade_fraction * sw_down / (self.RHO_AIR * self.CP_AIR * 50)
        dT_et = -latent_flux / (self.RHO_AIR * self.CP_AIR * 100)
        return {
            "dT_shade": dT_shade, "dT_evapotranspiration": dT_et,
            "dT_total": dT_shade + dT_et,
            "sw_reduction": sw_reduction, "latent_flux": latent_flux,
            "shade_fraction": shade_fraction,
        }

    def green_roof_cooling(self, sw_down, T_air):
        coverage = self.green.green_roof_coverage
        if coverage == 0:
            return {"dT_roof_surface": 0.0, "dT_air": 0.0, "load_reduction_pct": 0.0}
        dT_roof = -25.0 * coverage
        dT_air = dT_roof * 0.08
        return {"dT_roof_surface": dT_roof, "dT_air": dT_air,
                "load_reduction_pct": coverage * 25.0}

    def cool_material_cooling(self, sw_down):
        delta_alpha_roof = self.cool.cool_roof_albedo - self.cool.original_roof_albedo
        sw_saved_roof = sw_down * delta_alpha_roof * self.cool.cool_roof_fraction
        dT_roof = -sw_saved_roof * 0.007
        delta_alpha_pave = self.cool.cool_pavement_albedo - self.cool.original_pavement_albedo
        sw_saved_pave = sw_down * delta_alpha_pave * self.cool.cool_pavement_fraction
        dT_pave = -sw_saved_pave * 0.005
        return {
            "dT_cool_roof": dT_roof, "dT_cool_pavement": dT_pave,
            "dT_total": dT_roof + dT_pave,
            "sw_saved_roof": sw_saved_roof, "sw_saved_pave": sw_saved_pave,
            "total_sw_reduction": sw_saved_roof + sw_saved_pave,
        }

    def water_feature_cooling(self, T_air, rh=0.60):
        if self.water.misting_coverage > 0:
            rate = self.water.misting_rate / 3600
            latent = rate * self.water.misting_coverage * self.LATENT_HEAT
            dT_mist = -latent / (self.RHO_AIR * self.CP_AIR * 200)
        else:
            dT_mist = 0.0
            latent = 0.0
        dT_pond = -2.0 * self.water.retention_pond_fraction
        return {"dT_misting": dT_mist, "dT_pond": dT_pond,
                "dT_total": dT_mist + dT_pond, "latent_flux": latent}

    def total_cooling_effect(self, sw_down, T_air, rh=0.60):
        tree = self.tree_shade_cooling(sw_down, T_air)
        green_roof = self.green_roof_cooling(sw_down, T_air)
        cool_mat = self.cool_material_cooling(sw_down)
        water = self.water_feature_cooling(T_air, rh)
        raw_sum = (tree["dT_total"] + green_roof["dT_air"]
                   + cool_mat["dT_total"] + water["dT_total"])
        dT_net = raw_sum * 0.85
        return {
            "dT_trees": tree["dT_total"], "dT_green_roof": green_roof["dT_air"],
            "dT_cool_materials": cool_mat["dT_total"],
            "dT_water_features": water["dT_total"],
            "dT_raw_sum": raw_sum, "dT_net": dT_net,
            "sw_total_reduction": tree["sw_reduction"] + cool_mat["total_sw_reduction"],
            "latent_total": tree["latent_flux"] + water.get("latent_flux", 0),
        }


MITIGATION_SCENARIOS = {
    "baseline": CoolingEffectModel(
        green=GreenInfrastructure(tree_coverage=0.10, grass_coverage=0.05,
                                   green_roof_coverage=0.0, park_area_fraction=0.05),
        cool=CoolMaterials(), water=WaterFeatures()),
    "moderate_greening": CoolingEffectModel(
        green=GreenInfrastructure(tree_coverage=0.20, grass_coverage=0.10,
                                   green_roof_coverage=0.15, park_area_fraction=0.08,
                                   leaf_area_index=4.0),
        cool=CoolMaterials(cool_roof_fraction=0.30, cool_pavement_fraction=0.20),
        water=WaterFeatures(misting_coverage=0.02, retention_pond_fraction=0.01)),
    "aggressive_mitigation": CoolingEffectModel(
        green=GreenInfrastructure(tree_coverage=0.30, grass_coverage=0.15,
                                   green_roof_coverage=0.40, park_area_fraction=0.12,
                                   leaf_area_index=5.0, stomatal_conductance=0.012),
        cool=CoolMaterials(cool_roof_fraction=0.60, cool_pavement_fraction=0.50,
                           cool_roof_albedo=0.65, cool_pavement_albedo=0.45),
        water=WaterFeatures(misting_coverage=0.05, misting_rate=0.8,
                            retention_pond_fraction=0.03, pervious_pavement_fraction=0.20)),
}
