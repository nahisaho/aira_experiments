"""
Anthropogenic Heat Emission Model

Spatiotemporal modeling of anthropogenic heat flux (QF) from three sectors:
  - Transportation (QF_traffic)
  - Building HVAC / air conditioning (QF_building)
  - Industrial processes (QF_industry)

Based on Sailor & Lu (2004) inventory approach with Tokyo-specific profiles.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional


@dataclass
class TrafficHeatParams:
    road_density: float = 5000.0
    vehicle_density_peak: float = 2000.0
    avg_heat_per_vehicle: float = 3.0
    ev_fraction: float = 0.10
    ev_heat_factor: float = 0.3


@dataclass
class BuildingHeatParams:
    floor_area_ratio: float = 3.0
    cooling_intensity: float = 80.0
    heating_intensity: float = 50.0
    cop_cooling: float = 3.5
    cop_heating: float = 4.0
    internal_gains: float = 25.0
    waste_heat_fraction: float = 0.3


@dataclass
class IndustryHeatParams:
    industry_fraction: float = 0.05
    heat_intensity: float = 50.0
    operating_hours: tuple = (8, 20)


class AnthropogenicHeatModel:
    """
    Computes spatiotemporal anthropogenic heat flux.
    QF = QF_traffic + QF_building + QF_industry
    """

    TRAFFIC_DIURNAL = np.array([
        0.15, 0.10, 0.08, 0.08, 0.12, 0.30,
        0.60, 1.20, 1.80, 1.50, 1.30, 1.20,
        1.25, 1.20, 1.30, 1.40, 1.60, 1.80,
        1.50, 1.10, 0.80, 0.60, 0.40, 0.25,
    ])
    TRAFFIC_DIURNAL = TRAFFIC_DIURNAL / TRAFFIC_DIURNAL.mean()

    BUILDING_DIURNAL_OFFICE = np.array([
        0.30, 0.25, 0.25, 0.25, 0.25, 0.30,
        0.50, 0.80, 1.20, 1.50, 1.60, 1.60,
        1.50, 1.60, 1.60, 1.50, 1.40, 1.20,
        0.90, 0.70, 0.50, 0.40, 0.35, 0.30,
    ])
    BUILDING_DIURNAL_OFFICE = BUILDING_DIURNAL_OFFICE / BUILDING_DIURNAL_OFFICE.mean()

    BUILDING_DIURNAL_RESIDENTIAL = np.array([
        0.40, 0.35, 0.30, 0.30, 0.30, 0.40,
        0.70, 1.00, 0.80, 0.60, 0.60, 0.70,
        0.80, 0.80, 0.80, 0.90, 1.10, 1.40,
        1.60, 1.70, 1.60, 1.30, 0.90, 0.60,
    ])
    BUILDING_DIURNAL_RESIDENTIAL = BUILDING_DIURNAL_RESIDENTIAL / BUILDING_DIURNAL_RESIDENTIAL.mean()

    DOW_FACTOR = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.7, 0.5])

    def __init__(self, traffic=None, building=None, industry=None, land_use="commercial"):
        self.traffic = traffic or TrafficHeatParams()
        self.building = building or BuildingHeatParams()
        self.industry = industry or IndustryHeatParams()
        self.land_use = land_use

    def compute_traffic_heat(self, hour, day_of_week=2, month=8):
        base = (self.traffic.vehicle_density_peak * 1e-6
                * self.traffic.avg_heat_per_vehicle * 1000)
        ice_fraction = 1 - self.traffic.ev_fraction
        effective = base * (ice_fraction + self.traffic.ev_fraction * self.traffic.ev_heat_factor)
        temporal = self.TRAFFIC_DIURNAL[hour % 24] * self.DOW_FACTOR[day_of_week % 7]
        return effective * temporal

    def compute_building_heat(self, hour, T_outdoor, day_of_week=2, month=8):
        profile = (self.BUILDING_DIURNAL_OFFICE if self.land_use == "commercial"
                   else self.BUILDING_DIURNAL_RESIDENTIAL)
        temporal = profile[hour % 24] * self.DOW_FACTOR[day_of_week % 7]
        T_outdoor_C = T_outdoor - 273.15 if T_outdoor > 200 else T_outdoor
        far = self.building.floor_area_ratio

        if T_outdoor_C > 26:
            cooling_load = self.building.cooling_intensity * (T_outdoor_C - 26) / 10
            cooling_load = min(cooling_load, self.building.cooling_intensity)
            waste_heat = cooling_load * (1 + 1/self.building.cop_cooling) * self.building.waste_heat_fraction
        elif T_outdoor_C < 10:
            heating_load = self.building.heating_intensity * (10 - T_outdoor_C) / 15
            waste_heat = heating_load / self.building.cop_heating * self.building.waste_heat_fraction
        else:
            waste_heat = self.building.internal_gains * 0.1

        return waste_heat * far * temporal

    def compute_industry_heat(self, hour, day_of_week=2):
        if day_of_week >= 5:
            return 0.0
        h_start, h_end = self.industry.operating_hours
        if h_start <= hour < h_end:
            return self.industry.heat_intensity * self.industry.industry_fraction
        return self.industry.heat_intensity * self.industry.industry_fraction * 0.1

    def total_anthropogenic_heat(self, hour, T_outdoor, day_of_week=2, month=8):
        qf_t = self.compute_traffic_heat(hour, day_of_week, month)
        qf_b = self.compute_building_heat(hour, T_outdoor, day_of_week, month)
        qf_i = self.compute_industry_heat(hour, day_of_week)
        return {
            "QF_traffic": qf_t, "QF_building": qf_b,
            "QF_industry": qf_i, "QF_total": qf_t + qf_b + qf_i,
        }


TOKYO_HEAT_PROFILES = {
    "marunouchi": AnthropogenicHeatModel(
        traffic=TrafficHeatParams(road_density=8000, vehicle_density_peak=3000),
        building=BuildingHeatParams(floor_area_ratio=8.0, cooling_intensity=100,
                                     cop_cooling=4.0, internal_gains=35),
        industry=IndustryHeatParams(industry_fraction=0.01, heat_intensity=20),
        land_use="commercial"),
    "shinjuku": AnthropogenicHeatModel(
        traffic=TrafficHeatParams(road_density=7000, vehicle_density_peak=2500),
        building=BuildingHeatParams(floor_area_ratio=6.0, cooling_intensity=90,
                                     cop_cooling=3.8, internal_gains=30),
        land_use="commercial"),
    "residential_23ku": AnthropogenicHeatModel(
        traffic=TrafficHeatParams(road_density=4000, vehicle_density_peak=1200),
        building=BuildingHeatParams(floor_area_ratio=1.5, cooling_intensity=60,
                                     cop_cooling=3.0, internal_gains=15),
        land_use="residential"),
    "suburban": AnthropogenicHeatModel(
        traffic=TrafficHeatParams(road_density=3000, vehicle_density_peak=800),
        building=BuildingHeatParams(floor_area_ratio=0.8, cooling_intensity=50,
                                     cop_cooling=3.0, internal_gains=10),
        land_use="residential"),
}


def project_anthropogenic_heat_2050(model, ev_penetration=0.80,
                                     cop_improvement=1.3, population_change=0.85):
    from copy import deepcopy
    projected = deepcopy(model)
    projected.traffic.ev_fraction = ev_penetration
    projected.building.cop_cooling *= cop_improvement
    projected.building.cop_heating *= cop_improvement
    projected.building.cooling_intensity *= 1.15
    projected.building.floor_area_ratio *= population_change
    return projected
