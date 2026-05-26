"""
IFC to Simulation Model Converter
Automated conversion of IFC/BIM data to thermal, CFD, and daylighting simulation models.
"""
import json
import math
from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class Material:
    name: str
    thermal_conductivity: float  # W/(m·K)
    density: float  # kg/m³
    specific_heat: float  # J/(kg·K)
    thickness: float  # m
    solar_absorptance: float = 0.7
    visible_absorptance: float = 0.7
    roughness: str = "MediumRough"

    @property
    def thermal_resistance(self) -> float:
        return self.thickness / self.thermal_conductivity

    @property
    def u_value(self) -> float:
        return 1.0 / (self.thermal_resistance + 0.13 + 0.04)


@dataclass
class Window:
    name: str
    width: float  # m
    height: float  # m
    u_value: float  # W/(m²·K)
    shgc: float  # Solar Heat Gain Coefficient
    vlt: float  # Visible Light Transmittance
    orientation: float  # degrees from north

    @property
    def area(self) -> float:
        return self.width * self.height


@dataclass
class ThermalZone:
    name: str
    floor_area: float  # m²
    volume: float  # m³
    height: float  # m
    walls: List[Dict] = field(default_factory=list)
    windows: List[Window] = field(default_factory=list)
    occupancy_density: float = 0.1  # people/m²
    lighting_density: float = 10.0  # W/m²
    equipment_density: float = 15.0  # W/m²
    ventilation_rate: float = 0.006  # m³/s per person

    @property
    def total_window_area(self) -> float:
        return sum(w.area for w in self.windows)

    @property
    def wwr(self) -> float:
        total_wall = sum(w.get('area', 0) for w in self.walls)
        return self.total_window_area / total_wall if total_wall > 0 else 0


@dataclass
class BuildingModel:
    name: str
    location: str
    latitude: float
    longitude: float
    climate_zone: str
    zones: List[ThermalZone] = field(default_factory=list)
    materials: List[Material] = field(default_factory=list)

    @property
    def total_floor_area(self) -> float:
        return sum(z.floor_area for z in self.zones)

    @property
    def total_volume(self) -> float:
        return sum(z.volume for z in self.zones)


class IFCConverter:
    """Converts IFC building data to simulation-ready models."""

    # Standard material library
    MATERIAL_LIBRARY = {
        "concrete_200": Material("Concrete200", 1.4, 2300, 880, 0.200),
        "concrete_150": Material("Concrete150", 1.4, 2300, 880, 0.150),
        "insulation_xps_50": Material("XPS50", 0.034, 35, 1400, 0.050),
        "insulation_xps_100": Material("XPS100", 0.034, 35, 1400, 0.100),
        "insulation_rockwool_100": Material("Rockwool100", 0.038, 100, 840, 0.100),
        "gypsum_board": Material("GypsumBoard", 0.16, 800, 1090, 0.0125),
        "plywood": Material("Plywood", 0.15, 560, 1210, 0.012),
        "air_gap": Material("AirGap", 0.025, 1.2, 1006, 0.05),
        "glass_double_lowe": Material("DoubleLowE", 1.1, 2500, 840, 0.024),
    }

    WINDOW_LIBRARY = {
        "standard_double": Window("StdDouble", 1.5, 1.8, 2.7, 0.76, 0.70, 0),
        "low_e_double": Window("LowEDouble", 1.5, 1.8, 1.6, 0.40, 0.65, 0),
        "triple_low_e": Window("TripleLowE", 1.5, 1.8, 0.8, 0.25, 0.55, 0),
        "high_perf": Window("HighPerf", 1.5, 1.8, 0.7, 0.22, 0.50, 0),
    }

    def __init__(self):
        self.model = None

    def create_reference_building(self, building_type: str = "office") -> BuildingModel:
        """Create a reference building model for ZEB case study."""
        model = BuildingModel(
            name="ZEB_Office_CaseStudy",
            location="Tokyo, Japan",
            latitude=35.68,
            longitude=139.77,
            climate_zone="4A"
        )

        # Standard materials
        model.materials = [
            self.MATERIAL_LIBRARY["concrete_200"],
            self.MATERIAL_LIBRARY["insulation_xps_100"],
            self.MATERIAL_LIBRARY["gypsum_board"],
        ]

        # Create zones for a 3-story office building
        floor_configs = [
            ("Ground_Floor", 500, 3.5),
            ("Second_Floor", 500, 3.2),
            ("Third_Floor", 500, 3.2),
        ]

        orientations = [0, 90, 180, 270]  # N, E, S, W
        wall_lengths = [25.0, 20.0, 25.0, 20.0]  # rectangular 25m x 20m

        for floor_name, area, height in floor_configs:
            zone = ThermalZone(
                name=floor_name,
                floor_area=area,
                volume=area * height,
                height=height,
            )

            for orient, length in zip(orientations, wall_lengths):
                wall_area = length * height
                zone.walls.append({
                    'orientation': orient,
                    'area': wall_area,
                    'u_value': 0.35,
                })
                # Windows: higher WWR on south, lower on north
                wwr = 0.45 if orient == 180 else (0.35 if orient in [90, 270] else 0.25)
                win_area = wall_area * wwr
                n_windows = max(1, int(win_area / (1.5 * 1.8)))
                for i in range(n_windows):
                    win = Window(
                        name=f"{floor_name}_Win_{orient}_{i}",
                        width=1.5,
                        height=1.8,
                        u_value=1.6,
                        shgc=0.40,
                        vlt=0.65,
                        orientation=orient
                    )
                    zone.windows.append(win)

            model.zones.append(zone)

        self.model = model
        return model

    def generate_energyplus_params(self) -> Dict:
        """Generate EnergyPlus simulation parameters from the building model."""
        if not self.model:
            raise ValueError("No building model loaded")

        params = {
            "simulation_control": {
                "run_period": "Annual",
                "timestep": 6,  # per hour
                "solar_distribution": "FullExteriorWithReflections",
            },
            "location": {
                "city": self.model.location,
                "latitude": self.model.latitude,
                "longitude": self.model.longitude,
                "timezone": 9,
                "elevation": 40,
            },
            "zones": [],
            "hvac": {
                "system_type": "VAV_with_Reheat",
                "cooling_cop": 3.5,
                "heating_cop": 4.0,
                "fan_efficiency": 0.7,
                "heat_recovery_effectiveness": 0.75,
            },
        }

        for zone in self.model.zones:
            zone_data = {
                "name": zone.name,
                "floor_area": zone.floor_area,
                "volume": zone.volume,
                "occupancy_density": zone.occupancy_density,
                "lighting_density": zone.lighting_density,
                "equipment_density": zone.equipment_density,
                "infiltration_rate": 0.3,  # ACH
                "ventilation_rate_per_person": zone.ventilation_rate,
                "heating_setpoint": 20.0,
                "cooling_setpoint": 26.0,
                "total_window_area": zone.total_window_area,
                "wwr": zone.wwr,
            }
            params["zones"].append(zone_data)

        return params

    def generate_cfd_params(self) -> Dict:
        """Generate CFD simulation parameters for natural ventilation analysis."""
        if not self.model:
            raise ValueError("No building model loaded")

        zone = self.model.zones[0]  # Use ground floor for CFD
        return {
            "domain": {
                "length": 25.0,
                "width": 20.0,
                "height": zone.height,
                "mesh_resolution": 0.25,
            },
            "boundary_conditions": {
                "inlet_velocity": 3.0,  # m/s
                "inlet_direction": 180,  # from south
                "outlet_pressure": 0,
                "temperature_outdoor": 28.0,
                "temperature_indoor": 26.0,
            },
            "openings": [
                {"wall": "south", "width": 1.5, "height": 1.2, "sill_height": 0.9,
                 "discharge_coeff": 0.65, "count": 4},
                {"wall": "north", "width": 1.5, "height": 1.2, "sill_height": 0.9,
                 "discharge_coeff": 0.65, "count": 3},
            ],
            "turbulence_model": "k-epsilon",
            "solver": "simpleFoam",
            "iterations": 2000,
            "convergence_criterion": 1e-5,
        }

    def generate_radiance_params(self) -> Dict:
        """Generate Radiance/Honeybee daylighting simulation parameters."""
        if not self.model:
            raise ValueError("No building model loaded")

        return {
            "location": {
                "latitude": self.model.latitude,
                "longitude": self.model.longitude,
                "timezone": 9,
            },
            "radiance_parameters": {
                "ab": 5,   # ambient bounces
                "ad": 2048,  # ambient divisions
                "as_": 1024,  # ambient super-samples
                "ar": 512,  # ambient resolution
                "aa": 0.1,  # ambient accuracy
            },
            "analysis_grid": {
                "height": 0.8,  # work plane height
                "spacing": 0.5,  # grid spacing in m
            },
            "materials": {
                "wall_reflectance": 0.5,
                "ceiling_reflectance": 0.8,
                "floor_reflectance": 0.2,
                "glass_transmittance": 0.65,
            },
            "metrics": ["sDA300/50", "ASE1000/250", "UDI100-3000", "DA300"],
            "schedules": {
                "occupancy_start": 8,
                "occupancy_end": 18,
            },
        }

    def export_model_summary(self) -> Dict:
        """Export a comprehensive model summary."""
        if not self.model:
            raise ValueError("No building model loaded")

        return {
            "building_name": self.model.name,
            "location": self.model.location,
            "climate_zone": self.model.climate_zone,
            "total_floor_area_m2": self.model.total_floor_area,
            "total_volume_m3": self.model.total_volume,
            "num_zones": len(self.model.zones),
            "zones": [
                {
                    "name": z.name,
                    "floor_area": z.floor_area,
                    "volume": z.volume,
                    "num_windows": len(z.windows),
                    "total_window_area": round(z.total_window_area, 2),
                    "wwr": round(z.wwr, 3),
                }
                for z in self.model.zones
            ],
        }


if __name__ == "__main__":
    converter = IFCConverter()
    model = converter.create_reference_building("office")
    summary = converter.export_model_summary()
    print(json.dumps(summary, indent=2))
    ep_params = converter.generate_energyplus_params()
    print(f"\nEnergyPlus zones: {len(ep_params['zones'])}")
    cfd_params = converter.generate_cfd_params()
    print(f"CFD domain: {cfd_params['domain']}")
    rad_params = converter.generate_radiance_params()
    print(f"Radiance metrics: {rad_params['metrics']}")
