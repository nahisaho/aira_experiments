"""
Thermal Load Simulation Engine
Simplified EnergyPlus-compatible thermal simulation for building energy analysis.
"""
import math
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple


@dataclass
class WeatherData:
    """Hourly weather data for Tokyo (simplified TMY)."""
    month: int
    hour: int
    dry_bulb_temp: float  # °C
    relative_humidity: float  # %
    global_horizontal_radiation: float  # W/m²
    direct_normal_radiation: float  # W/m²
    diffuse_horizontal_radiation: float  # W/m²
    wind_speed: float  # m/s
    wind_direction: float  # degrees


class ThermalSimulation:
    """Simplified thermal load calculation engine compatible with EnergyPlus methods."""

    def __init__(self, building_params: Dict):
        self.params = building_params
        self.results = {}

    def _generate_tokyo_weather(self) -> List[WeatherData]:
        """Generate simplified Tokyo weather data (TMY-based)."""
        monthly_temp = [5.2, 5.7, 8.7, 13.9, 18.2, 21.4, 25.0, 26.4, 22.8, 17.5, 12.1, 7.6]
        monthly_rh = [52, 56, 62, 67, 72, 80, 82, 78, 75, 68, 60, 55]
        monthly_solar_peak = [350, 420, 500, 550, 580, 520, 600, 580, 450, 380, 330, 300]
        daily_range = [7.5, 8.0, 8.5, 9.0, 8.5, 7.0, 7.5, 8.0, 7.5, 8.0, 8.0, 7.5]

        weather = []
        for m in range(12):
            for h in range(24):
                # Diurnal temperature variation
                hour_factor = math.sin(math.pi * (h - 6) / 12) if 6 <= h <= 18 else -0.5
                temp = monthly_temp[m] + daily_range[m] * hour_factor * 0.5

                # Solar radiation (simplified)
                if 6 <= h <= 18:
                    solar_factor = max(0, math.sin(math.pi * (h - 6) / 12))
                    ghr = monthly_solar_peak[m] * solar_factor
                    dnr = ghr * 0.6
                    dhr = ghr * 0.4
                else:
                    ghr = dnr = dhr = 0

                weather.append(WeatherData(
                    month=m + 1, hour=h,
                    dry_bulb_temp=temp,
                    relative_humidity=monthly_rh[m],
                    global_horizontal_radiation=ghr,
                    direct_normal_radiation=dnr,
                    diffuse_horizontal_radiation=dhr,
                    wind_speed=3.5,
                    wind_direction=180
                ))
        return weather

    def calculate_envelope_heat_transfer(self, zone: Dict, outdoor_temp: float) -> float:
        """Calculate conduction heat transfer through building envelope."""
        indoor_temp = (zone['heating_setpoint'] + zone['cooling_setpoint']) / 2
        delta_t = outdoor_temp - indoor_temp

        # Wall conduction (U-value based)
        wall_area = zone['floor_area'] * 0.8  # approximate wall area
        window_area = zone['total_window_area']
        opaque_area = wall_area - window_area
        wall_u = 0.35  # W/(m²·K) - insulated wall

        q_wall = opaque_area * wall_u * delta_t
        q_window = window_area * 1.6 * delta_t  # Window U-value
        q_roof = zone['floor_area'] * 0.2 * delta_t * 0.3  # top floor factor

        return q_wall + q_window + q_roof

    def calculate_solar_gains(self, zone: Dict, solar_radiation: float,
                               month: int) -> float:
        """Calculate solar heat gains through windows."""
        shgc = 0.40
        window_area = zone['total_window_area']
        solar_angle_factor = self._solar_orientation_factor(month)
        return window_area * shgc * solar_radiation * solar_angle_factor * 0.3

    def _solar_orientation_factor(self, month: int) -> float:
        """Simplified solar orientation factor based on season."""
        summer_months = [6, 7, 8]
        winter_months = [12, 1, 2]
        if month in summer_months:
            return 0.7
        elif month in winter_months:
            return 1.2
        return 0.9

    def calculate_internal_gains(self, zone: Dict, hour: int) -> float:
        """Calculate internal heat gains from occupants, lighting, equipment."""
        # Occupancy schedule
        if 8 <= hour <= 18:
            occ_factor = 1.0 if 9 <= hour <= 17 else 0.5
        else:
            occ_factor = 0.05

        people_heat = zone['floor_area'] * zone['occupancy_density'] * 120 * occ_factor  # 120 W/person
        lighting_heat = zone['floor_area'] * zone['lighting_density'] * occ_factor
        equipment_heat = zone['floor_area'] * zone['equipment_density'] * occ_factor

        return people_heat + lighting_heat + equipment_heat

    def calculate_ventilation_load(self, zone: Dict, outdoor_temp: float) -> float:
        """Calculate ventilation/infiltration thermal load."""
        indoor_temp = (zone['heating_setpoint'] + zone['cooling_setpoint']) / 2
        delta_t = outdoor_temp - indoor_temp

        # Mechanical ventilation
        n_people = zone['floor_area'] * zone['occupancy_density']
        vent_flow = n_people * zone['ventilation_rate_per_person']
        q_vent = vent_flow * 1.2 * 1006 * delta_t

        # Infiltration
        infiltration_ach = zone['infiltration_rate']
        inf_flow = zone['volume'] * infiltration_ach / 3600
        q_inf = inf_flow * 1.2 * 1006 * delta_t

        # Heat recovery
        hr_eff = self.params['hvac']['heat_recovery_effectiveness']
        q_vent_recovered = q_vent * hr_eff

        return (q_vent - q_vent_recovered) + q_inf

    def run_annual_simulation(self) -> Dict:
        """Run annual hourly thermal simulation."""
        weather = self._generate_tokyo_weather()
        zones = self.params['zones']

        monthly_heating = np.zeros(12)
        monthly_cooling = np.zeros(12)
        monthly_total_energy = np.zeros(12)
        hourly_loads = {z['name']: {'heating': [], 'cooling': []} for z in zones}

        for wd in weather:
            for zone in zones:
                q_envelope = self.calculate_envelope_heat_transfer(zone, wd.dry_bulb_temp)
                q_solar = self.calculate_solar_gains(zone, wd.global_horizontal_radiation, wd.month)
                q_internal = self.calculate_internal_gains(zone, wd.hour)
                q_vent = self.calculate_ventilation_load(zone, wd.dry_bulb_temp)

                total_load = q_envelope + q_solar + q_internal + q_vent

                if total_load > 0:
                    # Cooling needed
                    cooling_energy = total_load / (self.params['hvac']['cooling_cop'] * 1000)
                    monthly_cooling[wd.month - 1] += cooling_energy
                    hourly_loads[zone['name']]['cooling'].append(cooling_energy)
                    hourly_loads[zone['name']]['heating'].append(0)
                else:
                    # Heating needed
                    heating_energy = abs(total_load) / (self.params['hvac']['heating_cop'] * 1000)
                    monthly_heating[wd.month - 1] += heating_energy
                    hourly_loads[zone['name']]['heating'].append(heating_energy)
                    hourly_loads[zone['name']]['cooling'].append(0)

        total_floor_area = sum(z['floor_area'] for z in zones)
        annual_heating = monthly_heating.sum()
        annual_cooling = monthly_cooling.sum()

        # Lighting and equipment energy
        annual_lighting = sum(z['floor_area'] * z['lighting_density'] * 2500 / 1000 for z in zones)
        annual_equipment = sum(z['floor_area'] * z['equipment_density'] * 2500 / 1000 for z in zones)
        # Fan energy
        annual_fan = total_floor_area * 15 * 2500 / 1000

        annual_total = annual_heating + annual_cooling + annual_lighting + annual_equipment + annual_fan

        self.results = {
            'monthly_heating_kWh': monthly_heating.tolist(),
            'monthly_cooling_kWh': monthly_cooling.tolist(),
            'annual_heating_kWh': round(annual_heating, 1),
            'annual_cooling_kWh': round(annual_cooling, 1),
            'annual_lighting_kWh': round(annual_lighting, 1),
            'annual_equipment_kWh': round(annual_equipment, 1),
            'annual_fan_kWh': round(annual_fan, 1),
            'annual_total_kWh': round(annual_total, 1),
            'eui_kWh_m2': round(annual_total / total_floor_area, 1),
            'heating_eui': round(annual_heating / total_floor_area, 1),
            'cooling_eui': round(annual_cooling / total_floor_area, 1),
            'peak_heating_kW': round(max(max(hourly_loads[z['name']]['heating']) for z in zones), 1),
            'peak_cooling_kW': round(max(max(hourly_loads[z['name']]['cooling']) for z in zones), 1),
            'total_floor_area': total_floor_area,
        }

        return self.results

    def run_zeb_optimization(self) -> Dict:
        """Run ZEB optimization with improved envelope and renewable energy."""
        # Baseline
        baseline = self.run_annual_simulation()
        baseline_eui = baseline['eui_kWh_m2']

        # ZEB improvements
        improvements = {}

        # 1. Enhanced insulation (reduce wall/roof U-value)
        insulation_reduction = 0.30
        improvements['enhanced_insulation'] = {
            'description': 'High-performance insulation (U=0.2 walls, U=0.15 roof)',
            'energy_reduction_pct': round(insulation_reduction * 100 / baseline_eui * baseline['heating_eui'], 1),
        }

        # 2. High-performance glazing
        glazing_reduction = 0.15
        improvements['high_perf_glazing'] = {
            'description': 'Triple low-E glazing (U=0.8, SHGC=0.25)',
            'energy_reduction_pct': round(glazing_reduction * 100 * (baseline['heating_eui'] + baseline['cooling_eui']) / baseline_eui, 1),
        }

        # 3. LED lighting
        led_reduction = 0.50
        improvements['led_lighting'] = {
            'description': 'LED lighting with daylight dimming (5 W/m²)',
            'energy_reduction_pct': round(led_reduction * baseline['annual_lighting_kWh'] / baseline['annual_total_kWh'] * 100, 1),
        }

        # 4. Heat recovery ventilation
        hr_improvement = 0.15
        improvements['heat_recovery'] = {
            'description': 'Enhanced HRV (90% effectiveness)',
            'energy_reduction_pct': round(hr_improvement * 100 * (baseline['heating_eui'] + baseline['cooling_eui']) / baseline_eui, 1),
        }

        # 5. Natural ventilation
        nat_vent_saving = 0.10
        improvements['natural_ventilation'] = {
            'description': 'Cross-ventilation with automated controls',
            'energy_reduction_pct': round(nat_vent_saving * baseline['annual_cooling_kWh'] / baseline['annual_total_kWh'] * 100, 1),
        }

        total_reduction_pct = sum(imp['energy_reduction_pct'] for imp in improvements.values())
        total_reduction_pct = min(total_reduction_pct, 65)  # cap at realistic maximum

        optimized_eui = baseline_eui * (1 - total_reduction_pct / 100)
        remaining_energy = baseline['annual_total_kWh'] * (1 - total_reduction_pct / 100)

        # PV generation
        pv_area = baseline['total_floor_area'] * 0.6  # 60% of roof area usable
        pv_efficiency = 0.20
        tokyo_solar = 1400  # kWh/m²/year
        pv_generation = pv_area * pv_efficiency * tokyo_solar
        pv_capacity_kw = pv_area * 0.20  # 200 W/m²

        zeb_balance = remaining_energy - pv_generation
        zeb_ratio = pv_generation / remaining_energy

        return {
            'baseline': baseline,
            'improvements': improvements,
            'total_reduction_pct': round(total_reduction_pct, 1),
            'optimized_eui_kWh_m2': round(optimized_eui, 1),
            'remaining_energy_kWh': round(remaining_energy, 1),
            'pv_area_m2': round(pv_area, 1),
            'pv_capacity_kW': round(pv_capacity_kw, 1),
            'pv_generation_kWh': round(pv_generation, 1),
            'zeb_balance_kWh': round(zeb_balance, 1),
            'zeb_ratio': round(zeb_ratio, 2),
            'is_zeb': zeb_ratio >= 1.0,
        }


if __name__ == "__main__":
    from ifc_converter import IFCConverter
    import json

    converter = IFCConverter()
    model = converter.create_reference_building()
    ep_params = converter.generate_energyplus_params()

    sim = ThermalSimulation(ep_params)
    results = sim.run_annual_simulation()
    print("=== Annual Thermal Simulation Results ===")
    print(json.dumps(results, indent=2))

    zeb = sim.run_zeb_optimization()
    print("\n=== ZEB Optimization Results ===")
    print(json.dumps(zeb, indent=2))
