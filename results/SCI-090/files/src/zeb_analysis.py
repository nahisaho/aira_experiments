"""
ZEB (Net Zero Energy Building) Design Analysis
Evaluates energy balance, renewable energy generation, and ZEB compliance.
"""
import numpy as np

class ZEBAnalysis:
    """Net Zero Energy Building evaluation and optimization."""
    
    def __init__(self, thermal_results, building_config):
        self.thermal = thermal_results
        self.config = building_config
        self.floor_area = building_config["total_floor_area"]
        
    def calculate_primary_energy(self):
        """Calculate primary energy consumption with system efficiencies."""
        # HVAC system efficiencies
        cop_cooling = 4.5  # High-efficiency VRF
        cop_heating = 3.8  # Heat pump COP
        
        # Electricity for HVAC
        e_cooling = self.thermal["annual_cooling_kWh"] / cop_cooling
        e_heating = self.thermal["annual_heating_kWh"] / cop_heating
        
        # Lighting (with daylight-responsive control, 30% savings)
        e_lighting = self.thermal["annual_lighting_kWh"] * 0.70
        
        # Equipment
        e_equipment = self.thermal["annual_equipment_kWh"]
        
        # Ventilation fans, pumps
        e_aux = self.floor_area * 5.0  # 5 kWh/m2/yr
        
        # DHW (domestic hot water via heat pump, COP=3.0)
        dhw_demand = self.floor_area * 10.0  # 10 kWh/m2/yr
        e_dhw = dhw_demand / 3.0
        
        total_electricity = e_cooling + e_heating + e_lighting + e_equipment + e_aux + e_dhw
        
        # Primary energy factor (grid electricity)
        pef = 2.0  # Japanese grid
        primary_energy = total_electricity * pef
        
        self.energy_breakdown = {
            "cooling_kWh": e_cooling,
            "heating_kWh": e_heating,
            "lighting_kWh": e_lighting,
            "equipment_kWh": e_equipment,
            "auxiliary_kWh": e_aux,
            "dhw_kWh": e_dhw,
            "total_electricity_kWh": total_electricity,
            "primary_energy_kWh": primary_energy,
            "EUI_electricity": total_electricity / self.floor_area,
            "EUI_primary": primary_energy / self.floor_area,
        }
        return self.energy_breakdown
    
    def design_pv_system(self, weather_ghi):
        """Design rooftop PV system for ZEB compliance."""
        roof_area = 800  # m2 (top floor area)
        usable_ratio = 0.70  # accounting for equipment, access
        pv_area = roof_area * usable_ratio
        
        # PV parameters
        pv_efficiency = 0.22  # High-efficiency monocrystalline
        performance_ratio = 0.80  # system losses, inverter, wiring
        tilt = 30  # degrees
        
        # Hourly PV generation
        hours = len(weather_ghi)
        pv_generation = np.zeros(hours)
        
        for h in range(hours):
            # Tilt factor approximation
            hour = h % 24
            day = h / 24.0
            solar_altitude = max(0, np.sin(2 * np.pi * (hour - 6) / 24))
            tilt_factor = 1.0 + 0.3 * np.sin(np.radians(tilt))
            
            irradiance = weather_ghi[h] * tilt_factor
            pv_generation[h] = pv_area * pv_efficiency * performance_ratio * irradiance / 1000
        
        annual_generation = np.sum(pv_generation)
        
        # Monthly generation
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        monthly_gen = []
        start = 0
        for days in days_in_month:
            end = start + days * 24
            monthly_gen.append(np.sum(pv_generation[start:end]))
            start = end
        
        self.pv_results = {
            "pv_area_m2": pv_area,
            "pv_efficiency": pv_efficiency,
            "annual_generation_kWh": annual_generation,
            "monthly_generation": monthly_gen,
            "hourly_generation": pv_generation,
            "capacity_kWp": pv_area * pv_efficiency,
        }
        return self.pv_results
    
    def evaluate_zeb_compliance(self):
        """Evaluate ZEB compliance (energy balance)."""
        total_consumption = self.energy_breakdown["total_electricity_kWh"]
        total_generation = self.pv_results["annual_generation_kWh"]
        
        net_energy = total_consumption - total_generation
        zeb_ratio = total_generation / total_consumption
        
        # Monthly balance
        monthly_consumption = []
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        daily_consumption = total_consumption / 365
        for days in days_in_month:
            monthly_consumption.append(daily_consumption * days)
        
        monthly_balance = [
            g - c for g, c in zip(self.pv_results["monthly_generation"], monthly_consumption)
        ]
        
        # ZEB classification
        if zeb_ratio >= 1.0:
            zeb_class = "ZEB (Net Zero)"
        elif zeb_ratio >= 0.75:
            zeb_class = "Nearly ZEB"
        elif zeb_ratio >= 0.50:
            zeb_class = "ZEB Ready"
        else:
            zeb_class = "Below ZEB Ready"
        
        self.zeb_results = {
            "total_consumption_kWh": total_consumption,
            "total_generation_kWh": total_generation,
            "net_energy_kWh": net_energy,
            "zeb_ratio": zeb_ratio,
            "zeb_classification": zeb_class,
            "monthly_consumption": monthly_consumption,
            "monthly_generation": self.pv_results["monthly_generation"],
            "monthly_balance": monthly_balance,
            "EUI_net": net_energy / self.floor_area,
        }
        return self.zeb_results
    
    def optimization_scenarios(self, weather_ghi):
        """Evaluate different design scenarios for ZEB achievement."""
        scenarios = [
            {"name": "Baseline", "cop_c": 3.5, "cop_h": 3.0, "pv_eff": 0.18, "daylight_saving": 0.0},
            {"name": "High-COP HVAC", "cop_c": 5.0, "cop_h": 4.5, "pv_eff": 0.18, "daylight_saving": 0.0},
            {"name": "Premium PV", "cop_c": 4.5, "cop_h": 3.8, "pv_eff": 0.24, "daylight_saving": 0.0},
            {"name": "Daylight Control", "cop_c": 4.5, "cop_h": 3.8, "pv_eff": 0.22, "daylight_saving": 0.30},
            {"name": "Full Optimization", "cop_c": 5.0, "cop_h": 4.5, "pv_eff": 0.24, "daylight_saving": 0.35},
        ]
        
        results = []
        for sc in scenarios:
            e_cooling = self.thermal["annual_cooling_kWh"] / sc["cop_c"]
            e_heating = self.thermal["annual_heating_kWh"] / sc["cop_h"]
            e_lighting = self.thermal["annual_lighting_kWh"] * (1 - sc["daylight_saving"])
            e_total = e_cooling + e_heating + e_lighting + self.thermal["annual_equipment_kWh"] + self.floor_area * 5.0 + self.floor_area * 10.0 / 3.0
            
            pv_area = 800 * 0.70
            pv_gen = np.sum(weather_ghi) * pv_area * sc["pv_eff"] * 0.80 * 1.3 / 1000
            
            results.append({
                "scenario": sc["name"],
                "consumption_kWh": e_total,
                "generation_kWh": pv_gen,
                "zeb_ratio": pv_gen / e_total,
                "EUI_net": (e_total - pv_gen) / self.floor_area,
            })
        
        return results
