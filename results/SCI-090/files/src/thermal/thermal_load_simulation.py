"""
熱負荷シミュレーションモジュール
EnergyPlus連携による年間熱負荷計算・HVAC性能評価

機能:
- 年間冷暖房負荷計算
- 月別・時間別エネルギー消費プロファイル
- HVAC機器容量算定
- ZEB達成度評価
"""

import numpy as np
import json
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class WeatherData:
    """気象データ（東京の代表的な年間データ）"""
    location: str = "Tokyo"
    latitude: float = 35.68
    longitude: float = 139.77

    # 月別平均外気温 [°C]
    monthly_temp: list[float] = field(default_factory=lambda: [
        5.2, 5.7, 8.7, 13.9, 18.2, 21.4,
        25.0, 26.4, 22.8, 17.5, 12.1, 7.6
    ])

    # 月別平均日射量 [kWh/m²/day]
    monthly_solar: list[float] = field(default_factory=lambda: [
        2.8, 3.2, 3.6, 4.2, 4.5, 3.8,
        4.1, 4.6, 3.4, 2.9, 2.6, 2.5
    ])

    # 月別平均湿度 [%]
    monthly_humidity: list[float] = field(default_factory=lambda: [
        52, 53, 56, 62, 69, 75,
        77, 73, 75, 68, 62, 56
    ])

    # 月間日数
    days_per_month: list[int] = field(default_factory=lambda: [
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    ])


@dataclass
class HVACSystem:
    """HVAC機器スペック"""
    name: str
    cooling_cop: float = 4.0
    heating_cop: float = 3.5
    fan_power: float = 0.5  # W/(m³/h)
    pump_power: float = 0.3  # W/(m³/h)
    heat_recovery_efficiency: float = 0.75
    supply_air_temp_cooling: float = 15.0  # °C
    supply_air_temp_heating: float = 35.0  # °C


@dataclass
class ThermalLoadResult:
    """熱負荷計算結果"""
    monthly_heating_load: list[float]  # kWh
    monthly_cooling_load: list[float]  # kWh
    monthly_lighting_energy: list[float]  # kWh
    monthly_equipment_energy: list[float]  # kWh
    monthly_fan_energy: list[float]  # kWh
    monthly_pump_energy: list[float]  # kWh
    peak_heating_load: float  # kW
    peak_cooling_load: float  # kW
    annual_primary_energy: float  # kWh/m²/yr
    zeb_score: float  # %


class ThermalLoadSimulator:
    """
    簡易熱負荷シミュレーター
    EnergyPlusの代理モデルとして月別定常計算を実施
    """

    # 建物使用スケジュール（平日）
    OCCUPANCY_SCHEDULE = [
        0, 0, 0, 0, 0, 0, 0, 0.5,
        1.0, 1.0, 1.0, 0.8, 0.5, 1.0, 1.0, 1.0,
        1.0, 0.5, 0.2, 0, 0, 0, 0, 0
    ]

    def __init__(self, weather: WeatherData, hvac: HVACSystem):
        self.weather = weather
        self.hvac = hvac

    def calculate_zone_load(
        self,
        floor_area: float,
        volume: float,
        wall_area: float,
        window_area: float,
        wall_u: float,
        window_u: float,
        window_shgc: float,
        occupancy_density: float,
        lighting_density: float,
        equipment_density: float,
        infiltration_ach: float,
        setpoint_heating: float,
        setpoint_cooling: float,
    ) -> ThermalLoadResult:
        """ゾーン熱負荷計算"""
        monthly_heating = []
        monthly_cooling = []
        monthly_lighting = []
        monthly_equipment = []
        monthly_fan = []
        monthly_pump = []

        peak_heating = 0
        peak_cooling = 0

        for month in range(12):
            t_out = self.weather.monthly_temp[month]
            solar = self.weather.monthly_solar[month]
            days = self.weather.days_per_month[month]
            weekdays = int(days * 5 / 7)

            # 内部発熱
            people_heat = occupancy_density * floor_area * 120  # W (120W/person)
            lighting_heat = lighting_density * floor_area  # W
            equipment_heat = equipment_density * floor_area  # W
            total_internal = people_heat + lighting_heat + equipment_heat

            # 日射取得
            solar_gain = window_area * solar * 1000 / 24 * window_shgc  # W (avg)

            # 外皮負荷（冷房）
            dt_cooling = max(0, t_out - setpoint_cooling)
            envelope_cooling = (
                wall_area * wall_u * dt_cooling
                + window_area * window_u * dt_cooling
            )

            # 外皮負荷（暖房）
            dt_heating = max(0, setpoint_heating - t_out)
            envelope_heating = (
                wall_area * wall_u * dt_heating
                + window_area * window_u * dt_heating
            )

            # 換気負荷
            air_density = 1.2  # kg/m³
            air_cp = 1005  # J/(kg·K)
            ventilation_flow = volume * infiltration_ach / 3600  # m³/s
            vent_cooling = max(0, ventilation_flow * air_density * air_cp * dt_cooling)
            vent_heating = max(0, ventilation_flow * air_density * air_cp * dt_heating)

            # 全熱換気なので回収分を差し引く
            vent_cooling *= (1 - self.hvac.heat_recovery_efficiency)
            vent_heating *= (1 - self.hvac.heat_recovery_efficiency)

            # 時間稼働率を考慮した月間負荷
            avg_occupancy = np.mean(self.OCCUPANCY_SCHEDULE)
            operating_hours = weekdays * 10  # hours/month

            # 冷房負荷
            cooling_load_w = max(0,
                envelope_cooling + solar_gain + total_internal * avg_occupancy
                + vent_cooling - envelope_heating
            )
            monthly_cooling_kwh = cooling_load_w * operating_hours / 1000
            monthly_cooling.append(monthly_cooling_kwh)

            # 暖房負荷
            heating_load_w = max(0,
                envelope_heating + vent_heating
                - total_internal * avg_occupancy * 0.5
                - solar_gain * 0.3
            )
            monthly_heating_kwh = heating_load_w * operating_hours / 1000
            monthly_heating.append(monthly_heating_kwh)

            # ピーク負荷更新
            peak_cooling = max(peak_cooling, cooling_load_w / 1000)
            peak_heating = max(peak_heating, heating_load_w / 1000)

            # 照明・機器エネルギー
            monthly_lighting.append(lighting_density * floor_area * operating_hours / 1000)
            monthly_equipment.append(equipment_density * floor_area * operating_hours / 1000)

            # 搬送エネルギー
            air_flow_rate = volume * 3  # m³/h
            fan_energy = self.hvac.fan_power * air_flow_rate * operating_hours / 1000
            pump_energy = self.hvac.pump_power * air_flow_rate * 0.3 * operating_hours / 1000
            monthly_fan.append(fan_energy)
            monthly_pump.append(pump_energy)

        # 年間一次エネルギー消費量
        total_cooling = sum(monthly_cooling) / self.hvac.cooling_cop
        total_heating = sum(monthly_heating) / self.hvac.heating_cop
        total_lighting = sum(monthly_lighting)
        total_equipment = sum(monthly_equipment)
        total_fan = sum(monthly_fan)
        total_pump = sum(monthly_pump)

        annual_total = (total_cooling + total_heating + total_lighting
                       + total_equipment + total_fan + total_pump)
        annual_per_area = annual_total / floor_area

        # ZEBスコア（BEI基準：基準一次エネルギーに対する削減率）
        baseline_energy = 300  # kWh/m²/yr (事務所標準)
        zeb_score = (1 - annual_per_area / baseline_energy) * 100

        return ThermalLoadResult(
            monthly_heating_load=monthly_heating,
            monthly_cooling_load=monthly_cooling,
            monthly_lighting_energy=monthly_lighting,
            monthly_equipment_energy=monthly_equipment,
            monthly_fan_energy=monthly_fan,
            monthly_pump_energy=monthly_pump,
            peak_heating_load=peak_heating,
            peak_cooling_load=peak_cooling,
            annual_primary_energy=annual_per_area,
            zeb_score=zeb_score,
        )


def run_thermal_simulation():
    """5階建ZEBオフィスの熱負荷シミュレーション"""
    weather = WeatherData()
    hvac = HVACSystem(
        name="VRF_with_DOAS",
        cooling_cop=5.5,
        heating_cop=4.2,
        fan_power=0.4,
        pump_power=0.2,
        heat_recovery_efficiency=0.80,
    )

    simulator = ThermalLoadSimulator(weather, hvac)

    # ビルパラメータ
    floor_area = 5000.0
    volume = 17500.0
    wall_area = 2215.0
    window_area = 886.0
    wall_u = 0.28
    window_u = 1.6
    window_shgc = 0.40

    result = simulator.calculate_zone_load(
        floor_area=floor_area,
        volume=volume,
        wall_area=wall_area,
        window_area=window_area,
        wall_u=wall_u,
        window_u=window_u,
        window_shgc=window_shgc,
        occupancy_density=0.1,
        lighting_density=8.0,  # LED照明
        equipment_density=12.0,  # 省エネ機器
        infiltration_ach=0.3,
        setpoint_heating=20.0,
        setpoint_cooling=26.0,
    )

    # 結果保存
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    output = {
        "simulation_type": "Annual Thermal Load Simulation",
        "building": {
            "name": "ZEB_Office_Demo",
            "floor_area_m2": floor_area,
            "volume_m3": volume,
            "hvac_system": hvac.name,
            "cooling_cop": hvac.cooling_cop,
            "heating_cop": hvac.heating_cop,
        },
        "results": {
            "peak_heating_kW": round(result.peak_heating_load, 1),
            "peak_cooling_kW": round(result.peak_cooling_load, 1),
            "annual_primary_energy_kWh_m2": round(result.annual_primary_energy, 1),
            "zeb_score_percent": round(result.zeb_score, 1),
            "monthly_data": {},
        },
    }

    for i, m in enumerate(months):
        output["results"]["monthly_data"][m] = {
            "heating_kWh": round(result.monthly_heating_load[i], 1),
            "cooling_kWh": round(result.monthly_cooling_load[i], 1),
            "lighting_kWh": round(result.monthly_lighting_energy[i], 1),
            "equipment_kWh": round(result.monthly_equipment_energy[i], 1),
            "fan_kWh": round(result.monthly_fan_energy[i], 1),
            "pump_kWh": round(result.monthly_pump_energy[i], 1),
        }

    with open("results/thermal_simulation_results.json", "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"=== 熱負荷シミュレーション結果 ===")
    print(f"ピーク暖房負荷: {result.peak_heating_load:.1f} kW")
    print(f"ピーク冷房負荷: {result.peak_cooling_load:.1f} kW")
    print(f"年間一次エネルギー消費量: {result.annual_primary_energy:.1f} kWh/m²/yr")
    print(f"ZEBスコア: {result.zeb_score:.1f}%")

    return result, output


if __name__ == "__main__":
    run_thermal_simulation()
