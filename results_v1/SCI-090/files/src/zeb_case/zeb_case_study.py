"""
ZEB（ネットゼロエネルギービル）設計ケーススタディ
建築・設備・再エネの統合評価

機能:
- 基準建物 vs ZEB設計の比較
- 省エネ技術の個別・複合効果評価
- 太陽光発電システム設計
- ZEB達成度判定（ZEB, Nearly ZEB, ZEB Ready, ZEB Oriented）
"""

import numpy as np
import json
from dataclasses import dataclass, field


@dataclass
class BuildingDesign:
    """建物設計パラメータ"""
    name: str
    floor_area: float  # m²
    num_floors: int
    wall_u: float  # W/(m²·K)
    roof_u: float  # W/(m²·K)
    window_u: float  # W/(m²·K)
    window_shgc: float
    wwr: float  # Window-to-Wall Ratio
    infiltration_ach: float
    lighting_density: float  # W/m²
    equipment_density: float  # W/m²
    hvac_cooling_cop: float
    hvac_heating_cop: float
    heat_recovery_eff: float
    hot_water_efficiency: float
    elevator_efficiency: float  # kWh/m²/yr
    pv_area: float  # m²
    pv_efficiency: float
    natural_ventilation_hours: float  # hours/yr of free cooling


@dataclass
class EnergyBreakdown:
    """エネルギー内訳"""
    heating: float  # kWh/m²/yr
    cooling: float  # kWh/m²/yr
    lighting: float  # kWh/m²/yr
    equipment: float  # kWh/m²/yr
    ventilation: float  # kWh/m²/yr
    hot_water: float  # kWh/m²/yr
    elevator: float  # kWh/m²/yr
    pv_generation: float  # kWh/m²/yr

    @property
    def total_consumption(self) -> float:
        return (self.heating + self.cooling + self.lighting +
                self.equipment + self.ventilation + self.hot_water +
                self.elevator)

    @property
    def net_energy(self) -> float:
        return self.total_consumption - self.pv_generation

    @property
    def zeb_ratio(self) -> float:
        """ZEB達成率（消費に対する創エネ比率）"""
        if self.total_consumption == 0:
            return 100
        return self.pv_generation / self.total_consumption * 100


class ZEBEvaluator:
    """ZEB評価エンジン"""

    # 東京の月別水平面全天日射量 [kWh/m²/month]
    MONTHLY_SOLAR_RADIATION = [
        86.8, 93.0, 111.6, 126.0, 139.5, 114.0,
        127.1, 142.6, 102.0, 89.9, 78.0, 77.5
    ]

    # 基準一次エネルギー消費量 (事務所ビル, 東京, kWh/m²/yr)
    BASELINE_ENERGY = {
        "heating": 15.0,
        "cooling": 45.0,
        "lighting": 30.0,
        "equipment": 40.0,
        "ventilation": 20.0,
        "hot_water": 12.0,
        "elevator": 8.0,
    }

    def calculate_energy(self, design: BuildingDesign) -> EnergyBreakdown:
        """建物設計に基づくエネルギー消費量の算定"""

        # 暖房エネルギー
        envelope_factor = design.wall_u / 0.5  # 基準U値=0.5に対する比
        window_factor = design.window_u / 3.0  # 基準U値=3.0に対する比
        infiltration_factor = design.infiltration_ach / 0.7

        heating = self.BASELINE_ENERGY["heating"] * (
            0.4 * envelope_factor + 0.3 * window_factor +
            0.3 * infiltration_factor
        ) / design.hvac_heating_cop * 3.5  # COP基準=3.5

        # 冷房エネルギー
        solar_gain_factor = (design.window_shgc * design.wwr) / (0.7 * 0.4)
        nat_vent_saving = design.natural_ventilation_hours / 8760 * 0.3

        cooling = self.BASELINE_ENERGY["cooling"] * (
            0.3 * envelope_factor + 0.3 * solar_gain_factor +
            0.2 * window_factor + 0.2
        ) * (1 - nat_vent_saving) / design.hvac_cooling_cop * 3.0

        # 照明エネルギー
        lighting = design.lighting_density / 20.0 * self.BASELINE_ENERGY["lighting"]

        # 機器エネルギー
        equipment = design.equipment_density / 20.0 * self.BASELINE_ENERGY["equipment"]

        # 換気エネルギー
        ventilation = self.BASELINE_ENERGY["ventilation"] * (
            1 - design.heat_recovery_eff * 0.6
        )

        # 給湯エネルギー
        hot_water = self.BASELINE_ENERGY["hot_water"] * (
            1 / design.hot_water_efficiency
        )

        # エレベーターエネルギー
        elevator = design.elevator_efficiency

        # PV発電量
        annual_solar = sum(self.MONTHLY_SOLAR_RADIATION)  # kWh/m²/yr
        pv_system_efficiency = design.pv_efficiency * 0.85  # インバータ損失等
        total_pv_generation = design.pv_area * annual_solar * pv_system_efficiency
        pv_per_area = total_pv_generation / design.floor_area

        return EnergyBreakdown(
            heating=round(heating, 2),
            cooling=round(cooling, 2),
            lighting=round(lighting, 2),
            equipment=round(equipment, 2),
            ventilation=round(ventilation, 2),
            hot_water=round(hot_water, 2),
            elevator=round(elevator, 2),
            pv_generation=round(pv_per_area, 2),
        )

    @staticmethod
    def classify_zeb(energy: EnergyBreakdown, baseline_total: float) -> str:
        """ZEB分類判定"""
        reduction_rate = (1 - energy.total_consumption / baseline_total) * 100
        renewable_coverage = energy.zeb_ratio  # PV / consumption * 100

        if energy.net_energy <= 0:
            return "ZEB"
        elif reduction_rate >= 50 and renewable_coverage >= 75:
            return "Nearly ZEB"
        elif reduction_rate >= 50:
            return "ZEB Ready"
        elif reduction_rate >= 40:
            return "ZEB Oriented"
        else:
            return "Non-ZEB"


def run_zeb_case_study():
    """ZEBケーススタディの実行"""
    evaluator = ZEBEvaluator()

    # ケース1: 基準建物（省エネ基準適合レベル）
    baseline = BuildingDesign(
        name="Baseline Office (H28 Standard)",
        floor_area=5000, num_floors=5,
        wall_u=0.53, roof_u=0.24, window_u=3.49,
        window_shgc=0.70, wwr=0.40,
        infiltration_ach=0.7,
        lighting_density=20.0, equipment_density=20.0,
        hvac_cooling_cop=3.0, hvac_heating_cop=3.0,
        heat_recovery_eff=0.0,
        hot_water_efficiency=0.80,
        elevator_efficiency=8.0,
        pv_area=0, pv_efficiency=0,
        natural_ventilation_hours=0,
    )

    # ケース2: ZEB Ready（高断熱・高効率設備）
    zeb_ready = BuildingDesign(
        name="ZEB Ready Office",
        floor_area=5000, num_floors=5,
        wall_u=0.28, roof_u=0.15, window_u=1.6,
        window_shgc=0.40, wwr=0.40,
        infiltration_ach=0.3,
        lighting_density=8.0, equipment_density=12.0,
        hvac_cooling_cop=5.5, hvac_heating_cop=4.2,
        heat_recovery_eff=0.80,
        hot_water_efficiency=3.5,
        elevator_efficiency=5.0,
        pv_area=0, pv_efficiency=0,
        natural_ventilation_hours=800,
    )

    # ケース3: Nearly ZEB（ZEB Ready + PV）
    nearly_zeb = BuildingDesign(
        name="Nearly ZEB Office",
        floor_area=5000, num_floors=5,
        wall_u=0.28, roof_u=0.15, window_u=1.6,
        window_shgc=0.40, wwr=0.40,
        infiltration_ach=0.3,
        lighting_density=8.0, equipment_density=12.0,
        hvac_cooling_cop=5.5, hvac_heating_cop=4.2,
        heat_recovery_eff=0.80,
        hot_water_efficiency=3.5,
        elevator_efficiency=5.0,
        pv_area=600, pv_efficiency=0.20,
        natural_ventilation_hours=800,
    )

    # ケース4: ZEB（最高性能 + 大容量PV）
    zeb = BuildingDesign(
        name="ZEB Office",
        floor_area=5000, num_floors=5,
        wall_u=0.20, roof_u=0.10, window_u=0.8,
        window_shgc=0.25, wwr=0.35,
        infiltration_ach=0.2,
        lighting_density=6.0, equipment_density=10.0,
        hvac_cooling_cop=7.0, hvac_heating_cop=5.5,
        heat_recovery_eff=0.90,
        hot_water_efficiency=4.0,
        elevator_efficiency=4.0,
        pv_area=1000, pv_efficiency=0.22,
        natural_ventilation_hours=1200,
    )

    cases = [baseline, zeb_ready, nearly_zeb, zeb]
    results = []

    baseline_energy = evaluator.calculate_energy(baseline)
    baseline_total = baseline_energy.total_consumption

    output = {
        "study_type": "ZEB Design Case Study",
        "location": "Tokyo, Japan",
        "building_type": "Office",
        "floor_area_m2": 5000,
        "baseline_energy_kWh_m2_yr": round(baseline_total, 1),
        "cases": [],
    }

    for design in cases:
        energy = evaluator.calculate_energy(design)
        zeb_class = evaluator.classify_zeb(energy, baseline_total)
        reduction = (1 - energy.total_consumption / baseline_total) * 100

        case_result = {
            "name": design.name,
            "envelope": {
                "wall_u": design.wall_u,
                "roof_u": design.roof_u,
                "window_u": design.window_u,
                "window_shgc": design.window_shgc,
                "wwr": design.wwr,
            },
            "systems": {
                "cooling_cop": design.hvac_cooling_cop,
                "heating_cop": design.hvac_heating_cop,
                "heat_recovery": design.heat_recovery_eff,
                "lighting_density_W_m2": design.lighting_density,
            },
            "energy_breakdown_kWh_m2_yr": {
                "heating": energy.heating,
                "cooling": energy.cooling,
                "lighting": energy.lighting,
                "equipment": energy.equipment,
                "ventilation": energy.ventilation,
                "hot_water": energy.hot_water,
                "elevator": energy.elevator,
                "total_consumption": round(energy.total_consumption, 2),
                "pv_generation": energy.pv_generation,
                "net_energy": round(energy.net_energy, 2),
            },
            "evaluation": {
                "zeb_classification": zeb_class,
                "energy_reduction_pct": round(reduction, 1),
                "renewable_ratio_pct": round(energy.zeb_ratio, 1),
            },
        }
        output["cases"].append(case_result)
        results.append((design, energy, zeb_class))

    # 省エネ技術別の効果分析
    tech_effects = []
    tech_configs = [
        ("High-Performance Envelope", {"wall_u": 0.28, "window_u": 1.6}),
        ("LED Lighting", {"lighting_density": 8.0}),
        ("High-COP HVAC", {"hvac_cooling_cop": 5.5, "hvac_heating_cop": 4.2}),
        ("Heat Recovery Ventilation", {"heat_recovery_eff": 0.80}),
        ("Natural Ventilation", {"natural_ventilation_hours": 800}),
        ("HP Hot Water", {"hot_water_efficiency": 3.5}),
    ]

    for tech_name, params in tech_configs:
        # 基準建物に1つだけ技術を適用
        modified = BuildingDesign(
            name=f"Baseline + {tech_name}",
            floor_area=5000, num_floors=5,
            wall_u=params.get("wall_u", baseline.wall_u),
            roof_u=params.get("roof_u", baseline.roof_u),
            window_u=params.get("window_u", baseline.window_u),
            window_shgc=params.get("window_shgc", baseline.window_shgc),
            wwr=baseline.wwr,
            infiltration_ach=params.get("infiltration_ach", baseline.infiltration_ach),
            lighting_density=params.get("lighting_density", baseline.lighting_density),
            equipment_density=baseline.equipment_density,
            hvac_cooling_cop=params.get("hvac_cooling_cop", baseline.hvac_cooling_cop),
            hvac_heating_cop=params.get("hvac_heating_cop", baseline.hvac_heating_cop),
            heat_recovery_eff=params.get("heat_recovery_eff", baseline.heat_recovery_eff),
            hot_water_efficiency=params.get("hot_water_efficiency", baseline.hot_water_efficiency),
            elevator_efficiency=baseline.elevator_efficiency,
            pv_area=0, pv_efficiency=0,
            natural_ventilation_hours=params.get("natural_ventilation_hours", 0),
        )
        mod_energy = evaluator.calculate_energy(modified)
        saving = baseline_total - mod_energy.total_consumption
        saving_pct = saving / baseline_total * 100

        tech_effects.append({
            "technology": tech_name,
            "energy_saving_kWh_m2_yr": round(saving, 2),
            "saving_percentage": round(saving_pct, 1),
        })

    output["technology_effects"] = tech_effects

    with open("results/zeb_case_study_results.json", "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    # サマリー出力
    print("=== ZEBケーススタディ結果 ===")
    print(f"基準エネルギー消費量: {baseline_total:.1f} kWh/m²/yr")
    print()
    for design, energy, zeb_class in results:
        print(f"【{design.name}】")
        print(f"  総消費: {energy.total_consumption:.1f} kWh/m²/yr")
        print(f"  PV発電: {energy.pv_generation:.1f} kWh/m²/yr")
        print(f"  Net: {energy.net_energy:.1f} kWh/m²/yr")
        print(f"  判定: {zeb_class}")
        print()

    print("=== 技術別省エネ効果 ===")
    for te in tech_effects:
        print(f"  {te['technology']}: {te['saving_percentage']:.1f}% "
              f"({te['energy_saving_kWh_m2_yr']:.1f} kWh/m²/yr)")

    return results, output


if __name__ == "__main__":
    run_zeb_case_study()
