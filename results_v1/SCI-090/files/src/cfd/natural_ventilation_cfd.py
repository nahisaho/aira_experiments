"""
自然換気CFD解析モジュール
クロスベンチレーション性能評価・風環境シミュレーション

機能:
- 風向・風速条件別の自然換気量算定
- クロスベンチレーション有効性評価
- 室内気流分布の簡易計算
- 通風性能指標の算出
"""

import numpy as np
import json
from dataclasses import dataclass, field


@dataclass
class Opening:
    """開口部定義"""
    name: str
    area: float  # m²
    azimuth: float  # degrees from north
    height: float  # 開口中心高さ m
    discharge_coeff: float = 0.6
    opening_factor: float = 0.5  # 開口割合
    zone: str = ""


@dataclass
class WindCondition:
    """外部風条件"""
    direction: float  # degrees from north
    speed: float  # m/s
    temperature: float  # °C
    profile_exponent: float = 0.25  # 地表面粗度係数


@dataclass
class CFDResult:
    """CFD解析結果"""
    wind_direction: float
    wind_speed: float
    ventilation_rate: float  # m³/s
    ach: float  # air changes per hour
    avg_indoor_velocity: float  # m/s
    max_indoor_velocity: float  # m/s
    pressure_diff: float  # Pa
    effectiveness: float  # 換気効率 (0-1)
    comfort_zone_ratio: float  # 快適域割合


class NaturalVentilationCFD:
    """簡易CFDベースの自然換気解析"""

    AIR_DENSITY = 1.2  # kg/m³

    # 風圧係数（方位別、矩形建物）
    WIND_PRESSURE_COEFFICIENTS = {
        0: {"windward": 0.7, "leeward": -0.4, "side": -0.6},
        45: {"windward": 0.4, "leeward": -0.3, "side": -0.5},
        90: {"windward": 0.7, "leeward": -0.4, "side": -0.6},
    }

    def __init__(self, openings: list[Opening], building_volume: float,
                 building_height: float = 17.5):
        self.openings = openings
        self.building_volume = building_volume
        self.building_height = building_height

    def _get_pressure_coefficient(self, surface_azimuth: float,
                                  wind_direction: float) -> float:
        """面の風圧係数を計算"""
        relative_angle = abs(surface_azimuth - wind_direction) % 360
        if relative_angle > 180:
            relative_angle = 360 - relative_angle

        if relative_angle < 45:
            return 0.7 - 0.3 * (relative_angle / 45)
        elif relative_angle < 135:
            return -0.6
        else:
            return -0.4 + 0.1 * ((relative_angle - 135) / 45)

    def _wind_velocity_at_height(self, v_ref: float, height: float,
                                 ref_height: float = 10.0,
                                 alpha: float = 0.25) -> float:
        """べき乗則による高さ補正"""
        return v_ref * (height / ref_height) ** alpha

    def calculate_ventilation(self, wind: WindCondition) -> CFDResult:
        """風力換気量の計算（オリフィス方程式ベース）"""

        # 各開口部の風圧係数と動圧
        opening_data = []
        for op in self.openings:
            cp = self._get_pressure_coefficient(op.azimuth, wind.direction)
            v_local = self._wind_velocity_at_height(
                wind.speed, op.height, alpha=wind.profile_exponent
            )
            dynamic_pressure = 0.5 * self.AIR_DENSITY * v_local ** 2
            wind_pressure = cp * dynamic_pressure

            effective_area = op.area * op.opening_factor * op.discharge_coeff
            opening_data.append({
                "opening": op,
                "cp": cp,
                "pressure": wind_pressure,
                "effective_area": effective_area,
            })

        # 風上と風下の開口を分類
        windward = [d for d in opening_data if d["cp"] > 0]
        leeward = [d for d in opening_data if d["cp"] <= 0]

        if not windward or not leeward:
            return CFDResult(
                wind_direction=wind.direction,
                wind_speed=wind.speed,
                ventilation_rate=0,
                ach=0,
                avg_indoor_velocity=0,
                max_indoor_velocity=0,
                pressure_diff=0,
                effectiveness=0,
                comfort_zone_ratio=0,
            )

        # 平均圧力差
        avg_p_windward = np.mean([d["pressure"] for d in windward])
        avg_p_leeward = np.mean([d["pressure"] for d in leeward])
        delta_p = avg_p_windward - avg_p_leeward

        # 等価開口面積（直列抵抗モデル）
        a_windward = sum(d["effective_area"] for d in windward)
        a_leeward = sum(d["effective_area"] for d in leeward)

        if a_windward + a_leeward == 0:
            a_eq = 0
        else:
            a_eq = (a_windward * a_leeward) / np.sqrt(a_windward**2 + a_leeward**2)

        # 換気量
        if delta_p > 0:
            q = a_eq * np.sqrt(2 * delta_p / self.AIR_DENSITY)
        else:
            q = 0

        ach = q * 3600 / self.building_volume if self.building_volume > 0 else 0

        # 室内気流速度推定
        total_floor_area = self.building_volume / self.building_height * 5  # 概算
        cross_section = total_floor_area / 5  # 1フロア断面
        avg_velocity = q / max(cross_section * 0.3, 0.01)  # 有効断面
        max_velocity = avg_velocity * 2.5  # 最大速度はジェット効果

        # 換気効率
        if ach > 0:
            # クロスベンチレーションの効率は対面開口配置で高い
            cross_vent_ratio = min(a_windward, a_leeward) / max(a_windward, a_leeward)
            effectiveness = 0.5 + 0.3 * cross_vent_ratio
        else:
            effectiveness = 0

        # 快適域割合（風速0.3-1.5m/sの範囲）
        if avg_velocity < 0.1:
            comfort_ratio = 0.2
        elif avg_velocity < 0.3:
            comfort_ratio = 0.5
        elif avg_velocity <= 1.5:
            comfort_ratio = 0.85
        else:
            comfort_ratio = max(0.3, 0.85 - (avg_velocity - 1.5) * 0.2)

        return CFDResult(
            wind_direction=wind.direction,
            wind_speed=wind.speed,
            ventilation_rate=round(q, 3),
            ach=round(ach, 2),
            avg_indoor_velocity=round(avg_velocity, 3),
            max_indoor_velocity=round(max_velocity, 3),
            pressure_diff=round(delta_p, 2),
            effectiveness=round(effectiveness, 3),
            comfort_zone_ratio=round(comfort_ratio, 3),
        )

    def parametric_study(self, wind_speeds: list[float],
                         wind_directions: list[float],
                         temperature: float = 28.0) -> list[CFDResult]:
        """パラメトリック解析：風向・風速を変えて換気性能を評価"""
        results = []
        for speed in wind_speeds:
            for direction in wind_directions:
                wind = WindCondition(
                    direction=direction,
                    speed=speed,
                    temperature=temperature,
                )
                result = self.calculate_ventilation(wind)
                results.append(result)
        return results


def create_building_openings() -> list[Opening]:
    """ZEBオフィスの開口部定義"""
    openings = []
    floor_height = 3.5

    for floor in range(1, 6):
        h = (floor - 0.5) * floor_height
        for azimuth, label in [(0, "N"), (90, "E"), (180, "S"), (270, "W")]:
            openings.append(Opening(
                name=f"Window_{label}_F{floor}",
                area=44.3,  # m² (WWR 40%)
                azimuth=azimuth,
                height=h,
                discharge_coeff=0.6,
                opening_factor=0.5,
                zone=f"Floor_{floor}_Zone",
            ))

    return openings


def run_cfd_analysis():
    """自然換気CFD解析の実行"""
    openings = create_building_openings()
    building_volume = 17500  # m³

    cfd = NaturalVentilationCFD(openings, building_volume)

    # パラメトリック解析
    wind_speeds = [1.0, 2.0, 3.0, 4.0, 5.0]
    wind_directions = [0, 45, 90, 135, 180, 225, 270, 315]

    results = cfd.parametric_study(wind_speeds, wind_directions)

    # 結果整理
    output = {
        "analysis_type": "Natural Ventilation CFD Analysis",
        "building_volume_m3": building_volume,
        "num_openings": len(openings),
        "parametric_results": [],
        "summary": {},
    }

    ach_matrix = np.zeros((len(wind_speeds), len(wind_directions)))

    for i, result in enumerate(results):
        speed_idx = i // len(wind_directions)
        dir_idx = i % len(wind_directions)
        ach_matrix[speed_idx, dir_idx] = result.ach

        output["parametric_results"].append({
            "wind_speed_m_s": result.wind_speed,
            "wind_direction_deg": result.wind_direction,
            "ventilation_rate_m3_s": result.ventilation_rate,
            "ach": result.ach,
            "avg_velocity_m_s": result.avg_indoor_velocity,
            "pressure_diff_Pa": result.pressure_diff,
            "effectiveness": result.effectiveness,
            "comfort_ratio": result.comfort_zone_ratio,
        })

    # サマリー統計
    all_ach = [r.ach for r in results]
    all_eff = [r.effectiveness for r in results]
    all_comfort = [r.comfort_zone_ratio for r in results]

    output["summary"] = {
        "avg_ach": round(np.mean(all_ach), 2),
        "max_ach": round(np.max(all_ach), 2),
        "min_ach": round(np.min(all_ach), 2),
        "avg_effectiveness": round(np.mean(all_eff), 3),
        "avg_comfort_ratio": round(np.mean(all_comfort), 3),
        "optimal_wind_direction_deg": wind_directions[
            np.argmax(np.mean(ach_matrix, axis=0))
        ],
        "cross_ventilation_viable": bool(np.mean(all_ach) > 3.0),
        "ach_matrix": ach_matrix.tolist(),
        "wind_speeds": wind_speeds,
        "wind_directions": wind_directions,
    }

    with open("results/cfd_ventilation_results.json", "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print("=== 自然換気CFD解析結果 ===")
    print(f"平均換気回数: {output['summary']['avg_ach']} ACH")
    print(f"最大換気回数: {output['summary']['max_ach']} ACH")
    print(f"平均換気効率: {output['summary']['avg_effectiveness']}")
    print(f"平均快適域割合: {output['summary']['avg_comfort_ratio']}")
    print(f"最適風向: {output['summary']['optimal_wind_direction_deg']}°")
    print(f"クロスベンチレーション可否: {output['summary']['cross_ventilation_viable']}")

    return results, output


if __name__ == "__main__":
    run_cfd_analysis()
