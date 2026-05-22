"""
昼光シミュレーションモジュール
Radiance/Honeybeeベースの昼光環境評価

機能:
- 空間昼光自律性 (sDA) 計算
- 年間日射グレア確率 (ASE) 計算
- 昼光率 (DF) 分布
- LEED v4 Daylight Credit 適合判定
- 照明制御による省エネ効果推定
"""

import numpy as np
import json
from dataclasses import dataclass, field


@dataclass
class RoomGeometry:
    """室空間ジオメトリ"""
    name: str
    width: float  # m
    depth: float  # m
    height: float  # m
    window_width: float  # m
    window_height: float  # m
    sill_height: float  # m
    orientation: float  # degrees from north
    glazing_vlt: float = 0.40  # 可視光透過率
    wall_reflectance: float = 0.50
    floor_reflectance: float = 0.20
    ceiling_reflectance: float = 0.80
    shading_type: str = "none"  # none, overhang, louver, blind


@dataclass
class DaylightResult:
    """昼光シミュレーション結果"""
    room_name: str
    daylight_factor_avg: float  # %
    daylight_factor_min: float  # %
    daylight_factor_uniformity: float  # min/avg
    sda_300_50: float  # Spatial Daylight Autonomy (% area ≥300lux for ≥50% hours)
    ase_1000_250: float  # Annual Sunlight Exposure (% area ≥1000lux for ≥250 hours)
    leed_compliant: bool
    lighting_energy_saving: float  # % reduction
    illuminance_grid: list[list[float]]  # lux grid


class DaylightSimulator:
    """簡易昼光シミュレーター（Radiance代理モデル）"""

    # 月別CIE曇天空照度 [klux] at Tokyo
    MONTHLY_DIFFUSE_ILLUMINANCE = [
        12.5, 15.0, 18.0, 22.0, 25.0, 20.0,
        22.0, 26.0, 18.0, 15.0, 13.0, 11.5,
    ]

    # 月別直達日射照度 [klux] at Tokyo (南面垂直面)
    MONTHLY_DIRECT_ILLUMINANCE = [
        35.0, 40.0, 45.0, 50.0, 55.0, 40.0,
        48.0, 58.0, 38.0, 32.0, 28.0, 30.0,
    ]

    def __init__(self, grid_resolution: float = 0.5):
        self.grid_resolution = grid_resolution

    def simulate_room(self, room: RoomGeometry) -> DaylightResult:
        """室内昼光環境のシミュレーション"""
        # グリッド作成
        nx = int(room.width / self.grid_resolution)
        ny = int(room.depth / self.grid_resolution)
        nx = max(nx, 2)
        ny = max(ny, 2)

        # BRS daylight factor formula (simplified)
        window_area = room.window_width * room.window_height
        floor_area = room.width * room.depth
        total_surface_area = 2 * (
            room.width * room.depth
            + room.width * room.height
            + room.depth * room.height
        )

        # 平均反射率
        avg_reflectance = (
            room.wall_reflectance * 2 * (room.width + room.depth) * room.height
            + room.floor_reflectance * floor_area
            + room.ceiling_reflectance * floor_area
        ) / total_surface_area

        # 昼光率グリッド計算
        df_grid = np.zeros((ny, nx))
        illuminance_grid = np.zeros((ny, nx))

        for iy in range(ny):
            for ix in range(nx):
                x = (ix + 0.5) * self.grid_resolution
                y = (iy + 0.5) * self.grid_resolution

                # 窓からの距離に基づく減衰
                dist_from_window = y  # 窓面からの奥行き
                window_center_x = room.width / 2
                lateral_offset = abs(x - window_center_x)

                # 立体角因子（simplified view factor）
                # 窓の見かけの立体角
                h_angle = np.arctan2(room.window_width / 2 - lateral_offset,
                                     max(dist_from_window, 0.1))
                h_angle = max(h_angle, 0.05)

                v_top = room.sill_height + room.window_height - 0.8  # work plane at 0.8m
                v_bot = max(room.sill_height - 0.8, 0)
                v_angle = np.arctan2(v_top, max(dist_from_window, 0.1)) - \
                          np.arctan2(v_bot, max(dist_from_window, 0.1))
                v_angle = max(v_angle, 0.02)

                # Sky component
                sky_component = (h_angle * v_angle) / (2 * np.pi) * 100
                sky_component *= room.glazing_vlt

                # Internally reflected component
                irc = window_area * room.glazing_vlt * avg_reflectance / \
                      (total_surface_area * (1 - avg_reflectance)) * 100 * 0.85

                # Externally reflected component (簡易)
                erc = 0.5

                # シェーディング補正
                shading_factor = 1.0
                if room.shading_type == "overhang":
                    shading_factor = 0.75
                elif room.shading_type == "louver":
                    shading_factor = 0.60
                elif room.shading_type == "blind":
                    shading_factor = 0.50

                df = (sky_component + irc + erc) * shading_factor
                df = max(df, 0.1)

                df_grid[iy, ix] = df

                # 年間平均照度（曇天空ベース）
                avg_sky_illuminance = np.mean(self.MONTHLY_DIFFUSE_ILLUMINANCE) * 1000
                illuminance_grid[iy, ix] = df / 100 * avg_sky_illuminance

        # 統計量算出
        df_avg = float(np.mean(df_grid))
        df_min = float(np.min(df_grid))
        df_uniformity = df_min / df_avg if df_avg > 0 else 0

        # sDA計算（300lux以上を50%以上の時間達成するグリッドポイントの割合）
        annual_hours = 3650  # 年間稼働時間概算
        points_meeting_sda = 0
        total_points = nx * ny

        for iy in range(ny):
            for ix in range(nx):
                hours_above_300 = 0
                for month in range(12):
                    sky_ill = self.MONTHLY_DIFFUSE_ILLUMINANCE[month] * 1000
                    direct_ill = self.MONTHLY_DIRECT_ILLUMINANCE[month] * 1000
                    total_ill = (df_grid[iy, ix] / 100 * sky_ill +
                                df_grid[iy, ix] / 100 * direct_ill * 0.3)

                    monthly_hours = annual_hours / 12
                    if total_ill >= 300:
                        hours_above_300 += monthly_hours

                if hours_above_300 / annual_hours >= 0.50:
                    points_meeting_sda += 1

        sda = (points_meeting_sda / total_points) * 100

        # ASE計算（1000lux以上を年間250時間以上受けるポイントの割合）
        points_exceeding_ase = 0
        for iy in range(ny):
            for ix in range(nx):
                hours_above_1000 = 0
                dist_from_win = (iy + 0.5) * self.grid_resolution
                for month in range(12):
                    direct_ill = self.MONTHLY_DIRECT_ILLUMINANCE[month] * 1000
                    # 直達日射は窓近傍のみ影響、奥行に応じて急減
                    depth_factor = max(0, 1.0 - dist_from_win / (room.depth * 0.5))
                    point_ill = df_grid[iy, ix] / 100 * direct_ill * depth_factor
                    if point_ill >= 1000:
                        hours_above_1000 += 30 * 4  # 直射の影響時間は限定的

                if hours_above_1000 >= 250:
                    points_exceeding_ase += 1

        ase = (points_exceeding_ase / total_points) * 100

        # LEED v4 判定
        leed_compliant = sda >= 55.0 and ase <= 10.0

        # 照明省エネ効果推定
        # sDAが高い = 昼光利用可能時間が長い = 照明削減可能
        if sda >= 70:
            lighting_saving = 45.0
        elif sda >= 55:
            lighting_saving = 35.0
        elif sda >= 40:
            lighting_saving = 25.0
        else:
            lighting_saving = 15.0

        return DaylightResult(
            room_name=room.name,
            daylight_factor_avg=round(df_avg, 2),
            daylight_factor_min=round(df_min, 2),
            daylight_factor_uniformity=round(df_uniformity, 3),
            sda_300_50=round(sda, 1),
            ase_1000_250=round(ase, 1),
            leed_compliant=leed_compliant,
            lighting_energy_saving=lighting_saving,
            illuminance_grid=illuminance_grid.tolist(),
        )


def run_daylight_simulation():
    """各フロア・方位の昼光シミュレーション実行"""
    simulator = DaylightSimulator(grid_resolution=1.0)

    rooms = []
    orientations = [
        (0, "North"), (90, "East"), (180, "South"), (270, "West")
    ]

    for floor in range(1, 6):
        for az, label in orientations:
            room = RoomGeometry(
                name=f"Floor{floor}_{label}_Zone",
                width=15.0,
                depth=10.0,
                height=3.5,
                window_width=12.0,
                window_height=2.1,
                sill_height=0.9,
                orientation=az,
                glazing_vlt=0.40,
                wall_reflectance=0.50,
                floor_reflectance=0.20,
                ceiling_reflectance=0.80,
                shading_type="overhang" if az == 180 else "none",
            )
            rooms.append(room)

    results = []
    for room in rooms:
        result = simulator.simulate_room(room)
        results.append(result)

    # 結果集計
    output = {
        "analysis_type": "Daylight Simulation (Radiance/Honeybee proxy)",
        "grid_resolution_m": 1.0,
        "rooms": [],
        "summary": {},
    }

    sda_values = []
    ase_values = []
    df_values = []
    savings = []

    for r in results:
        room_data = {
            "name": r.room_name,
            "daylight_factor_avg": r.daylight_factor_avg,
            "daylight_factor_min": r.daylight_factor_min,
            "uniformity": r.daylight_factor_uniformity,
            "sDA_300_50": r.sda_300_50,
            "ASE_1000_250": r.ase_1000_250,
            "leed_compliant": r.leed_compliant,
            "lighting_energy_saving_pct": r.lighting_energy_saving,
        }
        output["rooms"].append(room_data)
        sda_values.append(r.sda_300_50)
        ase_values.append(r.ase_1000_250)
        df_values.append(r.daylight_factor_avg)
        savings.append(r.lighting_energy_saving)

    output["summary"] = {
        "avg_sDA": round(np.mean(sda_values), 1),
        "avg_ASE": round(np.mean(ase_values), 1),
        "avg_daylight_factor": round(np.mean(df_values), 2),
        "leed_compliance_rate": round(
            sum(1 for r in results if r.leed_compliant) / len(results) * 100, 1
        ),
        "avg_lighting_saving_pct": round(np.mean(savings), 1),
        "total_rooms_analyzed": len(results),
    }

    with open("results/daylight_simulation_results.json", "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    # 代表室の照度グリッドを保存
    rep_result = results[4]  # Floor2_North
    np.savetxt("results/illuminance_grid_sample.csv",
               np.array(rep_result.illuminance_grid),
               delimiter=",", fmt="%.1f",
               header="Illuminance grid (lux) - " + rep_result.room_name)

    print("=== 昼光シミュレーション結果 ===")
    print(f"解析室数: {len(results)}")
    print(f"平均sDA(300/50%): {output['summary']['avg_sDA']}%")
    print(f"平均ASE(1000/250): {output['summary']['avg_ASE']}%")
    print(f"平均昼光率: {output['summary']['avg_daylight_factor']}%")
    print(f"LEED適合率: {output['summary']['leed_compliance_rate']}%")
    print(f"照明省エネ効果: {output['summary']['avg_lighting_saving_pct']}%")

    return results, output


if __name__ == "__main__":
    run_daylight_simulation()
