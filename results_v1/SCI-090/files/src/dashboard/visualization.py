"""
統合ダッシュボード・可視化モジュール
構造・設備・環境シミュレーション結果の統合表示

機能:
- 月別エネルギー消費・発電バランス
- 自然換気ACHヒートマップ
- 昼光率分布マップ
- ZEB達成度比較チャート
- 省エネ技術効果ウォーターフォール
"""

import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.dpi'] = 150


def load_results():
    """全シミュレーション結果の読み込み"""
    with open("results/thermal_simulation_results.json") as f:
        thermal = json.load(f)
    with open("results/cfd_ventilation_results.json") as f:
        cfd = json.load(f)
    with open("results/daylight_simulation_results.json") as f:
        daylight = json.load(f)
    with open("results/zeb_case_study_results.json") as f:
        zeb = json.load(f)
    return thermal, cfd, daylight, zeb


def plot_monthly_energy(thermal_data: dict, save_path: str):
    """月別エネルギー消費プロファイル"""
    months = list(thermal_data["results"]["monthly_data"].keys())
    data = thermal_data["results"]["monthly_data"]

    heating = [data[m]["heating_kWh"] for m in months]
    cooling = [data[m]["cooling_kWh"] for m in months]
    lighting = [data[m]["lighting_kWh"] for m in months]
    equipment = [data[m]["equipment_kWh"] for m in months]
    fan = [data[m]["fan_kWh"] for m in months]

    x = np.arange(len(months))
    width = 0.65

    fig, ax = plt.subplots(figsize=(12, 6))

    colors = ['#d62728', '#2196F3', '#FFC107', '#4CAF50', '#9C27B0']
    labels = ['Heating', 'Cooling', 'Lighting', 'Equipment', 'Fan/Pump']
    datasets = [heating, cooling, lighting, equipment, fan]

    bottom = np.zeros(len(months))
    for i, (d, c, l) in enumerate(zip(datasets, colors, labels)):
        ax.bar(x, d, width, bottom=bottom, color=c, label=l, alpha=0.85)
        bottom += np.array(d)

    ax.set_xlabel('Month')
    ax.set_ylabel('Energy Consumption (kWh)')
    ax.set_title('Monthly Energy Consumption Profile - ZEB Office (5,000 m²)')
    ax.set_xticks(x)
    ax.set_xticklabels(months)
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.3)

    # 注釈
    peak = thermal_data["results"]["peak_cooling_kW"]
    annual = thermal_data["results"]["annual_primary_energy_kWh_m2"]
    ax.text(0.02, 0.95, f'Peak Cooling: {peak} kW\nAnnual: {annual} kWh/m²/yr',
            transform=ax.transAxes, fontsize=9,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_ventilation_heatmap(cfd_data: dict, save_path: str):
    """風向×風速による換気回数ヒートマップ"""
    summary = cfd_data["summary"]
    ach_matrix = np.array(summary["ach_matrix"])
    speeds = summary["wind_speeds"]
    directions = summary["wind_directions"]

    fig, ax = plt.subplots(figsize=(10, 6))

    dir_labels = [f"{d}°" for d in directions]

    im = ax.imshow(ach_matrix, cmap='YlGnBu', aspect='auto',
                   interpolation='nearest')

    ax.set_xticks(range(len(directions)))
    ax.set_xticklabels(dir_labels)
    ax.set_yticks(range(len(speeds)))
    ax.set_yticklabels([f"{s} m/s" for s in speeds])
    ax.set_xlabel('Wind Direction')
    ax.set_ylabel('Wind Speed')
    ax.set_title('Natural Ventilation Performance (ACH) - Wind Direction vs Speed')

    # 値をセルに表示
    for i in range(len(speeds)):
        for j in range(len(directions)):
            val = ach_matrix[i, j]
            color = 'white' if val > np.max(ach_matrix) * 0.6 else 'black'
            ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                    color=color, fontsize=9)

    cbar = plt.colorbar(im, ax=ax, label='ACH (Air Changes per Hour)')

    # 3ACH基準線
    ax.text(0.02, 0.02, f'Cross-ventilation viable: ACH ≥ 3.0\n'
            f'Avg ACH: {summary["avg_ach"]}\n'
            f'Optimal direction: {summary["optimal_wind_direction_deg"]}°',
            transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_daylight_performance(daylight_data: dict, save_path: str):
    """昼光性能比較チャート"""
    rooms = daylight_data["rooms"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # sDA分布
    sda_values = [r["sDA_300_50"] for r in rooms]
    ase_values = [r["ASE_1000_250"] for r in rooms]
    df_values = [r["daylight_factor_avg"] for r in rooms]

    # 方位別集計
    orientations = ["North", "East", "South", "West"]
    orient_sda = {o: [] for o in orientations}
    orient_ase = {o: [] for o in orientations}
    orient_df = {o: [] for o in orientations}

    for r in rooms:
        for o in orientations:
            if o in r["name"]:
                orient_sda[o].append(r["sDA_300_50"])
                orient_ase[o].append(r["ASE_1000_250"])
                orient_df[o].append(r["daylight_factor_avg"])

    x = np.arange(len(orientations))
    width = 0.35

    # sDA by orientation
    ax = axes[0]
    sda_means = [np.mean(orient_sda[o]) for o in orientations]
    bars = ax.bar(x, sda_means, width, color=['#1565C0', '#43A047', '#E65100', '#6A1B9A'],
                  alpha=0.8)
    ax.axhline(y=55, color='red', linestyle='--', label='LEED threshold (55%)')
    ax.set_xlabel('Orientation')
    ax.set_ylabel('sDA300/50% (%)')
    ax.set_title('Spatial Daylight Autonomy by Orientation')
    ax.set_xticks(x)
    ax.set_xticklabels(orientations)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 100)

    for bar, val in zip(bars, sda_means):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.0f}%', ha='center', va='bottom', fontsize=9)

    # ASE by orientation
    ax = axes[1]
    ase_means = [np.mean(orient_ase[o]) for o in orientations]
    bars = ax.bar(x, ase_means, width, color=['#1565C0', '#43A047', '#E65100', '#6A1B9A'],
                  alpha=0.8)
    ax.axhline(y=10, color='red', linestyle='--', label='LEED limit (10%)')
    ax.set_xlabel('Orientation')
    ax.set_ylabel('ASE1000/250h (%)')
    ax.set_title('Annual Sunlight Exposure by Orientation')
    ax.set_xticks(x)
    ax.set_xticklabels(orientations)
    ax.legend(fontsize=8)

    for bar, val in zip(bars, ase_means):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{val:.0f}%', ha='center', va='bottom', fontsize=9)

    # Daylight Factor distribution
    ax = axes[2]
    df_means = [np.mean(orient_df[o]) for o in orientations]
    bars = ax.bar(x, df_means, width, color=['#1565C0', '#43A047', '#E65100', '#6A1B9A'],
                  alpha=0.8)
    ax.axhline(y=2.0, color='green', linestyle='--', label='Target DF (2%)')
    ax.set_xlabel('Orientation')
    ax.set_ylabel('Avg Daylight Factor (%)')
    ax.set_title('Average Daylight Factor by Orientation')
    ax.set_xticks(x)
    ax.set_xticklabels(orientations)
    ax.legend(fontsize=8)

    for bar, val in zip(bars, df_means):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=9)

    plt.suptitle('Daylight Performance Analysis', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_zeb_comparison(zeb_data: dict, save_path: str):
    """ZEBケース比較チャート"""
    cases = zeb_data["cases"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # エネルギー消費内訳の積み上げ棒グラフ
    ax = axes[0]
    categories = ['heating', 'cooling', 'lighting', 'equipment',
                   'ventilation', 'hot_water', 'elevator']
    cat_labels = ['Heating', 'Cooling', 'Lighting', 'Equipment',
                  'Ventilation', 'Hot Water', 'Elevator']
    colors = ['#d62728', '#2196F3', '#FFC107', '#4CAF50',
              '#9C27B0', '#FF9800', '#607D8B']

    case_names = [c["name"].replace(" Office", "").replace(" (H28 Standard)", "")
                  for c in cases]
    x = np.arange(len(cases))
    width = 0.6

    bottom = np.zeros(len(cases))
    for cat, color, label in zip(categories, colors, cat_labels):
        vals = [c["energy_breakdown_kWh_m2_yr"][cat] for c in cases]
        ax.bar(x, vals, width, bottom=bottom, color=color, label=label, alpha=0.85)
        bottom += np.array(vals)

    # PV発電を負の値で表示
    pv_vals = [-c["energy_breakdown_kWh_m2_yr"]["pv_generation"] for c in cases]
    ax.bar(x, pv_vals, width, color='#FFD700', edgecolor='#B8860B',
           linewidth=1.5, label='PV Generation', alpha=0.9, hatch='///')

    ax.set_xlabel('Design Case')
    ax.set_ylabel('Energy (kWh/m²/yr)')
    ax.set_title('Energy Balance Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(case_names, rotation=15, ha='right', fontsize=9)
    ax.legend(loc='upper right', fontsize=7, ncol=2)
    ax.grid(axis='y', alpha=0.3)
    ax.axhline(y=0, color='black', linewidth=0.8)

    # ZEB達成度レーダーチャート風の棒グラフ
    ax = axes[1]
    metrics = ['energy_reduction_pct', 'renewable_ratio_pct']
    metric_labels = ['Energy Reduction (%)', 'Renewable Ratio (%)']

    x2 = np.arange(len(cases))
    width2 = 0.35

    reductions = [c["evaluation"]["energy_reduction_pct"] for c in cases]
    renewables = [c["evaluation"]["renewable_ratio_pct"] for c in cases]

    bars1 = ax.bar(x2 - width2/2, reductions, width2, color='#1976D2',
                   label='Energy Reduction', alpha=0.85)
    bars2 = ax.bar(x2 + width2/2, renewables, width2, color='#388E3C',
                   label='Renewable Ratio', alpha=0.85)

    # ZEB分類表示
    for i, case in enumerate(cases):
        zeb_class = case["evaluation"]["zeb_classification"]
        color_map = {
            "Non-ZEB": "#F44336", "ZEB Oriented": "#FF9800",
            "ZEB Ready": "#FFC107", "Nearly ZEB": "#8BC34A", "ZEB": "#4CAF50"
        }
        ax.text(i, max(reductions[i], renewables[i]) + 3,
                zeb_class, ha='center', fontsize=9, fontweight='bold',
                color=color_map.get(zeb_class, 'black'),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor=color_map.get(zeb_class, 'gray'), alpha=0.8))

    ax.set_xlabel('Design Case')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('ZEB Achievement Evaluation')
    ax.set_xticks(x2)
    ax.set_xticklabels(case_names, rotation=15, ha='right', fontsize=9)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # ZEB基準線
    ax.axhline(y=100, color='green', linestyle=':', alpha=0.5, label='ZEB (100%)')
    ax.axhline(y=75, color='#8BC34A', linestyle='--', alpha=0.4)
    ax.axhline(y=50, color='#FFC107', linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_technology_waterfall(zeb_data: dict, save_path: str):
    """省エネ技術効果ウォーターフォールチャート"""
    baseline = zeb_data["baseline_energy_kWh_m2_yr"]
    techs = zeb_data["technology_effects"]

    fig, ax = plt.subplots(figsize=(12, 6))

    labels = ["Baseline"] + [t["technology"] for t in techs] + ["Combined\n(ZEB Ready)"]
    savings = [0] + [t["energy_saving_kWh_m2_yr"] for t in techs]

    # ZEB Ready の消費量を取得
    zeb_ready_case = zeb_data["cases"][1]
    zeb_ready_total = zeb_ready_case["energy_breakdown_kWh_m2_yr"]["total_consumption"]
    combined_saving = baseline - zeb_ready_total
    savings.append(combined_saving)

    # ウォーターフォール構築
    running_total = [baseline]
    for s in savings[1:-1]:
        running_total.append(running_total[-1] - s)
    running_total.append(zeb_ready_total)

    colors_wf = ['#1565C0']  # baseline
    for s in savings[1:-1]:
        colors_wf.append('#4CAF50' if s > 0 else '#F44336')
    colors_wf.append('#E65100')  # combined

    # 底部（floating bars）
    bottoms = [0]
    for i in range(1, len(running_total) - 1):
        bottoms.append(running_total[i])
    bottoms.append(0)

    heights = [baseline]
    for s in savings[1:-1]:
        heights.append(s)
    heights.append(zeb_ready_total)

    x = np.arange(len(labels))
    bars = ax.bar(x, heights, bottom=bottoms, color=colors_wf, alpha=0.85,
                  edgecolor='gray', linewidth=0.5)

    # 接続線
    for i in range(len(labels) - 2):
        ax.plot([i + 0.4, i + 0.6], [running_total[i+1], running_total[i+1]],
                color='gray', linewidth=0.8, linestyle='-')

    # 値ラベル
    for i, (bar, h) in enumerate(zip(bars, heights)):
        if i == 0:
            text = f'{baseline:.0f}'
            y_pos = baseline / 2
        elif i == len(labels) - 1:
            text = f'{zeb_ready_total:.0f}'
            y_pos = zeb_ready_total / 2
        else:
            text = f'-{h:.1f}'
            y_pos = bottoms[i] + h / 2

        ax.text(bar.get_x() + bar.get_width()/2, y_pos,
                text, ha='center', va='center', fontsize=8, fontweight='bold')

    ax.set_xlabel('Energy Conservation Measures')
    ax.set_ylabel('Energy (kWh/m²/yr)')
    ax.set_title('Energy Conservation Technology Waterfall - Individual Effects')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    # 注釈
    reduction_pct = (1 - zeb_ready_total / baseline) * 100
    ax.text(0.98, 0.95, f'Total Reduction: {reduction_pct:.0f}%\n'
            f'Baseline: {baseline:.0f} kWh/m²/yr\n'
            f'ZEB Ready: {zeb_ready_total:.0f} kWh/m²/yr',
            transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_integrated_dashboard(thermal, cfd, daylight, zeb, save_path: str):
    """統合ダッシュボード"""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # Panel 1: 月別エネルギー
    ax1 = fig.add_subplot(gs[0, 0])
    months = list(thermal["results"]["monthly_data"].keys())
    data = thermal["results"]["monthly_data"]
    heating = [data[m]["heating_kWh"] / 1000 for m in months]
    cooling = [data[m]["cooling_kWh"] / 1000 for m in months]
    x = np.arange(len(months))
    ax1.bar(x - 0.2, heating, 0.35, color='#d62728', label='Heating', alpha=0.8)
    ax1.bar(x + 0.2, cooling, 0.35, color='#2196F3', label='Cooling', alpha=0.8)
    ax1.set_title('Monthly Thermal Load (MWh)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(months, fontsize=7)
    ax1.legend(fontsize=7)
    ax1.grid(axis='y', alpha=0.3)

    # Panel 2: 換気ヒートマップ
    ax2 = fig.add_subplot(gs[0, 1])
    ach_matrix = np.array(cfd["summary"]["ach_matrix"])
    im = ax2.imshow(ach_matrix, cmap='YlGnBu', aspect='auto')
    ax2.set_title('Ventilation ACH (Direction × Speed)')
    ax2.set_xlabel('Wind Direction')
    ax2.set_ylabel('Wind Speed')
    ax2.set_xticks(range(len(cfd["summary"]["wind_directions"])))
    ax2.set_xticklabels([f'{d}°' for d in cfd["summary"]["wind_directions"]], fontsize=7)
    ax2.set_yticks(range(len(cfd["summary"]["wind_speeds"])))
    ax2.set_yticklabels([f'{s}m/s' for s in cfd["summary"]["wind_speeds"]], fontsize=7)
    plt.colorbar(im, ax=ax2, shrink=0.8)

    # Panel 3: 昼光sDA
    ax3 = fig.add_subplot(gs[0, 2])
    orientations = ["North", "East", "South", "West"]
    orient_sda = {o: [] for o in orientations}
    for r in daylight["rooms"]:
        for o in orientations:
            if o in r["name"]:
                orient_sda[o].append(r["sDA_300_50"])
    sda_means = [np.mean(orient_sda[o]) for o in orientations]
    bars = ax3.bar(orientations, sda_means,
                   color=['#1565C0', '#43A047', '#E65100', '#6A1B9A'], alpha=0.8)
    ax3.axhline(y=55, color='red', linestyle='--', alpha=0.7, label='LEED 55%')
    ax3.set_title('Daylight Autonomy (sDA) by Orientation')
    ax3.set_ylabel('sDA (%)')
    ax3.set_ylim(0, 100)
    ax3.legend(fontsize=7)

    # Panel 4: ZEB比較
    ax4 = fig.add_subplot(gs[1, 0])
    cases = zeb["cases"]
    case_names = [c["name"].split(" Office")[0].split("(")[0].strip() for c in cases]
    consumption = [c["energy_breakdown_kWh_m2_yr"]["total_consumption"] for c in cases]
    pv = [c["energy_breakdown_kWh_m2_yr"]["pv_generation"] for c in cases]
    net = [c["energy_breakdown_kWh_m2_yr"]["net_energy"] for c in cases]

    x4 = np.arange(len(cases))
    ax4.bar(x4, consumption, 0.35, color='#F44336', alpha=0.8, label='Consumption')
    ax4.bar(x4, [-p for p in pv], 0.35, color='#4CAF50', alpha=0.8, label='PV Gen')
    ax4.plot(x4, net, 'ko-', markersize=6, label='Net Energy')
    ax4.axhline(y=0, color='black', linewidth=0.8)
    ax4.set_title('ZEB Energy Balance')
    ax4.set_xticks(x4)
    ax4.set_xticklabels(case_names, fontsize=7, rotation=15)
    ax4.set_ylabel('kWh/m²/yr')
    ax4.legend(fontsize=7)

    # Panel 5: 省エネ技術効果
    ax5 = fig.add_subplot(gs[1, 1])
    techs = zeb["technology_effects"]
    tech_names = [t["technology"].replace("High-Performance ", "HP\n")
                  .replace("Heat Recovery ", "HR\n")
                  .replace("Natural ", "Nat.\n")
                  for t in techs]
    savings_pct = [t["saving_percentage"] for t in techs]
    colors_bar = plt.cm.viridis(np.linspace(0.3, 0.9, len(techs)))
    ax5.barh(range(len(techs)), savings_pct, color=colors_bar, alpha=0.85)
    ax5.set_yticks(range(len(techs)))
    ax5.set_yticklabels(tech_names, fontsize=7)
    ax5.set_xlabel('Energy Saving (%)')
    ax5.set_title('Technology Impact Analysis')
    for i, v in enumerate(savings_pct):
        ax5.text(v + 0.3, i, f'{v:.1f}%', va='center', fontsize=8)

    # Panel 6: KPI Summary
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')

    kpis = [
        ("Building", "ZEB Office Demo (5,000 m²)"),
        ("Location", "Tokyo, Japan (35.68°N)"),
        ("", ""),
        ("Peak Cooling", f"{thermal['results']['peak_cooling_kW']} kW"),
        ("Annual Energy", f"{thermal['results']['annual_primary_energy_kWh_m2']} kWh/m²/yr"),
        ("ZEB Score", f"{thermal['results']['zeb_score_percent']}%"),
        ("", ""),
        ("Avg ACH", f"{cfd['summary']['avg_ach']}"),
        ("Cross-Vent", f"{'Viable' if cfd['summary']['cross_ventilation_viable'] else 'Limited'}"),
        ("", ""),
        ("Avg sDA", f"{daylight['summary']['avg_sDA']}%"),
        ("LEED Rate", f"{daylight['summary']['leed_compliance_rate']}%"),
        ("Lighting Save", f"{daylight['summary']['avg_lighting_saving_pct']}%"),
        ("", ""),
        ("Best ZEB", f"{cases[-1]['evaluation']['zeb_classification']}"),
    ]

    y_pos = 0.95
    for label, value in kpis:
        if label == "":
            y_pos -= 0.03
            continue
        ax6.text(0.05, y_pos, label + ":", fontsize=9, fontweight='bold',
                transform=ax6.transAxes, verticalalignment='top')
        ax6.text(0.55, y_pos, value, fontsize=9,
                transform=ax6.transAxes, verticalalignment='top')
        y_pos -= 0.065

    ax6.set_title('Key Performance Indicators', fontsize=11)
    ax6.add_patch(mpatches.FancyBboxPatch(
        (0.02, 0.02), 0.96, 0.96, transform=ax6.transAxes,
        boxstyle="round,pad=0.02", facecolor='#F5F5F5',
        edgecolor='gray', linewidth=1))

    fig.suptitle('BIM-Integrated Environmental Performance Simulation Dashboard',
                 fontsize=16, fontweight='bold', y=0.98)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def generate_all_figures():
    """全図表の生成"""
    print("=== 図表生成開始 ===")
    thermal, cfd, daylight, zeb = load_results()

    plot_monthly_energy(thermal, "figures/fig1_monthly_energy.png")
    plot_ventilation_heatmap(cfd, "figures/fig2_ventilation_heatmap.png")
    plot_daylight_performance(daylight, "figures/fig3_daylight_performance.png")
    plot_zeb_comparison(zeb, "figures/fig4_zeb_comparison.png")
    plot_technology_waterfall(zeb, "figures/fig5_technology_waterfall.png")
    plot_integrated_dashboard(thermal, cfd, daylight, zeb,
                              "figures/fig6_integrated_dashboard.png")

    print("=== 全図表生成完了 ===")


if __name__ == "__main__":
    generate_all_figures()
