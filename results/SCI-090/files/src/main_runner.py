#!/usr/bin/env python3
"""
BIM統合環境性能シミュレーションシステム - メインランナー
全モジュールを順次実行し、統合結果を生成
"""

import sys
import os
import json
from datetime import datetime

src_dir = os.path.dirname(os.path.abspath(__file__))
workspace_dir = os.path.dirname(src_dir)
sys.path.insert(0, src_dir)
os.chdir(workspace_dir)


def main():
    timestamp = datetime.now().isoformat()
    print("=" * 60)
    print("BIM-Integrated Environmental Performance Simulation System")
    print(f"Execution: {timestamp}")
    print("=" * 60)

    log_entries = []

    def log(phase, event, **kwargs):
        entry = {
            "timestamp": datetime.now().isoformat(),
            "phase": phase,
            "event_type": event,
            "actor": "co-scientist",
            **kwargs,
        }
        log_entries.append(entry)

    log("init", "run_started", status="ok")

    # Phase 1: IFC変換
    print("\n[Phase 1] IFC Model Conversion")
    print("-" * 40)
    from ifc_converter.ifc_to_simulation import run_conversion_pipeline
    model, conv_summary = run_conversion_pipeline()
    log("ifc_conversion", "skill_selected",
        skill_or_tool="ifc_converter",
        files_written=["results/building_model.idf", "results/building_model.rad",
                       "results/cfd_config.json", "results/conversion_summary.json"],
        status="ok")

    # Phase 2: 熱負荷シミュレーション
    print("\n[Phase 2] Thermal Load Simulation")
    print("-" * 40)
    from thermal.thermal_load_simulation import run_thermal_simulation
    thermal_result, thermal_output = run_thermal_simulation()
    log("thermal_simulation", "skill_selected",
        skill_or_tool="thermal_load_simulator",
        files_written=["results/thermal_simulation_results.json"],
        status="ok")

    # Phase 3: 自然換気CFD解析
    print("\n[Phase 3] Natural Ventilation CFD Analysis")
    print("-" * 40)
    from cfd.natural_ventilation_cfd import run_cfd_analysis
    cfd_results, cfd_output = run_cfd_analysis()
    log("cfd_analysis", "skill_selected",
        skill_or_tool="natural_ventilation_cfd",
        files_written=["results/cfd_ventilation_results.json"],
        status="ok")

    # Phase 4: 昼光シミュレーション
    print("\n[Phase 4] Daylight Simulation")
    print("-" * 40)
    from daylight.daylight_simulation import run_daylight_simulation
    daylight_results, daylight_output = run_daylight_simulation()
    log("daylight_simulation", "skill_selected",
        skill_or_tool="daylight_simulator",
        files_written=["results/daylight_simulation_results.json",
                       "results/illuminance_grid_sample.csv"],
        status="ok")

    # Phase 5: ZEBケーススタディ
    print("\n[Phase 5] ZEB Case Study")
    print("-" * 40)
    from zeb_case.zeb_case_study import run_zeb_case_study
    zeb_results, zeb_output = run_zeb_case_study()
    log("zeb_case_study", "skill_selected",
        skill_or_tool="zeb_evaluator",
        files_written=["results/zeb_case_study_results.json"],
        status="ok")

    # Phase 6: 可視化・ダッシュボード生成
    print("\n[Phase 6] Dashboard & Visualization")
    print("-" * 40)
    from dashboard.visualization import generate_all_figures
    generate_all_figures()
    log("visualization", "skill_selected",
        skill_or_tool="dashboard_visualization",
        files_written=[
            "figures/fig1_monthly_energy.png",
            "figures/fig2_ventilation_heatmap.png",
            "figures/fig3_daylight_performance.png",
            "figures/fig4_zeb_comparison.png",
            "figures/fig5_technology_waterfall.png",
            "figures/fig6_integrated_dashboard.png",
        ],
        status="ok")

    log("finalize", "run_completed", status="ok")

    # ログ保存
    os.makedirs("logs", exist_ok=True)
    with open("logs/process-log.jsonl", "w") as f:
        for entry in log_entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print("\n" + "=" * 60)
    print("All simulations completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
