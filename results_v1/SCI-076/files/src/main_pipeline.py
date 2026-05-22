"""
Main pipeline runner: マルチモーダル作物生育予測・収量推定システム
Executes all modules and collects results.
"""
import sys
import json
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))

LOGS_DIR = Path(__file__).parent.parent / "logs"
LOGS_DIR.mkdir(exist_ok=True)


def log_event(phase, event_type, skill="main_pipeline", **kwargs):
    entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        **kwargs
    }
    with open(LOGS_DIR / "process-log.jsonl", "a") as f:
        f.write(json.dumps(entry, default=str) + "\n")


def main():
    log_event("init", "run_started")
    
    print("=" * 70)
    print("  Multimodal Crop Growth Prediction & Yield Estimation System")
    print("  ケーススタディ: 新潟県 水稲栽培")
    print("=" * 70)
    
    # Module 1: Vegetation Indices
    print("\n[1/6] Computing vegetation indices...")
    log_event("module1", "skill_selected", skill="vegetation_indices")
    from vegetation_indices import run_vegetation_analysis
    vi_stats, doy, ndvi_spatial = run_vegetation_analysis()
    log_event("module1", "handoff_completed",
              files_written=["fig01_vegetation_indices_timeseries.png",
                             "fig02_spatial_vegetation_map.png",
                             "vegetation_indices.npz"])
    
    # Module 2: Weather + Crop Model
    print("\n[2/6] Running crop growth model...")
    log_event("module2", "skill_selected", skill="weather_crop_model")
    from weather_crop_model import run_crop_model_analysis
    weather_df, model_results, stage_dates = run_crop_model_analysis()
    log_event("module2", "handoff_completed",
              files_written=["fig03_weather_data.png", "fig04_crop_model_results.png",
                             "weather_niigata_2025.csv", "crop_model_output.csv"])
    
    # Module 3: Soil Interpolation
    print("\n[3/6] Performing soil data interpolation...")
    log_event("module3", "skill_selected", skill="soil_interpolation")
    from soil_interpolation import run_soil_interpolation
    sensor_df, soil_results = run_soil_interpolation()
    log_event("module3", "handoff_completed",
              files_written=["fig05_soil_interpolation.png", "fig06_variograms.png",
                             "soil_sensor_data.csv", "soil_interpolation.npz"])
    
    # Module 4: Deep Learning Yield Prediction
    print("\n[4/6] Training CNN-LSTM yield model...")
    log_event("module4", "skill_selected", skill="deep_learning_yield")
    from deep_learning_yield import run_deep_learning_yield
    dl_metrics, yield_map = run_deep_learning_yield()
    log_event("module4", "handoff_completed",
              files_written=["fig07_dl_model_performance.png", "fig08_yield_map.png",
                             "dl_model_metrics.json", "predicted_yield_map.npy"])
    
    # Module 5: Variable Rate Fertilization
    print("\n[5/6] Generating VRA prescription maps...")
    log_event("module5", "skill_selected", skill="variable_rate_fertilization")
    from variable_rate_fertilization import run_vra_analysis
    vra_results = run_vra_analysis()
    log_event("module5", "handoff_completed",
              files_written=["fig09_vra_prescription.png", "fig10_economic_analysis.png",
                             "vra_maps.npz"])
    
    # Module 6: GEE/GeoPandas Pipeline
    print("\n[6/6] Building GEE/GeoPandas pipeline...")
    log_event("module6", "skill_selected", skill="gee_pipeline")
    from gee_pipeline import run_geopandas_analysis
    gdf, field_summary = run_geopandas_analysis()
    log_event("module6", "handoff_completed",
              files_written=["fig11_pipeline_architecture.png", "fig12_field_analysis.png",
                             "rice_fields.geojson", "field_summary.json",
                             "gee_pipeline_template.py"])
    
    log_event("final", "run_completed", status="ok")
    
    print("\n" + "=" * 70)
    print("  All modules completed successfully!")
    print("=" * 70)
    
    return {
        'vi_stats': {k: {'peak': float(v['mean'].max())} for k, v in vi_stats.items()},
        'crop_model': {
            'final_yield': float(model_results['yield_tha'].iloc[-1]),
            'peak_lai': float(model_results['lai'].max()),
            'stage_dates': {k: v.strftime('%Y-%m-%d') for k, v in stage_dates.items()},
        },
        'dl_metrics': dl_metrics,
        'vra_results': vra_results,
        'field_summary': field_summary,
    }


if __name__ == "__main__":
    results = main()
    
    with open(Path(__file__).parent.parent / "results" / "all_results_summary.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    
    print("\nAll results saved to results/all_results_summary.json")
