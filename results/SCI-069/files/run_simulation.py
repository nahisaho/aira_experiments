"""
Tokyo 2050 UHI Prediction — Main Simulation Runner

Integrates all components:
  1. UCM for 4 Tokyo districts
  2. Anthropogenic heat with 2050 projections
  3. Mitigation scenarios (baseline / moderate / aggressive)
  4. WBGT risk assessment
  5. Visualization
"""

import sys
import os
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.urban_canopy_model import (
    UrbanCanopyModel, BuildingMorphology, TOKYO_MORPHOLOGY,
    SURFACE_MATERIALS, compute_local_climate_zone
)
from src.anthropogenic_heat import (
    AnthropogenicHeatModel, TOKYO_HEAT_PROFILES, project_anthropogenic_heat_2050
)
from src.cooling_effects import MITIGATION_SCENARIOS, CoolingEffectModel
from src.wrf_ucm_coupling import (
    OfflineCouplingEngine, WRFNamelistGenerator, WRFUrbanConfig, TOKYO_DOMAINS
)
from src.wbgt_risk import (
    WBGTCalculator, HeatStrokeRiskAssessor, PopulationExposure
)


def generate_tokyo_summer_forcing(climate_scenario="current"):
    """
    Generate synthetic atmospheric forcing for a typical Tokyo August day.
    Climate scenarios: 'current' (2020s), 'rcp45_2050', 'rcp85_2050'
    """
    hours = np.arange(24)

    # Temperature: sinusoidal with realistic Tokyo summer values
    T_base = 28.0  # daily mean [°C]
    T_amp = 5.0    # diurnal amplitude

    if climate_scenario == "rcp45_2050":
        T_base += 1.8  # CMIP6 median warming for RCP4.5
    elif climate_scenario == "rcp85_2050":
        T_base += 3.2  # CMIP6 median for RCP8.5

    T_air = T_base + T_amp * np.sin(2 * np.pi * (hours - 6) / 24)
    T_air_K = T_air + 273.15

    # Solar radiation
    solar_max = 850  # W/m² peak
    sw_down = np.maximum(0, solar_max * np.sin(np.pi * (hours - 5) / 14))
    sw_down[hours < 5] = 0
    sw_down[hours > 19] = 0

    # Longwave: from air temperature
    sigma = 5.67e-8
    emiss_atm = 0.85
    lw_down = emiss_atm * sigma * T_air_K**4

    # Wind
    u_star = 0.3 + 0.15 * np.sin(2 * np.pi * (hours - 14) / 24)

    # Humidity
    q_air = 0.014 + 0.003 * np.sin(2 * np.pi * (hours + 3) / 24)

    # RH for WBGT calculation
    rh = 65 + 15 * np.cos(2 * np.pi * (hours - 6) / 24)

    # Wind speed at 10m
    wind_10m = 2.5 + 1.0 * np.sin(2 * np.pi * (hours - 14) / 24)

    return {
        "sw_down": sw_down,
        "lw_down": lw_down,
        "T_air": T_air_K,
        "T_air_C": T_air,
        "u_star": u_star,
        "q_air": q_air,
        "rh": rh,
        "wind_10m": wind_10m,
        "hours": hours,
    }


def run_district_simulation(district_name, morphology, heat_profile,
                             forcing, cooling_model=None):
    """Run UCM simulation for one district."""
    ucm = UrbanCanopyModel(morphology)
    engine = OfflineCouplingEngine(ucm, heat_profile, cooling_model)
    results = engine.run_diurnal_cycle(forcing)
    return results


def run_full_simulation():
    """Execute the complete Tokyo UHI simulation."""
    print("=" * 70)
    print("  Tokyo Urban Heat Island Prediction System")
    print("  UHI-Predict v1.0 — Simulation Run")
    print("=" * 70)

    districts = ["marunouchi", "shinjuku", "residential_23ku", "suburban"]
    scenarios = ["current", "rcp45_2050", "rcp85_2050"]

    all_results = {}
    summary_data = []

    # ── 1. District Classification ──
    print("\n[1/6] District Classification (LCZ)")
    print("-" * 40)
    lcz_data = []
    for d in districts:
        morph = TOKYO_MORPHOLOGY[d]
        lcz = compute_local_climate_zone(morph)
        print(f"  {d:20s}: {lcz}  (H/W={morph.hw_ratio:.1f}, SVF={morph.sky_view_factor:.3f})")
        lcz_data.append({
            "district": d, "lcz": lcz,
            "hw_ratio": round(morph.hw_ratio, 2),
            "svf": round(morph.sky_view_factor, 3),
            "building_height": morph.building_height_mean,
            "building_fraction": morph.building_fraction,
        })

    # ── 2. Current Climate Baseline ──
    print("\n[2/6] Current Climate Baseline Simulation")
    print("-" * 40)
    forcing_current = generate_tokyo_summer_forcing("current")

    for d in districts:
        morph = TOKYO_MORPHOLOGY[d]
        heat = TOKYO_HEAT_PROFILES[d]
        results = run_district_simulation(d, morph, heat, forcing_current)
        all_results[f"{d}_current_baseline"] = results

        uhi_max = np.max(results['UHI_intensity'])
        uhi_mean = np.mean(results['UHI_intensity'])
        qf_max = np.max(results['QF'])
        print(f"  {d:20s}: UHI_max={uhi_max:+.2f}°C  UHI_mean={uhi_mean:+.2f}°C  QF_max={qf_max:.1f} W/m²")

        summary_data.append({
            "district": d, "scenario": "current_baseline",
            "UHI_max": round(uhi_max, 2), "UHI_mean": round(uhi_mean, 2),
            "QF_max": round(qf_max, 1),
            "T_canyon_max": round(np.max(results['T_canyon']) - 273.15, 1),
        })

    # ── 3. Mitigation Scenarios ──
    print("\n[3/6] Mitigation Scenario Analysis")
    print("-" * 40)
    mitigation_results = {}

    for scenario_name, cooling_model in MITIGATION_SCENARIOS.items():
        print(f"\n  Scenario: {scenario_name}")
        for d in districts:
            morph = TOKYO_MORPHOLOGY[d]
            heat = TOKYO_HEAT_PROFILES[d]
            results = run_district_simulation(d, morph, heat, forcing_current, cooling_model)
            key = f"{d}_{scenario_name}"
            mitigation_results[key] = results

            uhi_max = np.max(results['UHI_intensity'])
            cooling = np.mean(results['dT_cooling'])
            print(f"    {d:20s}: UHI_max={uhi_max:+.2f}°C  cooling={cooling:+.2f}°C")

            summary_data.append({
                "district": d, "scenario": scenario_name,
                "UHI_max": round(uhi_max, 2),
                "UHI_mean": round(np.mean(results['UHI_intensity']), 2),
                "QF_max": round(np.max(results['QF']), 1),
                "dT_cooling_mean": round(cooling, 2),
            })

    # ── 4. 2050 Projections ──
    print("\n[4/6] 2050 Climate Projections")
    print("-" * 40)
    projection_results = {}

    for climate in ["rcp45_2050", "rcp85_2050"]:
        forcing_future = generate_tokyo_summer_forcing(climate)
        print(f"\n  Climate: {climate} (T_base +{'1.8' if '45' in climate else '3.2'}°C)")

        for d in districts:
            morph = TOKYO_MORPHOLOGY[d]
            heat_2050 = project_anthropogenic_heat_2050(TOKYO_HEAT_PROFILES[d])

            # Baseline (no mitigation)
            results = run_district_simulation(d, morph, heat_2050, forcing_future)
            key = f"{d}_{climate}_baseline"
            projection_results[key] = results
            uhi_max = np.max(results['UHI_intensity'])
            T_max = np.max(results['T_canyon']) - 273.15

            # With aggressive mitigation
            results_mit = run_district_simulation(
                d, morph, heat_2050, forcing_future,
                MITIGATION_SCENARIOS["aggressive_mitigation"])
            key_mit = f"{d}_{climate}_mitigated"
            projection_results[key_mit] = results_mit
            uhi_mit = np.max(results_mit['UHI_intensity'])
            T_mit = np.max(results_mit['T_canyon']) - 273.15

            print(f"    {d:20s}: T_max={T_max:.1f}°C (baseline) → {T_mit:.1f}°C (mitigated)  "
                  f"UHI={uhi_max:+.1f} → {uhi_mit:+.1f}°C")

            summary_data.append({
                "district": d, "scenario": f"{climate}_baseline",
                "UHI_max": round(uhi_max, 2),
                "T_canyon_max": round(T_max, 1),
            })
            summary_data.append({
                "district": d, "scenario": f"{climate}_mitigated",
                "UHI_max": round(uhi_mit, 2),
                "T_canyon_max": round(T_mit, 1),
            })

    # ── 5. WBGT Heat Stroke Risk ──
    print("\n[5/6] Heat Stroke Risk Assessment (WBGT)")
    print("-" * 40)

    wbgt_calc = WBGTCalculator()
    ward_populations = {
        "marunouchi": PopulationExposure(total_population=67_000, outdoor_fraction=0.25,
                                          elderly_fraction=0.22, outdoor_worker_fraction=0.08),
        "shinjuku": PopulationExposure(total_population=350_000, outdoor_fraction=0.20,
                                        elderly_fraction=0.22, outdoor_worker_fraction=0.04),
        "residential_23ku": PopulationExposure(total_population=740_000, outdoor_fraction=0.12,
                                                elderly_fraction=0.25, outdoor_worker_fraction=0.03),
        "suburban": PopulationExposure(total_population=920_000, outdoor_fraction=0.10,
                                        elderly_fraction=0.28, outdoor_worker_fraction=0.02),
    }

    wbgt_results = {}
    for climate_label, forcing in [("current", forcing_current),
                                     ("rcp85_2050", generate_tokyo_summer_forcing("rcp85_2050"))]:
        print(f"\n  Climate: {climate_label}")
        for d in districts:
            T_C = forcing['T_air_C']
            # Add UHI effect from simulation
            key = f"{d}_{climate_label}_baseline" if climate_label != "current" else f"{d}_current_baseline"
            if key in all_results:
                uhi = all_results[key]['UHI_intensity']
            elif key in projection_results:
                uhi = projection_results[key]['UHI_intensity']
            else:
                uhi = np.zeros(24)

            T_effective = T_C + uhi
            assessor = HeatStrokeRiskAssessor(ward_populations[d])
            risk = assessor.daily_risk_profile(
                T_effective, forcing['rh'], forcing['wind_10m'], forcing['sw_down'])

            wbgt_results[f"{d}_{climate_label}"] = risk
            print(f"    {d:20s}: WBGT_peak={risk['peak_wbgt']:.1f}°C ({risk['hourly_risk_level'][risk['peak_hour']]}) "
                  f"at {risk['peak_hour']:02d}:00  "
                  f"danger_hrs={risk['danger_hours']}  "
                  f"patients={risk['daily_total_patients']:.0f}/day")

    # ── 6. WRF Configuration ──
    print("\n[6/6] WRF-UCM Configuration Generated")
    print("-" * 40)
    namelist_gen = WRFNamelistGenerator()
    config = WRFUrbanConfig()
    physics_nml = namelist_gen.generate_physics_namelist(config, TOKYO_DOMAINS)
    domains_nml = namelist_gen.generate_domains_namelist(TOKYO_DOMAINS)

    with open("results/wrf_namelist_physics.txt", "w") as f:
        f.write(physics_nml)
    with open("results/wrf_namelist_domains.txt", "w") as f:
        f.write(domains_nml)
    print("  → results/wrf_namelist_physics.txt")
    print("  → results/wrf_namelist_domains.txt")

    for dom in TOKYO_DOMAINS:
        print(f"    {dom.name}: {dom.dx:.0f}m resolution, {dom.nx}×{dom.ny} grid, dt={dom.dt:.0f}s")

    # ── Save Results ──
    print("\n" + "=" * 70)
    print("  Saving Results")
    print("=" * 70)

    # Summary CSV
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv("results/simulation_summary.csv", index=False)
    print("  → results/simulation_summary.csv")

    # LCZ classification
    df_lcz = pd.DataFrame(lcz_data)
    df_lcz.to_csv("results/lcz_classification.csv", index=False)
    print("  → results/lcz_classification.csv")

    # Detailed hourly data for key scenarios
    for key_name in ["marunouchi_current_baseline", "shinjuku_current_baseline"]:
        if key_name in all_results:
            df = pd.DataFrame(all_results[key_name])
            df['hour'] = np.arange(24)
            df.to_csv(f"results/{key_name}_hourly.csv", index=False)
            print(f"  → results/{key_name}_hourly.csv")

    # WBGT results
    wbgt_summary = []
    for key, risk in wbgt_results.items():
        wbgt_summary.append({
            "scenario": key,
            "peak_wbgt": round(risk['peak_wbgt'], 1),
            "peak_hour": risk['peak_hour'],
            "danger_hours": risk['danger_hours'],
            "daily_patients": round(risk['daily_total_patients'], 0),
        })
    df_wbgt = pd.DataFrame(wbgt_summary)
    df_wbgt.to_csv("results/wbgt_risk_summary.csv", index=False)
    print("  → results/wbgt_risk_summary.csv")

    # Return all data for visualization
    return {
        "all_results": all_results,
        "mitigation_results": mitigation_results,
        "projection_results": projection_results,
        "wbgt_results": wbgt_results,
        "forcing_current": forcing_current,
        "summary": summary_data,
        "lcz": lcz_data,
    }


if __name__ == "__main__":
    data = run_full_simulation()
    print("\n✓ Simulation complete. Running visualization...")

    # Run visualization
    exec(open("run_visualization.py").read())
