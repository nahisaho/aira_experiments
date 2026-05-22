"""Run the end-to-end synthetic InSAR analysis workflow."""

from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from insar_pipeline import (
    classify_alert_level,
    detect_acceleration,
    detect_sse,
    detect_strain_anomaly,
    era5_correction,
    estimate_aps,
    estimate_coupling,
    estimate_deformation,
    estimate_ps_coherence,
    estimate_velocity,
    extract_transient,
    fit_linear_trend,
    fit_seasonal,
    gacos_correction,
    generate_alert,
    integrate_gps,
    kalman_filter_decompose,
    los_to_3d,
    monitor_stress_accumulation,
    select_pairs,
    select_psc,
    spatial_clustering,
    spatial_filter,
    unwrap_phase,
)
from insar_pipeline.displacement_3d import propagate_errors

DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
LOG_PATH = ROOT / "logs" / "process-log.jsonl"



def log_event(phase: str, event_type: str, files_written: list[str], status: str = "ok", skill_or_tool: str = "run_analysis") -> None:
    record = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "files_written": files_written,
        "status": status,
    }
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")



def linear_slope(cube: np.ndarray, time_years: np.ndarray) -> np.ndarray:
    t = time_years - np.mean(time_years)
    denom = np.sum(t**2)
    centered = cube.reshape(cube.shape[0], -1)
    slope = (t[:, None] * centered).sum(axis=0) / denom
    return slope.reshape(cube.shape[1:])



def fdr_bh(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty_like(ranked)
    n = len(p)
    prev = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        value = min(prev, ranked[idx] * n / rank)
        adjusted[idx] = value
        prev = value
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result.tolist()



def compare_groups(name: str, inside: np.ndarray, outside: np.ndarray) -> dict[str, float | str]:
    inside = np.asarray(inside, dtype=float)
    outside = np.asarray(outside, dtype=float)
    inside = inside[np.isfinite(inside)]
    outside = outside[np.isfinite(outside)]
    sample_inside = inside if inside.size <= 500 else inside[:: max(inside.size // 500, 1)]
    sample_outside = outside if outside.size <= 500 else outside[:: max(outside.size // 500, 1)]
    normal_inside = stats.normaltest(sample_inside).pvalue > 0.05 if sample_inside.size >= 8 else False
    normal_outside = stats.normaltest(sample_outside).pvalue > 0.05 if sample_outside.size >= 8 else False
    variance_ratio = np.var(sample_inside, ddof=1) / max(np.var(sample_outside, ddof=1), 1.0e-9)
    if normal_inside and normal_outside and 0.25 <= variance_ratio <= 4.0:
        test_name = "Welch t-test"
        pvalue = float(stats.ttest_ind(inside, outside, equal_var=False).pvalue)
        assumptions = "approx_normal"
    else:
        test_name = "Mann-Whitney U"
        pvalue = float(stats.mannwhitneyu(inside, outside, alternative="two-sided").pvalue)
        assumptions = "nonparametric"
    mean_inside = float(np.mean(inside))
    mean_outside = float(np.mean(outside))
    diff = mean_inside - mean_outside
    pooled = np.sqrt(((inside.size - 1) * np.var(inside, ddof=1) + (outside.size - 1) * np.var(outside, ddof=1)) / max(inside.size + outside.size - 2, 1))
    effect_size = float(diff / max(pooled, 1.0e-9))
    se = np.sqrt(np.var(inside, ddof=1) / inside.size + np.var(outside, ddof=1) / outside.size)
    ci = [float(diff - 1.96 * se), float(diff + 1.96 * se)]
    return {
        "name": name,
        "test": test_name,
        "assumptions": assumptions,
        "pvalue": pvalue,
        "effect_size_cohens_d": effect_size,
        "mean_inside": mean_inside,
        "mean_outside": mean_outside,
        "difference": float(diff),
        "ci95_low": ci[0],
        "ci95_high": ci[1],
        "n_inside": int(inside.size),
        "n_outside": int(outside.size),
    }



def write_csv(path: Path, header: list[str], rows: list[list[float | int | str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)



def apply_atmospheric_pipeline(cube: np.ndarray, elevation: np.ndarray, incidence: float) -> tuple[np.ndarray, dict[str, object]]:
    era_corrected, era_delay = era5_correction(cube, elevation, incidence)
    gacos_corrected, gacos_delay = gacos_correction(era_corrected, elevation)
    filtered = spatial_filter(gacos_corrected)
    aps, variogram = estimate_aps(filtered)
    corrected = filtered - aps
    meta = {
        "era_delay_mean": float(np.mean(era_delay)),
        "gacos_delay_mean": float(np.mean(gacos_delay)),
        "aps_range_pixels": float(variogram["range_pixels"]),
    }
    return corrected, meta



def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    log_event("execute", "run_started", [])

    amplitude_stack = np.load(DATA_DIR / "amplitude_stack.npy")
    elevation = np.load(DATA_DIR / "elevation.npy")
    time_days = np.load(DATA_DIR / "times_days.npy")
    time_years = np.load(DATA_DIR / "times_years.npy")
    perp_baselines = np.load(DATA_DIR / "perp_baselines.npy")
    wrapped_phase_asc = np.load(DATA_DIR / "wrapped_phase_asc.npy")
    wrapped_phase_desc = np.load(DATA_DIR / "wrapped_phase_desc.npy")
    los_asc_observed = np.load(DATA_DIR / "los_asc_observed.npy")
    los_desc_observed = np.load(DATA_DIR / "los_desc_observed.npy")
    true_east = np.load(DATA_DIR / "true_east.npy")
    true_north = np.load(DATA_DIR / "true_north.npy")
    true_up = np.load(DATA_DIR / "true_up.npy")
    gps_points = np.load(DATA_DIR / "gps_points.npy")
    sse_mask = np.load(DATA_DIR / "sse_mask.npy").astype(bool)

    rows, cols = elevation.shape
    wavelength = 0.056
    asc_incidence = np.deg2rad(34.0)
    desc_incidence = np.deg2rad(37.0)
    asc_azimuth = np.deg2rad(100.0)
    desc_azimuth = np.deg2rad(280.0)

    psc_mask, da = select_psc(amplitude_stack)
    coherence = estimate_ps_coherence(wrapped_phase_asc, psc_mask)
    unwrapped_asc = unwrap_phase(wrapped_phase_asc)
    ps_results = estimate_velocity(unwrapped_asc, time_years, wavelength=wavelength, dem_sensitivity=perp_baselines / np.std(perp_baselines))
    ps_velocity = ps_results["velocity"]
    ps_velocity_masked = np.where(psc_mask, ps_velocity, np.nan)

    corrected_asc, aps_meta_asc = apply_atmospheric_pipeline(los_asc_observed, elevation, asc_incidence)
    corrected_desc, aps_meta_desc = apply_atmospheric_pipeline(los_desc_observed, elevation, desc_incidence)

    pairs = select_pairs(time_days, perp_baselines)
    pair_stack = np.stack([corrected_asc[j] - corrected_asc[i] for i, j in pairs])
    sbas = estimate_deformation(pair_stack.reshape(len(pairs), -1), pairs, len(time_days))
    sbas_deformation = sbas["deformation"].reshape(len(time_days), rows, cols)
    sbas_velocity = linear_slope(sbas_deformation, time_years)

    representative_series = np.mean(corrected_asc[:, sse_mask], axis=1)
    linear = fit_linear_trend(time_years, representative_series)
    seasonal = fit_seasonal(time_years, representative_series - linear["fitted"])
    robust = extract_transient(time_years, representative_series)
    kalman = kalman_filter_decompose(time_years, representative_series)

    strain_series = np.gradient(np.mean(corrected_asc[:, sse_mask], axis=1), time_years)
    cusum = detect_strain_anomaly(strain_series)
    acceleration = detect_acceleration(representative_series, time_years)
    transient_map = np.mean(corrected_asc[26:34], axis=0) - np.mean(corrected_asc[18:26], axis=0)
    anomaly_mask = np.abs(transient_map) > (np.mean(np.abs(transient_map)) + 2.0 * np.std(transient_map))
    clusters = spatial_clustering(anomaly_mask)
    precursor_alert = classify_alert_level(cusum, acceleration, clusters)

    north_prior = np.full((rows, cols), np.mean(gps_points[:, 3]))
    east, north, up = los_to_3d(corrected_asc[-1], corrected_desc[-1], asc_incidence, desc_incidence, asc_azimuth, desc_azimuth, north_prior=north_prior, north_sigma=0.02)
    (east_ref, north_ref, up_ref), gps_stats = integrate_gps((east, north, up), gps_points)
    enu_errors = propagate_errors(0.003, np.array([asc_incidence, desc_incidence]), np.array([asc_azimuth, desc_azimuth]), north_sigma=0.02)

    coupling_map = estimate_coupling(np.nan_to_num(sbas_velocity))
    stress = monitor_stress_accumulation(coupling_map)
    sse_events = detect_sse(robust["transient"])
    nankai_alert = generate_alert(coupling_map, sse_events, precursor_alert["level"])

    true_up_velocity = linear_slope(true_up, time_years)
    rmse_velocity = float(np.sqrt(np.nanmean((sbas_velocity - true_up_velocity) ** 2)))
    rmse_up_final = float(np.sqrt(np.nanmean((up_ref - true_up[-1]) ** 2)))
    rmse_east_final = float(np.sqrt(np.nanmean((east_ref - true_east[-1]) ** 2)))
    rmse_north_final = float(np.sqrt(np.nanmean((north_ref - true_north[-1]) ** 2)))

    comparisons = [
        compare_groups("velocity_ps", ps_velocity_masked[sse_mask], ps_velocity_masked[~sse_mask]),
        compare_groups("transient_map", transient_map[sse_mask], transient_map[~sse_mask]),
        compare_groups("coupling_map", coupling_map[sse_mask], coupling_map[~sse_mask]),
    ]
    adjusted = fdr_bh([item["pvalue"] for item in comparisons])
    for item, adj in zip(comparisons, adjusted):
        item["fdr_pvalue"] = float(adj)

    summary = {
        "grid_shape": [rows, cols],
        "n_times": int(len(time_days)),
        "psc_count": int(np.sum(psc_mask)),
        "psc_fraction": float(np.mean(psc_mask)),
        "mean_temporal_coherence": float(np.nanmean(coherence)),
        "pair_count": int(len(pairs)),
        "network_rank": int(sbas["inversion"]["rank"]),
        "aps_range_pixels_asc": aps_meta_asc["aps_range_pixels"],
        "aps_range_pixels_desc": aps_meta_desc["aps_range_pixels"],
        "linear_velocity_m_per_yr": float(linear["slope"]),
        "linear_velocity_ci95": [float(linear["slope_ci"][0]), float(linear["slope_ci"][1])],
        "annual_amplitude_m": float(seasonal["annual_amplitude"]),
        "semi_annual_amplitude_m": float(seasonal["semi_annual_amplitude"]),
        "transient_peak_m": float(np.max(np.abs(robust["transient"]))),
        "precursor_alert": precursor_alert,
        "nankai_alert": nankai_alert,
        "sse_events": sse_events,
        "mean_coupling": float(np.mean(coupling_map)),
        "max_stress": float(stress["summary"]["max_stress"]),
        "rmse_velocity_m_per_yr": rmse_velocity,
        "rmse_east_final_m": rmse_east_final,
        "rmse_north_final_m": rmse_north_final,
        "rmse_up_final_m": rmse_up_final,
        "gps_rms_m": float(gps_stats["gps_rms"]),
        "enu_std_m": [float(value) for value in enu_errors["std_enu"]],
    }

    files_written: list[str] = []
    np.save(RESULTS_DIR / "psc_mask.npy", psc_mask.astype(np.uint8))
    np.save(RESULTS_DIR / "amplitude_dispersion.npy", da)
    np.save(RESULTS_DIR / "coherence.npy", np.nan_to_num(coherence, nan=0.0))
    np.save(RESULTS_DIR / "ps_velocity.npy", np.nan_to_num(ps_velocity_masked, nan=0.0))
    np.save(RESULTS_DIR / "sbas_velocity.npy", sbas_velocity)
    np.save(RESULTS_DIR / "transient_map.npy", transient_map)
    np.save(RESULTS_DIR / "coupling_map.npy", coupling_map)
    np.save(RESULTS_DIR / "stress_map.npy", stress["stress"])
    np.save(RESULTS_DIR / "east_final.npy", east_ref)
    np.save(RESULTS_DIR / "north_final.npy", north_ref)
    np.save(RESULTS_DIR / "up_final.npy", up_ref)
    np.save(RESULTS_DIR / "cluster_labels.npy", clusters["labels"])
    files_written.extend([
        "results/psc_mask.npy",
        "results/amplitude_dispersion.npy",
        "results/coherence.npy",
        "results/ps_velocity.npy",
        "results/sbas_velocity.npy",
        "results/transient_map.npy",
        "results/coupling_map.npy",
        "results/stress_map.npy",
        "results/east_final.npy",
        "results/north_final.npy",
        "results/up_final.npy",
        "results/cluster_labels.npy",
    ])

    velocity_rows = []
    for row in range(rows):
        for col in range(cols):
            velocity_rows.append([
                row,
                col,
                int(psc_mask[row, col]),
                float(da[row, col]),
                float(np.nan_to_num(coherence[row, col], nan=0.0)),
                float(np.nan_to_num(ps_velocity_masked[row, col], nan=0.0)),
                float(sbas_velocity[row, col]),
                float(coupling_map[row, col]),
                float(stress["stress"][row, col]),
                float(east_ref[row, col]),
                float(north_ref[row, col]),
                float(up_ref[row, col]),
            ])
    write_csv(
        RESULTS_DIR / "velocity_field.csv",
        ["row", "col", "psc", "amplitude_dispersion", "coherence", "ps_velocity_m_per_yr", "sbas_velocity_m_per_yr", "coupling", "stress", "east_final_m", "north_final_m", "up_final_m"],
        velocity_rows,
    )
    files_written.append("results/velocity_field.csv")

    ts_rows = []
    cumulative_alert = np.full(len(time_days), precursor_alert["score"])
    for idx in range(len(time_days)):
        ts_rows.append([
            idx,
            float(time_days[idx]),
            float(time_years[idx]),
            float(representative_series[idx]),
            float(linear["fitted"][idx]),
            float(seasonal["seasonal"][idx]),
            float(robust["transient"][idx]),
            float(kalman["level"][idx]),
            float(kalman["slope"][idx]),
            float(strain_series[idx]),
            float(cusum["cusum_positive"][idx]),
            float(cusum["cusum_negative"][idx]),
            float(acceleration["acceleration"][idx]),
            float(acceleration["zscore"][idx]),
            float(cumulative_alert[idx]),
        ])
    write_csv(
        RESULTS_DIR / "time_series_decomposition.csv",
        ["time_index", "days", "years", "observed_m", "linear_m", "seasonal_m", "transient_m", "kalman_level_m", "kalman_slope_m_per_yr", "strain_rate_m_per_yr", "cusum_positive", "cusum_negative", "acceleration_m_per_yr2", "acceleration_zscore", "alert_score"],
        ts_rows,
    )
    files_written.append("results/time_series_decomposition.csv")

    with (RESULTS_DIR / "analysis_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    files_written.append("results/analysis_summary.json")

    with (RESULTS_DIR / "precursor_alert.json").open("w", encoding="utf-8") as handle:
        json.dump({"precursor_alert": precursor_alert, "nankai_alert": nankai_alert, "sse_events": sse_events, "clusters": {"n_clusters": clusters["n_clusters"], "cluster_sizes": clusters["cluster_sizes"]}}, handle, indent=2, ensure_ascii=False)
    files_written.append("results/precursor_alert.json")

    with (RESULTS_DIR / "statistics.json").open("w", encoding="utf-8") as handle:
        json.dump(comparisons, handle, indent=2, ensure_ascii=False)
    files_written.append("results/statistics.json")

    stats_md = ["# Statistical Summary", "", "Assumption checks were applied before significance testing.", ""]
    for item in comparisons:
        stats_md.extend(
            [
                f"## {item['name']}",
                f"- Test: {item['test']} ({item['assumptions']})",
                f"- Mean difference: {item['difference']:.6f} m or m/yr",
                f"- 95% CI: [{item['ci95_low']:.6f}, {item['ci95_high']:.6f}]",
                f"- Cohen's d: {item['effect_size_cohens_d']:.3f}",
                f"- Raw p-value: {item['pvalue']:.3e}",
                f"- FDR-adjusted p-value: {item['fdr_pvalue']:.3e}",
                "",
            ]
        )
    (RESULTS_DIR / "statistical-summary.md").write_text("\n".join(stats_md), encoding="utf-8")
    files_written.append("results/statistical-summary.md")

    print(f"PSC selected: {summary['psc_count']} ({summary['psc_fraction']:.2%})")
    print(f"SBAS pairs: {summary['pair_count']} | network rank: {summary['network_rank']}")
    print(f"Representative LOS velocity: {summary['linear_velocity_m_per_yr']:.5f} m/yr")
    print(f"Precursor alert: {precursor_alert['level']} | Nankai alert: {nankai_alert['level']}")
    print(f"Velocity RMSE: {rmse_velocity:.5f} m/yr | 3D Up RMSE: {rmse_up_final:.5f} m")

    log_event("verify", "file_written", files_written)
    log_event("verify", "run_completed", files_written)


if __name__ == "__main__":
    main()
