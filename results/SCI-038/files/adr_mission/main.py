from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from adr_mission import LOGS_DIR, RESULTS_DIR
from adr_mission.capture_mechanism import analyze_capture_mechanisms
from adr_mission.debris_catalog import generate_debris_catalog
from adr_mission.debris_rotation import analyze_rotation
from adr_mission.mission_optimizer import optimize_mission
from adr_mission.orbit_transition import build_delta_v_matrix
from adr_mission.rendezvous import simulate_rendezvous_scenarios
from adr_mission.target_selection import score_targets

LOG_PATH = LOGS_DIR / "process-log.jsonl"
SUMMARY_PATH = RESULTS_DIR / "mission_summary.json"
STAT_SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"


def log_event(phase: str, event_type: str, skill_or_tool: str, handoff_in: dict, handoff_out: dict, files_written: list[str], status: str = "ok") -> None:
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in,
        "handoff_out": handoff_out,
        "files_written": files_written,
        "status": status,
    }
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, ensure_ascii=False) + "\n")


def main() -> None:
    log_event("EXECUTE", "skill_selected", "co-scientist-data-analysis", {"pipeline": "ADR mission design"}, {"modules": 7}, [])

    log_event("EXECUTE", "handoff_started", "debris_catalog", {"objects_requested": 50}, {}, [])
    catalog = generate_debris_catalog()
    log_event("EXECUTE", "handoff_completed", "debris_catalog", {"objects_requested": 50}, {"objects_generated": int(len(catalog))}, ["data/debris_catalog.csv"])
    log_event("EXECUTE", "file_written", "debris_catalog", {}, {"files": ["data/debris_catalog.csv"]}, ["data/debris_catalog.csv"])

    log_event("EXECUTE", "handoff_started", "target_selection", {"catalog_size": int(len(catalog))}, {}, [])
    scored, top_targets = score_targets(catalog)
    log_event("EXECUTE", "handoff_completed", "target_selection", {"catalog_size": int(len(catalog))}, {"top_target": top_targets.iloc[0]["debris_id"]}, ["results/target_scores.csv", "figures/target_selection.png"])
    log_event("EXECUTE", "file_written", "target_selection", {}, {"files": ["results/target_scores.csv", "figures/target_selection.png", "results/target_selection_summary.json"]}, ["results/target_scores.csv", "figures/target_selection.png", "results/target_selection_summary.json"])

    log_event("EXECUTE", "handoff_started", "orbit_transition", {"targets": len(top_targets)}, {}, [])
    dv_df, time_df = build_delta_v_matrix(top_targets)
    log_event("EXECUTE", "handoff_completed", "orbit_transition", {"targets": len(top_targets)}, {"max_delta_v_m_s": float(dv_df.to_numpy().max())}, ["results/delta_v_matrix.csv", "results/transfer_time_matrix_days.csv", "figures/delta_v_heatmap.png"])
    log_event("EXECUTE", "file_written", "orbit_transition", {}, {"files": ["results/delta_v_matrix.csv", "results/transfer_time_matrix_days.csv", "figures/delta_v_heatmap.png", "results/orbit_transition_summary.json"]}, ["results/delta_v_matrix.csv", "results/transfer_time_matrix_days.csv", "figures/delta_v_heatmap.png", "results/orbit_transition_summary.json"])

    log_event("EXECUTE", "handoff_started", "rendezvous", {"scenarios": 3}, {}, [])
    rendezvous_df = simulate_rendezvous_scenarios()
    log_event("EXECUTE", "handoff_completed", "rendezvous", {"scenarios": 3}, {"rows": int(len(rendezvous_df))}, ["results/rendezvous_trajectories.csv", "figures/rendezvous_trajectory.png"])
    log_event("EXECUTE", "file_written", "rendezvous", {}, {"files": ["results/rendezvous_trajectories.csv", "figures/rendezvous_trajectory.png", "results/rendezvous_summary.json"]}, ["results/rendezvous_trajectories.csv", "figures/rendezvous_trajectory.png", "results/rendezvous_summary.json"])

    log_event("EXECUTE", "handoff_started", "debris_rotation", {"integration_points": 2400}, {}, [])
    rotation_df, rotation_summary = analyze_rotation()
    log_event("EXECUTE", "handoff_completed", "debris_rotation", {"integration_points": int(len(rotation_df))}, rotation_summary, ["results/rotation_analysis.csv", "figures/debris_rotation.png"])
    log_event("EXECUTE", "file_written", "debris_rotation", {}, {"files": ["results/rotation_analysis.csv", "figures/debris_rotation.png", "results/rotation_summary.json"]}, ["results/rotation_analysis.csv", "figures/debris_rotation.png", "results/rotation_summary.json"])

    log_event("EXECUTE", "handoff_started", "capture_mechanism", {"samples": 121}, {}, [])
    capture_df, capture_summary = analyze_capture_mechanisms()
    log_event("EXECUTE", "handoff_completed", "capture_mechanism", {"samples": int(len(capture_df))}, capture_summary, ["results/capture_analysis.csv", "figures/capture_mechanisms.png"])
    log_event("EXECUTE", "file_written", "capture_mechanism", {}, {"files": ["results/capture_analysis.csv", "figures/capture_mechanisms.png", "results/capture_summary.json"]}, ["results/capture_analysis.csv", "figures/capture_mechanisms.png", "results/capture_summary.json"])

    log_event("EXECUTE", "handoff_started", "mission_optimizer", {"candidate_targets": 10}, {}, [])
    mission_result = optimize_mission(top_targets, dv_df, time_df)
    log_event("EXECUTE", "handoff_completed", "mission_optimizer", {"candidate_targets": 10}, {"selected_route_length": len(mission_result["selected_route"]), "feasible": mission_result["feasible"]}, ["results/optimal_mission_sequence.json", "figures/mission_optimization.png"])
    log_event("EXECUTE", "file_written", "mission_optimizer", {}, {"files": ["results/optimal_mission_sequence.json", "figures/mission_optimization.png"]}, ["results/optimal_mission_sequence.json", "figures/mission_optimization.png"])

    summary = {
        "catalog_size": int(len(catalog)),
        "top_target": top_targets.iloc[0]["debris_id"],
        "top_target_score": float(top_targets.iloc[0]["combined_score"]),
        "mean_decay_lifetime_days": float(catalog["decay_lifetime_days"].mean()),
        "max_pairwise_delta_v_m_s": float(dv_df.to_numpy().max()),
        "mean_pairwise_delta_v_m_s": float(dv_df.replace(0.0, np.nan).stack().mean()),
        "best_rendezvous_delta_v_m_s": float(rendezvous_df.groupby("scenario")["delta_v_total_m_s"].first().min()),
        "estimated_rotation_period_s": float(rotation_summary["estimated_period_s"]),
        "selected_mission_delta_v_m_s": float(mission_result["total_delta_v_m_s"]),
        "selected_mission_time_days": float(mission_result["total_time_days"]),
        "selected_mission_targets_serviced": int(mission_result["targets_serviced"]),
        "selected_mission_feasible": bool(mission_result["feasible"]),
        "selected_route": mission_result["selected_route"],
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    stats_md = f"""# Statistical Summary\n\n- Catalog size: {summary['catalog_size']} debris objects\n- Top target: {summary['top_target']} (score={summary['top_target_score']:.4f})\n- Mean decay lifetime: {summary['mean_decay_lifetime_days']:.2f} days\n- Mean pairwise ΔV among top targets: {summary['mean_pairwise_delta_v_m_s']:.2f} m/s\n- Best rendezvous ΔV: {summary['best_rendezvous_delta_v_m_s']:.2f} m/s\n- Estimated debris rotation period: {summary['estimated_rotation_period_s']:.2f} s\n- Selected mission route ΔV: {summary['selected_mission_delta_v_m_s']:.2f} m/s\n- Selected mission duration: {summary['selected_mission_time_days']:.2f} days\n- Targets serviced within budget: {summary['selected_mission_targets_serviced']}\n- Fuel-budget feasibility: {summary['selected_mission_feasible']}\n"""
    STAT_SUMMARY_PATH.write_text(stats_md, encoding="utf-8")
    log_event("REPORT", "report_finalized", "main_summary", {"summary_keys": list(summary.keys())}, summary, ["results/mission_summary.json", "results/statistical-summary.md"])
    log_event("LOG", "run_completed", "adr_mission.main", {"status": "complete"}, {"summary_file": "results/mission_summary.json"}, ["logs/process-log.jsonl"], status="ok")

    print("ADR mission design pipeline completed.")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
