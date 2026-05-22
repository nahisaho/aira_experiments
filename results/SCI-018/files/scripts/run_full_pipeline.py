from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from amr_framework import RESULTS_DIR, ROOT as PACKAGE_ROOT, ensure_output_dirs, log_event, save_json, seed_everything
from amr_framework.arg_detection import run_component as run_arg_detection
from amr_framework.evolutionary_paths import run_component as run_evolutionary_paths
from amr_framework.fitness_landscape import run_component as run_fitness_landscape
from amr_framework.hgt_network import run_component as run_hgt_network
from amr_framework.spatiotemporal import run_component as run_spatiotemporal
from amr_framework.treatment_optimizer import run_component as run_treatment_optimizer


def main() -> None:
    ensure_output_dirs()
    seed_everything(42)
    log_event(
        phase="pipeline",
        event_type="run_started",
        skill_or_tool="run_full_pipeline.py",
        handoff_in={"cwd": str(PACKAGE_ROOT), "seed": 42},
    )
    log_event(
        phase="pipeline",
        event_type="prompt_received",
        skill_or_tool="user_request",
        handoff_in={"objective": "Build a comprehensive AMR evolution prediction framework with six functional components and reproducible outputs."},
    )
    log_event(
        phase="pipeline",
        event_type="skill_selected",
        skill_or_tool="co-scientist-infectious-disease",
        handoff_out={"reason": "AMR prediction, pathogen genomics, epidemiological modeling, and treatment optimization are in scope."},
    )

    components = [
        ("Component 1/6 - ARG detection", "component_1", run_arg_detection),
        ("Component 2/6 - Fitness landscape", "component_2", run_fitness_landscape),
        ("Component 3/6 - Evolutionary paths", "component_3", run_evolutionary_paths),
        ("Component 4/6 - HGT network", "component_4", run_hgt_network),
        ("Component 5/6 - Spatiotemporal dynamics", "component_5", run_spatiotemporal),
        ("Component 6/6 - Treatment optimization", "component_6", run_treatment_optimizer),
    ]

    summary = {}
    for label, phase, runner in components:
        print(f"[AMR] {label} starting...")
        log_event(phase=phase, event_type="handoff_started", skill_or_tool=label)
        result = runner(seed=42)
        summary[phase] = result
        print(f"[AMR] {label} completed.")

    summary_path = RESULTS_DIR / "pipeline_summary.json"
    save_json(summary_path, summary)
    log_event(
        phase="pipeline",
        event_type="run_completed",
        skill_or_tool="run_full_pipeline.py",
        handoff_out={"summary_path": str(summary_path.relative_to(PACKAGE_ROOT)), "components": list(summary.keys())},
        files_written=[str(summary_path.relative_to(PACKAGE_ROOT))],
    )
    print("[AMR] Pipeline summary saved to results/pipeline_summary.json")
    print("[AMR] Full pipeline completed successfully.")


if __name__ == "__main__":
    main()
