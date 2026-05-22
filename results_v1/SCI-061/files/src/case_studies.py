import json
import os
import random
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - import guard
    raise RuntimeError(f"matplotlib with Agg backend is required: {exc}") from exc

try:
    from circuit_spec import CircuitSpec
    from pipeline import AutoDesignPipeline
    from stochastic_sim import GillespieSimulator, TauLeapingSimulator, run_ensemble
except Exception as exc:  # pragma: no cover - import guard
    raise RuntimeError(f"Case study imports failed: {exc}") from exc


WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES_DIR = os.path.join(WORKSPACE_ROOT, "figures")
RESULTS_DIR = os.path.join(WORKSPACE_ROOT, "results")
LOGS_DIR = os.path.join(WORKSPACE_ROOT, "logs")
DATA_DIR = os.path.join(WORKSPACE_ROOT, "data")
PROCESS_LOG = os.path.join(LOGS_DIR, "process-log.jsonl")
REPORT_PATH = os.path.join(WORKSPACE_ROOT, "report.md")
PREPROCESS_LOG = os.path.join(DATA_DIR, "preprocessing-log.md")


def _ensure_directories():
    for path in (FIGURES_DIR, RESULTS_DIR, LOGS_DIR, DATA_DIR):
        os.makedirs(path, exist_ok=True)


def _timestamp():
    return datetime.now(timezone.utc).isoformat()


def _log_event(phase, event_type, skill_or_tool, handoff_in=None, handoff_out=None, files_written=None, status="ok"):
    _ensure_directories()
    payload = {
        "timestamp": _timestamp(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files_written or [],
        "status": status,
    }
    with open(PROCESS_LOG, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload) + "\n")


def _write_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(_jsonify(payload), handle, indent=2, sort_keys=True)


def _jsonify(value):
    if isinstance(value, dict):
        return {str(key): _jsonify(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if hasattr(value, "to_sbol_dict"):
        return value.to_sbol_dict()
    if hasattr(value, "to_dict"):
        return value.to_dict()
    return value


def _record_preprocessing():
    _ensure_directories()
    text = (
        "# Preprocessing Log\n\n"
        "- No external datasets were imported.\n"
        "- Random seeds were fixed for numpy and random.\n"
        "- Figures were rendered with the Agg backend at 300 DPI.\n"
        "- Numerical outputs were exported as JSON summaries only.\n"
    )
    with open(PREPROCESS_LOG, "w", encoding="utf-8") as handle:
        handle.write(text)


def _design_names(design_candidate):
    return [part.name for part in getattr(design_candidate, "parts", [])]


def _build_toggle_spec():
    return CircuitSpec.toggle_switch()


def _build_repressilator_spec():
    spec = CircuitSpec(name="repressilator")
    for signal_name in ("TetR", "CI", "LacI"):
        spec.add_signal(signal_name)
    spec.add_signal("GFP", signal_type="output")
    spec.add_gate("g_tetR", "NOT", [], "TetR")
    spec.add_gate("g_ci", "NOT", [], "CI")
    spec.add_gate("g_lacI", "NOT", [], "LacI")
    spec.add_gate("g_gfp", "BUFFER", ["LacI"], "GFP")
    spec.add_feedback("LacI", "g_tetR")
    spec.add_feedback("TetR", "g_ci")
    spec.add_feedback("CI", "g_lacI")
    return spec


def _toggle_bistability_analysis(results, species_a="LacI", species_b="TetR"):
    final_a = np.array([run.trajectories[species_a][-1] for run in results], dtype=float)
    final_b = np.array([run.trajectories[species_b][-1] for run in results], dtype=float)
    diff = final_a - final_b
    threshold = max(5.0, 0.2 * float(np.mean(final_a + final_b + 1.0)))
    laci_high = diff > threshold
    tetr_high = diff < -threshold
    centers = {
        "LacI_high": {
            "LacI": float(np.mean(final_a[laci_high])) if np.any(laci_high) else 0.0,
            "TetR": float(np.mean(final_b[laci_high])) if np.any(laci_high) else 0.0,
        },
        "TetR_high": {
            "LacI": float(np.mean(final_a[tetr_high])) if np.any(tetr_high) else 0.0,
            "TetR": float(np.mean(final_b[tetr_high])) if np.any(tetr_high) else 0.0,
        },
    }
    separation = float(np.linalg.norm([
        centers["LacI_high"]["LacI"] - centers["TetR_high"]["LacI"],
        centers["LacI_high"]["TetR"] - centers["TetR_high"]["TetR"],
    ]))
    is_bistable = bool(np.mean(laci_high) >= 0.15 and np.mean(tetr_high) >= 0.15 and separation > threshold)
    return {
        "is_bistable": is_bistable,
        "threshold": float(threshold),
        "cluster_fractions": {
            "LacI_high": float(np.mean(laci_high)),
            "TetR_high": float(np.mean(tetr_high)),
        },
        "cluster_centers": centers,
        "separation": separation,
        "steady_states": [
            {species_a: float(a), species_b: float(b)}
            for a, b in zip(final_a.tolist(), final_b.tolist())
        ],
    }


def _oscillation_analysis(results, species="LacI"):
    periods = []
    amplitudes = []
    phase_portrait = None
    for run in results:
        grid = np.linspace(float(run.time_points[0]), float(run.time_points[-1]), 600)
        signal = np.interp(grid, run.time_points, run.trajectories[species])
        centered = signal - float(np.mean(signal))
        crossings = np.where(np.diff(np.signbit(centered)))[0]
        if len(crossings) >= 3:
            half_periods = np.diff(grid[crossings])
            periods.append(2.0 * float(np.mean(half_periods)))
        amplitudes.append(0.5 * float(np.quantile(signal, 0.95) - np.quantile(signal, 0.05)))
        if phase_portrait is None and "TetR" in run.trajectories and "CI" in run.trajectories:
            phase_portrait = {
                "TetR": np.interp(grid, run.time_points, run.trajectories["TetR"]).tolist(),
                "CI": np.interp(grid, run.time_points, run.trajectories["CI"]).tolist(),
            }
    return {
        "species": species,
        "mean_period": float(np.mean(periods)) if periods else None,
        "periods": [float(val) for val in periods],
        "mean_amplitude": float(np.mean(amplitudes)) if amplitudes else 0.0,
        "amplitudes": [float(val) for val in amplitudes],
        "phase_portrait": phase_portrait or {"TetR": [], "CI": []},
    }


def _plot_toggle_trajectories(results, path):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for run in results:
        axes[0].plot(run.time_points, run.trajectories["LacI"], alpha=0.35, color="#440154")
        axes[1].plot(run.time_points, run.trajectories["TetR"], alpha=0.35, color="#21918c")
    axes[0].set_ylabel("LacI copies")
    axes[1].set_ylabel("TetR copies")
    axes[1].set_xlabel("Time")
    axes[0].set_title("Toggle switch trajectories")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_toggle_steady_states(analysis, path):
    states = analysis["steady_states"]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        [item["LacI"] for item in states],
        [item["TetR"] for item in states],
        c=np.linspace(0, 1, len(states)),
        cmap="viridis",
        edgecolor="black",
        linewidth=0.3,
    )
    ax.set_xlabel("LacI steady state")
    ax.set_ylabel("TetR steady state")
    ax.set_title("Toggle switch steady states")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_repressilator_oscillations(results, path):
    reference = results[0]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for species, color in (("TetR", "#440154"), ("CI", "#3b528b"), ("LacI", "#5ec962")):
        ax.plot(reference.time_points, reference.trajectories[species], label=species, color=color, linewidth=1.5)
    ax.set_xlabel("Time")
    ax.set_ylabel("Protein copies")
    ax.set_title("Repressilator oscillations")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_phase_portrait(phase_portrait, path):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot(phase_portrait.get("TetR", []), phase_portrait.get("CI", []), color="#21918c", linewidth=1.4)
    ax.set_xlabel("TetR")
    ax.set_ylabel("CI")
    ax.set_title("Repressilator phase portrait")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_design_comparison(toggle_summary, repressilator_summary, path):
    labels = ["Combined", "Robustness", "Context", "Dynamic"]
    toggle_values = [
        toggle_summary["pipeline_combined_score"],
        toggle_summary["robustness_score"],
        toggle_summary["context_score"],
        toggle_summary["dynamic_score"],
    ]
    repressilator_values = [
        repressilator_summary["pipeline_combined_score"],
        repressilator_summary["robustness_score"],
        repressilator_summary["context_score"],
        repressilator_summary["dynamic_score"],
    ]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(x - width / 2, toggle_values, width, label="Toggle switch", color="#440154")
    ax.bar(x + width / 2, repressilator_values, width, label="Repressilator", color="#35b779")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title("Best design comparison")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _write_report(toggle_summary, repressilator_summary, summary_path):
    text = f"""# DRAFT — NOT FOR DISTRIBUTION

Generated: {_timestamp()}

## Methods

- AutoDesignPipeline ranked candidate designs using simulated performance, robustness, and context quality.
- Toggle switch redesign used Gillespie ensemble simulation (100 time units, 20 runs).
- Repressilator redesign used tau-leaping ensemble simulation (200 time units, 10 runs).
- Random seeds were fixed for reproducibility.

## Results

### Toggle switch
- Combined score: {toggle_summary['pipeline_combined_score']:.3f}
- Robustness score: {toggle_summary['robustness_score']:.3f}
- Context score: {toggle_summary['context_score']:.3f}
- Bistable: {toggle_summary['bistability']['is_bistable']}

### Repressilator
- Combined score: {repressilator_summary['pipeline_combined_score']:.3f}
- Robustness score: {repressilator_summary['robustness_score']:.3f}
- Context score: {repressilator_summary['context_score']:.3f}
- Mean oscillation period: {repressilator_summary['oscillation'].get('mean_period')}

## Discussion

These case studies provide lightweight computational redesign examples suitable for rapid screening. The scores are heuristic and should be interpreted as design-stage priorities rather than definitive biological predictions.

## File inventory

- figures/toggle_switch_trajectories.png
- figures/toggle_switch_steady_states.png
- figures/repressilator_oscillations.png
- figures/repressilator_phase_portrait.png
- figures/design_comparison.png
- results/toggle_switch_results.json
- results/repressilator_results.json
- results/pipeline_summary.json
- data/preprocessing-log.md
- logs/process-log.jsonl
"""
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(text)


def run_toggle_switch_study():
    """
    Toggle Switch Redesign:
    1. Define toggle switch circuit spec (two mutually repressing genes: LacI ⊣ TetR, TetR ⊣ LacI)
    2. Run the auto-design pipeline
    3. Simulate best design with Gillespie (short run, ~100 time units, 20 runs)
    4. Analyze bistability (check for two distinct steady states)
    5. Compute robustness score
    6. Apply context corrections
    7. Save results to results/ and figures to figures/
    """
    _ensure_directories()
    np.random.seed(7)
    random.seed(7)
    _log_event("PLAN", "handoff_started", "toggle_case_study", {"study": "toggle_switch"})
    spec = _build_toggle_spec()
    pipeline = AutoDesignPipeline()
    pipeline_result = pipeline.run(spec, target_behavior={"mode": "toggle", "expression": 1.0}, n_candidates=4, n_sim_runs=4)
    best_design, pipeline_sim, robustness_score, context_score = pipeline_result.ranked_designs[0]
    model = pipeline.build_simulation_model(spec, best_design)
    gillespie = GillespieSimulator(max_steps=150000)
    simulation_runs, ensemble_stats = run_ensemble(gillespie, model, t_end=100.0, n_runs=20, seed=11)
    bistability = _toggle_bistability_analysis(simulation_runs)
    figure_paths = {
        "trajectories": os.path.join(FIGURES_DIR, "toggle_switch_trajectories.png"),
        "steady_states": os.path.join(FIGURES_DIR, "toggle_switch_steady_states.png"),
    }
    _plot_toggle_trajectories(simulation_runs, figure_paths["trajectories"])
    _plot_toggle_steady_states(bistability, figure_paths["steady_states"])
    result_path = os.path.join(RESULTS_DIR, "toggle_switch_results.json")
    payload = {
        "study": "toggle_switch",
        "pipeline_execution_time": pipeline_result.execution_time,
        "pipeline_summary": pipeline_result.summary_statistics,
        "pipeline_combined_score": pipeline_sim["combined_score"],
        "best_design_parts": _design_names(best_design),
        "best_design_sbol": best_design.to_sbol_dict(),
        "dynamic_score": pipeline_sim["dynamic_metrics"]["performance_score"],
        "bistability": bistability,
        "ensemble_statistics": ensemble_stats,
        "robustness_score": robustness_score,
        "robustness_summary": pipeline_sim["robustness_summary"],
        "context_score": context_score,
        "context_summary": pipeline_sim["context_summary"],
        "figure_paths": figure_paths,
    }
    _write_json(result_path, payload)
    _log_event(
        "EXECUTE",
        "handoff_completed",
        "toggle_case_study",
        {"study": "toggle_switch"},
        {"best_design_parts": _design_names(best_design), "bistable": bistability["is_bistable"]},
        [result_path, figure_paths["trajectories"], figure_paths["steady_states"]],
    )
    return payload


def run_repressilator_study():
    """
    Repressilator Redesign:
    1. Define repressilator circuit (3-node ring: TetR ⊣ CI ⊣ LacI ⊣ TetR)
    2. Run auto-design pipeline
    3. Simulate with tau-leaping (faster for oscillatory dynamics, ~200 time units, 10 runs)
    4. Analyze oscillation period and amplitude
    5. Compute robustness
    6. Context effect analysis
    7. Save results and figures
    """
    _ensure_directories()
    np.random.seed(13)
    random.seed(13)
    _log_event("PLAN", "handoff_started", "repressilator_case_study", {"study": "repressilator"})
    spec = _build_repressilator_spec()
    pipeline = AutoDesignPipeline()
    pipeline_result = pipeline.run(spec, target_behavior={"mode": "repressilator", "expression": 1.0}, n_candidates=4, n_sim_runs=4)
    best_design, pipeline_sim, robustness_score, context_score = pipeline_result.ranked_designs[0]
    model = pipeline.build_simulation_model(spec, best_design)
    tau_leaping = TauLeapingSimulator(max_steps=200000, epsilon=0.04, exact_threshold=8)
    simulation_runs, ensemble_stats = run_ensemble(tau_leaping, model, t_end=200.0, n_runs=10, seed=19)
    oscillation = _oscillation_analysis(simulation_runs, species="LacI")
    figure_paths = {
        "oscillations": os.path.join(FIGURES_DIR, "repressilator_oscillations.png"),
        "phase_portrait": os.path.join(FIGURES_DIR, "repressilator_phase_portrait.png"),
    }
    _plot_repressilator_oscillations(simulation_runs, figure_paths["oscillations"])
    _plot_phase_portrait(oscillation["phase_portrait"], figure_paths["phase_portrait"])
    result_path = os.path.join(RESULTS_DIR, "repressilator_results.json")
    payload = {
        "study": "repressilator",
        "pipeline_execution_time": pipeline_result.execution_time,
        "pipeline_summary": pipeline_result.summary_statistics,
        "pipeline_combined_score": pipeline_sim["combined_score"],
        "best_design_parts": _design_names(best_design),
        "best_design_sbol": best_design.to_sbol_dict(),
        "dynamic_score": pipeline_sim["dynamic_metrics"]["performance_score"],
        "oscillation": oscillation,
        "ensemble_statistics": ensemble_stats,
        "robustness_score": robustness_score,
        "robustness_summary": pipeline_sim["robustness_summary"],
        "context_score": context_score,
        "context_summary": pipeline_sim["context_summary"],
        "figure_paths": figure_paths,
    }
    _write_json(result_path, payload)
    _log_event(
        "EXECUTE",
        "handoff_completed",
        "repressilator_case_study",
        {"study": "repressilator"},
        {"best_design_parts": _design_names(best_design), "mean_period": oscillation.get("mean_period")},
        [result_path, figure_paths["oscillations"], figure_paths["phase_portrait"]],
    )
    return payload


def run_all_case_studies():
    """Run both and save summary"""
    _ensure_directories()
    _record_preprocessing()
    _log_event("PLAN", "run_started", "case_studies", {"workspace": WORKSPACE_ROOT})
    _log_event("PLAN", "prompt_received", "case_studies", {"requested_outputs": [
        "figures/toggle_switch_trajectories.png",
        "figures/toggle_switch_steady_states.png",
        "figures/repressilator_oscillations.png",
        "figures/repressilator_phase_portrait.png",
        "figures/design_comparison.png",
        "results/toggle_switch_results.json",
        "results/repressilator_results.json",
        "results/pipeline_summary.json",
    ]})
    _log_event("PLAN", "skill_selected", "co-scientist-data-analysis", {"scope": "case studies"})
    toggle_summary = run_toggle_switch_study()
    repressilator_summary = run_repressilator_study()
    comparison_path = os.path.join(FIGURES_DIR, "design_comparison.png")
    _plot_design_comparison(toggle_summary, repressilator_summary, comparison_path)
    summary_path = os.path.join(RESULTS_DIR, "pipeline_summary.json")
    summary_payload = {
        "generated_at": _timestamp(),
        "toggle_switch": toggle_summary,
        "repressilator": repressilator_summary,
        "comparison_figure": comparison_path,
        "files_generated": [
            "figures/toggle_switch_trajectories.png",
            "figures/toggle_switch_steady_states.png",
            "figures/repressilator_oscillations.png",
            "figures/repressilator_phase_portrait.png",
            "figures/design_comparison.png",
            "results/toggle_switch_results.json",
            "results/repressilator_results.json",
            "results/pipeline_summary.json",
        ],
    }
    _write_json(summary_path, summary_payload)
    _write_report(toggle_summary, repressilator_summary, REPORT_PATH)
    _log_event("REPORT", "report_finalized", "case_studies", files_written=[REPORT_PATH, summary_path, comparison_path])
    _log_event("LOG", "run_completed", "case_studies", handoff_out={"summary": summary_path}, files_written=[PROCESS_LOG], status="ok")
    return summary_payload


if __name__ == '__main__':
    run_all_case_studies()
