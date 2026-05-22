"""
Main runner: Orchestrates all simulation modules and generates report data.
"""

import sys
import os
import json
import time

sys.path.insert(0, os.path.dirname(__file__))

from src.network_analysis import run_network_analysis
from src.metrics import run_metrics_analysis
from src.funding_mechanisms import run_funding_simulation
from src.abm_model import run_abm_simulation
from src.diversity_optimization import run_diversity_optimization
from src.kakenhi_case_study import run_kakenhi_case_study


def main():
    start_time = time.time()
    os.makedirs("logs", exist_ok=True)

    log_entries = []

    def log_event(phase, event_type, skill="", files=None, status="ok"):
        entry = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "phase": phase,
            "event_type": event_type,
            "actor": "co-scientist",
            "skill_or_tool": skill,
            "files_written": files or [],
            "status": status,
        }
        log_entries.append(entry)

    log_event("init", "run_started", "run_simulation.py")

    # === Module 1: Network Analysis ===
    print("=" * 60)
    print("Module 1: Research Network Analysis")
    print("=" * 60)
    log_event("network", "module_started", "network_analysis")
    coauth_G, citation_G, network_results = run_network_analysis()
    log_event("network", "module_completed", "network_analysis",
              ["figures/fig1_coauthorship_network.png",
               "figures/fig2_citation_network.png",
               "figures/fig3_coauth_degree_dist.png",
               "figures/fig4_citation_degree_dist.png",
               "results/network_metrics.json"])

    # === Module 2: Metrics Analysis ===
    print("\n" + "=" * 60)
    print("Module 2: Research Output Metrics")
    print("=" * 60)
    log_event("metrics", "module_started", "metrics")
    df_profiles, biases = run_metrics_analysis()
    log_event("metrics", "module_completed", "metrics",
              ["figures/fig5_metric_comparison.png",
               "data/researcher_profiles.csv",
               "results/metric_biases.json"])

    # === Module 3: Funding Mechanisms ===
    print("\n" + "=" * 60)
    print("Module 3: Funding Mechanism Simulation")
    print("=" * 60)
    log_event("funding", "module_started", "funding_mechanisms")
    funding_results = run_funding_simulation(df_profiles)
    log_event("funding", "module_completed", "funding_mechanisms",
              ["figures/fig6_mechanism_comparison.png",
               "results/funding_mechanisms.json"])

    # === Module 4: ABM Simulation ===
    print("\n" + "=" * 60)
    print("Module 4: Agent-Based Model Simulation")
    print("=" * 60)
    log_event("abm", "module_started", "abm_model")
    abm_results = run_abm_simulation(n_researchers=200, n_steps=20)
    log_event("abm", "module_completed", "abm_model",
              ["figures/fig7_abm_career_simulation.png",
               "results/abm_results.json"])

    # === Module 5: Diversity Optimization ===
    print("\n" + "=" * 60)
    print("Module 5: Diversity-Constrained Optimization")
    print("=" * 60)
    log_event("optimization", "module_started", "diversity_optimization")
    opt_results = run_diversity_optimization(df_profiles)
    log_event("optimization", "module_completed", "diversity_optimization",
              ["figures/fig8_pareto_frontier.png",
               "results/diversity_optimization.json",
               "results/pareto_points.csv"])

    # === Module 6: KAKENHI Case Study ===
    print("\n" + "=" * 60)
    print("Module 6: KAKENHI Case Study")
    print("=" * 60)
    log_event("kakenhi", "module_started", "kakenhi_case_study")
    kakenhi_results = run_kakenhi_case_study()
    log_event("kakenhi", "module_completed", "kakenhi_case_study",
              ["figures/fig9_kakenhi_case_study.png",
               "data/kakenhi_applicants.csv",
               "results/kakenhi_case_study.json"])

    # === Summary ===
    elapsed = time.time() - start_time
    print("\n" + "=" * 60)
    print(f"All modules complete. Total time: {elapsed:.1f}s")
    print("=" * 60)

    log_event("final", "run_completed", "run_simulation.py")

    # Write process log
    with open("logs/process-log.jsonl", "w") as f:
        for entry in log_entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print("Process log written to logs/process-log.jsonl")


if __name__ == "__main__":
    main()
