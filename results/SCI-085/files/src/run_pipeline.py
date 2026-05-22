"""
Perturb-seq Analysis Pipeline — Main Orchestrator
===================================================
Runs all 6 modules sequentially with logging.
"""

import os
import sys
import json
from datetime import datetime

# Ensure src is in path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LOG_DIR = "logs"
os.makedirs(LOG_DIR, exist_ok=True)

LOG_FILE = os.path.join(LOG_DIR, "process-log.jsonl")


def log_event(phase, event_type, skill="pipeline", status="ok", **kwargs):
    entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "status": status,
        **kwargs,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


def main():
    print("=" * 70)
    print("  Perturb-seq Analysis Pipeline")
    print("=" * 70)
    log_event("pipeline", "run_started")

    # Module 0: Data Simulation
    print("\n[0/6] Generating simulated data...")
    log_event("setup", "skill_selected", skill="data-simulation")
    import importlib
    mod0 = importlib.import_module("src.00_setup")
    mod0.simulate_perturbseq()
    log_event("setup", "handoff_completed", files_written=["data/perturbseq_simulated.h5ad"])

    # Module 1: QC & Guide Detection
    print("\n[1/6] Quality control & guide detection...")
    log_event("qc", "skill_selected", skill="perturbation-analysis")
    mod1 = importlib.import_module("src.01_qc_guide_detection")
    mod1.run_qc_pipeline()
    log_event("qc", "handoff_completed",
              files_written=["data/perturbseq_qc_filtered.h5ad", "results/01_qc_stats.json"])

    # Module 2: Differential Expression
    print("\n[2/6] Differential expression & co-expression modules...")
    log_event("de", "skill_selected", skill="perturbation-analysis")
    mod2 = importlib.import_module("src.02_differential_expression")
    mod2.run_de_pipeline()
    log_event("de", "handoff_completed",
              files_written=["results/02_de_results.csv", "results/02_de_summary.json"])

    # Module 3: Causal Graph
    print("\n[3/6] Causal graph estimation...")
    log_event("causal", "skill_selected", skill="causal-inference")
    mod3 = importlib.import_module("src.03_causal_graph")
    mod3.run_causal_pipeline()
    log_event("causal", "handoff_completed",
              files_written=["results/03_causal_edges.csv", "results/03_causal_summary.json"])

    # Module 4: Epistasis
    print("\n[4/6] Epistasis detection...")
    log_event("epistasis", "skill_selected", skill="perturbation-analysis")
    mod4 = importlib.import_module("src.04_epistasis")
    mod4.run_epistasis_pipeline()
    log_event("epistasis", "handoff_completed",
              files_written=["results/04_epistasis_results.csv", "results/04_epistasis_summary.json"])

    # Module 5: Latent Representation
    print("\n[5/6] Latent representation learning...")
    log_event("latent", "skill_selected", skill="deep-learning")
    mod5 = importlib.import_module("src.05_latent_representation")
    mod5.run_latent_pipeline()
    log_event("latent", "handoff_completed",
              files_written=["results/05_perturbation_similarity.csv", "results/05_latent_summary.json"])

    # Module 6: Essential Gene Network
    print("\n[6/6] Essential gene network case study...")
    log_event("essential", "skill_selected", skill="network-analysis")
    mod6 = importlib.import_module("src.06_essential_gene_network")
    mod6.run_essential_pipeline()
    log_event("essential", "handoff_completed",
              files_written=["results/06_essential_summary.json", "results/06_network_edges.csv"])

    # Final
    log_event("pipeline", "run_completed")
    print("\n" + "=" * 70)
    print("  Pipeline complete! Check results/ and figures/ for outputs.")
    print("=" * 70)


if __name__ == "__main__":
    main()
