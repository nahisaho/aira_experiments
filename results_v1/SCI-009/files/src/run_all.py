"""
PROTAC Computational Design Framework — Master Runner
Executes all 6 modules in sequence and collects results for the final report.
"""
import sys, os, json, datetime
sys.path.insert(0, os.path.dirname(__file__))

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)
os.makedirs("logs", exist_ok=True)

from protac_utils import log_event

def main():
    start = datetime.datetime.utcnow()
    log_event("pipeline", "run_started", "co-scientist",
              {"modules": 6, "target": "BRD4 PROTAC design"})

    print("=" * 65)
    print("  PROTAC Computational Design Framework")
    print("  Rosetta/AmberTools-inspired workflow (RDKit/MMFF94)")
    print("=" * 65)

    results = {}

    # ---- Module 1: Ternary complex modeling ----
    print("\n" + "─" * 65)
    from importlib import import_module
    m1 = import_module("01_ternary_complex_modeling")
    results["ternary"] = m1.run_ternary_modeling()

    # ---- Module 2: Linker optimization ----
    print("\n" + "─" * 65)
    import matplotlib.patches as mpatches  # ensure imported
    m2 = import_module("02_linker_optimization")
    results["linker"] = m2.run_linker_optimization()

    # ---- Module 3: E3 selectivity ----
    print("\n" + "─" * 65)
    m3 = import_module("03_e3_selectivity_prediction")
    results["e3"] = m3.run_e3_selectivity()

    # ---- Module 4: ADMET ----
    print("\n" + "─" * 65)
    m4 = import_module("04_admet_prediction")
    admet_df, admet_models = m4.run_admet_prediction()
    results["admet"] = (admet_df, admet_models)

    # ---- Module 5: SAR analysis ----
    print("\n" + "─" * 65)
    m5 = import_module("05_sar_analysis")
    sar_df, sar_models, cliff_df = m5.run_sar_analysis()
    results["sar"] = (sar_df, sar_models, cliff_df)

    # ---- Module 6: BRD4 case study ----
    print("\n" + "─" * 65)
    m6 = import_module("06_brd4_case_study")
    results["brd4"] = m6.run_brd4_case_study()

    elapsed = (datetime.datetime.utcnow() - start).total_seconds()
    print("\n" + "=" * 65)
    print(f"  All modules completed in {elapsed:.1f}s")
    print("=" * 65)

    log_event("pipeline", "run_completed", "co-scientist",
              {"elapsed_s": round(elapsed, 1), "status": "success"})

    # Collect summary metrics
    summary = {
        "ternary_best_linker": str(results["ternary"].iloc[0]["linker_name"]),
        "linker_opt_best": str(results["linker"].iloc[0]["linker_name"]),
        "e3_roc_auc": float(results["e3"]["roc_auc"]),
        "brd4_best_DC50_nM": float(results["brd4"]["DC50_nM"].min()),
        "brd4_best_Dmax_pct": float(results["brd4"]["Dmax_pct"].max()),
        "n_activity_cliffs": len(cliff_df),
        "sar_pDC50_R2": float(sar_models["pDC50"]["r2_cv"]),
        "sar_Dmax_R2": float(sar_models["Dmax_pct"]["r2_cv"]),
    }
    with open("results/pipeline_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\nKey Results:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    return summary, results


if __name__ == "__main__":
    main()
