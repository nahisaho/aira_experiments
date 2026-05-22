#!/usr/bin/env python3
"""GEM制約条件ベースフラックス解析 — 全モジュール実行"""
import os, sys, json
from datetime import datetime

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ".")

log = []
def log_event(phase, evt, skill="", files=None, status="ok"):
    log.append({"timestamp": datetime.now().isoformat(), "phase": phase, "event_type": evt,
                "actor": "co-scientist", "skill_or_tool": skill, "files_written": files or [], "status": status})

log_event("init", "run_started", "metabolic-modeling")

print("=" * 70)
print("  GEM Constraint-Based Flux Analysis Framework — All Modules")
print("=" * 70)

from scripts.s01_fba_optimization import run_fba_optimization
from scripts.s02_13c_mfa_integration import run_13c_mfa_integration
from scripts.s03_dynamic_fba import run_dynamic_fba
from scripts.s04_enzyme_constraints import run_enzyme_constraints
from scripts.s05_condition_specific import run_condition_specific
from scripts.s06_lysine_optimization import run_lysine_optimization

print("\n" + "▶"*20 + " MODULE 1: FBA Optimization " + "▶"*20)
r1 = run_fba_optimization()
log_event("m1", "completed", "fba", ["figures/01_*", "results/01_*"])

print("\n" + "▶"*20 + " MODULE 2: 13C-MFA Integration " + "▶"*20)
r2 = run_13c_mfa_integration()
log_event("m2", "completed", "13c-mfa", ["figures/02_*", "results/02_*"])

print("\n" + "▶"*20 + " MODULE 3: Dynamic FBA " + "▶"*20)
r3 = run_dynamic_fba()
log_event("m3", "completed", "dfba", ["figures/03_*", "results/03_*"])

print("\n" + "▶"*20 + " MODULE 4: Enzyme Constraints " + "▶"*20)
r4 = run_enzyme_constraints()
log_event("m4", "completed", "gecko", ["figures/04_*", "results/04_*"])

print("\n" + "▶"*20 + " MODULE 5: Condition-Specific " + "▶"*20)
r5 = run_condition_specific()
log_event("m5", "completed", "gimme", ["figures/05_*", "results/05_*"])

print("\n" + "▶"*20 + " MODULE 6: Lysine Optimization " + "▶"*20)
r6 = run_lysine_optimization()
log_event("m6", "completed", "lysine", ["figures/06_*", "results/06_*"])

log_event("final", "run_completed")

os.makedirs("logs", exist_ok=True)
with open("logs/process-log.jsonl", "w") as f:
    for e in log: f.write(json.dumps(e) + "\n")

print("\n" + "=" * 70)
print("  ALL MODULES COMPLETED SUCCESSFULLY")
print("=" * 70)

# Summary
print(json.dumps({
    "m1_growth": r1.get("standard_fba",{}).get("objective_value"),
    "m2_chi2": r2.get("chi_square",{}).get("reduced"),
    "m3_max_biomass": r3.get("summary",{}).get("max_biomass"),
    "m4_budgets": len(r4.get("budget_analysis",[])),
    "m5_conditions": {c: r5.get(c,{}).get("growth_rate") for c in ["aerobic","anaerobic","stress"]},
    "m6_max_lysine": r6.get("max_theoretical_lysine"),
    "m6_strategies": r6.get("optknock_strategies"),
}, indent=2))
