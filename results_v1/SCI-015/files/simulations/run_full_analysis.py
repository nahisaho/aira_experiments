"""
Main NCC Framework Simulation
==============================
Runs the complete analysis pipeline:
1. IIT Φ calculation across consciousness levels
2. PCI simulation (TMS-EEG)
3. Global Workspace Theory metrics
4. Clinical consciousness classification
5. Multi-theory comparison
6. AI consciousness implications

Run:
    python simulations/run_full_analysis.py
"""
import sys
import os
import json
import time
import numpy as np

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.ncc_framework import (
    PhiCalculator,
    PCISimulator,
    GlobalWorkspaceAnalyzer,
    ConsciousnessClassifier,
    generate_anesthesia_data,
)
from src.ncc_framework.clinical import (
    generate_clinical_dataset,
    CLASS_NAMES,
    STATE_CONSCIOUSNESS,
)
from src.ncc_framework.visualization import (
    plot_phi_vs_consciousness,
    plot_pci_spectrum,
    plot_gwt_metrics,
    plot_clinical_features,
    plot_multi_index_comparison,
    plot_confusion_matrix,
    plot_eeg_examples,
    plot_phi_heatmap,
)

from sklearn.metrics import confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier

import warnings
warnings.filterwarnings("ignore")


# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES_DIR = os.path.join(ROOT, "figures")
RESULTS_DIR = os.path.join(ROOT, "results")
DATA_DIR    = os.path.join(ROOT, "data")
LOGS_DIR    = os.path.join(ROOT, "logs")

for d in [FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR]:
    os.makedirs(d, exist_ok=True)

LOG_FILE = os.path.join(LOGS_DIR, "process-log.jsonl")


def log_event(phase, event_type, skill_or_tool, handoff_in=None,
              handoff_out=None, files_written=None, status="ok", **kwargs):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files_written or [],
        "status": status,
        **kwargs,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


def save_results(data, filename):
    path = os.path.join(RESULTS_DIR, filename)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=float)
    return path


# ─── Consciousness levels to sweep ───────────────────────────────────────────
LEVELS = np.linspace(0.05, 1.0, 16)

print("=" * 70)
print("NCC FRAMEWORK — Full Simulation Pipeline")
print("=" * 70)

log_event("INIT", "run_started", "ncc_framework", status="ok")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. IIT Φ CALCULATION
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[1/5] Computing IIT Φ across consciousness levels...")
log_event("IIT", "handoff_started", "PhiCalculator",
          handoff_in={"n_levels": len(LEVELS), "max_nodes": 4})

phi_calc = PhiCalculator(max_nodes=4, n_states=2)
phi_values = []
phi_stds = []
phi_matrices_by_state = {}

for i, level in enumerate(LEVELS):
    level_phi = []
    # Average over 3 seeds for stability
    for seed in [42, 43, 44]:
        data = generate_anesthesia_data(
            n_channels=8, n_samples=512,
            consciousness_level=float(level), seed=seed
        )
        p = phi_calc.compute_phi(data)
        level_phi.append(p)

    phi_values.append(np.mean(level_phi))
    phi_stds.append(np.std(level_phi))
    print(f"  Level {level:.2f} → Φ = {phi_values[-1]:.4f} ± {phi_stds[-1]:.4f}")

phi_values = np.array(phi_values)
phi_stds = np.array(phi_stds)

# Pairwise Φ matrices for representative states
print("\n  Computing pairwise Φ matrices for representative states...")
for state_label, cl_level in [("VS", 0.15), ("MCS", 0.45), ("CTRL", 0.95)]:
    data = generate_anesthesia_data(n_channels=6, n_samples=512,
                                     consciousness_level=cl_level, seed=42)
    phi_mat = phi_calc.integrated_information_matrix(data)
    phi_matrices_by_state[state_label] = phi_mat.tolist()

    plot_phi_heatmap(
        phi_mat,
        save_path=os.path.join(FIGURES_DIR, f"phi_matrix_{state_label}.png"),
        title=f"Pairwise Φ Matrix — {state_label} (level={cl_level})"
    )

# Plot
plot_phi_vs_consciousness(
    phi_values, LEVELS,
    save_path=os.path.join(FIGURES_DIR, "phi_vs_consciousness.png"),
    phi_std=phi_stds,
)

iit_results = {
    "consciousness_levels": LEVELS.tolist(),
    "phi_values": phi_values.tolist(),
    "phi_stds": phi_stds.tolist(),
    "correlation_phi_level": float(np.corrcoef(LEVELS, phi_values)[0, 1]),
    "phi_matrices": phi_matrices_by_state,
}
save_results(iit_results, "iit_phi_results.json")

log_event("IIT", "handoff_completed", "PhiCalculator",
          handoff_out={"correlation": iit_results["correlation_phi_level"]},
          files_written=["results/iit_phi_results.json",
                         "figures/phi_vs_consciousness.png"])

print(f"  Φ–consciousness correlation: r = {iit_results['correlation_phi_level']:.3f}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. PCI SIMULATION
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[2/5] Simulating Perturbational Complexity Index (PCI)...")
log_event("PCI", "handoff_started", "PCISimulator",
          handoff_in={"n_levels": len(LEVELS), "n_trials": 5})

pci_sim = PCISimulator(n_channels=20, fs=1000.0, alpha=0.01)
pci_results = pci_sim.pci_across_levels(LEVELS, n_trials=3, seed=100)

pci_values = np.array([r["pci"] for r in pci_results])

for res in pci_results:
    print(f"  Level {res['consciousness_level']:.2f} → PCI = {res['pci']:.4f} ± {res['pci_std']:.4f}")

# Plot
plot_pci_spectrum(
    pci_results,
    save_path=os.path.join(FIGURES_DIR, "pci_spectrum.png"),
)

save_results(pci_results, "pci_results.json")

pci_corr = float(np.corrcoef(LEVELS, pci_values)[0, 1])
log_event("PCI", "handoff_completed", "PCISimulator",
          handoff_out={"correlation": pci_corr},
          files_written=["results/pci_results.json", "figures/pci_spectrum.png"])

print(f"  PCI–consciousness correlation: r = {pci_corr:.3f}")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. GLOBAL WORKSPACE THEORY METRICS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[3/5] Computing Global Workspace Theory metrics...")
log_event("GWT", "handoff_started", "GlobalWorkspaceAnalyzer",
          handoff_in={"n_levels": len(LEVELS)})

gwt_analyzer = GlobalWorkspaceAnalyzer(threshold=0.25, fs=256.0)
gwt_results = []

for level in LEVELS:
    data = generate_anesthesia_data(
        n_channels=8, n_samples=512,
        consciousness_level=float(level), seed=42
    )
    metrics = gwt_analyzer.analyze(data)
    metrics.pop("connectivity_matrix", None)  # remove large arrays
    metrics.pop("adjacency_matrix", None)
    gci = gwt_analyzer.gwt_consciousness_index(data)
    metrics["gwt_index"] = gci
    metrics["consciousness_level"] = float(level)
    gwt_results.append(metrics)
    print(f"  Level {level:.2f} → GEff={metrics['global_efficiency']:.3f} "
          f"Ignition={metrics['ignition_index']:.3f} GCI={gci:.3f}")

gwt_indices = np.array([r["gwt_index"] for r in gwt_results])

plot_gwt_metrics(
    gwt_results, LEVELS,
    save_path=os.path.join(FIGURES_DIR, "gwt_metrics.png"),
)

save_results(gwt_results, "gwt_results.json")

gwt_corr = float(np.corrcoef(LEVELS, gwt_indices)[0, 1])
log_event("GWT", "handoff_completed", "GlobalWorkspaceAnalyzer",
          handoff_out={"correlation": gwt_corr},
          files_written=["results/gwt_results.json", "figures/gwt_metrics.png"])

print(f"  GWT Index–consciousness correlation: r = {gwt_corr:.3f}")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. CLINICAL CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[4/5] Training and evaluating clinical consciousness classifier...")
log_event("CLINICAL", "handoff_started", "ConsciousnessClassifier",
          handoff_in={"n_per_class": 30})

X, y, feature_names = generate_clinical_dataset(
    n_samples_per_class=30,
    n_channels=8,
    n_time=1024,
    seed=42,
)

classifier = ConsciousnessClassifier(n_channels=8, fs=256.0)
cv_results = classifier.cross_validate(X, y, cv=5)

print(f"  RF accuracy: {cv_results['rf_accuracy']:.3f} ± {cv_results['rf_accuracy_std']:.3f}")
print(f"  LDA accuracy: {cv_results['lda_accuracy']:.3f} ± {cv_results['lda_accuracy_std']:.3f}")

# Full fit for confusion matrix
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
rf.fit(X_scaled, y)
y_pred = rf.predict(X_scaled)

# Confusion matrix
labels_order = sorted(CLASS_NAMES.keys())
class_names_ordered = [CLASS_NAMES[i] for i in labels_order]
cm = confusion_matrix(y, y_pred, labels=labels_order, normalize="true")

plot_confusion_matrix(
    cm, class_names_ordered,
    save_path=os.path.join(FIGURES_DIR, "confusion_matrix.png"),
    title="Clinical Consciousness Classification (RF, Training Set)"
)

# LDA projection + feature importance
plot_clinical_features(
    X, y, feature_names, CLASS_NAMES,
    save_path=os.path.join(FIGURES_DIR, "clinical_lda_features.png"),
)

# Feature importances
feat_imp = dict(sorted(
    zip(feature_names, rf.feature_importances_.tolist()),
    key=lambda x: -x[1]
))

clf_results = {
    "cross_validation": cv_results,
    "feature_importances": feat_imp,
    "top_5_features": list(feat_imp.keys())[:5],
    "n_classes": len(CLASS_NAMES),
    "class_names": CLASS_NAMES,
}
save_results(clf_results, "classification_results.json")

print(f"  Top features: {clf_results['top_5_features']}")

log_event("CLINICAL", "handoff_completed", "ConsciousnessClassifier",
          handoff_out=cv_results,
          files_written=["results/classification_results.json"])


# ═══════════════════════════════════════════════════════════════════════════════
# 5. MULTI-THEORY COMPARISON + EEG EXAMPLES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[5/5] Generating comparative figures...")

plot_multi_index_comparison(
    LEVELS, phi_values, pci_values, gwt_indices,
    save_path=os.path.join(FIGURES_DIR, "multi_index_comparison.png"),
)

plot_eeg_examples(
    consciousness_levels=[1.0, 0.55, 0.35, 0.15],
    n_channels=4, n_samples=512, fs=256.0,
    save_path=os.path.join(FIGURES_DIR, "eeg_examples.png"),
    seed=42,
)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. AI CONSCIOUSNESS IMPLICATIONS — METRIC SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
# Compute reference values for clinical states
print("\n[+] Computing reference values per clinical state...")

ai_summary = {}
for state, cl_level in STATE_CONSCIOUSNESS.items():
    data = generate_anesthesia_data(n_channels=8, n_samples=512,
                                     consciousness_level=cl_level, seed=42)
    phi_v = phi_calc.compute_phi(data)
    gwt_v = gwt_analyzer.gwt_consciousness_index(data)

    # PCI at this level
    pci_r = pci_sim.simulate_and_compute(cl_level, n_trials=2, seed=42)

    ai_summary[state] = {
        "consciousness_level": cl_level,
        "phi": phi_v,
        "pci": pci_r["pci"],
        "gwt_index": gwt_v,
        "composite_ncc_index": (phi_v + pci_r["pci"] + gwt_v) / 3,
    }
    print(f"  {state}: Φ={phi_v:.4f} PCI={pci_r['pci']:.4f} GWT={gwt_v:.4f}")

save_results(ai_summary, "clinical_state_reference.json")


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("RESULTS SUMMARY")
print("=" * 70)
print(f"  IIT Φ–level correlation:        r = {iit_results['correlation_phi_level']:.3f}")
print(f"  PCI–level correlation:           r = {pci_corr:.3f}")
print(f"  GWT Index–level correlation:     r = {gwt_corr:.3f}")
print(f"  Clinical RF 5-fold accuracy:     {cv_results['rf_accuracy']:.3f} ± {cv_results['rf_accuracy_std']:.3f}")
print(f"  Clinical LDA 5-fold accuracy:    {cv_results['lda_accuracy']:.3f} ± {cv_results['lda_accuracy_std']:.3f}")

# Save composite summary
summary = {
    "iit": {"phi_level_correlation": iit_results["correlation_phi_level"],
            "max_phi": float(phi_values.max()), "min_phi": float(phi_values.min())},
    "pci": {"pci_level_correlation": pci_corr,
            "max_pci": float(pci_values.max()), "min_pci": float(pci_values.min())},
    "gwt": {"gwt_level_correlation": gwt_corr,
            "max_gwt": float(gwt_indices.max()), "min_gwt": float(gwt_indices.min())},
    "clinical_classification": cv_results,
    "clinical_state_reference": ai_summary,
    "figures_generated": [
        "phi_vs_consciousness.png",
        "pci_spectrum.png",
        "gwt_metrics.png",
        "clinical_lda_features.png",
        "confusion_matrix.png",
        "multi_index_comparison.png",
        "eeg_examples.png",
        "phi_matrix_VS.png",
        "phi_matrix_MCS.png",
        "phi_matrix_CTRL.png",
    ]
}
save_results(summary, "full_summary.json")

log_event("REPORT", "run_completed", "ncc_framework",
          handoff_out=summary,
          files_written=["results/full_summary.json"],
          status="ok")

print("\n[DONE] All results saved to results/ and figures/")
print("=" * 70)
