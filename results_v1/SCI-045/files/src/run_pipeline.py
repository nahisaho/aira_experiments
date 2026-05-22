#!/usr/bin/env python3
"""
Main Pipeline: Epigenetic Clock Development and Evaluation
Runs all analyses and generates figures and results.
"""
import sys, json, time
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from data_simulator import generate_all_datasets
from traditional_clocks import (
    HorvathStyleClock, GrimAgeStyleClock, ImprovedElasticNetClock,
    evaluate_clock, cross_validate_clock
)
from deep_clock import DeepClockTrainer
from analysis import (
    age_acceleration_analysis, tissue_specificity_analysis,
    intervention_sensitivity, longevity_validation, compute_age_acceleration
)
from visualization import (
    plot_prediction_scatter, plot_model_comparison, plot_tissue_performance,
    plot_age_acceleration, plot_intervention_effects, plot_training_history,
    plot_longevity_comparison
)

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)
DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)
LOGS_DIR = Path("logs")
LOGS_DIR.mkdir(exist_ok=True)

all_results = {}

def log_event(phase, event_type, **kwargs):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "phase": phase, "event_type": event_type,
        "actor": "co-scientist", **kwargs
    }
    with open(LOGS_DIR / "process-log.jsonl", "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")

# ── Phase 1: Data Generation ──
print("=" * 60)
print("PHASE 1: Data Generation")
print("=" * 60)
log_event("data_generation", "run_started")

datasets = generate_all_datasets()
for name, df in datasets.items():
    df.to_csv(DATA_DIR / f"{name}.csv", index=False)
    print(f"  {name}: {df.shape}")

log_event("data_generation", "run_completed", files_written=[
    f"data/{n}.csv" for n in datasets.keys()])

# ── Phase 2: Traditional Clock Training & Evaluation ──
print("\n" + "=" * 60)
print("PHASE 2: Traditional Clock Models")
print("=" * 60)
log_event("traditional_clocks", "run_started")

blood = datasets["blood_train"]
cpg_cols = [c for c in blood.columns if c.startswith("cg")]
X_blood = blood[cpg_cols].values
y_blood = blood["true_bio_age"].values

n_train = int(len(blood) * 0.8)
idx = np.random.permutation(len(blood))
train_idx, test_idx = idx[:n_train], idx[n_train:]

# 2a. Horvath-style clock
print("\n--- Horvath-Style Clock ---")
horvath = HorvathStyleClock(alpha=0.1, l1_ratio=0.5)
horvath.fit(X_blood[train_idx], y_blood[train_idx])
horvath_pred = horvath.predict(X_blood[test_idx])
horvath_metrics = evaluate_clock(y_blood[test_idx], horvath_pred, "Horvath-Style")
print(f"  MAE={horvath_metrics['MAE']}, R²={horvath_metrics['R2']}, "
      f"Non-zero CpGs={horvath.n_nonzero_}")
plot_prediction_scatter(y_blood[test_idx], horvath_pred,
                        "Horvath-Style Clock", "horvath_scatter.png")

# 2b. GrimAge-style clock
print("\n--- GrimAge-Style Clock ---")
grim = GrimAgeStyleClock()
grim.fit(X_blood[train_idx], y_blood[train_idx])
grim_pred = grim.predict(X_blood[test_idx])
grim_metrics = evaluate_clock(y_blood[test_idx], grim_pred, "GrimAge-Style")
print(f"  MAE={grim_metrics['MAE']}, R²={grim_metrics['R2']}")
plot_prediction_scatter(y_blood[test_idx], grim_pred,
                        "GrimAge-Style Clock", "grimage_scatter.png")

# 2c. Improved ElasticNet clock
print("\n--- Improved ElasticNet Clock ---")
improved = ImprovedElasticNetClock()
improved.fit(X_blood[train_idx], y_blood[train_idx])
improved_pred = improved.predict(X_blood[test_idx])
improved_metrics = evaluate_clock(y_blood[test_idx], improved_pred, "Improved-ElasticNet")
print(f"  MAE={improved_metrics['MAE']}, R²={improved_metrics['R2']}, "
      f"Non-zero={improved.n_nonzero_}")
plot_prediction_scatter(y_blood[test_idx], improved_pred,
                        "Improved ElasticNet Clock", "improved_elasticnet_scatter.png")

# Cross-validation for Horvath
print("\n--- 5-Fold Cross-Validation (Horvath) ---")
cv_preds, cv_overall, cv_folds = cross_validate_clock(
    HorvathStyleClock, X_blood, y_blood, n_folds=5)
print(f"  CV MAE={cv_overall['MAE']}, R²={cv_overall['R2']}")

all_results["traditional_clocks"] = {
    "horvath": horvath_metrics,
    "grimage": grim_metrics,
    "improved_elasticnet": improved_metrics,
    "horvath_cv": cv_overall,
}

log_event("traditional_clocks", "run_completed")

# ── Phase 3: Deep Learning Clock ──
print("\n" + "=" * 60)
print("PHASE 3: Deep Learning Clock")
print("=" * 60)
log_event("deep_clock", "run_started")

# Prepare multi-tissue training data
train_tissues = pd.concat([
    datasets["blood_train"].iloc[train_idx],
    datasets["brain"].iloc[:200],
    datasets["liver"].iloc[:200],
], ignore_index=True)

val_tissues = pd.concat([
    datasets["blood_train"].iloc[test_idx],
    datasets["brain"].iloc[200:],
    datasets["liver"].iloc[200:],
], ignore_index=True)

n_cpg = len(cpg_cols)
trainer = DeepClockTrainer(n_cpg=n_cpg, n_tissues=5, epochs=80, patience=12,
                           lr=5e-4, batch_size=64)
print("  Training deep clock...")
trainer.fit(train_tissues, val_tissues)

# Predict on blood test set
blood_test = datasets["blood_train"].iloc[test_idx].copy()
deep_pred = trainer.predict(blood_test)
deep_metrics = evaluate_clock(blood_test["true_bio_age"].values, deep_pred, "Deep-Clock")
print(f"  MAE={deep_metrics['MAE']}, R²={deep_metrics['R2']}")

plot_prediction_scatter(blood_test["true_bio_age"].values, deep_pred,
                        "Deep Learning Clock", "deep_clock_scatter.png")
plot_training_history(trainer.history)

all_results["deep_clock"] = deep_metrics

log_event("deep_clock", "run_completed")

# ── Phase 4: Model Comparison ──
print("\n" + "=" * 60)
print("PHASE 4: Model Comparison")
print("=" * 60)

all_metrics = [horvath_metrics, grim_metrics, improved_metrics, deep_metrics]
plot_model_comparison(all_metrics)

comparison_df = pd.DataFrame(all_metrics)
comparison_df.to_csv(RESULTS_DIR / "model_comparison.csv", index=False)
print(comparison_df.to_string(index=False))

# ── Phase 5: Tissue Specificity ──
print("\n" + "=" * 60)
print("PHASE 5: Tissue-Specific Analysis")
print("=" * 60)
log_event("tissue_analysis", "run_started")

tissue_preds = {}
for tissue in ["blood", "brain", "liver", "skin", "muscle"]:
    df_tissue = datasets[tissue]
    X_t = df_tissue[cpg_cols].values
    y_t = df_tissue["true_bio_age"].values
    # Use Horvath clock trained on blood
    pred = horvath.predict(X_t)
    tissue_preds[tissue] = {"true": y_t, "pred": pred}

tissue_results = tissue_specificity_analysis(tissue_preds)
for t, r in tissue_results.items():
    print(f"  {t}: MAE={r['MAE']}, R²={r['R2']}")

plot_tissue_performance(tissue_results)
all_results["tissue_specificity"] = tissue_results

pd.DataFrame(tissue_results).T.to_csv(RESULTS_DIR / "tissue_performance.csv")
log_event("tissue_analysis", "run_completed")

# ── Phase 6: Age Acceleration Analysis ──
print("\n" + "=" * 60)
print("PHASE 6: Age Acceleration Biomarker")
print("=" * 60)
log_event("age_acceleration", "run_started")

blood_full = datasets["blood_train"].copy()
blood_full["predicted_age"] = horvath.predict(X_blood)
blood_aa, aa_results = age_acceleration_analysis(blood_full)
print(f"  Mean acceleration: {aa_results['mean_acceleration']}")
print(f"  Correlation with true offset: {aa_results.get('correlation_with_true_offset', 'N/A')}")

plot_age_acceleration(
    blood_aa["age_acceleration"].values,
    groups=blood_aa["sex"].map({0: "Female", 1: "Male"}),
    filename="age_acceleration_blood.png"
)

all_results["age_acceleration"] = aa_results
log_event("age_acceleration", "run_completed")

# ── Phase 7: Intervention Sensitivity ──
print("\n" + "=" * 60)
print("PHASE 7: Intervention Sensitivity")
print("=" * 60)
log_event("intervention", "run_started")

intv_df = datasets["intervention"].copy()
X_intv = intv_df[cpg_cols].values
intv_df["predicted_age"] = horvath.predict(X_intv)
intv_results = intervention_sensitivity(intv_df)

for k, v in intv_results.items():
    if isinstance(v, dict):
        print(f"  {k}: Δ={v['delta_vs_control']}, d={v['cohens_d']}, "
              f"p={v['p_value']}, detect={v['detectable']}")

plot_intervention_effects(intv_results)
all_results["intervention_sensitivity"] = intv_results

log_event("intervention", "run_completed")

# ── Phase 8: Longevity Validation ──
print("\n" + "=" * 60)
print("PHASE 8: Longevity Cohort Validation")
print("=" * 60)
log_event("longevity", "run_started")

long_df = datasets["longevity"].copy()
X_long = long_df[cpg_cols].values
long_df["predicted_age"] = horvath.predict(X_long)

normal_df = datasets["blood"].copy()
X_norm = normal_df[cpg_cols].values
normal_df["predicted_age"] = horvath.predict(X_norm)

long_results = longevity_validation(long_df, normal_df)
print(f"  Longevity accel: {long_results['longevity_mean_accel']}")
print(f"  Normal accel: {long_results['normal_mean_accel']}")
print(f"  Cohen's d: {long_results['cohens_d']}, p={long_results['p_value']}")

accel_long = compute_age_acceleration(
    long_df["chronological_age"].values, long_df["predicted_age"].values)
accel_norm = compute_age_acceleration(
    normal_df["chronological_age"].values, normal_df["predicted_age"].values)
plot_longevity_comparison(accel_long, accel_norm)

all_results["longevity_validation"] = long_results
log_event("longevity", "run_completed")

# ── Save All Results ──
print("\n" + "=" * 60)
print("SAVING FINAL RESULTS")
print("=" * 60)

with open(RESULTS_DIR / "all_results.json", "w") as f:
    json.dump(all_results, f, indent=2, ensure_ascii=False, default=str)
print(f"  Results saved to {RESULTS_DIR / 'all_results.json'}")

log_event("pipeline", "run_completed", status="ok",
          files_written=["results/all_results.json", "results/model_comparison.csv",
                         "results/tissue_performance.csv"])

print("\n✅ Pipeline complete!")
