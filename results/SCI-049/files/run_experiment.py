"""
Main Experiment: End-to-end demonstration of the anomaly detection pipeline.
Generates synthetic CERN/LIGO-like data and runs all 6 modules.
"""
import numpy as np
import pandas as pd
import json
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.changepoint_detection import PELTDetector, BOCPDDetector
from src.multivariate_outlier import IsolationForestDetector, DeepSVDDDetector
from src.physics_constraints import (
    PhysicsConstrainedScorer, energy_conservation_constraint,
    range_constraint, positive_definite_constraint
)
from src.drift_detection import ADWINDetector, PageHinkleyDetector, RetrainingTrigger
from src.explainable_anomaly import ExplainableAnomalyDetector
from src.streaming_pipeline import (
    StreamingAnomalyPipeline, mahalanobis_detector, zscore_detector,
    ARCHITECTURE_SPEC
)

np.random.seed(42)

# ═══════════════════════════════════════════
# Data Generation: Synthetic CERN-like events
# ═══════════════════════════════════════════

def generate_experiment_data(n_samples=5000, n_features=8, anomaly_fraction=0.05):
    """Generate synthetic particle physics-like data with known anomalies."""
    n_normal = int(n_samples * (1 - anomaly_fraction))
    n_anomaly = n_samples - n_normal

    feature_names = [
        "energy", "transverse_momentum", "pseudorapidity", "azimuthal_angle",
        "invariant_mass", "missing_ET", "jet_multiplicity", "track_isolation"
    ][:n_features]

    # Normal events (Standard Model background)
    normal = np.column_stack([
        np.random.exponential(50, n_normal),      # energy [GeV]
        np.random.exponential(20, n_normal),      # pT [GeV]
        np.random.normal(0, 2, n_normal),         # eta
        np.random.uniform(-np.pi, np.pi, n_normal),  # phi
        np.random.exponential(80, n_normal),      # invariant mass [GeV]
        np.random.exponential(10, n_normal),      # missing ET [GeV]
        np.random.poisson(3, n_normal).astype(float),  # jet count
        np.random.exponential(0.5, n_normal),     # isolation
    ])[:, :n_features]

    # Anomalous events (BSM physics signals + detector artifacts)
    anomalies = np.column_stack([
        np.random.exponential(200, n_anomaly),    # higher energy
        np.random.exponential(80, n_anomaly),     # higher pT
        np.random.normal(0, 0.5, n_anomaly),     # more central
        np.random.uniform(-np.pi, np.pi, n_anomaly),
        np.random.normal(125, 5, n_anomaly),      # Higgs-like mass peak
        np.random.exponential(50, n_anomaly),     # higher missing ET
        np.random.poisson(6, n_anomaly).astype(float),
        np.random.exponential(2, n_anomaly),
    ])[:, :n_features]

    X = np.vstack([normal, anomalies])
    y_true = np.array([0] * n_normal + [1] * n_anomaly)

    # Shuffle
    idx = np.random.permutation(n_samples)
    X, y_true = X[idx], y_true[idx]

    return X, y_true, feature_names


def generate_timeseries_with_changepoints(n=2000):
    """Generate time series with known changepoints and drift."""
    segments = [
        np.random.normal(10, 1, 500),
        np.random.normal(15, 1.5, 400),
        np.random.normal(10, 2, 300),
        np.random.normal(20, 1, 400),
        np.random.normal(12, 1, 400),
    ]
    true_cps = [500, 900, 1200, 1600]
    ts = np.concatenate(segments)[:n]
    # Add trend/drift in last segment
    drift = np.linspace(0, 5, 400)
    ts[1600:2000] += drift[:min(400, n-1600)]
    return ts, true_cps


# ═══════════════════════════════════════════
# Run All Experiments
# ═══════════════════════════════════════════

def run_experiments():
    results = {}
    log_entries = []

    def log(phase, event, **kwargs):
        entry = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "phase": phase,
            "event_type": event,
            "actor": "co-scientist",
            **kwargs
        }
        log_entries.append(entry)
        print(f"  [{phase}] {event}")

    print("=" * 70)
    print("Large-Scale Scientific Data QC & Anomaly Detection System")
    print("=" * 70)

    # ── Generate Data ──
    log("data", "generating_data")
    X, y_true, feat_names = generate_experiment_data(5000, 8)
    ts, true_cps = generate_timeseries_with_changepoints(2000)

    pd.DataFrame(X, columns=feat_names).to_csv("data/synthetic_events.csv", index=False)
    pd.DataFrame({"value": ts}).to_csv("data/timeseries.csv", index=False)
    log("data", "data_saved", files=["data/synthetic_events.csv", "data/timeseries.csv"])

    # ═══════════════════════════════════════
    # Experiment 1: Change Point Detection
    # ═══════════════════════════════════════
    print("\n── Experiment 1: Change Point Detection ──")

    # PELT
    log("changepoint", "pelt_start")
    pelt = PELTDetector(model="rbf", min_size=30)
    pelt_result = pelt.fit_predict(ts)
    detected_cps = [bp for bp in pelt_result["breakpoints"] if bp < len(ts)]
    print(f"  PELT: Detected {pelt_result['n_changepoints']} changepoints at {detected_cps}")
    results["pelt"] = {
        "detected_changepoints": detected_cps,
        "true_changepoints": true_cps,
        "n_detected": pelt_result["n_changepoints"],
        "n_true": len(true_cps),
    }
    log("changepoint", "pelt_complete")

    # BOCPD
    log("changepoint", "bocpd_start")
    bocpd = BOCPDDetector(hazard_rate=250)
    bocpd_result = bocpd.detect(ts[:1000], threshold=0.3)
    print(f"  BOCPD: Detected {bocpd_result['n_changepoints']} changepoints")
    results["bocpd"] = {
        "detected_changepoints": bocpd_result["changepoints"][:20],
        "n_detected": bocpd_result["n_changepoints"],
    }
    log("changepoint", "bocpd_complete")

    # ═══════════════════════════════════════
    # Experiment 2: Multivariate Outlier Detection
    # ═══════════════════════════════════════
    print("\n── Experiment 2: Multivariate Outlier Detection ──")

    # Isolation Forest
    log("outlier", "iforest_start")
    iforest = IsolationForestDetector(contamination=0.05, n_estimators=200)
    iforest_result = iforest.fit_predict(X, feat_names)
    iforest_labels = iforest_result["labels"]
    tp = int(np.sum((iforest_labels == -1) & (y_true == 1)))
    fp = int(np.sum((iforest_labels == -1) & (y_true == 0)))
    fn = int(np.sum((iforest_labels == 1) & (y_true == 1)))
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"  Isolation Forest: {iforest_result['n_anomalies']} anomalies detected")
    print(f"    Precision={precision:.3f}, Recall={recall:.3f}, F1={f1:.3f}")
    results["isolation_forest"] = {
        "n_anomalies": iforest_result["n_anomalies"],
        "anomaly_rate": iforest_result["anomaly_rate"],
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "feature_importance": iforest_result["feature_importance"],
        "score_stats": iforest_result["score_stats"],
    }
    log("outlier", "iforest_complete")

    # Deep SVDD
    log("outlier", "deep_svdd_start")
    dsvdd = DeepSVDDDetector(encoding_dim=4, nu=0.05)
    dsvdd.fit(X[y_true == 0][:3000], n_epochs=50, lr=0.005)
    dsvdd_result = dsvdd.predict(X)
    dsvdd_labels = dsvdd_result["labels"]
    tp_d = int(np.sum((dsvdd_labels == -1) & (y_true == 1)))
    fp_d = int(np.sum((dsvdd_labels == -1) & (y_true == 0)))
    fn_d = int(np.sum((dsvdd_labels == 1) & (y_true == 1)))
    prec_d = tp_d / (tp_d + fp_d) if (tp_d + fp_d) > 0 else 0
    rec_d = tp_d / (tp_d + fn_d) if (tp_d + fn_d) > 0 else 0
    f1_d = 2 * prec_d * rec_d / (prec_d + rec_d) if (prec_d + rec_d) > 0 else 0

    print(f"  Deep SVDD: {dsvdd_result['n_anomalies']} anomalies, "
          f"Precision={prec_d:.3f}, Recall={rec_d:.3f}, F1={f1_d:.3f}")
    results["deep_svdd"] = {
        "n_anomalies": dsvdd_result["n_anomalies"],
        "precision": round(prec_d, 4),
        "recall": round(rec_d, 4),
        "f1": round(f1_d, 4),
        "radius": dsvdd_result["radius"],
    }
    log("outlier", "deep_svdd_complete")

    # ═══════════════════════════════════════
    # Experiment 3: Physics-Constrained Scoring
    # ═══════════════════════════════════════
    print("\n── Experiment 3: Physics-Constrained Anomaly Scoring ──")
    log("physics", "scoring_start")

    scorer = PhysicsConstrainedScorer(stat_weight=0.5, phys_weight=0.5)
    scorer.add_constraint("energy_positive", lambda d: np.clip(-d[:, 0], 0, None),
                          weight=2.0, description="Energy > 0")
    scorer.add_constraint("pT_range", lambda d: np.clip(d[:, 1] - 500, 0, None),
                          weight=1.5, description="pT < 500 GeV")
    scorer.add_constraint("mass_positive", lambda d: np.clip(-d[:, 4], 0, None),
                          weight=2.0, description="Invariant mass > 0")
    scorer.add_constraint("eta_range", lambda d: np.clip(np.abs(d[:, 2]) - 5, 0, None),
                          weight=1.0, description="|η| < 5")

    stat_scores = -iforest_result["scores"]  # higher = more anomalous
    phys_result = scorer.score(X, stat_scores)

    n_phys_violations = sum(c["n_violations"] for c in phys_result["constraint_details"])
    print(f"  Physics violations: {n_phys_violations} total across {len(phys_result['constraint_details'])} constraints")
    for c in phys_result["constraint_details"]:
        print(f"    {c['name']}: {c['n_violations']} violations ({c['violation_rate']:.4f})")

    results["physics_constrained"] = {
        "constraint_details": phys_result["constraint_details"],
        "total_violations": n_phys_violations,
        "score_correlation": float(np.corrcoef(
            phys_result["statistical_scores"], phys_result["physics_scores"]
        )[0, 1]) if np.std(phys_result["physics_scores"]) > 0 else 0.0,
    }
    log("physics", "scoring_complete")

    # ═══════════════════════════════════════
    # Experiment 4: Drift Detection
    # ═══════════════════════════════════════
    print("\n── Experiment 4: Concept Drift Detection ──")
    log("drift", "detection_start")

    # Generate drifting data
    drift_data = np.concatenate([
        np.random.normal(0, 1, 500),
        np.random.normal(0.5, 1, 300),
        np.random.normal(0, 1, 400),
        np.random.normal(2.0, 1.5, 300),
        np.random.normal(0, 1, 500),
    ])

    adwin = ADWINDetector(delta=0.01)
    adwin_result = adwin.detect_batch(drift_data)
    print(f"  ADWIN: {adwin_result['n_drifts']} drifts detected")

    ph = PageHinkleyDetector(delta=0.01, threshold=30)
    ph_result = ph.detect_batch(drift_data)
    print(f"  Page-Hinkley: {ph_result['n_drifts']} drifts detected")

    # Retraining trigger simulation
    trigger = RetrainingTrigger(performance_threshold=0.05, drift_patience=3)
    retrain_events = []
    baseline_perf = 0.95
    for step in range(0, len(drift_data), 50):
        current_perf = baseline_perf - 0.02 * np.sin(step / 200) - (step > 1200) * 0.1
        decision = trigger.evaluate(step, current_perf, baseline_perf,
                                    step in adwin_result.get("drift_points", []))
        if decision.should_retrain:
            retrain_events.append({
                "step": step, "reason": decision.reason,
                "type": decision.drift_type, "severity": round(decision.severity, 4)
            })
    print(f"  Retrain triggers: {len(retrain_events)}")

    results["drift_detection"] = {
        "adwin_drifts": adwin_result["n_drifts"],
        "page_hinkley_drifts": ph_result["n_drifts"],
        "retrain_events": retrain_events,
        "n_retrain_triggers": len(retrain_events),
    }
    log("drift", "detection_complete")

    # ═══════════════════════════════════════
    # Experiment 5: Explainable Anomaly Detection
    # ═══════════════════════════════════════
    print("\n── Experiment 5: Explainable Anomaly Detection ──")
    log("explain", "analysis_start")

    explainer = ExplainableAnomalyDetector(contamination=0.05)
    explain_result = explainer.fit_predict_explain(X, feat_names)

    print(f"  Anomalies with explanations: {explain_result['n_anomalies']}")
    print(f"  Top global features: ", end="")
    top3 = list(explain_result["global_feature_importance"].items())[:3]
    print(", ".join(f"{k}={v:.3f}" for k, v in top3))
    print(f"  Decision rules extracted: {len(explain_result['decision_rules'])}")
    print(f"  Root cause clusters: {len(explain_result['root_causes'])}")

    if explain_result["root_causes"]:
        for rc in explain_result["root_causes"][:3]:
            print(f"    {rc['feature']}: {rc['count']} anomalies "
                  f"({rc['fraction']:.1%}), z={rc['mean_z_score']:.2f} ({rc['direction']})")

    results["explainable"] = {
        "n_anomalies": explain_result["n_anomalies"],
        "global_feature_importance": explain_result["global_feature_importance"],
        "n_rules": len(explain_result["decision_rules"]),
        "top_rules": explain_result["decision_rules"][:3],
        "root_causes": explain_result["root_causes"][:5],
        "n_local_explanations": len(explain_result["local_explanations"]),
    }
    log("explain", "analysis_complete")

    # ═══════════════════════════════════════
    # Experiment 6: Streaming Pipeline
    # ═══════════════════════════════════════
    print("\n── Experiment 6: Streaming Pipeline (CERN/LIGO-scale) ──")
    log("streaming", "pipeline_start")

    pipeline = StreamingAnomalyPipeline(
        window_size=500, alert_threshold=0.7, n_features=8
    )
    pipeline.add_detector(mahalanobis_detector)
    pipeline.add_detector(zscore_detector)

    # Simulate streaming
    stream_results = pipeline.process_batch(X[:2000])
    summary = pipeline.get_summary()

    print(f"  Processed: {stream_results['n_processed']} events")
    print(f"  Anomalies: {stream_results['n_anomalies']} ({stream_results['anomaly_rate']:.3f})")
    print(f"  Throughput: {stream_results['metrics']['throughput_per_sec']:.0f} events/sec")
    print(f"  Alerts: {summary['n_alerts']} "
          f"(Critical={summary['alert_severity_counts']['critical']}, "
          f"Warning={summary['alert_severity_counts']['warning']})")

    results["streaming"] = {
        "events_processed": stream_results["n_processed"],
        "anomalies_detected": stream_results["n_anomalies"],
        "anomaly_rate": round(stream_results["anomaly_rate"], 4),
        "mean_score": round(stream_results["mean_score"], 4),
        "throughput_per_sec": round(stream_results["metrics"]["throughput_per_sec"], 1),
        "alert_summary": summary["alert_severity_counts"],
        "architecture": ARCHITECTURE_SPEC,
    }
    log("streaming", "pipeline_complete")

    # ── Save Results ──
    print("\n── Saving Results ──")

    # Make results JSON serializable
    def make_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, dict):
            return {k: make_serializable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [make_serializable(v) for v in obj]
        return obj

    results_serial = make_serializable(results)
    with open("results/experiment_results.json", "w") as f:
        json.dump(results_serial, f, indent=2)

    with open("logs/process-log.jsonl", "w") as f:
        for entry in log_entries:
            f.write(json.dumps(entry) + "\n")

    print("  Saved: results/experiment_results.json")
    print("  Saved: logs/process-log.jsonl")

    return results, X, y_true, feat_names, ts, true_cps, iforest_result, \
           explain_result, phys_result, drift_data, adwin_result, bocpd


if __name__ == "__main__":
    run_experiments()
