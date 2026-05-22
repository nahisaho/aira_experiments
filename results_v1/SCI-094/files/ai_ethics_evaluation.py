#!/usr/bin/env python3
"""
AI Ethics Quantitative Evaluation Framework
============================================
Integrates fairness, explainability, privacy, robustness,
and environmental impact metrics into a unified scoring pipeline.
"""

import json
import os
import warnings
import time
from datetime import datetime, timezone, timedelta

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.calibration import calibration_curve
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
np.random.seed(42)

JST = timezone(timedelta(hours=9))

# ============================================================
# 0. Synthetic Medical-AI Dataset
# ============================================================

def generate_medical_dataset(n=3000):
    """Generate synthetic medical diagnosis dataset with protected attributes."""
    X, y = make_classification(
        n_samples=n, n_features=20, n_informative=12,
        n_redundant=4, n_classes=2, weights=[0.6, 0.4],
        flip_y=0.05, random_state=42,
    )
    df = pd.DataFrame(X, columns=[f"feature_{i}" for i in range(20)])
    # Protected attributes
    df["age_group"] = np.random.choice([0, 1], size=n, p=[0.45, 0.55])  # 0=young, 1=old
    df["sex"] = np.random.choice([0, 1], size=n, p=[0.5, 0.5])
    # Inject bias: shift predictions for age_group=1
    y_biased = y.copy()
    bias_mask = (df["age_group"] == 1) & (y == 0)
    flip_idx = np.where(bias_mask)[0][:int(bias_mask.sum() * 0.15)]
    y_biased[flip_idx] = 1
    df["diagnosis"] = y_biased
    df.to_csv("data/medical_dataset.csv", index=False)
    return df


# ============================================================
# 1. Fairness Evaluation Module
# ============================================================

class FairnessEvaluator:
    """Evaluate Statistical Parity, Equalized Odds, and Calibration across groups."""

    def __init__(self, y_true, y_pred, y_prob, sensitive_attr, group_names=None):
        self.y_true = np.array(y_true)
        self.y_pred = np.array(y_pred)
        self.y_prob = np.array(y_prob)
        self.sensitive = np.array(sensitive_attr)
        self.groups = np.unique(self.sensitive)
        self.group_names = group_names or {g: f"Group_{g}" for g in self.groups}

    def statistical_parity_difference(self):
        rates = {}
        for g in self.groups:
            mask = self.sensitive == g
            rates[g] = self.y_pred[mask].mean()
        vals = list(rates.values())
        return {"group_positive_rates": {self.group_names[g]: float(rates[g]) for g in self.groups},
                "spd": float(max(vals) - min(vals))}

    def equalized_odds_difference(self):
        results = {}
        for label in [0, 1]:
            rates = {}
            for g in self.groups:
                mask = (self.sensitive == g) & (self.y_true == label)
                if mask.sum() > 0:
                    rates[g] = self.y_pred[mask].mean()
                else:
                    rates[g] = 0.0
            vals = list(rates.values())
            results[f"label_{label}"] = {
                "group_rates": {self.group_names[g]: float(rates[g]) for g in self.groups},
                "difference": float(max(vals) - min(vals)),
            }
        eo_diff = max(results["label_0"]["difference"], results["label_1"]["difference"])
        results["equalized_odds_difference"] = float(eo_diff)
        return results

    def calibration_difference(self, n_bins=10):
        cal_errors = {}
        for g in self.groups:
            mask = self.sensitive == g
            if mask.sum() < n_bins:
                cal_errors[g] = float("nan")
                continue
            prob_true, prob_pred = calibration_curve(
                self.y_true[mask], self.y_prob[mask], n_bins=n_bins, strategy="uniform"
            )
            cal_errors[g] = float(np.mean(np.abs(prob_true - prob_pred)))
        return {"group_calibration_errors": {self.group_names[g]: cal_errors[g] for g in self.groups},
                "max_calibration_gap": float(max(cal_errors.values()) - min(cal_errors.values()))}

    def integrated_fairness_score(self):
        spd = self.statistical_parity_difference()["spd"]
        eod = self.equalized_odds_difference()["equalized_odds_difference"]
        cal = self.calibration_difference()["max_calibration_gap"]
        # Weighted composite (lower is fairer)
        score = 1.0 - (0.4 * min(spd, 1) + 0.35 * min(eod, 1) + 0.25 * min(cal, 1))
        return {"spd": spd, "eod": eod, "calibration_gap": cal,
                "integrated_fairness_score": float(np.clip(score, 0, 1))}

    def evaluate_all(self):
        return {
            "statistical_parity": self.statistical_parity_difference(),
            "equalized_odds": self.equalized_odds_difference(),
            "calibration": self.calibration_difference(),
            "integrated": self.integrated_fairness_score(),
        }


# ============================================================
# 2. Explainability Evaluation Module
# ============================================================

class ExplainabilityEvaluator:
    """Quantify SHAP consistency and explanation stability."""

    def __init__(self, model, X, feature_names=None):
        self.model = model
        self.X = np.array(X)
        self.feature_names = feature_names or [f"f{i}" for i in range(self.X.shape[1])]

    def _compute_shap(self, X_subset):
        import shap
        explainer = shap.TreeExplainer(self.model)
        shap_values = explainer.shap_values(X_subset)
        if isinstance(shap_values, list):
            shap_values = shap_values[1]
        return shap_values

    def shap_consistency(self, n_runs=5, sample_size=200):
        """Measure rank correlation of feature importance across bootstrap samples."""
        from scipy.stats import spearmanr
        rankings = []
        for i in range(n_runs):
            idx = np.random.choice(len(self.X), size=min(sample_size, len(self.X)), replace=True)
            sv = self._compute_shap(self.X[idx])
            importance = np.abs(sv).mean(axis=0)
            rank = np.argsort(-importance)
            rankings.append(rank)
        correlations = []
        for i in range(len(rankings)):
            for j in range(i + 1, len(rankings)):
                corr, _ = spearmanr(rankings[i], rankings[j])
                correlations.append(corr)
        return {"mean_rank_correlation": float(np.mean(correlations)),
                "std_rank_correlation": float(np.std(correlations)),
                "n_runs": n_runs}

    def explanation_stability(self, n_perturbations=50, noise_scale=0.01, sample_size=100):
        """Measure how stable SHAP explanations are under small input perturbations."""
        idx = np.random.choice(len(self.X), size=min(sample_size, len(self.X)), replace=False)
        X_sample = self.X[idx]
        base_shap = self._compute_shap(X_sample)

        stability_scores = []
        for _ in range(n_perturbations):
            noise = np.random.normal(0, noise_scale, X_sample.shape)
            X_perturbed = X_sample + noise
            perturbed_shap = self._compute_shap(X_perturbed)
            cosine_sims = []
            for k in range(len(X_sample)):
                a, b = base_shap[k], perturbed_shap[k]
                norm_prod = np.linalg.norm(a) * np.linalg.norm(b)
                if norm_prod > 0:
                    cosine_sims.append(float(np.dot(a, b) / norm_prod))
            stability_scores.append(np.mean(cosine_sims))

        return {"mean_stability": float(np.mean(stability_scores)),
                "std_stability": float(np.std(stability_scores)),
                "noise_scale": noise_scale}

    def explainability_score(self):
        consistency = self.shap_consistency()
        stability = self.explanation_stability()
        score = 0.5 * max(0, consistency["mean_rank_correlation"]) + 0.5 * max(0, stability["mean_stability"])
        return {"shap_consistency": consistency,
                "explanation_stability": stability,
                "integrated_explainability_score": float(np.clip(score, 0, 1))}


# ============================================================
# 3. Privacy Risk Evaluation Module
# ============================================================

class PrivacyEvaluator:
    """Evaluate privacy risk via membership inference attack simulation."""

    def __init__(self, model, X_train, y_train, X_test, y_test):
        self.model = model
        self.X_train = np.array(X_train)
        self.y_train = np.array(y_train)
        self.X_test = np.array(X_test)
        self.y_test = np.array(y_test)

    def membership_inference_attack(self, n_shadow=500):
        """Simulate a threshold-based membership inference attack."""
        # Get prediction confidence for train (members) and test (non-members)
        if hasattr(self.model, "predict_proba"):
            train_conf = np.max(self.model.predict_proba(self.X_train[:n_shadow]), axis=1)
            test_conf = np.max(self.model.predict_proba(self.X_test[:n_shadow]), axis=1)
        else:
            train_conf = np.abs(self.model.decision_function(self.X_train[:n_shadow]))
            test_conf = np.abs(self.model.decision_function(self.X_test[:n_shadow]))

        # Labels: 1=member, 0=non-member
        labels = np.concatenate([np.ones(len(train_conf)), np.zeros(len(test_conf))])
        confs = np.concatenate([train_conf, test_conf])

        # Threshold-based attack: predict member if confidence > threshold
        best_acc = 0.5
        best_threshold = 0.5
        for threshold in np.linspace(confs.min(), confs.max(), 100):
            preds = (confs >= threshold).astype(int)
            acc = accuracy_score(labels, preds)
            if acc > best_acc:
                best_acc = acc
                best_threshold = threshold

        attack_auc = roc_auc_score(labels, confs)
        # Privacy score: higher is better (more resistant)
        privacy_score = 1.0 - (attack_auc - 0.5) * 2  # Map AUC 0.5-1.0 → score 1.0-0.0
        privacy_score = float(np.clip(privacy_score, 0, 1))

        return {
            "attack_accuracy": float(best_acc),
            "attack_auc": float(attack_auc),
            "best_threshold": float(best_threshold),
            "membership_inference_resistance": privacy_score,
            "risk_level": "LOW" if privacy_score > 0.7 else ("MEDIUM" if privacy_score > 0.4 else "HIGH"),
        }

    def overfitting_gap(self):
        train_acc = accuracy_score(self.y_train, self.model.predict(self.X_train))
        test_acc = accuracy_score(self.y_test, self.model.predict(self.X_test))
        return {"train_accuracy": float(train_acc),
                "test_accuracy": float(test_acc),
                "gap": float(train_acc - test_acc)}

    def evaluate_all(self):
        mia = self.membership_inference_attack()
        gap = self.overfitting_gap()
        return {"membership_inference": mia, "overfitting_gap": gap,
                "privacy_risk_score": mia["membership_inference_resistance"]}


# ============================================================
# 4. Robustness Evaluation Module
# ============================================================

class RobustnessEvaluator:
    """Evaluate model robustness against adversarial perturbations and distribution shift."""

    def __init__(self, model, X_test, y_test):
        self.model = model
        self.X_test = np.array(X_test)
        self.y_test = np.array(y_test)

    def adversarial_perturbation(self, epsilons=None):
        """FGSM-like perturbation for tree models (gradient-free approximation)."""
        if epsilons is None:
            epsilons = [0.01, 0.05, 0.1, 0.2, 0.5]

        base_acc = accuracy_score(self.y_test, self.model.predict(self.X_test))
        results = {"base_accuracy": float(base_acc), "perturbation_results": []}

        for eps in epsilons:
            # Random directional perturbation
            perturbation = np.random.uniform(-eps, eps, self.X_test.shape)
            X_perturbed = self.X_test + perturbation
            perturbed_acc = accuracy_score(self.y_test, self.model.predict(X_perturbed))
            results["perturbation_results"].append({
                "epsilon": eps,
                "perturbed_accuracy": float(perturbed_acc),
                "accuracy_drop": float(base_acc - perturbed_acc),
            })
        return results

    def distribution_shift(self, shift_magnitudes=None):
        """Evaluate robustness under covariate shift."""
        if shift_magnitudes is None:
            shift_magnitudes = [0.0, 0.5, 1.0, 2.0, 3.0]

        base_acc = accuracy_score(self.y_test, self.model.predict(self.X_test))
        results = {"base_accuracy": float(base_acc), "shift_results": []}

        for mag in shift_magnitudes:
            X_shifted = self.X_test + mag * np.ones_like(self.X_test)
            shifted_acc = accuracy_score(self.y_test, self.model.predict(X_shifted))
            results["shift_results"].append({
                "shift_magnitude": mag,
                "shifted_accuracy": float(shifted_acc),
                "accuracy_drop": float(base_acc - shifted_acc),
            })
        return results

    def robustness_score(self):
        adv = self.adversarial_perturbation()
        dist = self.distribution_shift()
        # Average retention at moderate perturbation
        adv_retention = 1.0 - np.mean([r["accuracy_drop"] for r in adv["perturbation_results"][:3]])
        dist_retention = 1.0 - np.mean([r["accuracy_drop"] for r in dist["shift_results"][:3]])
        score = 0.5 * np.clip(adv_retention, 0, 1) + 0.5 * np.clip(dist_retention, 0, 1)
        return {"adversarial": adv, "distribution_shift": dist,
                "integrated_robustness_score": float(np.clip(score, 0, 1))}


# ============================================================
# 5. Environmental Impact Module
# ============================================================

class EnvironmentalEvaluator:
    """Estimate computational CO2 emissions."""

    def __init__(self):
        self.measurements = []

    def estimate_emissions(self, training_time_seconds, n_params, hardware="cpu",
                           power_draw_watts=None, carbon_intensity_gco2_kwh=400):
        """Estimate CO2 emissions from training."""
        if power_draw_watts is None:
            power_draw_watts = {"cpu": 65, "gpu_t4": 70, "gpu_v100": 300, "gpu_a100": 400}.get(hardware, 65)

        energy_kwh = (power_draw_watts * training_time_seconds) / (3600 * 1000)
        co2_grams = energy_kwh * carbon_intensity_gco2_kwh
        # Efficiency: params per gram CO2
        efficiency = n_params / max(co2_grams, 1e-6)

        result = {
            "training_time_seconds": training_time_seconds,
            "hardware": hardware,
            "power_draw_watts": power_draw_watts,
            "energy_kwh": float(energy_kwh),
            "co2_grams": float(co2_grams),
            "co2_kg": float(co2_grams / 1000),
            "carbon_intensity_gco2_kwh": carbon_intensity_gco2_kwh,
            "params_per_gco2": float(efficiency),
        }
        self.measurements.append(result)
        return result

    def environmental_score(self):
        """Score based on CO2 efficiency (normalized, higher is better)."""
        if not self.measurements:
            return {"environmental_score": 1.0, "total_co2_grams": 0.0}
        total_co2 = sum(m["co2_grams"] for m in self.measurements)
        # Benchmark: 100g CO2 → score 0.5; 0g → 1.0; 1000g → ~0
        score = np.exp(-total_co2 / 200)
        return {"total_co2_grams": float(total_co2),
                "total_co2_kg": float(total_co2 / 1000),
                "environmental_score": float(np.clip(score, 0, 1)),
                "measurements": self.measurements}


# ============================================================
# 6. Integrated Ethics Dashboard
# ============================================================

class EthicsDashboard:
    """Aggregate all ethics dimensions into a unified dashboard."""

    def __init__(self):
        self.scores = {}

    def add_dimension(self, name, score, details):
        self.scores[name] = {"score": float(score), "details": details}

    def overall_ethics_score(self, weights=None):
        if weights is None:
            weights = {
                "fairness": 0.25,
                "explainability": 0.20,
                "privacy": 0.20,
                "robustness": 0.20,
                "environment": 0.15,
            }
        total = 0.0
        for dim, w in weights.items():
            if dim in self.scores:
                total += w * self.scores[dim]["score"]
        return float(np.clip(total, 0, 1))

    def generate_report(self):
        report = {
            "timestamp": datetime.now(JST).isoformat(),
            "dimensions": self.scores,
            "overall_ethics_score": self.overall_ethics_score(),
            "grade": self._grade(self.overall_ethics_score()),
        }
        return report

    @staticmethod
    def _grade(score):
        if score >= 0.9: return "A"
        if score >= 0.8: return "B"
        if score >= 0.7: return "C"
        if score >= 0.6: return "D"
        return "F"

    def plot_radar(self, save_path="figures/ethics_radar.png"):
        dims = list(self.scores.keys())
        vals = [self.scores[d]["score"] for d in dims]
        N = len(dims)
        angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
        vals_plot = vals + vals[:1]
        angles += angles[:1]

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
        ax.fill(angles, vals_plot, alpha=0.25, color="#2196F3")
        ax.plot(angles, vals_plot, "o-", linewidth=2, color="#1565C0")
        ax.set_xticks(angles[:-1])
        labels = [d.replace("_", "\n").title() for d in dims]
        ax.set_xticklabels(labels, fontsize=12, fontweight="bold")
        ax.set_ylim(0, 1)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=9)
        ax.set_title(f"AI Ethics Evaluation Radar\nOverall Score: {self.overall_ethics_score():.3f}",
                      fontsize=14, fontweight="bold", pad=20)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    def plot_bar(self, save_path="figures/ethics_bar.png"):
        dims = list(self.scores.keys())
        vals = [self.scores[d]["score"] for d in dims]
        colors = ["#4CAF50" if v >= 0.7 else ("#FF9800" if v >= 0.5 else "#F44336") for v in vals]

        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(dims, vals, color=colors, edgecolor="white", height=0.6)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}", va="center", fontsize=12, fontweight="bold")
        ax.set_xlim(0, 1.15)
        ax.set_xlabel("Score (0-1)", fontsize=12)
        ax.set_title("AI Ethics Dimension Scores", fontsize=14, fontweight="bold")
        ax.axvline(x=0.7, color="gray", linestyle="--", alpha=0.5, label="Threshold (0.7)")
        ax.legend()
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path


# ============================================================
# Main Pipeline Execution
# ============================================================

def run_pipeline():
    log_entries = []
    def log(phase, event, **kwargs):
        entry = {"timestamp": datetime.now(JST).isoformat(), "phase": phase,
                 "event_type": event, "actor": "co-scientist",
                 "skill_or_tool": "ai-ethics-evaluation", **kwargs}
        log_entries.append(entry)

    log("pipeline", "run_started")
    print("=" * 60)
    print("AI Ethics Quantitative Evaluation Framework")
    print("=" * 60)

    # --- Data Generation ---
    print("\n[1/6] Generating synthetic medical dataset...")
    df = generate_medical_dataset(n=3000)
    feature_cols = [c for c in df.columns if c.startswith("feature_")]
    X = df[feature_cols].values
    y = df["diagnosis"].values
    sensitive = df["age_group"].values
    group_names = {0: "Young (<50)", 1: "Old (≥50)"}

    X_train, X_test, y_train, y_test, s_train, s_test = train_test_split(
        X, y, sensitive, test_size=0.3, random_state=42, stratify=y
    )
    log("pipeline", "data_generated", files_written=["data/medical_dataset.csv"])

    # --- Model Training ---
    print("[2/6] Training models...")
    t0 = time.time()
    model = GradientBoostingClassifier(n_estimators=150, max_depth=5, random_state=42)
    model.fit(X_train, y_train)
    train_time = time.time() - t0

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    print(f"  Model accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print(f"  Model AUC: {roc_auc_score(y_test, y_prob):.4f}")
    print(f"  Training time: {train_time:.2f}s")

    # --- 1. Fairness Evaluation ---
    print("\n[3/6] Evaluating Fairness...")
    fairness_eval = FairnessEvaluator(y_test, y_pred, y_prob, s_test, group_names)
    fairness_results = fairness_eval.evaluate_all()
    integrated_fairness = fairness_results["integrated"]["integrated_fairness_score"]
    print(f"  Statistical Parity Difference: {fairness_results['integrated']['spd']:.4f}")
    print(f"  Equalized Odds Difference: {fairness_results['integrated']['eod']:.4f}")
    print(f"  Calibration Gap: {fairness_results['integrated']['calibration_gap']:.4f}")
    print(f"  → Integrated Fairness Score: {integrated_fairness:.4f}")
    log("fairness", "evaluation_completed")

    # Fairness figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    # SP bar
    sp = fairness_results["statistical_parity"]["group_positive_rates"]
    axes[0].bar(sp.keys(), sp.values(), color=["#2196F3", "#FF9800"])
    axes[0].set_title("Statistical Parity\n(Positive Rate by Group)", fontweight="bold")
    axes[0].set_ylabel("Positive Rate")
    axes[0].set_ylim(0, 1)
    # EO
    eo = fairness_results["equalized_odds"]
    for lbl in ["label_0", "label_1"]:
        rates = eo[lbl]["group_rates"]
        x = list(rates.keys())
        axes[1].bar([f"{xi}\n(y={lbl[-1]})" for xi in x],
                    rates.values(), alpha=0.7, label=f"True Label={lbl[-1]}")
    axes[1].set_title("Equalized Odds\n(TPR/FPR by Group)", fontweight="bold")
    axes[1].legend()
    axes[1].set_ylim(0, 1)
    # Calibration
    cal = fairness_results["calibration"]["group_calibration_errors"]
    axes[2].bar(cal.keys(), cal.values(), color=["#4CAF50", "#F44336"])
    axes[2].set_title("Calibration Error\nby Group", fontweight="bold")
    axes[2].set_ylabel("Mean Calibration Error")
    plt.tight_layout()
    plt.savefig("figures/fairness_metrics.png", dpi=300, bbox_inches="tight")
    plt.close()

    # --- 2. Explainability Evaluation ---
    print("\n[4/6] Evaluating Explainability...")
    explain_eval = ExplainabilityEvaluator(model, X_test, feature_cols)
    explain_results = explain_eval.explainability_score()
    explain_score = explain_results["integrated_explainability_score"]
    print(f"  SHAP Consistency (rank corr): {explain_results['shap_consistency']['mean_rank_correlation']:.4f}")
    print(f"  Explanation Stability: {explain_results['explanation_stability']['mean_stability']:.4f}")
    print(f"  → Integrated Explainability Score: {explain_score:.4f}")
    log("explainability", "evaluation_completed")

    # SHAP importance plot
    import shap
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test[:200])
    if isinstance(shap_values, list):
        shap_values = shap_values[1]
    fig, ax = plt.subplots(figsize=(10, 8))
    mean_abs = np.abs(shap_values).mean(axis=0)
    sorted_idx = np.argsort(mean_abs)[-15:]
    ax.barh([feature_cols[i] for i in sorted_idx], mean_abs[sorted_idx], color="#7B1FA2")
    ax.set_xlabel("Mean |SHAP Value|", fontsize=12)
    ax.set_title("Top 15 Feature Importances (SHAP)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/shap_importance.png", dpi=300, bbox_inches="tight")
    plt.close()

    # --- 3. Privacy Evaluation ---
    print("\n[5/6] Evaluating Privacy Risk...")
    privacy_eval = PrivacyEvaluator(model, X_train, y_train, X_test, y_test)
    privacy_results = privacy_eval.evaluate_all()
    privacy_score = privacy_results["privacy_risk_score"]
    print(f"  MIA Attack Accuracy: {privacy_results['membership_inference']['attack_accuracy']:.4f}")
    print(f"  MIA Attack AUC: {privacy_results['membership_inference']['attack_auc']:.4f}")
    print(f"  Overfitting Gap: {privacy_results['overfitting_gap']['gap']:.4f}")
    print(f"  → Privacy Resistance Score: {privacy_score:.4f}")
    print(f"  Risk Level: {privacy_results['membership_inference']['risk_level']}")
    log("privacy", "evaluation_completed")

    # Privacy figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # MIA confidence distributions
    train_conf = np.max(model.predict_proba(X_train[:500]), axis=1)
    test_conf = np.max(model.predict_proba(X_test[:500]), axis=1)
    axes[0].hist(train_conf, bins=30, alpha=0.6, label="Train (Members)", color="#F44336")
    axes[0].hist(test_conf, bins=30, alpha=0.6, label="Test (Non-Members)", color="#2196F3")
    axes[0].set_xlabel("Prediction Confidence", fontsize=11)
    axes[0].set_ylabel("Count", fontsize=11)
    axes[0].set_title("Membership Inference\nConfidence Distribution", fontweight="bold")
    axes[0].legend()
    # Overfitting gap
    gap = privacy_results["overfitting_gap"]
    axes[1].bar(["Train Accuracy", "Test Accuracy"], [gap["train_accuracy"], gap["test_accuracy"]],
                color=["#FF9800", "#4CAF50"])
    axes[1].set_ylim(0, 1)
    axes[1].set_title(f"Overfitting Gap: {gap['gap']:.4f}", fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/privacy_metrics.png", dpi=300, bbox_inches="tight")
    plt.close()

    # --- 4. Robustness Evaluation ---
    print("\n[6/6] Evaluating Robustness...")
    robust_eval = RobustnessEvaluator(model, X_test, y_test)
    robust_results = robust_eval.robustness_score()
    robust_score = robust_results["integrated_robustness_score"]
    print(f"  Adversarial perturbation results:")
    for r in robust_results["adversarial"]["perturbation_results"]:
        print(f"    ε={r['epsilon']:.2f}: accuracy={r['perturbed_accuracy']:.4f} (drop={r['accuracy_drop']:.4f})")
    print(f"  → Integrated Robustness Score: {robust_score:.4f}")
    log("robustness", "evaluation_completed")

    # Robustness figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    adv = robust_results["adversarial"]
    eps_vals = [r["epsilon"] for r in adv["perturbation_results"]]
    acc_vals = [r["perturbed_accuracy"] for r in adv["perturbation_results"]]
    axes[0].plot(eps_vals, acc_vals, "o-", linewidth=2, color="#D32F2F", markersize=8)
    axes[0].axhline(y=adv["base_accuracy"], linestyle="--", color="gray", alpha=0.7, label="Baseline")
    axes[0].set_xlabel("Perturbation Magnitude (ε)", fontsize=11)
    axes[0].set_ylabel("Accuracy", fontsize=11)
    axes[0].set_title("Adversarial Robustness", fontweight="bold")
    axes[0].legend()

    dist = robust_results["distribution_shift"]
    shift_vals = [r["shift_magnitude"] for r in dist["shift_results"]]
    shift_acc = [r["shifted_accuracy"] for r in dist["shift_results"]]
    axes[1].plot(shift_vals, shift_acc, "s-", linewidth=2, color="#1565C0", markersize=8)
    axes[1].set_xlabel("Distribution Shift Magnitude", fontsize=11)
    axes[1].set_ylabel("Accuracy", fontsize=11)
    axes[1].set_title("Distribution Shift Robustness", fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/robustness_metrics.png", dpi=300, bbox_inches="tight")
    plt.close()

    # --- 5. Environmental Impact ---
    env_eval = EnvironmentalEvaluator()
    env_result = env_eval.estimate_emissions(
        training_time_seconds=train_time,
        n_params=model.n_estimators * (2 ** model.max_depth),
        hardware="cpu",
    )
    env_scores = env_eval.environmental_score()
    env_score = env_scores["environmental_score"]
    print(f"\n  Environmental Impact:")
    print(f"    Energy: {env_result['energy_kwh']:.6f} kWh")
    print(f"    CO2: {env_result['co2_grams']:.4f} g")
    print(f"    → Environmental Score: {env_score:.4f}")
    log("environment", "evaluation_completed")

    # --- Integrated Dashboard ---
    dashboard = EthicsDashboard()
    dashboard.add_dimension("fairness", integrated_fairness, fairness_results)
    dashboard.add_dimension("explainability", explain_score, explain_results)
    dashboard.add_dimension("privacy", privacy_score, privacy_results)
    dashboard.add_dimension("robustness", robust_score, robust_results)
    dashboard.add_dimension("environment", env_score, env_scores)

    final_report = dashboard.generate_report()
    dashboard.plot_radar()
    dashboard.plot_bar()

    # Save results
    with open("results/ethics_evaluation_results.json", "w") as f:
        json.dump(final_report, f, indent=2, default=str)

    # Summary table
    summary_df = pd.DataFrame([
        {"Dimension": dim.title(), "Score": info["score"],
         "Grade": EthicsDashboard._grade(info["score"])}
        for dim, info in final_report["dimensions"].items()
    ])
    summary_df.loc[len(summary_df)] = {"Dimension": "OVERALL", "Score": final_report["overall_ethics_score"],
                                        "Grade": final_report["grade"]}
    summary_df.to_csv("results/ethics_summary.csv", index=False)

    print("\n" + "=" * 60)
    print("INTEGRATED ETHICS EVALUATION RESULTS")
    print("=" * 60)
    print(summary_df.to_string(index=False))
    print(f"\nOverall Ethics Score: {final_report['overall_ethics_score']:.4f} (Grade: {final_report['grade']})")

    # Heatmap figure
    fig, ax = plt.subplots(figsize=(8, 2))
    score_data = [[info["score"] for info in final_report["dimensions"].values()]]
    sns.heatmap(score_data, annot=True, fmt=".3f", cmap="RdYlGn",
                xticklabels=[d.title() for d in final_report["dimensions"].keys()],
                yticklabels=["Score"], vmin=0, vmax=1, ax=ax,
                annot_kws={"fontsize": 14, "fontweight": "bold"})
    ax.set_title("Ethics Evaluation Heatmap", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/ethics_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()

    log("pipeline", "run_completed", status="ok",
        files_written=["results/ethics_evaluation_results.json", "results/ethics_summary.csv",
                        "figures/ethics_radar.png", "figures/ethics_bar.png",
                        "figures/fairness_metrics.png", "figures/shap_importance.png",
                        "figures/privacy_metrics.png", "figures/robustness_metrics.png",
                        "figures/ethics_heatmap.png"])

    # Write process log
    with open("logs/process-log.jsonl", "w") as f:
        for entry in log_entries:
            f.write(json.dumps(entry, default=str) + "\n")

    return final_report


if __name__ == "__main__":
    run_pipeline()
