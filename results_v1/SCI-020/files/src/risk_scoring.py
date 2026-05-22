"""
Risk Scoring and Alert Threshold Optimization Module.
Weighted composite scoring with Bayesian updates and ROC-based threshold optimization.
"""

import numpy as np
import pandas as pd
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")

try:
    from sklearn.metrics import roc_curve, auc
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


# Default component weights
DEFAULT_WEIGHTS = {
    "genomic": 0.30,
    "epidemiology": 0.35,
    "nlp": 0.20,
    "rt": 0.15,
}

RISK_LEVELS = [
    (90, "CRITICAL"),
    (75, "HIGH"),
    (50, "MEDIUM"),
    (25, "LOW"),
    (0,  "MINIMAL"),
]


def classify_risk_level(score: float) -> str:
    """Map a [0-100] score to a risk level label."""
    for threshold, label in RISK_LEVELS:
        if score >= threshold:
            return label
    return "MINIMAL"


class RiskScoreCalculator:
    """Computes composite risk scores from multiple signal sources."""

    def __init__(self, weights: Optional[Dict[str, float]] = None):
        self.weights = weights or DEFAULT_WEIGHTS
        # Normalize weights
        total = sum(self.weights.values())
        self.weights = {k: v / total for k, v in self.weights.items()}

    def _normalize_rt_score(self, latest_rt: float) -> float:
        """Map Rt to [0-100] risk score."""
        if latest_rt <= 0.5:
            return 0.0
        elif latest_rt >= 3.0:
            return 100.0
        else:
            return round((latest_rt - 0.5) / 2.5 * 100, 2)

    def compute_composite_risk(self,
                                genomic_score: float,
                                epi_score: float,
                                nlp_score: float,
                                rt_score: float) -> Dict:
        """
        Weighted composite risk score [0-100].
        Each input should already be in [0-100] range.
        """
        components = {
            "genomic": np.clip(genomic_score, 0, 100),
            "epidemiology": np.clip(epi_score, 0, 100),
            "nlp": np.clip(nlp_score, 0, 100),
            "rt": np.clip(rt_score, 0, 100),
        }
        composite = sum(self.weights[k] * components[k] for k in components)
        composite = round(float(composite), 2)

        return {
            "composite_score": composite,
            "risk_level": classify_risk_level(composite),
            "component_scores": {k: round(float(v), 2) for k, v in components.items()},
            "component_weights": {k: round(v, 3) for k, v in self.weights.items()},
            "timestamp": datetime.now().isoformat(),
        }

    def apply_bayesian_update(self, prior_score: float,
                               new_evidence_score: float,
                               evidence_weight: float = 0.3) -> float:
        """
        Bayesian-like update: posterior = (1-w)*prior + w*evidence.
        """
        posterior = (1 - evidence_weight) * prior_score + evidence_weight * new_evidence_score
        return round(float(np.clip(posterior, 0, 100)), 2)

    def compute_uncertainty(self, component_scores: Dict[str, float],
                             n_samples: int = 1000, seed: int = 42) -> Dict:
        """Monte Carlo uncertainty quantification for composite score."""
        rng = np.random.default_rng(seed)
        samples = []
        for _ in range(n_samples):
            # Add noise to each component
            noisy = {
                k: np.clip(v + rng.normal(0, 5), 0, 100)
                for k, v in component_scores.items()
            }
            composite = sum(self.weights.get(k, 0.25) * noisy[k] for k in noisy)
            samples.append(float(composite))

        samples = np.array(samples)
        return {
            "mean": round(float(samples.mean()), 2),
            "std": round(float(samples.std()), 2),
            "ci_lower_95": round(float(np.percentile(samples, 2.5)), 2),
            "ci_upper_95": round(float(np.percentile(samples, 97.5)), 2),
            "ci_width": round(float(np.percentile(samples, 97.5) - np.percentile(samples, 2.5)), 2),
        }


class AlertThresholdOptimizer:
    """
    Optimizes alert thresholds using ROC analysis.
    Simulates historical outbreak data to find optimal thresholds.
    """

    def generate_historical_data(self, n_events: int = 500,
                                  seed: int = 42) -> pd.DataFrame:
        """
        Generate synthetic historical risk scores with ground truth labels.
        Labels: 1 = true outbreak, 0 = no outbreak.
        """
        rng = np.random.default_rng(seed)
        # Positive cases (true outbreaks): higher risk scores
        n_pos = n_events // 3
        n_neg = n_events - n_pos

        pos_scores = np.clip(rng.normal(68, 15, n_pos), 0, 100)
        neg_scores = np.clip(rng.normal(32, 18, n_neg), 0, 100)

        scores = np.concatenate([pos_scores, neg_scores])
        labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
        idx = rng.permutation(len(scores))

        return pd.DataFrame({
            "risk_score": np.round(scores[idx], 2),
            "outbreak_label": labels[idx].astype(int),
        })

    def optimize_thresholds(self, historical_df: pd.DataFrame,
                             target_sensitivity: float = 0.95) -> Dict:
        """
        ROC-based threshold optimization.
        Finds threshold achieving target_sensitivity with maximum specificity.
        """
        scores = historical_df["risk_score"].values
        labels = historical_df["outbreak_label"].values

        if HAS_SKLEARN:
            fpr, tpr, thresholds = roc_curve(labels, scores)
            roc_auc = auc(fpr, tpr)

            # Find threshold closest to target sensitivity
            target_idx = np.argmin(np.abs(tpr - target_sensitivity))
            optimal_threshold = float(thresholds[target_idx])
            achieved_sensitivity = float(tpr[target_idx])
            achieved_specificity = float(1 - fpr[target_idx])
        else:
            # Manual ROC
            thresholds = np.linspace(0, 100, 100)
            tpr_list, fpr_list = [], []
            for t in thresholds:
                pred = (scores >= t).astype(int)
                tp = np.sum((pred == 1) & (labels == 1))
                fn = np.sum((pred == 0) & (labels == 1))
                fp = np.sum((pred == 1) & (labels == 0))
                tn = np.sum((pred == 0) & (labels == 0))
                tpr_list.append(tp / max(tp + fn, 1))
                fpr_list.append(fp / max(fp + tn, 1))
            tpr = np.array(tpr_list)
            fpr_arr = np.array(fpr_list)
            roc_auc = float(np.trapz(tpr[::-1], fpr_arr[::-1]))
            target_idx = np.argmin(np.abs(tpr - target_sensitivity))
            optimal_threshold = float(thresholds[target_idx])
            achieved_sensitivity = float(tpr[target_idx])
            achieved_specificity = float(1 - fpr_arr[target_idx])

        # Youden's J statistic optimal point
        if HAS_SKLEARN:
            j_scores = tpr - fpr
            j_optimal_idx = np.argmax(j_scores)
            youdens_threshold = float(thresholds[j_optimal_idx])
        else:
            youdens_threshold = optimal_threshold

        return {
            "roc_auc": round(float(roc_auc), 4),
            "optimal_threshold_95sens": round(optimal_threshold, 2),
            "achieved_sensitivity": round(achieved_sensitivity, 4),
            "achieved_specificity": round(achieved_specificity, 4),
            "youdens_j_threshold": round(youdens_threshold, 2),
            "target_sensitivity": target_sensitivity,
            "n_historical_events": len(historical_df),
        }


class AlertGenerator:
    """Generates structured alerts from risk scores."""

    def generate_alert(self, risk_result: Dict,
                        genomic_result: Dict,
                        epi_result: Dict,
                        rt_result: Dict,
                        nlp_result: Dict) -> Dict:
        """Generate a structured outbreak alert."""
        score = risk_result["composite_score"]
        level = risk_result["risk_level"]
        ts = datetime.now().isoformat()

        # Build evidence summary
        evidence = []
        if rt_result.get("epidemic_growing"):
            evidence.append(f"Rt > 1 in {rt_result.get('n_countries_growing', 0)} countries "
                            f"(global mean Rt={rt_result.get('global_rt_mean', 'N/A')})")
        if genomic_result.get("n_novel_emerging", 0) > 0:
            evidence.append(f"{genomic_result['n_novel_emerging']} emerging lineage(s) detected")
        if epi_result.get("n_anomalies_detected", 0) > 0:
            evidence.append(f"{epi_result['n_anomalies_detected']} epidemiological anomalies")
        if nlp_result.get("n_critical", 0) > 0:
            evidence.append(f"{nlp_result['n_critical']} CRITICAL-severity alerts in surveillance feeds")

        return {
            "alert_id": f"ALERT-{ts[:10].replace('-', '')}-{int(score):03d}",
            "timestamp": ts,
            "risk_score": score,
            "risk_level": level,
            "action_required": level in ["HIGH", "CRITICAL"],
            "evidence_summary": evidence,
            "recommended_actions": _get_recommended_actions(level),
        }


def _get_recommended_actions(risk_level: str) -> List[str]:
    """Return recommended actions based on risk level."""
    base = ["Continue routine surveillance", "Monitor situation closely"]
    if risk_level == "MINIMAL":
        return base
    elif risk_level == "LOW":
        return base + ["Increase sequencing frequency", "Alert regional WHO offices"]
    elif risk_level == "MEDIUM":
        return base + ["Activate enhanced surveillance", "Notify national health authorities",
                       "Prepare response protocols"]
    elif risk_level == "HIGH":
        return ["Activate emergency operations center", "Issue public health advisory",
                "Mobilize response teams", "Increase testing capacity",
                "Consider travel advisories"]
    else:  # CRITICAL
        return ["Declare public health emergency", "Activate PHEIC protocols",
                "Immediate WHO notification", "Deploy rapid response teams",
                "Implement containment measures", "Emergency resource mobilization"]


def run_risk_scoring(genomic_result: Dict, epi_result: Dict,
                     rt_result: Dict, nlp_result: Dict,
                     config: Optional[Dict] = None) -> Dict:
    """Run the full risk scoring pipeline."""
    config = config or {}
    weights = config.get("weights", DEFAULT_WEIGHTS)

    calculator = RiskScoreCalculator(weights=weights)
    optimizer = AlertThresholdOptimizer()
    alert_gen = AlertGenerator()

    # Extract component scores
    genomic_score = min(100, genomic_result.get("n_novel_emerging", 0) * 20
                        + genomic_result.get("diversity_metrics", {}).get("shannon_diversity", 0) * 10)

    epi_score = epi_result.get("integrated_signals", {}).get("composite_epi_score", 30.0)

    nlp_score = nlp_result.get("signal_score", 40.0)

    # Rt score from global mean
    global_rt = rt_result.get("global_rt_mean", 1.0)
    rt_score = calculator._normalize_rt_score(global_rt)

    # Composite risk
    risk_result = calculator.compute_composite_risk(
        genomic_score, epi_score, nlp_score, rt_score
    )

    # Uncertainty quantification
    uncertainty = calculator.compute_uncertainty(risk_result["component_scores"])
    risk_result["uncertainty"] = uncertainty

    # Threshold optimization
    historical_df = optimizer.generate_historical_data()
    threshold_result = optimizer.optimize_thresholds(historical_df)

    # Generate alert
    alert = alert_gen.generate_alert(risk_result, genomic_result, epi_result,
                                      rt_result, nlp_result)

    return {
        "risk_result": risk_result,
        "threshold_optimization": threshold_result,
        "alert": alert,
        "component_scores": {
            "genomic_score": round(float(genomic_score), 2),
            "epi_score": round(float(epi_score), 2),
            "nlp_score": round(float(nlp_score), 2),
            "rt_score": round(float(rt_score), 2),
        },
    }


if __name__ == "__main__":
    # Demo run
    dummy_genomic = {"n_novel_emerging": 2, "diversity_metrics": {"shannon_diversity": 1.8}}
    dummy_epi = {"integrated_signals": {"composite_epi_score": 45.0}, "n_anomalies_detected": 3}
    dummy_rt = {"global_rt_mean": 1.35, "n_countries_growing": 3, "epidemic_growing": True}
    dummy_nlp = {"signal_score": 55.0, "n_critical": 1, "n_high": 2}
    results = run_risk_scoring(dummy_genomic, dummy_epi, dummy_rt, dummy_nlp)
    print(f"Composite risk: {results['risk_result']['composite_score']}")
    print(f"Risk level: {results['risk_result']['risk_level']}")
    print(f"ROC AUC: {results['threshold_optimization']['roc_auc']}")
    print(f"Alert: {results['alert']['alert_id']}")
