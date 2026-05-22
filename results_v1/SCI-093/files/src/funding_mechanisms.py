"""
Module 3: Funding Allocation Mechanism Simulation
Simulates peer review, lottery, and automated allocation mechanisms.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import json
import os

np.random.seed(42)


class FundingMechanism:
    def __init__(self, name):
        self.name = name

    def allocate(self, applicants_df, budget, grant_size):
        raise NotImplementedError


class PeerReviewMechanism(FundingMechanism):
    def __init__(self, n_reviewers=3, bias_strength=0.1):
        super().__init__("Peer Review")
        self.n_reviewers = n_reviewers
        self.bias_strength = bias_strength

    def allocate(self, applicants_df, budget, grant_size):
        n_grants = int(budget / grant_size)
        scores = []
        for _, app in applicants_df.iterrows():
            true_quality = app.get("score_hybrid", app.get("h_index", 5) / 50)
            reviewer_scores = []
            for _ in range(self.n_reviewers):
                bias = 0
                if app.get("career_stage") == "senior":
                    bias += self.bias_strength * 0.5
                if app.get("h_index", 0) > 20:
                    bias += self.bias_strength * 0.3
                noise = np.random.normal(0, 0.15)
                reviewer_scores.append(np.clip(true_quality + bias + noise, 0, 1))
            scores.append(np.mean(reviewer_scores))
        applicants_df = applicants_df.copy()
        applicants_df["review_score"] = scores
        applicants_df = applicants_df.sort_values("review_score", ascending=False)
        funded = applicants_df.head(n_grants).index.tolist()
        return funded, applicants_df


class LotteryMechanism(FundingMechanism):
    def __init__(self, threshold_percentile=50):
        super().__init__("Lottery")
        self.threshold_percentile = threshold_percentile

    def allocate(self, applicants_df, budget, grant_size):
        n_grants = int(budget / grant_size)
        quality = applicants_df.get("score_hybrid",
                                    applicants_df.get("h_index", pd.Series(np.zeros(len(applicants_df)))) / 50)
        threshold = np.percentile(quality, self.threshold_percentile)
        qualified = applicants_df[quality >= threshold]
        if len(qualified) <= n_grants:
            funded = qualified.index.tolist()
        else:
            funded = list(np.random.choice(qualified.index, size=n_grants, replace=False))
        return funded, applicants_df


class AutomatedMechanism(FundingMechanism):
    def __init__(self, diversity_weight=0.3):
        super().__init__("Automated")
        self.diversity_weight = diversity_weight

    def allocate(self, applicants_df, budget, grant_size):
        n_grants = int(budget / grant_size)
        applicants_df = applicants_df.copy()
        quality = applicants_df.get("score_hybrid",
                                    applicants_df.get("h_index", pd.Series(np.zeros(len(applicants_df)))) / 50)
        diversity_bonus = np.zeros(len(applicants_df))
        for i, (_, row) in enumerate(applicants_df.iterrows()):
            bonus = 0
            if row.get("career_stage") == "early":
                bonus += 0.1
            if row.get("gender") == "F":
                bonus += 0.05
            if row.get("region") in ["LatinAmerica", "Africa"]:
                bonus += 0.1
            diversity_bonus[i] = bonus
        composite = (1 - self.diversity_weight) * quality.values + self.diversity_weight * diversity_bonus
        applicants_df["auto_score"] = composite
        applicants_df = applicants_df.sort_values("auto_score", ascending=False)
        funded = applicants_df.head(n_grants).index.tolist()
        return funded, applicants_df


def evaluate_allocation(applicants_df, funded_ids, label=""):
    """Evaluate funding allocation outcomes."""
    funded = applicants_df.loc[funded_ids]
    total = applicants_df
    results = {"mechanism": label, "n_funded": len(funded)}

    if "score_hybrid" in total.columns:
        results["avg_quality_funded"] = round(funded["score_hybrid"].mean(), 4)
        results["avg_quality_all"] = round(total["score_hybrid"].mean(), 4)
        results["quality_ratio"] = round(
            funded["score_hybrid"].mean() / max(0.001, total["score_hybrid"].mean()), 3)

    for g in ["M", "F", "Other"]:
        results[f"pct_funded_{g}"] = round(
            (funded["gender"] == g).mean() * 100, 1) if "gender" in funded.columns else 0
        results[f"pct_total_{g}"] = round(
            (total["gender"] == g).mean() * 100, 1) if "gender" in total.columns else 0

    for stage in ["early", "mid", "senior"]:
        results[f"pct_funded_{stage}"] = round(
            (funded["career_stage"] == stage).mean() * 100, 1) if "career_stage" in funded.columns else 0

    if "region" in funded.columns:
        region_counts = funded["region"].value_counts(normalize=True)
        entropy = -np.sum(region_counts * np.log2(region_counts + 1e-10))
        results["region_entropy_funded"] = round(entropy, 3)
        region_counts_all = total["region"].value_counts(normalize=True)
        entropy_all = -np.sum(region_counts_all * np.log2(region_counts_all + 1e-10))
        results["region_entropy_all"] = round(entropy_all, 3)

    if "field" in funded.columns:
        field_counts = funded["field"].value_counts(normalize=True)
        field_entropy = -np.sum(field_counts * np.log2(field_counts + 1e-10))
        results["field_entropy_funded"] = round(field_entropy, 3)

    return results


def plot_mechanism_comparison(all_results, output_dir="figures"):
    """Plot comparison of allocation mechanisms."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    df = pd.DataFrame(all_results)

    ax = axes[0, 0]
    mechanisms = df["mechanism"]
    quality = df["avg_quality_funded"]
    bars = ax.bar(mechanisms, quality, color=["#2196F3", "#4CAF50", "#FF9800"], width=0.6)
    ax.axhline(y=df["avg_quality_all"].iloc[0], color="red", linestyle="--",
               label="Population Mean", alpha=0.7)
    ax.set_title("Average Quality of Funded Researchers", fontsize=13)
    ax.set_ylabel("Hybrid Score", fontsize=11)
    ax.legend()
    for bar, val in zip(bars, quality):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f"{val:.3f}", ha="center", fontsize=10)

    ax = axes[0, 1]
    x = np.arange(len(mechanisms))
    w = 0.25
    ax.bar(x - w, df["pct_funded_M"], w, label="Male", color="#2196F3", alpha=0.8)
    ax.bar(x, df["pct_funded_F"], w, label="Female", color="#E91E63", alpha=0.8)
    ax.bar(x + w, df["pct_funded_Other"], w, label="Other", color="#9C27B0", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(mechanisms)
    ax.set_title("Gender Distribution of Funded Researchers (%)", fontsize=13)
    ax.set_ylabel("Percentage", fontsize=11)
    ax.legend()

    ax = axes[1, 0]
    x = np.arange(len(mechanisms))
    w = 0.25
    ax.bar(x - w, df["pct_funded_early"], w, label="Early", color="#4CAF50", alpha=0.8)
    ax.bar(x, df["pct_funded_mid"], w, label="Mid", color="#FF9800", alpha=0.8)
    ax.bar(x + w, df["pct_funded_senior"], w, label="Senior", color="#F44336", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(mechanisms)
    ax.set_title("Career Stage Distribution of Funded (%)", fontsize=13)
    ax.set_ylabel("Percentage", fontsize=11)
    ax.legend()

    ax = axes[1, 1]
    entropy_funded = df["region_entropy_funded"]
    entropy_all = df["region_entropy_all"]
    x = np.arange(len(mechanisms))
    ax.bar(x - 0.15, entropy_funded, 0.3, label="Funded", color="#2196F3")
    ax.bar(x + 0.15, entropy_all, 0.3, label="All Applicants", color="#BDBDBD")
    ax.set_xticks(x)
    ax.set_xticklabels(mechanisms)
    ax.set_title("Region Diversity (Shannon Entropy)", fontsize=13)
    ax.set_ylabel("Entropy (bits)", fontsize=11)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"{output_dir}/fig6_mechanism_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_funding_simulation(applicants_df, output_dir="figures", results_dir="results"):
    """Run funding mechanism simulation."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    budget = 1_000_000
    grant_size = 50_000

    mechanisms = [
        PeerReviewMechanism(n_reviewers=3, bias_strength=0.15),
        LotteryMechanism(threshold_percentile=40),
        AutomatedMechanism(diversity_weight=0.3),
    ]

    all_results = []
    for mech in mechanisms:
        funded, scored_df = mech.allocate(applicants_df, budget, grant_size)
        result = evaluate_allocation(applicants_df, funded, mech.name)
        all_results.append(result)
        print(f"  {mech.name}: {result['n_funded']} funded, "
              f"quality={result.get('avg_quality_funded', 'N/A'):.3f}, "
              f"region_entropy={result.get('region_entropy_funded', 'N/A'):.3f}")

    plot_mechanism_comparison(all_results, output_dir)

    with open(f"{results_dir}/funding_mechanisms.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print("Funding mechanism simulation complete.")
    return all_results
