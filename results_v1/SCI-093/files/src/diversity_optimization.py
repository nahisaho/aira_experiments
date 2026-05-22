"""
Module 5: Diversity-Constrained Optimization
Multi-objective optimization of funding allocation under diversity constraints.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import linprog, minimize
import json
import os

np.random.seed(42)


def build_optimization_problem(df, budget=1_000_000, grant_size=50_000):
    """
    Formulate the funding allocation as a multi-objective optimization:
    - Maximize: total research quality
    - Subject to: diversity constraints (gender, region, career stage, field)
    """
    n = len(df)
    n_grants = int(budget / grant_size)

    quality = df["score_hybrid"].values if "score_hybrid" in df.columns else np.random.random(n)

    constraints_info = {
        "n_applicants": n,
        "n_grants": n_grants,
        "budget": budget,
        "grant_size": grant_size,
    }

    return quality, n_grants, constraints_info


def greedy_diverse_allocation(df, quality, n_grants,
                              gender_targets=None, region_targets=None,
                              early_career_min=0.2):
    """
    Greedy algorithm with diversity constraints.
    Iteratively select highest quality researcher that doesn't violate constraints.
    """
    if gender_targets is None:
        gender_targets = {"M": 0.50, "F": 0.40, "Other": 0.10}
    if region_targets is None:
        region_targets = {"Asia": 0.30, "Europe": 0.25, "NorthAmerica": 0.25,
                          "LatinAmerica": 0.12, "Africa": 0.08}

    selected = []
    remaining = list(range(len(df)))

    gender_counts = {g: 0 for g in gender_targets}
    region_counts = {r: 0 for r in region_targets}
    early_count = 0

    sorted_idx = sorted(remaining, key=lambda i: quality[i], reverse=True)

    for idx in sorted_idx:
        if len(selected) >= n_grants:
            break

        row = df.iloc[idx]
        g = row["gender"]
        r = row["region"]
        stage = row["career_stage"]

        # Check constraints (soft: allow exceeding by small margin)
        g_limit = int(n_grants * gender_targets.get(g, 0.5) * 1.3) + 1
        r_limit = int(n_grants * region_targets.get(r, 0.3) * 1.5) + 1

        if gender_counts.get(g, 0) < g_limit and region_counts.get(r, 0) < r_limit:
            selected.append(idx)
            gender_counts[g] = gender_counts.get(g, 0) + 1
            region_counts[r] = region_counts.get(r, 0) + 1
            if stage == "early":
                early_count += 1

    # If we haven't filled all slots, fill remaining without constraints
    if len(selected) < n_grants:
        for idx in sorted_idx:
            if idx not in selected and len(selected) < n_grants:
                selected.append(idx)

    return selected


def pareto_frontier_simulation(df, quality, n_grants, n_samples=500):
    """
    Monte Carlo simulation of Pareto frontier between efficiency and diversity.
    """
    pareto_points = []

    for _ in range(n_samples):
        # Random diversity weight
        div_weight = np.random.random()

        # Score = (1 - w) * quality + w * diversity_bonus
        diversity_bonus = np.zeros(len(df))
        for i, (_, row) in enumerate(df.iterrows()):
            bonus = 0
            if row.get("career_stage") == "early":
                bonus += 0.15
            if row.get("gender") == "F":
                bonus += 0.1
            elif row.get("gender") == "Other":
                bonus += 0.15
            if row.get("region") in ["LatinAmerica", "Africa"]:
                bonus += 0.15
            elif row.get("region") == "Asia":
                bonus += 0.05
            diversity_bonus[i] = bonus

        composite = (1 - div_weight) * quality + div_weight * diversity_bonus
        top_idx = np.argsort(composite)[-n_grants:]

        # Compute metrics for this allocation
        avg_quality = quality[top_idx].mean()

        funded = df.iloc[top_idx]
        # Diversity: Shannon entropy of region distribution
        region_counts = funded["region"].value_counts(normalize=True)
        region_entropy = -np.sum(region_counts * np.log2(region_counts + 1e-10))

        # Gender balance (distance from equal representation)
        gender_counts = funded["gender"].value_counts(normalize=True)
        gender_balance = 1 - abs(gender_counts.get("M", 0) - gender_counts.get("F", 0))

        # Early career proportion
        early_pct = (funded["career_stage"] == "early").mean()

        pareto_points.append({
            "diversity_weight": round(div_weight, 3),
            "avg_quality": round(avg_quality, 4),
            "region_entropy": round(region_entropy, 3),
            "gender_balance": round(gender_balance, 3),
            "early_career_pct": round(early_pct, 3),
            "composite_diversity": round(
                0.4 * region_entropy / 2.5 + 0.3 * gender_balance + 0.3 * early_pct, 3
            ),
        })

    return pd.DataFrame(pareto_points)


def plot_pareto_frontier(pareto_df, output_dir="figures"):
    """Plot Pareto frontier of efficiency vs diversity."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Quality vs Region Diversity
    ax = axes[0]
    scatter = ax.scatter(pareto_df["avg_quality"], pareto_df["region_entropy"],
                         c=pareto_df["diversity_weight"], cmap="viridis", alpha=0.5, s=15)
    ax.set_xlabel("Average Quality (Efficiency)", fontsize=12)
    ax.set_ylabel("Region Entropy (Diversity)", fontsize=12)
    ax.set_title("Efficiency vs Region Diversity", fontsize=13)
    plt.colorbar(scatter, ax=ax, label="Diversity Weight")

    # Quality vs Gender Balance
    ax = axes[1]
    scatter = ax.scatter(pareto_df["avg_quality"], pareto_df["gender_balance"],
                         c=pareto_df["diversity_weight"], cmap="viridis", alpha=0.5, s=15)
    ax.set_xlabel("Average Quality (Efficiency)", fontsize=12)
    ax.set_ylabel("Gender Balance", fontsize=12)
    ax.set_title("Efficiency vs Gender Balance", fontsize=13)
    plt.colorbar(scatter, ax=ax, label="Diversity Weight")

    # Quality vs Composite Diversity
    ax = axes[2]
    scatter = ax.scatter(pareto_df["avg_quality"], pareto_df["composite_diversity"],
                         c=pareto_df["diversity_weight"], cmap="viridis", alpha=0.5, s=15)
    ax.set_xlabel("Average Quality (Efficiency)", fontsize=12)
    ax.set_ylabel("Composite Diversity Score", fontsize=12)
    ax.set_title("Pareto Frontier: Efficiency vs Diversity", fontsize=13)
    plt.colorbar(scatter, ax=ax, label="Diversity Weight")

    # Find and mark Pareto-optimal points
    qualities = pareto_df["avg_quality"].values
    diversities = pareto_df["composite_diversity"].values
    pareto_mask = np.ones(len(qualities), dtype=bool)
    for i in range(len(qualities)):
        for j in range(len(qualities)):
            if i != j and qualities[j] >= qualities[i] and diversities[j] >= diversities[i]:
                if qualities[j] > qualities[i] or diversities[j] > diversities[i]:
                    pareto_mask[i] = False
                    break
    pareto_optimal = pareto_df[pareto_mask].sort_values("avg_quality")
    ax.plot(pareto_optimal["avg_quality"], pareto_optimal["composite_diversity"],
            "r-", linewidth=2, label="Pareto Frontier", alpha=0.7)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"{output_dir}/fig8_pareto_frontier.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    return pareto_mask


def run_diversity_optimization(df, output_dir="figures", results_dir="results"):
    """Run diversity-constrained optimization."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    quality, n_grants, info = build_optimization_problem(df)

    # Greedy diverse allocation
    selected = greedy_diverse_allocation(df, quality, n_grants)
    funded_df = df.iloc[selected]

    greedy_result = {
        "avg_quality": round(quality[selected].mean(), 4),
        "gender_dist": funded_df["gender"].value_counts().to_dict(),
        "region_dist": funded_df["region"].value_counts().to_dict(),
        "career_stage_dist": funded_df["career_stage"].value_counts().to_dict(),
        "field_dist": funded_df["field"].value_counts().to_dict(),
    }

    # Pareto frontier
    pareto_df = pareto_frontier_simulation(df, quality, n_grants, n_samples=500)
    pareto_mask = plot_pareto_frontier(pareto_df, output_dir)

    # Summary statistics
    opt_results = {
        "greedy_allocation": greedy_result,
        "pareto_analysis": {
            "n_samples": len(pareto_df),
            "n_pareto_optimal": int(pareto_mask.sum()),
            "quality_range": [round(pareto_df["avg_quality"].min(), 4),
                              round(pareto_df["avg_quality"].max(), 4)],
            "diversity_range": [round(pareto_df["composite_diversity"].min(), 3),
                                round(pareto_df["composite_diversity"].max(), 3)],
            "optimal_tradeoff_point": {
                "avg_quality": round(pareto_df.loc[pareto_mask, "avg_quality"].median(), 4),
                "composite_diversity": round(
                    pareto_df.loc[pareto_mask, "composite_diversity"].median(), 3),
            },
        },
        "constraints": info,
    }

    with open(f"{results_dir}/diversity_optimization.json", "w") as f:
        json.dump(opt_results, f, indent=2, default=str)

    pareto_df.to_csv(f"{results_dir}/pareto_points.csv", index=False)

    print("Diversity optimization complete.")
    print(f"  Greedy avg quality: {greedy_result['avg_quality']:.4f}")
    print(f"  Pareto optimal points: {int(pareto_mask.sum())}/{len(pareto_df)}")

    return opt_results


if __name__ == "__main__":
    from metrics import generate_researcher_profiles, compute_composite_scores
    df = generate_researcher_profiles(200)
    df = compute_composite_scores(df)
    run_diversity_optimization(df)
