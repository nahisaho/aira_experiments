"""
Module 6: KAKENHI (科研費) Funding Allocation Case Study
Simulates the Japanese research funding system (KAKENHI) with
realistic parameters for efficiency evaluation.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import json
import os

np.random.seed(42)


# KAKENHI grant categories (simplified)
KAKENHI_CATEGORIES = {
    "S": {"budget_per_grant": 50_000_000, "n_grants": 5, "success_rate": 0.05,
           "label": "基盤研究(S)"},
    "A": {"budget_per_grant": 20_000_000, "n_grants": 50, "success_rate": 0.15,
           "label": "基盤研究(A)"},
    "B": {"budget_per_grant": 8_000_000, "n_grants": 200, "success_rate": 0.25,
           "label": "基盤研究(B)"},
    "C": {"budget_per_grant": 3_000_000, "n_grants": 500, "success_rate": 0.30,
           "label": "基盤研究(C)"},
    "Wakate": {"budget_per_grant": 3_000_000, "n_grants": 300, "success_rate": 0.35,
               "label": "若手研究"},
    "Start": {"budget_per_grant": 1_500_000, "n_grants": 200, "success_rate": 0.40,
              "label": "研究活動スタート支援"},
}

# Japanese university categories
UNIVERSITY_TYPES = {
    "RU11": {"weight": 1.3, "pct": 0.15, "label": "RU11 (旧帝大等)"},
    "national_large": {"weight": 1.1, "pct": 0.20, "label": "大規模国立大"},
    "national_small": {"weight": 0.9, "pct": 0.25, "label": "中小規模国立大"},
    "private_large": {"weight": 0.85, "pct": 0.20, "label": "大規模私大"},
    "private_small": {"weight": 0.7, "pct": 0.15, "label": "中小規模私大"},
    "other": {"weight": 0.6, "pct": 0.05, "label": "その他機関"},
}


def generate_kakenhi_applicants(n=1000):
    """Generate synthetic KAKENHI applicant pool."""
    applicants = []

    for i in range(n):
        # University type (weighted sampling)
        uni_types = list(UNIVERSITY_TYPES.keys())
        uni_probs = [UNIVERSITY_TYPES[u]["pct"] for u in uni_types]
        uni_type = np.random.choice(uni_types, p=uni_probs)
        uni_weight = UNIVERSITY_TYPES[uni_type]["weight"]

        career_years = np.random.randint(1, 40)
        career_stage = "early" if career_years < 7 else ("mid" if career_years < 20 else "senior")

        # Research quality (influenced by university resources)
        base_talent = np.random.beta(2, 5)
        effective_quality = np.clip(base_talent * uni_weight + np.random.normal(0, 0.05), 0, 1)

        h_index = max(0, int(effective_quality * career_years * 2))
        papers = max(1, int(effective_quality * career_years * 3))
        citations = max(0, int(papers * np.random.lognormal(1.5, 1.0)))

        gender = np.random.choice(["M", "F"], p=[0.78, 0.22])  # Japan ratio
        field = np.random.choice(
            ["理工系", "生命科学", "人文社会", "医学", "学際"],
            p=[0.30, 0.25, 0.20, 0.15, 0.10]
        )

        applicants.append({
            "id": i,
            "university_type": uni_type,
            "university_label": UNIVERSITY_TYPES[uni_type]["label"],
            "career_years": career_years,
            "career_stage": career_stage,
            "gender": gender,
            "field": field,
            "talent": round(base_talent, 4),
            "effective_quality": round(effective_quality, 4),
            "h_index": h_index,
            "total_papers": papers,
            "total_citations": citations,
        })

    return pd.DataFrame(applicants)


def simulate_kakenhi_review(applicants_df, category="B", n_runs=10):
    """Simulate KAKENHI peer review process with reviewer variability."""
    cat_info = KAKENHI_CATEGORIES[category]
    n_grants = cat_info["n_grants"]
    n_applicants = len(applicants_df)

    results_over_runs = []

    for run in range(n_runs):
        scores = []
        for _, app in applicants_df.iterrows():
            true_q = app["effective_quality"]

            # 3 reviewers with noise and bias
            reviewer_scores = []
            for _ in range(3):
                # Bias: institutional prestige effect
                inst_bias = 0
                if app["university_type"] == "RU11":
                    inst_bias = 0.05
                elif app["university_type"] in ["national_large"]:
                    inst_bias = 0.02

                # Senior researcher bias
                seniority_bias = min(0.05, app["career_years"] * 0.002)

                noise = np.random.normal(0, 0.12)
                score = np.clip(true_q + inst_bias + seniority_bias + noise, 0, 1)
                reviewer_scores.append(score)

            scores.append(np.mean(reviewer_scores))

        applicants_df_run = applicants_df.copy()
        applicants_df_run["review_score"] = scores
        applicants_df_run = applicants_df_run.sort_values("review_score", ascending=False)

        funded = applicants_df_run.head(min(n_grants, n_applicants))

        run_result = {
            "run": run,
            "avg_quality_funded": round(funded["effective_quality"].mean(), 4),
            "avg_quality_all": round(applicants_df["effective_quality"].mean(), 4),
            "pct_RU11": round((funded["university_type"] == "RU11").mean() * 100, 1),
            "pct_female": round((funded["gender"] == "F").mean() * 100, 1),
            "pct_early": round((funded["career_stage"] == "early").mean() * 100, 1),
            "avg_hindex_funded": round(funded["h_index"].mean(), 1),
            "score_variance": round(np.var(scores), 4),
        }

        # Reviewer agreement (correlation between runs)
        run_result["funded_ids"] = funded["id"].tolist()[:n_grants]
        results_over_runs.append(run_result)

    # Inter-run agreement (Jaccard similarity)
    agreements = []
    for i in range(len(results_over_runs)):
        for j in range(i + 1, len(results_over_runs)):
            set_i = set(results_over_runs[i]["funded_ids"])
            set_j = set(results_over_runs[j]["funded_ids"])
            jaccard = len(set_i & set_j) / max(1, len(set_i | set_j))
            agreements.append(jaccard)

    return results_over_runs, np.mean(agreements)


def plot_kakenhi_results(applicants_df, review_results, agreement,
                         output_dir="figures"):
    """Plot KAKENHI case study results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # University type distribution
    ax = axes[0, 0]
    uni_dist_all = applicants_df["university_type"].value_counts(normalize=True)
    funded_ids = set(review_results[0]["funded_ids"])
    funded_df = applicants_df[applicants_df["id"].isin(funded_ids)]
    uni_dist_funded = funded_df["university_type"].value_counts(normalize=True)

    x = np.arange(len(uni_dist_all))
    w = 0.35
    label_map = {"RU11": "RU11", "national_large": "Nat-L", "national_small": "Nat-S",
                 "private_large": "Prv-L", "private_small": "Prv-S", "other": "Other"}
    labels = [label_map.get(u, u) for u in uni_dist_all.index]
    ax.bar(x - w/2, uni_dist_all.values * 100, w, label="All Applicants", color="#BDBDBD")
    funded_vals = [uni_dist_funded.get(u, 0) * 100 for u in uni_dist_all.index]
    ax.bar(x + w/2, funded_vals, w, label="Funded", color="#2196F3")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_title("University Type Distribution", fontsize=13)
    ax.set_ylabel("Percentage (%)", fontsize=11)
    ax.legend()

    # Quality distribution
    ax = axes[0, 1]
    qualities = [r["avg_quality_funded"] for r in review_results]
    ax.hist(qualities, bins=15, color="#4CAF50", alpha=0.7, edgecolor="black")
    ax.axvline(x=np.mean(qualities), color="red", linestyle="--",
               label=f"Mean = {np.mean(qualities):.3f}")
    ax.set_title(f"Funded Quality Distribution (N={len(review_results)} runs)", fontsize=13)
    ax.set_xlabel("Average Quality of Funded Researchers", fontsize=11)
    ax.set_ylabel("Frequency", fontsize=11)
    ax.legend()

    # Reviewer agreement
    ax = axes[1, 0]
    categories = list(KAKENHI_CATEGORIES.keys())
    # Simulate agreement for each category
    agreements_by_cat = {}
    for cat in categories:
        n_g = min(KAKENHI_CATEGORIES[cat]["n_grants"], len(applicants_df))
        _, agr = simulate_kakenhi_review(applicants_df, cat, n_runs=5)
        agreements_by_cat[cat] = agr

    cat_labels = [c for c in categories]
    agr_values = [agreements_by_cat[c] for c in categories]
    bars = ax.bar(cat_labels, agr_values, color="#FF9800", alpha=0.8)
    ax.set_title("Reviewer Agreement (Jaccard) by Category", fontsize=13)
    ax.set_ylabel("Jaccard Similarity", fontsize=11)
    ax.tick_params(axis="x", rotation=30)
    for bar, val in zip(bars, agr_values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{val:.2f}", ha="center", fontsize=9)

    # Gender and career stage
    ax = axes[1, 1]
    pct_female = [r["pct_female"] for r in review_results]
    pct_early = [r["pct_early"] for r in review_results]
    runs = range(len(review_results))
    ax.plot(list(runs), pct_female, "o-", label="Female (%)", color="#E91E63", markersize=5)
    ax.plot(list(runs), pct_early, "s-", label="Early Career (%)", color="#4CAF50", markersize=5)
    ax.axhline(y=22, color="#E91E63", linestyle="--", alpha=0.5, label="Population Female %")
    ax.set_title("Diversity in Funded Researchers Across Runs", fontsize=13)
    ax.set_xlabel("Simulation Run", fontsize=11)
    ax.set_ylabel("Percentage (%)", fontsize=11)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(f"{output_dir}/fig9_kakenhi_case_study.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    return agreements_by_cat


def run_kakenhi_case_study(output_dir="figures", results_dir="results", data_dir="data"):
    """Run full KAKENHI case study."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    applicants = generate_kakenhi_applicants(n=1000)
    applicants.to_csv(f"{data_dir}/kakenhi_applicants.csv", index=False)

    # Simulate category B (most common)
    review_results, agreement = simulate_kakenhi_review(applicants, category="B", n_runs=10)

    agreements_by_cat = plot_kakenhi_results(applicants, review_results, agreement, output_dir)

    # Compile results
    case_results = {
        "n_applicants": len(applicants),
        "category_B_results": {
            "n_runs": len(review_results),
            "mean_quality_funded": round(np.mean([r["avg_quality_funded"] for r in review_results]), 4),
            "std_quality_funded": round(np.std([r["avg_quality_funded"] for r in review_results]), 4),
            "mean_pct_RU11": round(np.mean([r["pct_RU11"] for r in review_results]), 1),
            "mean_pct_female": round(np.mean([r["pct_female"] for r in review_results]), 1),
            "mean_pct_early": round(np.mean([r["pct_early"] for r in review_results]), 1),
            "reviewer_agreement_jaccard": round(agreement, 3),
        },
        "agreements_by_category": {k: round(v, 3) for k, v in agreements_by_cat.items()},
        "institutional_bias": {
            "RU11_pct_applicants": round((applicants["university_type"] == "RU11").mean() * 100, 1),
            "RU11_pct_funded_mean": round(
                np.mean([r["pct_RU11"] for r in review_results]), 1),
        },
        "gender_gap": {
            "female_pct_applicants": round((applicants["gender"] == "F").mean() * 100, 1),
            "female_pct_funded_mean": round(
                np.mean([r["pct_female"] for r in review_results]), 1),
        },
    }

    with open(f"{results_dir}/kakenhi_case_study.json", "w") as f:
        json.dump(case_results, f, indent=2, default=str)

    print("KAKENHI case study complete.")
    print(f"  Applicants: {case_results['n_applicants']}")
    print(f"  Category B agreement: {case_results['category_B_results']['reviewer_agreement_jaccard']:.3f}")
    print(f"  RU11 representation: {case_results['institutional_bias']['RU11_pct_applicants']:.1f}% → "
          f"{case_results['institutional_bias']['RU11_pct_funded_mean']:.1f}% (funded)")
    print(f"  Female representation: {case_results['gender_gap']['female_pct_applicants']:.1f}% → "
          f"{case_results['gender_gap']['female_pct_funded_mean']:.1f}% (funded)")

    return case_results


if __name__ == "__main__":
    run_kakenhi_case_study()
