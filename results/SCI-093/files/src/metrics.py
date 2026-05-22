"""
Module 2: Research Output Measurement Metrics
Compares traditional metrics with alternative metrics and analyzes their limitations.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import json
import os

np.random.seed(42)


def generate_researcher_profiles(n=200):
    """Generate synthetic researcher profiles with various metrics."""
    profiles = []
    for i in range(n):
        career_years = np.random.randint(1, 35)
        base_productivity = np.random.lognormal(1.5, 0.8)

        total_papers = max(1, int(base_productivity * career_years * 0.3))
        total_citations = max(0, int(total_papers * np.random.lognormal(2.0, 1.0)))
        h_index = min(total_papers, int(np.sqrt(total_citations) * 0.5))

        avg_jif = max(0.5, np.random.lognormal(0.8, 0.6))
        jif_weighted = total_papers * avg_jif

        altmetric = max(0, int(np.random.exponential(20) + total_citations * 0.01))

        field_avg_citations = max(1, np.random.lognormal(2.5, 0.5))
        fwci = (total_citations / max(1, total_papers)) / field_avg_citations

        n_fields = min(5, max(1, int(np.random.exponential(1.5))))
        interdisciplinary = n_fields / 5.0

        mentoring = max(0, int(career_years * np.random.exponential(0.3)))

        open_access_ratio = np.random.beta(2, 3)
        data_sharing_ratio = np.random.beta(1.5, 4)
        open_science = (open_access_ratio + data_sharing_ratio) / 2

        patent_count = max(0, int(np.random.exponential(0.5)))
        policy_citations = max(0, int(np.random.exponential(0.3)))
        societal_impact = (patent_count * 2 + policy_citations * 3) / 10.0

        n_countries = min(20, max(1, int(np.random.exponential(2))))
        collab_diversity = n_countries / 20.0

        gender = np.random.choice(["M", "F", "Other"], p=[0.60, 0.35, 0.05])
        region = np.random.choice(
            ["Asia", "Europe", "NorthAmerica", "LatinAmerica", "Africa"],
            p=[0.30, 0.30, 0.25, 0.10, 0.05]
        )
        field = np.random.choice(["Physics", "Biology", "CS", "Chemistry", "Medicine"])
        career_stage = "early" if career_years < 7 else ("mid" if career_years < 20 else "senior")

        profiles.append({
            "id": i, "career_years": career_years, "career_stage": career_stage,
            "gender": gender, "region": region, "field": field,
            "total_papers": total_papers, "total_citations": total_citations,
            "h_index": h_index, "avg_jif": round(avg_jif, 2),
            "jif_weighted": round(jif_weighted, 2),
            "altmetric": altmetric, "fwci": round(fwci, 3),
            "interdisciplinary": round(interdisciplinary, 3),
            "mentoring": mentoring, "open_science": round(open_science, 3),
            "societal_impact": round(societal_impact, 3),
            "collab_diversity": round(collab_diversity, 3),
        })

    return pd.DataFrame(profiles)


def compute_composite_scores(df):
    """Compute composite scores combining traditional and alternative metrics."""
    metrics_trad = ["h_index", "total_citations", "jif_weighted"]
    metrics_alt = ["fwci", "altmetric", "interdisciplinary", "mentoring",
                   "open_science", "societal_impact", "collab_diversity"]

    df_norm = df.copy()
    for col in metrics_trad + metrics_alt:
        vmin, vmax = df[col].min(), df[col].max()
        if vmax > vmin:
            df_norm[f"{col}_norm"] = (df[col] - vmin) / (vmax - vmin)
        else:
            df_norm[f"{col}_norm"] = 0.0

    df_norm["score_traditional"] = (
        df_norm["h_index_norm"] * 0.4 +
        df_norm["total_citations_norm"] * 0.4 +
        df_norm["jif_weighted_norm"] * 0.2
    )

    df_norm["score_alternative"] = (
        df_norm["fwci_norm"] * 0.20 +
        df_norm["altmetric_norm"] * 0.10 +
        df_norm["interdisciplinary_norm"] * 0.15 +
        df_norm["mentoring_norm"] * 0.15 +
        df_norm["open_science_norm"] * 0.15 +
        df_norm["societal_impact_norm"] * 0.10 +
        df_norm["collab_diversity_norm"] * 0.15
    )

    df_norm["score_hybrid"] = (
        df_norm["score_traditional"] * 0.4 +
        df_norm["score_alternative"] * 0.6
    )

    return df_norm


def analyze_metric_biases(df):
    """Analyze biases in traditional vs alternative metrics."""
    results = {}

    for metric in ["score_traditional", "score_alternative", "score_hybrid"]:
        bias_analysis = {}
        bias_analysis["gender_means"] = df.groupby("gender")[metric].mean().to_dict()
        bias_analysis["region_means"] = df.groupby("region")[metric].mean().to_dict()
        bias_analysis["career_stage_means"] = df.groupby("career_stage")[metric].mean().to_dict()
        bias_analysis["field_means"] = df.groupby("field")[metric].mean().to_dict()
        corr = df[metric].corr(df["career_years"])
        bias_analysis["career_years_correlation"] = round(corr, 3)
        results[metric] = bias_analysis

    return results


def plot_metric_comparison(df, output_dir="figures"):
    """Plot comparison of traditional vs alternative scores."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    ax = axes[0, 0]
    scatter = ax.scatter(df["score_traditional"], df["score_alternative"],
                         c=df["career_years"], cmap="viridis", alpha=0.6, s=20)
    ax.set_xlabel("Traditional Score", fontsize=11)
    ax.set_ylabel("Alternative Score", fontsize=11)
    ax.set_title("Traditional vs Alternative Metrics", fontsize=13)
    plt.colorbar(scatter, ax=ax, label="Career Years")
    corr = df["score_traditional"].corr(df["score_alternative"])
    ax.text(0.05, 0.95, f"r = {corr:.3f}", transform=ax.transAxes,
            fontsize=11, verticalalignment="top")

    ax = axes[0, 1]
    metrics = ["score_traditional", "score_alternative", "score_hybrid"]
    gender_data = df.groupby("gender")[metrics].mean()
    gender_data.plot(kind="bar", ax=ax, width=0.7)
    ax.set_title("Scores by Gender", fontsize=13)
    ax.set_ylabel("Mean Score", fontsize=11)
    ax.tick_params(axis="x", rotation=0)
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    stage_order = ["early", "mid", "senior"]
    stage_data = df.groupby("career_stage")[metrics].mean().reindex(stage_order)
    stage_data.plot(kind="bar", ax=ax, width=0.7)
    ax.set_title("Scores by Career Stage", fontsize=13)
    ax.set_ylabel("Mean Score", fontsize=11)
    ax.tick_params(axis="x", rotation=0)
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    region_data = df.groupby("region")[metrics].mean()
    region_data.plot(kind="bar", ax=ax, width=0.7)
    ax.set_title("Scores by Region", fontsize=13)
    ax.set_ylabel("Mean Score", fontsize=11)
    ax.tick_params(axis="x", rotation=15)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(f"{output_dir}/fig5_metric_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_metrics_analysis(output_dir="figures", results_dir="results", data_dir="data"):
    """Run full metrics analysis pipeline."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    df = generate_researcher_profiles(n=200)
    df = compute_composite_scores(df)
    biases = analyze_metric_biases(df)

    df.to_csv(f"{data_dir}/researcher_profiles.csv", index=False)

    with open(f"{results_dir}/metric_biases.json", "w") as f:
        json.dump(biases, f, indent=2, default=str)

    plot_metric_comparison(df, output_dir)

    corr_trad_alt = df["score_traditional"].corr(df["score_alternative"])
    print("Metrics analysis complete.")
    print(f"  Correlation (traditional vs alternative): {corr_trad_alt:.3f}")
    print(f"  Traditional career_years corr: {biases['score_traditional']['career_years_correlation']}")
    print(f"  Alternative career_years corr: {biases['score_alternative']['career_years_correlation']}")

    return df, biases


if __name__ == "__main__":
    run_metrics_analysis()
