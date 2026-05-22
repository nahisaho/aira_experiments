#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
N_PROJECTS = 300
N_VOLUNTEER_RECORDS = 5000
DPI = 300
PALETTE = [
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
    "#F0E442",
    "#000000",
]

PLATFORMS = ["Zooniverse", "iNaturalist", "Galaxy_Zoo", "Foldit", "SciStarter", "eBird"]
FIELDS = ["Ecology", "Astronomy", "Biochemistry", "Climate", "Public_Health", "Biodiversity"]
EDUCATION_LEVELS = ["Secondary", "Undergraduate", "Graduate", "Self-taught", "Teacher"]

PLATFORM_PROFILES = {
    "Zooniverse": {"volunteer_scale": 3200, "retention": 0.43, "media": 8, "education": 10, "geo": 42},
    "iNaturalist": {"volunteer_scale": 4100, "retention": 0.52, "media": 11, "education": 12, "geo": 58},
    "Galaxy_Zoo": {"volunteer_scale": 2800, "retention": 0.38, "media": 13, "education": 9, "geo": 35},
    "Foldit": {"volunteer_scale": 1400, "retention": 0.47, "media": 10, "education": 8, "geo": 26},
    "SciStarter": {"volunteer_scale": 1900, "retention": 0.41, "media": 7, "education": 14, "geo": 31},
    "eBird": {"volunteer_scale": 5200, "retention": 0.63, "media": 9, "education": 11, "geo": 72},
}

FIELD_PROFILES = {
    "Ecology": {"contrib": 1.10, "quality": 0.82, "pubs": 1.15, "cost": 0.26, "pro_quality": 0.90},
    "Astronomy": {"contrib": 1.25, "quality": 0.84, "pubs": 1.05, "cost": 0.24, "pro_quality": 0.91},
    "Biochemistry": {"contrib": 0.85, "quality": 0.79, "pubs": 0.92, "cost": 0.34, "pro_quality": 0.93},
    "Climate": {"contrib": 1.05, "quality": 0.81, "pubs": 1.02, "cost": 0.29, "pro_quality": 0.89},
    "Public_Health": {"contrib": 0.78, "quality": 0.77, "pubs": 0.88, "cost": 0.37, "pro_quality": 0.92},
    "Biodiversity": {"contrib": 1.18, "quality": 0.83, "pubs": 1.10, "cost": 0.22, "pro_quality": 0.90},
}

EDUCATION_BASELINE = {
    "Secondary": 49,
    "Undergraduate": 60,
    "Graduate": 72,
    "Self-taught": 56,
    "Teacher": 68,
}


def ensure_output_dirs(root: Path) -> dict[str, Path]:
    figures_dir = root / "figures"
    results_dir = root / "results"
    figures_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    return {
        "figures": figures_dir,
        "results": results_dir,
    }



def mean_ci(series: pd.Series, confidence: float = 0.95) -> tuple[float, float, float]:
    values = pd.Series(series, dtype=float).dropna()
    if len(values) == 0:
        return float("nan"), float("nan"), float("nan")
    mean = float(values.mean())
    if len(values) == 1:
        return mean, mean, mean
    sem = stats.sem(values)
    margin = float(stats.t.ppf((1 + confidence) / 2, len(values) - 1) * sem)
    return mean, mean - margin, mean + margin



def benjamini_hochberg(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty_like(ranked)
    n = len(ranked)
    for i, value in enumerate(ranked, start=1):
        adjusted[i - 1] = value * n / i
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out.tolist()



def cohens_dz(differences: pd.Series) -> float:
    diffs = pd.Series(differences, dtype=float).dropna()
    if len(diffs) < 2:
        return 0.0
    sd = float(diffs.std(ddof=1))
    return 0.0 if sd == 0 else float(diffs.mean() / sd)



def paired_mean_difference_summary(before: pd.Series, after: pd.Series, label: str) -> dict[str, float | str]:
    diffs = pd.Series(after, dtype=float) - pd.Series(before, dtype=float)
    ci_mean, ci_low, ci_high = mean_ci(diffs)
    sample = diffs.iloc[: min(500, len(diffs))]
    shapiro_stat, shapiro_p = stats.shapiro(sample) if len(sample) >= 3 else (np.nan, np.nan)
    t_stat, t_p = stats.ttest_rel(after, before)
    wilcoxon_stat, wilcoxon_p = stats.wilcoxon(after, before, zero_method="wilcox")
    return {
        "comparison": label,
        "n": int(len(diffs)),
        "mean_difference": float(diffs.mean()),
        "ci95_low": float(ci_low),
        "ci95_high": float(ci_high),
        "cohens_dz": cohens_dz(diffs),
        "paired_t_stat": float(t_stat),
        "paired_t_p": float(t_p),
        "wilcoxon_stat": float(wilcoxon_stat),
        "wilcoxon_p": float(wilcoxon_p),
        "shapiro_w": float(shapiro_stat),
        "shapiro_p": float(shapiro_p),
    }



def simulate_projects(rng: np.random.Generator, n_projects: int) -> pd.DataFrame:
    rows: list[dict[str, int | float | str]] = []
    for idx in range(1, n_projects + 1):
        platform = rng.choice(PLATFORMS)
        field = rng.choice(FIELDS)
        profile = PLATFORM_PROFILES[platform]
        field_profile = FIELD_PROFILES[field]

        start_year = int(rng.integers(2007, 2025))
        maturity = 1.0 + max(0, 2024 - start_year) / 25
        volunteers = int(np.clip(rng.lognormal(np.log(profile["volunteer_scale"]), 0.55) / maturity, 150, 22000))
        contributions = int(
            np.clip(
                volunteers * rng.uniform(3.5, 11.5) * field_profile["contrib"] * rng.uniform(0.9, 1.2),
                1000,
                350000,
            )
        )
        geographic_reach = int(np.clip(rng.normal(profile["geo"], 10), 5, 120))
        retention = float(
            np.clip(
                rng.normal(profile["retention"] + (2024 - start_year) * 0.003, 0.08),
                0.12,
                0.88,
            )
        )
        data_quality = float(
            np.clip(
                rng.normal(
                    field_profile["quality"]
                    + retention * 0.10
                    + np.log10(volunteers) * 0.01,
                    0.04,
                ),
                0.55,
                0.99,
            )
        )
        publications = int(
            np.clip(
                rng.poisson(np.log1p(contributions) * field_profile["pubs"] * data_quality / 1.8),
                0,
                40,
            )
        )
        media_coverage = float(
            np.clip(
                rng.normal(35 + profile["media"] + publications * 2.3 + geographic_reach * 0.18, 8),
                5,
                100,
            )
        )
        education_impact = float(
            np.clip(
                rng.normal(42 + profile["education"] + retention * 24 + publications * 0.5, 7),
                10,
                100,
            )
        )
        cost_ratio = float(
            np.clip(
                field_profile["cost"]
                * (1.16 - data_quality * 0.30)
                * rng.normal(1.0, 0.12)
                / (1 + np.log10(contributions) * 0.14),
                0.03,
                0.80,
            )
        )

        rows.append(
            {
                "project_id": f"P{idx:03d}",
                "platform": platform,
                "field": field,
                "start_year": start_year,
                "num_volunteers": volunteers,
                "num_contributions": contributions,
                "data_quality_score": round(data_quality, 3),
                "publications_count": publications,
                "media_coverage_score": round(media_coverage, 2),
                "education_impact_score": round(education_impact, 2),
                "geographic_reach": geographic_reach,
                "volunteer_retention_rate": round(retention, 3),
                "cost_per_datapoint_vs_professional": round(cost_ratio, 3),
            }
        )
    return pd.DataFrame(rows)



def simulate_volunteer_records(
    rng: np.random.Generator,
    projects: pd.DataFrame,
    n_records: int,
) -> pd.DataFrame:
    project_weights = projects["num_volunteers"].to_numpy(dtype=float)
    project_weights = project_weights / project_weights.sum()
    chosen_idx = rng.choice(projects.index.to_numpy(), size=n_records, p=project_weights)

    records: list[dict[str, int | float | str]] = []
    for idx, project_idx in enumerate(chosen_idx, start=1):
        project = projects.loc[project_idx]
        avg_contrib = project["num_contributions"] / project["num_volunteers"]
        duration = int(
            np.clip(
                rng.gamma(shape=1.6 + project["volunteer_retention_rate"] * 5, scale=2.6),
                1,
                36,
            )
        )
        contributions = int(np.clip(rng.gamma(shape=1.8, scale=max(avg_contrib, 1.0)), 1, 1500))
        education = rng.choice(EDUCATION_LEVELS, p=[0.26, 0.32, 0.17, 0.16, 0.09])
        literacy_pre = float(np.clip(rng.normal(EDUCATION_BASELINE[education], 8), 20, 95))
        literacy_gain = (
            1.4
            + duration * 0.32
            + np.sqrt(contributions) * 0.14
            + (project["education_impact_score"] - 50) * 0.08
            + rng.normal(0, 2.6)
        )
        literacy_post = float(np.clip(literacy_pre + literacy_gain, 25, 100))
        records.append(
            {
                "volunteer_id": f"V{idx:05d}",
                "project_id": project["project_id"],
                "contributions": contributions,
                "duration_months": duration,
                "prior_education": education,
                "science_literacy_pre": round(literacy_pre, 2),
                "science_literacy_post": round(literacy_post, 2),
            }
        )
    return pd.DataFrame(records)



def add_derived_metrics(projects: pd.DataFrame, volunteer_records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    projects = projects.copy()
    volunteer_records = volunteer_records.copy()
    projects["professional_quality_benchmark"] = projects["field"].map(
        {field: profile["pro_quality"] for field, profile in FIELD_PROFILES.items()}
    )
    z_media = stats.zscore(projects["media_coverage_score"], ddof=1)
    z_pubs = stats.zscore(projects["publications_count"], ddof=1)
    z_geo = stats.zscore(projects["geographic_reach"], ddof=1)
    z_edu = stats.zscore(projects["education_impact_score"], ddof=1)
    projects["altmetric_attention_score"] = np.clip(50 + 12 * (0.45 * z_media + 0.25 * z_pubs + 0.15 * z_geo + 0.15 * z_edu), 0, 100)
    volunteer_records = volunteer_records.merge(
        projects[["project_id", "platform", "field", "volunteer_retention_rate"]],
        on="project_id",
        how="left",
    )
    volunteer_records["literacy_gain"] = (
        volunteer_records["science_literacy_post"] - volunteer_records["science_literacy_pre"]
    )
    volunteer_records["retained_12m"] = volunteer_records["duration_months"] >= 12
    return projects, volunteer_records



def analyze_platforms(projects: pd.DataFrame) -> dict[str, object]:
    grouped = projects.groupby("platform")
    summary = grouped.agg(
        projects=("project_id", "count"),
        mean_volunteers=("num_volunteers", "mean"),
        mean_contributions=("num_contributions", "mean"),
        mean_retention=("volunteer_retention_rate", "mean"),
        mean_quality=("data_quality_score", "mean"),
        mean_publications=("publications_count", "mean"),
    ).round(3)

    ci_records = []
    for platform, group in grouped:
        v_mean, v_low, v_high = mean_ci(group["num_volunteers"])
        c_mean, c_low, c_high = mean_ci(group["num_contributions"])
        r_mean, r_low, r_high = mean_ci(group["volunteer_retention_rate"])
        ci_records.append(
            {
                "platform": platform,
                "volunteers_ci95": [v_low, v_high],
                "contributions_ci95": [c_low, c_high],
                "retention_ci95": [r_low, r_high],
                "mean_volunteers": v_mean,
                "mean_contributions": c_mean,
                "mean_retention": r_mean,
            }
        )

    tests = []
    raw_p = []
    metrics = {
        "num_volunteers": "Participation volume",
        "num_contributions": "Contribution output",
        "volunteer_retention_rate": "Retention rate",
    }
    for metric, label in metrics.items():
        groups = [projects.loc[projects["platform"] == platform, metric].to_numpy() for platform in PLATFORMS]
        f_stat, p_value = stats.f_oneway(*groups)
        overall = np.concatenate(groups)
        ss_between = sum(len(group) * (group.mean() - overall.mean()) ** 2 for group in groups)
        ss_total = float(((overall - overall.mean()) ** 2).sum())
        eta_squared = 0.0 if ss_total == 0 else float(ss_between / ss_total)
        raw_p.append(float(p_value))
        tests.append({
            "metric": metric,
            "label": label,
            "f_stat": float(f_stat),
            "p_value": float(p_value),
            "eta_squared": eta_squared,
        })
    adjusted = benjamini_hochberg(raw_p)
    for test, adj in zip(tests, adjusted):
        test["fdr_p_value"] = float(adj)

    return {
        "summary": summary.reset_index().to_dict(orient="records"),
        "confidence_intervals": ci_records,
        "anova_tests": tests,
    }



def analyze_volunteers(volunteer_records: pd.DataFrame) -> dict[str, object]:
    education_summary = volunteer_records.groupby("prior_education").agg(
        volunteers=("volunteer_id", "count"),
        mean_duration_months=("duration_months", "mean"),
        mean_contributions=("contributions", "mean"),
        mean_literacy_gain=("literacy_gain", "mean"),
        retained_12m_rate=("retained_12m", "mean"),
    ).round(3)

    retention_by_platform = volunteer_records.groupby("platform").agg(
        mean_duration_months=("duration_months", "mean"),
        retained_12m_rate=("retained_12m", "mean"),
        mean_contributions=("contributions", "mean"),
    ).round(3)

    rho_duration, p_duration = stats.spearmanr(volunteer_records["duration_months"], volunteer_records["contributions"])

    return {
        "education_summary": education_summary.reset_index().to_dict(orient="records"),
        "retention_by_platform": retention_by_platform.reset_index().to_dict(orient="records"),
        "duration_contribution_correlation": {
            "spearman_rho": float(rho_duration),
            "p_value": float(p_duration),
        },
    }



def analyze_scientific_output(projects: pd.DataFrame) -> dict[str, object]:
    comparison = paired_mean_difference_summary(
        projects["professional_quality_benchmark"],
        projects["data_quality_score"],
        "Citizen science quality minus professional benchmark",
    )
    comparison["citizen_quality_mean"] = float(projects["data_quality_score"].mean())
    comparison["professional_quality_mean"] = float(projects["professional_quality_benchmark"].mean())
    comparison["mean_quality_gap"] = float(
        (projects["data_quality_score"] - projects["professional_quality_benchmark"]).mean()
    )

    field_gaps = projects.assign(
        quality_gap=projects["data_quality_score"] - projects["professional_quality_benchmark"]
    ).groupby("field").agg(
        mean_quality_gap=("quality_gap", "mean"),
        mean_publications=("publications_count", "mean"),
    ).round(3)

    return {
        "overall": comparison,
        "field_quality_gaps": field_gaps.reset_index().to_dict(orient="records"),
    }



def analyze_cost_effectiveness(projects: pd.DataFrame) -> dict[str, object]:
    rho_quality, p_quality = stats.spearmanr(
        projects["cost_per_datapoint_vs_professional"], projects["data_quality_score"]
    )
    rho_output, p_output = stats.spearmanr(
        projects["cost_per_datapoint_vs_professional"], projects["publications_count"]
    )
    by_field = projects.groupby("field").agg(
        mean_cost_ratio=("cost_per_datapoint_vs_professional", "mean"),
        mean_quality=("data_quality_score", "mean"),
        mean_publications=("publications_count", "mean"),
    ).round(3)
    return {
        "cost_quality_correlation": {"spearman_rho": float(rho_quality), "p_value": float(p_quality)},
        "cost_output_correlation": {"spearman_rho": float(rho_output), "p_value": float(p_output)},
        "by_field": by_field.reset_index().to_dict(orient="records"),
    }



def analyze_literacy(volunteer_records: pd.DataFrame) -> dict[str, object]:
    overall = paired_mean_difference_summary(
        volunteer_records["science_literacy_pre"],
        volunteer_records["science_literacy_post"],
        "Pre/post science literacy",
    )
    by_education = volunteer_records.groupby("prior_education").agg(
        mean_pre=("science_literacy_pre", "mean"),
        mean_post=("science_literacy_post", "mean"),
        mean_gain=("literacy_gain", "mean"),
    ).round(3)
    return {
        "overall": overall,
        "by_education": by_education.reset_index().to_dict(orient="records"),
    }



def analyze_geography(projects: pd.DataFrame) -> dict[str, object]:
    diversity = projects.groupby("platform").agg(
        mean_countries=("geographic_reach", "mean"),
        median_countries=("geographic_reach", "median"),
        max_countries=("geographic_reach", "max"),
    ).round(3)
    rho, p_value = stats.spearmanr(projects["geographic_reach"], projects["num_volunteers"])
    return {
        "platform_diversity": diversity.reset_index().to_dict(orient="records"),
        "reach_participation_correlation": {
            "spearman_rho": float(rho),
            "p_value": float(p_value),
        },
    }



def analyze_altmetrics(projects: pd.DataFrame) -> dict[str, object]:
    rho_media, p_media = stats.spearmanr(projects["altmetric_attention_score"], projects["media_coverage_score"])
    rho_pubs, p_pubs = stats.spearmanr(projects["altmetric_attention_score"], projects["publications_count"])
    return {
        "altmetric_media_correlation": {"spearman_rho": float(rho_media), "p_value": float(p_media)},
        "altmetric_publication_correlation": {"spearman_rho": float(rho_pubs), "p_value": float(p_pubs)},
        "altmetric_score_mean": float(projects["altmetric_attention_score"].mean()),
    }



def plot_platform_comparison(projects: pd.DataFrame, figure_path: Path) -> None:
    summary = projects.groupby("platform").agg(
        Volunteers=("num_volunteers", "mean"),
        Contributions=("num_contributions", "mean"),
        Retention=("volunteer_retention_rate", "mean"),
    ).reindex(PLATFORMS)
    normalized = summary.copy()
    normalized["Volunteers"] = 100 * normalized["Volunteers"] / normalized["Volunteers"].max()
    normalized["Contributions"] = 100 * normalized["Contributions"] / normalized["Contributions"].max()
    normalized["Retention"] = 100 * normalized["Retention"] / normalized["Retention"].max()

    x = np.arange(len(PLATFORMS))
    width = 0.24

    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    ax.bar(x - width, normalized["Volunteers"], width, label="Participation", color=PALETTE[0])
    ax.bar(x, normalized["Contributions"], width, label="Output", color=PALETTE[2])
    ax.bar(x + width, normalized["Retention"], width, label="Retention", color=PALETTE[1])
    ax.set_xticks(x)
    ax.set_xticklabels(PLATFORMS, rotation=20, ha="right")
    ax.set_ylabel("Relative platform metric (% of max)")
    ax.set_title("Citizen science platform comparison")
    ax.set_ylim(0, 110)
    ax.legend(frameon=False)
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)



def plot_retention_curves(volunteer_records: pd.DataFrame, figure_path: Path) -> None:
    max_month = 24
    months = np.arange(1, max_month + 1)

    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    for color, platform in zip(PALETTE, PLATFORMS):
        platform_records = volunteer_records.loc[volunteer_records["platform"] == platform, "duration_months"]
        retention_curve = [float((platform_records >= month).mean()) for month in months]
        ax.plot(months, retention_curve, label=platform, color=color, linewidth=2)

    ax.set_xlabel("Participation duration (months)")
    ax.set_ylabel("Share of volunteers retained")
    ax.set_title("Volunteer retention curves by platform")
    ax.set_ylim(0, 1.02)
    ax.set_xlim(1, max_month)
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)



def plot_literacy_impact(volunteer_records: pd.DataFrame, figure_path: Path) -> None:
    pre = volunteer_records["science_literacy_pre"].to_numpy(dtype=float)
    post = volunteer_records["science_literacy_post"].to_numpy(dtype=float)
    sample = volunteer_records.sample(n=min(250, len(volunteer_records)), random_state=SEED)

    fig, ax = plt.subplots(figsize=(7.8, 5.8))
    violin = ax.violinplot([pre, post], positions=[1, 2], widths=0.7, showmeans=True, showextrema=False)
    for body, color in zip(violin["bodies"], [PALETTE[0], PALETTE[2]]):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.55)
    violin["cmeans"].set_color(PALETTE[5])
    violin["cmeans"].set_linewidth(2)

    ax.scatter(np.full(len(sample), 1), sample["science_literacy_pre"], color=PALETTE[0], alpha=0.08, s=9)
    ax.scatter(np.full(len(sample), 2), sample["science_literacy_post"], color=PALETTE[2], alpha=0.08, s=9)
    for _, row in sample.iterrows():
        ax.plot([1, 2], [row["science_literacy_pre"], row["science_literacy_post"]], color="gray", alpha=0.05)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Pre participation", "Post participation"])
    ax.set_ylabel("Science literacy score")
    ax.set_title("Science literacy impact of participation")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)



def plot_cost_effectiveness(projects: pd.DataFrame, figure_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for idx, field in enumerate(FIELDS):
        subset = projects.loc[projects["field"] == field]
        ax.scatter(
            subset["cost_per_datapoint_vs_professional"],
            subset["data_quality_score"],
            s=36 + subset["publications_count"] * 2,
            alpha=0.75,
            color=PALETTE[idx],
            label=field.replace("_", " "),
            edgecolor="white",
            linewidth=0.4,
        )
    ax.set_xlabel("Cost per data point vs professional (ratio)")
    ax.set_ylabel("Data quality score")
    ax.set_title("Cost-effectiveness by field")
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)



def json_ready(obj):
    if isinstance(obj, dict):
        return {str(key): json_ready(value) for key, value in obj.items()}
    if isinstance(obj, list):
        return [json_ready(value) for value in obj]
    if isinstance(obj, tuple):
        return [json_ready(value) for value in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (pd.Timestamp, datetime)):
        return obj.isoformat()
    return obj



def main() -> None:
    rng = np.random.default_rng(SEED)
    root = Path(__file__).resolve().parent.parent
    output_dirs = ensure_output_dirs(root)

    projects = simulate_projects(rng, N_PROJECTS)
    volunteer_records = simulate_volunteer_records(rng, projects, N_VOLUNTEER_RECORDS)
    projects, volunteer_records = add_derived_metrics(projects, volunteer_records)

    figures = {
        "citizen_platform_comparison": output_dirs["figures"] / "citizen_platform_comparison.png",
        "citizen_retention_curve": output_dirs["figures"] / "citizen_retention_curve.png",
        "citizen_literacy_impact": output_dirs["figures"] / "citizen_literacy_impact.png",
        "citizen_cost_effectiveness": output_dirs["figures"] / "citizen_cost_effectiveness.png",
    }
    plot_platform_comparison(projects, figures["citizen_platform_comparison"])
    plot_retention_curves(volunteer_records, figures["citizen_retention_curve"])
    plot_literacy_impact(volunteer_records, figures["citizen_literacy_impact"])
    plot_cost_effectiveness(projects, figures["citizen_cost_effectiveness"])

    results = {
        "metadata": {
            "random_seed": SEED,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "n_projects": int(len(projects)),
            "n_volunteer_records": int(len(volunteer_records)),
            "platforms": PLATFORMS,
            "fields": FIELDS,
        },
        "platform_comparison": analyze_platforms(projects),
        "volunteer_demographics_and_retention": analyze_volunteers(volunteer_records),
        "scientific_output_quality": analyze_scientific_output(projects),
        "cost_effectiveness": analyze_cost_effectiveness(projects),
        "science_literacy_impact": analyze_literacy(volunteer_records),
        "geographic_diversity": analyze_geography(projects),
        "altmetrics_media_correlation": analyze_altmetrics(projects),
        "figure_files": {name: str(path) for name, path in figures.items()},
    }

    results_path = output_dirs["results"] / "citizen_science_results.json"
    with results_path.open("w", encoding="utf-8") as handle:
        json.dump(json_ready(results), handle, indent=2)

    print(f"Saved results to {results_path}")
    for figure_name, path in figures.items():
        print(f"Saved {figure_name} -> {path}")



if __name__ == "__main__":
    main()
