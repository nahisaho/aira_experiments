from __future__ import annotations

import json
import random
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
N_PAPERS = 1200

ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"
LOGS_DIR = ROOT / "logs"
REPORT_PATH = ROOT / "report.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"
STAT_SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
PREPROCESSING_LOG_PATH = DATA_DIR / "preprocessing-log.md"

OKABE_ITO = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#7F7F7F",
    "sky": "#56B4E9",
}

FIELD_REVIEW_MEAN = {
    "Biology": 112,
    "Neuroscience": 116,
    "Medicine": 132,
    "Public Health": 126,
    "Physics": 96,
    "Computer Science": 84,
    "Mathematics": 102,
    "Economics": 124,
    "Social Science": 120,
}

FIELD_ROUND_MEAN = {
    "Biology": 2.4,
    "Neuroscience": 2.5,
    "Medicine": 2.7,
    "Public Health": 2.5,
    "Physics": 2.2,
    "Computer Science": 2.0,
    "Mathematics": 2.1,
    "Economics": 2.3,
    "Social Science": 2.4,
}

FIELD_CITATION_1Y = {
    "Biology": 10.5,
    "Neuroscience": 12.0,
    "Medicine": 11.5,
    "Public Health": 9.5,
    "Physics": 8.0,
    "Computer Science": 9.0,
    "Mathematics": 6.5,
    "Economics": 7.0,
    "Social Science": 6.8,
}

SERVER_FIELD_PROBS = {
    "bioRxiv": {
        "Biology": 0.38,
        "Neuroscience": 0.24,
        "Medicine": 0.04,
        "Public Health": 0.03,
        "Physics": 0.03,
        "Computer Science": 0.06,
        "Mathematics": 0.02,
        "Economics": 0.02,
        "Social Science": 0.18,
    },
    "medRxiv": {
        "Biology": 0.04,
        "Neuroscience": 0.04,
        "Medicine": 0.47,
        "Public Health": 0.24,
        "Physics": 0.02,
        "Computer Science": 0.03,
        "Mathematics": 0.02,
        "Economics": 0.02,
        "Social Science": 0.12,
    },
    "arXiv": {
        "Biology": 0.02,
        "Neuroscience": 0.03,
        "Medicine": 0.01,
        "Public Health": 0.01,
        "Physics": 0.34,
        "Computer Science": 0.34,
        "Mathematics": 0.20,
        "Economics": 0.02,
        "Social Science": 0.03,
    },
    "SSRN": {
        "Biology": 0.01,
        "Neuroscience": 0.01,
        "Medicine": 0.03,
        "Public Health": 0.02,
        "Physics": 0.01,
        "Computer Science": 0.06,
        "Mathematics": 0.03,
        "Economics": 0.39,
        "Social Science": 0.44,
    },
    "none": {
        "Biology": 0.16,
        "Neuroscience": 0.10,
        "Medicine": 0.16,
        "Public Health": 0.10,
        "Physics": 0.10,
        "Computer Science": 0.12,
        "Mathematics": 0.08,
        "Economics": 0.08,
        "Social Science": 0.10,
    },
}

SERVER_PROBS = {
    "bioRxiv": 0.22,
    "medRxiv": 0.12,
    "arXiv": 0.24,
    "SSRN": 0.10,
    "none": 0.32,
}

SERVER_DOWNLOAD_BASE = {
    "bioRxiv": 1250,
    "medRxiv": 1550,
    "arXiv": 1850,
    "SSRN": 900,
}

SERVER_COMMENT_BASE = {
    "bioRxiv": 2.2,
    "medRxiv": 2.6,
    "arXiv": 1.8,
    "SSRN": 1.5,
}


def ensure_dirs() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)


def confidence_interval(values_a: pd.Series, values_b: pd.Series, rng: np.random.Generator, iterations: int = 2000) -> tuple[float, float]:
    a = values_a.to_numpy(dtype=float)
    b = values_b.to_numpy(dtype=float)
    boot = []
    for _ in range(iterations):
        sample_a = rng.choice(a, size=len(a), replace=True)
        sample_b = rng.choice(b, size=len(b), replace=True)
        boot.append(sample_a.mean() - sample_b.mean())
    return tuple(np.percentile(boot, [2.5, 97.5]))


def cohens_d(values_a: pd.Series, values_b: pd.Series) -> float:
    a = values_a.to_numpy(dtype=float)
    b = values_b.to_numpy(dtype=float)
    pooled_sd = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled_sd)


def shapiro_safe(values: pd.Series) -> dict[str, float | None]:
    sample = values.to_numpy(dtype=float)
    if len(sample) > 5000:
        sample = sample[:5000]
    if len(sample) < 3:
        return {"statistic": None, "p_value": None}
    stat, p = stats.shapiro(sample)
    return {"statistic": float(stat), "p_value": float(p)}


def choose_field(server: str, rng: np.random.Generator) -> str:
    field_probs = SERVER_FIELD_PROBS[server]
    fields = list(field_probs.keys())
    probs = list(field_probs.values())
    return str(rng.choice(fields, p=probs))


def simulate_dataset(n_papers: int = N_PAPERS, seed: int = SEED) -> pd.DataFrame:
    random.seed(seed)
    rng = np.random.default_rng(seed)
    servers = list(SERVER_PROBS.keys())
    server_probs = list(SERVER_PROBS.values())
    start_date = pd.Timestamp("2018-01-01")
    end_date = pd.Timestamp("2021-12-31")
    span_days = (end_date - start_date).days

    records: list[dict[str, object]] = []
    for idx in range(1, n_papers + 1):
        server = str(rng.choice(servers, p=server_probs))
        is_preprint = server != "none"
        field = choose_field(server, rng)
        quality = rng.normal(0, 0.9)
        visibility = rng.normal(0, 0.7)

        submission_date = start_date + pd.Timedelta(days=int(rng.integers(0, span_days + 1)))

        if is_preprint:
            num_versions = int(np.clip(1 + rng.poisson(1.35 + max(quality, 0) * 0.3), 1, 6))
            days_before_submission = int(np.clip(rng.normal(48 + 11 * (num_versions - 1), 18), 7, 210))
            preprint_date = submission_date - pd.Timedelta(days=days_before_submission)

            engagement_multiplier = np.exp(0.35 * quality + 0.25 * visibility)
            downloads = int(
                np.clip(
                    rng.lognormal(np.log(SERVER_DOWNLOAD_BASE[server]), 0.45)
                    * (1 + 0.12 * (num_versions - 1))
                    * engagement_multiplier,
                    120,
                    25000,
                )
            )
            comments_rate = SERVER_COMMENT_BASE[server] + 0.65 * (num_versions - 1) + 0.45 * max(quality, 0)
            comments = int(np.clip(rng.poisson(max(comments_rate, 0.2)), 0, 28))
        else:
            num_versions = 0
            preprint_date = pd.NaT
            downloads = 0
            comments = 0

        base_review = max(35, rng.normal(FIELD_REVIEW_MEAN[field], 14))
        review_multiplier = 1.0
        if is_preprint:
            feedback_gain = min(0.10, 0.008 * comments + 0.018 * max(num_versions - 1, 0))
            review_multiplier -= 0.10 + feedback_gain
        review_multiplier *= 1.0 - 0.02 * max(quality, 0)
        review_multiplier *= 1.0 + 0.03 * max(-quality, 0)
        review_duration_days = int(np.clip(rng.normal(base_review * review_multiplier, 8), 28, 240))

        review_rounds_signal = (
            FIELD_ROUND_MEAN[field]
            + rng.normal(0, 0.35)
            - (0.22 if is_preprint else 0.0)
            - (0.045 * comments if is_preprint else 0.0)
            - (0.08 * max(num_versions - 1, 0) if is_preprint else 0.0)
            - 0.10 * max(quality, 0)
            + 0.14 * max(-quality, 0)
        )
        peer_review_rounds = int(np.clip(np.rint(review_rounds_signal), 1, 5))

        revision_overhead = max(8, rng.normal(16 + 7 * (peer_review_rounds - 1), 4))
        acceptance_date = submission_date + pd.Timedelta(days=int(review_duration_days + revision_overhead))
        publication_lag = int(np.clip(rng.normal(38, 10), 12, 95))
        publication_date = acceptance_date + pd.Timedelta(days=publication_lag)

        citation_baseline = FIELD_CITATION_1Y[field] * np.exp(0.32 * quality + rng.normal(0, 0.26))
        early_citation_multiplier = 1.0
        if is_preprint:
            early_citation_multiplier += 0.18 + min(0.10, 0.008 * comments + 0.015 * (num_versions - 1))
            early_citation_multiplier += min(0.06, np.log1p(downloads) / 140)
        citation_count_1yr = int(np.clip(np.rint(citation_baseline * early_citation_multiplier), 0, 140))

        long_term_growth = 2.25 + 0.14 * max(quality, 0) + rng.normal(0, 0.18)
        citation_count_3yr = int(
            np.clip(
                np.rint(citation_count_1yr * max(long_term_growth, 1.55) + rng.normal(2.0, 3.5)),
                citation_count_1yr,
                360,
            )
        )

        altmetric_base = 4.0 + 1.15 * citation_count_1yr + (0.013 * downloads if is_preprint else 0.0) + 1.2 * comments
        altmetric_score = float(np.clip(rng.normal(altmetric_base, 7.5), 1, 320))
        media_rate = max(0.15, 0.12 * citation_count_1yr + 0.05 * comments + 0.018 * altmetric_score)
        media_mentions = int(np.clip(rng.poisson(media_rate), 0, 40))

        records.append(
            {
                "paper_id": f"P{idx:04d}",
                "preprint_server": server,
                "field": field,
                "submission_date": submission_date,
                "preprint_date": preprint_date,
                "acceptance_date": acceptance_date,
                "publication_date": publication_date,
                "num_preprint_versions": num_versions,
                "preprint_downloads": downloads,
                "preprint_comments": comments,
                "peer_review_rounds": peer_review_rounds,
                "review_duration_days": review_duration_days,
                "citation_count_1yr": citation_count_1yr,
                "citation_count_3yr": citation_count_3yr,
                "altmetric_score": round(altmetric_score, 2),
                "media_mentions": media_mentions,
            }
        )

    df = pd.DataFrame.from_records(records)
    for col in ["submission_date", "preprint_date", "acceptance_date", "publication_date"]:
        df[col] = pd.to_datetime(df[col])
    df["has_preprint"] = df["preprint_server"] != "none"
    df["submission_to_publication_days"] = (df["publication_date"] - df["submission_date"]).dt.days
    df["submission_to_acceptance_days"] = (df["acceptance_date"] - df["submission_date"]).dt.days
    df["days_preprint_before_submission"] = (df["submission_date"] - df["preprint_date"]).dt.days.fillna(0).astype(int)
    df["citation_acceleration_ratio"] = np.where(
        df["citation_count_3yr"] > 0,
        df["citation_count_1yr"] / df["citation_count_3yr"],
        0.0,
    )
    return df


def build_citation_trajectory(df: pd.DataFrame, months: int = 36) -> pd.DataFrame:
    rows = []
    for label, group in [("Preprint", df[df["has_preprint"]]), ("No preprint", df[~df["has_preprint"]])]:
        for month in range(1, months + 1):
            month_values = []
            for _, paper in group.iterrows():
                first_year = paper["citation_count_1yr"]
                third_year = paper["citation_count_3yr"]
                early_shape = 0.82 if paper["has_preprint"] else 1.02
                late_shape = 0.88 if paper["has_preprint"] else 1.0
                if month <= 12:
                    cum_cites = first_year * (month / 12) ** early_shape
                else:
                    added = (third_year - first_year) * ((month - 12) / 24) ** late_shape
                    cum_cites = first_year + added
                month_values.append(cum_cites)
            rows.append({"group": label, "month": month, "mean_citations": float(np.mean(month_values))})
    return pd.DataFrame(rows)


def analyze(df: pd.DataFrame) -> dict[str, object]:
    rng = np.random.default_rng(SEED)
    preprint = df[df["has_preprint"]]
    non_preprint = df[~df["has_preprint"]]

    duration_pre = preprint["submission_to_publication_days"]
    duration_non = non_preprint["submission_to_publication_days"]
    review_pre = preprint["review_duration_days"]
    review_non = non_preprint["review_duration_days"]
    citation_pre = preprint["citation_count_1yr"]
    citation_non = non_preprint["citation_count_1yr"]

    t_pub = stats.ttest_ind(duration_pre, duration_non, equal_var=False)
    t_cites = stats.ttest_ind(citation_pre, citation_non, equal_var=False)
    levene_pub = stats.levene(duration_pre, duration_non)

    high_feedback_threshold = float(preprint["preprint_comments"].median())
    high_feedback = preprint[preprint["preprint_comments"] >= high_feedback_threshold]
    low_feedback = preprint[preprint["preprint_comments"] < high_feedback_threshold]
    t_rounds = stats.ttest_ind(high_feedback["peer_review_rounds"], low_feedback["peer_review_rounds"], equal_var=False)
    comments_rounds_corr = stats.spearmanr(preprint["preprint_comments"], preprint["peer_review_rounds"])
    versions_rounds_corr = stats.spearmanr(preprint["num_preprint_versions"], preprint["peer_review_rounds"])
    versions_duration_corr = stats.spearmanr(preprint["num_preprint_versions"], preprint["review_duration_days"])
    comments_cites_corr = stats.spearmanr(preprint["preprint_comments"], preprint["citation_count_3yr"])
    downloads_cites_corr = stats.spearmanr(preprint["preprint_downloads"], preprint["citation_count_3yr"])

    server_df = preprint.groupby("preprint_server", as_index=False).agg(
        papers=("paper_id", "count"),
        mean_review_duration_days=("review_duration_days", "mean"),
        mean_review_rounds=("peer_review_rounds", "mean"),
        mean_citation_count_1yr=("citation_count_1yr", "mean"),
        mean_citation_count_3yr=("citation_count_3yr", "mean"),
        mean_downloads=("preprint_downloads", "mean"),
        mean_comments=("preprint_comments", "mean"),
        mean_altmetric_score=("altmetric_score", "mean"),
    )

    review_groups = [group["review_duration_days"].to_numpy() for _, group in preprint.groupby("preprint_server")]
    citation_groups = [group["citation_count_1yr"].to_numpy() for _, group in preprint.groupby("preprint_server")]
    anova_review = stats.f_oneway(*review_groups)
    anova_citations = stats.f_oneway(*citation_groups)
    corrected_pvals = {
        "review_duration_days": float(min(anova_review.pvalue * 2, 1.0)),
        "citation_count_1yr": float(min(anova_citations.pvalue * 2, 1.0)),
    }

    ci_pub = confidence_interval(duration_pre, duration_non, rng)
    ci_cites = confidence_interval(citation_pre, citation_non, rng)
    ci_rounds = confidence_interval(high_feedback["peer_review_rounds"], low_feedback["peer_review_rounds"], rng)

    review_reduction_pct = (1 - review_pre.mean() / review_non.mean()) * 100
    early_citation_advantage_pct = (citation_pre.mean() / citation_non.mean() - 1) * 100

    version_summary = (
        preprint.assign(
            version_bucket=pd.cut(
                preprint["num_preprint_versions"],
                bins=[0, 1, 2, 6],
                labels=["1 version", "2 versions", "3+ versions"],
            )
        )
        .groupby("version_bucket", observed=False)
        .agg(
            mean_comments=("preprint_comments", "mean"),
            mean_review_rounds=("peer_review_rounds", "mean"),
            mean_review_duration_days=("review_duration_days", "mean"),
            mean_citation_count_1yr=("citation_count_1yr", "mean"),
        )
        .round(2)
    )

    return {
        "seed": SEED,
        "n_papers": int(len(df)),
        "preprint_share": float(df["has_preprint"].mean()),
        "realized_effects": {
            "mean_review_duration_preprint_days": round(float(review_pre.mean()), 2),
            "mean_review_duration_non_preprint_days": round(float(review_non.mean()), 2),
            "review_time_reduction_pct": round(float(review_reduction_pct), 2),
            "mean_1yr_citations_preprint": round(float(citation_pre.mean()), 2),
            "mean_1yr_citations_non_preprint": round(float(citation_non.mean()), 2),
            "early_citation_advantage_pct": round(float(early_citation_advantage_pct), 2),
        },
        "assumption_checks": {
            "submission_to_publication_shapiro_preprint": shapiro_safe(duration_pre),
            "submission_to_publication_shapiro_non_preprint": shapiro_safe(duration_non),
            "submission_to_publication_levene": {
                "statistic": float(levene_pub.statistic),
                "p_value": float(levene_pub.pvalue),
            },
        },
        "submission_to_publication_analysis": {
            "preprint_mean_days": round(float(duration_pre.mean()), 2),
            "non_preprint_mean_days": round(float(duration_non.mean()), 2),
            "mean_difference_days": round(float(duration_pre.mean() - duration_non.mean()), 2),
            "mean_difference_95ci": [round(float(ci_pub[0]), 2), round(float(ci_pub[1]), 2)],
            "welch_t_statistic": round(float(t_pub.statistic), 4),
            "welch_p_value": float(t_pub.pvalue),
            "cohens_d": round(float(cohens_d(duration_pre, duration_non)), 4),
        },
        "feedback_efficiency_analysis": {
            "high_feedback_threshold_comments": high_feedback_threshold,
            "high_feedback_mean_rounds": round(float(high_feedback["peer_review_rounds"].mean()), 2),
            "low_feedback_mean_rounds": round(float(low_feedback["peer_review_rounds"].mean()), 2),
            "mean_difference_rounds": round(float(high_feedback["peer_review_rounds"].mean() - low_feedback["peer_review_rounds"].mean()), 2),
            "mean_difference_95ci": [round(float(ci_rounds[0]), 2), round(float(ci_rounds[1]), 2)],
            "welch_t_statistic": round(float(t_rounds.statistic), 4),
            "welch_p_value": float(t_rounds.pvalue),
            "spearman_comments_vs_rounds": {
                "rho": round(float(comments_rounds_corr.correlation), 4),
                "p_value": float(comments_rounds_corr.pvalue),
            },
        },
        "early_citation_analysis": {
            "preprint_mean_1yr_citations": round(float(citation_pre.mean()), 2),
            "non_preprint_mean_1yr_citations": round(float(citation_non.mean()), 2),
            "mean_difference": round(float(citation_pre.mean() - citation_non.mean()), 2),
            "mean_difference_95ci": [round(float(ci_cites[0]), 2), round(float(ci_cites[1]), 2)],
            "welch_t_statistic": round(float(t_cites.statistic), 4),
            "welch_p_value": float(t_cites.pvalue),
            "cohens_d": round(float(cohens_d(citation_pre, citation_non)), 4),
        },
        "version_evolution_analysis": {
            "spearman_versions_vs_review_rounds": {
                "rho": round(float(versions_rounds_corr.correlation), 4),
                "p_value": float(versions_rounds_corr.pvalue),
            },
            "spearman_versions_vs_review_duration": {
                "rho": round(float(versions_duration_corr.correlation), 4),
                "p_value": float(versions_duration_corr.pvalue),
            },
            "version_bucket_summary": version_summary.reset_index().to_dict(orient="records"),
        },
        "server_specific_trends": {
            "anova_review_duration_days": {
                "f_statistic": round(float(anova_review.statistic), 4),
                "p_value": float(anova_review.pvalue),
                "bonferroni_corrected_p_value": corrected_pvals["review_duration_days"],
            },
            "anova_citation_count_1yr": {
                "f_statistic": round(float(anova_citations.statistic), 4),
                "p_value": float(anova_citations.pvalue),
                "bonferroni_corrected_p_value": corrected_pvals["citation_count_1yr"],
            },
            "server_metrics": server_df.round(2).to_dict(orient="records"),
        },
        "engagement_correlation_analysis": {
            "spearman_downloads_vs_citations_3yr": {
                "rho": round(float(downloads_cites_corr.correlation), 4),
                "p_value": float(downloads_cites_corr.pvalue),
            },
            "spearman_comments_vs_citations_3yr": {
                "rho": round(float(comments_cites_corr.correlation), 4),
                "p_value": float(comments_cites_corr.pvalue),
            },
        },
    }


def save_figure(fig: plt.Figure, name: str) -> None:
    png_path = FIGURES_DIR / f"{name}.png"
    svg_path = FIGURES_DIR / f"{name}.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")


def save_figures(df: pd.DataFrame, trajectory_df: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, ax = plt.subplots(figsize=(8, 6))
    box_data = [
        df.loc[df["has_preprint"], "review_duration_days"],
        df.loc[~df["has_preprint"], "review_duration_days"],
    ]
    box = ax.boxplot(box_data, patch_artist=True, tick_labels=["Preprint", "No preprint"])
    for patch, color in zip(box["boxes"], [OKABE_ITO["blue"], OKABE_ITO["grey"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_title("Review duration by preprint status")
    ax.set_ylabel("Review duration (days)")
    ax.set_xlabel("Paper group")
    fig.tight_layout()
    save_figure(fig, "preprint_timeline")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    for label, color in [("Preprint", OKABE_ITO["blue"]), ("No preprint", OKABE_ITO["grey"] )]:
        subset = trajectory_df[trajectory_df["group"] == label]
        ax.plot(subset["month"], subset["mean_citations"], label=label, color=color, linewidth=2.5)
    ax.set_title("Citation accumulation over 36 months")
    ax.set_xlabel("Months since publication")
    ax.set_ylabel("Mean cumulative citations")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "preprint_citation_trajectory")
    plt.close(fig)

    server_metrics = (
        df[df["has_preprint"]]
        .groupby("preprint_server")
        .agg(
            review_duration_days=("review_duration_days", "mean"),
            citation_count_1yr=("citation_count_1yr", "mean"),
            engagement_index=("preprint_downloads", lambda x: np.mean(np.log1p(x)))
            
        )
        .loc[["bioRxiv", "arXiv", "medRxiv", "SSRN"]]
    )
    indexed = server_metrics.copy()
    for column in indexed.columns:
        indexed[column] = indexed[column] / indexed[column].mean() * 100

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(indexed.index))
    width = 0.24
    metric_colors = [OKABE_ITO["orange"], OKABE_ITO["blue"], OKABE_ITO["green"]]
    for i, column in enumerate(indexed.columns):
        ax.bar(x + (i - 1) * width, indexed[column], width=width, label=column.replace("_", " ").title(), color=metric_colors[i])
    ax.set_xticks(x)
    ax.set_xticklabels(indexed.index)
    ax.set_ylabel("Index (overall preprint server mean = 100)")
    ax.set_title("Server comparison across review, citation, and engagement metrics")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "preprint_server_comparison")
    plt.close(fig)

    scatter_df = df[df["has_preprint"]].copy()
    server_colors = {
        "bioRxiv": OKABE_ITO["green"],
        "arXiv": OKABE_ITO["blue"],
        "medRxiv": OKABE_ITO["vermillion"],
        "SSRN": OKABE_ITO["purple"],
    }
    fig, ax = plt.subplots(figsize=(8, 6))
    for server, group in scatter_df.groupby("preprint_server"):
        ax.scatter(
            group["preprint_comments"],
            group["peer_review_rounds"],
            alpha=0.55,
            s=28,
            label=server,
            color=server_colors[server],
        )
    slope, intercept, _, _, _ = stats.linregress(scatter_df["preprint_comments"], scatter_df["peer_review_rounds"])
    xline = np.linspace(0, scatter_df["preprint_comments"].max(), 100)
    ax.plot(xline, intercept + slope * xline, color="black", linestyle="--", linewidth=1.5)
    ax.set_title("Preprint feedback and peer review rounds")
    ax.set_xlabel("Preprint comments")
    ax.set_ylabel("Peer review rounds")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    save_figure(fig, "preprint_feedback_impact")
    plt.close(fig)


def write_preprocessing_log() -> None:
    PREPROCESSING_LOG_PATH.write_text(
        "# Preprocessing log\n\n"
        "- Random seeds set: numpy=42, random=42.\n"
        "- Generated 1,200 simulated papers with date, engagement, review, and citation fields.\n"
        "- Derived indicators: has_preprint, submission_to_publication_days, submission_to_acceptance_days, days_preprint_before_submission, citation_acceleration_ratio.\n"
        "- Missing date values are only permitted for preprint_date when preprint_server='none'.\n"
        "- Dates exported as ISO-8601 strings in the processed CSV.\n",
        encoding="utf-8",
    )



def write_statistical_summary(results: dict[str, object]) -> None:
    realized = results["realized_effects"]
    pub = results["submission_to_publication_analysis"]
    feedback = results["feedback_efficiency_analysis"]
    cites = results["early_citation_analysis"]
    engage = results["engagement_correlation_analysis"]
    server = results["server_specific_trends"]
    text = f"""# Statistical summary\n\n- Review-time reduction for preprints: {realized['review_time_reduction_pct']}% (mean {realized['mean_review_duration_preprint_days']} vs {realized['mean_review_duration_non_preprint_days']} days).\n- Submission-to-publication difference: {pub['mean_difference_days']} days, 95% CI {tuple(pub['mean_difference_95ci'])}, Welch p={pub['welch_p_value']:.4g}, Cohen's d={pub['cohens_d']}.\n- Feedback effect on review rounds: {feedback['mean_difference_rounds']} rounds, 95% CI {tuple(feedback['mean_difference_95ci'])}, Welch p={feedback['welch_p_value']:.4g}.\n- Early citation advantage: {realized['early_citation_advantage_pct']}% (mean difference {cites['mean_difference']} citations, 95% CI {tuple(cites['mean_difference_95ci'])}, Welch p={cites['welch_p_value']:.4g}, Cohen's d={cites['cohens_d']}).\n- Engagement correlations with 3-year citations: downloads rho={engage['spearman_downloads_vs_citations_3yr']['rho']} (p={engage['spearman_downloads_vs_citations_3yr']['p_value']:.4g}), comments rho={engage['spearman_comments_vs_citations_3yr']['rho']} (p={engage['spearman_comments_vs_citations_3yr']['p_value']:.4g}).\n- Server comparisons used Bonferroni correction for two ANOVA families: review p={server['anova_review_duration_days']['bonferroni_corrected_p_value']:.4g}, 1-year citations p={server['anova_citation_count_1yr']['bonferroni_corrected_p_value']:.4g}.\n"""
    STAT_SUMMARY_PATH.write_text(text, encoding="utf-8")



def write_report(results: dict[str, object]) -> None:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    realized = results["realized_effects"]
    feedback = results["feedback_efficiency_analysis"]
    versions = results["version_evolution_analysis"]
    engage = results["engagement_correlation_analysis"]
    server_rows = results["server_specific_trends"]["server_metrics"]
    top_server = max(server_rows, key=lambda row: row["mean_citation_count_1yr"])
    report = f"""# DRAFT — NOT FOR DISTRIBUTION\n\n## Preprint Role Evaluation Report\n\n- Timestamp: {timestamp}\n- Objective: Evaluate how simulated preprint posting relates to research communication speed, peer-review efficiency, and citation accumulation.\n\n## Methods\n\nA reproducible simulation generated {results['n_papers']} papers across bioRxiv, medRxiv, arXiv, SSRN, and non-preprint workflows using numpy/pandas with seed 42. The data model imposed two designed effects: roughly 15% shorter review duration for preprints and higher early citations for preprinted papers. Analyses used Welch t-tests, Spearman correlations, ANOVA for server-level contrasts, bootstrap 95% confidence intervals, and Bonferroni correction for the two ANOVA families. Figures use colorblind-friendly palettes and English-only labels.\n\n## Results\n\n- Preprints represented {results['preprint_share']:.1%} of papers in the simulated cohort.\n- Mean review duration was {realized['mean_review_duration_preprint_days']} days for preprints versus {realized['mean_review_duration_non_preprint_days']} days for non-preprints, a {realized['review_time_reduction_pct']}% reduction.\n- Preprints gained an early citation advantage of {realized['early_citation_advantage_pct']}% at 1 year.\n- Higher-feedback preprints had {feedback['high_feedback_mean_rounds']} review rounds on average versus {feedback['low_feedback_mean_rounds']} for lower-feedback preprints.\n- Version evolution showed comments and review efficiency shifting with additional preprint versions; see results JSON for bucketed summaries.\n- Among servers, {top_server['preprint_server']} showed the highest mean 1-year citation count ({top_server['mean_citation_count_1yr']}).\n- Engagement correlated positively with 3-year citations: downloads rho={engage['spearman_downloads_vs_citations_3yr']['rho']}, comments rho={engage['spearman_comments_vs_citations_3yr']['rho']}.\n\n## Discussion\n\nThe simulation supports a plausible narrative that preprints accelerate dissemination and modestly streamline peer review when public feedback accumulates before formal review. These results are synthetic rather than empirical, so effect sizes should be interpreted as scenario-based illustrations instead of real-world causal estimates.\n\n## Figure inventory\n\n- `figures/preprint_timeline.png` and `.svg`: box plot of review duration by preprint status.\n- `figures/preprint_citation_trajectory.png` and `.svg`: citation accumulation curves over 36 months.\n- `figures/preprint_server_comparison.png` and `.svg`: grouped server comparison chart.\n- `figures/preprint_feedback_impact.png` and `.svg`: comments versus review rounds scatter plot.\n\n## File inventory\n\n- `src/preprint_role_evaluation.py`\n- `data/preprint_simulated_dataset.csv`\n- `data/preprocessing-log.md`\n- `results/preprint_results.json`\n- `results/statistical-summary.md`\n- `report.md`\n- `logs/process-log.jsonl`\n"""
    REPORT_PATH.write_text(report, encoding="utf-8")



def write_process_log() -> None:
    timestamp = datetime.now(timezone.utc).isoformat()
    entries = [
        {"timestamp": timestamp, "phase": "PLAN", "event_type": "run_started", "actor": "co-scientist", "skill_or_tool": "co-scientist-data-analysis", "handoff_in": {"request": "simulate preprint server evaluation"}, "handoff_out": {}, "files_written": [], "status": "ok"},
        {"timestamp": timestamp, "phase": "PLAN", "event_type": "prompt_received", "actor": "co-scientist", "skill_or_tool": "user_prompt", "handoff_in": {"n_papers": 1200, "seed": 42}, "handoff_out": {}, "files_written": [], "status": "ok"},
        {"timestamp": timestamp, "phase": "PLAN", "event_type": "skill_selected", "actor": "co-scientist", "skill_or_tool": "co-scientist-data-analysis", "handoff_in": {}, "handoff_out": {"reason": "simulation, statistics, and visualization task"}, "files_written": [], "status": "ok"},
        {"timestamp": timestamp, "phase": "EXECUTE", "event_type": "handoff_started", "actor": "co-scientist", "skill_or_tool": "simulate_dataset", "handoff_in": {"fields": 16}, "handoff_out": {}, "files_written": [], "status": "ok"},
        {"timestamp": timestamp, "phase": "EXECUTE", "event_type": "handoff_completed", "actor": "co-scientist", "skill_or_tool": "analyze", "handoff_in": {}, "handoff_out": {"outputs": ["json", "figures", "csv", "report"]}, "files_written": [], "status": "ok"},
        {"timestamp": timestamp, "phase": "EXECUTE", "event_type": "file_written", "actor": "co-scientist", "skill_or_tool": "python", "handoff_in": {}, "handoff_out": {}, "files_written": ["data/preprint_simulated_dataset.csv", "data/preprocessing-log.md", "results/preprint_results.json", "results/statistical-summary.md", "figures/preprint_timeline.png", "figures/preprint_timeline.svg", "figures/preprint_citation_trajectory.png", "figures/preprint_citation_trajectory.svg", "figures/preprint_server_comparison.png", "figures/preprint_server_comparison.svg", "figures/preprint_feedback_impact.png", "figures/preprint_feedback_impact.svg", "report.md"], "status": "ok"},
        {"timestamp": timestamp, "phase": "REPORT", "event_type": "report_finalized", "actor": "co-scientist", "skill_or_tool": "write_report", "handoff_in": {}, "handoff_out": {"report": "report.md"}, "files_written": ["report.md"], "status": "ok"},
        {"timestamp": timestamp, "phase": "LOG", "event_type": "run_completed", "actor": "co-scientist", "skill_or_tool": "process_log", "handoff_in": {}, "handoff_out": {"status": "completed"}, "files_written": ["logs/process-log.jsonl"], "status": "ok"},
    ]
    with PROCESS_LOG_PATH.open("w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(json.dumps(entry) + "\n")



def main() -> None:
    ensure_dirs()
    df = simulate_dataset()
    trajectory_df = build_citation_trajectory(df)
    results = analyze(df)

    df_to_save = df.copy()
    for col in ["submission_date", "preprint_date", "acceptance_date", "publication_date"]:
        df_to_save[col] = df_to_save[col].dt.strftime("%Y-%m-%d")
    df_to_save.to_csv(DATA_DIR / "preprint_simulated_dataset.csv", index=False)
    trajectory_df.to_csv(RESULTS_DIR / "preprint_citation_trajectory.csv", index=False)

    save_figures(df, trajectory_df)

    with open(RESULTS_DIR / "preprint_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    write_preprocessing_log()
    write_statistical_summary(results)
    write_report(results)
    write_process_log()

    print(f"Saved simulated dataset to {DATA_DIR / 'preprint_simulated_dataset.csv'}")
    print(f"Saved figures to {FIGURES_DIR}")
    print(f"Saved results to {RESULTS_DIR / 'preprint_results.json'}")
    print(f"Saved report to {REPORT_PATH}")
    print(f"Saved process log to {PROCESS_LOG_PATH}")


if __name__ == "__main__":
    main()
