#!/usr/bin/env python3
"""Open Access Citation Advantage (OACA) causal estimation pipeline.

This script simulates bibliometric data with realistic confounding, estimates the
open access citation advantage with multiple causal inference methods, and saves
figures plus structured results.

Outputs
-------
- figures/oaca_propensity_distribution.png
- figures/oaca_causal_estimates.png
- figures/oaca_sensitivity_analysis.png
- results/oaca_results.json
- data/oaca_simulated_bibliometric_data.csv

Run with:
    python src/oaca_causal_estimation.py
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy import stats
from sklearn.linear_model import LinearRegression, LogisticRegression

SEED = 42
N_PAPERS = 2000
BOOTSTRAP_SAMPLES = 200
OA_MULTIPLIER = 1.3
CALIPER = 0.08

ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"

FIELD_CONFIG = {
    "Biomedicine": {"prob": 0.23, "citation_effect": 0.42, "altmetric_effect": 0.35, "author_rate": 7.0},
    "Physics": {"prob": 0.16, "citation_effect": 0.28, "altmetric_effect": 0.12, "author_rate": 5.4},
    "Computer Science": {"prob": 0.18, "citation_effect": 0.18, "altmetric_effect": 0.25, "author_rate": 4.1},
    "Economics": {"prob": 0.12, "citation_effect": 0.05, "altmetric_effect": 0.08, "author_rate": 2.8},
    "Psychology": {"prob": 0.14, "citation_effect": 0.12, "altmetric_effect": 0.16, "author_rate": 4.2},
    "Environmental Science": {"prob": 0.17, "citation_effect": 0.22, "altmetric_effect": 0.18, "author_rate": 5.1},
}


@dataclass
class MethodEstimate:
    name: str
    ate: float
    ci_lower: float
    ci_upper: float
    effect_size_d: float


def set_global_seed(seed: int = SEED) -> np.random.Generator:
    random.seed(seed)
    np.random.seed(seed)
    return np.random.default_rng(seed)


def ensure_directories() -> None:
    for directory in (FIGURES_DIR, RESULTS_DIR, DATA_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def build_journal_table(rng: np.random.Generator) -> pd.DataFrame:
    journals: List[Dict[str, object]] = []
    journal_id = 1
    switch_year = 2020
    for field, config in FIELD_CONFIG.items():
        for within_field in range(8):
            quality = rng.normal(0.0, 0.35)
            base_if = np.clip(rng.lognormal(mean=1.25 + config["citation_effect"] / 3, sigma=0.35), 1.0, 25.0)
            journals.append(
                {
                    "journal_id": f"J{journal_id:02d}",
                    "field": field,
                    "journal_name": f"{field.split()[0]} Journal {within_field + 1}",
                    "journal_base_if": float(base_if),
                    "journal_quality": float(quality),
                    "is_switch_journal": within_field < 2,
                    "switch_year": switch_year if within_field < 2 else np.nan,
                    "oa_policy_bias": 1.35 if within_field < 2 else rng.normal(0.0, 0.2),
                }
            )
            journal_id += 1
    return pd.DataFrame(journals)


def simulate_bibliometric_data(n_papers: int = N_PAPERS, seed: int = SEED) -> Tuple[pd.DataFrame, float]:
    rng = set_global_seed(seed)
    journals = build_journal_table(rng)

    fields = list(FIELD_CONFIG.keys())
    field_probs = np.array([FIELD_CONFIG[field]["prob"] for field in fields], dtype=float)
    field_probs /= field_probs.sum()

    sampled_fields = rng.choice(fields, size=n_papers, p=field_probs)
    years = rng.integers(2015, 2025, size=n_papers)

    records: List[Dict[str, object]] = []
    true_effects: List[float] = []

    for i in range(n_papers):
        field = sampled_fields[i]
        year = int(years[i])
        field_journals = journals.loc[journals["field"] == field].reset_index(drop=True)
        journal_row = field_journals.iloc[rng.integers(0, len(field_journals))]

        institution_rank = int(np.clip(np.round(rng.beta(2.8, 2.0) * 499) + 1, 1, 500))
        prestige = (501 - institution_rank) / 500
        author_h_index = int(
            np.clip(
                np.round(
                    rng.gamma(shape=4.2, scale=4.3)
                    + 0.55 * journal_row["journal_base_if"]
                    + 7.5 * prestige
                ),
                3,
                95,
            )
        )
        num_authors = int(
            np.clip(
                rng.poisson(lam=FIELD_CONFIG[field]["author_rate"] + 1.0 * prestige + 0.05 * journal_row["journal_base_if"]),
                1,
                20,
            )
        )
        journal_if = float(np.clip(journal_row["journal_base_if"] + rng.normal(0, 0.35), 0.5, 30.0))

        post_switch = bool(journal_row["is_switch_journal"] and year >= int(journal_row["switch_year"]))
        oa_logit = (
            -1.55
            + 0.18 * journal_if
            + 0.028 * author_h_index
            + 0.08 * num_authors
            + 0.85 * prestige
            + 0.10 * (year - 2015)
            + FIELD_CONFIG[field]["altmetric_effect"] * 0.45
            + float(journal_row["oa_policy_bias"])
            - 0.004 * institution_rank
        )
        oa_probability = 1.0 / (1.0 + np.exp(-oa_logit))
        oa_probability = float(np.clip(oa_probability, 0.02, 0.98))
        is_oa = 1 if post_switch else int(rng.binomial(1, oa_probability))

        age = 2025 - year
        log_mu_control = (
            1.55
            + FIELD_CONFIG[field]["citation_effect"]
            + 0.10 * age
            + 0.06 * journal_if
            + 0.018 * author_h_index
            + 0.055 * np.sqrt(num_authors)
            + 0.50 * prestige
            + 0.20 * float(journal_row["journal_quality"])
            - 0.0012 * institution_rank
        )
        mu_control = float(np.exp(log_mu_control))
        mu_treated = float(mu_control * OA_MULTIPLIER)
        mu_observed = mu_treated if is_oa else mu_control

        citation_lambda = rng.gamma(shape=4.0, scale=mu_observed / 4.0)
        citation_count = int(rng.poisson(lam=max(citation_lambda, 0.1)))

        altmetric_mean = np.exp(
            1.25
            + 0.11 * (year - 2015)
            + FIELD_CONFIG[field]["altmetric_effect"]
            + 0.05 * journal_if
            + 0.010 * author_h_index
            + 0.16 * is_oa
        )
        altmetric_score = int(np.round(rng.gamma(shape=3.0, scale=altmetric_mean / 3.0)))

        true_effects.append(mu_treated - mu_control)
        records.append(
            {
                "paper_id": f"P{i + 1:04d}",
                "year": year,
                "field": field,
                "journal_id": journal_row["journal_id"],
                "journal_name": journal_row["journal_name"],
                "journal_impact_factor": round(journal_if, 3),
                "author_h_index": author_h_index,
                "num_authors": num_authors,
                "institution_rank": institution_rank,
                "is_switch_journal": int(bool(journal_row["is_switch_journal"])),
                "switch_year": int(journal_row["switch_year"]) if not pd.isna(journal_row["switch_year"]) else None,
                "post_switch": int(post_switch),
                "oa_probability": round(oa_probability, 6),
                "is_oa": is_oa,
                "citation_count": citation_count,
                "altmetric_score": altmetric_score,
                "true_mu_control": round(mu_control, 6),
                "true_mu_treated": round(mu_treated, 6),
                "true_individual_effect": round(mu_treated - mu_control, 6),
            }
        )

    df = pd.DataFrame(records)
    return df, float(np.mean(true_effects))


COVARIATE_COLUMNS = [
    "year",
    "field",
    "journal_impact_factor",
    "author_h_index",
    "num_authors",
    "institution_rank",
]


def prepare_design_matrix(df: pd.DataFrame, columns: List[str], exclude_from_scaling: List[str] | None = None) -> pd.DataFrame:
    exclude_from_scaling = exclude_from_scaling or []
    design = pd.get_dummies(df[columns].copy(), columns=[col for col in columns if df[col].dtype == "object"], drop_first=True)
    for col in design.select_dtypes(include=[np.number]).columns:
        if col in exclude_from_scaling:
            continue
        std = design[col].std(ddof=0)
        if std > 0:
            design[col] = (design[col] - design[col].mean()) / std
        else:
            design[col] = 0.0
    return design.astype(float)


def estimate_propensity_scores(df: pd.DataFrame) -> np.ndarray:
    x = prepare_design_matrix(df, COVARIATE_COLUMNS)
    model = LogisticRegression(max_iter=2000, solver="lbfgs")
    model.fit(x, df["is_oa"].to_numpy())
    scores = model.predict_proba(x)[:, 1]
    return np.clip(scores, 0.01, 0.99)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights * values) / np.sum(weights))


def weighted_var(values: np.ndarray, weights: np.ndarray) -> float:
    mean_val = weighted_mean(values, weights)
    return float(np.sum(weights * (values - mean_val) ** 2) / np.sum(weights))


def cohens_d(values_t: np.ndarray, values_c: np.ndarray, weights_t: np.ndarray | None = None, weights_c: np.ndarray | None = None) -> float:
    if weights_t is None:
        weights_t = np.ones_like(values_t, dtype=float)
    if weights_c is None:
        weights_c = np.ones_like(values_c, dtype=float)
    mean_t = weighted_mean(values_t, weights_t)
    mean_c = weighted_mean(values_c, weights_c)
    var_t = weighted_var(values_t, weights_t)
    var_c = weighted_var(values_c, weights_c)
    pooled_sd = np.sqrt(max((var_t + var_c) / 2.0, 1e-9))
    return float((mean_t - mean_c) / pooled_sd)


def naive_estimator(df: pd.DataFrame) -> Dict[str, object]:
    treated = df.loc[df["is_oa"] == 1, "citation_count"].to_numpy(dtype=float)
    control = df.loc[df["is_oa"] == 0, "citation_count"].to_numpy(dtype=float)
    ate = float(treated.mean() - control.mean())
    return {"ate": ate, "effect_size_d": cohens_d(treated, control)}


def psm_estimator(df: pd.DataFrame, propensity_scores: np.ndarray | None = None) -> Dict[str, object]:
    if propensity_scores is None:
        propensity_scores = estimate_propensity_scores(df)

    treated_mask = df["is_oa"].to_numpy(dtype=int) == 1
    control_mask = ~treated_mask
    y = df["citation_count"].to_numpy(dtype=float)

    treated_idx = np.where(treated_mask)[0]
    control_idx = np.where(control_mask)[0]
    treated_scores = propensity_scores[treated_idx]
    control_scores = propensity_scores[control_idx]

    matched_treated_outcomes: List[float] = []
    matched_control_outcomes: List[float] = []
    pair_differences: List[float] = []
    unit_effects: List[float] = []

    for idx, score in zip(treated_idx, treated_scores):
        distances = np.abs(control_scores - score)
        best = np.argmin(distances)
        if distances[best] <= CALIPER:
            matched_control = control_idx[best]
            diff = y[idx] - y[matched_control]
            matched_treated_outcomes.append(y[idx])
            matched_control_outcomes.append(y[matched_control])
            pair_differences.append(diff)
            unit_effects.append(diff)

    for idx, score in zip(control_idx, control_scores):
        distances = np.abs(treated_scores - score)
        best = np.argmin(distances)
        if distances[best] <= CALIPER:
            matched_treated = treated_idx[best]
            diff = y[matched_treated] - y[idx]
            matched_treated_outcomes.append(y[matched_treated])
            matched_control_outcomes.append(y[idx])
            unit_effects.append(diff)

    if not unit_effects:
        raise RuntimeError("No propensity-score matches found within the caliper.")

    return {
        "ate": float(np.mean(unit_effects)),
        "effect_size_d": cohens_d(np.array(matched_treated_outcomes), np.array(matched_control_outcomes)),
        "pair_differences": pair_differences,
        "matched_fraction": float(len(unit_effects) / len(df)),
    }


def ipw_estimator(df: pd.DataFrame, propensity_scores: np.ndarray | None = None) -> Dict[str, object]:
    if propensity_scores is None:
        propensity_scores = estimate_propensity_scores(df)

    treatment = df["is_oa"].to_numpy(dtype=float)
    outcome = df["citation_count"].to_numpy(dtype=float)
    weights_t = treatment / propensity_scores
    weights_c = (1.0 - treatment) / (1.0 - propensity_scores)

    treated_values = outcome[treatment == 1]
    control_values = outcome[treatment == 0]
    treated_weights = weights_t[treatment == 1]
    control_weights = weights_c[treatment == 0]

    ate = weighted_mean(treated_values, treated_weights) - weighted_mean(control_values, control_weights)
    effect_size = cohens_d(treated_values, control_values, treated_weights, control_weights)
    return {"ate": float(ate), "effect_size_d": effect_size}


def fit_outcome_model(df: pd.DataFrame) -> Tuple[LinearRegression, pd.DataFrame]:
    x = prepare_design_matrix(df, COVARIATE_COLUMNS + ["is_oa"])
    model = LinearRegression()
    model.fit(x, df["citation_count"].to_numpy(dtype=float))
    return model, x


def doubly_robust_estimator(df: pd.DataFrame, propensity_scores: np.ndarray | None = None) -> Dict[str, object]:
    if propensity_scores is None:
        propensity_scores = estimate_propensity_scores(df)

    x_base = prepare_design_matrix(df, COVARIATE_COLUMNS)
    x_treated = x_base.copy()
    x_control = x_base.copy()
    x_treated["is_oa"] = 1.0
    x_control["is_oa"] = 0.0

    x_full = x_base.copy()
    x_full["is_oa"] = df["is_oa"].to_numpy(dtype=float)

    outcome_model = LinearRegression()
    outcome_model.fit(x_full, df["citation_count"].to_numpy(dtype=float))

    m1 = outcome_model.predict(x_treated)
    m0 = outcome_model.predict(x_control)
    treatment = df["is_oa"].to_numpy(dtype=float)
    outcome = df["citation_count"].to_numpy(dtype=float)

    pseudo_outcome = m1 - m0 + treatment * (outcome - m1) / propensity_scores - (1.0 - treatment) * (outcome - m0) / (1.0 - propensity_scores)
    ate = float(np.mean(pseudo_outcome))
    effect_size = float(ate / np.sqrt(max((np.var(m1, ddof=0) + np.var(m0, ddof=0)) / 2.0, 1e-9)))
    return {"ate": ate, "effect_size_d": effect_size}


def did_estimator(df: pd.DataFrame) -> Dict[str, object]:
    did_df = df.copy()
    did_df["treated_group"] = did_df["is_switch_journal"].astype(int)
    did_df["post_period"] = (did_df["year"] >= 2020).astype(int)
    did_df["treat_post"] = did_df["treated_group"] * did_df["post_period"]

    covariates = [
        "treated_group",
        "post_period",
        "treat_post",
        "journal_impact_factor",
        "author_h_index",
        "num_authors",
        "institution_rank",
        "field",
    ]
    x = prepare_design_matrix(did_df, covariates, exclude_from_scaling=["treated_group", "post_period", "treat_post"])
    model = LinearRegression()
    model.fit(x, did_df["citation_count"].to_numpy(dtype=float))
    coefficient = float(model.coef_[list(x.columns).index("treat_post")])

    post_treated = did_df.loc[(did_df["treated_group"] == 1) & (did_df["post_period"] == 1), "citation_count"].to_numpy(dtype=float)
    post_control = did_df.loc[(did_df["treated_group"] == 0) & (did_df["post_period"] == 1), "citation_count"].to_numpy(dtype=float)
    effect_size = float(coefficient / np.sqrt(max((np.var(post_treated, ddof=0) + np.var(post_control, ddof=0)) / 2.0, 1e-9)))
    return {"ate": coefficient, "effect_size_d": effect_size}


def bootstrap_confidence_interval(
    df: pd.DataFrame,
    estimator: Callable[[pd.DataFrame], Dict[str, object]],
    seed: int,
    n_bootstrap: int = BOOTSTRAP_SAMPLES,
) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates: List[float] = []
    n = len(df)

    for _ in range(n_bootstrap):
        sample_idx = rng.integers(0, n, size=n)
        sample = df.iloc[sample_idx].reset_index(drop=True)
        try:
            estimates.append(float(estimator(sample)["ate"]))
        except Exception:
            continue

    if len(estimates) < max(50, n_bootstrap // 2):
        raise RuntimeError("Bootstrap produced too few successful replicates.")

    lower, upper = np.percentile(estimates, [2.5, 97.5])
    return float(lower), float(upper)


def rosenbaum_sensitivity(pair_differences: np.ndarray, gamma_values: np.ndarray) -> pd.DataFrame:
    non_zero_diffs = pair_differences[np.abs(pair_differences) > 1e-9]
    positive = int(np.sum(non_zero_diffs > 0))
    total = int(len(non_zero_diffs))
    if total == 0:
        raise RuntimeError("Rosenbaum sensitivity analysis requires non-zero matched pair differences.")

    rows = []
    for gamma in gamma_values:
        p_low = 1.0 / (1.0 + gamma)
        p_high = gamma / (1.0 + gamma)
        lower_p_value = float(stats.binom.sf(positive - 1, total, p_low))
        upper_p_value = float(stats.binom.sf(positive - 1, total, p_high))
        rows.append({"gamma": float(gamma), "lower_p_value": lower_p_value, "upper_p_value": upper_p_value})
    return pd.DataFrame(rows)


def make_propensity_plot(df: pd.DataFrame, propensity_scores: np.ndarray) -> Path:
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, 2))
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.linspace(0, 1, 25)
    ax.hist(propensity_scores[df["is_oa"] == 1], bins=bins, alpha=0.70, color=colors[1], label="OA", density=True)
    ax.hist(propensity_scores[df["is_oa"] == 0], bins=bins, alpha=0.60, color=colors[0], label="Non-OA", density=True)
    ax.set_xlabel("Estimated propensity score")
    ax.set_ylabel("Density")
    ax.set_title("Propensity score distributions")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    output_path = FIGURES_DIR / "oaca_propensity_distribution.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def make_forest_plot(estimates: List[MethodEstimate], true_ate: float) -> Path:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(estimates)))
    y_positions = np.arange(len(estimates))[::-1]

    for y, estimate, color in zip(y_positions, estimates, colors):
        ax.plot([estimate.ci_lower, estimate.ci_upper], [y, y], color=color, lw=2.5)
        ax.scatter(estimate.ate, y, s=70, color=color, zorder=3)

    ax.axvline(true_ate, color="black", linestyle="--", linewidth=1.5, label=f"True ATE = {true_ate:.2f}")
    ax.set_yticks(y_positions)
    ax.set_yticklabels([estimate.name for estimate in estimates])
    ax.set_xlabel("Average treatment effect on citations")
    ax.set_title("OACA estimates across causal methods")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, loc="lower right")

    output_path = FIGURES_DIR / "oaca_causal_estimates.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def make_sensitivity_plot(sensitivity_df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.viridis(np.linspace(0.25, 0.75, 2))
    ax.plot(sensitivity_df["gamma"], sensitivity_df["lower_p_value"], color=colors[0], linewidth=2.2, label="Lower bound p-value")
    ax.plot(sensitivity_df["gamma"], sensitivity_df["upper_p_value"], color=colors[1], linewidth=2.2, label="Upper bound p-value")
    ax.axhline(0.05, color="black", linestyle="--", linewidth=1.2, label="p = 0.05")
    ax.set_xlabel("Rosenbaum Γ")
    ax.set_ylabel("One-sided p-value bound")
    ax.set_title("Rosenbaum sensitivity analysis")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    output_path = FIGURES_DIR / "oaca_sensitivity_analysis.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def summarize_data(df: pd.DataFrame) -> Dict[str, object]:
    oa_mask = df["is_oa"] == 1
    return {
        "n_papers": int(len(df)),
        "oa_share": float(df["is_oa"].mean()),
        "mean_citations_oa": float(df.loc[oa_mask, "citation_count"].mean()),
        "mean_citations_non_oa": float(df.loc[~oa_mask, "citation_count"].mean()),
        "mean_jif_oa": float(df.loc[oa_mask, "journal_impact_factor"].mean()),
        "mean_jif_non_oa": float(df.loc[~oa_mask, "journal_impact_factor"].mean()),
        "mean_hindex_oa": float(df.loc[oa_mask, "author_h_index"].mean()),
        "mean_hindex_non_oa": float(df.loc[~oa_mask, "author_h_index"].mean()),
        "switch_journal_share": float(df["is_switch_journal"].mean()),
        "years": {str(year): int(count) for year, count in df["year"].value_counts().sort_index().items()},
        "fields": {field: int(count) for field, count in df["field"].value_counts().sort_index().items()},
    }


def main() -> None:
    ensure_directories()
    df, true_ate = simulate_bibliometric_data()
    propensity_scores = estimate_propensity_scores(df)

    naive = naive_estimator(df)
    psm = psm_estimator(df, propensity_scores)
    ipw = ipw_estimator(df, propensity_scores)
    dr = doubly_robust_estimator(df, propensity_scores)
    did = did_estimator(df)

    estimates = {
        "Naive comparison": naive,
        "Propensity score matching": psm,
        "Inverse probability weighting": ipw,
        "Doubly robust": dr,
        "Difference-in-differences": did,
    }

    ci_builders = {
        "Naive comparison": lambda sample: naive_estimator(sample),
        "Propensity score matching": lambda sample: psm_estimator(sample),
        "Inverse probability weighting": lambda sample: ipw_estimator(sample),
        "Doubly robust": lambda sample: doubly_robust_estimator(sample),
        "Difference-in-differences": lambda sample: did_estimator(sample),
    }

    ordered_method_names = list(estimates.keys())
    method_estimates: List[MethodEstimate] = []
    for idx, name in enumerate(ordered_method_names):
        lower, upper = bootstrap_confidence_interval(df, ci_builders[name], seed=SEED + 100 + idx)
        result = estimates[name]
        method_estimates.append(
            MethodEstimate(
                name=name,
                ate=float(result["ate"]),
                ci_lower=lower,
                ci_upper=upper,
                effect_size_d=float(result["effect_size_d"]),
            )
        )

    gamma_values = np.linspace(1.0, 3.0, 11)
    sensitivity_df = rosenbaum_sensitivity(np.array(psm["pair_differences"], dtype=float), gamma_values)

    propensity_plot = make_propensity_plot(df, propensity_scores)
    forest_plot = make_forest_plot(method_estimates, true_ate)
    sensitivity_plot = make_sensitivity_plot(sensitivity_df)

    data_path = DATA_DIR / "oaca_simulated_bibliometric_data.csv"
    df.to_csv(data_path, index=False)

    results = {
        "metadata": {
            "seed": SEED,
            "n_papers": N_PAPERS,
            "bootstrap_samples": BOOTSTRAP_SAMPLES,
            "oa_multiplier": OA_MULTIPLIER,
            "caliper": CALIPER,
            "output_files": {
                "data": str(data_path),
                "propensity_plot": str(propensity_plot),
                "forest_plot": str(forest_plot),
                "sensitivity_plot": str(sensitivity_plot),
            },
            "package_versions": {
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
                "scikit_learn": sklearn.__version__,
                "matplotlib": matplotlib.__version__,
            },
        },
        "data_summary": summarize_data(df),
        "ground_truth": {
            "embedded_oa_multiplier": OA_MULTIPLIER,
            "true_average_treatment_effect": true_ate,
        },
        "methods": {
            estimate.name: {
                "ate": estimate.ate,
                "ci_95": [estimate.ci_lower, estimate.ci_upper],
                "effect_size_d": estimate.effect_size_d,
            }
            for estimate in method_estimates
        },
        "matching_diagnostics": {
            "psm_matched_fraction": float(psm.get("matched_fraction", np.nan)),
            "mean_propensity_score_oa": float(propensity_scores[df["is_oa"] == 1].mean()),
            "mean_propensity_score_non_oa": float(propensity_scores[df["is_oa"] == 0].mean()),
        },
        "sensitivity_analysis": {
            "method": "Rosenbaum sign-test bounds on matched pairs",
            "gamma_values": sensitivity_df["gamma"].round(3).tolist(),
            "lower_p_values": sensitivity_df["lower_p_value"].round(6).tolist(),
            "upper_p_values": sensitivity_df["upper_p_value"].round(6).tolist(),
        },
    }

    results_path = RESULTS_DIR / "oaca_results.json"
    with results_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)

    print("OACA causal estimation pipeline complete.")
    print(f"Saved data: {data_path}")
    print(f"Saved results: {results_path}")
    print(f"Saved figures: {propensity_plot.name}, {forest_plot.name}, {sensitivity_plot.name}")


if __name__ == "__main__":
    main()
