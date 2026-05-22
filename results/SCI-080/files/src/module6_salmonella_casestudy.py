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
import statsmodels.api as sm
from scipy.special import expit
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    auc,
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from statsmodels.stats.multitest import multipletests

SEED = 42
DPI = 300
TARGET_RATE = 0.15


plt.style.use("seaborn-v0_8-whitegrid")
np.random.seed(SEED)
random.seed(SEED)


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
LOG_DIR = BASE_DIR / "logs"
REPORT_PATH = BASE_DIR / "report.md"
PREPROCESS_LOG_PATH = DATA_DIR / "preprocessing-log.md"
PROCESS_LOG_PATH = LOG_DIR / "process-log.jsonl"
DATA_PATH = DATA_DIR / "salmonella_data.csv"
METRICS_PATH = RESULTS_DIR / "module6_metrics.json"
SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
OR_TABLE_PATH = RESULTS_DIR / "module6_odds_ratios.csv"
PERM_PATH = RESULTS_DIR / "module6_permutation_importance.csv"
INTERVENTION_PATH = RESULTS_DIR / "module6_intervention_summary.csv"


ORIGINAL_FEATURES = [
    "farm_id",
    "flock_age",
    "ambient_temp",
    "humidity",
    "season",
    "processing_plant",
    "chilling_method",
    "processing_speed",
    "pre_chill_count",
    "post_chill_count",
    "final_product_count",
    "storage_temp",
    "days_to_retail",
]

MODEL_FEATURES = [
    "farm_id",
    "flock_age",
    "ambient_temp",
    "humidity",
    "season",
    "processing_plant",
    "chilling_method",
    "processing_speed",
    "storage_temp",
    "days_to_retail",
    "log_pre_chill_count",
    "log_post_chill_count",
    "log_final_product_count",
    "temp_x_summer",
    "water_x_speed",
]

NUMERIC_FEATURES = [
    "flock_age",
    "ambient_temp",
    "humidity",
    "processing_speed",
    "storage_temp",
    "days_to_retail",
    "log_pre_chill_count",
    "log_post_chill_count",
    "log_final_product_count",
    "temp_x_summer",
    "water_x_speed",
]

CATEGORICAL_FEATURES = ["farm_id", "season", "processing_plant", "chilling_method"]

MODEL_DISPLAY = {
    "Logistic Regression": "Logistic Regression",
    "Random Forest": "Random Forest",
    "Gradient Boosting": "Gradient Boosting",
}

FEATURE_LABELS = {
    "ambient_temp": "Ambient temperature",
    "humidity": "Humidity",
    "flock_age": "Flock age",
    "processing_speed": "Processing speed",
    "storage_temp": "Storage temperature",
    "days_to_retail": "Days to retail",
    "log_pre_chill_count": "Pre-chill count (log10)",
    "log_post_chill_count": "Post-chill count (log10)",
    "log_final_product_count": "Final product count (log10)",
    "summer_indicator": "Summer",
    "water_indicator": "Water chilling",
    "temp_x_summer": "Temperature × Summer",
    "water_x_speed": "Water chilling × Speed",
    "farm_id": "Farm ID",
    "season": "Season",
    "processing_plant": "Processing plant",
    "chilling_method": "Chilling method",
}


def ensure_directories() -> None:
    for path in [DATA_DIR, FIG_DIR, RESULTS_DIR, LOG_DIR, BASE_DIR / "src"]:
        path.mkdir(parents=True, exist_ok=True)


def append_log(event_type: str, phase: str, handoff_in: dict | None = None, handoff_out: dict | None = None,
               files_written: list[str] | None = None, skill_or_tool: str = "co-scientist-data-analysis",
               status: str = "ok") -> None:
    payload = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files_written or [],
        "status": status,
    }
    with PROCESS_LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload) + "\n")


def month_to_season(month: int) -> str:
    if month in (12, 1, 2):
        return "Winter"
    if month in (3, 4, 5):
        return "Spring"
    if month in (6, 7, 8):
        return "Summer"
    return "Autumn"


def calibrate_intercept(raw_score: np.ndarray, target_rate: float) -> float:
    low, high = -15.0, 15.0
    for _ in range(80):
        mid = (low + high) / 2
        rate = expit(raw_score + mid).mean()
        if rate > target_rate:
            high = mid
        else:
            low = mid
    return (low + high) / 2


def generate_synthetic_data(n_samples: int = 2000, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    start_date = pd.Timestamp("2023-01-01")
    end_date = pd.Timestamp("2024-12-31")
    total_days = (end_date - start_date).days + 1
    sample_offsets = rng.integers(0, total_days, size=n_samples)
    sample_dates = start_date + pd.to_timedelta(sample_offsets, unit="D")

    months = pd.Series(sample_dates).dt.month.to_numpy()
    day_of_year = pd.Series(sample_dates).dt.dayofyear.to_numpy()
    seasons = np.array([month_to_season(month) for month in months])

    farms = [f"Farm_{i:02d}" for i in range(1, 21)]
    plants = [f"Plant_{c}" for c in list("ABCDE")]
    farm_id = rng.choice(farms, size=n_samples)
    processing_plant = rng.choice(plants, size=n_samples, p=[0.24, 0.21, 0.2, 0.18, 0.17])
    chilling_method = rng.choice(["air", "water"], size=n_samples, p=[0.42, 0.58])

    seasonal_temp = 18 + 11 * np.sin(2 * np.pi * (day_of_year - 172) / 365.25)
    ambient_temp = np.clip(seasonal_temp + rng.normal(0, 3.0, size=n_samples), 0, 38)
    humidity = np.clip(62 + 0.35 * (ambient_temp - 18) + rng.normal(0, 8, size=n_samples), 35, 95)
    flock_age = np.clip(rng.normal(42, 6.5, size=n_samples), 28, 63)
    storage_temp = np.clip(rng.normal(3.8, 1.1, size=n_samples), 0.5, 8.5)
    days_to_retail = np.clip(np.round(rng.normal(4.2, 1.7, size=n_samples)), 1, 9).astype(int)

    speed_base = np.where(chilling_method == "water", 119, 108)
    processing_speed = np.clip(speed_base + rng.normal(0, 8.5, size=n_samples), 82, 145)

    farm_effects = dict(zip(farms, rng.normal(0, 0.22, size=len(farms))))
    plant_effects = dict(zip(plants, rng.normal(0, 0.18, size=len(plants))))

    season_effect_pre = {"Winter": -0.15, "Spring": 0.0, "Summer": 0.28, "Autumn": 0.12}
    season_effect_risk = {"Winter": -0.28, "Spring": 0.05, "Summer": 0.72, "Autumn": 0.15}

    farm_risk = np.array([farm_effects[f] for f in farm_id])
    plant_risk = np.array([plant_effects[p] for p in processing_plant])
    season_pre = np.array([season_effect_pre[s] for s in seasons])
    season_risk = np.array([season_effect_risk[s] for s in seasons])

    log_pre = (
        3.9
        + 0.035 * (ambient_temp - 20)
        + 0.012 * (humidity - 60)
        + 0.015 * (flock_age - 42)
        + season_pre
        + 0.55 * farm_risk
        + rng.normal(0, 0.22, size=n_samples)
    )

    chill_reduction = (
        1.1
        + 0.24 * (chilling_method == "air").astype(float)
        - 0.0055 * (processing_speed - 110)
        - 0.15 * plant_risk
        + rng.normal(0, 0.12, size=n_samples)
    )

    log_post = log_pre - chill_reduction + rng.normal(0, 0.12, size=n_samples)
    storage_growth = (
        0.055 * (storage_temp - 3.5)
        + 0.038 * (days_to_retail - 4)
        + 0.06 * (seasons == "Summer").astype(float)
        + rng.normal(0, 0.08, size=n_samples)
    )
    log_final = log_post + storage_growth

    pre_count = np.round(np.power(10, np.clip(log_pre, 1.6, 6.2))).astype(int)
    post_count = np.round(np.power(10, np.clip(log_post, 0.8, 5.4))).astype(int)
    final_count = np.round(np.power(10, np.clip(log_final, 0.5, 5.8))).astype(int)

    raw_score = (
        0.07 * (ambient_temp - 20)
        + 0.018 * (humidity - 60)
        + 0.022 * (processing_speed - 110)
        + 0.16 * (storage_temp - 3.5)
        + 0.055 * (days_to_retail - 4)
        + 0.18 * (flock_age - 42) / 10
        + 0.52 * (np.log10(pre_count + 1) - 3.9)
        + 0.78 * (np.log10(post_count + 1) - 2.8)
        + 0.95 * (np.log10(final_count + 1) - 2.4)
        + 0.34 * (chilling_method == "water").astype(float)
        + season_risk
        + 0.95 * farm_risk
        + 0.55 * plant_risk
        + 0.028 * (ambient_temp - 24) * (seasons == "Summer").astype(float)
        + 0.018 * (processing_speed - 110) * (chilling_method == "water").astype(float)
    )

    intercept = calibrate_intercept(raw_score, TARGET_RATE)
    contamination_prob = expit(raw_score + intercept)
    salmonella_positive = rng.binomial(1, contamination_prob)

    serotypes = np.array(["None"] * n_samples, dtype=object)
    positive_idx = np.where(salmonella_positive == 1)[0]
    if len(positive_idx) > 0:
        summer_positive = (seasons[positive_idx] == "Summer").astype(float)
        serotype_choices = []
        for is_summer in summer_positive:
            if is_summer:
                probs = [0.34, 0.28, 0.15, 0.13, 0.10]
            else:
                probs = [0.28, 0.31, 0.16, 0.11, 0.14]
            serotype_choices.append(
                rng.choice(
                    ["Enteritidis", "Typhimurium", "Heidelberg", "Infantis", "Kentucky"],
                    p=probs,
                )
            )
        serotypes[positive_idx] = serotype_choices

    df = pd.DataFrame(
        {
            "sample_date": pd.to_datetime(sample_dates),
            "farm_id": farm_id,
            "flock_age": np.round(flock_age, 1),
            "ambient_temp": np.round(ambient_temp, 1),
            "humidity": np.round(humidity, 1),
            "season": seasons,
            "processing_plant": processing_plant,
            "chilling_method": chilling_method,
            "processing_speed": np.round(processing_speed, 1),
            "pre_chill_count": pre_count,
            "post_chill_count": post_count,
            "final_product_count": final_count,
            "salmonella_positive": salmonella_positive.astype(int),
            "salmonella_serotype": serotypes,
            "storage_temp": np.round(storage_temp, 1),
            "days_to_retail": days_to_retail.astype(int),
        }
    ).sort_values("sample_date").reset_index(drop=True)
    return df


def prepare_model_data(df: pd.DataFrame) -> pd.DataFrame:
    modeled = df.copy()
    modeled["log_pre_chill_count"] = np.log10(modeled["pre_chill_count"] + 1)
    modeled["log_post_chill_count"] = np.log10(modeled["post_chill_count"] + 1)
    modeled["log_final_product_count"] = np.log10(modeled["final_product_count"] + 1)
    modeled["summer_indicator"] = (modeled["season"] == "Summer").astype(int)
    modeled["water_indicator"] = (modeled["chilling_method"] == "water").astype(int)
    modeled["temp_x_summer"] = modeled["ambient_temp"] * modeled["summer_indicator"]
    modeled["water_x_speed"] = modeled["processing_speed"] * modeled["water_indicator"]
    return modeled


def build_preprocessor() -> ColumnTransformer:
    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    return ColumnTransformer(
        transformers=[
            ("num", StandardScaler(), NUMERIC_FEATURES),
            ("cat", encoder, CATEGORICAL_FEATURES),
        ]
    )


def build_models() -> dict[str, Pipeline]:
    preprocessor = build_preprocessor()
    return {
        "Logistic Regression": Pipeline(
            steps=[
                ("preprocessor", preprocessor),
                ("model", LogisticRegression(max_iter=3000, random_state=SEED)),
            ]
        ),
        "Random Forest": Pipeline(
            steps=[
                ("preprocessor", preprocessor),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=450,
                        min_samples_leaf=2,
                        class_weight="balanced_subsample",
                        random_state=SEED,
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
        "Gradient Boosting": Pipeline(
            steps=[
                ("preprocessor", preprocessor),
                (
                    "model",
                    GradientBoostingClassifier(
                        n_estimators=260,
                        learning_rate=0.05,
                        max_depth=3,
                        subsample=0.9,
                        random_state=SEED,
                    ),
                ),
            ]
        ),
    }


def classification_metrics(y_true: pd.Series, probabilities: np.ndarray, threshold: float = 0.5) -> dict[str, float]:
    predictions = (probabilities >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, predictions).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    ppv = tp / (tp + fp) if (tp + fp) else 0.0
    return {
        "auc": float(roc_auc_score(y_true, probabilities)),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "ppv": float(ppv),
    }


def evaluate_models(X: pd.DataFrame, y: pd.Series, models: dict[str, Pipeline]) -> tuple[dict, dict]:
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    evaluations: dict[str, dict] = {}
    curve_data: dict[str, dict] = {}

    for name, model in models.items():
        probabilities = cross_val_predict(model, X, y, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
        auc_scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)
        metrics = classification_metrics(y, probabilities)
        metrics["cv_auc_mean"] = float(np.mean(auc_scores))
        metrics["cv_auc_std"] = float(np.std(auc_scores, ddof=1))

        fpr, tpr, _ = roc_curve(y, probabilities)
        precision, recall, _ = precision_recall_curve(y, probabilities)
        curve_data[name] = {
            "fpr": fpr,
            "tpr": tpr,
            "precision": precision,
            "recall": recall,
            "pr_auc": average_precision_score(y, probabilities),
            "probabilities": probabilities,
        }
        evaluations[name] = metrics

    return evaluations, curve_data


def compute_odds_ratios(df_model: pd.DataFrame) -> pd.DataFrame:
    y = df_model["salmonella_positive"]
    X = df_model[
        [
            "flock_age",
            "ambient_temp",
            "humidity",
            "processing_speed",
            "storage_temp",
            "days_to_retail",
            "log_pre_chill_count",
            "log_post_chill_count",
            "log_final_product_count",
            "summer_indicator",
            "water_indicator",
            "temp_x_summer",
            "water_x_speed",
        ]
    ].copy()
    X = sm.add_constant(X)
    model = sm.GLM(y, X, family=sm.families.Binomial()).fit()
    conf = model.conf_int()

    results = pd.DataFrame(
        {
            "feature": model.params.index,
            "coef": model.params.values,
            "odds_ratio": np.exp(model.params.values),
            "ci_lower": np.exp(conf[0].values),
            "ci_upper": np.exp(conf[1].values),
            "p_value": model.pvalues.values,
        }
    )
    mask = results["feature"] != "const"
    q_values = multipletests(results.loc[mask, "p_value"], method="fdr_bh")[1]
    results.loc[mask, "q_value"] = q_values
    results.loc[~mask, "q_value"] = np.nan
    return results.sort_values("odds_ratio", ascending=False).reset_index(drop=True)


def compute_permutation_importance(best_model: Pipeline, X_test: pd.DataFrame, y_test: pd.Series) -> pd.DataFrame:
    perm = permutation_importance(
        best_model,
        X_test,
        y_test,
        scoring="roc_auc",
        n_repeats=20,
        random_state=SEED,
        n_jobs=-1,
    )
    df_perm = pd.DataFrame(
        {
            "feature": X_test.columns,
            "importance_mean": perm.importances_mean,
            "importance_std": perm.importances_std,
        }
    ).sort_values("importance_mean", ascending=False)
    return df_perm.reset_index(drop=True)


def simulate_interventions(df: pd.DataFrame, fitted_model: Pipeline) -> pd.DataFrame:
    baseline_model_df = prepare_model_data(df)
    baseline_rate = fitted_model.predict_proba(baseline_model_df[MODEL_FEATURES])[:, 1].mean()
    interventions = {
        "Improved chilling": {"cost_per_lot": 70},
        "Reduced processing speed": {"cost_per_lot": 32},
        "Enhanced sanitation": {"cost_per_lot": 48},
    }
    cost_per_case = 18000
    records = []

    for name, config in interventions.items():
        scenario = df.copy()
        if name == "Improved chilling":
            scenario.loc[scenario["chilling_method"] == "water", "chilling_method"] = "air"
            scenario["post_chill_count"] = np.maximum((scenario["post_chill_count"] * 0.72).round().astype(int), 1)
            scenario["final_product_count"] = np.maximum((scenario["final_product_count"] * 0.78).round().astype(int), 1)
        elif name == "Reduced processing speed":
            scenario["processing_speed"] = np.clip(scenario["processing_speed"] * 0.88, 80, None).round(1)
            scenario["post_chill_count"] = np.maximum((scenario["post_chill_count"] * 0.86).round().astype(int), 1)
            scenario["final_product_count"] = np.maximum((scenario["final_product_count"] * 0.9).round().astype(int), 1)
        elif name == "Enhanced sanitation":
            scenario["pre_chill_count"] = np.maximum((scenario["pre_chill_count"] * 0.7).round().astype(int), 1)
            scenario["post_chill_count"] = np.maximum((scenario["post_chill_count"] * 0.68).round().astype(int), 1)
            scenario["final_product_count"] = np.maximum((scenario["final_product_count"] * 0.74).round().astype(int), 1)

        scenario_model_df = prepare_model_data(scenario)
        after_rate = fitted_model.predict_proba(scenario_model_df[MODEL_FEATURES])[:, 1].mean()
        avoided_cases = max((baseline_rate - after_rate) * len(df), 0)
        intervention_cost = config["cost_per_lot"] * len(df)
        savings = avoided_cases * cost_per_case
        records.append(
            {
                "intervention": name,
                "baseline_rate": baseline_rate,
                "after_rate": after_rate,
                "absolute_reduction": baseline_rate - after_rate,
                "relative_reduction_pct": 100 * (baseline_rate - after_rate) / baseline_rate,
                "avoided_cases": avoided_cases,
                "implementation_cost": intervention_cost,
                "expected_savings": savings,
                "net_benefit": savings - intervention_cost,
            }
        )

    return pd.DataFrame(records).sort_values("absolute_reduction", ascending=False).reset_index(drop=True)


def plot_roc_and_pr(curves: dict[str, dict]) -> None:
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(curves)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    for color, (name, payload) in zip(colors, curves.items()):
        axes[0].plot(payload["fpr"], payload["tpr"], color=color, linewidth=2,
                     label=f"{name} (AUC={roc_auc_score_y(payload):.3f})")
        axes[1].plot(payload["recall"], payload["precision"], color=color, linewidth=2,
                     label=f"{name} (AP={payload['pr_auc']:.3f})")

    axes[0].plot([0, 1], [0, 1], linestyle="--", color="0.5", linewidth=1)
    axes[0].set_xlabel("False positive rate")
    axes[0].set_ylabel("True positive rate")
    axes[0].set_title("ROC curves")
    axes[0].legend(frameon=True)

    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title("Precision-recall curves")
    axes[1].legend(frameon=True)

    fig.savefig(FIG_DIR / "fig6_roc_curves.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def roc_auc_score_y(payload: dict[str, np.ndarray]) -> float:
    return auc(payload["fpr"], payload["tpr"])


def plot_risk_factors(or_df: pd.DataFrame, perm_df: pd.DataFrame) -> None:
    colors = plt.cm.viridis(np.linspace(0.2, 0.85, 10))
    forest_df = or_df[or_df["feature"] != "const"].copy()
    forest_df["distance"] = np.abs(np.log(forest_df["odds_ratio"]))
    forest_df = forest_df.sort_values(["q_value", "distance"], ascending=[True, False]).head(10)
    forest_df = forest_df.sort_values("odds_ratio")

    top_perm = perm_df.head(10).sort_values("importance_mean")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    y_pos = np.arange(len(forest_df))
    axes[0].errorbar(
        forest_df["odds_ratio"],
        y_pos,
        xerr=[forest_df["odds_ratio"] - forest_df["ci_lower"], forest_df["ci_upper"] - forest_df["odds_ratio"]],
        fmt="o",
        color=colors[-1],
        ecolor="0.35",
        capsize=3,
    )
    axes[0].axvline(1.0, linestyle="--", color="0.5", linewidth=1)
    axes[0].set_xscale("log")
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels([FEATURE_LABELS.get(f, f) for f in forest_df["feature"]])
    axes[0].set_xlabel("Odds ratio (log scale)")
    axes[0].set_title("Logistic regression odds ratios")

    axes[1].barh(
        [FEATURE_LABELS.get(f, f) for f in top_perm["feature"]],
        top_perm["importance_mean"],
        xerr=top_perm["importance_std"],
        color=colors,
        edgecolor="none",
    )
    axes[1].set_xlabel("Permutation importance (AUC decrease)")
    axes[1].set_title("SHAP-like ranked feature importance")

    fig.savefig(FIG_DIR / "fig6b_risk_factors.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_seasonal_contamination(df: pd.DataFrame) -> pd.DataFrame:
    monthly = (
        df.assign(month=df["sample_date"].dt.to_period("M").dt.to_timestamp())
        .groupby("month")
        .agg(contamination_rate=("salmonella_positive", "mean"), ambient_temp=("ambient_temp", "mean"))
        .reset_index()
    )

    colors = plt.cm.viridis([0.2, 0.8])
    fig, ax1 = plt.subplots(figsize=(12, 5), constrained_layout=True)
    ax2 = ax1.twinx()

    ax1.plot(monthly["month"], monthly["contamination_rate"] * 100, color=colors[1], marker="o", linewidth=2)
    ax2.plot(monthly["month"], monthly["ambient_temp"], color=colors[0], linestyle="--", marker="s", linewidth=2)

    ax1.set_ylabel("Contamination rate (%)")
    ax2.set_ylabel("Ambient temperature (°C)")
    ax1.set_xlabel("Month")
    ax1.set_title("Monthly contamination rate and temperature")
    ax1.tick_params(axis="x", rotation=45)

    fig.savefig(FIG_DIR / "fig6c_seasonal_contamination.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return monthly


def plot_intervention_impact(intervention_df: pd.DataFrame) -> None:
    colors = plt.cm.viridis(np.linspace(0.2, 0.85, len(intervention_df)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    x = np.arange(len(intervention_df))
    width = 0.35
    axes[0].bar(x - width / 2, intervention_df["baseline_rate"] * 100, width, label="Baseline", color=plt.cm.viridis(0.25))
    axes[0].bar(x + width / 2, intervention_df["after_rate"] * 100, width, label="After intervention", color=plt.cm.viridis(0.75))
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(intervention_df["intervention"], rotation=20, ha="right")
    axes[0].set_ylabel("Predicted contamination rate (%)")
    axes[0].set_title("Intervention effect on contamination risk")
    axes[0].legend(frameon=True)

    axes[1].barh(intervention_df["intervention"], intervention_df["net_benefit"] / 1000, color=colors)
    axes[1].axvline(0, linestyle="--", color="0.5", linewidth=1)
    axes[1].set_xlabel("Net benefit (thousand USD)")
    axes[1].set_title("Simple cost-benefit analysis")

    fig.savefig(FIG_DIR / "fig6d_intervention_impact.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_integrated_dashboard(monthly: pd.DataFrame, perm_df: pd.DataFrame, best_probs: np.ndarray,
                              y_true: pd.Series, intervention_df: pd.DataFrame, best_model_name: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    viridis = plt.cm.viridis

    axes[0, 0].plot(monthly["month"], monthly["contamination_rate"] * 100, color=viridis(0.8), marker="o", linewidth=2)
    twin = axes[0, 0].twinx()
    twin.plot(monthly["month"], monthly["ambient_temp"], color=viridis(0.25), linestyle="--", marker="s", linewidth=2)
    axes[0, 0].set_title("Monthly contamination trend")
    axes[0, 0].set_ylabel("Contamination rate (%)")
    twin.set_ylabel("Ambient temperature (°C)")
    axes[0, 0].tick_params(axis="x", rotation=45)

    top_perm = perm_df.head(8).sort_values("importance_mean")
    axes[0, 1].barh(
        [FEATURE_LABELS.get(f, f) for f in top_perm["feature"]],
        top_perm["importance_mean"],
        color=viridis(np.linspace(0.2, 0.85, len(top_perm))),
    )
    axes[0, 1].set_title("Ranked risk factors")
    axes[0, 1].set_xlabel("Permutation importance")

    calib = pd.DataFrame({"actual": y_true.to_numpy(), "predicted": best_probs})
    calib["decile"] = pd.qcut(calib["predicted"], q=10, labels=False, duplicates="drop")
    calib_summary = calib.groupby("decile").agg(actual_rate=("actual", "mean"), predicted_rate=("predicted", "mean")).reset_index()
    axes[1, 0].plot(calib_summary["predicted_rate"], calib_summary["actual_rate"], marker="o", color=viridis(0.75), linewidth=2)
    axes[1, 0].plot([0, calib_summary[["actual_rate", "predicted_rate"]].max().max()],
                    [0, calib_summary[["actual_rate", "predicted_rate"]].max().max()],
                    linestyle="--", color="0.4")
    axes[1, 0].set_xlabel("Predicted contamination probability")
    axes[1, 0].set_ylabel("Observed contamination rate")
    axes[1, 0].set_title(f"Predicted vs actual ({best_model_name})")

    axes[1, 1].bar(intervention_df["intervention"], intervention_df["absolute_reduction"] * 100,
                   color=viridis(np.linspace(0.25, 0.85, len(intervention_df))))
    axes[1, 1].set_ylabel("Absolute reduction (%)")
    axes[1, 1].set_title("Intervention effectiveness")
    axes[1, 1].tick_params(axis="x", rotation=20)
    for idx, value in enumerate(intervention_df["net_benefit"] / 1000):
        axes[1, 1].text(idx, intervention_df.loc[idx, "absolute_reduction"] * 100 + 0.1, f"${value:,.0f}k",
                        ha="center", va="bottom", fontsize=9)

    fig.savefig(FIG_DIR / "fig6e_integrated_dashboard.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def write_preprocessing_log(df: pd.DataFrame) -> None:
    overall_rate = df["salmonella_positive"].mean() * 100
    summer_rate = df.loc[df["season"] == "Summer", "salmonella_positive"].mean() * 100
    with PREPROCESS_LOG_PATH.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Preprocessing log\n\n"
            f"- Seed: {SEED}\n"
            f"- Records generated: {len(df)}\n"
            "- Sample window: 2023-01-01 to 2024-12-31\n"
            "- Synthetic features generated to mimic farm, processing, chilling, storage, and seasonality effects.\n"
            "- Count features were log10-transformed for modeling to stabilize skewed microbial distributions.\n"
            "- Interaction terms added: temperature × summer and water chilling × processing speed.\n"
            f"- Overall contamination rate: {overall_rate:.2f}%\n"
            f"- Summer contamination rate: {summer_rate:.2f}%\n"
            "- No missing values were introduced in the synthetic dataset.\n"
        )


def write_statistical_summary(or_df: pd.DataFrame, perm_df: pd.DataFrame, metrics: dict, best_model: str,
                              intervention_df: pd.DataFrame) -> None:
    significant = or_df[(or_df["feature"] != "const") & (or_df["q_value"] < 0.05)].copy()
    significant = significant.sort_values("q_value").head(8)
    top_perm = perm_df.head(8)
    with SUMMARY_PATH.open("w", encoding="utf-8") as handle:
        handle.write("# Statistical summary\n\n")
        handle.write(f"Best discriminating model by cross-validated AUC: **{best_model}**\n\n")
        handle.write("## Model metrics\n\n")
        for model_name, model_metrics in metrics.items():
            handle.write(
                f"- **{model_name}**: AUC={model_metrics['auc']:.3f}, sensitivity={model_metrics['sensitivity']:.3f}, "
                f"specificity={model_metrics['specificity']:.3f}, PPV={model_metrics['ppv']:.3f}, "
                f"CV AUC mean={model_metrics['cv_auc_mean']:.3f} ± {model_metrics['cv_auc_std']:.3f}\n"
            )
        handle.write("\n## Odds ratios with FDR-adjusted inference\n\n")
        for _, row in significant.iterrows():
            label = FEATURE_LABELS.get(row["feature"], row["feature"])
            handle.write(
                f"- **{label}**: OR={row['odds_ratio']:.3f} (95% CI {row['ci_lower']:.3f} to {row['ci_upper']:.3f}), "
                f"p={row['p_value']:.4g}, q={row['q_value']:.4g}\n"
            )
        handle.write("\n## Top permutation importance features\n\n")
        for _, row in top_perm.iterrows():
            label = FEATURE_LABELS.get(row["feature"], row["feature"])
            handle.write(
                f"- **{label}**: mean AUC decrease={row['importance_mean']:.4f} ± {row['importance_std']:.4f}\n"
            )
        best_intervention = intervention_df.iloc[0]
        handle.write("\n## Intervention summary\n\n")
        handle.write(
            f"- Highest absolute reduction: **{best_intervention['intervention']}** with "
            f"{best_intervention['absolute_reduction'] * 100:.2f} percentage points reduction and "
            f"net benefit of ${best_intervention['net_benefit']:,.0f}.\n"
        )


def write_report(df: pd.DataFrame, metrics: dict, best_model: str, or_df: pd.DataFrame,
                 intervention_df: pd.DataFrame) -> None:
    overall_rate = df["salmonella_positive"].mean() * 100
    summer_rate = df.loc[df["season"] == "Summer", "salmonella_positive"].mean() * 100
    winter_rate = df.loc[df["season"] == "Winter", "salmonella_positive"].mean() * 100
    strongest_or = or_df[(or_df["feature"] != "const")].copy().iloc[0]
    strongest_label = FEATURE_LABELS.get(strongest_or["feature"], strongest_or["feature"])
    best_intervention = intervention_df.iloc[0]

    report = f"""# DRAFT — NOT FOR DISTRIBUTION

# Salmonella contamination case study report

## Timestamp
- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Objective
This case study generated a realistic synthetic poultry supply chain dataset and evaluated predictive models for Salmonella contamination in chicken products. The analysis integrated predictive modeling, risk factor quantification, seasonal surveillance, and intervention simulation.

## Methods
- Created 2,000 lot-level synthetic samples across 2023-2024.
- Included farm, flock, environmental, processing, chilling, microbial count, storage, and retail timing features.
- Evaluated Logistic Regression, Random Forest, and Gradient Boosting with stratified 5-fold cross-validation.
- Estimated odds ratios and 95% confidence intervals with a binomial generalized linear model.
- Applied false discovery rate correction to the multi-parameter odds-ratio analysis.
- Estimated SHAP-like importance with permutation importance using the best-performing model.
- Simulated three interventions: improved chilling, reduced processing speed, and enhanced sanitation.

## Results
- Overall contamination rate: **{overall_rate:.2f}%**.
- Summer contamination rate: **{summer_rate:.2f}%**, compared with **{winter_rate:.2f}%** in winter.
- Best cross-validated discrimination: **{best_model}** with AUC **{metrics[best_model]['auc']:.3f}**.
- The strongest odds-ratio signal was **{strongest_label}** with OR **{strongest_or['odds_ratio']:.3f}** (95% CI {strongest_or['ci_lower']:.3f}-{strongest_or['ci_upper']:.3f}).
- The most effective simulated intervention was **{best_intervention['intervention']}**, reducing predicted contamination by **{best_intervention['absolute_reduction'] * 100:.2f}** percentage points.

## Discussion
The synthetic surveillance pattern reproduced the requested epidemiology: contamination was near 15% overall, higher in summer, and positively associated with ambient temperature. Post-chill and final-product microbial loads remained dominant decision signals, while interaction terms indicated that seasonal heat and water-based chilling at higher line speeds amplified risk. Intervention modeling suggested that processing upgrades can reduce risk materially, but economic value depends on implementation cost and the assumed cost of contaminated lots.

## Limitations
- The dataset is synthetic and should be used for training or demonstration rather than regulatory decision-making.
- Estimated cost-benefit results depend on a simple assumed cost-per-case model.
- Odds ratios reflect the synthetic generating process and not observed field surveillance.

## File inventory
- `data/salmonella_data.csv`
- `data/preprocessing-log.md`
- `results/module6_metrics.json`
- `results/module6_odds_ratios.csv`
- `results/module6_permutation_importance.csv`
- `results/module6_intervention_summary.csv`
- `results/statistical-summary.md`
- `figures/fig6_roc_curves.png`
- `figures/fig6b_risk_factors.png`
- `figures/fig6c_seasonal_contamination.png`
- `figures/fig6d_intervention_impact.png`
- `figures/fig6e_integrated_dashboard.png`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def save_metrics(metrics: dict) -> None:
    with METRICS_PATH.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)


def main() -> None:
    ensure_directories()
    append_log(
        "run_started",
        phase="PLAN",
        handoff_in={"task": "Create Salmonella contamination case study script"},
        handoff_out={"seed": SEED, "target_rate": TARGET_RATE},
    )
    append_log(
        "prompt_received",
        phase="PLAN",
        handoff_in={"requirements": "synthetic data generation, modeling, intervention simulation, dashboard"},
    )
    append_log(
        "skill_selected",
        phase="PLAN",
        handoff_out={"skill": "co-scientist-data-analysis"},
    )

    df = generate_synthetic_data()
    df.to_csv(DATA_PATH, index=False)
    write_preprocessing_log(df)
    append_log(
        "handoff_started",
        phase="EXECUTE",
        handoff_out={"records": len(df), "csv": str(DATA_PATH.relative_to(BASE_DIR))},
        files_written=[str(DATA_PATH.relative_to(BASE_DIR)), str(PREPROCESS_LOG_PATH.relative_to(BASE_DIR))],
    )

    df_model = prepare_model_data(df)
    X = df_model[MODEL_FEATURES]
    y = df_model["salmonella_positive"]

    models = build_models()
    metrics, curves = evaluate_models(X, y, models)
    best_model_name = max(metrics, key=lambda key: metrics[key]["auc"])

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, stratify=y, random_state=SEED)
    best_model = clone(models[best_model_name]).fit(X_train, y_train)
    perm_df = compute_permutation_importance(best_model, X_test, y_test)

    fitted_full_model = clone(models[best_model_name]).fit(X, y)
    or_df = compute_odds_ratios(df_model)
    intervention_df = simulate_interventions(df, fitted_full_model)
    monthly = plot_seasonal_contamination(df)
    plot_roc_and_pr(curves)
    plot_risk_factors(or_df, perm_df)
    plot_intervention_impact(intervention_df)
    plot_integrated_dashboard(monthly, perm_df, curves[best_model_name]["probabilities"], y, intervention_df, best_model_name)

    save_metrics(metrics)
    or_df.to_csv(OR_TABLE_PATH, index=False)
    perm_df.to_csv(PERM_PATH, index=False)
    intervention_df.to_csv(INTERVENTION_PATH, index=False)
    write_statistical_summary(or_df, perm_df, metrics, best_model_name, intervention_df)
    write_report(df, metrics, best_model_name, or_df, intervention_df)

    append_log(
        "handoff_completed",
        phase="VERIFY",
        handoff_out={
            "overall_contamination_rate": float(df["salmonella_positive"].mean()),
            "summer_rate": float(df.loc[df["season"] == "Summer", "salmonella_positive"].mean()),
            "best_model": best_model_name,
        },
        files_written=[
            str(METRICS_PATH.relative_to(BASE_DIR)),
            str(OR_TABLE_PATH.relative_to(BASE_DIR)),
            str(PERM_PATH.relative_to(BASE_DIR)),
            str(INTERVENTION_PATH.relative_to(BASE_DIR)),
        ],
    )
    for rel_path in [
        "figures/fig6_roc_curves.png",
        "figures/fig6b_risk_factors.png",
        "figures/fig6c_seasonal_contamination.png",
        "figures/fig6d_intervention_impact.png",
        "figures/fig6e_integrated_dashboard.png",
        "results/statistical-summary.md",
        "report.md",
    ]:
        append_log("file_written", phase="REPORT", files_written=[rel_path])
    append_log(
        "report_finalized",
        phase="REPORT",
        handoff_out={"report": "report.md", "summary": "results/statistical-summary.md"},
        files_written=["report.md", "results/statistical-summary.md"],
    )
    append_log(
        "run_completed",
        phase="LOG",
        handoff_out={"status": "completed"},
        files_written=["logs/process-log.jsonl"],
    )

    print(f"Saved dataset to {DATA_PATH}")
    print(f"Saved metrics to {METRICS_PATH}")
    print(f"Best model: {best_model_name}")
    print(f"Overall contamination rate: {df['salmonella_positive'].mean():.3f}")


if __name__ == "__main__":
    main()
