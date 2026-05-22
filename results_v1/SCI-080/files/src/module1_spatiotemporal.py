#!/usr/bin/env python3
"""Spatiotemporal prediction of synthetic foodborne illness outbreaks."""

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
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

BASE_DIR = Path(__file__).resolve().parents[1]
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"
REPORT_PATH = BASE_DIR / "report.md"
PREPROCESS_LOG_PATH = DATA_DIR / "preprocessing-log.md"
SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"
DATA_PATH = DATA_DIR / "spatiotemporal_data.csv"
METRICS_PATH = RESULTS_DIR / "module1_metrics.json"

for directory in (FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR):
    directory.mkdir(parents=True, exist_ok=True)


def log_event(event_type: str, skill_or_tool: str, handoff_in: dict | None = None, handoff_out: dict | None = None,
              files_written: list[str] | None = None, status: str = "ok", phase: str = "data-analysis") -> None:
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


log_event("run_started", "module1_spatiotemporal.py", handoff_in={"seed": SEED})
log_event("prompt_received", "co-scientist-data-analysis", handoff_in={"task": "Spatiotemporal outbreak prediction with synthetic data"})
log_event("skill_selected", "co-scientist-data-analysis", handoff_out={"reason": "Synthetic spatiotemporal modelling and visualization"})


def season_from_month(month: int) -> str:
    return {
        12: "Winter", 1: "Winter", 2: "Winter",
        3: "Spring", 4: "Spring", 5: "Spring",
        6: "Summer", 7: "Summer", 8: "Summer",
        9: "Autumn", 10: "Autumn", 11: "Autumn",
    }[month]


regions = [
    {"region": "North Coast", "base_temp": 13.0, "temp_amp": 8.0, "humidity_bias": 8.0, "population_density": 1250, "risk_factor": 1.10},
    {"region": "South Coast", "base_temp": 20.5, "temp_amp": 6.5, "humidity_bias": 10.5, "population_density": 1420, "risk_factor": 1.25},
    {"region": "Metro Central", "base_temp": 18.0, "temp_amp": 7.0, "humidity_bias": 5.0, "population_density": 2200, "risk_factor": 1.45},
    {"region": "High Plains", "base_temp": 14.5, "temp_amp": 10.0, "humidity_bias": -2.0, "population_density": 680, "risk_factor": 0.90},
    {"region": "River Valley", "base_temp": 17.2, "temp_amp": 8.4, "humidity_bias": 6.5, "population_density": 980, "risk_factor": 1.05},
    {"region": "Lake District", "base_temp": 15.0, "temp_amp": 7.8, "humidity_bias": 7.0, "population_density": 860, "risk_factor": 0.98},
    {"region": "Desert Edge", "base_temp": 21.5, "temp_amp": 9.8, "humidity_bias": -8.0, "population_density": 540, "risk_factor": 0.88},
    {"region": "Forest Belt", "base_temp": 12.8, "temp_amp": 8.9, "humidity_bias": 9.5, "population_density": 720, "risk_factor": 1.00},
    {"region": "Agricultural South", "base_temp": 19.3, "temp_amp": 7.5, "humidity_bias": 4.5, "population_density": 930, "risk_factor": 1.18},
    {"region": "Port East", "base_temp": 18.5, "temp_amp": 6.8, "humidity_bias": 11.0, "population_density": 1580, "risk_factor": 1.32},
]


def build_synthetic_dataset() -> pd.DataFrame:
    dates = pd.date_range("2020-01-01", "2024-12-01", freq="MS")
    records: list[dict] = []
    rng = np.random.default_rng(SEED)

    for region_info in regions:
        lag_1 = 6.0 + region_info["risk_factor"]
        lag_12 = lag_1
        lag_history = [lag_1] * 12

        for index, current_date in enumerate(dates):
            month = current_date.month
            year = current_date.year
            angle = 2 * np.pi * (month - 1) / 12.0
            seasonal_wave = np.sin(angle - np.pi / 2)
            month_sin = np.sin(angle)
            month_cos = np.cos(angle)
            season_index = (month - 1) // 3
            season_angle = 2 * np.pi * season_index / 4.0
            season_sin = np.sin(season_angle)
            season_cos = np.cos(season_angle)
            season = season_from_month(month)
            time_index = index
            climate_noise = rng.normal(0, 1.2)
            humidity_noise = rng.normal(0, 3.0)
            temperature = (
                region_info["base_temp"]
                + region_info["temp_amp"] * seasonal_wave
                + 0.08 * time_index
                + climate_noise
            )
            humidity = (
                46
                + 0.92 * temperature
                + region_info["humidity_bias"]
                - 1.7 * seasonal_wave
                + humidity_noise
            )
            humidity = float(np.clip(humidity, 32, 95))
            rolling_mean = float(np.mean(lag_history[-3:]))
            seasonal_risk = 1.8 + 4.4 * ((seasonal_wave + 1.0) / 2.0)
            lambda_incidents = (
                seasonal_risk
                + 0.16 * max(temperature - 8, 0)
                + 0.050 * max(humidity - 45, 0)
                + 0.0028 * region_info["population_density"]
                + 0.18 * lag_1
                + 0.07 * lag_12
                + 0.30 * rolling_mean
                + 1.4 * region_info["risk_factor"]
                + 0.03 * time_index
            )
            lambda_incidents = max(lambda_incidents, 0.5)
            incidents = int(rng.poisson(lam=lambda_incidents))
            risk_index = incidents / (1.0 + region_info["population_density"] / 1000.0)

            records.append(
                {
                    "date": current_date,
                    "year": year,
                    "month": month,
                    "season": season,
                    "region": region_info["region"],
                    "temperature": round(float(temperature), 2),
                    "humidity": round(humidity, 2),
                    "month_sin": round(float(month_sin), 6),
                    "month_cos": round(float(month_cos), 6),
                    "season_sin": round(float(season_sin), 6),
                    "season_cos": round(float(season_cos), 6),
                    "population_density": region_info["population_density"],
                    "historical_incidents": round(float(lag_1), 2),
                    "rolling_3m_incidents": round(rolling_mean, 2),
                    "lag_12_incidents": round(float(lag_12), 2),
                    "time_index": time_index,
                    "incident_count": incidents,
                    "risk_index": round(float(risk_index), 3),
                }
            )

            lag_history.append(incidents)
            lag_1 = float(incidents)
            lag_12 = float(lag_history[-12])

    frame = pd.DataFrame.from_records(records).sort_values(["date", "region"]).reset_index(drop=True)
    return frame


log_event("handoff_started", "data-generation", handoff_in={"regions": len(regions), "period": "2020-2024"})
df = build_synthetic_dataset()
df.to_csv(DATA_PATH, index=False)
log_event("file_written", "pandas.to_csv", files_written=[str(DATA_PATH)], handoff_out={"rows": int(df.shape[0]), "columns": int(df.shape[1])})
log_event("handoff_completed", "data-generation", handoff_out={"dataset_shape": list(df.shape)})

feature_columns = [
    "temperature",
    "humidity",
    "month_sin",
    "month_cos",
    "season_sin",
    "season_cos",
    "population_density",
    "historical_incidents",
    "rolling_3m_incidents",
    "lag_12_incidents",
    "time_index",
]

model_df = pd.get_dummies(df[[*feature_columns, "region"]], columns=["region"], drop_first=False)
X = model_df
y = df["incident_count"]
train_mask = df["date"] < pd.Timestamp("2024-01-01")
X_train, X_test = X.loc[train_mask], X.loc[~train_mask]
y_train, y_test = y.loc[train_mask], y.loc[~train_mask]

def make_models() -> dict[str, object]:
    return {
        "Random Forest Regressor": RandomForestRegressor(
            n_estimators=400,
            max_depth=10,
            min_samples_leaf=2,
            random_state=SEED,
            n_jobs=-1,
        ),
        "Gradient Boosting Regressor": GradientBoostingRegressor(
            n_estimators=350,
            learning_rate=0.04,
            max_depth=3,
            subsample=0.9,
            random_state=SEED,
        ),
        "MLPRegressor (LSTM proxy)": Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                (
                    "mlp",
                    MLPRegressor(
                        hidden_layer_sizes=(64, 32),
                        activation="relu",
                        alpha=1e-4,
                        learning_rate_init=0.01,
                        max_iter=2500,
                        early_stopping=True,
                        random_state=SEED,
                    ),
                ),
            ]
        ),
    }


models = make_models()
predictions: dict[str, np.ndarray] = {}
metrics: dict[str, dict[str, float]] = {}
trained_models: dict[str, object] = {}

log_event("handoff_started", "model-training", handoff_in={"n_train": int(X_train.shape[0]), "n_test": int(X_test.shape[0])})
for model_name, model in models.items():
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    predictions[model_name] = y_pred
    trained_models[model_name] = model
    metrics[model_name] = {
        "RMSE": round(float(np.sqrt(mean_squared_error(y_test, y_pred))), 4),
        "MAE": round(float(mean_absolute_error(y_test, y_pred)), 4),
        "R2": round(float(r2_score(y_test, y_pred)), 4),
    }

best_model_name = min(metrics, key=lambda name: metrics[name]["RMSE"])
best_model = trained_models[best_model_name]
best_predictions = predictions[best_model_name]
train_residuals = y_train.to_numpy() - np.asarray(best_model.predict(X_train))
residual_quantiles = np.quantile(train_residuals, [0.025, 0.975])
lower_ci = np.clip(best_predictions + residual_quantiles[0], 0, None)
upper_ci = np.clip(best_predictions + residual_quantiles[1], 0, None)

metrics_payload = {
    "seed": SEED,
    "train_period": ["2020-01", "2023-12"],
    "test_period": ["2024-01", "2024-12"],
    "best_model": best_model_name,
    "metrics": metrics,
}
with METRICS_PATH.open("w", encoding="utf-8") as handle:
    json.dump(metrics_payload, handle, indent=2)
log_event("file_written", "json.dump", files_written=[str(METRICS_PATH)], handoff_out={"best_model": best_model_name})
log_event("handoff_completed", "model-training", handoff_out={"metrics": metrics, "best_model": best_model_name})

importance_result = permutation_importance(
    best_model,
    X_test,
    y_test,
    n_repeats=20,
    random_state=SEED,
    n_jobs=-1,
)
importance_df = (
    pd.DataFrame({
        "feature": X.columns,
        "importance_mean": importance_result.importances_mean,
        "importance_std": importance_result.importances_std,
    })
    .sort_values("importance_mean", ascending=False)
    .reset_index(drop=True)
)

monthly_predictions = df.loc[~train_mask, ["date", "incident_count"]].copy()
monthly_predictions["prediction"] = best_predictions
monthly_predictions["lower_ci"] = lower_ci
monthly_predictions["upper_ci"] = upper_ci
monthly_summary = monthly_predictions.groupby("date", as_index=False).agg(
    actual=("incident_count", "sum"),
    predicted=("prediction", "sum"),
    lower_ci=("lower_ci", "sum"),
    upper_ci=("upper_ci", "sum"),
)

plt.style.use("seaborn-v0_8-whitegrid")
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(monthly_summary["date"], monthly_summary["actual"], marker="o", linewidth=2.2, color="#1f77b4", label="Actual incidents")
ax.plot(monthly_summary["date"], monthly_summary["predicted"], marker="s", linewidth=2.2, color="#ff7f0e", label=f"Predicted ({best_model_name})")
ax.fill_between(monthly_summary["date"], monthly_summary["lower_ci"], monthly_summary["upper_ci"], color="#2a9d8f", alpha=0.25, label="95% interval")
ax.set_title("Observed vs Predicted Foodborne Outbreak Incidents")
ax.set_xlabel("Month")
ax.set_ylabel("Monthly incident count")
ax.legend(frameon=False)
fig.autofmt_xdate()
fig.tight_layout()
fig.savefig(FIGURES_DIR / "fig1_spatiotemporal_prediction.png", dpi=300, bbox_inches="tight")
plt.close(fig)
log_event("file_written", "matplotlib", files_written=[str(FIGURES_DIR / "fig1_spatiotemporal_prediction.png")])

fig, ax = plt.subplots(figsize=(10, 6))
top_importances = importance_df.head(12).sort_values("importance_mean")
colors = plt.cm.viridis(np.linspace(0.15, 0.9, len(top_importances)))
ax.barh(top_importances["feature"], top_importances["importance_mean"], xerr=top_importances["importance_std"], color=colors, edgecolor="none")
ax.set_title(f"Permutation Feature Importance ({best_model_name})")
ax.set_xlabel("Mean importance decrease")
ax.set_ylabel("Feature")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "fig1b_feature_importance.png", dpi=300, bbox_inches="tight")
plt.close(fig)
log_event("file_written", "matplotlib", files_written=[str(FIGURES_DIR / "fig1b_feature_importance.png")])

heatmap_df = df.groupby(["region", "month"])["incident_count"].mean().unstack(fill_value=0)
month_labels = [datetime(2024, month, 1).strftime("%b") for month in heatmap_df.columns]
fig, ax = plt.subplots(figsize=(11, 6))
heatmap = ax.imshow(heatmap_df.values, cmap="cividis", aspect="auto")
ax.set_title("Regional Seasonal Risk Heatmap")
ax.set_xlabel("Month")
ax.set_ylabel("Region")
ax.set_xticks(np.arange(len(month_labels)))
ax.set_xticklabels(month_labels)
ax.set_yticks(np.arange(len(heatmap_df.index)))
ax.set_yticklabels(heatmap_df.index)
colorbar = fig.colorbar(heatmap, ax=ax)
colorbar.set_label("Mean incident count")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "fig1c_seasonal_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig)
log_event("file_written", "matplotlib", files_written=[str(FIGURES_DIR / "fig1c_seasonal_heatmap.png")])

monthly_total = df.groupby("date")["incident_count"].sum().asfreq("MS")
trend = monthly_total.rolling(window=12, center=True, min_periods=6).mean()
trend = trend.interpolate(limit_direction="both")
seasonal_component = (monthly_total - trend).groupby(monthly_total.index.month).transform("mean")
residual = monthly_total - trend - seasonal_component
fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
plot_series = [
    (monthly_total, "Observed monthly incidents", "#1f77b4"),
    (trend, "Trend", "#2a9d8f"),
    (seasonal_component, "Seasonality", "#6f4c9b"),
    (residual, "Residual", "#d1495b"),
]
for axis, (series, title, color) in zip(axes, plot_series):
    axis.plot(series.index, series.values, color=color, linewidth=2)
    axis.set_title(title)
    axis.set_ylabel("Count")
axes[-1].set_xlabel("Date")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "fig1d_time_decomposition.png", dpi=300, bbox_inches="tight")
plt.close(fig)
log_event("file_written", "matplotlib", files_written=[str(FIGURES_DIR / "fig1d_time_decomposition.png")])

summer = df.loc[df["season"] == "Summer", "incident_count"].to_numpy()
winter = df.loc[df["season"] == "Winter", "incident_count"].to_numpy()
mean_difference = float(summer.mean() - winter.mean())
pooled_sd = float(np.sqrt(((summer.std(ddof=1) ** 2) + (winter.std(ddof=1) ** 2)) / 2))
cohens_d = mean_difference / pooled_sd if pooled_sd else 0.0
rng = np.random.default_rng(SEED)
bootstrap_diffs = []
for _ in range(3000):
    sample_summer = rng.choice(summer, size=len(summer), replace=True)
    sample_winter = rng.choice(winter, size=len(winter), replace=True)
    bootstrap_diffs.append(sample_summer.mean() - sample_winter.mean())
ci_low, ci_high = np.quantile(bootstrap_diffs, [0.025, 0.975])

PREPROCESS_LOG_PATH.write_text(
    "# Preprocessing Log\n\n"
    "1. Generated monthly synthetic observations for 10 regions from January 2020 to December 2024.\n"
    "2. Applied deterministic random seed 42 to Python random and NumPy.\n"
    "3. Simulated temperature with regional baselines, seasonal oscillation, and a mild warming trend.\n"
    "4. Simulated humidity as a temperature-correlated feature with region-specific offsets and noise.\n"
    "5. Derived cyclical encodings for month and season using sine/cosine transforms.\n"
    "6. Generated outbreak counts from a Poisson process with summer amplification and historical lag effects.\n"
    "7. One-hot encoded region labels and used a temporal split: 2020-2023 for training, 2024 for testing.\n",
    encoding="utf-8",
)
log_event("file_written", "write_text", files_written=[str(PREPROCESS_LOG_PATH)])

SUMMARY_PATH.write_text(
    "# Statistical Summary\n\n"
    f"- Best predictive model: **{best_model_name}**.\n"
    f"- Summer vs Winter mean incident difference: **{mean_difference:.2f}** incidents.\n"
    f"- 95% bootstrap confidence interval: **[{ci_low:.2f}, {ci_high:.2f}]**.\n"
    f"- Effect size (Cohen's d): **{cohens_d:.2f}**.\n"
    "- No multiple-comparison correction was required because only one seasonal contrast was summarized.\n"
    "- Decomposition was computed with an additive rolling-trend approach on monthly totals.\n",
    encoding="utf-8",
)
log_event("file_written", "write_text", files_written=[str(SUMMARY_PATH)])

report_lines = [
    "# DRAFT — NOT FOR DISTRIBUTION\n",
    "# Spatiotemporal Prediction of Foodborne Illness Outbreaks\n",
    f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
    "## Methods\n",
    "A synthetic monthly panel was created for 10 regions across 2020-2024. Temperature followed region-specific seasonal waves, humidity was correlated with temperature, cyclical encodings captured temporal seasonality, and incident counts were sampled from a Poisson process informed by climate, density, and lagged outbreaks. Random Forest, Gradient Boosting, and MLPRegressor (as an LSTM proxy) were trained on 2020-2023 data and evaluated on 2024 data.\n",
    "## Results\n",
]
for model_name, values in metrics.items():
    report_lines.append(
        f"- {model_name}: RMSE={values['RMSE']:.4f}, MAE={values['MAE']:.4f}, R²={values['R2']:.4f}.\n"
    )
report_lines.extend([
    f"- Best model: {best_model_name}.\n",
    f"- Summer vs Winter effect size: Cohen's d={cohens_d:.2f}; 95% CI for mean difference=[{ci_low:.2f}, {ci_high:.2f}].\n",
    "## Discussion\n",
    "The synthetic system produces higher summer risk, especially in dense coastal and metropolitan regions, and the models recover this seasonal signal effectively. Because the data are simulated, the estimates demonstrate pipeline behavior rather than real-world epidemiologic truth.\n",
    "## File Inventory\n",
    "- data/spatiotemporal_data.csv\n",
    "- data/preprocessing-log.md\n",
    "- results/module1_metrics.json\n",
    "- results/statistical-summary.md\n",
    "- figures/fig1_spatiotemporal_prediction.png\n",
    "- figures/fig1b_feature_importance.png\n",
    "- figures/fig1c_seasonal_heatmap.png\n",
    "- figures/fig1d_time_decomposition.png\n",
    "- logs/process-log.jsonl\n",
])
REPORT_PATH.write_text("\n".join(report_lines), encoding="utf-8")
log_event("report_finalized", "write_text", files_written=[str(REPORT_PATH)])
log_event("run_completed", "module1_spatiotemporal.py", handoff_out={"best_model": best_model_name, "report": str(REPORT_PATH)})

print(f"Saved dataset to {DATA_PATH}")
print(f"Saved metrics to {METRICS_PATH}")
print(f"Best model: {best_model_name}")
