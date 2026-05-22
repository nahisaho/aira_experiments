"""
Analysis modules: age acceleration, tissue specificity,
intervention sensitivity, and longevity validation.
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import mean_absolute_error, r2_score


def compute_age_acceleration(chronological_age, predicted_age):
    slope, intercept, _, _, _ = stats.linregress(chronological_age, predicted_age)
    expected = slope * chronological_age + intercept
    return predicted_age - expected


def age_acceleration_analysis(df, pred_col="predicted_age"):
    accel = compute_age_acceleration(df["chronological_age"].values, df[pred_col].values)
    df = df.copy()
    df["age_acceleration"] = accel
    results = {
        "mean_acceleration": round(float(np.mean(accel)), 3),
        "std_acceleration": round(float(np.std(accel)), 3),
        "range": [round(float(np.min(accel)), 3), round(float(np.max(accel)), 3)],
    }
    if "biological_age_offset" in df.columns:
        r, p = stats.pearsonr(accel, df["biological_age_offset"].values)
        results["correlation_with_true_offset"] = round(r, 4)
        results["correlation_p_value"] = float(f"{p:.2e}")
    if "sex" in df.columns:
        male_accel = accel[df["sex"] == 1]
        female_accel = accel[df["sex"] == 0]
        t, p = stats.ttest_ind(male_accel, female_accel)
        results["sex_difference"] = {
            "male_mean": round(float(np.mean(male_accel)), 3),
            "female_mean": round(float(np.mean(female_accel)), 3),
            "t_statistic": round(float(t), 3),
            "p_value": float(f"{p:.4f}"),
        }
    return df, results


def tissue_specificity_analysis(tissue_predictions: dict):
    results = {}
    for tissue, data in tissue_predictions.items():
        y_true, y_pred = data["true"], data["pred"]
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)
        r, _ = stats.pearsonr(y_true, y_pred)
        results[tissue] = {"MAE": round(mae, 3), "R2": round(r2, 4),
                           "Pearson_r": round(r, 4), "n_samples": len(y_true)}
    return results


def intervention_sensitivity(df, pred_col="predicted_age"):
    df = df.copy()
    accel = compute_age_acceleration(df["chronological_age"].values, df[pred_col].values)
    df["age_acceleration"] = accel
    control = df[df["intervention"] == "none"]["age_acceleration"]
    results = {"control_mean_accel": round(float(control.mean()), 3)}
    for intv in df["intervention"].unique():
        if intv == "none":
            continue
        treated = df[df["intervention"] == intv]["age_acceleration"]
        t, p = stats.ttest_ind(control, treated)
        cohens_d = (control.mean() - treated.mean()) / np.sqrt(
            (control.std()**2 + treated.std()**2) / 2)
        results[intv] = {
            "mean_accel": round(float(treated.mean()), 3),
            "delta_vs_control": round(float(control.mean() - treated.mean()), 3),
            "cohens_d": round(float(cohens_d), 3),
            "t_statistic": round(float(t), 3),
            "p_value": float(f"{p:.4f}"),
            "n_samples": len(treated),
            "detectable": bool(p < 0.05),
        }
    return results


def longevity_validation(df_longevity, df_normal, pred_col="predicted_age"):
    accel_long = compute_age_acceleration(
        df_longevity["chronological_age"].values, df_longevity[pred_col].values)
    accel_norm = compute_age_acceleration(
        df_normal["chronological_age"].values, df_normal[pred_col].values)
    t, p = stats.ttest_ind(accel_long, accel_norm)
    cohens_d = (np.mean(accel_norm) - np.mean(accel_long)) / np.sqrt(
        (np.std(accel_norm)**2 + np.std(accel_long)**2) / 2)
    return {
        "longevity_mean_accel": round(float(np.mean(accel_long)), 3),
        "normal_mean_accel": round(float(np.mean(accel_norm)), 3),
        "t_statistic": round(float(t), 3),
        "p_value": float(f"{p:.4f}"),
        "cohens_d": round(float(cohens_d), 3),
        "longevity_n": len(accel_long),
        "normal_n": len(accel_norm),
    }
