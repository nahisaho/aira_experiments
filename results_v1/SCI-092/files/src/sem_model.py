"""
Component 5: Trust-Acceptance causal model using SEM path analysis.
"""

import numpy as np
import pandas as pd
import semopy
import json

np.random.seed(42)

def generate_sem_data(n=800):
    scientific_literacy = np.random.normal(0, 1, n)
    media_exposure = np.random.normal(0, 1, n)
    institutional_trust = 0.45 * scientific_literacy + np.random.normal(0, 0.7, n)
    perceived_risk = -0.35 * institutional_trust + 0.30 * media_exposure + np.random.normal(0, 0.6, n)
    perceived_benefit = 0.40 * institutional_trust + 0.25 * scientific_literacy - 0.15 * media_exposure + np.random.normal(0, 0.6, n)
    acceptance = 0.20 * institutional_trust - 0.35 * perceived_risk + 0.45 * perceived_benefit + np.random.normal(0, 0.5, n)

    data = pd.DataFrame()
    data["trust1"] = institutional_trust + np.random.normal(0, 0.4, n)
    data["trust2"] = 0.9 * institutional_trust + np.random.normal(0, 0.4, n)
    data["trust3"] = 0.85 * institutional_trust + np.random.normal(0, 0.45, n)
    data["risk1"] = perceived_risk + np.random.normal(0, 0.35, n)
    data["risk2"] = 0.92 * perceived_risk + np.random.normal(0, 0.38, n)
    data["risk3"] = 0.88 * perceived_risk + np.random.normal(0, 0.4, n)
    data["benefit1"] = perceived_benefit + np.random.normal(0, 0.35, n)
    data["benefit2"] = 0.93 * perceived_benefit + np.random.normal(0, 0.36, n)
    data["benefit3"] = 0.87 * perceived_benefit + np.random.normal(0, 0.42, n)
    data["accept1"] = acceptance + np.random.normal(0, 0.3, n)
    data["accept2"] = 0.91 * acceptance + np.random.normal(0, 0.35, n)
    data["accept3"] = 0.86 * acceptance + np.random.normal(0, 0.38, n)
    data["literacy1"] = scientific_literacy + np.random.normal(0, 0.3, n)
    data["literacy2"] = 0.9 * scientific_literacy + np.random.normal(0, 0.35, n)
    data["media1"] = media_exposure + np.random.normal(0, 0.35, n)
    data["media2"] = 0.88 * media_exposure + np.random.normal(0, 0.4, n)
    return data

def fit_sem_model(data):
    model_spec = """
    Trust =~ trust1 + trust2 + trust3
    Risk =~ risk1 + risk2 + risk3
    Benefit =~ benefit1 + benefit2 + benefit3
    Acceptance =~ accept1 + accept2 + accept3
    Literacy =~ literacy1 + literacy2
    Media =~ media1 + media2
    Trust ~ Literacy
    Risk ~ Trust + Media
    Benefit ~ Trust + Literacy + Media
    Acceptance ~ Trust + Risk + Benefit
    """
    model = semopy.Model(model_spec)
    model.fit(data)
    estimates = model.inspect()
    fit_stats = semopy.calc_stats(model)
    return model, estimates, fit_stats

def _safe_float(val, default=0.0):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default

def extract_sem_results(estimates, fit_stats):
    structural_paths = estimates[estimates["op"] == "~"].copy()
    path_results = []
    for _, row in structural_paths.iterrows():
        path_results.append({
            "from": row["rval"], "to": row["lval"],
            "estimate": round(_safe_float(row["Estimate"]), 4),
            "std_error": round(_safe_float(row["Std. Err"]), 4),
            "z_value": round(_safe_float(row["z-value"]), 4),
            "p_value": round(_safe_float(row["p-value"], 1.0), 6),
            "significant": _safe_float(row["p-value"], 1.0) < 0.05
        })
    loadings = estimates[estimates["op"] == "=~"].copy()
    loading_results = []
    for _, row in loadings.iterrows():
        loading_results.append({
            "latent": row["lval"], "indicator": row["rval"],
            "loading": round(_safe_float(row["Estimate"]), 4),
            "std_error": round(_safe_float(row["Std. Err"]), 4),
            "p_value": round(_safe_float(row["p-value"], 1.0), 6)
        })
    fit_dict = {}
    for col in fit_stats.columns:
        val = fit_stats[col].values[0]
        if isinstance(val, (int, float, np.integer, np.floating)):
            fit_dict[col] = round(float(val), 4)
    fit_interpretation = {
        "CFI": "good" if fit_dict.get("CFI", 0) > 0.95 else "acceptable" if fit_dict.get("CFI", 0) > 0.90 else "poor",
        "RMSEA": "good" if fit_dict.get("RMSEA", 1) < 0.06 else "acceptable" if fit_dict.get("RMSEA", 1) < 0.08 else "poor",
        "GFI": "good" if fit_dict.get("GFI", 0) > 0.95 else "acceptable" if fit_dict.get("GFI", 0) > 0.90 else "poor"
    }
    return {"structural_paths": path_results, "factor_loadings": loading_results,
            "fit_indices": fit_dict, "fit_interpretation": fit_interpretation}

def run_sem_analysis():
    data = generate_sem_data()
    data.to_csv("data/sem_data.csv", index=False)
    model, estimates, fit_stats = fit_sem_model(data)
    results = extract_sem_results(estimates, fit_stats)
    with open("results/sem_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return data, model, results

if __name__ == "__main__":
    data, model, results = run_sem_analysis()
    print("SEM analysis completed.")
    for p in results["structural_paths"]:
        sig = "***" if p["p_value"] < 0.001 else "**" if p["p_value"] < 0.01 else "*" if p["p_value"] < 0.05 else "ns"
        print(f"  {p['from']} → {p['to']}: β={p['estimate']:.3f} ({sig})")
    print(f"\nFit: CFI={results['fit_indices'].get('CFI', 'N/A')}, RMSEA={results['fit_indices'].get('RMSEA', 'N/A')}")
