"""
Component 3: Psychometric paradigm model for risk perception (Slovic framework).
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import FactorAnalysis
from sklearn.preprocessing import StandardScaler
import json

np.random.seed(42)

RISK_DIMENSIONS = {
    "dread": ["controllability", "dread", "catastrophic_potential", "fatal_consequences",
              "inequitable_risk", "involuntary_exposure", "future_generation_risk"],
    "unknown": ["observability", "knowledge_exposed", "knowledge_science", "novelty", "delayed_effects"]
}

def generate_risk_perception_data(n_respondents=500):
    technologies = {
        "gene_editing": {"dread_center": [4.5, 4.8, 4.2, 3.5, 4.0, 3.8, 5.5],
                         "unknown_center": [4.0, 4.5, 3.8, 5.2, 4.5]},
        "AI": {"dread_center": [3.8, 4.0, 4.5, 3.0, 3.5, 3.2, 4.0],
               "unknown_center": [3.5, 4.0, 3.5, 4.8, 4.0]},
        "nuclear_fusion": {"dread_center": [4.0, 4.5, 5.0, 4.0, 3.8, 4.2, 4.5],
                           "unknown_center": [4.5, 5.0, 4.2, 5.5, 4.8]}
    }
    all_data = []
    for tech, centers in technologies.items():
        for i in range(n_respondents):
            row = {"respondent_id": f"R{i:04d}", "technology": tech,
                   "age_group": np.random.choice(["18-29", "30-44", "45-59", "60+"], p=[0.25, 0.30, 0.25, 0.20]),
                   "education": np.random.choice(["high_school", "bachelor", "master", "doctorate"], p=[0.30, 0.35, 0.25, 0.10]),
                   "gender": np.random.choice(["male", "female", "other"], p=[0.48, 0.48, 0.04])}
            for j, dim in enumerate(RISK_DIMENSIONS["dread"]):
                row[dim] = int(np.clip(np.random.normal(centers["dread_center"][j], 1.0), 1, 7))
            for j, dim in enumerate(RISK_DIMENSIONS["unknown"]):
                row[dim] = int(np.clip(np.random.normal(centers["unknown_center"][j], 1.2), 1, 7))
            dread_avg = np.mean([row[d] for d in RISK_DIMENSIONS["dread"]])
            unknown_avg = np.mean([row[d] for d in RISK_DIMENSIONS["unknown"]])
            row["overall_risk_perception"] = int(np.clip(0.6 * dread_avg + 0.4 * unknown_avg + np.random.normal(0, 0.5), 1, 7))
            row["acceptance"] = int(np.clip(8 - row["overall_risk_perception"] + np.random.normal(0, 0.8), 1, 7))
            all_data.append(row)
    return pd.DataFrame(all_data)

def factor_analysis_psychometric(df):
    all_items = RISK_DIMENSIONS["dread"] + RISK_DIMENSIONS["unknown"]
    X = df[all_items].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    fa = FactorAnalysis(n_components=2, random_state=42)
    fa.fit(X_scaled)
    loadings = pd.DataFrame(fa.components_.T, columns=["Factor1_Dread", "Factor2_Unknown"], index=all_items)
    factor_scores = fa.transform(X_scaled)
    df = df.copy()
    df["dread_factor"] = factor_scores[:, 0]
    df["unknown_factor"] = factor_scores[:, 1]
    return loadings, df

def run_psychometric_analysis():
    df = generate_risk_perception_data()
    df.to_csv("data/psychometric_risk_data.csv", index=False)
    loadings, df_with_factors = factor_analysis_psychometric(df)
    loadings.to_csv("results/factor_loadings.csv")
    results = {}
    for tech in df["technology"].unique():
        sub = df_with_factors[df_with_factors["technology"] == tech]
        r, p = stats.pearsonr(sub["overall_risk_perception"], sub["acceptance"])
        demo_effects = {}
        for demo in ["age_group", "education", "gender"]:
            groups = sub.groupby(demo)["overall_risk_perception"].mean()
            f_stat, f_p = stats.f_oneway(*[sub[sub[demo] == g]["overall_risk_perception"].values for g in sub[demo].unique()])
            demo_effects[demo] = {
                "group_means": {k: round(v, 3) for k, v in groups.items()},
                "F_statistic": round(float(f_stat), 3), "p_value": round(float(f_p), 4)
            }
        results[tech] = {
            "dread_factor_mean": round(float(sub["dread_factor"].mean()), 4),
            "unknown_factor_mean": round(float(sub["unknown_factor"].mean()), 4),
            "risk_acceptance_correlation": round(float(r), 4),
            "risk_acceptance_p_value": round(float(p), 6),
            "mean_risk_perception": round(float(sub["overall_risk_perception"].mean()), 3),
            "mean_acceptance": round(float(sub["acceptance"].mean()), 3),
            "demographic_effects": demo_effects
        }
    with open("results/psychometric_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return df_with_factors, loadings, results

if __name__ == "__main__":
    df, loadings, results = run_psychometric_analysis()
    print("Psychometric analysis completed.")
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = results[tech]
        print(f"  {tech}: dread={r['dread_factor_mean']:.3f}, unknown={r['unknown_factor_mean']:.3f}, r={r['risk_acceptance_correlation']:.3f}")
