"""
Component 1: Meta-analysis framework for public opinion survey data
on emerging technologies (gene editing, AI, nuclear fusion).
"""

import numpy as np
import pandas as pd
from scipy import stats
import json

np.random.seed(42)

def generate_synthetic_studies():
    technologies = ["gene_editing", "AI", "nuclear_fusion"]
    countries = ["Japan", "USA", "UK", "Germany", "China", "South Korea"]
    years = list(range(2015, 2026))
    studies = []
    base_acceptance = {"gene_editing": 0.38, "AI": 0.55, "nuclear_fusion": 0.50}
    country_offset = {"Japan": -0.05, "USA": 0.05, "UK": 0.02,
                      "Germany": -0.02, "China": 0.08, "South Korea": 0.0}
    year_trend = {t: 0.012 for t in technologies}
    year_trend["AI"] = 0.018
    study_id = 1
    for tech in technologies:
        for country in countries:
            for year in np.random.choice(years, size=np.random.randint(2, 5), replace=False):
                n = np.random.randint(300, 3000)
                true_p = np.clip(
                    base_acceptance[tech] + country_offset[country] + year_trend[tech] * (year - 2020)
                    + np.random.normal(0, 0.04), 0.1, 0.9)
                observed_p = np.clip(np.random.binomial(n, true_p) / n, 0.05, 0.95)
                se = np.sqrt(observed_p * (1 - observed_p) / n)
                studies.append({
                    "study_id": f"S{study_id:03d}", "technology": tech,
                    "country": country, "year": int(year), "n": n,
                    "acceptance_rate": round(observed_p, 4),
                    "se": round(se, 4),
                    "quality_score": round(np.random.uniform(0.5, 1.0), 2)
                })
                study_id += 1
    return pd.DataFrame(studies)

def random_effects_meta_analysis(effects, variances):
    k = len(effects)
    effects = np.array(effects)
    variances = np.array(variances)
    w_fixed = 1.0 / variances
    Q = np.sum(w_fixed * (effects - np.average(effects, weights=w_fixed))**2)
    df = k - 1
    C = np.sum(w_fixed) - np.sum(w_fixed**2) / np.sum(w_fixed)
    tau2 = max(0, (Q - df) / C)
    w_random = 1.0 / (variances + tau2)
    pooled = np.average(effects, weights=w_random)
    se_pooled = np.sqrt(1.0 / np.sum(w_random))
    ci_lower = pooled - 1.96 * se_pooled
    ci_upper = pooled + 1.96 * se_pooled
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    return {
        "pooled_effect": round(pooled, 4), "se": round(se_pooled, 4),
        "ci_95": [round(ci_lower, 4), round(ci_upper, 4)],
        "tau2": round(tau2, 6), "I2": round(I2, 2),
        "Q": round(Q, 2), "Q_p_value": round(1 - stats.chi2.cdf(Q, df), 4), "k": k
    }

def run_meta_analysis():
    df = generate_synthetic_studies()
    df.to_csv("data/meta_analysis_studies.csv", index=False)
    results = {}
    for tech in df["technology"].unique():
        sub = df[df["technology"] == tech]
        ma = random_effects_meta_analysis(sub["acceptance_rate"].values, sub["se"].values**2)
        results[tech] = ma
        subgroup = {}
        for country in sub["country"].unique():
            csub = sub[sub["country"] == country]
            if len(csub) >= 2:
                subgroup[country] = random_effects_meta_analysis(
                    csub["acceptance_rate"].values, csub["se"].values**2)
        results[tech + "_by_country"] = subgroup

    for tech in df["technology"].unique():
        sub = df[df["technology"] == tech]
        slope, intercept, r_value, p_value, std_err = stats.linregress(sub["year"], sub["acceptance_rate"])
        results[tech + "_year_trend"] = {
            "slope": round(slope, 5), "intercept": round(intercept, 4),
            "r_squared": round(r_value**2, 4), "p_value": round(p_value, 4)
        }
    with open("results/meta_analysis_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return df, results

if __name__ == "__main__":
    df, results = run_meta_analysis()
    print("Meta-analysis completed.")
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = results[tech]
        print(f"  {tech}: pooled={r['pooled_effect']:.3f} 95%CI=[{r['ci_95'][0]:.3f}, {r['ci_95'][1]:.3f}] I²={r['I2']:.1f}%")
