"""
Component 4: Quantitative evaluation of framing effects on technology acceptance.
"""

import numpy as np
import pandas as pd
from scipy import stats
import json

np.random.seed(42)

FRAMING_CONDITIONS = {
    "benefit_frame": {"effect_size": 0.35},
    "risk_frame": {"effect_size": -0.40},
    "balanced_frame": {"effect_size": 0.05},
    "scientific_frame": {"effect_size": 0.15},
    "ethical_frame": {"effect_size": -0.20}
}

def generate_framing_experiment_data(n_per_condition=120):
    technologies = ["gene_editing", "AI", "nuclear_fusion"]
    base_acceptance = {"gene_editing": 4.0, "AI": 4.5, "nuclear_fusion": 4.2}
    all_data = []
    for tech in technologies:
        for frame, config in FRAMING_CONDITIONS.items():
            for i in range(n_per_condition):
                acceptance = np.clip(base_acceptance[tech] + config["effect_size"] * 1.5 + np.random.normal(0, 1.2), 1, 7)
                risk = np.clip(7 - acceptance + np.random.normal(0, 0.8), 1, 7)
                benefit = np.clip(acceptance + np.random.normal(0, 0.7), 1, 7)
                credibility = np.clip(4.5 + (0.3 if frame == "scientific_frame" else -0.2 if frame == "ethical_frame" else 0) + np.random.normal(0, 1.0), 1, 7)
                all_data.append({
                    "participant_id": f"P{len(all_data):05d}", "technology": tech,
                    "framing_condition": frame,
                    "acceptance": round(float(acceptance), 2),
                    "perceived_risk": round(float(risk), 2),
                    "perceived_benefit": round(float(benefit), 2),
                    "information_credibility": round(float(credibility), 2),
                    "prior_knowledge": int(np.random.randint(1, 8)),
                    "age": int(np.clip(np.random.normal(40, 12), 18, 80))
                })
    return pd.DataFrame(all_data)

def analyze_framing_effects(df):
    results = {}
    for tech in df["technology"].unique():
        sub = df[df["technology"] == tech]
        groups = [sub[sub["framing_condition"] == f]["acceptance"].values for f in FRAMING_CONDITIONS.keys()]
        f_stat, p_value = stats.f_oneway(*groups)
        grand_mean = sub["acceptance"].mean()
        ss_between = sum(len(g) * (g.mean() - grand_mean)**2 for g in groups)
        ss_total = sum((sub["acceptance"] - grand_mean)**2)
        eta_squared = ss_between / ss_total
        pairwise = {}
        frames = list(FRAMING_CONDITIONS.keys())
        n_comp = len(frames) * (len(frames) - 1) // 2
        for i in range(len(frames)):
            for j in range(i + 1, len(frames)):
                g1 = sub[sub["framing_condition"] == frames[i]]["acceptance"]
                g2 = sub[sub["framing_condition"] == frames[j]]["acceptance"]
                t, p = stats.ttest_ind(g1, g2)
                d = (g1.mean() - g2.mean()) / np.sqrt(((len(g1)-1)*g1.std()**2 + (len(g2)-1)*g2.std()**2) / (len(g1) + len(g2) - 2))
                pairwise[f"{frames[i]}_vs_{frames[j]}"] = {
                    "t_statistic": round(float(t), 3), "p_value": round(float(p), 4),
                    "p_adjusted": round(float(min(p * n_comp, 1.0)), 4), "cohens_d": round(float(d), 3)
                }
        condition_means = sub.groupby("framing_condition")["acceptance"].agg(["mean", "std", "count"]).round(3).to_dict("index")
        results[tech] = {
            "anova": {"F_statistic": round(float(f_stat), 3), "p_value": round(float(p_value), 6),
                      "eta_squared": round(float(eta_squared), 4),
                      "effect_interpretation": "large" if eta_squared > 0.14 else "medium" if eta_squared > 0.06 else "small"},
            "condition_means": condition_means, "pairwise_comparisons": pairwise, "n_per_tech": len(sub)
        }
    with open("results/framing_effects_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return results

def run_framing_analysis():
    df = generate_framing_experiment_data()
    df.to_csv("data/framing_experiment_data.csv", index=False)
    results = analyze_framing_effects(df)
    return df, results

if __name__ == "__main__":
    df, results = run_framing_analysis()
    print("Framing effects analysis completed.")
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = results[tech]["anova"]
        print(f"  {tech}: F={r['F_statistic']:.2f}, p={r['p_value']:.4f}, η²={r['eta_squared']:.4f} ({r['effect_interpretation']})")
