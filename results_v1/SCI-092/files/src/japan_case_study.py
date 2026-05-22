"""
Component 6: Case study — Genome-edited food acceptance in Japan.
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import json

np.random.seed(42)

def generate_japan_genome_food_data(n=600):
    data = []
    for i in range(n):
        age = int(np.clip(np.random.normal(45, 15), 18, 85))
        age_group = "18-29" if age < 30 else "30-44" if age < 45 else "45-59" if age < 60 else "60+"
        gender = np.random.choice(["male", "female"], p=[0.48, 0.52])
        education = np.random.choice(["high_school", "vocational", "bachelor", "master_doctorate"], p=[0.25, 0.15, 0.40, 0.20])
        region = np.random.choice(["kanto", "kansai", "chubu", "tohoku_hokkaido", "kyushu_okinawa", "chugoku_shikoku"], p=[0.35, 0.20, 0.15, 0.10, 0.12, 0.08])
        edu_bonus = {"high_school": 0, "vocational": 0.3, "bachelor": 0.6, "master_doctorate": 1.2}
        sci_literacy = np.clip(np.random.normal(3.5 + edu_bonus[education], 1.2), 1, 7)
        media_exp = np.clip(np.random.normal(4.5, 1.3), 1, 7)
        trust_gov = np.clip(np.random.normal(3.2, 1.3), 1, 7)
        trust_sci = np.clip(np.random.normal(4.5, 1.1), 1, 7)
        trust_corp = np.clip(np.random.normal(3.0, 1.2), 1, 7)
        trust_composite = (trust_gov + trust_sci + trust_corp) / 3
        knowledge = np.clip(0.5 * sci_literacy + np.random.normal(0, 1.0), 1, 7)
        gmo_distinction = np.clip(0.3 * knowledge + np.random.normal(2.5, 1.5), 1, 7)
        perceived_risk = np.clip(5.0 - 0.2 * trust_composite - 0.15 * knowledge + 0.1 * media_exp + np.random.normal(0, 1.0), 1, 7)
        perceived_benefit = np.clip(3.0 + 0.25 * trust_composite + 0.2 * knowledge - 0.05 * media_exp + np.random.normal(0, 1.0), 1, 7)
        naturalness_concern = np.clip(np.random.normal(5.0, 1.2), 1, 7)
        labeling_pref = np.clip(np.random.normal(5.8, 1.0), 1, 7)
        acceptance = np.clip(0.25 * trust_composite + 0.30 * perceived_benefit - 0.25 * perceived_risk - 0.15 * naturalness_concern + 0.10 * gmo_distinction + np.random.normal(2.5, 0.8), 1, 7)
        purchase_intention = np.clip(acceptance * 0.8 + 0.2 * perceived_benefit + np.random.normal(-0.5, 0.7), 1, 7)
        data.append({
            "respondent_id": f"JP{i:04d}", "age": age, "age_group": age_group, "gender": gender,
            "education": education, "region": region,
            "scientific_literacy": round(float(sci_literacy), 2),
            "media_exposure": round(float(media_exp), 2),
            "trust_government": round(float(trust_gov), 2),
            "trust_scientists": round(float(trust_sci), 2),
            "trust_corporations": round(float(trust_corp), 2),
            "trust_composite": round(float(trust_composite), 2),
            "knowledge_genome_editing": round(float(knowledge), 2),
            "gmo_distinction": round(float(gmo_distinction), 2),
            "perceived_risk": round(float(perceived_risk), 2),
            "perceived_benefit": round(float(perceived_benefit), 2),
            "naturalness_concern": round(float(naturalness_concern), 2),
            "labeling_preference": round(float(labeling_pref), 2),
            "acceptance": round(float(acceptance), 2),
            "purchase_intention": round(float(purchase_intention), 2)
        })
    return pd.DataFrame(data)

def analyze_japan_case(df):
    results = {}
    results["overall"] = {
        "n": len(df),
        "acceptance_mean": round(float(df["acceptance"].mean()), 3),
        "acceptance_std": round(float(df["acceptance"].std()), 3),
        "acceptance_median": round(float(df["acceptance"].median()), 3),
        "pct_accepting": round(float((df["acceptance"] >= 4).mean() * 100), 1),
        "pct_rejecting": round(float((df["acceptance"] < 3).mean() * 100), 1),
        "purchase_intention_mean": round(float(df["purchase_intention"].mean()), 3),
        "labeling_preference_mean": round(float(df["labeling_preference"].mean()), 3)
    }
    for demo in ["age_group", "gender", "education", "region"]:
        group_stats = df.groupby(demo).agg({"acceptance": ["mean", "std"], "perceived_risk": "mean", "perceived_benefit": "mean", "trust_composite": "mean"}).round(3)
        group_stats.columns = ["_".join(c) for c in group_stats.columns]
        results[f"by_{demo}"] = group_stats.to_dict("index")
        groups = [df[df[demo] == g]["acceptance"].values for g in df[demo].unique()]
        f_stat, p_val = stats.f_oneway(*groups)
        results[f"{demo}_anova"] = {"F": round(float(f_stat), 3), "p": round(float(p_val), 4)}

    key_vars = ["scientific_literacy", "media_exposure", "trust_composite", "knowledge_genome_editing",
                "perceived_risk", "perceived_benefit", "naturalness_concern", "acceptance", "purchase_intention"]
    corr = df[key_vars].corr().round(3)
    results["correlation_matrix"] = corr.to_dict()

    predictors = ["trust_composite", "perceived_risk", "perceived_benefit", "naturalness_concern",
                  "knowledge_genome_editing", "scientific_literacy", "media_exposure", "gmo_distinction"]
    X = df[predictors].values
    y = df["acceptance"].values
    reg = LinearRegression().fit(X, y)
    results["regression"] = {
        "r_squared": round(float(reg.score(X, y)), 4),
        "coefficients": {p: round(float(c), 4) for p, c in zip(predictors, reg.coef_)},
        "intercept": round(float(reg.intercept_), 4)
    }
    results["japan_specific"] = {
        "naturalness_concern_mean": round(float(df["naturalness_concern"].mean()), 3),
        "naturalness_acceptance_corr": round(float(df["naturalness_concern"].corr(df["acceptance"])), 4),
        "labeling_pref_mean": round(float(df["labeling_preference"].mean()), 3),
        "gmo_distinction_mean": round(float(df["gmo_distinction"].mean()), 3),
        "trust_scientists_vs_gov_gap": round(float(df["trust_scientists"].mean() - df["trust_government"].mean()), 3),
        "trust_scientists_vs_corp_gap": round(float(df["trust_scientists"].mean() - df["trust_corporations"].mean()), 3)
    }
    return results

def run_japan_case_study():
    df = generate_japan_genome_food_data()
    df.to_csv("data/japan_genome_food_data.csv", index=False)
    results = analyze_japan_case(df)
    with open("results/japan_case_study_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return df, results

if __name__ == "__main__":
    df, results = run_japan_case_study()
    o = results["overall"]
    print(f"Japan GE Food Case Study (n={o['n']}):")
    print(f"  Acceptance: M={o['acceptance_mean']:.2f} (SD={o['acceptance_std']:.2f})")
    print(f"  % Accepting: {o['pct_accepting']:.1f}%, % Rejecting: {o['pct_rejecting']:.1f}%")
    r = results["regression"]
    print(f"  Regression R²={r['r_squared']:.3f}")
