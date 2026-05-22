"""
Module 5: Degradation Activity (DC50/Dmax) SAR Analysis Automation
Builds quantitative SAR models for DC50 and Dmax from PROTAC structural features.
Implements: RF/GBM regression, SHAP-like feature contributions, activity cliff analysis.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score

from protac_utils import log_event, compute_descriptors, morgan_fingerprint

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)

# --------------------------------------------------------------------------
# 5.1  Synthetic PROTAC SAR dataset (BRD4 PROTACs)
# --------------------------------------------------------------------------

def generate_sar_dataset(n: int = 400, seed: int = 42) -> pd.DataFrame:
    """
    Synthetic SAR dataset for BRD4-targeting PROTACs.
    DC50 (nM): concentration for 50% degradation.
    Dmax (%): maximum degradation achieved.

    Literature-based correlates:
      - DC50 ↓ (more potent) with: high ternary complex stability, optimal linker length
      - Dmax ↑ with: VHL > CRBN for BRD4, Fsp3 > 0.35, low TPSA
    """
    rng = np.random.RandomState(seed)
    records = []
    e3_classes = ["VHL", "CRBN", "IAP"]

    for _ in range(n):
        # Structural features
        linker_len  = rng.randint(3, 18)   # number of atoms
        linker_type = rng.choice(["PEG", "Alkyl", "Amide", "Pip"])
        e3 = rng.choice(e3_classes, p=[0.4, 0.4, 0.2])
        mw   = 650 + linker_len * 18 + rng.normal(0, 80)
        logp = rng.uniform(1.5, 5.5)
        tpsa = rng.uniform(100, 260)
        hbd  = rng.randint(2, 9)
        hba  = rng.randint(8, 18)
        rot  = linker_len + rng.randint(3, 8)
        fsp3 = rng.uniform(0.15, 0.65)

        # Ternary complex stability proxy
        geom_score = np.exp(-((linker_len - 10) ** 2) / 18.0)
        coop = 2.5 * geom_score * (1.0 / (1 + 0.15 * rot)) + 0.5

        # DC50 model (nM) — lower is better
        # Best: linker_len 8-12, VHL, moderate LogP, low TPSA
        dc50_log = (2.2                          # baseline ~160 nM
                    - 0.8 * geom_score           # geometry contribution
                    - 0.3 * (1 if e3 == "VHL" else 0.1 if e3 == "CRBN" else -0.3)
                    - 0.05 * logp + 0.003 * tpsa
                    + 0.04 * hbd + 0.002 * mw * 0.001
                    + rng.normal(0, 0.3))
        dc50 = max(0.5, 10 ** dc50_log)

        # Dmax model (%) — higher is better
        # Best: VHL, Fsp3 > 0.4, optimal linker
        dmax = (75
                + 15 * geom_score
                + (10 if e3 == "VHL" else 5 if e3 == "CRBN" else -5)
                + 20 * (fsp3 - 0.35)
                - 0.1 * tpsa + 0.5
                + rng.normal(0, 8))
        dmax = np.clip(dmax, 5, 99)

        # Binary activity labels
        active_dc50 = dc50 < 100  # DC50 < 100 nM = active
        active_dmax = dmax > 60   # Dmax > 60% = good degrader

        records.append({
            "linker_len":   linker_len,
            "linker_type":  linker_type,
            "E3_ligase":    e3,
            "MW":           round(mw, 1),
            "LogP":         round(logp, 2),
            "TPSA":         round(tpsa, 1),
            "HBD":          hbd,
            "HBA":          hba,
            "RotBonds":     rot,
            "Fsp3":         round(fsp3, 3),
            "geom_score":   round(geom_score, 4),
            "cooperativity":round(coop, 3),
            "pDC50":        round(-np.log10(dc50 * 1e-9), 4),  # -log[M]
            "DC50_nM":      round(dc50, 2),
            "Dmax_pct":     round(dmax, 1),
            "active_dc50":  int(active_dc50),
            "active_dmax":  int(active_dmax),
        })

    return pd.DataFrame(records)


# --------------------------------------------------------------------------
# 5.2  QSAR models for DC50 and Dmax
# --------------------------------------------------------------------------

def train_sar_models(df: pd.DataFrame) -> dict:
    num_features = ["linker_len", "MW", "LogP", "TPSA", "HBD", "HBA",
                    "RotBonds", "Fsp3", "geom_score", "cooperativity"]
    cat_features = {"linker_type": ["PEG", "Alkyl", "Amide", "Pip"],
                    "E3_ligase":   ["VHL", "CRBN", "IAP"]}

    # One-hot encode categorical
    X_num = df[num_features].values.astype(np.float32)
    X_cat_parts = []
    for col, cats in cat_features.items():
        for cat in cats:
            X_cat_parts.append((df[col] == cat).astype(np.float32).values.reshape(-1, 1))
    X = np.hstack([X_num] + X_cat_parts)
    feature_names = num_features + [f"{c}_{v}" for c, vals in cat_features.items()
                                    for v in vals]

    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X)
    cv = KFold(n_splits=5, shuffle=True, random_state=42)

    results = {}
    for target, target_label in [("pDC50", "pDC50 (-log[M])"),
                                  ("Dmax_pct", "Dmax (%)")]:
        y = df[target].values
        rf  = RandomForestRegressor(n_estimators=200, max_depth=8,
                                    min_samples_leaf=3, random_state=42)
        gbm = GradientBoostingRegressor(n_estimators=150, max_depth=4,
                                        learning_rate=0.05, random_state=42)
        r2_rf  = cross_val_score(rf,  X_sc, y, cv=cv, scoring="r2").mean()
        r2_gbm = cross_val_score(gbm, X_sc, y, cv=cv, scoring="r2").mean()
        best = rf if r2_rf >= r2_gbm else gbm
        best.fit(X_sc, y)
        y_pred_cv = cross_val_predict(best, X_sc, y, cv=cv)

        results[target] = {
            "model": best, "scaler": scaler,
            "X": X_sc, "y": y, "y_pred": y_pred_cv,
            "r2_cv": round(max(r2_rf, r2_gbm), 4),
            "rmse_cv": round(np.sqrt(mean_squared_error(y, y_pred_cv)), 4),
            "feature_names": feature_names,
            "importances": best.feature_importances_,
            "model_name": "RF" if r2_rf >= r2_gbm else "GBM",
        }
        print(f"  {target:10s}: R²={results[target]['r2_cv']:.3f}  "
              f"RMSE={results[target]['rmse_cv']:.3f}  "
              f"({results[target]['model_name']})")
    return results


# --------------------------------------------------------------------------
# 5.3  Activity cliff analysis
# --------------------------------------------------------------------------

def activity_cliff_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify activity cliffs: pairs of structurally similar compounds
    with large DC50 differences (>10x).
    """
    from itertools import combinations

    # Use continuous features for similarity
    feats = ["linker_len", "MW", "LogP", "TPSA", "Fsp3", "geom_score"]
    X = StandardScaler().fit_transform(df[feats].values)

    cliff_records = []
    n = min(len(df), 100)  # limit for speed
    for i, j in combinations(range(n), 2):
        # Euclidean similarity in feature space
        dist = np.linalg.norm(X[i] - X[j])
        sim = 1.0 / (1.0 + dist)
        dc50_ratio = max(df.iloc[i]["DC50_nM"], df.iloc[j]["DC50_nM"]) / \
                     (min(df.iloc[i]["DC50_nM"], df.iloc[j]["DC50_nM"]) + 1e-6)
        dmax_diff = abs(df.iloc[i]["Dmax_pct"] - df.iloc[j]["Dmax_pct"])
        if sim > 0.6 and (dc50_ratio > 10 or dmax_diff > 30):
            cliff_records.append({
                "idx_A": i, "idx_B": j,
                "similarity": round(sim, 4),
                "DC50_A": df.iloc[i]["DC50_nM"],
                "DC50_B": df.iloc[j]["DC50_nM"],
                "DC50_ratio": round(dc50_ratio, 2),
                "Dmax_A": df.iloc[i]["Dmax_pct"],
                "Dmax_B": df.iloc[j]["Dmax_pct"],
                "Dmax_diff": round(dmax_diff, 1),
                "E3_A": df.iloc[i]["E3_ligase"],
                "E3_B": df.iloc[j]["E3_ligase"],
                "linker_len_A": df.iloc[i]["linker_len"],
                "linker_len_B": df.iloc[j]["linker_len"],
            })

    cliff_df = pd.DataFrame(cliff_records)
    return cliff_df.sort_values("DC50_ratio", ascending=False) if len(cliff_df) else cliff_df


# --------------------------------------------------------------------------
# 5.4  Main + plots
# --------------------------------------------------------------------------

def run_sar_analysis():
    print("[Module 5] DC50/Dmax SAR analysis automation ...")
    log_event("sar_analysis", "handoff_started", "co-scientist-molecular-docking",
              {"method": "RF/GBM QSAR + activity cliff", "targets": ["pDC50", "Dmax_pct"]})

    df = generate_sar_dataset(400, seed=42)
    df.to_csv("data/brd4_sar_dataset.csv", index=False)

    sar_models = train_sar_models(df)
    cliff_df   = activity_cliff_analysis(df)
    cliff_df.to_csv("results/activity_cliffs.csv", index=False)

    # ---- Figure: SAR analysis ----
    fig, axes = plt.subplots(2, 3, figsize=(17, 11))

    # A: pDC50 predicted vs observed
    ax = axes[0, 0]
    r = sar_models["pDC50"]
    ax.scatter(r["y"], r["y_pred"], alpha=0.4, s=15, c="steelblue")
    lims = [min(r["y"].min(), r["y_pred"].min()) - 0.2,
            max(r["y"].max(), r["y_pred"].max()) + 0.2]
    ax.plot(lims, lims, "r--", lw=1)
    ax.set_xlabel("Observed pDC50")
    ax.set_ylabel("Predicted pDC50")
    ax.set_title(f"A. pDC50 Model (R²={r['r2_cv']:.3f}, RMSE={r['rmse_cv']:.3f})")

    # B: Dmax predicted vs observed
    ax = axes[0, 1]
    r = sar_models["Dmax_pct"]
    ax.scatter(r["y"], r["y_pred"], alpha=0.4, s=15, c="darkorange")
    lims = [min(r["y"].min(), r["y_pred"].min()) - 2,
            max(r["y"].max(), r["y_pred"].max()) + 2]
    ax.plot(lims, lims, "r--", lw=1)
    ax.set_xlabel("Observed Dmax (%)")
    ax.set_ylabel("Predicted Dmax (%)")
    ax.set_title(f"B. Dmax Model (R²={r['r2_cv']:.3f}, RMSE={r['rmse_cv']:.3f})")

    # C: Feature importance for pDC50
    ax = axes[0, 2]
    r = sar_models["pDC50"]
    fn = r["feature_names"]
    imp = r["importances"]
    idx = np.argsort(imp)[-12:]
    ax.barh([fn[i] for i in idx], imp[idx], color="steelblue")
    ax.set_xlabel("Feature Importance")
    ax.set_title("C. pDC50 Feature Importances")

    # D: DC50 by E3 ligase (violin)
    ax = axes[1, 0]
    e3_groups = [df[df["E3_ligase"] == e]["DC50_nM"].clip(upper=500) for e in ["VHL", "CRBN", "IAP"]]
    parts = ax.violinplot(e3_groups, positions=[1, 2, 3], showmedians=True)
    for pc in parts["bodies"]:
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["VHL", "CRBN", "IAP"])
    ax.set_ylabel("DC50 (nM, clipped at 500)")
    ax.set_title("D. DC50 Distribution by E3 Ligase")

    # E: Linker length vs DC50
    ax = axes[1, 1]
    df_grp = df.groupby("linker_len")["DC50_nM"].agg(["median", "std"])
    ax.errorbar(df_grp.index, df_grp["median"], yerr=df_grp["std"].clip(upper=200),
                fmt="o-", color="steelblue", capsize=4, lw=2)
    ax.axvline(x=10, color="red", ls="--", lw=1, label="Optimal linker (10 atoms)")
    ax.set_xlabel("Linker Length (heavy atoms)")
    ax.set_ylabel("Median DC50 (nM)")
    ax.set_title("E. Linker Length vs DC50")
    ax.legend(fontsize=8)
    ax.set_ylim(0, None)

    # F: Dmax vs cooperativity
    ax = axes[1, 2]
    scatter = ax.scatter(df["cooperativity"], df["Dmax_pct"],
                         c=df["DC50_nM"].clip(upper=500), cmap="RdYlGn_r",
                         s=15, alpha=0.5, vmin=0, vmax=500)
    plt.colorbar(scatter, ax=ax, label="DC50 (nM)")
    z = np.polyfit(df["cooperativity"], df["Dmax_pct"], 1)
    xfit = np.linspace(df["cooperativity"].min(), df["cooperativity"].max(), 50)
    ax.plot(xfit, np.polyval(z, xfit), "r--", lw=2, label="linear fit")
    r_val, p_val = pearsonr(df["cooperativity"], df["Dmax_pct"])
    ax.set_xlabel("Predicted Cooperativity (α)")
    ax.set_ylabel("Dmax (%)")
    ax.set_title(f"F. Cooperativity vs Dmax (r={r_val:.3f}, p={p_val:.1e})")
    ax.legend(fontsize=8)

    plt.suptitle("BRD4 PROTAC SAR Analysis — DC50 and Dmax QSAR Models",
                 fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig("figures/05_sar_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Heatmap: E3 × linker type × Dmax
    fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))
    pivot_dc50 = df.pivot_table(values="DC50_nM", index="E3_ligase",
                                columns="linker_type", aggfunc="median")
    pivot_dmax = df.pivot_table(values="Dmax_pct", index="E3_ligase",
                                columns="linker_type", aggfunc="median")
    sns.heatmap(pivot_dc50, annot=True, fmt=".0f", cmap="RdYlGn_r",
                ax=axes2[0], cbar_kws={"label": "Median DC50 (nM)"})
    axes2[0].set_title("Median DC50 (nM): E3 × Linker Type")
    sns.heatmap(pivot_dmax, annot=True, fmt=".1f", cmap="RdYlGn",
                ax=axes2[1], cbar_kws={"label": "Median Dmax (%)"})
    axes2[1].set_title("Median Dmax (%): E3 × Linker Type")
    plt.suptitle("SAR Heatmap — PROTAC Degradation Activity", fontsize=12)
    plt.tight_layout()
    plt.savefig("figures/05_sar_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()

    model_summary = {k: {"r2_cv": v["r2_cv"], "rmse_cv": v["rmse_cv"],
                         "model": v["model_name"]}
                     for k, v in sar_models.items()}
    pd.DataFrame(model_summary).T.to_csv("results/sar_model_summary.csv")

    log_event("sar_analysis", "handoff_completed", "co-scientist-molecular-docking",
              {"n_compounds": len(df), "n_cliffs": len(cliff_df),
               "models": model_summary},
              files_written=["data/brd4_sar_dataset.csv",
                             "results/activity_cliffs.csv",
                             "results/sar_model_summary.csv",
                             "figures/05_sar_analysis.png",
                             "figures/05_sar_heatmap.png"])
    print(f"  Activity cliffs identified: {len(cliff_df)}")
    return df, sar_models, cliff_df


if __name__ == "__main__":
    run_sar_analysis()
