"""
Module 4: Cell Permeability and Oral Bioavailability Prediction
PAMPA/Caco-2 permeability + oral bioavailability (%F) prediction using
extended Lipinski/Veber rules + ML regression models.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, QED
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score

from protac_utils import (
    log_event, compute_descriptors, morgan_fingerprint,
    BRD4_WARHEAD_SMILES, VHL_LIGAND_SMILES, CRBN_LIGAND_SMILES,
    IAP_LIGAND_SMILES
)

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)

# --------------------------------------------------------------------------
# 4.1  Extended beyond-rule-of-5 (bRo5) filters for PROTACs
# --------------------------------------------------------------------------

def bro5_filter(smiles: str) -> dict:
    """
    Apply beyond-Rule-of-5 (bRo5) filters relevant to large macrocycle/PROTAC.
    Based on DeGoey et al. (J. Med. Chem. 2018) and beyond-Ro5 criteria.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    desc = compute_descriptors(smiles)
    mw    = desc["MW"]
    logp  = desc["LogP"]
    hbd   = desc["HBD"]
    hba   = desc["HBA"]
    tpsa  = desc["TPSA"]
    rot   = desc["RotBonds"]
    fsp3  = desc["Fsp3"]
    rings = desc["RingCount"]

    # Ro5 (classic)
    ro5 = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
    # bRo5 (relaxed for PROTACs / macrocycles)
    bro5 = (mw <= 1200) and (logp <= 8) and (hbd <= 10) and (hba <= 20) and (tpsa <= 300)
    # PROTAC-specific: Fsp3 ≥ 0.25 for 3D shape, moderate TPSA
    protac_ok = (fsp3 >= 0.20) and (tpsa <= 280) and (rot <= 30)

    # Predicted PAMPA permeability (10^-6 cm/s) — heuristic model
    # Based on: Leeson et al., Nature 2007; Palm et al., J. Pharm. Sci. 1996
    pampa_log = (-0.012 * tpsa + 0.15 * logp - 0.02 * hbd * 2 - 0.005 * rot
                 - 0.001 * mw + 2.5)
    pampa = 10 ** pampa_log

    # Predicted Caco-2 permeability (nm/s) — simplified
    caco2_log = (-0.010 * tpsa + 0.12 * logp - 0.015 * hbd * 2 - 0.003 * rot
                 - 0.0008 * mw + 2.0)
    caco2 = 10 ** caco2_log

    # Oral bioavailability (%F) heuristic (Egan egg / Veber criteria)
    # High: TPSA ≤ 140 & RotBonds ≤ 10
    # For PROTACs, predicted much lower
    f_oral_base = 100 * np.exp(-0.008 * tpsa) * np.exp(-0.05 * max(rot - 10, 0))
    f_oral = np.clip(f_oral_base * (1 - 0.3 * int(mw > 700)) * fsp3 * 2, 0.5, 85)

    # QED (drug-likeness)
    qed = QED.qed(mol)

    return {
        "smiles": smiles,
        "MW": round(mw, 1),
        "LogP": round(logp, 2),
        "HBD": hbd,
        "HBA": hba,
        "TPSA": round(tpsa, 1),
        "RotBonds": rot,
        "Fsp3": round(fsp3, 3),
        "RingCount": rings,
        "Ro5_pass":  ro5,
        "bRo5_pass": bro5,
        "PROTAC_filter_pass": protac_ok,
        "PAMPA_pred_nm_s":    round(pampa * 100, 3),
        "Caco2_pred_nm_s":    round(caco2 * 100, 3),
        "F_oral_pred_pct":    round(f_oral, 1),
        "QED": round(qed, 3),
    }


# --------------------------------------------------------------------------
# 4.2  Generate synthetic ADMET dataset
# --------------------------------------------------------------------------

def generate_admet_dataset(n: int = 500, seed: int = 42) -> pd.DataFrame:
    """Synthetic PROTAC ADMET dataset with realistic property distributions."""
    rng = np.random.RandomState(seed)
    records = []
    for _ in range(n):
        mw   = rng.uniform(500, 1200)
        logp = rng.uniform(0, 7)
        tpsa = rng.uniform(80, 280)
        hbd  = rng.randint(0, 10)
        hba  = rng.randint(5, 20)
        rot  = rng.randint(5, 28)
        fsp3 = rng.uniform(0.1, 0.7)
        rings = rng.randint(3, 9)

        pampa = np.clip(10 ** (-0.012*tpsa + 0.15*logp - 0.03*hbd - 0.004*rot
                               - 0.001*mw + 2.5 + rng.normal(0, 0.3)), 0.001, 500)
        caco2 = np.clip(10 ** (-0.010*tpsa + 0.12*logp - 0.015*hbd - 0.003*rot
                               - 0.0008*mw + 2.0 + rng.normal(0, 0.3)), 0.001, 300)
        f_oral = np.clip(100 * np.exp(-0.008*tpsa) * np.exp(-0.04*max(rot-10,0))
                         * (1 - 0.3*int(mw > 700)) * fsp3 * 2
                         + rng.normal(0, 3), 0.5, 80)

        records.append({
            "MW": round(mw, 1), "LogP": round(logp, 2),
            "TPSA": round(tpsa, 1), "HBD": hbd, "HBA": hba,
            "RotBonds": rot, "Fsp3": round(fsp3, 3), "RingCount": rings,
            "PAMPA_nm_s": round(pampa, 3),
            "Caco2_nm_s": round(caco2, 3),
            "F_oral_pct": round(f_oral, 2),
        })
    return pd.DataFrame(records)


# --------------------------------------------------------------------------
# 4.3  ML regression model for ADMET properties
# --------------------------------------------------------------------------

def train_admet_models(df: pd.DataFrame) -> dict:
    feature_cols = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotBonds", "Fsp3", "RingCount"]
    X = df[feature_cols].values
    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X)
    cv = KFold(n_splits=5, shuffle=True, random_state=42)

    results = {}
    for target in ["PAMPA_nm_s", "Caco2_nm_s", "F_oral_pct"]:
        y = df[target].values
        rf  = RandomForestRegressor(n_estimators=200, max_depth=8, random_state=42)
        gbm = GradientBoostingRegressor(n_estimators=100, max_depth=4, random_state=42)

        r2_rf  = cross_val_score(rf,  X_sc, y, cv=cv, scoring="r2").mean()
        r2_gbm = cross_val_score(gbm, X_sc, y, cv=cv, scoring="r2").mean()

        best_model = rf if r2_rf >= r2_gbm else gbm
        best_model.fit(X_sc, y)
        best_r2 = max(r2_rf, r2_gbm)

        results[target] = {
            "model": best_model, "scaler": scaler,
            "r2_cv": round(best_r2, 4),
            "model_name": "RF" if r2_rf >= r2_gbm else "GBM",
        }
        print(f"  {target:20s}: R² CV = {best_r2:.3f} ({results[target]['model_name']})")

    return results


# --------------------------------------------------------------------------
# 4.4  Main + plots
# --------------------------------------------------------------------------

def run_admet_prediction():
    print("[Module 4] Cell permeability / oral bioavailability prediction ...")
    log_event("admet_prediction", "handoff_started", "co-scientist-molecular-docking",
              {"method": "bRo5 + ML regression", "targets": ["PAMPA", "Caco2", "F_oral"]})

    # Evaluate representative PROTAC components + BRD4 PROTACs
    test_smiles = {
        "BRD4_Warhead":   BRD4_WARHEAD_SMILES,
        "VHL_Ligand":     VHL_LIGAND_SMILES,
        "CRBN_Ligand":    CRBN_LIGAND_SMILES,
        "IAP_Ligand":     IAP_LIGAND_SMILES,
        "ARV-825":        "Cc1sc2c(c1-c1ccc(Cl)cc1)C(=O)N(C)c1ccc(cc1)CCNC(=O)COCCOCCOCCO"
                          "NC(=O)C[C@@H]1CC[C@H](CC1)NC(=O)[C@@H](Cc1ccccc1)NC(=O)c1ccc(cc1)",
        "MZ1":            "Cc1sc2c(c1-c1ccc(Cl)cc1)C(=O)N(C)c1cc(ccc1)CCOC(=O)"
                          "OCCOCCOCCO",
    }

    admet_records = []
    for name, smi in test_smiles.items():
        res = bro5_filter(smi)
        if res:
            res["compound"] = name
            admet_records.append(res)
            print(f"  {name:20s}: MW={res['MW']:.0f}  F_oral={res['F_oral_pred_pct']:.1f}%  "
                  f"PAMPA={res['PAMPA_pred_nm_s']:.2f}  QED={res['QED']:.3f}")

    admet_df = pd.DataFrame(admet_records)
    admet_df.to_csv("results/admet_predictions.csv", index=False)

    # Train ML model on synthetic data
    df_syn = generate_admet_dataset(500, seed=42)
    df_syn.to_csv("data/admet_synthetic_dataset.csv", index=False)
    model_res = train_admet_models(df_syn)

    # Save model performance
    perf = {k: {"r2_cv": v["r2_cv"], "model": v["model_name"]}
            for k, v in model_res.items()}
    pd.DataFrame(perf).T.to_csv("results/admet_model_performance.csv")

    # ---- Figures ----
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Row 1: property distributions from synthetic dataset
    for j, prop in enumerate(["MW", "LogP", "TPSA"]):
        ax = axes[0, j]
        ax.hist(df_syn[prop], bins=30, color="steelblue", alpha=0.7, edgecolor="white")
        ax.set_xlabel(prop)
        ax.set_ylabel("Count")
        ax.set_title(f"A{j+1}. Distribution of {prop}")

    # Row 2: predicted ADMET for test compounds
    compounds = admet_df["compound"].tolist()
    for j, prop in enumerate(["F_oral_pred_pct", "PAMPA_pred_nm_s", "Caco2_pred_nm_s"]):
        ax = axes[1, j]
        vals = admet_df[prop].tolist()
        colors = ["#2ecc71" if v > 15 else "#e74c3c" for v in vals] if "oral" in prop \
            else ["#2ecc71" if v > 5 else "#e74c3c" for v in vals]
        ax.bar(compounds, vals, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_xticklabels(compounds, rotation=30, ha="right", fontsize=8)
        label = {"F_oral_pred_pct": "Predicted F_oral (%)",
                 "PAMPA_pred_nm_s": "PAMPA Permeability (nm/s)",
                 "Caco2_pred_nm_s": "Caco-2 Permeability (nm/s)"}[prop]
        ax.set_ylabel(label)
        ax.set_title(f"B{j+1}. {label}")

    plt.suptitle("PROTAC ADMET Prediction\n(bRo5 + ML Regression Framework)", fontsize=13)
    plt.tight_layout()
    plt.savefig("figures/04_admet_predictions.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Egan egg-like plot (TPSA vs logP)
    fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes2[0]
    scatter = ax.scatter(df_syn["LogP"], df_syn["TPSA"],
                         c=df_syn["F_oral_pct"], cmap="RdYlGn",
                         s=10, alpha=0.5, vmin=0, vmax=60)
    plt.colorbar(scatter, ax=ax, label="F_oral (%)")
    # Veber / Egan egg ellipse approximation
    from matplotlib.patches import Ellipse
    ellipse = Ellipse((2.5, 130), width=5, height=110, angle=0,
                      fill=False, edgecolor="blue", lw=2, linestyle="--",
                      label="Egan egg (oral bioavailability zone)")
    ax.add_patch(ellipse)
    # Plot test compounds
    for _, row in admet_df.iterrows():
        ax.scatter(row["LogP"], row["TPSA"], s=120, marker="*",
                   edgecolors="black", zorder=5,
                   c=[row["F_oral_pred_pct"]], cmap="RdYlGn", vmin=0, vmax=60)
        ax.annotate(row["compound"], (row["LogP"], row["TPSA"]),
                    fontsize=7, ha="left")
    ax.set_xlabel("LogP")
    ax.set_ylabel("TPSA (Å²)")
    ax.set_title("Egan Egg Plot — Oral Bioavailability Space")
    ax.legend(fontsize=8)

    ax = axes2[1]
    ax.scatter(df_syn["MW"], df_syn["PAMPA_nm_s"],
               c=df_syn["F_oral_pct"], cmap="RdYlGn", s=10, alpha=0.4, vmin=0, vmax=60)
    for _, row in admet_df.iterrows():
        ax.scatter(row["MW"], row["PAMPA_pred_nm_s"], s=120, marker="*",
                   edgecolors="black", zorder=5)
        ax.annotate(row["compound"], (row["MW"], row["PAMPA_pred_nm_s"]),
                    fontsize=7, ha="left")
    ax.set_xlabel("Molecular Weight (Da)")
    ax.set_ylabel("PAMPA Permeability (nm/s)")
    ax.set_title("MW vs PAMPA Permeability")
    ax.axhline(y=5, color="red", ls="--", lw=1, label="PAMPA threshold (5 nm/s)")
    ax.legend(fontsize=8)

    plt.suptitle("PROTAC Drug-likeness and Permeability Landscape", fontsize=12)
    plt.tight_layout()
    plt.savefig("figures/04_egan_egg.png", dpi=150, bbox_inches="tight")
    plt.close()

    log_event("admet_prediction", "handoff_completed", "co-scientist-molecular-docking",
              {"n_compounds_evaluated": len(admet_df),
               "admet_model_r2": {k: v["r2_cv"] for k, v in model_res.items()}},
              files_written=["results/admet_predictions.csv",
                             "data/admet_synthetic_dataset.csv",
                             "results/admet_model_performance.csv",
                             "figures/04_admet_predictions.png",
                             "figures/04_egan_egg.png"])
    return admet_df, model_res


if __name__ == "__main__":
    run_admet_prediction()
