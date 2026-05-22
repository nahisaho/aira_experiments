"""
Module 3: E3 Ligase (VHL/CRBN/IAP) Selectivity Prediction Model
Builds an ML classifier (Random Forest + SVM ensemble) to predict E3 selectivity
from PROTAC molecular fingerprints and physicochemical descriptors.
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
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict
from sklearn.metrics import (classification_report, confusion_matrix,
                             roc_auc_score, roc_curve, auc)
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA

from protac_utils import (
    log_event, compute_descriptors, morgan_fingerprint,
    BRD4_WARHEAD_SMILES, VHL_LIGAND_SMILES, CRBN_LIGAND_SMILES, IAP_LIGAND_SMILES
)

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

# --------------------------------------------------------------------------
# 3.1  Synthetic PROTAC dataset with E3-selectivity labels
# --------------------------------------------------------------------------

def generate_protac_dataset(n_samples: int = 300, seed: int = 42) -> pd.DataFrame:
    """
    Generate a synthetic PROTAC dataset with known E3-selectivity labels.
    Feature patterns are based on literature SAR for VHL/CRBN/IAP PROTACs.
    """
    rng = np.random.RandomState(seed)

    # E3-selective SMARTS patterns (simplified)
    e3_profiles = {
        "VHL":  {"MW_center": 950, "MW_std": 120, "LogP_center": 3.5, "LogP_std": 1.2,
                 "TPSA_center": 180, "TPSA_std": 30, "HBD_center": 4, "label": 0},
        "CRBN": {"MW_center": 850, "MW_std": 100, "LogP_center": 2.5, "LogP_std": 1.0,
                 "TPSA_center": 150, "TPSA_std": 25, "HBD_center": 3, "label": 1},
        "IAP":  {"MW_center": 800, "MW_std": 110, "LogP_center": 4.0, "LogP_std": 1.5,
                 "TPSA_center": 130, "TPSA_std": 20, "HBD_center": 2, "label": 2},
    }

    records = []
    for e3, profile in e3_profiles.items():
        n = n_samples // 3
        for _ in range(n):
            mw   = rng.normal(profile["MW_center"],   profile["MW_std"])
            logp = rng.normal(profile["LogP_center"],  profile["LogP_std"])
            tpsa = rng.normal(profile["TPSA_center"],  profile["TPSA_std"])
            hbd  = int(np.clip(rng.normal(profile["HBD_center"], 1.5), 0, 10))
            hba  = int(np.clip(rng.normal(hbd + 5, 2), 0, 20))
            rot  = int(np.clip(rng.normal(12, 3), 3, 25))
            mw   = max(mw, 400)

            # Create a pseudo-SMILES-fingerprint as Gaussian noise with E3 signature
            fp_base = np.zeros(128)
            if profile["label"] == 0:  # VHL: pyridinyl, hydroxy features
                fp_base[0:30]  = rng.normal(0.6, 0.2, 30)
            elif profile["label"] == 1:  # CRBN: glutarimide, phthalimide
                fp_base[30:70] = rng.normal(0.6, 0.2, 40)
            else:  # IAP: bivalent SMAC mimetic
                fp_base[70:110] = rng.normal(0.6, 0.2, 40)
            fp_base += rng.normal(0, 0.1, 128)
            fp_base = np.clip(fp_base, 0, 1)

            records.append({
                "E3_ligase": e3,
                "label": profile["label"],
                "MW":     round(mw, 1),
                "LogP":   round(logp, 2),
                "TPSA":   round(np.clip(tpsa, 50, 300), 1),
                "HBD":    hbd,
                "HBA":    hba,
                "RotBonds": rot,
                **{f"fp_{i}": round(fp_base[i], 4) for i in range(128)},
            })

    df = pd.DataFrame(records).sample(frac=1, random_state=seed).reset_index(drop=True)
    return df


# --------------------------------------------------------------------------
# 3.2  Build and evaluate selectivity model
# --------------------------------------------------------------------------

def build_selectivity_model(df: pd.DataFrame):
    fp_cols = [c for c in df.columns if c.startswith("fp_")]
    phys_cols = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotBonds"]
    X = df[phys_cols + fp_cols].values.astype(np.float32)
    y = df["label"].values

    le = LabelEncoder()
    y = le.fit_transform(df["E3_ligase"])
    class_names = le.classes_

    # Ensemble: RF + GBM
    rf = RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=-1)
    gbm = GradientBoostingClassifier(n_estimators=100, max_depth=4, random_state=42)
    scaler = StandardScaler()

    X_sc = scaler.fit_transform(X)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # CV scores
    rf_scores  = cross_val_score(rf,  X_sc, y, cv=cv, scoring="f1_macro")
    gbm_scores = cross_val_score(gbm, X_sc, y, cv=cv, scoring="f1_macro")

    print(f"  RF  CV F1-macro: {rf_scores.mean():.3f} ± {rf_scores.std():.3f}")
    print(f"  GBM CV F1-macro: {gbm_scores.mean():.3f} ± {gbm_scores.std():.3f}")

    # Fit on all data for analysis
    rf.fit(X_sc, y)
    gbm.fit(X_sc, y)
    y_pred_rf  = cross_val_predict(rf,  X_sc, y, cv=cv)
    y_pred_gbm = cross_val_predict(gbm, X_sc, y, cv=cv)
    y_prob_rf  = cross_val_predict(rf,  X_sc, y, cv=cv, method="predict_proba")
    y_prob_gbm = cross_val_predict(gbm, X_sc, y, cv=cv, method="predict_proba")

    # Ensemble probability
    y_prob = (y_prob_rf + y_prob_gbm) / 2
    y_pred_ens = y_prob.argmax(axis=1)

    report = classification_report(y, y_pred_ens, target_names=class_names, output_dict=True)
    cm = confusion_matrix(y, y_pred_ens)

    # Feature importances
    feat_names = phys_cols + fp_cols
    importances = rf.feature_importances_
    feat_df = pd.DataFrame({"feature": feat_names, "importance": importances})
    feat_df = feat_df.sort_values("importance", ascending=False).head(20)

    # ROC AUC (one-vs-rest)
    roc_auc = roc_auc_score(y, y_prob, multi_class="ovr", average="macro")
    print(f"  Ensemble ROC-AUC (macro): {roc_auc:.3f}")

    return {
        "rf": rf, "gbm": gbm, "scaler": scaler,
        "class_names": class_names, "le": le,
        "y": y, "y_pred": y_pred_ens, "y_prob": y_prob,
        "report": report, "cm": cm, "feat_df": feat_df,
        "roc_auc": roc_auc,
        "rf_cv_f1": rf_scores.mean(), "gbm_cv_f1": gbm_scores.mean(),
    }


# --------------------------------------------------------------------------
# 3.3  Visualizations
# --------------------------------------------------------------------------

def plot_e3_selectivity(df: pd.DataFrame, model_results: dict):
    class_names = model_results["class_names"]
    cm = model_results["cm"]
    report = model_results["report"]
    feat_df = model_results["feat_df"]
    y = model_results["y"]
    y_prob = model_results["y_prob"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # A: Confusion matrix
    ax = axes[0, 0]
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=class_names, yticklabels=class_names, ax=ax)
    ax.set_xlabel("Predicted E3 Ligase")
    ax.set_ylabel("True E3 Ligase")
    ax.set_title("A. Confusion Matrix (5-fold CV)")

    # B: Feature importance (top 20)
    ax = axes[0, 1]
    # Show only physico-chemical features for interpretability
    phys_feats = feat_df[~feat_df["feature"].str.startswith("fp_")]
    fp_sum = feat_df[feat_df["feature"].str.startswith("fp_")]["importance"].sum()
    disp_df = pd.concat([
        phys_feats,
        pd.DataFrame([{"feature": "Fingerprint (sum)", "importance": fp_sum}])
    ]).sort_values("importance", ascending=True)
    ax.barh(disp_df["feature"], disp_df["importance"], color="steelblue")
    ax.set_xlabel("Feature Importance")
    ax.set_title("B. Top Feature Importances (Random Forest)")

    # C: ROC curves
    ax = axes[1, 0]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for i, (cls, color) in enumerate(zip(class_names, colors)):
        fpr, tpr, _ = roc_curve((y == i).astype(int), y_prob[:, i])
        roc_auc_cls = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=color, lw=2,
                label=f"{cls} (AUC={roc_auc_cls:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("C. ROC Curves (One-vs-Rest)")
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)

    # D: PCA of chemical space colored by E3 selectivity
    ax = axes[1, 1]
    phys_cols = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotBonds"]
    X_phys = df[phys_cols].values
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(StandardScaler().fit_transform(X_phys))
    for i, (cls, color) in enumerate(zip(class_names, colors)):
        mask = df["label"].values == i
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1], c=color, alpha=0.5,
                   s=20, label=cls)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("D. Chemical Space PCA by E3 Selectivity")
    ax.legend(fontsize=9)

    plt.suptitle("E3 Ligase (VHL/CRBN/IAP) Selectivity Prediction Model",
                 fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig("figures/03_e3_selectivity_model.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Per-class metrics bar chart
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    metrics = ["precision", "recall", "f1-score"]
    x = np.arange(len(class_names))
    w = 0.25
    for j, metric in enumerate(metrics):
        vals = [report[cls][metric] for cls in class_names]
        ax2.bar(x + j * w, vals, w, label=metric.capitalize())
    ax2.set_xticks(x + w)
    ax2.set_xticklabels(class_names)
    ax2.set_ylim(0, 1.1)
    ax2.set_ylabel("Score")
    ax2.set_title("E3 Ligase Selectivity — Per-Class Metrics")
    ax2.legend()
    plt.tight_layout()
    plt.savefig("figures/03_e3_per_class_metrics.png", dpi=150, bbox_inches="tight")
    plt.close()


# --------------------------------------------------------------------------
# 3.4  Main
# --------------------------------------------------------------------------

def run_e3_selectivity():
    print("[Module 3] E3 ligase selectivity prediction model ...")
    log_event("e3_selectivity", "handoff_started", "co-scientist-molecular-docking",
              {"method": "RF+GBM ensemble", "classes": ["VHL", "CRBN", "IAP"]})

    df = generate_protac_dataset(n_samples=300, seed=42)
    df.to_csv("data/e3_selectivity_dataset.csv", index=False)

    model_results = build_selectivity_model(df)
    plot_e3_selectivity(df, model_results)

    # Save classification report
    rpt_df = pd.DataFrame(model_results["report"]).T
    rpt_df.to_csv("results/e3_selectivity_report.csv")

    summary = {
        "roc_auc_macro": round(model_results["roc_auc"], 4),
        "rf_cv_f1_mean": round(model_results["rf_cv_f1"], 4),
        "gbm_cv_f1_mean": round(model_results["gbm_cv_f1"], 4),
        "classes": list(model_results["class_names"]),
    }

    log_event("e3_selectivity", "handoff_completed", "co-scientist-molecular-docking",
              summary,
              files_written=["data/e3_selectivity_dataset.csv",
                             "results/e3_selectivity_report.csv",
                             "figures/03_e3_selectivity_model.png",
                             "figures/03_e3_per_class_metrics.png"])

    print(f"  ROC-AUC (macro): {model_results['roc_auc']:.3f}")
    return model_results


if __name__ == "__main__":
    run_e3_selectivity()
