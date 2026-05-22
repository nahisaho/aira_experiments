"""
Module 4: Combinatorial Perturbation Interaction (Epistasis) Detection
=======================================================================
- Expected vs. observed combinatorial effects
- Epistasis scoring (synergy / buffering / suppression)
- Interaction significance testing
"""

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os
import json
import warnings
warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)


def compute_perturbation_signatures(adata, control_key="non-targeting"):
    """
    Compute mean expression shift for each single and combo perturbation
    relative to controls. Returns dict of perturbation → gene signature vector.
    """
    if "log_normalized" in adata.layers:
        X = adata.layers["log_normalized"]
    else:
        X = adata.X

    if hasattr(X, "toarray"):
        X = X.toarray()

    ctrl_mask = adata.obs["perturbation"] == control_key
    ctrl_mean = X[ctrl_mask].mean(axis=0)

    signatures = {}
    for pert in adata.obs["perturbation"].unique():
        if pert == control_key:
            continue
        mask = adata.obs["perturbation"] == pert
        if mask.sum() < 5:
            continue
        signatures[pert] = X[mask].mean(axis=0) - ctrl_mean

    return signatures, ctrl_mean


def detect_epistasis(adata, signatures, alpha=0.05):
    """
    For each combinatorial perturbation A|B:
      - Expected effect = effect(A) + effect(B)  (additive model)
      - Observed effect = measured combo effect
      - Epistasis = observed - expected

    Tests:
      - Synergy: |observed| > |expected| (stronger than additive)
      - Buffering: |observed| < |expected| (weaker than additive)
      - Suppression: sign(observed) ≠ sign(expected)
    """
    combo_perts = [p for p in signatures.keys() if "|" in p]
    single_perts = [p for p in signatures.keys() if "|" not in p]

    results = []

    for combo in combo_perts:
        guides = combo.split("|")
        if len(guides) != 2:
            continue

        g1, g2 = guides[0], guides[1]
        if g1 not in signatures or g2 not in signatures:
            continue

        observed = signatures[combo]
        expected = signatures[g1] + signatures[g2]
        epistasis = observed - expected

        # Global metrics
        obs_magnitude = np.linalg.norm(observed)
        exp_magnitude = np.linalg.norm(expected)
        epi_magnitude = np.linalg.norm(epistasis)

        # Cosine similarity between observed and expected
        cos_sim = np.dot(observed, expected) / (obs_magnitude * exp_magnitude + 1e-10)

        # Classification
        if obs_magnitude > exp_magnitude * 1.2:
            interaction_type = "synergy"
        elif obs_magnitude < exp_magnitude * 0.8:
            interaction_type = "buffering"
        elif cos_sim < 0:
            interaction_type = "suppression"
        else:
            interaction_type = "additive"

        # Per-gene epistasis scores (top genes)
        top_epi_genes = np.argsort(np.abs(epistasis))[-10:]
        gene_details = {adata.var_names[g]: float(epistasis[g]) for g in top_epi_genes}

        results.append({
            "combo": combo,
            "guide_1": g1,
            "guide_2": g2,
            "observed_magnitude": float(obs_magnitude),
            "expected_magnitude": float(exp_magnitude),
            "epistasis_magnitude": float(epi_magnitude),
            "cosine_similarity": float(cos_sim),
            "interaction_type": interaction_type,
            "top_epistasis_genes": gene_details,
        })

    return pd.DataFrame(results)


def permutation_test_epistasis(adata, combo, g1, g2, control_key="non-targeting",
                                n_perms=1000):
    """
    Permutation test for epistasis significance.
    Null: random reassignment of cells to perturbation groups.
    """
    if "log_normalized" in adata.layers:
        X = adata.layers["log_normalized"]
    else:
        X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    ctrl_mask = adata.obs["perturbation"] == control_key
    combo_mask = adata.obs["perturbation"] == combo
    g1_mask = adata.obs["perturbation"] == g1
    g2_mask = adata.obs["perturbation"] == g2

    ctrl_mean = X[ctrl_mask].mean(axis=0)
    obs_combo = X[combo_mask].mean(axis=0) - ctrl_mean
    obs_g1 = X[g1_mask].mean(axis=0) - ctrl_mean
    obs_g2 = X[g2_mask].mean(axis=0) - ctrl_mean

    observed_epi = np.linalg.norm(obs_combo - (obs_g1 + obs_g2))

    # Pool all relevant cells
    all_mask = ctrl_mask | combo_mask | g1_mask | g2_mask
    pooled = X[all_mask]
    n_ctrl = ctrl_mask.sum()
    n_combo_cells = combo_mask.sum()
    n_g1_cells = g1_mask.sum()
    n_g2_cells = g2_mask.sum()

    null_epi = []
    for _ in range(n_perms):
        perm_idx = np.random.permutation(len(pooled))
        perm_ctrl = pooled[perm_idx[:n_ctrl]].mean(axis=0)
        perm_combo = pooled[perm_idx[n_ctrl:n_ctrl + n_combo_cells]].mean(axis=0) - perm_ctrl
        perm_g1 = pooled[perm_idx[n_ctrl + n_combo_cells:n_ctrl + n_combo_cells + n_g1_cells]].mean(axis=0) - perm_ctrl
        perm_g2 = pooled[perm_idx[-n_g2_cells:]].mean(axis=0) - perm_ctrl
        null_epi.append(np.linalg.norm(perm_combo - (perm_g1 + perm_g2)))

    p_value = np.mean(np.array(null_epi) >= observed_epi)
    return p_value


def plot_epistasis(epi_df, out_dir="figures"):
    """Visualize epistasis results."""
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # 1. Observed vs Expected magnitude
    ax = axes[0]
    colors = {"synergy": "red", "buffering": "blue", "suppression": "orange", "additive": "grey"}
    for itype, color in colors.items():
        mask = epi_df["interaction_type"] == itype
        if mask.sum() > 0:
            ax.scatter(epi_df.loc[mask, "expected_magnitude"],
                      epi_df.loc[mask, "observed_magnitude"],
                      c=color, label=itype, s=80, alpha=0.7, edgecolors="white")
    lims = [0, max(epi_df["expected_magnitude"].max(), epi_df["observed_magnitude"].max()) * 1.1]
    ax.plot(lims, lims, "k--", alpha=0.3, label="Additive expectation")
    ax.set_xlabel("Expected Effect Magnitude")
    ax.set_ylabel("Observed Effect Magnitude")
    ax.set_title("Epistasis: Observed vs Expected")
    ax.legend(fontsize=8)

    # 2. Interaction type distribution
    ax = axes[1]
    type_counts = epi_df["interaction_type"].value_counts()
    ax.bar(type_counts.index, type_counts.values,
           color=[colors.get(t, "grey") for t in type_counts.index],
           edgecolor="white")
    ax.set_ylabel("N Combinations")
    ax.set_title("Interaction Types")

    # 3. Epistasis magnitude ranking
    ax = axes[2]
    epi_sorted = epi_df.sort_values("epistasis_magnitude", ascending=True)
    bar_colors = [colors.get(t, "grey") for t in epi_sorted["interaction_type"]]
    ax.barh(range(len(epi_sorted)), epi_sorted["epistasis_magnitude"],
            color=bar_colors, edgecolor="white")
    ax.set_yticks(range(len(epi_sorted)))
    ax.set_yticklabels(epi_sorted["combo"], fontsize=7)
    ax.set_xlabel("Epistasis Magnitude")
    ax.set_title("Epistasis Score Ranking")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "04_epistasis.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, "04_epistasis.svg"), bbox_inches="tight")
    plt.close()
    print(f"✓ Epistasis plots saved")


def run_epistasis_pipeline(adata_path="data/perturbseq_processed.h5ad"):
    """Main epistasis detection pipeline."""
    adata = ad.read_h5ad(adata_path)
    print(f"Loaded: {adata.shape}")

    # Compute signatures
    print("Computing perturbation signatures...")
    signatures, ctrl_mean = compute_perturbation_signatures(adata)
    print(f"  Single perturbations: {sum(1 for k in signatures if '|' not in k)}")
    print(f"  Combo perturbations: {sum(1 for k in signatures if '|' in k)}")

    # Detect epistasis
    print("Detecting epistasis...")
    epi_df = detect_epistasis(adata, signatures)

    if len(epi_df) == 0:
        print("⚠ No combinatorial perturbations found for epistasis analysis")
        summary = {"n_combos": 0, "message": "No combinatorial perturbations detected"}
        with open("results/04_epistasis_summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        return None, summary

    # Permutation tests for top interactions
    print("Running permutation tests...")
    p_values = []
    for _, row in epi_df.iterrows():
        try:
            p = permutation_test_epistasis(adata, row["combo"], row["guide_1"],
                                           row["guide_2"], n_perms=500)
        except Exception:
            p = 1.0
        p_values.append(p)

    epi_df["p_value"] = p_values
    if len(p_values) > 1:
        _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")
        epi_df["p_adj"] = p_adj
    else:
        epi_df["p_adj"] = p_values
    epi_df["significant"] = epi_df["p_adj"] < 0.05

    # Plot
    plot_epistasis(epi_df)

    # Save
    os.makedirs("results", exist_ok=True)
    epi_save = epi_df.drop(columns=["top_epistasis_genes"])
    epi_save.to_csv("results/04_epistasis_results.csv", index=False)

    summary = {
        "n_combos_tested": len(epi_df),
        "n_significant": int(epi_df["significant"].sum()),
        "interaction_types": epi_df["interaction_type"].value_counts().to_dict(),
        "mean_epistasis_magnitude": float(epi_df["epistasis_magnitude"].mean()),
        "top_interaction": epi_df.sort_values("epistasis_magnitude", ascending=False).iloc[0][
            ["combo", "interaction_type", "epistasis_magnitude"]].to_dict() if len(epi_df) > 0 else {},
    }
    with open("results/04_epistasis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✓ Epistasis analysis complete:")
    print(f"  Combos tested: {summary['n_combos_tested']}")
    print(f"  Significant: {summary['n_significant']}")
    print(f"  Types: {summary['interaction_types']}")

    return epi_df, summary


if __name__ == "__main__":
    run_epistasis_pipeline()
