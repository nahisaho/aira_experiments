"""
Module 2: Gene Program Variation Detection
===========================================
- Differential expression per perturbation (Wilcoxon rank-sum)
- Co-expression module detection (Hotspot / NMF)
- Gene program activity scoring
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import os
import json
import warnings
warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)


def preprocess(adata):
    """Standard preprocessing: normalize, log-transform, HVG, scale."""
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log_normalized"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor="seurat_v3",
                                 layer="raw_counts")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, random_state=SEED)
    sc.pp.neighbors(adata, n_pcs=20, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    return adata


def differential_expression_per_perturbation(adata, control_key="non-targeting",
                                              min_cells=20, alpha=0.05):
    """
    Wilcoxon rank-sum test: each perturbation vs. non-targeting controls.
    Returns DataFrame with log2FC, p-value, adjusted p-value.
    """
    control_mask = adata.obs["perturbation"] == control_key
    perturbations = [p for p in adata.obs["perturbation"].unique()
                     if p != control_key and "|" not in p]

    # Use log-normalized data
    if "log_normalized" in adata.layers:
        X = adata.layers["log_normalized"]
    else:
        X = adata.X

    if hasattr(X, "toarray"):
        X = X.toarray()

    ctrl_expr = X[control_mask]
    results = []

    for pert in perturbations:
        pert_mask = adata.obs["perturbation"] == pert
        n_pert = pert_mask.sum()
        if n_pert < min_cells:
            continue

        pert_expr = X[pert_mask]

        for g_idx in range(adata.n_vars):
            ctrl_vals = ctrl_expr[:, g_idx]
            pert_vals = pert_expr[:, g_idx]

            mean_ctrl = np.mean(ctrl_vals)
            mean_pert = np.mean(pert_vals)
            log2fc = np.log2(mean_pert + 1) - np.log2(mean_ctrl + 1)

            stat, pval = ranksums(pert_vals, ctrl_vals)

            results.append({
                "perturbation": pert,
                "gene": adata.var_names[g_idx],
                "log2FC": log2fc,
                "mean_ctrl": mean_ctrl,
                "mean_pert": mean_pert,
                "p_value": pval,
                "n_cells": int(n_pert),
            })

    df = pd.DataFrame(results)
    if len(df) == 0:
        return df

    # Multiple testing correction (BH-FDR per perturbation)
    adj_pvals = []
    for pert in df["perturbation"].unique():
        mask = df["perturbation"] == pert
        _, pvals_adj, _, _ = multipletests(df.loc[mask, "p_value"], method="fdr_bh")
        adj_pvals.extend(pvals_adj)
    df["p_adj"] = adj_pvals
    df["significant"] = df["p_adj"] < alpha

    return df


def nmf_coexpression_modules(adata, n_modules=8, n_top_genes=500):
    """
    NMF-based co-expression module detection.
    Returns module membership and gene loading matrix.
    """
    # Use HVGs for module detection
    hvg_mask = adata.var["highly_variable"].values if "highly_variable" in adata.var else \
        np.ones(adata.n_vars, dtype=bool)

    if "log_normalized" in adata.layers:
        X = adata.layers["log_normalized"][:, hvg_mask]
    else:
        X = adata.X[:, hvg_mask]

    if hasattr(X, "toarray"):
        X = X.toarray()

    # Ensure non-negative
    X = np.maximum(X, 0)

    # Fit NMF
    model = NMF(n_components=n_modules, init="nndsvda", random_state=SEED,
                max_iter=500, l1_ratio=0.5)
    W = model.fit_transform(X)  # cells × modules
    H = model.components_        # modules × genes

    gene_names = adata.var_names[hvg_mask]

    # Module membership: assign each gene to its strongest module
    module_assignments = {}
    for mod_idx in range(n_modules):
        top_genes_idx = np.argsort(H[mod_idx])[-n_top_genes // n_modules:]
        module_assignments[f"Module_{mod_idx}"] = gene_names[top_genes_idx].tolist()

    # Store cell-level module scores
    for mod_idx in range(n_modules):
        adata.obs[f"NMF_module_{mod_idx}"] = W[:, mod_idx]

    return adata, module_assignments, H, gene_names


def plot_de_results(de_df, module_assignments, adata, out_dir="figures"):
    """Visualization of DE and co-expression results."""
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # 1. Volcano plot (top perturbation)
    ax = axes[0, 0]
    top_pert = de_df.groupby("perturbation")["significant"].sum().idxmax()
    sub = de_df[de_df["perturbation"] == top_pert].copy()
    sub["-log10(p_adj)"] = -np.log10(sub["p_adj"].clip(1e-300))

    sig_mask = sub["significant"]
    ax.scatter(sub.loc[~sig_mask, "log2FC"], sub.loc[~sig_mask, "-log10(p_adj)"],
               c="grey", alpha=0.3, s=8, label="NS")
    ax.scatter(sub.loc[sig_mask, "log2FC"], sub.loc[sig_mask, "-log10(p_adj)"],
               c="crimson", alpha=0.6, s=12, label="Significant")
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(f"Volcano Plot: {top_pert}")
    ax.legend(fontsize=8)

    # 2. DE gene counts per perturbation
    ax = axes[0, 1]
    de_counts = de_df[de_df["significant"]].groupby("perturbation").size().sort_values(ascending=True)
    if len(de_counts) > 0:
        de_counts.plot.barh(ax=ax, color="steelblue", edgecolor="white")
    ax.set_xlabel("N Significant DEGs (FDR < 0.05)")
    ax.set_title("DEGs per Perturbation")

    # 3. Module sizes
    ax = axes[1, 0]
    mod_sizes = {k: len(v) for k, v in module_assignments.items()}
    ax.bar(mod_sizes.keys(), mod_sizes.values(), color="coral", edgecolor="white")
    ax.set_ylabel("N Genes")
    ax.set_title("Co-expression Module Sizes (NMF)")
    ax.tick_params(axis="x", rotation=45)

    # 4. UMAP colored by top module
    ax = axes[1, 1]
    if "X_umap" in adata.obsm:
        umap = adata.obsm["X_umap"]
        scores = adata.obs["NMF_module_0"].values
        sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=scores, cmap="viridis",
                             s=3, alpha=0.5)
        plt.colorbar(sc_plot, ax=ax, label="Module 0 Score")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    ax.set_title("UMAP: NMF Module 0 Activity")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "02_de_coexpression.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, "02_de_coexpression.svg"), bbox_inches="tight")
    plt.close()
    print(f"✓ DE/co-expression plots saved")


def run_de_pipeline(adata_path="data/perturbseq_qc_filtered.h5ad"):
    """Main DE + co-expression pipeline."""
    adata = ad.read_h5ad(adata_path)
    print(f"Loaded: {adata.shape}")

    # Preprocess
    adata = preprocess(adata)

    # Differential expression
    print("Running differential expression analysis...")
    de_df = differential_expression_per_perturbation(adata)

    # Co-expression modules
    print("Detecting co-expression modules (NMF)...")
    adata, modules, H, gene_names = nmf_coexpression_modules(adata)

    # Plots
    plot_de_results(de_df, modules, adata)

    # Save results
    os.makedirs("results", exist_ok=True)
    de_df.to_csv("results/02_de_results.csv", index=False)

    de_summary = {
        "total_tests": len(de_df),
        "significant_genes": int(de_df["significant"].sum()),
        "perturbations_tested": int(de_df["perturbation"].nunique()),
        "de_genes_per_perturbation": de_df[de_df["significant"]].groupby("perturbation").size().to_dict(),
        "n_modules": len(modules),
        "module_sizes": {k: len(v) for k, v in modules.items()},
    }
    with open("results/02_de_summary.json", "w") as f:
        json.dump(de_summary, f, indent=2)

    with open("results/02_modules.json", "w") as f:
        json.dump(modules, f, indent=2)

    adata.write_h5ad("data/perturbseq_processed.h5ad")

    print(f"✓ DE analysis complete:")
    print(f"  Tests: {de_summary['total_tests']}, Significant: {de_summary['significant_genes']}")
    print(f"  Modules detected: {de_summary['n_modules']}")

    return adata, de_df, modules


if __name__ == "__main__":
    run_de_pipeline()
