"""
Module 1: Perturbation Assignment QC & Guide Detection
=======================================================
- Guide UMI thresholding (mixture model)
- Multiplet detection
- Assignment confidence scoring
- QC visualization
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
import os
import json

SEED = 42
np.random.seed(SEED)


def guide_umi_thresholding(adata, min_umi=5, gmm_components=2):
    """
    Assign guides to cells using Gaussian Mixture Model on log(UMI+1).
    Returns per-cell guide assignment with confidence scores.
    """
    guide_umi = adata.obsm["guide_umi"]
    guide_names = adata.uns["guide_names"]
    n_cells, n_guides = guide_umi.shape

    # Fit GMM per guide to separate signal from noise
    thresholds = {}
    for g_idx in range(n_guides):
        counts = guide_umi[:, g_idx]
        log_counts = np.log1p(counts).reshape(-1, 1)

        if counts.max() < min_umi:
            thresholds[guide_names[g_idx]] = min_umi
            continue

        gmm = GaussianMixture(n_components=gmm_components, random_state=SEED)
        gmm.fit(log_counts)

        # Threshold = midpoint between component means
        means = sorted(gmm.means_.flatten())
        if len(means) == 2:
            threshold = np.exp((means[0] + means[1]) / 2) - 1
        else:
            threshold = min_umi
        thresholds[guide_names[g_idx]] = max(threshold, min_umi)

    # Assign guides
    assigned_guides = []
    confidence_scores = []
    n_detected = []

    for cell_idx in range(n_cells):
        detected = []
        max_score = 0
        for g_idx in range(n_guides):
            umi = guide_umi[cell_idx, g_idx]
            thr = thresholds[guide_names[g_idx]]
            if umi >= thr:
                detected.append(guide_names[g_idx])
                score = min(umi / thr, 10) / 10  # normalized confidence
                max_score = max(max_score, score)

        if len(detected) == 0:
            assigned_guides.append("unassigned")
            confidence_scores.append(0.0)
        elif len(detected) == 1:
            assigned_guides.append(detected[0])
            confidence_scores.append(max_score)
        else:
            assigned_guides.append("|".join(sorted(detected)))
            confidence_scores.append(max_score * 0.8)  # penalize multiplets

        n_detected.append(len(detected))

    adata.obs["guide_assignment"] = assigned_guides
    adata.obs["guide_confidence"] = confidence_scores
    adata.obs["n_guides_detected"] = n_detected

    return adata, thresholds


def cell_level_qc(adata):
    """Standard scRNA-seq QC metrics."""
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Flag low-quality cells
    min_genes = 200
    max_genes = np.percentile(adata.obs["n_genes_by_counts"], 99)
    min_counts = 500

    adata.obs["pass_qc"] = (
        (adata.obs["n_genes_by_counts"] >= min_genes) &
        (adata.obs["n_genes_by_counts"] <= max_genes) &
        (adata.obs["total_counts"] >= min_counts)
    )
    return adata


def plot_qc(adata, thresholds, out_dir="figures"):
    """Generate QC plots."""
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 1. Guide UMI distribution
    guide_umi = adata.obsm["guide_umi"]
    ax = axes[0, 0]
    ax.hist(np.log1p(guide_umi.sum(axis=1)), bins=50, color="steelblue", edgecolor="white")
    ax.set_xlabel("log(Total Guide UMI + 1)")
    ax.set_ylabel("Cells")
    ax.set_title("Guide UMI Distribution")

    # 2. Guides per cell
    ax = axes[0, 1]
    n_det = adata.obs["n_guides_detected"]
    ax.bar(*np.unique(n_det, return_counts=True), color="coral", edgecolor="white")
    ax.set_xlabel("N Guides Detected")
    ax.set_ylabel("Cells")
    ax.set_title("Guides per Cell")

    # 3. Confidence score distribution
    ax = axes[0, 2]
    ax.hist(adata.obs["guide_confidence"], bins=50, color="seagreen", edgecolor="white")
    ax.set_xlabel("Assignment Confidence")
    ax.set_ylabel("Cells")
    ax.set_title("Guide Assignment Confidence")

    # 4. Genes per cell
    ax = axes[1, 0]
    ax.hist(adata.obs["n_genes_by_counts"], bins=50, color="mediumpurple", edgecolor="white")
    ax.set_xlabel("N Genes")
    ax.set_ylabel("Cells")
    ax.set_title("Genes per Cell")

    # 5. UMI counts per cell
    ax = axes[1, 1]
    ax.hist(np.log1p(adata.obs["total_counts"]), bins=50, color="goldenrod", edgecolor="white")
    ax.set_xlabel("log(Total UMI + 1)")
    ax.set_ylabel("Cells")
    ax.set_title("UMI Counts per Cell")

    # 6. Assignment pie chart
    ax = axes[1, 2]
    assignment_cats = pd.Series(adata.obs["guide_assignment"].values)
    cats = {
        "Control (NT)": (assignment_cats == "non-targeting").sum(),
        "Single Guide": ((assignment_cats != "non-targeting") &
                         (assignment_cats != "unassigned") &
                         (~assignment_cats.str.contains(r"\|", regex=True))).sum(),
        "Combo": assignment_cats.str.contains(r"\|", regex=True).sum(),
        "Unassigned": (assignment_cats == "unassigned").sum(),
    }
    ax.pie(cats.values(), labels=cats.keys(), autopct="%1.1f%%",
           colors=["#4CAF50", "#2196F3", "#FF9800", "#9E9E9E"])
    ax.set_title("Cell Assignment Categories")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "01_qc_guide_detection.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, "01_qc_guide_detection.svg"), bbox_inches="tight")
    plt.close()
    print(f"✓ QC plots saved to {out_dir}/01_qc_guide_detection.png")


def run_qc_pipeline(adata_path="data/perturbseq_simulated.h5ad"):
    """Main QC pipeline."""
    adata = ad.read_h5ad(adata_path)
    print(f"Loaded: {adata.shape}")

    # Step 1: Guide detection
    adata, thresholds = guide_umi_thresholding(adata)

    # Step 2: Cell QC
    adata = cell_level_qc(adata)

    # Step 3: Plots
    plot_qc(adata, thresholds)

    # Step 4: Summary stats
    stats = {
        "total_cells": adata.n_obs,
        "pass_qc": int(adata.obs["pass_qc"].sum()),
        "fail_qc": int((~adata.obs["pass_qc"]).sum()),
        "unassigned": int((adata.obs["guide_assignment"] == "unassigned").sum()),
        "single_guide": int(adata.obs["n_guides_detected"].eq(1).sum()),
        "multi_guide": int(adata.obs["n_guides_detected"].gt(1).sum()),
        "mean_confidence": float(adata.obs["guide_confidence"].mean()),
        "median_guide_umi": float(np.median(adata.obsm["guide_umi"].sum(axis=1))),
        "thresholds": {k: float(v) for k, v in thresholds.items()},
    }

    os.makedirs("results", exist_ok=True)
    with open("results/01_qc_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Save filtered data
    adata_filtered = adata[adata.obs["pass_qc"] & (adata.obs["guide_assignment"] != "unassigned")].copy()
    adata_filtered.write_h5ad("data/perturbseq_qc_filtered.h5ad")

    print(f"✓ QC complete: {stats['pass_qc']}/{stats['total_cells']} cells pass QC")
    print(f"  Unassigned: {stats['unassigned']}, Multi-guide: {stats['multi_guide']}")
    print(f"  Mean confidence: {stats['mean_confidence']:.3f}")
    print(f"✓ Filtered data saved: {adata_filtered.shape}")

    return adata_filtered, stats


if __name__ == "__main__":
    run_qc_pipeline()
