"""
Module 2 — Spatially Variable Gene (SVG) Detection
Combines Squidpy (Moran's I / Sepal) and SpatialDE for robust SVG identification.
"""
from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

logger = logging.getLogger(__name__)


# ── Moran's I (Squidpy) ─────────────────────────────────────────────────────
def compute_morans_i(
    adata: ad.AnnData,
    n_perms: int = 1000,
    n_jobs: int = 4,
) -> pd.DataFrame:
    """
    Compute Moran's I spatial autocorrelation for all genes.

    Returns a DataFrame with columns: I, pval_norm, var_norm, pval_sim.
    """
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        n_perms=n_perms,
        n_jobs=n_jobs,
    )
    df = adata.uns["moranI"].copy()
    df = df.sort_values("I", ascending=False)
    logger.info(
        "Moran's I computed for %d genes; top gene I=%.4f",
        len(df),
        df["I"].iloc[0],
    )
    return df


# ── SpatialDE ────────────────────────────────────────────────────────────────
def run_spatialde(
    adata: ad.AnnData,
    layer: str | None = None,
) -> pd.DataFrame:
    """
    Run SpatialDE Gaussian-process-based SVG detection.

    Returns DataFrame with columns: g (gene), pval, qval, l (length scale),
    FSV (fraction spatial variance).
    """
    import SpatialDE

    counts = (
        pd.DataFrame(
            adata.layers[layer].toarray()
            if hasattr(adata.layers.get(layer, adata.X), "toarray")
            else (adata.layers[layer] if layer else adata.X),
            index=adata.obs_names,
            columns=adata.var_names,
        )
    )
    coords = pd.DataFrame(
        adata.obsm["spatial"],
        index=adata.obs_names,
        columns=["x", "y"],
    )

    # Stabilize variance
    norm_expr = SpatialDE.stabilize(counts.T).T
    resid_expr = SpatialDE.regress_out(coords, norm_expr.T, "1 + x + y").T

    results = SpatialDE.run(coords, resid_expr)
    results = results.sort_values("qval")

    logger.info(
        "SpatialDE: %d SVGs at q < 0.05",
        (results["qval"] < 0.05).sum(),
    )
    return results


# ── Sepal score (Squidpy diffusion-based) ────────────────────────────────────
def compute_sepal_score(
    adata: ad.AnnData,
    n_jobs: int = 4,
) -> pd.DataFrame:
    """
    Compute Sepal diffusion-based spatial variability score.
    Higher score → more spatially structured expression.
    """
    sq.gr.spatial_autocorr(adata, mode="geary", n_jobs=n_jobs)
    df = adata.uns["gearyC"].copy()
    df = df.sort_values("C")  # lower Geary's C = stronger clustering
    logger.info("Geary's C (proxy Sepal) computed for %d genes", len(df))
    return df


# ── Consensus ranking ───────────────────────────────────────────────────────
def consensus_svg(
    moran_df: pd.DataFrame,
    spatialde_df: pd.DataFrame,
    alpha: float = 0.05,
    correction: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Combine Moran's I and SpatialDE results into a consensus SVG list.
    A gene is 'consensus SVG' if significant in BOTH methods.
    """
    from statsmodels.stats.multitest import multipletests

    # Moran significant genes
    _, moran_qval, _, _ = multipletests(moran_df["pval_sim"], method=correction)
    moran_sig = set(moran_df.index[moran_qval < alpha])

    # SpatialDE significant genes
    spde_sig = set(spatialde_df.loc[spatialde_df["qval"] < alpha, "g"])

    consensus = moran_sig & spde_sig
    logger.info(
        "Consensus SVGs: %d (Moran sig=%d, SpatialDE sig=%d)",
        len(consensus),
        len(moran_sig),
        len(spde_sig),
    )

    # Build combined table
    records = []
    for gene in consensus:
        rec = {"gene": gene}
        if gene in moran_df.index:
            rec["moran_I"] = moran_df.loc[gene, "I"]
            rec["moran_qval"] = moran_qval[list(moran_df.index).index(gene)]
        row = spatialde_df.loc[spatialde_df["g"] == gene]
        if len(row):
            rec["spatialde_qval"] = row["qval"].values[0]
            rec["FSV"] = row["FSV"].values[0]
            rec["length_scale"] = row["l"].values[0]
        records.append(rec)

    df = pd.DataFrame(records).sort_values("moran_I", ascending=False)
    return df


# ── Visualization ────────────────────────────────────────────────────────────
def plot_svg_maps(
    adata: ad.AnnData,
    genes: list[str],
    save_path: str = "figures/svg_expression_maps.png",
    n_cols: int = 4,
) -> None:
    """Plot spatial expression maps for top SVGs."""
    n_genes = len(genes)
    n_rows = int(np.ceil(n_genes / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = np.atleast_2d(axes)

    coords = adata.obsm["spatial"]
    expr_mat = adata[:, genes].X
    if hasattr(expr_mat, "toarray"):
        expr_mat = expr_mat.toarray()

    for idx, gene in enumerate(genes):
        ax = axes[idx // n_cols, idx % n_cols]
        sc_ = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=expr_mat[:, idx],
            s=4,
            cmap="viridis",
            edgecolors="none",
        )
        ax.set_title(gene, fontsize=10)
        ax.axis("off")
        plt.colorbar(sc_, ax=ax, fraction=0.046, pad=0.04)

    for idx in range(n_genes, n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].axis("off")

    plt.suptitle("Spatially Variable Genes", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("SVG maps saved → %s", save_path)


def save_svg_results(
    df: pd.DataFrame,
    save_path: str = "results/spatially_variable_genes.csv",
) -> None:
    df.to_csv(save_path, index=False)
    logger.info("SVG table saved → %s (%d genes)", save_path, len(df))
