"""
Module 1 — Spot Deconvolution with cell2location
Estimates cell-type composition per spatial spot using a reference scRNA-seq atlas.
"""
from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


# ── Reference model (NB regression) ─────────────────────────────────────────
def prepare_reference(
    adata_ref: ad.AnnData,
    cell_type_key: str = "cell_type",
    batch_key: str | None = "sample",
    n_epochs: int = 250,
    lr: float = 0.002,
) -> dict:
    """
    Train a Negative Binomial regression on reference scRNA-seq data
    to extract per-gene, per-cell-type expression signatures.

    Returns
    -------
    dict with keys:
        'inf_aver'   — pd.DataFrame (genes × cell_types) of inferred averages
        'model_ref'  — trained reference model object
    """
    import cell2location
    from cell2location.models import RegressionModel

    # Subset to protein-coding if available
    adata_ref = adata_ref.copy()
    if "counts" in adata_ref.layers:
        adata_ref.X = adata_ref.layers["counts"].copy()

    RegressionModel.setup_anndata(
        adata_ref,
        labels_key=cell_type_key,
        batch_key=batch_key,
    )
    model_ref = RegressionModel(adata_ref)
    model_ref.train(max_epochs=n_epochs, lr=lr, use_gpu=True)

    # Export estimated cell-type signatures
    adata_ref = model_ref.export_posterior(
        adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
    )
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"].copy()
    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]

    logger.info(
        "Reference model trained: %d genes × %d cell types",
        inf_aver.shape[0],
        inf_aver.shape[1],
    )
    return {"inf_aver": inf_aver, "model_ref": model_ref}


# ── Spatial mapping ──────────────────────────────────────────────────────────
def run_cell2location(
    adata_sp: ad.AnnData,
    inf_aver: pd.DataFrame,
    n_cells_per_location: int = 8,
    detection_alpha: int = 20,
    n_epochs: int = 30_000,
) -> ad.AnnData:
    """
    Map cell-type signatures onto spatial data using cell2location.

    Adds `.obsm['q05_cell_abundance_w_sf']` — posterior 5th-percentile
    cell-type abundances (recommended for downstream use).
    """
    import cell2location
    from cell2location.models import Cell2location

    adata_sp = adata_sp.copy()
    if "counts" in adata_sp.layers:
        adata_sp.X = adata_sp.layers["counts"].copy()

    # Intersect genes
    shared_genes = adata_sp.var_names.intersection(inf_aver.index)
    adata_sp = adata_sp[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes, :]

    Cell2location.setup_anndata(adata_sp)
    model_sp = Cell2location(
        adata_sp,
        cell_state_df=inf_aver,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
    )
    model_sp.train(
        max_epochs=n_epochs,
        batch_size=None,
        train_size=1,
        use_gpu=True,
    )

    adata_sp = model_sp.export_posterior(
        adata_sp, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
    )

    logger.info(
        "cell2location mapping complete: %d spots, %d cell types",
        adata_sp.n_obs,
        inf_aver.shape[1],
    )
    return adata_sp


# ── Convenience: full pipeline ───────────────────────────────────────────────
def deconvolve(
    adata_sp: ad.AnnData,
    adata_ref: ad.AnnData,
    cell_type_key: str = "cell_type",
    batch_key: str | None = "sample",
    n_cells_per_location: int = 8,
    detection_alpha: int = 20,
    n_epochs_ref: int = 250,
    n_epochs_spatial: int = 30_000,
) -> ad.AnnData:
    """End-to-end deconvolution: train reference → map spatial."""
    ref_out = prepare_reference(
        adata_ref,
        cell_type_key=cell_type_key,
        batch_key=batch_key,
        n_epochs=n_epochs_ref,
    )
    adata_sp = run_cell2location(
        adata_sp,
        inf_aver=ref_out["inf_aver"],
        n_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
        n_epochs=n_epochs_spatial,
    )
    return adata_sp


# ── Visualization ────────────────────────────────────────────────────────────
def plot_deconvolution(
    adata: ad.AnnData,
    cell_types: list[str] | None = None,
    save_path: str = "figures/deconvolution_map.png",
    n_cols: int = 4,
) -> None:
    """Plot spatial maps of estimated cell-type abundances."""
    abundance_key = "q05_cell_abundance_w_sf"
    if abundance_key not in adata.obsm:
        logger.warning("Abundance key '%s' not found. Skipping plot.", abundance_key)
        return

    df = adata.obsm[abundance_key]
    if cell_types is None:
        cell_types = list(df.columns)

    n_types = len(cell_types)
    n_rows = int(np.ceil(n_types / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = np.atleast_2d(axes)

    coords = adata.obsm["spatial"]
    for idx, ct in enumerate(cell_types):
        ax = axes[idx // n_cols, idx % n_cols]
        sc_ = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=df[ct].values,
            s=4,
            cmap="viridis",
            edgecolors="none",
        )
        ax.set_title(ct, fontsize=10)
        ax.axis("off")
        plt.colorbar(sc_, ax=ax, fraction=0.046, pad=0.04)

    # Hide unused axes
    for idx in range(n_types, n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].axis("off")

    plt.suptitle("Cell-type abundance (cell2location q05)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Deconvolution map saved → %s", save_path)


def export_abundances(
    adata: ad.AnnData,
    save_path: str = "results/cell_type_abundances.csv",
) -> pd.DataFrame:
    """Export cell-type abundance matrix as CSV."""
    abundance_key = "q05_cell_abundance_w_sf"
    if abundance_key not in adata.obsm:
        raise KeyError(f"'{abundance_key}' not in adata.obsm")
    df = adata.obsm[abundance_key].copy()
    df.index = adata.obs_names
    df.to_csv(save_path)
    logger.info("Abundances exported → %s (%d spots × %d types)", save_path, *df.shape)
    return df
