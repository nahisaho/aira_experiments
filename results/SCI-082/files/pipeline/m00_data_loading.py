"""
Module 0 — Data Loading, QC, and Preprocessing
Supports Visium (10x), MERFISH, and Slide-seq platforms.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq

logger = logging.getLogger(__name__)


# ── Platform-specific loaders ───────────────────────────────────────────────
def load_visium(input_dir: str | Path) -> ad.AnnData:
    """Load 10x Visium data from Space Ranger output."""
    adata = sc.read_visium(input_dir)
    adata.var_names_make_unique()
    adata.obs["platform"] = "visium"
    logger.info("Loaded Visium data: %s spots × %s genes", *adata.shape)
    return adata


def load_merfish(input_dir: str | Path) -> ad.AnnData:
    """Load MERFISH data (cell-level count matrix + spatial coords)."""
    adata = sc.read_h5ad(Path(input_dir) / "merfish.h5ad")
    adata.obs["platform"] = "merfish"
    logger.info("Loaded MERFISH data: %s cells × %s genes", *adata.shape)
    return adata


def load_slideseq(input_dir: str | Path) -> ad.AnnData:
    """Load Slide-seq data (bead-level)."""
    adata = sc.read_h5ad(Path(input_dir) / "slideseq.h5ad")
    adata.obs["platform"] = "slideseq"
    logger.info("Loaded Slide-seq data: %s beads × %s genes", *adata.shape)
    return adata


LOADERS = {
    "visium": load_visium,
    "merfish": load_merfish,
    "slideseq": load_slideseq,
}


def load_spatial_data(platform: str, input_dir: str | Path) -> ad.AnnData:
    loader = LOADERS.get(platform)
    if loader is None:
        raise ValueError(f"Unsupported platform: {platform}")
    return loader(input_dir)


# ── Quality control ─────────────────────────────────────────────────────────
def run_qc(
    adata: ad.AnnData,
    min_counts: int = 500,
    max_counts: int = 50_000,
    min_genes: int = 200,
    max_genes: int = 8_000,
    max_pct_mito: float = 20.0,
    min_cells_per_gene: int = 10,
) -> ad.AnnData:
    """Apply QC filters and annotate mitochondrial / ribosomal fractions."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, inplace=True
    )

    n_before = adata.n_obs
    adata = adata[
        (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["total_counts"] <= max_counts)
        & (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito),
        :,
    ].copy()

    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    n_after = adata.n_obs
    logger.info("QC: %d → %d spots/cells retained (removed %d)", n_before, n_after, n_before - n_after)
    return adata


# ── Normalization ────────────────────────────────────────────────────────────
def normalize(adata: ad.AnnData, method: str = "log_normalize") -> ad.AnnData:
    """Normalize counts. Keeps raw layer for downstream methods."""
    adata.layers["counts"] = adata.X.copy()

    if method == "log_normalize":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif method == "sctransform":
        # Placeholder: scvi-tools or R-based SCTransform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        raise ValueError(f"Unknown normalization: {method}")

    logger.info("Normalization complete: method=%s", method)
    return adata


# ── Dimensionality reduction & clustering ────────────────────────────────────
def reduce_and_cluster(
    adata: ad.AnnData,
    n_hvg: int = 3000,
    n_pcs: int = 30,
    resolution: float = 0.8,
) -> ad.AnnData:
    """HVG selection → PCA → neighbors → UMAP → Leiden clustering."""
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts")
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    logger.info("Clustering complete: %d clusters found", adata.obs["leiden"].nunique())
    return adata


# ── Spatial graph ────────────────────────────────────────────────────────────
def build_spatial_graph(
    adata: ad.AnnData,
    coord_type: str = "generic",
    n_neighs: int = 6,
    radius: Optional[float] = None,
) -> ad.AnnData:
    """Build spatial connectivity graph with Squidpy."""
    if radius:
        sq.gr.spatial_neighbors(adata, coord_type=coord_type, radius=radius)
    else:
        sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)
    logger.info("Spatial graph built: coord_type=%s, n_neighs=%s", coord_type, n_neighs)
    return adata


# ── QC plots ─────────────────────────────────────────────────────────────────
def plot_qc(adata: ad.AnnData, save_path: str = "figures/qc_violin.png") -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    sc.pl.violin(adata, keys="total_counts", ax=axes[0], show=False)
    sc.pl.violin(adata, keys="n_genes_by_counts", ax=axes[1], show=False)
    sc.pl.violin(adata, keys="pct_counts_mt", ax=axes[2], show=False)
    sc.pl.violin(adata, keys="pct_counts_ribo", ax=axes[3], show=False)
    for ax, title in zip(axes, ["Total counts", "Genes detected", "% Mito", "% Ribo"]):
        ax.set_title(title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("QC violin plot saved → %s", save_path)
