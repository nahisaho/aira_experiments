"""
Module 6 — Tumor–Immune Microenvironment (TIME) Case Study
Comprehensive analysis of immune infiltration, exhaustion, and
spatial organisation around tumor boundaries.
"""
from __future__ import annotations

import logging
from typing import Optional

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


# ── Tumor boundary detection ────────────────────────────────────────────────
def define_tumor_boundary(
    adata: ad.AnnData,
    tumor_label: str = "Tumor",
    cell_type_key: str = "cell_type",
    boundary_distance_um: float = 200.0,
) -> pd.DataFrame:
    """
    Classify spots into spatial zones relative to the tumor boundary:
        - tumor_core     : tumor spots far from boundary
        - tumor_edge     : tumor spots near boundary
        - stroma_near    : stroma spots near tumor
        - stroma_far     : stroma spots far from tumor

    Returns DataFrame with zone assignments and distance to tumor boundary.
    """
    is_tumor = adata.obs[cell_type_key] == tumor_label
    coords = adata.obsm["spatial"]
    tumor_coords = coords[is_tumor]
    non_tumor_coords = coords[~is_tumor]

    # Distance from every spot to nearest tumor spot
    if len(tumor_coords) == 0:
        logger.warning("No tumor spots found. Cannot define boundary.")
        return pd.DataFrame()

    dist_to_tumor = cdist(coords, tumor_coords).min(axis=1)

    # Distance from every tumor spot to nearest non-tumor spot
    dist_tumor_to_stroma = cdist(tumor_coords, non_tumor_coords).min(axis=1)

    zones = pd.Series("unknown", index=adata.obs_names, name="tumor_zone")
    tumor_mask = is_tumor.values

    # Tumor spots
    tumor_boundary_dist = np.full(len(adata), np.nan)
    tumor_boundary_dist[tumor_mask] = dist_tumor_to_stroma
    tumor_boundary_dist[~tumor_mask] = dist_to_tumor[~tumor_mask]

    zones[tumor_mask & (dist_tumor_to_stroma <= boundary_distance_um)] = "tumor_edge"
    zones[tumor_mask & (dist_tumor_to_stroma > boundary_distance_um)] = "tumor_core"
    zones[~tumor_mask & (dist_to_tumor <= boundary_distance_um)] = "stroma_near"
    zones[~tumor_mask & (dist_to_tumor > boundary_distance_um)] = "stroma_far"

    result = pd.DataFrame(
        {
            "tumor_zone": zones,
            "dist_to_boundary_um": tumor_boundary_dist,
        },
        index=adata.obs_names,
    )
    adata.obs["tumor_zone"] = zones.values
    adata.obs["dist_to_boundary_um"] = tumor_boundary_dist

    for z in ["tumor_core", "tumor_edge", "stroma_near", "stroma_far"]:
        n = (zones == z).sum()
        logger.info("Zone '%s': %d spots", z, n)

    return result


# ── Immune infiltration gradient ─────────────────────────────────────────────
def compute_immune_gradient(
    adata: ad.AnnData,
    immune_cell_types: list[str],
    cell_type_key: str = "cell_type",
    distance_bins: list[float] | None = None,
) -> pd.DataFrame:
    """
    Quantify immune cell density as a function of distance from tumor boundary.
    """
    if distance_bins is None:
        distance_bins = [0, 50, 100, 200, 500, 1000]

    dist = adata.obs["dist_to_boundary_um"].values
    ct = adata.obs[cell_type_key].values
    is_tumor_side = adata.obs["tumor_zone"].isin(["tumor_core", "tumor_edge"]).values

    # Focus on stroma side distances
    records = []
    for i in range(len(distance_bins) - 1):
        lo, hi = distance_bins[i], distance_bins[i + 1]
        mask = (~is_tumor_side) & (dist >= lo) & (dist < hi)
        n_total = mask.sum()
        if n_total == 0:
            continue
        for ict in immune_cell_types:
            n_immune = ((ct == ict) & mask).sum()
            records.append(
                {
                    "distance_bin": f"{lo}-{hi}",
                    "cell_type": ict,
                    "count": n_immune,
                    "total_spots": n_total,
                    "density": n_immune / n_total,
                }
            )

    df = pd.DataFrame(records)
    logger.info("Immune gradient: %d distance bins × %d immune types", len(distance_bins) - 1, len(immune_cell_types))
    return df


# ── Exhaustion scoring ───────────────────────────────────────────────────────
def score_exhaustion(
    adata: ad.AnnData,
    exhaustion_markers: list[str],
    score_name: str = "exhaustion_score",
) -> None:
    """
    Compute exhaustion module score (scanpy.tl.score_genes) for
    T-cell exhaustion markers (PDCD1, LAG3, HAVCR2, TIGIT, TOX).
    """
    available = [g for g in exhaustion_markers if g in adata.var_names]
    if len(available) < 2:
        logger.warning("Too few exhaustion markers found (%d). Skipping.", len(available))
        return

    sc.tl.score_genes(adata, gene_list=available, score_name=score_name)
    logger.info(
        "Exhaustion score computed using %d/%d markers (mean=%.3f)",
        len(available),
        len(exhaustion_markers),
        adata.obs[score_name].mean(),
    )


def score_checkpoint_ligands(
    adata: ad.AnnData,
    checkpoint_ligands: list[str],
    score_name: str = "checkpoint_ligand_score",
) -> None:
    """Score checkpoint ligand expression (CD274/PD-L1, PDCD1LG2/PD-L2, etc.)."""
    available = [g for g in checkpoint_ligands if g in adata.var_names]
    if len(available) < 1:
        logger.warning("No checkpoint ligands found. Skipping.")
        return
    sc.tl.score_genes(adata, gene_list=available, score_name=score_name)
    logger.info("Checkpoint ligand score computed using %d markers", len(available))


# ── Spatial immune–tumor interaction hotspots ────────────────────────────────
def identify_interaction_hotspots(
    adata: ad.AnnData,
    tumor_label: str = "Tumor",
    immune_labels: list[str] | None = None,
    cell_type_key: str = "cell_type",
    n_neighs: int = 10,
) -> pd.DataFrame:
    """
    Identify hotspots where tumor and immune cells co-localize,
    based on spatial neighbor composition.
    """
    if immune_labels is None:
        immune_labels = ["CD8_T", "CD4_T", "Macrophage", "NK"]

    sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic")

    adj = adata.obsp["spatial_connectivities"]
    ct = adata.obs[cell_type_key].values

    records = []
    tumor_idx = np.where(ct == tumor_label)[0]
    for idx in tumor_idx:
        neighbors = adj[idx].nonzero()[1]
        neighbor_types = ct[neighbors]
        for ict in immune_labels:
            n_immune_neighbors = (neighbor_types == ict).sum()
            if n_immune_neighbors > 0:
                records.append(
                    {
                        "tumor_spot": adata.obs_names[idx],
                        "immune_type": ict,
                        "n_immune_neighbors": n_immune_neighbors,
                        "x": adata.obsm["spatial"][idx, 0],
                        "y": adata.obsm["spatial"][idx, 1],
                    }
                )

    df = pd.DataFrame(records)
    logger.info("Interaction hotspots: %d tumor–immune contacts identified", len(df))
    return df


# ── Visualization ────────────────────────────────────────────────────────────
def plot_tumor_zones(
    adata: ad.AnnData,
    save_path: str = "figures/tumor_zones.png",
) -> None:
    """Spatial map colored by tumor zone."""
    coords = adata.obsm["spatial"]
    zones = adata.obs["tumor_zone"]
    zone_colors = {
        "tumor_core": "#d62728",
        "tumor_edge": "#ff7f0e",
        "stroma_near": "#2ca02c",
        "stroma_far": "#1f77b4",
        "unknown": "#7f7f7f",
    }

    fig, ax = plt.subplots(figsize=(10, 10))
    for zone, color in zone_colors.items():
        mask = zones == zone
        if mask.sum() == 0:
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=color,
            s=4,
            label=zone,
            edgecolors="none",
        )
    ax.legend(markerscale=4, fontsize=10)
    ax.set_title("Tumor Microenvironment Zones")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Tumor zone map saved → %s", save_path)


def plot_immune_gradient(
    gradient_df: pd.DataFrame,
    save_path: str = "figures/immune_gradient.png",
) -> None:
    """Bar plot of immune cell density vs. distance from tumor."""
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(
        data=gradient_df,
        x="distance_bin",
        y="density",
        hue="cell_type",
        ax=ax,
        palette="Set2",
    )
    ax.set_xlabel("Distance from tumor boundary (μm)")
    ax.set_ylabel("Immune cell density")
    ax.set_title("Immune Infiltration Gradient")
    ax.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Immune gradient plot saved → %s", save_path)


def plot_exhaustion_spatial(
    adata: ad.AnnData,
    score_name: str = "exhaustion_score",
    save_path: str = "figures/exhaustion_spatial.png",
) -> None:
    """Spatial map of exhaustion score."""
    if score_name not in adata.obs.columns:
        logger.warning("Score '%s' not found. Skipping plot.", score_name)
        return

    coords = adata.obsm["spatial"]
    fig, ax = plt.subplots(figsize=(10, 10))
    sc_ = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=adata.obs[score_name].values,
        s=4,
        cmap="inferno",
        edgecolors="none",
    )
    ax.set_title(f"T-cell Exhaustion Score (spatial)")
    ax.axis("off")
    plt.colorbar(sc_, ax=ax, label=score_name, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Exhaustion spatial map saved → %s", save_path)


def plot_tumor_immune_landscape(
    adata: ad.AnnData,
    save_path: str = "figures/tumor_immune_landscape.png",
) -> None:
    """Multi-panel overview of the tumor-immune landscape."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    coords = adata.obsm["spatial"]

    # Panel 1: Tumor zones
    zone_colors = {
        "tumor_core": "#d62728",
        "tumor_edge": "#ff7f0e",
        "stroma_near": "#2ca02c",
        "stroma_far": "#1f77b4",
    }
    for zone, color in zone_colors.items():
        mask = adata.obs["tumor_zone"] == zone
        if mask.sum():
            axes[0, 0].scatter(coords[mask, 0], coords[mask, 1], c=color, s=2, label=zone)
    axes[0, 0].legend(markerscale=4, fontsize=8)
    axes[0, 0].set_title("Tumor Zones")
    axes[0, 0].axis("off")

    # Panel 2: Distance to boundary
    sc1 = axes[0, 1].scatter(
        coords[:, 0], coords[:, 1],
        c=adata.obs.get("dist_to_boundary_um", np.zeros(len(adata))),
        s=2, cmap="coolwarm",
    )
    axes[0, 1].set_title("Distance to Tumor Boundary (μm)")
    axes[0, 1].axis("off")
    plt.colorbar(sc1, ax=axes[0, 1], fraction=0.046)

    # Panel 3: Exhaustion score
    if "exhaustion_score" in adata.obs:
        sc2 = axes[1, 0].scatter(
            coords[:, 0], coords[:, 1],
            c=adata.obs["exhaustion_score"], s=2, cmap="inferno",
        )
        plt.colorbar(sc2, ax=axes[1, 0], fraction=0.046)
    axes[1, 0].set_title("T-cell Exhaustion Score")
    axes[1, 0].axis("off")

    # Panel 4: Checkpoint ligand score
    if "checkpoint_ligand_score" in adata.obs:
        sc3 = axes[1, 1].scatter(
            coords[:, 0], coords[:, 1],
            c=adata.obs["checkpoint_ligand_score"], s=2, cmap="magma",
        )
        plt.colorbar(sc3, ax=axes[1, 1], fraction=0.046)
    axes[1, 1].set_title("Checkpoint Ligand Score")
    axes[1, 1].axis("off")

    plt.suptitle("Tumor–Immune Microenvironment Landscape", fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Tumor–immune landscape saved → %s", save_path)


def save_time_report(
    zone_df: pd.DataFrame,
    gradient_df: pd.DataFrame,
    hotspots_df: pd.DataFrame,
    save_path: str = "results/tumor_immune_report.csv",
) -> None:
    """Save combined TIME analysis results."""
    zone_df.to_csv(save_path.replace(".csv", "_zones.csv"))
    gradient_df.to_csv(save_path.replace(".csv", "_gradient.csv"), index=False)
    hotspots_df.to_csv(save_path.replace(".csv", "_hotspots.csv"), index=False)
    logger.info("TIME report saved → %s*", save_path)
