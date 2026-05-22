"""
Module 5 — 3D Spatial Reconstruction from Serial Sections
Aligns consecutive tissue sections and integrates into a 3D coordinate system.
"""
from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import procrustes

logger = logging.getLogger(__name__)


# ── Section alignment ───────────────────────────────────────────────────────
def pairwise_icp(
    source: np.ndarray,
    target: np.ndarray,
    max_iterations: int = 200,
    tolerance: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Iterative Closest Point (ICP) alignment of 2D point clouds.

    Parameters
    ----------
    source, target : ndarray of shape (N, 2) and (M, 2)
    max_iterations : int
    tolerance : float – convergence threshold on RMSE delta.

    Returns
    -------
    T           : (3, 3) homogeneous transformation matrix
    aligned     : (N, 2) transformed source coordinates
    final_error : float – final RMSE
    """
    from sklearn.neighbors import NearestNeighbors

    src = source.copy()
    prev_error = np.inf

    T_total = np.eye(3)

    for i in range(max_iterations):
        # Find nearest neighbors in target
        nn = NearestNeighbors(n_neighbors=1, algorithm="kd_tree")
        nn.fit(target)
        distances, indices = nn.kneighbors(src)
        matched_target = target[indices.ravel()]

        # Compute centroids
        centroid_src = src.mean(axis=0)
        centroid_tgt = matched_target.mean(axis=0)

        # Center
        src_c = src - centroid_src
        tgt_c = matched_target - centroid_tgt

        # SVD for optimal rotation
        H = src_c.T @ tgt_c
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Correct reflection
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        t = centroid_tgt - R @ centroid_src

        # Build homogeneous transform
        T_step = np.eye(3)
        T_step[:2, :2] = R
        T_step[:2, 2] = t

        T_total = T_step @ T_total
        src = (R @ src.T).T + t

        # Check convergence
        mean_error = distances.mean()
        if abs(prev_error - mean_error) < tolerance:
            logger.debug("ICP converged at iteration %d (error=%.6f)", i, mean_error)
            break
        prev_error = mean_error

    return T_total, src, prev_error


def align_sections(
    adatas: dict[str, ad.AnnData],
    section_order: list[str],
    max_iterations: int = 200,
    tolerance: float = 1e-6,
) -> dict[str, np.ndarray]:
    """
    Sequentially align serial sections using ICP.
    The first section serves as the reference frame.

    Returns dict mapping section_id → aligned (N, 2) coordinates.
    """
    aligned_coords = {}
    ref_id = section_order[0]
    aligned_coords[ref_id] = adatas[ref_id].obsm["spatial"].copy()

    for i in range(1, len(section_order)):
        prev_id = section_order[i - 1]
        curr_id = section_order[i]

        target_pts = aligned_coords[prev_id]
        source_pts = adatas[curr_id].obsm["spatial"].copy()

        _, aligned, error = pairwise_icp(
            source_pts, target_pts, max_iterations=max_iterations, tolerance=tolerance
        )
        aligned_coords[curr_id] = aligned
        logger.info(
            "Aligned %s → %s: RMSE=%.4f",
            curr_id,
            prev_id,
            error,
        )

    return aligned_coords


# ── 3D coordinate assembly ───────────────────────────────────────────────────
def build_3d_coordinates(
    adatas: dict[str, ad.AnnData],
    aligned_coords: dict[str, np.ndarray],
    section_order: list[str],
    z_spacing_um: float = 10.0,
) -> pd.DataFrame:
    """
    Combine aligned 2D coordinates with z-axis (section depth).

    Returns DataFrame with columns: barcode, section, x, y, z.
    """
    records = []
    for z_idx, sec_id in enumerate(section_order):
        coords = aligned_coords[sec_id]
        barcodes = adatas[sec_id].obs_names
        z = z_idx * z_spacing_um
        for j, bc in enumerate(barcodes):
            records.append(
                {
                    "barcode": bc,
                    "section": sec_id,
                    "x": coords[j, 0],
                    "y": coords[j, 1],
                    "z": z,
                }
            )
    df = pd.DataFrame(records)
    logger.info(
        "3D coordinates assembled: %d points across %d sections",
        len(df),
        len(section_order),
    )
    return df


# ── Concatenate sections into single AnnData ─────────────────────────────────
def merge_sections(
    adatas: dict[str, ad.AnnData],
    coords_3d: pd.DataFrame,
    section_order: list[str],
) -> ad.AnnData:
    """
    Merge all sections into a single AnnData with 3D spatial coordinates.
    """
    adata_list = []
    for sec_id in section_order:
        a = adatas[sec_id].copy()
        a.obs["section"] = sec_id
        sec_coords = coords_3d.loc[coords_3d["section"] == sec_id, ["x", "y", "z"]].values
        a.obsm["spatial_3d"] = sec_coords
        adata_list.append(a)

    adata_merged = ad.concat(adata_list, join="inner")
    adata_merged.obs_names_make_unique()
    logger.info("Merged AnnData: %d spots × %d genes", *adata_merged.shape)
    return adata_merged


# ── Visualization ────────────────────────────────────────────────────────────
def plot_3d_reconstruction(
    coords_3d: pd.DataFrame,
    color_by: str = "section",
    save_path: str = "figures/3d_reconstruction.png",
) -> None:
    """3D scatter plot of reconstructed tissue coordinates."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    groups = coords_3d.groupby(color_by)
    cmap = plt.cm.get_cmap("tab10", len(groups))

    for i, (name, group) in enumerate(groups):
        ax.scatter(
            group["x"],
            group["y"],
            group["z"],
            s=1,
            c=[cmap(i)],
            label=str(name),
            alpha=0.6,
        )

    ax.set_xlabel("X (μm)")
    ax.set_ylabel("Y (μm)")
    ax.set_zlabel("Z (μm)")
    ax.set_title("3D Spatial Reconstruction")
    ax.legend(markerscale=6, fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("3D reconstruction plot saved → %s", save_path)


def plot_3d_gene_expression(
    coords_3d: pd.DataFrame,
    adata_merged: ad.AnnData,
    gene: str,
    save_path: str = "figures/3d_gene_expression.png",
) -> None:
    """3D scatter colored by gene expression."""
    if gene not in adata_merged.var_names:
        logger.warning("Gene '%s' not in merged AnnData. Skipping.", gene)
        return

    expr = adata_merged[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = expr.ravel()

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")
    sc_ = ax.scatter(
        coords_3d["x"],
        coords_3d["y"],
        coords_3d["z"],
        c=expr,
        s=1,
        cmap="viridis",
        alpha=0.7,
    )
    ax.set_xlabel("X (μm)")
    ax.set_ylabel("Y (μm)")
    ax.set_zlabel("Z (μm)")
    ax.set_title(f"3D expression: {gene}")
    plt.colorbar(sc_, ax=ax, label="Expression", shrink=0.5)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("3D gene expression plot saved → %s", save_path)


def save_3d_coordinates(
    coords_3d: pd.DataFrame,
    save_path: str = "results/coordinates_3d.csv",
) -> None:
    coords_3d.to_csv(save_path, index=False)
    logger.info("3D coordinates saved → %s (%d points)", save_path, len(coords_3d))
