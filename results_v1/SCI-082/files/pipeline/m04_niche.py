"""
Module 4 — Tissue Niche Identification
Identifies microenvironmental niches by integrating cell-type composition
and spatially-aware clustering.
"""
from __future__ import annotations

import logging

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


# ── Cell-type neighborhood profile ──────────────────────────────────────────
def compute_neighborhood_profile(
    adata: ad.AnnData,
    cell_type_key: str = "cell_type",
    n_neighs: int = 15,
    coord_type: str = "generic",
) -> pd.DataFrame:
    """
    For each spot/cell, compute the proportion of each cell type
    among its spatial neighbors → neighborhood composition profile.
    """
    # Ensure spatial graph
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type=coord_type)

    # One-hot encode cell types
    ct_dummies = pd.get_dummies(adata.obs[cell_type_key])
    ct_matrix = ct_dummies.values.astype(float)

    # Adjacency-weighted average of neighbor cell types
    from scipy.sparse import issparse

    adj = adata.obsp["spatial_connectivities"]
    if issparse(adj):
        neighbor_profile = adj.dot(ct_matrix)
    else:
        neighbor_profile = np.dot(adj, ct_matrix)

    # Normalize to proportions
    row_sums = neighbor_profile.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    neighbor_profile = neighbor_profile / row_sums

    profile_df = pd.DataFrame(
        neighbor_profile,
        index=adata.obs_names,
        columns=[f"nhood_{c}" for c in ct_dummies.columns],
    )
    logger.info(
        "Neighborhood profile: %d spots × %d cell-type features",
        *profile_df.shape,
    )
    return profile_df


# ── Cell-type abundance profile (from deconvolution) ─────────────────────────
def get_deconv_profile(adata: ad.AnnData) -> pd.DataFrame | None:
    """Extract cell2location abundance as niche features."""
    key = "q05_cell_abundance_w_sf"
    if key in adata.obsm:
        df = adata.obsm[key].copy()
        df.columns = [f"deconv_{c}" for c in df.columns]
        logger.info("Deconvolution profile: %d spots × %d types", *df.shape)
        return df
    logger.warning("No deconvolution abundances found; skipping.")
    return None


# ── Niche clustering ─────────────────────────────────────────────────────────
def identify_niches(
    adata: ad.AnnData,
    features_df: pd.DataFrame,
    n_clusters: int | None = None,
    max_k: int = 10,
    algorithm: str = "leiden",
    resolution: float = 0.8,
    random_state: int = 42,
) -> pd.Series:
    """
    Cluster spots into tissue niches based on local composition features.

    Parameters
    ----------
    n_clusters : int or None
        If None, auto-select via silhouette score (KMeans path).
    algorithm : str
        'leiden' (graph-based) or 'kmeans'.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(features_df.values)

    if algorithm == "kmeans":
        if n_clusters is None:
            # Auto-select k via silhouette
            scores = {}
            for k in range(2, max_k + 1):
                km = KMeans(n_clusters=k, random_state=random_state, n_init=10)
                labels = km.fit_predict(X_scaled)
                scores[k] = silhouette_score(X_scaled, labels)
            n_clusters = max(scores, key=scores.get)
            logger.info("Auto-selected k=%d (silhouette=%.3f)", n_clusters, scores[n_clusters])

        km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
        niche_labels = pd.Series(
            km.fit_predict(X_scaled).astype(str),
            index=features_df.index,
            name="niche",
        )
    elif algorithm == "leiden":
        # PCA → neighbors → Leiden on the feature space
        adata_tmp = ad.AnnData(X=X_scaled, obs=pd.DataFrame(index=features_df.index))
        sc.tl.pca(adata_tmp, n_comps=min(20, X_scaled.shape[1]))
        sc.pp.neighbors(adata_tmp, n_pcs=min(20, X_scaled.shape[1]))
        sc.tl.leiden(adata_tmp, resolution=resolution, key_added="niche")
        niche_labels = adata_tmp.obs["niche"].copy()
    else:
        raise ValueError(f"Unknown clustering algorithm: {algorithm}")

    adata.obs["niche"] = niche_labels.values
    n_niches = niche_labels.nunique()
    logger.info("Identified %d tissue niches (algorithm=%s)", n_niches, algorithm)
    return niche_labels


# ── Niche characterization ───────────────────────────────────────────────────
def characterize_niches(
    adata: ad.AnnData,
    features_df: pd.DataFrame,
    niche_key: str = "niche",
) -> pd.DataFrame:
    """
    Summarize mean composition per niche → identify dominant cell types.
    """
    features_df = features_df.copy()
    features_df["niche"] = adata.obs[niche_key].values
    summary = features_df.groupby("niche").mean()
    logger.info("Niche characterization: %d niches × %d features", *summary.shape)
    return summary


# ── Visualization ────────────────────────────────────────────────────────────
def plot_niche_map(
    adata: ad.AnnData,
    niche_key: str = "niche",
    save_path: str = "figures/niche_map.png",
) -> None:
    """Spatial scatter of niche assignments."""
    coords = adata.obsm["spatial"]
    niches = adata.obs[niche_key].astype("category")
    n_niches = niches.cat.categories.size

    fig, ax = plt.subplots(figsize=(10, 10))
    cmap = plt.cm.get_cmap("tab20", n_niches)
    for i, niche in enumerate(niches.cat.categories):
        mask = niches == niche
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[cmap(i)],
            s=4,
            label=f"Niche {niche}",
            edgecolors="none",
        )
    ax.legend(markerscale=4, fontsize=8, loc="upper right")
    ax.set_title("Tissue Microenvironment Niches", fontsize=14)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Niche map saved → %s", save_path)


def plot_niche_composition(
    summary_df: pd.DataFrame,
    save_path: str = "figures/niche_composition.png",
) -> None:
    """Stacked bar chart of cell-type composition per niche."""
    fig, ax = plt.subplots(figsize=(12, 6))
    summary_df.plot(kind="bar", stacked=True, ax=ax, cmap="tab20", width=0.8)
    ax.set_xlabel("Niche")
    ax.set_ylabel("Mean proportion")
    ax.set_title("Cell-type composition per niche")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Niche composition plot saved → %s", save_path)


def save_niche_results(
    niche_labels: pd.Series,
    summary_df: pd.DataFrame,
    labels_path: str = "results/niche_assignments.csv",
    summary_path: str = "results/niche_summary.csv",
) -> None:
    niche_labels.to_csv(labels_path)
    summary_df.to_csv(summary_path)
    logger.info("Niche results saved → %s, %s", labels_path, summary_path)
