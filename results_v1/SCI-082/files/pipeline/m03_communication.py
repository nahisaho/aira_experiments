"""
Module 3 — Cell–Cell Communication Analysis
Integrates LIANA (consensus LR scoring) with Squidpy spatial context.
"""
from __future__ import annotations

import logging

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq

logger = logging.getLogger(__name__)


# ── LIANA ligand-receptor inference ──────────────────────────────────────────
def run_liana(
    adata: ad.AnnData,
    groupby: str = "cell_type",
    resource_name: str = "consensus",
    min_expr_prop: float = 0.1,
    use_raw: bool = False,
) -> pd.DataFrame:
    """
    Run LIANA multi-method consensus scoring for ligand–receptor interactions.

    Parameters
    ----------
    groupby : str
        Column in adata.obs with cell-type labels (or deconvolution-derived labels).
    resource_name : str
        LR database: 'consensus', 'cellphonedb', 'cellchat', 'connectomedb2020'.
    min_expr_prop : float
        Minimum proportion of cells expressing ligand/receptor.

    Returns
    -------
    pd.DataFrame with columns including:
        source, target, ligand_complex, receptor_complex,
        magnitude_rank, specificity_rank, ...
    """
    import liana as li

    li.mt.rank_aggregate(
        adata,
        groupby=groupby,
        resource_name=resource_name,
        expr_prop=min_expr_prop,
        use_raw=use_raw,
        verbose=True,
    )
    lr_res = adata.uns["liana_res"].copy()
    logger.info(
        "LIANA: %d LR interactions detected across %d cell-type pairs",
        len(lr_res),
        lr_res.groupby(["source", "target"]).ngroups,
    )
    return lr_res


# ── Squidpy spatial LR co-expression ────────────────────────────────────────
def spatial_lr_coexpression(
    adata: ad.AnnData,
    lr_pairs: list[tuple[str, str]],
    n_perms: int = 1000,
    n_jobs: int = 4,
) -> pd.DataFrame:
    """
    Evaluate spatial co-expression of LR pairs using Squidpy's
    co-occurrence or interaction matrix.
    """
    results = []
    for ligand, receptor in lr_pairs:
        if ligand not in adata.var_names or receptor not in adata.var_names:
            continue
        # Bivariate Moran's I (spatial cross-correlation)
        sq.gr.spatial_autocorr(
            adata[:, [ligand, receptor]],
            mode="moran",
            n_perms=n_perms,
            n_jobs=n_jobs,
        )
        mi = adata[:, [ligand, receptor]].uns.get("moranI", None)
        if mi is not None:
            results.append(
                {
                    "ligand": ligand,
                    "receptor": receptor,
                    "ligand_moran_I": mi.loc[ligand, "I"] if ligand in mi.index else np.nan,
                    "receptor_moran_I": mi.loc[receptor, "I"] if receptor in mi.index else np.nan,
                }
            )
    df = pd.DataFrame(results)
    logger.info("Spatial LR co-expression: %d pairs evaluated", len(df))
    return df


# ── Squidpy interaction matrix ───────────────────────────────────────────────
def compute_interaction_matrix(
    adata: ad.AnnData,
    cluster_key: str = "cell_type",
    n_perms: int = 1000,
) -> None:
    """Compute neighborhood enrichment (interaction) scores."""
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    logger.info("Neighborhood enrichment computed for key '%s'", cluster_key)


# ── Prioritize spatially-resolved interactions ───────────────────────────────
def prioritize_interactions(
    liana_df: pd.DataFrame,
    spatial_lr_df: pd.DataFrame | None = None,
    top_n: int = 50,
    rank_col: str = "magnitude_rank",
) -> pd.DataFrame:
    """
    Rank interactions by LIANA magnitude_rank, optionally boosted
    by spatial LR co-expression evidence.
    """
    df = liana_df.sort_values(rank_col, ascending=True).head(top_n).copy()
    if spatial_lr_df is not None and not spatial_lr_df.empty:
        df = df.merge(
            spatial_lr_df,
            left_on=["ligand_complex", "receptor_complex"],
            right_on=["ligand", "receptor"],
            how="left",
        )
    logger.info("Top %d interactions prioritized", len(df))
    return df


# ── Visualization ────────────────────────────────────────────────────────────
def plot_communication_network(
    lr_df: pd.DataFrame,
    top_n: int = 30,
    save_path: str = "figures/communication_network.png",
) -> None:
    """Plot a heatmap of top LR interaction scores between cell types."""
    top = lr_df.sort_values("magnitude_rank").head(top_n)
    pivot = top.pivot_table(
        index="source",
        columns="target",
        values="magnitude_rank",
        aggfunc="count",
        fill_value=0,
    )

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("Target cell type")
    ax.set_ylabel("Source cell type")
    ax.set_title(f"Cell–cell communication (top {top_n} LR pairs)")
    plt.colorbar(im, ax=ax, label="Number of interactions")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Communication network plot saved → %s", save_path)


def plot_nhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str = "cell_type",
    save_path: str = "figures/nhood_enrichment.png",
) -> None:
    """Plot neighborhood enrichment heatmap."""
    fig, ax = plt.subplots(figsize=(8, 8))
    sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, ax=ax)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Neighborhood enrichment plot saved → %s", save_path)


def save_lr_results(
    df: pd.DataFrame,
    save_path: str = "results/ligand_receptor_results.csv",
) -> None:
    df.to_csv(save_path, index=False)
    logger.info("LR results saved → %s (%d interactions)", save_path, len(df))
