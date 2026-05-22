"""
beta_diversity.py — PCoA + PERMANOVA from MetaPhlAn 4 profiles.

Distances: Bray-Curtis, Jaccard, Aitchison (CLR + Euclidean).
"""
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

profiles = pd.read_csv(snakemake.input.profiles, sep="\t", comment="#", index_col=0)
metadata = pd.read_csv(snakemake.input.metadata, sep="\t").set_index("sample")
dist_metrics = snakemake.params.metrics

species = profiles[profiles.index.str.contains(r"s__") &
                   ~profiles.index.str.contains(r"t__")]
abundance = species.T
abundance = abundance.loc[abundance.index.isin(metadata.index)]

def clr_transform(df):
    pseudo = df.replace(0, 1e-6)
    log_data = np.log(pseudo)
    return log_data.subtract(log_data.mean(axis=1), axis=0)

dist_matrices = {}
for metric in dist_metrics:
    if metric == "aitchison":
        clr = clr_transform(abundance)
        dm_array = pdist(clr.values, metric="euclidean")
    else:
        dm_array = pdist(abundance.values, metric=metric.replace("_", ""))
    dist_matrices[metric] = DistanceMatrix(squareform(dm_array), ids=abundance.index.tolist())

# PERMANOVA
permanova_results = []
for metric, dm in dist_matrices.items():
    grouping = metadata.loc[dm.ids, "group"]
    result = permanova(dm, grouping, permutations=999)
    permanova_results.append({
        "metric": metric,
        "pseudo_F": round(result["test statistic"], 4),
        "p_value": result["p-value"],
        "permutations": 999,
        "sample_size": result["sample size"],
    })

pd.DataFrame(permanova_results).to_csv(snakemake.output.permanova, sep="\t", index=False)

# PCoA on first metric (Bray-Curtis)
primary_metric = dist_metrics[0]
dm = dist_matrices[primary_metric]
pcoa_result = pcoa(dm)
pcoa_df = pcoa_result.samples[["PC1", "PC2"]].copy()
pcoa_df["sample"] = dm.ids
pcoa_df["group"] = [metadata.loc[s, "group"] for s in dm.ids]
pcoa_df.to_csv(snakemake.output.pcoa, sep="\t", index=False)

var_explained = pcoa_result.proportion_explained

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
colors = {"healthy": "#2ca02c", "disease": "#d62728"}
for grp, color in colors.items():
    mask = pcoa_df["group"] == grp
    ax.scatter(pcoa_df.loc[mask, "PC1"], pcoa_df.loc[mask, "PC2"],
               c=color, label=grp, s=100, edgecolors="black", linewidth=0.5, alpha=0.8)
ax.set_xlabel(f"PC1 ({var_explained['PC1']*100:.1f}%)", fontsize=12)
ax.set_ylabel(f"PC2 ({var_explained['PC2']*100:.1f}%)", fontsize=12)
ax.set_title(f"PCoA — {primary_metric.replace('_', ' ').title()} Distance", fontsize=14)
ax.legend(title="Group", fontsize=11)

perm_res = permanova_results[0]
ax.annotate(f"PERMANOVA: F={perm_res['pseudo_F']:.2f}, p={perm_res['p_value']:.3f}",
            xy=(0.02, 0.02), xycoords="axes fraction", fontsize=10, fontstyle="italic")
plt.tight_layout()
fig.savefig(snakemake.output.figure, dpi=300, bbox_inches="tight")
plt.close()
