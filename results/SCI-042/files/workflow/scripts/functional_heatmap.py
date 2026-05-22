"""
functional_heatmap.py — Heatmap of differentially abundant metabolic pathways.

Displays CLR-transformed pathway abundances with sample and group annotations.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

pathways = pd.read_csv(snakemake.input.pathways, sep="\t", index_col=0)
da = pd.read_csv(snakemake.input.da, sep="\t")
metadata = pd.read_csv(snakemake.input.metadata, sep="\t").set_index("sample")

# Filter to significant features if available, else top variable
sig_features = da[da["significant"] == True]["feature"].tolist() if "significant" in da.columns else []
if len(sig_features) < 5:
    variance = pathways.var(axis=1).nlargest(30)
    plot_features = variance.index.tolist()
else:
    plot_features = sig_features[:30]

# Subset and CLR-transform
sub = pathways.loc[pathways.index.isin(plot_features)]
if sub.empty:
    sub = pathways.head(20)

pseudo = sub.replace(0, 1e-6)
log_data = np.log(pseudo)
clr = log_data.subtract(log_data.mean(axis=0), axis=1)

# Annotation colors
samples_ordered = metadata.sort_values("group").index.tolist()
samples_in_data = [s for s in samples_ordered if s in clr.columns]
clr = clr[samples_in_data]

col_colors = pd.Series(
    {s: "#2ca02c" if metadata.loc[s, "group"] == "healthy" else "#d62728"
     for s in samples_in_data},
    name="Group"
)

g = sns.clustermap(
    clr,
    cmap="RdBu_r",
    center=0,
    col_colors=col_colors,
    figsize=(12, max(8, len(clr) * 0.35)),
    yticklabels=True,
    xticklabels=True,
    col_cluster=False,
    dendrogram_ratio=(0.1, 0.08),
    cbar_pos=(0.02, 0.8, 0.03, 0.15),
    linewidths=0.5,
)
g.ax_heatmap.set_xlabel("Samples", fontsize=11)
g.ax_heatmap.set_ylabel("Pathways", fontsize=11)
g.fig.suptitle("Differentially Abundant Metabolic Pathways (CLR)", y=1.02, fontsize=14)

plt.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
plt.close()
