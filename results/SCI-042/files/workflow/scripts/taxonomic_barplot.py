"""
taxonomic_barplot.py — Stacked barplot of taxonomic composition.

Top 15 taxa at phylum & genus level; remaining grouped as "Other".
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

profiles = pd.read_csv(snakemake.input.profiles, sep="\t", comment="#", index_col=0)
metadata = pd.read_csv(snakemake.input.metadata, sep="\t")

def filter_level(df, level_prefix, exclude_prefix=None):
    mask = df.index.str.contains(f"{level_prefix}__")
    if exclude_prefix:
        mask = mask & ~df.index.str.contains(f"{exclude_prefix}__")
    sub = df[mask].copy()
    sub.index = sub.index.map(lambda x: x.split("|")[-1].replace(f"{level_prefix}__", ""))
    return sub

top_n = 15
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

for ax, (level, prefix, excl) in zip(axes, [("Phylum", "p", "c"), ("Genus", "g", "s")]):
    data = filter_level(profiles, prefix, excl).T
    data = data.div(data.sum(axis=1), axis=0) * 100

    top_taxa = data.mean().nlargest(top_n).index.tolist()
    other = data.drop(columns=top_taxa, errors="ignore").sum(axis=1)
    plot_df = data[top_taxa].copy()
    plot_df["Other"] = other

    # Sort by group
    order = metadata.sort_values("group")["sample"].tolist()
    plot_df = plot_df.reindex([s for s in order if s in plot_df.index])

    plot_df.plot.bar(stacked=True, ax=ax, width=0.85,
                     colormap="tab20", edgecolor="none")
    ax.set_title(f"{level}-level Composition", fontsize=13)
    ax.set_ylabel("Relative Abundance (%)")
    ax.set_xlabel("")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, ncol=1)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)

plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
plt.close()
