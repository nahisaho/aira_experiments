"""
alpha_diversity.py — Compute alpha diversity metrics from MetaPhlAn 4 profiles.

Metrics: Shannon, Simpson, observed features, Chao1.
Outputs: TSV table + grouped boxplot (SVG).
"""
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

profiles = pd.read_csv(snakemake.input.profiles, sep="\t", comment="#", index_col=0)
metadata = pd.read_csv(snakemake.input.metadata, sep="\t")
metrics_list = snakemake.params.metrics

# Keep species-level rows
species = profiles[profiles.index.str.contains(r"s__") &
                   ~profiles.index.str.contains(r"t__")]
abundance = species.T  # samples × species

def shannon(x):
    x = x[x > 0]
    p = x / x.sum()
    return -np.sum(p * np.log(p))

def simpson(x):
    x = x[x > 0]
    p = x / x.sum()
    return 1 - np.sum(p ** 2)

def observed(x):
    return (x > 0).sum()

def chao1(x):
    n = x.sum()
    f1 = (x == 1).sum()
    f2 = max((x == 2).sum(), 1)
    return observed(x) + (f1 * (f1 - 1)) / (2 * (f2 + 1))

metric_funcs = {"shannon": shannon, "simpson": simpson,
                "observed_features": observed, "chao1": chao1}

results = []
for sample in abundance.index:
    row = {"sample": sample}
    for m in metrics_list:
        row[m] = metric_funcs[m](abundance.loc[sample].values)
    results.append(row)

df = pd.DataFrame(results)
df = df.merge(metadata[["sample", "group"]], on="sample")
df.to_csv(snakemake.output.table, sep="\t", index=False)

# Boxplot
fig, axes = plt.subplots(1, len(metrics_list), figsize=(4*len(metrics_list), 5))
if len(metrics_list) == 1:
    axes = [axes]
palette = {"healthy": "#2ca02c", "disease": "#d62728"}
for ax, m in zip(axes, metrics_list):
    sns.boxplot(data=df, x="group", y=m, ax=ax, palette=palette, width=0.5)
    sns.stripplot(data=df, x="group", y=m, ax=ax, color="black", size=6, alpha=0.7)
    ax.set_title(m.replace("_", " ").title(), fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel(m)
    # Mann-Whitney U test
    g1 = df[df["group"] == "healthy"][m].dropna()
    g2 = df[df["group"] == "disease"][m].dropna()
    if len(g1) >= 2 and len(g2) >= 2:
        stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")
        ax.annotate(f"p = {pval:.3f}", xy=(0.5, 0.95), xycoords="axes fraction",
                    ha="center", fontsize=10, fontstyle="italic")
plt.tight_layout()
fig.savefig(snakemake.output.figure, dpi=300, bbox_inches="tight")
plt.close()
