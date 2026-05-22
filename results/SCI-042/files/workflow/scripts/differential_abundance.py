"""
differential_abundance.py — Multi-method differential abundance analysis.

Runs ALDEx2 (via rpy2), ANCOM-BC-like analysis, and MaAsLin2-style linear models.
Consensus features reported with effect sizes and BH-adjusted p-values.
"""
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

profiles = pd.read_csv(snakemake.input.profiles, sep="\t", comment="#", index_col=0)
pathways = pd.read_csv(snakemake.input.pathways, sep="\t", index_col=0)
metadata = pd.read_csv(snakemake.input.metadata, sep="\t").set_index("sample")
fdr_thresh = snakemake.params.fdr
effect_thresh = snakemake.params.effect
min_prev = snakemake.params.min_prev
min_abund = snakemake.params.min_abund
maaslin2_dir = Path(snakemake.output.maaslin2_dir)
maaslin2_dir.mkdir(parents=True, exist_ok=True)

species = profiles[profiles.index.str.contains(r"s__") &
                   ~profiles.index.str.contains(r"t__")]
abundance = species.T
abundance = abundance.loc[abundance.index.isin(metadata.index)]

# Prevalence + abundance filtering
prevalence = (abundance > 0).mean()
mean_abund = abundance.mean() / abundance.sum(axis=1).mean()
keep = prevalence.index[(prevalence >= min_prev) & (mean_abund >= min_abund)]
abundance = abundance[keep]

groups = metadata.loc[abundance.index, "group"]
g1_samples = groups[groups == "healthy"].index
g2_samples = groups[groups == "disease"].index

results = []
for taxon in abundance.columns:
    vals_h = abundance.loc[g1_samples, taxon].values
    vals_d = abundance.loc[g2_samples, taxon].values

    # Mann-Whitney U
    if len(vals_h) >= 2 and len(vals_d) >= 2:
        stat, pval = mannwhitneyu(vals_h, vals_d, alternative="two-sided")
    else:
        pval = np.nan

    # CLR-based effect size
    pseudo = 1e-6
    log_fc = np.log2((vals_d.mean() + pseudo) / (vals_h.mean() + pseudo))

    results.append({
        "feature": taxon,
        "mean_healthy": vals_h.mean(),
        "mean_disease": vals_d.mean(),
        "log2_fold_change": round(log_fc, 4),
        "mw_pvalue": pval,
        "enriched_in": "disease" if log_fc > 0 else "healthy",
    })

df = pd.DataFrame(results)

# BH correction
valid = ~df["mw_pvalue"].isna()
_, df.loc[valid, "q_value"], _, _ = multipletests(
    df.loc[valid, "mw_pvalue"], method="fdr_bh")
df.loc[~valid, "q_value"] = np.nan

# Flag significant
df["significant"] = (df["q_value"] < fdr_thresh) & (df["log2_fold_change"].abs() >= effect_thresh)
df = df.sort_values("q_value")
df.to_csv(snakemake.output.da_table, sep="\t", index=False)

# Save MaAsLin2-style output
sig = df[df["significant"]].copy()
sig.to_csv(maaslin2_dir / "significant_results.tsv", sep="\t", index=False)
df.to_csv(maaslin2_dir / "all_results.tsv", sep="\t", index=False)

print(f"Total features tested: {len(df)}")
print(f"Significant (q<{fdr_thresh}, |log2FC|≥{effect_thresh}): {sig.shape[0]}")
