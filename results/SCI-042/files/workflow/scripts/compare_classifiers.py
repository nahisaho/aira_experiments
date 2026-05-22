"""
compare_classifiers.py — Concordance analysis: Kraken2/Bracken vs MetaPhlAn 4

Computes Spearman correlation and Jaccard index of detected species between
two classification approaches per sample.
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

bracken_files = snakemake.input.bracken
metaphlan_files = snakemake.input.metaphlan
output_file = snakemake.output[0]

results = []

for bf, mf in zip(sorted(bracken_files), sorted(metaphlan_files)):
    sample = Path(bf).stem.replace(".bracken", "")

    # Parse Bracken
    bdf = pd.read_csv(bf, sep="\t")
    bracken_species = dict(zip(
        bdf["name"].str.strip(),
        bdf["fraction_total_reads"]
    ))

    # Parse MetaPhlAn 4 (species-level rows only)
    mdf = pd.read_csv(mf, sep="\t", comment="#")
    mdf.columns = ["clade_name", "clade_taxid", "relative_abundance",
                    "additional_species"][:len(mdf.columns)]
    species_rows = mdf[mdf["clade_name"].str.contains(r"s__") &
                       ~mdf["clade_name"].str.contains(r"t__")]
    metaphlan_species = {}
    for _, row in species_rows.iterrows():
        sp_name = row["clade_name"].split("|")[-1].replace("s__", "").replace("_", " ")
        metaphlan_species[sp_name] = float(row["relative_abundance"])

    # Jaccard index of detected species
    set_b = set(bracken_species.keys())
    set_m = set(metaphlan_species.keys())
    jaccard = len(set_b & set_m) / len(set_b | set_m) if (set_b | set_m) else 0

    # Spearman correlation on shared species
    shared = sorted(set_b & set_m)
    if len(shared) >= 3:
        vals_b = [bracken_species[s] for s in shared]
        vals_m = [metaphlan_species[s] for s in shared]
        rho, pval = spearmanr(vals_b, vals_m)
    else:
        rho, pval = np.nan, np.nan

    results.append({
        "sample": sample,
        "kraken2_species_count": len(set_b),
        "metaphlan4_species_count": len(set_m),
        "shared_species_count": len(shared),
        "jaccard_index": round(jaccard, 4),
        "spearman_rho": round(rho, 4) if not np.isnan(rho) else np.nan,
        "spearman_pvalue": pval,
    })

pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)
