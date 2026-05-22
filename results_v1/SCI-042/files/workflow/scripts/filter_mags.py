"""
filter_mags.py — Filter MAGs by MIMAG quality standards

Reads CheckM2 quality report, applies completeness/contamination thresholds,
classifies MAGs as high-quality or medium-quality, and copies passing bins.
"""
import pandas as pd
import shutil
from pathlib import Path

report = pd.read_csv(snakemake.input.report, sep="\t")
bins_dir = Path(snakemake.input.bins)
out_dir = Path(snakemake.output[0])
summary_file = snakemake.output[1]

comp_thresh = snakemake.params.comp
cont_thresh = snakemake.params.cont
hq_comp = snakemake.params.hq_comp
hq_cont = snakemake.params.hq_cont

out_dir.mkdir(parents=True, exist_ok=True)

# Filter
passing = report[
    (report["Completeness"] >= comp_thresh) &
    (report["Contamination"] <= cont_thresh)
].copy()

# Quality tier
passing["quality_tier"] = passing.apply(
    lambda r: "high" if r["Completeness"] >= hq_comp and r["Contamination"] <= hq_cont
              else "medium",
    axis=1
)

# Copy passing bins
for _, row in passing.iterrows():
    bin_name = row["Name"]
    src = bins_dir / f"{bin_name}.fa"
    if src.exists():
        shutil.copy2(src, out_dir / f"{bin_name}.fa")

passing.to_csv(summary_file, sep="\t", index=False)
print(f"Filtered MAGs: {len(passing)} passed "
      f"(comp≥{comp_thresh}%, cont≤{cont_thresh}%)")
print(f"  High-quality: {(passing['quality_tier']=='high').sum()}")
print(f"  Medium-quality: {(passing['quality_tier']=='medium').sum()}")
