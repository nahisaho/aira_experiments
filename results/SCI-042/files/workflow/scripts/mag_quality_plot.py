"""
mag_quality_plot.py — Scatter plot of MAG completeness vs contamination.

Color-coded by MIMAG quality tier; size proportional to genome size.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

report = pd.read_csv(snakemake.input.report, sep="\t")

fig, ax = plt.subplots(figsize=(9, 7))

# MIMAG quality tiers
def classify_mag(row):
    if row["Completeness"] >= 90 and row["Contamination"] <= 5:
        return "High Quality"
    elif row["Completeness"] >= 50 and row["Contamination"] <= 10:
        return "Medium Quality"
    else:
        return "Low Quality"

report["quality_tier"] = report.apply(classify_mag, axis=1)

colors = {"High Quality": "#2ca02c", "Medium Quality": "#ff7f0e", "Low Quality": "#d62728"}
for tier, color in colors.items():
    mask = report["quality_tier"] == tier
    sub = report[mask]
    size = sub.get("Genome_Size", pd.Series([50]*len(sub)))
    ax.scatter(sub["Completeness"], sub["Contamination"],
               c=color, label=f"{tier} (n={len(sub)})",
               s=80, alpha=0.7, edgecolors="black", linewidth=0.5)

# Reference lines for MIMAG thresholds
ax.axvline(x=50, color="gray", linestyle="--", alpha=0.5, label="MIMAG 50% comp.")
ax.axvline(x=90, color="gray", linestyle=":", alpha=0.5, label="MIMAG 90% comp.")
ax.axhline(y=5, color="lightblue", linestyle="--", alpha=0.5, label="MIMAG 5% cont.")
ax.axhline(y=10, color="lightcoral", linestyle="--", alpha=0.5, label="MIMAG 10% cont.")

ax.set_xlabel("Completeness (%)", fontsize=12)
ax.set_ylabel("Contamination (%)", fontsize=12)
ax.set_title("MAG Quality Assessment (CheckM2)", fontsize=14)
ax.legend(loc="upper left", fontsize=9)
ax.set_xlim(0, 105)
ax.set_ylim(-1, max(report["Contamination"].max() * 1.1, 15))
ax.invert_yaxis()

plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
plt.close()
