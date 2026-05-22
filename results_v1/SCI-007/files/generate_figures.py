"""Generate all figures for the antibody design experiment."""
import os, sys, json, random
import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import Counter

random.seed(42); np.random.seed(42); torch.manual_seed(42)
sys.path.insert(0, "/app/projects/031250d9-fdbc-4fbc-8aec-563fa17e5354/workspace")
from antibody_model import AMINO_ACIDS

BASE    = "/app/projects/031250d9-fdbc-4fbc-8aec-563fa17e5354/workspace"
FIGURES = BASE + "/figures"
RESULTS = BASE + "/results"

with open(RESULTS + "/training_history.json") as f:       history   = json.load(f)
with open(RESULTS + "/pdl1_top_candidates.json") as f:    top_cands = json.load(f)
with open(RESULTS + "/pdl1_summary_statistics.json") as f: summary  = json.load(f)
with open(RESULTS + "/optimization_history.json") as f:   opt_hist  = json.load(f)
df = pd.read_csv(RESULTS + "/pdl1_candidate_table.csv")
df["type"] = df["is_benchmark"].map({True: "Benchmark", False: "Generated"})

sns.set_theme(style="whitegrid", palette="colorblind")
p = sns.color_palette("viridis", 6)

# ── Figure 1: Training Curves ──────────────────────────────────────
epochs = range(1, len(history["train_loss"]) + 1)
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0,0].plot(epochs, history["train_loss"], color=p[0], lw=2)
axes[0,0].set_xlabel("Epoch"); axes[0,0].set_ylabel("Loss")
axes[0,0].set_title("Total Training Loss")

axes[0,1].plot(epochs, history["val_kd_pearson"], color=p[1], lw=2, label="Pearson r")
axes[0,1].plot(epochs, history["val_kd_spearman"], color=p[2], lw=2, ls="--", label="Spearman ρ")
axes[0,1].axhline(0, color="gray", lw=0.8, ls=":")
axes[0,1].set_xlabel("Epoch"); axes[0,1].set_ylabel("Correlation")
axes[0,1].set_title("Binding Affinity (log Kd) Prediction"); axes[0,1].set_ylim(-0.3, 0.5); axes[0,1].legend()

axes[1,0].plot(epochs, history["val_kd_rmse"], color=p[3], lw=2)
axes[1,0].set_xlabel("Epoch"); axes[1,0].set_ylabel("RMSE")
axes[1,0].set_title("Log Kd RMSE")

axes[1,1].plot(epochs, history["val_tm_pearson"], color=p[4], lw=2, label="Pearson r")
axes[1,1].plot(epochs, history["val_tm_rmse"], color=p[5], lw=2, ls="--", label="RMSE")
axes[1,1].set_xlabel("Epoch"); axes[1,1].set_ylabel("Value")
axes[1,1].set_title("Tm Prediction"); axes[1,1].legend()

fig.suptitle("Antibody Design Model — Training Curves", fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGURES + "/fig1_training_curves.png", dpi=150, bbox_inches="tight"); plt.close()
print("fig1 saved")

# ── Figure 2: Property Distributions ──────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
props = [
    ("log_kd", "Predicted log Kd", False),
    ("tm", "Predicted Tm (arb)", True),
    ("humanization_score", "Humanization Score", True),
    ("immunogenicity_risk", "Immunogenicity Risk", False),
    ("aggregation_score", "Aggregation Score", False),
    ("developability_index", "Developability Index", True),
]
for ax, (col, title, hi_better) in zip(axes.flat, props):
    sns.histplot(data=df, x=col, hue="type", ax=ax, bins=15,
                 palette={"Benchmark": p[0], "Generated": p[3]},
                 alpha=0.7, stat="density")
    ax.set_title(title); ax.set_xlabel("")
    arrow = "↑ better" if hi_better else "↓ better"
    ax.text(0.97, 0.95, arrow, transform=ax.transAxes, ha="right", va="top", fontsize=9, color="gray")
fig.suptitle("Candidate Property Distributions: Generated vs Benchmark", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGURES + "/fig2_property_distributions.png", dpi=150, bbox_inches="tight"); plt.close()
print("fig2 saved")

# ── Figure 3: Multi-Objective Scatter ─────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
gen_df = df[df["type"] == "Generated"]
ben_df = df[df["type"] == "Benchmark"]

def scatter_2d(ax, xc, yc, xl, yl, title):
    sc = ax.scatter(gen_df[xc], gen_df[yc], c=gen_df["developability_index"],
                    cmap="viridis", alpha=0.7, s=60, label="Generated", zorder=3)
    if len(ben_df) > 0:
        ax.scatter(ben_df[xc], ben_df[yc], marker="*", s=200, color="red", zorder=5, label="Benchmark")
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title)
    plt.colorbar(sc, ax=ax, label="Developability"); ax.legend(fontsize=8)

scatter_2d(axes[0], "log_kd", "humanization_score", "Predicted log Kd", "Humanization Score", "Affinity vs Humanization")
scatter_2d(axes[1], "aggregation_score", "tm", "Aggregation Score", "Predicted Tm", "Aggregation vs Stability")
scatter_2d(axes[2], "immunogenicity_risk", "developability_index", "Immunogenicity Risk", "Developability Index", "Immunogenicity vs Developability")
fig.suptitle("Multi-Objective Property Space — PD-L1 CDR-H3 Candidates", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGURES + "/fig3_multi_objective_scatter.png", dpi=150, bbox_inches="tight"); plt.close()
print("fig3 saved")

# ── Figure 4: Optimization Convergence ────────────────────────────
oh = opt_hist["best_scores_history"]
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(range(1, len(oh)+1), oh, color=p[0], lw=2)
ax.fill_between(range(1, len(oh)+1), oh, alpha=0.2, color=p[0])
ax.set_xlabel("Generation"); ax.set_ylabel("Best Composite Score")
ax.set_title("Genetic Algorithm Convergence (PD-L1 CDR-H3)")
ax.grid(True, alpha=0.4)
fig.tight_layout()
fig.savefig(FIGURES + "/fig4_optimization_convergence.png", dpi=150, bbox_inches="tight"); plt.close()
print("fig4 saved")

# ── Figure 5: Length Distribution + Top-10 Heatmap ────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
lengths = df["length"].tolist()
axes[0].hist(lengths, bins=range(5, 26), color=p[2], edgecolor="white", alpha=0.8)
axes[0].axvline(np.mean(lengths), color="red", ls="--", label=f"Mean={np.mean(lengths):.1f}")
axes[0].set_xlabel("CDR-H3 Length (AA)"); axes[0].set_ylabel("Count")
axes[0].set_title("CDR-H3 Length Distribution"); axes[0].legend()

top10_df = pd.DataFrame([{
    "Seq": c["label"][:14],
    "log_Kd": c["log_kd"],
    "Tm": c["tm"],
    "Human.": c["humanization_score"],
    "Immuno.": c["immunogenicity_risk"],
    "Agg.": c["aggregation_score"],
    "Dev.": c["developability_index"],
} for c in top_cands[:10]]).set_index("Seq")

sns.heatmap(top10_df, ax=axes[1], cmap="viridis", annot=True, fmt=".2f",
            linewidths=0.5, annot_kws={"size": 8})
axes[1].set_title("Top-10 PD-L1 Candidates Property Heatmap")
axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=8, rotation=0)
fig.tight_layout()
fig.savefig(FIGURES + "/fig5_cdrh3_analysis.png", dpi=150, bbox_inches="tight"); plt.close()
print("fig5 saved")

# ── Figure 6: Amino Acid Composition ──────────────────────────────
gen_seqs = df[df["type"] == "Generated"]["sequence"].tolist()
all_aas = "".join(gen_seqs)
counts = Counter(all_aas)
aa_order = sorted(AMINO_ACIDS)
freqs = [counts.get(aa, 0) / max(len(all_aas), 1) for aa in aa_order]

fig, ax = plt.subplots(figsize=(10, 5))
x = np.arange(len(aa_order)); w = 0.35
ax.bar(x - w/2, freqs, w, label="Generated CDRs", color=p[0], alpha=0.8)
ax.bar(x + w/2, [1/20]*20, w, label="Uniform Background", color=p[3], alpha=0.6)
ax.set_xticks(x); ax.set_xticklabels(aa_order)
ax.set_xlabel("Amino Acid"); ax.set_ylabel("Frequency")
ax.set_title("Amino Acid Composition of Generated CDR-H3 Sequences")
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES + "/fig6_aa_composition.png", dpi=150, bbox_inches="tight"); plt.close()
print("fig6 saved")

print("\nAll figures saved to", FIGURES)
print("Files:", sorted(os.listdir(FIGURES)))
