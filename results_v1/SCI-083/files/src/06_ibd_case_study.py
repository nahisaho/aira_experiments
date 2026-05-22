#!/usr/bin/env python3
"""
Module 6: IBD ケーススタディ — 統合解析パイプライン
Inflammatory Bowel Disease (IBD) Case Study

Integrates all modules:
  1. Peak annotation (Module 1)
  2. Correlation network (Module 2)
  3. Causal inference (Module 3)
  4. Pathway enrichment (Module 4)
  5. Biomarker scoring (Module 5)

Cohort: Control (n=50), UC (n=50), CD (n=50)
Omics: 16S rRNA + untargeted LC-MS/MS metabolomics
"""

import os
import sys
import json
import logging
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Resolve imports from src
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# mixOmics Integration Pipeline (R script generator)
# ---------------------------------------------------------------------------

def generate_mixomics_r_script(output_dir: str) -> str:
    """
    mixOmics R パイプライン生成
    DIABLO (Data Integration Analysis for Biomarker discovery using Latent cOmponents)
    """
    r_script = f"""
# ============================================================
# mixOmics DIABLO Integration Pipeline for IBD Case Study
# ============================================================
library(mixOmics)

# --- Data loading ---
taxa <- read.csv("{output_dir}/taxa_abundance.csv", row.names=1)
metabolites <- read.csv("{output_dir}/metabolite_abundance.csv", row.names=1)
metadata <- read.csv("{output_dir}/metadata.csv")

Y <- factor(metadata$group, levels = c("Control", "UC", "CD"))

# --- Pre-processing ---
# CLR transform for taxa
taxa.clr <- logratio.transfo(as.matrix(taxa), logratio = "CLR", offset = 1)

# Log transform for metabolites
metabolites.log <- log2(as.matrix(metabolites) + 1)

# --- Design matrix ---
# Specify expected correlation between blocks
design <- matrix(c(0,   1,   0.1,
                   1,   0,   0.1,
                   0.1, 0.1, 0  ), ncol=3, nrow=3,
                 dimnames = list(c("taxa", "metabolites", "clinical"),
                                 c("taxa", "metabolites", "clinical")))

# --- DIABLO model (block.splsda) ---
X <- list(taxa = taxa.clr, metabolites = metabolites.log)

# Tune keepX (number of features per component)
# tune.diablo <- tune.block.splsda(
#   X, Y, ncomp = 2,
#   test.keepX = list(taxa = c(5, 10, 15), metabolites = c(10, 15, 20)),
#   design = design[1:2, 1:2],
#   validation = "Mfold", folds = 5, nrepeat = 10
# )

# Fit final model
diablo.model <- block.splsda(
  X, Y, ncomp = 2,
  keepX = list(taxa = c(10, 10), metabolites = c(15, 15)),
  design = design[1:2, 1:2]
)

# --- Performance evaluation ---
perf.diablo <- perf(diablo.model, validation = "Mfold", folds = 5,
                     nrepeat = 10, progressBar = TRUE)

# --- Visualization ---
pdf("{output_dir}/../figures/diablo_plotIndiv.pdf")
plotIndiv(diablo.model, ind.names = FALSE, legend = TRUE,
          title = "DIABLO Sample Plot (IBD)")
dev.off()

pdf("{output_dir}/../figures/diablo_circosPlot.pdf")
circosPlot(diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks = c("steelblue", "darkorange"),
           color.cor = c("red", "blue"))
dev.off()

pdf("{output_dir}/../figures/diablo_loadings.pdf")
plotLoadings(diablo.model, comp = 1, contrib = "max",
             method = "median", legend.color = c("blue", "red", "green"))
dev.off()

pdf("{output_dir}/../figures/diablo_network.pdf")
network(diablo.model, blocks = c(1, 2), cutoff = 0.4,
        color.node = c("steelblue", "darkorange"))
dev.off()

# --- Export selected features ---
selected.taxa <- selectVar(diablo.model, block = "taxa", comp = 1)$taxa$name
selected.met <- selectVar(diablo.model, block = "metabolites", comp = 1)$metabolites$name

write.csv(data.frame(feature = selected.taxa), "{output_dir}/diablo_selected_taxa.csv")
write.csv(data.frame(feature = selected.met), "{output_dir}/diablo_selected_metabolites.csv")

# --- Save performance ---
sink("{output_dir}/diablo_performance.txt")
print(perf.diablo)
sink()

cat("DIABLO pipeline completed successfully.\\n")
"""
    return r_script


# ---------------------------------------------------------------------------
# MelonnPan Pipeline (R script generator)
# ---------------------------------------------------------------------------

def generate_melonnpan_r_script(output_dir: str) -> str:
    """
    MelonnPan: Metagenomic prediction of community metabolomes
    菌叢データから代謝物プロファイルを予測
    """
    r_script = f"""
# ============================================================
# MelonnPan: Metabolite Prediction from Microbiome Data
# ============================================================
library(melonnpan)

# --- Data loading ---
taxa <- read.csv("{output_dir}/taxa_abundance.csv", row.names=1)
metabolites <- read.csv("{output_dir}/metabolite_abundance.csv", row.names=1)

# --- MelonnPan training ---
# Train the model on paired microbiome-metabolome data
melonnpan.train(
  metab = metabolites,
  taxa = taxa,
  output = "{output_dir}/melonnpan_output"
)

# --- Prediction on new samples ---
melonnpan.predict(
  taxa = taxa,
  weight.matrix = "{output_dir}/melonnpan_output/MelonnPan_Trained_Weights.txt",
  output = "{output_dir}/melonnpan_predictions"
)

# --- Evaluation ---
# Compare predicted vs measured metabolites
predicted <- read.table(
  "{output_dir}/melonnpan_predictions/MelonnPan_Predicted_Metabolites.txt",
  header=TRUE, row.names=1
)

# Correlation analysis
cors <- sapply(1:ncol(metabolites), function(i) {{
  if(colnames(metabolites)[i] %in% colnames(predicted)) {{
    cor(metabolites[,i], predicted[,colnames(metabolites)[i]], method="spearman")
  }} else {{ NA }}
}})

cor.df <- data.frame(
  metabolite = colnames(metabolites),
  spearman_r = cors,
  well_predicted = cors > 0.3
)
write.csv(cor.df, "{output_dir}/melonnpan_prediction_quality.csv")

cat("MelonnPan pipeline completed.\\n")
"""
    return r_script


# ---------------------------------------------------------------------------
# IBD-specific differential analysis
# ---------------------------------------------------------------------------

def ibd_differential_analysis(taxa_df: pd.DataFrame, met_df: pd.DataFrame,
                               metadata: pd.DataFrame) -> dict:
    """IBD 特異的差分解析 (Control vs UC, Control vs CD, UC vs CD)"""
    logger.info("Running IBD differential analysis")

    comparisons = [
        ("Control", "UC"),
        ("Control", "CD"),
        ("UC", "CD"),
    ]

    results = {}
    for grp1, grp2 in comparisons:
        comp_name = f"{grp1}_vs_{grp2}"
        mask1 = metadata["group"] == grp1
        mask2 = metadata["group"] == grp2

        taxa_results = []
        taxa_num = taxa_df.select_dtypes(include=[np.number])
        for col in taxa_num.columns[:20]:
            stat, pval = stats.mannwhitneyu(
                taxa_num.loc[mask1, col],
                taxa_num.loc[mask2, col],
                alternative="two-sided"
            )
            fc = np.median(taxa_num.loc[mask2, col]) / (np.median(taxa_num.loc[mask1, col]) + 1e-10)
            taxa_results.append({
                "feature": col, "omic": "taxa",
                "pvalue": pval, "fold_change": round(fc, 3),
            })

        met_results = []
        met_num = met_df.select_dtypes(include=[np.number])
        for col in met_num.columns[:30]:
            stat, pval = stats.mannwhitneyu(
                met_num.loc[mask1, col],
                met_num.loc[mask2, col],
                alternative="two-sided"
            )
            fc = np.median(met_num.loc[mask2, col]) / (np.median(met_num.loc[mask1, col]) + 1e-10)
            met_results.append({
                "feature": col, "omic": "metabolome",
                "pvalue": pval, "fold_change": round(fc, 3),
            })

        all_res = pd.DataFrame(taxa_results + met_results)
        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(all_res["pvalue"], method="fdr_bh")
        all_res["fdr_qvalue"] = fdr

        sig = all_res[all_res["fdr_qvalue"] < 0.05]
        results[comp_name] = {
            "n_tested": len(all_res),
            "n_significant": len(sig),
            "top_features": sig.sort_values("pvalue").head(10).to_dict("records"),
        }

    return results


# ---------------------------------------------------------------------------
# IBD Disease Activity Score
# ---------------------------------------------------------------------------

def compute_ibd_activity_score(taxa_df: pd.DataFrame, met_df: pd.DataFrame,
                                metadata: pd.DataFrame) -> pd.DataFrame:
    """
    IBD activity composite score:
      Score = w1 * (Faecalibacterium depletion) + w2 * (butyrate depletion)
            + w3 * (E.coli enrichment) + w4 * (calprotectin proxy)
    """
    logger.info("Computing IBD activity composite score")

    taxa_num = taxa_df.select_dtypes(include=[np.number])
    met_num = met_df.select_dtypes(include=[np.number])

    n = len(metadata)

    # Standardize key features
    def zscore(x):
        return (x - x.mean()) / (x.std() + 1e-10)

    # Use first few columns as proxies
    faecal_z = -zscore(taxa_num.iloc[:, 1])  # Faecalibacterium depletion
    ecoli_z = zscore(taxa_num.iloc[:, 14] if taxa_num.shape[1] > 14 else taxa_num.iloc[:, -1])
    butyrate_z = -zscore(met_num.iloc[:, 0])  # butyrate depletion
    tmao_z = zscore(met_num.iloc[:, 5] if met_num.shape[1] > 5 else met_num.iloc[:, -1])

    composite = 0.3 * faecal_z + 0.2 * ecoli_z + 0.3 * butyrate_z + 0.2 * tmao_z

    scores_df = metadata[["sample_id", "group"]].copy()
    scores_df["ibd_activity_score"] = composite.values.round(3)
    scores_df["risk_category"] = pd.cut(
        composite, bins=[-np.inf, -0.5, 0.5, np.inf],
        labels=["Low", "Moderate", "High"]
    )

    return scores_df


# ---------------------------------------------------------------------------
# Main IBD case study pipeline
# ---------------------------------------------------------------------------

def run_ibd_case_study(output_dir: str = "results") -> dict:
    os.makedirs(output_dir, exist_ok=True)
    figures_dir = os.path.join(os.path.dirname(output_dir), "figures")
    os.makedirs(figures_dir, exist_ok=True)

    logger.info("=" * 60)
    logger.info("IBD Case Study: Integrated Multi-omic Analysis")
    logger.info("=" * 60)

    # ---- Module 2: Generate data & correlation network ----
    from importlib import import_module
    try:
        mod2 = import_module("02_correlation_network")
        data = mod2.generate_simulated_data()
    except ImportError:
        logger.info("Generating data locally")
        np.random.seed(42)
        n = 150
        groups = ["Control"]*50 + ["UC"]*50 + ["CD"]*50
        metadata = pd.DataFrame({
            "sample_id": [f"S{i:03d}" for i in range(n)],
            "group": groups,
            "age": np.random.normal(45, 12, n).astype(int).clip(18, 80),
            "sex": np.random.choice(["M", "F"], n),
            "bmi": np.random.normal(25, 4, n).round(1).clip(16, 45),
        })

        taxa_data = np.random.lognormal(2, 1.5, (n, 80))
        taxa_names = [f"g__genus_{i}" for i in range(80)]
        for i in range(50, 150):
            taxa_data[i, 1] *= 0.3
            taxa_data[i, 2] *= 0.4
            taxa_data[i, 14] *= 2.5

        taxa_df = pd.DataFrame(taxa_data, columns=taxa_names)
        taxa_df.insert(0, "sample_id", metadata["sample_id"])

        met_data = np.random.lognormal(8, 2, (n, 200))
        met_names = [f"met_{i:04d}" for i in range(200)]
        for i in range(50, 150):
            met_data[i, 0] *= 0.3
            met_data[i, 5] *= 2.0
            met_data[i, 10] *= 0.4

        met_df = pd.DataFrame(met_data, columns=met_names)
        met_df.insert(0, "sample_id", metadata["sample_id"])
        data = {"metadata": metadata, "taxa": taxa_df, "metabolites": met_df}

    metadata = data["metadata"]
    taxa_df = data["taxa"]
    met_df = data["metabolites"]

    metadata.to_csv(os.path.join(output_dir, "metadata.csv"), index=False)
    taxa_df.to_csv(os.path.join(output_dir, "taxa_abundance.csv"), index=False)
    met_df.to_csv(os.path.join(output_dir, "metabolite_abundance.csv"), index=False)

    # ---- IBD Differential Analysis ----
    diff_results = ibd_differential_analysis(taxa_df, met_df, metadata)

    # ---- IBD Activity Score ----
    activity_scores = compute_ibd_activity_score(taxa_df, met_df, metadata)
    activity_scores.to_csv(os.path.join(output_dir, "ibd_activity_scores.csv"), index=False)

    # Score statistics by group
    score_stats = activity_scores.groupby("group")["ibd_activity_score"].agg(
        ["mean", "std", "median"]
    ).round(3).to_dict()

    # ---- Generate R scripts ----
    mixomics_script = generate_mixomics_r_script(output_dir)
    melonnpan_script = generate_melonnpan_r_script(output_dir)

    with open(os.path.join(output_dir, "run_mixomics_diablo.R"), "w") as f:
        f.write(mixomics_script)
    with open(os.path.join(output_dir, "run_melonnpan.R"), "w") as f:
        f.write(melonnpan_script)

    # ---- Summary ----
    summary = {
        "cohort": {
            "total_samples": len(metadata),
            "groups": metadata["group"].value_counts().to_dict(),
            "mean_age": round(metadata["age"].mean(), 1),
            "sex_ratio": metadata["sex"].value_counts().to_dict(),
        },
        "differential_analysis": {
            comp: {
                "n_significant": res["n_significant"],
                "n_tested": res["n_tested"],
            }
            for comp, res in diff_results.items()
        },
        "activity_score_stats": score_stats,
        "risk_distribution": activity_scores["risk_category"].value_counts().to_dict(),
    }

    with open(os.path.join(output_dir, "ibd_case_study_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    with open(os.path.join(output_dir, "differential_results.json"), "w") as f:
        json.dump(diff_results, f, indent=2, default=str)

    logger.info("IBD Case Study pipeline completed")
    return summary


if __name__ == "__main__":
    summary = run_ibd_case_study(output_dir="../results")
    print(json.dumps(summary, indent=2, default=str))
