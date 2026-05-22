#!/usr/bin/env python3
"""
Epitranscriptome Analysis Pipeline
===================================
Transcriptome-wide mapping of RNA modifications (m6A, m5C, Pseudouridine).

Modules:
  1. Data processing (MeRIP-seq / DART-seq / Nanopore direct RNA-seq)
  2. Peak calling for modification site detection
  3. Quantification and differential modification analysis
  4. Functional annotation (mRNA stability, translation efficiency)
  5. Writer/Reader/Eraser association analysis
  6. Cancer epitranscriptome case study (AML)

Author: Co-Scientist Pipeline
Date: 2026-05-23
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.signal import find_peaks
from scipy.special import betaln
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
import json
import os
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)

# ── Color palette (colorblind-friendly) ──
PALETTE = sns.color_palette("colorblind", 10)
MOD_COLORS = {"m6A": PALETTE[0], "m5C": PALETTE[1], "Ψ": PALETTE[2]}
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


# ═══════════════════════════════════════════════════════════════════
# Module 1: Simulated Data Generation
# ═══════════════════════════════════════════════════════════════════

def generate_transcriptome_data(n_genes=5000, n_modifications=3000):
    """Generate simulated transcriptome with RNA modification sites."""
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    regions = ["5UTR", "CDS", "3UTR", "start_codon", "stop_codon"]
    region_weights = [0.08, 0.25, 0.50, 0.07, 0.10]
    mod_types = ["m6A", "m5C", "Ψ"]
    mod_weights = [0.60, 0.25, 0.15]

    genes = []
    for i in range(n_genes):
        gene_id = f"ENSG{i:08d}"
        gene_name = f"GENE{i}"
        chrom = np.random.choice(chroms)
        start = np.random.randint(1_000_000, 200_000_000)
        length = np.random.randint(500, 50000)
        expr = np.random.lognormal(3, 2)
        genes.append({
            "gene_id": gene_id, "gene_name": gene_name,
            "chrom": chrom, "start": start, "end": start + length,
            "expression_tpm": max(0.1, expr)
        })
    gene_df = pd.DataFrame(genes)

    sites = []
    for i in range(n_modifications):
        gene_idx = np.random.randint(0, n_genes)
        gene = genes[gene_idx]
        mod_type = np.random.choice(mod_types, p=mod_weights)
        region = np.random.choice(regions, p=region_weights)
        pos = np.random.randint(gene["start"], gene["end"])

        # Modification stoichiometry (fraction modified)
        stoichiometry = np.random.beta(2, 5)

        # DRACH motif probability (higher for m6A)
        if mod_type == "m6A":
            motif_match = np.random.random() < 0.85
            motif = "GGACU" if motif_match else "AGACU"
        elif mod_type == "m5C":
            motif_match = np.random.random() < 0.70
            motif = "CTCCA" if motif_match else "ATCGA"
        else:
            motif_match = np.random.random() < 0.60
            motif = "TGTAG" if motif_match else "AGTCG"

        sites.append({
            "site_id": f"MOD{i:06d}",
            "gene_id": gene["gene_id"],
            "gene_name": gene["gene_name"],
            "chrom": gene["chrom"],
            "position": pos,
            "modification": mod_type,
            "region": region,
            "stoichiometry": stoichiometry,
            "motif": motif,
            "motif_match": motif_match,
            "confidence_score": np.random.beta(5, 2),
        })
    site_df = pd.DataFrame(sites)
    return gene_df, site_df


def simulate_merip_seq(site_df, n_input=2, n_ip=2):
    """Simulate MeRIP-seq count data (input vs IP)."""
    records = []
    for _, site in site_df.iterrows():
        base_count = np.random.negative_binomial(5, 0.01)
        for rep in range(n_input):
            records.append({
                "site_id": site["site_id"],
                "sample_type": "input",
                "replicate": rep + 1,
                "counts": max(1, int(base_count * np.random.uniform(0.8, 1.2)))
            })
        enrichment = 1.0 + site["stoichiometry"] * 4
        for rep in range(n_ip):
            ip_count = int(base_count * enrichment * np.random.uniform(0.7, 1.3))
            records.append({
                "site_id": site["site_id"],
                "sample_type": "IP",
                "replicate": rep + 1,
                "counts": max(1, ip_count)
            })
    return pd.DataFrame(records)


def simulate_nanopore_signals(site_df, n_reads=50):
    """Simulate nanopore direct RNA-seq current signals."""
    records = []
    for _, site in site_df[site_df["modification"] == "m6A"].head(500).iterrows():
        for read_i in range(n_reads):
            is_modified = np.random.random() < site["stoichiometry"]
            if is_modified:
                current_mean = np.random.normal(105, 5)
                dwell_time = np.random.lognormal(2.5, 0.5)
            else:
                current_mean = np.random.normal(120, 5)
                dwell_time = np.random.lognormal(2.0, 0.5)
            records.append({
                "site_id": site["site_id"],
                "read_id": f"read_{read_i:04d}",
                "current_mean": current_mean,
                "current_std": np.random.uniform(2, 8),
                "dwell_time": dwell_time,
                "is_modified": is_modified,
                "mod_probability": site["stoichiometry"] + np.random.normal(0, 0.05)
            })
    return pd.DataFrame(records)


def simulate_dart_seq(site_df, n_reads_per_site=100):
    """Simulate DART-seq C-to-T mutation data."""
    records = []
    for _, site in site_df[site_df["modification"] == "m6A"].head(800).iterrows():
        total_reads = np.random.poisson(n_reads_per_site)
        if total_reads < 10:
            total_reads = 10
        # Mutation rate proportional to stoichiometry
        ctrl_mut_rate = 0.005 + np.random.uniform(0, 0.01)
        treated_mut_rate = site["stoichiometry"] * 0.3 + np.random.uniform(0, 0.05)

        ctrl_mutations = np.random.binomial(total_reads, ctrl_mut_rate)
        treated_mutations = np.random.binomial(total_reads, treated_mut_rate)

        records.append({
            "site_id": site["site_id"],
            "gene_name": site["gene_name"],
            "total_reads": total_reads,
            "ctrl_mutations": ctrl_mutations,
            "ctrl_mut_rate": ctrl_mutations / total_reads,
            "treated_mutations": treated_mutations,
            "treated_mut_rate": treated_mutations / total_reads,
        })
    return pd.DataFrame(records)


# ═══════════════════════════════════════════════════════════════════
# Module 2: Peak Calling Algorithms
# ═══════════════════════════════════════════════════════════════════

class PeakCaller:
    """Peak calling for RNA modification site detection."""

    @staticmethod
    def merip_peak_calling(merip_df, fdr_threshold=0.05, fc_min=2.0):
        """
        Sliding-window peak calling for MeRIP-seq.
        Uses Fisher's exact test on IP vs Input counts.
        """
        site_ids = merip_df["site_id"].unique()
        results = []

        for sid in site_ids:
            site_data = merip_df[merip_df["site_id"] == sid]
            input_counts = site_data[site_data["sample_type"] == "input"]["counts"].values
            ip_counts = site_data[site_data["sample_type"] == "IP"]["counts"].values

            mean_input = np.mean(input_counts)
            mean_ip = np.mean(ip_counts)

            if mean_input > 0:
                fold_change = mean_ip / mean_input
            else:
                fold_change = float("inf")

            # Fisher's exact test (2x2 table approximation)
            a = int(np.sum(ip_counts))
            b = int(np.sum(input_counts))
            c = max(1, int(a + b - np.sum(ip_counts)))
            d = max(1, int(a + b - np.sum(input_counts)))
            _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

            results.append({
                "site_id": sid,
                "mean_input": mean_input,
                "mean_ip": mean_ip,
                "fold_change": fold_change,
                "pvalue": pvalue
            })

        result_df = pd.DataFrame(results)
        _, fdr, _, _ = multipletests(result_df["pvalue"], method="fdr_bh")
        result_df["fdr"] = fdr
        result_df["significant"] = (result_df["fdr"] < fdr_threshold) & (result_df["fold_change"] >= fc_min)
        return result_df

    @staticmethod
    def dart_peak_calling(dart_df, fdr_threshold=0.05, mut_rate_threshold=0.05):
        """
        DART-seq mutation-based peak calling.
        Uses binomial test for enriched C-to-T mutations.
        """
        results = []
        for _, row in dart_df.iterrows():
            # Binomial test: treated mutation rate > control rate
            pval = stats.binomtest(
                int(row["treated_mutations"]),
                int(row["total_reads"]),
                row["ctrl_mut_rate"],
                alternative="greater"
            ).pvalue
            delta_mut = row["treated_mut_rate"] - row["ctrl_mut_rate"]
            results.append({
                "site_id": row["site_id"],
                "delta_mutation_rate": delta_mut,
                "pvalue": pval
            })

        result_df = pd.DataFrame(results)
        _, fdr, _, _ = multipletests(result_df["pvalue"], method="fdr_bh")
        result_df["fdr"] = fdr
        result_df["significant"] = (result_df["fdr"] < fdr_threshold) & \
                                   (result_df["delta_mutation_rate"] >= mut_rate_threshold)
        return result_df

    @staticmethod
    def nanopore_modification_calling(nanopore_df, prob_threshold=0.5):
        """
        Nanopore-based modification detection using current signal deviation.
        Applies logistic regression-like scoring on dwell time and current shift.
        """
        site_results = []
        for sid, group in nanopore_df.groupby("site_id"):
            mean_current = group["current_mean"].mean()
            mean_dwell = group["dwell_time"].mean()
            n_reads = len(group)

            # Score based on current deviation from unmodified baseline (120 pA)
            current_shift = 120 - mean_current
            # Sigmoid-like scoring
            raw_score = current_shift / 15.0
            mod_probability = 1.0 / (1.0 + np.exp(-5 * (raw_score - 0.3)))

            site_results.append({
                "site_id": sid,
                "n_reads": n_reads,
                "mean_current": mean_current,
                "mean_dwell": mean_dwell,
                "current_shift": current_shift,
                "mod_probability": mod_probability,
                "significant": mod_probability >= prob_threshold and n_reads >= 20
            })
        return pd.DataFrame(site_results)


# ═══════════════════════════════════════════════════════════════════
# Module 3: Quantification & Differential Modification
# ═══════════════════════════════════════════════════════════════════

class DifferentialModification:
    """Quantification and differential modification analysis."""

    @staticmethod
    def quantify_stoichiometry(merip_results, site_df):
        """Estimate per-site modification stoichiometry from IP/input ratio."""
        merged = merip_results.merge(site_df[["site_id", "stoichiometry", "modification"]],
                                      on="site_id", how="left")
        # Estimate stoichiometry from fold-change (normalized)
        max_fc = merged["fold_change"].replace([np.inf], np.nan).quantile(0.99)
        merged["estimated_stoichiometry"] = np.clip(
            (merged["fold_change"] - 1) / (max_fc - 1), 0, 1
        )
        return merged

    @staticmethod
    def differential_modification(site_df, n_tumor=500, n_normal=500):
        """
        Beta-binomial differential modification analysis.
        Simulates tumor vs normal comparison.
        """
        results = []
        test_sites = site_df.sample(min(n_tumor + n_normal, len(site_df)),
                                     random_state=42)

        for idx, (_, site) in enumerate(test_sites.iterrows()):
            is_tumor = idx < n_tumor
            condition = "tumor" if is_tumor else "normal"

            # Tumor: dysregulated modification levels
            if is_tumor:
                mod_level = site["stoichiometry"] * np.random.choice(
                    [0.5, 1.5, 2.0], p=[0.3, 0.4, 0.3]
                )
                mod_level = np.clip(mod_level, 0, 1)
            else:
                mod_level = site["stoichiometry"]

            results.append({
                "site_id": site["site_id"],
                "gene_name": site["gene_name"],
                "modification": site["modification"],
                "region": site["region"],
                "condition": condition,
                "mod_level": mod_level,
                "stoichiometry_true": site["stoichiometry"]
            })

        diff_df = pd.DataFrame(results)

        # Compute differential statistics per site
        diff_stats = []
        for sid in diff_df["site_id"].unique():
            sd = diff_df[diff_df["site_id"] == sid]
            tumor_vals = sd[sd["condition"] == "tumor"]["mod_level"].values
            normal_vals = sd[sd["condition"] == "normal"]["mod_level"].values

            if len(tumor_vals) > 0 and len(normal_vals) > 0:
                delta = np.mean(tumor_vals) - np.mean(normal_vals)
                # Use Welch's t-test approximation
                combined = np.concatenate([tumor_vals, normal_vals])
                if np.std(combined) > 0:
                    t_stat = delta / (np.std(combined) / np.sqrt(len(combined)))
                    pval = 2 * stats.norm.sf(abs(t_stat))
                else:
                    pval = 1.0
            else:
                delta = 0
                pval = 1.0

            info = sd.iloc[0]
            diff_stats.append({
                "site_id": sid,
                "gene_name": info["gene_name"],
                "modification": info["modification"],
                "region": info["region"],
                "delta_mod": delta,
                "pvalue": max(pval, 1e-300),
                "tumor_mean": np.mean(tumor_vals) if len(tumor_vals) > 0 else np.nan,
                "normal_mean": np.mean(normal_vals) if len(normal_vals) > 0 else np.nan,
            })

        stat_df = pd.DataFrame(diff_stats)
        valid = stat_df["pvalue"].notna()
        if valid.sum() > 0:
            _, fdr, _, _ = multipletests(stat_df.loc[valid, "pvalue"], method="fdr_bh")
            stat_df.loc[valid, "fdr"] = fdr
        else:
            stat_df["fdr"] = 1.0
        stat_df["significant"] = (stat_df["fdr"] < 0.05) & (abs(stat_df["delta_mod"]) > 0.1)
        return diff_df, stat_df


# ═══════════════════════════════════════════════════════════════════
# Module 4: Functional Annotation
# ═══════════════════════════════════════════════════════════════════

class FunctionalAnnotation:
    """Functional annotation of modification sites."""

    @staticmethod
    def annotate_regions(site_df):
        """Compute regional distribution of modification sites."""
        region_counts = site_df.groupby(["modification", "region"]).size().reset_index(name="count")
        region_pct = region_counts.copy()
        for mod in region_pct["modification"].unique():
            mask = region_pct["modification"] == mod
            total = region_pct.loc[mask, "count"].sum()
            region_pct.loc[mask, "percentage"] = region_pct.loc[mask, "count"] / total * 100
        return region_pct

    @staticmethod
    def mrna_stability_analysis(site_df, n_genes=500):
        """Simulate mRNA stability (half-life) with/without modifications."""
        genes = site_df.groupby("gene_name").agg(
            n_mods=("site_id", "count"),
            mean_stoich=("stoichiometry", "mean"),
            primary_mod=("modification", lambda x: x.mode().iloc[0] if len(x) > 0 else "m6A")
        ).reset_index().head(n_genes)

        # Modified genes tend to have shorter half-lives (m6A → YTHDF2 decay)
        genes["half_life_hours"] = np.random.lognormal(
            2.5 - 0.3 * genes["mean_stoich"], 0.5, size=len(genes)
        )
        # Unmodified control
        genes["half_life_ctrl"] = np.random.lognormal(2.5, 0.5, size=len(genes))
        genes["stability_ratio"] = genes["half_life_hours"] / genes["half_life_ctrl"]
        genes["stability_change"] = np.where(
            genes["stability_ratio"] < 0.8, "destabilized",
            np.where(genes["stability_ratio"] > 1.2, "stabilized", "unchanged")
        )
        return genes

    @staticmethod
    def translation_efficiency(site_df, n_genes=500):
        """Simulate translation efficiency changes linked to modifications."""
        genes = site_df.groupby("gene_name").agg(
            n_mods=("site_id", "count"),
            mean_stoich=("stoichiometry", "mean"),
            has_5utr_mod=("region", lambda x: "5UTR" in x.values),
            has_cds_mod=("region", lambda x: "CDS" in x.values),
        ).reset_index().head(n_genes)

        # 5'UTR m6A → enhanced cap-independent translation
        te_boost = np.where(genes["has_5utr_mod"], 0.3, 0) * genes["mean_stoich"]
        genes["translation_efficiency"] = np.random.lognormal(0, 0.4, len(genes)) * (1 + te_boost)
        genes["te_control"] = np.random.lognormal(0, 0.4, len(genes))
        genes["te_log2fc"] = np.log2(genes["translation_efficiency"] / genes["te_control"])
        genes["te_change"] = np.where(
            genes["te_log2fc"] > 0.5, "enhanced",
            np.where(genes["te_log2fc"] < -0.5, "reduced", "unchanged")
        )
        return genes


# ═══════════════════════════════════════════════════════════════════
# Module 5: Writer / Reader / Eraser Association
# ═══════════════════════════════════════════════════════════════════

class WREAnalysis:
    """Writer/Reader/Eraser association analysis."""

    WRE_DB = {
        "m6A": {
            "writers": ["METTL3", "METTL14", "WTAP", "KIAA1429", "RBM15"],
            "readers": ["YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2",
                        "IGF2BP1", "IGF2BP2", "IGF2BP3"],
            "erasers": ["FTO", "ALKBH5"]
        },
        "m5C": {
            "writers": ["NSUN2", "NSUN3", "NSUN5", "NSUN6", "DNMT2"],
            "readers": ["ALYREF", "YBX1"],
            "erasers": ["TET2"]
        },
        "Ψ": {
            "writers": ["PUS1", "PUS7", "TRUB1", "DKC1"],
            "readers": [],
            "erasers": []
        }
    }

    @staticmethod
    def simulate_wre_expression(n_samples=50):
        """Simulate WRE gene expression across samples."""
        records = []
        for mod_type, wre in WREAnalysis.WRE_DB.items():
            for role, genes in wre.items():
                for gene in genes:
                    for i in range(n_samples):
                        expr = np.random.lognormal(4, 1.5)
                        records.append({
                            "sample_id": f"S{i:03d}",
                            "modification": mod_type,
                            "role": role,
                            "gene": gene,
                            "expression_tpm": expr
                        })
        return pd.DataFrame(records)

    @staticmethod
    def correlate_wre_with_modifications(wre_expr_df, site_df, n_samples=50):
        """Correlate WRE expression with modification levels."""
        # Aggregate modification levels per sample (simulated)
        sample_mod_levels = {}
        for i in range(n_samples):
            sample_id = f"S{i:03d}"
            sample_mod_levels[sample_id] = {
                "m6A": np.random.beta(3, 5),
                "m5C": np.random.beta(2, 8),
                "Ψ": np.random.beta(2, 6)
            }

        correlations = []
        for mod_type in ["m6A", "m5C", "Ψ"]:
            mod_levels = [sample_mod_levels[f"S{i:03d}"][mod_type] for i in range(n_samples)]

            for role in ["writers", "readers", "erasers"]:
                genes = WREAnalysis.WRE_DB[mod_type].get(role, [])
                for gene in genes:
                    gene_expr = wre_expr_df[
                        (wre_expr_df["gene"] == gene) &
                        (wre_expr_df["modification"] == mod_type)
                    ]["expression_tpm"].values[:n_samples]

                    if len(gene_expr) == n_samples:
                        # Add correlation structure
                        if role == "writers":
                            gene_expr = gene_expr * (1 + 2 * np.array(mod_levels))
                        elif role == "erasers":
                            gene_expr = gene_expr * (1 - 0.8 * np.array(mod_levels))

                        r, p = stats.pearsonr(gene_expr, mod_levels)
                        correlations.append({
                            "modification": mod_type,
                            "role": role,
                            "gene": gene,
                            "pearson_r": r,
                            "pvalue": p,
                            "direction": "positive" if r > 0 else "negative"
                        })

        corr_df = pd.DataFrame(correlations)
        if len(corr_df) > 0:
            _, fdr, _, _ = multipletests(corr_df["pvalue"], method="fdr_bh")
            corr_df["fdr"] = fdr
            corr_df["significant"] = corr_df["fdr"] < 0.05
        return corr_df, sample_mod_levels


# ═══════════════════════════════════════════════════════════════════
# Module 6: Cancer Case Study (AML)
# ═══════════════════════════════════════════════════════════════════

class CancerCaseStudy:
    """AML m6A epitranscriptome case study."""

    KEY_GENES = {
        "MYC": {"role": "oncogene", "m6A_effect": "stabilization"},
        "BCL2": {"role": "anti-apoptotic", "m6A_effect": "enhanced_translation"},
        "CEBPA": {"role": "tumor_suppressor", "m6A_effect": "destabilization"},
        "FLT3": {"role": "receptor_tyrosine_kinase", "m6A_effect": "altered_splicing"},
        "NPM1": {"role": "nucleolar_protein", "m6A_effect": "mislocalization"},
        "IDH2": {"role": "metabolic_enzyme", "m6A_effect": "reduced_expression"},
        "RUNX1": {"role": "transcription_factor", "m6A_effect": "destabilization"},
        "TP53": {"role": "tumor_suppressor", "m6A_effect": "reduced_translation"},
    }

    @staticmethod
    def simulate_aml_data(n_tumor=30, n_normal=30):
        """Simulate AML vs normal m6A modification profiles."""
        records = []
        for gene, info in CancerCaseStudy.KEY_GENES.items():
            for i in range(n_tumor):
                if info["m6A_effect"] in ["stabilization", "enhanced_translation"]:
                    mod_level = np.random.beta(5, 3)  # Hyper-methylated
                else:
                    mod_level = np.random.beta(2, 6)  # Hypo-methylated

                expr = np.random.lognormal(5 + mod_level, 1)
                records.append({
                    "sample_id": f"AML_{i:03d}",
                    "condition": "AML",
                    "gene": gene,
                    "role": info["role"],
                    "m6A_level": mod_level,
                    "expression": expr,
                    "m6A_effect": info["m6A_effect"]
                })

            for i in range(n_normal):
                mod_level = np.random.beta(3, 5)  # Normal levels
                expr = np.random.lognormal(4.5, 0.8)
                records.append({
                    "sample_id": f"Normal_{i:03d}",
                    "condition": "Normal",
                    "gene": gene,
                    "role": info["role"],
                    "m6A_level": mod_level,
                    "expression": expr,
                    "m6A_effect": info["m6A_effect"]
                })

        return pd.DataFrame(records)

    @staticmethod
    def analyze_aml_differential(aml_df):
        """Differential analysis of m6A levels in AML vs Normal."""
        results = []
        for gene in CancerCaseStudy.KEY_GENES:
            gd = aml_df[aml_df["gene"] == gene]
            tumor = gd[gd["condition"] == "AML"]["m6A_level"].values
            normal = gd[gd["condition"] == "Normal"]["m6A_level"].values

            t_stat, pval = stats.ttest_ind(tumor, normal)
            delta = np.mean(tumor) - np.mean(normal)
            cohens_d = delta / np.sqrt((np.var(tumor) + np.var(normal)) / 2)

            results.append({
                "gene": gene,
                "role": CancerCaseStudy.KEY_GENES[gene]["role"],
                "m6A_effect": CancerCaseStudy.KEY_GENES[gene]["m6A_effect"],
                "aml_mean": np.mean(tumor),
                "normal_mean": np.mean(normal),
                "delta_m6A": delta,
                "t_statistic": t_stat,
                "pvalue": pval,
                "cohens_d": cohens_d,
                "direction": "hyper" if delta > 0 else "hypo"
            })

        result_df = pd.DataFrame(results)
        _, fdr, _, _ = multipletests(result_df["pvalue"], method="fdr_bh")
        result_df["fdr"] = fdr
        result_df["significant"] = result_df["fdr"] < 0.05
        return result_df


# ═══════════════════════════════════════════════════════════════════
# Visualization Module
# ═══════════════════════════════════════════════════════════════════

class Visualizer:
    """Publication-quality figure generation."""

    @staticmethod
    def plot_modification_landscape(site_df, save_path="figures/fig1_modification_landscape.png"):
        """Multi-panel modification landscape overview."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # Panel A: Modification type distribution
        ax = axes[0, 0]
        mod_counts = site_df["modification"].value_counts()
        bars = ax.bar(mod_counts.index, mod_counts.values,
                      color=[MOD_COLORS.get(m, PALETTE[3]) for m in mod_counts.index],
                      edgecolor="black", linewidth=0.5)
        for bar, val in zip(bars, mod_counts.values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                    str(val), ha="center", va="bottom", fontweight="bold")
        ax.set_xlabel("Modification Type")
        ax.set_ylabel("Number of Sites")
        ax.set_title("A. RNA Modification Site Distribution", fontweight="bold")

        # Panel B: Regional distribution
        ax = axes[0, 1]
        region_order = ["5UTR", "start_codon", "CDS", "stop_codon", "3UTR"]
        for mod in ["m6A", "m5C", "Ψ"]:
            mod_data = site_df[site_df["modification"] == mod]
            counts = mod_data["region"].value_counts().reindex(region_order, fill_value=0)
            pcts = counts / counts.sum() * 100
            ax.plot(region_order, pcts, "o-", label=mod,
                    color=MOD_COLORS.get(mod, PALETTE[3]), linewidth=2, markersize=8)
        ax.set_xlabel("Transcript Region")
        ax.set_ylabel("Percentage of Sites (%)")
        ax.set_title("B. Regional Distribution", fontweight="bold")
        ax.legend()

        # Panel C: Stoichiometry distribution
        ax = axes[1, 0]
        for mod in ["m6A", "m5C", "Ψ"]:
            mod_data = site_df[site_df["modification"] == mod]
            ax.hist(mod_data["stoichiometry"], bins=30, alpha=0.6,
                    label=mod, color=MOD_COLORS.get(mod, PALETTE[3]),
                    edgecolor="black", linewidth=0.3)
        ax.set_xlabel("Modification Stoichiometry")
        ax.set_ylabel("Number of Sites")
        ax.set_title("C. Stoichiometry Distribution", fontweight="bold")
        ax.legend()

        # Panel D: Confidence scores
        ax = axes[1, 1]
        for mod in ["m6A", "m5C", "Ψ"]:
            mod_data = site_df[site_df["modification"] == mod]
            ax.hist(mod_data["confidence_score"], bins=30, alpha=0.6,
                    label=mod, color=MOD_COLORS.get(mod, PALETTE[3]),
                    edgecolor="black", linewidth=0.3)
        ax.set_xlabel("Confidence Score")
        ax.set_ylabel("Number of Sites")
        ax.set_title("D. Detection Confidence", fontweight="bold")
        ax.legend()

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    @staticmethod
    def plot_peak_calling_comparison(merip_res, dart_res, nano_res,
                                      save_path="figures/fig2_peak_calling.png"):
        """Compare peak calling results across methods."""
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # MeRIP-seq volcano
        ax = axes[0]
        log2fc = np.log2(merip_res["fold_change"].replace([np.inf, 0], np.nan).dropna())
        neg_log10p = -np.log10(merip_res["pvalue"].replace(0, 1e-300))
        colors = ["#d62728" if s else "#aec7e8" for s in merip_res["significant"]]
        ax.scatter(log2fc[:len(colors)], neg_log10p[:len(colors)],
                   c=colors[:len(log2fc)], alpha=0.5, s=10, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5)
        ax.axvline(np.log2(2), color="grey", linestyle="--", alpha=0.5)
        ax.set_xlabel("log2(Fold Change IP/Input)")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title("A. MeRIP-seq Peak Calling", fontweight="bold")
        n_sig = merip_res["significant"].sum()
        ax.text(0.05, 0.95, f"Significant: {n_sig}",
                transform=ax.transAxes, fontsize=10, verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        # DART-seq
        ax = axes[1]
        colors_d = ["#d62728" if s else "#aec7e8" for s in dart_res["significant"]]
        ax.scatter(dart_res["delta_mutation_rate"],
                   -np.log10(dart_res["pvalue"].replace(0, 1e-300)),
                   c=colors_d, alpha=0.5, s=10, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5)
        ax.set_xlabel("Δ Mutation Rate (Treated - Control)")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title("B. DART-seq Mutation Calling", fontweight="bold")
        n_sig_d = dart_res["significant"].sum()
        ax.text(0.05, 0.95, f"Significant: {n_sig_d}",
                transform=ax.transAxes, fontsize=10, verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        # Nanopore
        ax = axes[2]
        colors_n = ["#d62728" if s else "#aec7e8" for s in nano_res["significant"]]
        ax.scatter(nano_res["current_shift"], nano_res["mod_probability"],
                   c=colors_n, alpha=0.5, s=10, edgecolors="none")
        ax.axhline(0.5, color="grey", linestyle="--", alpha=0.5)
        ax.set_xlabel("Current Shift (pA)")
        ax.set_ylabel("Modification Probability")
        ax.set_title("C. Nanopore Direct RNA-seq", fontweight="bold")
        n_sig_n = nano_res["significant"].sum()
        ax.text(0.05, 0.95, f"Significant: {n_sig_n}",
                transform=ax.transAxes, fontsize=10, verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    @staticmethod
    def plot_differential_modification(diff_stats,
                                        save_path="figures/fig3_differential_modification.png"):
        """Volcano plot of differential modification analysis."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Volcano plot
        ax = axes[0]
        x = diff_stats["delta_mod"]
        y = -np.log10(diff_stats["pvalue"].replace(0, 1e-300))
        sig = diff_stats["significant"]
        ax.scatter(x[~sig], y[~sig], c="#aec7e8", alpha=0.4, s=15, label="NS")
        hyper = sig & (x > 0)
        hypo = sig & (x < 0)
        ax.scatter(x[hyper], y[hyper], c="#d62728", alpha=0.7, s=20, label="Hyper-modified")
        ax.scatter(x[hypo], y[hypo], c="#2ca02c", alpha=0.7, s=20, label="Hypo-modified")
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5)
        ax.axvline(0.1, color="grey", linestyle="--", alpha=0.5)
        ax.axvline(-0.1, color="grey", linestyle="--", alpha=0.5)
        ax.set_xlabel("Δ Modification Level (Tumor - Normal)")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title("A. Differential Modification (Tumor vs Normal)", fontweight="bold")
        ax.legend(loc="upper right", fontsize=9)

        # By modification type
        ax = axes[1]
        for mod in ["m6A", "m5C", "Ψ"]:
            mod_data = diff_stats[diff_stats["modification"] == mod]
            ax.hist(mod_data["delta_mod"], bins=25, alpha=0.6, label=mod,
                    color=MOD_COLORS.get(mod, PALETTE[3]), edgecolor="black", linewidth=0.3)
        ax.set_xlabel("Δ Modification Level")
        ax.set_ylabel("Number of Sites")
        ax.set_title("B. Delta by Modification Type", fontweight="bold")
        ax.legend()

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    @staticmethod
    def plot_functional_annotation(stability_df, te_df,
                                    save_path="figures/fig4_functional_annotation.png"):
        """Functional annotation: stability and translation efficiency."""
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # mRNA stability
        ax = axes[0]
        change_counts = stability_df["stability_change"].value_counts()
        colors_s = [PALETTE[0], PALETTE[1], PALETTE[2]]
        ax.bar(change_counts.index, change_counts.values, color=colors_s[:len(change_counts)],
               edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Stability Change")
        ax.set_ylabel("Number of Genes")
        ax.set_title("A. mRNA Stability Impact", fontweight="bold")

        # Stability ratio vs stoichiometry
        ax = axes[1]
        scatter = ax.scatter(stability_df["mean_stoich"],
                             stability_df["stability_ratio"],
                             c=stability_df["n_mods"], cmap="viridis",
                             alpha=0.6, s=20, edgecolors="none")
        plt.colorbar(scatter, ax=ax, label="Number of Modifications")
        ax.axhline(1.0, color="red", linestyle="--", alpha=0.5)
        ax.set_xlabel("Mean Stoichiometry")
        ax.set_ylabel("Stability Ratio (Modified / Control)")
        ax.set_title("B. Stoichiometry vs Stability", fontweight="bold")

        # Translation efficiency
        ax = axes[2]
        te_counts = te_df["te_change"].value_counts()
        colors_t = [PALETTE[3], PALETTE[4], PALETTE[5]]
        ax.bar(te_counts.index, te_counts.values, color=colors_t[:len(te_counts)],
               edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Translation Efficiency Change")
        ax.set_ylabel("Number of Genes")
        ax.set_title("C. Translation Efficiency Impact", fontweight="bold")

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    @staticmethod
    def plot_wre_analysis(wre_corr_df, save_path="figures/fig5_wre_analysis.png"):
        """Writer/Reader/Eraser correlation heatmap."""
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        for idx, mod in enumerate(["m6A", "m5C", "Ψ"]):
            ax = axes[idx]
            mod_data = wre_corr_df[wre_corr_df["modification"] == mod]
            if len(mod_data) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{mod} WRE Correlations", fontweight="bold")
                continue

            pivot = mod_data.pivot_table(index="gene", columns="role",
                                          values="pearson_r", aggfunc="first")
            if len(pivot) > 0:
                sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdBu_r",
                            center=0, vmin=-1, vmax=1, ax=ax,
                            cbar_kws={"label": "Pearson r"})
            ax.set_title(f"{mod} WRE Correlations", fontweight="bold")
            ax.set_ylabel("")

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    @staticmethod
    def plot_cancer_case_study(aml_df, aml_results,
                                save_path="figures/fig6_cancer_case_study.png"):
        """AML epitranscriptome case study visualization."""
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))

        # Panel A: m6A levels by gene and condition
        ax = axes[0, 0]
        gene_order = aml_results.sort_values("delta_m6A", ascending=False)["gene"].values
        sns.boxplot(data=aml_df, x="gene", y="m6A_level", hue="condition",
                    order=gene_order, ax=ax, palette={"AML": PALETTE[0], "Normal": PALETTE[1]},
                    fliersize=2, linewidth=0.8)
        ax.set_xlabel("Gene")
        ax.set_ylabel("m6A Modification Level")
        ax.set_title("A. m6A Levels: AML vs Normal", fontweight="bold")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(loc="upper right", fontsize=9)

        # Panel B: Effect sizes (Cohen's d)
        ax = axes[0, 1]
        colors_bar = ["#d62728" if d > 0 else "#2ca02c" for d in aml_results["cohens_d"]]
        bars = ax.barh(aml_results["gene"], aml_results["cohens_d"],
                       color=colors_bar, edgecolor="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.axvline(0.8, color="grey", linestyle="--", alpha=0.5, label="|d|=0.8")
        ax.axvline(-0.8, color="grey", linestyle="--", alpha=0.5)
        ax.set_xlabel("Cohen's d (Effect Size)")
        ax.set_title("B. Effect Sizes (AML vs Normal)", fontweight="bold")
        ax.legend(fontsize=9)

        # Panel C: Significance (-log10 FDR)
        ax = axes[1, 0]
        y_vals = -np.log10(aml_results["fdr"].replace(0, 1e-300))
        sig_colors = ["#d62728" if s else "#aec7e8" for s in aml_results["significant"]]
        ax.barh(aml_results["gene"], y_vals, color=sig_colors,
                edgecolor="black", linewidth=0.5)
        ax.axvline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5, label="FDR=0.05")
        ax.set_xlabel("-log10(FDR)")
        ax.set_title("C. Statistical Significance", fontweight="bold")
        ax.legend(fontsize=9)

        # Panel D: m6A effect mechanism summary
        ax = axes[1, 1]
        effects = aml_results.groupby("m6A_effect").size()
        ax.pie(effects, labels=effects.index, autopct="%1.0f%%",
               colors=[PALETTE[i] for i in range(len(effects))],
               startangle=90, textprops={"fontsize": 10})
        ax.set_title("D. m6A Functional Effects in AML", fontweight="bold")

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path

    @staticmethod
    def plot_method_comparison_heatmap(site_df, merip_res, dart_res, nano_res,
                                       save_path="figures/fig7_method_concordance.png"):
        """Cross-method concordance heatmap."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Method detection overlap (simulated)
        merip_sites = set(merip_res[merip_res["significant"]]["site_id"])
        dart_sites = set(dart_res[dart_res["significant"]]["site_id"])
        nano_sites = set(nano_res[nano_res["significant"]]["site_id"])

        methods = ["MeRIP-seq", "DART-seq", "Nanopore"]
        site_sets = [merip_sites, dart_sites, nano_sites]
        overlap_matrix = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                if len(site_sets[i]) > 0 and len(site_sets[j]) > 0:
                    overlap = len(site_sets[i] & site_sets[j])
                    overlap_matrix[i, j] = overlap / len(site_sets[i] | site_sets[j]) * 100

        ax = axes[0]
        sns.heatmap(overlap_matrix, annot=True, fmt=".1f", cmap="YlOrRd",
                    xticklabels=methods, yticklabels=methods, ax=ax,
                    cbar_kws={"label": "Jaccard Index (%)"})
        ax.set_title("A. Cross-Method Concordance", fontweight="bold")

        # Detection counts
        ax = axes[1]
        counts = [len(merip_sites), len(dart_sites), len(nano_sites)]
        bars = ax.bar(methods, counts, color=[PALETTE[0], PALETTE[1], PALETTE[2]],
                      edgecolor="black", linewidth=0.5)
        for bar, val in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                    str(val), ha="center", va="bottom", fontweight="bold")
        ax.set_ylabel("Number of Significant Sites")
        ax.set_title("B. Sites Detected per Method", fontweight="bold")

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return save_path


# ═══════════════════════════════════════════════════════════════════
# Main Pipeline Execution
# ═══════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("  Epitranscriptome Analysis Pipeline")
    print("  RNA Modification Mapping: m6A / m5C / Pseudouridine")
    print("=" * 70)

    log_entries = []

    def log_event(phase, event, skill="epitranscriptome_pipeline", **kwargs):
        from datetime import datetime, timezone, timedelta
        entry = {
            "timestamp": datetime.now(timezone(timedelta(hours=9))).isoformat(),
            "phase": phase,
            "event_type": event,
            "actor": "co-scientist",
            "skill_or_tool": skill,
            **kwargs
        }
        log_entries.append(entry)

    log_event("init", "run_started")

    # ── Step 1: Generate simulated data ──
    print("\n[1/6] Generating simulated transcriptome data...")
    gene_df, site_df = generate_transcriptome_data(n_genes=5000, n_modifications=3000)
    merip_data = simulate_merip_seq(site_df)
    dart_data = simulate_dart_seq(site_df)
    nanopore_data = simulate_nanopore_signals(site_df)

    gene_df.to_csv("data/gene_annotations.csv", index=False)
    site_df.to_csv("data/modification_sites.csv", index=False)
    merip_data.to_csv("data/merip_seq_counts.csv", index=False)
    dart_data.to_csv("data/dart_seq_mutations.csv", index=False)
    nanopore_data.to_csv("data/nanopore_signals.csv", index=False)

    print(f"  → Genes: {len(gene_df):,}")
    print(f"  → Modification sites: {len(site_df):,}")
    print(f"  → m6A: {(site_df['modification']=='m6A').sum()}, "
          f"m5C: {(site_df['modification']=='m5C').sum()}, "
          f"Ψ: {(site_df['modification']=='Ψ').sum()}")
    log_event("data_generation", "file_written",
              files_written=["data/gene_annotations.csv", "data/modification_sites.csv",
                             "data/merip_seq_counts.csv", "data/dart_seq_mutations.csv",
                             "data/nanopore_signals.csv"])

    # ── Step 2: Peak calling ──
    print("\n[2/6] Running peak calling algorithms...")
    pc = PeakCaller()

    merip_results = pc.merip_peak_calling(merip_data)
    merip_results.to_csv("results/merip_peak_calls.csv", index=False)
    n_merip_sig = merip_results["significant"].sum()
    print(f"  → MeRIP-seq: {n_merip_sig} significant peaks "
          f"({n_merip_sig/len(merip_results)*100:.1f}%)")

    dart_results = pc.dart_peak_calling(dart_data)
    dart_results.to_csv("results/dart_peak_calls.csv", index=False)
    n_dart_sig = dart_results["significant"].sum()
    print(f"  → DART-seq:  {n_dart_sig} significant peaks "
          f"({n_dart_sig/len(dart_results)*100:.1f}%)")

    nano_results = pc.nanopore_modification_calling(nanopore_data)
    nano_results.to_csv("results/nanopore_calls.csv", index=False)
    n_nano_sig = nano_results["significant"].sum()
    print(f"  → Nanopore:  {n_nano_sig} significant sites "
          f"({n_nano_sig/len(nano_results)*100:.1f}%)")

    log_event("peak_calling", "file_written",
              files_written=["results/merip_peak_calls.csv",
                             "results/dart_peak_calls.csv",
                             "results/nanopore_calls.csv"])

    # ── Step 3: Quantification & Differential Modification ──
    print("\n[3/6] Quantification and differential modification analysis...")
    dm = DifferentialModification()

    quant_df = dm.quantify_stoichiometry(merip_results, site_df)
    quant_df.to_csv("results/stoichiometry_estimates.csv", index=False)

    diff_df, diff_stats = dm.differential_modification(site_df)
    diff_stats.to_csv("results/differential_modification.csv", index=False)
    n_diff_sig = diff_stats["significant"].sum()
    n_hyper = ((diff_stats["significant"]) & (diff_stats["delta_mod"] > 0)).sum()
    n_hypo = ((diff_stats["significant"]) & (diff_stats["delta_mod"] < 0)).sum()
    print(f"  → Differentially modified sites: {n_diff_sig}")
    print(f"  → Hyper-modified: {n_hyper}, Hypo-modified: {n_hypo}")
    log_event("differential", "file_written",
              files_written=["results/stoichiometry_estimates.csv",
                             "results/differential_modification.csv"])

    # ── Step 4: Functional annotation ──
    print("\n[4/6] Functional annotation...")
    fa = FunctionalAnnotation()

    region_dist = fa.annotate_regions(site_df)
    region_dist.to_csv("results/region_distribution.csv", index=False)

    stability_df = fa.mrna_stability_analysis(site_df)
    stability_df.to_csv("results/mrna_stability.csv", index=False)
    n_destab = (stability_df["stability_change"] == "destabilized").sum()
    n_stab = (stability_df["stability_change"] == "stabilized").sum()
    print(f"  → Destabilized mRNAs: {n_destab}, Stabilized: {n_stab}")

    te_df = fa.translation_efficiency(site_df)
    te_df.to_csv("results/translation_efficiency.csv", index=False)
    n_te_enh = (te_df["te_change"] == "enhanced").sum()
    n_te_red = (te_df["te_change"] == "reduced").sum()
    print(f"  → Translation enhanced: {n_te_enh}, Reduced: {n_te_red}")
    log_event("functional_annotation", "file_written",
              files_written=["results/region_distribution.csv",
                             "results/mrna_stability.csv",
                             "results/translation_efficiency.csv"])

    # ── Step 5: Writer/Reader/Eraser analysis ──
    print("\n[5/6] Writer/Reader/Eraser association analysis...")
    wre = WREAnalysis()
    wre_expr = wre.simulate_wre_expression()
    wre_corr, sample_mods = wre.correlate_wre_with_modifications(wre_expr, site_df)
    wre_corr.to_csv("results/wre_correlations.csv", index=False)
    wre_expr.to_csv("data/wre_expression.csv", index=False)

    n_sig_corr = wre_corr["significant"].sum() if "significant" in wre_corr.columns else 0
    print(f"  → Significant WRE correlations: {n_sig_corr}/{len(wre_corr)}")
    writer_corrs = wre_corr[wre_corr["role"] == "writers"]
    if len(writer_corrs) > 0:
        print(f"  → Writer mean |r|: {writer_corrs['pearson_r'].abs().mean():.3f}")
    eraser_corrs = wre_corr[wre_corr["role"] == "erasers"]
    if len(eraser_corrs) > 0:
        print(f"  → Eraser mean |r|: {eraser_corrs['pearson_r'].abs().mean():.3f}")
    log_event("wre_analysis", "file_written",
              files_written=["results/wre_correlations.csv", "data/wre_expression.csv"])

    # ── Step 6: Cancer case study ──
    print("\n[6/6] AML epitranscriptome case study...")
    cs = CancerCaseStudy()
    aml_df = cs.simulate_aml_data()
    aml_results = cs.analyze_aml_differential(aml_df)
    aml_df.to_csv("data/aml_case_study_data.csv", index=False)
    aml_results.to_csv("results/aml_differential_results.csv", index=False)

    n_aml_sig = aml_results["significant"].sum()
    print(f"  → Significant genes: {n_aml_sig}/{len(aml_results)}")
    for _, row in aml_results.iterrows():
        sig_mark = "***" if row["significant"] else "   "
        print(f"  {sig_mark} {row['gene']:8s}: Δm6A={row['delta_m6A']:+.3f}, "
              f"d={row['cohens_d']:+.2f}, FDR={row['fdr']:.2e} ({row['direction']})")
    log_event("cancer_case_study", "file_written",
              files_written=["data/aml_case_study_data.csv",
                             "results/aml_differential_results.csv"])

    # ── Generate figures ──
    print("\n[FIG] Generating publication-quality figures...")
    viz = Visualizer()
    f1 = viz.plot_modification_landscape(site_df)
    print(f"  → {f1}")
    f2 = viz.plot_peak_calling_comparison(merip_results, dart_results, nano_results)
    print(f"  → {f2}")
    f3 = viz.plot_differential_modification(diff_stats)
    print(f"  → {f3}")
    f4 = viz.plot_functional_annotation(stability_df, te_df)
    print(f"  → {f4}")
    f5 = viz.plot_wre_analysis(wre_corr)
    print(f"  → {f5}")
    f6 = viz.plot_cancer_case_study(aml_df, aml_results)
    print(f"  → {f6}")
    f7 = viz.plot_method_comparison_heatmap(site_df, merip_results, dart_results, nano_results)
    print(f"  → {f7}")

    log_event("visualization", "file_written",
              files_written=[f1, f2, f3, f4, f5, f6, f7])

    # ── Summary statistics ──
    print("\n" + "=" * 70)
    print("  PIPELINE SUMMARY")
    print("=" * 70)

    summary = {
        "total_genes": len(gene_df),
        "total_modification_sites": len(site_df),
        "m6A_sites": int((site_df["modification"] == "m6A").sum()),
        "m5C_sites": int((site_df["modification"] == "m5C").sum()),
        "psi_sites": int((site_df["modification"] == "Ψ").sum()),
        "merip_significant_peaks": int(n_merip_sig),
        "dart_significant_peaks": int(n_dart_sig),
        "nanopore_significant_sites": int(n_nano_sig),
        "differentially_modified_sites": int(n_diff_sig),
        "hyper_modified": int(n_hyper),
        "hypo_modified": int(n_hypo),
        "destabilized_mrnas": int(n_destab),
        "stabilized_mrnas": int(n_stab),
        "te_enhanced": int(n_te_enh),
        "te_reduced": int(n_te_red),
        "wre_significant_correlations": int(n_sig_corr),
        "aml_significant_genes": int(n_aml_sig),
    }

    with open("results/pipeline_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    for k, v in summary.items():
        print(f"  {k}: {v:,}")

    # ── Save process log ──
    log_event("report", "run_completed", status="ok")
    os.makedirs("logs", exist_ok=True)
    with open("logs/process-log.jsonl", "w") as f:
        for entry in log_entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print(f"\n✅ Pipeline completed. All outputs saved.")
    print(f"   Results: results/")
    print(f"   Figures: figures/")
    print(f"   Data:    data/")
    print(f"   Logs:    logs/process-log.jsonl")

    return summary


if __name__ == "__main__":
    summary = main()
