#!/usr/bin/env python3
"""
Module 4: 代謝パスウェイ富化解析（微生物代謝 + 宿主代謝の統合）
Integrated Metabolic Pathway Enrichment Analysis

Components:
  1. KEGG pathway enrichment (ORA + GSEA)
  2. Microbial metabolic pathway mapping (MetaCyc, HUMAnN3)
  3. Host metabolic pathway mapping (Reactome, KEGG)
  4. Joint pathway topology analysis
  5. Pathway-level integration score
"""

import os
import json
import logging
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# KEGG Pathway Database (curated subset)
# ---------------------------------------------------------------------------

MICROBIAL_PATHWAYS = {
    "ko00650": {
        "name": "Butanoate metabolism",
        "compounds": ["C00246", "C00042", "C00024", "C01412", "C00332"],
        "source": "microbial",
        "category": "SCFA biosynthesis",
    },
    "ko00640": {
        "name": "Propanoate metabolism",
        "compounds": ["C00163", "C00024", "C00100", "C05668"],
        "source": "microbial",
        "category": "SCFA biosynthesis",
    },
    "ko00380": {
        "name": "Tryptophan metabolism",
        "compounds": ["C00078", "C02693", "C00331", "C00643", "C01598"],
        "source": "both",
        "category": "Amino acid metabolism",
    },
    "ko00120": {
        "name": "Primary bile acid biosynthesis",
        "compounds": ["C02691", "C00695", "C05122", "C01921"],
        "source": "both",
        "category": "Bile acid metabolism",
    },
    "ko00121": {
        "name": "Secondary bile acid biosynthesis",
        "compounds": ["C02691", "C03990", "C05465", "C05466"],
        "source": "microbial",
        "category": "Bile acid metabolism",
    },
    "ko00680": {
        "name": "Methane metabolism",
        "compounds": ["C00067", "C00132", "C01438"],
        "source": "microbial",
        "category": "Energy metabolism",
    },
    "ko00910": {
        "name": "Nitrogen metabolism",
        "compounds": ["C00014", "C00088", "C00244", "C00169"],
        "source": "microbial",
        "category": "Nitrogen cycling",
    },
}

HOST_PATHWAYS = {
    "hsa00590": {
        "name": "Arachidonic acid metabolism",
        "compounds": ["C00219", "C00427", "C05356", "C05966", "C14768"],
        "source": "host",
        "category": "Lipid metabolism",
    },
    "hsa04080": {
        "name": "Neuroactive ligand-receptor interaction",
        "compounds": ["C00078", "C00187", "C00334", "C00547"],
        "source": "host",
        "category": "Signaling",
    },
    "hsa00760": {
        "name": "Nicotinate and nicotinamide metabolism",
        "compounds": ["C00153", "C00253", "C03150"],
        "source": "host",
        "category": "Vitamin metabolism",
    },
    "hsa00340": {
        "name": "Histidine metabolism",
        "compounds": ["C00135", "C00388", "C05130", "C00025"],
        "source": "host",
        "category": "Amino acid metabolism",
    },
    "hsa04979": {
        "name": "Cholesterol metabolism",
        "compounds": ["C00187", "C02530", "C05122"],
        "source": "host",
        "category": "Lipid metabolism",
    },
}

ALL_PATHWAYS = {**MICROBIAL_PATHWAYS, **HOST_PATHWAYS}


# ---------------------------------------------------------------------------
# Over-Representation Analysis (ORA)
# ---------------------------------------------------------------------------

def pathway_ora(significant_compounds: list, background_compounds: list,
                pathway_db: dict = None) -> pd.DataFrame:
    """
    Fisher's exact test によるパスウェイ ORA
    """
    logger.info("Running pathway Over-Representation Analysis (ORA)")
    if pathway_db is None:
        pathway_db = ALL_PATHWAYS

    sig_set = set(significant_compounds)
    bg_set = set(background_compounds)
    n_bg = len(bg_set)
    n_sig = len(sig_set)

    results = []
    for pw_id, pw_info in pathway_db.items():
        pw_compounds = set(pw_info["compounds"])
        pw_in_bg = pw_compounds & bg_set
        pw_in_sig = pw_compounds & sig_set

        a = len(pw_in_sig)
        b = len(pw_in_bg) - a
        c = n_sig - a
        d = n_bg - a - b - c

        if a == 0:
            continue

        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        results.append({
            "pathway_id": pw_id,
            "pathway_name": pw_info["name"],
            "source": pw_info["source"],
            "category": pw_info["category"],
            "hits": a,
            "pathway_size": len(pw_in_bg),
            "odds_ratio": round(odds_ratio, 3),
            "pvalue": round(p_value, 6),
            "hit_compounds": list(pw_in_sig),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
        df["fdr_qvalue"] = fdr
        df = df.sort_values("pvalue")

    logger.info(f"ORA: {len(df)} pathways tested, {(df['pvalue'] < 0.05).sum()} significant")
    return df


# ---------------------------------------------------------------------------
# Gene Set Enrichment Analysis (GSEA) style
# ---------------------------------------------------------------------------

def pathway_gsea(compound_scores: dict, pathway_db: dict = None,
                 n_perm: int = 1000) -> pd.DataFrame:
    """
    GSEA-style enrichment using compound-level scores (e.g., fold change)
    """
    logger.info("Running pathway GSEA")
    if pathway_db is None:
        pathway_db = ALL_PATHWAYS

    sorted_compounds = sorted(compound_scores.keys(), key=lambda x: compound_scores[x], reverse=True)
    n_total = len(sorted_compounds)

    results = []
    for pw_id, pw_info in pathway_db.items():
        pw_set = set(pw_info["compounds"])
        pw_in_data = pw_set & set(sorted_compounds)

        if len(pw_in_data) < 2:
            continue

        # Calculate enrichment score (running sum)
        n_hit = len(pw_in_data)
        n_miss = n_total - n_hit
        es_max = 0
        running_sum = 0

        for comp in sorted_compounds:
            if comp in pw_in_data:
                running_sum += 1.0 / n_hit
            else:
                running_sum -= 1.0 / max(n_miss, 1)
            if abs(running_sum) > abs(es_max):
                es_max = running_sum

        # Permutation test
        np.random.seed(42)
        null_es = []
        for _ in range(n_perm):
            perm = np.random.permutation(sorted_compounds)
            rs = 0
            es_null = 0
            for comp in perm:
                if comp in pw_in_data:
                    rs += 1.0 / n_hit
                else:
                    rs -= 1.0 / max(n_miss, 1)
                if abs(rs) > abs(es_null):
                    es_null = rs
            null_es.append(es_null)

        p_value = np.mean(np.abs(null_es) >= abs(es_max))

        results.append({
            "pathway_id": pw_id,
            "pathway_name": pw_info["name"],
            "source": pw_info["source"],
            "enrichment_score": round(es_max, 4),
            "normalized_es": round(es_max / (np.std(null_es) + 1e-10), 4),
            "pvalue": round(p_value, 4),
            "n_hits": n_hit,
        })

    df = pd.DataFrame(results).sort_values("pvalue")
    logger.info(f"GSEA: {len(df)} pathways, {(df['pvalue'] < 0.05).sum()} enriched")
    return df


# ---------------------------------------------------------------------------
# Joint Pathway Topology Analysis
# ---------------------------------------------------------------------------

def joint_pathway_analysis(microbial_enrichment: pd.DataFrame,
                           host_enrichment: pd.DataFrame) -> pd.DataFrame:
    """
    微生物–宿主 統合パスウェイ解析:
    Combined p-values (Fisher's method) for pathways present in both domains.
    """
    logger.info("Running joint pathway topology analysis")

    micro = microbial_enrichment.set_index("pathway_name")
    host = host_enrichment.set_index("pathway_name")

    shared_categories = set()
    for _, row in microbial_enrichment.iterrows():
        for _, hrow in host_enrichment.iterrows():
            if row.get("category") == hrow.get("category"):
                shared_categories.add(row["category"])

    joint_results = []
    all_pathways_combined = pd.concat([microbial_enrichment, host_enrichment])

    for category in shared_categories:
        cat_pathways = all_pathways_combined[all_pathways_combined["category"] == category]
        pvals = cat_pathways["pvalue"].values
        pvals = pvals[pvals > 0]

        if len(pvals) >= 2:
            chi2_stat = -2 * np.sum(np.log(pvals))
            combined_p = 1 - stats.chi2.cdf(chi2_stat, 2 * len(pvals))

            joint_results.append({
                "category": category,
                "n_pathways": len(pvals),
                "fisher_chi2": round(chi2_stat, 2),
                "combined_pvalue": round(combined_p, 6),
                "pathways_involved": list(cat_pathways["pathway_name"].unique()),
            })

    df = pd.DataFrame(joint_results).sort_values("combined_pvalue")
    logger.info(f"Joint analysis: {len(df)} shared categories identified")
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pathway_enrichment_pipeline(output_dir: str = "results") -> dict:
    os.makedirs(output_dir, exist_ok=True)

    # Simulated significant compounds (from differential analysis)
    np.random.seed(42)
    all_kegg_compounds = []
    for pw in ALL_PATHWAYS.values():
        all_kegg_compounds.extend(pw["compounds"])
    all_kegg_compounds = list(set(all_kegg_compounds))

    # Significant compounds (IBD vs Control)
    sig_compounds = list(np.random.choice(
        all_kegg_compounds,
        size=min(15, len(all_kegg_compounds)),
        replace=False
    ))

    # Compound scores for GSEA
    compound_scores = {c: np.random.normal(0, 1) for c in all_kegg_compounds}
    for c in sig_compounds:
        compound_scores[c] = abs(np.random.normal(2, 0.5))

    # ORA - microbial pathways
    ora_micro = pathway_ora(sig_compounds, all_kegg_compounds, MICROBIAL_PATHWAYS)
    ora_host = pathway_ora(sig_compounds, all_kegg_compounds, HOST_PATHWAYS)

    # GSEA
    gsea_results = pathway_gsea(compound_scores)

    # Joint analysis
    ora_all = pathway_ora(sig_compounds, all_kegg_compounds, ALL_PATHWAYS)
    joint = joint_pathway_analysis(
        ora_all[ora_all["source"].isin(["microbial", "both"])],
        ora_all[ora_all["source"].isin(["host", "both"])],
    )

    # Save results
    ora_all.to_csv(os.path.join(output_dir, "pathway_ora_results.csv"), index=False)
    gsea_results.to_csv(os.path.join(output_dir, "pathway_gsea_results.csv"), index=False)
    joint.to_csv(os.path.join(output_dir, "joint_pathway_results.csv"), index=False)

    summary = {
        "total_pathways_tested": len(ALL_PATHWAYS),
        "significant_ora": int((ora_all["pvalue"] < 0.05).sum()) if len(ora_all) > 0 else 0,
        "significant_gsea": int((gsea_results["pvalue"] < 0.05).sum()) if len(gsea_results) > 0 else 0,
        "shared_categories": len(joint),
        "top_pathway": ora_all.iloc[0]["pathway_name"] if len(ora_all) > 0 else "N/A",
        "microbial_pathways": len(ora_micro),
        "host_pathways": len(ora_host),
    }

    with open(os.path.join(output_dir, "pathway_enrichment_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    return summary


if __name__ == "__main__":
    summary = run_pathway_enrichment_pipeline(output_dir="../results")
    print(json.dumps(summary, indent=2))
