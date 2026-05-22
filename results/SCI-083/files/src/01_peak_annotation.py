#!/usr/bin/env python3
"""
Module 1: 非標的メタボロミクス ピーク同定・アノテーション自動化パイプライン
Untargeted Metabolomics Peak Identification & Annotation Pipeline

Workflow:
  1. mzML/mzXML raw data → feature extraction (pyOpenMS / XCMS via rpy2)
  2. Peak grouping & retention time alignment
  3. Adduct/isotope deconvolution
  4. Database matching (HMDB, KEGG, MassBank, METLIN)
  5. MS2 spectral matching (cosine similarity)
  6. Confidence scoring (Metabolomics Standards Initiative levels 1-4)
"""

import os
import json
import logging
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PeakFeature:
    feature_id: str
    mz: float
    rt: float  # retention time (seconds)
    intensity: float
    adduct: str = "[M+H]+"
    isotope_pattern: list = field(default_factory=list)
    ms2_spectrum: Optional[dict] = None


@dataclass
class Annotation:
    feature_id: str
    compound_name: str
    formula: str
    inchikey: str
    hmdb_id: str
    kegg_id: str
    mass_error_ppm: float
    ms2_cosine_score: float
    msi_level: int  # Metabolomics Standards Initiative confidence level (1-4)


# ---------------------------------------------------------------------------
# Feature extraction (pyOpenMS wrapper)
# ---------------------------------------------------------------------------

def extract_features_pyopenms(mzml_path: str, noise_threshold: float = 1e4) -> pd.DataFrame:
    """
    pyOpenMS を用いた特徴量抽出。
    本番では FeatureFinderMetabo を使用。
    ここではシミュレーションデータで代替。
    """
    logger.info(f"Feature extraction from: {mzml_path}")

    np.random.seed(42)
    n_features = 2500
    features = pd.DataFrame({
        "feature_id": [f"F{i:05d}" for i in range(n_features)],
        "mz": np.random.uniform(80, 1200, n_features).round(4),
        "rt": np.random.uniform(30, 1800, n_features).round(2),
        "intensity": np.random.lognormal(mean=12, sigma=2, size=n_features).round(0),
        "adduct": np.random.choice(
            ["[M+H]+", "[M+Na]+", "[M-H]-", "[M+NH4]+", "[M+K]+"],
            n_features, p=[0.45, 0.15, 0.25, 0.10, 0.05]
        ),
    })
    features = features[features["intensity"] > noise_threshold].reset_index(drop=True)
    logger.info(f"Extracted {len(features)} features above noise threshold")
    return features


# ---------------------------------------------------------------------------
# XCMS-based extraction (R via rpy2)
# ---------------------------------------------------------------------------

def extract_features_xcms_r(mzml_path: str) -> pd.DataFrame:
    """
    XCMS (R) 経由のピーク検出 — rpy2 が利用可能な場合のみ。
    CentWave → groupChromPeaks → adjustRtime → featureDefinitions
    """
    r_script = """
    library(xcms)
    library(MSnbase)

    # --- CentWave peak detection ---
    raw <- readMSData(files = "{mzml_path}", mode = "onDisk")
    cwp <- CentWaveParam(peakwidth = c(5, 30), noise = 1e4, snthresh = 5)
    xdata <- findChromPeaks(raw, param = cwp)

    # --- RT alignment (Obiwarp) ---
    xdata <- adjustRtime(xdata, param = ObiwarpParam())

    # --- Peak grouping ---
    pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(raw))),
                            bw = 5, minFraction = 0.5)
    xdata <- groupChromPeaks(xdata, param = pdp)

    # --- Feature table ---
    ft <- featureDefinitions(xdata)
    fv <- featureValues(xdata, value = "into")
    result <- cbind(ft, fv)
    write.csv(result, "xcms_features.csv")
    """.format(mzml_path=mzml_path)

    logger.info("XCMS R script prepared (requires rpy2 or direct R execution)")
    return r_script


# ---------------------------------------------------------------------------
# Adduct & isotope deconvolution
# ---------------------------------------------------------------------------

ADDUCT_MASSES = {
    "[M+H]+": 1.007276,
    "[M+Na]+": 22.989218,
    "[M+K]+": 38.963158,
    "[M+NH4]+": 18.034164,
    "[M-H]-": -1.007276,
    "[M+Cl]-": 34.969402,
    "[M-H2O+H]+": -17.002740,
}

def deconvolve_adducts(features: pd.DataFrame, mz_tol_ppm: float = 10.0) -> pd.DataFrame:
    """アダクト・同位体デコンボリューション"""
    logger.info("Running adduct/isotope deconvolution")
    features = features.copy()

    neutral_masses = []
    for _, row in features.iterrows():
        adduct_shift = ADDUCT_MASSES.get(row["adduct"], 1.007276)
        neutral = row["mz"] - adduct_shift
        neutral_masses.append(round(neutral, 4))

    features["neutral_mass"] = neutral_masses

    # Group by neutral mass within tolerance
    features = features.sort_values("neutral_mass").reset_index(drop=True)
    group_id = 0
    groups = [0]
    for i in range(1, len(features)):
        mass_diff_ppm = abs(features.loc[i, "neutral_mass"] - features.loc[i-1, "neutral_mass"]) / features.loc[i, "neutral_mass"] * 1e6
        if mass_diff_ppm > mz_tol_ppm:
            group_id += 1
        groups.append(group_id)
    features["deconv_group"] = groups

    n_groups = features["deconv_group"].nunique()
    logger.info(f"Deconvolved into {n_groups} unique compound groups")
    return features


# ---------------------------------------------------------------------------
# Database matching (HMDB, KEGG, MassBank)
# ---------------------------------------------------------------------------

def match_databases(features: pd.DataFrame, mz_tol_ppm: float = 5.0) -> pd.DataFrame:
    """
    質量データベースマッチング (シミュレーション)
    本番では HMDB REST API / KEGG API / local SQLite DB を使用
    """
    logger.info("Matching against metabolite databases (HMDB, KEGG, MassBank)")

    np.random.seed(123)
    n = len(features)

    known_metabolites = [
        ("Butyrate", "C4H8O2", "INCHI_BUTYRATE", "HMDB0000039", "C00246"),
        ("Propionate", "C3H6O2", "INCHI_PROPIONATE", "HMDB0000237", "C00163"),
        ("Tryptophan", "C11H12N2O2", "INCHI_TRP", "HMDB0000929", "C00078"),
        ("Indole-3-propionic acid", "C11H11NO2", "INCHI_IPA", "HMDB0002302", "C02693"),
        ("Bile acid (CDCA)", "C24H40O4", "INCHI_CDCA", "HMDB0000518", "C02691"),
        ("Trimethylamine N-oxide", "C3H9NO", "INCHI_TMAO", "HMDB0000925", "C01104"),
        ("p-Cresol sulfate", "C7H8O4S", "INCHI_PCS", "HMDB0011635", "C01468"),
        ("Hippuric acid", "C9H9NO3", "INCHI_HIP", "HMDB0000714", "C01586"),
        ("Phenylacetic acid", "C8H8O2", "INCHI_PAA", "HMDB0000209", "C07086"),
        ("Succinate", "C4H6O4", "INCHI_SUC", "HMDB0000254", "C00042"),
    ]

    annotations = []
    for _, row in features.iterrows():
        if np.random.random() < 0.35:
            met = known_metabolites[np.random.randint(0, len(known_metabolites))]
            ann = {
                "feature_id": row["feature_id"],
                "compound_name": met[0],
                "formula": met[1],
                "inchikey": met[2],
                "hmdb_id": met[3],
                "kegg_id": met[4],
                "mass_error_ppm": round(np.random.uniform(0.1, mz_tol_ppm), 2),
                "ms2_cosine_score": round(np.random.uniform(0.5, 0.99), 3),
                "msi_level": np.random.choice([1, 2, 3, 4], p=[0.05, 0.25, 0.40, 0.30]),
            }
            annotations.append(ann)

    ann_df = pd.DataFrame(annotations)
    logger.info(f"Annotated {len(ann_df)} features ({len(ann_df)/n*100:.1f}%)")
    return ann_df


# ---------------------------------------------------------------------------
# MS2 spectral matching
# ---------------------------------------------------------------------------

def cosine_similarity(spec_a: np.ndarray, spec_b: np.ndarray) -> float:
    """Modified cosine similarity for MS2 spectra"""
    if spec_a.size == 0 or spec_b.size == 0:
        return 0.0
    dot = np.dot(spec_a, spec_b)
    norm_a = np.linalg.norm(spec_a)
    norm_b = np.linalg.norm(spec_b)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return float(dot / (norm_a * norm_b))


def assign_msi_levels(annotations: pd.DataFrame) -> pd.DataFrame:
    """
    MSI confidence levels:
      Level 1: RT + MS2 match to authentic standard
      Level 2: MS2 match to spectral library
      Level 3: Putative class (accurate mass + chemical class)
      Level 4: Unknown (accurate mass only)
    """
    logger.info("Assigning MSI confidence levels")
    ann = annotations.copy()
    ann.loc[ann["ms2_cosine_score"] >= 0.9, "msi_level"] = 1
    ann.loc[(ann["ms2_cosine_score"] >= 0.7) & (ann["ms2_cosine_score"] < 0.9), "msi_level"] = 2
    ann.loc[(ann["ms2_cosine_score"] >= 0.5) & (ann["ms2_cosine_score"] < 0.7), "msi_level"] = 3
    ann.loc[ann["ms2_cosine_score"] < 0.5, "msi_level"] = 4
    return ann


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_peak_annotation_pipeline(mzml_path: str = "simulated", output_dir: str = "results") -> dict:
    """メインパイプライン実行"""
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Feature extraction
    features = extract_features_pyopenms(mzml_path)

    # Step 2: Adduct/isotope deconvolution
    features = deconvolve_adducts(features)

    # Step 3: Database matching
    annotations = match_databases(features)

    # Step 4: MSI level assignment
    annotations = assign_msi_levels(annotations)

    # Save results
    features.to_csv(os.path.join(output_dir, "peak_features.csv"), index=False)
    annotations.to_csv(os.path.join(output_dir, "annotations.csv"), index=False)

    summary = {
        "total_features": len(features),
        "annotated_features": len(annotations),
        "annotation_rate": round(len(annotations) / len(features) * 100, 1),
        "msi_level_distribution": annotations["msi_level"].value_counts().to_dict(),
        "unique_compounds": annotations["compound_name"].nunique(),
        "deconv_groups": features["deconv_group"].nunique(),
    }

    with open(os.path.join(output_dir, "peak_annotation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    logger.info(f"Pipeline complete: {summary}")
    return summary


if __name__ == "__main__":
    summary = run_peak_annotation_pipeline(
        output_dir="../results"
    )
    print(json.dumps(summary, indent=2))
