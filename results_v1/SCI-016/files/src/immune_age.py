from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _minmax(series: pd.Series) -> pd.Series:
    denom = series.max() - series.min()
    if denom == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - series.min()) / denom


def estimate_immune_age(
    diversity: pd.DataFrame,
    clonotypes: pd.DataFrame,
    public_tcrs: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    singleton_counts = clonotypes.assign(singleton=(clonotypes["clone_count"] == 1).astype(int)).groupby("sample_id")["singleton"].mean()
    public_counts = public_tcrs.groupby("sample_id").size() if not public_tcrs.empty else pd.Series(dtype=float)
    mean_cdr3 = clonotypes.groupby("sample_id")["cdr3_length"].mean()

    immune = diversity.copy()
    immune["singleton_ratio"] = immune["sample_id"].map(singleton_counts).fillna(0)
    immune["public_tcr_count"] = immune["sample_id"].map(public_counts).fillna(0).astype(int)
    immune["mean_cdr3_length"] = immune["sample_id"].map(mean_cdr3).fillna(immune["mean_cdr3_length"])

    score_components = pd.DataFrame(
        {
            "age_component": _minmax(immune["chronological_age"]),
            "diversity_loss": 1 - _minmax(immune["shannon_entropy"]),
            "memory_bias": 1 - immune["singleton_ratio"].clip(0, 1),
            "clonal_expansion": _minmax(immune["top10_clone_frequency"]),
            "public_burden": _minmax(immune["public_tcr_count"]),
            "cdr3_shift": _minmax((immune["mean_cdr3_length"] - 14).abs()),
        }
    )
    weights = {
        "age_component": 0.30,
        "diversity_loss": 0.20,
        "memory_bias": 0.20,
        "clonal_expansion": 0.20,
        "public_burden": 0.05,
        "cdr3_shift": 0.05,
    }
    raw_score = sum(score_components[col] * weight for col, weight in weights.items())
    immune["immune_age_score"] = (20 + 75 * raw_score).round(2)
    immune["naive_memory_ratio_proxy"] = ((immune["singleton_ratio"] + 1e-3) / (1 - immune["singleton_ratio"] + 1e-3)).round(3)
    immune["immunologically_aged"] = immune["immune_age_score"] > (immune["chronological_age"] + 10)

    cols = [
        "sample_id",
        "sample_type",
        "chronological_age",
        "immune_age_score",
        "immunologically_aged",
        "mean_cdr3_length",
        "singleton_ratio",
        "naive_memory_ratio_proxy",
        "shannon_entropy",
        "top10_clone_frequency",
        "public_tcr_count",
    ]
    immune[cols].to_csv(output_path, sep="\t", index=False)
    return immune[cols]
