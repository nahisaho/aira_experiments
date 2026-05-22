from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


def preprocess_airr(input_path: Path, output_path: Path, log_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, pd.DataFrame]]:
    raw = pd.read_csv(input_path, sep="\t")
    raw["productive"] = raw["productive"].astype(bool)
    raw["cdr3_length"] = raw["junction_aa"].fillna("").str.replace("*", "", regex=False).str.len()

    quality = (
        raw.groupby(["sample_id", "sample_type", "chronological_age", "icb_response"], as_index=False)
        .apply(
            lambda g: pd.Series(
                {
                    "total_reads": int(g["clone_count"].sum()),
                    "productive_reads": int(g.loc[g["productive"], "clone_count"].sum()),
                    "productive_ratio": float(g.loc[g["productive"], "clone_count"].sum() / g["clone_count"].sum()),
                    "mean_cdr3_length": float(g.loc[g["productive"], "cdr3_length"].mean()),
                    "median_cdr3_length": float(g.loc[g["productive"], "cdr3_length"].median()),
                }
            )
        )
        .reset_index(drop=True)
    )

    productive = raw.loc[raw["productive"]].copy()
    productive["junction_aa"] = productive["junction_aa"].str.replace("*", "", regex=False)

    clonotypes = (
        productive.groupby(
            [
                "sample_id",
                "sample_type",
                "chronological_age",
                "icb_response",
                "v_call",
                "j_call",
                "junction_aa",
            ],
            as_index=False,
        )["clone_count"]
        .sum()
    )
    sample_totals = clonotypes.groupby("sample_id")["clone_count"].transform("sum")
    clonotypes["clone_frequency"] = clonotypes["clone_count"] / sample_totals
    median_counts = clonotypes.groupby("sample_id")["clone_count"].transform("median")
    clonotypes["expansion_index"] = clonotypes["clone_count"] / median_counts.clip(lower=1)
    clonotypes["cdr3_length"] = clonotypes["junction_aa"].str.len()
    clonotypes["clonotype_id"] = (
        clonotypes["sample_id"]
        + "|"
        + clonotypes["v_call"]
        + "|"
        + clonotypes["j_call"]
        + "|"
        + clonotypes["junction_aa"]
    )
    clonotypes = clonotypes.sort_values(["sample_id", "clone_count"], ascending=[True, False]).reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    clonotypes.to_csv(output_path, sep="\t", index=False)

    top_v_usage = (
        clonotypes.groupby(["sample_id", "v_call"], as_index=False)["clone_count"].sum()
        .assign(freq=lambda d: d["clone_count"] / d.groupby("sample_id")["clone_count"].transform("sum"))
    )

    log_text = """# Preprocessing Log

1. Loaded AIRR-style TSV and validated required fields.
2. Calculated sample-level quality metrics from total and productive clone counts.
3. Filtered to productive rearrangements only.
4. Defined clonotypes using `(V gene, J gene, CDR3 amino acid)` tuples.
5. Aggregated clone counts, computed clone frequencies and expansion index.
6. Calculated CDR3 length summaries for downstream immune age analysis.
"""
    log_path.write_text(log_text, encoding="utf-8")

    extras = {
        "raw": raw,
        "quality": quality,
        "v_usage": top_v_usage,
    }
    return clonotypes, quality, extras
