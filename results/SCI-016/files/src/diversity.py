from __future__ import annotations

from math import exp, log
from pathlib import Path

import numpy as np
import pandas as pd


def _shannon(counts: np.ndarray) -> tuple[float, float]:
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum()), float(-(p * np.log(p)).sum())


def _chao1(counts: np.ndarray) -> float:
    s_obs = float((counts > 0).sum())
    n1 = float((counts == 1).sum())
    n2 = float((counts == 2).sum())
    if n2 == 0:
        return s_obs + (n1 * (n1 - 1)) / 2
    return s_obs + (n1**2) / (2 * n2)


def _d50(counts: np.ndarray) -> int:
    ordered = np.sort(counts)[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    return int(np.argmax(cumulative >= 0.5) + 1)


def calculate_diversity_metrics(clonotypes: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    records = []
    for sample_id, group in clonotypes.groupby("sample_id"):
        counts = group["clone_count"].to_numpy(dtype=float)
        shannon_log2, shannon_nat = _shannon(counts)
        richness = int((counts > 0).sum())
        simpson = float((counts / counts.sum()) ** 2).sum()
        gini_simpson = 1.0 - simpson
        pielou = float(shannon_nat / log(richness)) if richness > 1 else 0.0
        top10_freq = float(group.nlargest(10, "clone_frequency")["clone_frequency"].sum())
        records.append(
            {
                "sample_id": sample_id,
                "sample_type": group["sample_type"].iloc[0],
                "chronological_age": int(group["chronological_age"].iloc[0]),
                "icb_response": int(group["icb_response"].iloc[0]),
                "richness": richness,
                "shannon_entropy": shannon_log2,
                "chao1": _chao1(counts),
                "hill_q0": float(richness),
                "hill_q1": float(exp(shannon_nat)),
                "hill_q2": float(1.0 / simpson) if simpson > 0 else 0.0,
                "gini_simpson": gini_simpson,
                "pielou_evenness": pielou,
                "d50": _d50(counts),
                "top10_clone_frequency": top10_freq,
                "top1_clone_frequency": float(group["clone_frequency"].max()),
                "singleton_ratio": float((group["clone_count"] == 1).mean()),
                "mean_cdr3_length": float(group["cdr3_length"].mean()),
                "clone_expansion_index": float(group["expansion_index"].mean()),
            }
        )
    diversity = pd.DataFrame(records).sort_values("sample_id").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    diversity.to_csv(output_path, sep="\t", index=False)
    return diversity
