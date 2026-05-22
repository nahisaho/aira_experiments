from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd


BASE_REFERENCE = [
    ("TRBV20-1", "CASSIRSSYEQYF", "CMV pp65", "HLA-A*02:01"),
    ("TRBV20-1", "CASSPPTGELF", "CMV pp65", "HLA-B*07:02"),
    ("TRBV27", "CASSLEGQYF", "Flu M1", "HLA-A*02:01"),
    ("TRBV20-1", "CATSDRLAGGYNEQYF", "EBV BMLF", "HLA-A*24:02"),
    ("TRBV12-3", "CASSIGTGELF", "COVID Spike", "HLA-A*02:01"),
    ("TRBV7-2", "CASSLGQNTLYF", "MART1", "HLA-A*02:01"),
    ("TRBV5-1", "CASSYVGNTIYF", "NYESO1", "HLA-A*24:02"),
    ("TRBV6-5", "CASSQETQYF", "MAGEA3", "HLA-B*07:02"),
    ("TRBV28", "CASSQGANTEAFF", "CMV IE1", "HLA-B*07:02"),
    ("TRBV29-1", "CASSLGVNEKLFF", "Flu NP", "HLA-A*24:02"),
]

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
RNG = np.random.default_rng(42)


def _mutate(seq: str) -> str:
    if len(seq) < 5:
        return seq
    seq_list = list(seq)
    pos = int(RNG.integers(1, len(seq_list) - 1))
    seq_list[pos] = RNG.choice(AMINO_ACIDS)
    return "".join(seq_list)


def build_vdjdb_reference() -> pd.DataFrame:
    entries = []
    for v_call, cdr3, antigen, hla in BASE_REFERENCE:
        entries.append({"v_call": v_call, "junction_aa": cdr3, "antigen": antigen, "hla_restriction": hla})
        for _ in range(5):
            entries.append(
                {
                    "v_call": v_call,
                    "junction_aa": _mutate(cdr3),
                    "antigen": antigen,
                    "hla_restriction": hla,
                }
            )
    ref = pd.DataFrame(entries).drop_duplicates().reset_index(drop=True)
    ref["reference_id"] = [f"VDJDB_{i:03d}" for i in range(1, len(ref) + 1)]
    return ref


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        curr = [i]
        for j, cb in enumerate(b, start=1):
            ins = curr[j - 1] + 1
            delete = prev[j] + 1
            sub = prev[j - 1] + (ca != cb)
            curr.append(min(ins, delete, sub))
        prev = curr
    return prev[-1]


def _identity(a: str, b: str) -> float:
    return 1.0 - levenshtein(a, b) / max(len(a), len(b))


def identify_public_tcrs(clonotypes: pd.DataFrame, output_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    reference = build_vdjdb_reference()
    matches: List[dict] = []
    for row in clonotypes.itertuples(index=False):
        candidates = reference.loc[reference["v_call"] == row.v_call]
        best = None
        best_identity = 0.0
        for ref in candidates.itertuples(index=False):
            ident = _identity(row.junction_aa, ref.junction_aa)
            if ident > best_identity:
                best_identity = ident
                best = ref
        if best is not None and best_identity >= 0.8:
            confidence = round(min(0.99, 0.65 * best_identity + 0.35 * min(row.clone_frequency * 6, 1.0)), 3)
            matches.append(
                {
                    "sample_id": row.sample_id,
                    "sample_type": row.sample_type,
                    "v_call": row.v_call,
                    "j_call": row.j_call,
                    "junction_aa": row.junction_aa,
                    "clone_count": int(row.clone_count),
                    "clone_frequency": float(row.clone_frequency),
                    "matched_reference": best.reference_id,
                    "matched_cdr3": best.junction_aa,
                    "identity": round(best_identity, 3),
                    "antigen": best.antigen,
                    "hla_restriction": best.hla_restriction,
                    "confidence_score": confidence,
                    "tumor_reactive": int(best.antigen in {"MART1", "NYESO1", "MAGEA3"}),
                }
            )
    public_df = pd.DataFrame(matches).sort_values(["sample_id", "confidence_score"], ascending=[True, False]).reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    public_df.to_csv(output_path, sep="\t", index=False)
    return public_df, reference
