from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

SEED = 42
RNG = np.random.default_rng(SEED)

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
CODON_MAP = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"
}

TRBV_FREQS = {
    "TRBV20-1": 0.12,
    "TRBV28": 0.08,
    "TRBV5-1": 0.07,
    "TRBV7-2": 0.06,
    "TRBV12-3": 0.05,
    "TRBV27": 0.05,
    "TRBV6-5": 0.04,
    "TRBV29-1": 0.04,
    "TRBV18": 0.03,
    "TRBV4-1": 0.05,
    "TRBV6-1": 0.04,
    "TRBV7-9": 0.04,
    "TRBV9": 0.04,
    "TRBV14": 0.04,
    "TRBV19": 0.04,
    "TRBV24-1": 0.04,
    "TRBV30": 0.03,
    "TRBV10-3": 0.03,
    "TRBV11-2": 0.03,
    "TRBV15": 0.02,
}

TRBJ_FREQS = {
    "TRBJ1-1": 0.12, "TRBJ1-2": 0.10, "TRBJ1-4": 0.09, "TRBJ1-5": 0.08,
    "TRBJ1-6": 0.09, "TRBJ2-1": 0.11, "TRBJ2-3": 0.10, "TRBJ2-5": 0.12,
    "TRBJ2-7": 0.10, "TRBJ2-2": 0.09,
}

TRBD_GENES = ["TRBD1", "TRBD2", "TRBD2*02"]

SAMPLES = [
    {"sample_id": "HD01", "sample_type": "healthy", "chronological_age": 24, "icb_response": 0},
    {"sample_id": "HD02", "sample_type": "healthy", "chronological_age": 31, "icb_response": 0},
    {"sample_id": "HD03", "sample_type": "healthy", "chronological_age": 39, "icb_response": 0},
    {"sample_id": "HD04", "sample_type": "healthy", "chronological_age": 47, "icb_response": 0},
    {"sample_id": "HD05", "sample_type": "healthy", "chronological_age": 55, "icb_response": 0},
    {"sample_id": "CA01", "sample_type": "cancer", "chronological_age": 58, "icb_response": 0},
    {"sample_id": "CA02", "sample_type": "cancer", "chronological_age": 64, "icb_response": 0},
    {"sample_id": "CA03", "sample_type": "cancer", "chronological_age": 69, "icb_response": 0},
    {"sample_id": "IR01", "sample_type": "icb_responder", "chronological_age": 49, "icb_response": 1},
    {"sample_id": "IR02", "sample_type": "icb_responder", "chronological_age": 61, "icb_response": 1},
]

PUBLIC_TCRS = [
    {"v_call": "TRBV20-1", "j_call": "TRBJ1-2", "junction_aa": "CASSIRSSYEQYF", "antigen": "CMV pp65", "hla": "HLA-A*02:01"},
    {"v_call": "TRBV20-1", "j_call": "TRBJ1-2", "junction_aa": "CASSPPTGELF", "antigen": "CMV pp65", "hla": "HLA-B*07:02"},
    {"v_call": "TRBV27", "j_call": "TRBJ1-1", "junction_aa": "CASSLEGQYF", "antigen": "Flu M1", "hla": "HLA-A*02:01"},
    {"v_call": "TRBV20-1", "j_call": "TRBJ2-7", "junction_aa": "CATSDRLAGGYNEQYF", "antigen": "EBV BMLF", "hla": "HLA-A*24:02"},
    {"v_call": "TRBV12-3", "j_call": "TRBJ2-1", "junction_aa": "CASSIGTGELF", "antigen": "COVID Spike", "hla": "HLA-A*02:01"},
    {"v_call": "TRBV7-2", "j_call": "TRBJ2-5", "junction_aa": "CASSLGQNTLYF", "antigen": "MART1", "hla": "HLA-A*02:01"},
    {"v_call": "TRBV5-1", "j_call": "TRBJ2-3", "junction_aa": "CASSYVGNTIYF", "antigen": "NYESO1", "hla": "HLA-A*24:02"},
    {"v_call": "TRBV6-5", "j_call": "TRBJ2-5", "junction_aa": "CASSQETQYF", "antigen": "MAGEA3", "hla": "HLA-B*07:02"},
]


def _weighted_choice(weights: Dict[str, float], size: int) -> np.ndarray:
    keys = np.array(list(weights.keys()))
    probs = np.array(list(weights.values()), dtype=float)
    probs = probs / probs.sum()
    return RNG.choice(keys, size=size, p=probs)


def _random_cdr3(length: int) -> str:
    core = "".join(RNG.choice(AMINO_ACIDS, size=max(length - 2, 1)))
    tail = "F" if RNG.random() < 0.75 else "YF"
    seq = f"C{core}{tail}"
    return seq[: max(length, 8)]


def _aa_to_nt(seq: str) -> str:
    return "".join(CODON_MAP.get(aa, "NNN") for aa in seq)


def _sample_reads(sample_type: str) -> int:
    if sample_type == "healthy":
        return int(RNG.integers(9000, 15001))
    if sample_type == "cancer":
        return int(RNG.integers(6000, 12001))
    return int(RNG.integers(12000, 19001))


def _sample_clones(sample_type: str) -> int:
    if sample_type == "healthy":
        return int(RNG.integers(650, 1000))
    if sample_type == "cancer":
        return int(RNG.integers(280, 520))
    return int(RNG.integers(720, 1100))


def _generate_clone_counts(n_clones: int, total_reads: int, sample_type: str) -> np.ndarray:
    alpha = {"healthy": 1.3, "cancer": 0.45, "icb_responder": 1.1}[sample_type]
    weights = RNG.dirichlet(np.full(n_clones, alpha))
    if sample_type == "cancer":
        boosted = RNG.choice(np.arange(n_clones), size=8, replace=False)
        weights[boosted] *= RNG.uniform(3.0, 8.0, size=8)
    elif sample_type == "icb_responder":
        boosted = RNG.choice(np.arange(n_clones), size=10, replace=False)
        weights[boosted] *= RNG.uniform(1.5, 3.0, size=10)
    weights = weights / weights.sum()
    counts = np.maximum(1, np.floor(weights * total_reads).astype(int))
    diff = total_reads - counts.sum()
    while diff != 0:
        idx = int(RNG.integers(0, n_clones))
        if diff > 0:
            counts[idx] += 1
            diff -= 1
        elif counts[idx] > 1:
            counts[idx] -= 1
            diff += 1
    return counts


def _inject_public_tcrs(records: List[dict], sample: dict, total_reads: int) -> None:
    n_public = 2 if sample["sample_type"] == "healthy" else 4 if sample["sample_type"] == "cancer" else 5
    chosen = RNG.choice(len(PUBLIC_TCRS), size=n_public, replace=False)
    for idx in chosen:
        entry = PUBLIC_TCRS[idx]
        base_count = int(total_reads * RNG.uniform(0.005, 0.03))
        if sample["sample_type"] == "cancer":
            base_count = int(base_count * RNG.uniform(1.5, 3.0))
        if sample["sample_type"] == "icb_responder" and entry["antigen"] in {"MART1", "NYESO1", "MAGEA3"}:
            base_count = int(base_count * RNG.uniform(2.0, 3.5))
        records.append(
            {
                "sequence_id": f"{sample['sample_id']}_PUB_{entry['junction_aa']}",
                "sample_id": sample["sample_id"],
                "sample_type": sample["sample_type"],
                "chronological_age": sample["chronological_age"],
                "icb_response": sample["icb_response"],
                "v_call": entry["v_call"],
                "d_call": RNG.choice(TRBD_GENES),
                "j_call": entry["j_call"],
                "junction": _aa_to_nt(entry["junction_aa"]),
                "junction_aa": entry["junction_aa"],
                "productive": True,
                "clone_count": max(base_count, 10),
            }
        )


def simulate_airr_dataset(output_path: Path) -> pd.DataFrame:
    rows: List[dict] = []
    for sample in SAMPLES:
        total_reads = _sample_reads(sample["sample_type"])
        n_clones = _sample_clones(sample["sample_type"])
        v_calls = _weighted_choice(TRBV_FREQS, n_clones)
        j_calls = _weighted_choice(TRBJ_FREQS, n_clones)
        d_calls = RNG.choice(TRBD_GENES, size=n_clones, replace=True)
        counts = _generate_clone_counts(n_clones, total_reads, sample["sample_type"])
        productive_prob = {"healthy": 0.92, "cancer": 0.86, "icb_responder": 0.9}[sample["sample_type"]]

        for i in range(n_clones):
            cdr3_len = int(RNG.integers(10, 18))
            junction_aa = _random_cdr3(cdr3_len)
            productive = bool(RNG.random() < productive_prob)
            if not productive and not junction_aa.endswith("*"):
                junction_aa = junction_aa[:-1] + "*"
            rows.append(
                {
                    "sequence_id": f"{sample['sample_id']}_{i:05d}",
                    "sample_id": sample["sample_id"],
                    "sample_type": sample["sample_type"],
                    "chronological_age": sample["chronological_age"],
                    "icb_response": sample["icb_response"],
                    "v_call": v_calls[i],
                    "d_call": d_calls[i],
                    "j_call": j_calls[i],
                    "junction": _aa_to_nt(junction_aa.replace("*", "F")),
                    "junction_aa": junction_aa,
                    "productive": productive,
                    "clone_count": int(counts[i]),
                }
            )
        _inject_public_tcrs(rows, sample, total_reads)

    df = pd.DataFrame(rows)
    df = df.sort_values(["sample_id", "clone_count"], ascending=[True, False]).reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


if __name__ == "__main__":
    simulate_airr_dataset(Path(__file__).resolve().parent / "simulated_tcr_seq.tsv")
