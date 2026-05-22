"""
CRISPR-Cas9 Off-Target Prediction — Preprocessing Pipeline
Supports GUIDE-seq and CIRCLE-seq data formats.
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Optional
import re
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# ─── Constants ────────────────────────────────────────────────────────────────
NUC_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4, "-": 5}
SEQ_LEN = 23  # 20-nt protospacer + 3-nt PAM
EPIGENETIC_DIM = 8  # ATAC-seq bins (4) + CpG methylation bins (4)

MISMATCH_TYPES = {
    ("A", "C"): 0, ("A", "G"): 1, ("A", "T"): 2,
    ("C", "A"): 3, ("C", "G"): 4, ("C", "T"): 5,
    ("G", "A"): 6, ("G", "C"): 7, ("G", "T"): 8,
    ("T", "A"): 9, ("T", "C"): 10, ("T", "G"): 11,
    ("match", "match"): 12,
    ("DNA_bulge", "DNA_bulge"): 13,
    ("RNA_bulge", "RNA_bulge"): 14,
}


# ─── Sequence Encoding ────────────────────────────────────────────────────────

def one_hot_encode(seq: str, length: int = SEQ_LEN) -> np.ndarray:
    """One-hot encode a nucleotide sequence → (length, 4) float32 array."""
    seq = seq.upper().ljust(length, "N")[:length]
    arr = np.zeros((length, 4), dtype=np.float32)
    for i, nt in enumerate(seq):
        idx = NUC_MAP.get(nt, 4)
        if idx < 4:
            arr[i, idx] = 1.0
    return arr


def encode_mismatch_pattern(
    guide: str, target: str, length: int = SEQ_LEN
) -> np.ndarray:
    """
    Encode mismatch pattern between guide RNA and genomic target sequence.
    Returns (length, 15) float32 array: one-hot over 15 mismatch/match types.
    """
    guide = guide.upper().ljust(length, "N")[:length]
    target = target.upper().ljust(length, "N")[:length]
    arr = np.zeros((length, len(MISMATCH_TYPES)), dtype=np.float32)

    for i, (g, t) in enumerate(zip(guide, target)):
        if g == "-":
            key = ("RNA_bulge", "RNA_bulge")
        elif t == "-":
            key = ("DNA_bulge", "DNA_bulge")
        elif g == t:
            key = ("match", "match")
        else:
            key = (g, t)
        idx = MISMATCH_TYPES.get(key, 12)
        arr[i, idx] = 1.0
    return arr


def positional_mismatch_vector(guide: str, target: str) -> np.ndarray:
    """
    Binary mismatch indicator per position (length=23).
    Positions closer to PAM (seed region) are weighted ×2.
    """
    vec = np.zeros(SEQ_LEN, dtype=np.float32)
    guide = guide.upper().ljust(SEQ_LEN, "N")[:SEQ_LEN]
    target = target.upper().ljust(SEQ_LEN, "N")[:SEQ_LEN]
    for i, (g, t) in enumerate(zip(guide, target)):
        if g != t and g != "N" and t != "N":
            # seed region = last 12 nt before PAM (positions 8–19)
            weight = 2.0 if 8 <= i <= 19 else 1.0
            vec[i] = weight
    return vec


# ─── Epigenetic Feature Encoding ─────────────────────────────────────────────

def encode_epigenetics(
    atac_signal: Optional[np.ndarray] = None,
    methylation: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Encode epigenetic features into a fixed-length (EPIGENETIC_DIM,) vector.
    atac_signal: log-normalised ATAC-seq read counts over site window.
    methylation: CpG methylation fraction (0–1) per position.
    Returns zeros when data is absent (missing-data safe).
    """
    feat = np.zeros(EPIGENETIC_DIM, dtype=np.float32)
    if atac_signal is not None:
        sig = np.array(atac_signal, dtype=np.float32)
        sig = np.log1p(sig)
        # 4 percentile bins: min, 33rd, 66th, max
        feat[0] = sig.min()
        feat[1] = np.percentile(sig, 33)
        feat[2] = np.percentile(sig, 66)
        feat[3] = sig.max()
    if methylation is not None:
        meth = np.array(methylation, dtype=np.float32)
        feat[4] = meth.mean()
        feat[5] = meth.std()
        feat[6] = (meth > 0.5).mean()   # fraction of hypermethylated CpGs
        feat[7] = (meth < 0.1).mean()   # fraction of unmethylated CpGs
    return feat


# ─── GUIDE-seq / CIRCLE-seq Parsers ──────────────────────────────────────────

def parse_guide_seq(filepath: str) -> pd.DataFrame:
    """
    Parse a GUIDE-seq output table (BED-like TSV).
    Expected columns: chr, start, end, name, reads, strand,
                      guide_seq, target_seq, mismatches, label
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"GUIDE-seq file not found: {filepath}")

    df = pd.read_csv(filepath, sep="\t", comment="#")
    required = {"guide_seq", "target_seq", "reads", "mismatches"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GUIDE-seq file missing columns: {missing}")

    # Normalise sequence length
    df["guide_seq"] = df["guide_seq"].str.upper().str[:SEQ_LEN]
    df["target_seq"] = df["target_seq"].str.upper().str[:SEQ_LEN]

    # Binary label: reads > threshold → positive off-target site
    if "label" not in df.columns:
        threshold = df["reads"].quantile(0.25)
        df["label"] = (df["reads"] > threshold).astype(int)

    logger.info("GUIDE-seq: loaded %d sites from %s", len(df), filepath)
    return df


def parse_circle_seq(filepath: str) -> pd.DataFrame:
    """
    Parse a CIRCLE-seq output table.
    Expected columns: chr, start, end, guide_seq, target_seq,
                      read_count, mismatches, label
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"CIRCLE-seq file not found: {filepath}")

    df = pd.read_csv(filepath, sep="\t", comment="#")
    required = {"guide_seq", "target_seq", "read_count", "mismatches"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CIRCLE-seq file missing columns: {missing}")

    df["guide_seq"] = df["guide_seq"].str.upper().str[:SEQ_LEN]
    df["target_seq"] = df["target_seq"].str.upper().str[:SEQ_LEN]

    if "label" not in df.columns:
        threshold = df["read_count"].quantile(0.25)
        df["label"] = (df["read_count"] > threshold).astype(int)

    # Rename to unified column names
    df = df.rename(columns={"read_count": "reads"})
    logger.info("CIRCLE-seq: loaded %d sites from %s", len(df), filepath)
    return df


# ─── Feature Matrix Builder ───────────────────────────────────────────────────

class CRISPRFeatureBuilder:
    """
    Converts a DataFrame of off-target candidate sites into
    multi-channel tensors ready for CNN input.

    Output tensor shape per sample:
        sequence_channels : (SEQ_LEN, 4+4+15) = (23, 23) — guide + target OH + mismatch
        scalar_features   : (SEQ_LEN + EPIGENETIC_DIM,)  — positional + epigenetic
    """

    def __init__(self, include_epigenetics: bool = True):
        self.include_epigenetics = include_epigenetics

    def build_sequence_tensor(self, row: pd.Series) -> np.ndarray:
        guide_oh   = one_hot_encode(row["guide_seq"])        # (23, 4)
        target_oh  = one_hot_encode(row["target_seq"])       # (23, 4)
        mismatch   = encode_mismatch_pattern(row["guide_seq"], row["target_seq"])  # (23, 15)
        return np.concatenate([guide_oh, target_oh, mismatch], axis=1)  # (23, 23)

    def build_scalar_vector(self, row: pd.Series) -> np.ndarray:
        pos_mm = positional_mismatch_vector(row["guide_seq"], row["target_seq"])  # (23,)
        if self.include_epigenetics:
            atac  = row.get("atac_signal", None)
            meth  = row.get("methylation", None)
            epi   = encode_epigenetics(atac, meth)  # (8,)
            return np.concatenate([pos_mm, epi])   # (31,)
        return pos_mm  # (23,)

    def transform(
        self, df: pd.DataFrame
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns:
            X_seq    : (N, 23, 23) float32
            X_scalar : (N, 31) float32
            y        : (N,) int32
        """
        seq_list, scalar_list = [], []
        for _, row in df.iterrows():
            seq_list.append(self.build_sequence_tensor(row))
            scalar_list.append(self.build_scalar_vector(row))

        X_seq    = np.stack(seq_list).astype(np.float32)
        X_scalar = np.stack(scalar_list).astype(np.float32)
        y        = df["label"].values.astype(np.int32)

        logger.info(
            "Feature matrix built: X_seq=%s X_scalar=%s y=%s",
            X_seq.shape, X_scalar.shape, y.shape,
        )
        return X_seq, X_scalar, y


# ─── Synthetic Data Generator (for testing / benchmarking) ───────────────────

NUC = list("ACGT")

def _random_seq(length: int) -> str:
    return "".join(np.random.choice(NUC, length))

def _mutate_seq(seq: str, n_mismatches: int) -> str:
    positions = np.random.choice(len(seq) - 3, n_mismatches, replace=False)
    seq = list(seq)
    for p in positions:
        others = [n for n in NUC if n != seq[p]]
        seq[p] = np.random.choice(others)
    return "".join(seq)


def generate_synthetic_dataset(
    n_guides: int = 20,
    sites_per_guide: int = 50,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generate a synthetic CRISPR off-target dataset for unit tests and demos.
    Label logic: 0–1 mismatches → positive (1); ≥4 mismatches → negative (0).
    """
    np.random.seed(seed)
    rows = []
    for _ in range(n_guides):
        guide = _random_seq(SEQ_LEN)
        for _ in range(sites_per_guide):
            n_mm = np.random.choice([0, 1, 2, 3, 4, 5], p=[0.05, 0.15, 0.25, 0.25, 0.2, 0.1])
            target = _mutate_seq(guide, n_mm)
            label  = 1 if n_mm <= 2 else 0
            reads  = int(np.random.lognormal(6, 1.5)) if label else int(np.random.lognormal(2, 1))
            rows.append({
                "guide_seq":    guide,
                "target_seq":   target,
                "mismatches":   n_mm,
                "reads":        reads,
                "label":        label,
                "atac_signal":  None,
                "methylation":  None,
            })
    df = pd.DataFrame(rows)
    logger.info("Generated synthetic dataset: %d rows, pos=%d neg=%d",
                len(df), df["label"].sum(), (df["label"] == 0).sum())
    return df


if __name__ == "__main__":
    df = generate_synthetic_dataset()
    builder = CRISPRFeatureBuilder(include_epigenetics=True)
    X_seq, X_scalar, y = builder.transform(df)
    print(f"X_seq: {X_seq.shape}, X_scalar: {X_scalar.shape}, y: {y.shape}")
