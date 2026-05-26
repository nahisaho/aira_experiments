"""
CRISPR-Cas9 Off-Target Prediction: Data Preprocessing Pipeline
Handles GUIDE-seq and CIRCLE-seq data preprocessing, feature encoding,
and epigenetic feature integration.
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict, List, Optional
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
import os

# Nucleotide encoding
NUC_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
NUC_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

# Mismatch type encoding (12 possible + match)
MISMATCH_TYPES = [
    'AA', 'AC', 'AG', 'AT',
    'CA', 'CC', 'CG', 'CT',
    'GA', 'GC', 'GG', 'GT',
    'TA', 'TC', 'TG', 'TT'
]


def one_hot_encode_sequence(seq: str, max_len: int = 23) -> np.ndarray:
    """One-hot encode a DNA sequence."""
    encoding = np.zeros((max_len, 4), dtype=np.float32)
    for i, nuc in enumerate(seq[:max_len]):
        if nuc in NUC_MAP and NUC_MAP[nuc] < 4:
            encoding[i, NUC_MAP[nuc]] = 1.0
    return encoding


def encode_mismatch_matrix(grna: str, target: str, max_len: int = 23) -> np.ndarray:
    """
    Create mismatch encoding matrix between gRNA and target.
    Returns a (max_len, 16) binary matrix for mismatch types.
    """
    matrix = np.zeros((max_len, 16), dtype=np.float32)
    for i in range(min(len(grna), len(target), max_len)):
        pair = grna[i] + target[i]
        if pair in MISMATCH_TYPES:
            idx = MISMATCH_TYPES.index(pair)
            matrix[i, idx] = 1.0
    return matrix


def encode_position_features(grna: str, target: str, max_len: int = 23) -> np.ndarray:
    """
    Encode positional mismatch features:
    - Binary mismatch indicator per position
    - Distance from PAM (normalized)
    - Consecutive mismatch count
    """
    features = np.zeros((max_len, 3), dtype=np.float32)
    min_len = min(len(grna), len(target), max_len)

    # Binary mismatch indicator
    for i in range(min_len):
        if grna[i] != target[i]:
            features[i, 0] = 1.0

    # Normalized distance from PAM (PAM is at 3' end, positions 21-23)
    for i in range(max_len):
        features[i, 1] = 1.0 - (i / max_len)

    # Consecutive mismatch count
    consec = 0
    for i in range(min_len):
        if grna[i] != target[i]:
            consec += 1
        else:
            consec = 0
        features[i, 2] = min(consec / 5.0, 1.0)

    return features


def encode_epigenetic_features(
    chromatin_accessibility: Optional[np.ndarray] = None,
    methylation: Optional[np.ndarray] = None,
    histone_h3k4me3: Optional[np.ndarray] = None,
    histone_h3k27ac: Optional[np.ndarray] = None,
    max_len: int = 23
) -> np.ndarray:
    """
    Encode epigenetic features for each position.
    Each feature is a signal value at the corresponding genomic position.
    Returns: (max_len, n_epi_features) array
    """
    epi_features = []

    if chromatin_accessibility is not None:
        ca = np.zeros(max_len, dtype=np.float32)
        ca[:len(chromatin_accessibility)] = chromatin_accessibility[:max_len]
        epi_features.append(ca)
    
    if methylation is not None:
        meth = np.zeros(max_len, dtype=np.float32)
        meth[:len(methylation)] = methylation[:max_len]
        epi_features.append(meth)
    
    if histone_h3k4me3 is not None:
        h3k4 = np.zeros(max_len, dtype=np.float32)
        h3k4[:len(histone_h3k4me3)] = histone_h3k4me3[:max_len]
        epi_features.append(h3k4)
    
    if histone_h3k27ac is not None:
        h3k27 = np.zeros(max_len, dtype=np.float32)
        h3k27[:len(histone_h3k27ac)] = histone_h3k27ac[:max_len]
        epi_features.append(h3k27)

    if not epi_features:
        return np.zeros((max_len, 4), dtype=np.float32)

    return np.stack(epi_features, axis=-1)


def generate_synthetic_dataset(
    n_samples: int = 10000,
    seq_len: int = 23,
    positive_ratio: float = 0.1,
    seed: int = 42
) -> Dict:
    """
    Generate synthetic CRISPR off-target dataset mimicking GUIDE-seq/CIRCLE-seq.
    Positive samples have mismatches correlated with off-target activity.
    """
    rng = np.random.RandomState(seed)
    nucs = list('ACGT')

    grna_seqs = []
    target_seqs = []
    labels = []
    epi_data = {
        'chromatin_accessibility': [],
        'methylation': [],
        'h3k4me3': [],
        'h3k27ac': []
    }

    n_positive = int(n_samples * positive_ratio)
    n_negative = n_samples - n_positive

    for i in range(n_samples):
        grna = ''.join(rng.choice(nucs, seq_len))
        target = list(grna)

        if i < n_positive:
            # Positive: introduce 1-6 mismatches, biased toward seed region
            n_mm = rng.randint(1, 7)
            positions = rng.choice(seq_len, n_mm, replace=False)
            for pos in positions:
                alt_nucs = [n for n in nucs if n != target[pos]]
                target[pos] = rng.choice(alt_nucs)
            labels.append(1)
            # Higher chromatin accessibility for positive samples
            ca = rng.beta(5, 2, seq_len).astype(np.float32)
            meth = rng.beta(2, 5, seq_len).astype(np.float32)
        else:
            # Negative: more mismatches, especially in seed region
            n_mm = rng.randint(4, 10)
            positions = rng.choice(seq_len, min(n_mm, seq_len), replace=False)
            for pos in positions:
                alt_nucs = [n for n in nucs if n != target[pos]]
                target[pos] = rng.choice(alt_nucs)
            labels.append(0)
            ca = rng.beta(2, 5, seq_len).astype(np.float32)
            meth = rng.beta(5, 2, seq_len).astype(np.float32)

        target = ''.join(target)
        grna_seqs.append(grna)
        target_seqs.append(target)
        epi_data['chromatin_accessibility'].append(ca)
        epi_data['methylation'].append(meth)
        epi_data['h3k4me3'].append(rng.beta(3, 3, seq_len).astype(np.float32))
        epi_data['h3k27ac'].append(rng.beta(3, 3, seq_len).astype(np.float32))

    return {
        'grna_seqs': grna_seqs,
        'target_seqs': target_seqs,
        'labels': np.array(labels, dtype=np.float32),
        'chromatin_accessibility': np.array(epi_data['chromatin_accessibility']),
        'methylation': np.array(epi_data['methylation']),
        'h3k4me3': np.array(epi_data['h3k4me3']),
        'h3k27ac': np.array(epi_data['h3k27ac'])
    }


def preprocess_dataset(data: Dict, max_len: int = 23) -> Tuple[np.ndarray, np.ndarray]:
    """
    Full preprocessing pipeline: encode sequences + mismatches + epigenetics.
    Returns: (X, y) where X has shape (n_samples, max_len, n_channels)
    """
    n = len(data['grna_seqs'])
    all_features = []

    for i in range(n):
        grna = data['grna_seqs'][i]
        target = data['target_seqs'][i]

        # Sequence encoding (4 channels each = 8)
        grna_enc = one_hot_encode_sequence(grna, max_len)
        target_enc = one_hot_encode_sequence(target, max_len)

        # Mismatch encoding (16 channels)
        mm_enc = encode_mismatch_matrix(grna, target, max_len)

        # Positional features (3 channels)
        pos_enc = encode_position_features(grna, target, max_len)

        # Epigenetic features (4 channels)
        epi_enc = encode_epigenetic_features(
            data['chromatin_accessibility'][i],
            data['methylation'][i],
            data['h3k4me3'][i],
            data['h3k27ac'][i],
            max_len
        )

        # Concatenate: (max_len, 4+4+16+3+4 = 31 channels)
        features = np.concatenate([grna_enc, target_enc, mm_enc, pos_enc, epi_enc], axis=-1)
        all_features.append(features)

    X = np.array(all_features, dtype=np.float32)
    y = data['labels']

    return X, y


def create_cv_splits(
    X: np.ndarray,
    y: np.ndarray,
    n_folds: int = 5,
    seed: int = 42
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Create stratified k-fold cross-validation splits."""
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    splits = [(train_idx, val_idx) for train_idx, val_idx in skf.split(X, y)]
    return splits


if __name__ == '__main__':
    print("Generating synthetic dataset...")
    data = generate_synthetic_dataset(n_samples=5000)
    print(f"  Samples: {len(data['labels'])}")
    print(f"  Positive: {int(data['labels'].sum())}")
    print(f"  Negative: {int(len(data['labels']) - data['labels'].sum())}")

    print("\nPreprocessing...")
    X, y = preprocess_dataset(data)
    print(f"  X shape: {X.shape}")
    print(f"  y shape: {y.shape}")
    print(f"  Feature channels: {X.shape[2]}")

    print("\nCreating CV splits...")
    splits = create_cv_splits(X, y)
    for i, (train_idx, val_idx) in enumerate(splits):
        print(f"  Fold {i+1}: train={len(train_idx)}, val={len(val_idx)}")
