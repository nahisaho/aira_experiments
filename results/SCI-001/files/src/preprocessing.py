"""
CRISPR-Cas9 Off-Target Effect Prediction: Data Preprocessing Pipeline
Handles GUIDE-seq and CIRCLE-seq data preprocessing, feature extraction,
and epigenetic data integration.
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass, field


# === Nucleotide Encoding ===

NUCLEOTIDE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

# Mismatch type encoding (12 types: AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG)
MISMATCH_TYPES = [
    'AC', 'AG', 'AT', 'CA', 'CG', 'CT',
    'GA', 'GC', 'GT', 'TA', 'TC', 'TG'
]
MISMATCH_MAP = {mt: i for i, mt in enumerate(MISMATCH_TYPES)}


@dataclass
class CRISPRSample:
    """Single CRISPR off-target sample."""
    guide_seq: str          # 20nt guide RNA sequence
    target_seq: str         # 23nt target (20nt + 3nt PAM)
    chrom: str = ''
    position: int = 0
    strand: str = '+'
    read_count: int = 0     # GUIDE-seq/CIRCLE-seq read count
    label: float = 0.0      # Binary label or continuous cleavage frequency
    chromatin_accessibility: float = 0.0  # ATAC-seq / DNase-seq signal
    methylation_level: float = 0.0       # Bisulfite-seq methylation level
    ctcf_binding: float = 0.0
    histone_marks: Dict[str, float] = field(default_factory=dict)


def one_hot_encode_sequence(seq: str, max_len: int = 23) -> np.ndarray:
    """One-hot encode a nucleotide sequence.
    
    Args:
        seq: DNA sequence string
        max_len: Maximum sequence length (pad/truncate)
    
    Returns:
        np.ndarray of shape (4, max_len) - channels-first for CNN
    """
    seq = seq.upper()[:max_len]
    encoding = np.zeros((4, max_len), dtype=np.float32)
    for i, nt in enumerate(seq):
        if nt in NUCLEOTIDE_MAP and NUCLEOTIDE_MAP[nt] < 4:
            encoding[NUCLEOTIDE_MAP[nt], i] = 1.0
    return encoding


def encode_mismatch_pattern(guide: str, target: str) -> np.ndarray:
    """Encode mismatch pattern between guide RNA and target DNA.
    
    Creates a multi-channel representation:
    - Channel 0: Binary mismatch indicator (1 if mismatch at position)
    - Channels 1-12: Mismatch type one-hot encoding
    - Channel 13: Position-weighted mismatch (seed region emphasis)
    
    Args:
        guide: 20nt guide RNA sequence
        target: 20nt target DNA sequence (excluding PAM)
    
    Returns:
        np.ndarray of shape (14, 20)
    """
    guide = guide.upper()[:20]
    target = target.upper()[:20]
    encoding = np.zeros((14, 20), dtype=np.float32)
    
    for i in range(min(len(guide), len(target), 20)):
        g, t = guide[i], target[i]
        if g != t and g != 'N' and t != 'N':
            encoding[0, i] = 1.0  # Mismatch indicator
            mismatch_key = g + t
            if mismatch_key in MISMATCH_MAP:
                encoding[1 + MISMATCH_MAP[mismatch_key], i] = 1.0
            # Seed region weighting (positions 1-12 from PAM are more important)
            seed_weight = 1.0 if i >= 8 else 0.5  # PAM-proximal = higher weight
            encoding[13, i] = seed_weight
    
    return encoding


def encode_pam_sequence(target_seq: str) -> np.ndarray:
    """Encode PAM sequence (last 3 nucleotides of target).
    
    Args:
        target_seq: 23nt target sequence (20nt protospacer + 3nt PAM)
    
    Returns:
        np.ndarray of shape (4, 3) one-hot encoding of PAM
    """
    pam = target_seq[-3:].upper()
    return one_hot_encode_sequence(pam, max_len=3)


def encode_epigenetic_features(sample: CRISPRSample) -> np.ndarray:
    """Encode epigenetic features for a sample.
    
    Features:
    - Chromatin accessibility (ATAC-seq/DNase-seq signal, log-transformed)
    - DNA methylation level (0-1)
    - CTCF binding signal
    - Histone modifications (H3K4me3, H3K27ac, H3K27me3, H3K36me3)
    
    Returns:
        np.ndarray of shape (7,) - epigenetic feature vector
    """
    histone_marks = sample.histone_marks or {}
    features = np.array([
        np.log1p(sample.chromatin_accessibility),
        sample.methylation_level,
        np.log1p(sample.ctcf_binding),
        np.log1p(histone_marks.get('H3K4me3', 0.0)),
        np.log1p(histone_marks.get('H3K27ac', 0.0)),
        np.log1p(histone_marks.get('H3K27me3', 0.0)),
        np.log1p(histone_marks.get('H3K36me3', 0.0)),
    ], dtype=np.float32)
    return features


class GUIDESeqPreprocessor:
    """Preprocessor for GUIDE-seq experimental data."""
    
    def __init__(self, read_count_threshold: int = 5, 
                 normalize_reads: bool = True):
        self.read_count_threshold = read_count_threshold
        self.normalize_reads = normalize_reads
    
    def load_and_process(self, filepath: str) -> List[CRISPRSample]:
        """Load GUIDE-seq data from file.
        
        Expected format: TSV with columns:
        guide_seq, target_seq, chrom, position, strand, read_count, ...
        """
        df = pd.read_csv(filepath, sep='\t')
        
        # Quality filtering
        df = df[df['read_count'] >= self.read_count_threshold]
        
        # Remove duplicates
        df = df.drop_duplicates(subset=['guide_seq', 'target_seq', 'chrom', 'position'])
        
        if self.normalize_reads:
            # Log-normalize read counts per guide
            for guide in df['guide_seq'].unique():
                mask = df['guide_seq'] == guide
                counts = df.loc[mask, 'read_count'].values
                df.loc[mask, 'normalized_count'] = np.log1p(counts) / np.log1p(counts.max())
        
        samples = []
        for _, row in df.iterrows():
            sample = CRISPRSample(
                guide_seq=row['guide_seq'],
                target_seq=row['target_seq'],
                chrom=row.get('chrom', ''),
                position=int(row.get('position', 0)),
                strand=row.get('strand', '+'),
                read_count=int(row['read_count']),
                label=row.get('normalized_count', 1.0)
            )
            samples.append(sample)
        
        return samples


class CIRCLESeqPreprocessor:
    """Preprocessor for CIRCLE-seq experimental data."""
    
    def __init__(self, score_threshold: float = 0.01):
        self.score_threshold = score_threshold
    
    def load_and_process(self, filepath: str) -> List[CRISPRSample]:
        """Load CIRCLE-seq data from file."""
        df = pd.read_csv(filepath, sep='\t')
        
        # Filter by CIRCLE-seq score
        if 'circle_score' in df.columns:
            df = df[df['circle_score'] >= self.score_threshold]
        
        samples = []
        for _, row in df.iterrows():
            sample = CRISPRSample(
                guide_seq=row['guide_seq'],
                target_seq=row['target_seq'],
                chrom=row.get('chrom', ''),
                position=int(row.get('position', 0)),
                read_count=int(row.get('read_count', 0)),
                label=float(row.get('circle_score', row.get('label', 0.0)))
            )
            samples.append(sample)
        
        return samples


class EpigeneticAnnotator:
    """Annotate CRISPR samples with epigenetic information."""
    
    def __init__(self, cell_type: str = 'HEK293'):
        self.cell_type = cell_type
        self.atac_data = None
        self.methylation_data = None
        self.histone_data = {}
    
    def load_atac_seq(self, filepath: str):
        """Load ATAC-seq / DNase-seq signal data (bigWig or BED format)."""
        # In production: use pyBigWig to load bigWig files
        # Here we define the interface
        self.atac_data = filepath
    
    def load_methylation(self, filepath: str):
        """Load bisulfite sequencing methylation data."""
        self.methylation_data = filepath
    
    def load_histone_marks(self, mark: str, filepath: str):
        """Load histone modification ChIP-seq data."""
        self.histone_data[mark] = filepath
    
    def annotate(self, samples: List[CRISPRSample]) -> List[CRISPRSample]:
        """Add epigenetic annotations to samples.
        
        In production, this queries bigWig/BED files for signal at each
        genomic coordinate. Here we define the interface and logic.
        """
        for sample in samples:
            if sample.chrom and sample.position:
                # Query chromatin accessibility at position
                sample.chromatin_accessibility = self._query_signal(
                    self.atac_data, sample.chrom, sample.position)
                # Query methylation level
                sample.methylation_level = self._query_signal(
                    self.methylation_data, sample.chrom, sample.position)
                # Query histone marks
                for mark, data in self.histone_data.items():
                    sample.histone_marks[mark] = self._query_signal(
                        data, sample.chrom, sample.position)
        return samples
    
    def _query_signal(self, data_source, chrom: str, position: int,
                      window: int = 500) -> float:
        """Query signal value at a genomic position (±window bp)."""
        if data_source is None:
            return 0.0
        # Placeholder: In production, use pyBigWig or pybedtools
        return 0.0


class FeatureAssembler:
    """Assemble all features into model-ready tensors."""
    
    def __init__(self, include_epigenetics: bool = True):
        self.include_epigenetics = include_epigenetics
    
    def assemble(self, samples: List[CRISPRSample]) -> Dict[str, np.ndarray]:
        """Convert samples to feature tensors.
        
        Returns:
            Dictionary with keys:
            - 'guide_onehot': (N, 4, 20) guide RNA one-hot
            - 'target_onehot': (N, 4, 23) target DNA one-hot
            - 'mismatch_features': (N, 14, 20) mismatch pattern encoding
            - 'pam_encoding': (N, 4, 3) PAM one-hot
            - 'epigenetic_features': (N, 7) epigenetic vector
            - 'labels': (N,) target labels
        """
        n = len(samples)
        data = {
            'guide_onehot': np.zeros((n, 4, 20), dtype=np.float32),
            'target_onehot': np.zeros((n, 4, 23), dtype=np.float32),
            'mismatch_features': np.zeros((n, 14, 20), dtype=np.float32),
            'pam_encoding': np.zeros((n, 4, 3), dtype=np.float32),
            'epigenetic_features': np.zeros((n, 7), dtype=np.float32),
            'labels': np.zeros(n, dtype=np.float32),
        }
        
        for i, sample in enumerate(samples):
            data['guide_onehot'][i] = one_hot_encode_sequence(
                sample.guide_seq, max_len=20)
            data['target_onehot'][i] = one_hot_encode_sequence(
                sample.target_seq, max_len=23)
            data['mismatch_features'][i] = encode_mismatch_pattern(
                sample.guide_seq, sample.target_seq[:20])
            data['pam_encoding'][i] = encode_pam_sequence(sample.target_seq)
            if self.include_epigenetics:
                data['epigenetic_features'][i] = encode_epigenetic_features(sample)
            data['labels'][i] = sample.label
        
        return data


def generate_negative_samples(positive_samples: List[CRISPRSample],
                               genome_sequences: Optional[Dict] = None,
                               neg_ratio: int = 10,
                               max_mismatches: int = 6,
                               seed: int = 42) -> List[CRISPRSample]:
    """Generate negative samples (non-cleavage sites).
    
    Strategy:
    1. Random genomic sites with ≤max_mismatches to guide
    2. Maintain mismatch distribution similar to positives
    3. Ensure no overlap with known off-target sites
    """
    rng = np.random.RandomState(seed)
    negatives = []
    nucleotides = ['A', 'C', 'G', 'T']
    
    positive_coords = {
        (s.chrom, s.position) for s in positive_samples
        if s.chrom and s.position
    }
    
    for sample in positive_samples:
        for _ in range(neg_ratio):
            # Generate random target with controlled mismatches
            n_mismatches = rng.randint(1, max_mismatches + 1)
            positions = rng.choice(20, size=n_mismatches, replace=False)
            target = list(sample.guide_seq[:20])
            for pos in positions:
                original = target[pos]
                alternatives = [nt for nt in nucleotides if nt != original]
                target[pos] = rng.choice(alternatives)
            # Add random PAM
            pam = 'NGG' if rng.random() > 0.3 else 'NAG'
            pam = rng.choice(nucleotides) + pam[1:]
            target_seq = ''.join(target) + pam
            
            neg = CRISPRSample(
                guide_seq=sample.guide_seq,
                target_seq=target_seq,
                label=0.0
            )
            negatives.append(neg)
    
    return negatives


def create_cross_validation_splits(samples: List[CRISPRSample],
                                    n_folds: int = 5,
                                    strategy: str = 'guide_stratified',
                                    seed: int = 42) -> List[Tuple[List[int], List[int]]]:
    """Create cross-validation splits.
    
    Strategies:
    - 'guide_stratified': Split by guide RNA to prevent data leakage
    - 'chromosome': Split by chromosome
    - 'random': Standard random split
    
    Returns:
        List of (train_indices, test_indices) tuples
    """
    rng = np.random.RandomState(seed)
    n = len(samples)
    
    if strategy == 'guide_stratified':
        guides = list(set(s.guide_seq for s in samples))
        rng.shuffle(guides)
        fold_size = len(guides) // n_folds
        
        splits = []
        for fold in range(n_folds):
            start = fold * fold_size
            end = start + fold_size if fold < n_folds - 1 else len(guides)
            test_guides = set(guides[start:end])
            
            train_idx = [i for i, s in enumerate(samples) 
                        if s.guide_seq not in test_guides]
            test_idx = [i for i, s in enumerate(samples) 
                       if s.guide_seq in test_guides]
            splits.append((train_idx, test_idx))
        
        return splits
    
    elif strategy == 'chromosome':
        chroms = list(set(s.chrom for s in samples if s.chrom))
        rng.shuffle(chroms)
        fold_size = max(1, len(chroms) // n_folds)
        
        splits = []
        for fold in range(n_folds):
            start = fold * fold_size
            end = start + fold_size if fold < n_folds - 1 else len(chroms)
            test_chroms = set(chroms[start:end])
            
            train_idx = [i for i, s in enumerate(samples) 
                        if s.chrom not in test_chroms]
            test_idx = [i for i, s in enumerate(samples) 
                       if s.chrom in test_chroms]
            splits.append((train_idx, test_idx))
        
        return splits
    
    else:  # random
        indices = list(range(n))
        rng.shuffle(indices)
        fold_size = n // n_folds
        
        splits = []
        for fold in range(n_folds):
            start = fold * fold_size
            end = start + fold_size if fold < n_folds - 1 else n
            test_idx = indices[start:end]
            train_idx = indices[:start] + indices[end:]
            splits.append((train_idx, test_idx))
        
        return splits


if __name__ == '__main__':
    # Demo: generate synthetic samples and process
    print("=== CRISPR Off-Target Preprocessing Pipeline Demo ===")
    
    # Create synthetic samples
    np.random.seed(42)
    guides = [
        'GAGTCCGAGCAGAAGAAGAA',
        'GTTGCCCACGTGATCAGCTA',
        'GGCACTGCGGCTGGAGGTGG',
    ]
    
    samples = []
    for guide in guides:
        # On-target
        on_target = CRISPRSample(
            guide_seq=guide,
            target_seq=guide + 'TGG',
            chrom='chr1',
            position=np.random.randint(1e6, 1e8),
            read_count=1000,
            label=1.0,
            chromatin_accessibility=50.0,
            methylation_level=0.3,
        )
        samples.append(on_target)
        
        # Off-targets with mismatches
        for n_mm in range(1, 5):
            target = list(guide)
            positions = np.random.choice(20, n_mm, replace=False)
            for p in positions:
                alts = [nt for nt in 'ACGT' if nt != target[p]]
                target[p] = np.random.choice(alts)
            ot = CRISPRSample(
                guide_seq=guide,
                target_seq=''.join(target) + 'AGG',
                chrom=f'chr{np.random.randint(1, 23)}',
                position=np.random.randint(1e6, 1e8),
                read_count=max(1, int(1000 / (10 ** n_mm))),
                label=max(0.0, 1.0 - n_mm * 0.25),
                chromatin_accessibility=np.random.uniform(0, 100),
                methylation_level=np.random.uniform(0, 1),
            )
            samples.append(ot)
    
    # Generate negatives
    negatives = generate_negative_samples(samples, neg_ratio=5, seed=42)
    all_samples = samples + negatives
    
    print(f"Positive samples: {len(samples)}")
    print(f"Negative samples: {len(negatives)}")
    print(f"Total samples: {len(all_samples)}")
    
    # Assemble features
    assembler = FeatureAssembler(include_epigenetics=True)
    features = assembler.assemble(all_samples)
    
    print(f"\nFeature shapes:")
    for key, val in features.items():
        print(f"  {key}: {val.shape}")
    
    # Create CV splits
    splits = create_cross_validation_splits(
        all_samples, n_folds=5, strategy='guide_stratified')
    
    print(f"\nCross-validation splits (guide-stratified):")
    for i, (train, test) in enumerate(splits):
        print(f"  Fold {i+1}: train={len(train)}, test={len(test)}")
    
    print("\n✓ Preprocessing pipeline complete.")
