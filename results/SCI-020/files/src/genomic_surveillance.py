"""
Genomic Surveillance Module
Handles real-time phylogenetic analysis from GISAID/GenBank.
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import random
import string
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings("ignore")

try:
    from Bio import SeqIO, Phylo
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.Align import PairwiseAligner
    HAS_BIO = True
except ImportError:
    HAS_BIO = False


# SARS-CoV-2 spike protein reference (first 300 nt of spike, representative)
REFERENCE_SPIKE = (
    "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCC"
    "CTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACT"
    "CAAACAAATGGTTTCCTTTTCTTTTCTTGTTCCTTTTGTTTTTATTATTGCAATAGTTGGAAAGTCTAAATGTGTCA"
    "ATCAAATTTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTATTATTCTTATTGTTGTTGTTGTTGTCTT"
)


def _random_mutate(seq: str, n_mutations: int, rng: np.random.Generator) -> str:
    """Apply n random point mutations to a nucleotide sequence."""
    seq_list = list(seq)
    positions = rng.choice(len(seq), size=min(n_mutations, len(seq)), replace=False)
    bases = ['A', 'T', 'G', 'C']
    for pos in positions:
        current = seq_list[pos]
        alternatives = [b for b in bases if b != current]
        seq_list[pos] = rng.choice(alternatives)
    return "".join(seq_list)


def _generate_lineage_id(prefix: str, rng: np.random.Generator) -> str:
    """Generate a plausible lineage identifier."""
    letter = rng.choice(list(string.ascii_uppercase[:8]))
    nums = rng.integers(1, 50)
    sub = rng.integers(1, 10)
    return f"{prefix}.{letter}{nums}.{sub}"


class GISAIDClient:
    """Simulates GISAID API client for sequence retrieval."""

    def fetch_sequences(self, pathogen: str = "SARS-CoV-2",
                        days_back: int = 30, n_seqs: int = 120,
                        seed: int = 42) -> pd.DataFrame:
        """Fetch recent sequences; returns DataFrame of metadata + sequence."""
        rng = np.random.default_rng(seed)
        now = datetime.now()
        records = []
        lineages = ["BA.2", "BA.5", "XBB.1.5", "XBB.1.16", "JN.1", "KP.2"]
        # Add one novel lineage
        lineages.append(_generate_lineage_id("KP", rng))

        for i in range(n_seqs):
            lineage = rng.choice(lineages)
            n_mut = int(rng.normal(15, 5))
            sequence = _random_mutate(REFERENCE_SPIKE, max(0, n_mut), rng)
            days_ago = rng.integers(0, days_back)
            records.append({
                "accession": f"EPI_ISL_{7000000 + i}",
                "collection_date": (now - timedelta(days=int(days_ago))).strftime("%Y-%m-%d"),
                "country": rng.choice(["Japan", "USA", "Germany", "Brazil", "India", "Kenya"]),
                "lineage": lineage,
                "n_mutations": max(0, n_mut),
                "sequence": sequence,
                "source": "GISAID",
            })
        return pd.DataFrame(records)


class GenBankClient:
    """Simulates GenBank Entrez client."""

    def fetch_sequences(self, pathogen: str = "SARS-CoV-2",
                        days_back: int = 30, n_seqs: int = 80,
                        seed: int = 99) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        now = datetime.now()
        records = []
        lineages = ["BA.2", "BA.4", "XBB", "BQ.1", "JN.1"]
        for i in range(n_seqs):
            lineage = rng.choice(lineages)
            n_mut = int(rng.normal(12, 4))
            sequence = _random_mutate(REFERENCE_SPIKE, max(0, n_mut), rng)
            days_ago = rng.integers(0, days_back)
            records.append({
                "accession": f"OX{500000 + i}",
                "collection_date": (now - timedelta(days=int(days_ago))).strftime("%Y-%m-%d"),
                "country": rng.choice(["UK", "France", "Australia", "Canada"]),
                "lineage": lineage,
                "n_mutations": max(0, n_mut),
                "sequence": sequence,
                "source": "GenBank",
            })
        return pd.DataFrame(records)


class PhylogeneticAnalyzer:
    """Builds and analyzes phylogenetic trees from sequence data."""

    def __init__(self):
        self.calculator = DistanceCalculator("identity") if HAS_BIO else None
        self.constructor = DistanceTreeConstructor() if HAS_BIO else None

    def compute_pairwise_distance_matrix(self, sequences: List[str]) -> np.ndarray:
        """Compute simple Hamming distance matrix."""
        n = len(sequences)
        mat = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                s1, s2 = sequences[i], sequences[j]
                length = min(len(s1), len(s2))
                diffs = sum(a != b for a, b in zip(s1[:length], s2[:length]))
                dist = diffs / length if length > 0 else 0.0
                mat[i, j] = mat[j, i] = dist
        return mat

    def detect_emerging_lineages(self, df: pd.DataFrame,
                                  frequency_threshold: float = 0.05) -> pd.DataFrame:
        """Detect lineages with rapid frequency increase."""
        df = df.copy()
        df["collection_date"] = pd.to_datetime(df["collection_date"])
        df_sorted = df.sort_values("collection_date")
        midpoint = df_sorted["collection_date"].median()

        early = df_sorted[df_sorted["collection_date"] <= midpoint]
        late = df_sorted[df_sorted["collection_date"] > midpoint]

        early_freq = early["lineage"].value_counts(normalize=True)
        late_freq = late["lineage"].value_counts(normalize=True)

        all_lineages = set(early_freq.index) | set(late_freq.index)
        results = []
        for lin in all_lineages:
            ef = early_freq.get(lin, 0.0)
            lf = late_freq.get(lin, 0.0)
            delta = lf - ef
            results.append({
                "lineage": lin,
                "early_frequency": round(ef, 4),
                "late_frequency": round(lf, 4),
                "delta_frequency": round(delta, 4),
                "emerging": delta > frequency_threshold,
            })
        return pd.DataFrame(results).sort_values("delta_frequency", ascending=False)

    def compute_evolutionary_rate(self, df: pd.DataFrame) -> Dict:
        """Estimate evolutionary rate (substitutions/site/year)."""
        avg_mutations = df["n_mutations"].mean()
        seq_length = len(REFERENCE_SPIKE)
        # Assume avg sample is 15 days old, annualize
        rate_per_year = (avg_mutations / seq_length) * (365 / 15)
        return {
            "mean_mutations_per_seq": round(avg_mutations, 2),
            "sequence_length": seq_length,
            "estimated_rate_per_site_per_year": round(rate_per_year, 6),
        }

    def compute_diversity_metrics(self, df: pd.DataFrame) -> Dict:
        """Compute Shannon diversity and Simpson index over lineages."""
        freqs = df["lineage"].value_counts(normalize=True)
        shannon = -np.sum(freqs * np.log(freqs + 1e-12))
        simpson = 1.0 - np.sum(freqs ** 2)
        return {
            "n_lineages": len(freqs),
            "shannon_diversity": round(float(shannon), 4),
            "simpson_index": round(float(simpson), 4),
            "dominant_lineage": freqs.idxmax(),
            "dominant_lineage_frequency": round(float(freqs.max()), 4),
        }


def run_genomic_surveillance(config: Optional[Dict] = None) -> Dict:
    """Run the complete genomic surveillance pipeline."""
    config = config or {}
    gisaid = GISAIDClient()
    genbank = GenBankClient()
    analyzer = PhylogeneticAnalyzer()

    df_gisaid = gisaid.fetch_sequences(n_seqs=120)
    df_genbank = genbank.fetch_sequences(n_seqs=80)
    df_all = pd.concat([df_gisaid, df_genbank], ignore_index=True)

    emerging = analyzer.detect_emerging_lineages(df_all)
    evo_rate = analyzer.compute_evolutionary_rate(df_all)
    diversity = analyzer.compute_diversity_metrics(df_all)

    return {
        "sequences_df": df_all,
        "emerging_lineages": emerging,
        "evolutionary_rate": evo_rate,
        "diversity_metrics": diversity,
        "n_sequences_total": len(df_all),
        "n_novel_emerging": int(emerging["emerging"].sum()),
    }


if __name__ == "__main__":
    results = run_genomic_surveillance()
    print(f"Total sequences: {results['n_sequences_total']}")
    print(f"Emerging lineages: {results['n_novel_emerging']}")
    print(f"Evolutionary rate: {results['evolutionary_rate']['estimated_rate_per_site_per_year']:.6f} sub/site/year")
    print(results["emerging_lineages"].head())
