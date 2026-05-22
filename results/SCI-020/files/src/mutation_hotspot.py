"""
Mutation Hotspot Prediction and Functional Impact Assessment.
Uses Shannon entropy for hotspot identification and BLOSUM62 for impact scoring.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings("ignore")


# BLOSUM62 substitution matrix (subset of common amino acids)
BLOSUM62 = {
    ('A','A'):4,('A','R'):-1,('A','N'):-2,('A','D'):-2,('A','C'):0,
    ('A','Q'):-1,('A','E'):-1,('A','G'):0,('A','H'):-2,('A','I'):-1,
    ('A','L'):-1,('A','K'):-1,('A','M'):-1,('A','F'):-2,('A','P'):-1,
    ('A','S'):1,('A','T'):0,('A','W'):-3,('A','Y'):-2,('A','V'):0,
    ('R','R'):5,('R','N'):0,('R','D'):-2,('R','C'):-3,('R','Q'):1,
    ('R','E'):0,('R','G'):-2,('R','H'):0,('R','I'):-3,('R','L'):-2,
    ('R','K'):2,('R','M'):-1,('R','F'):-3,('R','P'):-2,('R','S'):-1,
    ('R','T'):-1,('R','W'):-3,('R','Y'):-2,('R','V'):-3,
    ('N','N'):6,('N','D'):1,('N','C'):-3,('N','Q'):0,('N','E'):0,
    ('N','G'):0,('N','H'):1,('N','I'):-3,('N','L'):-3,('N','K'):0,
    ('N','M'):-2,('N','F'):-3,('N','P'):-2,('N','S'):1,('N','T'):0,
    ('N','W'):-4,('N','Y'):-2,('N','V'):-3,
    ('D','D'):6,('D','C'):-3,('D','Q'):0,('D','E'):2,('D','G'):-1,
    ('D','H'):-1,('D','I'):-3,('D','L'):-4,('D','K'):-1,('D','M'):-3,
    ('D','F'):-3,('D','P'):-1,('D','S'):0,('D','T'):-1,('D','W'):-4,
    ('D','Y'):-3,('D','V'):-3,
    ('C','C'):9,('C','Q'):-3,('C','E'):-4,('C','G'):-3,('C','H'):-3,
    ('C','I'):-1,('C','L'):-1,('C','K'):-3,('C','M'):-1,('C','F'):-2,
    ('C','P'):-3,('C','S'):-1,('C','T'):-1,('C','W'):-2,('C','Y'):-2,('C','V'):-1,
    ('Q','Q'):5,('Q','E'):2,('Q','G'):-2,('Q','H'):0,('Q','I'):-3,
    ('Q','L'):-2,('Q','K'):1,('Q','M'):0,('Q','F'):-3,('Q','P'):-1,
    ('Q','S'):0,('Q','T'):-1,('Q','W'):-2,('Q','Y'):-1,('Q','V'):-2,
    ('E','E'):5,('E','G'):-2,('E','H'):0,('E','I'):-3,('E','L'):-3,
    ('E','K'):1,('E','M'):-2,('E','F'):-3,('E','P'):-1,('E','S'):0,
    ('E','T'):-1,('E','W'):-3,('E','Y'):-2,('E','V'):-2,
    ('G','G'):6,('G','H'):-2,('G','I'):-4,('G','L'):-4,('G','K'):-2,
    ('G','M'):-3,('G','F'):-3,('G','P'):-2,('G','S'):0,('G','T'):-2,
    ('G','W'):-2,('G','Y'):-3,('G','V'):-3,
    ('H','H'):8,('H','I'):-3,('H','L'):-3,('H','K'):-1,('H','M'):-2,
    ('H','F'):-1,('H','P'):-2,('H','S'):-1,('H','T'):-2,('H','W'):-2,
    ('H','Y'):2,('H','V'):-3,
    ('I','I'):4,('I','L'):2,('I','K'):-1,('I','M'):1,('I','F'):0,
    ('I','P'):-3,('I','S'):-2,('I','T'):-1,('I','W'):-3,('I','Y'):-1,('I','V'):3,
    ('L','L'):4,('L','K'):-2,('L','M'):2,('L','F'):0,('L','P'):-3,
    ('L','S'):-2,('L','T'):-1,('L','W'):-2,('L','Y'):-1,('L','V'):1,
    ('K','K'):5,('K','M'):-1,('K','F'):-3,('K','P'):-1,('K','S'):0,
    ('K','T'):-1,('K','W'):-3,('K','Y'):-2,('K','V'):-2,
    ('M','M'):5,('M','F'):0,('M','P'):-2,('M','S'):-1,('M','T'):-1,
    ('M','W'):-1,('M','Y'):-1,('M','V'):1,
    ('F','F'):6,('F','P'):-4,('F','S'):-2,('F','T'):-2,('F','W'):1,
    ('F','Y'):3,('F','V'):-1,
    ('P','P'):7,('P','S'):-1,('P','T'):-1,('P','W'):-4,('P','Y'):-3,('P','V'):-2,
    ('S','S'):4,('S','T'):1,('S','W'):-3,('S','Y'):-2,('S','V'):-2,
    ('T','T'):5,('T','W'):-2,('T','Y'):-2,('T','V'):0,
    ('W','W'):11,('W','Y'):2,('W','V'):-3,
    ('Y','Y'):7,('Y','V'):-1,
    ('V','V'):4,
}


def blosum62_score(aa1: str, aa2: str) -> int:
    """Return BLOSUM62 score for an amino acid substitution."""
    if aa1 == aa2:
        return BLOSUM62.get((aa1, aa1), 0)
    key = (aa1, aa2) if (aa1, aa2) in BLOSUM62 else (aa2, aa1)
    return BLOSUM62.get(key, -4)


# Known spike protein ACE2-binding and antibody epitope positions (representative)
KNOWN_EPITOPE_POSITIONS = list(range(417, 508)) + list(range(484, 506))
RBD_POSITIONS = list(range(319, 541))


class MutationHotspotPredictor:
    """Identifies mutation hotspots using Shannon entropy."""

    def identify_hotspots(self, sequences: List[str],
                          window_size: int = 10) -> pd.DataFrame:
        """
        Compute per-position entropy across alignment; return hotspot DataFrame.
        """
        if not sequences:
            return pd.DataFrame()

        min_len = min(len(s) for s in sequences)
        sequences = [s[:min_len] for s in sequences]
        n_seqs = len(sequences)

        results = []
        for pos in range(min_len):
            bases = [seq[pos] for seq in sequences]
            counts = {}
            for b in bases:
                counts[b] = counts.get(b, 0) + 1
            entropy = 0.0
            for cnt in counts.values():
                p = cnt / n_seqs
                if p > 0:
                    entropy -= p * np.log2(p)
            results.append({
                "position": pos,
                "entropy": round(entropy, 4),
                "dominant_base": max(counts, key=counts.get),
                "n_variants": len(counts),
            })

        df = pd.DataFrame(results)
        # Sliding-window mean entropy
        df["window_entropy"] = df["entropy"].rolling(window=window_size,
                                                      center=True).mean().fillna(df["entropy"])
        threshold = df["window_entropy"].mean() + df["window_entropy"].std()
        df["is_hotspot"] = df["window_entropy"] > threshold
        return df

    def compute_mutation_frequency(self, df_seqs: pd.DataFrame) -> pd.DataFrame:
        """Count per-lineage mutation frequencies."""
        records = []
        for lineage, grp in df_seqs.groupby("lineage"):
            records.append({
                "lineage": lineage,
                "n_sequences": len(grp),
                "mean_mutations": round(grp["n_mutations"].mean(), 2),
                "std_mutations": round(grp["n_mutations"].std(), 2),
                "max_mutations": int(grp["n_mutations"].max()),
            })
        return pd.DataFrame(records).sort_values("mean_mutations", ascending=False)


class FunctionalImpactAssessor:
    """Assesses functional impact of mutations using BLOSUM62 and epitope overlap."""

    def generate_spike_mutations(self, n: int = 30, seed: int = 42) -> pd.DataFrame:
        """Generate representative spike protein mutations for analysis."""
        rng = np.random.default_rng(seed)
        amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
        positions = rng.integers(1, 1274, size=n)  # Spike protein length
        from_aa = rng.choice(amino_acids, size=n)
        to_aa = rng.choice(amino_acids, size=n)
        # Ensure they differ
        for i in range(n):
            while to_aa[i] == from_aa[i]:
                to_aa[i] = rng.choice(amino_acids)

        records = []
        for i in range(n):
            pos = int(positions[i])
            fa, ta = from_aa[i], to_aa[i]
            score = blosum62_score(fa, ta)
            in_rbd = pos in RBD_POSITIONS
            in_epitope = pos in KNOWN_EPITOPE_POSITIONS
            impact = "benign" if score >= 1 else ("moderate" if score >= -2 else "deleterious")
            records.append({
                "mutation": f"{fa}{pos}{ta}",
                "position": pos,
                "from_aa": fa,
                "to_aa": ta,
                "blosum62_score": score,
                "impact_class": impact,
                "in_rbd": in_rbd,
                "in_epitope": in_epitope,
                "immune_escape_risk": in_epitope and score < 0,
            })
        return pd.DataFrame(records)

    def rank_mutations_by_risk(self, mutations_df: pd.DataFrame) -> pd.DataFrame:
        """Compute composite risk score and rank mutations."""
        df = mutations_df.copy()
        # Risk = -blosum62 (lower is worse) + RBD bonus + epitope bonus
        df["risk_score"] = (
            (-df["blosum62_score"]).clip(lower=0) * 10
            + df["in_rbd"].astype(int) * 20
            + df["in_epitope"].astype(int) * 25
            + df["immune_escape_risk"].astype(int) * 15
        )
        df["risk_score"] = df["risk_score"].clip(0, 100)
        return df.sort_values("risk_score", ascending=False).reset_index(drop=True)

    def assess_immune_escape(self, mutations_df: pd.DataFrame) -> Dict:
        """Summary statistics on immune escape potential."""
        n_rbd = int(mutations_df["in_rbd"].sum())
        n_epitope = int(mutations_df["in_epitope"].sum())
        n_escape = int(mutations_df["immune_escape_risk"].sum())
        high_risk = mutations_df[mutations_df["risk_score"] > 50] if "risk_score" in mutations_df.columns else pd.DataFrame()
        return {
            "n_mutations_analyzed": len(mutations_df),
            "n_rbd_mutations": n_rbd,
            "n_epitope_mutations": n_epitope,
            "n_immune_escape_risk": n_escape,
            "n_high_risk_mutations": len(high_risk),
            "immune_escape_rate": round(n_escape / max(len(mutations_df), 1), 3),
        }


def run_mutation_analysis(sequences_df: Optional[pd.DataFrame] = None) -> Dict:
    """Run the complete mutation hotspot and impact pipeline."""
    predictor = MutationHotspotPredictor()
    assessor = FunctionalImpactAssessor()

    if sequences_df is not None:
        seqs = sequences_df["sequence"].tolist()
        hotspots_df = predictor.identify_hotspots(seqs[:50])  # limit for speed
        mutation_freq = predictor.compute_mutation_frequency(sequences_df)
    else:
        hotspots_df = pd.DataFrame()
        mutation_freq = pd.DataFrame()

    mutations_df = assessor.generate_spike_mutations(n=30)
    mutations_df = assessor.rank_mutations_by_risk(mutations_df)
    escape_stats = assessor.assess_immune_escape(mutations_df)

    return {
        "hotspots_df": hotspots_df,
        "mutations_df": mutations_df,
        "mutation_freq": mutation_freq,
        "escape_stats": escape_stats,
        "n_hotspot_positions": int(hotspots_df["is_hotspot"].sum()) if len(hotspots_df) else 0,
        "top_mutations": mutations_df.head(5)["mutation"].tolist() if len(mutations_df) else [],
    }


if __name__ == "__main__":
    results = run_mutation_analysis()
    print(f"Hotspot positions: {results['n_hotspot_positions']}")
    print(f"Top high-risk mutations: {results['top_mutations']}")
    print(f"Immune escape stats: {results['escape_stats']}")
