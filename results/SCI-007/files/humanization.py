"""
Humanization Score Prediction & Immunogenicity Risk Assessment
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import List, Tuple, Optional
from antibody_model import CDRStructureEncoder, VOCAB_SIZE, PAD_IDX, AMINO_ACIDS


# ─────────────────────────────────────────
# Human germline sequence database (simplified)
# In production, use IMGT or OAS database
# ─────────────────────────────────────────
HUMAN_VH_GERMLINES = {
    # Representative IGHV germline CDR-H3 anchor sequences
    "IGHV1-2": "EVQLVESGGGLVQPGGSLRLSCAASGFTFS",
    "IGHV1-69": "EVQLVESGGGLVQPGGSLRLSCAASGGTFS",
    "IGHV3-23": "EVQLVESGGGLVQPGRSLRLSCAASGFTFS",
    "IGHV4-34": "QVQLQESGPGLVKPSQTLSLTCTVSGGSIS",
    "IGHV5-51": "EVQLVQSGAEVKKPGASVKVSCKASGYTFT",
}

# MHC-II binding motifs associated with immunogenicity (simplified)
IMMUNOGENIC_MOTIFS = [
    "FVNQHLCG",  # insulin-like
    "WGQGTLVT",  # common VH framework
    "YYCARSGYD",  # example CDR-H3
]


# ─────────────────────────────────────────
# 1. Humanization Score Predictor
# ─────────────────────────────────────────
class HumanizationScorePredictor(nn.Module):
    """
    Predicts humanization score (0–1) for an antibody sequence.
    Based on:
      - Sequence similarity to human germline (Hu score)
      - Z-score from human antibody repertoire
      - Deep learning component trained on OAS data
    """

    def __init__(self, d_model: int = 256, dropout: float = 0.1):
        super().__init__()
        self.encoder = CDRStructureEncoder(d_model=d_model)
        self.classifier = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(128, 64),
            nn.GELU(),
            nn.Linear(64, 1),
            nn.Sigmoid(),
        )

    def forward(self, seq_tokens: torch.Tensor) -> torch.Tensor:
        """Returns humanization probability (B,)"""
        enc = self.encoder(seq_tokens)          # (B, L, d)
        pooled = enc.mean(dim=1)                # (B, d)
        return self.classifier(pooled).squeeze(-1)

    @staticmethod
    def compute_sequence_similarity(seq1: str, seq2: str) -> float:
        """Pairwise sequence identity (simple)."""
        min_len = min(len(seq1), len(seq2))
        if min_len == 0:
            return 0.0
        matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
        return matches / max(len(seq1), len(seq2))

    @staticmethod
    def germline_similarity_score(cdr_seq: str) -> dict:
        """
        Computes max similarity to known human IGHV germlines.
        Returns dict with best match and score.
        """
        best_score = 0.0
        best_germline = None
        for germline_id, germline_seq in HUMAN_VH_GERMLINES.items():
            score = HumanizationScorePredictor.compute_sequence_similarity(
                cdr_seq, germline_seq
            )
            if score > best_score:
                best_score = score
                best_germline = germline_id
        return {"best_germline": best_germline, "similarity": best_score}

    @staticmethod
    def hu_score(cdr_seq: str, framework_seq: str = "") -> float:
        """
        Simplified Hu score calculation.
        Based on percentage of human-like positions.
        """
        full_seq = framework_seq + cdr_seq
        human_aa_freq = {
            "A": 0.07, "R": 0.06, "N": 0.04, "D": 0.05, "C": 0.02,
            "Q": 0.04, "E": 0.06, "G": 0.07, "H": 0.02, "I": 0.06,
            "L": 0.10, "K": 0.06, "M": 0.02, "F": 0.04, "P": 0.05,
            "S": 0.07, "T": 0.06, "W": 0.01, "Y": 0.03, "V": 0.07,
        }
        score = sum(human_aa_freq.get(aa, 0) for aa in full_seq)
        return min(score / len(full_seq) if full_seq else 0.0, 1.0) * 10.0


# ─────────────────────────────────────────
# 2. Immunogenicity Risk Predictor
# ─────────────────────────────────────────
class ImmunogenicityPredictor(nn.Module):
    """
    Predicts T-cell epitope risk (immunogenicity) based on:
      - MHC-II binding score approximation
      - Sequence-based risk features
      - Deep learning classifier
    """

    def __init__(self, d_model: int = 256, n_hla_alleles: int = 8, dropout: float = 0.1):
        super().__init__()
        self.n_hla_alleles = n_hla_alleles
        self.encoder = CDRStructureEncoder(d_model=d_model)
        # Per-allele binding score head
        self.mhc_head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(128, n_hla_alleles),
            nn.Sigmoid(),
        )
        # Overall immunogenicity risk
        self.risk_head = nn.Sequential(
            nn.Linear(d_model + n_hla_alleles, 64),
            nn.GELU(),
            nn.Linear(64, 1),
            nn.Sigmoid(),
        )

    def forward(
        self, seq_tokens: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Returns:
          mhc_scores: (B, n_hla_alleles) — MHC-II binding probability per allele
          risk_score: (B,) — overall immunogenicity risk
        """
        enc = self.encoder(seq_tokens)          # (B, L, d)
        pooled = enc.mean(dim=1)                # (B, d)
        mhc_scores = self.mhc_head(pooled)      # (B, n_hla)
        feat = torch.cat([pooled, mhc_scores], dim=-1)
        risk_score = self.risk_head(feat).squeeze(-1)
        return mhc_scores, risk_score

    @staticmethod
    def check_known_motifs(seq: str) -> List[str]:
        """Scan for known immunogenic motifs."""
        found = []
        for motif in IMMUNOGENIC_MOTIFS:
            if motif in seq:
                found.append(motif)
        return found

    @staticmethod
    def estimate_t_cell_epitope_score(seq: str, window: int = 9) -> float:
        """
        Simplified T-cell epitope scoring using BLOSUM-based hydrophobicity.
        In production, use NetMHCIIpan.
        """
        hydrophobicity = {
            "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
            "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
            "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
            "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
        }
        if len(seq) < window:
            return 0.0
        scores = []
        for i in range(len(seq) - window + 1):
            kmer = seq[i : i + window]
            score = sum(hydrophobicity.get(aa, 0) for aa in kmer) / window
            scores.append(score)
        return max(scores) if scores else 0.0


# ─────────────────────────────────────────
# 3. Combined Humanization + Immunogenicity Assessment
# ─────────────────────────────────────────
def assess_humanization_immunogenicity(
    sequences: List[str],
    model_human: HumanizationScorePredictor,
    model_immuno: ImmunogenicityPredictor,
    device: str = "cpu",
) -> List[dict]:
    """
    Full assessment pipeline for a list of CDR-H3 sequences.
    Returns list of dicts with all scores.
    """
    from antibody_model import encode_sequence, PAD_IDX

    results = []
    for seq in sequences:
        tok = encode_sequence(seq)
        tok_padded = F.pad(tok, (0, max(0, 25 - len(tok))), value=PAD_IDX)
        tok_batch = tok_padded.unsqueeze(0).to(device)

        with torch.no_grad():
            hu_score_dl = model_human(tok_batch).item()
            mhc_scores, risk = model_immuno(tok_batch)
            mhc_scores = mhc_scores.squeeze(0).cpu().numpy()
            risk_val = risk.item()

        germline = HumanizationScorePredictor.germline_similarity_score(seq)
        hu_rule = HumanizationScorePredictor.hu_score(seq)
        t_cell = ImmunogenicityPredictor.estimate_t_cell_epitope_score(seq)
        motifs = ImmunogenicityPredictor.check_known_motifs(seq)

        results.append({
            "sequence": seq,
            "length": len(seq),
            "humanization_score_dl": round(hu_score_dl, 4),
            "hu_score_rule": round(hu_rule, 4),
            "best_germline": germline["best_germline"],
            "germline_similarity": round(germline["similarity"], 4),
            "immunogenicity_risk": round(risk_val, 4),
            "t_cell_epitope_score": round(t_cell, 4),
            "mhc_binding_mean": round(float(mhc_scores.mean()), 4),
            "known_immunogenic_motifs": motifs,
        })
    return results


if __name__ == "__main__":
    print("=== Humanization & Immunogenicity Module Test ===")
    device = "cpu"
    human_model = HumanizationScorePredictor(d_model=128)
    immuno_model = ImmunogenicityPredictor(d_model=128, n_hla_alleles=8)

    test_seqs = [
        "ARSGYDGFDY",    # CDR-H3 from existing antibody
        "YWYCARDLGYYY",  # synthetic CDR-H3
        "GQGTTLTVSS",    # CDR-L3 like
    ]

    results = assess_humanization_immunogenicity(test_seqs, human_model, immuno_model, device)
    for r in results:
        print(f"\n  Seq: {r['sequence']}")
        print(f"    Humanization (DL):   {r['humanization_score_dl']:.3f}")
        print(f"    Hu-score (rule):     {r['hu_score_rule']:.3f}")
        print(f"    Germline similarity: {r['germline_similarity']:.3f} ({r['best_germline']})")
        print(f"    Immunogenicity risk: {r['immunogenicity_risk']:.3f}")
        print(f"    T-cell epitope:      {r['t_cell_epitope_score']:.3f}")
