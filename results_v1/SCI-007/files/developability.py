"""
Developability Prediction Module
Predicts expression yield and aggregation propensity for antibody candidates.
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import List, Optional, Tuple
from antibody_model import CDRStructureEncoder, VOCAB_SIZE, PAD_IDX, AMINO_ACIDS, encode_sequence


# ─────────────────────────────────────────
# Biophysical feature calculators
# ─────────────────────────────────────────
AA_PROPERTIES = {
    # aa: (hydrophobicity, charge, size, flexibility)
    "A": (1.8,  0.0, 1.0, 0.36),
    "R": (-4.5, 1.0, 4.0, 0.53),
    "N": (-3.5, 0.0, 2.0, 0.46),
    "D": (-3.5, -1.0, 2.0, 0.51),
    "C": (2.5,  0.0, 1.5, 0.35),
    "Q": (-3.5, 0.0, 3.0, 0.49),
    "E": (-3.5, -1.0, 3.0, 0.50),
    "G": (-0.4, 0.0, 0.5, 0.54),
    "H": (-3.2, 0.5, 3.0, 0.32),
    "I": (4.5,  0.0, 2.0, 0.46),
    "L": (3.8,  0.0, 2.0, 0.59),
    "K": (-3.9, 1.0, 3.5, 0.47),
    "M": (1.9,  0.0, 2.5, 0.60),
    "F": (2.8,  0.0, 3.0, 0.31),
    "P": (-1.6, 0.0, 1.5, 0.51),
    "S": (-0.8, 0.0, 1.0, 0.51),
    "T": (-0.7, 0.0, 1.5, 0.44),
    "W": (-0.9, 0.0, 4.0, 0.31),
    "Y": (-1.3, 0.0, 3.5, 0.42),
    "V": (4.2,  0.0, 1.5, 0.39),
}

# Aggregation-prone motifs
AGG_PRONE_MOTIFS = [
    "VH",  "LY", "YL", "VL", "FW", "WF",
    "II",  "LL", "VV", "IV", "VI",
    "VHVH", "FLLL", "WWWW",
]


def compute_biophysical_features(seq: str) -> np.ndarray:
    """
    Compute per-sequence biophysical features:
    [mean_hydrophobicity, net_charge, mean_size, mean_flexibility,
     fraction_hydrophobic, fraction_charged, GRAVY, instability_index,
     patches_hydrophobic, patches_charged]
    """
    if not seq:
        return np.zeros(10)

    hydro = [AA_PROPERTIES.get(aa, (0, 0, 0, 0))[0] for aa in seq]
    charge = [AA_PROPERTIES.get(aa, (0, 0, 0, 0))[1] for aa in seq]
    size = [AA_PROPERTIES.get(aa, (0, 0, 0, 0))[2] for aa in seq]
    flex = [AA_PROPERTIES.get(aa, (0, 0, 0, 0))[3] for aa in seq]

    mean_hydro = np.mean(hydro)
    net_charge = sum(charge)
    mean_size = np.mean(size)
    mean_flex = np.mean(flex)
    frac_hydro = sum(1 for h in hydro if h > 0) / len(seq)
    frac_charged = sum(1 for c in charge if c != 0) / len(seq)
    gravy = sum(hydro) / len(seq)

    # Simplified instability index (based on DIWV dipeptide weights)
    instability = sum(
        1.0 if seq[i] in "RFWY" else 0.0
        for i in range(len(seq))
    ) / len(seq) * 100

    # Count hydrophobic/charged patches (runs of 3+)
    hydro_patches = sum(
        1 for i in range(len(seq) - 2)
        if all(hydro[i + j] > 1.5 for j in range(3))
    )
    charge_patches = sum(
        1 for i in range(len(seq) - 2)
        if abs(sum(charge[i + j] for j in range(3))) >= 2.5
    )

    return np.array([
        mean_hydro, net_charge, mean_size, mean_flex,
        frac_hydro, frac_charged, gravy, instability,
        hydro_patches, charge_patches,
    ], dtype=np.float32)


def compute_aggregation_score(seq: str) -> float:
    """
    Rule-based aggregation propensity score (0–1).
    Higher = more aggregation prone.
    """
    score = 0.0
    # Check aggregation-prone motifs
    for motif in AGG_PRONE_MOTIFS:
        if motif in seq:
            score += 0.1 * len(motif)

    # High hydrophobicity penalty
    bp = compute_biophysical_features(seq)
    if bp[0] > 2.0:  # mean_hydro
        score += (bp[0] - 2.0) * 0.2
    if bp[8] > 0:    # hydrophobic patches
        score += bp[8] * 0.15

    return min(score, 1.0)


# ─────────────────────────────────────────
# 1. Expression Yield Predictor
# ─────────────────────────────────────────
class ExpressionYieldPredictor(nn.Module):
    """
    Predicts relative expression yield (0–1, normalized titer) from sequence.
    Uses sequence encoder + biophysical features.
    """

    def __init__(self, d_model: int = 256, n_biophys: int = 10, dropout: float = 0.1):
        super().__init__()
        self.encoder = CDRStructureEncoder(d_model=d_model)
        self.biophys_proj = nn.Linear(n_biophys, d_model // 4)
        self.predictor = nn.Sequential(
            nn.Linear(d_model + d_model // 4, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, 64),
            nn.GELU(),
            nn.Linear(64, 1),
            nn.Sigmoid(),
        )

    def forward(
        self,
        seq_tokens: torch.Tensor,       # (B, L)
        biophys_feats: torch.Tensor,    # (B, n_biophys)
    ) -> torch.Tensor:
        enc = self.encoder(seq_tokens)          # (B, L, d)
        seq_feat = enc.mean(dim=1)              # (B, d)
        phys_feat = F.gelu(self.biophys_proj(biophys_feats))   # (B, d/4)
        combined = torch.cat([seq_feat, phys_feat], dim=-1)
        return self.predictor(combined).squeeze(-1)  # (B,)


# ─────────────────────────────────────────
# 2. Aggregation Propensity Predictor
# ─────────────────────────────────────────
class AggregationPredictor(nn.Module):
    """
    Predicts aggregation propensity and colloidal stability (B22 proxy).
    Outputs: [aggregation_score, b22_proxy]
    """

    def __init__(self, d_model: int = 256, n_biophys: int = 10, dropout: float = 0.1):
        super().__init__()
        self.encoder = CDRStructureEncoder(d_model=d_model)
        self.biophys_proj = nn.Linear(n_biophys, d_model // 4)
        self.predictor = nn.Sequential(
            nn.Linear(d_model + d_model // 4, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, 64),
            nn.GELU(),
            nn.Linear(64, 2),  # [aggregation_score, b22_proxy]
        )
        self.agg_sigmoid = nn.Sigmoid()

    def forward(
        self,
        seq_tokens: torch.Tensor,
        biophys_feats: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        enc = self.encoder(seq_tokens)
        seq_feat = enc.mean(dim=1)
        phys_feat = F.gelu(self.biophys_proj(biophys_feats))
        combined = torch.cat([seq_feat, phys_feat], dim=-1)
        out = self.predictor(combined)          # (B, 2)
        agg_score = self.agg_sigmoid(out[:, 0])
        b22 = out[:, 1]                         # not bounded (can be negative)
        return agg_score, b22


# ─────────────────────────────────────────
# 3. Polyreactivity Predictor (PSR score)
# ─────────────────────────────────────────
class PolyreactivityPredictor(nn.Module):
    """
    Predicts polyspecificity reagent (PSR) score — proxy for non-specific binding.
    Low PSR is desired for developability.
    """

    def __init__(self, d_model: int = 256, dropout: float = 0.1):
        super().__init__()
        self.encoder = CDRStructureEncoder(d_model=d_model)
        self.head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(128, 1),
            nn.Sigmoid(),
        )

    def forward(self, seq_tokens: torch.Tensor) -> torch.Tensor:
        enc = self.encoder(seq_tokens)
        return self.head(enc.mean(dim=1)).squeeze(-1)


# ─────────────────────────────────────────
# 4. Developability Index (composite)
# ─────────────────────────────────────────
def compute_developability_index(
    expression: float,
    aggregation: float,
    polyreactivity: float,
    humanization: float,
    immunogenicity: float,
    weights: Optional[dict] = None,
) -> float:
    """
    Composite developability index (0–1).
    Higher = better developability.
    """
    if weights is None:
        weights = {
            "expression": 0.25,
            "aggregation": 0.25,       # inverted
            "polyreactivity": 0.20,    # inverted
            "humanization": 0.20,
            "immunogenicity": 0.10,    # inverted
        }

    score = (
        weights["expression"] * expression
        + weights["aggregation"] * (1.0 - aggregation)
        + weights["polyreactivity"] * (1.0 - polyreactivity)
        + weights["humanization"] * humanization
        + weights["immunogenicity"] * (1.0 - immunogenicity)
    )
    return float(np.clip(score, 0, 1))


# ─────────────────────────────────────────
# 5. Full Developability Assessment
# ─────────────────────────────────────────
def assess_developability(
    sequences: List[str],
    expr_model: ExpressionYieldPredictor,
    agg_model: AggregationPredictor,
    psr_model: PolyreactivityPredictor,
    device: str = "cpu",
) -> List[dict]:
    """Batch developability assessment for a list of antibody sequences."""
    results = []
    for seq in sequences:
        tok = encode_sequence(seq)
        tok_padded = F.pad(tok, (0, max(0, 25 - len(tok))), value=PAD_IDX)
        tok_batch = tok_padded.unsqueeze(0).to(device)

        biophys = compute_biophysical_features(seq)
        biophys_tensor = torch.tensor(biophys, dtype=torch.float32).unsqueeze(0).to(device)

        with torch.no_grad():
            expr = expr_model(tok_batch, biophys_tensor).item()
            agg_score, b22 = agg_model(tok_batch, biophys_tensor)
            agg_val = agg_score.item()
            b22_val = b22.item()
            psr_val = psr_model(tok_batch).item()

        rule_agg = compute_aggregation_score(seq)
        dev_idx = compute_developability_index(expr, agg_val, psr_val, 0.8, 0.2)

        results.append({
            "sequence": seq,
            "expression_yield": round(expr, 4),
            "aggregation_score_dl": round(agg_val, 4),
            "aggregation_score_rule": round(rule_agg, 4),
            "b22_proxy": round(b22_val, 4),
            "polyreactivity_psr": round(psr_val, 4),
            "gravy_score": round(float(biophys[6]), 4),
            "instability_index": round(float(biophys[7]), 4),
            "net_charge": round(float(biophys[1]), 1),
            "frac_hydrophobic": round(float(biophys[4]), 4),
            "developability_index": round(dev_idx, 4),
        })
    return results


if __name__ == "__main__":
    print("=== Developability Module Test ===")
    device = "cpu"
    expr_model = ExpressionYieldPredictor(d_model=128)
    agg_model = AggregationPredictor(d_model=128)
    psr_model = PolyreactivityPredictor(d_model=128)

    test_seqs = [
        "ARSGYDGFDY",
        "YWYCARDLGYYY",
        "GQGTTLTVSS",
        "VVVLLLFFFY",  # high aggregation
    ]

    results = assess_developability(test_seqs, expr_model, agg_model, psr_model, device)
    for r in results:
        print(f"\n  Seq: {r['sequence']}")
        print(f"    Expression:          {r['expression_yield']:.3f}")
        print(f"    Aggregation (DL):    {r['aggregation_score_dl']:.3f}")
        print(f"    Aggregation (rule):  {r['aggregation_score_rule']:.3f}")
        print(f"    Polyreactivity PSR:  {r['polyreactivity_psr']:.3f}")
        print(f"    Developability Idx:  {r['developability_index']:.3f}")
