"""
PD-L1 Target Antibody Case Study
In silico design and validation of anti-PD-L1 antibody CDR-H3 sequences.
"""

import torch
import torch.nn.functional as F
import numpy as np
import json
from typing import List, Dict

from antibody_model import (
    AntibodyDesignModel, VOCAB_SIZE, PAD_IDX,
    encode_sequence, decode_sequence, CDR_H3_MAX_LEN
)
from humanization import HumanizationScorePredictor, ImmunogenicityPredictor, assess_humanization_immunogenicity
from developability import (
    ExpressionYieldPredictor, AggregationPredictor, PolyreactivityPredictor,
    assess_developability, compute_developability_index
)
from optimization import run_multi_attribute_optimization, OptimizationWeights


# ─────────────────────────────────────────
# PD-L1 Context (representative sequence segment)
# PD-L1 (CD274) extracellular domain residues 19–239
# Using simplified token representation of key binding epitope
# ─────────────────────────────────────────
PDL1_EPITOPE = (
    "FTIVNPEDSSQIVILNGSQHSLTFQNLTVNRQGLSTATEIEAFEKETGFLLNKVSDGFYPEPVTVSWNSGALTSGVHTFPAVLQSSG"
)

# Known benchmark anti-PD-L1 CDR-H3 sequences
# Based on published structures (simplified/anonymized for in silico demo)
BENCHMARK_CDRS = {
    "atezolizumab_CDR-H3": "SSYSGFFDYWGQGT",   # approximate
    "durvalumab_CDR-H3":   "GGYYDFWSGPFDH",     # approximate
    "avelumab_CDR-H3":     "SGYYVDHYGMDV",      # approximate
    "reference_binder":    "ARSGYDGFAMDY",
}

# Framework region (VH3 germline-based)
FRAMEWORK_VH3 = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDT"
    "AVYYCAK"
)


def encode_protein_segment(seq: str, max_len: int = 120) -> torch.Tensor:
    """Encode a protein sequence to token tensor, truncated to max_len."""
    tok = encode_sequence(seq[:max_len])
    padded = F.pad(tok, (0, max(0, max_len - len(tok))), value=PAD_IDX)
    return padded.unsqueeze(0)


def run_pdl1_case_study(
    model: AntibodyDesignModel,
    human_model: HumanizationScorePredictor,
    immuno_model: ImmunogenicityPredictor,
    expr_model: ExpressionYieldPredictor,
    agg_model: AggregationPredictor,
    psr_model: PolyreactivityPredictor,
    device: str = "cpu",
    n_generated: int = 50,
    n_generations: int = 30,
) -> Dict:
    """
    Full PD-L1 antibody design case study.
    Returns comprehensive results dict.
    """
    print("\n[PD-L1 Case Study] Encoding antigen and framework regions...")

    ag_tokens = encode_protein_segment(PDL1_EPITOPE, max_len=80).to(device)
    fw_tokens = encode_protein_segment(FRAMEWORK_VH3, max_len=100).to(device)

    # ── Step 1: Generate novel CDR-H3 sequences via diffusion ──
    print("[PD-L1 Case Study] Generating novel CDR-H3 sequences via diffusion model...")
    model.eval()
    with torch.no_grad():
        ag_enc = model.encoder(ag_tokens)
        fw_enc = model.encoder(fw_tokens)
        ag_enc_exp = ag_enc.expand(n_generated, -1, -1)
        fw_enc_exp = fw_enc.expand(n_generated, -1, -1)
        generated_tokens = model.diffusion.sample(
            ag_enc_exp, fw_enc_exp,
            cdr_length=12,
            n_samples=n_generated,
            temperature=0.8,
            device=str(device),
        )

    generated_seqs = [decode_sequence(t) for t in generated_tokens]
    print(f"  Generated {len(generated_seqs)} novel CDR-H3 sequences")

    # ── Step 2: Benchmark CDR scoring ──
    print("[PD-L1 Case Study] Scoring benchmark antibodies...")
    all_seqs = list(BENCHMARK_CDRS.values()) + generated_seqs[:20]
    all_labels = list(BENCHMARK_CDRS.keys()) + [f"gen_{i+1:03d}" for i in range(20)]

    # ── Step 3: Humanization & immunogenicity ──
    print("[PD-L1 Case Study] Assessing humanization and immunogenicity...")
    human_results = assess_humanization_immunogenicity(
        all_seqs, human_model, immuno_model, device
    )

    # ── Step 4: Developability ──
    print("[PD-L1 Case Study] Assessing developability...")
    dev_results = assess_developability(all_seqs, expr_model, agg_model, psr_model, device)

    # ── Step 5: In silico affinity & stability ──
    print("[PD-L1 Case Study] Predicting affinity and stability...")
    affinity_results = []
    for seq in all_seqs:
        tok = encode_sequence(seq)
        tok_p = F.pad(tok, (0, max(0, CDR_H3_MAX_LEN - len(tok))), value=PAD_IDX)
        tok_b = tok_p.unsqueeze(0).to(device)
        with torch.no_grad():
            props = model.predict_properties(tok_b, ag_tokens.expand(1, -1))
        affinity_results.append({
            "log_kd": round(props["log_kd"].item(), 4),
            "delta_delta_G": round(props["delta_delta_G"].item(), 4),
            "tm": round(props["Tm"].item() * 90.0, 2),  # denormalize
        })

    # ── Step 6: Multi-attribute optimization ──
    print("[PD-L1 Case Study] Running multi-attribute optimization...")
    seed_seqs = list(BENCHMARK_CDRS.values()) + generated_seqs[:5]
    opt_result = run_multi_attribute_optimization(
        model, ag_tokens.squeeze(0).unsqueeze(0),
        fw_tokens.squeeze(0).unsqueeze(0),
        seed_seqs,
        n_generations=n_generations,
        pop_size=30,
        device=str(device),
    )

    # ── Step 7: Compile combined results ──
    print("[PD-L1 Case Study] Compiling results...")
    combined = []
    for i, (seq, label) in enumerate(zip(all_seqs, all_labels)):
        hr = human_results[i]
        dr = dev_results[i]
        ar = affinity_results[i]
        is_benchmark = label in BENCHMARK_CDRS

        dev_idx = compute_developability_index(
            expression=dr["expression_yield"],
            aggregation=dr["aggregation_score_dl"],
            polyreactivity=dr["polyreactivity_psr"],
            humanization=hr["humanization_score_dl"],
            immunogenicity=hr["immunogenicity_risk"],
        )

        combined.append({
            "label": label,
            "sequence": seq,
            "length": len(seq),
            "is_benchmark": is_benchmark,
            "log_kd": ar["log_kd"],
            "tm": ar["tm"],
            "delta_delta_G": ar["delta_delta_G"],
            "humanization_score": hr["humanization_score_dl"],
            "germline_similarity": hr["germline_similarity"],
            "immunogenicity_risk": hr["immunogenicity_risk"],
            "t_cell_epitope_score": hr["t_cell_epitope_score"],
            "expression_yield": dr["expression_yield"],
            "aggregation_score": dr["aggregation_score_dl"],
            "polyreactivity_psr": dr["polyreactivity_psr"],
            "developability_index": dev_idx,
        })

    # Sort by developability + predicted affinity composite
    combined.sort(
        key=lambda x: x["developability_index"] + (1 / (1 + np.exp(x["log_kd"]))),
        reverse=True,
    )

    return {
        "all_candidates": combined,
        "optimization_result": opt_result,
        "n_generated": n_generated,
        "n_benchmark": len(BENCHMARK_CDRS),
        "pareto_size": len(opt_result["pareto_front"]),
        "top_candidates": combined[:10],
    }


def compute_summary_statistics(candidates: List[Dict]) -> Dict:
    """Compute summary statistics across candidate pool."""
    generated = [c for c in candidates if not c["is_benchmark"]]
    benchmarks = [c for c in candidates if c["is_benchmark"]]

    def stats(vals):
        a = np.array(vals)
        return {"mean": round(float(np.mean(a)), 4), "std": round(float(np.std(a)), 4),
                "min": round(float(np.min(a)), 4), "max": round(float(np.max(a)), 4)}

    return {
        "generated": {
            "n": len(generated),
            "log_kd": stats([c["log_kd"] for c in generated]),
            "tm": stats([c["tm"] for c in generated]),
            "humanization": stats([c["humanization_score"] for c in generated]),
            "immunogenicity": stats([c["immunogenicity_risk"] for c in generated]),
            "developability": stats([c["developability_index"] for c in generated]),
            "aggregation": stats([c["aggregation_score"] for c in generated]),
        },
        "benchmark": {
            "n": len(benchmarks),
            "log_kd": stats([c["log_kd"] for c in benchmarks]) if benchmarks else {},
            "developability": stats([c["developability_index"] for c in benchmarks]) if benchmarks else {},
        },
    }


if __name__ == "__main__":
    print("=== PD-L1 Case Study Module Test ===")
    device = "cpu"
    model = AntibodyDesignModel(d_model=128, T=50)
    human_m = HumanizationScorePredictor(d_model=128)
    immuno_m = ImmunogenicityPredictor(d_model=128, n_hla_alleles=8)
    expr_m = ExpressionYieldPredictor(d_model=128)
    agg_m = AggregationPredictor(d_model=128)
    psr_m = PolyreactivityPredictor(d_model=128)

    result = run_pdl1_case_study(
        model, human_m, immuno_m, expr_m, agg_m, psr_m,
        device=device, n_generated=10, n_generations=5
    )
    print(f"\nTop candidate: {result['top_candidates'][0]['label']}")
    print(f"  Sequence:         {result['top_candidates'][0]['sequence']}")
    print(f"  Dev. index:       {result['top_candidates'][0]['developability_index']:.3f}")
