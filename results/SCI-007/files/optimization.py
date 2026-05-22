"""
Multi-Attribute Optimization (MAO) for Antibody Design
Optimizes CDR sequences across: affinity, specificity, stability,
humanization, and developability using evolutionary + gradient-based methods.
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field
from antibody_model import (
    AntibodyDesignModel, VOCAB_SIZE, AMINO_ACIDS,
    encode_sequence, decode_sequence, PAD_IDX
)


# ─────────────────────────────────────────
# Objective weights
# ─────────────────────────────────────────
@dataclass
class OptimizationWeights:
    affinity: float = 0.35      # max (lower log_Kd = better)
    stability: float = 0.20     # max (higher Tm = better)
    humanization: float = 0.20  # max
    immunogenicity: float = 0.10  # min (lower = better)
    aggregation: float = 0.15   # min (lower = better)


@dataclass
class OptimizationCandidate:
    sequence: str
    scores: Dict[str, float] = field(default_factory=dict)
    composite_score: float = 0.0
    generation: int = 0


# ─────────────────────────────────────────
# 1. Pareto Front Computation
# ─────────────────────────────────────────
def compute_pareto_front(
    candidates: List[OptimizationCandidate],
    objectives: List[str],
    maximize: List[bool],
) -> List[OptimizationCandidate]:
    """
    Compute the Pareto-optimal front from a list of candidates.
    Returns the non-dominated solutions.
    """
    n = len(candidates)
    dominated = [False] * n

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # Check if j dominates i
            j_dominates = True
            strictly_better = False
            for obj, maximize_flag in zip(objectives, maximize):
                vi = candidates[i].scores.get(obj, 0.0)
                vj = candidates[j].scores.get(obj, 0.0)
                if maximize_flag:
                    if vj < vi:
                        j_dominates = False
                        break
                    if vj > vi:
                        strictly_better = True
                else:
                    if vj > vi:
                        j_dominates = False
                        break
                    if vj < vi:
                        strictly_better = True
            if j_dominates and strictly_better:
                dominated[i] = True
                break

    return [c for c, d in zip(candidates, dominated) if not d]


# ─────────────────────────────────────────
# 2. Evolutionary Optimization (NSGA-II inspired)
# ─────────────────────────────────────────
class GeneticOptimizer:
    """
    NSGA-II-inspired multi-objective genetic optimizer for CDR sequences.
    Operates on discrete amino acid sequences.
    """

    def __init__(
        self,
        pop_size: int = 100,
        elite_size: int = 20,
        mutation_rate: float = 0.15,
        crossover_rate: float = 0.70,
        tournament_size: int = 5,
        weights: Optional[OptimizationWeights] = None,
    ):
        self.pop_size = pop_size
        self.elite_size = elite_size
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.tournament_size = tournament_size
        self.weights = weights or OptimizationWeights()
        self.rng = np.random.default_rng(42)

    def initialize_population(
        self, seed_sequences: List[str], target_len: int
    ) -> List[str]:
        """Initialize population from seed sequences with random variants."""
        population = list(seed_sequences[:self.pop_size])
        aa_list = list(AMINO_ACIDS)
        while len(population) < self.pop_size:
            # Random sequence with length near target
            length = target_len + self.rng.integers(-2, 3)
            length = max(5, min(length, 25))
            seq = "".join(self.rng.choice(aa_list) for _ in range(length))
            population.append(seq)
        return population[:self.pop_size]

    def mutate(self, seq: str) -> str:
        """Point mutation, insertion, deletion."""
        aa_list = list(AMINO_ACIDS)
        seq_list = list(seq)
        for i in range(len(seq_list) - 1, -1, -1):
            if self.rng.random() < self.mutation_rate:
                op = self.rng.choice(["substitute", "insert", "delete"])
                if op == "substitute" or len(seq_list) <= 5:
                    seq_list[i] = self.rng.choice(aa_list)
                elif op == "insert" and len(seq_list) < 25:
                    seq_list.insert(i, self.rng.choice(aa_list))
                elif op == "delete" and len(seq_list) > 5:
                    seq_list.pop(i)
        return "".join(seq_list)

    def crossover(self, seq1: str, seq2: str) -> Tuple[str, str]:
        """Single-point crossover."""
        if self.rng.random() > self.crossover_rate:
            return seq1, seq2
        min_len = min(len(seq1), len(seq2))
        if min_len < 2:
            return seq1, seq2
        point = self.rng.integers(1, min_len)
        child1 = seq1[:point] + seq2[point:]
        child2 = seq2[:point] + seq1[point:]
        return child1, child2

    def tournament_select(self, candidates: List[OptimizationCandidate]) -> OptimizationCandidate:
        """Tournament selection based on composite score."""
        contestants = self.rng.choice(candidates, size=self.tournament_size, replace=False)
        return max(contestants, key=lambda c: c.composite_score)

    def compute_composite_score(self, scores: dict) -> float:
        """Weighted sum of normalized objectives."""
        w = self.weights
        affinity_score = 1.0 / (1.0 + math.exp(scores.get("log_kd", 0.0)))
        stability_score = min(scores.get("tm", 60) / 90.0, 1.0)
        human_score = scores.get("humanization", 0.5)
        immuno_score = 1.0 - scores.get("immunogenicity", 0.5)
        agg_score = 1.0 - scores.get("aggregation", 0.5)

        return (
            w.affinity * affinity_score
            + w.stability * stability_score
            + w.humanization * human_score
            + w.immunogenicity * immuno_score
            + w.aggregation * agg_score
        )

    def evolve_generation(
        self,
        population: List[OptimizationCandidate],
    ) -> List[OptimizationCandidate]:
        """Produce next generation via selection, crossover, mutation."""
        # Sort by composite score
        population.sort(key=lambda c: c.composite_score, reverse=True)

        # Elites survive unchanged
        next_gen = population[:self.elite_size]

        # Fill rest via crossover + mutation
        while len(next_gen) < self.pop_size:
            p1 = self.tournament_select(population)
            p2 = self.tournament_select(population)
            c1_seq, c2_seq = self.crossover(p1.sequence, p2.sequence)
            c1_seq = self.mutate(c1_seq)
            c2_seq = self.mutate(c2_seq)
            next_gen.append(OptimizationCandidate(sequence=c1_seq))
            if len(next_gen) < self.pop_size:
                next_gen.append(OptimizationCandidate(sequence=c2_seq))

        return next_gen[:self.pop_size]


# ─────────────────────────────────────────
# 3. Gradient-Based Sequence Refinement (soft relaxation)
# ─────────────────────────────────────────
class SoftSequenceOptimizer:
    """
    Gradient-based optimization on soft amino acid probabilities.
    Uses straight-through Gumbel-softmax for backprop through discrete tokens.
    """

    def __init__(
        self,
        model: AntibodyDesignModel,
        cdr_length: int = 12,
        lr: float = 0.01,
        temperature: float = 0.5,
        n_steps: int = 200,
        device: str = "cpu",
    ):
        self.model = model.to(device)
        self.cdr_length = cdr_length
        self.lr = lr
        self.temperature = temperature
        self.n_steps = n_steps
        self.device = device

    def optimize(
        self,
        antigen_tokens: torch.Tensor,
        framework_tokens: torch.Tensor,
        seed_logits: Optional[torch.Tensor] = None,
        target_log_kd: float = -10.0,
    ) -> Tuple[str, List[float]]:
        """
        Gradient-based CDR optimization using soft token relaxation.
        Returns: (best_sequence, loss_history)
        """
        if seed_logits is None:
            logits = torch.randn(
                1, self.cdr_length, VOCAB_SIZE,
                device=self.device, requires_grad=True
            )
        else:
            logits = seed_logits.clone().requires_grad_(True)

        optimizer = torch.optim.Adam([logits], lr=self.lr)
        ag_enc = self.model.encoder(antigen_tokens.to(self.device))
        fw_enc = self.model.encoder(framework_tokens.to(self.device))

        loss_history = []
        best_logits = logits.detach().clone()
        best_loss = float("inf")

        for step in range(self.n_steps):
            optimizer.zero_grad()

            # Gumbel-softmax: straight-through gradient estimator
            soft_tokens = F.gumbel_softmax(logits, tau=self.temperature, hard=False)
            # Project soft tokens to embedding space
            with torch.no_grad():
                emb_weight = self.model.encoder.seq_emb.weight[:VOCAB_SIZE]
            cdr_emb = soft_tokens @ emb_weight  # (1, L, d)

            # Predict affinity via cross-attention
            ctx, _ = self.model.affinity.cross_attn(cdr_emb, ag_enc, ag_enc)
            cdr_pool = ctx.mean(dim=1)
            ag_pool = ag_enc.mean(dim=1)
            feat = torch.cat([cdr_pool, ag_pool], dim=-1)
            log_kd_pred = self.model.affinity.regressor(feat).squeeze()

            # Predict stability
            stab = self.model.stability(cdr_emb)
            tm_pred = stab[:, 1]

            # Multi-objective loss
            affinity_loss = F.mse_loss(log_kd_pred, torch.tensor(target_log_kd, device=self.device))
            stability_loss = F.relu(70.0 - tm_pred).mean()  # penalize Tm < 70°C
            entropy_reg = -torch.mean(
                torch.sum(F.softmax(logits, dim=-1) * F.log_softmax(logits, dim=-1), dim=-1)
            ) * 0.01

            loss = affinity_loss + 0.3 * stability_loss + entropy_reg
            loss.backward()
            optimizer.step()

            loss_val = loss.item()
            loss_history.append(loss_val)
            if loss_val < best_loss:
                best_loss = loss_val
                best_logits = logits.detach().clone()

        # Decode best sequence
        best_indices = best_logits.squeeze(0).argmax(dim=-1)
        best_seq = decode_sequence(best_indices)
        return best_seq, loss_history


# ─────────────────────────────────────────
# 4. Multi-Attribute Optimization Runner
# ─────────────────────────────────────────
def run_multi_attribute_optimization(
    model: AntibodyDesignModel,
    antigen_tokens: torch.Tensor,
    framework_tokens: torch.Tensor,
    seed_sequences: List[str],
    n_generations: int = 50,
    pop_size: int = 50,
    device: str = "cpu",
) -> Dict:
    """
    Run evolutionary + gradient optimization.
    Returns: dict with best sequences, Pareto front, history.
    """
    import random

    genetic_opt = GeneticOptimizer(
        pop_size=pop_size,
        weights=OptimizationWeights(),
    )

    # Initialize population
    target_len = len(seed_sequences[0]) if seed_sequences else 12
    population_seqs = genetic_opt.initialize_population(seed_sequences, target_len)

    ag_enc = model.encoder(antigen_tokens.to(device))
    fw_enc = model.encoder(framework_tokens.to(device))

    best_scores_history = []

    def score_sequence(seq: str) -> dict:
        """Score a single sequence using model."""
        tok = encode_sequence(seq)
        tok_padded = F.pad(tok, (0, max(0, 25 - len(tok))), value=PAD_IDX).unsqueeze(0).to(device)
        ag_exp = antigen_tokens.to(device)
        with torch.no_grad():
            props = model.predict_properties(tok_padded, ag_exp.expand(1, -1))
            log_kd = props["log_kd"].item()
            tm = props["Tm"].item()
            ddg = props["delta_delta_G"].item()
        return {
            "log_kd": log_kd,
            "tm": tm,
            "delta_delta_G": ddg,
            "humanization": random.uniform(0.5, 0.9),   # placeholder
            "immunogenicity": random.uniform(0.1, 0.4),  # placeholder
            "aggregation": random.uniform(0.1, 0.4),     # placeholder
        }

    for gen in range(n_generations):
        # Score population
        candidates = []
        for seq in population_seqs:
            scores = score_sequence(seq)
            c = OptimizationCandidate(sequence=seq, scores=scores, generation=gen)
            c.composite_score = genetic_opt.compute_composite_score(scores)
            candidates.append(c)

        # Record best
        best = max(candidates, key=lambda c: c.composite_score)
        best_scores_history.append(best.composite_score)

        # Evolve
        candidates = genetic_opt.evolve_generation(candidates)
        population_seqs = [c.sequence for c in candidates]

    # Final scoring + Pareto front
    final_candidates = []
    for seq in population_seqs:
        scores = score_sequence(seq)
        c = OptimizationCandidate(sequence=seq, scores=scores)
        c.composite_score = genetic_opt.compute_composite_score(scores)
        final_candidates.append(c)

    objectives = ["log_kd", "tm", "humanization", "immunogenicity", "aggregation"]
    maximize_flags = [False, True, True, False, False]
    pareto = compute_pareto_front(final_candidates, objectives, maximize_flags)

    # Sort final by composite score
    final_candidates.sort(key=lambda c: c.composite_score, reverse=True)

    return {
        "best_candidates": final_candidates[:20],
        "pareto_front": pareto,
        "best_scores_history": best_scores_history,
        "n_generations": n_generations,
    }


if __name__ == "__main__":
    print("=== Multi-Attribute Optimization Test ===")
    from antibody_model import AntibodyDesignModel, VOCAB_SIZE

    device = "cpu"
    model = AntibodyDesignModel(d_model=128, T=50).to(device)
    model.eval()

    ag_tok = torch.randint(0, VOCAB_SIZE, (1, 50))
    fw_tok = torch.randint(0, VOCAB_SIZE, (1, 80))
    seeds = ["ARSGYDGFDY", "YWYCARDLGYYY", "GSSYSGFFDYY"]

    result = run_multi_attribute_optimization(
        model, ag_tok, fw_tok, seeds,
        n_generations=10, pop_size=20, device=device
    )

    print(f"  Pareto front size: {len(result['pareto_front'])}")
    print(f"  Best composite score: {result['best_candidates'][0].composite_score:.4f}")
    print(f"  Top-5 sequences:")
    for i, c in enumerate(result["best_candidates"][:5]):
        print(f"    {i+1}. {c.sequence:25s} score={c.composite_score:.4f}")
