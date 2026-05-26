"""Multi-objective optimization helpers for diffusion-based antibody design."""
from __future__ import annotations

from typing import Dict, List, Optional

import torch

from data_utils import decode_sequence
from property_predictor import MultiPropertyOptimizer


def pareto_dominates(candidate: torch.Tensor, reference: torch.Tensor) -> bool:
    """Return True if candidate dominates reference under maximize-all objectives."""
    return bool(torch.all(candidate >= reference) and torch.any(candidate > reference))


def pareto_front(points: torch.Tensor) -> torch.Tensor:
    """Return indices of the Pareto-optimal set."""
    num_points = points.size(0)
    keep = torch.ones(num_points, dtype=torch.bool, device=points.device)
    for i in range(num_points):
        if not keep[i]:
            continue
        for j in range(num_points):
            if i != j and pareto_dominates(points[j], points[i]):
                keep[i] = False
                break
    return torch.where(keep)[0]


def weighted_sum_scalarization(points: torch.Tensor, weights: torch.Tensor) -> torch.Tensor:
    """Compute weighted linear scalarization scores."""
    weights = weights.to(points.device, dtype=points.dtype)
    return (points * weights).sum(dim=-1)


def non_dominated_sort(points: torch.Tensor) -> List[torch.Tensor]:
    """Simplified NSGA-II non-dominated sorting."""
    num_points = points.size(0)
    dominates = [[] for _ in range(num_points)]
    dominated_count = [0 for _ in range(num_points)]
    fronts: List[List[int]] = [[]]
    for i in range(num_points):
        for j in range(num_points):
            if i == j:
                continue
            if pareto_dominates(points[i], points[j]):
                dominates[i].append(j)
            elif pareto_dominates(points[j], points[i]):
                dominated_count[i] += 1
        if dominated_count[i] == 0:
            fronts[0].append(i)
    current = 0
    while current < len(fronts) and fronts[current]:
        next_front: List[int] = []
        for index in fronts[current]:
            for dominated in dominates[index]:
                dominated_count[dominated] -= 1
                if dominated_count[dominated] == 0:
                    next_front.append(dominated)
        if next_front:
            fronts.append(next_front)
        current += 1
    return [torch.tensor(front, dtype=torch.long, device=points.device) for front in fronts if front]


def crowding_distance(points: torch.Tensor) -> torch.Tensor:
    """Estimate crowding distance within a front."""
    if points.size(0) <= 2:
        return torch.full((points.size(0),), float("inf"), device=points.device)
    distance = torch.zeros(points.size(0), device=points.device)
    for dim in range(points.size(1)):
        values = points[:, dim]
        order = torch.argsort(values)
        distance[order[0]] = float("inf")
        distance[order[-1]] = float("inf")
        span = (values[order[-1]] - values[order[0]]).clamp_min(1e-6)
        for idx in range(1, points.size(0) - 1):
            distance[order[idx]] += (values[order[idx + 1]] - values[order[idx - 1]]) / span
    return distance


def guided_sampling(
    diffusion_model,
    property_optimizer: MultiPropertyOptimizer,
    batch_size: int = 32,
    seq_length: Optional[int] = None,
    condition: Optional[torch.Tensor] = None,
    guidance_scale: float = 1.5,
    temperature: float = 1.0,
    device: Optional[torch.device] = None,
) -> torch.Tensor:
    """Guide diffusion sampling using property gradients on denoiser logits."""
    if device is None:
        device = next(diffusion_model.parameters()).device
    property_optimizer = property_optimizer.to(device)

    def modify_logits(logits: torch.Tensor, _x_t: torch.Tensor, _timestep: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        step = int(_timestep[0].item())
        if step > 40 or step % 5 != 0:
            return logits
        with torch.enable_grad():
            guided_logits = logits.detach().requires_grad_(True)
            score = property_optimizer.weighted_score(guided_logits).mean()
            gradient = torch.autograd.grad(score, guided_logits)[0]
        return (guided_logits + guidance_scale * gradient * mask.unsqueeze(-1)).detach()

    return diffusion_model.sample(
        batch_size=batch_size,
        seq_length=seq_length,
        condition=condition,
        device=device,
        temperature=temperature,
        return_tokens=True,
        logit_modifier=modify_logits,
    )


def generate_pareto_candidates(
    diffusion_model,
    property_optimizer: MultiPropertyOptimizer,
    num_candidates: int = 64,
    num_return: int = 10,
    seq_length: Optional[int] = None,
    condition: Optional[torch.Tensor] = None,
    guidance_scale: float = 1.5,
    device: Optional[torch.device] = None,
) -> List[Dict[str, object]]:
    """Generate candidates and return the top Pareto-optimal subset."""
    if device is None:
        device = next(diffusion_model.parameters()).device
    tokens = guided_sampling(
        diffusion_model=diffusion_model,
        property_optimizer=property_optimizer,
        batch_size=num_candidates,
        seq_length=seq_length,
        condition=condition,
        guidance_scale=guidance_scale,
        device=device,
    )
    metrics = property_optimizer.evaluate(tokens)
    objectives = property_optimizer.objective_vector(tokens)
    fronts = non_dominated_sort(objectives)
    candidates: List[Dict[str, object]] = []
    seen = set()
    for front in fronts:
        scores = property_optimizer.weighted_score(tokens[front])
        ranked = front[torch.argsort(scores, descending=True)]
        for index in ranked.tolist():
            sequence = decode_sequence(tokens[index])
            if sequence in seen:
                continue
            seen.add(sequence)
            candidates.append(
                {
                    "sequence": sequence,
                    "objectives": objectives[index].detach().cpu(),
                    "binding_affinity": float(metrics["binding_affinity"][index].item()),
                    "stability_tm": float(metrics["stability_tm"][index].item()),
                    "humanness": float(metrics["humanness"][index].item()),
                    "immunogenicity": float(metrics["immunogenicity"][index].item()),
                    "expression_level": float(metrics["expression_level"][index].item()),
                    "aggregation_propensity": float(metrics["aggregation_propensity"][index].item()),
                    "weighted_score": float(property_optimizer.weighted_score(tokens[index : index + 1])[0].item()),
                }
            )
            if len(candidates) >= num_return:
                return candidates
    return candidates


__all__ = [
    "crowding_distance",
    "generate_pareto_candidates",
    "guided_sampling",
    "non_dominated_sort",
    "pareto_dominates",
    "pareto_front",
    "weighted_sum_scalarization",
]
