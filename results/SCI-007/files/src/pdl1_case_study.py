"""Synthetic PD-L1 antibody design case study using the diffusion toolkit."""
from __future__ import annotations

from pprint import pprint
from typing import Dict, List

import torch

from data_utils import AntibodyCDRDataset, compute_sequence_properties, simulate_structure_features
from diffusion_model import DiscreteDiffusionModel
from optimizer import generate_pareto_candidates
from property_predictor import MultiPropertyOptimizer

REFERENCE_ANTIBODIES = {
    "atezolizumab": "ARYFDYWGQGTL",
    "durvalumab": "ARGYGSYFDY",
}


def simulate_pdl1_interface_features(batch_size: int = 1, feature_dim: int = 32, seed: int = 17) -> torch.Tensor:
    """Create synthetic PD-L1 interface features capturing hydrophobic and electrostatic patches."""
    torch.manual_seed(seed)
    base = torch.zeros(feature_dim, dtype=torch.float32)
    positions = torch.arange(feature_dim, dtype=torch.float32)
    base += 0.35 * torch.sin(positions / 2.3)
    base += 0.25 * torch.cos(positions / 3.7)
    base[:6] += torch.tensor([0.8, 0.7, 0.6, -0.2, 0.4, 0.5])
    base[6:12] += torch.tensor([-0.5, -0.4, -0.3, 0.1, 0.2, -0.1])
    base[12:18] += torch.tensor([0.45, 0.35, 0.25, 0.2, 0.15, 0.1])
    return base.unsqueeze(0).repeat(batch_size, 1)


def _sequence_identity(query: str, reference: str) -> float:
    length = max(len(query), len(reference), 1)
    matches = sum(1 for left, right in zip(query, reference) if left == right)
    return matches / length


def _kmer_overlap(query: str, reference: str, k: int = 2) -> float:
    query_kmers = {query[i : i + k] for i in range(max(len(query) - k + 1, 1))}
    ref_kmers = {reference[i : i + k] for i in range(max(len(reference) - k + 1, 1))}
    union = len(query_kmers | ref_kmers) or 1
    return len(query_kmers & ref_kmers) / union


def estimate_binding_energy(sequence: str, interface_features: torch.Tensor) -> float:
    """Estimate a synthetic binding energy; more negative is better."""
    props = compute_sequence_properties(sequence)
    interface_hydrophobic = float(interface_features[:6].mean().item())
    interface_charge = float(interface_features[6:12].mean().item())
    aromatic_bonus = 7.5 * props["aromatic_fraction"]
    charge_match = -abs(props["charge"] + 3.0 * interface_charge)
    length_term = -0.25 * abs(len(sequence) - 13)
    hydrophobic_match = 2.4 * props["hydrophobic_fraction"] * max(interface_hydrophobic, 0.0)
    return -(5.0 + aromatic_bonus + hydrophobic_match + charge_match + length_term)


def structural_compatibility(sequence: str) -> float:
    """Estimate loop compatibility with a flat protein-protein interface."""
    structure = simulate_structure_features(sequence)
    coords = structure["relative_positions"]
    if coords.numel() == 0:
        return 0.0
    compactness = torch.norm(coords - coords.mean(dim=0, keepdim=True), dim=-1).mean().item()
    phi_smoothness = structure["phi"].diff().abs().mean().item() if len(sequence) > 1 else 0.0
    length_score = max(0.0, 1.0 - abs(len(sequence) - 13) / 10.0)
    compactness_score = max(0.0, 1.0 - compactness / 10.0)
    torsion_score = max(0.0, 1.0 - phi_smoothness / 90.0)
    return 0.4 * length_score + 0.35 * compactness_score + 0.25 * torsion_score


def assess_developability(sequence: str, property_optimizer: MultiPropertyOptimizer) -> Dict[str, float]:
    """Evaluate developability-related metrics for a single candidate."""
    metrics = property_optimizer.evaluate([sequence])
    expression = float(metrics["expression_level"][0].item())
    aggregation = float(metrics["aggregation_propensity"][0].item())
    humanness = float(metrics["humanness"][0].item())
    composite = 0.45 * expression + 0.35 * (1.0 - aggregation) + 0.20 * humanness
    return {
        "expression_level": expression,
        "aggregation_propensity": aggregation,
        "humanness": humanness,
        "developability_score": composite,
    }


def compare_with_references(sequence: str) -> Dict[str, Dict[str, float]]:
    """Compare a candidate with surrogate PD-L1 therapeutic reference loops."""
    return {
        name: {
            "sequence_identity": _sequence_identity(sequence, reference),
            "motif_overlap": _kmer_overlap(sequence, reference),
        }
        for name, reference in REFERENCE_ANTIBODIES.items()
    }


def _warm_start_diffusion(model: DiscreteDiffusionModel, steps: int = 25, batch_size: int = 32, seed: int = 17) -> None:
    dataset = AntibodyCDRDataset(num_samples=max(batch_size * 4, 128), seed=seed)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    device = next(model.parameters()).device
    generator = torch.Generator(device="cpu").manual_seed(seed)
    for _ in range(steps):
        indices = torch.randint(0, len(dataset), (batch_size,), generator=generator)
        items = [dataset[int(index)] for index in indices]
        batch = {
            "tokens": torch.stack([item["tokens"] for item in items]).to(device),
            "mask": torch.stack([item["mask"] for item in items]).to(device),
            "coords": torch.stack([item["coords"] for item in items]).to(device),
            "antigen_features": torch.stack([item["antigen_features"] for item in items]).to(device),
        }
        model.training_step(batch, optimizer)


def validate_candidates(
    candidates: List[Dict[str, object]],
    interface_features: torch.Tensor,
    property_optimizer: MultiPropertyOptimizer,
) -> List[Dict[str, object]]:
    """Run in silico validation on generated candidates."""
    validated: List[Dict[str, object]] = []
    for candidate in candidates:
        sequence = str(candidate["sequence"])
        binding_energy = estimate_binding_energy(sequence, interface_features[0])
        compatibility = structural_compatibility(sequence)
        developability = assess_developability(sequence, property_optimizer)
        reference_similarity = compare_with_references(sequence)
        aggregate_similarity = sum(metric["motif_overlap"] for metric in reference_similarity.values()) / len(reference_similarity)
        final_score = (
            0.35 * candidate["weighted_score"]
            + 0.25 * developability["developability_score"]
            + 0.20 * compatibility
            + 0.10 * (-binding_energy / 15.0)
            + 0.10 * aggregate_similarity
        )
        validated.append(
            {
                **candidate,
                "binding_energy": binding_energy,
                "structural_compatibility": compatibility,
                "developability": developability,
                "reference_similarity": reference_similarity,
                "validation_score": final_score,
            }
        )
    validated.sort(key=lambda item: item["validation_score"], reverse=True)
    return validated


def generate_pdl1_candidates(num_candidates: int = 48, train_steps: int = 25, seed: int = 17) -> List[Dict[str, object]]:
    """Generate Pareto-optimized antibody candidates for a simulated PD-L1 interface."""
    device = torch.device("cpu")
    interface_features = simulate_pdl1_interface_features(batch_size=num_candidates, seed=seed)
    model = DiscreteDiffusionModel(condition_dim=interface_features.size(-1)).to(device)
    property_optimizer = MultiPropertyOptimizer().to(device)
    _warm_start_diffusion(model, steps=train_steps, seed=seed)
    return generate_pareto_candidates(
        diffusion_model=model,
        property_optimizer=property_optimizer,
        num_candidates=num_candidates,
        num_return=min(12, num_candidates),
        condition=interface_features.to(device),
        guidance_scale=1.6,
        device=device,
    )


def run_pdl1_case_study(num_candidates: int = 48, train_steps: int = 25, seed: int = 17) -> Dict[str, object]:
    """Run a fully synthetic PD-L1 antibody discovery workflow."""
    device = torch.device("cpu")
    interface_features = simulate_pdl1_interface_features(batch_size=num_candidates, seed=seed).to(device)
    model = DiscreteDiffusionModel(condition_dim=interface_features.size(-1)).to(device)
    property_optimizer = MultiPropertyOptimizer().to(device)
    _warm_start_diffusion(model, steps=train_steps, seed=seed)
    candidates = generate_pareto_candidates(
        diffusion_model=model,
        property_optimizer=property_optimizer,
        num_candidates=num_candidates,
        num_return=min(12, num_candidates),
        condition=interface_features,
        guidance_scale=1.6,
        device=device,
    )
    validated = validate_candidates(candidates, interface_features, property_optimizer)
    top_binding = sum(item["binding_energy"] for item in validated[:5]) / max(len(validated[:5]), 1)
    top_compatibility = sum(item["structural_compatibility"] for item in validated[:5]) / max(len(validated[:5]), 1)
    return {
        "target": "PD-L1",
        "num_candidates_generated": num_candidates,
        "num_pareto_selected": len(candidates),
        "reference_antibodies": REFERENCE_ANTIBODIES,
        "average_top5_binding_energy": top_binding,
        "average_top5_structural_compatibility": top_compatibility,
        "top_candidates": validated[:5],
    }


if __name__ == "__main__":
    pprint(run_pdl1_case_study(num_candidates=24, train_steps=10))
