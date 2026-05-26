"""Utilities for synthetic antibody CDR-H3 data generation and analysis."""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence

import torch
from torch.utils.data import Dataset

AMINO_ACIDS: List[str] = list("ACDEFGHIKLMNPQRSTVWY")
PAD_TOKEN = "-"
UNK_TOKEN = "X"
SPECIAL_TOKENS = [PAD_TOKEN, UNK_TOKEN]
AA_VOCAB: List[str] = AMINO_ACIDS + SPECIAL_TOKENS
AA_TO_IDX: Dict[str, int] = {aa: idx for idx, aa in enumerate(AA_VOCAB)}
IDX_TO_AA: Dict[int, str] = {idx: aa for aa, idx in AA_TO_IDX.items()}
AA_VOCAB_SIZE = len(AA_VOCAB)
PAD_IDX = AA_TO_IDX[PAD_TOKEN]
UNK_IDX = AA_TO_IDX[UNK_TOKEN]
MIN_CDR_H3_LEN = 8
MAX_CDR_H3_LEN = 20

CHARGE_TABLE: Dict[str, float] = {
    "A": 0.0,
    "C": 0.0,
    "D": -1.0,
    "E": -1.0,
    "F": 0.0,
    "G": 0.0,
    "H": 0.1,
    "I": 0.0,
    "K": 1.0,
    "L": 0.0,
    "M": 0.0,
    "N": 0.0,
    "P": 0.0,
    "Q": 0.0,
    "R": 1.0,
    "S": 0.0,
    "T": 0.0,
    "V": 0.0,
    "W": 0.0,
    "Y": 0.0,
    PAD_TOKEN: 0.0,
    UNK_TOKEN: 0.0,
}
HYDROPHOBICITY_TABLE: Dict[str, float] = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
    PAD_TOKEN: 0.0,
    UNK_TOKEN: 0.0,
}
MOLECULAR_WEIGHT_TABLE: Dict[str, float] = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "Y": 181.19,
    PAD_TOKEN: 0.0,
    UNK_TOKEN: 136.90,
}

_POSITION_WEIGHTS = {
    "start": {
        "A": 1.5,
        "R": 1.4,
        "G": 1.1,
        "S": 1.0,
        "T": 0.8,
        "D": 0.7,
        "Y": 0.5,
        "all": 0.25,
    },
    "middle": {
        "Y": 1.6,
        "G": 1.6,
        "D": 1.1,
        "S": 1.1,
        "W": 1.0,
        "F": 1.0,
        "N": 0.9,
        "R": 0.8,
        "A": 0.7,
        "all": 0.35,
    },
    "end": {
        "Y": 1.6,
        "D": 1.5,
        "F": 1.3,
        "W": 0.8,
        "G": 0.7,
        "R": 0.5,
        "all": 0.2,
    },
}
_COMMON_HUMAN_LIKE = set("AGSDTYV")
_AROMATIC = set("FWY")
_HYDROPHOBIC = set("AILMFWVY")


def _weighted_choice(weights: Dict[str, float], rng: random.Random) -> str:
    expanded = {aa: weights.get("all", 0.2) for aa in AMINO_ACIDS}
    for aa, weight in weights.items():
        if aa != "all":
            expanded[aa] = expanded.get(aa, 0.0) + weight
    total = sum(expanded.values())
    threshold = rng.random() * total
    running = 0.0
    for aa in AMINO_ACIDS:
        running += expanded[aa]
        if running >= threshold:
            return aa
    return "Y"


def encode_sequence(sequence: str, max_length: int = MAX_CDR_H3_LEN) -> torch.Tensor:
    """Encode a CDR-H3 sequence into token indices with right padding."""
    sequence = sequence.upper().strip()
    tokens = [AA_TO_IDX.get(residue, UNK_IDX) for residue in sequence[:max_length]]
    if len(tokens) < max_length:
        tokens.extend([PAD_IDX] * (max_length - len(tokens)))
    return torch.tensor(tokens, dtype=torch.long)


def decode_sequence(tokens: Sequence[int] | torch.Tensor, strip_pad: bool = True) -> str:
    """Decode token indices back into an amino-acid string."""
    if isinstance(tokens, torch.Tensor):
        tokens = tokens.detach().cpu().tolist()
    residues = [IDX_TO_AA.get(int(token), UNK_TOKEN) for token in tokens]
    sequence = "".join(residues)
    return sequence.replace(PAD_TOKEN, "") if strip_pad else sequence


def compute_sequence_charge(sequence: str) -> float:
    """Compute an approximate net charge at neutral pH."""
    side_chain_charge = sum(CHARGE_TABLE.get(residue, 0.0) for residue in sequence)
    terminal_charge = 0.1 if sequence else 0.0
    return side_chain_charge + terminal_charge


def compute_hydrophobicity(sequence: str) -> float:
    """Compute average Kyte-Doolittle hydrophobicity."""
    if not sequence:
        return 0.0
    return sum(HYDROPHOBICITY_TABLE.get(residue, 0.0) for residue in sequence) / len(sequence)


def compute_molecular_weight(sequence: str) -> float:
    """Compute an approximate molecular weight in Daltons."""
    if not sequence:
        return 0.0
    weight = sum(MOLECULAR_WEIGHT_TABLE.get(residue, MOLECULAR_WEIGHT_TABLE[UNK_TOKEN]) for residue in sequence)
    water_loss = 18.015 * max(len(sequence) - 1, 0)
    return weight - water_loss


def compute_sequence_properties(sequence: str) -> Dict[str, float]:
    """Compute simple physicochemical properties for a sequence."""
    length = max(len(sequence), 1)
    aromatic_fraction = sum(residue in _AROMATIC for residue in sequence) / length
    hydrophobic_fraction = sum(residue in _HYDROPHOBIC for residue in sequence) / length
    return {
        "length": float(len(sequence)),
        "charge": compute_sequence_charge(sequence),
        "hydrophobicity": compute_hydrophobicity(sequence),
        "molecular_weight": compute_molecular_weight(sequence),
        "aromatic_fraction": aromatic_fraction,
        "hydrophobic_fraction": hydrophobic_fraction,
    }


def generate_synthetic_cdr_h3(length: Optional[int] = None, rng: Optional[random.Random] = None) -> str:
    """Generate a realistic-looking synthetic CDR-H3 sequence."""
    rng = rng or random.Random()
    length = length or rng.randint(MIN_CDR_H3_LEN, MAX_CDR_H3_LEN)
    residues: List[str] = []
    for index in range(length):
        if index == 0:
            residue = _weighted_choice(_POSITION_WEIGHTS["start"], rng)
        elif index == length - 1:
            residue = _weighted_choice(_POSITION_WEIGHTS["end"], rng)
        else:
            residue = _weighted_choice(_POSITION_WEIGHTS["middle"], rng)
        residues.append(residue)
    if length >= 11 and rng.random() < 0.65:
        motif = rng.choice(["GYG", "YFD", "GWD", "DYG"])
        start = rng.randint(2, max(2, length - len(motif) - 2))
        residues[start : start + len(motif)] = list(motif)
    cysteine_positions = [i for i, aa in enumerate(residues) if aa == "C"]
    for pos in cysteine_positions[1:]:
        residues[pos] = rng.choice(["S", "A", "Y"])
    return "".join(residues)


def simulate_structure_features(sequence: str) -> Dict[str, torch.Tensor]:
    """Simulate coarse torsion and coordinate features for a loop sequence."""
    seed = sum((index + 1) * ord(residue) for index, residue in enumerate(sequence))
    rng = random.Random(seed)
    phi: List[float] = []
    psi: List[float] = []
    positions: List[List[float]] = []
    current = [0.0, 0.0, 0.0]
    direction = 0.0
    for index, residue in enumerate(sequence):
        flexibility = 1.4 if residue in "GSTPN" else 0.8
        phi_angle = -65.0 + 30.0 * math.sin(0.31 * (index + 1) + rng.random()) + 10.0 * flexibility
        psi_angle = 135.0 * math.cos(0.23 * (index + 1) + rng.random())
        phi.append(phi_angle)
        psi.append(psi_angle)
        step = 1.2 + 0.15 * flexibility + 0.05 * rng.random()
        direction += math.radians(phi_angle / 3.0)
        current = [
            current[0] + step * math.cos(direction),
            current[1] + 0.35 * math.sin(math.radians(psi_angle)),
            current[2] + step * math.sin(direction),
        ]
        positions.append(current.copy())
    coords = torch.tensor(positions, dtype=torch.float32)
    if coords.numel() > 0:
        coords = coords - coords[0]
    return {
        "phi": torch.tensor(phi, dtype=torch.float32),
        "psi": torch.tensor(psi, dtype=torch.float32),
        "relative_positions": coords,
    }


def _simulate_affinity(sequence: str, antigen_features: torch.Tensor) -> float:
    props = compute_sequence_properties(sequence)
    antigen_bias = float(antigen_features[:4].mean().item())
    return 6.0 + 2.4 * props["aromatic_fraction"] + 0.6 * props["hydrophobic_fraction"] - 0.2 * abs(props["charge"]) + antigen_bias


def _simulate_stability(sequence: str) -> float:
    props = compute_sequence_properties(sequence)
    gly_pro_penalty = 7.0 * (sequence.count("G") + 1.5 * sequence.count("P")) / max(len(sequence), 1)
    return 62.0 + 7.5 * props["hydrophobic_fraction"] - gly_pro_penalty


def _simulate_humanness(sequence: str) -> float:
    human_like_fraction = sum(residue in _COMMON_HUMAN_LIKE for residue in sequence) / max(len(sequence), 1)
    cysteine_penalty = 0.25 * max(sequence.count("C") - 1, 0)
    return max(0.0, min(1.0, 0.45 + 0.7 * human_like_fraction - cysteine_penalty))


def _simulate_immunogenicity(sequence: str) -> float:
    rare_fraction = sum(residue in {"C", "M", "W"} for residue in sequence) / max(len(sequence), 1)
    return max(0.0, min(1.0, 0.2 + 0.9 * rare_fraction + 0.05 * abs(compute_sequence_charge(sequence))))


def _simulate_developability(sequence: str) -> Dict[str, float]:
    props = compute_sequence_properties(sequence)
    expression = max(0.0, min(1.0, 0.8 - 0.08 * abs(props["charge"]) + 0.05 * props["hydrophobic_fraction"]))
    aggregation = max(0.0, min(1.0, 0.15 + 0.8 * props["hydrophobic_fraction"] + 0.2 * props["aromatic_fraction"]))
    return {"expression_level": expression, "aggregation_propensity": aggregation}


def generate_synthetic_antibody_dataset(
    num_samples: int = 1000,
    seed: int = 7,
    antigen_feature_dim: int = 32,
) -> List[Dict[str, object]]:
    """Generate a synthetic dataset of antibody loop sequences and labels."""
    rng = random.Random(seed)
    dataset: List[Dict[str, object]] = []
    for _ in range(num_samples):
        sequence = generate_synthetic_cdr_h3(rng=rng)
        structure_features = simulate_structure_features(sequence)
        antigen_features = torch.tensor(
            [rng.uniform(-1.0, 1.0) for _ in range(antigen_feature_dim)],
            dtype=torch.float32,
        )
        affinity = _simulate_affinity(sequence, antigen_features)
        stability = _simulate_stability(sequence)
        humanness = _simulate_humanness(sequence)
        immunogenicity = _simulate_immunogenicity(sequence)
        developability = _simulate_developability(sequence)
        dataset.append(
            {
                "sequence": sequence,
                "length": len(sequence),
                "encoded": encode_sequence(sequence),
                "structure_features": structure_features,
                "antigen_features": antigen_features,
                "binding_affinity": affinity,
                "stability_tm": stability,
                "humanness": humanness,
                "immunogenicity": immunogenicity,
                **developability,
            }
        )
    return dataset


class AntibodyCDRDataset(Dataset):
    """PyTorch dataset for synthetic antibody CDR-H3 sequences."""

    def __init__(
        self,
        samples: Optional[List[Dict[str, object]]] = None,
        num_samples: int = 1000,
        max_length: int = MAX_CDR_H3_LEN,
        seed: int = 7,
    ) -> None:
        self.max_length = max_length
        self.samples = samples or generate_synthetic_antibody_dataset(num_samples=num_samples, seed=seed)

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, index: int) -> Dict[str, torch.Tensor | str]:
        sample = self.samples[index]
        sequence = str(sample["sequence"])
        structure = sample["structure_features"]
        length = min(len(sequence), self.max_length)
        coords = torch.zeros(self.max_length, 3, dtype=torch.float32)
        phi = torch.zeros(self.max_length, dtype=torch.float32)
        psi = torch.zeros(self.max_length, dtype=torch.float32)
        coords[:length] = structure["relative_positions"][:length]
        phi[:length] = structure["phi"][:length]
        psi[:length] = structure["psi"][:length]
        mask = torch.zeros(self.max_length, dtype=torch.float32)
        mask[:length] = 1.0
        return {
            "tokens": encode_sequence(sequence, max_length=self.max_length),
            "mask": mask,
            "coords": coords,
            "phi": phi,
            "psi": psi,
            "length": torch.tensor(length, dtype=torch.long),
            "antigen_features": sample["antigen_features"].float(),
            "binding_affinity": torch.tensor(float(sample["binding_affinity"]), dtype=torch.float32),
            "stability_tm": torch.tensor(float(sample["stability_tm"]), dtype=torch.float32),
            "humanness": torch.tensor(float(sample["humanness"]), dtype=torch.float32),
            "immunogenicity": torch.tensor(float(sample["immunogenicity"]), dtype=torch.float32),
            "expression_level": torch.tensor(float(sample["expression_level"]), dtype=torch.float32),
            "aggregation_propensity": torch.tensor(float(sample["aggregation_propensity"]), dtype=torch.float32),
            "sequence": sequence,
        }


__all__ = [
    "AA_TO_IDX",
    "AA_VOCAB",
    "AA_VOCAB_SIZE",
    "AMINO_ACIDS",
    "AntibodyCDRDataset",
    "CHARGE_TABLE",
    "HYDROPHOBICITY_TABLE",
    "IDX_TO_AA",
    "MAX_CDR_H3_LEN",
    "MIN_CDR_H3_LEN",
    "MOLECULAR_WEIGHT_TABLE",
    "PAD_IDX",
    "PAD_TOKEN",
    "UNK_IDX",
    "UNK_TOKEN",
    "compute_hydrophobicity",
    "compute_molecular_weight",
    "compute_sequence_charge",
    "compute_sequence_properties",
    "decode_sequence",
    "encode_sequence",
    "generate_synthetic_antibody_dataset",
    "generate_synthetic_cdr_h3",
    "simulate_structure_features",
]
