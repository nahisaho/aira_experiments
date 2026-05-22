"""Improved Synthetic Accessibility (SA) Score.

Extends the original SA Score by Ertl & Schuffenhauer with:
- Reaction feasibility scoring
- Starting material availability heuristic
- Stereochemistry complexity penalty
- Ring system complexity analysis
"""

import math
import numpy as np
from typing import Dict, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
from rdkit.Chem import rdchem


def calculate_fragment_score(mol: Chem.Mol) -> float:
    """Fragment contribution to SA score based on ECFP-like fragments."""
    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    bits = fp.GetNonzeroElements()
    n_frags = len(bits)
    if n_frags == 0:
        return 0.0
    freq_sum = sum(min(v, 4) for v in bits.values())
    return freq_sum / n_frags


def calculate_complexity_penalty(mol: Chem.Mol) -> float:
    """Structural complexity penalty based on multiple factors."""
    penalties = 0.0

    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings > 4:
        penalties += (n_rings - 4) * 0.5

    # Spiro and bridged ring penalty
    n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    n_bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    penalties += n_spiro * 0.8
    penalties += n_bridgehead * 1.0

    # Stereocenters penalty
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    n_chiral = len(chiral_centers)
    if n_chiral > 2:
        penalties += (n_chiral - 2) * 0.6

    # Macrocycle penalty
    ring_sizes = [len(r) for r in ring_info.AtomRings()]
    n_macrocycles = sum(1 for s in ring_sizes if s > 8)
    penalties += n_macrocycles * 1.5

    return penalties


def calculate_reaction_feasibility(mol: Chem.Mol) -> float:
    """Score based on presence of well-known synthetic handles."""
    score = 0.0
    smiles = Chem.MolToSmiles(mol)

    # Functional groups that are easy synthetic handles
    easy_handles = [
        ("[NH2]", 0.3),
        ("[OH]", 0.2),
        ("[C](=O)[OH]", 0.3),
        ("[Br]", 0.4),
        ("[Cl]", 0.3),
        ("c1ccccc1", 0.2),  # phenyl
    ]

    for smarts, bonus in easy_handles:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            score += bonus

    return min(score, 2.0)


def calculate_starting_material_score(mol: Chem.Mol) -> float:
    """Heuristic score for starting material availability.

    Small, simple molecules with common functional groups
    are more likely to be commercially available.
    """
    n_atoms = mol.GetNumHeavyAtoms()
    mw = Descriptors.ExactMolWt(mol)

    # Simple molecules are more available
    if n_atoms <= 6:
        size_score = 1.0
    elif n_atoms <= 12:
        size_score = 0.7
    elif n_atoms <= 20:
        size_score = 0.4
    else:
        size_score = 0.1

    # Common element bonus
    atom_types = set(atom.GetSymbol() for atom in mol.GetAtoms())
    common_elements = {"C", "N", "O", "S", "F", "Cl", "Br"}
    if atom_types.issubset(common_elements):
        element_score = 1.0
    else:
        element_score = 0.5

    return (size_score + element_score) / 2.0


def improved_sa_score(mol: Chem.Mol) -> Dict[str, float]:
    """Calculate improved Synthetic Accessibility score.

    Returns a dictionary with component scores and the final SA score.
    Score range: 1 (easy) to 10 (hard)
    """
    if mol is None:
        return {"sa_score": 10.0, "components": {}}

    fragment_score = calculate_fragment_score(mol)
    complexity_penalty = calculate_complexity_penalty(mol)
    reaction_feasibility = calculate_reaction_feasibility(mol)
    starting_material = calculate_starting_material_score(mol)

    n_atoms = mol.GetNumHeavyAtoms()
    size_penalty = max(0, (n_atoms - 20) * 0.1)

    # Molecular weight contribution
    mw = Descriptors.ExactMolWt(mol)
    mw_penalty = max(0, (mw - 500) * 0.005)

    # Rotatable bonds – flexibility penalty
    n_rotatable = Descriptors.NumRotatableBonds(mol)
    rot_penalty = max(0, (n_rotatable - 10) * 0.15)

    # Combine scores
    raw_score = (
        3.0
        - fragment_score * 0.5
        + complexity_penalty
        - reaction_feasibility * 0.8
        - starting_material * 0.5
        + size_penalty
        + mw_penalty
        + rot_penalty
    )

    sa_score = max(1.0, min(10.0, raw_score))

    return {
        "sa_score": round(sa_score, 3),
        "components": {
            "fragment_score": round(fragment_score, 3),
            "complexity_penalty": round(complexity_penalty, 3),
            "reaction_feasibility": round(reaction_feasibility, 3),
            "starting_material_score": round(starting_material, 3),
            "size_penalty": round(size_penalty, 3),
            "mw_penalty": round(mw_penalty, 3),
            "rotatable_bond_penalty": round(rot_penalty, 3),
        },
        "molecular_properties": {
            "heavy_atoms": n_atoms,
            "molecular_weight": round(mw, 2),
            "num_rings": mol.GetRingInfo().NumRings(),
            "num_stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            "num_rotatable_bonds": n_rotatable,
            "logP": round(Crippen.MolLogP(mol), 2),
        },
    }


def sa_score_from_smiles(smiles: str) -> Dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    return improved_sa_score(mol)
