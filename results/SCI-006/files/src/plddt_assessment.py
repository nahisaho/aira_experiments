"""
Module 1: AlphaFold2 pLDDT-based Docking Suitability Assessment

Evaluates AlphaFold2 predicted structures for molecular docking suitability
based on per-residue pLDDT confidence scores, with binding site quality
metrics and adaptive docking strategy selection.
"""

import json
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Tuple, Optional
from enum import Enum


class DockingSuitability(Enum):
    """Docking suitability classification based on pLDDT scores."""
    HIGH = "high"        # pLDDT >= 90: direct docking
    MODERATE = "moderate" # 70 <= pLDDT < 90: docking with caution
    LOW = "low"          # 50 <= pLDDT < 70: requires MD refinement
    UNSUITABLE = "unsuitable"  # pLDDT < 50: not recommended


@dataclass
class ResidueConfidence:
    """Per-residue confidence assessment."""
    residue_id: int
    residue_name: str
    chain_id: str
    plddt: float
    is_binding_site: bool = False
    suitability: str = ""

    def __post_init__(self):
        if not self.suitability:
            self.suitability = classify_plddt(self.plddt).value


@dataclass
class BindingSiteAssessment:
    """Comprehensive binding site quality assessment."""
    site_residues: List[ResidueConfidence]
    mean_plddt: float = 0.0
    min_plddt: float = 0.0
    max_plddt: float = 0.0
    std_plddt: float = 0.0
    fraction_high_confidence: float = 0.0
    fraction_disordered: float = 0.0
    overall_suitability: str = ""
    recommended_strategy: str = ""
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self):
        if self.site_residues:
            self._compute_metrics()

    def _compute_metrics(self):
        plddts = [r.plddt for r in self.site_residues]
        self.mean_plddt = float(np.mean(plddts))
        self.min_plddt = float(np.min(plddts))
        self.max_plddt = float(np.max(plddts))
        self.std_plddt = float(np.std(plddts))
        self.fraction_high_confidence = sum(1 for p in plddts if p >= 90) / len(plddts)
        self.fraction_disordered = sum(1 for p in plddts if p < 50) / len(plddts)
        self.overall_suitability = classify_plddt(self.mean_plddt).value
        self.recommended_strategy = self._determine_strategy()
        self._generate_warnings()

    def _determine_strategy(self) -> str:
        if self.mean_plddt >= 90 and self.min_plddt >= 70:
            return "rigid_docking"
        elif self.mean_plddt >= 70:
            if self.std_plddt > 15:
                return "flexible_docking_with_ensemble"
            return "flexible_docking"
        elif self.mean_plddt >= 50:
            return "md_refinement_then_docking"
        else:
            return "homology_modeling_recommended"

    def _generate_warnings(self):
        if self.fraction_disordered > 0.3:
            self.warnings.append(
                f"WARNING: {self.fraction_disordered:.0%} of binding site residues "
                "have pLDDT < 50 (likely disordered)"
            )
        if self.std_plddt > 20:
            self.warnings.append(
                f"WARNING: High pLDDT variance (σ={self.std_plddt:.1f}) suggests "
                "mixed confidence in binding site"
            )
        if self.min_plddt < 30:
            self.warnings.append(
                "WARNING: Extremely low confidence residues detected in binding site. "
                "Consider experimental structure if available."
            )


def classify_plddt(score: float) -> DockingSuitability:
    """Classify pLDDT score into docking suitability category."""
    if score >= 90:
        return DockingSuitability.HIGH
    elif score >= 70:
        return DockingSuitability.MODERATE
    elif score >= 50:
        return DockingSuitability.LOW
    else:
        return DockingSuitability.UNSUITABLE


def parse_plddt_from_bfactor(pdb_path: str) -> List[ResidueConfidence]:
    """
    Parse pLDDT scores from B-factor column of AlphaFold2 PDB file.
    AlphaFold2 stores pLDDT in the B-factor column for CA atoms.
    """
    residues = []
    seen_residues = set()

    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    chain_id = line[21]
                    res_id = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    plddt = float(line[60:66].strip())

                    key = (chain_id, res_id)
                    if key not in seen_residues:
                        seen_residues.add(key)
                        residues.append(ResidueConfidence(
                            residue_id=res_id,
                            residue_name=res_name,
                            chain_id=chain_id,
                            plddt=plddt
                        ))
    except FileNotFoundError:
        pass

    return residues


def identify_binding_site(
    residues: List[ResidueConfidence],
    binding_site_residue_ids: List[int],
    chain_id: str = "A"
) -> List[ResidueConfidence]:
    """Mark and extract binding site residues."""
    binding_residues = []
    for res in residues:
        if res.residue_id in binding_site_residue_ids and res.chain_id == chain_id:
            res.is_binding_site = True
            binding_residues.append(res)
    return binding_residues


def assess_binding_site(
    residues: List[ResidueConfidence],
    binding_site_residue_ids: List[int],
    chain_id: str = "A"
) -> BindingSiteAssessment:
    """
    Perform comprehensive binding site assessment based on pLDDT scores.
    """
    site_residues = identify_binding_site(residues, binding_site_residue_ids, chain_id)
    return BindingSiteAssessment(site_residues=site_residues)


def compute_local_confidence_map(
    residues: List[ResidueConfidence],
    window_size: int = 5
) -> Dict[int, float]:
    """
    Compute smoothed local confidence scores using sliding window average.
    Useful for identifying confident structural regions.
    """
    plddts = [(r.residue_id, r.plddt) for r in residues]
    plddts.sort(key=lambda x: x[0])

    local_scores = {}
    values = [p[1] for p in plddts]

    for i, (res_id, _) in enumerate(plddts):
        start = max(0, i - window_size // 2)
        end = min(len(values), i + window_size // 2 + 1)
        local_scores[res_id] = float(np.mean(values[start:end]))

    return local_scores


def generate_docking_recommendation(assessment: BindingSiteAssessment) -> Dict:
    """Generate comprehensive docking recommendation report."""
    return {
        "overall_suitability": assessment.overall_suitability,
        "recommended_strategy": assessment.recommended_strategy,
        "binding_site_metrics": {
            "mean_plddt": round(assessment.mean_plddt, 2),
            "min_plddt": round(assessment.min_plddt, 2),
            "max_plddt": round(assessment.max_plddt, 2),
            "std_plddt": round(assessment.std_plddt, 2),
            "fraction_high_confidence": round(assessment.fraction_high_confidence, 3),
            "fraction_disordered": round(assessment.fraction_disordered, 3),
        },
        "strategy_parameters": _get_strategy_parameters(assessment),
        "warnings": assessment.warnings,
    }


def _get_strategy_parameters(assessment: BindingSiteAssessment) -> Dict:
    """Get recommended parameters for the selected docking strategy."""
    strategy = assessment.recommended_strategy
    params = {
        "rigid_docking": {
            "exhaustiveness": 32,
            "num_poses": 9,
            "energy_range": 3.0,
            "flexible_residues": [],
        },
        "flexible_docking": {
            "exhaustiveness": 64,
            "num_poses": 20,
            "energy_range": 5.0,
            "flexible_residues": [
                r.residue_id for r in assessment.site_residues
                if r.plddt < 80
            ],
        },
        "flexible_docking_with_ensemble": {
            "exhaustiveness": 64,
            "num_poses": 50,
            "energy_range": 5.0,
            "num_ensemble_structures": 5,
            "md_equilibration_ns": 10,
            "flexible_residues": [
                r.residue_id for r in assessment.site_residues
                if r.plddt < 80
            ],
        },
        "md_refinement_then_docking": {
            "md_production_ns": 100,
            "clustering_method": "RMSD",
            "num_representative_structures": 10,
            "exhaustiveness": 64,
            "num_poses": 50,
        },
        "homology_modeling_recommended": {
            "note": "AlphaFold2 structure unreliable for this binding site. "
                    "Use experimental structure or template-based homology model.",
        },
    }
    return params.get(strategy, {})


# --- Synthetic data generation for demonstration ---

def generate_synthetic_plddt_profile(
    n_residues: int = 300,
    seed: int = 42
) -> List[ResidueConfidence]:
    """Generate a realistic synthetic pLDDT profile for demonstration."""
    rng = np.random.RandomState(seed)

    plddts = np.zeros(n_residues)
    # Core structured regions (high confidence)
    for start, end in [(20, 80), (100, 180), (200, 270)]:
        plddts[start:end] = rng.normal(92, 4, end - start)
    # Loop regions (moderate confidence)
    for start, end in [(80, 100), (180, 200)]:
        plddts[start:end] = rng.normal(72, 8, end - start)
    # Terminal regions (low confidence)
    plddts[:20] = rng.normal(45, 12, 20)
    plddts[270:] = rng.normal(40, 15, n_residues - 270)

    plddts = np.clip(plddts, 0, 100)

    residues = []
    amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                   "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                   "THR", "TRP", "TYR", "VAL"]
    for i in range(n_residues):
        residues.append(ResidueConfidence(
            residue_id=i + 1,
            residue_name=rng.choice(amino_acids),
            chain_id="A",
            plddt=float(plddts[i])
        ))
    return residues
