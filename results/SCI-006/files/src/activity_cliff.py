"""
Module 5: Activity Cliff Detection and Chemical Space Exploration

Implements methods for identifying activity cliffs (structurally similar
compounds with large activity differences) and strategies for exploring
chemical space around lead compounds.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from itertools import combinations


@dataclass
class Compound:
    """Compound with structural and activity information."""
    compound_id: str
    smiles: str
    pki: float  # -log10(Ki)
    fingerprint: Optional[np.ndarray] = None
    mw: float = 0.0
    logp: float = 0.0
    hbd: int = 0
    hba: int = 0
    tpsa: float = 0.0
    rotatable_bonds: int = 0
    scaffold: str = ""


@dataclass
class ActivityCliff:
    """An activity cliff pair."""
    compound_a: str
    compound_b: str
    similarity: float
    activity_diff: float  # |pKi_a - pKi_b|
    sali: float  # Structure-Activity Landscape Index
    cliff_type: str = ""  # "potency_gain", "potency_loss", "scaffold_hop"
    
    structural_difference: str = ""  # Human-readable description
    
    def __post_init__(self):
        if not self.cliff_type:
            if self.activity_diff > 2:
                self.cliff_type = "major_cliff"
            elif self.activity_diff > 1:
                self.cliff_type = "moderate_cliff"
            else:
                self.cliff_type = "minor_cliff"


@dataclass
class ChemicalSpaceAnalysis:
    """Analysis of chemical space coverage."""
    n_compounds: int
    n_clusters: int
    coverage_score: float  # 0-1
    diversity_score: float  # 0-1
    
    # Principal component analysis
    pca_variance_explained: List[float] = field(default_factory=list)
    pca_coordinates: Optional[np.ndarray] = None
    
    # Cluster information
    cluster_sizes: List[int] = field(default_factory=list)
    cluster_centroids: Optional[np.ndarray] = None
    
    # Exploration recommendations
    underexplored_regions: List[Dict] = field(default_factory=list)


def compute_tanimoto_similarity(fp_a: np.ndarray, fp_b: np.ndarray) -> float:
    """Compute Tanimoto similarity between two binary fingerprints."""
    intersection = np.sum(np.logical_and(fp_a, fp_b))
    union = np.sum(np.logical_or(fp_a, fp_b))
    if union == 0:
        return 0.0
    return float(intersection / union)


def compute_sali(similarity: float, activity_diff: float) -> float:
    """
    Compute Structure-Activity Landscape Index (SALI).
    SALI = |pKi_a - pKi_b| / (1 - Tanimoto)
    Higher SALI indicates a more dramatic activity cliff.
    """
    if similarity >= 1.0:
        return float('inf')
    return activity_diff / (1.0 - similarity)


def detect_activity_cliffs(
    compounds: List[Compound],
    similarity_threshold: float = 0.7,
    activity_threshold: float = 1.0,
    sali_threshold: float = 5.0
) -> List[ActivityCliff]:
    """
    Detect activity cliffs in a compound set.
    
    An activity cliff is defined as a pair of compounds with:
    - Tanimoto similarity >= similarity_threshold
    - |ΔpKi| >= activity_threshold
    - SALI >= sali_threshold
    """
    cliffs = []
    
    for i, j in combinations(range(len(compounds)), 2):
        comp_a = compounds[i]
        comp_b = compounds[j]
        
        if comp_a.fingerprint is None or comp_b.fingerprint is None:
            continue
        
        sim = compute_tanimoto_similarity(comp_a.fingerprint, comp_b.fingerprint)
        
        if sim < similarity_threshold:
            continue
        
        activity_diff = abs(comp_a.pki - comp_b.pki)
        
        if activity_diff < activity_threshold:
            continue
        
        sali = compute_sali(sim, activity_diff)
        
        if sali >= sali_threshold:
            cliffs.append(ActivityCliff(
                compound_a=comp_a.compound_id,
                compound_b=comp_b.compound_id,
                similarity=sim,
                activity_diff=activity_diff,
                sali=sali,
            ))
    
    cliffs.sort(key=lambda c: c.sali, reverse=True)
    return cliffs


def analyze_chemical_space(
    compounds: List[Compound],
    n_components: int = 2
) -> ChemicalSpaceAnalysis:
    """Analyze chemical space coverage and diversity."""
    fps = np.array([c.fingerprint for c in compounds if c.fingerprint is not None])
    
    if len(fps) == 0:
        return ChemicalSpaceAnalysis(n_compounds=0, n_clusters=0,
                                      coverage_score=0, diversity_score=0)
    
    # PCA
    fps_centered = fps - fps.mean(axis=0)
    try:
        U, S, Vt = np.linalg.svd(fps_centered, full_matrices=False)
        variance_explained = (S ** 2) / np.sum(S ** 2)
        pca_coords = U[:, :n_components] * S[:n_components]
    except np.linalg.LinAlgError:
        pca_coords = fps_centered[:, :n_components]
        variance_explained = np.ones(n_components) / n_components
    
    # Simple clustering (k-means style)
    n_clusters = min(10, len(fps) // 5 + 1)
    rng = np.random.RandomState(42)
    centroids = fps[rng.choice(len(fps), n_clusters, replace=False)]
    
    for _ in range(20):
        distances = np.array([
            [np.sum((fp - c) ** 2) for c in centroids]
            for fp in fps
        ])
        labels = distances.argmin(axis=1)
        for k in range(n_clusters):
            mask = labels == k
            if mask.any():
                centroids[k] = fps[mask].mean(axis=0)
    
    cluster_sizes = [int(np.sum(labels == k)) for k in range(n_clusters)]
    
    # Diversity: average pairwise distance
    n_sample = min(100, len(fps))
    sample_idx = rng.choice(len(fps), n_sample, replace=False)
    dists = []
    for i, j in combinations(sample_idx, 2):
        dists.append(1 - compute_tanimoto_similarity(fps[i], fps[j]))
    diversity = float(np.mean(dists)) if dists else 0.0
    
    # Coverage: fraction of space covered by clusters
    coverage = min(1.0, n_clusters / 20 * diversity)
    
    return ChemicalSpaceAnalysis(
        n_compounds=len(compounds),
        n_clusters=n_clusters,
        coverage_score=float(coverage),
        diversity_score=float(diversity),
        pca_variance_explained=variance_explained[:n_components].tolist(),
        pca_coordinates=pca_coords,
        cluster_sizes=cluster_sizes,
        underexplored_regions=_find_underexplored(pca_coords, labels, n_clusters),
    )


def _find_underexplored(coords, labels, n_clusters):
    """Identify underexplored regions in chemical space."""
    regions = []
    grid_size = 10
    x_range = (coords[:, 0].min(), coords[:, 0].max())
    y_range = (coords[:, 1].min(), coords[:, 1].max())
    
    x_bins = np.linspace(x_range[0], x_range[1], grid_size + 1)
    y_bins = np.linspace(y_range[0], y_range[1], grid_size + 1)
    
    for i in range(grid_size):
        for j in range(grid_size):
            mask = (
                (coords[:, 0] >= x_bins[i]) & (coords[:, 0] < x_bins[i+1]) &
                (coords[:, 1] >= y_bins[j]) & (coords[:, 1] < y_bins[j+1])
            )
            count = mask.sum()
            if count == 0:
                regions.append({
                    "center_x": float((x_bins[i] + x_bins[i+1]) / 2),
                    "center_y": float((y_bins[j] + y_bins[j+1]) / 2),
                    "density": 0,
                })
    
    return regions[:5]


def generate_exploration_strategy(
    cliffs: List[ActivityCliff],
    space_analysis: ChemicalSpaceAnalysis,
    compounds: List[Compound]
) -> Dict:
    """Generate chemical space exploration strategy based on cliff analysis."""
    strategies = []
    
    # Strategy 1: Cliff investigation
    if cliffs:
        top_cliffs = cliffs[:5]
        strategies.append({
            "name": "Activity Cliff Investigation",
            "priority": "high",
            "description": "Synthesize analogs that interpolate between cliff pairs",
            "targets": [
                {
                    "pair": (c.compound_a, c.compound_b),
                    "sali": c.sali,
                    "expected_insight": "SAR discontinuity analysis"
                }
                for c in top_cliffs
            ],
        })
    
    # Strategy 2: Underexplored regions
    if space_analysis.underexplored_regions:
        strategies.append({
            "name": "Chemical Space Expansion",
            "priority": "medium",
            "description": "Design compounds targeting underexplored regions",
            "n_target_regions": len(space_analysis.underexplored_regions),
        })
    
    # Strategy 3: Scaffold hopping
    strategies.append({
        "name": "Scaffold Hopping",
        "priority": "medium",
        "description": "Explore alternative scaffolds maintaining key pharmacophores",
        "approach": "bioisostere_replacement",
    })
    
    return {
        "n_strategies": len(strategies),
        "strategies": strategies,
        "estimated_compounds_to_synthesize": sum(
            len(s.get("targets", [])) * 3 for s in strategies
        ) + 20,
        "chemical_space_coverage_current": space_analysis.coverage_score,
        "chemical_space_coverage_target": min(1.0, space_analysis.coverage_score + 0.2),
    }


def generate_synthetic_compounds(
    n_compounds: int = 100,
    seed: int = 42
) -> List[Compound]:
    """Generate synthetic compound dataset for demonstration."""
    rng = np.random.RandomState(seed)
    fp_size = 1024
    
    compounds = []
    # Generate several scaffold classes
    n_scaffolds = 5
    scaffold_fps = [rng.randint(0, 2, fp_size) for _ in range(n_scaffolds)]
    scaffold_base_pki = rng.uniform(5, 8, n_scaffolds)
    
    smiles_templates = [
        "c1ccc(NC(=O)c2ccccc2)cc1",
        "c1ccc(-c2cnc3ccccc3n2)cc1",
        "O=C(Nc1ccccc1)c1ccc(O)cc1",
        "c1ccc(CSc2nnc(-c3ccccc3)o2)cc1",
        "c1ccc(-n2c(=O)c3ccccc3nc2=O)cc1",
    ]
    
    for i in range(n_compounds):
        scaffold_idx = i % n_scaffolds
        
        # Fingerprint: scaffold + random variation
        fp = scaffold_fps[scaffold_idx].copy()
        n_mutations = rng.randint(10, 50)
        mutation_pos = rng.choice(fp_size, n_mutations, replace=False)
        fp[mutation_pos] = 1 - fp[mutation_pos]
        
        # Activity: base + noise, with occasional cliffs
        pki = scaffold_base_pki[scaffold_idx] + rng.normal(0, 0.5)
        if rng.random() < 0.05:  # 5% chance of activity cliff
            pki += rng.choice([-2.5, 2.5])
        pki = np.clip(pki, 3, 11)
        
        compounds.append(Compound(
            compound_id=f"CMPD_{i+1:04d}",
            smiles=smiles_templates[scaffold_idx],
            pki=float(pki),
            fingerprint=fp.astype(float),
            mw=float(rng.uniform(250, 600)),
            logp=float(rng.uniform(0, 5)),
            hbd=int(rng.randint(0, 5)),
            hba=int(rng.randint(2, 10)),
            tpsa=float(rng.uniform(40, 140)),
            rotatable_bonds=int(rng.randint(1, 10)),
            scaffold=f"Scaffold_{scaffold_idx + 1}",
        ))
    
    return compounds
