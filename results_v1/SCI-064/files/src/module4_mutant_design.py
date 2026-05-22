"""
Module 4: Mutant Library Computational Design (Binding Affinity Tuning)
=======================================================================
Rational design of mutation libraries to tune biosensor properties:
- Saturation mutagenesis scoring
- Rosetta-like energy function for stability prediction
- Binding affinity perturbation (ΔΔG) estimation
- Multi-objective optimization (affinity vs. stability vs. selectivity)
- Directed evolution simulation
"""

import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import pearsonr
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import json
import os
import itertools


# -- Amino acid properties --
AA_PROPERTIES = {
    'A': {'mw': 89.1, 'hydro': 1.8, 'charge': 0, 'volume': 88.6, 'flexibility': 0.36},
    'R': {'mw': 174.2, 'hydro': -4.5, 'charge': 1, 'volume': 173.4, 'flexibility': 0.53},
    'N': {'mw': 132.1, 'hydro': -3.5, 'charge': 0, 'volume': 114.1, 'flexibility': 0.46},
    'D': {'mw': 133.1, 'hydro': -3.5, 'charge': -1, 'volume': 111.1, 'flexibility': 0.51},
    'C': {'mw': 121.2, 'hydro': 2.5, 'charge': 0, 'volume': 108.5, 'flexibility': 0.35},
    'Q': {'mw': 146.2, 'hydro': -3.5, 'charge': 0, 'volume': 143.8, 'flexibility': 0.49},
    'E': {'mw': 147.1, 'hydro': -3.5, 'charge': -1, 'volume': 138.4, 'flexibility': 0.50},
    'G': {'mw': 75.0, 'hydro': -0.4, 'charge': 0, 'volume': 60.1, 'flexibility': 0.54},
    'H': {'mw': 155.2, 'hydro': -3.2, 'charge': 0.5, 'volume': 153.2, 'flexibility': 0.32},
    'I': {'mw': 131.2, 'hydro': 4.5, 'charge': 0, 'volume': 166.7, 'flexibility': 0.46},
    'L': {'mw': 131.2, 'hydro': 3.8, 'charge': 0, 'volume': 166.7, 'flexibility': 0.40},
    'K': {'mw': 146.2, 'hydro': -3.9, 'charge': 1, 'volume': 168.6, 'flexibility': 0.53},
    'M': {'mw': 149.2, 'hydro': 1.9, 'charge': 0, 'volume': 162.9, 'flexibility': 0.42},
    'F': {'mw': 165.2, 'hydro': 2.8, 'charge': 0, 'volume': 189.9, 'flexibility': 0.31},
    'P': {'mw': 115.1, 'hydro': -1.6, 'charge': 0, 'volume': 112.7, 'flexibility': 0.51},
    'S': {'mw': 105.1, 'hydro': -0.8, 'charge': 0, 'volume': 89.0, 'flexibility': 0.51},
    'T': {'mw': 119.1, 'hydro': -0.7, 'charge': 0, 'volume': 116.1, 'flexibility': 0.44},
    'W': {'mw': 204.2, 'hydro': -0.9, 'charge': 0, 'volume': 227.8, 'flexibility': 0.31},
    'Y': {'mw': 181.2, 'hydro': -1.3, 'charge': 0, 'volume': 193.6, 'flexibility': 0.42},
    'V': {'mw': 117.1, 'hydro': 4.2, 'charge': 0, 'volume': 140.0, 'flexibility': 0.39},
}

# Blosum62-derived substitution energies (simplified)
SUBSTITUTION_ENERGY = {}
aa_list = list(AA_PROPERTIES.keys())
np.random.seed(42)
for a1 in aa_list:
    for a2 in aa_list:
        if a1 == a2:
            SUBSTITUTION_ENERGY[(a1, a2)] = 0.0
        else:
            # Based on property differences
            dh = abs(AA_PROPERTIES[a1]['hydro'] - AA_PROPERTIES[a2]['hydro'])
            dv = abs(AA_PROPERTIES[a1]['volume'] - AA_PROPERTIES[a2]['volume']) / 50.0
            dc = abs(AA_PROPERTIES[a1]['charge'] - AA_PROPERTIES[a2]['charge'])
            SUBSTITUTION_ENERGY[(a1, a2)] = 0.5 * dh + 0.3 * dv + 2.0 * dc


@dataclass
class Mutation:
    """Represents a point mutation."""
    position: int
    wild_type: str
    mutant: str
    ddG_stability: float  # kcal/mol
    ddG_binding: float    # kcal/mol
    effect_on_Kd: float   # fold change
    effect_on_hill: float # change in Hill coefficient
    confidence: float     # prediction confidence


@dataclass
class MutantDesign:
    """A designed mutant with multiple mutations."""
    mutations: List[Mutation]
    total_ddG_stability: float
    total_ddG_binding: float
    predicted_Kd: float
    predicted_hill: float
    predicted_dynamic_range: float
    fitness_score: float


def compute_ddG_stability(position: int, wt_aa: str, mut_aa: str,
                           burial_fraction: float = 0.5,
                           secondary_structure: str = "helix") -> float:
    """
    Estimate stability change (ΔΔG) for a point mutation.
    Simplified Rosetta-like energy function.
    """
    # Base substitution energy
    base_energy = SUBSTITUTION_ENERGY.get((wt_aa, mut_aa), 2.0)
    
    # Burial effect (buried mutations have larger impact)
    burial_factor = 1.0 + burial_fraction * 1.5
    
    # Secondary structure context
    ss_factors = {"helix": 1.0, "sheet": 0.8, "coil": 0.6}
    ss_factor = ss_factors.get(secondary_structure, 0.7)
    
    # Proline in helix is destabilizing
    if mut_aa == 'P' and secondary_structure == "helix":
        base_energy += 3.0
    
    # Glycine in sheet is destabilizing
    if mut_aa == 'G' and secondary_structure == "sheet":
        base_energy += 1.5
    
    # Charge introduction/removal in buried region
    if burial_fraction > 0.7:
        charge_change = abs(AA_PROPERTIES[mut_aa]['charge'] - AA_PROPERTIES[wt_aa]['charge'])
        base_energy += charge_change * 2.5
    
    return base_energy * burial_factor * ss_factor


def compute_ddG_binding(position: int, wt_aa: str, mut_aa: str,
                         is_binding_site: bool = False,
                         is_metal_coordinating: bool = False,
                         distance_to_ligand: float = 10.0) -> float:
    """
    Estimate binding affinity change (ΔΔG_binding) for a mutation.
    """
    if not is_binding_site and distance_to_ligand > 8.0:
        # Indirect effect only
        return SUBSTITUTION_ENERGY.get((wt_aa, mut_aa), 0.5) * 0.1
    
    base_effect = SUBSTITUTION_ENERGY.get((wt_aa, mut_aa), 1.0)
    
    # Distance-dependent decay
    distance_factor = np.exp(-distance_to_ligand / 5.0)
    
    # Metal coordination is critical
    if is_metal_coordinating:
        if wt_aa == 'C' and mut_aa != 'C':
            return 5.0 + base_effect  # Loss of Cys coordination
        if wt_aa == 'H' and mut_aa != 'H':
            return 3.0 + base_effect  # Loss of His coordination
        if wt_aa in ['D', 'E'] and mut_aa not in ['D', 'E']:
            return 2.5 + base_effect
    
    # Hydrophobic pocket mutations
    if is_binding_site:
        vol_change = abs(AA_PROPERTIES[mut_aa]['volume'] - AA_PROPERTIES[wt_aa]['volume'])
        base_effect += vol_change / 30.0  # Cavity/clash penalty
    
    return base_effect * distance_factor * (2.0 if is_binding_site else 0.5)


def design_saturation_mutagenesis(tf_type: str,
                                    n_residues: int,
                                    binding_residues: List[int],
                                    metal_residues: List[int],
                                    wt_sequence: str) -> List[Mutation]:
    """
    Design comprehensive saturation mutagenesis library.
    Score all possible single mutations at key positions.
    """
    np.random.seed(hash(tf_type) % 2**31)
    
    mutations = []
    target_positions = set()
    
    # Focus positions: binding site ± 2 residues
    for r in binding_residues:
        for offset in range(-2, 3):
            pos = r + offset
            if 0 <= pos < n_residues:
                target_positions.add(pos)
    
    # Also include allosteric pathway residues (linker region)
    linker_start = n_residues // 3
    linker_end = 2 * n_residues // 3
    for pos in range(linker_start, min(linker_end, linker_start + 15)):
        target_positions.add(pos)
    
    for pos in sorted(target_positions):
        wt_aa = wt_sequence[pos % len(wt_sequence)]
        
        is_binding = pos in binding_residues
        is_metal = pos in metal_residues
        burial = 0.7 if is_binding else 0.4
        dist_ligand = 3.0 if is_binding else 8.0 + np.random.uniform(0, 5)
        
        ss = "helix" if pos % 7 < 4 else ("sheet" if pos % 7 < 6 else "coil")
        
        for mut_aa in AA_PROPERTIES.keys():
            if mut_aa == wt_aa:
                continue
            
            ddG_stab = compute_ddG_stability(pos, wt_aa, mut_aa, burial, ss)
            ddG_bind = compute_ddG_binding(pos, wt_aa, mut_aa,
                                            is_binding, is_metal, dist_ligand)
            
            # Effect on Kd
            kd_fold = np.exp(ddG_bind / 0.593)  # RT at 298K
            
            # Effect on Hill coefficient
            hill_change = 0.0
            if is_binding:
                hill_change = -0.1 * ddG_bind / 2.0
            elif linker_start <= pos <= linker_end:
                hill_change = -0.05 * abs(ddG_stab) / 2.0
            
            confidence = max(0.1, 1.0 - 0.1 * max(abs(ddG_stab), abs(ddG_bind)))
            
            mutations.append(Mutation(
                position=pos,
                wild_type=wt_aa,
                mutant=mut_aa,
                ddG_stability=round(ddG_stab, 3),
                ddG_binding=round(ddG_bind, 3),
                effect_on_Kd=round(kd_fold, 4),
                effect_on_hill=round(hill_change, 4),
                confidence=round(min(1.0, max(0.0, confidence)), 3)
            ))
    
    return mutations


def optimize_mutant_combination(mutations: List[Mutation],
                                  target_Kd: float,
                                  wt_Kd: float,
                                  wt_hill: float,
                                  max_mutations: int = 3,
                                  stability_threshold: float = 3.0) -> List[MutantDesign]:
    """
    Find optimal mutation combinations using multi-objective optimization.
    Objectives: minimize |Kd - target_Kd|, maximize stability, maximize Hill.
    """
    # Filter destabilizing mutations
    viable = [m for m in mutations if m.ddG_stability < stability_threshold]
    
    # Group by position
    by_position = {}
    for m in viable:
        if m.position not in by_position:
            by_position[m.position] = []
        by_position[m.position].append(m)
    
    # Score top mutations per position
    top_per_position = {}
    for pos, muts in by_position.items():
        scored = sorted(muts, key=lambda m: abs(np.log10(max(m.effect_on_Kd, 1e-10)) -
                                                  np.log10(target_Kd / wt_Kd)))
        top_per_position[pos] = scored[:3]
    
    # Generate combinations
    positions = sorted(top_per_position.keys())[:10]
    designs = []
    
    for n_mut in range(1, min(max_mutations + 1, len(positions) + 1)):
        for pos_combo in itertools.combinations(positions, n_mut):
            # Try best mutation at each position
            combo_mutations = []
            for pos in pos_combo:
                combo_mutations.append(top_per_position[pos][0])
            
            total_ddG_stab = sum(m.ddG_stability for m in combo_mutations)
            total_ddG_bind = sum(m.ddG_binding for m in combo_mutations)
            
            if total_ddG_stab > stability_threshold * n_mut:
                continue
            
            pred_Kd = wt_Kd * np.prod([m.effect_on_Kd for m in combo_mutations])
            pred_hill = wt_hill + sum(m.effect_on_hill for m in combo_mutations)
            pred_hill = max(0.5, min(4.0, pred_hill))
            
            # Dynamic range estimation
            pred_dr = max(1.0, 100.0 * np.exp(-total_ddG_stab / 3.0))
            
            # Fitness: weighted multi-objective
            kd_score = -abs(np.log10(max(pred_Kd, 1e-10)) - np.log10(target_Kd))
            stab_score = -total_ddG_stab / 5.0
            hill_score = pred_hill / 3.0
            dr_score = np.log10(max(pred_dr, 1.0)) / 2.0
            
            fitness = 0.4 * kd_score + 0.3 * stab_score + 0.15 * hill_score + 0.15 * dr_score
            
            designs.append(MutantDesign(
                mutations=combo_mutations,
                total_ddG_stability=round(total_ddG_stab, 3),
                total_ddG_binding=round(total_ddG_bind, 3),
                predicted_Kd=round(pred_Kd, 6),
                predicted_hill=round(pred_hill, 3),
                predicted_dynamic_range=round(pred_dr, 1),
                fitness_score=round(fitness, 4)
            ))
    
    designs.sort(key=lambda d: d.fitness_score, reverse=True)
    return designs[:20]


def run_mutant_design(output_dir: str = "results") -> Dict:
    """Run mutant library design pipeline."""
    os.makedirs(output_dir, exist_ok=True)
    
    tf_configs = {
        "MerR": {"n_residues": 144, "binding": [82, 117, 126],
                 "metal": [82, 117, 126], "wt_Kd": 0.1, "wt_hill": 1.8,
                 "targets": {"lower_Kd": 0.01, "higher_Kd": 1.0, "broader_range": 0.1}},
        "ArsR": {"n_residues": 117, "binding": [32, 34, 37],
                 "metal": [32, 34, 37], "wt_Kd": 1.0, "wt_hill": 1.5,
                 "targets": {"lower_Kd": 0.1, "higher_Kd": 10.0, "broader_range": 1.0}},
        "CadC": {"n_residues": 122, "binding": [7, 11, 58, 60],
                 "metal": [7, 11, 58, 60], "wt_Kd": 0.5, "wt_hill": 2.0,
                 "targets": {"lower_Kd": 0.05, "higher_Kd": 5.0, "broader_range": 0.5}},
    }
    
    aa_seq = "MAKVLQFHDIGKQYLSETLKAVCDNFPDIPYRGYTREEIQAELVFPGSYKVLHKIDEQFER" * 3
    
    all_results = {}
    
    for tf, config in tf_configs.items():
        print(f"  Designing mutants for {tf}...")
        
        mutations = design_saturation_mutagenesis(
            tf, config["n_residues"],
            config["binding"], config["metal"],
            aa_seq[:config["n_residues"]]
        )
        
        designs_per_target = {}
        for target_name, target_Kd in config["targets"].items():
            designs = optimize_mutant_combination(
                mutations, target_Kd,
                config["wt_Kd"], config["wt_hill"]
            )
            
            designs_per_target[target_name] = [{
                "mutations": [f"{m.wild_type}{m.position}{m.mutant}" for m in d.mutations],
                "ddG_stability": d.total_ddG_stability,
                "ddG_binding": d.total_ddG_binding,
                "predicted_Kd": d.predicted_Kd,
                "predicted_hill": d.predicted_hill,
                "dynamic_range": d.predicted_dynamic_range,
                "fitness": d.fitness_score,
            } for d in designs[:5]]
        
        # Statistics
        beneficial = [m for m in mutations if m.ddG_binding < 0]
        neutral = [m for m in mutations if -0.5 <= m.ddG_binding <= 0.5]
        deleterious = [m for m in mutations if m.ddG_binding > 2.0]
        
        all_results[tf] = {
            "total_mutations_scored": len(mutations),
            "beneficial_count": len(beneficial),
            "neutral_count": len(neutral),
            "deleterious_count": len(deleterious),
            "positions_analyzed": len(set(m.position for m in mutations)),
            "designs": designs_per_target,
            "mutation_landscape": {
                "mean_ddG_stability": round(np.mean([m.ddG_stability for m in mutations]), 3),
                "mean_ddG_binding": round(np.mean([m.ddG_binding for m in mutations]), 3),
                "std_ddG_stability": round(np.std([m.ddG_stability for m in mutations]), 3),
                "std_ddG_binding": round(np.std([m.ddG_binding for m in mutations]), 3),
            }
        }
    
    with open(os.path.join(output_dir, "mutant_design.json"), 'w') as f:
        json.dump(all_results, f, indent=2)
    
    return all_results


if __name__ == "__main__":
    results = run_mutant_design()
    for tf, data in results.items():
        print(f"\n=== {tf} ===")
        print(f"  Total mutations scored: {data['total_mutations_scored']}")
        print(f"  Beneficial: {data['beneficial_count']}")
        for target, designs in data['designs'].items():
            if designs:
                print(f"  Best for {target}: {designs[0]['mutations']} "
                      f"(Kd={designs[0]['predicted_Kd']:.4f})")
