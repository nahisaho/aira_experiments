"""
Module 1: Ligand Binding Pocket Structural Analysis & Docking
=============================================================
Analyzes allosteric transcription factor binding pockets using
structural bioinformatics approaches including:
- Binding pocket identification (fpocket-like geometric analysis)
- Molecular property calculation for ligands
- Scoring function for docking pose evaluation
- Binding energy estimation
"""

import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.distance import cdist
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
import json
import os


@dataclass
class Atom:
    """Represents an atom with 3D coordinates and properties."""
    name: str
    element: str
    x: float
    y: float
    z: float
    residue: str = ""
    chain: str = "A"
    bfactor: float = 0.0
    charge: float = 0.0
    radius: float = 1.7  # van der Waals radius

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


@dataclass
class BindingPocket:
    """Represents an identified binding pocket."""
    center: np.ndarray
    volume: float
    surface_area: float
    residues: List[str]
    druggability_score: float
    hydrophobicity: float
    atoms: List[Atom] = field(default_factory=list)


@dataclass
class Ligand:
    """Represents a ligand molecule for docking."""
    name: str
    smiles: str
    molecular_weight: float
    logP: float
    hbd: int  # H-bond donors
    hba: int  # H-bond acceptors
    rotatable_bonds: int
    atoms: List[Atom] = field(default_factory=list)
    binding_energy: float = 0.0


# -- VdW radii for common elements --
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    'Fe': 1.80, 'Zn': 1.39, 'Cu': 1.40, 'Mn': 1.61, 'Mg': 1.73,
    'Ca': 1.97, 'Hg': 1.55, 'Pb': 2.02, 'Cd': 1.58, 'As': 1.85,
    'Cr': 1.66, 'Ni': 1.63, 'Co': 1.64,
}

# -- Hydrophobicity scale (Kyte-Doolittle) --
HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2,
}


def generate_model_aTF_structure(tf_type: str = "MerR") -> List[Atom]:
    """
    Generate a model allosteric transcription factor structure.
    Supports: MerR (mercury), ArsR (arsenic), CadC (cadmium),
    CueR (copper), SmtB (zinc), NahR (naphthalene), DmpR (phenol).
    """
    np.random.seed(42)
    
    tf_configs = {
        "MerR": {"n_residues": 144, "binding_metals": ["Hg"],
                 "key_cysteines": [82, 117, 126], "dbd_range": (1, 80),
                 "linker_range": (81, 95), "mbd_range": (96, 144)},
        "ArsR": {"n_residues": 117, "binding_metals": ["As"],
                 "key_cysteines": [32, 34, 37], "dbd_range": (1, 70),
                 "linker_range": (71, 80), "mbd_range": (81, 117)},
        "CadC": {"n_residues": 122, "binding_metals": ["Cd", "Pb", "Zn"],
                 "key_cysteines": [7, 11, 58, 60], "dbd_range": (1, 65),
                 "linker_range": (66, 78), "mbd_range": (79, 122)},
        "CueR": {"n_residues": 135, "binding_metals": ["Cu"],
                 "key_cysteines": [112, 120], "dbd_range": (1, 75),
                 "linker_range": (76, 90), "mbd_range": (91, 135)},
        "SmtB": {"n_residues": 122, "binding_metals": ["Zn"],
                 "key_cysteines": [14, 61, 64, 97], "dbd_range": (1, 60),
                 "linker_range": (61, 75), "mbd_range": (76, 122)},
        "NahR": {"n_residues": 305, "binding_metals": [],
                 "key_cysteines": [], "dbd_range": (1, 120),
                 "linker_range": (121, 160), "mbd_range": (161, 305)},
        "DmpR": {"n_residues": 283, "binding_metals": [],
                 "key_cysteines": [], "dbd_range": (1, 100),
                 "linker_range": (101, 140), "mbd_range": (141, 283)},
    }
    
    config = tf_configs.get(tf_type, tf_configs["MerR"])
    residue_names = list(HYDROPHOBICITY.keys())
    atoms = []
    
    for i in range(1, config["n_residues"] + 1):
        if i in config.get("key_cysteines", []):
            res_name = "CYS"
        else:
            res_name = residue_names[np.random.randint(len(residue_names))]
        
        # Generate backbone with helical structure
        if config["dbd_range"][0] <= i <= config["dbd_range"][1]:
            # DNA-binding domain: helix-turn-helix
            theta = i * 100 * np.pi / 180
            r = 8.0
            z_offset = i * 1.5
        elif config["linker_range"][0] <= i <= config["linker_range"][1]:
            # Linker region: extended
            theta = i * 160 * np.pi / 180
            r = 12.0
            z_offset = config["dbd_range"][1] * 1.5 + (i - config["linker_range"][0]) * 3.0
        else:
            # Metal/ligand binding domain
            theta = i * 100 * np.pi / 180
            r = 10.0
            z_offset = (config["dbd_range"][1] * 1.5 +
                       (config["linker_range"][1] - config["linker_range"][0]) * 3.0 +
                       (i - config["mbd_range"][0]) * 1.5)
        
        # Backbone atoms (CA, N, C, O)
        for atom_name, offset in [("N", -0.5), ("CA", 0.0), ("C", 0.5), ("O", 1.0)]:
            element = atom_name[0]
            atoms.append(Atom(
                name=atom_name, element=element,
                x=r * np.cos(theta) + np.random.normal(0, 0.3) + offset * 0.3,
                y=r * np.sin(theta) + np.random.normal(0, 0.3),
                z=z_offset + offset * 0.5,
                residue=f"{res_name}{i}", chain="A",
                bfactor=20.0 + np.random.normal(0, 5),
                radius=VDW_RADII.get(element, 1.7)
            ))
        
        # Side chain CB
        atoms.append(Atom(
            name="CB", element="C",
            x=r * np.cos(theta) + 1.5 * np.cos(theta + np.pi/3),
            y=r * np.sin(theta) + 1.5 * np.sin(theta + np.pi/3),
            z=z_offset + 0.8,
            residue=f"{res_name}{i}", chain="A",
            radius=VDW_RADII.get("C", 1.7)
        ))
    
    return atoms


def identify_binding_pockets(atoms: List[Atom],
                              probe_radius: float = 1.4,
                              min_pocket_volume: float = 100.0,
                              grid_spacing: float = 1.0) -> List[BindingPocket]:
    """
    Identify binding pockets using alpha-sphere based method.
    """
    coords = np.array([a.coords for a in atoms])
    
    # Grid-based pocket detection
    margin = 5.0
    x_range = np.arange(coords[:, 0].min() - margin, coords[:, 0].max() + margin, grid_spacing)
    y_range = np.arange(coords[:, 1].min() - margin, coords[:, 1].max() + margin, grid_spacing)
    z_range = np.arange(coords[:, 2].min() - margin, coords[:, 2].max() + margin, grid_spacing)
    
    # Sample grid points for efficiency
    np.random.seed(123)
    n_samples = min(5000, len(x_range) * len(y_range) * len(z_range))
    grid_points = np.column_stack([
        np.random.uniform(coords[:, 0].min() - margin, coords[:, 0].max() + margin, n_samples),
        np.random.uniform(coords[:, 1].min() - margin, coords[:, 1].max() + margin, n_samples),
        np.random.uniform(coords[:, 2].min() - margin, coords[:, 2].max() + margin, n_samples),
    ])
    
    # Find cavity points
    distances = cdist(grid_points, coords)
    radii = np.array([a.radius for a in atoms])
    min_dist = distances.min(axis=1)
    
    # Pocket points: outside protein surface but within interaction distance
    cavity_mask = (min_dist > probe_radius + 1.0) & (min_dist < probe_radius + 5.0)
    cavity_points = grid_points[cavity_mask]
    
    if len(cavity_points) < 10:
        # Relax criteria
        cavity_mask = (min_dist > probe_radius) & (min_dist < probe_radius + 8.0)
        cavity_points = grid_points[cavity_mask]
    
    # Cluster cavity points into pockets
    pockets = []
    if len(cavity_points) > 0:
        from scipy.cluster.hierarchy import fcluster, linkage
        if len(cavity_points) > 3:
            Z = linkage(cavity_points, method='average')
            clusters = fcluster(Z, t=8.0, criterion='distance')
            
            for cluster_id in np.unique(clusters):
                cluster_mask = clusters == cluster_id
                cluster_points = cavity_points[cluster_mask]
                
                if len(cluster_points) < 5:
                    continue
                
                center = cluster_points.mean(axis=0)
                
                # Calculate pocket properties
                try:
                    hull = ConvexHull(cluster_points)
                    volume = hull.volume
                    surface_area = hull.area
                except:
                    volume = len(cluster_points) * grid_spacing**3
                    surface_area = volume ** (2/3) * 4.836
                
                if volume < min_pocket_volume:
                    continue
                
                # Find lining residues
                dist_to_center = cdist([center], coords)[0]
                nearby_mask = dist_to_center < 8.0
                nearby_atoms = [atoms[i] for i in range(len(atoms)) if nearby_mask[i]]
                residues = list(set(a.residue for a in nearby_atoms))
                
                # Calculate hydrophobicity
                hydro_scores = []
                for res in residues:
                    res_name = ''.join(c for c in res if c.isalpha())
                    if res_name in HYDROPHOBICITY:
                        hydro_scores.append(HYDROPHOBICITY[res_name])
                avg_hydro = np.mean(hydro_scores) if hydro_scores else 0.0
                
                # Druggability score (simplified Fpocket-like)
                druggability = min(1.0, max(0.0,
                    0.3 * min(1.0, volume / 500.0) +
                    0.25 * min(1.0, len(residues) / 20.0) +
                    0.25 * (1.0 / (1.0 + np.exp(-avg_hydro / 2.0))) +
                    0.2 * min(1.0, surface_area / 400.0)
                ))
                
                pockets.append(BindingPocket(
                    center=center,
                    volume=volume,
                    surface_area=surface_area,
                    residues=residues,
                    druggability_score=druggability,
                    hydrophobicity=avg_hydro,
                    atoms=nearby_atoms
                ))
    
    # Sort by druggability
    pockets.sort(key=lambda p: p.druggability_score, reverse=True)
    return pockets


def generate_environmental_ligands() -> Dict[str, List[Ligand]]:
    """Generate ligand libraries for environmental pollutant classes."""
    
    heavy_metals = [
        Ligand("Mercury(II)", "Hg", 200.59, -0.5, 0, 0, 0),
        Ligand("Arsenic(III)", "[As]", 74.92, -0.4, 0, 0, 0),
        Ligand("Cadmium(II)", "[Cd]", 112.41, -0.3, 0, 0, 0),
        Ligand("Lead(II)", "[Pb]", 207.2, -0.2, 0, 0, 0),
        Ligand("Chromium(VI)", "[Cr]", 52.0, -0.6, 0, 0, 0),
        Ligand("Copper(II)", "[Cu]", 63.55, -0.1, 0, 0, 0),
        Ligand("Zinc(II)", "[Zn]", 65.38, -0.3, 0, 0, 0),
        Ligand("Nickel(II)", "[Ni]", 58.69, -0.2, 0, 0, 0),
    ]
    
    organic_solvents = [
        Ligand("Benzene", "c1ccccc1", 78.11, 2.13, 0, 0, 0),
        Ligand("Toluene", "Cc1ccccc1", 92.14, 2.73, 0, 0, 0),
        Ligand("Naphthalene", "c1ccc2ccccc2c1", 128.17, 3.30, 0, 0, 0),
        Ligand("Phenol", "Oc1ccccc1", 94.11, 1.46, 1, 1, 0),
        Ligand("Xylene", "Cc1ccccc1C", 106.17, 3.12, 0, 0, 0),
        Ligand("Chloroform", "ClC(Cl)Cl", 119.38, 1.97, 0, 0, 0),
        Ligand("TCE", "ClC=C(Cl)Cl", 131.39, 2.42, 0, 0, 0),
        Ligand("PCB-1", "c1ccc(-c2ccccc2Cl)cc1", 188.65, 4.56, 0, 0, 0),
    ]
    
    pesticides = [
        Ligand("Atrazine", "CCNc1nc(Cl)nc(NC(C)C)n1", 215.68, 2.61, 2, 4, 3),
        Ligand("DDT", "ClC(Cl)=C(c1ccc(Cl)cc1)c1ccc(Cl)cc1", 354.49, 6.91, 0, 0, 2),
    ]
    
    return {
        "heavy_metals": heavy_metals,
        "organic_solvents": organic_solvents,
        "pesticides": pesticides,
    }


def score_docking_pose(pocket: BindingPocket, ligand: Ligand) -> Dict:
    """
    Score a docking pose using a simplified scoring function.
    Components: van der Waals, electrostatic, H-bond, desolvation.
    """
    np.random.seed(hash(ligand.name) % 2**31)
    
    # Shape complementarity score
    vol_ratio = min(ligand.molecular_weight / 5.0, pocket.volume) / max(pocket.volume, 1.0)
    shape_score = -2.0 * np.exp(-((vol_ratio - 0.4)**2) / 0.1)
    
    # Hydrophobic matching
    if ligand.logP > 2.0 and pocket.hydrophobicity > 0:
        hydro_score = -1.5 * min(1.0, ligand.logP / 5.0) * min(1.0, pocket.hydrophobicity / 3.0)
    elif ligand.logP < 0 and pocket.hydrophobicity < 0:
        hydro_score = -1.0
    else:
        hydro_score = 0.5  # mismatch penalty
    
    # H-bond score
    hbond_score = -0.5 * min(ligand.hbd + ligand.hba, 6)
    
    # Rotatable bond penalty (entropy)
    rot_penalty = 0.3 * ligand.rotatable_bonds
    
    # Metal coordination bonus (for metal ligands)
    metal_bonus = 0.0
    if ligand.molecular_weight > 50 and ligand.logP < 0:
        cys_count = sum(1 for r in pocket.residues if 'CYS' in r)
        his_count = sum(1 for r in pocket.residues if 'HIS' in r)
        asp_count = sum(1 for r in pocket.residues if 'ASP' in r or 'GLU' in r)
        metal_bonus = -2.0 * (cys_count * 0.8 + his_count * 0.5 + asp_count * 0.3)
    
    # Total binding energy (kcal/mol)
    total_energy = shape_score + hydro_score + hbond_score + rot_penalty + metal_bonus
    total_energy += np.random.normal(0, 0.5)  # noise
    
    # Dissociation constant estimation
    kd_uM = np.exp(total_energy / (0.593))  # RT at 298K ≈ 0.593 kcal/mol
    
    ligand.binding_energy = total_energy
    
    return {
        "ligand": ligand.name,
        "total_energy_kcal": round(total_energy, 2),
        "shape_score": round(shape_score, 3),
        "hydrophobic_score": round(hydro_score, 3),
        "hbond_score": round(hbond_score, 3),
        "rotation_penalty": round(rot_penalty, 3),
        "metal_coordination": round(metal_bonus, 3),
        "estimated_Kd_uM": round(kd_uM, 4),
        "pocket_volume_A3": round(pocket.volume, 1),
        "druggability": round(pocket.druggability_score, 3),
    }


def run_structural_analysis(output_dir: str = "results") -> Dict:
    """Run complete structural analysis pipeline for all TF types."""
    os.makedirs(output_dir, exist_ok=True)
    
    tf_types = ["MerR", "ArsR", "CadC", "CueR", "SmtB", "NahR", "DmpR"]
    ligand_library = generate_environmental_ligands()
    
    all_results = {}
    
    for tf in tf_types:
        print(f"  Analyzing {tf}...")
        atoms = generate_model_aTF_structure(tf)
        pockets = identify_binding_pockets(atoms)
        
        if not pockets:
            all_results[tf] = {"pockets_found": 0, "docking_results": []}
            continue
        
        # Dock all relevant ligands to the best pocket
        best_pocket = pockets[0]
        docking_results = []
        
        # Select appropriate ligands based on TF type
        if tf in ["MerR", "ArsR", "CadC", "CueR", "SmtB"]:
            ligands = ligand_library["heavy_metals"]
        else:
            ligands = ligand_library["organic_solvents"]
        
        for lig in ligands:
            result = score_docking_pose(best_pocket, lig)
            docking_results.append(result)
        
        docking_results.sort(key=lambda x: x["total_energy_kcal"])
        
        all_results[tf] = {
            "n_atoms": len(atoms),
            "pockets_found": len(pockets),
            "best_pocket": {
                "volume": round(best_pocket.volume, 1),
                "surface_area": round(best_pocket.surface_area, 1),
                "druggability": round(best_pocket.druggability_score, 3),
                "hydrophobicity": round(best_pocket.hydrophobicity, 2),
                "n_lining_residues": len(best_pocket.residues),
            },
            "docking_results": docking_results[:5],  # top 5
        }
    
    # Save results
    with open(os.path.join(output_dir, "structural_analysis.json"), 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    return all_results


if __name__ == "__main__":
    results = run_structural_analysis()
    for tf, data in results.items():
        print(f"\n=== {tf} ===")
        print(f"  Pockets: {data.get('pockets_found', 0)}")
        if data.get('docking_results'):
            best = data['docking_results'][0]
            print(f"  Best ligand: {best['ligand']} ({best['total_energy_kcal']} kcal/mol)")
