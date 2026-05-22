"""
Module 2: Allosteric Communication Pathway Analysis via Molecular Dynamics
==========================================================================
Analyzes allosteric signal transduction from ligand binding domain (LBD)
to DNA-binding domain (DBD) using:
- Normal mode analysis (elastic network model)
- Dynamic cross-correlation maps
- Perturbation response scanning
- Allosteric pathway identification
"""

import numpy as np
from scipy.linalg import eigh
from scipy.spatial.distance import cdist, squareform, pdist
from scipy.signal import correlate
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import json
import os


@dataclass
class AlloPathway:
    """Represents an allosteric communication pathway."""
    residues: List[int]
    correlation_strength: float
    pathway_length: float
    bottleneck_residue: int
    signal_transfer_efficiency: float


def build_elastic_network(ca_coords: np.ndarray,
                           cutoff: float = 12.0,
                           spring_constant: float = 1.0) -> np.ndarray:
    """
    Build Anisotropic Network Model (ANM) Hessian matrix.
    """
    n_atoms = len(ca_coords)
    n_dof = 3 * n_atoms
    hessian = np.zeros((n_dof, n_dof))
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            diff = ca_coords[j] - ca_coords[i]
            dist = np.linalg.norm(diff)
            
            if dist < cutoff:
                # Distance-dependent spring constant (Tirion)
                gamma = spring_constant * (cutoff / dist) ** 2
                
                # 3x3 super-element
                outer = np.outer(diff, diff) / (dist ** 2) * gamma
                
                # Fill Hessian
                for a in range(3):
                    for b in range(3):
                        hessian[3*i+a, 3*j+b] = -outer[a, b]
                        hessian[3*j+b, 3*i+a] = -outer[a, b]
                        hessian[3*i+a, 3*i+b] += outer[a, b]
                        hessian[3*j+a, 3*j+b] += outer[a, b]
    
    return hessian


def compute_normal_modes(hessian: np.ndarray, n_modes: int = 20) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute normal modes from the Hessian matrix.
    Returns eigenvalues and eigenvectors (skipping first 6 trivial modes).
    """
    eigenvalues, eigenvectors = eigh(hessian)
    
    # Skip first 6 trivial modes (translation + rotation)
    nontrivial_start = 6
    selected_eigenvalues = eigenvalues[nontrivial_start:nontrivial_start + n_modes]
    selected_eigenvectors = eigenvectors[:, nontrivial_start:nontrivial_start + n_modes]
    
    # Replace any near-zero eigenvalues
    selected_eigenvalues = np.maximum(selected_eigenvalues, 1e-10)
    
    return selected_eigenvalues, selected_eigenvectors


def compute_cross_correlation(eigenvectors: np.ndarray,
                                eigenvalues: np.ndarray,
                                n_atoms: int) -> np.ndarray:
    """
    Compute dynamic cross-correlation matrix (DCCM) from normal modes.
    C_ij = <Δr_i · Δr_j> / sqrt(<Δr_i²><Δr_j²>)
    """
    n_modes = len(eigenvalues)
    
    # Compute covariance matrix
    cov = np.zeros((n_atoms, n_atoms))
    
    for k in range(n_modes):
        mode = eigenvectors[:, k].reshape(n_atoms, 3)
        inv_eigenval = 1.0 / eigenvalues[k]
        
        for i in range(n_atoms):
            for j in range(i, n_atoms):
                dot_product = np.dot(mode[i], mode[j]) * inv_eigenval
                cov[i, j] += dot_product
                if i != j:
                    cov[j, i] += dot_product
    
    # Normalize to correlation
    diag = np.sqrt(np.diag(cov))
    diag[diag == 0] = 1e-10
    correlation = cov / np.outer(diag, diag)
    
    return np.clip(correlation, -1.0, 1.0)


def perturbation_response_scanning(hessian: np.ndarray,
                                     n_atoms: int,
                                     perturbed_residues: List[int]) -> np.ndarray:
    """
    Perturbation Response Scanning (PRS) to identify allosteric coupling.
    Apply unit force at perturbed residues and measure response everywhere.
    """
    n_dof = 3 * n_atoms
    
    # Pseudo-inverse of Hessian (compliance matrix)
    eigenvalues, eigenvectors = eigh(hessian)
    
    # Regularize
    eigenvalues = np.maximum(eigenvalues, 1e-8)
    inv_eigenvalues = 1.0 / eigenvalues
    inv_eigenvalues[:6] = 0  # zero out trivial modes
    
    compliance = eigenvectors @ np.diag(inv_eigenvalues) @ eigenvectors.T
    
    # Response matrix
    response = np.zeros((len(perturbed_residues), n_atoms))
    
    for idx, res in enumerate(perturbed_residues):
        # Apply random unit forces in all 3 directions
        total_response = np.zeros(n_atoms)
        for direction in range(3):
            force = np.zeros(n_dof)
            force[3 * res + direction] = 1.0
            
            displacement = compliance @ force
            
            # Compute magnitude of response per residue
            for j in range(n_atoms):
                total_response[j] += np.linalg.norm(displacement[3*j:3*j+3]) ** 2
        
        response[idx] = np.sqrt(total_response / 3.0)
    
    return response


def find_allosteric_pathways(correlation_matrix: np.ndarray,
                              source_residues: List[int],
                              target_residues: List[int],
                              threshold: float = 0.3,
                              max_pathways: int = 5) -> List[AlloPathway]:
    """
    Find allosteric communication pathways using correlation-weighted
    shortest path algorithm.
    """
    n = correlation_matrix.shape[0]
    
    # Convert correlation to distance
    abs_corr = np.abs(correlation_matrix)
    # Distance = -log(|correlation|), with threshold
    distance_matrix = np.full((n, n), np.inf)
    for i in range(n):
        for j in range(n):
            if abs_corr[i, j] > threshold and i != j:
                distance_matrix[i, j] = -np.log(abs_corr[i, j])
    
    # Dijkstra's algorithm for each source-target pair
    pathways = []
    
    for src in source_residues:
        for tgt in target_residues:
            # Dijkstra
            dist = np.full(n, np.inf)
            dist[src] = 0
            visited = np.zeros(n, dtype=bool)
            prev = np.full(n, -1, dtype=int)
            
            for _ in range(n):
                # Find unvisited node with minimum distance
                unvisited_dist = np.where(visited, np.inf, dist)
                u = np.argmin(unvisited_dist)
                
                if dist[u] == np.inf:
                    break
                if u == tgt:
                    break
                
                visited[u] = True
                
                for v in range(n):
                    if not visited[v] and distance_matrix[u, v] < np.inf:
                        alt = dist[u] + distance_matrix[u, v]
                        if alt < dist[v]:
                            dist[v] = alt
                            prev[v] = u
            
            if dist[tgt] < np.inf:
                # Reconstruct path
                path = []
                current = tgt
                while current != -1:
                    path.append(current)
                    current = prev[current]
                path.reverse()
                
                # Calculate pathway properties
                path_corrs = []
                for k in range(len(path) - 1):
                    path_corrs.append(abs_corr[path[k], path[k+1]])
                
                avg_corr = np.mean(path_corrs) if path_corrs else 0.0
                
                # Find bottleneck (weakest link)
                bottleneck_idx = np.argmin(path_corrs) if path_corrs else 0
                bottleneck = path[bottleneck_idx + 1] if path_corrs else path[0]
                
                # Signal transfer efficiency
                efficiency = np.prod(path_corrs) if path_corrs else 0.0
                
                pathways.append(AlloPathway(
                    residues=path,
                    correlation_strength=avg_corr,
                    pathway_length=dist[tgt],
                    bottleneck_residue=bottleneck,
                    signal_transfer_efficiency=efficiency
                ))
    
    # Sort by efficiency
    pathways.sort(key=lambda p: p.signal_transfer_efficiency, reverse=True)
    return pathways[:max_pathways]


def generate_md_trajectory(ca_coords: np.ndarray,
                            eigenvectors: np.ndarray,
                            eigenvalues: np.ndarray,
                            n_frames: int = 500,
                            temperature: float = 300.0) -> np.ndarray:
    """
    Generate pseudo-MD trajectory from normal mode analysis.
    Uses Boltzmann-weighted mode amplitudes.
    """
    n_atoms = len(ca_coords)
    n_modes = len(eigenvalues)
    kB = 0.001987  # kcal/(mol·K)
    
    trajectory = np.zeros((n_frames, n_atoms, 3))
    
    np.random.seed(42)
    
    for frame in range(n_frames):
        displacement = np.zeros((n_atoms, 3))
        
        for k in range(n_modes):
            # Boltzmann amplitude
            amplitude = np.sqrt(kB * temperature / eigenvalues[k])
            phase = 2 * np.pi * frame / (20 + k * 5)  # quasi-periodic
            noise = np.random.normal(0, 0.1)
            
            mode_shape = eigenvectors[:, k].reshape(n_atoms, 3)
            displacement += amplitude * np.sin(phase + noise) * mode_shape
        
        trajectory[frame] = ca_coords + displacement
    
    return trajectory


def compute_rmsf(trajectory: np.ndarray) -> np.ndarray:
    """Compute root mean square fluctuation per residue."""
    mean_coords = trajectory.mean(axis=0)
    deviations = trajectory - mean_coords[np.newaxis, :, :]
    msf = (deviations ** 2).sum(axis=2).mean(axis=0)
    return np.sqrt(msf)


def run_md_analysis(output_dir: str = "results") -> Dict:
    """Run complete MD analysis pipeline."""
    os.makedirs(output_dir, exist_ok=True)
    
    from src.module1_structural_analysis import generate_model_aTF_structure
    
    tf_types = ["MerR", "ArsR", "CadC", "CueR", "SmtB"]
    all_results = {}
    
    for tf in tf_types:
        print(f"  MD analysis for {tf}...")
        atoms = generate_model_aTF_structure(tf)
        
        # Extract CA coordinates
        ca_atoms = [a for a in atoms if a.name == "CA"]
        ca_coords = np.array([a.coords for a in ca_atoms])
        n_atoms = len(ca_coords)
        
        if n_atoms < 10:
            continue
        
        # Build ANM and compute modes
        hessian = build_elastic_network(ca_coords, cutoff=12.0)
        eigenvalues, eigenvectors = compute_normal_modes(hessian, n_modes=min(20, n_atoms - 7))
        
        # Dynamic cross-correlation
        dccm = compute_cross_correlation(eigenvectors, eigenvalues, n_atoms)
        
        # Generate trajectory and RMSF
        trajectory = generate_md_trajectory(ca_coords, eigenvectors, eigenvalues, n_frames=200)
        rmsf = compute_rmsf(trajectory)
        
        # Define source (MBD) and target (DBD) residues
        n_third = n_atoms // 3
        source_residues = list(range(2 * n_third, min(n_atoms - 1, 2 * n_third + 5)))
        target_residues = list(range(0, min(5, n_third)))
        
        # PRS analysis
        prs_response = perturbation_response_scanning(hessian, n_atoms, source_residues)
        
        # Find allosteric pathways
        pathways = find_allosteric_pathways(dccm, source_residues, target_residues)
        
        # Compile results
        all_results[tf] = {
            "n_residues": n_atoms,
            "n_modes_computed": len(eigenvalues),
            "eigenvalue_spectrum": eigenvalues[:10].tolist(),
            "collectivity": float(np.exp(-np.sum(
                (eigenvalues / eigenvalues.sum()) *
                np.log(eigenvalues / eigenvalues.sum() + 1e-20)
            ) / np.log(len(eigenvalues)))),
            "mean_rmsf_A": float(np.mean(rmsf)),
            "max_rmsf_A": float(np.max(rmsf)),
            "max_rmsf_residue": int(np.argmax(rmsf)),
            "dccm_stats": {
                "mean_abs_correlation": float(np.mean(np.abs(dccm))),
                "max_positive_correlation": float(np.max(dccm[~np.eye(n_atoms, dtype=bool)])),
                "max_negative_correlation": float(np.min(dccm)),
                "n_highly_correlated_pairs": int(np.sum(np.abs(dccm) > 0.5) // 2),
            },
            "prs_max_response_residue": int(np.argmax(prs_response.mean(axis=0))),
            "n_allosteric_pathways": len(pathways),
            "best_pathway": {
                "residues": pathways[0].residues if pathways else [],
                "length": len(pathways[0].residues) if pathways else 0,
                "efficiency": float(pathways[0].signal_transfer_efficiency) if pathways else 0,
                "bottleneck": int(pathways[0].bottleneck_residue) if pathways else -1,
                "avg_correlation": float(pathways[0].correlation_strength) if pathways else 0,
            } if pathways else None,
            "rmsf_profile": rmsf.tolist(),
        }
    
    # Save results
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)
    
    with open(os.path.join(output_dir, "md_analysis.json"), 'w') as f:
        json.dump(all_results, f, indent=2, cls=NumpyEncoder)
    
    return all_results


if __name__ == "__main__":
    results = run_md_analysis()
    for tf, data in results.items():
        print(f"\n=== {tf} ===")
        print(f"  Residues: {data['n_residues']}")
        print(f"  Mean RMSF: {data['mean_rmsf_A']:.2f} Å")
        if data.get('best_pathway'):
            print(f"  Best pathway efficiency: {data['best_pathway']['efficiency']:.4f}")
