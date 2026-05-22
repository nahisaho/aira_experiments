"""
Perturb-seq Analysis Framework — Environment Setup & Data Simulation
====================================================================
Generates synthetic Perturb-seq data for pipeline development and testing.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
import os
import json
from datetime import datetime

SEED = 42
np.random.seed(SEED)

# --- Configuration ---
N_CELLS = 5000
N_GENES = 2000
N_GUIDES = 20          # 20 perturbation targets
N_PROGRAMS = 5         # latent gene programs
COMBO_FRAC = 0.10      # 10% cells get combinatorial perturbations
NT_FRAC = 0.15          # 15% non-targeting controls
MOI = 0.3               # multiplicity of infection (guide detection noise)

OUT_DIR = "data"
os.makedirs(OUT_DIR, exist_ok=True)


def simulate_perturbseq():
    """Create synthetic Perturb-seq AnnData with realistic structure."""

    # 1. Gene programs: each program activates a subset of genes
    program_matrix = np.zeros((N_PROGRAMS, N_GENES))
    genes_per_program = N_GENES // (N_PROGRAMS * 2)
    for p in range(N_PROGRAMS):
        start = p * genes_per_program
        program_matrix[p, start:start + genes_per_program] = np.random.uniform(0.5, 2.0, genes_per_program)

    # 2. Assign perturbations
    guide_names = [f"gene_{i}_guide" for i in range(N_GUIDES)]
    guide_names_with_nt = guide_names + ["non-targeting"]

    assignments = []
    for i in range(N_CELLS):
        if np.random.rand() < NT_FRAC:
            assignments.append("non-targeting")
        elif np.random.rand() < COMBO_FRAC:
            g1, g2 = np.random.choice(guide_names, 2, replace=False)
            assignments.append(f"{g1}|{g2}")
        else:
            assignments.append(np.random.choice(guide_names))

    # 3. Build expression matrix
    # Baseline expression (log-normal)
    baseline = np.random.lognormal(mean=1.0, sigma=0.8, size=(N_CELLS, N_GENES))

    # Perturbation effects: each guide affects 1-2 gene programs
    guide_to_program = {}
    for i, g in enumerate(guide_names):
        affected = np.random.choice(N_PROGRAMS, size=np.random.randint(1, 3), replace=False)
        guide_to_program[g] = affected

    for cell_idx, assign in enumerate(assignments):
        if assign == "non-targeting":
            continue
        guides = assign.split("|")
        effect = np.zeros(N_GENES)
        for g in guides:
            for prog in guide_to_program[g]:
                direction = np.random.choice([-1, 1])
                effect += direction * program_matrix[prog] * np.random.uniform(0.3, 1.0)
        baseline[cell_idx] *= np.exp(effect * 0.5)

    # 4. Add technical noise & sparsity
    baseline = np.random.poisson(baseline).astype(np.float32)
    baseline[baseline < 0] = 0

    # 5. Guide UMI counts (detection quality)
    guide_umi = np.zeros((N_CELLS, len(guide_names_with_nt)))
    for cell_idx, assign in enumerate(assignments):
        guides = assign.split("|")
        for g in guides:
            gidx = guide_names_with_nt.index(g)
            detected = np.random.binomial(1, 1 - MOI)
            if detected:
                guide_umi[cell_idx, gidx] = np.random.poisson(50)
            else:
                guide_umi[cell_idx, gidx] = np.random.poisson(2)
        # Add ambient guide noise
        noise_idx = np.random.choice(len(guide_names_with_nt), size=2, replace=False)
        for ni in noise_idx:
            guide_umi[cell_idx, ni] += np.random.poisson(1)

    # 6. Construct AnnData
    gene_names = [f"Gene_{i}" for i in range(N_GENES)]
    cell_ids = [f"Cell_{i}" for i in range(N_CELLS)]

    adata = ad.AnnData(
        X=sparse.csr_matrix(baseline),
        obs=pd.DataFrame({
            "cell_id": cell_ids,
            "perturbation": assignments,
            "n_guides": [len(a.split("|")) for a in assignments],
            "is_control": [a == "non-targeting" for a in assignments],
        }, index=cell_ids),
        var=pd.DataFrame({"gene_name": gene_names}, index=gene_names),
    )

    # Store guide UMI as obsm
    adata.obsm["guide_umi"] = guide_umi
    adata.uns["guide_names"] = guide_names_with_nt
    adata.uns["guide_to_program"] = {k: v.tolist() for k, v in guide_to_program.items()}
    adata.uns["program_matrix"] = program_matrix
    adata.uns["seed"] = SEED

    # 7. Save
    adata.write_h5ad(os.path.join(OUT_DIR, "perturbseq_simulated.h5ad"))

    # Summary
    summary = {
        "n_cells": N_CELLS,
        "n_genes": N_GENES,
        "n_guides": N_GUIDES,
        "n_programs": N_PROGRAMS,
        "combo_fraction": COMBO_FRAC,
        "nt_fraction": NT_FRAC,
        "n_control_cells": sum(1 for a in assignments if a == "non-targeting"),
        "n_single_perturb": sum(1 for a in assignments if a != "non-targeting" and "|" not in a),
        "n_combo_perturb": sum(1 for a in assignments if "|" in a),
        "timestamp": datetime.now().isoformat(),
    }
    with open(os.path.join(OUT_DIR, "simulation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✓ Simulated Perturb-seq data saved: {adata.shape}")
    print(f"  Controls: {summary['n_control_cells']}, Single: {summary['n_single_perturb']}, Combo: {summary['n_combo_perturb']}")
    return adata


if __name__ == "__main__":
    simulate_perturbseq()
