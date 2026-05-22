"""
Module 1: Ternary Complex (POI–PROTAC–E3) Structural Modeling
Implements geometry scoring and interaction analysis for ternary complex formation.
Rosetta-inspired scoring terms are approximated with RDKit distance/angle geometry.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
try:
    from rdkit.Chem import Draw
    _DRAW_AVAILABLE = True
except ImportError:
    Draw = None
    _DRAW_AVAILABLE = False
from io import BytesIO
import json

from protac_utils import (
    log_event, smiles_to_mol, compute_descriptors,
    BRD4_WARHEAD_SMILES, VHL_LIGAND_SMILES, CRBN_LIGAND_SMILES,
    IAP_LIGAND_SMILES, LINKER_LIBRARY
)

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

# --------------------------------------------------------------------------
# 1.1  Geometric scoring model for ternary complex
# --------------------------------------------------------------------------

def compute_gyration_radius(mol) -> float:
    """Radius of gyration (Å) as proxy for molecular compactness."""
    conf = mol.GetConformer()
    coords = conf.GetPositions()
    center = coords.mean(axis=0)
    rg = np.sqrt(((coords - center) ** 2).sum(axis=1).mean())
    return rg

def linker_end_to_end(mol) -> float:
    """Approximate end-to-end distance by finding the two most distal atoms."""
    conf = mol.GetConformer()
    coords = conf.GetPositions()
    dists = np.linalg.norm(coords[:, None] - coords[None, :], axis=-1)
    return dists.max()

def ternary_complex_score(poi_smiles: str, e3_smiles: str,
                          linker_smiles: str, n_atoms_linker: int) -> dict:
    """
    Heuristic ternary complex formation score.

    Score components (Rosetta-inspired):
      - geometry_score   : linker length matching (0–1, higher=better)
      - flexibility_score: linker conformational flexibility (fewer rot bonds = higher)
      - cooperativity    : predicted cooperativity α (PPIs + geometry)
      - buried_surface   : approximate buried surface area proxy (Å²)
    """
    poi_mol  = smiles_to_mol(poi_smiles)
    e3_mol   = smiles_to_mol(e3_smiles)
    linker_mol = smiles_to_mol(linker_smiles)

    if not all([poi_mol, e3_mol, linker_mol]):
        return {}

    rg_poi  = compute_gyration_radius(poi_mol)
    rg_e3   = compute_gyration_radius(e3_mol)
    e2e     = linker_end_to_end(linker_mol)

    # Optimal linker bridges the two binding pockets (typically 12–20 Å apart)
    optimal_distance = (rg_poi + rg_e3) * 0.9  # heuristic
    geometry_score = np.exp(-((e2e - optimal_distance) ** 2) / (2 * 4.0 ** 2))

    rot_bonds_linker = rdMolDescriptors.CalcNumRotatableBonds(linker_mol)
    flexibility_score = 1.0 / (1.0 + 0.15 * rot_bonds_linker)

    # Cooperativity: empirically, medium-length linkers show highest cooperativity
    alpha_cooperativity = 2.5 * geometry_score * flexibility_score + 0.5

    # Buried surface area proxy (sum of accessible surface of components)
    poi_desc  = compute_descriptors(poi_smiles)
    e3_desc   = compute_descriptors(e3_smiles)
    bsa_proxy = (poi_desc.get("MW", 0) ** 0.67 + e3_desc.get("MW", 0) ** 0.67) * 2.5

    # Rosetta-like composite score (lower = better ternary complex)
    composite_score = -(geometry_score * 5 + flexibility_score * 3 +
                        alpha_cooperativity * 2 + bsa_proxy * 0.01)

    return {
        "linker_smiles":      linker_smiles,
        "n_atoms_linker":     n_atoms_linker,
        "rg_poi_A":           round(rg_poi, 2),
        "rg_e3_A":            round(rg_e3, 2),
        "e2e_linker_A":       round(e2e, 2),
        "optimal_dist_A":     round(optimal_distance, 2),
        "geometry_score":     round(geometry_score, 4),
        "flexibility_score":  round(flexibility_score, 4),
        "alpha_cooperativity":round(alpha_cooperativity, 3),
        "bsa_proxy_A2":       round(bsa_proxy, 1),
        "composite_score":    round(composite_score, 3),
    }


# --------------------------------------------------------------------------
# 1.2  Run ternary complex analysis for all linker lengths
# --------------------------------------------------------------------------

def run_ternary_modeling():
    print("[Module 1] Ternary complex structural modeling ...")
    log_event("ternary_complex", "handoff_started", "co-scientist-molecular-docking",
              {"poi": "BRD4", "e3": "VHL", "linkers": list(LINKER_LIBRARY.keys())})

    records = []
    for name, lsmi in LINKER_LIBRARY.items():
        n_atoms = Chem.MolFromSmiles(lsmi).GetNumHeavyAtoms()
        score = ternary_complex_score(BRD4_WARHEAD_SMILES, VHL_LIGAND_SMILES, lsmi, n_atoms)
        if score:
            score["linker_name"] = name
            records.append(score)
            print(f"  {name:12s}: geometry={score['geometry_score']:.3f}  "
                  f"alpha={score['alpha_cooperativity']:.3f}  "
                  f"composite={score['composite_score']:.3f}")

    df = pd.DataFrame(records)
    df = df.sort_values("composite_score")

    # Save CSV
    csv_path = "results/ternary_complex_scores.csv"
    df.to_csv(csv_path, index=False)

    # ---------- Figure 1: composite score by linker ----------
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(df)))

    # Panel A: geometry score
    ax = axes[0]
    ax.barh(df["linker_name"], df["geometry_score"], color=colors)
    ax.set_xlabel("Geometry Score")
    ax.set_title("A. Linker Geometry Matching")
    ax.set_xlim(0, 1)

    # Panel B: cooperativity
    ax = axes[1]
    ax.barh(df["linker_name"], df["alpha_cooperativity"], color=colors)
    ax.set_xlabel("Cooperativity α")
    ax.set_title("B. Predicted Cooperativity (α)")
    ax.axvline(x=1.0, color="red", ls="--", lw=1, label="α=1 (no cooperativity)")
    ax.legend(fontsize=8)

    # Panel C: end-to-end vs optimal
    ax = axes[2]
    ax.scatter(df["e2e_linker_A"], df["optimal_dist_A"], s=80,
               c=np.arange(len(df)), cmap="viridis")
    for _, row in df.iterrows():
        ax.annotate(row["linker_name"], (row["e2e_linker_A"], row["optimal_dist_A"]),
                    fontsize=7, ha="left", va="bottom")
    lim = max(df["e2e_linker_A"].max(), df["optimal_dist_A"].max()) + 2
    ax.plot([0, lim], [0, lim], "r--", lw=1, label="perfect match")
    ax.set_xlabel("Linker End-to-End Distance (Å)")
    ax.set_ylabel("Optimal POI–E3 Distance (Å)")
    ax.set_title("C. Linker Length vs Optimal Distance")
    ax.legend(fontsize=8)

    plt.suptitle("BRD4–PROTAC–VHL Ternary Complex Structural Scoring", fontsize=13, y=1.02)
    plt.tight_layout()
    fig_path = "figures/01_ternary_complex_scores.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    # ---------- Figure 2: 2D structure of best PROTAC scaffold ----------
    best_linker = df.iloc[0]["linker_name"]
    best_lsmi   = LINKER_LIBRARY[best_linker]
    # Rough concatenation for visualization (simplified PROTAC scaffold)
    full_protac_smiles = BRD4_WARHEAD_SMILES + "." + best_lsmi + "." + VHL_LIGAND_SMILES

    mols_to_draw = [
        Chem.MolFromSmiles(BRD4_WARHEAD_SMILES),
        Chem.MolFromSmiles(best_lsmi),
        Chem.MolFromSmiles(VHL_LIGAND_SMILES),
    ]
    labels = ["BRD4 Warhead (JQ1)", f"Best Linker ({best_linker})", "VHL Ligand (VH032)"]
    fig2_path = "figures/01_protac_fragments.png"
    if _DRAW_AVAILABLE:
        img = Draw.MolsToGridImage(mols_to_draw, molsPerRow=3,
                                   subImgSize=(400, 300), legends=labels)
        img.save(fig2_path)
    else:
        # Text fallback
        fig_fb, ax_fb = plt.subplots(figsize=(10, 3))
        ax_fb.axis("off")
        for idx, (lbl, mol) in enumerate(zip(labels, mols_to_draw)):
            smi_str = Chem.MolToSmiles(mol) if mol else "N/A"
            ax_fb.text(idx/3 + 0.05, 0.5,
                       f"{lbl}\n{smi_str[:40]}...",
                       transform=ax_fb.transAxes, fontsize=8, va="center")
        plt.title("PROTAC Fragment SMILES (2D rendering unavailable)")
        plt.tight_layout()
        plt.savefig(fig2_path, dpi=100, bbox_inches="tight")
        plt.close()

    log_event("ternary_complex", "handoff_completed", "co-scientist-molecular-docking",
              {"best_linker": best_linker, "best_score": float(df.iloc[0]["composite_score"])},
              files_written=[csv_path, fig_path, fig2_path])

    print(f"  Best linker: {best_linker} (composite={df.iloc[0]['composite_score']:.3f})")
    return df


if __name__ == "__main__":
    run_ternary_modeling()
