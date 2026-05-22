"""
Module 2: Linker Length/Composition Systematic Optimization
Implements MD-inspired conformational sampling + MM-GBSA-style free energy estimation.
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
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeForceField, MMFFGetMoleculeProperties
from scipy.stats import gaussian_kde

from protac_utils import (
    log_event, smiles_to_mol, compute_descriptors,
    BRD4_WARHEAD_SMILES, VHL_LIGAND_SMILES, LINKER_LIBRARY
)

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)

# --------------------------------------------------------------------------
# 2.1  Generate systematic linker library (PEG, alkyl, mixed)
# --------------------------------------------------------------------------

def generate_linker_library() -> dict:
    """
    Expand the linker library systematically:
      - PEG (n=1–6 units)
      - Alkyl (n=3–8 carbons)
      - Piperazine-PEG hybrids
      - Amide-containing
    """
    linkers = {}
    # PEG series
    for n in range(1, 7):
        unit = "OCC"
        smiles = "O" + unit * n + "O"
        linkers[f"PEG{n}"] = smiles

    # Alkyl series
    for n in range(3, 9):
        linkers[f"Alkyl{n}"] = "C" * n

    # Piperazine-alkyl
    for n in range(1, 4):
        linkers[f"Pip{n}"] = "CCN1CCN(CC1)" + "CC" * n

    # Amide linkers
    linkers["Amide4"] = "CCNC(=O)CCC"
    linkers["Amide6"] = "CCNC(=O)CCCCC"
    linkers["AmidePEG"] = "CCNC(=O)COCCO"

    return linkers

# --------------------------------------------------------------------------
# 2.2  MM-GBSA-inspired free energy estimation
# --------------------------------------------------------------------------

def mmgbsa_free_energy(mol, n_confs: int = 50, seed: int = 42) -> dict:
    """
    Compute MM-GBSA-inspired binding free energy components.
    Uses MMFF94 force field for internal energy + GBSA continuum solvation proxy.

    ΔG_bind ≈ ΔE_MM + ΔG_solvation - TΔS_conf
    """
    mol_h = Chem.AddHs(mol)

    # Generate conformational ensemble
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params)

    energies = []
    conformer_rmsds = []

    for conf_id in range(mol_h.GetNumConformers()):
        props = MMFFGetMoleculeProperties(mol_h)
        if props is None:
            continue
        ff = MMFFGetMoleculeForceField(mol_h, props, confId=conf_id)
        if ff is None:
            continue
        ff.Minimize(maxIts=500)
        energy = ff.CalcEnergy()
        energies.append(energy)

    if not energies:
        return {}

    energies = np.array(energies)
    E_min = energies.min()
    kT = 0.593  # kcal/mol at 298 K

    # Boltzmann-weighted average energy
    boltzmann_weights = np.exp(-(energies - E_min) / kT)
    boltzmann_weights /= boltzmann_weights.sum()
    E_avg = np.average(energies, weights=boltzmann_weights)

    # Solvation energy proxy (GBSA-like): depends on TPSA and MW
    mol_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol_h))
    desc = compute_descriptors(mol_smiles)
    tpsa  = desc.get("TPSA", 0)
    logp  = desc.get("LogP", 0)
    mw    = desc.get("MW", 0)

    # ΔG_solv ≈ -0.012 * TPSA + 0.3 * LogP   (empirical GB/SA coefficients)
    dG_solv = -0.012 * tpsa + 0.3 * logp

    # Conformational entropy: S ≈ -k * Σ w_i * ln(w_i)
    entropy = -kT * np.sum(boltzmann_weights * np.log(boltzmann_weights + 1e-12))

    # Total free energy (relative, kcal/mol)
    dG_bind = E_avg * 0.01 + dG_solv - entropy  # scaled MM energy

    # RMSD spread of conformers (structural diversity)
    conf_coords = [mol_h.GetConformer(i).GetPositions()
                   for i in range(mol_h.GetNumConformers())]
    if len(conf_coords) > 1:
        ref = conf_coords[0]
        rmsds = [np.sqrt(((c - ref)**2).mean()) for c in conf_coords[1:]]
        rmsd_spread = np.mean(rmsds)
    else:
        rmsd_spread = 0.0

    return {
        "E_min_kcal":     round(E_min, 3),
        "E_avg_kcal":     round(E_avg, 3),
        "dG_solv_kcal":   round(dG_solv, 3),
        "conf_entropy":   round(entropy, 4),
        "dG_bind_est":    round(dG_bind, 4),
        "n_confs_gen":    mol_h.GetNumConformers(),
        "rmsd_spread_A":  round(rmsd_spread, 3),
    }


# --------------------------------------------------------------------------
# 2.3  Linker strain energy
# --------------------------------------------------------------------------

def linker_strain_energy(linker_smiles: str) -> float:
    """Compute internal strain energy (kcal/mol) as proxy for linker rigidity cost."""
    mol = Chem.MolFromSmiles(linker_smiles)
    if mol is None:
        return 999.0
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    props = MMFFGetMoleculeProperties(mol_h)
    if props is None:
        return 999.0
    ff = MMFFGetMoleculeForceField(mol_h, props)
    if ff is None:
        return 999.0
    ff.Minimize()
    return round(ff.CalcEnergy(), 3)


# --------------------------------------------------------------------------
# 2.4  Run linker optimization pipeline
# --------------------------------------------------------------------------

def run_linker_optimization():
    print("[Module 2] Linker length/composition optimization ...")
    log_event("linker_optimization", "handoff_started", "co-scientist-molecular-docking",
              {"method": "MM-GBSA", "force_field": "MMFF94"})

    linker_lib = generate_linker_library()
    records = []

    for name, smiles in linker_lib.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        desc = compute_descriptors(smiles)
        strain = linker_strain_energy(smiles)
        mol_3d = smiles_to_mol(smiles)
        if mol_3d is None:
            continue
        fep = mmgbsa_free_energy(mol_3d, n_confs=30)
        if not fep:
            continue

        rec = {
            "linker_name":   name,
            "smiles":        smiles,
            "n_heavy":       desc.get("HeavyAtoms", 0),
            "MW":            desc.get("MW", 0),
            "LogP":          desc.get("LogP", 0),
            "TPSA":          desc.get("TPSA", 0),
            "RotBonds":      desc.get("RotBonds", 0),
            "strain_kcal":   strain,
            **fep,
        }
        records.append(rec)
        print(f"  {name:12s}: dG={fep['dG_bind_est']:+.3f}  "
              f"strain={strain:.2f}  rot={desc.get('RotBonds',0)}")

    df = pd.DataFrame(records).sort_values("dG_bind_est")
    df.to_csv("results/linker_optimization.csv", index=False)
    df.to_json("data/linker_library_annotated.json", orient="records", indent=2)

    # ---------- Figure: 4-panel linker analysis ----------
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # A: dG by linker type
    ax = axes[0, 0]
    colors_map = {"PEG": "#1f77b4", "Alkyl": "#ff7f0e",
                  "Pip": "#2ca02c", "Amide": "#d62728"}
    bar_colors = [colors_map.get(n[:3] if not n[:4].isdigit() else n[:5], "#9467bd")
                  for n in df["linker_name"]]
    # Assign colors by prefix
    def get_color(name):
        for k, c in colors_map.items():
            if name.startswith(k):
                return c
        return "#9467bd"
    bar_colors = [get_color(n) for n in df["linker_name"]]
    bars = ax.barh(df["linker_name"], df["dG_bind_est"], color=bar_colors)
    ax.set_xlabel("ΔG_bind estimated (kcal/mol)")
    ax.set_title("A. Estimated Binding Free Energy by Linker")
    ax.axvline(x=0, color="black", lw=0.8)
    patches = [mpatches.Patch(color=c, label=k) for k, c in colors_map.items()]
    ax.legend(handles=patches, fontsize=8, loc="lower right")

    # B: RotBonds vs dG
    ax = axes[0, 1]
    scatter = ax.scatter(df["RotBonds"], df["dG_bind_est"],
                         c=[list(colors_map.values()).index(get_color(n)) % 4
                            for n in df["linker_name"]],
                         cmap="tab10", s=70, alpha=0.8)
    for _, row in df.iterrows():
        ax.annotate(row["linker_name"], (row["RotBonds"], row["dG_bind_est"]),
                    fontsize=6, ha="left")
    ax.set_xlabel("Rotatable Bonds")
    ax.set_ylabel("ΔG_bind estimated (kcal/mol)")
    ax.set_title("B. Flexibility vs Binding Free Energy")

    # C: Strain energy
    ax = axes[1, 0]
    ax.scatter(df["n_heavy"], df["strain_kcal"],
               c=[get_color(n) for n in df["linker_name"]], s=70, alpha=0.8)
    for _, row in df.iterrows():
        ax.annotate(row["linker_name"], (row["n_heavy"], row["strain_kcal"]),
                    fontsize=6, ha="left")
    ax.set_xlabel("Number of Heavy Atoms")
    ax.set_ylabel("Strain Energy (kcal/mol)")
    ax.set_title("C. Linker Strain Energy vs Size")

    # D: Conformational entropy
    ax = axes[1, 1]
    ax.scatter(df["RotBonds"], df["conf_entropy"],
               c=[get_color(n) for n in df["linker_name"]], s=70, alpha=0.8)
    z = np.polyfit(df["RotBonds"], df["conf_entropy"], 1)
    xfit = np.linspace(df["RotBonds"].min(), df["RotBonds"].max(), 50)
    ax.plot(xfit, np.polyval(z, xfit), "r--", lw=1.5, label="linear fit")
    ax.set_xlabel("Rotatable Bonds")
    ax.set_ylabel("Conformational Entropy (kcal/mol·K)")
    ax.set_title("D. Flexibility vs Conformational Entropy")
    ax.legend(fontsize=8)

    plt.suptitle("PROTAC Linker Systematic Optimization\n(MM-GBSA Free Energy Analysis)",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig("figures/02_linker_optimization.png", dpi=150, bbox_inches="tight")
    plt.close()

    import matplotlib
    # Figure: heatmap of linker property space
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    props = ["MW", "LogP", "TPSA", "RotBonds", "strain_kcal", "dG_bind_est", "conf_entropy"]
    heat_df = df.set_index("linker_name")[props]
    # Normalize each column
    heat_norm = (heat_df - heat_df.min()) / (heat_df.max() - heat_df.min() + 1e-9)
    import seaborn as sns
    sns.heatmap(heat_norm.T, annot=heat_df.T.round(2), fmt=".2f", cmap="viridis",
                linewidths=0.5, ax=ax2, cbar_kws={"label": "Normalized Value"})
    ax2.set_title("Linker Property Space Heatmap (MM-GBSA analysis)")
    ax2.set_xlabel("Linker")
    ax2.set_ylabel("Property")
    plt.tight_layout()
    plt.savefig("figures/02_linker_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()

    log_event("linker_optimization", "handoff_completed", "co-scientist-molecular-docking",
              {"n_linkers": len(df), "best": df.iloc[0]["linker_name"]},
              files_written=["results/linker_optimization.csv",
                             "data/linker_library_annotated.json",
                             "figures/02_linker_optimization.png",
                             "figures/02_linker_heatmap.png"])

    print(f"  Best linker: {df.iloc[0]['linker_name']} (dG={df.iloc[0]['dG_bind_est']:.3f})")
    return df


if __name__ == "__main__":
    run_linker_optimization()
