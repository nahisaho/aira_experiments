from __future__ import annotations

import json
import random
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"
LOGS_DIR = ROOT / "logs"
PROCESS_LOG = LOGS_DIR / "process-log.jsonl"
DATASET_CSV = DATA_DIR / "synthetic_affinity_dataset.csv"

for _d in (FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR):
    _d.mkdir(parents=True, exist_ok=True)


def set_global_seed(seed: int = 42) -> None:
    np.random.seed(seed)
    random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def to_serializable(obj: Any) -> Any:
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {str(k): to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_serializable(v) for v in obj]
    return str(obj)


def save_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(to_serializable(payload), handle, indent=2)


def append_log(
    phase: str,
    event_type: str,
    skill_or_tool: str,
    handoff_in: Dict[str, Any] | None = None,
    handoff_out: Dict[str, Any] | None = None,
    files_written: Iterable[str] | None = None,
    status: str = "ok",
) -> None:
    record = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": to_serializable(handoff_in or {}),
        "handoff_out": to_serializable(handoff_out or {}),
        "files_written": [str(Path(f)) for f in (files_written or [])],
        "status": status,
    }
    with PROCESS_LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")


PERIODIC_TABLE = Chem.GetPeriodicTable()
HYBRIDIZATION_MAP = {
    Chem.HybridizationType.SP: 1,
    Chem.HybridizationType.SP2: 2,
    Chem.HybridizationType.SP3: 3,
    Chem.HybridizationType.SP3D: 4,
    Chem.HybridizationType.SP3D2: 5,
}
CHIRALITY_MAP = {
    Chem.rdchem.ChiralType.CHI_UNSPECIFIED: 0,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: 2,
    Chem.rdchem.ChiralType.CHI_OTHER: 3,
}


def approximate_group(atomic_num: int) -> int:
    groups = {1: 1, 5: 13, 6: 14, 7: 15, 8: 16, 9: 17, 15: 15, 16: 16, 17: 17, 35: 17, 53: 17}
    return groups.get(atomic_num, min(18, max(1, atomic_num % 18 or 18)))


def ring_size_for_atom(atom: Chem.Atom) -> int:
    ring_info = atom.GetOwningMol().GetRingInfo()
    for size in range(3, 9):
        if ring_info.IsAtomInRingOfSize(atom.GetIdx(), size):
            return size
    return 0


def get_base_smiles() -> List[str]:
    return [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
        "CC(=O)NC1=CC=C(O)C=C1",
        "CC(C)NCC(O)COc1ccccc1OCC=C",
        "CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl",
        "CCOC(=O)C1=CC=CC=C1NC(=O)N(CC)CC",
        "COC1=CC=C(C=C1)CCN",
        "CC1=CC(=O)NC(=O)N1",
        "CC(C)(C)NCC(O)COc1cccc2ccccc12",
        "CN1CCC(CC1)C2=CN=CC=C2",
        "CC(C)NCC(O)COc1ccc(cc1)CCO",
        "CC(C)(C)OC(=O)N1CCC(CC1)CN2C=NC3=C2N=CN=C3N",
        "CN(C)CCOC(C1=CC=CC=C1)C2=CC=CC=C2",
        "CCOC(=O)N1CCN(CC1)C2=NC3=CC=CC=C3N=C2N",
        "COC1=CC2=C(C=C1)N=C(NC(=O)OC(C)(C)C)N=C2N",
        "CC1=C(C(=O)NC(=O)N1)N",
        "CCOC1=CC=C(C=C1)NC(=O)C2=CC=C(OCC)C=C2",
        "CCN1C=NC2=C1N=C(NC3=CC=CC=C3)N=C2N",
        "CC(C)OC1=CC=C(C=C1)C(C)C(=O)O",
        "CCOC(=O)C1CN(C(=O)O1)C2=CC=CC=C2",
        "COC1=CC=CC=C1OCCNCC(O)CO",
        "CCN(CC)CCOC(=O)C1=CC=CC=C1",
        "CC(C)(C)OC(=O)NCC1=CC=CC=C1",
        "CCOC(=O)C1=CN=CC=C1",
        "CC(C)OC1=CC=C(C=C1)NCC(O)CO",
        "COC1=CC=C(C=C1)OC(=O)NCC2=CC=CC=C2",
        "CCN1N=NN=C1NCC2=CC=CC=C2",
        "CC(C)C1=CC=C(C=C1)C(C)C(=O)NCCO",
        "COC1=CC=C(C=C1)C2=NC3=CC=CC=C3N2",
    ]


def sanitize_molecule(mol: Chem.Mol | None) -> Chem.Mol | None:
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True)


def generate_diverse_molecules(target_n: int = 200, seed: int = 42) -> List[Chem.Mol]:
    set_global_seed(seed)
    rng = random.Random(seed)
    base_mols = [sanitize_molecule(Chem.MolFromSmiles(s)) for s in get_base_smiles()]
    base_mols = [m for m in base_mols if m is not None]
    molecules: Dict[str, Chem.Mol] = {canonical_smiles(m): m for m in base_mols}
    fragment_smiles: List[str] = []
    for mol in base_mols:
        try:
            fragment_smiles.extend(list(BRICS.BRICSDecompose(mol, keepNonLeafNodes=True)))
        except Exception:
            continue
    fragment_smiles = [s for s in fragment_smiles if Chem.MolFromSmiles(s) is not None]
    rng.shuffle(fragment_smiles)

    attempts = 0
    while len(molecules) < target_n and attempts < 6000:
        attempts += 1
        chosen = rng.sample(fragment_smiles, k=min(3, len(fragment_smiles)))
        try:
            prods = BRICS.BRICSBuild([Chem.MolFromSmiles(s) for s in chosen if Chem.MolFromSmiles(s) is not None])
            for prod in prods:
                mol = sanitize_molecule(prod)
                if mol is None:
                    continue
                mw = Descriptors.MolWt(mol)
                if not 120 <= mw <= 650:
                    continue
                molecules.setdefault(canonical_smiles(mol), mol)
                if len(molecules) >= target_n:
                    break
        except Exception:
            continue

    atom_choices = [6, 7, 8, 9, 17]
    while len(molecules) < target_n:
        parent = Chem.Mol(rng.choice(list(molecules.values())))
        rw = Chem.RWMol(parent)
        atom = rw.GetAtomWithIdx(rng.randrange(rw.GetNumAtoms()))
        if atom.GetAtomicNum() > 1:
            atom.SetAtomicNum(rng.choice(atom_choices))
        mol = sanitize_molecule(rw.GetMol())
        if mol is not None:
            molecules.setdefault(canonical_smiles(mol), mol)
    return list(molecules.values())[:target_n]


def descriptor_bundle(mol: Chem.Mol) -> Dict[str, float]:
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": float(Lipinski.NumHDonors(mol)),
        "hba": float(Lipinski.NumHAcceptors(mol)),
        "qed": float(QED.qed(mol)),
        "rings": float(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_rings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "fsp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
    }


def synthetic_pic50_from_mol(mol: Chem.Mol) -> float:
    desc = descriptor_bundle(mol)
    smi = canonical_smiles(mol)
    seed = sum(ord(c) for c in smi) % 100000
    rng = np.random.default_rng(seed)
    value = (
        7.4
        - 0.0065 * abs(desc["mw"] - 360)
        + 0.45 * (4.2 - abs(desc["logp"] - 3.0))
        - 0.011 * abs(desc["tpsa"] - 90)
        - 0.18 * desc["hbd"]
        - 0.12 * desc["hba"]
        + 0.22 * desc["aromatic_rings"]
        + 0.35 * desc["qed"]
        + 0.15 * desc["fsp3"]
        + rng.normal(0, 0.18)
    )
    return float(np.clip(value, 4.0, 10.0))


def load_or_create_affinity_dataset(target_n: int = 200, seed: int = 42) -> pd.DataFrame:
    set_global_seed(seed)
    if DATASET_CSV.exists():
        return pd.read_csv(DATASET_CSV)
    mols = generate_diverse_molecules(target_n=target_n, seed=seed)
    rows = []
    for mol in mols:
        desc = descriptor_bundle(mol)
        rows.append({"smiles": canonical_smiles(mol), "pIC50": synthetic_pic50_from_mol(mol), **desc})
    df = pd.DataFrame(rows).drop_duplicates(subset="smiles").reset_index(drop=True)
    df.to_csv(DATASET_CSV, index=False)
    return df


def ensure_3d_molecule(mol: Chem.Mol, seed: int = 42) -> Chem.Mol:
    mol = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        params.useRandomCoords = True
        AllChem.EmbedMolecule(mol, params)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
    return mol
