"""
PROTAC Design Framework - Shared Utilities
Utility functions used across all modules.
"""
import json
import datetime
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
try:
    from rdkit.Chem import Draw
    _DRAW_AVAILABLE = True
except ImportError:
    _DRAW_AVAILABLE = False

LOG_PATH = "logs/process-log.jsonl"

os.makedirs("logs", exist_ok=True)

def log_event(phase: str, event_type: str, skill: str, details: dict,
              files_written: list = None, status: str = "ok"):
    entry = {
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "details": details,
        "files_written": files_written or [],
        "status": status,
    }
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry) + "\n")

def smiles_to_mol(smiles: str):
    """Parse SMILES and add 3D coordinates."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        # Fallback: random coordinates
        AllChem.EmbedMolecule(mol, AllChem.ETDG())
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol

def compute_descriptors(smiles: str) -> dict:
    """Compute RDKit physicochemical descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "MW":    Descriptors.MolWt(mol),
        "LogP":  Descriptors.MolLogP(mol),
        "HBD":   rdMolDescriptors.CalcNumHBD(mol),
        "HBA":   rdMolDescriptors.CalcNumHBA(mol),
        "TPSA":  Descriptors.TPSA(mol),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "ArRings":  rdMolDescriptors.CalcNumAromaticRings(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
        "Fsp3":  rdMolDescriptors.CalcFractionCSP3(mol),
    }

def morgan_fingerprint(smiles: str, radius: int = 2, nbits: int = 2048) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nbits)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    fp = gen.GetFingerprintAsNumPy(mol)
    return fp.astype(np.float32)

# --- PROTAC-specific SMILES fragments ---

# BRD4 warhead: JQ1-based bromodomain inhibitor
BRD4_WARHEAD_SMILES = "Cc1sc2c(c1-c1ccc(Cl)cc1)C(=O)N(C)c1ccc(cc1)C2"

# VHL E3 ligase ligand (VH032 derivative)
VHL_LIGAND_SMILES = "CC(C)(C)C(=O)N[C@@H]1CC[C@H](CC1)C(=O)N[C@@H](Cc1ccccc1)C(=O)O"

# CRBN ligand (thalidomide derivative)
CRBN_LIGAND_SMILES = "O=C1CN(C(=O)c2ccccc21)C1CCC(=O)NC1=O"

# IAP ligand (LCL-161 fragment)
IAP_LIGAND_SMILES = "CC(C)(C)c1nc2cc(F)ccc2c(=O)n1CC(N)C1CCCC1"

# Representative linkers (PEG / alkyl)
LINKER_LIBRARY = {
    "PEG2":    "OCCOCCO",
    "PEG3":    "OCCOCCOCCO",
    "PEG4":    "OCCOCCOCCOCCO",
    "Alkyl4":  "CCCCCC",
    "Alkyl6":  "CCCCCCCC",
    "PipeAm":  "N1CCN(CC1)CC",
    "Hybrid3": "OCCNCCNCO",
}
