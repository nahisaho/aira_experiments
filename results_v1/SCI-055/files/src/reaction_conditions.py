"""Reaction condition prediction module.

Predicts optimal reaction conditions including:
- Solvent selection
- Temperature range
- Catalyst recommendation
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem import rdMolDescriptors


# Reaction condition knowledge base
REACTION_CONDITIONS_DB = {
    "Amide bond formation": {
        "solvents": [("DMF", 0.35), ("DCM", 0.30), ("THF", 0.20), ("DMSO", 0.15)],
        "temperature_range": (0, 25),
        "catalysts": [("EDC/HOBt", 0.40), ("HATU", 0.35), ("DCC", 0.25)],
        "time_hours": (2, 24),
        "typical_yield": (70, 90),
    },
    "Suzuki coupling": {
        "solvents": [("DMF/H2O", 0.30), ("Dioxane/H2O", 0.30), ("THF/H2O", 0.25), ("Toluene", 0.15)],
        "temperature_range": (80, 120),
        "catalysts": [("Pd(PPh3)4", 0.35), ("Pd(dppf)Cl2", 0.30), ("Pd(OAc)2/SPhos", 0.35)],
        "time_hours": (4, 18),
        "typical_yield": (60, 95),
        "base": [("K2CO3", 0.40), ("Cs2CO3", 0.30), ("K3PO4", 0.30)],
    },
    "Ester hydrolysis": {
        "solvents": [("THF/H2O", 0.35), ("MeOH/H2O", 0.30), ("Dioxane/H2O", 0.20), ("EtOH/H2O", 0.15)],
        "temperature_range": (20, 60),
        "catalysts": [("NaOH", 0.40), ("LiOH", 0.35), ("KOH", 0.25)],
        "time_hours": (1, 12),
        "typical_yield": (80, 98),
    },
    "Reductive amination": {
        "solvents": [("DCE", 0.30), ("MeOH", 0.30), ("DCM", 0.25), ("THF", 0.15)],
        "temperature_range": (20, 40),
        "catalysts": [("NaBH(OAc)3", 0.40), ("NaBH3CN", 0.35), ("NaBH4", 0.25)],
        "time_hours": (2, 18),
        "typical_yield": (60, 85),
    },
    "Williamson ether synthesis": {
        "solvents": [("DMF", 0.35), ("DMSO", 0.25), ("Acetone", 0.25), ("THF", 0.15)],
        "temperature_range": (50, 80),
        "catalysts": [("K2CO3", 0.40), ("NaH", 0.35), ("Cs2CO3", 0.25)],
        "time_hours": (4, 24),
        "typical_yield": (65, 90),
    },
    "Fischer esterification": {
        "solvents": [("Toluene", 0.35), ("Benzene", 0.25), ("neat", 0.25), ("DCM", 0.15)],
        "temperature_range": (60, 120),
        "catalysts": [("H2SO4", 0.40), ("p-TsOH", 0.35), ("DMAP/DCC", 0.25)],
        "time_hours": (6, 48),
        "typical_yield": (50, 85),
    },
    "Heck reaction": {
        "solvents": [("DMF", 0.35), ("DMA", 0.25), ("NMP", 0.25), ("Dioxane", 0.15)],
        "temperature_range": (100, 140),
        "catalysts": [("Pd(OAc)2/PPh3", 0.35), ("Pd(PPh3)4", 0.30), ("Pd(dba)2", 0.35)],
        "time_hours": (8, 24),
        "typical_yield": (55, 85),
        "base": [("Et3N", 0.40), ("K2CO3", 0.30), ("NaOAc", 0.30)],
    },
    "Buchwald-Hartwig amination": {
        "solvents": [("Toluene", 0.35), ("Dioxane", 0.30), ("THF", 0.20), ("DME", 0.15)],
        "temperature_range": (80, 110),
        "catalysts": [("Pd2(dba)3/BINAP", 0.35), ("Pd(OAc)2/XPhos", 0.35), ("Pd(dppf)Cl2", 0.30)],
        "time_hours": (6, 24),
        "typical_yield": (60, 90),
        "base": [("NaOtBu", 0.40), ("Cs2CO3", 0.35), ("K3PO4", 0.25)],
    },
    "Aldol condensation": {
        "solvents": [("THF", 0.35), ("Et2O", 0.25), ("DCM", 0.25), ("Toluene", 0.15)],
        "temperature_range": (-78, 0),
        "catalysts": [("LDA", 0.35), ("NaOH", 0.30), ("L-Proline", 0.35)],
        "time_hours": (1, 12),
        "typical_yield": (50, 80),
    },
    "Wittig reaction": {
        "solvents": [("THF", 0.40), ("DCM", 0.25), ("DMF", 0.20), ("Et2O", 0.15)],
        "temperature_range": (-78, 25),
        "catalysts": [("n-BuLi", 0.40), ("NaHMDS", 0.35), ("KHMDS", 0.25)],
        "time_hours": (1, 8),
        "typical_yield": (55, 85),
    },
}


class ReactionConditionPredictor:
    """Predicts optimal reaction conditions for retrosynthetic steps."""

    def __init__(self):
        self.conditions_db = REACTION_CONDITIONS_DB

    def predict_conditions(self, reaction_name: str, product_smiles: str,
                           reactant_smiles: str) -> Dict:
        """Predict optimal reaction conditions for a given transformation."""
        # Look up known conditions
        if reaction_name in self.conditions_db:
            conditions = self.conditions_db[reaction_name]
        else:
            conditions = self._infer_conditions(product_smiles, reactant_smiles)

        # Score solvents based on molecular properties
        solvent_scores = self._rank_solvents(
            conditions.get("solvents", []),
            product_smiles, reactant_smiles
        )

        # Optimal temperature prediction
        temp_range = conditions.get("temperature_range", (20, 80))
        optimal_temp = self._predict_temperature(
            temp_range, product_smiles, reactant_smiles
        )

        # Catalyst recommendation
        catalyst_scores = conditions.get("catalysts", [("None", 1.0)])

        # Estimate yield
        yield_range = conditions.get("typical_yield", (40, 70))
        estimated_yield = self._estimate_yield(
            yield_range, product_smiles, reactant_smiles
        )

        result = {
            "reaction_name": reaction_name,
            "recommended_solvent": solvent_scores[0] if solvent_scores else ("Unknown", 0.5),
            "all_solvents": solvent_scores,
            "optimal_temperature_C": optimal_temp,
            "temperature_range_C": temp_range,
            "recommended_catalyst": catalyst_scores[0] if catalyst_scores else ("None", 1.0),
            "all_catalysts": catalyst_scores,
            "estimated_time_hours": conditions.get("time_hours", (4, 24)),
            "estimated_yield_percent": estimated_yield,
        }

        if "base" in conditions:
            result["recommended_base"] = conditions["base"][0]
            result["all_bases"] = conditions["base"]

        return result

    def _rank_solvents(self, solvents: List[Tuple], product_smi: str,
                       reactant_smi: str) -> List[Tuple]:
        """Rank solvents considering molecular polarity."""
        if not solvents:
            return [("DCM", 0.5)]

        prod_mol = Chem.MolFromSmiles(product_smi)
        if prod_mol is None:
            return solvents

        logp = Crippen.MolLogP(prod_mol)
        tpsa = Descriptors.TPSA(prod_mol)

        # Adjust solvent scores based on polarity matching
        polar_solvents = {"DMF", "DMSO", "H2O", "MeOH", "EtOH", "NMP", "DMA"}
        adjusted = []
        for solvent, score in solvents:
            adj = score
            solvent_base = solvent.split("/")[0]
            is_polar = solvent_base in polar_solvents
            if tpsa > 80 and is_polar:
                adj *= 1.2
            elif tpsa < 40 and not is_polar:
                adj *= 1.2
            adjusted.append((solvent, round(adj, 3)))

        adjusted.sort(key=lambda x: x[1], reverse=True)
        return adjusted

    def _predict_temperature(self, temp_range: Tuple, product_smi: str,
                             reactant_smi: str) -> float:
        """Predict optimal temperature within the given range."""
        prod_mol = Chem.MolFromSmiles(product_smi)
        if prod_mol is None:
            return sum(temp_range) / 2

        mw = Descriptors.ExactMolWt(prod_mol)
        n_rings = prod_mol.GetRingInfo().NumRings()

        # Higher MW and more rings → higher temperature
        t_low, t_high = temp_range
        position = 0.5
        if mw > 300:
            position += 0.15
        if n_rings > 3:
            position += 0.1
        position = min(max(position, 0.0), 1.0)
        return round(t_low + (t_high - t_low) * position)

    def _estimate_yield(self, yield_range: Tuple, product_smi: str,
                        reactant_smi: str) -> Tuple[float, float]:
        """Estimate reaction yield range."""
        prod_mol = Chem.MolFromSmiles(product_smi)
        if prod_mol is None:
            return yield_range

        # Complex molecules tend to have lower yields
        n_atoms = prod_mol.GetNumHeavyAtoms()
        n_stereo = len(Chem.FindMolChiralCenters(prod_mol, includeUnassigned=True))

        y_low, y_high = yield_range
        if n_atoms > 25:
            y_high -= 10
        if n_stereo > 1:
            y_high -= 5 * n_stereo
            y_low -= 5

        return (max(10, round(y_low)), max(20, round(y_high)))

    def _infer_conditions(self, product_smi: str, reactant_smi: str) -> Dict:
        """Infer generic conditions when no template match is found."""
        return {
            "solvents": [("THF", 0.30), ("DCM", 0.25), ("DMF", 0.25), ("MeOH", 0.20)],
            "temperature_range": (20, 80),
            "catalysts": [("None specified", 1.0)],
            "time_hours": (4, 24),
            "typical_yield": (40, 70),
        }

    def predict_for_route(self, route_steps: List[Dict]) -> List[Dict]:
        """Predict conditions for each step in a multi-step route."""
        enriched_steps = []
        for step in route_steps:
            conditions = self.predict_conditions(
                step.get("reaction", "Unknown"),
                step.get("product", ""),
                step.get("reactants", ""),
            )
            enriched_step = {**step, "conditions": conditions}
            enriched_steps.append(enriched_step)
        return enriched_steps
