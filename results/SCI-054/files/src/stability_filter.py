"""
Module 5: Water stability and synthesizability prediction filters.

Classifies MOFs by:
1. Water stability (hydrothermal/hydrolytic)
2. Synthetic accessibility (known precursors, reaction conditions)
3. Chemical stability under DAC operating conditions
"""
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class StabilityScore:
    """Stability assessment for a single MOF."""
    mof_id: str = ""
    water_stability_score: float = 0.0     # 0-1 probability
    water_stability_class: str = "unknown" # stable / moderate / unstable
    synthesizability_score: float = 0.0    # 0-1 probability
    thermal_stability_K: float = 0.0       # decomposition temperature
    acid_stability: bool = False
    mechanical_stability_GPa: float = 0.0  # bulk modulus
    overall_dac_suitability: float = 0.0   # composite score for DAC
    failure_reasons: List[str] = field(default_factory=list)

    @property
    def passes_dac_filter(self) -> bool:
        return (self.water_stability_score >= 0.5 and
                self.synthesizability_score >= 0.3 and
                self.thermal_stability_K >= 373.0 and
                len(self.failure_reasons) == 0)


class WaterStabilityPredictor:
    """
    Predict MOF water stability using rule-based and ML approaches.

    Key factors (from literature):
    1. Metal-ligand bond strength (thermodynamic stability)
    2. Coordination geometry shielding (kinetic stability)
    3. Hydrophobicity of pore surface
    4. Framework flexibility/rigidity
    """

    # Metal-node water stability rankings (empirical)
    METAL_STABILITY = {
        "Zr": 0.95,   # Zr-oxo clusters: exceptionally water-stable
        "Hf": 0.93,
        "Al": 0.85,   # Al-carboxylate: good stability
        "Ti": 0.80,
        "Cr": 0.75,   # Cr(III): kinetically inert
        "Fe": 0.50,   # Fe(III): moderate, Fe(II): poor
        "In": 0.60,
        "Ni": 0.55,
        "Co": 0.45,
        "Cu": 0.40,   # Cu paddlewheel: water-sensitive
        "Zn": 0.35,   # Zn-carboxylate: generally water-unstable
        "Cd": 0.25,
        "Mn": 0.30,
        "Mg": 0.20,
    }

    # Linker stability contributions
    LINKER_STABILITY = {
        "azolate": 0.90,          # strong M-N bonds
        "carboxylate_zr": 0.85,   # Zr-carboxylate: very stable
        "phosphonate": 0.80,
        "carboxylate": 0.50,      # general carboxylate
        "mixed_N_O": 0.65,
        "thiolate": 0.30,
        "unknown": 0.40,
    }

    # Known water-stable MOF families
    STABLE_FAMILIES = {
        "UiO-66", "UiO-67", "UiO-68",    # Zr-based
        "MOF-808", "NU-1000", "NU-1200",  # Zr-based
        "MIL-53", "MIL-100", "MIL-101",   # Al/Cr-based
        "MIL-125",                          # Ti-based
        "ZIF-8", "ZIF-67", "ZIF-90",      # Zn/Co imidazolate
        "PCN-222", "PCN-224", "PCN-250",   # Zr-porphyrin
        "BUT-12", "BUT-13",                # Azolate-based
    }

    def predict(self, mof_id: str, metal_type: str, linker_type: str,
                has_oms: bool = False, porosity: float = 0.5,
                functional_groups: Optional[List[str]] = None) -> float:
        """Predict water stability score (0-1)."""
        metal_score = self.METAL_STABILITY.get(metal_type, 0.4)

        # Adjust linker classification for Zr carboxylates
        if metal_type in ("Zr", "Hf") and linker_type == "carboxylate":
            linker_type = "carboxylate_zr"
        linker_score = self.LINKER_STABILITY.get(linker_type, 0.4)

        # Composite score with weights
        score = 0.5 * metal_score + 0.3 * linker_score + 0.2 * 0.5

        # Adjustments
        if has_oms and metal_type not in ("Zr", "Hf", "Al", "Cr"):
            score -= 0.15  # OMS on less-stable metals → water attack site

        if functional_groups:
            hydrophobic = {"CF3", "CH3", "F", "alkyl"}
            if any(fg in hydrophobic for fg in functional_groups):
                score += 0.10  # hydrophobic groups improve stability

        # Check known stable families
        for family in self.STABLE_FAMILIES:
            if family.lower() in mof_id.lower():
                score = max(score, 0.80)
                break

        return min(max(score, 0.0), 1.0)


class SynthesizabilityPredictor:
    """
    Predict synthetic accessibility of MOFs.

    Considers:
    1. Commercial availability of precursors
    2. Reported synthesis conditions
    3. Structural complexity
    4. Number of unique building blocks
    """

    COMMON_METALS = {
        "Zn", "Cu", "Zr", "Al", "Fe", "Co", "Ni", "Cr", "Mg", "Ca",
        "In", "Cd", "Mn", "Ti",
    }

    COMMON_LINKERS = {
        "BDC": 0.95,         # terephthalic acid
        "BTC": 0.90,         # trimesic acid
        "NDC": 0.85,         # naphthalenedicarboxylic acid
        "BPDC": 0.80,        # biphenyl-4,4'-dicarboxylic acid
        "imidazole": 0.90,
        "2-methylimidazole": 0.90,
        "pyrazole": 0.85,
        "triazole": 0.80,
        "TCPP": 0.60,        # porphyrin linkers
        "TATB": 0.55,
    }

    def predict(self, mof_id: str, metal_type: str, n_atom_types: int,
                n_atoms_per_uc: int, linker_type: str,
                source_db: str = "CoRE") -> float:
        """Predict synthesizability score (0-1)."""
        score = 0.5  # base score

        # Metal availability
        if metal_type in self.COMMON_METALS:
            score += 0.15

        # Structural complexity penalty
        if n_atoms_per_uc > 500:
            score -= 0.15
        elif n_atoms_per_uc > 200:
            score -= 0.05

        # Number of unique elements
        if n_atom_types <= 4:
            score += 0.10
        elif n_atom_types > 6:
            score -= 0.10

        # CoRE MOF structures are experimentally reported
        if source_db == "CoRE":
            score += 0.20  # known to be synthesizable

        # hMOF are hypothetical
        elif source_db == "hMOF":
            score -= 0.10

        # Linker type bonus
        if linker_type in ("carboxylate", "azolate"):
            score += 0.10

        return min(max(score, 0.0), 1.0)


class ThermalStabilityEstimator:
    """Estimate thermal decomposition temperature."""

    # Approximate decomposition temperatures by metal-linker class
    DECOMP_TEMPS = {
        ("Zr", "carboxylate"): 773,   # K — UiO-66 type
        ("Zr", "mixed_N_O"): 723,
        ("Al", "carboxylate"): 673,
        ("Cr", "carboxylate"): 623,
        ("Fe", "carboxylate"): 573,
        ("Ti", "carboxylate"): 623,
        ("Zn", "carboxylate"): 523,
        ("Cu", "carboxylate"): 473,
        ("Zn", "azolate"): 673,       # ZIF-type
        ("Co", "azolate"): 623,
        ("Ni", "azolate"): 648,
    }

    def estimate(self, metal_type: str, linker_type: str) -> float:
        key = (metal_type, linker_type)
        if key in self.DECOMP_TEMPS:
            return self.DECOMP_TEMPS[key]
        # Default based on metal alone
        metal_defaults = {"Zr": 700, "Al": 623, "Cr": 573, "Fe": 523,
                          "Ti": 573, "Zn": 473, "Cu": 423, "Co": 523,
                          "Ni": 548, "Mn": 473, "Mg": 423}
        return metal_defaults.get(metal_type, 473)


class StabilityFilter:
    """Combined stability and synthesizability filter for MOF screening."""

    def __init__(self, water_threshold: float = 0.5,
                 synth_threshold: float = 0.3,
                 min_thermal_K: float = 373.0):
        self.water_predictor = WaterStabilityPredictor()
        self.synth_predictor = SynthesizabilityPredictor()
        self.thermal_estimator = ThermalStabilityEstimator()
        self.water_threshold = water_threshold
        self.synth_threshold = synth_threshold
        self.min_thermal_K = min_thermal_K

    def evaluate(self, mof_id: str, metal_type: str, linker_type: str,
                 has_oms: bool, porosity: float, n_atom_types: int,
                 n_atoms_per_uc: int, source_db: str,
                 functional_groups: Optional[List[str]] = None) -> StabilityScore:
        """Evaluate stability and synthesizability."""
        score = StabilityScore(mof_id=mof_id)

        score.water_stability_score = self.water_predictor.predict(
            mof_id, metal_type, linker_type, has_oms, porosity, functional_groups
        )

        if score.water_stability_score >= 0.7:
            score.water_stability_class = "stable"
        elif score.water_stability_score >= 0.4:
            score.water_stability_class = "moderate"
        else:
            score.water_stability_class = "unstable"

        score.synthesizability_score = self.synth_predictor.predict(
            mof_id, metal_type, n_atom_types, n_atoms_per_uc,
            linker_type, source_db
        )

        score.thermal_stability_K = self.thermal_estimator.estimate(
            metal_type, linker_type
        )

        # Failure reasons
        if score.water_stability_score < self.water_threshold:
            score.failure_reasons.append(
                f"Water stability too low ({score.water_stability_score:.2f} "
                f"< {self.water_threshold})"
            )
        if score.synthesizability_score < self.synth_threshold:
            score.failure_reasons.append(
                f"Synthesizability too low ({score.synthesizability_score:.2f} "
                f"< {self.synth_threshold})"
            )
        if score.thermal_stability_K < self.min_thermal_K:
            score.failure_reasons.append(
                f"Thermal stability insufficient ({score.thermal_stability_K:.0f} K "
                f"< {self.min_thermal_K:.0f} K)"
            )

        # Overall DAC suitability
        score.overall_dac_suitability = (
            0.4 * score.water_stability_score +
            0.3 * min(score.thermal_stability_K / 773.0, 1.0) +
            0.3 * score.synthesizability_score
        )

        return score

    def filter_batch(self, mof_list: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """Filter a batch of MOFs, returning (passed, failed)."""
        passed, failed = [], []
        for mof in mof_list:
            score = self.evaluate(
                mof_id=mof.get("mof_id", ""),
                metal_type=mof.get("metal_type", ""),
                linker_type=mof.get("linker_type", ""),
                has_oms=mof.get("has_oms", False),
                porosity=mof.get("porosity", 0.5),
                n_atom_types=mof.get("n_atom_types", 4),
                n_atoms_per_uc=mof.get("n_atoms_per_uc", 100),
                source_db=mof.get("source_db", "CoRE"),
            )
            mof["stability_score"] = score
            if score.passes_dac_filter:
                passed.append(mof)
            else:
                failed.append(mof)

        logger.info(f"Stability filter: {len(passed)}/{len(mof_list)} passed")
        return passed, failed
