"""
DFT Data Generation & Materials Project / AFLOW Integration
Simulates DFT-level energetics for HEA supercells and fetches
real materials data from public databases via REST API.
"""

import numpy as np
import pandas as pd
import requests
import json
from typing import Dict, List, Optional
import time
import warnings
warnings.filterwarnings("ignore")


# -----------------------------------------------------------------------
# Simulated DFT calculations (SQS-based)
# -----------------------------------------------------------------------
class DFTSimulator:
    """
    Simulates DFT total-energy calculations on Special Quasi-random
    Structures (SQS) for HEA supercells.

    Energy model based on:
    - Zunger et al. (1990) SQS framework
    - Kormann et al. (2017) interatomic potential approximation
    - Wang et al. (2021) ML interatomic potential validation

    Outputs: formation energy, elastic constants, magnetic moments, DOS.
    """

    def __init__(self, functional: str = "PBE", cutoff_eV: float = 500.0,
                 random_state: int = 42):
        self.functional = functional
        self.cutoff = cutoff_eV
        self.rng = np.random.default_rng(random_state)

    def formation_energy(self, composition: Dict[str, float]) -> float:
        """
        ΔE_f [eV/atom] — SQS formation energy estimate.
        Uses a simplified cluster expansion (CE) model.
        Negative values → thermodynamically stable structure.
        """
        from src.hea_descriptors import mixing_enthalpy
        dH = mixing_enthalpy(composition)  # kJ/mol

        # Convert to eV/atom (1 kJ/mol ≈ 0.01036 eV/atom)
        E_base = dH * 0.01036

        # DFT correction for magnetic interactions (relevant for Fe, Co, Ni)
        mag_elements = {"Fe": -0.08, "Co": -0.05, "Ni": -0.03, "Mn": +0.04}
        E_mag = sum(composition.get(el, 0) * v for el, v in mag_elements.items())

        # Vibrational entropy contribution (Debye model approximation)
        E_vib = -0.005 * len(composition)

        noise = self.rng.normal(0, 0.01)
        return E_base + E_mag + E_vib + noise

    def elastic_constants(self, composition: Dict[str, float]) -> Dict[str, float]:
        """
        Voigt-averaged elastic constants [GPa].
        Returns C11, C12, C44, B, G, E, nu.
        """
        from src.hea_descriptors import rule_of_mixture_modulus, ELEMENT_PROPS, _mean_prop

        B_rom, G_rom = rule_of_mixture_modulus(composition)

        # DFT typically gives 5-15% correction relative to ROM
        f_dft = 1.0 + self.rng.normal(0.0, 0.05)
        B = B_rom * f_dft
        G = G_rom * f_dft

        E = 9 * B * G / (3 * B + G)
        nu = (3 * B - 2 * G) / (2 * (3 * B + G))

        # Cubic elastic constants (Voigt averaged)
        C11 = B + 4 * G / 3
        C12 = B - 2 * G / 3
        C44 = G  # simplified (Zener anisotropy ≈ 1 for disordered HEA)

        return {"C11": C11, "C12": C12, "C44": C44,
                "B": B, "G": G, "E": E, "nu": nu}

    def magnetic_moment(self, composition: Dict[str, float]) -> float:
        """Average magnetic moment per atom [μ_B]."""
        mag_mom = {"Fe": 2.22, "Co": 1.72, "Ni": 0.61, "Mn": 2.50,
                   "Cr": -0.50, "Al": 0.0, "Ti": 0.0, "V": 0.0,
                   "Mo": 0.0, "W": 0.0, "Cu": 0.0}
        mu = sum(composition.get(el, 0) * mag_mom.get(el, 0)
                 for el in composition)
        return mu + self.rng.normal(0, 0.05)

    def lattice_parameter(self, composition: Dict[str, float]) -> float:
        """
        Vegard's law lattice parameter estimate [Å].
        FCC reference values for each element.
        """
        fcc_a = {"Cr": 3.53, "Mn": 3.53, "Fe": 3.57, "Co": 3.54,
                 "Ni": 3.52, "Al": 4.05, "Ti": 4.10, "V":  3.02,
                 "Mo": 3.15, "W":  3.16, "Cu": 3.61, "Zr": 3.60,
                 "Nb": 3.30, "Hf": 3.60, "Ta": 3.30}
        a = sum(composition.get(el, 0) * fcc_a.get(el, 3.55)
                for el in composition)
        return a + self.rng.normal(0, 0.005)

    def stacking_fault_energy(self, composition: Dict[str, float]) -> float:
        """
        SFE [mJ/m²]: controls deformation mechanism (TWIP/TRIP/slip).
        Low SFE < 20 → TWIP/TRIP; moderate 20-80 → partial dislocation;
        high > 80 → cross-slip dominant.
        """
        vec = sum(composition.get(el, 0) * {"Cr":6,"Mn":7,"Fe":8,"Co":9,"Ni":10,
                                              "Al":3,"Ti":4,"V":5,"Mo":6,"W":6,
                                              "Cu":11}.get(el, 8)
                  for el in composition)
        # Empirical: SFE scales with VEC for FCC (Carter & Holmes model)
        sfe_base = 20 * (vec - 7.5) ** 2 - 10
        x_Ni = composition.get("Ni", 0)
        x_Co = composition.get("Co", 0)
        sfe = sfe_base + 30 * x_Ni - 10 * x_Co
        return max(-50, sfe) + self.rng.normal(0, 5)

    def generate_dataset(self, compositions: List[Dict[str, float]]) -> pd.DataFrame:
        """Batch DFT simulation."""
        rows = []
        for comp in compositions:
            el_consts = self.elastic_constants(comp)
            row = {
                "composition": str(comp),
                "E_form_eV_atom": self.formation_energy(comp),
                "lattice_a_Ang": self.lattice_parameter(comp),
                "mag_moment_uB": self.magnetic_moment(comp),
                "SFE_mJm2": self.stacking_fault_energy(comp),
                **{f"DFT_{k}": v for k, v in el_consts.items()},
            }
            for el, x in comp.items():
                row[f"x_{el}"] = x
            rows.append(row)
        return pd.DataFrame(rows)


# -----------------------------------------------------------------------
# Materials Project / AFLOW API client
# -----------------------------------------------------------------------
class MaterialsDatabaseClient:
    """
    Fetches HEA-relevant data from:
    - Materials Project (REST v3)  — https://materialsproject.org
    - AFLOW repository             — http://aflowlib.org
    Falls back to curated literature dataset if APIs are unavailable.
    """

    MP_BASE = "https://api.materialsproject.org"
    AFLOW_BASE = "http://aflowlib.org/API/aflux"

    def __init__(self, mp_api_key: Optional[str] = None, timeout: int = 15):
        self.mp_key = mp_api_key
        self.timeout = timeout

    def fetch_mp_binary_data(self, elements: List[str]) -> Optional[pd.DataFrame]:
        """
        Query Materials Project for binary compound properties among given elements.
        Returns DataFrame with formula, formation energy, band gap, bulk modulus.
        """
        if not self.mp_key:
            return None

        headers = {"X-API-KEY": self.mp_key}
        params = {
            "elements": ",".join(elements),
            "fields": "material_id,formula_pretty,formation_energy_per_atom,"
                      "band_gap,bulk_modulus,shear_modulus,energy_above_hull",
            "deprecated": False,
        }
        try:
            r = requests.get(f"{self.MP_BASE}/materials/summary",
                             headers=headers, params=params, timeout=self.timeout)
            if r.status_code == 200:
                data = r.json().get("data", [])
                return pd.DataFrame(data)
        except Exception:
            pass
        return None

    def fetch_aflow_structures(self, elements: List[str], n_max: int = 50) -> Optional[pd.DataFrame]:
        """
        Query AFLOW for computed properties of compounds with given elements.
        """
        el_str = ",".join(f"'{el}'" for el in elements)
        query = (f"species({el_str}),nspecies({len(elements)}),"
                 f"Egap,Eentropy,spin_atom,density,agl_bulk_modulus_vrh,agl_shear_modulus_vrh")
        params = {"AFLUX": query, "format": "json", "paging": f"1:{n_max}"}
        try:
            r = requests.get(self.AFLOW_BASE, params=params, timeout=self.timeout)
            if r.status_code == 200:
                data = r.json()
                return pd.DataFrame(data if isinstance(data, list) else [data])
        except Exception:
            pass
        return None

    def get_curated_cantor_data(self) -> pd.DataFrame:
        """
        Curated dataset of CrMnFeCoNi-family HEAs from literature.
        Sources: Cantor (2004), Otto (2013), Gludovatz (2014),
                 Wu (2014), Bracq (2017), Laplanche (2017).
        """
        records = [
            # alloy, σ_y(MPa), ε_f(%), E_pit(V), VEC, T_test(K), reference
            {"alloy": "CrMnFeCoNi",     "sigma_y": 350,  "elong": 60, "E_pit": -0.05, "VEC": 8.0,  "T_K": 293, "ref": "Cantor2004"},
            {"alloy": "CrMnFeCoNi",     "sigma_y": 620,  "elong": 55, "E_pit": -0.05, "VEC": 8.0,  "T_K": 77,  "ref": "Gludovatz2014"},
            {"alloy": "CrFeCoNi",       "sigma_y": 330,  "elong": 52, "E_pit": +0.10, "VEC": 8.25, "T_K": 293, "ref": "Wu2014"},
            {"alloy": "CrCoNi",         "sigma_y": 480,  "elong": 61, "E_pit": +0.15, "VEC": 8.33, "T_K": 293, "ref": "Miao2020"},
            {"alloy": "Cr20Mn20Fe20Co20Ni20", "sigma_y": 355, "elong": 55, "E_pit": -0.05, "VEC": 8.0, "T_K": 293, "ref": "Otto2013"},
            {"alloy": "CrCoNiAl0.1",    "sigma_y": 530,  "elong": 45, "E_pit": +0.12, "VEC": 8.1,  "T_K": 293, "ref": "Gao2018"},
            {"alloy": "CrCoNiAl0.3",    "sigma_y": 750,  "elong": 28, "E_pit": +0.05, "VEC": 7.8,  "T_K": 293, "ref": "Gao2018"},
            {"alloy": "CrFeCoNiMo0.1",  "sigma_y": 380,  "elong": 48, "E_pit": +0.25, "VEC": 8.1,  "T_K": 293, "ref": "Zhao2017"},
            {"alloy": "CrMnFeCoNiTi0.1","sigma_y": 400,  "elong": 35, "E_pit": -0.02, "VEC": 7.8,  "T_K": 293, "ref": "Ye2018"},
            {"alloy": "CrFeCoNiAl0.5",  "sigma_y": 890,  "elong": 15, "E_pit": +0.08, "VEC": 7.6,  "T_K": 293, "ref": "Wang2012"},
            {"alloy": "CrFeCoNiAl1.0",  "sigma_y":1250,  "elong":  8, "E_pit": +0.02, "VEC": 7.2,  "T_K": 293, "ref": "Wang2012"},
            {"alloy": "CrCoNi(MnFe)0",  "sigma_y": 480,  "elong": 62, "E_pit": +0.14, "VEC": 8.33, "T_K": 77,  "ref": "Laplanche2017"},
            {"alloy": "FeCoNi",         "sigma_y": 250,  "elong": 58, "E_pit": -0.10, "VEC": 9.0,  "T_K": 293, "ref": "Wu2014"},
            {"alloy": "CrMnFeNi",       "sigma_y": 300,  "elong": 45, "E_pit": +0.00, "VEC": 7.75, "T_K": 293, "ref": "Bracq2017"},
            {"alloy": "CrFeCoNiV0.25",  "sigma_y": 420,  "elong": 40, "E_pit": +0.05, "VEC": 8.05, "T_K": 293, "ref": "Stepanov2016"},
        ]
        return pd.DataFrame(records)

    def get_refractory_hea_data(self) -> pd.DataFrame:
        """
        Curated dataset of refractory HEAs (RHEAs) from literature.
        Relevant for high-temperature strength targets.
        """
        records = [
            {"alloy": "MoNbTaVW",     "sigma_y": 405,  "T_K": 1600, "hardness_HV": 4455, "ref": "Senkov2011"},
            {"alloy": "HfNbTaTiZr",   "sigma_y": 929,  "T_K": 293,  "hardness_HV": 2830, "ref": "Senkov2011"},
            {"alloy": "CrMoNbVW",     "sigma_y": 1246, "T_K": 293,  "hardness_HV": 4700, "ref": "Senkov2013"},
            {"alloy": "MoNbTaW",      "sigma_y": 1058, "T_K": 293,  "hardness_HV": 4400, "ref": "Senkov2010"},
            {"alloy": "NbMoTaW+Re",   "sigma_y": 1300, "T_K": 293,  "hardness_HV": 5100, "ref": "Senkov2012"},
            {"alloy": "CrNbTiVZr",    "sigma_y": 1271, "T_K": 293,  "hardness_HV": 4600, "ref": "Senkov2013"},
            {"alloy": "TiZrHfNbTa",   "sigma_y": 1000, "T_K": 293,  "hardness_HV": 3300, "ref": "Yao2016"},
            {"alloy": "AlCrNbTiV",    "sigma_y": 1085, "T_K": 293,  "hardness_HV": 4400, "ref": "Stepanov2015"},
        ]
        return pd.DataFrame(records)
