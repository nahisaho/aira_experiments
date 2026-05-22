"""
Materials Database: Ionic radii, electronegativities, and known perovskite properties.
Sources: Shannon (1976) ionic radii; experimental/DFT literature values.
"""

import numpy as np

# Shannon ionic radii (Angstrom), coordination number in parentheses
# Format: {element: {oxidation_state: {CN: radius}}}
IONIC_RADII = {
    # A-site cations (12-coordinate for cubic perovskite)
    "Cs": {1: {12: 1.88, 6: 1.74}},
    "Rb": {1: {12: 1.72, 6: 1.52}},
    "K":  {1: {12: 1.64, 6: 1.38}},
    "MA": {1: {12: 2.17}},   # methylammonium (effective)
    "FA": {1: {12: 2.53}},   # formamidinium (effective)
    "EA": {1: {12: 2.74}},   # ethylammonium
    "GA": {1: {12: 2.78}},   # guanidinium
    "DMA": {1: {12: 2.72}},  # dimethylammonium
    # B-site cations (6-coordinate octahedral)
    "Sn": {2: {6: 1.10}, 4: {6: 0.69}},
    "Ge": {2: {6: 0.73}, 4: {6: 0.53}},
    "Bi": {3: {6: 1.03}},
    "Sb": {3: {6: 0.76}},
    "Ti": {4: {6: 0.605}},
    "Zr": {4: {6: 0.72}},
    "Hf": {4: {6: 0.71}},
    "In": {3: {6: 0.80}},
    "Ga": {3: {6: 0.62}},
    "Mn": {2: {6: 0.83}},
    "Fe": {2: {6: 0.78}, 3: {6: 0.645}},
    "Cu": {2: {6: 0.73}},
    "Ag": {1: {6: 1.15}},
    "Au": {3: {6: 0.85}},
    "Pb": {2: {6: 1.19}},   # reference (toxic)
    # X-site anions (6-coordinate)
    "I":  {-1: {6: 2.20}},
    "Br": {-1: {6: 1.96}},
    "Cl": {-1: {6: 1.81}},
    "F":  {-1: {6: 1.33}},
    "SCN": {-1: {6: 2.17}},
    "HCOO": {-1: {6: 2.08}},  # formate
    "BF4":  {-1: {6: 2.23}},
}

# Electronegativity (Pauling scale)
ELECTRONEGATIVITY = {
    "Cs": 0.79, "Rb": 0.82, "K": 0.82, "MA": 2.3, "FA": 2.3,
    "EA": 2.3,  "GA": 2.3,  "DMA": 2.3,
    "Sn": 1.96, "Ge": 2.01, "Bi": 2.02, "Sb": 2.05,
    "Ti": 1.54, "Zr": 1.33, "Hf": 1.30,
    "In": 1.78, "Ga": 1.81, "Mn": 1.55, "Fe": 1.83,
    "Cu": 1.90, "Ag": 1.93, "Au": 2.54, "Pb": 2.33,
    "I": 2.66,  "Br": 2.96, "Cl": 3.16, "F": 3.98,
    "SCN": 2.58, "HCOO": 2.2, "BF4": 3.0,
}

# Known experimental / high-quality DFT band gaps (eV)
KNOWN_BANDGAPS = {
    # Sn-based
    ("MA", "Sn", "I"):   {"Eg": 1.20, "source": "exp"},
    ("FA", "Sn", "I"):   {"Eg": 1.41, "source": "exp"},
    ("Cs", "Sn", "I"):   {"Eg": 1.31, "source": "exp"},
    ("MA", "Sn", "Br"):  {"Eg": 2.15, "source": "exp"},
    ("Cs", "Sn", "Br"):  {"Eg": 1.75, "source": "exp"},
    ("MA", "Sn", "Cl"):  {"Eg": 2.90, "source": "DFT"},
    # Ge-based
    ("MA", "Ge", "I"):   {"Eg": 1.90, "source": "DFT"},
    ("Cs", "Ge", "I"):   {"Eg": 1.63, "source": "exp"},
    ("Rb", "Ge", "I"):   {"Eg": 1.80, "source": "DFT"},
    ("MA", "Ge", "Br"):  {"Eg": 2.70, "source": "DFT"},
    ("Cs", "Ge", "Br"):  {"Eg": 2.30, "source": "exp"},
    # Bi-based (A3Bi2X9 structure – adjusted)
    ("MA", "Bi", "I"):   {"Eg": 2.10, "source": "exp"},
    ("Cs", "Bi", "I"):   {"Eg": 2.20, "source": "exp"},
    ("Rb", "Bi", "I"):   {"Eg": 2.30, "source": "DFT"},
    ("MA", "Bi", "Br"):  {"Eg": 2.60, "source": "DFT"},
    ("Cs", "Bi", "Br"):  {"Eg": 2.70, "source": "exp"},
    # Reference Pb (for comparison)
    ("MA", "Pb", "I"):   {"Eg": 1.55, "source": "exp"},
    ("FA", "Pb", "I"):   {"Eg": 1.48, "source": "exp"},
    ("Cs", "Pb", "I"):   {"Eg": 1.73, "source": "exp"},
    ("MA", "Pb", "Br"):  {"Eg": 2.29, "source": "exp"},
    ("Cs", "Pb", "Br"):  {"Eg": 2.36, "source": "exp"},
}

# PCE records (%) and stability notes from literature
DEVICE_RECORDS = {
    ("MA", "Sn", "I"):  {"PCE": 13.24, "Voc": 0.88, "Jsc": 21.0, "FF": 0.72, "stability": "moderate"},
    ("FA", "Sn", "I"):  {"PCE": 14.81, "Voc": 0.91, "Jsc": 23.0, "FF": 0.71, "stability": "moderate"},
    ("Cs", "Sn", "I"):  {"PCE":  8.5,  "Voc": 0.72, "Jsc": 16.8, "FF": 0.70, "stability": "low"},
    ("Cs", "Ge", "I"):  {"PCE":  4.94, "Voc": 0.74, "Jsc": 11.1, "FF": 0.60, "stability": "low"},
    ("MA", "Bi", "I"):  {"PCE":  3.17, "Voc": 0.81, "Jsc":  5.5, "FF": 0.72, "stability": "high"},
    ("Cs", "Bi", "I"):  {"PCE":  2.10, "Voc": 0.74, "Jsc":  4.2, "FF": 0.68, "stability": "high"},
    ("MA", "Pb", "I"):  {"PCE": 25.7,  "Voc": 1.18, "Jsc": 26.1, "FF": 0.83, "stability": "reference"},
}

def get_ionic_radius(element, ox_state, cn=6):
    """Return ionic radius with fallback for missing CN."""
    try:
        radii = IONIC_RADII[element][ox_state]
        if cn in radii:
            return radii[cn]
        return list(radii.values())[0]
    except KeyError:
        raise ValueError(f"Ionic radius not found: {element}({ox_state}+) CN={cn}")

def get_all_candidates():
    """Return list of all candidate ABX3 combinations."""
    A_sites = ["Cs", "Rb", "MA", "FA", "EA", "DMA"]
    B_sites_Sn = [("Sn", 2)]
    B_sites_Ge = [("Ge", 2)]
    B_sites_Bi = [("Bi", 3)]
    X_sites = ["I", "Br", "Cl"]

    candidates = []
    for A in A_sites:
        for B, ox in B_sites_Sn + B_sites_Ge + B_sites_Bi:
            for X in X_sites:
                candidates.append({
                    "A": A, "B": B, "X": X,
                    "B_ox": ox,
                    "formula": f"{A}{B}{X}3",
                    "system": B,
                })
    return candidates
