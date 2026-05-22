"""
RMG-style Automated VOC Oxidation Reaction Network Generator
Generates gas-phase chemical reaction pathways for SOA precursors.
"""
import numpy as np
import networkx as nx
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import json
import logging

logger = logging.getLogger(__name__)


@dataclass
class Species:
    name: str
    formula: str
    smiles: str
    molecular_weight: float  # g/mol
    n_carbons: int
    n_oxygens: int
    n_doubles: int          # double bonds
    vapor_pressure: float   # Pa at 298K (estimated)
    henry_const: float      # M/atm
    is_radical: bool = False
    generation: int = 0     # oxidation generation


@dataclass
class Reaction:
    reactants: List[str]
    products: List[str]
    rate_constant: float    # cm3/molecule/s or s-1
    reaction_type: str
    activation_energy: float = 0.0  # kJ/mol
    branching_ratio: float = 1.0


# ── Primary VOC precursors (terpenes + isoprene) ──────────────────────────────
PRIMARY_VOCS: Dict[str, Species] = {
    "alpha_pinene": Species(
        name="alpha-pinene", formula="C10H16", smiles="CC1=CCC2CC1CC2(C)C",
        molecular_weight=136.23, n_carbons=10, n_oxygens=0, n_doubles=1,
        vapor_pressure=632.0, henry_const=0.023, generation=0
    ),
    "beta_pinene": Species(
        name="beta-pinene", formula="C10H16", smiles="C=C1CCC2CC1CC2(C)C",
        molecular_weight=136.23, n_carbons=10, n_oxygens=0, n_doubles=1,
        vapor_pressure=852.0, henry_const=0.017, generation=0
    ),
    "limonene": Species(
        name="limonene", formula="C10H16", smiles="CC1=CCC(CC1)C(=C)C",
        molecular_weight=136.23, n_carbons=10, n_oxygens=0, n_doubles=2,
        vapor_pressure=201.0, henry_const=0.058, generation=0
    ),
    "isoprene": Species(
        name="isoprene", formula="C5H8", smiles="C=CC(=C)C",
        molecular_weight=68.12, n_carbons=5, n_oxygens=0, n_doubles=2,
        vapor_pressure=74000.0, henry_const=0.0013, generation=0
    ),
    "toluene": Species(
        name="toluene", formula="C7H8", smiles="Cc1ccccc1",
        molecular_weight=92.14, n_carbons=7, n_oxygens=0, n_doubles=4,
        vapor_pressure=3793.0, henry_const=0.15, generation=0
    ),
}

# ── Atmospheric oxidants ──────────────────────────────────────────────────────
OXIDANTS = {
    "OH":  {"conc_molec_cm3": 2.0e6,   "units": "molecule/cm3"},
    "O3":  {"conc_molec_cm3": 7.4e11,  "units": "molecule/cm3"},   # 30 ppb
    "NO3": {"conc_molec_cm3": 2.5e8,   "units": "molecule/cm3"},
    "NO":  {"conc_molec_cm3": 2.5e10,  "units": "molecule/cm3"},   # 1 ppb
    "NO2": {"conc_molec_cm3": 6.2e10,  "units": "molecule/cm3"},   # 2.5 ppb
    "HO2": {"conc_molec_cm3": 1.0e8,   "units": "molecule/cm3"},
}

# ── Experimental / literature rate constants [cm3 molecule-1 s-1] ─────────────
# Source: NIST Chemical Kinetics Database, Atkinson et al. (2006)
RATE_CONSTANTS = {
    # OH reactions
    ("alpha_pinene", "OH"):  5.33e-11,
    ("beta_pinene",  "OH"):  7.89e-11,
    ("limonene",     "OH"):  1.71e-10,
    ("isoprene",     "OH"):  1.00e-10,
    ("toluene",      "OH"):  5.63e-12,
    # O3 reactions
    ("alpha_pinene", "O3"):  8.66e-17,
    ("beta_pinene",  "O3"):  1.50e-17,
    ("limonene",     "O3"):  2.00e-16,
    ("isoprene",     "O3"):  1.27e-17,
    ("toluene",      "O3"):  0.0,         # negligible
    # NO3 reactions
    ("alpha_pinene", "NO3"): 6.16e-12,
    ("beta_pinene",  "NO3"): 2.51e-12,
    ("limonene",     "NO3"): 1.22e-11,
    ("isoprene",     "NO3"): 6.78e-13,
    ("toluene",      "NO3"): 0.0,
}

# ── Generation-1 oxidation products (simplified structural isomers) ───────────
GENERATION1_PRODUCTS = {
    "alpha_pinene": {
        "OH": [
            ("pinanediol",      "C10H18O2",  0.040, 0.30),  # (name, formula, Psat_Pa, yield)
            ("alpha_pin_OH",    "C10H16O",   1.20,  0.25),
            ("pinic_acid",      "C9H14O4",   1.2e-4,0.10),
            ("pinonic_acid",    "C10H16O3",  0.072, 0.15),
            ("norpinic_acid",   "C8H12O4",   2.0e-5,0.05),
        ],
        "O3": [
            ("pinic_acid",      "C9H14O4",   1.2e-4,0.20),
            ("pinaldehyde",     "C10H16O",   9.0,   0.30),
            ("norpinaldehyde",  "C9H14O",    15.0,  0.15),
            ("acetone",         "C3H6O",     30800, 0.25),
            ("cis_pinonic",     "C10H16O3",  0.072, 0.10),
        ],
        "NO3": [
            ("alpha_pin_nitrate","C10H17NO4", 0.01,  0.50),
            ("pinonaldehyde",   "C10H16O",   9.0,   0.30),
        ],
    },
    "beta_pinene": {
        "OH": [
            ("nopinaldehyde",   "C9H14O",    22.0,  0.35),
            ("beta_pin_OH",     "C10H16O",   2.5,   0.30),
            ("pinic_acid",      "C9H14O4",   1.2e-4,0.10),
        ],
        "O3": [
            ("nopinone",        "C9H14O",    18.0,  0.45),
            ("formaldehyde",    "CH2O",      1.8e5, 0.50),
            ("pinic_acid",      "C9H14O4",   1.2e-4,0.05),
        ],
        "NO3": [
            ("beta_pin_nitrate","C10H17NO4", 0.02,  0.55),
        ],
    },
    "limonene": {
        "OH": [
            ("limonene_OH",     "C10H16O",   0.80,  0.40),
            ("limonic_acid",    "C9H14O4",   5.0e-5,0.12),
            ("limonaketone",    "C9H14O3",   0.30,  0.18),
        ],
        "O3": [
            ("limonic_acid",    "C9H14O4",   5.0e-5,0.25),
            ("7_OH_lim",        "C10H18O2",  0.15,  0.20),
            ("LIMAL",           "C10H16O2",  2.0,   0.25),
        ],
        "NO3": [
            ("limonene_nitrate","C10H17NO4", 0.03,  0.60),
        ],
    },
    "isoprene": {
        "OH": [
            ("ISOPOOH",         "C5H10O3",   4.0,   0.25),
            ("methacrolein",    "C4H6O",     9050,  0.23),
            ("MVK",             "C4H6O",     12300, 0.32),
            ("ISOP_OOH",        "C5H10O3",   4.0,   0.08),
            ("2MGA",            "C4H8O4",    0.50,  0.06),
        ],
        "O3": [
            ("methacrolein",    "C4H6O",     9050,  0.20),
            ("formaldehyde",    "CH2O",      1.8e5, 0.60),
            ("MVK",             "C4H6O",     12300, 0.16),
        ],
        "NO3": [
            ("delta_ISOP_NO3",  "C5H9NO4",   0.12,  0.70),
            ("beta_ISOP_NO3",   "C5H9NO4",   0.12,  0.30),
        ],
    },
    "toluene": {
        "OH": [
            ("cresol",          "C7H8O",     165.0, 0.18),
            ("benzaldehyde",    "C7H6O",     170.0, 0.07),
            ("toluene_RO2",     "C7H7O5",    0.05,  0.30),
            ("DHBO",            "C7H8O2",    8.0,   0.15),
        ],
        "O3": [],
        "NO3": [],
    },
}

# ── Generation-2 products (further oxidation) ─────────────────────────────────
GENERATION2_PRODUCTS = {
    "pinonaldehyde":  [("pinic_acid", 0.40), ("norpinic_acid", 0.30)],
    "methacrolein":   [("2MGA", 0.12), ("methylglyoxal", 0.25)],
    "MVK":            [("methylglyoxal", 0.30), ("2MGA", 0.08)],
    "cresol":         [("methylnitrophenol", 0.15), ("ring_frag_products", 0.50)],
    "pinaldehyde":    [("pinic_acid", 0.35), ("norpinaldehyde", 0.25)],
}


class ReactionNetworkGenerator:
    """
    Automated chemical reaction network generator (RMG-inspired).
    Builds directed graph of VOC oxidation pathways.
    """

    def __init__(self, max_generations: int = 3):
        self.max_generations = max_generations
        self.graph = nx.DiGraph()
        self.species_db: Dict[str, dict] = {}
        self.reactions: List[Reaction] = []
        self.stats: Dict[str, int] = {}

    def generate_network(self, voc_list: List[str]) -> nx.DiGraph:
        """Generate the full reaction network for the given VOC list."""
        # Add primary VOC nodes
        for voc in voc_list:
            if voc in PRIMARY_VOCS:
                sp = PRIMARY_VOCS[voc]
                self.graph.add_node(
                    sp.name,
                    formula=sp.formula,
                    MW=sp.molecular_weight,
                    Psat=sp.vapor_pressure,
                    generation=0,
                    is_radical=False,
                    n_carbons=sp.n_carbons,
                    n_oxygens=sp.n_oxygens,
                )
                self.species_db[sp.name] = {
                    "formula": sp.formula,
                    "MW": sp.molecular_weight,
                    "Psat": sp.vapor_pressure,
                }

        # Generation 1
        self._add_generation1(voc_list)
        # Generation 2
        self._add_generation2()

        self.stats = {
            "n_species":   self.graph.number_of_nodes(),
            "n_reactions": self.graph.number_of_edges(),
            "n_primary":   len(voc_list),
        }
        logger.info(f"Network built: {self.stats}")
        return self.graph

    def _add_generation1(self, voc_list: List[str]):
        for voc in voc_list:
            if voc not in GENERATION1_PRODUCTS:
                continue
            sp = PRIMARY_VOCS[voc]
            for oxidant, products in GENERATION1_PRODUCTS[voc].items():
                k = RATE_CONSTANTS.get((voc, oxidant), 0.0)
                if k == 0.0:
                    continue
                for prod_name, prod_formula, psat, yield_frac in products:
                    if prod_name not in self.graph:
                        # estimate MW from formula
                        mw = self._estimate_mw(prod_formula)
                        n_c = prod_formula.count("C") if "C" in prod_formula else 0
                        n_o = prod_formula.count("O") if "O" in prod_formula else 0
                        self.graph.add_node(
                            prod_name,
                            formula=prod_formula,
                            MW=mw,
                            Psat=psat,
                            generation=1,
                            is_radical=False,
                            n_carbons=n_c,
                            n_oxygens=n_o,
                        )
                        self.species_db[prod_name] = {
                            "formula": prod_formula,
                            "MW": mw,
                            "Psat": psat,
                        }
                    self.graph.add_edge(
                        sp.name, prod_name,
                        oxidant=oxidant,
                        rate_constant=k,
                        yield_frac=yield_frac,
                        generation=1,
                        reaction_type=f"VOC+{oxidant}",
                    )
                    rxn = Reaction(
                        reactants=[sp.name, oxidant],
                        products=[prod_name],
                        rate_constant=k,
                        reaction_type=f"VOC+{oxidant}",
                        branching_ratio=yield_frac,
                    )
                    self.reactions.append(rxn)

    def _add_generation2(self):
        for parent, children in GENERATION2_PRODUCTS.items():
            if parent not in self.graph:
                continue
            for child_name, yield_frac in children:
                if child_name not in self.graph:
                    self.graph.add_node(
                        child_name,
                        formula="CxHyOz",
                        MW=150.0,
                        Psat=0.01,
                        generation=2,
                        is_radical=False,
                        n_carbons=6,
                        n_oxygens=3,
                    )
                k2 = 5e-12  # typical gen-2 OH rate constant
                self.graph.add_edge(
                    parent, child_name,
                    oxidant="OH",
                    rate_constant=k2,
                    yield_frac=yield_frac,
                    generation=2,
                    reaction_type="gen2_oxidation",
                )
                rxn = Reaction(
                    reactants=[parent, "OH"],
                    products=[child_name],
                    rate_constant=k2,
                    reaction_type="gen2_oxidation",
                    branching_ratio=yield_frac,
                )
                self.reactions.append(rxn)

    @staticmethod
    def _estimate_mw(formula: str) -> float:
        """Rough MW estimate from molecular formula string."""
        import re
        atoms = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}
        mw = 0.0
        for sym, mass in atoms.items():
            m = re.search(rf"{sym}(\d*)", formula)
            if m:
                n = int(m.group(1)) if m.group(1) else 1
                mw += n * mass
        return mw

    def get_soa_precursors(self, psat_threshold: float = 10.0) -> List[str]:
        """Return species with Psat < threshold (Pa) — potential SOA formers."""
        return [
            n for n, d in self.graph.nodes(data=True)
            if d.get("generation", 0) >= 1 and d.get("Psat", 1e9) < psat_threshold
        ]

    def export_json(self, path: str):
        data = {
            "nodes": [
                {"id": n, **d}
                for n, d in self.graph.nodes(data=True)
            ],
            "edges": [
                {"source": u, "target": v, **d}
                for u, v, d in self.graph.edges(data=True)
            ],
            "stats": self.stats,
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Network exported to {path}")
        return data
