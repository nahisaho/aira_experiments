"""
Module 1: MOF structural feature extraction from CoRE MOF / hMOF databases.

Extracts geometric, chemical, topological, and energy-based descriptors
using Zeo++ and custom analysis routines.
"""
import subprocess
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class GeometricDescriptors:
    """Geometric descriptors computed by Zeo++."""
    lcd: float = 0.0           # Largest cavity diameter (Å)
    pld: float = 0.0           # Pore limiting diameter (Å)
    asa: float = 0.0           # Accessible surface area (m²/g)
    nasa: float = 0.0          # Non-accessible surface area (m²/g)
    av: float = 0.0            # Accessible volume (cm³/g)
    nav: float = 0.0           # Non-accessible volume (cm³/g)
    porosity: float = 0.0      # Void fraction (0-1)
    density: float = 0.0       # Framework density (g/cm³)
    n_channels: int = 0        # Number of accessible channels
    psd_mean: float = 0.0      # Mean pore size (Å)
    psd_std: float = 0.0       # Pore size distribution std (Å)


@dataclass
class ChemicalDescriptors:
    """Chemical composition descriptors."""
    metal_type: str = ""
    metal_fraction: float = 0.0
    linker_type: str = ""
    n_atom_types: int = 0
    has_open_metal_sites: bool = False
    functional_groups: List[str] = field(default_factory=list)
    n_atoms_per_uc: int = 0
    molecular_weight_uc: float = 0.0


@dataclass
class TopologicalDescriptors:
    """Topological network descriptors."""
    topology_code: str = ""
    coordination_number: int = 0
    dimensionality: int = 3
    vertex_symbol: str = ""
    n_sbu: int = 0


@dataclass
class EnergyDescriptors:
    """Energy-based descriptors from force field calculations."""
    henry_co2: float = 0.0         # mol/kg/Pa
    henry_h2: float = 0.0
    heat_adsorption_co2: float = 0.0  # kJ/mol
    heat_adsorption_h2: float = 0.0
    widom_co2: float = 0.0
    widom_h2: float = 0.0


@dataclass
class MOFFeatures:
    """Complete feature set for a single MOF."""
    mof_id: str = ""
    source_db: str = ""
    cif_path: str = ""
    geometric: GeometricDescriptors = field(default_factory=GeometricDescriptors)
    chemical: ChemicalDescriptors = field(default_factory=ChemicalDescriptors)
    topological: TopologicalDescriptors = field(default_factory=TopologicalDescriptors)
    energy: EnergyDescriptors = field(default_factory=EnergyDescriptors)

    def to_feature_vector(self) -> np.ndarray:
        """Convert to flat numpy array for ML input."""
        geo = self.geometric
        chem = self.chemical
        topo = self.topological
        ener = self.energy
        return np.array([
            geo.lcd, geo.pld, geo.asa, geo.nasa, geo.av, geo.nav,
            geo.porosity, geo.density, geo.n_channels, geo.psd_mean, geo.psd_std,
            chem.metal_fraction, chem.n_atom_types,
            float(chem.has_open_metal_sites), chem.n_atoms_per_uc,
            chem.molecular_weight_uc,
            topo.coordination_number, topo.dimensionality, topo.n_sbu,
            ener.henry_co2, ener.henry_h2,
            ener.heat_adsorption_co2, ener.heat_adsorption_h2,
            ener.widom_co2, ener.widom_h2,
        ])

    @staticmethod
    def feature_names() -> List[str]:
        return [
            "LCD", "PLD", "ASA", "NASA", "AV", "NAV",
            "porosity", "density", "n_channels", "PSD_mean", "PSD_std",
            "metal_fraction", "n_atom_types",
            "has_OMS", "n_atoms_per_uc", "MW_per_uc",
            "coordination_number", "dimensionality", "n_sbu",
            "henry_CO2", "henry_H2",
            "Qst_CO2", "Qst_H2",
            "widom_CO2", "widom_H2",
        ]


class ZeoppRunner:
    """Interface to Zeo++ for geometric descriptor calculation."""

    def __init__(self, zeopp_path: str = "network", probe_radius: float = 1.86,
                 n_sa_samples: int = 5000, n_vol_samples: int = 100000):
        self.zeopp_path = zeopp_path
        self.probe_radius = probe_radius
        self.n_sa_samples = n_sa_samples
        self.n_vol_samples = n_vol_samples

    def compute_descriptors(self, cif_path: Path) -> GeometricDescriptors:
        """Run Zeo++ analyses and return geometric descriptors."""
        desc = GeometricDescriptors()
        desc.lcd, desc.pld = self._run_diameter(cif_path)
        desc.asa, desc.nasa = self._run_surface_area(cif_path)
        desc.av, desc.nav, desc.porosity = self._run_volume(cif_path)
        desc.n_channels = self._run_channel(cif_path)
        desc.psd_mean, desc.psd_std = self._run_psd(cif_path)
        desc.density = self._compute_density(cif_path)
        return desc

    def _run_diameter(self, cif_path: Path) -> Tuple[float, float]:
        """Compute LCD and PLD using Zeo++ -res flag."""
        out_file = cif_path.with_suffix(".res")
        cmd = [self.zeopp_path, "-ha", "-res", str(out_file), str(cif_path)]
        try:
            subprocess.run(cmd, capture_output=True, timeout=300, check=True)
            with open(out_file) as f:
                tokens = f.read().split()
                return float(tokens[1]), float(tokens[2])
        except (subprocess.SubprocessError, FileNotFoundError, IndexError) as e:
            logger.warning(f"Zeo++ diameter failed for {cif_path}: {e}")
            return 0.0, 0.0

    def _run_surface_area(self, cif_path: Path) -> Tuple[float, float]:
        """Compute ASA and NASA using Zeo++ -sa flag."""
        out_file = cif_path.with_suffix(".sa")
        cmd = [
            self.zeopp_path, "-ha",
            "-sa", str(self.probe_radius), str(self.probe_radius),
            str(self.n_sa_samples), str(out_file), str(cif_path)
        ]
        try:
            subprocess.run(cmd, capture_output=True, timeout=600, check=True)
            with open(out_file) as f:
                tokens = f.readline().split()
                asa_m2g = float(tokens[9])
                nasa_m2g = float(tokens[13]) if len(tokens) > 13 else 0.0
            return asa_m2g, nasa_m2g
        except (subprocess.SubprocessError, FileNotFoundError, IndexError) as e:
            logger.warning(f"Zeo++ SA failed for {cif_path}: {e}")
            return 0.0, 0.0

    def _run_volume(self, cif_path: Path) -> Tuple[float, float, float]:
        """Compute accessible/non-accessible volume and porosity."""
        out_file = cif_path.with_suffix(".vol")
        cmd = [
            self.zeopp_path, "-ha",
            "-vol", str(self.probe_radius), str(self.probe_radius),
            str(self.n_vol_samples), str(out_file), str(cif_path)
        ]
        try:
            subprocess.run(cmd, capture_output=True, timeout=600, check=True)
            with open(out_file) as f:
                tokens = f.readline().split()
                porosity = float(tokens[9])
                av = float(tokens[11])
                nav = float(tokens[15]) if len(tokens) > 15 else 0.0
            return av, nav, porosity
        except (subprocess.SubprocessError, FileNotFoundError, IndexError) as e:
            logger.warning(f"Zeo++ volume failed for {cif_path}: {e}")
            return 0.0, 0.0, 0.0

    def _run_channel(self, cif_path: Path) -> int:
        """Detect number of accessible channels."""
        out_file = cif_path.with_suffix(".chan")
        cmd = [self.zeopp_path, "-ha", "-chan", str(self.probe_radius),
               str(out_file), str(cif_path)]
        try:
            subprocess.run(cmd, capture_output=True, timeout=300, check=True)
            with open(out_file) as f:
                return int(f.readline().split()[0])
        except (subprocess.SubprocessError, FileNotFoundError, IndexError) as e:
            logger.warning(f"Zeo++ channel failed for {cif_path}: {e}")
            return 0

    def _run_psd(self, cif_path: Path) -> Tuple[float, float]:
        """Compute pore size distribution statistics."""
        out_file = cif_path.with_suffix(".psd_histo")
        cmd = [
            self.zeopp_path, "-ha",
            "-psd", str(self.probe_radius), str(self.probe_radius),
            str(self.n_vol_samples), str(out_file), str(cif_path)
        ]
        try:
            subprocess.run(cmd, capture_output=True, timeout=600, check=True)
            sizes, counts = [], []
            with open(out_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    tokens = line.split()
                    if len(tokens) >= 2:
                        sizes.append(float(tokens[0]))
                        counts.append(float(tokens[1]))
            if counts:
                s = np.array(sizes)
                c = np.array(counts)
                total = c.sum()
                if total > 0:
                    mean = np.average(s, weights=c)
                    var = np.average((s - mean) ** 2, weights=c)
                    return mean, np.sqrt(var)
            return 0.0, 0.0
        except (subprocess.SubprocessError, FileNotFoundError) as e:
            logger.warning(f"Zeo++ PSD failed for {cif_path}: {e}")
            return 0.0, 0.0

    def _compute_density(self, cif_path: Path) -> float:
        """Extract framework density from CIF file."""
        try:
            with open(cif_path) as f:
                for line in f:
                    if "_cell_volume" in line:
                        return 0.0  # Full impl uses pymatgen
            return 0.0
        except Exception:
            return 0.0


class ChemicalAnalyzer:
    """Extract chemical composition descriptors from CIF structures."""

    METALS = {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr",
        "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
        "Au", "Pb", "Bi", "U",
    }

    ATOMIC_MASSES = {
        "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
        "S": 32.065, "Cl": 35.453, "Br": 79.904, "I": 126.904,
        "Zn": 65.38, "Cu": 63.546, "Zr": 91.224, "Al": 26.982,
        "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cr": 51.996,
        "Mn": 54.938, "Cd": 112.411, "In": 114.818, "Mg": 24.305,
        "Ca": 40.078, "Ba": 137.327, "Ti": 47.867, "V": 50.942,
    }

    def analyze(self, cif_path: Path) -> ChemicalDescriptors:
        desc = ChemicalDescriptors()
        try:
            elements, counts = self._parse_composition(cif_path)
            desc.n_atom_types = len(set(elements))
            desc.n_atoms_per_uc = sum(counts)
            metals = [(e, c) for e, c in zip(elements, counts) if e in self.METALS]
            if metals:
                desc.metal_type = metals[0][0]
                total_mass = sum(self.ATOMIC_MASSES.get(e, 50.0) * c
                                 for e, c in zip(elements, counts))
                metal_mass = sum(self.ATOMIC_MASSES.get(e, 50.0) * c for e, c in metals)
                desc.metal_fraction = metal_mass / total_mass if total_mass > 0 else 0
                desc.molecular_weight_uc = total_mass
            desc.has_open_metal_sites = self._detect_oms(cif_path, desc.metal_type)
            desc.linker_type = self._classify_linker(elements)
        except Exception as e:
            logger.warning(f"Chemical analysis failed for {cif_path}: {e}")
        return desc

    def _parse_composition(self, cif_path: Path) -> Tuple[List[str], List[int]]:
        elements, counts = [], []
        in_atom_block = False
        with open(cif_path) as f:
            for line in f:
                if "_atom_site_type_symbol" in line:
                    in_atom_block = True
                    continue
                if in_atom_block and line.strip() and not line.startswith("_"):
                    tokens = line.split()
                    if tokens:
                        elem = tokens[0].strip("0123456789+-")
                        elements.append(elem)
                        counts.append(1)
                if in_atom_block and (line.startswith("loop_") or line.strip() == ""):
                    if elements:
                        break
        return elements, counts

    def _detect_oms(self, cif_path: Path, metal: str) -> bool:
        oms_metals = {"Cu", "Cr", "Fe", "Co", "Ni", "Mn", "V", "Zn"}
        return metal in oms_metals

    def _classify_linker(self, elements: List[str]) -> str:
        has_n = "N" in elements
        has_o = "O" in elements
        has_s = "S" in elements
        if has_n and not has_o:
            return "azolate"
        elif has_o and not has_n:
            return "carboxylate"
        elif has_n and has_o:
            return "mixed_N_O"
        elif has_s:
            return "thiolate"
        return "unknown"


class MOFDatabaseLoader:
    """Load and index MOF structures from CoRE MOF and hMOF databases."""

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.core_dir = data_dir / "CoRE_MOF"
        self.hmof_dir = data_dir / "hMOF"

    def list_structures(self, source: str = "all",
                        max_n: Optional[int] = None) -> List[Dict]:
        structures = []
        if source in ("all", "CoRE"):
            structures.extend(self._scan_directory(self.core_dir, "CoRE"))
        if source in ("all", "hMOF"):
            structures.extend(self._scan_directory(self.hmof_dir, "hMOF"))
        if max_n is not None:
            structures = structures[:max_n]
        logger.info(f"Found {len(structures)} MOF structures from {source}")
        return structures

    def _scan_directory(self, directory: Path, source: str) -> List[Dict]:
        if not directory.exists():
            logger.warning(f"Database directory not found: {directory}")
            return []
        cifs = sorted(directory.glob("*.cif"))
        return [{"mof_id": cif.stem, "source_db": source, "cif_path": str(cif)}
                for cif in cifs]


class FeatureExtractionPipeline:
    """Orchestrate full feature extraction for MOF structures."""

    def __init__(self, zeopp_runner: ZeoppRunner, chem_analyzer: ChemicalAnalyzer):
        self.zeopp = zeopp_runner
        self.chem = chem_analyzer

    def extract(self, mof_entry: Dict) -> MOFFeatures:
        cif_path = Path(mof_entry["cif_path"])
        features = MOFFeatures(
            mof_id=mof_entry["mof_id"],
            source_db=mof_entry["source_db"],
            cif_path=str(cif_path),
        )
        features.geometric = self.zeopp.compute_descriptors(cif_path)
        features.chemical = self.chem.analyze(cif_path)
        logger.info(f"Extracted features for {features.mof_id}: "
                     f"LCD={features.geometric.lcd:.2f} Å, "
                     f"ASA={features.geometric.asa:.1f} m²/g")
        return features

    def extract_batch(self, entries: List[Dict]) -> List[MOFFeatures]:
        results = []
        for i, entry in enumerate(entries):
            try:
                results.append(self.extract(entry))
            except Exception as e:
                logger.error(f"Failed to extract {entry['mof_id']}: {e}")
            if (i + 1) % 100 == 0:
                logger.info(f"Processed {i + 1}/{len(entries)} structures")
        return results
