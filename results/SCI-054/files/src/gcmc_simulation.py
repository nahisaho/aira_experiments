"""
Module 2: Grand Canonical Monte Carlo (GCMC) adsorption simulation via RASPA.

Generates RASPA input files, executes simulations, and parses output
for CO2 and H2 adsorption isotherms.
"""
import json
import logging
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class AdsorptionPoint:
    """Single adsorption data point."""
    pressure: float            # bar
    loading_abs: float         # mol/kg (absolute)
    loading_excess: float      # mol/kg (excess)
    loading_mmol_g: float      # mmol/g
    enthalpy: float            # kJ/mol
    error_loading: float
    error_enthalpy: float


@dataclass
class AdsorptionIsotherm:
    """Complete adsorption isotherm for one gas in one MOF."""
    mof_id: str
    gas: str
    temperature: float
    points: List[AdsorptionPoint] = field(default_factory=list)
    henry_coefficient: float = 0.0    # mol/kg/Pa
    heat_of_adsorption: float = 0.0   # kJ/mol at zero loading

    def pressures(self) -> np.ndarray:
        return np.array([p.pressure for p in self.points])

    def loadings(self) -> np.ndarray:
        return np.array([p.loading_mmol_g for p in self.points])

    def to_dict(self) -> Dict:
        return {
            "mof_id": self.mof_id, "gas": self.gas,
            "temperature": self.temperature,
            "henry_coefficient": self.henry_coefficient,
            "Qst_zero_loading": self.heat_of_adsorption,
            "isotherm": [
                {"P_bar": p.pressure, "loading_mmol_g": p.loading_mmol_g,
                 "loading_mol_kg": p.loading_abs, "enthalpy_kJ_mol": p.enthalpy}
                for p in self.points
            ],
        }


class RASPAInputGenerator:
    """Generate RASPA simulation input files."""

    GCMC_TEMPLATE = """\
SimulationType                MonteCarlo
NumberOfCycles                {n_cycles_prod}
NumberOfInitializationCycles  {n_cycles_init}
PrintEvery                    {print_every}
RestartFile                   no

Forcefield                    {force_field}
CutOffVDW                    {cutoff}
ChargeMethod                 Ewald
EwaldPrecision               1e-6

Framework 0
FrameworkName                {framework_name}
UnitCells                    {unit_cells}
ExternalTemperature          {temperature}
ExternalPressure             {pressure}

Component 0 MoleculeName             {molecule}
             MoleculeDefinition       {molecule_def}
             TranslationProbability   0.5
             RotationProbability       0.5
             ReinsertionProbability    0.5
             SwapProbability           1.0
             CreateNumberOfMolecules   0
"""

    WIDOM_TEMPLATE = """\
SimulationType                MonteCarlo
NumberOfCycles                {n_cycles}
PrintEvery                    {print_every}

Forcefield                    {force_field}
CutOffVDW                    {cutoff}

Framework 0
FrameworkName                {framework_name}
UnitCells                    {unit_cells}
ExternalTemperature          {temperature}

Component 0 MoleculeName             {molecule}
             MoleculeDefinition       {molecule_def}
             WidomProbability          1.0
             CreateNumberOfMolecules   0
"""

    def __init__(self, raspa_dir: str = ".", force_field: str = "GenericMOFs",
                 cutoff: float = 12.8):
        self.raspa_dir = Path(raspa_dir)
        self.force_field = force_field
        self.cutoff = cutoff

    def generate_gcmc_input(self, framework_name: str, molecule: str,
                             temperature: float, pressure: float,
                             n_cycles_init: int = 10000,
                             n_cycles_prod: int = 50000,
                             unit_cells: str = "2 2 2",
                             output_dir: Optional[Path] = None) -> Path:
        if output_dir is None:
            output_dir = self.raspa_dir / "simulations" / framework_name / molecule
        output_dir.mkdir(parents=True, exist_ok=True)
        sim_dir = output_dir / f"P_{pressure:.6f}"
        sim_dir.mkdir(parents=True, exist_ok=True)

        molecule_def = self._get_molecule_definition(molecule)
        content = self.GCMC_TEMPLATE.format(
            n_cycles_prod=n_cycles_prod, n_cycles_init=n_cycles_init,
            print_every=max(n_cycles_prod // 10, 1),
            force_field=self.force_field, cutoff=self.cutoff,
            framework_name=framework_name, unit_cells=unit_cells,
            temperature=temperature, pressure=pressure * 1e5,
            molecule=molecule, molecule_def=molecule_def,
        )
        (sim_dir / "simulation.input").write_text(content)
        return sim_dir

    def generate_widom_input(self, framework_name: str, molecule: str,
                              temperature: float, n_cycles: int = 50000,
                              unit_cells: str = "2 2 2",
                              output_dir: Optional[Path] = None) -> Path:
        if output_dir is None:
            output_dir = self.raspa_dir / "simulations" / framework_name / f"{molecule}_widom"
        output_dir.mkdir(parents=True, exist_ok=True)
        molecule_def = self._get_molecule_definition(molecule)
        content = self.WIDOM_TEMPLATE.format(
            n_cycles=n_cycles, print_every=max(n_cycles // 10, 1),
            force_field=self.force_field, cutoff=self.cutoff,
            framework_name=framework_name, unit_cells=unit_cells,
            temperature=temperature, molecule=molecule,
            molecule_def=molecule_def,
        )
        (output_dir / "simulation.input").write_text(content)
        return output_dir

    def _get_molecule_definition(self, molecule: str) -> str:
        return {"CO2": "TraPPE", "H2": "Darkrim-Levesque", "N2": "TraPPE",
                "H2O": "TIP4P", "CH4": "TraPPE"}.get(molecule, "ExampleDefinition")

    def determine_unit_cells(self, cell_lengths: Tuple[float, float, float],
                              cutoff: float = 12.8) -> str:
        return " ".join(str(max(1, int(np.ceil(2 * cutoff / l))))
                        for l in cell_lengths)


class RASPAOutputParser:
    """Parse RASPA simulation output files."""

    def parse_gcmc_output(self, sim_dir: Path) -> Optional[AdsorptionPoint]:
        output_dir = sim_dir / "Output" / "System_0"
        if not output_dir.exists():
            output_files = list(sim_dir.rglob("*.data"))
            if not output_files:
                return None
            output_file = output_files[0]
        else:
            output_files = list(output_dir.glob("*.data"))
            if not output_files:
                return None
            output_file = output_files[0]
        return self._parse_data_file(output_file)

    def _parse_data_file(self, filepath: Path) -> Optional[AdsorptionPoint]:
        loading_abs = loading_excess = enthalpy = 0.0
        error_loading = error_enthalpy = pressure = 0.0
        try:
            content = filepath.read_text()
            for line in content.split("\n"):
                if "Average loading absolute [mol/kg]" in line:
                    tokens = line.split()
                    loading_abs = float(tokens[-4])
                    error_loading = float(tokens[-2])
                elif "Average loading excess [mol/kg]" in line:
                    loading_excess = float(line.split()[-4])
                elif "Enthalpy of adsorption" in line and "[K]" not in line:
                    tokens = line.split()
                    for i, t in enumerate(tokens):
                        if t == "[kJ/mol]":
                            enthalpy = float(tokens[i - 1])
                            if i + 2 < len(tokens):
                                error_enthalpy = float(tokens[i + 2])
                            break
                elif "External pressure:" in line:
                    pressure = float(line.split()[-2]) / 1e5
            return AdsorptionPoint(
                pressure=pressure, loading_abs=loading_abs,
                loading_excess=loading_excess, loading_mmol_g=loading_abs,
                enthalpy=enthalpy, error_loading=error_loading,
                error_enthalpy=error_enthalpy,
            )
        except Exception as e:
            logger.error(f"Failed to parse {filepath}: {e}")
            return None

    def parse_widom_output(self, sim_dir: Path) -> Tuple[float, float]:
        output_files = list(sim_dir.rglob("*.data"))
        if not output_files:
            return 0.0, 0.0
        henry = qst = 0.0
        try:
            for line in output_files[0].read_text().split("\n"):
                if "Average Henry coefficient" in line and "[mol/kg/Pa]" in line:
                    henry = float(line.split()[-4])
                elif "Average  <U_gh>_1-<U_h>_0" in line and "[kJ/mol]" in line:
                    qst = abs(float(line.split()[-4]))
        except Exception as e:
            logger.error(f"Widom parse failed: {e}")
        return henry, qst


class GCMCSimulator:
    """Full GCMC simulation workflow manager."""

    def __init__(self, raspa_path: str = "simulate",
                 input_gen: Optional[RASPAInputGenerator] = None,
                 parser: Optional[RASPAOutputParser] = None):
        self.raspa_path = raspa_path
        self.input_gen = input_gen or RASPAInputGenerator()
        self.parser = parser or RASPAOutputParser()

    def run_isotherm(self, framework_name: str, gas: str, temperature: float,
                      pressures: List[float], n_cycles_init: int = 10000,
                      n_cycles_prod: int = 50000) -> AdsorptionIsotherm:
        isotherm = AdsorptionIsotherm(mof_id=framework_name, gas=gas,
                                       temperature=temperature)
        for pressure in pressures:
            sim_dir = self.input_gen.generate_gcmc_input(
                framework_name=framework_name, molecule=gas,
                temperature=temperature, pressure=pressure,
                n_cycles_init=n_cycles_init, n_cycles_prod=n_cycles_prod,
            )
            if self._execute_raspa(sim_dir):
                point = self.parser.parse_gcmc_output(sim_dir)
                if point:
                    point.pressure = pressure
                    isotherm.points.append(point)

        widom_dir = self.input_gen.generate_widom_input(
            framework_name=framework_name, molecule=gas, temperature=temperature)
        if self._execute_raspa(widom_dir):
            henry, qst = self.parser.parse_widom_output(widom_dir)
            isotherm.henry_coefficient = henry
            isotherm.heat_of_adsorption = qst
        return isotherm

    def _execute_raspa(self, sim_dir: Path) -> bool:
        try:
            result = subprocess.run([self.raspa_path], cwd=str(sim_dir),
                                     capture_output=True, timeout=3600)
            if result.returncode != 0:
                logger.warning(f"RASPA failed in {sim_dir}")
                return False
            return True
        except (subprocess.SubprocessError, FileNotFoundError) as e:
            logger.warning(f"RASPA execution error: {e}")
            return False

    def compute_selectivity(self, isotherm_co2: AdsorptionIsotherm,
                             isotherm_n2: AdsorptionIsotherm,
                             y_co2: float = 0.0004) -> float:
        if isotherm_n2.henry_coefficient > 0:
            return isotherm_co2.henry_coefficient / isotherm_n2.henry_coefficient
        return float("inf")

    def compute_working_capacity(self, isotherm: AdsorptionIsotherm,
                                  p_ads: float, p_des: float) -> float:
        pressures = isotherm.pressures()
        loadings = isotherm.loadings()
        if len(pressures) < 2:
            return 0.0
        return max(0.0, np.interp(p_ads, pressures, loadings) -
                   np.interp(p_des, pressures, loadings))
