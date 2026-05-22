"""
Pipeline configuration for MOF high-throughput screening.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


@dataclass
class GCMCConfig:
    """RASPA GCMC simulation parameters."""
    temperature: float = 298.0          # K
    pressure_points: List[float] = field(default_factory=lambda: [
        0.0004, 0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0  # bar
    ])
    n_cycles_init: int = 10000
    n_cycles_prod: int = 50000
    force_field: str = "GenericMOFs"
    co2_model: str = "TraPPE"
    h2_model: str = "Darkrim-Levesque"
    cutoff: float = 12.8               # Å
    charge_method: str = "EQeq"


@dataclass
class ZeoppConfig:
    """Zeo++ geometric analysis parameters."""
    probe_radius_n2: float = 1.86       # Å (N2 probe)
    probe_radius_co2: float = 1.72      # Å (CO2 kinetic diameter / 2)
    n_samples_sa: int = 5000
    n_samples_vol: int = 100000
    ha: bool = True
    channel_probe: float = 1.86


@dataclass
class MLConfig:
    """ML model configuration."""
    model_type: str = "gradient_boosting"
    n_estimators: int = 500
    max_depth: int = 8
    learning_rate: float = 0.05
    test_size: float = 0.2
    cv_folds: int = 5
    random_state: int = 42
    use_geometric: bool = True
    use_chemical: bool = True
    use_topological: bool = True
    use_energy: bool = True


@dataclass
class StabilityConfig:
    """Water stability and synthesizability filter."""
    water_stability_threshold: float = 0.5
    synthesizability_threshold: float = 0.3
    metal_node_blacklist: List[str] = field(default_factory=lambda: [
        "Cr2+", "Fe2+", "Mn2+"
    ])
    linker_stability_classes: List[str] = field(default_factory=lambda: [
        "azolate", "carboxylate_zr", "phosphonate"
    ])


@dataclass
class DACConfig:
    """DAC-specific ranking criteria."""
    co2_concentration_ppm: float = 420.0
    co2_partial_pressure_bar: float = 0.000420
    desorption_temperature: float = 373.0    # K (100°C)
    min_working_capacity: float = 1.0        # mmol/g
    min_selectivity_co2_n2: float = 50.0
    max_regeneration_energy: float = 60.0    # kJ/mol
    humidity_tolerance: bool = True


@dataclass
class PipelineConfig:
    """Main pipeline configuration."""
    workspace: Path = Path(".")
    data_dir: Path = Path("data")
    results_dir: Path = Path("results")
    figures_dir: Path = Path("figures")
    logs_dir: Path = Path("logs")
    core_mof_version: str = "2019-ASR"
    hmof_subset: str = "all"
    max_structures: Optional[int] = None
    gcmc: GCMCConfig = field(default_factory=GCMCConfig)
    zeopp: ZeoppConfig = field(default_factory=ZeoppConfig)
    ml: MLConfig = field(default_factory=MLConfig)
    stability: StabilityConfig = field(default_factory=StabilityConfig)
    dac: DACConfig = field(default_factory=DACConfig)
    n_workers: int = 8
    batch_size: int = 100

    def __post_init__(self):
        for d in [self.data_dir, self.results_dir, self.figures_dir, self.logs_dir]:
            d.mkdir(parents=True, exist_ok=True)
