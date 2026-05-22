"""
Module 3: Free Energy Perturbation (FEP) and Metadynamics Comparison

Implements and compares two free energy calculation approaches:
- FEP: Alchemical free energy perturbation for relative binding free energies
- Metadynamics: Enhanced sampling with collective variables for absolute binding free energies

Provides comparative analysis framework and recommendation engine.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from enum import Enum


class FreeEnergyMethod(Enum):
    FEP = "fep"
    METADYNAMICS = "metadynamics"
    TI = "thermodynamic_integration"
    ABFE = "absolute_binding_free_energy"


@dataclass
class AlchemicalLeg:
    """Single leg of an alchemical free energy calculation."""
    lambda_windows: List[float]
    delta_g: float  # kcal/mol
    uncertainty: float
    overlap_matrix: Optional[List[List[float]]] = None
    convergence_time_ns: float = 0.0


@dataclass
class FEPResult:
    """Results from a Free Energy Perturbation calculation."""
    ligand_a: str
    ligand_b: str
    
    # Thermodynamic cycle legs
    complex_leg: Optional[AlchemicalLeg] = None
    solvent_leg: Optional[AlchemicalLeg] = None
    
    # Final result
    ddg_fep: float = 0.0  # ΔΔG (kcal/mol)
    ddg_uncertainty: float = 0.0
    ddg_experimental: Optional[float] = None
    
    # Quality metrics
    hysteresis: float = 0.0
    convergence_ratio: float = 0.0
    
    @property
    def signed_error(self) -> Optional[float]:
        if self.ddg_experimental is not None:
            return self.ddg_fep - self.ddg_experimental
        return None


@dataclass
class MetadynamicsResult:
    """Results from a Metadynamics calculation."""
    ligand: str
    
    # Collective variables
    cv_names: List[str] = field(default_factory=list)
    cv_ranges: List[Tuple[float, float]] = field(default_factory=list)
    
    # Free energy surface
    fes_bins: int = 100
    fes_values: Optional[np.ndarray] = None
    
    # Binding free energy
    dg_binding: float = 0.0  # kcal/mol
    dg_uncertainty: float = 0.0
    dg_experimental: Optional[float] = None
    
    # Convergence
    hills_deposited: int = 0
    total_simulation_time_ns: float = 0.0
    is_converged: bool = False
    
    # Quality metrics
    diffusion_coefficient: float = 0.0
    recrossing_events: int = 0


@dataclass
class MethodComparison:
    """Comparative analysis of FEP vs Metadynamics."""
    fep_results: List[FEPResult]
    metadynamics_results: List[MetadynamicsResult]
    
    # Accuracy comparison
    fep_rmse: float = 0.0
    fep_mae: float = 0.0
    fep_r_squared: float = 0.0
    fep_kendall_tau: float = 0.0
    
    metad_rmse: float = 0.0
    metad_mae: float = 0.0
    metad_r_squared: float = 0.0
    metad_kendall_tau: float = 0.0
    
    # Computational cost
    fep_total_gpu_hours: float = 0.0
    metad_total_gpu_hours: float = 0.0
    
    # Recommendations
    recommended_method: str = ""
    reasoning: str = ""


def generate_lambda_schedule(
    n_windows: int = 12,
    schedule_type: str = "optimal"
) -> List[float]:
    """Generate lambda schedule for alchemical transformations."""
    if schedule_type == "linear":
        return np.linspace(0, 1, n_windows).tolist()
    elif schedule_type == "optimal":
        # Trapezoidal scheme with denser sampling near endpoints
        lambdas = []
        for i in range(n_windows):
            x = i / (n_windows - 1)
            # Sigmoid-like transformation for better overlap
            lam = x ** 2 * (3 - 2 * x) if x <= 0.5 else 1 - (1 - x) ** 2 * (3 - 2 * (1 - x))
            lambdas.append(round(lam, 4))
        return lambdas
    elif schedule_type == "geometric":
        half = n_windows // 2
        lower = (np.geomspace(0.01, 0.5, half)).tolist()
        upper = (1 - np.geomspace(0.01, 0.5, n_windows - half)[::-1]).tolist()
        return [round(x, 4) for x in lower + upper]
    else:
        return np.linspace(0, 1, n_windows).tolist()


def generate_fep_openmm_code(
    ligand_a_sdf: str,
    ligand_b_sdf: str,
    n_lambda: int = 12,
    per_window_ns: float = 5.0
) -> str:
    """Generate OpenMM-based FEP calculation code."""
    lambdas = generate_lambda_schedule(n_lambda)
    return f'''
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmtools import alchemy, states, mcmc
from openmmtools.multistate import MultiStateReporter, MultiStateSampler
from perses.annihilation import NCMCHybridSystemFactory

# Define lambda schedule
lambda_schedule = {lambdas}

# Create alchemical system
alchemical_region = alchemy.AlchemicalRegion(
    alchemical_atoms=ligand_atoms,
    annihilate_electrostatics=True,
    annihilate_sterics=False,
    softcore_alpha=0.5,
    softcore_a=1,
    softcore_b=1,
    softcore_c=6
)

factory = alchemy.AbsoluteAlchemicalFactory()
alchemical_system = factory.create_alchemical_system(
    reference_system, alchemical_region
)

# Setup replica exchange with solute tempering (REST2)
n_replicas = len(lambda_schedule)
thermodynamic_states = []

for lam in lambda_schedule:
    alch_state = alchemy.AlchemicalState.from_system(alchemical_system)
    alch_state.lambda_electrostatics = lam
    alch_state.lambda_sterics = lam
    thermo_state = states.ThermodynamicState(
        alchemical_system, temperature=300*unit.kelvin
    )
    thermodynamic_states.append(thermo_state)

# MCMC move
mcmc_move = mcmc.LangevinDynamicsMove(
    timestep=4.0*unit.femtoseconds,
    collision_rate=1.0/unit.picosecond,
    n_steps=1250,  # 5 ps per iteration
    reassign_velocities=False
)

# Multi-state sampler (MBAR)
sampler = MultiStateSampler(mcmc_moves=mcmc_move, number_of_iterations=1000)
reporter = MultiStateReporter("fep_complex.nc", checkpoint_interval=10)
sampler.create(
    thermodynamic_states=thermodynamic_states,
    sampler_states=sampler_states,
    storage=reporter
)

# Run
sampler.run()

# Analyze with MBAR
from openmmtools.multistate import MultiStateSamplerAnalyzer
analyzer = MultiStateSamplerAnalyzer(reporter)
delta_f, ddelta_f = analyzer.get_free_energy()
print(f"ΔG = {{delta_f[0,-1]:.2f}} ± {{ddelta_f[0,-1]:.2f}} kT")
'''


def generate_metadynamics_code(
    cv_type: str = "funnel",
    height: float = 1.0,  # kJ/mol
    sigma: float = 0.05,  # nm
    pace: int = 500,
    temperature: float = 300.0
) -> str:
    """Generate OpenMM metadynamics code with PLUMED-style CVs."""
    return f'''
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm.app import metadynamics

# Define collective variables
# CV1: Distance between protein binding site COM and ligand COM
cv1 = mm.CustomCentroidBondForce(2, "distance(g1, g2)")
cv1.addGroup(binding_site_atoms)  # protein binding site
cv1.addGroup(ligand_atoms)         # ligand
cv1.addBond([0, 1])

# CV2: Coordination number (number of contacts)
cv2 = mm.CustomNonbondedForce(
    "step(r_cut - r) * (1 - (r/r_cut)^6) / (1 - (r/r_cut)^12)"
)
cv2.addGlobalParameter("r_cut", 0.6)  # 6 Å cutoff

# Funnel restraint (prevents ligand from escaping sideways)
funnel = mm.CustomCentroidBondForce(2, 
    "0.5*k_funnel*(max(0, r_perp - R_cyl))^2; "
    "r_perp = sqrt(dx^2 + dy^2); "
    "dx = x2 - x1 - proj_x; dy = y2 - y1 - proj_y; "
    "proj_x = ((x2-x1)*ax + (y2-y1)*ay + (z2-z1)*az)*ax; "
    "proj_y = ((x2-x1)*ax + (y2-y1)*ay + (z2-z1)*az)*ay"
)

# Well-tempered metadynamics
meta = metadynamics.Metadynamics(
    system,
    [cv1_biasvar, cv2_biasvar],
    temperature={temperature}*unit.kelvin,
    biasFactor=10,  # well-tempered: ΔT = (γ-1)*T
    height={height}*unit.kilojoules_per_mole,
    frequency={pace},
    biasDir="metad_hills",
    saveFrequency={pace * 10}
)

# Run metadynamics
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

n_steps = 500000000  # 1 μs total
meta.step(simulation, n_steps)

# Reconstruct free energy surface
fes = meta.getFreeEnergy()
print("Free energy surface reconstructed")
print(f"Binding ΔG = {{compute_binding_dg(fes):.2f}} kcal/mol")
'''


def simulate_fep_results(
    n_perturbations: int = 10,
    seed: int = 42
) -> List[FEPResult]:
    """Generate realistic simulated FEP results for demonstration."""
    rng = np.random.RandomState(seed)
    results = []
    
    ligand_names = [f"LIG_{chr(65+i)}" for i in range(n_perturbations + 1)]
    
    for i in range(n_perturbations):
        ddg_exp = rng.uniform(-3, 3)
        ddg_calc = ddg_exp + rng.normal(0, 0.8)
        
        n_lambda = 12
        lambdas = generate_lambda_schedule(n_lambda)
        
        results.append(FEPResult(
            ligand_a=ligand_names[i],
            ligand_b=ligand_names[i + 1],
            complex_leg=AlchemicalLeg(
                lambda_windows=lambdas,
                delta_g=float(rng.normal(-15, 5)),
                uncertainty=float(rng.uniform(0.1, 0.4)),
                convergence_time_ns=float(rng.uniform(2, 8))
            ),
            solvent_leg=AlchemicalLeg(
                lambda_windows=lambdas,
                delta_g=float(rng.normal(-10, 3)),
                uncertainty=float(rng.uniform(0.1, 0.3)),
                convergence_time_ns=float(rng.uniform(1, 5))
            ),
            ddg_fep=float(ddg_calc),
            ddg_uncertainty=float(rng.uniform(0.2, 0.6)),
            ddg_experimental=float(ddg_exp),
            hysteresis=float(rng.uniform(0.1, 0.5)),
            convergence_ratio=float(rng.uniform(0.85, 0.99)),
        ))
    
    return results


def simulate_metadynamics_results(
    n_ligands: int = 10,
    seed: int = 42
) -> List[MetadynamicsResult]:
    """Generate realistic simulated metadynamics results."""
    rng = np.random.RandomState(seed)
    results = []
    
    for i in range(n_ligands):
        dg_exp = rng.uniform(-12, -4)
        dg_calc = dg_exp + rng.normal(0, 1.5)
        
        results.append(MetadynamicsResult(
            ligand=f"LIG_{chr(65+i)}",
            cv_names=["protein-ligand distance", "coordination number"],
            cv_ranges=[(0.2, 3.0), (0, 15)],
            dg_binding=float(dg_calc),
            dg_uncertainty=float(rng.uniform(0.5, 1.5)),
            dg_experimental=float(dg_exp),
            hills_deposited=int(rng.randint(50000, 200000)),
            total_simulation_time_ns=float(rng.uniform(500, 2000)),
            is_converged=bool(rng.random() > 0.2),
            diffusion_coefficient=float(rng.uniform(0.5, 2.0)),
            recrossing_events=int(rng.randint(5, 50)),
        ))
    
    return results


def compare_methods(
    fep_results: List[FEPResult],
    metadynamics_results: List[MetadynamicsResult]
) -> MethodComparison:
    """Perform comprehensive comparison of FEP and metadynamics."""
    # FEP statistics
    fep_errors = [r.signed_error for r in fep_results if r.signed_error is not None]
    fep_calc = [r.ddg_fep for r in fep_results if r.ddg_experimental is not None]
    fep_exp = [r.ddg_experimental for r in fep_results if r.ddg_experimental is not None]
    
    # Metadynamics statistics
    metad_errors = [r.dg_binding - r.dg_experimental
                    for r in metadynamics_results if r.dg_experimental is not None]
    metad_calc = [r.dg_binding for r in metadynamics_results if r.dg_experimental is not None]
    metad_exp = [r.dg_experimental for r in metadynamics_results if r.dg_experimental is not None]
    
    def compute_stats(errors, calc, exp):
        rmse = np.sqrt(np.mean(np.array(errors) ** 2))
        mae = np.mean(np.abs(errors))
        if len(calc) > 1:
            corr = np.corrcoef(calc, exp)[0, 1]
            r2 = corr ** 2
        else:
            r2 = 0.0
        # Kendall tau
        tau = _kendall_tau(calc, exp)
        return rmse, mae, r2, tau
    
    fep_rmse, fep_mae, fep_r2, fep_tau = compute_stats(fep_errors, fep_calc, fep_exp)
    metad_rmse, metad_mae, metad_r2, metad_tau = compute_stats(metad_errors, metad_calc, metad_exp)
    
    # Cost estimation
    fep_gpu_hours = sum(
        len(r.complex_leg.lambda_windows) * 5 * 2  # windows × ns × 2 legs
        for r in fep_results if r.complex_leg
    ) * 0.5  # GPU hours per ns
    
    metad_gpu_hours = sum(
        r.total_simulation_time_ns for r in metadynamics_results
    ) * 0.5
    
    # Recommendation
    if fep_rmse < metad_rmse and fep_gpu_hours < metad_gpu_hours:
        rec = "FEP"
        reason = "FEP shows better accuracy with lower computational cost"
    elif metad_rmse < fep_rmse:
        rec = "Metadynamics"
        reason = "Metadynamics shows better accuracy despite higher cost"
    else:
        rec = "FEP"
        reason = "FEP preferred for relative rankings; use metadynamics for absolute ΔG"
    
    return MethodComparison(
        fep_results=fep_results,
        metadynamics_results=metadynamics_results,
        fep_rmse=float(fep_rmse),
        fep_mae=float(fep_mae),
        fep_r_squared=float(fep_r2),
        fep_kendall_tau=float(fep_tau),
        metad_rmse=float(metad_rmse),
        metad_mae=float(metad_mae),
        metad_r_squared=float(metad_r2),
        metad_kendall_tau=float(metad_tau),
        fep_total_gpu_hours=float(fep_gpu_hours),
        metad_total_gpu_hours=float(metad_gpu_hours),
        recommended_method=rec,
        reasoning=reason,
    )


def _kendall_tau(x, y):
    """Compute Kendall's tau rank correlation."""
    n = len(x)
    if n < 2:
        return 0.0
    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            if dx * dy > 0:
                concordant += 1
            elif dx * dy < 0:
                discordant += 1
    denom = n * (n - 1) / 2
    return (concordant - discordant) / denom if denom > 0 else 0.0
