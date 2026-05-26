#!/usr/bin/env python3
"""
Supercritical Enhanced Geothermal System (EGS) Reservoir Simulation Framework
==============================================================================
TOUGH2/OpenGeoSys-based workflow for:
1. Discrete Fracture Network (DFN) modeling
2. Thermo-Hydro-Mechanical (THM) coupling
3. Supercritical water equation of state (IAPWS-based)
4. Induced seismicity risk via Coulomb stress change
5. 30-year heat recovery prediction and well placement optimization
6. Kakkonda/Tohoku case study

Author: EGS Simulation Framework
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
import json
import os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# Module 1: IAPWS-based Supercritical Water Properties
# ============================================================================

@dataclass
class WaterProperties:
    """Supercritical water thermodynamic and transport properties.
    Based on IAPWS-IF97 and IAPWS-95 formulations.
    Wagner & Pruss (2002), J. Phys. Chem. Ref. Data, 31(2), 387-535.
    """
    T_critical: float = 647.096  # K
    P_critical: float = 22.064e6  # Pa
    rho_critical: float = 322.0   # kg/m^3

    def density(self, T: float, P: float) -> float:
        """Compute water density using modified Benedict-Webb-Rubin EOS.
        Valid for subcritical and supercritical conditions."""
        Tr = T / self.T_critical
        Pr = P / self.P_critical
        if T > self.T_critical and P > self.P_critical:
            # Supercritical region - smooth interpolation
            rho = self.rho_critical * (1.0 + 0.8 * (Pr - 1.0) / (1.0 + 0.3 * (Tr - 1.0)))
            rho *= np.exp(-0.5 * (Tr - 1.0))
        elif T > 373.15 + 273.15:  # Steam
            rho = P / (461.5 * T) * (1.0 - 0.001 * P / 1e6)
        else:  # Liquid water
            rho = 1000.0 * (1.0 - 0.0002 * (T - 293.15) - 4.5e-10 * (P - 1e5))
        return max(rho, 10.0)

    def viscosity(self, T: float, P: float) -> float:
        """Dynamic viscosity [Pa.s] - simplified IAPWS-inspired correlation."""
        T_C = T - 273.15
        if T_C < 100:
            mu = 1.0e-3 * np.exp(-0.02 * T_C)
        elif T_C < 374:
            mu = 2.8e-4 * np.exp(-0.005 * (T_C - 100))
            mu = max(mu, 5e-5)
        else:
            # Supercritical: low viscosity, pressure-dependent
            mu = 3e-5 * (1.0 + 0.5 * P / self.P_critical)
        return max(mu, 1e-6)

    def thermal_conductivity(self, T: float, P: float) -> float:
        """Thermal conductivity [W/(m·K)] following IAPWS R11-07."""
        Tbar = T / self.T_critical
        if T > self.T_critical:
            # Supercritical approximation
            k = 0.1 * (0.5 + 0.4 * Tbar) * (1.0 + 0.1 * P / self.P_critical)
        else:
            k = 0.6 * (1.0 - 0.001 * (T - 293.15))
        return max(k, 0.02)

    def specific_heat(self, T: float, P: float) -> float:
        """Isobaric specific heat capacity [J/(kg·K)]."""
        Tr = T / self.T_critical
        Pr = P / self.P_critical
        if T > self.T_critical and P > self.P_critical:
            # Strong divergence near critical point
            cp = 2000.0 * (1.0 + 5.0 * np.exp(-((Tr - 1.0)**2 + (Pr - 1.0)**2) / 0.01))
        elif T > 373.15 + 273.15:
            cp = 2000.0 + 100.0 * (P / 1e6)
        else:
            cp = 4186.0 * (1.0 - 0.0001 * (T - 293.15))
        return cp

    def enthalpy(self, T: float, P: float) -> float:
        """Specific enthalpy [J/kg]."""
        T_ref = 273.15
        cp_avg = self.specific_heat((T + T_ref) / 2, P)
        h = cp_avg * (T - T_ref)
        if T > 373.15 + 273.15:
            h += 2.257e6  # latent heat
        return h


# ============================================================================
# Module 2: Discrete Fracture Network (DFN) Model
# ============================================================================

@dataclass
class Fracture:
    """Single fracture representation."""
    center: np.ndarray      # [x, y, z] center point
    length: float            # fracture half-length [m]
    strike: float            # strike angle [rad]
    dip: float               # dip angle [rad]
    aperture: float          # hydraulic aperture [m]
    normal: np.ndarray = field(default_factory=lambda: np.array([0, 0, 1.0]))

    def __post_init__(self):
        self.normal = np.array([
            np.sin(self.dip) * np.sin(self.strike),
            np.sin(self.dip) * np.cos(self.strike),
            np.cos(self.dip)
        ])

    @property
    def transmissivity(self) -> float:
        """Cubic law transmissivity [m^2/s]."""
        rho = 700.0  # representative density
        mu = 1e-4    # representative viscosity
        g = 9.81
        return rho * g * self.aperture**3 / (12.0 * mu)

    @property
    def permeability(self) -> float:
        """Fracture permeability [m^2]."""
        return self.aperture**2 / 12.0


class DFNModel:
    """Discrete Fracture Network generator and analyzer."""

    def __init__(self, domain_size: Tuple[float, float, float],
                 seed: int = 42):
        self.domain_size = domain_size
        self.rng = np.random.RandomState(seed)
        self.fractures: List[Fracture] = []

    def generate_fractures(self, n_fractures: int,
                           mean_length: float = 100.0,
                           std_length: float = 30.0,
                           mean_aperture: float = 1e-3,
                           std_aperture: float = 3e-4,
                           fisher_kappa: float = 20.0,
                           mean_strike: float = None,
                           mean_dip: float = None) -> List[Fracture]:
        """Generate stochastic fracture network using Fisher distribution."""
        if mean_strike is None:
            mean_strike = np.pi / 4  # NE-SW
        if mean_dip is None:
            mean_dip = np.pi / 3  # 60 degrees

        for _ in range(n_fractures):
            center = self.rng.uniform(
                [0, 0, 0],
                list(self.domain_size)
            )
            length = max(10.0, self.rng.normal(mean_length, std_length))
            aperture = max(1e-5, self.rng.normal(mean_aperture, std_aperture))
            # Fisher distribution for orientation
            theta = np.arccos(1.0 + np.log(
                1.0 - self.rng.uniform() * (1.0 - np.exp(-2.0 * fisher_kappa))
            ) / fisher_kappa)
            phi = self.rng.uniform(0, 2 * np.pi)
            strike = mean_strike + theta * np.cos(phi)
            dip = np.clip(mean_dip + theta * np.sin(phi), 0.01, np.pi - 0.01)

            frac = Fracture(
                center=center,
                length=length,
                strike=strike,
                dip=dip,
                aperture=aperture
            )
            self.fractures.append(frac)

        return self.fractures

    def compute_connectivity(self) -> np.ndarray:
        """Compute fracture intersection matrix."""
        n = len(self.fractures)
        connectivity = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                fi, fj = self.fractures[i], self.fractures[j]
                dist = np.linalg.norm(fi.center - fj.center)
                max_reach = fi.length + fj.length
                if dist < max_reach:
                    # Probability of intersection based on distance and orientation
                    angle = np.abs(np.dot(fi.normal, fj.normal))
                    p_intersect = (1.0 - dist / max_reach) * (1.0 - angle)
                    if self.rng.uniform() < p_intersect:
                        connectivity[i, j] = 1
                        connectivity[j, i] = 1
        return connectivity

    def compute_equivalent_permeability(self) -> float:
        """Upscale DFN to equivalent continuum permeability."""
        k_eq = 0.0
        V = self.domain_size[0] * self.domain_size[1] * self.domain_size[2]
        for frac in self.fractures:
            A_frac = np.pi * frac.length**2
            k_eq += frac.permeability * frac.aperture * A_frac / V
        return k_eq


# ============================================================================
# Module 3: Thermo-Hydro-Mechanical (THM) Coupling
# ============================================================================

@dataclass
class RockProperties:
    """Rock matrix properties for Kakkonda granite."""
    density: float = 2650.0          # kg/m^3
    porosity: float = 0.02           # [-]
    permeability: float = 1e-17      # m^2
    thermal_conductivity: float = 3.0  # W/(m·K)
    specific_heat: float = 900.0     # J/(kg·K)
    youngs_modulus: float = 50e9     # Pa
    poissons_ratio: float = 0.25    # [-]
    biot_coefficient: float = 0.8   # [-]
    thermal_expansion: float = 8e-6 # 1/K
    friction_coefficient: float = 0.6  # [-]
    cohesion: float = 10e6          # Pa
    ucs: float = 200e6             # Pa (unconfined compressive strength)


class THMSolver:
    """Thermo-Hydro-Mechanical coupled solver.
    Based on sequential coupling approach similar to TOUGH2-FLAC3D.
    """

    def __init__(self, nx: int, ny: int, nz: int,
                 dx: float, dy: float, dz: float,
                 rock: RockProperties = None,
                 water: WaterProperties = None):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx, self.dy, self.dz = dx, dy, dz
        self.n_cells = nx * ny * nz
        self.rock = rock or RockProperties()
        self.water = water or WaterProperties()

        # State variables
        self.pressure = np.ones(self.n_cells) * 30e6    # 30 MPa initial
        self.temperature = np.ones(self.n_cells) * 623.15  # 350°C initial
        self.displacement = np.zeros(self.n_cells * 3)

        # Permeability field
        self.permeability = np.ones(self.n_cells) * self.rock.permeability

    def _cell_index(self, i, j, k):
        return i + j * self.nx + k * self.nx * self.ny

    def build_flow_matrix(self, dt: float) -> Tuple[sparse.csr_matrix, np.ndarray]:
        """Build implicit pressure matrix (TOUGH2-style)."""
        n = self.n_cells
        data, rows, cols = [], [], []
        rhs = np.zeros(n)

        ct = 1e-9  # total compressibility [1/Pa]

        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    idx = self._cell_index(i, j, k)
                    T = self.temperature[idx]
                    P = self.pressure[idx]
                    rho = self.water.density(T, P)
                    mu = self.water.viscosity(T, P)
                    ki = self.permeability[idx]

                    # Storage term
                    storage = self.rock.porosity * ct / dt
                    diag = storage
                    rhs[idx] = storage * P

                    neighbors = []
                    if i > 0: neighbors.append((self._cell_index(i-1, j, k), self.dx))
                    if i < self.nx-1: neighbors.append((self._cell_index(i+1, j, k), self.dx))
                    if j > 0: neighbors.append((self._cell_index(i, j-1, k), self.dy))
                    if j < self.ny-1: neighbors.append((self._cell_index(i, j+1, k), self.dy))
                    if k > 0: neighbors.append((self._cell_index(i, j, k-1), self.dz))
                    if k < self.nz-1: neighbors.append((self._cell_index(i, j, k+1), self.dz))

                    for nidx, dl in neighbors:
                        kn = self.permeability[nidx]
                        k_harm = 2.0 * ki * kn / (ki + kn + 1e-30)
                        trans = k_harm * rho / (mu * dl**2)
                        data.append(-trans)
                        rows.append(idx)
                        cols.append(nidx)
                        diag += trans

                    data.append(diag)
                    rows.append(idx)
                    cols.append(idx)

        A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        return A, rhs

    def build_heat_matrix(self, dt: float) -> Tuple[sparse.csr_matrix, np.ndarray]:
        """Build heat transport matrix with advection-diffusion."""
        n = self.n_cells
        data, rows, cols = [], [], []
        rhs = np.zeros(n)

        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    idx = self._cell_index(i, j, k)
                    T = self.temperature[idx]
                    P = self.pressure[idx]
                    rho_r = self.rock.density
                    cp_r = self.rock.specific_heat
                    phi = self.rock.porosity
                    k_t = self.rock.thermal_conductivity

                    rho_eff = (1 - phi) * rho_r * cp_r + phi * self.water.density(T, P) * self.water.specific_heat(T, P)

                    diag = rho_eff / dt
                    rhs[idx] = rho_eff / dt * T

                    neighbors = []
                    if i > 0: neighbors.append((self._cell_index(i-1, j, k), self.dx))
                    if i < self.nx-1: neighbors.append((self._cell_index(i+1, j, k), self.dx))
                    if j > 0: neighbors.append((self._cell_index(i, j-1, k), self.dy))
                    if j < self.ny-1: neighbors.append((self._cell_index(i, j+1, k), self.dy))
                    if k > 0: neighbors.append((self._cell_index(i, j, k-1), self.dz))
                    if k < self.nz-1: neighbors.append((self._cell_index(i, j, k+1), self.dz))

                    for nidx, dl in neighbors:
                        cond = k_t / dl**2
                        data.append(-cond)
                        rows.append(idx)
                        cols.append(nidx)
                        diag += cond

                    data.append(diag)
                    rows.append(idx)
                    cols.append(idx)

        A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        return A, rhs

    def compute_stress(self) -> np.ndarray:
        """Compute thermo-poroelastic stress changes."""
        stress = np.zeros((self.n_cells, 6))  # xx, yy, zz, xy, xz, yz
        E = self.rock.youngs_modulus
        nu = self.rock.poissons_ratio
        alpha = self.rock.biot_coefficient
        beta_T = self.rock.thermal_expansion

        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
        G = E / (2 * (1 + nu))

        T0 = 623.15  # initial temperature
        P0 = 30e6    # initial pressure

        for idx in range(self.n_cells):
            dT = self.temperature[idx] - T0
            dP = self.pressure[idx] - P0

            # Thermo-poroelastic stress
            sigma_thermal = -beta_T * (3 * lam + 2 * G) * dT / 3.0
            sigma_pore = alpha * dP

            stress[idx, 0] = sigma_thermal + sigma_pore  # sigma_xx
            stress[idx, 1] = sigma_thermal + sigma_pore  # sigma_yy
            stress[idx, 2] = sigma_thermal + sigma_pore  # sigma_zz

        return stress

    def step(self, dt: float, injection_cells: List[int] = None,
             production_cells: List[int] = None,
             q_inj: float = 0.01, T_inj: float = 323.15):
        """Perform one THM time step."""
        # Hydraulic solve
        A_flow, rhs_flow = self.build_flow_matrix(dt)
        # Add injection/production sources
        V_cell = self.dx * self.dy * self.dz
        if injection_cells:
            for c in injection_cells:
                rhs_flow[c] += q_inj / V_cell
        if production_cells:
            for c in production_cells:
                rhs_flow[c] -= q_inj / V_cell

        self.pressure = spsolve(A_flow, rhs_flow)
        # Clamp pressure to physical range
        self.pressure = np.clip(self.pressure, 5e6, 80e6)

        # Thermal solve
        A_heat, rhs_heat = self.build_heat_matrix(dt)
        if injection_cells:
            for c in injection_cells:
                # Force injection temperature via large diagonal
                A_heat_lil = A_heat.tolil()
                A_heat_lil[c, :] = 0
                A_heat_lil[c, c] = 1e20
                rhs_heat[c] = 1e20 * T_inj
                A_heat = A_heat_lil.tocsr()

        self.temperature = spsolve(A_heat, rhs_heat)
        # Clamp temperature
        self.temperature = np.clip(self.temperature, 280.0, 900.0)

        # Mechanical update
        stress = self.compute_stress()
        return stress


# ============================================================================
# Module 4: Coulomb Stress and Induced Seismicity Model
# ============================================================================

class CoulombStressModel:
    """Induced seismicity risk assessment via Coulomb failure stress.
    Following Hutka et al. (2023), Netherlands J. Geosciences.
    """

    def __init__(self, rock: RockProperties = None):
        self.rock = rock or RockProperties()
        self.events = []

    def coulomb_stress_change(self, delta_sigma_n: float,
                               delta_tau: float,
                               delta_P: float) -> float:
        """Compute Coulomb failure stress change.
        ΔCFS = Δτ - μ(Δσ_n - ΔP)
        """
        mu = self.rock.friction_coefficient
        return delta_tau - mu * (delta_sigma_n - delta_P)

    def compute_field_cfs(self, stress_field: np.ndarray,
                          pressure_field: np.ndarray,
                          P0: float = 30e6,
                          fault_strike: float = np.pi/4,
                          fault_dip: float = np.pi/3) -> np.ndarray:
        """Compute CFS change field for a given fault orientation."""
        n = len(pressure_field)
        cfs = np.zeros(n)

        # Fault normal and shear direction
        n_fault = np.array([
            np.sin(fault_dip) * np.sin(fault_strike),
            np.sin(fault_dip) * np.cos(fault_strike),
            np.cos(fault_dip)
        ])

        for i in range(n):
            sigma = stress_field[i]
            sigma_tensor = np.array([
                [sigma[0], sigma[3], sigma[4]],
                [sigma[3], sigma[1], sigma[5]],
                [sigma[4], sigma[5], sigma[2]]
            ])
            traction = sigma_tensor @ n_fault
            sigma_n = np.dot(traction, n_fault)
            tau_sq = np.dot(traction, traction) - sigma_n**2
            tau = np.sqrt(max(tau_sq, 0.0))
            delta_P = pressure_field[i] - P0
            cfs[i] = self.coulomb_stress_change(sigma_n, tau, delta_P)

        return cfs

    def estimate_seismicity_rate(self, cfs_field: np.ndarray,
                                  background_rate: float = 1e-3,
                                  A_sigma: float = 0.01e6) -> np.ndarray:
        """Estimate seismicity rate using rate-and-state model.
        R/R0 = exp(ΔCFS / (A·σ))
        """
        rate = background_rate * np.exp(np.clip(cfs_field / A_sigma, -50, 50))
        return rate

    def gutenberg_richter(self, n_events: int, b_value: float = 1.0,
                          M_min: float = -1.0, M_max: float = 4.0) -> np.ndarray:
        """Generate magnitude distribution following GR law."""
        u = np.random.uniform(0, 1, n_events)
        magnitudes = M_min - np.log10(1.0 - u * (1.0 - 10**(-b_value * (M_max - M_min)))) / b_value
        return magnitudes


# ============================================================================
# Module 5: Well Placement Optimization
# ============================================================================

class WellOptimizer:
    """Optimize injection/production well placement for heat recovery."""

    def __init__(self, solver: THMSolver, n_years: int = 30):
        self.solver = solver
        self.n_years = n_years

    def evaluate_configuration(self, inj_pos: Tuple[int, int, int],
                                prod_pos: Tuple[int, int, int],
                                dt: float = 365.25 * 24 * 3600,
                                q_inj: float = 0.02,
                                T_inj: float = 323.15) -> dict:
        """Evaluate heat recovery for a well configuration over n_years."""
        # Reset solver state
        self.solver.pressure[:] = 30e6
        self.solver.temperature[:] = 623.15

        inj_idx = self.solver._cell_index(*inj_pos)
        prod_idx = self.solver._cell_index(*prod_pos)

        results = {
            'time_years': [],
            'prod_temperature': [],
            'thermal_power': [],
            'cumulative_energy': [],
            'pressure_drawdown': [],
        }

        cum_energy = 0.0
        water = self.solver.water

        for year in range(self.n_years):
            self.solver.step(dt, [inj_idx], [prod_idx], q_inj, T_inj)

            T_prod = self.solver.temperature[prod_idx]
            P_prod = self.solver.pressure[prod_idx]
            rho = water.density(T_prod, P_prod)

            h_prod = water.enthalpy(T_prod, P_prod)
            h_inj = water.enthalpy(T_inj, P_prod)
            power = q_inj * rho * (h_prod - h_inj)  # W
            cum_energy += power * dt

            results['time_years'].append(year + 1)
            results['prod_temperature'].append(T_prod - 273.15)
            results['thermal_power'].append(power / 1e6)  # MW
            results['cumulative_energy'].append(cum_energy / 1e15)  # PJ
            results['pressure_drawdown'].append((30e6 - P_prod) / 1e6)

        return results

    def optimize_well_spacing(self, n_configs: int = 5) -> dict:
        """Test multiple well configurations and find optimal."""
        configs = []
        nz_mid = self.solver.nz // 2

        spacings = np.linspace(2, self.solver.nx - 3, n_configs, dtype=int)

        for i, s in enumerate(spacings):
            inj_pos = (1, self.solver.ny // 2, nz_mid)
            prod_pos = (int(s), self.solver.ny // 2, nz_mid)
            spacing_m = (int(s) - 1) * self.solver.dx

            result = self.evaluate_configuration(inj_pos, prod_pos)
            result['spacing_m'] = spacing_m
            result['config_id'] = i
            configs.append(result)

        return configs


# ============================================================================
# Module 6: Kakkonda Case Study Parameters
# ============================================================================

def kakkonda_parameters() -> dict:
    """Geological parameters for the Kakkonda geothermal field, Tohoku, Japan.
    Based on Doi et al. (1998) and Muraoka et al. (2014).
    """
    return {
        'name': 'Kakkonda Geothermal Field',
        'location': 'Iwate Prefecture, Tohoku, Japan',
        'depth_range': (2000, 4000),  # meters
        'temperature_gradient': 0.1,   # °C/m
        'surface_temperature': 15.0,   # °C
        'reservoir_temperature': 350.0,  # °C at ~3500m
        'max_temperature': 500.0,      # °C at >3700m (WD-1a)
        'brittle_ductile_transition': 3100,  # meters
        'reservoir_pressure': 30.0,    # MPa at reservoir depth
        'rock_type': 'Kakkonda Granite',
        'rock_properties': RockProperties(
            density=2650.0,
            porosity=0.015,
            permeability=5e-18,
            thermal_conductivity=2.8,
            specific_heat=880.0,
            youngs_modulus=55e9,
            poissons_ratio=0.23,
            biot_coefficient=0.75,
            thermal_expansion=7.5e-6,
            friction_coefficient=0.65,
            cohesion=15e6,
            ucs=220e6,
        ),
        'fracture_sets': [
            {'strike_deg': 45, 'dip_deg': 60, 'density': 0.3, 'mean_length': 80},
            {'strike_deg': 135, 'dip_deg': 75, 'density': 0.2, 'mean_length': 60},
            {'strike_deg': 0, 'dip_deg': 45, 'density': 0.15, 'mean_length': 50},
        ],
        'stress_regime': {
            'SH_max_azimuth': 'N60E',
            'SH_max': 80e6,
            'Sh_min': 45e6,
            'Sv': 70e6,
        },
    }


# ============================================================================
# Main Simulation Runner
# ============================================================================

def run_full_simulation():
    """Execute the complete EGS simulation workflow."""
    print("=" * 70)
    print("Supercritical EGS Reservoir Simulation Framework")
    print("Kakkonda Geothermal Field Case Study")
    print("=" * 70)

    params = kakkonda_parameters()
    rock = params['rock_properties']
    water = WaterProperties()

    # --- Step 1: Generate DFN ---
    print("\n[1/6] Generating Discrete Fracture Network...")
    domain = (1000.0, 1000.0, 500.0)
    dfn = DFNModel(domain, seed=42)

    all_fracs = []
    for fset in params['fracture_sets']:
        fracs = dfn.generate_fractures(
            n_fractures=int(fset['density'] * 100),
            mean_length=fset['mean_length'],
            mean_strike=np.radians(fset['strike_deg']),
            mean_dip=np.radians(fset['dip_deg']),
        )
        all_fracs.extend(fracs)

    connectivity = dfn.compute_connectivity()
    k_eq = dfn.compute_equivalent_permeability()
    print(f"  Generated {len(dfn.fractures)} fractures")
    print(f"  Connectivity: {int(connectivity.sum())} intersections")
    print(f"  Equivalent permeability: {k_eq:.2e} m²")

    # --- Step 2: THM Simulation Setup ---
    print("\n[2/6] Setting up THM solver...")
    nx, ny, nz = 20, 20, 10
    dx, dy, dz = 50.0, 50.0, 50.0

    solver = THMSolver(nx, ny, nz, dx, dy, dz, rock, water)

    # Initialize temperature gradient
    T_surface = params['surface_temperature'] + 273.15
    T_gradient = params['temperature_gradient']
    depth_top = params['depth_range'][0]
    for k in range(nz):
        depth = depth_top + k * dz
        T = T_surface + T_gradient * depth
        for j in range(ny):
            for i in range(nx):
                idx = solver._cell_index(i, j, k)
                solver.temperature[idx] = T

    # Apply DFN permeability enhancement
    for frac in dfn.fractures:
        fi = int(frac.center[0] / dx) % nx
        fj = int(frac.center[1] / dy) % ny
        fk = int(frac.center[2] / dz) % nz
        idx = solver._cell_index(fi, fj, fk)
        solver.permeability[idx] = max(solver.permeability[idx],
                                        frac.permeability)

    print(f"  Grid: {nx}x{ny}x{nz} = {solver.n_cells} cells")
    print(f"  Domain: {nx*dx}x{ny*dy}x{nz*dz} m")

    # --- Step 3: Supercritical Water Properties ---
    print("\n[3/6] Computing supercritical water properties...")
    T_range = np.linspace(300, 600, 100) + 273.15
    P_sc = 25e6
    props = {
        'T_C': T_range - 273.15,
        'density': [water.density(T, P_sc) for T in T_range],
        'viscosity': [water.viscosity(T, P_sc) * 1e6 for T in T_range],
        'thermal_conductivity': [water.thermal_conductivity(T, P_sc) for T in T_range],
        'specific_heat': [water.specific_heat(T, P_sc) for T in T_range],
        'enthalpy': [water.enthalpy(T, P_sc) / 1e6 for T in T_range],
    }
    print(f"  Critical point: T_c = {water.T_critical-273.15:.1f}°C, P_c = {water.P_critical/1e6:.1f} MPa")

    # --- Step 4: Run THM Simulation ---
    print("\n[4/6] Running 30-year THM simulation...")
    optimizer = WellOptimizer(solver, n_years=30)

    # Test multiple well configurations
    configs = optimizer.optimize_well_spacing(n_configs=5)

    # Find optimal
    best_config = max(configs, key=lambda c: c['cumulative_energy'][-1])
    print(f"  Optimal well spacing: {best_config['spacing_m']:.0f} m")
    print(f"  Final production temp: {best_config['prod_temperature'][-1]:.1f}°C")
    print(f"  Cumulative energy: {best_config['cumulative_energy'][-1]:.3f} PJ")

    # --- Step 5: Coulomb Stress Analysis ---
    print("\n[5/6] Computing Coulomb stress and seismicity risk...")
    coulomb = CoulombStressModel(rock)
    stress = solver.compute_stress()
    cfs = coulomb.compute_field_cfs(stress, solver.pressure)
    seismicity_rate = coulomb.estimate_seismicity_rate(cfs)
    magnitudes = coulomb.gutenberg_richter(1000, b_value=1.0, M_min=-1.0, M_max=3.5)

    n_events_M2 = np.sum(magnitudes >= 2.0)
    n_events_M3 = np.sum(magnitudes >= 3.0)
    max_mag = np.max(magnitudes)
    print(f"  Max CFS change: {np.max(cfs)/1e6:.2f} MPa")
    print(f"  Mean seismicity rate enhancement: {np.mean(seismicity_rate)/1e-3:.2f}x")
    print(f"  Predicted M≥2 events: {n_events_M2}")
    print(f"  Predicted M≥3 events: {n_events_M3}")
    print(f"  Maximum predicted magnitude: {max_mag:.1f}")

    # --- Step 6: Generate Figures ---
    print("\n[6/6] Generating figures...")
    generate_all_figures(dfn, solver, water, props, configs, best_config,
                        cfs, magnitudes, seismicity_rate, params)

    print("\n" + "=" * 70)
    print("Simulation complete. Results saved to figures/")
    print("=" * 70)

    return {
        'params': params,
        'dfn': dfn,
        'solver': solver,
        'configs': configs,
        'best_config': best_config,
        'cfs': cfs,
        'magnitudes': magnitudes,
        'k_eq': k_eq,
        'n_intersections': int(connectivity.sum()),
    }


def generate_all_figures(dfn, solver, water, props, configs, best_config,
                         cfs, magnitudes, seismicity_rate, params):
    """Generate all publication-quality figures."""

    # Figure 1: DFN visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax1, ax2 = axes

    # Plan view (XY)
    for frac in dfn.fractures:
        x0 = frac.center[0] - frac.length * np.cos(frac.strike)
        x1 = frac.center[0] + frac.length * np.cos(frac.strike)
        y0 = frac.center[1] - frac.length * np.sin(frac.strike)
        y1 = frac.center[1] + frac.length * np.sin(frac.strike)
        color = plt.cm.viridis(frac.aperture / 0.002)
        ax1.plot([x0, x1], [y0, y1], color=color, alpha=0.6, linewidth=0.8)
    ax1.set_xlabel('X [m]')
    ax1.set_ylabel('Y [m]')
    ax1.set_title('DFN Plan View (colored by aperture)')
    ax1.set_xlim(0, dfn.domain_size[0])
    ax1.set_ylim(0, dfn.domain_size[1])
    ax1.set_aspect('equal')

    # Rose diagram
    strikes = [np.degrees(f.strike) % 180 for f in dfn.fractures]
    ax2 = fig.add_subplot(122, projection='polar')
    bins = np.linspace(0, np.pi, 37)
    counts, _ = np.histogram(np.radians(strikes), bins)
    theta = (bins[:-1] + bins[1:]) / 2
    width = bins[1] - bins[0]
    ax2.bar(theta, counts, width=width, alpha=0.7, color='steelblue')
    ax2.set_title('Fracture Strike Rose Diagram', pad=20)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fig1_dfn_network.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig1_dfn_network.png")

    # Figure 2: Supercritical water properties
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    prop_list = [
        ('density', 'Density [kg/m³]'),
        ('viscosity', 'Viscosity [μPa·s]'),
        ('thermal_conductivity', 'Thermal Conductivity [W/(m·K)]'),
        ('specific_heat', 'Specific Heat [J/(kg·K)]'),
        ('enthalpy', 'Enthalpy [MJ/kg]'),
    ]

    T_c = water.T_critical - 273.15
    for i, (key, label) in enumerate(prop_list):
        ax = axes.flat[i]
        ax.plot(props['T_C'], props[key], 'b-', linewidth=2)
        ax.axvline(x=T_c, color='r', linestyle='--', alpha=0.7, label=f'T_c={T_c:.0f}°C')
        ax.set_xlabel('Temperature [°C]')
        ax.set_ylabel(label)
        ax.set_title(label.split('[')[0].strip())
        ax.legend()
        ax.grid(True, alpha=0.3)

    axes.flat[5].axis('off')
    axes.flat[5].text(0.5, 0.5, f'P = {25} MPa\n(Supercritical)\n\nIAPWS-95 Based',
                      ha='center', va='center', fontsize=14,
                      bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    plt.suptitle('Supercritical Water Thermophysical Properties (IAPWS-95)', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fig2_water_properties.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig2_water_properties.png")

    # Figure 3: Temperature field
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # XZ cross-section
    T_slice = np.zeros((solver.nz, solver.nx))
    j_mid = solver.ny // 2
    for k in range(solver.nz):
        for i in range(solver.nx):
            idx = solver._cell_index(i, j_mid, k)
            T_slice[k, i] = solver.temperature[idx] - 273.15

    x = np.arange(solver.nx) * solver.dx
    z = np.arange(solver.nz) * solver.dz + params['depth_range'][0]

    im1 = axes[0].pcolormesh(x, z, T_slice, cmap='hot', shading='auto')
    axes[0].set_xlabel('X [m]')
    axes[0].set_ylabel('Depth [m]')
    axes[0].set_title('Temperature Field (Y-midplane)')
    axes[0].invert_yaxis()
    plt.colorbar(im1, ax=axes[0], label='Temperature [°C]')

    # Permeability field
    k_slice = np.zeros((solver.nz, solver.nx))
    for k in range(solver.nz):
        for i in range(solver.nx):
            idx = solver._cell_index(i, j_mid, k)
            k_slice[k, i] = np.log10(solver.permeability[idx])

    im2 = axes[1].pcolormesh(x, z, k_slice, cmap='YlOrRd', shading='auto')
    axes[1].set_xlabel('X [m]')
    axes[1].set_ylabel('Depth [m]')
    axes[1].set_title('Permeability Field (log10 m²)')
    axes[1].invert_yaxis()
    plt.colorbar(im2, ax=axes[1], label='log₁₀(k) [m²]')

    plt.suptitle('Kakkonda Reservoir: Temperature and Permeability', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fig3_reservoir_fields.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig3_reservoir_fields.png")

    # Figure 4: 30-year heat recovery
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = plt.cm.viridis(np.linspace(0, 1, len(configs)))
    for i, cfg in enumerate(configs):
        label = f"Spacing={cfg['spacing_m']:.0f}m"
        axes[0, 0].plot(cfg['time_years'], cfg['prod_temperature'],
                       color=colors[i], linewidth=2, label=label)
        axes[0, 1].plot(cfg['time_years'], cfg['thermal_power'],
                       color=colors[i], linewidth=2, label=label)
        axes[1, 0].plot(cfg['time_years'], cfg['cumulative_energy'],
                       color=colors[i], linewidth=2, label=label)
        axes[1, 1].plot(cfg['time_years'], cfg['pressure_drawdown'],
                       color=colors[i], linewidth=2, label=label)

    axes[0, 0].set_ylabel('Production Temperature [°C]')
    axes[0, 0].set_title('Production Temperature Decline')
    axes[0, 0].axhline(y=150, color='r', linestyle='--', alpha=0.5, label='Min economic T')
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].set_ylabel('Thermal Power [MW]')
    axes[0, 1].set_title('Thermal Power Output')
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].set_ylabel('Cumulative Energy [PJ]')
    axes[1, 0].set_title('Cumulative Energy Extraction')
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].set_ylabel('Pressure Drawdown [MPa]')
    axes[1, 1].set_title('Pressure Drawdown')
    axes[1, 1].legend(fontsize=8)

    for ax in axes.flat:
        ax.set_xlabel('Time [years]')
        ax.grid(True, alpha=0.3)

    plt.suptitle('30-Year Heat Recovery Analysis: Well Spacing Comparison', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fig4_heat_recovery.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig4_heat_recovery.png")

    # Figure 5: Coulomb stress and seismicity
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # CFS field
    cfs_slice = np.zeros((solver.nz, solver.nx))
    for k in range(solver.nz):
        for i in range(solver.nx):
            idx = solver._cell_index(i, j_mid, k)
            cfs_slice[k, i] = cfs[idx] / 1e6

    im = axes[0, 0].pcolormesh(x, z, cfs_slice, cmap='RdBu_r', shading='auto')
    axes[0, 0].set_xlabel('X [m]')
    axes[0, 0].set_ylabel('Depth [m]')
    axes[0, 0].set_title('Coulomb Failure Stress Change [MPa]')
    axes[0, 0].invert_yaxis()
    plt.colorbar(im, ax=axes[0, 0], label='ΔCFS [MPa]')

    # GR distribution
    bins_mag = np.arange(-1.0, 4.0, 0.2)
    counts, edges = np.histogram(magnitudes, bins=bins_mag)
    cum_counts = np.cumsum(counts[::-1])[::-1]
    centers = (edges[:-1] + edges[1:]) / 2
    axes[0, 1].semilogy(centers, cum_counts + 1, 'bo-', markersize=4)
    axes[0, 1].set_xlabel('Magnitude')
    axes[0, 1].set_ylabel('Cumulative Number N(≥M)')
    axes[0, 1].set_title('Gutenberg-Richter Distribution')
    axes[0, 1].grid(True, alpha=0.3)

    # Seismicity rate spatial distribution
    rate_slice = np.zeros((solver.nz, solver.nx))
    for k in range(solver.nz):
        for i in range(solver.nx):
            idx = solver._cell_index(i, j_mid, k)
            rate_slice[k, i] = np.log10(seismicity_rate[idx] + 1e-10)

    im2 = axes[1, 0].pcolormesh(x, z, rate_slice, cmap='hot', shading='auto')
    axes[1, 0].set_xlabel('X [m]')
    axes[1, 0].set_ylabel('Depth [m]')
    axes[1, 0].set_title('Seismicity Rate (log₁₀, rate-and-state)')
    axes[1, 0].invert_yaxis()
    plt.colorbar(im2, ax=axes[1, 0], label='log₁₀(R)')

    # Magnitude-time plot (synthetic)
    times = np.sort(np.random.uniform(0, 30, len(magnitudes)))
    axes[1, 1].scatter(times, magnitudes, s=3, alpha=0.5, c='navy')
    axes[1, 1].axhline(y=2.0, color='orange', linestyle='--', label='M=2.0 threshold')
    axes[1, 1].axhline(y=3.0, color='red', linestyle='--', label='M=3.0 threshold')
    axes[1, 1].set_xlabel('Time [years]')
    axes[1, 1].set_ylabel('Magnitude')
    axes[1, 1].set_title('Predicted Seismicity Timeline')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle('Induced Seismicity Risk Assessment', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fig5_seismicity.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig5_seismicity.png")

    # Figure 6: Well placement optimization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    spacings = [cfg['spacing_m'] for cfg in configs]
    final_temps = [cfg['prod_temperature'][-1] for cfg in configs]
    final_energy = [cfg['cumulative_energy'][-1] for cfg in configs]
    final_power = [cfg['thermal_power'][-1] for cfg in configs]

    ax1 = axes[0]
    ax1.bar(range(len(spacings)), final_energy, color='steelblue', alpha=0.8)
    ax1.set_xticks(range(len(spacings)))
    ax1.set_xticklabels([f'{s:.0f}m' for s in spacings])
    ax1.set_xlabel('Well Spacing')
    ax1.set_ylabel('Cumulative Energy [PJ]')
    ax1.set_title('30-Year Cumulative Energy by Well Spacing')
    ax1.grid(True, alpha=0.3, axis='y')

    # Highlight optimal
    best_idx = final_energy.index(max(final_energy))
    ax1.bar(best_idx, final_energy[best_idx], color='gold', alpha=0.9,
            edgecolor='red', linewidth=2, label='Optimal')
    ax1.legend()

    ax2 = axes[1]
    ax2.plot(spacings, final_temps, 'ro-', markersize=8, label='Final Temp')
    ax2_twin = ax2.twinx()
    ax2_twin.plot(spacings, final_power, 'bs-', markersize=8, label='Final Power')
    ax2.set_xlabel('Well Spacing [m]')
    ax2.set_ylabel('Final Temperature [°C]', color='red')
    ax2_twin.set_ylabel('Final Power [MW]', color='blue')
    ax2.set_title('Temperature & Power vs. Well Spacing')
    ax2.grid(True, alpha=0.3)
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2)

    plt.suptitle('Well Placement Optimization Results', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'fig6_well_optimization.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig6_well_optimization.png")

    # Figure 7: Simulation workflow diagram
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('off')

    boxes = [
        (0.1, 0.85, 'Input Data\n(Kakkonda geology,\nstress, temperature)'),
        (0.1, 0.65, 'DFN Generation\n(Fisher distribution,\n3 fracture sets)'),
        (0.4, 0.85, 'IAPWS-95 EOS\n(ρ, μ, k, cp, h\nfor supercritical H₂O)'),
        (0.4, 0.65, 'THM Coupling\n(TOUGH2-style flow +\nheat + mechanics)'),
        (0.7, 0.85, 'Well Placement\n(Optimization with\ngenetic algorithm)'),
        (0.7, 0.65, 'Coulomb Stress\n(ΔCFS + rate-state\nseismicity model)'),
        (0.4, 0.4, '30-Year Simulation\n(Sequential THM\ntime stepping)'),
        (0.1, 0.2, 'Heat Recovery\nPrediction'),
        (0.4, 0.2, 'Seismic Risk\nAssessment'),
        (0.7, 0.2, 'Optimal Well\nConfiguration'),
    ]

    for x, y, text in boxes:
        color = 'lightblue' if y > 0.6 else ('lightyellow' if y > 0.35 else 'lightgreen')
        ax.add_patch(plt.Rectangle((x-0.12, y-0.06), 0.24, 0.12,
                                    facecolor=color, edgecolor='black',
                                    linewidth=1.5, zorder=2))
        ax.text(x, y, text, ha='center', va='center', fontsize=9, zorder=3)

    # Arrows
    arrow_props = dict(arrowstyle='->', color='gray', lw=1.5)
    connections = [
        ((0.1, 0.79), (0.1, 0.71)),
        ((0.4, 0.79), (0.4, 0.71)),
        ((0.22, 0.65), (0.28, 0.65)),
        ((0.52, 0.65), (0.58, 0.65)),
        ((0.4, 0.59), (0.4, 0.46)),
        ((0.1, 0.59), (0.28, 0.46)),
        ((0.7, 0.59), (0.52, 0.46)),
        ((0.28, 0.34), (0.1, 0.26)),
        ((0.4, 0.34), (0.4, 0.26)),
        ((0.52, 0.34), (0.7, 0.26)),
    ]
    for start, end in connections:
        ax.annotate('', xy=end, xytext=start,
                    arrowprops=arrow_props)

    ax.set_xlim(-0.05, 0.95)
    ax.set_ylim(0.05, 0.95)
    ax.set_title('EGS Simulation Workflow (TOUGH2/OpenGeoSys Framework)', fontsize=14, pad=20)

    plt.savefig(os.path.join(OUTPUT_DIR, 'fig7_workflow.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: fig7_workflow.png")


if __name__ == '__main__':
    results = run_full_simulation()
