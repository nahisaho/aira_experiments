#!/usr/bin/env python3
"""
Molecular Simulation Protocol for Concentrated Electrolyte Solutions
=====================================================================
Implements:
  1. Force field parameter optimization (ion-water, ion-ion interactions)
  2. Activity coefficient / osmotic pressure via Kirkwood-Buff integrals
  3. Ion transport (diffusion, conductivity) via Green-Kubo relations
  4. Solvation structure analysis (coordination number, solvation free energy)
  5. Anomalous transport phenomena in concentrated electrolytes
  6. Case study: EC/DMC/LiPF6 lithium-ion battery electrolyte

This script generates GROMACS/LAMMPS input protocols and performs
numerical calculations using model data to demonstrate the methodology.
"""

import numpy as np
import os

# ============================================================
# Constants
# ============================================================
KB = 1.380649e-23       # Boltzmann constant (J/K)
NA = 6.02214076e23      # Avogadro's number
E_CHARGE = 1.602176634e-19  # Elementary charge (C)
EPSILON_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
KB_KCAL = 1.987204e-3   # Boltzmann constant (kcal/mol/K)

# ============================================================
# 1. Force Field Parameters
# ============================================================

class ForceFieldParameters:
    """Force field parameters for EC/DMC/LiPF6 system."""
    
    def __init__(self):
        # Lennard-Jones parameters (sigma in nm, epsilon in kJ/mol)
        self.lj_params = {
            'Li': {'sigma': 0.1506, 'epsilon': 0.6947, 'charge': 0.80},
            'P':  {'sigma': 0.3740, 'epsilon': 0.8368, 'charge': 1.34},
            'F':  {'sigma': 0.3118, 'epsilon': 0.2552, 'charge': -0.39},
            'O_EC':  {'sigma': 0.2960, 'epsilon': 0.8786, 'charge': -0.4684},
            'C_EC':  {'sigma': 0.3400, 'epsilon': 0.3598, 'charge': 0.7714},
            'O_DMC': {'sigma': 0.3000, 'epsilon': 0.7113, 'charge': -0.4374},
            'C_DMC': {'sigma': 0.3500, 'epsilon': 0.2761, 'charge': 0.6580},
            'OW':  {'sigma': 0.3166, 'epsilon': 0.6502, 'charge': -0.8476},
            'HW':  {'sigma': 0.0000, 'epsilon': 0.0000, 'charge': 0.4238},
        }
        # Scaled charges (0.8 scaling for concentrated electrolytes)
        self.charge_scale = 0.80
        
    def get_cross_params(self, atom1, atom2):
        """Lorentz-Berthelot combining rules."""
        p1 = self.lj_params[atom1]
        p2 = self.lj_params[atom2]
        sigma_ij = (p1['sigma'] + p2['sigma']) / 2.0
        epsilon_ij = np.sqrt(p1['epsilon'] * p2['epsilon'])
        return sigma_ij, epsilon_ij
    
    def optimize_parameters(self, target_density, target_diff_coeff, T=298.15):
        """
        Iterative force field optimization against experimental data.
        Uses a simplified gradient-descent on LJ parameters.
        """
        print("=" * 60)
        print("Force Field Parameter Optimization")
        print("=" * 60)
        
        # Target properties for 1M LiPF6 in EC:DMC (1:1)
        rho_exp = target_density       # g/cm^3
        D_exp = target_diff_coeff      # cm^2/s
        
        # Initial parameters
        sigma_Li = self.lj_params['Li']['sigma']
        eps_Li = self.lj_params['Li']['epsilon']
        
        results = []
        for iteration in range(5):
            # Simulate density with current parameters
            rho_sim = self._compute_density(sigma_Li, eps_Li, T)
            D_sim = self._compute_diffusion(sigma_Li, eps_Li, T)
            
            err_rho = abs(rho_sim - rho_exp) / rho_exp * 100
            err_D = abs(D_sim - D_exp) / D_exp * 100
            obj = np.sqrt(err_rho**2 + err_D**2)
            
            results.append({
                'iteration': iteration + 1,
                'sigma_Li': sigma_Li,
                'eps_Li': eps_Li,
                'rho_sim': rho_sim,
                'rho_exp': rho_exp,
                'D_sim': D_sim,
                'D_exp': D_exp,
                'err_rho': err_rho,
                'err_D': err_D,
                'objective': obj,
            })
            
            print(f"Iter {iteration+1}: σ_Li={sigma_Li:.4f} nm, ε_Li={eps_Li:.4f} kJ/mol")
            print(f"  ρ_sim={rho_sim:.4f} g/cm³ (exp: {rho_exp:.4f}), err={err_rho:.2f}%")
            print(f"  D_sim={D_sim:.2e} cm²/s (exp: {D_exp:.2e}), err={err_D:.2f}%")
            
            # Gradient step
            d_sigma = 0.001 * (rho_exp - rho_sim) / rho_exp
            d_eps = 0.02 * (D_exp - D_sim) / D_exp
            sigma_Li += d_sigma
            eps_Li += d_eps
        
        # Update optimized values
        self.lj_params['Li']['sigma'] = sigma_Li
        self.lj_params['Li']['epsilon'] = eps_Li
        
        return results
    
    def _compute_density(self, sigma_Li, eps_Li, T):
        """Model density as function of LJ params (simplified)."""
        base = 1.2050  # base density g/cm^3
        correction = -0.15 * (sigma_Li - 0.15) + 0.03 * (eps_Li - 0.70)
        noise = np.random.normal(0, 0.001)
        return base + correction + noise
    
    def _compute_diffusion(self, sigma_Li, eps_Li, T):
        """Model diffusion coefficient (simplified)."""
        base = 2.5e-6  # cm^2/s
        correction = base * (0.3 * (eps_Li - 0.70) - 0.5 * (sigma_Li - 0.15))
        noise = np.random.normal(0, base * 0.02)
        return base + correction + noise


# ============================================================
# 2. Kirkwood-Buff Integration
# ============================================================

class KirkwoodBuffAnalysis:
    """Kirkwood-Buff integral analysis for activity/osmotic coefficients."""
    
    def __init__(self, T=298.15, box_length=5.0):
        self.T = T
        self.box_length = box_length  # nm
        
    def compute_rdf(self, pair_type, concentration, n_points=500):
        """
        Compute radial distribution function g(r) for given ion pair.
        Uses analytical models calibrated to MD data.
        """
        r = np.linspace(0.15, self.box_length / 2.0, n_points)
        
        if pair_type == 'Li-OW':
            # Li+ - water oxygen: strong first shell
            r0 = 0.196
            sigma1 = 0.015
            A1 = 8.5 - 1.2 * concentration
            r1 = 0.395
            sigma2 = 0.04
            A2 = 1.8 + 0.3 * concentration
        elif pair_type == 'Li-Li':
            r0 = 0.42
            sigma1 = 0.06
            A1 = 1.5 + 2.0 * concentration
            r1 = 0.68
            sigma2 = 0.08
            A2 = 1.2 + 0.8 * concentration
        elif pair_type == 'Li-PF6':
            r0 = 0.35
            sigma1 = 0.04
            A1 = 4.0 + 3.5 * concentration
            r1 = 0.58
            sigma2 = 0.06
            A2 = 1.5 + 1.0 * concentration
        elif pair_type == 'Li-O_EC':
            r0 = 0.200
            sigma1 = 0.018
            A1 = 6.5 - 0.5 * concentration
            r1 = 0.410
            sigma2 = 0.045
            A2 = 1.5 + 0.2 * concentration
        elif pair_type == 'PF6-OW':
            r0 = 0.36
            sigma1 = 0.04
            A1 = 2.5 - 0.3 * concentration
            r1 = 0.60
            sigma2 = 0.06
            A2 = 1.2
        else:
            r0 = 0.30
            sigma1 = 0.03
            A1 = 2.0
            r1 = 0.50
            sigma2 = 0.05
            A2 = 1.3
        
        g = np.ones_like(r)
        g += A1 * np.exp(-((r - r0) / sigma1)**2)
        g += A2 * np.exp(-((r - r1) / sigma2)**2)
        # Depletion region
        g *= (1.0 - np.exp(-((r / (r0 * 0.7))**12)))
        # Long-range oscillation damping
        g += 0.1 * np.exp(-r / 0.5) * np.cos(2 * np.pi * r / 0.3)
        g = np.maximum(g, 0)
        
        return r, g
    
    def compute_kb_integral(self, r, g):
        """
        Compute Kirkwood-Buff integral:
          G_ij = 4π ∫₀^∞ [g_ij(r) - 1] r² dr
        """
        integrand = 4.0 * np.pi * (g - 1.0) * r**2
        G = np.cumsum(integrand) * (r[1] - r[0])
        return G
    
    def compute_activity_coefficient(self, concentrations):
        """
        Compute mean ionic activity coefficient γ± from KB integrals.
        ln(γ±) depends on KB integrals G++, G--, G+-, G+s, G-s, Gss
        """
        print("\n" + "=" * 60)
        print("Kirkwood-Buff Analysis: Activity & Osmotic Coefficients")
        print("=" * 60)
        
        results = []
        for c in concentrations:
            rho_s = 55.5 - 2.0 * c  # mol/L solvent (approximate)
            rho_i = c
            
            # Compute KB integrals for key pairs
            r, g_pp = self.compute_rdf('Li-Li', c)
            G_pp = self.compute_kb_integral(r, g_pp)[-1]  # nm^3
            
            r, g_mm = self.compute_rdf('PF6-OW', c)
            G_mm = self.compute_kb_integral(r, g_mm)[-1]
            
            r, g_pm = self.compute_rdf('Li-PF6', c)
            G_pm = self.compute_kb_integral(r, g_pm)[-1]
            
            r, g_ps = self.compute_rdf('Li-OW', c)
            G_ps = self.compute_kb_integral(r, g_ps)[-1]
            
            # Activity coefficient (simplified KB expression)
            # ln(γ±) ≈ -ρ_s * (G_ps - G_ss) / (1 + ρ_i * (G_pp + G_mm - 2*G_pm))
            G_ss = -0.02  # solvent-solvent KB integral (nm^3)
            
            denom = 1.0 + rho_i * (G_pp + G_mm - 2 * G_pm)
            if abs(denom) < 0.01:
                denom = 0.01
            ln_gamma = -rho_s * (G_ps - G_ss) / denom
            gamma = np.exp(ln_gamma)
            
            # Osmotic coefficient
            phi = 1.0 - ln_gamma / 2.0
            phi = max(0.5, min(phi, 1.5))
            
            results.append({
                'concentration': c,
                'G_pp': G_pp,
                'G_pm': G_pm,
                'G_ps': G_ps,
                'gamma': gamma,
                'phi': phi,
                'ln_gamma': ln_gamma,
            })
            
            print(f"c = {c:.1f} M: γ± = {gamma:.4f}, φ = {phi:.4f}")
            print(f"  G_++ = {G_pp:.3f} nm³, G_+- = {G_pm:.3f} nm³, G_+s = {G_ps:.3f} nm³")
        
        return results


# ============================================================
# 3. Green-Kubo Transport Properties
# ============================================================

class GreenKuboTransport:
    """Green-Kubo and Einstein relation transport calculations."""
    
    def __init__(self, T=298.15, dt=0.001, n_steps=100000):
        self.T = T
        self.dt = dt  # ps
        self.n_steps = n_steps
        
    def generate_velocity_trajectory(self, mass, D_target, n_particles=100):
        """
        Generate velocity autocorrelation-consistent trajectories.
        Uses Ornstein-Uhlenbeck process.
        """
        gamma = KB * self.T / (mass * 1.66054e-27 * D_target * 1e-4)  # friction
        sigma_v = np.sqrt(KB * self.T / (mass * 1.66054e-27))
        
        n_frames = min(self.n_steps, 50000)
        velocities = np.zeros((n_particles, n_frames, 3))
        
        for p in range(n_particles):
            v = np.random.normal(0, sigma_v, 3)
            for t in range(n_frames):
                velocities[p, t] = v
                noise = np.random.normal(0, 1, 3)
                v = v * np.exp(-gamma * self.dt * 1e-12) + \
                    sigma_v * np.sqrt(1 - np.exp(-2 * gamma * self.dt * 1e-12)) * noise
        
        return velocities
    
    def compute_vacf(self, velocities):
        """Velocity autocorrelation function."""
        n_particles, n_frames, _ = velocities.shape
        max_lag = min(5000, n_frames // 2)
        vacf = np.zeros(max_lag)
        
        for lag in range(max_lag):
            corr = 0.0
            count = 0
            for p in range(n_particles):
                for t in range(n_frames - lag):
                    corr += np.dot(velocities[p, t], velocities[p, t + lag])
                    count += 1
                if count > 100000:
                    break
            vacf[lag] = corr / count if count > 0 else 0.0
        
        vacf /= vacf[0] if vacf[0] != 0 else 1.0
        return vacf
    
    def compute_msd(self, D_target, n_particles=100, concentrations=None):
        """
        Compute mean squared displacement for multiple concentrations.
        MSD = 6Dt for 3D diffusion.
        """
        if concentrations is None:
            concentrations = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0]
        
        print("\n" + "=" * 60)
        print("Green-Kubo / Einstein Transport Analysis")
        print("=" * 60)
        
        t_max = 1000.0  # ps
        t = np.linspace(0, t_max, 5000)
        
        results = []
        for c in concentrations:
            # Concentration-dependent diffusion
            D_c = D_target * np.exp(-0.35 * c)  # Decreases with concentration
            
            # Add subdiffusive regime at short times
            alpha = 1.0 - 0.08 * c  # anomalous exponent
            alpha = max(alpha, 0.7)
            
            # MSD with crossover from subdiffusive to diffusive
            tau_cross = 10.0 + 5.0 * c  # crossover time (ps)
            msd = np.where(
                t < tau_cross,
                6 * D_c * 1e-4 * (t * 1e-12)**alpha / (1e-12)**alpha * 1e20,
                6 * D_c * 1e-4 * t * 1e-12 * 1e20 + \
                    6 * D_c * 1e-4 * tau_cross * 1e-12 * 1e20 * \
                    ((tau_cross * 1e-12)**(alpha-1) - 1)
            )
            # Simpler: use smooth interpolation
            msd = 6 * D_c * 1e-4 * ((t * 1e-12)**alpha + (t * 1e-12)) / 2.0 * 1e20
            
            # Compute D from long-time MSD slope
            idx_fit = t > 500
            if np.sum(idx_fit) > 10:
                coeffs = np.polyfit(t[idx_fit] * 1e-12, msd[idx_fit] * 1e-20, 1)
                D_computed = coeffs[0] / 6.0 * 1e4  # cm^2/s
            else:
                D_computed = D_c
            
            results.append({
                'concentration': c,
                'D_target': D_c,
                'D_computed': D_computed,
                'alpha': alpha,
                't': t,
                'msd': msd,
            })
            
            print(f"c = {c:.1f} M: D = {D_computed:.2e} cm²/s, α = {alpha:.3f}")
        
        return results
    
    def compute_conductivity(self, concentrations, D_cation, D_anion, 
                             z_cat=1, z_an=-1):
        """
        Compute ionic conductivity from Nernst-Einstein and Green-Kubo.
        σ_NE = (F²/RT) * Σ c_i z_i² D_i  (Nernst-Einstein, upper bound)
        σ_GK includes cross-correlations (Haven ratio)
        """
        F = 96485.0  # Faraday constant C/mol
        R = 8.314    # J/mol/K
        
        results = []
        for c in concentrations:
            D_cat = D_cation * np.exp(-0.35 * c)
            D_an = D_anion * np.exp(-0.30 * c)
            
            # Nernst-Einstein
            sigma_NE = (F**2 / (R * self.T)) * c * 1000 * \
                       (z_cat**2 * D_cat * 1e-4 + z_an**2 * D_an * 1e-4)
            sigma_NE *= 10  # S/m -> mS/cm conversion factor adjustment
            
            # Haven ratio (accounts for ion-ion correlations)
            # H < 1 for associated electrolytes
            H = 0.75 - 0.10 * c  # decreases at higher concentration
            H = max(H, 0.35)
            
            sigma_GK = sigma_NE * H
            
            # Transference number
            t_plus = z_cat**2 * D_cat / (z_cat**2 * D_cat + z_an**2 * D_an)
            
            results.append({
                'concentration': c,
                'sigma_NE': sigma_NE,
                'sigma_GK': sigma_GK,
                'Haven_ratio': H,
                't_plus': t_plus,
                'D_cat': D_cat,
                'D_an': D_an,
            })
            
            print(f"c = {c:.1f} M: σ_NE = {sigma_NE:.2f} mS/cm, σ_GK = {sigma_GK:.2f} mS/cm, "
                  f"H = {H:.3f}, t+ = {t_plus:.3f}")
        
        return results


# ============================================================
# 4. Solvation Structure Analysis
# ============================================================

class SolvationAnalysis:
    """Solvation structure: coordination number and free energy."""
    
    def __init__(self, T=298.15):
        self.T = T
    
    def compute_coordination_number(self, r, g, r_cutoff, rho_bulk):
        """
        Coordination number from RDF:
          n(r) = 4π ρ ∫₀^r_cut g(r') r'² dr'
        """
        mask = r <= r_cutoff
        integrand = 4.0 * np.pi * rho_bulk * g[mask] * r[mask]**2
        cn = np.trapz(integrand, r[mask])
        return cn
    
    def compute_solvation_free_energy(self, r, g, pair_type='Li-OW'):
        """
        Potential of mean force (PMF):
          w(r) = -k_B T ln[g(r)]
        """
        g_safe = np.maximum(g, 1e-10)
        pmf = -KB_KCAL * self.T * np.log(g_safe)
        return pmf
    
    def analyze_solvation_shells(self, kb_analysis, concentrations):
        """Analyze solvation structure at multiple concentrations."""
        print("\n" + "=" * 60)
        print("Solvation Structure Analysis")
        print("=" * 60)
        
        pairs = ['Li-OW', 'Li-O_EC', 'Li-PF6']
        cutoffs = {'Li-OW': 0.28, 'Li-O_EC': 0.30, 'Li-PF6': 0.45}
        rho_bulk = {'Li-OW': 33.3, 'Li-O_EC': 10.0, 'Li-PF6': 0.6}  # nm^-3
        
        results = {}
        for pair in pairs:
            results[pair] = []
            for c in concentrations:
                rho = rho_bulk[pair]
                if pair == 'Li-OW':
                    rho = 33.3 * (1.0 - 0.05 * c)
                elif pair == 'Li-O_EC':
                    rho = 10.0 * (1.0 + 0.02 * c)
                elif pair == 'Li-PF6':
                    rho = 0.6 * c
                
                r, g = kb_analysis.compute_rdf(pair, c)
                cn = self.compute_coordination_number(r, g, cutoffs[pair], rho)
                pmf = self.compute_solvation_free_energy(r, g, pair)
                pmf_min = np.min(pmf)
                
                results[pair].append({
                    'concentration': c,
                    'coord_number': cn,
                    'pmf_min': pmf_min,
                    'r': r,
                    'g': g,
                    'pmf': pmf,
                })
                
                print(f"{pair}, c={c:.1f} M: CN = {cn:.2f}, ΔG_solv = {pmf_min:.2f} kcal/mol")
        
        return results


# ============================================================
# 5. Anomalous Transport Analysis
# ============================================================

class AnomalousTransport:
    """Analysis of anomalous transport in concentrated electrolytes."""
    
    def __init__(self, T=298.15):
        self.T = T
    
    def compute_anomalous_exponent(self, t, msd):
        """
        Compute local anomalous exponent:
          α(t) = d[ln(MSD)] / d[ln(t)]
        """
        log_t = np.log(t[1:])
        log_msd = np.log(np.maximum(msd[1:], 1e-30))
        alpha = np.gradient(log_msd, log_t)
        return t[1:], alpha
    
    def analyze_anomalous_transport(self, concentrations):
        """Compute anomalous transport indicators."""
        print("\n" + "=" * 60)
        print("Anomalous Transport Phenomena")
        print("=" * 60)
        
        t = np.logspace(-1, 3, 1000)  # ps
        D0 = 2.5e-6  # cm^2/s reference
        
        results = []
        for c in concentrations:
            alpha_long = 1.0 - 0.05 * c
            alpha_long = max(alpha_long, 0.70)
            
            alpha_short = 0.5 + 0.1 * (1.0 - c / 5.0)
            alpha_short = max(alpha_short, 0.3)
            
            # Crossover timescale
            tau_c = 5.0 * np.exp(0.5 * c)  # ps
            
            # MSD with anomalous regimes
            D_c = D0 * np.exp(-0.35 * c)
            msd = np.zeros_like(t)
            for i, ti in enumerate(t):
                if ti < tau_c:
                    msd[i] = 6 * D_c * 1e-4 * (ti * 1e-12)**alpha_short * 1e20
                else:
                    msd_cross = 6 * D_c * 1e-4 * (tau_c * 1e-12)**alpha_short * 1e20
                    msd[i] = msd_cross * (ti / tau_c)**alpha_long
            
            # Compute local exponent
            t_alpha, alpha_t = self.compute_anomalous_exponent(t, msd)
            
            # Non-Gaussian parameter
            ngp = 0.2 * c * np.exp(-t / (2 * tau_c))
            
            results.append({
                'concentration': c,
                'alpha_short': alpha_short,
                'alpha_long': alpha_long,
                'tau_crossover': tau_c,
                't': t,
                'msd': msd,
                't_alpha': t_alpha,
                'alpha_t': alpha_t,
                'ngp': ngp,
            })
            
            print(f"c = {c:.1f} M: α_short = {alpha_short:.3f}, α_long = {alpha_long:.3f}, "
                  f"τ_cross = {tau_c:.1f} ps")
        
        return results


# ============================================================
# 6. GROMACS/LAMMPS Protocol Generator
# ============================================================

class SimulationProtocol:
    """Generate GROMACS and LAMMPS input files."""
    
    @staticmethod
    def generate_gromacs_mdp(filename, run_type='production'):
        """Generate GROMACS .mdp file."""
        params = {
            'em': {
                'integrator': 'steep',
                'nsteps': '50000',
                'emtol': '100.0',
                'emstep': '0.01',
            },
            'nvt': {
                'integrator': 'md',
                'nsteps': '500000',
                'dt': '0.001',
                'tcoupl': 'v-rescale',
                'ref_t': '298.15',
                'tau_t': '0.1',
                'gen_vel': 'yes',
                'gen_temp': '298.15',
            },
            'npt': {
                'integrator': 'md',
                'nsteps': '1000000',
                'dt': '0.001',
                'tcoupl': 'v-rescale',
                'ref_t': '298.15',
                'tau_t': '0.1',
                'pcoupl': 'Parrinello-Rahman',
                'ref_p': '1.0',
                'tau_p': '2.0',
                'compressibility': '4.5e-5',
            },
            'production': {
                'integrator': 'md',
                'nsteps': '50000000',
                'dt': '0.001',
                'tcoupl': 'Nose-Hoover',
                'ref_t': '298.15',
                'tau_t': '1.0',
                'pcoupl': 'Parrinello-Rahman',
                'ref_p': '1.0',
                'tau_p': '5.0',
                'compressibility': '4.5e-5',
                'nstxout': '5000',
                'nstvout': '5000',
                'nstenergy': '1000',
                'nstlog': '1000',
            },
        }
        
        common = {
            'cutoff-scheme': 'Verlet',
            'nstlist': '20',
            'ns_type': 'grid',
            'coulombtype': 'PME',
            'rcoulomb': '1.2',
            'rvdw': '1.2',
            'pbc': 'xyz',
            'DispCorr': 'EnerPres',
            'constraints': 'h-bonds',
            'constraint_algorithm': 'LINCS',
        }
        
        p = params.get(run_type, params['production'])
        p.update(common)
        
        content = f"; GROMACS MDP file for {run_type}\n"
        content += f"; EC/DMC/LiPF6 electrolyte simulation\n\n"
        for key, val in p.items():
            content += f"{key:<25s} = {val}\n"
        
        with open(filename, 'w') as f:
            f.write(content)
        print(f"Generated: {filename}")
        return content
    
    @staticmethod
    def generate_lammps_input(filename):
        """Generate LAMMPS input script for EC/DMC/LiPF6."""
        content = """# LAMMPS input script for EC/DMC/LiPF6 electrolyte simulation
# Force field: OPLS-AA with scaled charges (0.8)

units           real
atom_style      full
boundary        p p p
pair_style      lj/cut/coul/long 12.0
kspace_style    pppm 1.0e-5
bond_style      harmonic
angle_style     harmonic
dihedral_style  opls

# Read data
read_data       ec_dmc_lipf6.data

# Force field parameters (scaled charges)
# Li+
pair_coeff  1 1  0.1660  1.5060  # Li-Li (sigma in Angstrom, eps in kcal/mol)
# P (PF6-)
pair_coeff  2 2  0.2000  3.7400
# F (PF6-)
pair_coeff  3 3  0.0610  3.1180

# Mixing rules: arithmetic/geometric (Lorentz-Berthelot)
pair_modify     mix arithmetic

# Neighbor list
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# ---- Energy Minimization ----
minimize        1.0e-4 1.0e-6 10000 100000

# ---- NVT Equilibration ----
velocity        all create 298.15 12345 dist gaussian
fix             1 all nvt temp 298.15 298.15 100.0
timestep        1.0
thermo          1000
thermo_style    custom step temp press pe ke etotal density
run             500000
unfix           1

# ---- NPT Equilibration ----
fix             2 all npt temp 298.15 298.15 100.0 iso 1.0 1.0 1000.0
run             1000000
unfix           2

# ---- Production Run (NVT for transport) ----
reset_timestep  0
fix             3 all nvt temp 298.15 298.15 100.0

# Output trajectories for analysis
dump            1 all custom 5000 traj.lammpstrj id type x y z vx vy vz
dump_modify     1 sort id

# Compute MSD
group           Li type 1
group           PF6 type 2 3
compute         msd_Li Li msd
compute         msd_PF6 PF6 msd
fix             msd_out all ave/time 100 1 100 c_msd_Li[4] c_msd_PF6[4] &
                file msd_output.dat

# Compute VACF for Green-Kubo conductivity
compute         vacf_Li Li vacf
fix             vacf_out Li ave/correlate 1 10000 10000 c_vacf_Li[1] c_vacf_Li[2] &
                c_vacf_Li[3] c_vacf_Li[4] type auto file vacf_Li.dat ave running

# Compute charge flux for ionic conductivity
compute         flux all heat/flux ke_atom pe_atom stress/atom NULL
# Custom ionic current computation would use fix ave/correlate

thermo          5000
run             50000000

write_data      final.data
"""
        with open(filename, 'w') as f:
            f.write(content)
        print(f"Generated: {filename}")
        return content


# ============================================================
# Main Execution
# ============================================================

def main():
    np.random.seed(42)
    output_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(output_dir)
    fig_dir = os.path.join(base_dir, 'figures')
    scripts_dir = os.path.join(base_dir, 'scripts')
    os.makedirs(fig_dir, exist_ok=True)
    
    concentrations = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0]
    
    # ---- 1. Force Field Optimization ----
    ff = ForceFieldParameters()
    ff_results = ff.optimize_parameters(
        target_density=1.2050,
        target_diff_coeff=2.5e-6,
        T=298.15
    )
    
    # ---- 2. Kirkwood-Buff Analysis ----
    kb = KirkwoodBuffAnalysis(T=298.15, box_length=5.0)
    kb_results = kb.compute_activity_coefficient(concentrations)
    
    # ---- 3. Transport Properties ----
    gk = GreenKuboTransport(T=298.15)
    msd_results = gk.compute_msd(D_target=2.5e-6, concentrations=concentrations)
    cond_results = gk.compute_conductivity(
        concentrations, D_cation=2.5e-6, D_anion=3.0e-6
    )
    
    # ---- 4. Solvation Analysis ----
    solv = SolvationAnalysis(T=298.15)
    solv_results = solv.analyze_solvation_shells(kb, concentrations)
    
    # ---- 5. Anomalous Transport ----
    anom = AnomalousTransport(T=298.15)
    anom_results = anom.analyze_anomalous_transport(concentrations)
    
    # ---- 6. Generate Simulation Input Files ----
    proto = SimulationProtocol()
    for run_type in ['em', 'nvt', 'npt', 'production']:
        proto.generate_gromacs_mdp(
            os.path.join(scripts_dir, f'{run_type}.mdp'), run_type
        )
    proto.generate_lammps_input(os.path.join(scripts_dir, 'lammps_input.in'))
    
    # ---- Save numerical results ----
    save_results(base_dir, ff_results, kb_results, msd_results, 
                 cond_results, solv_results, anom_results, concentrations, kb)
    
    print("\n" + "=" * 60)
    print("All calculations complete. Results saved to figures/ and scripts/")
    print("=" * 60)
    
    return {
        'ff': ff_results,
        'kb': kb_results,
        'msd': msd_results,
        'cond': cond_results,
        'solv': solv_results,
        'anom': anom_results,
    }


def save_results(base_dir, ff_results, kb_results, msd_results,
                 cond_results, solv_results, anom_results, concentrations, kb):
    """Save all numerical results as CSV files."""
    fig_dir = os.path.join(base_dir, 'figures')
    
    # FF optimization convergence
    with open(os.path.join(fig_dir, 'ff_optimization.csv'), 'w') as f:
        f.write("iteration,sigma_Li,eps_Li,rho_sim,err_rho,D_sim,err_D,objective\n")
        for r in ff_results:
            f.write(f"{r['iteration']},{r['sigma_Li']:.6f},{r['eps_Li']:.6f},"
                    f"{r['rho_sim']:.6f},{r['err_rho']:.4f},"
                    f"{r['D_sim']:.4e},{r['err_D']:.4f},{r['objective']:.4f}\n")
    
    # Activity coefficients
    with open(os.path.join(fig_dir, 'activity_coefficients.csv'), 'w') as f:
        f.write("concentration,gamma,phi,G_pp,G_pm,G_ps,ln_gamma\n")
        for r in kb_results:
            f.write(f"{r['concentration']:.2f},{r['gamma']:.6f},{r['phi']:.6f},"
                    f"{r['G_pp']:.6f},{r['G_pm']:.6f},{r['G_ps']:.6f},{r['ln_gamma']:.6f}\n")
    
    # Conductivity
    with open(os.path.join(fig_dir, 'conductivity.csv'), 'w') as f:
        f.write("concentration,sigma_NE,sigma_GK,Haven_ratio,t_plus\n")
        for r in cond_results:
            f.write(f"{r['concentration']:.2f},{r['sigma_NE']:.4f},{r['sigma_GK']:.4f},"
                    f"{r['Haven_ratio']:.4f},{r['t_plus']:.4f}\n")
    
    # Coordination numbers
    with open(os.path.join(fig_dir, 'coordination_numbers.csv'), 'w') as f:
        f.write("pair,concentration,coord_number,pmf_min\n")
        for pair, pair_results in solv_results.items():
            for r in pair_results:
                f.write(f"{pair},{r['concentration']:.2f},{r['coord_number']:.4f},"
                        f"{r['pmf_min']:.4f}\n")
    
    # Anomalous exponents
    with open(os.path.join(fig_dir, 'anomalous_transport.csv'), 'w') as f:
        f.write("concentration,alpha_short,alpha_long,tau_crossover\n")
        for r in anom_results:
            f.write(f"{r['concentration']:.2f},{r['alpha_short']:.4f},"
                    f"{r['alpha_long']:.4f},{r['tau_crossover']:.2f}\n")


if __name__ == '__main__':
    main()
