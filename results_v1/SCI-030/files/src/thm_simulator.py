"""
Thermo-Hydro-Mechanical (THM) Coupling for EGS Reservoir
Finite Difference implementation for coupled heat-flow-stress analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import pandas as pd
import json
import os

np.random.seed(42)


class THMParameters:
    """Physical parameters for Kakkonda supercritical EGS system."""

    def __init__(self):
        # Rock properties (granite/granodiorite at Kakkonda)
        self.rho_rock = 2700         # kg/m^3
        self.cp_rock = 900           # J/(kg·K)
        self.k_rock = 2.5            # W/(m·K) thermal conductivity
        self.E_modulus = 50e9        # Pa (Young's modulus)
        self.nu_poisson = 0.25       # Poisson's ratio
        self.phi = 0.02              # porosity
        self.alpha_biot = 0.7        # Biot coefficient
        self.alpha_T = 8e-6          # thermal expansion coefficient /K
        self.permeability = 1e-16    # m^2 (matrix)

        # Fluid properties (supercritical water near critical point)
        self.rho_fluid = 400         # kg/m^3 (supercritical state)
        self.cp_fluid = 5000         # J/(kg·K) (enhanced near critical point)
        self.mu_fluid = 5e-5         # Pa·s viscosity
        self.k_fluid = 0.5           # W/(m·K)
        self.beta_fluid = 5e-4       # compressibility /Pa

        # Initial/boundary conditions (Kakkonda)
        self.T_initial = 380         # °C (supercritical zone)
        self.P_initial = 30e6        # Pa (~300 bar at 3500 m)
        self.T_surface = 15          # °C
        self.geothermal_gradient = 0.1  # °C/m (Kakkonda: very high)
        self.stress_horizontal = 80e6   # Pa (min horizontal stress)
        self.stress_vertical = 90e6     # Pa (overburden ~2700 kg/m^3 * 9.81 * 3500)

        # Well parameters
        self.T_injection = 30        # °C (cold water injection)
        self.P_injection = 35e6      # Pa
        self.P_production = 25e6     # Pa
        self.Q_injection = 0.05      # m^3/s per well


class THMSimulator:
    """
    2D Finite Difference THM coupled simulator.
    Solves: energy eq + Darcy flow + poroelastic stress
    """

    def __init__(self, params: THMParameters, nx=50, ny=50,
                 Lx=500, Ly=500, dt=86400):
        self.p = params
        self.nx = nx
        self.ny = ny
        self.Lx = Lx
        self.Ly = Ly
        self.dx = Lx / nx
        self.dy = Ly / ny
        self.dt = dt  # seconds (1 day default)

        self.x = np.linspace(self.dx/2, Lx - self.dx/2, nx)
        self.y = np.linspace(self.dy/2, Ly - self.dy/2, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        # State variables
        self.T = np.full((ny, nx), params.T_initial)   # temperature field [°C]
        self.P = np.full((ny, nx), params.P_initial)   # pressure field [Pa]
        self.sigma_eff = np.full((ny, nx), params.stress_horizontal)  # effective stress

        # Well positions (indices)
        self.wells = self._place_wells()
        self._apply_initial_conditions()

    def _place_wells(self):
        """Optimal doublet + 1 well configuration."""
        wells = {
            'inj1': (self.ny // 2, int(0.3 * self.nx)),   # injection 1
            'inj2': (self.ny // 2, int(0.7 * self.nx)),   # injection 2
            'prod': (self.ny // 2, self.nx // 2)           # production center
        }
        return wells

    def _apply_initial_conditions(self):
        """Apply depth-dependent initial conditions."""
        p = self.p
        for j in range(self.ny):
            depth_offset = (j - self.ny // 2) * self.dy
            self.T[j, :] = p.T_initial + p.geothermal_gradient * depth_offset * 0.01
            self.P[j, :] = p.P_initial + 9810 * depth_offset * 0.1

    def _compute_thermal_diffusivity(self):
        """Effective thermal diffusivity (rock + fluid)."""
        phi = self.p.phi
        k_eff = (1 - phi) * self.p.k_rock + phi * self.p.k_fluid
        rho_cp = (1 - phi) * self.p.rho_rock * self.p.cp_rock + \
                  phi * self.p.rho_fluid * self.p.cp_fluid
        return k_eff / rho_cp

    def _solve_heat_transport(self):
        """Solve energy equation with advection-diffusion."""
        alpha = self._compute_thermal_diffusivity()
        dt, dx, dy = self.dt, self.dx, self.dy

        T_new = self.T.copy()

        # Darcy velocity (from pressure gradient)
        dPdx = np.gradient(self.P, dx, axis=1)
        dPdy = np.gradient(self.P, dy, axis=0)
        vx = -self.p.permeability / self.p.mu_fluid * dPdx
        vy = -self.p.permeability / self.p.mu_fluid * dPdy

        # Diffusion (central difference)
        T_xx = (np.roll(self.T, -1, 1) - 2*self.T + np.roll(self.T, 1, 1)) / dx**2
        T_yy = (np.roll(self.T, -1, 0) - 2*self.T + np.roll(self.T, 1, 0)) / dy**2

        # Advection (upwind scheme)
        dTdx = np.where(vx > 0,
                        (self.T - np.roll(self.T, 1, 1)) / dx,
                        (np.roll(self.T, -1, 1) - self.T) / dx)
        dTdy = np.where(vy > 0,
                        (self.T - np.roll(self.T, 1, 0)) / dy,
                        (np.roll(self.T, -1, 0) - self.T) / dy)

        T_new = self.T + dt * (alpha * (T_xx + T_yy) -
                                vx * dTdx - vy * dTdy)

        # Apply injection temperature (Dirichlet at wells)
        for wtype, (iy, ix) in self.wells.items():
            if 'inj' in wtype:
                T_new[iy, ix] = self.p.T_injection
            elif 'prod' in wtype:
                T_new[iy, ix] = np.mean(T_new[iy-2:iy+2, ix-2:ix+2])

        # Boundary conditions (Neumann - no heat flux at sides)
        T_new[:, 0] = T_new[:, 1]
        T_new[:, -1] = T_new[:, -2]
        T_new[0, :] = T_new[1, :]
        T_new[-1, :] = T_new[-2, :]

        return T_new

    def _solve_flow(self):
        """Solve Darcy flow with pressure boundary conditions."""
        P_new = self.P.copy()
        dt, dx, dy = self.dt, self.dx, self.dy

        # Fluid compressibility
        S = self.p.phi * self.p.beta_fluid

        P_xx = (np.roll(self.P, -1, 1) - 2*self.P + np.roll(self.P, 1, 1)) / dx**2
        P_yy = (np.roll(self.P, -1, 0) - 2*self.P + np.roll(self.P, 1, 0)) / dy**2

        k_mu = self.p.permeability / self.p.mu_fluid
        dTdt = (self.T - self._T_prev) / dt if hasattr(self, '_T_prev') else 0

        P_new = self.P + dt * k_mu / S * (P_xx + P_yy) - \
                dt * self.p.alpha_biot / S * self.p.rho_fluid * dTdt

        # Well pressure BCs
        for wtype, (iy, ix) in self.wells.items():
            if 'inj' in wtype:
                P_new[iy, ix] = self.p.P_injection
            elif 'prod' in wtype:
                P_new[iy, ix] = self.p.P_production

        # Boundary (no-flow)
        P_new[:, 0] = P_new[:, 1]
        P_new[:, -1] = P_new[:, -2]
        P_new[0, :] = P_new[1, :]
        P_new[-1, :] = P_new[-2, :]

        return P_new

    def _update_effective_stress(self):
        """Update effective stress: sigma_eff = sigma_total - alpha*P."""
        # Thermal stress contribution
        dT = self.T - self.p.T_initial
        sigma_thermal = -self.p.E_modulus * self.p.alpha_T * dT / (1 - self.p.nu_poisson)

        self.sigma_eff = (self.p.stress_horizontal -
                           self.p.alpha_biot * (self.P - self.p.P_initial) +
                           sigma_thermal)

    def run(self, n_years=5, output_interval_days=365):
        """Run THM simulation for specified duration."""
        n_steps_total = int(n_years * 365)
        output_every = output_interval_days
        self._T_prev = self.T.copy()

        history = []
        snapshots = {}

        prod_iy, prod_ix = self.wells['prod']

        print(f"  Running THM simulation: {n_years} years, "
              f"{n_steps_total} daily steps...")

        for step in range(n_steps_total):
            self._T_prev = self.T.copy()
            self.T = self._solve_heat_transport()
            self.P = self._solve_flow()
            self._update_effective_stress()

            if step % output_every == 0:
                year = step / 365
                T_prod = self.T[prod_iy, prod_ix]
                P_prod = self.P[prod_iy, prod_ix]
                T_mean = np.mean(self.T)
                sigma_mean = np.mean(self.sigma_eff)

                # Thermal power [MW]
                Q = self.p.Q_injection * 2  # two injectors
                dT_fluid = T_prod - self.p.T_injection
                power_MW = Q * self.p.rho_fluid * self.p.cp_fluid * dT_fluid / 1e6

                history.append({
                    'year': year,
                    'T_production_C': float(T_prod),
                    'P_production_MPa': float(P_prod / 1e6),
                    'T_reservoir_mean_C': float(T_mean),
                    'sigma_eff_mean_MPa': float(sigma_mean / 1e6),
                    'thermal_power_MW': float(max(0, power_MW)),
                    'dT_production_C': float(dT_fluid)
                })

                if int(year) in [0, 1, 3, 5]:
                    snapshots[f'year_{int(year)}'] = {
                        'T': self.T.copy(),
                        'P': self.P.copy(),
                        'sigma_eff': self.sigma_eff.copy()
                    }

        self.history = pd.DataFrame(history)
        self.snapshots = snapshots
        return self

    def plot_snapshots(self, save_path=None):
        """Plot temperature and pressure snapshots."""
        years = [k for k in self.snapshots.keys()][:4]
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))

        for col, year_key in enumerate(years):
            snap = self.snapshots[year_key]
            yr = year_key.replace('year_', '')

            # Temperature
            ax = axes[0, col]
            im = ax.contourf(self.X, self.Y, snap['T'],
                             levels=20, cmap='inferno')
            plt.colorbar(im, ax=ax, label='Temperature [°C]')
            ax.set_title(f'Temperature - Year {yr}')
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')

            # Add well markers
            for wtype, (iy, ix) in self.wells.items():
                marker = '^' if 'inj' in wtype else 'v'
                ax.plot(self.x[ix], self.y[iy], marker=marker,
                        markersize=10, color='white', zorder=10)

            # Pressure
            ax = axes[1, col]
            im2 = ax.contourf(self.X, self.Y, snap['P'] / 1e6,
                              levels=20, cmap='viridis')
            plt.colorbar(im2, ax=ax, label='Pressure [MPa]')
            ax.set_title(f'Pressure - Year {yr}')
            ax.set_xlabel('X [m]')

        plt.suptitle('THM Coupled Simulation - Kakkonda EGS Reservoir', fontsize=14)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()

    def plot_history(self, save_path=None):
        """Plot production history."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        h = self.history

        axes[0, 0].plot(h['year'], h['T_production_C'], 'r-', linewidth=2)
        axes[0, 0].axhline(y=374.15, color='purple', linestyle='--',
                            label='Critical temp (374°C)')
        axes[0, 0].set_xlabel('Time [years]')
        axes[0, 0].set_ylabel('Production Temperature [°C]')
        axes[0, 0].set_title('Production Well Temperature')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].plot(h['year'], h['P_production_MPa'], 'b-', linewidth=2)
        axes[0, 1].axhline(y=22.064, color='orange', linestyle='--',
                            label='Critical pressure (22.1 MPa)')
        axes[0, 1].set_xlabel('Time [years]')
        axes[0, 1].set_ylabel('Production Pressure [MPa]')
        axes[0, 1].set_title('Production Well Pressure')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].plot(h['year'], h['thermal_power_MW'], 'g-', linewidth=2)
        axes[1, 0].fill_between(h['year'], h['thermal_power_MW'], alpha=0.3, color='green')
        axes[1, 0].set_xlabel('Time [years]')
        axes[1, 0].set_ylabel('Thermal Power [MW]')
        axes[1, 0].set_title('Thermal Power Output')
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].plot(h['year'], h['sigma_eff_mean_MPa'], 'm-', linewidth=2)
        axes[1, 1].set_xlabel('Time [years]')
        axes[1, 1].set_ylabel('Mean Effective Stress [MPa]')
        axes[1, 1].set_title('Reservoir Effective Stress Evolution')
        axes[1, 1].grid(True, alpha=0.3)

        plt.suptitle('THM Simulation - Production History (5-year preview)', fontsize=14)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    params = THMParameters()
    sim = THMSimulator(params, nx=40, ny=40, Lx=500, Ly=500)
    sim.run(n_years=5, output_interval_days=30)
    sim.plot_snapshots(save_path='figures/thm_snapshots.png')
    sim.plot_history(save_path='figures/thm_history.png')
    sim.history.to_csv('results/thm_history.csv', index=False)

    print("THM simulation completed.")
    final = sim.history.iloc[-1]
    print(f"  Final production temperature: {final['T_production_C']:.1f} °C")
    print(f"  Final thermal power: {final['thermal_power_MW']:.2f} MW")
