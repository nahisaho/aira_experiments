"""
30-Year Heat Recovery Simulation and Optimal Well Placement
for Kakkonda / Tohoku Supercritical EGS
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import differential_evolution
import pandas as pd
import json
import os

np.random.seed(42)


class HeatRecoverySimulator:
    """
    Long-term (30-year) heat recovery simulation for EGS.
    Uses semi-analytical doublet model with thermal drawdown.
    """

    def __init__(self):
        # Reservoir parameters (Kakkonda 3500m zone)
        self.T_reservoir = 380.0   # °C initial
        self.T_injection = 30.0    # °C injection temperature
        self.T_surface = 15.0      # °C surface

        # Rock properties
        self.rho_rock = 2700       # kg/m^3
        self.cp_rock = 900         # J/(kg·K)
        self.k_rock = 2.5          # W/(m·K)
        self.phi = 0.02            # porosity

        # Fluid (supercritical steam near wellbore)
        self.rho_fluid = 400       # kg/m^3
        self.cp_fluid = 4000       # J/(kg·K)
        self.mu_fluid = 5e-5       # Pa·s

        # System geometry
        self.H_reservoir = 200     # m (reservoir thickness)
        self.permeability = 1e-15  # m^2 (enhanced by EGS stimulation)

        # Power conversion
        self.eta_thermal = 0.18    # thermal-to-electric efficiency (supercritical)
        self.eta_pump = 0.75       # pump efficiency

        # Economic parameters
        self.drilling_cost_per_m = 8000   # JPY/m → $/m approx
        self.capex_per_MW = 5e6           # $/MW installed
        self.opex_per_MWh = 0.05          # $/kWh operation

    def thermal_drawdown_model(self, Q_m3s, L_well_m, t_years):
        """
        Schulz (1990) / Sanyal & Butler (2005) thermal drawdown model.
        Accounts for heat extraction vs. rock heat capacity.
        """
        t_sec = t_years * 365.25 * 86400

        # Thermal front velocity
        v_thermal = Q_m3s * self.rho_fluid * self.cp_fluid / \
                    (self.H_reservoir * self.rho_rock * self.cp_rock * L_well_m)

        # Breakthrough time (thermal)
        t_break = L_well_m / v_thermal
        t_break = max(t_break, 1e-10)

        # Temperature drawdown factor
        # After breakthrough: exponential decay
        if t_sec < t_break:
            eta_T = 1.0  # full reservoir temperature
        else:
            decay = np.exp(-0.15 * (t_sec - t_break) / t_break)
            eta_T = 0.85 + 0.15 * decay

        T_prod = self.T_injection + eta_T * (self.T_reservoir - self.T_injection)
        return T_prod, eta_T

    def compute_30year_production(self, well_config, dt_years=0.5):
        """
        Compute 30-year production history for a given well configuration.
        
        well_config: dict with keys:
          - 'n_doublets': int
          - 'well_spacing': float (m)
          - 'Q_per_well': float (m^3/s)
          - 'depth': float (m)
        """
        t_arr = np.arange(0, 30 + dt_years, dt_years)
        results = []

        Q_total = well_config['n_doublets'] * well_config['Q_per_well']
        L = well_config['well_spacing']
        n_d = well_config['n_doublets']

        for t in t_arr:
            T_prod, eta_T = self.thermal_drawdown_model(
                well_config['Q_per_well'], L, t + 0.001)

            # Thermal power [MW]
            dT = T_prod - self.T_injection
            P_thermal = Q_total * self.rho_fluid * self.cp_fluid * dT / 1e6

            # Electric power [MWe]
            # Supercritical steam → higher efficiency at T > 374°C
            eta = self.eta_thermal * (1 + 0.2 * (T_prod - 300) / 100) \
                  if T_prod > 300 else self.eta_thermal * 0.5
            eta = min(eta, 0.25)
            P_electric = P_thermal * eta

            # Pumping power
            dP_pump = 10e6  # Pa differential
            P_pump = Q_total * dP_pump / (1e6 * self.eta_pump)

            # Net electric power
            P_net = max(0, P_electric - P_pump)

            # Cumulative energy
            E_cumulative = P_net * t * 8760  # MWh (assuming 8760 h/yr)

            results.append({
                'year': t,
                'T_production_C': T_prod,
                'eta_thermal': eta_T,
                'P_thermal_MW': P_thermal,
                'P_electric_MW': P_electric,
                'P_pump_MW': P_pump,
                'P_net_MW': P_net,
                'E_cumulative_GWh': E_cumulative / 1000,
                'Q_total_m3s': Q_total
            })

        return pd.DataFrame(results)

    def optimize_well_placement(self, domain_size=(500, 500),
                                 n_doublets=2):
        """
        Optimize well positions to maximize 30-year energy recovery.
        Uses differential evolution optimization.
        """
        Lx, Ly = domain_size

        def objective(params):
            """Negative 30-year energy (to minimize = maximize energy)."""
            Q = params[0]
            spacing = params[1]
            n_d = int(params[2])

            if spacing < 100 or Q < 0.01:
                return 1e10

            config = {
                'n_doublets': n_d,
                'well_spacing': spacing,
                'Q_per_well': Q,
                'depth': 3500
            }
            hist = self.compute_30year_production(config)
            return -hist['E_cumulative_GWh'].iloc[-1]

        bounds = [
            (0.01, 0.15),    # Q per well [m^3/s]
            (100, 500),       # well spacing [m]
            (1, 4)            # n_doublets (integer)
        ]

        result = differential_evolution(objective, bounds, maxiter=50,
                                         seed=42, tol=0.01, workers=1)

        optimal = {
            'Q_per_well_m3s': float(result.x[0]),
            'well_spacing_m': float(result.x[1]),
            'n_doublets': int(result.x[2]),
            'max_energy_GWh': float(-result.fun),
            'optimization_success': bool(result.success)
        }
        return optimal

    def well_placement_grid_search(self, n_scenarios=20):
        """Grid search over well configurations."""
        scenarios = []
        for spacing in [150, 250, 350, 450]:
            for Q in [0.03, 0.06, 0.10]:
                for n_d in [1, 2, 3]:
                    config = {
                        'n_doublets': n_d,
                        'well_spacing': spacing,
                        'Q_per_well': Q,
                        'depth': 3500
                    }
                    hist = self.compute_30year_production(config)
                    final = hist.iloc[-1]
                    scenarios.append({
                        'n_doublets': n_d,
                        'spacing_m': spacing,
                        'Q_m3s': Q,
                        'E_30yr_GWh': final['E_cumulative_GWh'],
                        'P_net_final_MW': final['P_net_MW'],
                        'T_prod_final_C': final['T_production_C'],
                        'P_thermal_final_MW': final['P_thermal_MW']
                    })

        return pd.DataFrame(scenarios).sort_values('E_30yr_GWh', ascending=False)

    def plot_results(self, save_path=None):
        """Comprehensive 30-year production visualization."""
        fig = plt.figure(figsize=(20, 16))
        fig.suptitle('30-Year Heat Recovery Analysis - Kakkonda Supercritical EGS',
                     fontsize=14, fontweight='bold')

        # Scenario comparison
        scenarios_df = self.well_placement_grid_search()

        # Best 3 scenarios
        top3_configs = scenarios_df.head(3)

        ax1 = fig.add_subplot(3, 3, 1)
        colors_sc = ['#E74C3C', '#2ECC71', '#3498DB']
        for idx, (_, sc) in enumerate(top3_configs.iterrows()):
            config = {
                'n_doublets': sc['n_doublets'],
                'well_spacing': sc['spacing_m'],
                'Q_per_well': sc['Q_m3s'],
                'depth': 3500
            }
            hist = self.compute_30year_production(config)
            label = f"N={sc['n_doublets']}, L={sc['spacing_m']:.0f}m, Q={sc['Q_m3s']:.2f}"
            ax1.plot(hist['year'], hist['T_production_C'],
                     color=colors_sc[idx], linewidth=2, label=label)

        ax1.axhline(y=374.15, color='purple', linestyle='--',
                    label='Critical Temp (374°C)', alpha=0.7)
        ax1.set_xlabel('Time [years]')
        ax1.set_ylabel('Production Temperature [°C]')
        ax1.set_title('Temperature Drawdown - Top Scenarios')
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)

        # Net power output
        ax2 = fig.add_subplot(3, 3, 2)
        for idx, (_, sc) in enumerate(top3_configs.iterrows()):
            config = {
                'n_doublets': sc['n_doublets'],
                'well_spacing': sc['spacing_m'],
                'Q_per_well': sc['Q_m3s'],
                'depth': 3500
            }
            hist = self.compute_30year_production(config)
            label = f"Scenario {idx+1}"
            ax2.plot(hist['year'], hist['P_net_MW'],
                     color=colors_sc[idx], linewidth=2, label=label)

        ax2.set_xlabel('Time [years]')
        ax2.set_ylabel('Net Electric Power [MWe]')
        ax2.set_title('Net Electric Power Output')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        # Cumulative energy
        ax3 = fig.add_subplot(3, 3, 3)
        for idx, (_, sc) in enumerate(top3_configs.iterrows()):
            config = {
                'n_doublets': sc['n_doublets'],
                'well_spacing': sc['spacing_m'],
                'Q_per_well': sc['Q_m3s'],
                'depth': 3500
            }
            hist = self.compute_30year_production(config)
            label = f"Scenario {idx+1}"
            ax3.plot(hist['year'], hist['E_cumulative_GWh'],
                     color=colors_sc[idx], linewidth=2, label=label)

        ax3.set_xlabel('Time [years]')
        ax3.set_ylabel('Cumulative Energy [GWh]')
        ax3.set_title('Cumulative Energy Production')
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)

        # Grid search heatmap
        ax4 = fig.add_subplot(3, 3, 4)
        pivot = scenarios_df[scenarios_df['n_doublets'] == 2].pivot_table(
            values='E_30yr_GWh', index='spacing_m', columns='Q_m3s')
        im = ax4.imshow(pivot.values, aspect='auto', cmap='viridis',
                        origin='lower')
        plt.colorbar(im, ax=ax4, label='30-yr Energy [GWh]')
        ax4.set_xticks(range(len(pivot.columns)))
        ax4.set_xticklabels([f'{q:.2f}' for q in pivot.columns])
        ax4.set_yticks(range(len(pivot.index)))
        ax4.set_yticklabels(pivot.index.astype(int))
        ax4.set_xlabel('Flow Rate [m³/s]')
        ax4.set_ylabel('Well Spacing [m]')
        ax4.set_title('30-yr Energy Heatmap (N=2 doublets)')

        # Optimal well layout diagram
        ax5 = fig.add_subplot(3, 3, 5)
        best = top3_configs.iloc[0]
        spacing = best['spacing_m']
        cx, cy = 250, 250

        # Injection wells
        inj_positions = [
            (cx - spacing/2, cy - spacing/4),
            (cx + spacing/2, cy - spacing/4)
        ]
        # Production well
        prod_position = (cx, cy + spacing/4)

        domain = plt.Rectangle((0, 0), 500, 500,
                                 fill=False, edgecolor='black', linewidth=2)
        ax5.add_patch(domain)

        for i, (ix, iy) in enumerate(inj_positions):
            circle = plt.Circle((ix, iy), 15, color='#3498DB', alpha=0.8)
            ax5.add_patch(circle)
            ax5.annotate(f'INJ-{i+1}', (ix, iy), ha='center', va='center',
                         fontsize=8, color='white', fontweight='bold')

        circle_prod = plt.Circle((prod_position[0], prod_position[1]),
                                  20, color='#E74C3C', alpha=0.8)
        ax5.add_patch(circle_prod)
        ax5.annotate('PROD', prod_position, ha='center', va='center',
                     fontsize=8, color='white', fontweight='bold')

        # Flow lines
        for ix, iy in inj_positions:
            ax5.annotate('', xy=prod_position, xytext=(ix, iy),
                         arrowprops=dict(arrowstyle='->', color='blue',
                                         lw=1.5, alpha=0.5))

        ax5.set_xlim(0, 500)
        ax5.set_ylim(0, 500)
        ax5.set_xlabel('X [m]')
        ax5.set_ylabel('Y [m]')
        ax5.set_title(f'Optimal Well Layout\n(spacing={spacing:.0f}m, N={best["n_doublets"]})')
        ax5.set_aspect('equal')
        ax5.grid(True, alpha=0.3)

        # Heat recovery factor
        ax6 = fig.add_subplot(3, 3, 6)
        config_base = {
            'n_doublets': int(best['n_doublets']),
            'well_spacing': float(best['spacing_m']),
            'Q_per_well': float(best['Q_m3s']),
            'depth': 3500
        }
        hist_best = self.compute_30year_production(config_base)
        recovery_factor = hist_best['eta_thermal']
        ax6.plot(hist_best['year'], recovery_factor * 100, 'g-', linewidth=2)
        ax6.fill_between(hist_best['year'], recovery_factor * 100, alpha=0.3)
        ax6.set_xlabel('Time [years]')
        ax6.set_ylabel('Thermal Recovery Factor [%]')
        ax6.set_title('Heat Recovery Factor Evolution')
        ax6.grid(True, alpha=0.3)

        # Thermal vs electric breakdown
        ax7 = fig.add_subplot(3, 3, 7)
        ax7.stackplot(hist_best['year'],
                       hist_best['P_pump_MW'],
                       hist_best['P_electric_MW'] - hist_best['P_pump_MW'],
                       hist_best['P_thermal_MW'] - hist_best['P_electric_MW'],
                       labels=['Pump Load', 'Net Electric', 'Residual Thermal'],
                       colors=['#E74C3C', '#2ECC71', '#F39C12'], alpha=0.7)
        ax7.set_xlabel('Time [years]')
        ax7.set_ylabel('Power [MW]')
        ax7.set_title('Power Balance (Best Scenario)')
        ax7.legend(fontsize=8)
        ax7.grid(True, alpha=0.3)

        # Capacity factor vs depth
        ax8 = fig.add_subplot(3, 3, 8)
        depths = np.arange(2000, 5000, 200)
        T_res = 15 + 0.1 * depths  # 100°C/km gradient (Kakkonda)
        # Approximate capacity: higher T → more power, but deeper = higher cost
        P_capacity = np.clip((T_res - 200) * 0.05, 0, 20)
        cf = 0.92 - 0.01 * (depths - 3500) / 500  # capacity factor

        ax8.plot(depths, P_capacity, 'b-', linewidth=2, label='Thermal Power [MW]')
        ax8.set_xlabel('Depth [m]')
        ax8.set_ylabel('Estimated Capacity [MW]')
        ax8t = ax8.twinx()
        ax8t.plot(depths, cf, 'r--', linewidth=2, label='Capacity Factor')
        ax8t.set_ylabel('Capacity Factor')
        ax8.axvline(x=3500, color='green', linestyle='--', linewidth=2,
                    label='Kakkonda target depth')
        ax8.set_title('Depth vs Capacity Analysis')
        ax8.legend(loc='upper left', fontsize=8)
        ax8t.legend(loc='lower right', fontsize=8)
        ax8.grid(True, alpha=0.3)

        # LCOE comparison
        ax9 = fig.add_subplot(3, 3, 9)
        tech_names = ['Conv.\nGeothermal', 'EGS\nSubcritical', 'EGS\nSupercritical\n(Kakkonda)',
                       'Solar\nPV', 'Wind\nOnshore']
        lcoe_values = [80, 150, 120, 50, 45]
        colors_lcoe = ['#2ECC71', '#E67E22', '#E74C3C', '#F1C40F', '#3498DB']
        bars = ax9.bar(tech_names, lcoe_values, color=colors_lcoe, alpha=0.8)
        ax9.set_ylabel('LCOE [USD/MWh]')
        ax9.set_title('LCOE Comparison (2024 estimates)')
        ax9.bar_label(bars, fmt='%d', fontsize=9)
        ax9.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()

        return scenarios_df, hist_best


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    sim = HeatRecoverySimulator()
    scenarios_df, best_hist = sim.plot_results(save_path='figures/heat_recovery_30yr.png')

    scenarios_df.to_csv('results/well_placement_scenarios.csv', index=False)
    best_hist.to_csv('results/best_scenario_30yr_history.csv', index=False)

    print("30-year heat recovery simulation completed.")
    best = scenarios_df.iloc[0]
    print(f"  Optimal: N={best['n_doublets']}, spacing={best['spacing_m']:.0f}m, "
          f"Q={best['Q_m3s']:.2f} m³/s")
    print(f"  30-yr energy: {best['E_30yr_GWh']:.1f} GWh")
    print(f"  Final net power: {best['P_net_final_MW']:.1f} MWe")
