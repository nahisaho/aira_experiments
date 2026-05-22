"""
Kakkonda / Tohoku Case Study - Geological Model and Integration
Full EGS simulation workflow integration for Japanese supercritical geothermal system
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import pandas as pd
import json
import os
import sys

np.random.seed(42)

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


class KakkondaGeologicalModel:
    """
    Geological model for Kakkonda geothermal field, Iwate Prefecture, Japan.
    
    Key features:
    - WD-1a borehole reached 3729m (380°C, T > critical)
    - New Kakkonda Granite (NKG): young pluton (~0.2 Ma)
    - World's first supercritical geothermal well (1995)
    - Located in NE Japan arc, close to volcanic front
    - Stress: σH (E-W) > σh (N-S) ≈ σV
    """

    def __init__(self):
        self.layers = self._define_geology()
        self.stress_profile = self._define_stress()
        self.temperature_profile = self._define_temperature()

    def _define_geology(self):
        """Stratigraphic column for Kakkonda area."""
        return [
            {'unit': 'Quaternary Volcanics', 'top': 0, 'base': 300,
             'lithology': 'Volcanic tuff, lava', 'k': 1e-14, 'phi': 0.20,
             'rho': 2200, 'color': '#F0E68C'},
            {'unit': 'Tertiary Sediments', 'top': 300, 'base': 800,
             'lithology': 'Sandstone, mudstone', 'k': 1e-15, 'phi': 0.12,
             'rho': 2400, 'color': '#DEB887'},
            {'unit': 'Miocene Granodiorite', 'top': 800, 'base': 2000,
             'lithology': 'Granodiorite (altered)', 'k': 1e-16, 'phi': 0.05,
             'rho': 2650, 'color': '#CD853F'},
            {'unit': 'Old Kakkonda Granite', 'top': 2000, 'base': 3200,
             'lithology': 'Granite (fractured)', 'k': 5e-16, 'phi': 0.03,
             'rho': 2700, 'color': '#BC8F8F'},
            {'unit': 'New Kakkonda Granite (NKG)', 'top': 3200, 'base': 4500,
             'lithology': 'Young granite, supercritical zone', 'k': 1e-17, 'phi': 0.01,
             'rho': 2750, 'color': '#8B4513'},
        ]

    def _define_stress(self):
        """In-situ stress profile (MPa vs depth m)."""
        depths = np.arange(0, 5000, 100)
        # Vertical (lithostatic): ~27 MPa/km
        sv = 0.027 * depths
        # Max horizontal (E-W, tectonic): slightly > vertical
        sh_max = 0.028 * depths + 5
        # Min horizontal (N-S): slightly < vertical
        sh_min = 0.025 * depths

        return pd.DataFrame({
            'depth_m': depths,
            'Sv_MPa': sv,
            'SHmax_MPa': sh_max,
            'Shmin_MPa': sh_min,
            'PP_MPa': 0.0098 * depths  # hydrostatic pore pressure
        })

    def _define_temperature(self):
        """Temperature profile (WD-1a well data approximation)."""
        depths = np.arange(0, 4000, 50)
        T = np.zeros_like(depths, dtype=float)

        for i, d in enumerate(depths):
            if d < 1000:
                T[i] = 15 + 0.05 * d        # 50°C/km
            elif d < 2500:
                T[i] = 65 + 0.08 * (d - 1000)  # 80°C/km
            elif d < 3500:
                T[i] = 185 + 0.12 * (d - 2500)  # 120°C/km
            else:
                T[i] = 305 + 0.15 * (d - 3500)  # 150°C/km

        return pd.DataFrame({'depth_m': depths, 'T_C': T})

    def get_supercritical_zone_depth(self):
        """Find depth where T > 374.15°C (water critical temperature)."""
        temp = self.temperature_profile
        sc_rows = temp[temp['T_C'] >= 374.15]
        if len(sc_rows) > 0:
            return float(sc_rows.iloc[0]['depth_m'])
        return None

    def plot_geological_model(self, save_path=None):
        """Plot comprehensive geological model."""
        fig = plt.figure(figsize=(20, 14))
        fig.suptitle('Kakkonda Geothermal Field - Geological Model\n'
                     'Iwate Prefecture, NE Japan', fontsize=14, fontweight='bold')

        # Stratigraphic column
        ax1 = fig.add_subplot(1, 4, 1)
        for layer in self.layers:
            height = layer['base'] - layer['top']
            rect = mpatches.Rectangle(
                (0, -layer['base']), 1, height,
                facecolor=layer['color'], edgecolor='black', linewidth=1
            )
            ax1.add_patch(rect)
            mid_depth = -(layer['top'] + layer['base']) / 2
            ax1.text(0.5, mid_depth, layer['unit'], ha='center', va='center',
                     fontsize=7, fontweight='bold', wrap=True)

        ax1.set_xlim(0, 1)
        ax1.set_ylim(-4500, 0)
        ax1.set_yticks([-d for d in range(0, 5000, 500)])
        ax1.set_yticklabels([f'{d}m' for d in range(0, 5000, 500)])
        ax1.set_xticks([])
        ax1.set_title('Stratigraphy')
        ax1.axhline(y=-3729, color='red', linestyle='--', linewidth=2,
                    label='WD-1a TD (3729m)')
        sc_depth = self.get_supercritical_zone_depth()
        if sc_depth:
            ax1.axhline(y=-sc_depth, color='purple', linestyle=':',
                        linewidth=2, label=f'SC zone ({sc_depth:.0f}m)')
        ax1.legend(fontsize=7, loc='lower right')

        # Temperature profile
        ax2 = fig.add_subplot(1, 4, 2)
        temp = self.temperature_profile
        ax2.plot(temp['T_C'], -temp['depth_m'], 'r-', linewidth=2.5)
        ax2.axvline(x=374.15, color='purple', linestyle='--',
                    label='Tc = 374.15°C')
        ax2.axvline(x=100, color='orange', linestyle=':', label='100°C')
        ax2.axhline(y=-3729, color='red', linestyle='--', alpha=0.5)
        ax2.fill_betweenx(-temp['depth_m'],
                           np.where(temp['T_C'] >= 374.15, temp['T_C'], 374.15),
                           374.15,
                           where=temp['T_C'] >= 374.15,
                           alpha=0.3, color='purple', label='Supercritical zone')
        ax2.set_xlabel('Temperature [°C]')
        ax2.set_ylabel('Depth [m]')
        ax2.set_ylim(-4000, 0)
        ax2.set_title('Temperature Profile\n(WD-1a approximation)')
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)

        # Stress profile
        ax3 = fig.add_subplot(1, 4, 3)
        stress = self.stress_profile
        ax3.plot(stress['Sv_MPa'], -stress['depth_m'], 'k-',
                 linewidth=2, label='Sv (Vertical)')
        ax3.plot(stress['SHmax_MPa'], -stress['depth_m'], 'r-',
                 linewidth=2, label='SHmax (E-W)')
        ax3.plot(stress['Shmin_MPa'], -stress['depth_m'], 'b-',
                 linewidth=2, label='Shmin (N-S)')
        ax3.plot(stress['PP_MPa'], -stress['depth_m'], 'g--',
                 linewidth=1.5, label='Pore Pressure')
        ax3.axhline(y=-3500, color='purple', linestyle=':',
                    alpha=0.7, label='Target depth (3500m)')
        ax3.set_xlabel('Stress / Pressure [MPa]')
        ax3.set_ylabel('Depth [m]')
        ax3.set_ylim(-4000, 0)
        ax3.set_title('In-Situ Stress Profile')
        ax3.legend(fontsize=7)
        ax3.grid(True, alpha=0.3)

        # Permeability profile
        ax4 = fig.add_subplot(1, 4, 4)
        all_depths = []
        all_k = []
        for layer in self.layers:
            for d in range(layer['top'], min(layer['base'], 4001), 100):
                all_depths.append(d)
                # Add fracture-enhanced permeability in NKG
                k = layer['k']
                if layer['unit'] == 'New Kakkonda Granite (NKG)':
                    k *= np.random.lognormal(0, 0.5)  # fracture variability
                all_k.append(k)

        ax4.semilogx(all_k, [-d for d in all_depths], 'b-', linewidth=1.5)
        ax4.fill_betweenx([-d for d in all_depths], [1e-20] * len(all_k),
                           all_k, alpha=0.3, color='blue')
        ax4.set_xlabel('Permeability [m²]')
        ax4.set_ylabel('Depth [m]')
        ax4.set_xlim(1e-19, 1e-13)
        ax4.set_ylim(-4000, 0)
        ax4.set_title('Permeability Profile')
        ax4.axhline(y=-3500, color='purple', linestyle=':',
                    alpha=0.7, label='EGS target zone')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()

        return fig


class KakkondaIntegration:
    """
    Master integration class for the full Kakkonda EGS simulation.
    Runs all modules and compiles final results.
    """

    def __init__(self):
        self.geo_model = KakkondaGeologicalModel()
        self.results = {}

    def run_full_workflow(self, output_dir='.'):
        """Execute complete simulation workflow."""
        os.makedirs(f'{output_dir}/figures', exist_ok=True)
        os.makedirs(f'{output_dir}/results', exist_ok=True)
        os.makedirs(f'{output_dir}/data', exist_ok=True)
        os.makedirs(f'{output_dir}/logs', exist_ok=True)

        print("=" * 60)
        print("KAKKONDA EGS SUPERCRITICAL SIMULATION WORKFLOW")
        print("=" * 60)

        import datetime
        log_entries = []

        def log_event(phase, event, skill, files=None, status='ok'):
            entry = {
                'timestamp': datetime.datetime.utcnow().isoformat() + 'Z',
                'phase': phase,
                'event_type': event,
                'actor': 'co-scientist',
                'skill_or_tool': skill,
                'files_written': files or [],
                'status': status
            }
            log_entries.append(entry)

        log_event('SETUP', 'run_started', 'kakkonda_integration')

        # 1. Geological model
        print("\n[1/6] Building geological model...")
        self.geo_model.plot_geological_model(
            save_path=f'{output_dir}/figures/geological_model.png')
        stress = self.geo_model.stress_profile
        temp = self.geo_model.temperature_profile
        stress.to_csv(f'{output_dir}/data/stress_profile.csv', index=False)
        temp.to_csv(f'{output_dir}/data/temperature_profile.csv', index=False)
        sc_depth = self.geo_model.get_supercritical_zone_depth()
        if sc_depth is None:
            sc_depth = 3200.0  # default estimate
        print(f"  Supercritical zone depth: {sc_depth:.0f} m")
        log_event('GEOLOGY', 'handoff_completed', 'geological_model',
                  ['figures/geological_model.png', 'data/stress_profile.csv'])

        # 2. DFN model
        print("\n[2/6] Running DFN model...")
        from dfn_model import DFNModel
        dfn = DFNModel(domain_size=(500, 500), depth=3500)
        dfn.generate()
        wells = {
            'Injection_Well_1': (150, 250),
            'Injection_Well_2': (350, 250),
            'Production_Well': (250, 250)
        }
        dfn.plot(save_path=f'{output_dir}/figures/dfn_model.png', wells=wells)
        dfn_stats = dfn.compute_connectivity()
        k_tensor = dfn.get_permeability_tensor()
        dfn_stats['Kxx_m2'] = float(k_tensor[0, 0])
        dfn_stats['Kyy_m2'] = float(k_tensor[1, 1])
        dfn.fractures.to_csv(f'{output_dir}/results/dfn_fractures.csv', index=False)
        with open(f'{output_dir}/results/dfn_statistics.json', 'w') as f:
            json.dump(dfn_stats, f, indent=2)
        print(f"  {dfn_stats['n_fractures']} fractures, "
              f"Kxx={dfn_stats['Kxx_m2']:.2e} m²")
        log_event('DFN', 'handoff_completed', 'dfn_model',
                  ['figures/dfn_model.png', 'results/dfn_statistics.json'])
        self.results['dfn'] = dfn_stats

        # 3. EOS properties
        print("\n[3/6] Computing supercritical water EOS...")
        from supercritical_eos import SupercriticalWaterEOS
        eos = SupercriticalWaterEOS()
        eos.plot_properties(save_path=f'{output_dir}/figures/eos_properties.png')
        eos_summary = eos.summary_at_kakkonda()
        eos_summary.to_csv(f'{output_dir}/results/eos_kakkonda_summary.csv', index=False)
        kakkonda_state = eos_summary[eos_summary['T_C'] == 380].iloc[0]
        print(f"  At 380°C / 30 MPa: rho={kakkonda_state['density_kg_m3']:.1f} kg/m³, "
              f"cp={kakkonda_state['cp_kJ_kgK']:.2f} kJ/(kg·K)")
        log_event('EOS', 'handoff_completed', 'supercritical_eos',
                  ['figures/eos_properties.png', 'results/eos_kakkonda_summary.csv'])
        self.results['eos'] = eos_summary.to_dict('records')

        # 4. THM simulation
        print("\n[4/6] Running THM coupled simulation (5-year preview)...")
        from thm_simulator import THMParameters, THMSimulator
        params = THMParameters()
        thm = THMSimulator(params, nx=40, ny=40, Lx=500, Ly=500)
        thm.run(n_years=5, output_interval_days=30)
        thm.plot_snapshots(save_path=f'{output_dir}/figures/thm_snapshots.png')
        thm.plot_history(save_path=f'{output_dir}/figures/thm_history.png')
        thm.history.to_csv(f'{output_dir}/results/thm_history.csv', index=False)
        final_thm = thm.history.iloc[-1]
        print(f"  Final T_prod={final_thm['T_production_C']:.1f}°C, "
              f"P_thermal={final_thm['thermal_power_MW']:.1f} MW")
        log_event('THM', 'handoff_completed', 'thm_simulator',
                  ['figures/thm_snapshots.png', 'figures/thm_history.png'])
        self.results['thm'] = final_thm.to_dict()

        # 5. Coulomb stress / seismic risk
        print("\n[5/6] Computing Coulomb stress and seismic risk...")
        from coulomb_stress import CoulombStressModel
        coulomb = CoulombStressModel()
        seismic_results = coulomb.plot_results(
            save_path=f'{output_dir}/figures/coulomb_stress.png')
        with open(f'{output_dir}/results/seismic_risk_analysis.json', 'w') as f:
            json.dump(seismic_results, f, indent=2, default=float)
        tl = seismic_results['traffic_light']
        print(f"  Traffic Light: {tl['status']}, Mw_max={seismic_results['Mw_max_credible']:.1f}")
        log_event('SEISMIC', 'handoff_completed', 'coulomb_stress',
                  ['figures/coulomb_stress.png', 'results/seismic_risk_analysis.json'])
        self.results['seismic'] = seismic_results

        # 6. 30-year heat recovery
        print("\n[6/6] Running 30-year heat recovery optimization...")
        from heat_recovery import HeatRecoverySimulator
        hr = HeatRecoverySimulator()
        scenarios_df, best_hist = hr.plot_results(
            save_path=f'{output_dir}/figures/heat_recovery_30yr.png')
        scenarios_df.to_csv(f'{output_dir}/results/well_placement_scenarios.csv', index=False)
        best_hist.to_csv(f'{output_dir}/results/best_scenario_30yr_history.csv', index=False)
        best_sc = scenarios_df.iloc[0]
        final_30yr = best_hist.iloc[-1]
        print(f"  Best: {best_sc['E_30yr_GWh']:.0f} GWh (30yr), "
              f"{final_30yr['P_net_MW']:.1f} MWe final")
        log_event('RECOVERY', 'handoff_completed', 'heat_recovery',
                  ['figures/heat_recovery_30yr.png', 'results/well_placement_scenarios.csv'])
        self.results['heat_recovery'] = {
            'best_30yr_GWh': float(best_sc['E_30yr_GWh']),
            'n_doublets': int(best_sc['n_doublets']),
            'spacing_m': float(best_sc['spacing_m']),
            'Q_m3s': float(best_sc['Q_m3s']),
            'final_net_MW': float(final_30yr['P_net_MW'])
        }

        # Summary visualization
        print("\n  Generating summary figure...")
        self._plot_summary(output_dir, sc_depth)
        log_event('SUMMARY', 'report_finalized', 'kakkonda_integration',
                  ['figures/kakkonda_summary.png'])

        log_event('COMPLETE', 'run_completed', 'kakkonda_integration', status='ok')

        # Save log
        with open(f'{output_dir}/logs/process-log.jsonl', 'w') as f:
            for entry in log_entries:
                f.write(json.dumps(entry) + '\n')

        print("\n" + "=" * 60)
        print("SIMULATION WORKFLOW COMPLETED")
        print("=" * 60)

        return self.results

    def _plot_summary(self, output_dir, sc_depth):
        """Summary dashboard figure."""
        fig = plt.figure(figsize=(18, 10))
        fig.patch.set_facecolor('#1a1a2e')

        fig.suptitle('KAKKONDA SUPERCRITICAL EGS — SIMULATION SUMMARY DASHBOARD',
                     fontsize=16, fontweight='bold', color='white', y=0.98)

        metrics = [
            ('Supercritical\nZone Depth', f'{sc_depth:.0f} m', '#E74C3C'),
            ('Reservoir\nTemperature', '380 °C', '#E67E22'),
            ('30-yr Energy\n(Best Config)', f"{self.results['heat_recovery']['best_30yr_GWh']:.0f} GWh", '#2ECC71'),
            ('Final Net\nPower', f"{self.results['heat_recovery']['final_net_MW']:.1f} MWe", '#3498DB'),
            ('Total\nFractures', f"{self.results['dfn']['n_fractures']}", '#9B59B6'),
            ('Seismic Risk\nStatus', self.results['seismic']['traffic_light']['status'], '#1ABC9C'),
        ]

        for i, (label, value, color) in enumerate(metrics):
            ax = fig.add_axes([0.02 + i * 0.16, 0.72, 0.14, 0.22])
            ax.set_facecolor('#16213e')
            ax.add_patch(plt.Rectangle((0.05, 0.05), 0.9, 0.9,
                                        facecolor=color, alpha=0.2,
                                        transform=ax.transAxes))
            ax.text(0.5, 0.65, value, ha='center', va='center',
                    fontsize=18, fontweight='bold', color=color,
                    transform=ax.transAxes)
            ax.text(0.5, 0.2, label, ha='center', va='center',
                    fontsize=9, color='lightgray',
                    transform=ax.transAxes)
            ax.axis('off')

        # Mini temperature profile
        ax_t = fig.add_axes([0.02, 0.08, 0.28, 0.58])
        ax_t.set_facecolor('#16213e')
        temp = self.geo_model.temperature_profile
        ax_t.plot(temp['T_C'], -temp['depth_m'], 'r-', linewidth=2.5)
        ax_t.axvline(x=374.15, color='purple', linestyle='--',
                      label='Tc = 374.15°C', linewidth=1.5)
        ax_t.fill_betweenx(-temp['depth_m'],
                             np.where(temp['T_C'] >= 374.15, temp['T_C'], 374.15),
                             374.15,
                             where=temp['T_C'] >= 374.15,
                             alpha=0.4, color='purple')
        ax_t.set_xlabel('Temperature [°C]', color='white')
        ax_t.set_ylabel('Depth [m]', color='white')
        ax_t.set_title('Temperature Profile', color='white')
        ax_t.tick_params(colors='white')
        ax_t.legend(fontsize=8, facecolor='#16213e', labelcolor='white')
        ax_t.grid(True, alpha=0.2, color='gray')
        for spine in ax_t.spines.values():
            spine.set_color('gray')

        # Mini 30-yr production (if available)
        ax_p = fig.add_axes([0.36, 0.08, 0.28, 0.58])
        ax_p.set_facecolor('#16213e')
        try:
            best_hist = pd.read_csv(f'{output_dir}/results/best_scenario_30yr_history.csv')
            ax_p.fill_between(best_hist['year'], best_hist['P_net_MW'],
                               alpha=0.4, color='#2ECC71')
            ax_p.plot(best_hist['year'], best_hist['P_net_MW'],
                       'g-', linewidth=2.5)
            ax_p.set_xlabel('Time [years]', color='white')
            ax_p.set_ylabel('Net Electric Power [MWe]', color='white')
            ax_p.set_title('30-Year Power Production', color='white')
        except Exception:
            ax_p.text(0.5, 0.5, 'Production data\nnot available',
                       ha='center', va='center', transform=ax_p.transAxes,
                       color='white')
        ax_p.tick_params(colors='white')
        ax_p.grid(True, alpha=0.2, color='gray')
        for spine in ax_p.spines.values():
            spine.set_color('gray')

        # Mini CFF evolution
        ax_s = fig.add_axes([0.70, 0.08, 0.28, 0.58])
        ax_s.set_facecolor('#16213e')
        try:
            cev = pd.DataFrame(self.results['seismic']['cff_evolution'])
            ax_s.semilogy(cev['time_days'], cev['CFF_max_MPa'],
                           'r-o', linewidth=2.5, markersize=6)
            ax_s.axhline(y=0.1, color='#2ECC71', linestyle='--', label='Green')
            ax_s.axhline(y=0.5, color='orange', linestyle='--', label='Yellow')
            ax_s.axhline(y=1.0, color='red', linestyle='--', label='Orange')
            ax_s.set_xlabel('Time [days]', color='white')
            ax_s.set_ylabel('Max ΔCFF [MPa]', color='white')
            ax_s.set_title('Seismic Risk Evolution', color='white')
            ax_s.legend(fontsize=7, facecolor='#16213e', labelcolor='white')
        except Exception:
            pass
        ax_s.tick_params(colors='white')
        ax_s.grid(True, alpha=0.2, color='gray')
        for spine in ax_s.spines.values():
            spine.set_color('gray')

        plt.savefig(f'{output_dir}/figures/kakkonda_summary.png',
                    dpi=150, bbox_inches='tight', facecolor='#1a1a2e')
        plt.close()


if __name__ == '__main__':
    integration = KakkondaIntegration()
    results = integration.run_full_workflow(output_dir='.')
    print("\nAll results saved to workspace/")
