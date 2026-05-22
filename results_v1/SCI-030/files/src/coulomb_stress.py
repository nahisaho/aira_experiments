"""
Coulomb Stress Change Modeling for Induced Seismicity Risk Assessment
EGS-induced seismicity at Kakkonda / Tohoku region
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import pandas as pd
import json
import os

np.random.seed(42)


class CoulombStressModel:
    """
    Coulomb Failure Function (CFF) analysis for EGS-induced seismicity.
    
    CFF = tau + mu * (sigma_n - P)
    where:
      tau     = shear stress change [MPa]
      mu      = friction coefficient (~0.6 for granite)
      sigma_n = normal stress change [MPa]
      P       = pore pressure change [MPa]
    
    Positive CFF → increased failure potential → seismicity risk
    """

    def __init__(self, mu_friction=0.6, lambda_skempton=0.5):
        self.mu = mu_friction         # friction coefficient
        self.lambda_s = lambda_skempton  # Skempton's coefficient

        # Kakkonda regional stress field
        # NE Japan: compressional regime (E-W direction)
        self.sigma1 = 90.0   # MPa (vertical, overburden at 3500m)
        self.sigma2 = 85.0   # MPa (max horizontal: E-W)
        self.sigma3 = 80.0   # MPa (min horizontal: N-S)

        # Regional fault orientations (Kakkonda area)
        self.fault_sets = [
            {'name': 'NNW-SSE_strike_slip', 'strike': 160, 'dip': 80, 'rake': 0},
            {'name': 'ENE-WSW_normal', 'strike': 70, 'dip': 60, 'rake': -90},
            {'name': 'EW_thrust', 'strike': 90, 'dip': 30, 'rake': 90},
        ]

    def compute_stress_change_from_injection(self, nx=80, ny=80,
                                              Lx=1000, Ly=1000,
                                              well_positions=None,
                                              dP_injection=10e6,
                                              time_days=365):
        """
        Compute pore pressure-induced stress changes from fluid injection.
        Uses Eshelby-type analytical solution for pressure diffusion.
        """
        x = np.linspace(-Lx/2, Lx/2, nx)
        y = np.linspace(-Ly/2, Ly/2, ny)
        X, Y = np.meshgrid(x, y)

        if well_positions is None:
            well_positions = [(-150, 0), (150, 0)]  # doublet

        # Hydraulic diffusivity (granite with fractures)
        D_hydraulic = 0.1  # m^2/s

        dP = np.zeros((ny, nx))
        for wx, wy in well_positions:
            r = np.sqrt((X - wx)**2 + (Y - wy)**2)
            r = np.maximum(r, 5.0)  # avoid singularity

            # Point source solution: P(r,t) = Q/(4πDt) * exp(-r²/4Dt)
            t = time_days * 86400
            dP_well = (dP_injection / (4 * np.pi * D_hydraulic * t)) * \
                       np.exp(-r**2 / (4 * D_hydraulic * t))

            # Normalize to injection pressure at well
            dP_well = dP_well / (dP_well.max() + 1e-10) * (dP_injection / 1e6)
            dP += dP_well

        # Stress changes from pressure (poroelastic coupling)
        alpha_biot = 0.7
        nu = 0.25
        factor = alpha_biot * (1 - 2*nu) / (2 * (1 - nu))

        d_sigma_n = -factor * dP     # normal stress change [MPa]
        d_tau = factor * 0.3 * dP    # shear stress change [MPa] (simplified)

        return X, Y, dP, d_sigma_n, d_tau

    def compute_cff(self, d_sigma_n, d_tau, dP):
        """Compute Coulomb Failure Function change."""
        # Effective normal stress change
        d_sigma_eff = d_sigma_n - dP  # Terzaghi effective stress

        cff = d_tau + self.mu * d_sigma_eff
        return cff

    def estimate_seismicity_rate(self, cff, r_background=0.1):
        """
        Estimate seismicity rate using Dieterich (1994) rate-and-state friction.
        R = r_background * exp(CFF / (a*sigma_n))
        """
        a = 0.01   # rate-state parameter
        sigma_ref = 10.0  # MPa reference normal stress

        rate = r_background * np.exp(cff / (a * sigma_ref))
        return np.clip(rate, 0, 100)

    def magnitude_gutenberg_richter(self, b_value=1.0, M_min=0.5):
        """
        Gutenberg-Richter magnitude-frequency distribution.
        log10(N) = a - b*M
        """
        M = np.linspace(M_min, 5.0, 100)
        N_cum = 10 ** (3.0 - b_value * (M - M_min))

        return M, N_cum

    def compute_max_credible_magnitude(self, volume_km3):
        """
        McGarr (2014) maximum magnitude from injection volume.
        log10(M0_max) = log10(G * Delta_V)
        """
        G = 30e9  # Pa (shear modulus)
        Delta_V = volume_km3 * 1e9  # m^3
        M0_max = G * Delta_V  # N·m
        Mw_max = (2/3) * np.log10(M0_max) - 6.07  # moment magnitude
        return Mw_max

    def run_risk_analysis(self, injection_years=5,
                          Q_injection_m3_s=0.05):
        """Run comprehensive seismic risk analysis."""
        results = {}

        # Time series of CFF evolution
        times_days = [30, 90, 180, 365, 730, 1825]
        cff_evolution = []

        for t_days in times_days:
            X, Y, dP, d_sn, d_tau = self.compute_stress_change_from_injection(
                time_days=t_days)
            cff = self.compute_cff(d_sn, d_tau, dP)
            cff_evolution.append({
                'time_days': t_days,
                'CFF_max_MPa': float(cff.max()),
                'CFF_mean_MPa': float(cff.mean()),
                'area_positive_CFF_km2': float(np.sum(cff > 0.01) *
                                                (1000/80)**2 / 1e6),
                'seismicity_rate': float(self.estimate_seismicity_rate(cff).mean())
            })

        results['cff_evolution'] = cff_evolution

        # Magnitude-frequency
        M, N_cum = self.magnitude_gutenberg_richter()
        results['mag_freq'] = {'M': M.tolist(), 'N_cum': N_cum.tolist()}

        # Maximum credible magnitude
        V_total = Q_injection_m3_s * injection_years * 365 * 86400 / 1e9
        Mw_max = self.compute_max_credible_magnitude(V_total)
        results['Mw_max_credible'] = float(Mw_max)
        results['total_injection_volume_km3'] = float(V_total)

        # Traffic light protocol
        results['traffic_light'] = self._traffic_light_protocol(cff_evolution[-1])

        return results

    def _traffic_light_protocol(self, last_cff):
        """TRAFFIC LIGHT PROTOCOL for induced seismicity management."""
        Mw_max = 3.0  # simplified
        if last_cff['CFF_max_MPa'] < 0.1:
            status = 'GREEN'
            action = 'Continue normal operations'
        elif last_cff['CFF_max_MPa'] < 0.5:
            status = 'YELLOW'
            action = 'Reduce injection rate by 50%, increase monitoring'
        elif last_cff['CFF_max_MPa'] < 1.0:
            status = 'ORANGE'
            action = 'Suspend injection, assess fault stability'
        else:
            status = 'RED'
            action = 'Immediate shutdown, fault reactivation risk'

        return {'status': status, 'action': action,
                'CFF_max': last_cff['CFF_max_MPa']}

    def plot_results(self, save_path=None):
        """Comprehensive seismic risk visualization."""
        fig = plt.figure(figsize=(20, 14))
        fig.suptitle('Coulomb Stress Change & Induced Seismicity Risk\n'
                     'Kakkonda EGS - Injection Scenario', fontsize=14)

        # CFF maps at different times
        time_snapshots = [30, 365, 1825]
        axes_cff = [fig.add_subplot(3, 3, i+1) for i in range(3)]

        for ax, t_days in zip(axes_cff, time_snapshots):
            X, Y, dP, d_sn, d_tau = self.compute_stress_change_from_injection(
                time_days=t_days)
            cff = self.compute_cff(d_sn, d_tau, dP)

            vmax = max(abs(cff.min()), abs(cff.max()), 0.01)
            im = ax.contourf(X/1000, Y/1000, cff,
                             levels=np.linspace(-vmax, vmax, 21),
                             cmap='RdBu_r', extend='both')
            plt.colorbar(im, ax=ax, label='ΔCFF [MPa]', shrink=0.8)

            # Well locations
            ax.plot(-0.15, 0, 'bv', markersize=10, label='Inj. Wells')
            ax.plot(0.15, 0, 'bv', markersize=10)
            ax.plot(0, 0, 'r^', markersize=12, label='Prod. Well')

            ax.set_xlabel('X [km]')
            ax.set_ylabel('Y [km]')
            ax.set_title(f'ΔCFF at t = {t_days} days')
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)

        # CFF evolution
        ax4 = fig.add_subplot(3, 3, 4)
        results = self.run_risk_analysis()
        cev = pd.DataFrame(results['cff_evolution'])
        ax4.semilogy(cev['time_days'], cev['CFF_max_MPa'], 'r-o', linewidth=2)
        ax4.axhline(y=0.1, color='green', linestyle='--', label='Green threshold')
        ax4.axhline(y=0.5, color='orange', linestyle='--', label='Yellow threshold')
        ax4.axhline(y=1.0, color='red', linestyle='--', label='Orange threshold')
        ax4.set_xlabel('Time [days]')
        ax4.set_ylabel('Max ΔCFF [MPa]')
        ax4.set_title('CFF Evolution Over Time')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        # Seismicity rate
        ax5 = fig.add_subplot(3, 3, 5)
        ax5.plot(cev['time_days'], cev['seismicity_rate'], 'b-s', linewidth=2)
        ax5.set_xlabel('Time [days]')
        ax5.set_ylabel('Seismicity Rate [events/day]')
        ax5.set_title('Induced Seismicity Rate')
        ax5.grid(True, alpha=0.3)

        # Gutenberg-Richter
        ax6 = fig.add_subplot(3, 3, 6)
        M = np.array(results['mag_freq']['M'])
        N = np.array(results['mag_freq']['N_cum'])
        ax6.semilogy(M, N, 'k-', linewidth=2, label='b=1.0 (Kakkonda)')
        ax6.axvline(x=results['Mw_max_credible'], color='red',
                    linestyle='--', linewidth=2,
                    label=f'Mw_max = {results["Mw_max_credible"]:.1f}')
        ax6.axvline(x=2.5, color='orange', linestyle=':',
                    label='Regulatory limit (M2.5)')
        ax6.set_xlabel('Magnitude (Mw)')
        ax6.set_ylabel('Cumulative Frequency')
        ax6.set_title('Gutenberg-Richter Distribution')
        ax6.legend(fontsize=8)
        ax6.grid(True, alpha=0.3)

        # Pore pressure diffusion
        ax7 = fig.add_subplot(3, 3, 7)
        r_arr = np.linspace(10, 1000, 200)
        for t_days in [30, 90, 180, 365]:
            t = t_days * 86400
            D = 0.1
            dP_arr = np.exp(-r_arr**2 / (4 * D * t))
            dP_arr = dP_arr / dP_arr.max() * 10  # normalize to 10 MPa
            ax7.plot(r_arr, dP_arr, label=f't = {t_days}d')
        ax7.set_xlabel('Distance from Well [m]')
        ax7.set_ylabel('Pore Pressure Change [MPa]')
        ax7.set_title('Pressure Diffusion Front')
        ax7.legend(fontsize=8)
        ax7.grid(True, alpha=0.3)

        # Traffic light status
        ax8 = fig.add_subplot(3, 3, 8)
        tl = results['traffic_light']
        colors_tl = {'GREEN': 'green', 'YELLOW': 'gold',
                     'ORANGE': 'orange', 'RED': 'red'}
        color = colors_tl.get(tl['status'], 'gray')
        circle = plt.Circle((0.5, 0.5), 0.35, color=color, zorder=2)
        ax8.add_patch(circle)
        ax8.set_xlim(0, 1)
        ax8.set_ylim(0, 1)
        ax8.text(0.5, 0.05, tl['action'], ha='center', va='bottom',
                 fontsize=9, wrap=True, transform=ax8.transAxes)
        ax8.text(0.5, 0.5, tl['status'], ha='center', va='center',
                 fontsize=16, fontweight='bold', color='white', zorder=3)
        ax8.set_title('Traffic Light Protocol Status')
        ax8.axis('off')

        # Positive CFF area
        ax9 = fig.add_subplot(3, 3, 9)
        ax9.fill_between(cev['time_days'],
                          cev['area_positive_CFF_km2'],
                          alpha=0.5, color='red')
        ax9.plot(cev['time_days'], cev['area_positive_CFF_km2'],
                 'r-o', linewidth=2)
        ax9.set_xlabel('Time [days]')
        ax9.set_ylabel('Area with ΔCFF > 0.01 MPa [km²]')
        ax9.set_title('Seismogenic Zone Evolution')
        ax9.grid(True, alpha=0.3)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()

        return results


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    model = CoulombStressModel()
    results = model.plot_results(save_path='figures/coulomb_stress.png')

    with open('results/seismic_risk_analysis.json', 'w') as f:
        json.dump(results, f, indent=2, default=float)

    print("Coulomb stress analysis completed.")
    tl = results['traffic_light']
    print(f"  Traffic Light Status: {tl['status']}")
    print(f"  Max CFF: {tl['CFF_max']:.3f} MPa")
    print(f"  Max Credible Magnitude: Mw {results['Mw_max_credible']:.1f}")
