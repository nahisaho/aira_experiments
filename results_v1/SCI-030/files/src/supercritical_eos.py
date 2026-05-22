"""
Supercritical Water Equation of State and Transport Properties
Based on IAPWS-IF97 formulation with near-critical enhancement
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import pandas as pd
import json
import os

try:
    from iapws import IAPWS97, IAPWS95
    IAPWS_AVAILABLE = True
except ImportError:
    IAPWS_AVAILABLE = False
    print("iapws not available, using simplified EOS")


class SupercriticalWaterEOS:
    """
    Equation of State for supercritical water (H2O).
    Covers the supercritical domain relevant to EGS (T: 200-600°C, P: 10-100 MPa).
    
    Critical point: Tc = 374.15°C, Pc = 22.064 MPa, rhoc = 322 kg/m^3
    """

    Tc = 374.15    # °C
    Pc = 22.064e6  # Pa
    rhoc = 322.0   # kg/m^3
    Rc = 8.314     # J/(mol·K)
    Mw = 0.018015  # kg/mol

    def __init__(self):
        self.T_critical = self.Tc
        self.P_critical = self.Pc

    def _reduced_vars(self, T_C, P_Pa):
        """Compute reduced temperature and pressure."""
        tau = (T_C + 273.15) / (self.Tc + 273.15)
        pi = P_Pa / self.Pc
        return tau, pi

    def density(self, T_C, P_Pa):
        """
        Water density [kg/m^3] using IAPWS-IF97 or simplified formulation.
        """
        if IAPWS_AVAILABLE:
            try:
                state = IAPWS97(T=T_C + 273.15, P=P_Pa / 1e6)
                if state.phase in ['Supercritical', 'Gas', 'Liquid', 'SupercriticalPressure']:
                    return state.rho
            except Exception:
                pass

        # Simplified EOS for supercritical region
        T_K = T_C + 273.15
        P_MPa = P_Pa / 1e6
        # Approximation valid for T > 350°C, P > 15 MPa
        rho = (P_MPa * 1000) / (0.461522 * T_K) * \
              (1 + 0.1 * np.exp(-0.01 * (T_K - 647)) * np.sqrt(P_MPa / 22.064))
        return max(10.0, min(1000.0, rho))

    def enthalpy(self, T_C, P_Pa):
        """Specific enthalpy [J/kg] using IAPWS-IF97."""
        if IAPWS_AVAILABLE:
            try:
                state = IAPWS97(T=T_C + 273.15, P=P_Pa / 1e6)
                return state.h * 1000  # kJ/kg -> J/kg
            except Exception:
                pass

        # Simplified: cp * T (approximate)
        cp = self.heat_capacity(T_C, P_Pa)
        return cp * (T_C - 0)

    def heat_capacity(self, T_C, P_Pa):
        """
        Isobaric heat capacity [J/(kg·K)].
        Shows strong enhancement near critical point.
        """
        if IAPWS_AVAILABLE:
            try:
                state = IAPWS97(T=T_C + 273.15, P=P_Pa / 1e6)
                return state.cp * 1000  # kJ/(kg·K) -> J/(kg·K)
            except Exception:
                pass

        # Simplified with critical enhancement
        T_K = T_C + 273.15
        P_MPa = P_Pa / 1e6

        # Near-critical enhancement
        dT = abs(T_K - (self.Tc + 273.15))
        dP = abs(P_MPa - self.Pc / 1e6)
        enhancement = 1 + 5 * np.exp(-0.005 * dT**2 - 0.1 * dP**2)

        cp_base = 2000 + 3000 * np.exp(-0.003 * (T_K - 500)**2)
        return cp_base * enhancement

    def viscosity(self, T_C, P_Pa):
        """
        Dynamic viscosity [Pa·s] using IAPWS formulation.
        Decreases significantly with temperature in supercritical state.
        """
        if IAPWS_AVAILABLE:
            try:
                state = IAPWS97(T=T_C + 273.15, P=P_Pa / 1e6)
                return state.mu
            except Exception:
                pass

        # IAPWS-2008 simplified correlation
        T_K = T_C + 273.15
        # Dilute gas contribution
        mu0 = 1e-6 * 1.002 * np.sqrt(T_K / 647.096) / \
              (1 + 0.396 / (T_K / 647.096))
        # Density correction
        rho = self.density(T_C, P_Pa)
        mu = mu0 * (1 + 0.00000033 * rho)
        return max(1e-6, min(1e-3, mu))

    def thermal_conductivity(self, T_C, P_Pa):
        """
        Thermal conductivity [W/(m·K)].
        Shows enhancement near critical point (lambda divergence).
        """
        if IAPWS_AVAILABLE:
            try:
                state = IAPWS97(T=T_C + 273.15, P=P_Pa / 1e6)
                return state.k
            except Exception:
                pass

        T_K = T_C + 273.15
        # Near-critical enhancement
        dT = abs(T_K - 647.096)
        enhancement = 1 + 3 * np.exp(-0.02 * dT)
        k_base = 0.6 * (T_K / 300) ** (-0.3)
        return k_base * enhancement

    def compute_grid(self, T_range=(200, 600), P_range=(10e6, 100e6),
                     n_T=50, n_P=50):
        """Compute EOS properties over T-P grid."""
        T_arr = np.linspace(T_range[0], T_range[1], n_T)
        P_arr = np.linspace(P_range[0], P_range[1], n_P)
        T_grid, P_grid = np.meshgrid(T_arr, P_arr)

        props = {
            'T': T_grid, 'P': P_grid / 1e6,
            'density': np.zeros_like(T_grid),
            'enthalpy': np.zeros_like(T_grid),
            'heat_capacity': np.zeros_like(T_grid),
            'viscosity': np.zeros_like(T_grid),
            'thermal_conductivity': np.zeros_like(T_grid)
        }

        for i in range(n_P):
            for j in range(n_T):
                T = T_arr[j]
                P = P_arr[i]
                props['density'][i, j] = self.density(T, P)
                props['enthalpy'][i, j] = self.enthalpy(T, P) / 1e6  # MJ/kg
                props['heat_capacity'][i, j] = self.heat_capacity(T, P) / 1000  # kJ/(kg·K)
                props['viscosity'][i, j] = self.viscosity(T, P) * 1e6  # µPa·s
                props['thermal_conductivity'][i, j] = self.thermal_conductivity(T, P)

        return props, T_arr, P_arr

    def plot_properties(self, save_path=None):
        """Plot EOS properties in T-P space."""
        print("  Computing EOS grid...")
        props, T_arr, P_arr = self.compute_grid(n_T=40, n_P=40)

        fig = plt.figure(figsize=(18, 12))
        fig.suptitle('Supercritical Water EOS - IAPWS-IF97\n'
                     'Kakkonda EGS Operating Range', fontsize=14)

        prop_list = [
            ('density', 'Density [kg/m³]', 'plasma'),
            ('enthalpy', 'Specific Enthalpy [MJ/kg]', 'inferno'),
            ('heat_capacity', 'Cp [kJ/(kg·K)]', 'hot'),
            ('viscosity', 'Viscosity [µPa·s]', 'viridis'),
            ('thermal_conductivity', 'Thermal Conductivity [W/(m·K)]', 'cividis'),
        ]

        for idx, (pname, label, cmap) in enumerate(prop_list):
            ax = fig.add_subplot(2, 3, idx + 1)
            T_grid = props['T']
            P_grid = props['P']

            im = ax.contourf(T_grid, P_grid, props[pname], levels=20, cmap=cmap)
            plt.colorbar(im, ax=ax, label=label, shrink=0.8)

            # Critical point
            ax.plot(self.Tc, self.Pc / 1e6, 'w*', markersize=12,
                    label='Critical point')
            # Kakkonda operating point
            ax.plot(380, 30, 'r^', markersize=10, label='Kakkonda (~3500 m)')
            # Critical isochor
            ax.axvline(x=self.Tc, color='white', linestyle='--', alpha=0.5, lw=0.8)
            ax.axhline(y=self.Pc/1e6, color='white', linestyle='--', alpha=0.5, lw=0.8)

            ax.set_xlabel('Temperature [°C]')
            ax.set_ylabel('Pressure [MPa]')
            ax.set_title(label)
            if idx == 0:
                ax.legend(fontsize=8)

        # Phase diagram summary
        ax6 = fig.add_subplot(2, 3, 6)
        phase_map = np.zeros_like(props['density'])
        for i in range(phase_map.shape[0]):
            for j in range(phase_map.shape[1]):
                T = T_arr[j]
                P_MPa = P_arr[i] / 1e6
                if T > self.Tc and P_MPa > self.Pc / 1e6:
                    phase_map[i, j] = 3   # Supercritical
                elif T > self.Tc:
                    phase_map[i, j] = 2   # Superheated vapor
                elif P_MPa > self.Pc / 1e6:
                    phase_map[i, j] = 1   # Compressed liquid
                else:
                    phase_map[i, j] = 0   # Two-phase / subcritical

        cmap_phase = plt.get_cmap('Set2', 4)
        im = ax6.contourf(props['T'], props['P'], phase_map,
                          levels=[-0.5, 0.5, 1.5, 2.5, 3.5],
                          cmap=cmap_phase)
        cbar = plt.colorbar(im, ax=ax6)
        cbar.set_ticks([0, 1, 2, 3])
        cbar.set_ticklabels(['Subcritical', 'Compressed Liq.', 'Superheated Vap.', 'Supercritical'])
        ax6.plot(self.Tc, self.Pc / 1e6, 'k*', markersize=14, label='Critical Point')
        ax6.plot(380, 30, 'rv', markersize=12, label='Kakkonda EGS')
        ax6.set_xlabel('Temperature [°C]')
        ax6.set_ylabel('Pressure [MPa]')
        ax6.set_title('Phase Diagram')
        ax6.legend(fontsize=8)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()

    def summary_at_kakkonda(self):
        """Compute and return fluid properties at Kakkonda conditions."""
        conditions = [
            ('Surface (1 km)', 80, 10e6),
            ('Transition (2.5 km)', 250, 25e6),
            ('Supercritical (3.5 km)', 380, 30e6),
            ('Deep SC (4 km)', 450, 40e6),
        ]
        results = []
        for name, T, P in conditions:
            results.append({
                'location': name,
                'T_C': T,
                'P_MPa': P / 1e6,
                'density_kg_m3': self.density(T, P),
                'enthalpy_kJ_kg': self.enthalpy(T, P) / 1000,
                'cp_kJ_kgK': self.heat_capacity(T, P) / 1000,
                'viscosity_uPas': self.viscosity(T, P) * 1e6,
                'k_W_mK': self.thermal_conductivity(T, P)
            })
        return pd.DataFrame(results)


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    eos = SupercriticalWaterEOS()
    eos.plot_properties(save_path='figures/eos_properties.png')

    summary = eos.summary_at_kakkonda()
    summary.to_csv('results/eos_kakkonda_summary.csv', index=False)

    print("EOS analysis completed.")
    print(summary.to_string(index=False))
