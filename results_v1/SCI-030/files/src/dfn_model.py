"""
Discrete Fracture Network (DFN) Model for EGS Reservoir Simulation
Kakkonda / Tohoku region geological conditions
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import pandas as pd
from scipy.stats import lognorm, vonmises
import json
import os

np.random.seed(42)


class FractureSet:
    """Single fracture set with statistical distributions."""

    def __init__(self, name, density, length_mean, length_std,
                 orientation_mean, orientation_kappa, aperture_mean, aperture_std):
        self.name = name
        self.density = density           # fractures per m^2 (P21)
        self.length_mean = length_mean   # m (log-normal)
        self.length_std = length_std
        self.orientation_mean = orientation_mean  # degrees (von Mises)
        self.orientation_kappa = orientation_kappa
        self.aperture_mean = aperture_mean       # mm
        self.aperture_std = aperture_std


class DFNModel:
    """
    Stochastic Discrete Fracture Network model.
    Kakkonda geological context:
     - Depth: 2500-4000 m (supercritical zone ~3729 m)
     - Dominant strike: NNW-SSE (tectonic extension direction)
     - Secondary: ENE-WSW
    """

    def __init__(self, domain_size=(500, 500), depth=3500):
        self.domain_size = domain_size  # (x, y) meters
        self.depth = depth              # meters below surface
        self.fractures = []
        self.fracture_sets = self._define_kakkonda_fractures()

    def _define_kakkonda_fractures(self):
        """
        Fracture sets based on Kakkonda/Tohoku geological data.
        Reference: Kakkonda WD-1a borehole data (Doi et al., 1998)
        """
        sets = [
            FractureSet(
                name="NNW-SSE_tensile",
                density=0.008,        # P21 density
                length_mean=80,       # m
                length_std=0.6,       # log-normal sigma
                orientation_mean=160,  # degrees (NNW-SSE)
                orientation_kappa=8.0, # concentration
                aperture_mean=0.5,    # mm
                aperture_std=0.2
            ),
            FractureSet(
                name="ENE-WSW_shear",
                density=0.005,
                length_mean=50,
                length_std=0.8,
                orientation_mean=75,   # ENE-WSW
                orientation_kappa=5.0,
                aperture_mean=0.3,
                aperture_std=0.15
            ),
            FractureSet(
                name="Horizontal_bedding",
                density=0.003,
                length_mean=120,
                length_std=0.5,
                orientation_mean=0,   # horizontal
                orientation_kappa=12.0,
                aperture_mean=0.2,
                aperture_std=0.1
            )
        ]
        return sets

    def generate(self):
        """Generate stochastic DFN realization."""
        self.fractures = []
        Lx, Ly = self.domain_size
        area = Lx * Ly

        for fset in self.fracture_sets:
            n_fractures = int(fset.density * area)
            for _ in range(n_fractures):
                # Center position
                cx = np.random.uniform(0, Lx)
                cy = np.random.uniform(0, Ly)

                # Length (log-normal)
                length = lognorm.rvs(fset.length_std,
                                     scale=fset.length_mean)
                length = np.clip(length, 5, 500)

                # Orientation (von Mises)
                angle_rad = vonmises.rvs(fset.orientation_kappa,
                                          loc=np.radians(fset.orientation_mean))
                angle_deg = np.degrees(angle_rad) % 180

                # Aperture (log-normal)
                aperture = lognorm.rvs(fset.aperture_std / fset.aperture_mean,
                                        scale=fset.aperture_mean)
                aperture = np.clip(aperture, 0.01, 5.0)  # mm

                # Hydraulic aperture (cubic law)
                hydraulic_aperture = aperture * 1e-3  # convert to m
                transmissivity = (hydraulic_aperture ** 3) / (12 * 1e-4)  # m^2/s

                # Endpoints
                dx = (length / 2) * np.cos(np.radians(angle_deg))
                dy = (length / 2) * np.sin(np.radians(angle_deg))
                x1, y1 = cx - dx, cy - dy
                x2, y2 = cx + dx, cy + dy

                self.fractures.append({
                    'set': fset.name,
                    'cx': cx, 'cy': cy,
                    'x1': x1, 'y1': y1,
                    'x2': x2, 'y2': y2,
                    'length': length,
                    'angle': angle_deg,
                    'aperture_mm': aperture,
                    'transmissivity': transmissivity,
                    'depth': self.depth
                })

        self.fractures = pd.DataFrame(self.fractures)
        return self

    def compute_connectivity(self):
        """Compute fracture intersection network (simplified)."""
        n = len(self.fractures)
        intersections = 0

        # Check pairwise intersections (sample for large networks)
        sample_idx = np.random.choice(n, min(200, n), replace=False)
        for i in range(len(sample_idx)):
            for j in range(i + 1, len(sample_idx)):
                fi = self.fractures.iloc[sample_idx[i]]
                fj = self.fractures.iloc[sample_idx[j]]
                if self._segments_intersect(
                    (fi.x1, fi.y1), (fi.x2, fi.y2),
                    (fj.x1, fj.y1), (fj.x2, fj.y2)
                ):
                    intersections += 1

        connectivity = intersections / (len(sample_idx) * (len(sample_idx) - 1) / 2)
        p32 = len(self.fractures) / (self.domain_size[0] *
                                      self.domain_size[1] * 10)  # per m^3
        return {
            'n_fractures': n,
            'intersections_sampled': intersections,
            'connectivity_index': connectivity,
            'P21_density': len(self.fractures) / (self.domain_size[0] * self.domain_size[1]),
            'P32_estimate': p32
        }

    @staticmethod
    def _segments_intersect(p1, p2, p3, p4):
        """Check if line segment p1-p2 intersects p3-p4."""
        def cross(o, a, b):
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

        d1 = cross(p3, p4, p1)
        d2 = cross(p3, p4, p2)
        d3 = cross(p1, p2, p3)
        d4 = cross(p1, p2, p4)

        if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and \
           ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
            return True
        return False

    def get_permeability_tensor(self):
        """
        Compute bulk permeability tensor from DFN.
        Snow (1969) approach for fractured rock.
        """
        k_tensor = np.zeros((2, 2))
        for _, f in self.fractures.iterrows():
            theta = np.radians(f['angle'])
            e = f['aperture_mm'] * 1e-3  # m
            spacing = 1.0 / (f['length'] * self.fracture_sets[0].density)
            spacing = max(spacing, 0.5)

            # Permeability from cubic law
            k_frac = e ** 3 / (12 * spacing)

            n = np.array([np.cos(theta + np.pi/2), np.sin(theta + np.pi/2)])
            I = np.eye(2)
            nn = np.outer(n, n)
            k_tensor += k_frac * (I - nn)

        k_tensor /= len(self.fractures)
        return k_tensor

    def plot(self, save_path=None, wells=None):
        """Plot DFN with well locations."""
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        colors = {'NNW-SSE_tensile': '#E74C3C',
                  'ENE-WSW_shear': '#2ECC71',
                  'Horizontal_bedding': '#3498DB'}

        ax = axes[0]
        for _, f in self.fractures.iterrows():
            color = colors.get(f['set'], 'gray')
            alpha = min(0.7, f['aperture_mm'] / 2.0)
            lw = f['aperture_mm'] * 0.8
            ax.plot([f.x1, f.x2], [f.y1, f.y2],
                    color=color, alpha=alpha, linewidth=lw)

        if wells:
            for wname, wpos in wells.items():
                marker = '^' if 'inj' in wname.lower() else 'v'
                ax.plot(wpos[0], wpos[1], marker=marker, markersize=12,
                        color='black', zorder=10,
                        label=wname.replace('_', ' '))

        patches = [mpatches.Patch(color=c, label=k.replace('_', ' '))
                   for k, c in colors.items()]
        ax.legend(handles=patches, loc='upper right', fontsize=8)
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_title(f'DFN Realization - Kakkonda Reservoir ({self.depth} m depth)')
        ax.set_xlim(0, self.domain_size[0])
        ax.set_ylim(0, self.domain_size[1])
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        # Rose diagram (fracture orientations)
        ax2 = axes[1]
        angles = self.fractures['angle'].values
        bins = np.arange(0, 181, 10)
        hist, _ = np.histogram(angles, bins=bins)
        theta = np.radians(bins[:-1])
        width = np.radians(10)

        ax2 = plt.subplot(122, projection='polar')
        bars = ax2.bar(theta, hist, width=width, alpha=0.7,
                       color=['#E74C3C' if 140 <= b <= 180 or 0 <= b <= 20 else
                              '#2ECC71' if 60 <= b <= 90 else '#3498DB'
                              for b in bins[:-1]])
        ax2.set_title('Fracture Orientation Rose Diagram', pad=20)
        ax2.set_theta_zero_location('N')
        ax2.set_theta_direction(-1)
        ax2.set_xlabel('Frequency')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        return fig


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    dfn = DFNModel(domain_size=(500, 500), depth=3500)
    dfn.generate()

    wells = {
        'Injection_Well_1': (150, 250),
        'Injection_Well_2': (350, 250),
        'Production_Well': (250, 250)
    }
    dfn.plot(save_path='figures/dfn_model.png', wells=wells)

    stats = dfn.compute_connectivity()
    k_tensor = dfn.get_permeability_tensor()

    stats['permeability_kxx_m2'] = float(k_tensor[0, 0])
    stats['permeability_kyy_m2'] = float(k_tensor[1, 1])
    stats['permeability_kxy_m2'] = float(k_tensor[0, 1])

    with open('results/dfn_statistics.json', 'w') as f:
        json.dump(stats, f, indent=2)

    dfn.fractures.to_csv('results/dfn_fractures.csv', index=False)
    print("DFN Model completed.")
    print(f"  Total fractures: {stats['n_fractures']}")
    print(f"  P21 density: {stats['P21_density']:.6f} m^-1")
    print(f"  Kxx = {k_tensor[0,0]:.3e} m^2")
    print(f"  Kyy = {k_tensor[1,1]:.3e} m^2")
