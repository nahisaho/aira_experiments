"""
Module 2: 乳化系の微視的構造と巨視的レオロジーの関係
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class EmulsionPhase:
    G_storage: float
    G_loss: float
    eta: float

    @property
    def G_complex(self) -> complex:
        return complex(self.G_storage, self.G_loss)


@dataclass
class EmulsionSystem:
    continuous: EmulsionPhase
    dispersed: EmulsionPhase
    volume_fraction: float
    droplet_radii: np.ndarray
    droplet_fractions: np.ndarray
    interfacial_tension: float


def palierne_complex_modulus(system, omega):
    G_storage = np.zeros_like(omega)
    G_loss = np.zeros_like(omega)
    for w_idx, w in enumerate(omega):
        Gm = complex(system.continuous.eta * w * 0.1, system.continuous.eta * w)
        Gd = complex(system.dispersed.G_storage, system.dispersed.G_loss + system.dispersed.eta * w)
        Gamma = system.interfacial_tension
        H_total = complex(0, 0)
        for R, f in zip(system.droplet_radii, system.droplet_fractions):
            alpha = Gamma / R
            numer = 4 * alpha * (2 * Gm + 5 * Gd) + (Gd - Gm) * (16 * Gm + 19 * Gd)
            denom = 40 * alpha * (Gm + Gd) + (2 * Gd + 3 * Gm) * (16 * Gm + 19 * Gd)
            if abs(denom) > 1e-30:
                H_total += f * numer / denom
        phi = system.volume_fraction
        Gstar = Gm * (1 + 3 * phi * H_total) / (1 - 2 * phi * H_total)
        G_storage[w_idx] = Gstar.real
        G_loss[w_idx] = Gstar.imag
    return G_storage, G_loss


def krieger_dougherty(phi, eta_0=1e-3, phi_max=0.64, intrinsic_viscosity=2.5):
    eta_rel = (1 - phi / phi_max) ** (-intrinsic_viscosity * phi_max)
    return eta_0 * eta_rel


def cross_model(gamma_dot, eta_0, eta_inf, K, n):
    return eta_inf + (eta_0 - eta_inf) / (1 + (K * gamma_dot) ** n)


@dataclass
class CGBead:
    position: np.ndarray
    velocity: np.ndarray
    force: np.ndarray
    bead_type: str
    mass: float = 1.0


@dataclass
class CGMDSimulation:
    beads: List[CGBead] = field(default_factory=list)
    box_size: np.ndarray = field(default_factory=lambda: np.array([20.0, 20.0, 20.0]))
    dt: float = 0.005
    temperature: float = 1.0
    gamma_dpd: float = 4.5
    a_params: dict = field(default_factory=dict)

    def initialize_emulsion(self, n_oil=500, n_water=2000, n_surfactant=200, droplet_radius=4.0):
        rng = np.random.default_rng(42)
        self.beads = []
        center = self.box_size / 2
        for _ in range(n_oil):
            while True:
                pos = center + rng.normal(0, droplet_radius / 2, 3)
                if np.linalg.norm(pos - center) < droplet_radius:
                    break
            self.beads.append(CGBead(pos, rng.normal(0, 0.1, 3), np.zeros(3), 'oil'))
        for _ in range(n_water):
            while True:
                pos = rng.uniform(0, self.box_size)
                if np.linalg.norm(pos - center) > droplet_radius * 1.2:
                    break
            self.beads.append(CGBead(pos, rng.normal(0, 0.1, 3), np.zeros(3), 'water'))
        for _ in range(n_surfactant):
            theta = rng.uniform(0, 2 * np.pi)
            phi_a = rng.uniform(0, np.pi)
            r = droplet_radius + rng.normal(0, 0.3)
            pos = center + r * np.array([np.sin(phi_a)*np.cos(theta),
                                         np.sin(phi_a)*np.sin(theta), np.cos(phi_a)])
            self.beads.append(CGBead(pos % self.box_size, rng.normal(0, 0.1, 3),
                                     np.zeros(3), 'surfactant_head'))
        self.a_params = {
            ('oil', 'oil'): 25.0, ('water', 'water'): 25.0,
            ('oil', 'water'): 100.0,
            ('surfactant_head', 'water'): 25.0, ('surfactant_head', 'oil'): 80.0,
        }

    def run(self, n_steps=500, save_interval=100):
        rng = np.random.default_rng(42)
        trajectory = {'step': [], 'droplet_radius': [], 'pressure': []}
        n = len(self.beads)
        for step in range(n_steps):
            # simplified velocity-verlet with DPD forces on subset
            for b in self.beads[:200]:
                b.force = np.zeros(3)
            for i in range(min(n, 200)):
                for j in range(i+1, min(n, 200)):
                    rij = self.beads[j].position - self.beads[i].position
                    rij -= self.box_size * np.round(rij / self.box_size)
                    dist = np.linalg.norm(rij)
                    if 1e-6 < dist < 1.0:
                        eij = rij / dist
                        w = 1 - dist
                        pair = tuple(sorted([self.beads[i].bead_type, self.beads[j].bead_type]))
                        a_ij = self.a_params.get(pair, 25.0)
                        F = a_ij * w * eij
                        self.beads[i].force -= F
                        self.beads[j].force += F
            for b in self.beads[:200]:
                b.velocity += self.dt * b.force / b.mass
                b.position = (b.position + self.dt * b.velocity) % self.box_size

            if step % save_interval == 0:
                oil_pos = np.array([b.position for b in self.beads if b.bead_type == 'oil'])
                if len(oil_pos) > 0:
                    c = oil_pos.mean(axis=0)
                    mean_r = np.linalg.norm(oil_pos - c, axis=1).mean()
                else:
                    mean_r = 0
                KE = sum(0.5*b.mass*np.dot(b.velocity, b.velocity) for b in self.beads[:200])
                trajectory['step'].append(step)
                trajectory['droplet_radius'].append(mean_r)
                trajectory['pressure'].append(KE / (3 * np.prod(self.box_size)))
        return trajectory


def droplet_size_distribution(d_mean, d_std, n_bins=20):
    sigma = np.sqrt(np.log(1 + (d_std / d_mean)**2))
    mu = np.log(d_mean) - sigma**2 / 2
    d = np.exp(np.linspace(mu - 3*sigma, mu + 3*sigma, n_bins))
    pdf = (1/(d*sigma*np.sqrt(2*np.pi))) * np.exp(-(np.log(d)-mu)**2/(2*sigma**2))
    fractions = pdf / pdf.sum()
    return d/2, fractions


def microstructure_to_rheology(phi, d_mean, d_std, eta_continuous=1e-3,
                                eta_dispersed=0.05, gamma_if=0.01, G_dispersed=100.0):
    radii, fractions = droplet_size_distribution(d_mean, d_std)
    system = EmulsionSystem(
        continuous=EmulsionPhase(0.0, 0.0, eta_continuous),
        dispersed=EmulsionPhase(G_dispersed, G_dispersed*0.1, eta_dispersed),
        volume_fraction=phi, droplet_radii=radii,
        droplet_fractions=fractions, interfacial_tension=gamma_if
    )
    omega = np.logspace(-2, 3, 50)
    G_prime, G_double = palierne_complex_modulus(system, omega)
    eta_apparent = krieger_dougherty(np.array([phi]), eta_0=eta_continuous)[0]
    d32 = np.sum(fractions*(2*radii)**3) / np.sum(fractions*(2*radii)**2)
    return {'omega': omega, 'G_storage': G_prime, 'G_loss': G_double,
            'eta_apparent': eta_apparent, 'd32': d32, 'phi': phi}
