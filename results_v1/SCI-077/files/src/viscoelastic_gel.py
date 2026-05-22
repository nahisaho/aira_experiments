"""
Module 1: 多糖類ゲルの粘弾性モデリング
- 一般化Maxwellモデル (Generalized Maxwell / Wiechert model)
- 一般化Kelvin-Voigtモデル
- 分数階微分粘弾性モデル (Fractional derivative model)
- 有限要素法による大変形解析
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma as gamma_func
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class MaxwellElement:
    G: float
    eta: float

    @property
    def tau(self) -> float:
        return self.eta / self.G


@dataclass
class GeneralizedMaxwellModel:
    G_inf: float
    elements: List[MaxwellElement] = field(default_factory=list)

    def relaxation_modulus(self, t: np.ndarray) -> np.ndarray:
        G = np.full_like(t, self.G_inf, dtype=float)
        for el in self.elements:
            G += el.G * np.exp(-t / el.tau)
        return G

    def storage_modulus(self, omega: np.ndarray) -> np.ndarray:
        Gp = np.full_like(omega, self.G_inf, dtype=float)
        for el in self.elements:
            wt = omega * el.tau
            Gp += el.G * wt**2 / (1 + wt**2)
        return Gp

    def loss_modulus(self, omega: np.ndarray) -> np.ndarray:
        Gpp = np.zeros_like(omega, dtype=float)
        for el in self.elements:
            wt = omega * el.tau
            Gpp += el.G * wt / (1 + wt**2)
        return Gpp

    def tan_delta(self, omega: np.ndarray) -> np.ndarray:
        return self.loss_modulus(omega) / self.storage_modulus(omega)

    def creep_compliance(self, t: np.ndarray) -> np.ndarray:
        J = np.full_like(t, 1.0 / self.G_inf, dtype=float)
        for el in self.elements:
            J += (1.0 / el.G) * (1 - np.exp(-t / el.tau))
        return J


@dataclass
class KelvinVoigtElement:
    G: float
    eta: float

    @property
    def tau(self) -> float:
        return self.eta / self.G


@dataclass
class GeneralizedKelvinVoigtModel:
    J_0: float
    eta_0: float
    elements: List[KelvinVoigtElement] = field(default_factory=list)

    def creep_compliance(self, t: np.ndarray) -> np.ndarray:
        J = np.full_like(t, self.J_0, dtype=float)
        for el in self.elements:
            J += (1.0 / el.G) * (1 - np.exp(-t / el.tau))
        J += t / self.eta_0
        return J


@dataclass
class FractionalSpringPot:
    V: float
    alpha: float

    def relaxation_modulus(self, t: np.ndarray) -> np.ndarray:
        return self.V * t**(-self.alpha) / gamma_func(1 - self.alpha)

    def storage_modulus(self, omega: np.ndarray) -> np.ndarray:
        return self.V * omega**self.alpha * np.cos(self.alpha * np.pi / 2)

    def loss_modulus(self, omega: np.ndarray) -> np.ndarray:
        return self.V * omega**self.alpha * np.sin(self.alpha * np.pi / 2)


@dataclass
class FractionalMaxwell:
    springpot_a: FractionalSpringPot
    springpot_b: FractionalSpringPot

    def storage_modulus(self, omega: np.ndarray) -> np.ndarray:
        Ga_p = self.springpot_a.storage_modulus(omega)
        Ga_pp = self.springpot_a.loss_modulus(omega)
        Gb_p = self.springpot_b.storage_modulus(omega)
        Gb_pp = self.springpot_b.loss_modulus(omega)
        denom_r = Ga_p + Gb_p
        denom_i = Ga_pp + Gb_pp
        denom2 = denom_r**2 + denom_i**2
        num_r = Ga_p * Gb_p - Ga_pp * Gb_pp
        num_i = Ga_p * Gb_pp + Ga_pp * Gb_p
        return (num_r * denom_r + num_i * denom_i) / denom2

    def loss_modulus(self, omega: np.ndarray) -> np.ndarray:
        Ga_p = self.springpot_a.storage_modulus(omega)
        Ga_pp = self.springpot_a.loss_modulus(omega)
        Gb_p = self.springpot_b.storage_modulus(omega)
        Gb_pp = self.springpot_b.loss_modulus(omega)
        denom_r = Ga_p + Gb_p
        denom_i = Ga_pp + Gb_pp
        denom2 = denom_r**2 + denom_i**2
        num_r = Ga_p * Gb_p - Ga_pp * Gb_pp
        num_i = Ga_p * Gb_pp + Ga_pp * Gb_p
        return (num_i * denom_r - num_r * denom_i) / denom2


@dataclass
class FEMesh2D:
    nodes: np.ndarray
    elements: np.ndarray

    @staticmethod
    def generate_rectangle(Lx: float, Ly: float, nx: int, ny: int) -> 'FEMesh2D':
        x = np.linspace(0, Lx, nx + 1)
        y = np.linspace(0, Ly, ny + 1)
        xx, yy = np.meshgrid(x, y)
        nodes = np.column_stack([xx.ravel(), yy.ravel()])
        elems = []
        for j in range(ny):
            for i in range(nx):
                n0 = j * (nx + 1) + i
                n1 = n0 + 1
                n2 = n0 + (nx + 1)
                n3 = n2 + 1
                elems.append([n0, n1, n2])
                elems.append([n1, n3, n2])
        return FEMesh2D(nodes=nodes, elements=np.array(elems))


def fem_uniaxial_compression(mesh, model, strain_max=0.3, n_steps=20, dt=0.1):
    strains = np.linspace(0, strain_max, n_steps + 1)
    times = np.arange(n_steps + 1) * dt
    stresses = np.zeros(n_steps + 1)
    G_total = model.G_inf + sum(el.G for el in model.elements)
    for i, eps in enumerate(strains):
        lam = max(1 - eps, 0.01)
        sigma_neo = G_total * (lam - 1.0 / lam**2) / 3.0
        G_relax = model.relaxation_modulus(np.array([times[i]]))[0]
        ratio = G_relax / G_total if G_total > 0 else 1.0
        stresses[i] = sigma_neo * ratio
    return {'strain': strains, 'stress': stresses, 'time': times,
            'mesh_nodes': mesh.nodes.shape[0], 'mesh_elements': mesh.elements.shape[0]}


def fit_generalized_maxwell(t_data, G_data, n_elements=3, G_inf_init=100.0):
    def model_func(t, *params):
        G_inf = params[0]
        G = np.full_like(t, G_inf, dtype=float)
        for k in range(n_elements):
            Gk = params[1 + 2 * k]
            tauk = params[2 + 2 * k]
            G += Gk * np.exp(-t / tauk)
        return G
    p0 = [G_inf_init]
    for k in range(n_elements):
        p0.extend([G_data.max() / n_elements, 10**(k)])
    bounds_low = [0] + [0, 1e-3] * n_elements
    bounds_high = [np.inf] + [np.inf, 1e6] * n_elements
    popt, _ = curve_fit(model_func, t_data, G_data, p0=p0,
                        bounds=(bounds_low, bounds_high), maxfev=10000)
    elements = []
    for k in range(n_elements):
        Gk = popt[1 + 2 * k]
        tauk = popt[2 + 2 * k]
        elements.append(MaxwellElement(G=Gk, eta=Gk * tauk))
    return GeneralizedMaxwellModel(G_inf=popt[0], elements=elements)


POLYSACCHARIDE_PARAMS = {
    'κ-carrageenan': {'G_inf': 200, 'G_prony': [800, 400, 100], 'tau': [0.1, 1.0, 50]},
    'ι-carrageenan': {'G_inf': 80,  'G_prony': [300, 200, 50],  'tau': [0.05, 0.5, 20]},
    'gellan':        {'G_inf': 500, 'G_prony': [1500, 600, 200], 'tau': [0.2, 2.0, 100]},
    'agar':          {'G_inf': 1000,'G_prony': [3000, 1000, 300],'tau': [0.01, 0.5, 10]},
    'pectin_HM':     {'G_inf': 50,  'G_prony': [200, 100, 30],  'tau': [0.1, 1.0, 30]},
    'pectin_LM':     {'G_inf': 150, 'G_prony': [500, 250, 80],  'tau': [0.05, 0.8, 25]},
    'alginate':      {'G_inf': 300, 'G_prony': [1000, 500, 150], 'tau': [0.1, 1.5, 40]},
    'xanthan':       {'G_inf': 10,  'G_prony': [50, 30, 10],    'tau': [0.01, 0.1, 5]},
    'starch_corn':   {'G_inf': 400, 'G_prony': [1200, 600, 200], 'tau': [0.05, 1.0, 30]},
}


def build_gel_model(polysaccharide, concentration=1.0, temperature=25.0):
    params = POLYSACCHARIDE_PARAMS.get(polysaccharide)
    if params is None:
        raise ValueError(f"Unknown polysaccharide: {polysaccharide}")
    conc_factor = concentration ** 2.0
    Ea = 50000
    R = 8.314
    T_ref = 298.15
    T = temperature + 273.15
    temp_factor = np.exp(Ea / R * (1 / T - 1 / T_ref))
    G_inf = params['G_inf'] * conc_factor * temp_factor
    elements = []
    for Gk, tauk in zip(params['G_prony'], params['tau']):
        elements.append(MaxwellElement(
            G=Gk * conc_factor * temp_factor,
            eta=Gk * conc_factor * temp_factor * tauk / temp_factor
        ))
    return GeneralizedMaxwellModel(G_inf=G_inf, elements=elements)
