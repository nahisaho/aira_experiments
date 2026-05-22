"""
Module 5: 3Dフードプリンティングの印刷性予測
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Tuple
from scipy.optimize import minimize


@dataclass
class FoodInkRheology:
    yield_stress: float
    K: float
    n: float
    G_storage: float
    G_loss: float
    thixotropy_index: float = 1.0
    recovery_time: float = 10.0

    @property
    def tan_delta(self):
        return self.G_loss / self.G_storage if self.G_storage > 0 else float('inf')

    def apparent_viscosity(self, gamma_dot):
        return self.yield_stress / (gamma_dot + 1e-10) + self.K * gamma_dot**(self.n-1)

    def shear_stress(self, gamma_dot):
        return self.yield_stress + self.K * gamma_dot**self.n

    def structural_recovery(self, t):
        return self.G_storage * (1 - np.exp(-t / self.recovery_time))


@dataclass
class ExtrusionModel:
    nozzle_diameter: float = 1.0
    nozzle_length: float = 20.0
    reservoir_diameter: float = 25.0

    def extrusion_pressure(self, ink, flow_rate):
        R = self.nozzle_diameter / 2 * 1e-3
        L = self.nozzle_length * 1e-3
        Q = flow_rate * 1e-9
        gamma_dot_wall = (3*ink.n+1)/(4*ink.n) * (32*Q)/(np.pi*(2*R)**3)
        tau_wall = ink.yield_stress + ink.K * gamma_dot_wall**ink.n
        dP = 2 * tau_wall * L / R
        dP_entry = 0.5*ink.K*(Q/(np.pi*R**2))**ink.n*((self.reservoir_diameter/self.nozzle_diameter)**2-1)
        return dP + dP_entry

    def extrusion_force(self, ink, flow_rate):
        P = self.extrusion_pressure(ink, flow_rate)
        A = np.pi * (self.reservoir_diameter/2*1e-3)**2
        return P * A


@dataclass
class PrintabilityPredictor:
    layer_height: float = 1.0
    line_width: float = 1.5
    print_speed: float = 20.0

    def die_swell_ratio(self, ink):
        swell = 0.1 + (1 + (2*ink.K)**2)**(1/6) * ink.n**0.5
        return min(swell, 2.0)

    def shape_retention_index(self, ink, n_layers=10):
        rho, g = 1200, 9.81
        h = n_layers * self.layer_height * 1e-3
        sigma_g = rho * g * h
        if ink.yield_stress > sigma_g:
            delta_h = sigma_g / (2*ink.G_storage) * h
            sri = 1 - delta_h / h
        else:
            sri = 1 / (sigma_g / ink.yield_stress)
        return np.clip(sri, 0, 1)

    def layer_adhesion_index(self, ink, time_between_layers=5.0):
        G_rec = ink.structural_recovery(np.array([time_between_layers]))[0]
        return np.clip(1 - G_rec/ink.G_storage, 0, 1)

    def line_uniformity(self, ink, flow_rate):
        swell = self.die_swell_ratio(ink)
        return np.clip(1 - abs(swell-1)/swell, 0, 1)

    def overall_printability_score(self, ink, flow_rate=10.0, n_layers=10):
        sri = self.shape_retention_index(ink, n_layers)
        adhesion = self.layer_adhesion_index(ink)
        uniformity = self.line_uniformity(ink, flow_rate)
        yield_s = np.exp(-((np.log10(max(ink.yield_stress,1))-2.5)/0.5)**2)
        tan_d_s = np.exp(-((ink.tan_delta-0.2)/0.15)**2)
        n_s = np.exp(-((ink.n-0.4)/0.15)**2)
        thixo_s = np.exp(-((ink.recovery_time-5)/5)**2)
        scores = {'shape_retention': sri, 'yield_stress_fit': yield_s,
                  'tan_delta_fit': tan_d_s, 'n_index_fit': n_s,
                  'layer_adhesion': adhesion, 'line_uniformity': uniformity,
                  'thixotropy_fit': thixo_s}
        weights = {'shape_retention':0.25,'yield_stress_fit':0.20,'tan_delta_fit':0.15,
                   'n_index_fit':0.10,'layer_adhesion':0.10,'line_uniformity':0.10,'thixotropy_fit':0.10}
        scores['overall'] = sum(scores[k]*weights[k] for k in weights)
        return scores


def optimize_printing_parameters(ink, target_layers=10):
    def objective(params):
        nd, sp, fr, wt = params
        pred = PrintabilityPredictor(layer_height=abs(nd)*0.8, line_width=abs(nd)*1.2, print_speed=abs(sp))
        return -pred.overall_printability_score(ink, flow_rate=abs(fr), n_layers=target_layers)['overall']
    result = minimize(objective, [1.0,20.0,10.0,5.0], method='Nelder-Mead', options={'maxiter':500})
    opt = np.abs(result.x)
    pred = PrintabilityPredictor(layer_height=opt[0]*0.8, line_width=opt[0]*1.2, print_speed=opt[1])
    scores = pred.overall_printability_score(ink, flow_rate=opt[2], n_layers=target_layers)
    ext = ExtrusionModel(nozzle_diameter=opt[0])
    return {'optimal_nozzle_diameter': opt[0], 'optimal_print_speed': opt[1],
            'optimal_flow_rate': opt[2], 'optimal_wait_time': opt[3],
            'extrusion_force': ext.extrusion_force(ink, opt[2]), 'printability_scores': scores}


def rheology_suitability_map(n_points=25):
    tau_y = np.logspace(1, 3.3, n_points)
    G_prime = np.logspace(2, 4, n_points)
    TY, GP = np.meshgrid(tau_y, G_prime)
    scores = np.zeros_like(TY)
    pred = PrintabilityPredictor()
    for i in range(n_points):
        for j in range(n_points):
            ink = FoodInkRheology(TY[i,j], TY[i,j]/100, 0.4, GP[i,j], GP[i,j]*0.2, recovery_time=5.0)
            scores[i,j] = pred.overall_printability_score(ink, n_layers=10)['overall']
    return {'yield_stress': TY, 'G_storage': GP, 'printability_score': scores}
