"""
Module 3: Dose-Response Mathematical Modeling (Extended Hill Equation)
=====================================================================
Mathematical framework for biosensor dose-response characterization:
- Standard Hill equation fitting
- Extended Hill model with basal/maximal activity
- Cooperative binding models
- Biphasic response modeling
- Dynamic range and sensitivity analysis
- Noise modeling and limit of detection
"""

import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.stats import norm
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Callable
import json
import os


@dataclass
class HillParameters:
    """Parameters for extended Hill equation."""
    Vmin: float        # Basal expression (leakiness)
    Vmax: float        # Maximum expression
    Kd: float          # Half-maximal concentration
    n: float           # Hill coefficient (cooperativity)
    sigma_noise: float # Measurement noise
    
    @property
    def dynamic_range(self) -> float:
        return self.Vmax / max(self.Vmin, 1e-10)
    
    @property
    def fold_induction(self) -> float:
        return self.Vmax / max(self.Vmin, 1e-10)
    
    @property
    def sensitivity(self) -> float:
        """Maximum slope of the dose-response curve."""
        return self.n * (self.Vmax - self.Vmin) / (4 * self.Kd)


@dataclass
class BiosensorMetrics:
    """Performance metrics for a biosensor design."""
    LOD: float               # Limit of detection
    LOQ: float               # Limit of quantification
    linear_range: Tuple[float, float]
    dynamic_range_dB: float
    sensitivity: float
    selectivity_ratio: float
    response_time: float     # arbitrary units
    signal_to_noise: float


# ---- Hill equation variants ----

def hill_equation(x: np.ndarray, Vmin: float, Vmax: float,
                   Kd: float, n: float) -> np.ndarray:
    """Standard extended Hill equation: V = Vmin + (Vmax - Vmin) * x^n / (Kd^n + x^n)"""
    xn = np.power(np.maximum(x, 1e-20), n)
    kn = np.power(Kd, n)
    return Vmin + (Vmax - Vmin) * xn / (kn + xn)


def hill_repressor(x: np.ndarray, Vmin: float, Vmax: float,
                    Kd: float, n: float) -> np.ndarray:
    """Repressor Hill equation (inverse)."""
    xn = np.power(np.maximum(x, 1e-20), n)
    kn = np.power(Kd, n)
    return Vmax - (Vmax - Vmin) * xn / (kn + xn)


def biphasic_hill(x: np.ndarray, Vmin: float, Vmax: float,
                   Kd1: float, n1: float,
                   Kd2: float, n2: float, w: float) -> np.ndarray:
    """Biphasic Hill model (dual binding site)."""
    term1 = np.power(np.maximum(x, 1e-20), n1) / (np.power(Kd1, n1) + np.power(np.maximum(x, 1e-20), n1))
    term2 = np.power(np.maximum(x, 1e-20), n2) / (np.power(Kd2, n2) + np.power(np.maximum(x, 1e-20), n2))
    return Vmin + (Vmax - Vmin) * (w * term1 + (1 - w) * term2)


def thermodynamic_model(x: np.ndarray, K_LR: float, K_DNA_apo: float,
                          K_DNA_holo: float, n_coop: float,
                          alpha: float, P_total: float) -> np.ndarray:
    """
    Thermodynamic model of allosteric TF-based biosensor.
    Models the partition function of all TF states:
    - Apo TF (unbound)
    - Holo TF (ligand-bound)
    - TF-DNA complex (apo and holo)
    """
    # Fractional occupancy of ligand binding
    f_bound = np.power(np.maximum(x, 1e-20), n_coop) / (np.power(K_LR, n_coop) + np.power(np.maximum(x, 1e-20), n_coop))
    
    # DNA binding affinity shifts with ligand binding
    K_DNA_eff = K_DNA_apo * (1 - f_bound) + K_DNA_holo * f_bound
    
    # Fractional occupancy of DNA
    f_DNA = P_total / (K_DNA_eff + P_total)
    
    # Reporter expression (proportional to free promoter for de-repression)
    reporter = alpha * (1 - f_DNA)
    
    return reporter


# ---- Fitting and analysis ----

def fit_hill_equation(concentrations: np.ndarray,
                       responses: np.ndarray,
                       model: str = "activator") -> HillParameters:
    """Fit Hill equation to dose-response data."""
    func = hill_equation if model == "activator" else hill_repressor
    
    # Initial guesses
    p0 = [np.min(responses), np.max(responses),
          np.median(concentrations), 1.5]
    
    bounds = ([0, 0, 1e-10, 0.1], [np.inf, np.inf, np.inf, 10.0])
    
    try:
        popt, pcov = curve_fit(func, concentrations, responses,
                                p0=p0, bounds=bounds, maxfev=10000)
        residuals = responses - func(concentrations, *popt)
        sigma = np.std(residuals)
    except:
        popt = p0
        sigma = np.std(responses) * 0.1
    
    return HillParameters(
        Vmin=popt[0], Vmax=popt[1], Kd=popt[2], n=popt[3],
        sigma_noise=sigma
    )


def calculate_LOD(params: HillParameters, confidence: float = 0.99) -> float:
    """Calculate limit of detection (3σ method)."""
    z = norm.ppf(confidence)
    # LOD: concentration where signal exceeds Vmin + 3*sigma
    signal_threshold = params.Vmin + 3 * params.sigma_noise
    
    if signal_threshold >= params.Vmax:
        return np.inf
    
    # Inverse Hill equation
    ratio = (signal_threshold - params.Vmin) / (params.Vmax - signal_threshold)
    if ratio <= 0:
        return np.inf
    
    LOD = params.Kd * np.power(ratio, 1.0 / params.n)
    return LOD


def calculate_linear_range(params: HillParameters,
                            linearity_threshold: float = 0.95) -> Tuple[float, float]:
    """Calculate the linear dynamic range of the biosensor."""
    # Linear range: where Hill curve approximates a line
    # For Hill equation, approximately EC10 to EC90
    EC10 = params.Kd * np.power(0.1 / 0.9, 1.0 / params.n)
    EC90 = params.Kd * np.power(0.9 / 0.1, 1.0 / params.n)
    return (EC10, EC90)


def generate_synthetic_dose_response(params: HillParameters,
                                       n_concentrations: int = 12,
                                       n_replicates: int = 3,
                                       conc_range: Tuple[float, float] = (1e-3, 1e3)) -> Tuple[np.ndarray, np.ndarray]:
    """Generate synthetic dose-response data with noise."""
    np.random.seed(42)
    
    concentrations = np.logspace(np.log10(conc_range[0]),
                                  np.log10(conc_range[1]),
                                  n_concentrations)
    
    all_concs = []
    all_responses = []
    
    for c in concentrations:
        for _ in range(n_replicates):
            response = hill_equation(np.array([c]),
                                      params.Vmin, params.Vmax,
                                      params.Kd, params.n)[0]
            response += np.random.normal(0, params.sigma_noise)
            response = max(0, response)
            all_concs.append(c)
            all_responses.append(response)
    
    return np.array(all_concs), np.array(all_responses)


def compute_biosensor_metrics(params: HillParameters,
                                target_conc_range: Tuple[float, float] = (1e-3, 1e3)) -> BiosensorMetrics:
    """Compute comprehensive biosensor performance metrics."""
    LOD = calculate_LOD(params)
    LOQ = LOD * 3.33  # LOQ ≈ 10σ
    linear_range = calculate_linear_range(params)
    
    dynamic_range_dB = 20 * np.log10(max(params.dynamic_range, 1.0))
    sensitivity = params.sensitivity
    
    # Signal-to-noise at Kd
    signal_at_kd = hill_equation(np.array([params.Kd]),
                                  params.Vmin, params.Vmax,
                                  params.Kd, params.n)[0]
    snr = (signal_at_kd - params.Vmin) / max(params.sigma_noise, 1e-10)
    
    return BiosensorMetrics(
        LOD=LOD,
        LOQ=LOQ,
        linear_range=linear_range,
        dynamic_range_dB=dynamic_range_dB,
        sensitivity=sensitivity,
        selectivity_ratio=1.0,
        response_time=1.0,
        signal_to_noise=snr
    )


def run_dose_response_modeling(output_dir: str = "results") -> Dict:
    """Run dose-response modeling for all biosensor configurations."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Define biosensor configurations
    configs = {
        "MerR_Hg": HillParameters(Vmin=50, Vmax=5000, Kd=0.1, n=1.8, sigma_noise=80),
        "ArsR_As": HillParameters(Vmin=100, Vmax=8000, Kd=1.0, n=1.5, sigma_noise=120),
        "CadC_Cd": HillParameters(Vmin=80, Vmax=6000, Kd=0.5, n=2.0, sigma_noise=100),
        "CueR_Cu": HillParameters(Vmin=60, Vmax=4500, Kd=0.05, n=2.2, sigma_noise=70),
        "SmtB_Zn": HillParameters(Vmin=120, Vmax=7000, Kd=2.0, n=1.3, sigma_noise=150),
        "NahR_Naphthalene": HillParameters(Vmin=200, Vmax=10000, Kd=5.0, n=1.1, sigma_noise=200),
        "DmpR_Phenol": HillParameters(Vmin=150, Vmax=9000, Kd=10.0, n=1.4, sigma_noise=180),
    }
    
    all_results = {}
    
    for name, params in configs.items():
        print(f"  Modeling {name}...")
        
        # Generate synthetic data
        concs, responses = generate_synthetic_dose_response(params)
        
        # Fit model
        fitted = fit_hill_equation(concs, responses)
        
        # Calculate metrics
        metrics = compute_biosensor_metrics(fitted)
        
        # Thermodynamic model comparison
        thermo_response = thermodynamic_model(
            np.logspace(-3, 3, 100),
            K_LR=params.Kd, K_DNA_apo=0.01, K_DNA_holo=10.0,
            n_coop=params.n, alpha=params.Vmax, P_total=1.0
        )
        
        all_results[name] = {
            "fitted_parameters": {
                "Vmin": round(fitted.Vmin, 2),
                "Vmax": round(fitted.Vmax, 2),
                "Kd_uM": round(fitted.Kd, 4),
                "Hill_coefficient": round(fitted.n, 3),
                "noise_sigma": round(fitted.sigma_noise, 2),
            },
            "metrics": {
                "LOD_uM": round(metrics.LOD, 6),
                "LOQ_uM": round(metrics.LOQ, 6),
                "linear_range_uM": [round(metrics.linear_range[0], 6),
                                     round(metrics.linear_range[1], 4)],
                "dynamic_range_dB": round(metrics.dynamic_range_dB, 2),
                "dynamic_range_fold": round(fitted.dynamic_range, 1),
                "sensitivity": round(metrics.sensitivity, 4),
                "SNR_at_Kd": round(metrics.signal_to_noise, 2),
            },
            "true_parameters": {
                "Vmin": params.Vmin,
                "Vmax": params.Vmax,
                "Kd_uM": params.Kd,
                "Hill_coefficient": params.n,
            },
            "n_data_points": len(concs),
        }
    
    with open(os.path.join(output_dir, "dose_response_modeling.json"), 'w') as f:
        json.dump(all_results, f, indent=2)
    
    return all_results


if __name__ == "__main__":
    results = run_dose_response_modeling()
    for name, data in results.items():
        print(f"\n=== {name} ===")
        m = data["metrics"]
        print(f"  LOD: {m['LOD_uM']:.4f} µM")
        print(f"  Dynamic range: {m['dynamic_range_fold']:.1f}-fold ({m['dynamic_range_dB']:.1f} dB)")
        print(f"  Hill coefficient: {data['fitted_parameters']['Hill_coefficient']:.2f}")
