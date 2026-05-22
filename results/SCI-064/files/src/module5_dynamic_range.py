"""
Module 5: Reporter Output Dynamic Range Maximization
=====================================================
Optimization of genetic circuit design to maximize biosensor dynamic range:
- Promoter-RBS strength optimization
- Genetic circuit architecture comparison
- Noise analysis and signal processing
- Multi-layer amplification cascade design
- Feedback loop optimization
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import json
import os


@dataclass
class CircuitParameters:
    """Parameters for a genetic circuit."""
    promoter_strength: float     # Max transcription rate (au/min)
    rbs_strength: float          # Translation efficiency
    protein_degradation: float   # Degradation rate (1/min)
    mrna_degradation: float      # mRNA degradation rate (1/min)
    tf_concentration: float      # TF total concentration (nM)
    cooperativity: float         # Hill coefficient
    Kd: float                    # TF-DNA binding Kd (nM)
    basal_rate: float            # Leaky transcription


@dataclass
class CircuitPerformance:
    """Performance metrics of a circuit design."""
    fold_induction: float
    dynamic_range_dB: float
    response_time_min: float
    noise_coefficient_variation: float
    output_ON: float
    output_OFF: float
    linear_range: Tuple[float, float]


def simulate_simple_circuit(params: CircuitParameters,
                              ligand_conc: float,
                              t_max: float = 300,
                              dt: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate simple TF-promoter-reporter circuit.
    dM/dt = alpha * f(L) + basal - delta_m * M
    dP/dt = beta * M - delta_p * P
    """
    def f_ligand(L, Kd, n):
        """Fractional activation by ligand."""
        return L**n / (Kd**n + L**n)
    
    def odes(y, t, p, L):
        M, P = y
        f_L = f_ligand(L, p.Kd, p.cooperativity)
        
        dM = p.promoter_strength * f_L + p.basal_rate - p.mrna_degradation * M
        dP = p.rbs_strength * M - p.protein_degradation * P
        return [dM, dP]
    
    t = np.arange(0, t_max, dt)
    y0 = [0.0, 0.0]
    sol = odeint(odes, y0, t, args=(params, ligand_conc))
    
    return t, sol


def simulate_cascade_circuit(params: CircuitParameters,
                               ligand_conc: float,
                               n_layers: int = 2,
                               t_max: float = 600) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate multi-layer amplification cascade.
    Layer 1: TF → intermediate protein
    Layer 2: intermediate → reporter
    """
    def odes(y, t, p, L, n_layers):
        n_vars = 2 * n_layers  # mRNA + protein per layer
        dydt = np.zeros(n_vars)
        
        # First layer: regulated by ligand
        f_L = L**p.cooperativity / (p.Kd**p.cooperativity + L**p.cooperativity)
        dydt[0] = p.promoter_strength * f_L + p.basal_rate - p.mrna_degradation * y[0]
        dydt[1] = p.rbs_strength * y[0] - p.protein_degradation * y[1]
        
        # Subsequent layers: regulated by previous protein
        for layer in range(1, n_layers):
            prev_protein = y[2 * (layer - 1) + 1]
            Kd_cascade = 50.0  # nM
            n_cascade = 2.0
            f_cascade = prev_protein**n_cascade / (Kd_cascade**n_cascade + prev_protein**n_cascade)
            
            idx = 2 * layer
            dydt[idx] = p.promoter_strength * 0.8 * f_cascade + p.basal_rate * 0.5 - p.mrna_degradation * y[idx]
            dydt[idx + 1] = p.rbs_strength * y[idx] - p.protein_degradation * 0.8 * y[idx + 1]
        
        return dydt
    
    t = np.arange(0, t_max, 0.1)
    y0 = np.zeros(2 * n_layers)
    sol = odeint(odes, y0, t, args=(params, ligand_conc, n_layers))
    
    return t, sol


def simulate_feedback_circuit(params: CircuitParameters,
                                ligand_conc: float,
                                feedback_type: str = "positive",
                                feedback_strength: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate circuit with feedback loop.
    Positive feedback: reporter enhances its own transcription
    Negative feedback: reporter represses its transcription
    """
    def odes(y, t, p, L, fb_type, fb_str):
        M, P = y
        f_L = L**p.cooperativity / (p.Kd**p.cooperativity + L**p.cooperativity)
        
        if fb_type == "positive":
            fb_factor = 1.0 + fb_str * P / (100.0 + P)
        else:
            fb_factor = 1.0 / (1.0 + fb_str * P / 100.0)
        
        dM = p.promoter_strength * f_L * fb_factor + p.basal_rate - p.mrna_degradation * M
        dP = p.rbs_strength * M - p.protein_degradation * P
        return [dM, dP]
    
    t = np.arange(0, 500, 0.1)
    y0 = [0.0, 0.0]
    sol = odeint(odes, y0, t, args=(params, ligand_conc, feedback_type, feedback_strength))
    
    return t, sol


def evaluate_circuit(params: CircuitParameters,
                      circuit_type: str = "simple",
                      n_layers: int = 1) -> CircuitPerformance:
    """Evaluate circuit performance across ligand concentrations."""
    # Simulate OFF state
    if circuit_type == "simple":
        _, sol_off = simulate_simple_circuit(params, 0.0)
        _, sol_on = simulate_simple_circuit(params, params.Kd * 100)
    elif circuit_type == "cascade":
        _, sol_off = simulate_cascade_circuit(params, 0.0, n_layers)
        _, sol_on = simulate_cascade_circuit(params, params.Kd * 100, n_layers)
    elif circuit_type == "positive_feedback":
        _, sol_off = simulate_feedback_circuit(params, 0.0, "positive")
        _, sol_on = simulate_feedback_circuit(params, params.Kd * 100, "positive")
    elif circuit_type == "negative_feedback":
        _, sol_off = simulate_feedback_circuit(params, 0.0, "negative")
        _, sol_on = simulate_feedback_circuit(params, params.Kd * 100, "negative")
    else:
        raise ValueError(f"Unknown circuit type: {circuit_type}")
    
    # Get steady-state outputs
    if circuit_type == "cascade":
        output_off = sol_off[-1, -1]
        output_on = sol_on[-1, -1]
        # Response time (time to reach 90% of steady state)
        target = output_on * 0.9
        protein_trace = sol_on[:, -1]
    else:
        output_off = sol_off[-1, 1]
        output_on = sol_on[-1, 1]
        target = output_on * 0.9
        protein_trace = sol_on[:, 1]
    
    output_off = max(output_off, 1e-6)
    fold_induction = output_on / output_off
    
    # Response time
    above_target = np.where(protein_trace >= target)[0]
    response_time = above_target[0] * 0.1 if len(above_target) > 0 else 999
    
    # Noise (CV) - simplified stochastic estimate
    cv = 1.0 / np.sqrt(max(output_on, 1.0)) + 0.05
    
    # Dynamic range
    dr_dB = 20 * np.log10(max(fold_induction, 1.0))
    
    return CircuitPerformance(
        fold_induction=round(fold_induction, 1),
        dynamic_range_dB=round(dr_dB, 2),
        response_time_min=round(response_time, 1),
        noise_coefficient_variation=round(cv, 4),
        output_ON=round(output_on, 2),
        output_OFF=round(output_off, 4),
        linear_range=(params.Kd * 0.1, params.Kd * 10)
    )


def optimize_dynamic_range(base_params: CircuitParameters,
                             circuit_type: str = "simple") -> Dict:
    """
    Optimize circuit parameters to maximize dynamic range.
    Uses differential evolution for global optimization.
    """
    def objective(x):
        p = CircuitParameters(
            promoter_strength=x[0],
            rbs_strength=x[1],
            protein_degradation=x[2],
            mrna_degradation=x[3],
            tf_concentration=base_params.tf_concentration,
            cooperativity=base_params.cooperativity,
            Kd=base_params.Kd,
            basal_rate=x[4]
        )
        try:
            perf = evaluate_circuit(p, circuit_type)
            # Maximize fold induction while keeping response time < 120 min
            penalty = max(0, perf.response_time_min - 120) * 10
            return -(np.log10(max(perf.fold_induction, 1.0)) - penalty)
        except:
            return 0
    
    bounds = [
        (1.0, 100.0),     # promoter_strength
        (0.1, 10.0),      # rbs_strength
        (0.001, 0.1),     # protein_degradation
        (0.01, 0.5),      # mrna_degradation
        (0.01, 5.0),      # basal_rate
    ]
    
    result = differential_evolution(objective, bounds, seed=42,
                                      maxiter=50, popsize=15, tol=0.01)
    
    opt_params = CircuitParameters(
        promoter_strength=result.x[0],
        rbs_strength=result.x[1],
        protein_degradation=result.x[2],
        mrna_degradation=result.x[3],
        tf_concentration=base_params.tf_concentration,
        cooperativity=base_params.cooperativity,
        Kd=base_params.Kd,
        basal_rate=result.x[4]
    )
    
    opt_perf = evaluate_circuit(opt_params, circuit_type)
    
    return {
        "optimized_parameters": {
            "promoter_strength": round(opt_params.promoter_strength, 3),
            "rbs_strength": round(opt_params.rbs_strength, 3),
            "protein_degradation": round(opt_params.protein_degradation, 5),
            "mrna_degradation": round(opt_params.mrna_degradation, 5),
            "basal_rate": round(opt_params.basal_rate, 4),
        },
        "performance": {
            "fold_induction": opt_perf.fold_induction,
            "dynamic_range_dB": opt_perf.dynamic_range_dB,
            "response_time_min": opt_perf.response_time_min,
            "noise_CV": opt_perf.noise_coefficient_variation,
            "output_ON": opt_perf.output_ON,
            "output_OFF": opt_perf.output_OFF,
        }
    }


def run_dynamic_range_optimization(output_dir: str = "results") -> Dict:
    """Run complete dynamic range optimization pipeline."""
    os.makedirs(output_dir, exist_ok=True)
    
    biosensors = {
        "MerR_Hg": CircuitParameters(50, 5.0, 0.01, 0.1, 100, 1.8, 0.1, 0.5),
        "ArsR_As": CircuitParameters(40, 4.0, 0.01, 0.1, 150, 1.5, 1.0, 1.0),
        "CadC_Cd": CircuitParameters(60, 6.0, 0.01, 0.1, 120, 2.0, 0.5, 0.8),
        "NahR_Naph": CircuitParameters(30, 3.0, 0.01, 0.1, 200, 1.1, 5.0, 2.0),
        "DmpR_Phenol": CircuitParameters(35, 3.5, 0.01, 0.1, 180, 1.4, 10.0, 1.5),
    }
    
    circuit_types = ["simple", "cascade", "positive_feedback", "negative_feedback"]
    
    all_results = {}
    
    for sensor_name, base_params in biosensors.items():
        print(f"  Optimizing {sensor_name}...")
        sensor_results = {}
        
        # Evaluate all circuit architectures
        for ct in circuit_types:
            n_layers = 2 if ct == "cascade" else 1
            
            # Baseline performance
            baseline = evaluate_circuit(base_params, ct, n_layers)
            
            # Optimized performance
            optimized = optimize_dynamic_range(base_params, ct)
            
            sensor_results[ct] = {
                "baseline": {
                    "fold_induction": baseline.fold_induction,
                    "dynamic_range_dB": baseline.dynamic_range_dB,
                    "response_time_min": baseline.response_time_min,
                    "noise_CV": baseline.noise_coefficient_variation,
                },
                "optimized": optimized,
            }
        
        # Find best architecture
        best_ct = max(circuit_types,
                       key=lambda ct: sensor_results[ct]["optimized"]["performance"]["fold_induction"])
        
        sensor_results["best_architecture"] = best_ct
        sensor_results["best_fold_induction"] = sensor_results[best_ct]["optimized"]["performance"]["fold_induction"]
        
        all_results[sensor_name] = sensor_results
    
    with open(os.path.join(output_dir, "dynamic_range_optimization.json"), 'w') as f:
        json.dump(all_results, f, indent=2)
    
    return all_results


if __name__ == "__main__":
    results = run_dynamic_range_optimization()
    for sensor, data in results.items():
        print(f"\n=== {sensor} ===")
        print(f"  Best architecture: {data['best_architecture']}")
        print(f"  Best fold induction: {data['best_fold_induction']:.1f}")
