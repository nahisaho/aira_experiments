"""
Pharmaceutical Intermediate Case Study: Continuous Suzuki-Miyaura Coupling
===========================================================================
Demonstrates end-to-end continuous flow synthesis optimization
for a Pd-catalyzed cross-coupling reaction.
"""

import numpy as np
import json, os

np.random.seed(42)

CASE_STUDY = {
    "reaction": "Suzuki-Miyaura Cross-Coupling",
    "substrate_A": "4-Bromoanisole",
    "substrate_B": "Phenylboronic acid",
    "catalyst": "Pd(PPh3)4",
    "base": "K2CO3",
    "solvent": "THF/H2O (4:1)",
    "product": "4-Methoxybiphenyl",
    "target_yield_pct": 95,
    "target_purity_pct": 99.5,
    "target_throughput_g_h": 50,
}

def suzuki_kinetics(temp_C, conc_A, conc_B, cat_pct, flow_rate, base_equiv=2.0):
    """Suzuki coupling kinetics with catalyst deactivation."""
    T_K = temp_C + 273.15
    Ea = 55000  # J/mol
    R = 8.314

    k0 = 5e7
    k = k0 * np.exp(-Ea / (R * T_K))

    # Catalyst effect (Pd loading)
    cat_factor = cat_pct / (cat_pct + 0.5)  # Michaelis-Menten-like

    # Base effect
    base_factor = min(1.0, base_equiv / 1.5)

    # Effective rate constant
    k_eff = k * cat_factor * base_factor

    # Residence time
    V_reactor = 10e-3  # 10 mL reactor volume
    tau = V_reactor / (flow_rate * 1e-3 / 60)  # seconds

    # Second-order reaction: dX/dt = k_eff * C_A0 * (1-X) * (theta_B - X)
    theta_B = conc_B / conc_A  # stoichiometric ratio
    Da = k_eff * conc_A * tau

    # Numerical integration (simple Euler)
    X = 0.0
    dt_sim = tau / 1000
    for _ in range(1000):
        if theta_B - X > 0 and 1 - X > 0:
            dX = k_eff * conc_A * (1 - X) * (theta_B - X) * dt_sim
            X += dX
            X = min(X, min(1.0, theta_B))

    # Side reactions (homo-coupling, protodeboronation)
    T_excess = max(0, temp_C - 90)
    homo_coupling = 0.005 * T_excess / 10 * X
    protodeboronation = 0.003 * (1 + T_excess / 20)

    selectivity = max(0.9, 1.0 - homo_coupling - protodeboronation)
    yield_pct = X * selectivity * 100

    # Product purity
    impurity_homo = homo_coupling * X * 100
    impurity_proto = protodeboronation * 100
    impurity_Pd = 0.1 * cat_pct  # Pd leaching
    purity = max(90, 100 - impurity_homo - impurity_proto - impurity_Pd)

    mw_product = 184.24  # g/mol
    throughput = yield_pct / 100 * conc_A * flow_rate * 1e-3 / 60 * mw_product * 3600

    return {
        "conversion_pct": round(X * 100, 2),
        "selectivity": round(selectivity, 4),
        "yield_pct": round(yield_pct, 2),
        "purity_pct": round(purity, 2),
        "throughput_g_h": round(throughput, 3),
        "residence_time_s": round(tau, 1),
        "damkohler": round(Da, 2),
        "impurities": {
            "homo_coupling_pct": round(impurity_homo, 3),
            "protodeboronation_pct": round(impurity_proto, 3),
            "Pd_residual_ppm": round(impurity_Pd * 10, 1),
        },
    }

def batch_vs_continuous_comparison():
    """Compare batch and continuous flow for the Suzuki coupling."""
    batch = {
        "mode": "Batch",
        "reaction_time_h": 4,
        "yield_pct": 88,
        "purity_pct": 97.5,
        "throughput_g_day": 120,
        "Pd_loading_mol_pct": 5.0,
        "solvent_volume_mL": 500,
        "temperature_C": 80,
        "scalability": "Linear (vessel size)",
        "reproducibility_rsd_pct": 5.2,
        "PMI": 32,  # Process Mass Intensity
        "E_factor": 28,
    }

    # Optimized continuous conditions
    cont_result = suzuki_kinetics(
        temp_C=100, conc_A=0.5, conc_B=0.6, cat_pct=2.0, flow_rate=2.0
    )
    continuous = {
        "mode": "Continuous Flow",
        "reaction_time_s": cont_result["residence_time_s"],
        "yield_pct": cont_result["yield_pct"],
        "purity_pct": cont_result["purity_pct"],
        "throughput_g_day": round(cont_result["throughput_g_h"] * 24, 1),
        "Pd_loading_mol_pct": 2.0,
        "solvent_volume_mL": 10,  # reactor volume
        "temperature_C": 100,
        "scalability": "Numbering up + mild scale-up",
        "reproducibility_rsd_pct": 1.8,
        "PMI": 8,
        "E_factor": 5,
    }

    advantages = [
        "60% reduction in catalyst loading",
        f"Residence time: {cont_result['residence_time_s']:.0f}s vs 4h (batch)",
        f"Improved reproducibility: RSD {continuous['reproducibility_rsd_pct']}% vs {batch['reproducibility_rsd_pct']}%",
        f"PMI reduction: {batch['PMI']} → {continuous['PMI']} ({(1-continuous['PMI']/batch['PMI'])*100:.0f}%)",
        f"E-factor improvement: {batch['E_factor']} → {continuous['E_factor']}",
        "Superheated conditions accessible safely",
        "Reduced Pd contamination in product",
        "Real-time quality control via PAT",
    ]

    return {"batch": batch, "continuous": continuous, "advantages": advantages}

def process_control_integration():
    """Design process control software architecture."""
    return {
        "software_architecture": {
            "control_layer": {
                "platform": "Siemens SIMATIC / Allen-Bradley CompactLogix",
                "protocol": "OPC-UA",
                "scan_rate_ms": 100,
                "safety": "SIL-2 rated E-stop and pressure relief",
            },
            "scada_layer": {
                "platform": "Ignition SCADA / LabVIEW",
                "features": [
                    "Real-time process visualization",
                    "Alarm management (ISA-18.2)",
                    "Historical data logging",
                    "Recipe management (ISA-88)",
                ],
            },
            "analytics_layer": {
                "platform": "Python + Flask API / MATLAB",
                "features": [
                    "Bayesian optimization engine",
                    "Multivariate statistical process control (MSPC)",
                    "Digital twin synchronization",
                    "Predictive maintenance",
                ],
            },
            "pat_integration": {
                "hplc": {
                    "model": "Agilent 1260 Infinity II",
                    "interface": "ICF (Instrument Control Framework)",
                    "method_time_min": 3,
                    "auto_sampling": True,
                },
                "ftir": {
                    "model": "Mettler Toledo ReactIR 702L",
                    "interface": "iC IR software API",
                    "sampling_interval_s": 15,
                    "probe": "DiComp (diamond ATR)",
                },
                "flow_sensors": {
                    "type": "Coriolis mass flow",
                    "model": "Bronkhorst mini CORI-FLOW",
                    "accuracy_pct": 0.2,
                },
            },
        },
        "communication_flow": [
            "Sensors → PLC (OPC-UA, 100ms cycle)",
            "PLC → SCADA (OPC-UA, 500ms cycle)",
            "SCADA → Analytics (REST API, 1-5s cycle)",
            "Analytics → PLC (setpoint updates, event-driven)",
            "PAT → Analytics (method-dependent: 15s IR, 3min HPLC)",
        ],
        "data_pipeline": {
            "raw_data": "InfluxDB (time-series, 100ms resolution)",
            "processed_data": "PostgreSQL (experiment records, batch genealogy)",
            "ml_models": "MLflow (model versioning, A/B testing)",
            "reports": "Automated PDF/HTML via Jinja2 templates",
        },
    }

def run_case_study():
    # Parameter sweep
    temps = [60, 70, 80, 90, 100, 110, 120]
    flow_rates = [0.5, 1.0, 1.5, 2.0, 3.0]
    cat_loadings = [1.0, 2.0, 3.0, 5.0]

    sweep_results = []
    for T in temps:
        for Q in flow_rates:
            for cat in cat_loadings:
                r = suzuki_kinetics(T, 0.5, 0.6, cat, Q)
                r.update({"temperature_C": T, "flow_rate_mL_min": Q, "catalyst_mol_pct": cat})
                sweep_results.append(r)

    # Find optimal
    best = max(sweep_results, key=lambda x: x["yield_pct"] if x["purity_pct"] >= 99.0 else 0)

    comparison = batch_vs_continuous_comparison()
    control = process_control_integration()

    results = {
        "case_study_info": CASE_STUDY,
        "parameter_sweep": {
            "n_conditions": len(sweep_results),
            "top_5_conditions": sorted(
                [r for r in sweep_results if r["purity_pct"] >= 99.0],
                key=lambda x: -x["yield_pct"]
            )[:5],
        },
        "optimal_conditions": {
            "temperature_C": best["temperature_C"],
            "flow_rate_mL_min": best["flow_rate_mL_min"],
            "catalyst_mol_pct": best["catalyst_mol_pct"],
            "concentration_M": 0.5,
            "stoichiometry": 1.2,
            "yield_pct": best["yield_pct"],
            "purity_pct": best["purity_pct"],
            "throughput_g_h": best["throughput_g_h"],
            "residence_time_s": best["residence_time_s"],
        },
        "batch_vs_continuous": comparison,
        "process_control": control,
    }

    return results

if __name__ == "__main__":
    results = run_case_study()
    os.makedirs("results", exist_ok=True)

    with open("results/case_study_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("=== Pharmaceutical Case Study: Suzuki-Miyaura Coupling ===")
    opt = results["optimal_conditions"]
    print(f"\nOptimal Conditions:")
    print(f"  Temperature: {opt['temperature_C']}°C")
    print(f"  Flow Rate: {opt['flow_rate_mL_min']} mL/min")
    print(f"  Catalyst: {opt['catalyst_mol_pct']} mol%")
    print(f"  Yield: {opt['yield_pct']}%")
    print(f"  Purity: {opt['purity_pct']}%")
    print(f"  Throughput: {opt['throughput_g_h']} g/h")

    comp = results["batch_vs_continuous"]
    print(f"\n--- Batch vs Continuous ---")
    print(f"  Batch yield: {comp['batch']['yield_pct']}%")
    print(f"  Continuous yield: {comp['continuous']['yield_pct']}%")
    print(f"\nAdvantages of continuous flow:")
    for a in comp["advantages"]:
        print(f"  • {a}")
