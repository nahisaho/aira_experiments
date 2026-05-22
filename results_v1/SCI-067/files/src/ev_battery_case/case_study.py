"""
EV Battery Manufacturing LCA Case Study.

Complete LCA for a 75 kWh NMC811 lithium-ion battery pack,
demonstrating the full automation pipeline.
"""
from __future__ import annotations

import json
import math
import os
import sys

# Add parent to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from nlp_extractor.process_tree_builder import NLPProcessTreeBuilder, ProcessTree
from ecoinvent_matcher.matcher import EcoinventMatcher
from uncertainty.propagation import (
    MonteCarloSimulator, TaylorExpansion, UncertainParameter
)
from hotspot.analysis import (
    HotspotAnalyzer, ScenarioComparator, ScenarioDefinition
)
from scope3.estimator import Scope3Estimator


# ===========================================================================
# EV Battery Bill of Materials (NMC811 - 75 kWh pack)
# ===========================================================================
# Based on published LCA literature (Dai et al. 2019, Kelly et al. 2020,
# GREET model 2023, Ecoinvent 3.10)

BATTERY_BOM = [
    # Cathode materials (NMC811: 80% Ni, 10% Mn, 10% Co)
    {"component": "NMC811 cathode active material", "material": "nickel sulfate",
     "mass_kg": 72.0, "process": "Cathode material production", "supplier_location": "CN"},
    {"component": "NMC811 cathode active material", "material": "cobalt sulfate",
     "mass_kg": 9.0, "process": "Cathode material production", "supplier_location": "CD"},
    {"component": "NMC811 cathode active material", "material": "lithium carbonate",
     "mass_kg": 18.5, "process": "Cathode material production", "supplier_location": "CL"},
    {"component": "NMC811 cathode active material", "material": "manganese sulfate",
     "mass_kg": 8.5, "process": "Cathode material production", "supplier_location": "ZA"},
    {"component": "Cathode foil", "material": "aluminium",
     "mass_kg": 14.0, "process": "Aluminium foil production", "supplier_location": "CN"},

    # Anode materials
    {"component": "Anode active material", "material": "natural graphite",
     "mass_kg": 52.0, "process": "Graphite anode production", "supplier_location": "CN"},
    {"component": "Anode foil", "material": "copper",
     "mass_kg": 18.0, "process": "Copper foil production", "supplier_location": "CN"},

    # Electrolyte
    {"component": "Electrolyte", "material": "electrolyte (LiPF6)",
     "mass_kg": 28.0, "process": "Electrolyte production", "supplier_location": "JP"},

    # Separator
    {"component": "Separator", "material": "separator (PE/PP)",
     "mass_kg": 8.0, "process": "Separator production", "supplier_location": "KR"},

    # Cell housing and module components
    {"component": "Cell housing", "material": "aluminium",
     "mass_kg": 35.0, "process": "Cell housing manufacturing", "supplier_location": "CN"},
    {"component": "Module housing", "material": "aluminium",
     "mass_kg": 25.0, "process": "Module assembly", "supplier_location": "DE"},

    # Pack components
    {"component": "Battery pack enclosure", "material": "steel",
     "mass_kg": 45.0, "process": "Pack enclosure manufacturing", "supplier_location": "DE"},
    {"component": "BMS electronics", "material": "copper",
     "mass_kg": 3.5, "process": "BMS production", "supplier_location": "CN"},
    {"component": "Cooling system", "material": "aluminium",
     "mass_kg": 12.0, "process": "Cooling system manufacturing", "supplier_location": "DE"},
    {"component": "Wiring harness", "material": "copper",
     "mass_kg": 5.0, "process": "Wiring production", "supplier_location": "CN"},
    {"component": "Thermal interface material", "material": "plastics (average)",
     "mass_kg": 4.0, "process": "TIM application", "supplier_location": "US"},
]

# Energy consumption for manufacturing processes
ENERGY_CONSUMPTION = {
    "cell_manufacturing": {
        "electricity_kwh": 42.0,   # per kWh battery capacity
        "heat_mj": 25.0,          # per kWh battery capacity
        "location": "CN",          # China grid
    },
    "module_assembly": {
        "electricity_kwh": 3.5,
        "heat_mj": 1.2,
        "location": "DE",
    },
    "pack_assembly": {
        "electricity_kwh": 5.0,
        "heat_mj": 2.0,
        "location": "DE",
    },
    "dry_room_operation": {
        "electricity_kwh": 15.0,   # significant energy for dehumidification
        "heat_mj": 0,
        "location": "CN",
    },
}

# Transport distances (tkm = tonne-kilometres)
TRANSPORT_DATA = [
    {"route": "Raw materials → Cell plant", "mode": "ocean freight",
     "distance_km": 8000, "mass_t": 0.25},
    {"route": "Cell plant → Module plant", "mode": "truck transport",
     "distance_km": 500, "mass_t": 0.35},
    {"route": "Module plant → OEM", "mode": "truck transport",
     "distance_km": 200, "mass_t": 0.42},
]

# End-of-life assumptions
EOL_ASSUMPTIONS = {
    "collection_rate": 0.95,
    "recycling_rate": 0.70,
    "hydrometallurgical_fraction": 0.60,
    "pyrometallurgical_fraction": 0.40,
    "material_recovery": {
        "cobalt": 0.95,
        "nickel": 0.92,
        "lithium": 0.80,
        "copper": 0.98,
        "aluminium": 0.90,
    },
}


def build_process_tree() -> ProcessTree:
    """Build the EV battery process tree from BOM."""
    builder = NLPProcessTreeBuilder()
    return builder.build_from_bom(BATTERY_BOM, "NMC811 75kWh Battery Pack")


def compute_process_impacts() -> dict[str, dict]:
    """
    Compute per-process GHG impacts using emission factors.

    Returns dict of {process_id: {"name", "gwp", "ap", "ep", "ced", "category"}}
    """
    estimator = Scope3Estimator()

    impacts = {}

    # Material impacts (Category 1)
    for i, item in enumerate(BATTERY_BOM):
        pid = f"mat_{i:03d}"
        ef = estimator._find_factor(item["material"], item.get("supplier_location", "GLO"))
        factor = ef.factor if ef else 3.0  # default factor
        gwp = item["mass_kg"] * factor
        impacts[pid] = {
            "name": f"{item['component']} ({item['material']})",
            "gwp": round(gwp, 2),
            "ap": round(gwp * 0.008, 4),   # rough AP proxy
            "ep": round(gwp * 0.002, 4),   # rough EP proxy
            "ced": round(gwp * 12.0, 2),   # rough CED proxy (MJ)
            "unit": "kg CO2-eq",
            "category": "material",
        }

    # Energy impacts (Category 3)
    capacity_kwh = 75.0
    grid_factors = {
        "CN": 0.58,
        "DE": 0.35,
        "US": 0.42,
        "GLO": 0.45,
    }
    for stage, data in ENERGY_CONSUMPTION.items():
        pid = f"energy_{stage}"
        elec = data["electricity_kwh"] * capacity_kwh
        loc = data["location"]
        gf = grid_factors.get(loc, 0.45)
        gwp_elec = elec * gf
        gwp_heat = data["heat_mj"] * capacity_kwh * 0.056  # natural gas factor
        gwp_total = gwp_elec + gwp_heat
        impacts[pid] = {
            "name": f"Energy - {stage.replace('_', ' ').title()}",
            "gwp": round(gwp_total, 2),
            "ap": round(gwp_total * 0.006, 4),
            "ep": round(gwp_total * 0.001, 4),
            "ced": round(elec * 3.6 + data["heat_mj"] * capacity_kwh, 2),
            "unit": "kg CO2-eq",
            "category": "energy",
        }

    # Transport impacts (Category 4)
    for i, t in enumerate(TRANSPORT_DATA):
        pid = f"transport_{i:03d}"
        tkm = t["distance_km"] * t["mass_t"]
        ef = estimator._find_factor(t["mode"])
        factor = ef.factor if ef else 0.05
        gwp = tkm * factor
        impacts[pid] = {
            "name": f"Transport - {t['route']}",
            "gwp": round(gwp, 2),
            "ap": round(gwp * 0.01, 4),
            "ep": round(gwp * 0.003, 4),
            "ced": round(gwp * 15.0, 2),
            "unit": "kg CO2-eq",
            "category": "transport",
        }

    # End-of-life credits
    total_recyclable_mass = sum(
        item["mass_kg"]
        for item in BATTERY_BOM
        if item["material"] in ["copper", "aluminium", "nickel sulfate", "cobalt sulfate"]
    )
    eol_credit = total_recyclable_mass * (-0.8) * EOL_ASSUMPTIONS["recycling_rate"]
    impacts["eol_recycling"] = {
        "name": "End-of-Life Recycling Credits",
        "gwp": round(eol_credit, 2),
        "ap": round(eol_credit * 0.005, 4),
        "ep": round(eol_credit * 0.001, 4),
        "ced": round(eol_credit * 10.0, 2),
        "unit": "kg CO2-eq",
        "category": "end-of-life",
    }

    return impacts


def build_uncertainty_parameters() -> list[UncertainParameter]:
    """Build uncertain parameters for Monte Carlo simulation."""
    return [
        UncertainParameter("nickel_mass", 72.0, "lognormal",
                           {"mu": math.log(72.0), "sigma": 0.15}),
        UncertainParameter("cobalt_mass", 9.0, "lognormal",
                           {"mu": math.log(9.0), "sigma": 0.20}),
        UncertainParameter("lithium_mass", 18.5, "lognormal",
                           {"mu": math.log(18.5), "sigma": 0.15}),
        UncertainParameter("graphite_mass", 52.0, "normal",
                           {"mean": 52.0, "std": 5.0}),
        UncertainParameter("electricity_intensity", 42.0, "lognormal",
                           {"mu": math.log(42.0), "sigma": 0.25}),
        UncertainParameter("grid_carbon_intensity", 0.58, "normal",
                           {"mean": 0.58, "std": 0.05}),
        UncertainParameter("nickel_ef", 6.8, "lognormal",
                           {"mu": math.log(6.8), "sigma": 0.20}),
        UncertainParameter("cobalt_ef", 12.3, "lognormal",
                           {"mu": math.log(12.3), "sigma": 0.30}),
        UncertainParameter("recycling_rate", 0.70, "triangular",
                           {"low": 0.50, "mode": 0.70, "high": 0.90}),
        UncertainParameter("transport_distance", 8000, "uniform",
                           {"low": 5000, "high": 12000}),
    ]


def lca_model(**kwargs) -> float:
    """
    Simplified LCA model for uncertainty propagation.

    Computes total GWP (kg CO2-eq) for NMC811 75kWh battery pack.
    """
    # Material contributions
    ni_gwp = kwargs.get("nickel_mass", 72.0) * kwargs.get("nickel_ef", 6.8)
    co_gwp = kwargs.get("cobalt_mass", 9.0) * kwargs.get("cobalt_ef", 12.3)
    li_gwp = kwargs.get("lithium_mass", 18.5) * 7.5
    gr_gwp = kwargs.get("graphite_mass", 52.0) * 1.8
    al_gwp = (14.0 + 35.0 + 25.0 + 12.0) * 8.2
    cu_gwp = (18.0 + 3.5 + 5.0) * 3.5
    elec_gwp = 28.0 * 10.2
    sep_gwp = 8.0 * 4.5
    steel_gwp = 45.0 * 2.1
    plastic_gwp = 4.0 * 3.1

    material_total = (ni_gwp + co_gwp + li_gwp + gr_gwp + al_gwp + cu_gwp
                      + elec_gwp + sep_gwp + steel_gwp + plastic_gwp)

    # Energy contributions
    capacity = 75.0
    elec_int = kwargs.get("electricity_intensity", 42.0)
    grid_ci = kwargs.get("grid_carbon_intensity", 0.58)
    energy_gwp = elec_int * capacity * grid_ci

    # Additional energy stages
    energy_gwp += 3.5 * capacity * 0.35  # module assembly (DE grid)
    energy_gwp += 5.0 * capacity * 0.35  # pack assembly
    energy_gwp += 15.0 * capacity * 0.58  # dry room
    energy_gwp += (25.0 + 1.2 + 2.0) * capacity * 0.056  # heat

    # Transport
    transport_dist = kwargs.get("transport_distance", 8000)
    transport_gwp = transport_dist * 0.25 * 0.008  # ocean
    transport_gwp += 500 * 0.35 * 0.062  # truck 1
    transport_gwp += 200 * 0.42 * 0.062  # truck 2

    # End-of-life recycling credit
    recycling_rate = kwargs.get("recycling_rate", 0.70)
    recyclable_mass = 72.0 + 9.0 + 18.0 + 5.0 + 35.0 + 25.0 + 12.0
    eol_credit = recyclable_mass * (-0.8) * recycling_rate

    total = material_total + energy_gwp + transport_gwp + eol_credit
    return total


def define_scenarios() -> list[ScenarioDefinition]:
    """Define comparison scenarios for the EV battery LCA."""
    return [
        ScenarioDefinition(
            name="Renewable Energy Manufacturing",
            description="100% renewable electricity for cell manufacturing",
            parameter_changes={"energy": 0.15},  # 85% reduction in energy GWP
        ),
        ScenarioDefinition(
            name="LFP Chemistry Switch",
            description="Switch from NMC811 to LFP cathode chemistry",
            parameter_changes={
                "Cathode": 0.55,   # LFP has ~45% lower cathode impact
                "cobalt": 0.0,     # No cobalt needed
            },
        ),
        ScenarioDefinition(
            name="Closed-Loop Recycling",
            description="95% recycling rate with closed-loop material recovery",
            parameter_changes={"End-of-Life": 1.8},  # Larger recycling credits
        ),
        ScenarioDefinition(
            name="Localized Supply Chain",
            description="European supply chain (reduced transport, cleaner grid)",
            parameter_changes={
                "Transport": 0.30,  # 70% transport reduction
                "energy": 0.60,     # Cleaner EU grid
            },
        ),
        ScenarioDefinition(
            name="Best Case 2030",
            description="Combined improvements: renewable energy + recycling + local supply",
            parameter_changes={
                "energy": 0.20,
                "Transport": 0.25,
                "End-of-Life": 2.0,
            },
        ),
    ]


def run_case_study() -> dict:
    """Execute the complete EV battery LCA case study."""
    results = {}

    # Step 1: Build process tree
    tree = build_process_tree()
    results["process_tree"] = tree.to_dict()
    results["total_mass_kg"] = sum(item["mass_kg"] for item in BATTERY_BOM)

    # Step 2: Compute impacts
    impacts = compute_process_impacts()
    results["process_impacts"] = impacts

    total_gwp = sum(v["gwp"] for v in impacts.values())
    total_ced = sum(v["ced"] for v in impacts.values())
    results["total_gwp_kg_co2eq"] = round(total_gwp, 2)
    results["gwp_per_kwh"] = round(total_gwp / 75.0, 2)
    results["total_ced_mj"] = round(total_ced, 2)

    # Step 3: Uncertainty analysis
    params = build_uncertainty_parameters()

    mc = MonteCarloSimulator(seed=42)
    mc_result = mc.run(lca_model, params, n_iterations=10000)
    results["monte_carlo"] = mc_result.to_dict()

    taylor = TaylorExpansion()
    taylor_result = taylor.propagate(lca_model, params)
    results["taylor_expansion"] = taylor_result.to_dict()

    # Step 4: Hotspot analysis
    analyzer = HotspotAnalyzer()
    hotspots = analyzer.analyze(impacts, "GWP")
    results["hotspots"] = [
        {
            "process": h.process_name,
            "contribution_pct": h.contribution_pct,
            "absolute_gwp": h.absolute_value,
            "improvement": h.improvement_potential,
            "alternatives": h.suggested_alternatives,
        }
        for h in hotspots
    ]

    # Step 5: Scenario comparison
    comparator = ScenarioComparator()
    scenarios = define_scenarios()
    comparison = comparator.compare(impacts, scenarios)
    results["scenario_comparison"] = {
        "baseline_gwp": comparison.scenarios[0].total_gwp,
        "scenarios": [
            {
                "name": s.scenario_name,
                "gwp": s.total_gwp,
            }
            for s in comparison.scenarios
        ],
        "relative_changes": comparison.relative_changes,
    }

    # Step 6: Scope 3 estimation
    estimator = Scope3Estimator()
    scope3_data = {
        1: [{"source": item["material"], "amount": item["mass_kg"], "unit": "kg"}
            for item in BATTERY_BOM],
        3: [
            {"source": "electricity (China grid)",
             "amount": (42.0 + 15.0) * 75.0, "unit": "kWh"},
            {"source": "electricity (EU grid)",
             "amount": (3.5 + 5.0) * 75.0, "unit": "kWh"},
            {"source": "natural gas", "amount": 28.2 * 75.0 / 38.0, "unit": "m3"},
        ],
        4: [{"source": t["mode"], "amount": t["distance_km"] * t["mass_t"], "unit": "tkm"}
            for t in TRANSPORT_DATA],
        5: [{"source": "landfill (mixed waste)", "amount": 15.0, "unit": "kg"},
            {"source": "incineration", "amount": 5.0, "unit": "kg"}],
        12: [{"source": "battery recycling (hydromet)",
              "amount": results["total_mass_kg"] * 0.95 * 0.70 * 0.60, "unit": "kg"},
             {"source": "battery recycling (pyromet)",
              "amount": results["total_mass_kg"] * 0.95 * 0.70 * 0.40, "unit": "kg"}],
    }

    scope3_report = estimator.estimate_full_scope3(
        scope3_data, region="GLO",
        company_name="EV Battery Manufacturer (NMC811)",
        year=2024,
    )
    results["scope3"] = scope3_report.summary()

    return results


if __name__ == "__main__":
    results = run_case_study()
    output_path = os.path.join(os.path.dirname(__file__), "..", "..", "results", "ev_battery_lca_results.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results saved to {output_path}")
    print(f"Total GWP: {results['total_gwp_kg_co2eq']} kg CO2-eq")
    print(f"GWP per kWh: {results['gwp_per_kwh']} kg CO2-eq/kWh")
