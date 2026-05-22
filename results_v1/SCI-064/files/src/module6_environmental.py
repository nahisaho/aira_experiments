"""
Module 6: Environmental Pollutant Detection Application
========================================================
Integration of all modules for environmental monitoring biosensors:
- Multi-analyte detection panel design
- Cross-reactivity analysis
- Field deployment considerations
- Regulatory compliance (EPA/WHO limits)
- Performance benchmarking against standard methods
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import json
import os


@dataclass
class EnvironmentalStandard:
    """Regulatory standard for an environmental pollutant."""
    analyte: str
    epa_mcl_ppb: float          # EPA Maximum Contaminant Level
    who_guideline_ppb: float    # WHO guideline value
    eu_standard_ppb: float      # EU drinking water standard
    molecular_weight: float     # g/mol
    
    @property
    def epa_mcl_uM(self) -> float:
        return self.epa_mcl_ppb / self.molecular_weight * 1000
    
    @property
    def who_guideline_uM(self) -> float:
        return self.who_guideline_ppb / self.molecular_weight * 1000


# Environmental standards database
STANDARDS = {
    "Mercury": EnvironmentalStandard("Mercury(II)", 2.0, 6.0, 1.0, 200.59),
    "Arsenic": EnvironmentalStandard("Arsenic(III)", 10.0, 10.0, 10.0, 74.92),
    "Cadmium": EnvironmentalStandard("Cadmium(II)", 5.0, 3.0, 5.0, 112.41),
    "Lead": EnvironmentalStandard("Lead(II)", 15.0, 10.0, 10.0, 207.20),
    "Chromium": EnvironmentalStandard("Chromium(VI)", 100.0, 50.0, 50.0, 52.00),
    "Copper": EnvironmentalStandard("Copper(II)", 1300.0, 2000.0, 2000.0, 63.55),
    "Zinc": EnvironmentalStandard("Zinc(II)", 5000.0, 3000.0, None, 65.38),
    "Nickel": EnvironmentalStandard("Nickel(II)", 100.0, 70.0, 20.0, 58.69),
    "Naphthalene": EnvironmentalStandard("Naphthalene", 100.0, None, None, 128.17),
    "Phenol": EnvironmentalStandard("Phenol", 300.0, None, None, 94.11),
    "Benzene": EnvironmentalStandard("Benzene", 5.0, 10.0, 1.0, 78.11),
    "Toluene": EnvironmentalStandard("Toluene", 1000.0, 700.0, None, 92.14),
}


@dataclass
class BiosensorSpec:
    """Complete biosensor specification for field deployment."""
    name: str
    target_analyte: str
    tf_type: str
    LOD_uM: float
    LOD_ppb: float
    linear_range_uM: Tuple[float, float]
    dynamic_range_fold: float
    response_time_min: float
    selectivity: Dict[str, float]
    meets_epa: bool
    meets_who: bool
    circuit_type: str
    key_mutations: List[str]


def design_biosensor_panel() -> Dict[str, BiosensorSpec]:
    """Design complete biosensor panel for environmental monitoring."""
    np.random.seed(42)
    
    panel = {}
    
    # Heavy metal biosensors
    metal_configs = [
        ("Hg_Sensor", "Mercury", "MerR", 0.001, 1.8, 100,
         ["C82S_sensitivity", "L117F_selectivity"], "cascade"),
        ("As_Sensor", "Arsenic", "ArsR", 0.01, 1.5, 80,
         ["C32A_tuning", "D37E_stability"], "positive_feedback"),
        ("Cd_Sensor", "Cadmium", "CadC", 0.005, 2.0, 60,
         ["C7S_specificity", "C60A_range"], "simple"),
        ("Pb_Sensor", "Lead", "CadC", 0.008, 1.6, 50,
         ["D11N_Pb_preference", "H58A_selectivity"], "cascade"),
        ("Cu_Sensor", "Copper", "CueR", 0.05, 2.2, 120,
         ["C112S_tuning", "M120L_dynamic_range"], "simple"),
        ("Zn_Sensor", "Zinc", "SmtB", 0.1, 1.3, 70,
         ["H97A_sensitivity", "C61S_specificity"], "simple"),
        ("Cr_Sensor", "Chromium", "CadC", 0.02, 1.4, 40,
         ["E60D_Cr_binding", "C11A_selectivity"], "cascade"),
        ("Ni_Sensor", "Nickel", "SmtB", 0.03, 1.5, 55,
         ["H14A_Ni_preference", "D64E_tuning"], "positive_feedback"),
    ]
    
    # Organic pollutant biosensors
    organic_configs = [
        ("Naphthalene_Sensor", "Naphthalene", "NahR", 0.05, 1.1, 50,
         ["F165W_pocket_reshape", "L201V_affinity"], "cascade"),
        ("Phenol_Sensor", "Phenol", "DmpR", 0.1, 1.4, 45,
         ["W148F_specificity", "Y201A_sensitivity"], "positive_feedback"),
        ("Benzene_Sensor", "Benzene", "DmpR", 0.005, 1.3, 35,
         ["F160W_benzene_binding", "I180V_selectivity"], "cascade"),
        ("Toluene_Sensor", "Toluene", "DmpR", 0.2, 1.2, 40,
         ["L175V_toluene_fit", "M200I_range"], "simple"),
    ]
    
    for config in metal_configs + organic_configs:
        name, analyte, tf, lod_uM, hill, dr_fold, mutations, circuit = config
        
        std = STANDARDS.get(analyte)
        if std is None:
            continue
        
        lod_ppb = lod_uM * std.molecular_weight / 1000
        linear_max = lod_uM * dr_fold
        
        # Cross-reactivity analysis
        selectivity = {}
        for other_analyte in STANDARDS:
            if other_analyte == analyte:
                selectivity[other_analyte] = 1.0
            elif other_analyte in ["Mercury", "Arsenic", "Cadmium", "Lead", "Copper", "Zinc", "Chromium", "Nickel"]:
                if analyte in ["Mercury", "Arsenic", "Cadmium", "Lead", "Copper", "Zinc", "Chromium", "Nickel"]:
                    selectivity[other_analyte] = round(np.random.uniform(0.001, 0.15), 4)
                else:
                    selectivity[other_analyte] = round(np.random.uniform(0.0001, 0.01), 4)
            else:
                if analyte in ["Naphthalene", "Phenol", "Benzene", "Toluene"]:
                    selectivity[other_analyte] = round(np.random.uniform(0.01, 0.2), 4)
                else:
                    selectivity[other_analyte] = round(np.random.uniform(0.0001, 0.01), 4)
        
        response_time = 30 + np.random.uniform(0, 90) if circuit == "simple" else 60 + np.random.uniform(0, 120)
        
        meets_epa = lod_ppb < std.epa_mcl_ppb * 0.5
        meets_who = lod_ppb < (std.who_guideline_ppb * 0.5 if std.who_guideline_ppb else float('inf'))
        
        panel[name] = BiosensorSpec(
            name=name,
            target_analyte=analyte,
            tf_type=tf,
            LOD_uM=round(lod_uM, 6),
            LOD_ppb=round(lod_ppb, 4),
            linear_range_uM=(round(lod_uM * 3, 6), round(linear_max, 4)),
            dynamic_range_fold=dr_fold,
            response_time_min=round(response_time, 1),
            selectivity=selectivity,
            meets_epa=meets_epa,
            meets_who=meets_who,
            circuit_type=circuit,
            key_mutations=mutations,
        )
    
    return panel


def benchmark_against_standards(panel: Dict[str, BiosensorSpec]) -> Dict:
    """Compare biosensor performance against gold-standard methods."""
    standard_methods = {
        "ICP-MS": {"LOD_ppb": 0.001, "cost_per_sample": 50, "time_hours": 4,
                    "portability": "lab_only", "throughput": "medium"},
        "AAS": {"LOD_ppb": 0.1, "cost_per_sample": 20, "time_hours": 2,
                "portability": "lab_only", "throughput": "medium"},
        "XRF": {"LOD_ppb": 100, "cost_per_sample": 5, "time_hours": 0.1,
                "portability": "portable", "throughput": "high"},
        "Colorimetric_kit": {"LOD_ppb": 10, "cost_per_sample": 2, "time_hours": 0.5,
                              "portability": "field", "throughput": "high"},
        "GC-MS": {"LOD_ppb": 0.01, "cost_per_sample": 100, "time_hours": 6,
                   "portability": "lab_only", "throughput": "low"},
    }
    
    benchmarks = {}
    for sensor_name, spec in panel.items():
        comparisons = {}
        for method_name, method in standard_methods.items():
            lod_ratio = spec.LOD_ppb / method["LOD_ppb"] if method["LOD_ppb"] > 0 else float('inf')
            comparisons[method_name] = {
                "LOD_ratio": round(lod_ratio, 2),
                "cost_advantage": f"${method['cost_per_sample']} vs ~$0.10",
                "time_advantage": f"{method['time_hours']}h vs {spec.response_time_min/60:.1f}h",
                "portability_advantage": method["portability"] != "field",
            }
        benchmarks[sensor_name] = comparisons
    
    return benchmarks


def field_deployment_analysis(panel: Dict[str, BiosensorSpec]) -> Dict:
    """Analyze field deployment considerations."""
    deployment = {}
    
    for name, spec in panel.items():
        # Environmental matrix effects
        matrix_factors = {
            "pH_sensitivity": "moderate" if spec.tf_type in ["MerR", "CueR"] else "high",
            "temperature_range_C": (4, 42),
            "ionic_strength_tolerance_mM": 500 if "metal" in name.lower() or spec.target_analyte in
                ["Mercury","Arsenic","Cadmium","Lead","Copper","Zinc","Chromium","Nickel"] else 200,
            "shelf_life_days": 30 if spec.circuit_type == "simple" else 14,
            "sample_preparation": "minimal" if spec.tf_type in ["MerR", "CueR"] else "filtration_required",
        }
        
        # Multiplexing compatibility
        multiplex_compatible = []
        for other_name, other_spec in panel.items():
            if other_name != name and other_spec.target_analyte != spec.target_analyte:
                max_cross = max(spec.selectivity.get(other_spec.target_analyte, 0),
                               other_spec.selectivity.get(spec.target_analyte, 0))
                if max_cross < 0.05:
                    multiplex_compatible.append(other_name)
        
        deployment[name] = {
            "matrix_factors": matrix_factors,
            "multiplex_compatible_sensors": multiplex_compatible[:5],
            "recommended_applications": _get_applications(spec),
            "regulatory_compliance": {
                "EPA_compliant": spec.meets_epa,
                "WHO_compliant": spec.meets_who,
                "required_LOD_ppb": STANDARDS.get(spec.target_analyte, 
                    EnvironmentalStandard("", 100, 100, 100, 100)).epa_mcl_ppb,
                "achieved_LOD_ppb": spec.LOD_ppb,
            }
        }
    
    return deployment


def _get_applications(spec: BiosensorSpec) -> List[str]:
    """Determine recommended applications based on sensor specs."""
    apps = []
    std = STANDARDS.get(spec.target_analyte)
    if std is None:
        return ["General screening"]
    
    if spec.meets_epa:
        apps.append("Drinking water monitoring")
    if spec.LOD_ppb < std.epa_mcl_ppb:
        apps.append("Wastewater effluent testing")
    if spec.response_time_min < 60:
        apps.append("Real-time field monitoring")
    if spec.dynamic_range_fold > 50:
        apps.append("Industrial process control")
    apps.append("Environmental screening")
    if spec.target_analyte in ["Mercury", "Lead", "Arsenic"]:
        apps.append("Contaminated site assessment")
    
    return apps


def run_environmental_application(output_dir: str = "results") -> Dict:
    """Run complete environmental application analysis."""
    os.makedirs(output_dir, exist_ok=True)
    
    print("  Designing biosensor panel...")
    panel = design_biosensor_panel()
    
    print("  Benchmarking against standards...")
    benchmarks = benchmark_against_standards(panel)
    
    print("  Field deployment analysis...")
    deployment = field_deployment_analysis(panel)
    
    # Compile results
    results = {
        "panel_summary": {
            "total_sensors": len(panel),
            "heavy_metal_sensors": sum(1 for s in panel.values()
                                        if s.target_analyte in ["Mercury","Arsenic","Cadmium","Lead","Copper","Zinc","Chromium","Nickel"]),
            "organic_sensors": sum(1 for s in panel.values()
                                    if s.target_analyte in ["Naphthalene","Phenol","Benzene","Toluene"]),
            "epa_compliant": sum(1 for s in panel.values() if s.meets_epa),
            "who_compliant": sum(1 for s in panel.values() if s.meets_who),
        },
        "sensors": {},
        "benchmarks": benchmarks,
        "deployment": deployment,
    }
    
    for name, spec in panel.items():
        results["sensors"][name] = {
            "target": spec.target_analyte,
            "tf_type": spec.tf_type,
            "LOD_uM": spec.LOD_uM,
            "LOD_ppb": spec.LOD_ppb,
            "linear_range_uM": list(spec.linear_range_uM),
            "dynamic_range_fold": spec.dynamic_range_fold,
            "response_time_min": spec.response_time_min,
            "circuit_type": spec.circuit_type,
            "key_mutations": spec.key_mutations,
            "meets_epa": spec.meets_epa,
            "meets_who": spec.meets_who,
            "max_cross_reactivity": round(max(
                v for k, v in spec.selectivity.items() if k != spec.target_analyte
            ), 4) if spec.selectivity else 0,
        }
    
    with open(os.path.join(output_dir, "environmental_application.json"), 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    return results


if __name__ == "__main__":
    results = run_environmental_application()
    print(f"\nPanel: {results['panel_summary']['total_sensors']} sensors")
    print(f"EPA compliant: {results['panel_summary']['epa_compliant']}")
    print(f"WHO compliant: {results['panel_summary']['who_compliant']}")
    
    for name, data in results["sensors"].items():
        print(f"\n  {name}: LOD = {data['LOD_ppb']:.2f} ppb, "
              f"DR = {data['dynamic_range_fold']:.0f}x, "
              f"EPA: {'✓' if data['meets_epa'] else '✗'}")
