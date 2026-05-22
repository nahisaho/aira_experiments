#!/usr/bin/env python3
"""
InVEST/ARIES-based Ecosystem Service Evaluation Pipeline
=========================================================
里山生態系サービスの空間的定量化パイプライン

Pipeline stages:
  1. Data preparation (land cover, DEM, climate, soil)
  2. InVEST model execution (Carbon, SDR, NDR, Water Yield, Habitat Quality, Pollination, Recreation)
  3. ARIES k.LAB integration (flood regulation, cultural services)
  4. Result aggregation and monetary valuation
  5. SEEA-EA account compilation
  6. Sensitivity analysis and uncertainty quantification
"""

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import math


# ==============================================================================
# 1. Configuration & Data Classes
# ==============================================================================

@dataclass
class StudyArea:
    """里山ケーススタディ対象地域"""
    name: str = "Noto Peninsula Satoyama (GIAHS)"
    total_area_ha: float = 5000.0
    crs: str = "EPSG:6674"  # JGD2011 / Japan Plane Rectangular CS VII
    bbox: Tuple[float, float, float, float] = (136.6, 37.0, 137.3, 37.5)
    dem_resolution_m: float = 30.0
    land_cover_classes: Dict[str, float] = field(default_factory=lambda: {
        "mixed_forest": 1500.0,
        "paddy_field": 1200.0,
        "upland_crop": 500.0,
        "bamboo_grove": 300.0,
        "irrigation_pond": 200.0,
        "grassland_meadow": 400.0,
        "settlement": 600.0,
        "other": 300.0,
    })


@dataclass
class CarbonPools:
    """Carbon density by land cover type (tC/ha)"""
    above_ground: float
    below_ground: float
    soil: float
    dead_matter: float

    @property
    def total(self) -> float:
        return self.above_ground + self.below_ground + self.soil + self.dead_matter


# Literature-based carbon pool values for Satoyama landscape
CARBON_POOLS: Dict[str, CarbonPools] = {
    "mixed_forest":     CarbonPools(82.0, 26.0, 95.0, 12.0),
    "paddy_field":      CarbonPools(3.5,  0.5,  65.0,  1.0),
    "upland_crop":      CarbonPools(5.0,  1.5,  55.0,  0.5),
    "bamboo_grove":     CarbonPools(45.0, 15.0, 70.0,  5.0),
    "irrigation_pond":  CarbonPools(0.0,  0.0,  40.0,  3.0),
    "grassland_meadow": CarbonPools(8.0,  12.0, 80.0,  2.0),
    "settlement":       CarbonPools(15.0,  3.0, 30.0,  1.0),
    "other":            CarbonPools(5.0,  2.0,  40.0,  1.0),
}

# Biophysical parameters
USLE_PARAMS = {
    "mixed_forest":     {"C": 0.003, "P": 1.0, "K": 0.028},
    "paddy_field":      {"C": 0.15,  "P": 0.10, "K": 0.035},
    "upland_crop":      {"C": 0.30,  "P": 0.50, "K": 0.035},
    "bamboo_grove":     {"C": 0.01,  "P": 1.0, "K": 0.028},
    "irrigation_pond":  {"C": 0.0,   "P": 1.0, "K": 0.020},
    "grassland_meadow": {"C": 0.01,  "P": 1.0, "K": 0.030},
    "settlement":       {"C": 0.0,   "P": 1.0, "K": 0.025},
    "other":            {"C": 0.10,  "P": 1.0, "K": 0.030},
}

# Water yield parameters (Budyko curve)
WATER_PARAMS = {
    "precipitation_mm": 2100,  # Annual precipitation (Noto region)
    "pet_mm": 800,  # Potential evapotranspiration
    "Kc": {  # Crop coefficients
        "mixed_forest": 1.0, "paddy_field": 1.2, "upland_crop": 0.9,
        "bamboo_grove": 1.0, "irrigation_pond": 1.1, "grassland_meadow": 0.8,
        "settlement": 0.3, "other": 0.6,
    },
    "Z_parameter": 11.0,  # Zhang parameter (seasonality)
}

# Nutrient loading/retention parameters (kg/ha/year)
NUTRIENT_PARAMS = {
    "mixed_forest":     {"N_load": 2.0,  "N_eff": 0.90, "P_load": 0.2, "P_eff": 0.85},
    "paddy_field":      {"N_load": 35.0, "N_eff": 0.40, "P_load": 5.0, "P_eff": 0.35},
    "upland_crop":      {"N_load": 45.0, "N_eff": 0.25, "P_load": 8.0, "P_eff": 0.20},
    "bamboo_grove":     {"N_load": 3.0,  "N_eff": 0.85, "P_load": 0.3, "P_eff": 0.80},
    "irrigation_pond":  {"N_load": 5.0,  "N_eff": 0.70, "P_load": 1.0, "P_eff": 0.65},
    "grassland_meadow": {"N_load": 8.0,  "N_eff": 0.75, "P_load": 1.0, "P_eff": 0.70},
    "settlement":       {"N_load": 25.0, "N_eff": 0.10, "P_load": 4.0, "P_eff": 0.08},
    "other":            {"N_load": 10.0, "N_eff": 0.50, "P_load": 2.0, "P_eff": 0.45},
}

# Habitat quality parameters
HABITAT_PARAMS = {
    "sensitivity": {
        "mixed_forest": {"habitat": True, "quality": 1.0,
                         "threats": {"settlement": 0.8, "upland_crop": 0.5, "paddy_field": 0.3}},
        "paddy_field": {"habitat": True, "quality": 0.5,
                        "threats": {"settlement": 0.7, "upland_crop": 0.3}},
        "bamboo_grove": {"habitat": True, "quality": 0.6,
                         "threats": {"settlement": 0.6}},
        "grassland_meadow": {"habitat": True, "quality": 0.7,
                             "threats": {"settlement": 0.7, "upland_crop": 0.4}},
        "irrigation_pond": {"habitat": True, "quality": 0.8,
                            "threats": {"settlement": 0.9, "upland_crop": 0.5, "paddy_field": 0.2}},
    },
    "threats": {
        "settlement": {"max_dist_km": 5.0, "weight": 1.0, "decay": "exponential"},
        "upland_crop": {"max_dist_km": 3.0, "weight": 0.7, "decay": "linear"},
        "paddy_field": {"max_dist_km": 2.0, "weight": 0.3, "decay": "linear"},
    },
    "half_saturation": 0.5,
}

# Monetary valuation unit prices
UNIT_PRICES = {
    "carbon_social_cost": 12000,    # JPY/tCO2 (Japan Carbon Tax 2025 level)
    "carbon_market_price": 3000,    # JPY/tCO2 (J-Credit average)
    "water_supply_price": 200,      # JPY/m3
    "erosion_damage_cost": 5000,    # JPY/ton sediment
    "N_treatment_cost": 2500,       # JPY/kg-N removed
    "P_treatment_cost": 12000,      # JPY/kg-P removed
    "recreation_value": 3500,       # JPY/person-day
    "pollination_value": 150000,    # JPY/ha cropland/year
    "flood_damage_avoided": 50000,  # JPY/ha/year (expected annual damage)
}

# Declining Discount Rate schedule
DDR_SCHEDULE = {
    (0, 30): 0.035,
    (31, 75): 0.025,
    (76, 125): 0.020,
    (126, 200): 0.015,
    (201, 300): 0.010,
    (301, 999): 0.005,
}


# ==============================================================================
# 2. InVEST Model Simulation Functions
# ==============================================================================

def compute_carbon_storage(study_area: StudyArea) -> Dict:
    """InVEST Carbon Storage and Sequestration Model (simplified)"""
    results = {}
    total_storage_tC = 0.0

    for lc, area_ha in study_area.land_cover_classes.items():
        pools = CARBON_POOLS.get(lc)
        if pools is None:
            continue
        storage_tC = pools.total * area_ha
        total_storage_tC += storage_tC
        results[lc] = {
            "area_ha": area_ha,
            "density_tC_per_ha": pools.total,
            "above_ground_tC": pools.above_ground * area_ha,
            "below_ground_tC": pools.below_ground * area_ha,
            "soil_tC": pools.soil * area_ha,
            "dead_matter_tC": pools.dead_matter * area_ha,
            "total_tC": storage_tC,
        }

    total_tCO2 = total_storage_tC * 44 / 12
    results["summary"] = {
        "total_storage_tC": round(total_storage_tC, 1),
        "total_storage_tCO2": round(total_tCO2, 1),
        "mean_density_tC_per_ha": round(total_storage_tC / study_area.total_area_ha, 1),
        "monetary_value_social_JPY": round(total_tCO2 * UNIT_PRICES["carbon_social_cost"]),
        "monetary_value_market_JPY": round(total_tCO2 * UNIT_PRICES["carbon_market_price"]),
    }
    return results


def compute_sediment_delivery(study_area: StudyArea) -> Dict:
    """InVEST SDR Model (simplified USLE-based)"""
    R_factor = 5500  # Rainfall erosivity (MJ·mm/ha/h/year) for Noto
    LS_factor = 2.5  # Average slope-length factor (moderate terrain)

    results = {}
    total_potential = 0.0
    total_actual = 0.0

    for lc, area_ha in study_area.land_cover_classes.items():
        params = USLE_PARAMS.get(lc, {"C": 0.1, "P": 1.0, "K": 0.03})
        potential = R_factor * params["K"] * LS_factor * area_ha  # C=1, P=1
        actual = potential * params["C"] * params["P"]
        avoided = potential - actual
        total_potential += potential
        total_actual += actual
        results[lc] = {
            "potential_erosion_ton": round(potential, 1),
            "actual_erosion_ton": round(actual, 1),
            "avoided_erosion_ton": round(avoided, 1),
            "retention_pct": round((1 - params["C"] * params["P"]) * 100, 1),
        }

    total_avoided = total_potential - total_actual
    results["summary"] = {
        "total_potential_erosion_ton": round(total_potential, 1),
        "total_actual_erosion_ton": round(total_actual, 1),
        "total_avoided_erosion_ton": round(total_avoided, 1),
        "overall_retention_pct": round((1 - total_actual / total_potential) * 100, 1),
        "monetary_value_JPY": round(total_avoided * UNIT_PRICES["erosion_damage_cost"]),
    }
    return results


def compute_water_yield(study_area: StudyArea) -> Dict:
    """InVEST Annual Water Yield (Budyko curve)"""
    P = WATER_PARAMS["precipitation_mm"]
    PET = WATER_PARAMS["pet_mm"]
    Z = WATER_PARAMS["Z_parameter"]

    results = {}
    total_yield_m3 = 0.0

    for lc, area_ha in study_area.land_cover_classes.items():
        Kc = WATER_PARAMS["Kc"].get(lc, 0.6)
        AET_PET_ratio = Kc * PET / P
        # Budyko-Zhang equation: AET/P = 1 + w - [1 + w^ω]^(1/ω)
        # Simplified: AET/P ≈ (1 + AET_PET_ratio - (1 + AET_PET_ratio**Z)**(1/Z))
        omega = Z
        budyko = 1 + AET_PET_ratio - (1 + AET_PET_ratio**omega)**(1.0/omega)
        aet_frac = min(max(budyko, 0.0), 1.0)
        water_yield_mm = P * (1 - aet_frac)
        water_yield_m3 = water_yield_mm * area_ha * 10  # mm*ha -> m3
        total_yield_m3 += water_yield_m3
        results[lc] = {
            "aet_fraction": round(aet_frac, 3),
            "water_yield_mm": round(water_yield_mm, 1),
            "water_yield_m3": round(water_yield_m3, 0),
        }

    results["summary"] = {
        "total_water_yield_m3": round(total_yield_m3, 0),
        "mean_yield_mm": round(total_yield_m3 / (study_area.total_area_ha * 10), 1),
        "monetary_value_JPY": round(total_yield_m3 * UNIT_PRICES["water_supply_price"]),
    }
    return results


def compute_nutrient_retention(study_area: StudyArea) -> Dict:
    """InVEST NDR Model (simplified)"""
    results = {}
    total_N_load = total_N_retained = total_P_load = total_P_retained = 0.0

    for lc, area_ha in study_area.land_cover_classes.items():
        params = NUTRIENT_PARAMS.get(lc, {"N_load": 5, "N_eff": 0.5, "P_load": 1, "P_eff": 0.5})
        n_load = params["N_load"] * area_ha
        n_retained = n_load * params["N_eff"]
        p_load = params["P_load"] * area_ha
        p_retained = p_load * params["P_eff"]

        total_N_load += n_load
        total_N_retained += n_retained
        total_P_load += p_load
        total_P_retained += p_retained

        results[lc] = {
            "N_load_kg": round(n_load, 1),
            "N_retained_kg": round(n_retained, 1),
            "N_export_kg": round(n_load - n_retained, 1),
            "P_load_kg": round(p_load, 1),
            "P_retained_kg": round(p_retained, 1),
            "P_export_kg": round(p_load - p_retained, 1),
        }

    results["summary"] = {
        "total_N_load_kg": round(total_N_load, 1),
        "total_N_retained_kg": round(total_N_retained, 1),
        "total_N_export_kg": round(total_N_load - total_N_retained, 1),
        "N_retention_pct": round(total_N_retained / total_N_load * 100, 1),
        "total_P_load_kg": round(total_P_load, 1),
        "total_P_retained_kg": round(total_P_retained, 1),
        "total_P_export_kg": round(total_P_load - total_P_retained, 1),
        "P_retention_pct": round(total_P_retained / total_P_load * 100, 1),
        "monetary_value_N_JPY": round(total_N_retained * UNIT_PRICES["N_treatment_cost"]),
        "monetary_value_P_JPY": round(total_P_retained * UNIT_PRICES["P_treatment_cost"]),
        "monetary_value_total_JPY": round(
            total_N_retained * UNIT_PRICES["N_treatment_cost"]
            + total_P_retained * UNIT_PRICES["P_treatment_cost"]
        ),
    }
    return results


def compute_habitat_quality(study_area: StudyArea) -> Dict:
    """InVEST Habitat Quality Model (simplified)"""
    results = {}
    total_quality_ha = 0.0
    total_habitat_ha = 0.0

    for lc, area_ha in study_area.land_cover_classes.items():
        hparams = HABITAT_PARAMS["sensitivity"].get(lc)
        if hparams and hparams["habitat"]:
            # Simplified degradation based on proximity to threats
            threat_score = 0.0
            for threat, sensitivity in hparams["threats"].items():
                threat_area = study_area.land_cover_classes.get(threat, 0)
                threat_proportion = threat_area / study_area.total_area_ha
                threat_score += sensitivity * threat_proportion
            degradation = min(threat_score, 1.0)
            k = HABITAT_PARAMS["half_saturation"]
            quality = hparams["quality"] * (1 - (degradation**2.5 / (degradation**2.5 + k**2.5)))
            total_quality_ha += quality * area_ha
            total_habitat_ha += area_ha
            results[lc] = {
                "is_habitat": True,
                "base_quality": hparams["quality"],
                "degradation_score": round(degradation, 3),
                "quality_index": round(quality, 3),
                "quality_weighted_ha": round(quality * area_ha, 1),
            }
        else:
            results[lc] = {
                "is_habitat": False,
                "quality_index": 0.0,
                "quality_weighted_ha": 0.0,
            }

    mean_quality = total_quality_ha / total_habitat_ha if total_habitat_ha > 0 else 0.0
    results["summary"] = {
        "total_habitat_area_ha": round(total_habitat_ha, 1),
        "total_quality_weighted_ha": round(total_quality_ha, 1),
        "mean_habitat_quality": round(mean_quality, 3),
        "habitat_coverage_pct": round(total_habitat_ha / study_area.total_area_ha * 100, 1),
    }
    return results


def compute_pollination(study_area: StudyArea) -> Dict:
    """InVEST Crop Pollination Model (simplified)"""
    nesting_suitability = {
        "mixed_forest": 0.9, "bamboo_grove": 0.6, "grassland_meadow": 0.8,
        "paddy_field": 0.1, "upland_crop": 0.1, "irrigation_pond": 0.0,
        "settlement": 0.2, "other": 0.3,
    }
    floral_resources = {
        "mixed_forest": 0.7, "bamboo_grove": 0.3, "grassland_meadow": 0.9,
        "paddy_field": 0.2, "upland_crop": 0.5, "irrigation_pond": 0.1,
        "settlement": 0.3, "other": 0.4,
    }

    total_abundance = 0.0
    total_area = 0.0
    results = {}
    for lc, area_ha in study_area.land_cover_classes.items():
        ns = nesting_suitability.get(lc, 0.3)
        fr = floral_resources.get(lc, 0.3)
        abundance = (ns * fr) ** 0.5  # geometric mean
        total_abundance += abundance * area_ha
        total_area += area_ha
        results[lc] = {
            "nesting_suitability": ns,
            "floral_resources": fr,
            "pollinator_abundance_index": round(abundance, 3),
        }

    crop_area = study_area.land_cover_classes.get("upland_crop", 0)
    mean_abundance = total_abundance / total_area if total_area > 0 else 0.0
    results["summary"] = {
        "mean_pollinator_abundance": round(mean_abundance, 3),
        "pollination_dependent_crop_ha": crop_area,
        "monetary_value_JPY": round(crop_area * UNIT_PRICES["pollination_value"] * mean_abundance),
    }
    return results


def compute_recreation(study_area: StudyArea) -> Dict:
    """InVEST Recreation Model (simplified)"""
    attractiveness = {
        "mixed_forest": 0.8, "paddy_field": 0.5, "bamboo_grove": 0.6,
        "irrigation_pond": 0.7, "grassland_meadow": 0.6, "settlement": 0.1,
        "upland_crop": 0.3, "other": 0.2,
    }
    base_visits_per_ha = 5.0  # person-days/ha/year (regional average)

    total_visits = 0.0
    results = {}
    for lc, area_ha in study_area.land_cover_classes.items():
        attr = attractiveness.get(lc, 0.3)
        visits = base_visits_per_ha * attr * area_ha
        total_visits += visits
        results[lc] = {
            "attractiveness": attr,
            "estimated_visits": round(visits, 0),
        }

    results["summary"] = {
        "total_annual_visits": round(total_visits, 0),
        "visits_per_ha": round(total_visits / study_area.total_area_ha, 1),
        "monetary_value_JPY": round(total_visits * UNIT_PRICES["recreation_value"]),
    }
    return results


# ==============================================================================
# 3. Discount Rate & NPV Calculations
# ==============================================================================

def get_ddr(year: int) -> float:
    """Get declining discount rate for a given year"""
    for (y_start, y_end), rate in DDR_SCHEDULE.items():
        if y_start <= year <= y_end:
            return rate
    return 0.005


def compute_npv(annual_flow: float, horizon_years: int = 25,
                rate_type: str = "ddr") -> Dict:
    """Compute NPV of ecosystem service flows"""
    rates = {"constant_low": 0.014, "constant_mid": 0.035, "constant_high": 0.05}
    npv = 0.0
    discount_factors = []

    for t in range(horizon_years):
        if rate_type == "ddr":
            r = get_ddr(t)
        else:
            r = rates.get(rate_type, 0.035)
        df = 1.0 / (1.0 + r) ** t
        npv += annual_flow * df
        if t % 5 == 0:
            discount_factors.append({"year": t, "rate": r, "df": round(df, 4)})

    return {
        "annual_flow_JPY": annual_flow,
        "horizon_years": horizon_years,
        "rate_type": rate_type,
        "npv_JPY": round(npv),
        "sample_discount_factors": discount_factors,
    }


# ==============================================================================
# 4. SEEA-EA Account Compilation
# ==============================================================================

def compile_seea_accounts(study_area: StudyArea, service_results: Dict) -> Dict:
    """Compile SEEA-EA compatible ecosystem accounts"""

    # Extent account
    extent = {lc: {"area_ha": area, "pct": round(area / study_area.total_area_ha * 100, 1)}
              for lc, area in study_area.land_cover_classes.items()}

    # Condition account (composite index)
    condition_indicators = {
        "vegetation_cover": 0.75,
        "species_richness": 0.65,
        "soil_organic_carbon": 0.80,
        "water_quality": 0.60,
        "landscape_connectivity": service_results["habitat_quality"]["summary"]["mean_habitat_quality"],
    }
    n = len(condition_indicators)
    composite = math.exp(sum(math.log(v) for v in condition_indicators.values()) / n)

    # Service flow account (physical units)
    service_flows = {
        "carbon_storage_tC": service_results["carbon"]["summary"]["total_storage_tC"],
        "erosion_avoided_ton": service_results["sediment"]["summary"]["total_avoided_erosion_ton"],
        "water_yield_m3": service_results["water"]["summary"]["total_water_yield_m3"],
        "N_retained_kg": service_results["nutrient"]["summary"]["total_N_retained_kg"],
        "P_retained_kg": service_results["nutrient"]["summary"]["total_P_retained_kg"],
        "recreation_visits": service_results["recreation"]["summary"]["total_annual_visits"],
        "mean_pollinator_index": service_results["pollination"]["summary"]["mean_pollinator_abundance"],
    }

    # Monetary account
    annual_value = sum([
        service_results["carbon"]["summary"]["monetary_value_social_JPY"] * 0.02,  # annualized
        service_results["sediment"]["summary"]["monetary_value_JPY"],
        service_results["water"]["summary"]["monetary_value_JPY"],
        service_results["nutrient"]["summary"]["monetary_value_total_JPY"],
        service_results["recreation"]["summary"]["monetary_value_JPY"],
        service_results["pollination"]["summary"]["monetary_value_JPY"],
    ])

    npv_results = {}
    for rate_type in ["ddr", "constant_low", "constant_mid", "constant_high"]:
        for horizon in [25, 50, 100]:
            key = f"{rate_type}_{horizon}y"
            npv_results[key] = compute_npv(annual_value, horizon, rate_type)

    return {
        "extent_account": extent,
        "condition_account": {
            "indicators": condition_indicators,
            "composite_condition_index": round(composite, 3),
        },
        "service_flow_account": service_flows,
        "monetary_account": {
            "annual_total_value_JPY": round(annual_value),
            "annual_value_per_ha_JPY": round(annual_value / study_area.total_area_ha),
            "npv_scenarios": npv_results,
        },
    }


# ==============================================================================
# 5. WTP Integration (Choice Experiment Results)
# ==============================================================================

def simulate_wtp_estimates() -> Dict:
    """Simulated WTP estimates from choice experiment (Mixed Logit)"""
    return {
        "model": "Mixed Logit (simulated exchange values)",
        "n_respondents": 600,
        "n_observations": 2400,
        "log_likelihood": -2156.3,
        "pseudo_r2": 0.32,
        "coefficients": {
            "ASC_sq": {"mean": 0.45, "se": 0.12, "p": 0.001, "distribution": "fixed"},
            "biodiversity": {"mean": 0.82, "se": 0.15, "p": 0.001, "sd": 0.45, "distribution": "normal"},
            "water_quality": {"mean": 0.68, "se": 0.13, "p": 0.001, "sd": 0.38, "distribution": "normal"},
            "landscape": {"mean": 0.95, "se": 0.18, "p": 0.001, "sd": 0.52, "distribution": "normal"},
            "carbon": {"mean": 0.35, "se": 0.11, "p": 0.002, "sd": 0.28, "distribution": "normal"},
            "cost": {"mean": -0.00042, "se": 0.00005, "p": 0.001, "distribution": "fixed"},
        },
        "wtp_estimates_jpy_household_year": {
            "biodiversity_20pct": {"mean": 1952, "ci95_low": 1285, "ci95_high": 2619},
            "biodiversity_40pct": {"mean": 3905, "ci95_low": 2571, "ci95_high": 5239},
            "water_quality_standard": {"mean": 1619, "ci95_low": 1023, "ci95_high": 2215},
            "water_quality_swimmable": {"mean": 3238, "ci95_low": 2047, "ci95_high": 4430},
            "landscape_partial": {"mean": 2262, "ci95_low": 1452, "ci95_high": 3071},
            "landscape_traditional": {"mean": 4524, "ci95_low": 2905, "ci95_high": 6143},
            "carbon_increase_30pct": {"mean": 833, "ci95_low": 417, "ci95_high": 1250},
        },
        "aggregate_wtp": {
            "n_households": 15000,
            "total_annual_wtp_JPY": {
                "all_attributes_max": round(15000 * (3905 + 3238 + 4524 + 833)),
                "conservative_mean": round(15000 * (1952 + 1619 + 2262 + 833)),
            }
        },
        "latent_classes": {
            "class_1_conservationist": {"share": 0.35, "highest_wtp": "landscape_traditional"},
            "class_2_pragmatist": {"share": 0.45, "highest_wtp": "water_quality_swimmable"},
            "class_3_indifferent": {"share": 0.20, "highest_wtp": "none (protest responses)"},
        }
    }


# ==============================================================================
# 6. Sensitivity Analysis
# ==============================================================================

def run_sensitivity_analysis(study_area: StudyArea) -> Dict:
    """Run sensitivity analysis on key parameters"""
    results = {}

    # Carbon price sensitivity
    carbon_prices = [3000, 6000, 12000, 20000, 50000]
    carbon_tCO2 = sum(
        CARBON_POOLS[lc].total * area * 44 / 12
        for lc, area in study_area.land_cover_classes.items()
        if lc in CARBON_POOLS
    )
    results["carbon_price_sensitivity"] = {
        str(p): {"price_JPY_tCO2": p, "value_JPY": round(carbon_tCO2 * p)}
        for p in carbon_prices
    }

    # Discount rate sensitivity on 50-year NPV
    base_annual = 500_000_000  # approximate annual flow
    rates = [0.01, 0.02, 0.035, 0.05, 0.07]
    results["discount_rate_sensitivity"] = {}
    for r in rates:
        npv = sum(base_annual / (1 + r)**t for t in range(50))
        results["discount_rate_sensitivity"][str(r)] = {
            "rate": r,
            "npv_50y_JPY": round(npv),
            "ratio_to_baseline": round(npv / sum(base_annual / 1.035**t for t in range(50)), 2),
        }

    # Land use change scenarios
    scenarios = {
        "baseline": dict(study_area.land_cover_classes),
        "abandonment": {
            "mixed_forest": 1200, "paddy_field": 800, "upland_crop": 300,
            "bamboo_grove": 700, "irrigation_pond": 100, "grassland_meadow": 200,
            "settlement": 600, "other": 1100,
        },
        "conservation": {
            "mixed_forest": 2000, "paddy_field": 1000, "upland_crop": 400,
            "bamboo_grove": 200, "irrigation_pond": 300, "grassland_meadow": 600,
            "settlement": 400, "other": 100,
        },
        "urbanization": {
            "mixed_forest": 1000, "paddy_field": 800, "upland_crop": 400,
            "bamboo_grove": 200, "irrigation_pond": 100, "grassland_meadow": 200,
            "settlement": 1800, "other": 500,
        },
    }

    results["land_use_scenarios"] = {}
    for name, lc_map in scenarios.items():
        sa = StudyArea()
        sa.land_cover_classes = lc_map
        carbon = compute_carbon_storage(sa)
        sediment = compute_sediment_delivery(sa)
        water = compute_water_yield(sa)
        nutrient = compute_nutrient_retention(sa)
        results["land_use_scenarios"][name] = {
            "carbon_tC": carbon["summary"]["total_storage_tC"],
            "erosion_avoided_ton": sediment["summary"]["total_avoided_erosion_ton"],
            "water_yield_m3": water["summary"]["total_water_yield_m3"],
            "N_retained_kg": nutrient["summary"]["total_N_retained_kg"],
        }

    return results


# ==============================================================================
# 7. Main Pipeline Execution
# ==============================================================================

def run_pipeline():
    """Execute the full evaluation pipeline"""
    print("=" * 70)
    print("Ecosystem Service Valuation Pipeline - Satoyama Case Study")
    print("=" * 70)

    study_area = StudyArea()

    # Run all InVEST models
    print("\n[1/7] Computing Carbon Storage and Sequestration...")
    carbon = compute_carbon_storage(study_area)

    print("[2/7] Computing Sediment Delivery Ratio...")
    sediment = compute_sediment_delivery(study_area)

    print("[3/7] Computing Annual Water Yield...")
    water = compute_water_yield(study_area)

    print("[4/7] Computing Nutrient Delivery Ratio...")
    nutrient = compute_nutrient_retention(study_area)

    print("[5/7] Computing Habitat Quality...")
    habitat = compute_habitat_quality(study_area)

    print("[6/7] Computing Pollination Services...")
    pollination = compute_pollination(study_area)

    print("[7/7] Computing Recreation Value...")
    recreation = compute_recreation(study_area)

    service_results = {
        "carbon": carbon,
        "sediment": sediment,
        "water": water,
        "nutrient": nutrient,
        "habitat_quality": habitat,
        "pollination": pollination,
        "recreation": recreation,
    }

    # Compile SEEA-EA accounts
    print("\n[SEEA-EA] Compiling ecosystem accounts...")
    seea = compile_seea_accounts(study_area, service_results)

    # WTP estimates
    print("[WTP] Integrating choice experiment results...")
    wtp = simulate_wtp_estimates()

    # Sensitivity analysis
    print("[SA] Running sensitivity analysis...")
    sensitivity = run_sensitivity_analysis(study_area)

    # Compile full results
    full_results = {
        "metadata": {
            "study_area": study_area.name,
            "total_area_ha": study_area.total_area_ha,
            "analysis_date": "2026-05-23",
            "framework_version": "1.0.0",
        },
        "invest_results": service_results,
        "seea_ea_accounts": seea,
        "wtp_estimates": wtp,
        "sensitivity_analysis": sensitivity,
    }

    # Save results
    output_dir = Path("results")
    output_dir.mkdir(exist_ok=True)
    with open(output_dir / "pipeline_results.json", "w", encoding="utf-8") as f:
        json.dump(full_results, f, ensure_ascii=False, indent=2)

    # Print summary
    print("\n" + "=" * 70)
    print("PIPELINE RESULTS SUMMARY")
    print("=" * 70)
    print(f"\nStudy Area: {study_area.name}")
    print(f"Total Area: {study_area.total_area_ha} ha")

    print("\n--- InVEST Model Results ---")
    print(f"Carbon Storage: {carbon['summary']['total_storage_tC']:,.1f} tC "
          f"({carbon['summary']['total_storage_tCO2']:,.1f} tCO2)")
    print(f"  Social Value: ¥{carbon['summary']['monetary_value_social_JPY']:,.0f}")
    print(f"Erosion Avoided: {sediment['summary']['total_avoided_erosion_ton']:,.1f} ton/yr "
          f"(retention {sediment['summary']['overall_retention_pct']:.1f}%)")
    print(f"  Value: ¥{sediment['summary']['monetary_value_JPY']:,.0f}")
    print(f"Water Yield: {water['summary']['total_water_yield_m3']:,.0f} m³/yr")
    print(f"  Value: ¥{water['summary']['monetary_value_JPY']:,.0f}")
    print(f"N Retained: {nutrient['summary']['total_N_retained_kg']:,.1f} kg/yr "
          f"({nutrient['summary']['N_retention_pct']:.1f}%)")
    print(f"P Retained: {nutrient['summary']['total_P_retained_kg']:,.1f} kg/yr "
          f"({nutrient['summary']['P_retention_pct']:.1f}%)")
    print(f"  Value: ¥{nutrient['summary']['monetary_value_total_JPY']:,.0f}")
    print(f"Habitat Quality: {habitat['summary']['mean_habitat_quality']:.3f} "
          f"(coverage {habitat['summary']['habitat_coverage_pct']:.1f}%)")
    print(f"Pollinator Abundance: {pollination['summary']['mean_pollinator_abundance']:.3f}")
    print(f"  Value: ¥{pollination['summary']['monetary_value_JPY']:,.0f}")
    print(f"Recreation: {recreation['summary']['total_annual_visits']:,.0f} person-days/yr")
    print(f"  Value: ¥{recreation['summary']['monetary_value_JPY']:,.0f}")

    print(f"\n--- SEEA-EA Monetary Account ---")
    print(f"Annual Total Value: ¥{seea['monetary_account']['annual_total_value_JPY']:,.0f}")
    print(f"Annual Value per ha: ¥{seea['monetary_account']['annual_value_per_ha_JPY']:,.0f}/ha")
    print(f"Condition Index: {seea['condition_account']['composite_condition_index']:.3f}")

    npv_25 = seea['monetary_account']['npv_scenarios']['ddr_25y']['npv_JPY']
    npv_50 = seea['monetary_account']['npv_scenarios']['ddr_50y']['npv_JPY']
    npv_100 = seea['monetary_account']['npv_scenarios']['ddr_100y']['npv_JPY']
    print(f"\nNPV (DDR): 25yr=¥{npv_25:,.0f} / 50yr=¥{npv_50:,.0f} / 100yr=¥{npv_100:,.0f}")

    print(f"\n--- WTP Estimates (Choice Experiment) ---")
    for attr, vals in wtp['wtp_estimates_jpy_household_year'].items():
        print(f"  {attr}: ¥{vals['mean']:,.0f}/世帯/年 "
              f"[95%CI: ¥{vals['ci95_low']:,.0f}-¥{vals['ci95_high']:,.0f}]")

    print(f"\nAggregate WTP (conservative): "
          f"¥{wtp['aggregate_wtp']['total_annual_wtp_JPY']['conservative_mean']:,.0f}/yr")

    print(f"\n--- Land Use Scenario Comparison ---")
    for scenario, vals in sensitivity['land_use_scenarios'].items():
        print(f"  {scenario}: Carbon={vals['carbon_tC']:,.0f}tC, "
              f"Erosion avoided={vals['erosion_avoided_ton']:,.0f}t, "
              f"Water={vals['water_yield_m3']:,.0f}m³")

    print(f"\nResults saved to: results/pipeline_results.json")
    print("=" * 70)
    return full_results


if __name__ == "__main__":
    run_pipeline()
