"""
Probiotic and Prebiotic Effect Prediction Model
=================================================
Models the impact of probiotic supplementation and prebiotic
substrates on gut microbial community dynamics.

Probiotics: Direct addition of specific microbial species
Prebiotics: Selective substrates that promote beneficial microbes
Synbiotics: Combination of probiotics + prebiotics
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from src.glv_community_model import (
    GLVParameters, build_interaction_matrix, run_glv_simulation,
    compute_steady_state_composition, SPECIES_NAMES, SPECIES_SHORT
)
from src.scfa_flux_model import compute_scfa_production_rates, compute_scfa_accumulation


# Probiotic strains and their mapping to model species
PROBIOTIC_STRAINS = {
    'lactobacillus_rhamnosus_gg': {
        'species_idx': 6,  # L. rhamnosus
        'dose_cfu': 1e10,
        'model_abundance_per_dose': 50.0,
        'survival_rate': 0.01,  # fraction surviving gastric passage
        'colonization_half_life': 48,  # hours
        'description': 'LGG - immune modulation, pathogen exclusion',
    },
    'bifidobacterium_longum_bb536': {
        'species_idx': 3,  # B. longum
        'dose_cfu': 5e9,
        'model_abundance_per_dose': 30.0,
        'survival_rate': 0.05,
        'colonization_half_life': 72,
        'description': 'BB536 - allergy reduction, fiber fermentation',
    },
    'akkermansia_muciniphila': {
        'species_idx': 4,  # A. muciniphila
        'dose_cfu': 1e9,
        'model_abundance_per_dose': 15.0,
        'survival_rate': 0.001,  # strict anaerobe, low survival
        'colonization_half_life': 24,
        'description': 'Mucin-degrading, metabolic health',
    },
}

# Prebiotic substrates
PREBIOTIC_SUBSTRATES = {
    'inulin': {
        'description': 'Fructooligosaccharide from chicory root',
        'dose_g_per_day': 10.0,
        # Selective growth enhancement factors per species
        'growth_boost': np.array([
            1.1, 1.3, 1.2, 1.8, 1.0, 0.9, 1.1, 0.8, 1.1, 1.1
        ]),
        # Additional substrate supply [fiber, starch, protein, sugars, mucin]
        'substrate_addition': np.array([0.4, 0.0, 0.0, 0.1, 0.0]),
    },
    'galactooligosaccharides': {
        'description': 'GOS from lactose conversion',
        'dose_g_per_day': 5.0,
        'growth_boost': np.array([
            1.0, 1.1, 1.1, 2.0, 1.0, 0.9, 1.3, 0.8, 1.0, 1.0
        ]),
        'substrate_addition': np.array([0.2, 0.0, 0.0, 0.2, 0.0]),
    },
    'resistant_starch': {
        'description': 'Type 2/3 resistant starch',
        'dose_g_per_day': 15.0,
        'growth_boost': np.array([
            1.2, 1.4, 1.3, 1.1, 1.0, 0.9, 1.0, 0.8, 1.1, 2.0
        ]),
        'substrate_addition': np.array([0.1, 0.6, 0.0, 0.0, 0.0]),
    },
    'pectin': {
        'description': 'Apple/citrus pectin',
        'dose_g_per_day': 8.0,
        'growth_boost': np.array([
            1.5, 1.2, 1.1, 1.3, 1.0, 0.9, 1.0, 0.8, 1.3, 1.0
        ]),
        'substrate_addition': np.array([0.3, 0.0, 0.0, 0.1, 0.0]),
    },
}


def simulate_probiotic_intervention(
    probiotic_key: str,
    duration_days: int = 28,
    dosing_frequency_h: float = 24.0,
    baseline_abundances: np.ndarray = None,
    diet_substrate_func=None,
) -> dict:
    """Simulate the effect of probiotic supplementation."""
    probiotic = PROBIOTIC_STRAINS[probiotic_key]
    species_idx = probiotic['species_idx']
    effective_dose = (
        probiotic['model_abundance_per_dose'] *
        probiotic['survival_rate']
    )

    if baseline_abundances is None:
        baseline_abundances = np.array([
            100, 80, 60, 70, 30, 20, 40, 5, 50, 45
        ], dtype=float)

    # Modify initial abundances with probiotic addition
    probiotic_abundances = baseline_abundances.copy()
    probiotic_abundances[species_idx] += effective_dose

    # Create substrate function with probiotic's colonization dynamics
    half_life = probiotic['colonization_half_life']
    decay_rate = np.log(2) / half_life

    def probiotic_substrate_func(t):
        base = np.array([0.5, 0.3, 0.2, 0.1, 0.2])
        if diet_substrate_func is not None:
            base = diet_substrate_func(t)
        # Periodic dosing adds transient boost
        dose_effect = effective_dose * np.exp(-decay_rate * (t % dosing_frequency_h))
        return base

    # Run baseline (no probiotic)
    baseline_results = run_glv_simulation(
        initial_abundances=baseline_abundances.copy(),
        substrate_func=diet_substrate_func,
        t_span=(0, duration_days * 24),
        t_points=duration_days * 12,
    )

    # Run with probiotic
    probiotic_results = run_glv_simulation(
        initial_abundances=probiotic_abundances,
        substrate_func=probiotic_substrate_func,
        t_span=(0, duration_days * 24),
        t_points=duration_days * 12,
    )

    # Compute SCFA changes
    baseline_scfa = compute_scfa_production_rates(
        baseline_results['abundances'], baseline_results['substrates']
    )
    probiotic_scfa = compute_scfa_production_rates(
        probiotic_results['abundances'], probiotic_results['substrates']
    )

    return {
        'probiotic': probiotic_key,
        'probiotic_info': probiotic,
        'baseline': baseline_results,
        'intervention': probiotic_results,
        'baseline_scfa': baseline_scfa,
        'intervention_scfa': probiotic_scfa,
        'duration_days': duration_days,
    }


def simulate_prebiotic_intervention(
    prebiotic_key: str,
    duration_days: int = 28,
    baseline_abundances: np.ndarray = None,
) -> dict:
    """Simulate the effect of prebiotic supplementation."""
    prebiotic = PREBIOTIC_SUBSTRATES[prebiotic_key]

    if baseline_abundances is None:
        baseline_abundances = np.array([
            100, 80, 60, 70, 30, 20, 40, 5, 50, 45
        ], dtype=float)

    # Modify growth parameters with prebiotic boost
    params_baseline = GLVParameters()
    params_prebiotic = GLVParameters()
    params_prebiotic.mu_max = params_baseline.mu_max * prebiotic['growth_boost']

    # Create substrate function with prebiotic addition
    substrate_addition = prebiotic['substrate_addition'] / 24.0  # per hour

    def prebiotic_substrate_func(t):
        base = np.array([0.5, 0.3, 0.2, 0.1, 0.2])
        return base + substrate_addition

    # Run baseline
    baseline_results = run_glv_simulation(
        initial_abundances=baseline_abundances.copy(),
        params=params_baseline,
        t_span=(0, duration_days * 24),
        t_points=duration_days * 12,
    )

    # Run with prebiotic
    prebiotic_results = run_glv_simulation(
        initial_abundances=baseline_abundances.copy(),
        params=params_prebiotic,
        substrate_func=prebiotic_substrate_func,
        t_span=(0, duration_days * 24),
        t_points=duration_days * 12,
    )

    return {
        'prebiotic': prebiotic_key,
        'prebiotic_info': prebiotic,
        'baseline': baseline_results,
        'intervention': prebiotic_results,
        'duration_days': duration_days,
    }


def simulate_synbiotic(
    probiotic_key: str,
    prebiotic_key: str,
    duration_days: int = 28,
    baseline_abundances: np.ndarray = None,
) -> dict:
    """Simulate synbiotic (probiotic + prebiotic) intervention."""
    probiotic = PROBIOTIC_STRAINS[probiotic_key]
    prebiotic = PREBIOTIC_SUBSTRATES[prebiotic_key]

    if baseline_abundances is None:
        baseline_abundances = np.array([
            100, 80, 60, 70, 30, 20, 40, 5, 50, 45
        ], dtype=float)

    # Combine probiotic dose + prebiotic growth boost
    synbiotic_abundances = baseline_abundances.copy()
    synbiotic_abundances[probiotic['species_idx']] += (
        probiotic['model_abundance_per_dose'] * probiotic['survival_rate']
    )

    params = GLVParameters()
    params.mu_max = params.mu_max * prebiotic['growth_boost']

    substrate_addition = prebiotic['substrate_addition'] / 24.0

    def synbiotic_substrate_func(t):
        base = np.array([0.5, 0.3, 0.2, 0.1, 0.2])
        return base + substrate_addition

    # Run baseline
    baseline_results = run_glv_simulation(
        initial_abundances=baseline_abundances.copy(),
        t_span=(0, duration_days * 24),
        t_points=duration_days * 12,
    )

    # Run synbiotic
    synbiotic_results = run_glv_simulation(
        initial_abundances=synbiotic_abundances,
        params=params,
        substrate_func=synbiotic_substrate_func,
        t_span=(0, duration_days * 24),
        t_points=duration_days * 12,
    )

    return {
        'probiotic': probiotic_key,
        'prebiotic': prebiotic_key,
        'baseline': baseline_results,
        'intervention': synbiotic_results,
        'duration_days': duration_days,
    }


def compare_interventions(duration_days: int = 28) -> dict:
    """Compare all intervention types."""
    baseline_abundances = np.array([
        100, 80, 60, 70, 30, 20, 40, 5, 50, 45
    ], dtype=float)

    results = {
        'probiotic': {},
        'prebiotic': {},
        'synbiotic': None,
    }

    for pk in PROBIOTIC_STRAINS:
        results['probiotic'][pk] = simulate_probiotic_intervention(
            pk, duration_days=duration_days,
            baseline_abundances=baseline_abundances.copy()
        )

    for pk in PREBIOTIC_SUBSTRATES:
        results['prebiotic'][pk] = simulate_prebiotic_intervention(
            pk, duration_days=duration_days,
            baseline_abundances=baseline_abundances.copy()
        )

    # Example synbiotic
    results['synbiotic'] = simulate_synbiotic(
        'bifidobacterium_longum_bb536', 'inulin',
        duration_days=duration_days,
        baseline_abundances=baseline_abundances.copy()
    )

    return results


if __name__ == "__main__":
    print("Running probiotic intervention (LGG)...")
    prob_results = simulate_probiotic_intervention('lactobacillus_rhamnosus_gg')
    baseline_ss = compute_steady_state_composition(prob_results['baseline'])
    interv_ss = compute_steady_state_composition(prob_results['intervention'])

    print(f"\nBaseline Shannon: {baseline_ss['shannon_diversity']:.3f}")
    print(f"Probiotic Shannon: {interv_ss['shannon_diversity']:.3f}")
    print(f"Change: {interv_ss['shannon_diversity'] - baseline_ss['shannon_diversity']:+.3f}")
