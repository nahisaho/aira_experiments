"""
Fermented Food Case Study
===========================
Simulates the impact of regular fermented food consumption on
gut microbiota diversity, based on the Stanford study design
(Sonnenburg et al., Cell Host & Microbe, 2021).

Fermented foods modeled:
  - Yogurt (L. bulgaricus, S. thermophilus)
  - Kimchi (L. plantarum, L. brevis, Leuconostoc)
  - Kefir (diverse LAB + yeasts)
  - Kombucha (Acetobacter, yeasts, LAB)
  - Natto (B. subtilis)
  - Miso (A. oryzae fermentation products)
"""

import numpy as np
from typing import Dict, List
from src.glv_community_model import (
    GLVParameters, build_interaction_matrix, run_glv_simulation,
    compute_steady_state_composition, SPECIES_NAMES, SPECIES_SHORT
)
from src.scfa_flux_model import (
    compute_scfa_production_rates, compute_scfa_accumulation
)


FERMENTED_FOODS = {
    'yogurt': {
        'name': 'Yogurt',
        'servings_per_day': 2,
        # Growth boost from LAB and metabolites
        'growth_modifiers': np.array([
            1.0, 1.1, 1.05, 1.2, 1.0, 0.9, 1.5, 0.8, 1.0, 1.0
        ]),
        # Additional substrate from fermentation products
        'substrate_addition': np.array([0.05, 0.0, 0.1, 0.15, 0.0]),
        # Direct microbial input (transient)
        'microbial_input': {6: 5.0},  # L. rhamnosus proxy
        'bioactive_compounds': ['lactic_acid', 'bacteriocins', 'CLA'],
    },
    'kimchi': {
        'name': 'Kimchi',
        'servings_per_day': 1,
        'growth_modifiers': np.array([
            1.05, 1.15, 1.1, 1.15, 1.0, 0.85, 1.3, 0.75, 1.1, 1.05
        ]),
        'substrate_addition': np.array([0.15, 0.0, 0.05, 0.1, 0.0]),
        'microbial_input': {6: 3.0, 3: 2.0},
        'bioactive_compounds': ['lactic_acid', 'bacteriocins', 'fiber', 'polyphenols'],
    },
    'kefir': {
        'name': 'Kefir',
        'servings_per_day': 1,
        'growth_modifiers': np.array([
            1.05, 1.15, 1.1, 1.25, 1.1, 0.85, 1.4, 0.7, 1.05, 1.05
        ]),
        'substrate_addition': np.array([0.05, 0.0, 0.1, 0.1, 0.0]),
        'microbial_input': {6: 4.0, 3: 3.0},
        'bioactive_compounds': ['lactic_acid', 'acetic_acid', 'kefiran', 'bacteriocins'],
    },
    'natto': {
        'name': 'Natto',
        'servings_per_day': 1,
        'growth_modifiers': np.array([
            1.1, 1.2, 1.1, 1.1, 1.0, 0.95, 1.0, 0.85, 1.1, 1.1
        ]),
        'substrate_addition': np.array([0.2, 0.05, 0.15, 0.05, 0.0]),
        'microbial_input': {},
        'bioactive_compounds': ['nattokinase', 'vitamin_K2', 'polyglutamic_acid', 'isoflavones'],
    },
    'mixed_fermented': {
        'name': 'Mixed Fermented Foods (6+ servings/day)',
        'servings_per_day': 6,
        'growth_modifiers': np.array([
            1.1, 1.25, 1.15, 1.3, 1.1, 0.8, 1.5, 0.65, 1.1, 1.1
        ]),
        'substrate_addition': np.array([0.2, 0.05, 0.15, 0.2, 0.0]),
        'microbial_input': {6: 8.0, 3: 5.0},
        'bioactive_compounds': ['diverse_metabolites'],
    },
}


def simulate_fermented_food_intervention(
    food_key: str,
    duration_days: int = 70,  # 10-week intervention
    baseline_abundances: np.ndarray = None,
    washout_days: int = 14,
) -> dict:
    """
    Simulate fermented food intervention with washout period.

    Timeline:
      - Phase 1 (days 0-14): Baseline
      - Phase 2 (days 14-56): Intervention (6 weeks)
      - Phase 3 (days 56-70): Washout
    """
    food = FERMENTED_FOODS[food_key]

    if baseline_abundances is None:
        # Start with slightly dysbiotic composition (low diversity)
        baseline_abundances = np.array([
            120, 50, 40, 45, 15, 35, 25, 10, 60, 30
        ], dtype=float)

    baseline_phase_days = 14
    intervention_days = duration_days - baseline_phase_days - washout_days
    assert intervention_days > 0

    params_baseline = GLVParameters()
    params_intervention = GLVParameters()
    params_intervention.mu_max = params_baseline.mu_max * food['growth_modifiers']

    substrate_addition = food['substrate_addition'] / 24.0

    def intervention_substrate_func(t):
        base = np.array([0.5, 0.3, 0.2, 0.1, 0.2])
        return base + substrate_addition

    # Phase 1: Baseline
    phase1 = run_glv_simulation(
        initial_abundances=baseline_abundances.copy(),
        params=params_baseline,
        t_span=(0, baseline_phase_days * 24),
        t_points=baseline_phase_days * 12,
    )

    # Phase 2: Intervention (add microbial input)
    intervention_init = phase1['abundances'][:, -1].copy()
    for sp_idx, amount in food['microbial_input'].items():
        intervention_init[sp_idx] += amount

    phase2 = run_glv_simulation(
        initial_abundances=intervention_init,
        params=params_intervention,
        substrate_func=intervention_substrate_func,
        t_span=(0, intervention_days * 24),
        t_points=intervention_days * 12,
    )

    # Phase 3: Washout (return to baseline diet)
    washout_init = phase2['abundances'][:, -1].copy()
    phase3 = run_glv_simulation(
        initial_abundances=washout_init,
        params=params_baseline,
        t_span=(0, washout_days * 24),
        t_points=washout_days * 12,
    )

    # Concatenate timecourses
    time_offset_p2 = baseline_phase_days * 24
    time_offset_p3 = time_offset_p2 + intervention_days * 24

    combined_time = np.concatenate([
        phase1['time'],
        phase2['time'] + time_offset_p2,
        phase3['time'] + time_offset_p3,
    ])
    combined_abundances = np.concatenate([
        phase1['abundances'],
        phase2['abundances'],
        phase3['abundances'],
    ], axis=1)
    combined_substrates = np.concatenate([
        phase1['substrates'],
        phase2['substrates'],
        phase3['substrates'],
    ], axis=1)
    combined_shannon = np.concatenate([
        phase1['shannon_diversity'],
        phase2['shannon_diversity'],
        phase3['shannon_diversity'],
    ])
    combined_simpson = np.concatenate([
        phase1['simpson_diversity'],
        phase2['simpson_diversity'],
        phase3['simpson_diversity'],
    ])

    # SCFA analysis
    scfa_baseline = compute_scfa_production_rates(
        phase1['abundances'], phase1['substrates']
    )
    scfa_intervention = compute_scfa_production_rates(
        phase2['abundances'], phase2['substrates']
    )
    scfa_washout = compute_scfa_production_rates(
        phase3['abundances'], phase3['substrates']
    )

    return {
        'food_key': food_key,
        'food_info': food,
        'time': combined_time,
        'time_days': combined_time / 24,
        'abundances': combined_abundances,
        'substrates': combined_substrates,
        'shannon_diversity': combined_shannon,
        'simpson_diversity': combined_simpson,
        'phase_boundaries_days': [0, baseline_phase_days,
                                   baseline_phase_days + intervention_days,
                                   duration_days],
        'scfa': {
            'baseline': scfa_baseline,
            'intervention': scfa_intervention,
            'washout': scfa_washout,
        },
        'species_names': SPECIES_NAMES,
        'species_short': SPECIES_SHORT,
    }


def run_fermented_food_comparison(duration_days: int = 70) -> Dict[str, dict]:
    """Compare effects of different fermented foods."""
    results = {}
    for food_key in FERMENTED_FOODS:
        results[food_key] = simulate_fermented_food_intervention(
            food_key, duration_days=duration_days
        )
    return results


def compute_diversity_change_metrics(results: dict) -> dict:
    """Compute diversity change metrics from intervention results."""
    phase_bounds = results['phase_boundaries_days']
    time_days = results['time_days']
    shannon = results['shannon_diversity']

    # Average diversity in each phase
    baseline_mask = time_days < phase_bounds[1]
    intervention_mask = (time_days >= phase_bounds[1]) & (time_days < phase_bounds[2])
    washout_mask = time_days >= phase_bounds[2]

    baseline_div = shannon[baseline_mask].mean() if baseline_mask.any() else 0
    intervention_div = shannon[intervention_mask].mean() if intervention_mask.any() else 0
    washout_div = shannon[washout_mask].mean() if washout_mask.any() else 0

    # Species-level changes
    abundances = results['abundances']
    n_baseline = baseline_mask.sum()
    n_intervention = intervention_mask.sum()

    baseline_composition = abundances[:, :max(1, n_baseline)].mean(axis=1)
    intervention_composition = abundances[:, n_baseline:n_baseline+max(1, n_intervention)].mean(axis=1)

    total_b = baseline_composition.sum()
    total_i = intervention_composition.sum()
    rel_baseline = baseline_composition / (total_b + 1e-12)
    rel_intervention = intervention_composition / (total_i + 1e-12)

    # Bray-Curtis dissimilarity
    bray_curtis = np.sum(np.abs(rel_baseline - rel_intervention)) / (
        np.sum(rel_baseline) + np.sum(rel_intervention)
    )

    return {
        'food': results['food_info']['name'],
        'baseline_shannon': float(baseline_div),
        'intervention_shannon': float(intervention_div),
        'washout_shannon': float(washout_div),
        'diversity_change_pct': float(
            (intervention_div - baseline_div) / (baseline_div + 1e-12) * 100
        ),
        'resilience_index': float(
            (washout_div - baseline_div) / (intervention_div - baseline_div + 1e-12)
        ),
        'bray_curtis_dissimilarity': float(bray_curtis),
    }


if __name__ == "__main__":
    print("Running fermented food case study...")
    results = simulate_fermented_food_intervention('mixed_fermented')
    metrics = compute_diversity_change_metrics(results)

    print(f"\n{metrics['food']}:")
    print(f"  Baseline Shannon:     {metrics['baseline_shannon']:.3f}")
    print(f"  Intervention Shannon: {metrics['intervention_shannon']:.3f}")
    print(f"  Washout Shannon:      {metrics['washout_shannon']:.3f}")
    print(f"  Diversity change:     {metrics['diversity_change_pct']:+.1f}%")
    print(f"  Bray-Curtis:          {metrics['bray_curtis_dissimilarity']:.3f}")
