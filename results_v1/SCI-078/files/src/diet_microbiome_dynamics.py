"""
Diet-Microbiome Long-term Dynamics Simulator
==============================================
Simulates the impact of different dietary patterns on gut microbiota
composition over weeks to months.

Diet patterns modeled:
  1. Western diet (high fat, low fiber)
  2. Mediterranean diet (balanced, high fiber)
  3. Plant-based/vegan diet (very high fiber)
  4. High-protein diet (ketogenic-like)
"""

import numpy as np
from typing import Dict, List, Callable
from src.glv_community_model import (
    GLVParameters, build_interaction_matrix, run_glv_simulation,
    compute_steady_state_composition, SPECIES_NAMES, SPECIES_SHORT
)


DIET_PATTERNS = {
    'western': {
        'name': 'Western Diet',
        'description': 'High fat, high sugar, low fiber',
        # [fiber, starch, protein, simple_sugars, mucin]
        'substrate_supply': np.array([0.15, 0.4, 0.3, 0.6, 0.2]),
        'color': '#e74c3c',
    },
    'mediterranean': {
        'name': 'Mediterranean Diet',
        'description': 'Balanced, high fiber, moderate protein',
        'substrate_supply': np.array([0.8, 0.5, 0.25, 0.2, 0.2]),
        'color': '#2ecc71',
    },
    'plant_based': {
        'name': 'Plant-based Diet',
        'description': 'Very high fiber, low protein',
        'substrate_supply': np.array([1.2, 0.6, 0.1, 0.15, 0.15]),
        'color': '#27ae60',
    },
    'high_protein': {
        'name': 'High-protein Diet',
        'description': 'High protein, very low carb',
        'substrate_supply': np.array([0.1, 0.1, 0.8, 0.05, 0.2]),
        'color': '#3498db',
    },
}


def create_diet_substrate_func(diet_key: str, meal_frequency: int = 3) -> Callable:
    """Create time-dependent substrate supply function for a diet pattern."""
    diet = DIET_PATTERNS[diet_key]
    base_supply = diet['substrate_supply']

    def substrate_func(t):
        # Simulate meal pulses (3 meals/day = every 8 hours)
        hour_of_day = t % 24
        meal_times = np.linspace(7, 19, meal_frequency)  # 7am, 1pm, 7pm
        meal_pulse = sum(
            np.exp(-0.5 * ((hour_of_day - mt) / 1.0) ** 2)
            for mt in meal_times
        )
        # Basal supply + meal pulses
        return base_supply * (0.3 + 0.7 * meal_pulse / meal_frequency)

    return substrate_func


def simulate_diet_pattern(
    diet_key: str,
    duration_days: int = 30,
    initial_abundances: np.ndarray = None,
    params: GLVParameters = None,
) -> dict:
    """Simulate microbiota response to a specific diet pattern."""
    t_span = (0, duration_days * 24)
    t_points = duration_days * 24  # hourly resolution

    substrate_func = create_diet_substrate_func(diet_key)

    results = run_glv_simulation(
        initial_abundances=initial_abundances,
        params=params,
        substrate_func=substrate_func,
        t_span=t_span,
        t_points=t_points,
    )

    results['diet'] = diet_key
    results['diet_info'] = DIET_PATTERNS[diet_key]
    results['duration_days'] = duration_days

    return results


def simulate_diet_comparison(
    duration_days: int = 30,
    initial_abundances: np.ndarray = None,
) -> Dict[str, dict]:
    """Compare microbiota response across all diet patterns."""
    if initial_abundances is None:
        initial_abundances = np.array([
            100, 80, 60, 70, 30, 20, 40, 5, 50, 45
        ], dtype=float)

    results = {}
    for diet_key in DIET_PATTERNS:
        results[diet_key] = simulate_diet_pattern(
            diet_key,
            duration_days=duration_days,
            initial_abundances=initial_abundances.copy(),
        )

    return results


def simulate_diet_switch(
    diet_sequence: List[str],
    phase_durations: List[int],
    initial_abundances: np.ndarray = None,
) -> dict:
    """
    Simulate sequential diet changes.
    E.g., Western → Mediterranean transition over multiple phases.
    """
    if initial_abundances is None:
        initial_abundances = np.array([
            100, 80, 60, 70, 30, 20, 40, 5, 50, 45
        ], dtype=float)

    all_time = []
    all_abundances = []
    all_substrates = []
    all_shannon = []
    all_simpson = []
    phase_boundaries = [0]

    current_abundances = initial_abundances.copy()
    current_substrates = np.array([5.0, 3.0, 2.0, 1.0, 2.0])
    time_offset = 0

    for diet_key, duration in zip(diet_sequence, phase_durations):
        results = simulate_diet_pattern(
            diet_key,
            duration_days=duration,
            initial_abundances=current_abundances,
        )

        # Append results with time offset
        all_time.append(results['time'] + time_offset)
        all_abundances.append(results['abundances'])
        all_substrates.append(results['substrates'])
        all_shannon.append(results['shannon_diversity'])
        all_simpson.append(results['simpson_diversity'])

        # Update for next phase
        current_abundances = results['abundances'][:, -1]
        current_substrates = results['substrates'][:, -1]
        time_offset += duration * 24
        phase_boundaries.append(time_offset)

    return {
        'time': np.concatenate(all_time),
        'abundances': np.concatenate(all_abundances, axis=1),
        'substrates': np.concatenate(all_substrates, axis=1),
        'shannon_diversity': np.concatenate(all_shannon),
        'simpson_diversity': np.concatenate(all_simpson),
        'diet_sequence': diet_sequence,
        'phase_durations': phase_durations,
        'phase_boundaries': phase_boundaries,
        'species_names': SPECIES_NAMES,
        'species_short': SPECIES_SHORT,
    }


def compute_diet_impact_metrics(comparison_results: Dict[str, dict]) -> dict:
    """Compute summary metrics comparing diet impacts."""
    metrics = {}
    for diet_key, results in comparison_results.items():
        ss = compute_steady_state_composition(results)
        metrics[diet_key] = {
            'diet_name': DIET_PATTERNS[diet_key]['name'],
            'shannon_diversity': float(ss['shannon_diversity']),
            'simpson_diversity': float(ss['simpson_diversity']),
            'dominant_species': ss['dominant_species'],
            'butyrate_producers_fraction': float(
                ss['relative_abundance'][1] + ss['relative_abundance'][2]  # F.praus + R.intest
            ),
            'pathobiont_fraction': float(
                ss['relative_abundance'][5] + ss['relative_abundance'][7]  # E.coli + C.diff
            ),
            'firmicutes_bacteroidetes_ratio': float(
                (ss['relative_abundance'][1] + ss['relative_abundance'][2] +
                 ss['relative_abundance'][6] + ss['relative_abundance'][9]) /
                (ss['relative_abundance'][0] + ss['relative_abundance'][8] + 1e-12)
            ),
        }
    return metrics


if __name__ == "__main__":
    print("Running diet comparison simulation...")
    comparison = simulate_diet_comparison(duration_days=14)
    metrics = compute_diet_impact_metrics(comparison)

    for diet_key, m in metrics.items():
        print(f"\n{m['diet_name']}:")
        print(f"  Shannon diversity: {m['shannon_diversity']:.3f}")
        print(f"  Dominant species: {m['dominant_species']}")
        print(f"  Butyrate producers: {m['butyrate_producers_fraction']:.1%}")
        print(f"  Pathobionts: {m['pathobiont_fraction']:.1%}")
