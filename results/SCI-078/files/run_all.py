#!/usr/bin/env python3
"""
Main Runner: Gut Microbiota-Diet Interaction Systems Biology Framework
=======================================================================
Runs all simulations and generates figures and results.
"""

import sys
import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from datetime import datetime

# Set style
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'figure.facecolor': 'white',
})
sns.set_palette('colorblind')

from src.shime_digestion_model import (
    SHIMEParameters, FoodComposition, run_shime_simulation,
    compute_absorption_efficiency
)
from src.glv_community_model import (
    GLVParameters, build_interaction_matrix, run_glv_simulation,
    compute_steady_state_composition, SPECIES_NAMES, SPECIES_SHORT
)
from src.scfa_flux_model import (
    SCFAParameters, compute_scfa_production_rates,
    compute_scfa_accumulation, compute_scfa_timecourse
)
from src.diet_microbiome_dynamics import (
    DIET_PATTERNS, simulate_diet_comparison, simulate_diet_switch,
    compute_diet_impact_metrics
)
from src.probiotic_prebiotic_model import (
    PROBIOTIC_STRAINS, PREBIOTIC_SUBSTRATES,
    simulate_probiotic_intervention, simulate_prebiotic_intervention,
    simulate_synbiotic, compare_interventions
)
from src.fermented_food_casestudy import (
    FERMENTED_FOODS, simulate_fermented_food_intervention,
    run_fermented_food_comparison, compute_diversity_change_metrics
)
from src.micom_community_model import (
    run_community_metabolic_analysis, generate_gapseq_report
)


def ensure_dirs():
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)
    os.makedirs('data', exist_ok=True)
    os.makedirs('logs', exist_ok=True)


# ── Figure 1: SHIME Digestion Model ──────────────────────────────────────

def generate_shime_figures():
    print("[1/6] Running SHIME digestion simulation...")

    food = FoodComposition(
        name="high_fiber_meal",
        protein=25.0, starch=40.0, simple_sugars=8.0,
        dietary_fiber=20.0, soluble_fiber=8.0, insoluble_fiber=12.0,
        lipids=15.0, polyphenols=1.0, resistant_starch=8.0
    )
    results = run_shime_simulation(food, t_span=(0, 72), t_points=600)
    efficiency = compute_absorption_efficiency(results)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('SHIME Multi-Compartment Digestion Model', fontsize=14, fontweight='bold')

    compartment_names = ['Stomach', 'Small Intestine', 'Ascending Colon',
                         'Transverse Colon', 'Descending Colon']
    nutrient_names = ['Protein', 'Starch', 'Lipid', 'Fiber', 'Polyphenols']
    colors = ['#e74c3c', '#3498db', '#f39c12', '#2ecc71', '#9b59b6']

    # Panel A-E: Nutrient transit through compartments
    for comp_idx in range(5):
        ax = axes.flat[comp_idx]
        for nut_idx, (nut_name, color) in enumerate(zip(nutrient_names, colors)):
            y_idx = comp_idx * 5 + nut_idx
            ax.plot(results['time'], results['solution'][y_idx],
                    label=nut_name, color=color, linewidth=1.5)
        ax.set_title(compartment_names[comp_idx])
        ax.set_xlabel('Time (h)')
        ax.set_ylabel('Concentration (g)')
        ax.legend(loc='upper right', fontsize=7)
        ax.set_xlim(0, 72)

    # Panel F: Absorption efficiency bar chart
    ax = axes.flat[5]
    nutrients = list(efficiency['absorption_efficiency'].keys())
    values = [efficiency['absorption_efficiency'][n] * 100 for n in nutrients]
    bars = ax.bar(range(len(nutrients)), values, color=colors)
    ax.set_xticks(range(len(nutrients)))
    ax.set_xticklabels([n.capitalize() for n in nutrients], rotation=30)
    ax.set_ylabel('Absorption Efficiency (%)')
    ax.set_title('Nutrient Absorption Efficiency (72h)')
    ax.set_ylim(0, 105)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig('figures/fig1_shime_digestion.png', bbox_inches='tight')
    plt.savefig('figures/fig1_shime_digestion.svg', bbox_inches='tight')
    plt.close()

    # Save results
    with open('results/shime_absorption_efficiency.json', 'w') as f:
        json.dump({k: float(v) for k, v in efficiency['absorption_efficiency'].items()}, f, indent=2)

    print(f"  → Absorption efficiency: {efficiency['absorption_efficiency']}")
    return results, efficiency


# ── Figure 2: gLV Community Dynamics ─────────────────────────────────────

def generate_glv_figures():
    print("[2/6] Running gLV community simulation...")

    results = run_glv_simulation(t_span=(0, 720), t_points=1200)
    ss = compute_steady_state_composition(results)

    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)
    fig.suptitle('Generalized Lotka-Volterra (gLV) Community Dynamics', fontsize=14, fontweight='bold')

    # Panel A: Species abundance over time
    ax1 = fig.add_subplot(gs[0, 0:2])
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    for i in range(10):
        ax1.plot(results['time'] / 24, results['abundances'][i],
                 label=SPECIES_SHORT[i], color=colors[i], linewidth=1.5)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Abundance (a.u.)')
    ax1.set_title('(A) Species Abundance Dynamics')
    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)

    # Panel B: Relative abundance stacked area
    ax2 = fig.add_subplot(gs[1, 0:2])
    time_days = results['time'] / 24
    ax2.stackplot(time_days, results['relative_abundances'],
                  labels=SPECIES_SHORT, colors=colors, alpha=0.8)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Relative Abundance')
    ax2.set_title('(B) Community Composition Over Time')
    ax2.set_ylim(0, 1)
    ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)

    # Panel C: Diversity indices
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(time_days, results['shannon_diversity'], 'b-', label='Shannon', linewidth=2)
    ax3_twin = ax3.twinx()
    ax3_twin.plot(time_days, results['simpson_diversity'], 'r--', label='Simpson', linewidth=2)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Shannon Index', color='b')
    ax3_twin.set_ylabel('Simpson Index', color='r')
    ax3.set_title("(C) α-Diversity Indices")
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc='lower right')

    # Panel D: Interaction matrix heatmap
    ax4 = fig.add_subplot(gs[1, 2])
    A = results['interaction_matrix']
    im = ax4.imshow(A, cmap='RdBu_r', aspect='auto',
                     vmin=-0.001, vmax=0.001)
    ax4.set_xticks(range(10))
    ax4.set_yticks(range(10))
    ax4.set_xticklabels(SPECIES_SHORT, rotation=90, fontsize=6)
    ax4.set_yticklabels(SPECIES_SHORT, fontsize=6)
    ax4.set_title('(D) Interaction Matrix')
    plt.colorbar(im, ax=ax4, shrink=0.8, label='Interaction strength')

    plt.savefig('figures/fig2_glv_community.png', bbox_inches='tight')
    plt.savefig('figures/fig2_glv_community.svg', bbox_inches='tight')
    plt.close()

    # Save results
    ss_serializable = {
        'species': ss['species'],
        'relative_abundance': ss['relative_abundance'].tolist(),
        'shannon_diversity': float(ss['shannon_diversity']),
        'simpson_diversity': float(ss['simpson_diversity']),
        'dominant_species': ss['dominant_species'],
    }
    with open('results/glv_steady_state.json', 'w') as f:
        json.dump(ss_serializable, f, indent=2)

    print(f"  → Steady-state Shannon: {ss['shannon_diversity']:.3f}")
    print(f"  → Dominant species: {ss['dominant_species']}")
    return results, ss


# ── Figure 3: SCFA Flux Prediction ───────────────────────────────────────

def generate_scfa_figures(glv_results):
    print("[3/6] Running SCFA flux prediction...")

    scfa_results = compute_scfa_timecourse(glv_results)
    rates = scfa_results['production_rates']
    conc = scfa_results['concentrations']
    time_days = scfa_results['time'] / 24

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Short-Chain Fatty Acid (SCFA) Production Flux', fontsize=14, fontweight='bold')

    # Panel A: SCFA production rates over time
    ax = axes[0, 0]
    ax.plot(time_days, rates['acetate_rate'], label='Acetate', color='#e74c3c', linewidth=2)
    ax.plot(time_days, rates['propionate_rate'], label='Propionate', color='#3498db', linewidth=2)
    ax.plot(time_days, rates['butyrate_rate'], label='Butyrate', color='#2ecc71', linewidth=2)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Production Rate (mmol/h)')
    ax.set_title('(A) SCFA Production Rates')
    ax.legend()

    # Panel B: SCFA concentrations
    ax = axes[0, 1]
    ax.plot(time_days, conc['acetate_mM'], label='Acetate', color='#e74c3c', linewidth=2)
    ax.plot(time_days, conc['propionate_mM'], label='Propionate', color='#3498db', linewidth=2)
    ax.plot(time_days, conc['butyrate_total_mM'], label='Butyrate', color='#2ecc71', linewidth=2)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Concentration (mM)')
    ax.set_title('(B) SCFA Concentrations')
    ax.legend()

    # Panel C: SCFA ratios (stacked area)
    ax = axes[1, 0]
    ratios = rates['scfa_ratios']
    ax.stackplot(time_days,
                 ratios['acetate_fraction'],
                 ratios['propionate_fraction'],
                 ratios['butyrate_fraction'],
                 labels=['Acetate', 'Propionate', 'Butyrate'],
                 colors=['#e74c3c', '#3498db', '#2ecc71'], alpha=0.8)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Fraction')
    ax.set_title('(C) SCFA Molar Ratios')
    ax.set_ylim(0, 1)
    ax.legend(loc='center right')

    # Panel D: Butyrate utilization
    ax = axes[1, 1]
    ax.plot(time_days, conc['butyrate_colonocyte_mM'],
            label='Colonocyte utilization', color='#27ae60', linewidth=2)
    ax.plot(time_days, conc['butyrate_systemic_mM'],
            label='Systemic absorption', color='#f39c12', linewidth=2)
    ax.fill_between(time_days, 0, conc['butyrate_colonocyte_mM'],
                     alpha=0.2, color='#27ae60')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Butyrate (mM)')
    ax.set_title('(D) Butyrate Partitioning')
    ax.legend()

    plt.tight_layout()
    plt.savefig('figures/fig3_scfa_flux.png', bbox_inches='tight')
    plt.savefig('figures/fig3_scfa_flux.svg', bbox_inches='tight')
    plt.close()

    # Save results
    final_scfa = {
        'acetate_mM': float(np.mean(conc['acetate_mM'][-100:])),
        'propionate_mM': float(np.mean(conc['propionate_mM'][-100:])),
        'butyrate_total_mM': float(np.mean(conc['butyrate_total_mM'][-100:])),
        'total_scfa_mM': float(np.mean(conc['total_scfa_mM'][-100:])),
        'acetate_fraction': float(np.mean(ratios['acetate_fraction'][-100:])),
        'propionate_fraction': float(np.mean(ratios['propionate_fraction'][-100:])),
        'butyrate_fraction': float(np.mean(ratios['butyrate_fraction'][-100:])),
    }
    with open('results/scfa_steady_state.json', 'w') as f:
        json.dump(final_scfa, f, indent=2)

    print(f"  → Steady-state SCFA (mM): A={final_scfa['acetate_mM']:.2f}, "
          f"P={final_scfa['propionate_mM']:.2f}, B={final_scfa['butyrate_total_mM']:.2f}")
    return scfa_results


# ── Figure 4: Diet Pattern Comparison ─────────────────────────────────────

def generate_diet_figures():
    print("[4/6] Running diet comparison simulation...")

    comparison = simulate_diet_comparison(duration_days=30)
    metrics = compute_diet_impact_metrics(comparison)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Diet Pattern Impact on Gut Microbiota', fontsize=14, fontweight='bold')

    # Panel A: Shannon diversity comparison
    ax = axes[0, 0]
    diet_names = []
    shannon_values = []
    diet_colors = []
    for dk, m in metrics.items():
        diet_names.append(m['diet_name'])
        shannon_values.append(m['shannon_diversity'])
        diet_colors.append(DIET_PATTERNS[dk]['color'])
    bars = ax.bar(range(len(diet_names)), shannon_values, color=diet_colors)
    ax.set_xticks(range(len(diet_names)))
    ax.set_xticklabels(diet_names, rotation=30, ha='right')
    ax.set_ylabel('Shannon Diversity Index')
    ax.set_title('(A) α-Diversity by Diet Pattern')
    for bar, val in zip(bars, shannon_values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    # Panel B: Butyrate producers fraction
    ax = axes[0, 1]
    butyrate_fracs = [metrics[dk]['butyrate_producers_fraction'] * 100 for dk in metrics]
    pathobiont_fracs = [metrics[dk]['pathobiont_fraction'] * 100 for dk in metrics]
    x = np.arange(len(diet_names))
    width = 0.35
    ax.bar(x - width/2, butyrate_fracs, width, label='Butyrate Producers', color='#2ecc71')
    ax.bar(x + width/2, pathobiont_fracs, width, label='Pathobionts', color='#e74c3c')
    ax.set_xticks(x)
    ax.set_xticklabels(diet_names, rotation=30, ha='right')
    ax.set_ylabel('Relative Abundance (%)')
    ax.set_title('(B) Beneficial vs. Pathobiont Species')
    ax.legend()

    # Panel C: Diversity dynamics over time
    ax = axes[1, 0]
    for dk in comparison:
        r = comparison[dk]
        ax.plot(r['time'] / 24, r['shannon_diversity'],
                label=DIET_PATTERNS[dk]['name'],
                color=DIET_PATTERNS[dk]['color'], linewidth=2)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Shannon Diversity')
    ax.set_title('(C) Diversity Dynamics Over 30 Days')
    ax.legend()

    # Panel D: Diet switch simulation (Western → Mediterranean)
    ax = axes[1, 1]
    switch_results = simulate_diet_switch(
        diet_sequence=['western', 'mediterranean', 'western'],
        phase_durations=[14, 28, 14],
    )
    time_days = switch_results['time'] / 24
    ax.plot(time_days, switch_results['shannon_diversity'], 'k-', linewidth=2)
    for boundary in switch_results['phase_boundaries'][1:-1]:
        ax.axvline(x=boundary / 24, color='gray', linestyle='--', alpha=0.7)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Shannon Diversity')
    ax.set_title('(D) Diet Switch: Western → Mediterranean → Western')
    ax.annotate('Western', xy=(7, ax.get_ylim()[0]),
                ha='center', fontsize=9, color='#e74c3c')
    ax.annotate('Mediterranean', xy=(28, ax.get_ylim()[0]),
                ha='center', fontsize=9, color='#2ecc71')
    ax.annotate('Western', xy=(49, ax.get_ylim()[0]),
                ha='center', fontsize=9, color='#e74c3c')

    plt.tight_layout()
    plt.savefig('figures/fig4_diet_comparison.png', bbox_inches='tight')
    plt.savefig('figures/fig4_diet_comparison.svg', bbox_inches='tight')
    plt.close()

    # Save metrics
    with open('results/diet_impact_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    for dk, m in metrics.items():
        print(f"  → {m['diet_name']}: Shannon={m['shannon_diversity']:.3f}, "
              f"Butyrate producers={m['butyrate_producers_fraction']:.1%}")
    return comparison, metrics


# ── Figure 5: Probiotic/Prebiotic Effects ─────────────────────────────────

def generate_probiotic_figures():
    print("[5/6] Running probiotic/prebiotic simulations...")

    # Run prebiotic comparisons
    prebiotic_results = {}
    for pk in ['inulin', 'galactooligosaccharides', 'resistant_starch', 'pectin']:
        prebiotic_results[pk] = simulate_prebiotic_intervention(pk, duration_days=28)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Probiotic/Prebiotic Intervention Effects', fontsize=14, fontweight='bold')

    # Panel A: Prebiotic effect on Shannon diversity
    ax = axes[0, 0]
    prebiotic_names = []
    baseline_divs = []
    intervention_divs = []
    for pk, result in prebiotic_results.items():
        ss_base = compute_steady_state_composition(result['baseline'])
        ss_int = compute_steady_state_composition(result['intervention'])
        prebiotic_names.append(PREBIOTIC_SUBSTRATES[pk]['description'][:15])
        baseline_divs.append(ss_base['shannon_diversity'])
        intervention_divs.append(ss_int['shannon_diversity'])

    x = np.arange(len(prebiotic_names))
    width = 0.35
    ax.bar(x - width/2, baseline_divs, width, label='Baseline', color='#95a5a6')
    ax.bar(x + width/2, intervention_divs, width, label='With Prebiotic', color='#2ecc71')
    ax.set_xticks(x)
    ax.set_xticklabels(prebiotic_names, rotation=30, ha='right')
    ax.set_ylabel('Shannon Diversity')
    ax.set_title('(A) Prebiotic Effect on Diversity')
    ax.legend()

    # Panel B: Probiotic LGG effect over time
    ax = axes[0, 1]
    lgg = simulate_probiotic_intervention('lactobacillus_rhamnosus_gg', duration_days=28)
    time_days = lgg['baseline']['time'] / 24
    ax.plot(time_days, lgg['baseline']['shannon_diversity'],
            'k--', label='Baseline', linewidth=2)
    ax.plot(time_days, lgg['intervention']['shannon_diversity'],
            '#2ecc71', label='+ L. rhamnosus GG', linewidth=2)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Shannon Diversity')
    ax.set_title('(B) LGG Probiotic Effect on Diversity')
    ax.legend()

    # Panel C: Synbiotic vs individual interventions
    ax = axes[1, 0]
    synbiotic = simulate_synbiotic('bifidobacterium_longum_bb536', 'inulin', duration_days=28)
    bb536 = simulate_probiotic_intervention('bifidobacterium_longum_bb536', duration_days=28)
    inulin = prebiotic_results['inulin']

    ss_base = compute_steady_state_composition(synbiotic['baseline'])
    ss_probiotic = compute_steady_state_composition(bb536['intervention'])
    ss_prebiotic = compute_steady_state_composition(inulin['intervention'])
    ss_synbiotic = compute_steady_state_composition(synbiotic['intervention'])

    conditions = ['Baseline', 'Probiotic\n(BB536)', 'Prebiotic\n(Inulin)', 'Synbiotic\n(BB536+Inulin)']
    values = [ss_base['shannon_diversity'], ss_probiotic['shannon_diversity'],
              ss_prebiotic['shannon_diversity'], ss_synbiotic['shannon_diversity']]
    bar_colors = ['#95a5a6', '#3498db', '#2ecc71', '#9b59b6']
    bars = ax.bar(conditions, values, color=bar_colors)
    ax.set_ylabel('Shannon Diversity')
    ax.set_title('(C) Synbiotic Synergy Effect')
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    # Panel D: Species-level response to inulin
    ax = axes[1, 1]
    base_ss = compute_steady_state_composition(inulin['baseline'])
    int_ss = compute_steady_state_composition(inulin['intervention'])
    fold_change = np.log2(
        (int_ss['relative_abundance'] + 1e-6) /
        (base_ss['relative_abundance'] + 1e-6)
    )
    colors_fc = ['#2ecc71' if fc > 0 else '#e74c3c' for fc in fold_change]
    ax.barh(range(10), fold_change, color=colors_fc)
    ax.set_yticks(range(10))
    ax.set_yticklabels(SPECIES_SHORT, fontsize=8)
    ax.set_xlabel('Log₂ Fold Change')
    ax.set_title('(D) Species Response to Inulin')
    ax.axvline(x=0, color='k', linewidth=0.5)

    plt.tight_layout()
    plt.savefig('figures/fig5_probiotic_prebiotic.png', bbox_inches='tight')
    plt.savefig('figures/fig5_probiotic_prebiotic.svg', bbox_inches='tight')
    plt.close()

    # Save results
    probiotic_metrics = {
        'synbiotic_effect': {
            'baseline_shannon': float(ss_base['shannon_diversity']),
            'probiotic_only': float(ss_probiotic['shannon_diversity']),
            'prebiotic_only': float(ss_prebiotic['shannon_diversity']),
            'synbiotic': float(ss_synbiotic['shannon_diversity']),
        },
        'inulin_fold_changes': {
            name: float(fc) for name, fc in zip(SPECIES_SHORT, fold_change)
        }
    }
    with open('results/probiotic_prebiotic_metrics.json', 'w') as f:
        json.dump(probiotic_metrics, f, indent=2)

    print(f"  → Synbiotic Shannon: {ss_synbiotic['shannon_diversity']:.3f}")
    return prebiotic_results, probiotic_metrics


# ── Figure 6: Fermented Food Case Study ───────────────────────────────────

def generate_fermented_food_figures():
    print("[6/6] Running fermented food case study...")

    comparison = run_fermented_food_comparison(duration_days=70)
    all_metrics = {}
    for fk, results in comparison.items():
        all_metrics[fk] = compute_diversity_change_metrics(results)

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Fermented Food Intake: Impact on Gut Microbiota Diversity',
                 fontsize=14, fontweight='bold')

    food_colors = {
        'yogurt': '#f1c40f', 'kimchi': '#e74c3c', 'kefir': '#3498db',
        'natto': '#8e44ad', 'mixed_fermented': '#2ecc71'
    }

    # Panel A: Shannon diversity timecourse for mixed fermented foods
    ax = axes[0, 0]
    mixed = comparison['mixed_fermented']
    ax.plot(mixed['time_days'], mixed['shannon_diversity'], 'k-', linewidth=2)
    for bd in mixed['phase_boundaries_days'][1:-1]:
        ax.axvline(x=bd, color='gray', linestyle='--', alpha=0.7)
    ax.fill_between([0, mixed['phase_boundaries_days'][1]],
                     *ax.get_ylim(), alpha=0.1, color='gray', label='Baseline')
    ax.fill_between([mixed['phase_boundaries_days'][1], mixed['phase_boundaries_days'][2]],
                     *ax.get_ylim(), alpha=0.1, color='green', label='Intervention')
    ax.fill_between([mixed['phase_boundaries_days'][2], mixed['phase_boundaries_days'][3]],
                     *ax.get_ylim(), alpha=0.1, color='orange', label='Washout')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Shannon Diversity')
    ax.set_title('(A) Mixed Fermented Foods (6+ servings/day)')
    ax.legend(fontsize=7)

    # Panel B: Diversity change by food type
    ax = axes[0, 1]
    food_names = [all_metrics[fk]['food'] for fk in all_metrics]
    div_changes = [all_metrics[fk]['diversity_change_pct'] for fk in all_metrics]
    fc_colors = [food_colors.get(fk, '#95a5a6') for fk in all_metrics]
    bars = ax.bar(range(len(food_names)), div_changes, color=fc_colors)
    ax.set_xticks(range(len(food_names)))
    ax.set_xticklabels(food_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Diversity Change (%)')
    ax.set_title('(B) Diversity Change by Fermented Food')
    ax.axhline(y=0, color='k', linewidth=0.5)

    # Panel C: Bray-Curtis dissimilarity
    ax = axes[0, 2]
    bc_values = [all_metrics[fk]['bray_curtis_dissimilarity'] for fk in all_metrics]
    ax.bar(range(len(food_names)), bc_values, color=fc_colors)
    ax.set_xticks(range(len(food_names)))
    ax.set_xticklabels(food_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Bray-Curtis Dissimilarity')
    ax.set_title('(C) Community Composition Shift')

    # Panel D: Species composition - baseline vs intervention (mixed fermented)
    ax = axes[1, 0]
    mixed_metrics = all_metrics['mixed_fermented']
    time_days = mixed['time_days']
    for i in range(10):
        ax.plot(time_days, mixed['abundances'][i],
                label=SPECIES_SHORT[i], linewidth=1.2)
    for bd in mixed['phase_boundaries_days'][1:-1]:
        ax.axvline(x=bd, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Abundance (a.u.)')
    ax.set_title('(D) Species Dynamics - Mixed Fermented')
    ax.legend(fontsize=6, ncol=2)

    # Panel E: SCFA changes
    ax = axes[1, 1]
    scfa_baseline = mixed['scfa']['baseline']
    scfa_interv = mixed['scfa']['intervention']
    scfa_names = ['Acetate', 'Propionate', 'Butyrate']
    baseline_means = [
        np.mean(scfa_baseline['acetate_rate'][-50:]),
        np.mean(scfa_baseline['propionate_rate'][-50:]),
        np.mean(scfa_baseline['butyrate_rate'][-50:]),
    ]
    interv_means = [
        np.mean(scfa_interv['acetate_rate'][-50:]),
        np.mean(scfa_interv['propionate_rate'][-50:]),
        np.mean(scfa_interv['butyrate_rate'][-50:]),
    ]
    x = np.arange(3)
    width = 0.35
    ax.bar(x - width/2, baseline_means, width, label='Baseline', color='#95a5a6')
    ax.bar(x + width/2, interv_means, width, label='Fermented Foods', color='#2ecc71')
    ax.set_xticks(x)
    ax.set_xticklabels(scfa_names)
    ax.set_ylabel('Production Rate (mmol/h)')
    ax.set_title('(E) SCFA Production: Baseline vs. Intervention')
    ax.legend()

    # Panel F: Summary heatmap
    ax = axes[1, 2]
    heatmap_data = np.array([
        [all_metrics[fk]['baseline_shannon'] for fk in all_metrics],
        [all_metrics[fk]['intervention_shannon'] for fk in all_metrics],
        [all_metrics[fk]['washout_shannon'] for fk in all_metrics],
        [all_metrics[fk]['diversity_change_pct'] for fk in all_metrics],
    ])
    im = ax.imshow(heatmap_data, cmap='YlGn', aspect='auto')
    ax.set_xticks(range(len(food_names)))
    ax.set_xticklabels(food_names, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(4))
    ax.set_yticklabels(['Baseline H\'', 'Intervention H\'', 'Washout H\'', 'Change (%)'], fontsize=8)
    ax.set_title('(F) Summary Heatmap')
    plt.colorbar(im, ax=ax, shrink=0.8)
    for i in range(4):
        for j in range(len(food_names)):
            ax.text(j, i, f'{heatmap_data[i, j]:.1f}',
                    ha='center', va='center', fontsize=7,
                    color='white' if heatmap_data[i, j] > heatmap_data.max() * 0.6 else 'black')

    plt.tight_layout()
    plt.savefig('figures/fig6_fermented_food.png', bbox_inches='tight')
    plt.savefig('figures/fig6_fermented_food.svg', bbox_inches='tight')
    plt.close()

    # Save metrics
    with open('results/fermented_food_metrics.json', 'w') as f:
        json.dump(all_metrics, f, indent=2)

    for fk, m in all_metrics.items():
        print(f"  → {m['food']}: ΔDiversity={m['diversity_change_pct']:+.1f}%, "
              f"BC={m['bray_curtis_dissimilarity']:.3f}")
    return comparison, all_metrics


# ── MICOM/gapseq Analysis ─────────────────────────────────────────────────

def generate_micom_results():
    print("[MICOM] Running community metabolic analysis...")

    analysis = run_community_metabolic_analysis()
    report = generate_gapseq_report(analysis)

    with open('results/micom_community_fba.txt', 'w') as f:
        f.write(report)

    # Save structured results
    serializable = {
        'models': analysis['models'],
        'n_reactions': analysis['n_reactions'],
        'n_metabolites': analysis['n_metabolites'],
        'species_growth_rates': analysis['community_fba']['species_growth_rates'],
        'community_growth_rate': analysis['community_fba']['community_growth_rate'],
        'cross_feeding_pairs': analysis['community_fba']['cross_feeding']['cross_feeding_pairs'],
        'exchange_summary': analysis['exchange_summary'],
    }
    with open('results/micom_analysis.json', 'w') as f:
        json.dump(serializable, f, indent=2)

    print(report)
    return analysis


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Gut Microbiota-Diet Interaction Systems Biology Framework")
    print("=" * 70)
    print(f"Start time: {datetime.now().isoformat()}\n")

    ensure_dirs()

    # Run all simulations and generate figures
    shime_results, shime_eff = generate_shime_figures()
    glv_results, glv_ss = generate_glv_figures()
    scfa_results = generate_scfa_figures(glv_results)
    diet_results, diet_metrics = generate_diet_figures()
    prob_results, prob_metrics = generate_probiotic_figures()
    ferm_results, ferm_metrics = generate_fermented_food_figures()
    micom_analysis = generate_micom_results()

    # Write process log
    log_entry = {
        'timestamp': datetime.now().isoformat(),
        'phase': 'complete',
        'event_type': 'run_completed',
        'actor': 'co-scientist',
        'skill_or_tool': 'gut-microbiome-framework',
        'files_written': [
            'figures/fig1_shime_digestion.png',
            'figures/fig1_shime_digestion.svg',
            'figures/fig2_glv_community.png',
            'figures/fig2_glv_community.svg',
            'figures/fig3_scfa_flux.png',
            'figures/fig3_scfa_flux.svg',
            'figures/fig4_diet_comparison.png',
            'figures/fig4_diet_comparison.svg',
            'figures/fig5_probiotic_prebiotic.png',
            'figures/fig5_probiotic_prebiotic.svg',
            'figures/fig6_fermented_food.png',
            'figures/fig6_fermented_food.svg',
            'results/shime_absorption_efficiency.json',
            'results/glv_steady_state.json',
            'results/scfa_steady_state.json',
            'results/diet_impact_metrics.json',
            'results/probiotic_prebiotic_metrics.json',
            'results/fermented_food_metrics.json',
            'results/micom_community_fba.txt',
            'results/micom_analysis.json',
        ],
        'status': 'ok',
    }
    with open('logs/process-log.jsonl', 'w') as f:
        f.write(json.dumps(log_entry) + '\n')

    print(f"\n{'=' * 70}")
    print("All simulations completed successfully!")
    print(f"End time: {datetime.now().isoformat()}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
