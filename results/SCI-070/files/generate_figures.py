#!/usr/bin/env python3
"""
Visualization module for Ecosystem Service Valuation Framework
Generates publication-quality figures for report.md
"""

import json
import math
import os

# Use non-interactive backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def setup_style():
    """Set publication-quality plot style"""
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 13,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.spines.top': False,
        'axes.spines.right': False,
    })


def fig1_service_valuation_summary(results: dict, output_dir: str):
    """Figure 1: Ecosystem service monetary valuation summary (bar chart)"""
    services = {
        'Carbon\nStorage': results['invest_results']['carbon']['summary']['monetary_value_social_JPY'] * 0.02,
        'Erosion\nControl': results['invest_results']['sediment']['summary']['monetary_value_JPY'],
        'Water\nYield': results['invest_results']['water']['summary']['monetary_value_JPY'],
        'Nutrient\nRetention': results['invest_results']['nutrient']['summary']['monetary_value_total_JPY'],
        'Pollination': results['invest_results']['pollination']['summary']['monetary_value_JPY'],
        'Recreation': results['invest_results']['recreation']['summary']['monetary_value_JPY'],
    }

    categories = {
        'Carbon\nStorage': 'Regulating',
        'Erosion\nControl': 'Regulating',
        'Water\nYield': 'Provisioning',
        'Nutrient\nRetention': 'Regulating',
        'Pollination': 'Regulating',
        'Recreation': 'Cultural',
    }

    colors_map = {'Provisioning': '#2196F3', 'Regulating': '#4CAF50', 'Cultural': '#FF9800'}
    colors = [colors_map[categories[s]] for s in services.keys()]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(services.keys(), [v / 1e6 for v in services.values()], color=colors, edgecolor='white', linewidth=0.5)

    for bar, val in zip(bars, services.values()):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'¥{val/1e6:.0f}M', ha='center', va='bottom', fontsize=9)

    ax.set_ylabel('Annual Value (Million JPY)')
    ax.set_title('Ecosystem Service Monetary Valuation — Satoyama Case Study')

    legend_patches = [mpatches.Patch(color=c, label=l) for l, c in colors_map.items()]
    ax.legend(handles=legend_patches, loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig1_service_valuation.png'))
    plt.savefig(os.path.join(output_dir, 'fig1_service_valuation.svg'))
    plt.close()
    print("  Saved fig1_service_valuation.png/svg")


def fig2_carbon_pools(results: dict, output_dir: str):
    """Figure 2: Carbon storage by land cover and pool type (stacked bar)"""
    carbon = results['invest_results']['carbon']
    lc_order = ['mixed_forest', 'paddy_field', 'upland_crop', 'bamboo_grove',
                'irrigation_pond', 'grassland_meadow', 'settlement', 'other']
    lc_labels = ['Mixed\nForest', 'Paddy\nField', 'Upland\nCrop', 'Bamboo\nGrove',
                 'Irrigation\nPond', 'Grassland', 'Settlement', 'Other']

    pools = ['above_ground_tC', 'below_ground_tC', 'soil_tC', 'dead_matter_tC']
    pool_labels = ['Above-ground', 'Below-ground', 'Soil', 'Dead matter']
    pool_colors = ['#66BB6A', '#43A047', '#8D6E63', '#A1887F']

    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(len(lc_order))

    for pool, label, color in zip(pools, pool_labels, pool_colors):
        values = [carbon[lc][pool] / 1000 for lc in lc_order]
        ax.bar(lc_labels, values, bottom=bottom, label=label, color=color, edgecolor='white', linewidth=0.5)
        bottom += np.array(values)

    ax.set_ylabel('Carbon Storage (×1,000 tC)')
    ax.set_title('Carbon Storage by Land Cover Type and Pool')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig2_carbon_pools.png'))
    plt.savefig(os.path.join(output_dir, 'fig2_carbon_pools.svg'))
    plt.close()
    print("  Saved fig2_carbon_pools.png/svg")


def fig3_wtp_forest_plot(results: dict, output_dir: str):
    """Figure 3: WTP estimates with confidence intervals (forest plot)"""
    wtp = results['wtp_estimates']['wtp_estimates_jpy_household_year']
    labels = list(wtp.keys())
    means = [wtp[k]['mean'] for k in labels]
    ci_low = [wtp[k]['ci95_low'] for k in labels]
    ci_high = [wtp[k]['ci95_high'] for k in labels]

    display_labels = [
        'Biodiversity +20%', 'Biodiversity +40%',
        'Water: Standard', 'Water: Swimmable',
        'Landscape: Partial', 'Landscape: Traditional',
        'Carbon +30%',
    ]

    fig, ax = plt.subplots(figsize=(9, 6))
    y_pos = range(len(labels))
    xerr_low = [m - l for m, l in zip(means, ci_low)]
    xerr_high = [h - m for m, h in zip(means, ci_high)]

    ax.errorbar(means, y_pos, xerr=[xerr_low, xerr_high],
                fmt='o', color='#1565C0', markersize=8, capsize=5, linewidth=2)
    ax.axvline(x=0, color='grey', linestyle='--', alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(display_labels)
    ax.set_xlabel('WTP (JPY / household / year)')
    ax.set_title('Willingness-to-Pay Estimates — Choice Experiment Results')
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig3_wtp_estimates.png'))
    plt.savefig(os.path.join(output_dir, 'fig3_wtp_estimates.svg'))
    plt.close()
    print("  Saved fig3_wtp_estimates.png/svg")


def fig4_discount_rate_sensitivity(results: dict, output_dir: str):
    """Figure 4: NPV sensitivity to discount rates"""
    sa = results['sensitivity_analysis']['discount_rate_sensitivity']
    rates = sorted([float(r) for r in sa.keys()])
    npvs = [sa[str(r)]['npv_50y_JPY'] / 1e9 for r in rates]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(rates, npvs, 'o-', color='#D32F2F', linewidth=2, markersize=8)
    ax.fill_between(rates, npvs, alpha=0.1, color='#D32F2F')

    for r, npv in zip(rates, npvs):
        ax.annotate(f'¥{npv:.1f}B', (r, npv), textcoords="offset points",
                    xytext=(0, 12), ha='center', fontsize=9)

    ax.set_xlabel('Discount Rate')
    ax.set_ylabel('NPV over 50 years (Billion JPY)')
    ax.set_title('Sensitivity of Ecosystem Asset Value to Discount Rate')
    ax.set_xticks(rates)
    ax.set_xticklabels([f'{r:.0%}' for r in rates])
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig4_discount_sensitivity.png'))
    plt.savefig(os.path.join(output_dir, 'fig4_discount_sensitivity.svg'))
    plt.close()
    print("  Saved fig4_discount_sensitivity.png/svg")


def fig5_land_use_scenarios(results: dict, output_dir: str):
    """Figure 5: Land use change scenario comparison (radar/spider chart)"""
    scenarios = results['sensitivity_analysis']['land_use_scenarios']
    baseline = scenarios['baseline']

    labels = ['Carbon\nStorage', 'Erosion\nControl', 'Water\nYield', 'Nitrogen\nRetention']
    metrics = ['carbon_tC', 'erosion_avoided_ton', 'water_yield_m3', 'N_retained_kg']
    n = len(labels)
    angles = [i * 2 * math.pi / n for i in range(n)]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    scenario_colors = {'baseline': '#4CAF50', 'abandonment': '#FF9800',
                       'conservation': '#2196F3', 'urbanization': '#F44336'}

    for scenario_name, scenario_data in scenarios.items():
        values = [scenario_data[m] / baseline[m] for m in metrics]
        values += values[:1]
        ax.plot(angles, values, 'o-', linewidth=2, label=scenario_name.capitalize(),
                color=scenario_colors.get(scenario_name, '#999'))
        ax.fill(angles, values, alpha=0.1, color=scenario_colors.get(scenario_name, '#999'))

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.set_title('Land Use Scenario Comparison\n(Ratio to Baseline)', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig5_scenario_comparison.png'))
    plt.savefig(os.path.join(output_dir, 'fig5_scenario_comparison.svg'))
    plt.close()
    print("  Saved fig5_scenario_comparison.png/svg")


def fig6_seea_condition(results: dict, output_dir: str):
    """Figure 6: SEEA-EA Ecosystem Condition indicators"""
    cond = results['seea_ea_accounts']['condition_account']['indicators']
    labels = ['Vegetation\nCover', 'Species\nRichness', 'Soil Organic\nCarbon',
              'Water\nQuality', 'Landscape\nConnectivity']
    values = list(cond.values())
    composite = results['seea_ea_accounts']['condition_account']['composite_condition_index']

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ['#66BB6A' if v >= 0.7 else '#FFA726' if v >= 0.5 else '#EF5350' for v in values]
    bars = ax.barh(labels, values, color=colors, edgecolor='white', linewidth=0.5)

    for bar, val in zip(bars, values):
        ax.text(val + 0.02, bar.get_y() + bar.get_height()/2,
                f'{val:.2f}', ha='left', va='center', fontsize=10)

    ax.axvline(x=composite, color='#1565C0', linestyle='--', linewidth=2,
               label=f'Composite Index = {composite:.3f}')
    ax.set_xlim(0, 1.1)
    ax.set_xlabel('Condition Score (0 = Degraded, 1 = Reference)')
    ax.set_title('SEEA-EA Ecosystem Condition Account — Satoyama')
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig6_seea_condition.png'))
    plt.savefig(os.path.join(output_dir, 'fig6_seea_condition.svg'))
    plt.close()
    print("  Saved fig6_seea_condition.png/svg")


def generate_all_figures(results_path: str = 'results/pipeline_results.json',
                          output_dir: str = 'figures'):
    """Generate all figures for the report"""
    setup_style()
    os.makedirs(output_dir, exist_ok=True)

    with open(results_path, 'r', encoding='utf-8') as f:
        results = json.load(f)

    print("Generating figures...")
    fig1_service_valuation_summary(results, output_dir)
    fig2_carbon_pools(results, output_dir)
    fig3_wtp_forest_plot(results, output_dir)
    fig4_discount_rate_sensitivity(results, output_dir)
    fig5_land_use_scenarios(results, output_dir)
    fig6_seea_condition(results, output_dir)
    print(f"\nAll figures saved to {output_dir}/")


if __name__ == "__main__":
    generate_all_figures()
