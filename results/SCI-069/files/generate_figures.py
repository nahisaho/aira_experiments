"""
Figure Generation for UHI Simulation Results
Creates all publication-quality figures.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
from wrf_ucm_simulator import run_all_scenarios, WRFUCMSimulator
from anthropogenic_heat import AnthropogenicHeatModel
from mitigation import MitigationScenario

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150
FIGDIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGDIR, exist_ok=True)


def plot_building_morphology(sim):
    """Figure 1: Building morphology maps."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    data = [
        (sim.ucm.building_height, 'Building Height [m]', 'YlOrRd'),
        (sim.ucm.building_fraction, 'Building Plan Area Fraction', 'Oranges'),
        (sim.ucm.canyon_aspect, 'Canyon Aspect Ratio (H/W)', 'RdPu'),
        (sim.ucm.sky_view_factor, 'Sky View Factor', 'Blues_r'),
    ]

    for ax, (field, title, cmap) in zip(axes.flatten(), data):
        im = ax.imshow(field, cmap=cmap, origin='lower')
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('Grid X (500m)')
        ax.set_ylabel('Grid Y (500m)')
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle('Tokyo Urban Morphology Parameters', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig1_morphology.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig1_morphology.png")


def plot_anthropogenic_heat(sim):
    """Figure 2: Anthropogenic heat components and diurnal profile."""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 3, figure=fig)

    # Spatial maps at 14:00
    components = sim.ah_model.get_component_breakdown(14)
    titles = ['Traffic', 'HVAC', 'Industry']
    keys = ['traffic', 'hvac', 'industry']

    for idx, (key, title) in enumerate(zip(keys, titles)):
        ax = fig.add_subplot(gs[0, idx])
        im = ax.imshow(components[key], cmap='hot', origin='lower')
        ax.set_title(f'{title} Heat (14:00) [W/m²]')
        plt.colorbar(im, ax=ax, shrink=0.8)

    # Diurnal profile
    ax_diurnal = fig.add_subplot(gs[1, :])
    hours = np.arange(0, 24, 0.5)
    traffic_vals, hvac_vals, industry_vals, total_vals = [], [], [], []

    for h in hours:
        comps = sim.ah_model.get_component_breakdown(h)
        traffic_vals.append(np.mean(comps['traffic']))
        hvac_vals.append(np.mean(comps['hvac']))
        industry_vals.append(np.mean(comps['industry']))
        total_vals.append(np.mean(sum(comps.values())))

    ax_diurnal.fill_between(hours, 0, traffic_vals, alpha=0.4, label='Traffic')
    ax_diurnal.fill_between(hours, traffic_vals,
                           np.array(traffic_vals) + np.array(hvac_vals),
                           alpha=0.4, label='HVAC')
    ax_diurnal.fill_between(hours, np.array(traffic_vals) + np.array(hvac_vals),
                           np.array(traffic_vals) + np.array(hvac_vals) + np.array(industry_vals),
                           alpha=0.4, label='Industry')
    ax_diurnal.plot(hours, total_vals, 'k-', linewidth=2, label='Total')
    ax_diurnal.set_xlabel('Hour of Day')
    ax_diurnal.set_ylabel('Mean Anthropogenic Heat [W/m²]')
    ax_diurnal.set_title('Diurnal Profile of Anthropogenic Heat Emissions')
    ax_diurnal.legend()
    ax_diurnal.set_xlim(0, 24)
    ax_diurnal.grid(True, alpha=0.3)

    fig.suptitle('Anthropogenic Heat Emission Components', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig2_anthropogenic_heat.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig2_anthropogenic_heat.png")


def plot_uhi_diurnal(scenarios):
    """Figure 3: Diurnal UHI intensity comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    colors = {
        'baseline_2020': '#1f77b4',
        'green': '#2ca02c',
        'cool_roof': '#ff7f0e',
        'combined': '#d62728',
        'baseline_2050': '#9467bd',
        'mitigated_2050': '#8c564b'
    }
    labels = {
        'baseline_2020': 'Baseline 2020',
        'green': 'Green Infrastructure',
        'cool_roof': 'Cool Roofs',
        'combined': 'Combined Mitigation',
        'baseline_2050': 'Baseline 2050',
        'mitigated_2050': '2050 + Mitigation'
    }

    for name, result in scenarios.items():
        ax1.plot(result['hours'], result['UHI_mean'],
                color=colors[name], label=labels[name], linewidth=2)
        ax2.plot(result['hours'], result['UHI_max'],
                color=colors[name], label=labels[name], linewidth=2)

    for ax, title in [(ax1, 'Mean UHI Intensity'), (ax2, 'Maximum UHI Intensity')]:
        ax.set_xlabel('Hour of Day')
        ax.set_ylabel('UHI Intensity [K]')
        ax.set_title(title)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 23)

    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig3_uhi_diurnal.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig3_uhi_diurnal.png")


def plot_spatial_uhi(scenarios):
    """Figure 4: Spatial UHI distribution at 14:00."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    names = ['baseline_2020', 'green', 'cool_roof', 'combined', 'baseline_2050', 'mitigated_2050']
    titles = ['Baseline 2020', 'Green Infrastructure', 'Cool Roofs',
              'Combined', 'Baseline 2050', '2050 + Mitigation']

    vmin = min(np.min(scenarios[n]['spatial_UHI'][14]) for n in names)
    vmax = max(np.max(scenarios[n]['spatial_UHI'][14]) for n in names)

    for ax, name, title in zip(axes.flatten(), names, titles):
        uhi = scenarios[name]['spatial_UHI'][14]
        im = ax.imshow(uhi, cmap='RdYlBu_r', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=11)
        plt.colorbar(im, ax=ax, shrink=0.8, label='UHI [K]')

    fig.suptitle('Spatial UHI Distribution at 14:00 JST', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig4_spatial_uhi.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig4_spatial_uhi.png")


def plot_wbgt_analysis(scenarios):
    """Figure 5: WBGT analysis and heat stress risk."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Diurnal WBGT
    ax = axes[0, 0]
    for name in ['baseline_2020', 'combined', 'baseline_2050', 'mitigated_2050']:
        labels_map = {
            'baseline_2020': 'Baseline 2020',
            'combined': 'Combined 2020',
            'baseline_2050': 'Baseline 2050',
            'mitigated_2050': '2050 + Mitigation'
        }
        ax.plot(scenarios[name]['hours'], scenarios[name]['WBGT_max'],
               label=labels_map[name], linewidth=2)
    ax.axhline(y=28, color='orange', linestyle='--', alpha=0.7, label='Caution (28°C)')
    ax.axhline(y=31, color='red', linestyle='--', alpha=0.7, label='Danger (31°C)')
    ax.set_xlabel('Hour of Day')
    ax.set_ylabel('Maximum WBGT [°C]')
    ax.set_title('Diurnal Maximum WBGT')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Risk map at 14:00 - Baseline 2020
    ax = axes[0, 1]
    risk_colors = ['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c', '#8e44ad']
    risk_cmap = mcolors.ListedColormap(risk_colors)
    risk = scenarios['baseline_2020']['risk_map'][14]
    im = ax.imshow(risk, cmap=risk_cmap, origin='lower', vmin=0, vmax=4)
    ax.set_title('Heat Stress Risk (Baseline 2020, 14:00)')
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(['Low', 'Moderate', 'High', 'Very High', 'Extreme'])

    # Risk map - 2050
    ax = axes[1, 0]
    risk_2050 = scenarios['baseline_2050']['risk_map'][14]
    im = ax.imshow(risk_2050, cmap=risk_cmap, origin='lower', vmin=0, vmax=4)
    ax.set_title('Heat Stress Risk (Baseline 2050, 14:00)')
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(['Low', 'Moderate', 'High', 'Very High', 'Extreme'])

    # Risk map - 2050 mitigated
    ax = axes[1, 1]
    risk_m = scenarios['mitigated_2050']['risk_map'][14]
    im = ax.imshow(risk_m, cmap=risk_cmap, origin='lower', vmin=0, vmax=4)
    ax.set_title('Heat Stress Risk (2050 + Mitigation, 14:00)')
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(['Low', 'Moderate', 'High', 'Very High', 'Extreme'])

    fig.suptitle('WBGT Heat Stress Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig5_wbgt_risk.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig5_wbgt_risk.png")


def plot_cooling_effectiveness(scenarios):
    """Figure 6: Cooling effectiveness comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Bar chart of mean/max cooling
    mit_scenarios = {
        'Green': ('baseline_2020', 'green'),
        'Cool Roof': ('baseline_2020', 'cool_roof'),
        'Combined': ('baseline_2020', 'combined'),
        'Combined\n(2050)': ('baseline_2050', 'mitigated_2050'),
    }

    names = list(mit_scenarios.keys())
    mean_cool = []
    max_cool = []

    for base_key, mit_key in mit_scenarios.values():
        base_T = scenarios[base_key]['spatial_T'][14]
        mit_T = scenarios[mit_key]['spatial_T'][14]
        delta = base_T - mit_T
        mean_cool.append(np.mean(delta))
        max_cool.append(np.max(delta))

    x = np.arange(len(names))
    width = 0.35
    ax1.bar(x - width/2, mean_cool, width, label='Mean Cooling', color='#3498db')
    ax1.bar(x + width/2, max_cool, width, label='Max Cooling', color='#e74c3c')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)
    ax1.set_ylabel('Temperature Reduction [K]')
    ax1.set_title('Cooling Effectiveness at 14:00')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    # Spatial cooling pattern for combined scenario
    base_T = scenarios['baseline_2020']['spatial_T'][14]
    comb_T = scenarios['combined']['spatial_T'][14]
    cooling = base_T - comb_T
    im = ax2.imshow(cooling, cmap='Blues', origin='lower')
    ax2.set_title('Spatial Cooling Pattern (Combined, 14:00)')
    plt.colorbar(im, ax=ax2, label='Cooling [K]')

    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig6_cooling_effectiveness.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig6_cooling_effectiveness.png")


def plot_2050_projection(scenarios):
    """Figure 7: 2050 projection comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Temperature difference 2050 vs 2020
    T_2020 = scenarios['baseline_2020']['spatial_T'][14]
    T_2050 = scenarios['baseline_2050']['spatial_T'][14]
    T_2050m = scenarios['mitigated_2050']['spatial_T'][14]

    im1 = axes[0].imshow(T_2020 - 273.15, cmap='hot', origin='lower')
    axes[0].set_title('Canyon Temperature 2020 [°C]')
    plt.colorbar(im1, ax=axes[0], shrink=0.8)

    im2 = axes[1].imshow(T_2050 - 273.15, cmap='hot', origin='lower')
    axes[1].set_title('Canyon Temperature 2050 [°C]')
    plt.colorbar(im2, ax=axes[1], shrink=0.8)

    warming = T_2050 - T_2020
    im3 = axes[2].imshow(warming, cmap='Reds', origin='lower')
    axes[2].set_title('Warming 2050−2020 [K]')
    plt.colorbar(im3, ax=axes[2], shrink=0.8)

    fig.suptitle('Tokyo 2050 UHI Projection', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig7_2050_projection.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig7_2050_projection.png")


def plot_summary_table(scenarios):
    """Figure 8: Summary statistics table."""
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis('off')

    headers = ['Scenario', 'Mean UHI\n[K]', 'Peak UHI\n[K]', 'Peak T\n[°C]',
               'Peak WBGT\n[°C]', 'Mean Cooling\n[K]']
    data = []

    for name, label in [('baseline_2020', 'Baseline 2020'),
                        ('green', 'Green Infra.'),
                        ('cool_roof', 'Cool Roofs'),
                        ('combined', 'Combined'),
                        ('baseline_2050', 'Baseline 2050'),
                        ('mitigated_2050', '2050+Mitig.')]:
        r = scenarios[name]
        mean_uhi = np.mean(r['UHI_mean'])
        peak_uhi = max(r['UHI_max'])
        peak_T = max(r['T_canyon_max'])
        peak_wbgt = max(r['WBGT_max'])

        if name in ('green', 'cool_roof', 'combined'):
            base_T = scenarios['baseline_2020']['spatial_T'][14]
            mit_T = r['spatial_T'][14]
            cooling = np.mean(base_T - mit_T)
        elif name == 'mitigated_2050':
            base_T = scenarios['baseline_2050']['spatial_T'][14]
            mit_T = r['spatial_T'][14]
            cooling = np.mean(base_T - mit_T)
        else:
            cooling = 0.0

        data.append([label, f'{mean_uhi:.2f}', f'{peak_uhi:.2f}',
                     f'{peak_T:.1f}', f'{peak_wbgt:.1f}', f'{cooling:.2f}'])

    table = ax.table(cellText=data, colLabels=headers, loc='center',
                     cellLoc='center', colColours=['#dce6f1']*6)
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    ax.set_title('Summary of Simulation Results', fontsize=14, fontweight='bold', pad=20)
    fig.savefig(os.path.join(FIGDIR, 'fig8_summary_table.png'), bbox_inches='tight')
    plt.close()
    print("Saved fig8_summary_table.png")


if __name__ == "__main__":
    print("=" * 60)
    print("Running UHI Simulation and Generating Figures")
    print("=" * 60)

    # Run all scenarios
    scenarios = run_all_scenarios()

    # Print summary stats
    print("\n" + "=" * 60)
    print("SIMULATION RESULTS SUMMARY")
    print("=" * 60)
    for name, result in scenarios.items():
        peak_hour = np.argmax(result['UHI_max'])
        print(f"\n--- {name} ({result['year']}) ---")
        print(f"  Peak UHI: {result['UHI_max'][peak_hour]:.2f} K at {peak_hour}:00")
        print(f"  Mean daily UHI: {np.mean(result['UHI_mean']):.2f} K")
        print(f"  Peak WBGT: {max(result['WBGT_max']):.1f} °C")
        print(f"  Peak canyon T: {max(result['T_canyon_max']):.1f} °C")

    # Generate all figures
    print("\nGenerating figures...")
    sim = WRFUCMSimulator((50, 50))
    sim.initialize()

    plot_building_morphology(sim)
    plot_anthropogenic_heat(sim)
    plot_uhi_diurnal(scenarios)
    plot_spatial_uhi(scenarios)
    plot_wbgt_analysis(scenarios)
    plot_cooling_effectiveness(scenarios)
    plot_2050_projection(scenarios)
    plot_summary_table(scenarios)

    print("\nAll figures generated successfully!")
    print(f"Figures saved to: {FIGDIR}")
