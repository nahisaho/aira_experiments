"""
Integrated Building Performance Simulation Dashboard
Main entry point for running all simulations and generating visualizations.
"""
import sys
import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(__file__))

from ifc_converter import IFCConverter
from thermal_simulation import ThermalSimulation
from cfd_simulation import CFDSimulation
from daylight_simulation import DaylightSimulation


def create_figures_dir():
    fig_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    return fig_dir


def plot_system_architecture(fig_dir):
    """Create system architecture diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')

    # Title
    ax.text(7, 9.5, 'BIM-Integrated Environmental Performance Simulation System',
            ha='center', va='center', fontsize=14, fontweight='bold')

    # IFC/BIM Layer
    rect = plt.Rectangle((0.5, 7.5), 3, 1.5, fill=True, facecolor='#4472C4', alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.add_patch(rect)
    ax.text(2, 8.25, 'IFC/BIM\nData Source', ha='center', va='center', fontsize=10, color='white', fontweight='bold')

    # IFC Converter
    rect = plt.Rectangle((5, 7.5), 4, 1.5, fill=True, facecolor='#5B9BD5', alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.add_patch(rect)
    ax.text(7, 8.25, 'IFC Model Converter\n(Geometry, Materials, Zones)', ha='center', va='center', fontsize=9, color='white', fontweight='bold')

    # Arrow from BIM to Converter
    ax.annotate('', xy=(5, 8.25), xytext=(3.5, 8.25),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))

    # Simulation Engines
    sim_boxes = [
        (0.5, 4.5, '#ED7D31', 'EnergyPlus\nThermal Load\nSimulation'),
        (4, 4.5, '#70AD47', 'OpenFOAM\nCFD / Natural\nVentilation'),
        (7.5, 4.5, '#FFC000', 'Radiance\nDaylighting\nSimulation'),
        (11, 4.5, '#7030A0', 'Structural\n& HVAC\nAnalysis'),
    ]

    for x, y, color, label in sim_boxes:
        rect = plt.Rectangle((x, y), 3, 2, fill=True, facecolor=color, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x + 1.5, y + 1, label, ha='center', va='center', fontsize=9, color='white', fontweight='bold')

    # Arrows from converter to simulations
    for x, _, _, _ in sim_boxes:
        ax.annotate('', xy=(x + 1.5, 6.5), xytext=(7, 7.5),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=1.5, ls='--'))

    # Integration Dashboard
    rect = plt.Rectangle((2, 1), 10, 2.5, fill=True, facecolor='#C00000', alpha=0.8, edgecolor='black', linewidth=2)
    ax.add_patch(rect)
    ax.text(7, 2.25, 'Integrated Performance Dashboard\nZEB Assessment | Energy Balance | Comfort Analysis | LEED Compliance',
            ha='center', va='center', fontsize=11, color='white', fontweight='bold')

    # Arrows from sims to dashboard
    for x, y, _, _ in sim_boxes:
        ax.annotate('', xy=(x + 1.5, 3.5), xytext=(x + 1.5, 4.5),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'system_architecture.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> system_architecture.png")


def plot_monthly_energy(results, fig_dir):
    """Plot monthly heating and cooling loads."""
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    # Monthly heating and cooling
    x = np.arange(12)
    width = 0.35
    heating = np.array(results['monthly_heating_kWh'])
    cooling = np.array(results['monthly_cooling_kWh'])

    bars1 = ax1.bar(x - width/2, heating, width, label='Heating', color='#ED7D31', edgecolor='black', linewidth=0.5)
    bars2 = ax1.bar(x + width/2, cooling, width, label='Cooling', color='#4472C4', edgecolor='black', linewidth=0.5)

    ax1.set_xlabel('Month', fontsize=11)
    ax1.set_ylabel('Energy Consumption (kWh)', fontsize=11)
    ax1.set_title('Monthly Heating and Cooling Energy Consumption', fontsize=13, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(months)
    ax1.legend(fontsize=10)
    ax1.grid(axis='y', alpha=0.3)

    # EUI breakdown
    categories = ['Heating', 'Cooling', 'Lighting', 'Equipment', 'Fan/Pump']
    values = [
        results['heating_eui'],
        results['cooling_eui'],
        results['annual_lighting_kWh'] / results['total_floor_area'],
        results['annual_equipment_kWh'] / results['total_floor_area'],
        results['annual_fan_kWh'] / results['total_floor_area'],
    ]
    colors = ['#ED7D31', '#4472C4', '#FFC000', '#70AD47', '#7030A0']

    bars = ax2.barh(categories, values, color=colors, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Energy Use Intensity (kWh/m²/year)', fontsize=11)
    ax2.set_title(f'Annual Energy Use Intensity Breakdown (Total EUI: {results["eui_kWh_m2"]} kWh/m²/yr)', fontsize=13, fontweight='bold')
    for bar, val in zip(bars, values):
        ax2.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                 f'{val:.1f}', va='center', fontsize=10)
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'monthly_energy.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> monthly_energy.png")


def plot_cfd_results(cfd_result, fig_dir):
    """Plot CFD velocity and temperature fields."""
    vel = np.array(cfd_result['velocity_field'])
    temp = np.array(cfd_result['temperature_field'])

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Velocity field
    im1 = axes[0].imshow(vel, cmap='jet', aspect='auto', origin='lower',
                         extent=[0, 25, 0, 20])
    axes[0].set_title('Indoor Air Velocity Distribution (m/s)', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Length (m)', fontsize=10)
    axes[0].set_ylabel('Width (m)', fontsize=10)
    plt.colorbar(im1, ax=axes[0], label='Velocity (m/s)')

    # Temperature field
    im2 = axes[1].imshow(temp, cmap='coolwarm', aspect='auto', origin='lower',
                         extent=[0, 25, 0, 20])
    axes[1].set_title('Indoor Temperature Distribution (°C)', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Length (m)', fontsize=10)
    axes[1].set_ylabel('Width (m)', fontsize=10)
    plt.colorbar(im2, ax=axes[1], label='Temperature (°C)')

    plt.suptitle('CFD Natural Ventilation Analysis - Cross Ventilation', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'cfd_results.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> cfd_results.png")


def plot_cross_ventilation_scenarios(scenarios, fig_dir):
    """Plot cross-ventilation scenario comparison."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    speeds = [s['wind_speed'] for s in scenarios]
    achs = [s['ach'] for s in scenarios]
    temps = [s['avg_temp'] for s in scenarios]
    comforts = [s['comfort_fraction'] * 100 for s in scenarios]
    vels = [s['avg_velocity'] for s in scenarios]

    # ACH vs wind speed
    axes[0, 0].plot(speeds, achs, 'o-', color='#4472C4', linewidth=2, markersize=8)
    axes[0, 0].set_xlabel('External Wind Speed (m/s)', fontsize=10)
    axes[0, 0].set_ylabel('Air Changes per Hour (ACH)', fontsize=10)
    axes[0, 0].set_title('Ventilation Rate vs Wind Speed', fontsize=11, fontweight='bold')
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].axhline(y=4, color='red', linestyle='--', alpha=0.7, label='Min. req. (4 ACH)')
    axes[0, 0].legend()

    # Indoor velocity
    axes[0, 1].plot(speeds, vels, 's-', color='#70AD47', linewidth=2, markersize=8)
    axes[0, 1].set_xlabel('External Wind Speed (m/s)', fontsize=10)
    axes[0, 1].set_ylabel('Avg. Indoor Velocity (m/s)', fontsize=10)
    axes[0, 1].set_title('Indoor Air Speed vs Wind Speed', fontsize=11, fontweight='bold')
    axes[0, 1].axhspan(0.15, 0.8, alpha=0.15, color='green', label='Comfort range')
    axes[0, 1].grid(alpha=0.3)
    axes[0, 1].legend()

    # Temperature
    axes[1, 0].plot(speeds, temps, '^-', color='#ED7D31', linewidth=2, markersize=8)
    axes[1, 0].set_xlabel('External Wind Speed (m/s)', fontsize=10)
    axes[1, 0].set_ylabel('Avg. Indoor Temperature (°C)', fontsize=10)
    axes[1, 0].set_title('Indoor Temperature vs Wind Speed', fontsize=11, fontweight='bold')
    axes[1, 0].axhspan(23, 28, alpha=0.15, color='orange', label='Comfort range')
    axes[1, 0].grid(alpha=0.3)
    axes[1, 0].legend()

    # Comfort fraction
    axes[1, 1].bar(speeds, comforts, width=0.6, color='#5B9BD5', edgecolor='black', linewidth=0.5)
    axes[1, 1].set_xlabel('External Wind Speed (m/s)', fontsize=10)
    axes[1, 1].set_ylabel('Comfort Zone Fraction (%)', fontsize=10)
    axes[1, 1].set_title('Thermal Comfort Coverage', fontsize=11, fontweight='bold')
    axes[1, 1].set_ylim(0, 100)
    axes[1, 1].grid(axis='y', alpha=0.3)

    plt.suptitle('Cross-Ventilation Performance Under Various Wind Conditions', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'cross_ventilation.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> cross_ventilation.png")


def plot_daylight_results(daylight_metrics, fig_dir):
    """Plot daylighting analysis results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Sample illuminance grid
    illum = np.array(daylight_metrics['illuminance_grid_sample'])
    im1 = axes[0, 0].imshow(illum, cmap='YlOrRd', aspect='auto', origin='lower',
                             extent=[0, 25, 0, 20])
    axes[0, 0].set_title('Illuminance Distribution (lux) - Sep 12:00', fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel('Length (m)')
    axes[0, 0].set_ylabel('Width (m)')
    plt.colorbar(im1, ax=axes[0, 0], label='Illuminance (lux)')

    # sDA grid
    sda = np.array(daylight_metrics['sda_grid'])
    im2 = axes[0, 1].imshow(sda, cmap='RdYlGn', aspect='auto', origin='lower',
                             extent=[0, 25, 0, 20], vmin=0, vmax=100)
    axes[0, 1].set_title('Spatial Daylight Autonomy sDA₃₀₀/₅₀ (%)', fontsize=11, fontweight='bold')
    axes[0, 1].set_xlabel('Length (m)')
    axes[0, 1].set_ylabel('Width (m)')
    plt.colorbar(im2, ax=axes[0, 1], label='sDA (%)')

    # UDI grid
    udi = np.array(daylight_metrics['udi_grid'])
    im3 = axes[1, 0].imshow(udi, cmap='RdYlGn', aspect='auto', origin='lower',
                             extent=[0, 25, 0, 20], vmin=0, vmax=100)
    axes[1, 0].set_title('Useful Daylight Illuminance UDI₁₀₀₋₃₀₀₀ (%)', fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel('Length (m)')
    axes[1, 0].set_ylabel('Width (m)')
    plt.colorbar(im3, ax=axes[1, 0], label='UDI (%)')

    # Summary metrics bar chart
    metrics_labels = ['sDA₃₀₀/₅₀', 'ASE₁₀₀₀/₂₅₀', 'Mean DA₃₀₀', 'Mean UDI']
    metrics_values = [
        daylight_metrics['sDA300_50'],
        daylight_metrics['ASE1000_250'],
        daylight_metrics['mean_DA300'],
        daylight_metrics['mean_UDI_100_3000'],
    ]
    colors = ['#70AD47', '#ED7D31', '#4472C4', '#FFC000']
    bars = axes[1, 1].bar(metrics_labels, metrics_values, color=colors, edgecolor='black', linewidth=0.5)
    for bar, val in zip(bars, metrics_values):
        axes[1, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f'{val:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    axes[1, 1].set_ylabel('Percentage (%)', fontsize=10)
    axes[1, 1].set_title('Daylighting Performance Metrics Summary', fontsize=11, fontweight='bold')
    axes[1, 1].set_ylim(0, 100)
    axes[1, 1].grid(axis='y', alpha=0.3)

    # LEED assessment
    leed_text = f"LEED v4.1 Daylight: {daylight_metrics['leed_daylight_points']} point(s)"
    axes[1, 1].text(0.5, 0.95, leed_text, transform=axes[1, 1].transAxes,
                    ha='center', fontsize=9, fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor='lightgreen' if daylight_metrics['leed_sda_pass'] else 'lightyellow', alpha=0.8))

    plt.suptitle('Daylighting Analysis Results (Radiance/Honeybee)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'daylight_results.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> daylight_results.png")


def plot_zeb_analysis(zeb_results, fig_dir):
    """Plot ZEB case study results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    baseline = zeb_results['baseline']

    # 1. Energy reduction waterfall
    improvements = zeb_results['improvements']
    labels = ['Baseline']
    values = [baseline['eui_kWh_m2']]
    for key, imp in improvements.items():
        labels.append(imp['description'].split('(')[0].strip())
        reduction = baseline['eui_kWh_m2'] * imp['energy_reduction_pct'] / 100
        values.append(-reduction)
    labels.append('Optimized')
    values.append(zeb_results['optimized_eui_kWh_m2'])

    cumulative = [values[0]]
    for v in values[1:-1]:
        cumulative.append(cumulative[-1] + v)
    cumulative.append(values[-1])

    colors_wf = ['#4472C4'] + ['#70AD47'] * (len(values) - 2) + ['#ED7D31']
    bottoms = [0] + [min(cumulative[i], cumulative[i] - values[i+1]) for i in range(len(values) - 2)] + [0]

    for i, (label, val) in enumerate(zip(labels, values)):
        if i == 0 or i == len(values) - 1:
            axes[0, 0].bar(i, abs(val), color=colors_wf[i], edgecolor='black', linewidth=0.5)
        else:
            axes[0, 0].bar(i, abs(val), bottom=cumulative[i], color=colors_wf[i],
                          edgecolor='black', linewidth=0.5)

    axes[0, 0].set_xticks(range(len(labels)))
    axes[0, 0].set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    axes[0, 0].set_ylabel('EUI (kWh/m²/yr)', fontsize=10)
    axes[0, 0].set_title('Energy Reduction Measures', fontsize=11, fontweight='bold')
    axes[0, 0].grid(axis='y', alpha=0.3)

    # 2. ZEB energy balance
    balance_labels = ['Energy\nDemand', 'PV\nGeneration']
    balance_values = [zeb_results['remaining_energy_kWh'] / 1000,
                      zeb_results['pv_generation_kWh'] / 1000]
    balance_colors = ['#ED7D31', '#70AD47']
    bars = axes[0, 1].bar(balance_labels, balance_values, color=balance_colors,
                          edgecolor='black', linewidth=1, width=0.5)
    for bar, val in zip(bars, balance_values):
        axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f'{val:.0f} MWh', ha='center', fontsize=11, fontweight='bold')
    axes[0, 1].set_ylabel('Energy (MWh/year)', fontsize=10)
    axes[0, 1].set_title(f'ZEB Energy Balance (Ratio: {zeb_results["zeb_ratio"]:.2f})',
                         fontsize=11, fontweight='bold')
    zeb_status = "✓ ZEB Achieved" if zeb_results['is_zeb'] else "△ Nearly ZEB"
    axes[0, 1].text(0.5, 0.95, zeb_status, transform=axes[0, 1].transAxes,
                    ha='center', fontsize=12, fontweight='bold',
                    color='green' if zeb_results['is_zeb'] else 'orange',
                    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    axes[0, 1].grid(axis='y', alpha=0.3)

    # 3. Monthly energy balance
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    monthly_demand = np.array(baseline['monthly_heating_kWh']) + np.array(baseline['monthly_cooling_kWh'])
    monthly_demand = monthly_demand * (1 - zeb_results['total_reduction_pct'] / 100)
    # Tokyo monthly solar variation
    solar_factors = [0.75, 0.80, 0.90, 1.00, 1.05, 0.85, 0.95, 1.05, 0.90, 0.85, 0.80, 0.70]
    monthly_pv = np.array(solar_factors) * zeb_results['pv_generation_kWh'] / 12

    x = np.arange(12)
    axes[1, 0].bar(x - 0.2, monthly_demand, 0.4, label='Energy Demand', color='#ED7D31', edgecolor='black', linewidth=0.5)
    axes[1, 0].bar(x + 0.2, monthly_pv, 0.4, label='PV Generation', color='#70AD47', edgecolor='black', linewidth=0.5)
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(months)
    axes[1, 0].set_ylabel('Energy (kWh)', fontsize=10)
    axes[1, 0].set_title('Monthly Energy Balance', fontsize=11, fontweight='bold')
    axes[1, 0].legend()
    axes[1, 0].grid(axis='y', alpha=0.3)

    # 4. Improvement breakdown pie
    imp_labels = [imp['description'].split('(')[0].strip() for imp in improvements.values()]
    imp_values = [imp['energy_reduction_pct'] for imp in improvements.values()]
    imp_colors = ['#4472C4', '#5B9BD5', '#FFC000', '#70AD47', '#ED7D31']
    wedges, texts, autotexts = axes[1, 1].pie(imp_values, labels=imp_labels, autopct='%1.1f%%',
                                               colors=imp_colors, startangle=90, textprops={'fontsize': 8})
    axes[1, 1].set_title('Energy Reduction Contribution', fontsize=11, fontweight='bold')

    plt.suptitle('Net Zero Energy Building (ZEB) Case Study Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'zeb_analysis.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> zeb_analysis.png")


def plot_integrated_dashboard(thermal, cfd, daylight, zeb, fig_dir):
    """Create integrated performance dashboard."""
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.35)

    # Title
    fig.suptitle('Integrated Building Environmental Performance Dashboard',
                 fontsize=16, fontweight='bold', y=0.98)

    # 1. Energy Summary (top-left)
    ax1 = fig.add_subplot(gs[0, 0:2])
    categories = ['Heating', 'Cooling', 'Lighting', 'Equipment', 'Fan']
    vals = [thermal['heating_eui'], thermal['cooling_eui'],
            thermal['annual_lighting_kWh'] / thermal['total_floor_area'],
            thermal['annual_equipment_kWh'] / thermal['total_floor_area'],
            thermal['annual_fan_kWh'] / thermal['total_floor_area']]
    colors = ['#ED7D31', '#4472C4', '#FFC000', '#70AD47', '#7030A0']
    ax1.bar(categories, vals, color=colors, edgecolor='black', linewidth=0.5)
    ax1.set_ylabel('EUI (kWh/m²/yr)')
    ax1.set_title(f'Energy Use Intensity (Total: {thermal["eui_kWh_m2"]} kWh/m²/yr)', fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)

    # 2. CFD Summary (top-right)
    ax2 = fig.add_subplot(gs[0, 2:4])
    vel = np.array(cfd['velocity_field'])
    im = ax2.imshow(vel, cmap='jet', aspect='auto', origin='lower', extent=[0, 25, 0, 20])
    ax2.set_title(f'Airflow (Avg: {cfd["avg_indoor_velocity_ms"]} m/s, ACH: {cfd["air_changes_per_hour"]})',
                  fontweight='bold')
    ax2.set_xlabel('Length (m)')
    ax2.set_ylabel('Width (m)')
    plt.colorbar(im, ax=ax2, label='m/s')

    # 3. Daylight metrics (middle-left)
    ax3 = fig.add_subplot(gs[1, 0:2])
    dl_labels = ['sDA₃₀₀/₅₀', 'ASE₁₀₀₀/₂₅₀', 'DA₃₀₀', 'UDI']
    dl_vals = [daylight['sDA300_50'], daylight['ASE1000_250'],
               daylight['mean_DA300'], daylight['mean_UDI_100_3000']]
    dl_colors = ['#70AD47', '#ED7D31', '#4472C4', '#FFC000']
    bars = ax3.bar(dl_labels, dl_vals, color=dl_colors, edgecolor='black', linewidth=0.5)
    for bar, val in zip(bars, dl_vals):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 f'{val:.1f}%', ha='center', fontsize=9, fontweight='bold')
    ax3.set_ylabel('Percentage (%)')
    ax3.set_title(f'Daylighting Metrics (LEED: {daylight["leed_daylight_points"]} pts)', fontweight='bold')
    ax3.set_ylim(0, 100)
    ax3.grid(axis='y', alpha=0.3)

    # 4. Daylight distribution (middle-right)
    ax4 = fig.add_subplot(gs[1, 2:4])
    sda_grid = np.array(daylight['sda_grid'])
    im4 = ax4.imshow(sda_grid, cmap='RdYlGn', aspect='auto', origin='lower',
                     extent=[0, 25, 0, 20], vmin=0, vmax=100)
    ax4.set_title('Spatial Daylight Autonomy Distribution', fontweight='bold')
    ax4.set_xlabel('Length (m)')
    ax4.set_ylabel('Width (m)')
    plt.colorbar(im4, ax=ax4, label='sDA (%)')

    # 5. ZEB Balance (bottom-left)
    ax5 = fig.add_subplot(gs[2, 0:2])
    zeb_labels = ['Demand\n(Optimized)', 'PV\nGeneration']
    zeb_vals = [zeb['remaining_energy_kWh'] / 1000, zeb['pv_generation_kWh'] / 1000]
    zeb_colors = ['#ED7D31', '#70AD47']
    bars = ax5.bar(zeb_labels, zeb_vals, color=zeb_colors, edgecolor='black', linewidth=1, width=0.5)
    for bar, val in zip(bars, zeb_vals):
        ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 f'{val:.0f} MWh', ha='center', fontsize=10, fontweight='bold')
    ax5.set_ylabel('Energy (MWh/yr)')
    status = "✓ ZEB" if zeb['is_zeb'] else f"Ratio: {zeb['zeb_ratio']:.2f}"
    ax5.set_title(f'ZEB Energy Balance ({status})', fontweight='bold')
    ax5.grid(axis='y', alpha=0.3)

    # 6. Performance scorecard (bottom-right)
    ax6 = fig.add_subplot(gs[2, 2:4])
    ax6.axis('off')
    scorecard = [
        ['Metric', 'Value', 'Target', 'Status'],
        ['EUI (kWh/m²/yr)', f'{thermal["eui_kWh_m2"]}', '≤200', '✓' if thermal['eui_kWh_m2'] <= 200 else '△'],
        ['sDA₃₀₀/₅₀ (%)', f'{daylight["sDA300_50"]:.1f}', '≥55', '✓' if daylight['sDA300_50'] >= 55 else '△'],
        ['ASE₁₀₀₀/₂₅₀ (%)', f'{daylight["ASE1000_250"]:.1f}', '≤10', '✓' if daylight['ASE1000_250'] <= 10 else '△'],
        ['ACH (nat. vent.)', f'{cfd["air_changes_per_hour"]}', '≥4', '✓' if cfd['air_changes_per_hour'] >= 4 else '△'],
        ['Comfort Zone (%)', f'{cfd["comfort_zone_fraction"]*100:.1f}', '≥80', '✓' if cfd['comfort_zone_fraction'] >= 0.8 else '△'],
        ['ZEB Ratio', f'{zeb["zeb_ratio"]:.2f}', '≥1.0', '✓' if zeb['is_zeb'] else '△'],
        ['PV Capacity (kW)', f'{zeb["pv_capacity_kW"]:.0f}', '-', '-'],
    ]

    table = ax6.table(cellText=scorecard, loc='center', cellLoc='center',
                      colWidths=[0.35, 0.2, 0.2, 0.1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    # Style header
    for j in range(4):
        table[0, j].set_facecolor('#4472C4')
        table[0, j].set_text_props(color='white', fontweight='bold')

    for i in range(1, len(scorecard)):
        status = scorecard[i][3]
        if status == '✓':
            table[i, 3].set_facecolor('#C6EFCE')
        elif status == '△':
            table[i, 3].set_facecolor('#FFEB9C')

    ax6.set_title('Performance Scorecard', fontweight='bold', fontsize=12, pad=20)

    plt.savefig(os.path.join(fig_dir, 'integrated_dashboard.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> integrated_dashboard.png")


def main():
    print("=" * 60)
    print("BIM-Integrated Environmental Performance Simulation System")
    print("=" * 60)

    fig_dir = create_figures_dir()

    # Step 1: IFC Model Conversion
    print("\n[1/6] IFC Model Conversion...")
    converter = IFCConverter()
    model = converter.create_reference_building("office")
    summary = converter.export_model_summary()
    print(f"  Building: {summary['building_name']}")
    print(f"  Location: {summary['location']} (Climate Zone {summary['climate_zone']})")
    print(f"  Total Floor Area: {summary['total_floor_area_m2']} m²")
    print(f"  Zones: {summary['num_zones']}")

    # Step 2: Thermal Simulation
    print("\n[2/6] Thermal Load Simulation (EnergyPlus)...")
    ep_params = converter.generate_energyplus_params()
    thermal_sim = ThermalSimulation(ep_params)
    thermal_results = thermal_sim.run_annual_simulation()
    print(f"  Annual EUI: {thermal_results['eui_kWh_m2']} kWh/m²/yr")
    print(f"  Heating: {thermal_results['annual_heating_kWh']} kWh | Cooling: {thermal_results['annual_cooling_kWh']} kWh")
    print(f"  Peak Heating: {thermal_results['peak_heating_kW']} kW | Peak Cooling: {thermal_results['peak_cooling_kW']} kW")

    # Step 3: CFD Simulation
    print("\n[3/6] CFD Natural Ventilation Simulation...")
    cfd_params = converter.generate_cfd_params()
    cfd_sim = CFDSimulation(cfd_params)
    cfd_results = cfd_sim.run_simulation()
    print(f"  Iterations: {cfd_results['iterations']} | Converged: {cfd_results['converged']}")
    print(f"  Avg Indoor Velocity: {cfd_results['avg_indoor_velocity_ms']} m/s")
    print(f"  ACH: {cfd_results['air_changes_per_hour']} | Comfort: {cfd_results['comfort_zone_fraction']*100:.1f}%")

    # Step 4: Cross-ventilation scenarios
    print("\n[4/6] Cross-Ventilation Scenario Analysis...")
    cfd_sim2 = CFDSimulation(cfd_params)
    scenarios = cfd_sim2.evaluate_cross_ventilation_scenarios()
    for s in scenarios:
        print(f"  Wind {s['wind_speed']}m/s -> ACH={s['ach']}, Comfort={s['comfort_fraction']*100:.0f}%")

    # Step 5: Daylighting Simulation
    print("\n[5/6] Daylighting Simulation (Radiance/Honeybee)...")
    rad_params = converter.generate_radiance_params()
    daylight_sim = DaylightSimulation(rad_params, ep_params['zones'][0])
    daylight_results = daylight_sim.calculate_annual_metrics()
    print(f"  sDA300/50: {daylight_results['sDA300_50']}%")
    print(f"  ASE1000/250: {daylight_results['ASE1000_250']}%")
    print(f"  Mean DA300: {daylight_results['mean_DA300']}%")
    print(f"  LEED Daylight: {daylight_results['leed_daylight_points']} point(s)")

    # Step 6: ZEB Case Study
    print("\n[6/6] ZEB Optimization Case Study...")
    thermal_sim2 = ThermalSimulation(ep_params)
    zeb_results = thermal_sim2.run_zeb_optimization()
    print(f"  Baseline EUI: {zeb_results['baseline']['eui_kWh_m2']} kWh/m²/yr")
    print(f"  Optimized EUI: {zeb_results['optimized_eui_kWh_m2']} kWh/m²/yr")
    print(f"  Total Reduction: {zeb_results['total_reduction_pct']}%")
    print(f"  PV Generation: {zeb_results['pv_generation_kWh']} kWh/yr")
    print(f"  ZEB Ratio: {zeb_results['zeb_ratio']} | ZEB Achieved: {zeb_results['is_zeb']}")

    # Generate all figures
    print("\n[Generating Figures]...")
    plot_system_architecture(fig_dir)
    plot_monthly_energy(thermal_results, fig_dir)
    plot_cfd_results(cfd_results, fig_dir)
    plot_cross_ventilation_scenarios(scenarios, fig_dir)
    plot_daylight_results(daylight_results, fig_dir)
    plot_zeb_analysis(zeb_results, fig_dir)
    plot_integrated_dashboard(thermal_results, cfd_results, daylight_results, zeb_results, fig_dir)

    # Save all results as JSON
    all_results = {
        'building_summary': summary,
        'thermal_results': {k: v for k, v in thermal_results.items()},
        'cfd_results': {k: v for k, v in cfd_results.items()
                       if k not in ['velocity_field', 'temperature_field', 'residuals']},
        'cross_ventilation_scenarios': scenarios,
        'daylight_results': {k: v for k, v in daylight_results.items()
                            if 'grid' not in k},
        'zeb_results': {k: v for k, v in zeb_results.items()
                       if k != 'baseline'},
        'zeb_baseline': {k: v for k, v in zeb_results['baseline'].items()
                        if 'monthly' not in k},
    }

    results_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'simulation_results.json')
    with open(results_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n  Results saved to: simulation_results.json")

    print("\n" + "=" * 60)
    print("All simulations completed successfully!")
    print("=" * 60)

    return all_results


if __name__ == "__main__":
    results = main()
