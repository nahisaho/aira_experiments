"""
Integrated Dashboard Visualization
Generates all figures for the BIM-integrated environmental performance simulation.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
from ifc_converter import IFCParser
from thermal_simulation import WeatherData, ThermalSimulation
from cfd_ventilation import VentilationCFD, run_multi_scenario
from daylight_simulation import DaylightSimulation
from zeb_analysis import ZEBAnalysis

plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.bbox'] = 'tight'

def main():
    print("=" * 60)
    print("BIM-Integrated Environmental Performance Simulation")
    print("=" * 60)
    
    # Step 1: IFC Model Conversion
    print("\n[1/6] IFC Model Conversion...")
    parser = IFCParser()
    parser.load_sample_building()
    zones = parser.to_energyplus_zones()
    config = parser.get_building_summary()
    config["SHGC"] = parser.windows["SHGC"]
    print(f"  Building: {config['total_floor_area']} m², {config['num_zones']} zones, {config['num_floors']} floors")
    
    # Step 2: Weather Data
    print("\n[2/6] Generating Weather Data (Tokyo)...")
    weather_gen = WeatherData("Tokyo")
    weather = weather_gen.generate()
    print(f"  Temperature range: {weather['temperature'].min():.1f}°C to {weather['temperature'].max():.1f}°C")
    
    # Step 3: Thermal Simulation
    print("\n[3/6] Running Thermal Load Simulation...")
    thermal_sim = ThermalSimulation(zones, weather, config)
    thermal_results = thermal_sim.run()
    print(f"  Annual Heating: {thermal_results['annual_heating_kWh']:.0f} kWh ({thermal_results['EUI_heating']:.1f} kWh/m²)")
    print(f"  Annual Cooling: {thermal_results['annual_cooling_kWh']:.0f} kWh ({thermal_results['EUI_cooling']:.1f} kWh/m²)")
    print(f"  EUI Total: {thermal_results['EUI_total']:.1f} kWh/m²/yr")
    
    # Step 4: CFD Ventilation
    print("\n[4/6] Running CFD Ventilation Analysis...")
    cfd_results = run_multi_scenario()
    for name, data in cfd_results.items():
        m = data["metrics"]
        print(f"  {name}: ACH={m['ACH']:.1f}, v_avg={m['avg_velocity_ms']:.3f} m/s, comfort={m['comfort_ratio']:.1%}")
    
    # Step 5: Daylight Simulation
    print("\n[5/6] Running Daylight Simulation...")
    daylight_sim = DaylightSimulation()
    daylight_sim.configure_room(glazing_vlt=0.50)
    df = daylight_sim.calculate_daylight_factor()
    annual_dl = daylight_sim.calculate_annual_metrics(weather["ghi"])
    print(f"  Mean DF: {annual_dl['mean_DF']:.2f}%")
    print(f"  sDA300/50%: {annual_dl['sDA300_50']:.1f}%")
    print(f"  UDI(100-2000): {annual_dl['UDI_100_2000']:.1f}%")
    
    # Step 6: ZEB Analysis
    print("\n[6/6] ZEB Analysis...")
    zeb = ZEBAnalysis(thermal_results, config)
    energy = zeb.calculate_primary_energy()
    pv = zeb.design_pv_system(weather["ghi"])
    zeb_result = zeb.evaluate_zeb_compliance()
    scenarios = zeb.optimization_scenarios(weather["ghi"])
    print(f"  Total Consumption: {zeb_result['total_consumption_kWh']:.0f} kWh")
    print(f"  PV Generation: {zeb_result['total_generation_kWh']:.0f} kWh")
    print(f"  ZEB Ratio: {zeb_result['zeb_ratio']:.2f}")
    print(f"  Classification: {zeb_result['zeb_classification']}")
    
    # Generate Figures
    print("\n" + "=" * 60)
    print("Generating Figures...")
    os.makedirs("figures", exist_ok=True)
    
    # Figure 1: System Architecture
    generate_architecture_diagram()
    
    # Figure 2: Monthly Energy Balance
    generate_monthly_energy(thermal_results)
    
    # Figure 3: CFD Velocity Fields
    generate_cfd_plots(cfd_results)
    
    # Figure 4: Daylight Factor Map
    generate_daylight_plots(daylight_sim, annual_dl)
    
    # Figure 5: ZEB Energy Balance
    generate_zeb_plots(zeb_result, energy, scenarios)
    
    # Figure 6: Integrated Dashboard
    generate_dashboard(thermal_results, cfd_results, annual_dl, zeb_result)
    
    # Figure 7: Weather Data
    generate_weather_plot(weather)
    
    # Figure 8: Parametric Study
    generate_parametric_study(daylight_sim, weather)
    
    print("\nAll figures generated successfully!")
    print("=" * 60)
    
    return {
        "thermal": thermal_results,
        "cfd": cfd_results,
        "daylight": annual_dl,
        "zeb": zeb_result,
        "energy": energy,
        "scenarios": scenarios,
        "config": config,
    }


def generate_architecture_diagram():
    """Figure 1: System architecture diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis('off')
    ax.set_title('BIM-Integrated Environmental Performance Simulation System Architecture', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Boxes
    boxes = [
        (1, 6.5, 3, 1, 'IFC/BIM Model\n(Building Data)', '#4ECDC4'),
        (5, 6.5, 3, 1, 'IFC Parser &\nModel Converter', '#45B7D1'),
        (9, 6.5, 2.5, 1, 'Geometry\nValidator', '#96CEB4'),
        (0.5, 4.5, 2.2, 1, 'EnergyPlus\nThermal Sim', '#FF6B6B'),
        (3.2, 4.5, 2.2, 1, 'CFD\nVentilation', '#C06C84'),
        (5.9, 4.5, 2.2, 1, 'Radiance\nDaylight', '#F8B500'),
        (8.6, 4.5, 2.2, 1, 'Renewable\nEnergy', '#6C5CE7'),
        (3, 2.5, 6, 1, 'Integration Engine & Data Aggregation', '#2D3436'),
        (3, 0.5, 6, 1, 'Unified Dashboard & ZEB Evaluation', '#0984E3'),
    ]
    
    for x, y, w, h, label, color in boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                              facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.85)
        ax.add_patch(rect)
        text_color = 'white' if color in ['#2D3436', '#0984E3', '#6C5CE7', '#C06C84'] else 'black'
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=9, fontweight='bold', color=text_color)
    
    # Arrows
    arrows = [
        (2.5, 6.5, 5, 6.5),
        (8, 6.5, 9, 6.5),
        (6.5, 6.5, 1.6, 5.5),
        (6.5, 6.5, 4.3, 5.5),
        (6.5, 6.5, 7.0, 5.5),
        (6.5, 6.5, 9.7, 5.5),
        (1.6, 4.5, 6, 3.5),
        (4.3, 4.5, 6, 3.5),
        (7.0, 4.5, 6, 3.5),
        (9.7, 4.5, 6, 3.5),
        (6, 2.5, 6, 1.5),
    ]
    
    for x1, y1, x2, y2 in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    
    plt.savefig('figures/fig1_system_architecture.png', dpi=150)
    plt.close()
    print("  [✓] Figure 1: System Architecture")


def generate_monthly_energy(thermal):
    """Figure 2: Monthly energy load profiles."""
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Monthly heating and cooling
    ax = axes[0]
    x = np.arange(12)
    w = 0.35
    ax.bar(x - w/2, thermal['monthly_heating'], w, label='Heating', color='#FF6B6B', alpha=0.8)
    ax.bar(x + w/2, thermal['monthly_cooling'], w, label='Cooling', color='#45B7D1', alpha=0.8)
    ax.set_xlabel('Month')
    ax.set_ylabel('Energy (kWh)')
    ax.set_title('Monthly Heating and Cooling Loads')
    ax.set_xticks(x)
    ax.set_xticklabels(months, rotation=45)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Energy breakdown pie
    ax = axes[1]
    labels = ['Heating', 'Cooling', 'Lighting', 'Equipment']
    values = [
        thermal['annual_heating_kWh'],
        thermal['annual_cooling_kWh'],
        thermal['annual_lighting_kWh'],
        thermal['annual_equipment_kWh'],
    ]
    colors = ['#FF6B6B', '#45B7D1', '#F8B500', '#96CEB4']
    ax.pie(values, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90, 
           pctdistance=0.85, wedgeprops=dict(width=0.4))
    ax.set_title('Annual Energy Breakdown')
    
    plt.tight_layout()
    plt.savefig('figures/fig2_monthly_energy.png', dpi=150)
    plt.close()
    print("  [✓] Figure 2: Monthly Energy Loads")


def generate_cfd_plots(cfd_results):
    """Figure 3: CFD velocity field visualizations."""
    scenarios = ["Baseline", "Large_Opening", "Offset_Opening"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for idx, name in enumerate(scenarios):
        ax = axes[idx]
        vel = cfd_results[name]["velocity_field"]
        im = ax.imshow(vel, cmap='jet', aspect='auto', origin='lower',
                       extent=[0, 10, 0, 8], vmin=0, vmax=vel.max())
        ax.set_title(f'{name}\nACH={cfd_results[name]["metrics"]["ACH"]:.1f}')
        ax.set_xlabel('Width (m)')
        ax.set_ylabel('Depth (m)')
        plt.colorbar(im, ax=ax, label='Velocity (m/s)', shrink=0.8)
    
    plt.suptitle('Cross-Ventilation CFD Analysis: Velocity Fields', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/fig3_cfd_velocity.png', dpi=150)
    plt.close()
    print("  [✓] Figure 3: CFD Velocity Fields")
    
    # CFD comparison bar chart
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    names = list(cfd_results.keys())
    achs = [cfd_results[n]["metrics"]["ACH"] for n in names]
    comforts = [cfd_results[n]["metrics"]["comfort_ratio"] * 100 for n in names]
    
    colors = ['#4ECDC4', '#FF6B6B', '#45B7D1', '#F8B500', '#6C5CE7']
    
    ax = axes[0]
    ax.barh(names, achs, color=colors, alpha=0.8)
    ax.set_xlabel('Air Changes per Hour (ACH)')
    ax.set_title('Ventilation Rate Comparison')
    ax.grid(axis='x', alpha=0.3)
    
    ax = axes[1]
    ax.barh(names, comforts, color=colors, alpha=0.8)
    ax.set_xlabel('Comfort Zone Coverage (%)')
    ax.set_title('Thermal Comfort Coverage')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/fig4_cfd_comparison.png', dpi=150)
    plt.close()
    print("  [✓] Figure 4: CFD Comparison")


def generate_daylight_plots(sim, metrics):
    """Figure 5: Daylight analysis results."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Daylight Factor map
    ax = axes[0]
    im = ax.imshow(sim.daylight_factor, cmap='YlOrRd', aspect='auto', origin='lower',
                   extent=[0, sim.W, 0, sim.D])
    ax.set_title(f'Daylight Factor (%)\nMean={metrics["mean_DF"]:.2f}%')
    ax.set_xlabel('Width (m)')
    ax.set_ylabel('Depth (m)')
    plt.colorbar(im, ax=ax, label='DF (%)', shrink=0.8)
    
    # Daylight Autonomy map
    ax = axes[1]
    im = ax.imshow(metrics["da_grid"], cmap='RdYlGn', aspect='auto', origin='lower',
                   extent=[0, sim.W, 0, sim.D], vmin=0, vmax=100)
    ax.set_title(f'Daylight Autonomy (%)\nMean={metrics["mean_DA"]:.1f}%')
    ax.set_xlabel('Width (m)')
    ax.set_ylabel('Depth (m)')
    plt.colorbar(im, ax=ax, label='DA (%)', shrink=0.8)
    
    # Metrics summary
    ax = axes[2]
    ax.axis('off')
    metric_data = [
        ['Metric', 'Value', 'Target'],
        ['Mean DF', f'{metrics["mean_DF"]:.2f}%', '≥ 2.0%'],
        ['sDA₃₀₀/₅₀%', f'{metrics["sDA300_50"]:.1f}%', '≥ 55%'],
        ['ASE₁₀₀₀/₂₅₀', f'{metrics["ASE1000_250"]:.1f}%', '≤ 10%'],
        ['UDI₁₀₀₋₂₀₀₀', f'{metrics["UDI_100_2000"]:.1f}%', '≥ 60%'],
    ]
    table = ax.table(cellText=metric_data[1:], colLabels=metric_data[0],
                     cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.8)
    for i in range(len(metric_data[0])):
        table[0, i].set_facecolor('#4ECDC4')
        table[0, i].set_text_props(fontweight='bold', color='white')
    ax.set_title('Daylight Performance Metrics', fontweight='bold')
    
    plt.suptitle('Daylight Simulation Results (Radiance/Honeybee)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/fig5_daylight.png', dpi=150)
    plt.close()
    print("  [✓] Figure 5: Daylight Analysis")


def generate_zeb_plots(zeb_result, energy, scenarios):
    """Figure 6: ZEB analysis results."""
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Monthly energy balance
    ax = axes[0, 0]
    x = np.arange(12)
    ax.bar(x, zeb_result['monthly_consumption'], 0.4, label='Consumption', color='#FF6B6B', alpha=0.8)
    ax.bar(x + 0.4, zeb_result['monthly_generation'], 0.4, label='PV Generation', color='#4ECDC4', alpha=0.8)
    ax.plot(x + 0.2, zeb_result['monthly_balance'], 'k-o', label='Net Balance', markersize=4)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Month')
    ax.set_ylabel('Energy (kWh)')
    ax.set_title('Monthly Energy Balance')
    ax.set_xticks(x + 0.2)
    ax.set_xticklabels(months, rotation=45)
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)
    
    # Energy breakdown waterfall
    ax = axes[0, 1]
    categories = ['Cooling', 'Heating', 'Lighting', 'Equipment', 'Auxiliary', 'DHW']
    values = [energy['cooling_kWh'], energy['heating_kWh'], energy['lighting_kWh'],
              energy['equipment_kWh'], energy['auxiliary_kWh'], energy['dhw_kWh']]
    colors = ['#45B7D1', '#FF6B6B', '#F8B500', '#96CEB4', '#C06C84', '#6C5CE7']
    ax.barh(categories, values, color=colors, alpha=0.85)
    ax.set_xlabel('Annual Electricity (kWh)')
    ax.set_title('Electricity Consumption Breakdown')
    ax.grid(axis='x', alpha=0.3)
    for i, v in enumerate(values):
        ax.text(v + 200, i, f'{v:.0f}', va='center', fontsize=9)
    
    # Scenario comparison
    ax = axes[1, 0]
    sc_names = [s['scenario'] for s in scenarios]
    sc_ratios = [s['zeb_ratio'] for s in scenarios]
    bar_colors = ['#FF6B6B' if r < 0.75 else '#F8B500' if r < 1.0 else '#4ECDC4' for r in sc_ratios]
    ax.barh(sc_names, sc_ratios, color=bar_colors, alpha=0.85)
    ax.axvline(x=1.0, color='green', linestyle='--', label='ZEB Threshold', linewidth=2)
    ax.axvline(x=0.75, color='orange', linestyle='--', label='Nearly ZEB', linewidth=1)
    ax.set_xlabel('ZEB Ratio (Generation/Consumption)')
    ax.set_title('ZEB Compliance by Scenario')
    ax.legend(fontsize=8)
    ax.grid(axis='x', alpha=0.3)
    
    # ZEB summary gauge
    ax = axes[1, 1]
    ax.axis('off')
    ratio = zeb_result['zeb_ratio']
    classification = zeb_result['zeb_classification']
    
    summary = [
        f"ZEB Classification: {classification}",
        f"",
        f"Total Consumption: {zeb_result['total_consumption_kWh']:,.0f} kWh/yr",
        f"PV Generation: {zeb_result['total_generation_kWh']:,.0f} kWh/yr",
        f"Net Energy: {zeb_result['net_energy_kWh']:,.0f} kWh/yr",
        f"ZEB Ratio: {ratio:.2f}",
        f"Net EUI: {zeb_result['EUI_net']:.1f} kWh/m²/yr",
        f"",
        f"PV Area: {560:.0f} m²",
        f"PV Capacity: {560*0.22:.0f} kWp",
    ]
    
    bg_color = '#4ECDC4' if ratio >= 1.0 else '#F8B500' if ratio >= 0.75 else '#FF6B6B'
    rect = FancyBboxPatch((0.05, 0.05), 0.9, 0.9, boxstyle="round,pad=0.05",
                          facecolor=bg_color, alpha=0.2, edgecolor=bg_color, linewidth=2)
    ax.add_patch(rect)
    ax.text(0.5, 0.95, 'ZEB Evaluation Summary', ha='center', va='top',
            fontsize=14, fontweight='bold', transform=ax.transAxes)
    for i, line in enumerate(summary):
        weight = 'bold' if i == 0 else 'normal'
        ax.text(0.5, 0.82 - i * 0.08, line, ha='center', va='top',
                fontsize=11, fontweight=weight, transform=ax.transAxes)
    
    plt.suptitle('Net Zero Energy Building (ZEB) Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/fig6_zeb_analysis.png', dpi=150)
    plt.close()
    print("  [✓] Figure 6: ZEB Analysis")


def generate_dashboard(thermal, cfd, daylight, zeb):
    """Figure 7: Integrated performance dashboard."""
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(3, 4, hspace=0.4, wspace=0.35)
    
    fig.suptitle('Integrated Building Performance Dashboard', fontsize=16, fontweight='bold', y=0.98)
    
    # Panel 1: EUI summary
    ax = fig.add_subplot(gs[0, 0])
    ax.axis('off')
    eui_data = [
        ('Heating EUI', f"{thermal['EUI_heating']:.1f}", 'kWh/m²'),
        ('Cooling EUI', f"{thermal['EUI_cooling']:.1f}", 'kWh/m²'),
        ('Total EUI', f"{thermal['EUI_total']:.1f}", 'kWh/m²'),
    ]
    ax.text(0.5, 0.95, 'Energy Use Intensity', ha='center', va='top', fontsize=12, fontweight='bold',
            transform=ax.transAxes)
    for i, (label, val, unit) in enumerate(eui_data):
        ax.text(0.5, 0.7 - i*0.25, f"{label}\n{val} {unit}", ha='center', va='top',
                fontsize=11, transform=ax.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # Panel 2: Monthly loads
    ax = fig.add_subplot(gs[0, 1:3])
    months = np.arange(12)
    ax.fill_between(months, thermal['monthly_heating'], alpha=0.4, color='#FF6B6B', label='Heating')
    ax.fill_between(months, thermal['monthly_cooling'], alpha=0.4, color='#45B7D1', label='Cooling')
    ax.set_title('Monthly Thermal Loads', fontsize=10)
    ax.set_xticks(months)
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    ax.legend(fontsize=8)
    ax.set_ylabel('kWh')
    
    # Panel 3: ZEB status
    ax = fig.add_subplot(gs[0, 3])
    ax.axis('off')
    ratio = zeb['zeb_ratio']
    color = '#4ECDC4' if ratio >= 1.0 else '#F8B500' if ratio >= 0.75 else '#FF6B6B'
    circle = plt.Circle((0.5, 0.5), 0.35, color=color, alpha=0.3, transform=ax.transAxes)
    ax.add_patch(circle)
    ax.text(0.5, 0.55, f"{ratio:.0%}", ha='center', va='center', fontsize=24, fontweight='bold',
            transform=ax.transAxes, color=color)
    ax.text(0.5, 0.35, 'ZEB Ratio', ha='center', va='center', fontsize=10, transform=ax.transAxes)
    ax.text(0.5, 0.9, zeb['zeb_classification'], ha='center', va='top', fontsize=11, fontweight='bold',
            transform=ax.transAxes)
    
    # Panel 4: CFD summary
    ax = fig.add_subplot(gs[1, 0:2])
    baseline = cfd["Baseline"]
    vel = baseline["velocity_field"]
    im = ax.imshow(vel, cmap='jet', aspect='auto', origin='lower', extent=[0, 10, 0, 8])
    ax.set_title(f'Cross-Ventilation (Baseline)\nACH={baseline["metrics"]["ACH"]:.1f}', fontsize=10)
    ax.set_xlabel('Width (m)')
    ax.set_ylabel('Depth (m)')
    plt.colorbar(im, ax=ax, label='v (m/s)', shrink=0.7)
    
    # Panel 5: Daylight
    ax = fig.add_subplot(gs[1, 2:4])
    ax.axis('off')
    dl_items = [
        ('Mean DF', f"{daylight['mean_DF']:.2f}%", '≥2.0%'),
        ('sDA₃₀₀', f"{daylight['sDA300_50']:.1f}%", '≥55%'),
        ('ASE₁₀₀₀', f"{daylight['ASE1000_250']:.1f}%", '≤10%'),
        ('UDI', f"{daylight['UDI_100_2000']:.1f}%", '≥60%'),
    ]
    ax.text(0.5, 0.95, 'Daylight Performance', ha='center', va='top', fontsize=12, fontweight='bold',
            transform=ax.transAxes)
    for i, (name, val, target) in enumerate(dl_items):
        y = 0.75 - i * 0.2
        ax.text(0.2, y, name, ha='left', va='center', fontsize=11, transform=ax.transAxes)
        ax.text(0.55, y, val, ha='center', va='center', fontsize=11, fontweight='bold', transform=ax.transAxes)
        ax.text(0.8, y, f'Target: {target}', ha='center', va='center', fontsize=9, color='gray',
                transform=ax.transAxes)
    
    # Panel 6: ZEB balance
    ax = fig.add_subplot(gs[2, :])
    months_lbl = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    x = np.arange(12)
    ax.bar(x - 0.2, zeb['monthly_consumption'], 0.4, label='Consumption', color='#FF6B6B', alpha=0.7)
    ax.bar(x + 0.2, zeb['monthly_generation'], 0.4, label='PV Generation', color='#4ECDC4', alpha=0.7)
    ax.plot(x, zeb['monthly_balance'], 'k-o', label='Net Balance', markersize=4)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(months_lbl)
    ax.set_ylabel('Energy (kWh)')
    ax.set_title('Monthly Energy Balance (ZEB Evaluation)', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)
    
    plt.savefig('figures/fig7_dashboard.png', dpi=150)
    plt.close()
    print("  [✓] Figure 7: Integrated Dashboard")


def generate_weather_plot(weather):
    """Figure 8: Weather data visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    
    hours = np.arange(8760) / 24
    
    # Temperature
    ax = axes[0, 0]
    # Daily averages
    daily_temp = weather['temperature'].reshape(-1, 24).mean(axis=1)
    ax.plot(np.arange(365), daily_temp, color='#FF6B6B', alpha=0.8, linewidth=0.5)
    # Moving average
    window = 14
    ma = np.convolve(daily_temp, np.ones(window)/window, mode='valid')
    ax.plot(np.arange(len(ma)) + window//2, ma, color='darkred', linewidth=2)
    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title('Outdoor Air Temperature (Tokyo)')
    ax.grid(alpha=0.3)
    
    # Solar radiation
    ax = axes[0, 1]
    daily_ghi = weather['ghi'].reshape(-1, 24).sum(axis=1) / 1000
    ax.fill_between(np.arange(365), daily_ghi, alpha=0.4, color='#F8B500')
    ma_ghi = np.convolve(daily_ghi, np.ones(window)/window, mode='valid')
    ax.plot(np.arange(len(ma_ghi)) + window//2, ma_ghi, color='darkorange', linewidth=2)
    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Daily GHI (kWh/m²)')
    ax.set_title('Global Horizontal Irradiance')
    ax.grid(alpha=0.3)
    
    # Temperature histogram
    ax = axes[1, 0]
    ax.hist(weather['temperature'], bins=50, color='#45B7D1', alpha=0.7, edgecolor='white')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Frequency (hours)')
    ax.set_title('Temperature Distribution')
    ax.axvline(x=20, color='green', linestyle='--', label='Heating SP', alpha=0.7)
    ax.axvline(x=26, color='red', linestyle='--', label='Cooling SP', alpha=0.7)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    
    # Wind rose (simplified)
    ax = axes[1, 1]
    ax.hist(weather['wind_speed'], bins=40, color='#96CEB4', alpha=0.7, edgecolor='white')
    ax.set_xlabel('Wind Speed (m/s)')
    ax.set_ylabel('Frequency (hours)')
    ax.set_title('Wind Speed Distribution')
    ax.grid(alpha=0.3)
    
    plt.suptitle('Weather Data: Tokyo, Japan (TMY)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/fig8_weather.png', dpi=150)
    plt.close()
    print("  [✓] Figure 8: Weather Data")


def generate_parametric_study(daylight_sim, weather):
    """Figure 9: Parametric glazing study."""
    results = daylight_sim.parametric_glazing_study(weather["ghi"])
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    vlts = [r['VLT'] for r in results]
    sdas = [r['sDA'] for r in results]
    ases = [r['ASE'] for r in results]
    udis = [r['UDI'] for r in results]
    dfs = [r['mean_DF'] for r in results]
    
    ax = axes[0]
    ax.plot(vlts, sdas, 'o-', color='#4ECDC4', linewidth=2, markersize=8, label='sDA₃₀₀/₅₀%')
    ax.axhline(y=55, color='green', linestyle='--', alpha=0.5, label='Target (55%)')
    ax.set_xlabel('Visible Light Transmittance (VLT)')
    ax.set_ylabel('sDA (%)')
    ax.set_title('Spatial Daylight Autonomy')
    ax.legend()
    ax.grid(alpha=0.3)
    
    ax = axes[1]
    ax.plot(vlts, ases, 's-', color='#FF6B6B', linewidth=2, markersize=8, label='ASE₁₀₀₀/₂₅₀')
    ax.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='Limit (10%)')
    ax.set_xlabel('Visible Light Transmittance (VLT)')
    ax.set_ylabel('ASE (%)')
    ax.set_title('Annual Sunlight Exposure')
    ax.legend()
    ax.grid(alpha=0.3)
    
    ax = axes[2]
    ax.plot(vlts, dfs, 'D-', color='#F8B500', linewidth=2, markersize=8, label='Mean DF')
    ax.axhline(y=2.0, color='green', linestyle='--', alpha=0.5, label='Target (2.0%)')
    ax.set_xlabel('Visible Light Transmittance (VLT)')
    ax.set_ylabel('Daylight Factor (%)')
    ax.set_title('Mean Daylight Factor')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.suptitle('Parametric Study: Impact of Glazing VLT on Daylight Performance', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('figures/fig9_parametric.png', dpi=150)
    plt.close()
    print("  [✓] Figure 9: Parametric Study")


if __name__ == "__main__":
    all_results = main()
