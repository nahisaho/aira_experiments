"""
Visualization module — generates all figures for the UHI report.
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.urban_canopy_model import TOKYO_MORPHOLOGY, compute_local_climate_zone, UrbanCanopyModel
from src.anthropogenic_heat import TOKYO_HEAT_PROFILES, project_anthropogenic_heat_2050
from src.cooling_effects import MITIGATION_SCENARIOS
from src.wrf_ucm_coupling import OfflineCouplingEngine
from src.wbgt_risk import WBGTCalculator, HeatStrokeRiskAssessor, PopulationExposure

COLORS = {'marunouchi': '#d62728', 'shinjuku': '#ff7f0e',
          'residential_23ku': '#2ca02c', 'suburban': '#1f77b4'}
DISTRICTS = ["marunouchi", "shinjuku", "residential_23ku", "suburban"]


def generate_forcing(climate="current"):
    hours = np.arange(24)
    T_base = 28.0
    if climate == "rcp45_2050": T_base += 1.8
    elif climate == "rcp85_2050": T_base += 3.2
    T_air = T_base + 5.0 * np.sin(2 * np.pi * (hours - 6) / 24)
    T_air_K = T_air + 273.15
    solar_max = 850
    sw_down = np.maximum(0, solar_max * np.sin(np.pi * (hours - 5) / 14))
    sw_down[hours < 5] = 0; sw_down[hours > 19] = 0
    sigma = 5.67e-8
    lw_down = 0.85 * sigma * T_air_K**4
    u_star = 0.3 + 0.15 * np.sin(2 * np.pi * (hours - 14) / 24)
    q_air = 0.014 + 0.003 * np.sin(2 * np.pi * (hours + 3) / 24)
    rh = 65 + 15 * np.cos(2 * np.pi * (hours - 6) / 24)
    wind_10m = 2.5 + 1.0 * np.sin(2 * np.pi * (hours - 14) / 24)
    return {"sw_down": sw_down, "lw_down": lw_down, "T_air": T_air_K,
            "T_air_C": T_air, "u_star": u_star, "q_air": q_air,
            "rh": rh, "wind_10m": wind_10m, "hours": hours}


def run_sim(district, forcing, cooling=None, heat_profile=None):
    morph = TOKYO_MORPHOLOGY[district]
    heat = heat_profile or TOKYO_HEAT_PROFILES[district]
    ucm = UrbanCanopyModel(morph)
    engine = OfflineCouplingEngine(ucm, heat, cooling)
    return engine.run_diurnal_cycle(forcing)


def fig1_energy_balance():
    """Fig 1: Diurnal energy balance for Marunouchi."""
    forcing = generate_forcing()
    results = run_sim("marunouchi", forcing)
    hours = np.arange(24)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.plot(hours, results['Q_star'], 'k-', lw=2, label='Q* (Net Radiation)')
    ax1.plot(hours, results['QH'], 'r-', lw=2, label='QH (Sensible)')
    ax1.plot(hours, results['QE'], 'b-', lw=2, label='QE (Latent)')
    ax1.plot(hours, results['dQS'], 'g--', lw=2, label='ΔQS (Storage)')
    ax1.plot(hours, results['QF'], 'm:', lw=2, label='QF (Anthropogenic)')
    ax1.axhline(y=0, color='gray', ls='-', alpha=0.3)
    ax1.set_ylabel('Energy Flux [W/m²]')
    ax1.set_title('Diurnal Energy Balance — Marunouchi (August)')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3)

    ax2.stackplot(hours, results['QF_traffic'], results['QF_building'],
                  results['QF_industry'],
                  labels=['Traffic', 'Building HVAC', 'Industry'],
                  colors=['#ff9999', '#ff6666', '#cc3333'], alpha=0.8)
    ax2.set_xlabel('Hour (JST)')
    ax2.set_ylabel('Anthropogenic Heat [W/m²]')
    ax2.set_title('Anthropogenic Heat Flux Components')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 23)

    plt.tight_layout()
    plt.savefig('figures/fig1_energy_balance.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → figures/fig1_energy_balance.png")


def fig2_uhi_intensity():
    """Fig 2: UHI intensity comparison across districts."""
    forcing = generate_forcing()
    hours = np.arange(24)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for d in DISTRICTS:
        results = run_sim(d, forcing)
        ax1.plot(hours, results['UHI_intensity'], color=COLORS[d], lw=2, label=d)

    ax1.axhline(y=0, color='gray', ls='-', alpha=0.3)
    ax1.set_xlabel('Hour (JST)')
    ax1.set_ylabel('UHI Intensity [°C]')
    ax1.set_title('Diurnal UHI Intensity — Current Climate')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Bar chart of peak UHI by district
    peak_uhis = []
    for d in DISTRICTS:
        results = run_sim(d, forcing)
        peak_uhis.append(np.max(results['UHI_intensity']))

    bars = ax2.bar(DISTRICTS, peak_uhis, color=[COLORS[d] for d in DISTRICTS])
    ax2.set_ylabel('Peak UHI Intensity [°C]')
    ax2.set_title('Peak UHI Intensity by District')
    ax2.grid(True, alpha=0.3, axis='y')
    for bar, val in zip(bars, peak_uhis):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{val:.1f}°C', ha='center', fontsize=10)

    plt.tight_layout()
    plt.savefig('figures/fig2_uhi_intensity.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → figures/fig2_uhi_intensity.png")


def fig3_mitigation_scenarios():
    """Fig 3: Cooling effect of mitigation strategies."""
    forcing = generate_forcing()
    hours = np.arange(24)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    scenario_colors = {'baseline': '#d62728', 'moderate_greening': '#ff7f0e',
                       'aggressive_mitigation': '#2ca02c'}

    for idx, d in enumerate(DISTRICTS):
        ax = axes[idx // 2][idx % 2]
        for sname, cooling in MITIGATION_SCENARIOS.items():
            results = run_sim(d, forcing, cooling)
            ax.plot(hours, results['UHI_intensity'], color=scenario_colors[sname],
                    lw=2, label=sname.replace('_', ' ').title())
        ax.axhline(y=0, color='gray', ls='-', alpha=0.3)
        ax.set_title(f'{d.replace("_", " ").title()}')
        ax.set_xlabel('Hour (JST)')
        ax.set_ylabel('UHI Intensity [°C]')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('UHI Intensity Under Mitigation Scenarios', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('figures/fig3_mitigation_scenarios.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → figures/fig3_mitigation_scenarios.png")


def fig4_2050_projection():
    """Fig 4: 2050 temperature projections."""
    scenarios = ["current", "rcp45_2050", "rcp85_2050"]
    hours = np.arange(24)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    scenario_colors = {"current": '#1f77b4', "rcp45_2050": '#ff7f0e', "rcp85_2050": '#d62728'}

    for idx, d in enumerate(DISTRICTS):
        ax = axes[idx // 2][idx % 2]
        for climate in scenarios:
            forcing = generate_forcing(climate)
            heat = TOKYO_HEAT_PROFILES[d]
            if "2050" in climate:
                heat = project_anthropogenic_heat_2050(heat)
            results = run_sim(d, forcing, heat_profile=heat)
            T_canyon_C = results['T_canyon'] - 273.15
            label = climate.replace('_', ' ').upper()
            ax.plot(hours, T_canyon_C, color=scenario_colors[climate], lw=2, label=label)

        ax.set_title(f'{d.replace("_", " ").title()}')
        ax.set_xlabel('Hour (JST)')
        ax.set_ylabel('Canyon Air Temperature [°C]')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Canyon Temperature: Current vs 2050 Projections', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('figures/fig4_2050_projection.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → figures/fig4_2050_projection.png")


def fig5_wbgt_risk():
    """Fig 5: WBGT and heat stroke risk."""
    wbgt_calc = WBGTCalculator()
    hours = np.arange(24)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for idx, (climate, label) in enumerate([("current", "Current"),
                                             ("rcp85_2050", "RCP8.5 2050")]):
        forcing = generate_forcing(climate)
        ax_wbgt = axes[0][idx]
        ax_risk = axes[1][idx]

        for d in DISTRICTS:
            results = run_sim(d, forcing)
            uhi = results['UHI_intensity']
            T_eff = forcing['T_air_C'] + uhi
            wbgt = np.array([wbgt_calc.compute_wbgt_outdoor(
                T_eff[h], forcing['rh'][h], forcing['wind_10m'][h],
                forcing['sw_down'][h]) for h in range(24)])
            ax_wbgt.plot(hours, wbgt, color=COLORS[d], lw=2, label=d)

        # Risk thresholds
        ax_wbgt.axhspan(31, 40, alpha=0.1, color='red', label='Danger (≥31)')
        ax_wbgt.axhspan(28, 31, alpha=0.1, color='orange', label='Severe Warning')
        ax_wbgt.axhspan(25, 28, alpha=0.1, color='yellow', label='Warning')
        ax_wbgt.axhline(y=31, color='red', ls='--', alpha=0.5)
        ax_wbgt.axhline(y=28, color='orange', ls='--', alpha=0.5)
        ax_wbgt.set_title(f'WBGT — {label}')
        ax_wbgt.set_ylabel('WBGT [°C]')
        ax_wbgt.legend(fontsize=7, loc='upper left')
        ax_wbgt.grid(True, alpha=0.3)
        ax_wbgt.set_ylim(15, 38)

        # Patient estimation
        pops = {
            "marunouchi": PopulationExposure(67_000, 0.25, 0.22, outdoor_worker_fraction=0.08),
            "shinjuku": PopulationExposure(350_000, 0.20, 0.22, outdoor_worker_fraction=0.04),
            "residential_23ku": PopulationExposure(740_000, 0.12, 0.25, outdoor_worker_fraction=0.03),
            "suburban": PopulationExposure(920_000, 0.10, 0.28, outdoor_worker_fraction=0.02),
        }
        for d in DISTRICTS:
            results = run_sim(d, forcing)
            uhi = results['UHI_intensity']
            T_eff = forcing['T_air_C'] + uhi
            assessor = HeatStrokeRiskAssessor(pops[d])
            risk = assessor.daily_risk_profile(T_eff, forcing['rh'],
                                                forcing['wind_10m'], forcing['sw_down'])
            ax_risk.plot(hours, risk['hourly_patients'], color=COLORS[d], lw=2, label=d)

        ax_risk.set_title(f'Estimated Heat Stroke Patients — {label}')
        ax_risk.set_xlabel('Hour (JST)')
        ax_risk.set_ylabel('Patients / hour')
        ax_risk.legend(fontsize=7)
        ax_risk.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/fig5_wbgt_risk.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → figures/fig5_wbgt_risk.png")


def fig6_summary_heatmap():
    """Fig 6: Summary heatmap of UHI intensity across scenarios."""
    scenarios = ["current", "rcp45_2050", "rcp85_2050"]
    mitigations = ["baseline", "moderate_greening", "aggressive_mitigation"]

    data = np.zeros((len(DISTRICTS), len(scenarios) * len(mitigations)))
    col_labels = []

    for j, climate in enumerate(scenarios):
        forcing = generate_forcing(climate)
        for k, mit_name in enumerate(mitigations):
            col_idx = j * len(mitigations) + k
            col_labels.append(f"{climate.split('_')[0]}\n{mit_name.split('_')[0]}")
            cooling = MITIGATION_SCENARIOS[mit_name]

            for i, d in enumerate(DISTRICTS):
                heat = TOKYO_HEAT_PROFILES[d]
                if "2050" in climate:
                    heat = project_anthropogenic_heat_2050(heat)
                results = run_sim(d, forcing, cooling, heat)
                data[i, col_idx] = np.max(results['UHI_intensity'])

    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(data, aspect='auto', cmap='YlOrRd', vmin=0)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(range(len(DISTRICTS)))
    ax.set_yticklabels([d.replace('_', ' ').title() for d in DISTRICTS])

    for i in range(len(DISTRICTS)):
        for j in range(len(col_labels)):
            ax.text(j, i, f'{data[i,j]:.1f}', ha='center', va='center', fontsize=9,
                    color='white' if data[i,j] > 3 else 'black')

    ax.set_title('Peak UHI Intensity [°C] — Scenario × Mitigation Matrix')
    plt.colorbar(im, ax=ax, label='UHI Intensity [°C]')

    # Vertical separators between climate scenarios
    for sep in [2.5, 5.5]:
        ax.axvline(sep, color='white', lw=2)

    plt.tight_layout()
    plt.savefig('figures/fig6_summary_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  → figures/fig6_summary_heatmap.png")


if __name__ == "__main__":
    print("\nGenerating figures...")
    print("-" * 40)
    fig1_energy_balance()
    fig2_uhi_intensity()
    fig3_mitigation_scenarios()
    fig4_2050_projection()
    fig5_wbgt_risk()
    fig6_summary_heatmap()
    print("\n✓ All figures generated successfully.")
