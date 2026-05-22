"""
Visualization module for dark matter simulation results.
Generates publication-quality figures.
"""
import numpy as np
import json
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import matplotlib.colors as mcolors

# Colorblind-friendly palette
COLORS = {
    'Xe': '#0072B2',
    'Ar': '#E69F00',
    'Ge': '#009E73',
    'NaI': '#CC79A7',
    'combined': '#D55E00',
    'nu_floor': '#999999',
    'axion': '#56B4E9',
    'dark_photon': '#F0E442',
    'pbh': '#000000',
}

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
})


def load_results(results_dir: str) -> dict:
    """Load all result JSON files."""
    data = {}
    for fname in os.listdir(results_dir):
        if fname.endswith('.json') and fname != 'summary.json':
            key = fname.replace('.json', '')
            with open(os.path.join(results_dir, fname)) as f:
                data[key] = json.load(f)
    return data


def plot_sensitivity_curves(data: dict, output_dir: str):
    """Plot WIMP sensitivity curves for all detectors."""
    fig, ax = plt.subplots(figsize=(10, 8))

    m_dm = np.array(data['sensitivity']['m_dm_gev'])
    colors_list = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']

    for i, (name, det_data) in enumerate(data['sensitivity']['detectors'].items()):
        limits = np.array(det_data['exclusion_90cl'])
        ax.loglog(m_dm, limits, '-', color=colors_list[i % len(colors_list)],
                  linewidth=2, label=f'{name} (90% CL)')
        disc = np.array(det_data['discovery_3sigma'])
        ax.loglog(m_dm, disc, '--', color=colors_list[i % len(colors_list)],
                  linewidth=1.5, alpha=0.6, label=f'{name} (3σ disc.)')

    # Neutrino floor
    if 'neutrino_floor' in data:
        nf = data['neutrino_floor']
        m_nf = np.array(nf['m_dm_gev'])
        if 'Xe131_next_gen' in nf:
            floor = np.array(nf['Xe131_next_gen'])
            ax.fill_between(m_nf, floor, 1e-55, alpha=0.15, color='gray',
                           label='Neutrino fog (Xe)')

    ax.set_xlabel(r'Dark Matter Mass $m_\chi$ [GeV/$c^2$]')
    ax.set_ylabel(r'SI Cross Section $\sigma_{SI}$ [cm$^2$]')
    ax.set_title('WIMP-Nucleon SI Cross Section Sensitivity')
    ax.set_xlim(1, 1e4)
    ax.set_ylim(1e-50, 1e-42)
    ax.legend(loc='upper right', fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3, which='both')

    plt.savefig(os.path.join(output_dir, 'fig1_sensitivity_curves.png'))
    plt.savefig(os.path.join(output_dir, 'fig1_sensitivity_curves.svg'))
    plt.close()
    print("  📈 fig1_sensitivity_curves.png/svg")


def plot_non_wimp_candidates(data: dict, output_dir: str):
    """Plot non-WIMP candidate detection feasibility."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Axion panel
    ax = axes[0]
    axion = data['non_wimp']['axion']
    couplings = []
    events = []
    for key, val in axion.items():
        couplings.append(val['coupling'])
        events.append(val['total_events'])
    ax.barh(range(len(couplings)), events,
            color='#56B4E9', edgecolor='#0072B2')
    ax.set_yticks(range(len(couplings)))
    ax.set_yticklabels([f"$g_{{ae}}$={c:.0e}" for c in couplings])
    ax.set_xlabel('Expected Events')
    ax.set_title('Solar Axion (Axioelectric)')
    ax.set_xscale('log')
    ax.axvline(x=10, color='red', linestyle='--', alpha=0.5, label='Discovery threshold')
    ax.legend(fontsize=8)

    # Dark photon panel
    ax = axes[1]
    dp = data['non_wimp']['dark_photon']
    masses = []
    rates = []
    for key, val in dp.items():
        masses.append(val['mass_kev'])
        rates.append(val['rate_per_kg_day'])
    ax.semilogy(masses, rates, 'o-', color='#F0E442',
                markeredgecolor='#D55E00', markersize=8, linewidth=2)
    ax.set_xlabel('Dark Photon Mass [keV]')
    ax.set_ylabel('Rate [events/kg/day]')
    ax.set_title(r"Dark Photon ($\kappa = 10^{-15}$)")
    ax.grid(True, alpha=0.3)

    # PBH panel
    ax = axes[2]
    pbh = data['non_wimp']['pbh']
    pbh_masses = []
    pbh_events = []
    pbh_temps = []
    for key, val in pbh.items():
        pbh_masses.append(val['mass_g'])
        pbh_events.append(val['total_recoils'])
        pbh_temps.append(val['hawking_temp_gev'])

    ax2 = ax.twinx()
    ax.loglog(pbh_masses, pbh_events, 's-', color='#000000',
              markersize=8, linewidth=2, label='Recoil events')
    ax2.loglog(pbh_masses, pbh_temps, '^--', color='#D55E00',
               markersize=8, linewidth=1.5, label='$T_H$')
    ax.set_xlabel('PBH Mass [g]')
    ax.set_ylabel('Expected Recoils', color='#000000')
    ax2.set_ylabel('Hawking Temp [GeV]', color='#D55E00')
    ax.set_title('Primordial Black Holes')
    ax.legend(loc='upper left', fontsize=8)
    ax2.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig2_non_wimp_candidates.png'))
    plt.savefig(os.path.join(output_dir, 'fig2_non_wimp_candidates.svg'))
    plt.close()
    print("  📈 fig2_non_wimp_candidates.png/svg")


def plot_directional_detection(data: dict, output_dir: str):
    """Plot directional detector analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Angular distribution
    ax = axes[0]
    dir_data = data['directional']
    if 'angular_distributions' in dir_data:
        cos_theta = np.array(dir_data['angular_distributions']['cos_theta'])
        best = np.array(dir_data['angular_distributions']['best_case'])
        poor = np.array(dir_data['angular_distributions']['poor_case'])

        theta_deg = np.degrees(np.arccos(cos_theta))
        ax.plot(theta_deg, best, '-', color='#0072B2', linewidth=2,
                label='15° + H/T (CYGNUS-HD)')
        ax.plot(theta_deg, poor, '--', color='#CC79A7', linewidth=2,
                label='60° no-H/T (basic)')
        # Isotropic reference
        ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5,
                   label='Isotropic background')
        ax.set_xlabel(r'Angle from Cygnus direction $\theta$ [deg]')
        ax.set_ylabel('Relative Rate')
        ax.set_title('Directional Recoil Distribution ($m_\\chi$ = 50 GeV)')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Discovery reach vs angular resolution
    ax = axes[1]
    configs = ['ang15_htTrue', 'ang30_htTrue', 'ang60_htTrue',
               'ang15_htFalse', 'ang30_htFalse', 'ang60_htFalse']

    masses_to_plot = ['m=10', 'm=50', 'm=100', 'm=500']
    ang_res_values = [15, 30, 60]

    for ht in [True, False]:
        n_events_3s = []
        for ang in ang_res_values:
            key = f'ang{ang}_ht{ht}'
            if key in dir_data:
                n = dir_data[key].get('m=50', {}).get('n_events_3sigma', 100)
                n_events_3s.append(n)
            else:
                n_events_3s.append(100)

        style = '-o' if ht else '--s'
        label = 'With H/T' if ht else 'Without H/T'
        color = '#0072B2' if ht else '#E69F00'
        ax.semilogy(ang_res_values, n_events_3s, style, color=color,
                   markersize=8, linewidth=2, label=label)

    ax.set_xlabel('Angular Resolution [deg]')
    ax.set_ylabel(r'Events for 3$\sigma$ Detection')
    ax.set_title('Directional Discovery Reach ($m_\\chi$ = 50 GeV)')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(ang_res_values)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig3_directional_detection.png'))
    plt.savefig(os.path.join(output_dir, 'fig3_directional_detection.svg'))
    plt.close()
    print("  📈 fig3_directional_detection.png/svg")


def plot_neutrino_floor(data: dict, output_dir: str):
    """Plot neutrino floor predictions."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    nf = data['neutrino_floor']
    m_dm = np.array(nf['m_dm_gev'])

    # Left: Neutrino floor for different targets
    ax = axes[0]
    target_colors = {'Xe131': '#0072B2', 'Ar40': '#E69F00', 'Ge76': '#009E73'}
    for tname, color in target_colors.items():
        for exp_name, style in [('current', ':'), ('next_gen', '--'), ('ultimate', '-')]:
            key = f'{tname}_{exp_name}'
            if key in nf:
                floor = np.array(nf[key])
                valid = floor > 1e-55
                if np.any(valid):
                    ax.loglog(m_dm[valid], floor[valid], style,
                             color=color, linewidth=2,
                             label=f'{tname} ({exp_name})')

    ax.set_xlabel(r'Dark Matter Mass $m_\chi$ [GeV/$c^2$]')
    ax.set_ylabel(r'Neutrino Floor $\sigma_{SI}$ [cm$^2$]')
    ax.set_title('Neutrino Floor vs Target & Exposure')
    ax.set_xlim(1, 1e4)
    ax.set_ylim(1e-52, 1e-44)
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')

    # Right: Neutrino recoil spectra
    ax = axes[1]
    Er = np.array(nf['Er_kev'])
    spectra = nf.get('neutrino_spectra', {})

    source_styles = {
        'pp': ('-', '#0072B2'), '7Be_862': ('-', '#E69F00'),
        '8B': ('-', '#009E73'), 'hep': ('--', '#CC79A7'),
        'atm': ('--', '#D55E00')
    }

    for src, (style, color) in source_styles.items():
        key = f'Xe131_{src}'
        if key in spectra:
            rate = np.array(spectra[key])
            valid = rate > 0
            if np.any(valid):
                ax.semilogy(Er[valid], rate[valid], style, color=color,
                           linewidth=2, label=f'{src} (Xe)')

    ax.set_xlabel('Recoil Energy [keV]')
    ax.set_ylabel('Rate [events/keV/kg/day]')
    ax.set_title(r'CE$\nu$NS Recoil Spectra (Xe target)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 50)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig4_neutrino_floor.png'))
    plt.savefig(os.path.join(output_dir, 'fig4_neutrino_floor.svg'))
    plt.close()
    print("  📈 fig4_neutrino_floor.png/svg")


def plot_background_strategies(data: dict, output_dir: str):
    """Plot background reduction strategy comparison."""
    fig, ax = plt.subplots(figsize=(10, 6))

    bg = data['backgrounds']
    targets = list(bg.keys())
    strategies = list(bg[targets[0]]['strategies'].keys())

    x = np.arange(len(strategies))
    width = 0.25
    target_colors = {'Xe131': '#0072B2', 'Ar40': '#E69F00', 'Ge76': '#009E73'}

    for i, tname in enumerate(targets):
        reductions = []
        for s in strategies:
            rf = bg[tname]['strategies'][s].get('reduction_factor', 1.0)
            reductions.append(rf if rf is not None else 1.0)
        ax.bar(x + i * width, reductions, width,
               color=target_colors.get(tname, '#999999'),
               label=tname, edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Background Reduction Strategy')
    ax.set_ylabel('Reduction Factor (relative to baseline)')
    ax.set_title('Systematic Background Reduction Evaluation')
    ax.set_xticks(x + width)
    ax.set_xticklabels(strategies, rotation=30, ha='right')
    ax.legend()
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig5_background_strategies.png'))
    plt.savefig(os.path.join(output_dir, 'fig5_background_strategies.svg'))
    plt.close()
    print("  📈 fig5_background_strategies.png/svg")


def plot_multi_target(data: dict, output_dir: str):
    """Plot multi-target complementarity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    mt = data['multi_target']
    m_dm = np.array(mt['m_dm_gev'])

    # Left: Sensitivity comparison
    ax = axes[0]
    det_colors = {
        'Xe (DARWIN)': '#0072B2',
        'Ar (DS-20k)': '#E69F00',
        'Ge (SuperCDMS)': '#009E73',
        'NaI (COSINE)': '#CC79A7',
        'combined': '#D55E00',
    }

    for name, limits in mt['sensitivities'].items():
        limits = np.array(limits)
        color = det_colors.get(name, '#999999')
        style = '-' if name != 'combined' else '-'
        lw = 3 if name == 'combined' else 1.5
        ax.loglog(m_dm, limits, style, color=color, linewidth=lw, label=name)

    ax.set_xlabel(r'Dark Matter Mass $m_\chi$ [GeV/$c^2$]')
    ax.set_ylabel(r'$\sigma_{SI}$ [cm$^2$] (90% CL)')
    ax.set_title('Multi-Target Sensitivity Comparison')
    ax.set_xlim(1, 1e4)
    ax.set_ylim(1e-50, 1e-40)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # Right: Response at m=50 GeV
    ax = axes[1]
    resp = mt.get('response_m50', {})
    if resp:
        names = list(resp.keys())
        events = [resp[n]['total_events'] for n in names]
        peak_E = [resp[n]['peak_energy_kev'] for n in names]

        bar_colors = [det_colors.get(n, '#999999') for n in names]
        bars = ax.bar(range(len(names)), events, color=bar_colors,
                     edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels([n.split('(')[0].strip() for n in names],
                          rotation=15, ha='right')
        ax.set_ylabel('Expected Events')
        ax.set_title(r'Response at $m_\chi=50$ GeV, $\sigma=10^{-46}$ cm$^2$')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, axis='y')

        for bar, e_peak in zip(bars, peak_E):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                   f'E_peak={e_peak:.0f} keV', ha='center', va='bottom',
                   fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig6_multi_target.png'))
    plt.savefig(os.path.join(output_dir, 'fig6_multi_target.svg'))
    plt.close()
    print("  📈 fig6_multi_target.png/svg")


def plot_annual_modulation(data: dict, output_dir: str):
    """Plot annual modulation analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    am = data['annual_modulation']

    # Left: Time-dependent rate
    ax = axes[0]
    colors = {'NaI_DAMA': '#CC79A7', 'Xe_DARWIN': '#0072B2', 'Ar_DS20k': '#E69F00'}

    for name, mod_data in am.items():
        t = np.array(mod_data['t_days'])
        rates = np.array(mod_data['daily_rates'])
        if np.max(rates) > 0:
            norm_rates = rates / np.mean(rates)
            ax.plot(t, norm_rates, '-', color=colors.get(name, '#999999'),
                   linewidth=2, label=name)

    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Normalized Rate')
    ax.set_title('Annual Modulation Signal')
    ax.axvline(x=152, color='red', linestyle=':', alpha=0.3, label='June 2 (peak)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 365)

    # Right: Significance vs observation time
    ax = axes[1]
    for name, mod_data in am.items():
        years = []
        sigmas = []
        for yr_key, sig_data in mod_data['significance'].items():
            yr = float(yr_key.replace('yr', ''))
            years.append(yr)
            sigmas.append(sig_data['significance_sigma'])

        ax.plot(years, sigmas, 'o-', color=colors.get(name, '#999999'),
               markersize=8, linewidth=2, label=name)

    ax.axhline(y=3, color='green', linestyle='--', alpha=0.5, label=r'3$\sigma$')
    ax.axhline(y=5, color='red', linestyle='--', alpha=0.5, label=r'5$\sigma$')
    ax.set_xlabel('Observation Time [years]')
    ax.set_ylabel(r'Significance [$\sigma$]')
    ax.set_title('Modulation Detection Significance')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig7_annual_modulation.png'))
    plt.savefig(os.path.join(output_dir, 'fig7_annual_modulation.svg'))
    plt.close()
    print("  📈 fig7_annual_modulation.png/svg")


def generate_all_figures(results_dir: str = 'results',
                          figures_dir: str = 'figures'):
    """Generate all publication figures."""
    print("\n" + "=" * 50)
    print("🎨 Generating Publication Figures")
    print("=" * 50)

    os.makedirs(figures_dir, exist_ok=True)
    data = load_results(results_dir)

    if 'sensitivity' in data:
        plot_sensitivity_curves(data, figures_dir)
    if 'non_wimp' in data:
        plot_non_wimp_candidates(data, figures_dir)
    if 'directional' in data:
        plot_directional_detection(data, figures_dir)
    if 'neutrino_floor' in data:
        plot_neutrino_floor(data, figures_dir)
    if 'backgrounds' in data:
        plot_background_strategies(data, figures_dir)
    if 'multi_target' in data:
        plot_multi_target(data, figures_dir)
    if 'annual_modulation' in data:
        plot_annual_modulation(data, figures_dir)

    print("\n✅ All figures generated successfully.")


if __name__ == '__main__':
    generate_all_figures()
