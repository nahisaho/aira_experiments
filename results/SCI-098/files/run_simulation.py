"""
Monte Carlo simulation runner for dark matter direct detection.
GEANT4/ROOT-inspired architecture with event generation, propagation,
and analysis chain.
"""
import numpy as np
import json
import os
import sys
from datetime import datetime
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.core.constants import TARGETS, NuclearTarget, V_EARTH, V_0
from src.signals.dm_signals import (
    WIMPSignal, AxionSignal, DarkPhotonSignal, PrimordialBHSignal, SignalFactory
)
from src.detectors.detector_models import (
    DetectorConfig, DirectionalDetector, LiquidNobleDetector,
    get_xenon_nt, get_darwin, get_darkside20k, get_supercdms,
    get_cygnus, get_cosine100
)
from src.backgrounds.background_models import (
    NeutrinoFloor, BackgroundBudget
)
from src.analysis.statistics import (
    SensitivityCalculator, AnnualModulation, MultiTargetComplementarity,
    NeutrinoFloorCalculator
)


class MCEventGenerator:
    """Monte Carlo event generator (GEANT4-inspired)."""

    def __init__(self, seed: int = 42):
        self.rng = np.random.RandomState(seed)

    def generate_dm_events(self, signal_model, Er_range: Tuple[float, float],
                           n_events: int = 10000) -> np.ndarray:
        """Generate MC events from signal model using accept-reject."""
        Er_grid = np.linspace(Er_range[0], Er_range[1], 1000)
        dR = signal_model.differential_rate(Er_grid)

        if np.max(dR) == 0:
            return np.array([])

        dR_max = np.max(dR) * 1.1

        events = []
        n_generated = 0
        batch_size = n_events * 10

        while len(events) < n_events and n_generated < batch_size * 10:
            Er_try = self.rng.uniform(Er_range[0], Er_range[1], batch_size)
            u = self.rng.uniform(0, dR_max, batch_size)

            dR_try = np.interp(Er_try, Er_grid, dR)
            accepted = Er_try[u < dR_try]
            events.extend(accepted.tolist())
            n_generated += batch_size

        return np.array(events[:n_events])

    def generate_background_events(self, bg_rate: float,
                                     Er_range: Tuple[float, float],
                                     exposure_kg_day: float) -> np.ndarray:
        """Generate flat background events."""
        n_expected = bg_rate * (Er_range[1] - Er_range[0]) * exposure_kg_day
        n_events = self.rng.poisson(n_expected)
        return self.rng.uniform(Er_range[0], Er_range[1], n_events)


class SimulationRunner:
    """Main simulation orchestrator."""

    def __init__(self, output_dir: str = '.'):
        self.output_dir = output_dir
        self.results = {}
        self.mc = MCEventGenerator(seed=2024)

    def run_all(self):
        """Execute complete simulation suite."""
        print("=" * 70)
        print("Dark Matter Direct Detection Simulation Framework v1.0")
        print("GEANT4/ROOT-inspired Monte Carlo Framework")
        print(f"Run started: {datetime.now().isoformat()}")
        print("=" * 70)

        self.run_sensitivity_scan()
        self.run_non_wimp_candidates()
        self.run_directional_analysis()
        self.run_neutrino_floor()
        self.run_background_evaluation()
        self.run_multi_target_analysis()
        self.run_annual_modulation()

        self.save_results()
        print("\n✅ All simulations completed successfully.")

    def run_sensitivity_scan(self):
        """Task 1: WIMP sensitivity scan for multiple detectors."""
        print("\n" + "─" * 50)
        print("📊 [1/7] WIMP Sensitivity Scan")
        print("─" * 50)

        m_dm = np.logspace(np.log10(1), np.log10(1e4), 80)
        detectors = {
            'XENON-nT': get_xenon_nt(),
            'DARWIN': get_darwin(),
            'DarkSide-20k': get_darkside20k(),
            'SuperCDMS': get_supercdms(),
        }

        sensitivity_results = {}
        for name, det in detectors.items():
            calc = SensitivityCalculator(
                det.target, det.exposure_kg_day,
                det.threshold_kev, det.max_energy_kev,
                det.background_rate, det.efficiency
            )
            limits = calc.exclusion_limit_90cl(m_dm)
            disc = calc.discovery_reach_3sigma(m_dm)
            # Cap at physical bounds
            limits = np.clip(limits, 1e-50, 1e-38)
            disc = np.clip(disc, 1e-50, 1e-38)
            sensitivity_results[name] = {
                'exclusion_90cl': limits.tolist(),
                'discovery_3sigma': disc.tolist(),
            }
            valid = np.isfinite(limits) & (limits < 1e-38)
            if np.any(valid):
                min_limit = np.min(limits[valid])
                best_mass = m_dm[valid][np.argmin(limits[valid])]
            else:
                min_limit = limits[len(limits)//2]
                best_mass = m_dm[len(m_dm)//2]
            print(f"  {name}: σ_min = {min_limit:.2e} cm² at m_χ = {best_mass:.1f} GeV")

        self.results['sensitivity'] = {
            'm_dm_gev': m_dm.tolist(),
            'detectors': sensitivity_results
        }

    def run_non_wimp_candidates(self):
        """Task 2: Non-WIMP DM candidates detection feasibility."""
        print("\n" + "─" * 50)
        print("🔬 [2/7] Non-WIMP Dark Matter Candidates")
        print("─" * 50)

        target = TARGETS['Xe131']
        exposure = 40000 * 365 * 10  # DARWIN-scale, 10 year

        # Axion search
        print("  Axion (solar):")
        E_kev = np.linspace(0.5, 15, 200)
        axion_results = {}
        for g_ae in [1e-11, 1e-12, 1e-13]:
            ax = AxionSignal(1.0, g_ae, target, exposure)
            rate = ax.differential_rate(E_kev)
            total = ax.total_events((0.5, 15))
            axion_results[f'g_ae={g_ae:.0e}'] = {
                'coupling': g_ae,
                'total_events': total,
                'peak_rate': float(np.max(rate)),
            }
            print(f"    g_ae = {g_ae:.0e}: {total:.1f} events")

        # Dark photon
        print("  Dark Photon:")
        dp_results = {}
        for m_dp in [0.5, 1.0, 5.0, 10.0]:
            dp = DarkPhotonSignal(m_dp, 1e-15, target, exposure)
            rate = dp.dm_absorption_rate()
            dp_results[f'm_dp={m_dp}keV'] = {
                'mass_kev': m_dp,
                'kappa': 1e-15,
                'rate_per_kg_day': rate,
                'total_events': rate * exposure,
            }
            print(f"    m_A' = {m_dp} keV: rate = {rate:.2e} /kg/day")

        # Primordial Black Holes
        print("  Primordial Black Holes:")
        pbh_results = {}
        Er_kev = np.linspace(1, 50, 100)
        for m_pbh in [1e15, 1e16, 1e17]:
            pbh = PrimordialBHSignal(m_pbh, 1.0, target, exposure)
            T_H = pbh.hawking_temperature_gev()
            recoils = pbh.neutrino_induced_recoils(Er_kev)
            total = np.trapz(recoils, Er_kev) * exposure
            pbh_results[f'm_pbh={m_pbh:.0e}g'] = {
                'mass_g': m_pbh,
                'hawking_temp_gev': T_H,
                'total_recoils': total,
            }
            print(f"    M_PBH = {m_pbh:.0e} g: T_H = {T_H:.2e} GeV, events = {total:.2e}")

        self.results['non_wimp'] = {
            'axion': axion_results,
            'dark_photon': dp_results,
            'pbh': pbh_results,
        }

    def run_directional_analysis(self):
        """Task 3: Directional detector sensitivity (CYGNUS/MIMAC)."""
        print("\n" + "─" * 50)
        print("🧭 [3/7] Directional Detector Analysis (CYGNUS/MIMAC)")
        print("─" * 50)

        cygnus_config = get_cygnus()
        dir_results = {}

        for ang_res in [15, 30, 60]:
            for ht in [True, False]:
                det = DirectionalDetector(
                    cygnus_config,
                    angular_resolution_deg=ang_res,
                    head_tail_recognition=ht
                )

                m_dm_values = [10, 50, 100, 500]
                reach = {}
                for m_dm in m_dm_values:
                    r = det.discovery_reach_directional(m_dm)
                    reach[f'm={m_dm}'] = r

                key = f'ang{ang_res}_ht{ht}'
                dir_results[key] = reach
                n3 = reach['m=50']['n_events_3sigma']
                n5 = reach['m=50']['n_events_5sigma']
                ht_str = "H/T" if ht else "no-H/T"
                print(f"  {ang_res}° {ht_str}: 3σ={n3}, 5σ={n5} events (m_χ=50 GeV)")

        # Angular distribution
        cos_theta = np.linspace(-1, 1, 200)
        det_best = DirectionalDetector(cygnus_config, 15, True)
        det_poor = DirectionalDetector(cygnus_config, 60, False)

        ang_dist_best = det_best.recoil_direction_distribution(cos_theta, 50.0)
        ang_dist_poor = det_poor.recoil_direction_distribution(cos_theta, 50.0)

        dir_results['angular_distributions'] = {
            'cos_theta': cos_theta.tolist(),
            'best_case': ang_dist_best.tolist(),
            'poor_case': ang_dist_poor.tolist(),
        }

        self.results['directional'] = dir_results

    def run_neutrino_floor(self):
        """Task 4: Neutrino floor calculation."""
        print("\n" + "─" * 50)
        print("🌊 [4/7] Neutrino Floor (CEνNS) Prediction")
        print("─" * 50)

        m_dm = np.logspace(np.log10(1), np.log10(1e4), 60)
        floor_results = {}

        targets_for_floor = {
            'Xe131': TARGETS['Xe131'],
            'Ar40': TARGETS['Ar40'],
            'Ge76': TARGETS['Ge76'],
        }

        exposures = {
            'current': 1e3,       # ~ton-year
            'next_gen': 1e5,      # ~100 ton-year
            'ultimate': 1e7,      # ~10 kton-year
        }

        for tname, target in targets_for_floor.items():
            floor_calc = NeutrinoFloorCalculator(target)

            for exp_name, exp_val in exposures.items():
                floor = floor_calc.compute_floor(m_dm, exp_val)
                key = f'{tname}_{exp_name}'
                floor_results[key] = floor.tolist()
                min_floor = np.min(floor)
                print(f"  {tname} ({exp_name}): ν-floor min = {min_floor:.2e} cm²")

        # Neutrino recoil spectra
        Er = np.linspace(0.1, 50, 200)
        nu_spectra = {}
        for tname, target in targets_for_floor.items():
            nu = NeutrinoFloor(target)
            for src in ['pp', '7Be_862', '8B', 'hep', 'atm']:
                rate = nu.recoil_rate_source(src, Er, n_Enu=200)
                nu_spectra[f'{tname}_{src}'] = rate.tolist()

        floor_results['Er_kev'] = Er.tolist()
        floor_results['neutrino_spectra'] = nu_spectra
        floor_results['m_dm_gev'] = m_dm.tolist()
        self.results['neutrino_floor'] = floor_results

    def run_background_evaluation(self):
        """Task 5: Systematic background reduction evaluation."""
        print("\n" + "─" * 50)
        print("🛡️ [5/7] Background Reduction Strategy Evaluation")
        print("─" * 50)

        Er = np.linspace(0.5, 50, 200)
        bg_results = {}

        for tname in ['Xe131', 'Ar40', 'Ge76']:
            target = TARGETS[tname]
            budget = BackgroundBudget(target, shielding_factor=1e-6)

            # Baseline backgrounds
            bg_components = budget.total_background(Er, include_neutrinos=True)
            total = float(np.trapz(bg_components['total'], Er))

            # Reduction strategies
            strategies = budget.evaluate_reduction_strategies(Er)

            bg_results[tname] = {
                'baseline_total_rate': total,
                'strategies': {}
            }

            for sname, sdata in strategies.items():
                bg_results[tname]['strategies'][sname] = {
                    'params': sdata['params'],
                    'total_rate': sdata['total_rate'],
                    'reduction_factor': sdata['reduction_factor'],
                }
                rf = sdata['reduction_factor']
                if rf is not None:
                    print(f"  {tname} [{sname}]: reduction = {rf:.3f}×")

        self.results['backgrounds'] = bg_results

    def run_multi_target_analysis(self):
        """Task 6: Multi-target complementarity analysis."""
        print("\n" + "─" * 50)
        print("🎯 [6/7] Multi-Target Complementarity (Xe/Ar/Ge/NaI)")
        print("─" * 50)

        targets = {
            'Xe (DARWIN)': {
                'target': TARGETS['Xe131'],
                'exposure': 40000 * 3650,
                'threshold': 0.5,
                'max_energy': 70.0,
                'background': 5e-6,
                'efficiency': 0.90,
            },
            'Ar (DS-20k)': {
                'target': TARGETS['Ar40'],
                'exposure': 20000 * 3650,
                'threshold': 0.6,
                'max_energy': 200.0,
                'background': 1e-6,
                'efficiency': 0.80,
            },
            'Ge (SuperCDMS)': {
                'target': TARGETS['Ge76'],
                'exposure': 30 * 1825,
                'threshold': 0.04,
                'max_energy': 50.0,
                'background': 1e-4,
                'efficiency': 0.70,
            },
            'NaI (COSINE)': {
                'target': TARGETS['Na23'],
                'exposure': 106 * 1095,
                'threshold': 1.0,
                'max_energy': 20.0,
                'background': 2.7,
                'efficiency': 0.65,
            },
        }

        mtc = MultiTargetComplementarity(targets)
        m_dm = np.logspace(np.log10(1), np.log10(1e4), 80)

        # Sensitivity comparison
        sensitivities = mtc.combined_sensitivity(m_dm)
        # Cap values
        for k in sensitivities:
            sensitivities[k] = np.clip(sensitivities[k], 1e-50, 1e-38)
        mt_results = {
            'm_dm_gev': m_dm.tolist(),
            'sensitivities': {k: v.tolist() for k, v in sensitivities.items()},
        }

        # Response at benchmark point
        response = mtc.target_response_matrix(50.0, 1e-46)
        mt_results['response_m50'] = response

        # SI/SD discrimination
        si_sd = mtc.si_sd_discrimination(m_dm)
        mt_results['si_sd_discrimination'] = si_sd

        for name, resp in response.items():
            print(f"  {name}: {resp['total_events']:.1f} events (m_χ=50 GeV, σ=10⁻⁴⁶)")

        # Print combined improvement
        xe_min = np.min(sensitivities['Xe (DARWIN)'])
        comb_min = np.min(sensitivities['combined'])
        improvement = xe_min / comb_min if comb_min > 0 else 0
        print(f"  Combined improvement over best single: {improvement:.1f}×")

        self.results['multi_target'] = mt_results

    def run_annual_modulation(self):
        """Task 7: Annual modulation statistical power."""
        print("\n" + "─" * 50)
        print("📅 [7/7] Annual Modulation Signal Analysis")
        print("─" * 50)

        mod_results = {}

        configs = {
            'NaI_DAMA': {'target': 'Na23', 'm_dm': 10, 'sigma': 2e-41,
                         'exposure': 250 * 365, 'bg': 1.0, 'Er': (2, 6)},
            'Xe_DARWIN': {'target': 'Xe131', 'm_dm': 50, 'sigma': 1e-46,
                          'exposure': 40000 * 365, 'bg': 5e-6, 'Er': (1, 20)},
            'Ar_DS20k': {'target': 'Ar40', 'm_dm': 100, 'sigma': 1e-46,
                         'exposure': 20000 * 365, 'bg': 1e-6, 'Er': (5, 50)},
        }

        for name, cfg in configs.items():
            target = TARGETS[cfg['target']]
            am = AnnualModulation(target, cfg['m_dm'], cfg['sigma'], cfg['exposure'])

            # Modulation fraction
            mod_frac = am.modulation_fraction(cfg['Er'])

            # Detection significance for different observation periods
            sig_results = {}
            for n_yr in [1, 3, 5, 10]:
                sig = am.detection_significance(n_yr, cfg['Er'], cfg['bg'])
                sig_results[f'{n_yr}yr'] = sig
                print(f"  {name} ({n_yr}yr): σ = {sig['significance_sigma']:.2f}, "
                      f"N = {sig['total_events']:.0f}")

            # Time-dependent rate
            t_days = np.linspace(0, 365, 365)
            Er_grid = np.linspace(cfg['Er'][0], cfg['Er'][1], 50)
            daily_rates = []
            for t in t_days:
                R = am.modulated_rate(Er_grid, t)
                daily_rates.append(float(np.trapz(R, Er_grid)))

            mod_results[name] = {
                'config': {k: v for k, v in cfg.items() if k != 'target'},
                'target': cfg['target'],
                'modulation_fraction': mod_frac,
                'significance': sig_results,
                'daily_rates': daily_rates,
                't_days': t_days.tolist(),
            }

        self.results['annual_modulation'] = mod_results

    def save_results(self):
        """Save all results to JSON files."""
        results_dir = os.path.join(self.output_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)

        for key, data in self.results.items():
            filepath = os.path.join(results_dir, f'{key}.json')
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=2, default=str)
            print(f"  💾 Saved: results/{key}.json")

        # Save summary
        summary = {
            'framework': 'DMDDSF v1.0',
            'run_date': datetime.now().isoformat(),
            'modules': list(self.results.keys()),
            'total_results_files': len(self.results),
        }
        with open(os.path.join(results_dir, 'summary.json'), 'w') as f:
            json.dump(summary, f, indent=2)


if __name__ == '__main__':
    runner = SimulationRunner(output_dir='.')
    runner.run_all()
