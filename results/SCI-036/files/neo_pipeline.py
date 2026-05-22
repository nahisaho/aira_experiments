"""
NEO Risk Assessment Pipeline — Main Orchestrator
Integrates all six modules into a unified probabilistic risk assessment workflow.
"""

import numpy as np
import json
import time
import os
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# Import pipeline modules
from neo_orbital_propagation import (
    MonteCarloOrbitalPropagator, OrbitalElements, demo_neo
)
from neo_perturbations import (
    GravitationalPerturbationAnalyzer, YarkovskyParameters,
    compute_yarkovsky_uncertainty_band
)
from neo_keyhole import (
    BPlaneAnalyzer, compute_torino_palermo_scales
)
from neo_bayesian_update import (
    BayesianImpactUpdater, PriorBelief, AstrometricObservation
)
from neo_impact_model import (
    ImpactScenario, ImpactDamageModel, AtmosphericEntryModel,
    plot_damage_map, compare_scenarios
)
from neo_deflection import (
    NEOPhysicalParams, KineticImpactorMission, GravityTractorMission,
    simulate_full_deflection_campaign
)


LOG_PATH = "logs/process-log.jsonl"


def log_event(phase: str, event_type: str, skill: str,
               handoff_in: dict = None, handoff_out: dict = None,
               files_written: list = None, status: str = "ok") -> None:
    """Append an event to the process log."""
    entry = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files_written or [],
        "status": status,
    }
    os.makedirs("logs", exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry) + "\n")


def run_pipeline(
    n_clones: int = 200,
    t_years: float = 100.0,
    run_mcmc: bool = False,
    seed: int = 42
) -> dict:
    """
    Execute the complete NEO risk assessment pipeline.

    Parameters
    ----------
    n_clones  : Number of virtual asteroid clones for Monte Carlo
    t_years   : Integration time span [yr]
    run_mcmc  : Whether to run full PyMC MCMC (slow; ~3 min)
    seed      : Random seed
    """
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    log_event("pipeline", "run_started", "neo_pipeline",
               handoff_in={"n_clones": n_clones, "t_years": t_years})
    print("\n" + "="*60)
    print(" NEO RISK ASSESSMENT PIPELINE")
    print("="*60)

    # ----------------------------------------------------------------
    # 1. Define representative NEO (Apophis-like)
    # ----------------------------------------------------------------
    neo = demo_neo()
    print(f"\n[NEO] Target: Apophis-like NEO")
    print(f"  a = {neo.a:.4f} AU, e = {neo.e:.4f}, i = {np.degrees(neo.inc):.2f}°")
    print(f"  D = {neo.diameter:.2f} km, ρ = {neo.density:.0f} kg/m³")

    log_event("phase1", "prompt_received", "neo_orbital_propagation",
               handoff_in={"neo": str(neo)})

    # ----------------------------------------------------------------
    # 2. Monte Carlo Orbital Propagation
    # ----------------------------------------------------------------
    print("\n[Phase 1] Monte Carlo Orbital Propagation")
    t0 = time.time()
    propagator = MonteCarloOrbitalPropagator(neo, n_clones=n_clones, seed=seed)
    mc_results = propagator.propagate(t_years=t_years, n_outputs=100)
    propagator.plot_minimum_distances(mc_results, "figures/fig1_monte_carlo_orbits.png")
    elapsed_mc = time.time() - t0

    log_event("phase1", "file_written", "neo_orbital_propagation",
               files_written=["figures/fig1_monte_carlo_orbits.png"],
               handoff_out={"impact_prob": mc_results['impact_probability']})

    # ----------------------------------------------------------------
    # 3. Gravitational Perturbations + Yarkovsky
    # ----------------------------------------------------------------
    print("\n[Phase 2] Gravitational Perturbations & Yarkovsky Effect")
    t0 = time.time()
    yp = YarkovskyParameters(
        diameter=neo.diameter * 1000,  # km → m
        density=neo.density,
        albedo=0.23,
        thermal_conductivity=0.01,
        rotation_period=30.56 * 3600,   # 30.56 h for Apophis
        obliquity=np.radians(170),       # retrograde-like
    )
    perturb_analyzer = GravitationalPerturbationAnalyzer(neo, yarkovsky_params=yp)
    perturb_results = perturb_analyzer.analyze_secular_evolution(t_years=100.0, n_outputs=300)
    perturb_analyzer.plot_perturbations(perturb_results, "figures/fig2_perturbations.png")

    yark_band = compute_yarkovsky_uncertainty_band(yp, neo.a, 100.0, n_samples=500, seed=seed)
    elapsed_perturb = time.time() - t0

    log_event("phase2", "file_written", "neo_perturbations",
               files_written=["figures/fig2_perturbations.png"],
               handoff_out={"yark_drift_au": yark_band['mean']})

    # ----------------------------------------------------------------
    # 4. Keyhole Analysis
    # ----------------------------------------------------------------
    print("\n[Phase 3] Keyhole (Resonant Return) Search")
    t0 = time.time()
    # Apophis 2029 encounter: v_inf ≈ 5.87 km/s
    encounter_vel = 5.87
    # Uncertainty based on MC orbital spread (use a floor to avoid NaN/zero)
    mc_std = np.std(mc_results['min_distances_au'])
    sigma_b_km = mc_std * 1.496e8 if np.isfinite(mc_std) and mc_std > 0 else 2000.0
    sigma_b_km = max(sigma_b_km, 500.0)  # floor at 500 km

    bplane = BPlaneAnalyzer(neo, uncertainty_sigma_km=sigma_b_km)
    keyholes = bplane.search_keyholes(encounter_vel_kms=encounter_vel, n_resonances=7)
    bplane.plot_b_plane(keyholes, encounter_vel, "figures/fig3_b_plane_keyholes.png")

    # Torino/Palermo scales
    min_dist = mc_results['min_distances_au']
    best_approach = float(np.nanmin(min_dist[min_dist < 1e10])) if mc_results['n_impacts'] > 0 else 0.001
    scales = compute_torino_palermo_scales(
        mc_results['impact_probability'], neo.diameter, encounter_vel,
        dist_au=best_approach
    )
    elapsed_keyhole = time.time() - t0

    log_event("phase3", "file_written", "neo_keyhole",
               files_written=["figures/fig3_b_plane_keyholes.png"],
               handoff_out={"n_keyholes": len(keyholes),
                             "torino": scales['torino_scale'],
                             "palermo": scales['palermo_scale']})

    # ----------------------------------------------------------------
    # 5. Bayesian Probability Update
    # ----------------------------------------------------------------
    print("\n[Phase 4] Bayesian Probability Update")
    t0 = time.time()
    prior = PriorBelief(
        a_mean=neo.a, a_std=neo.da * 1000,
        e_mean=neo.e, e_std=neo.de * 1000,
        impact_prob_prior=max(mc_results['impact_probability'], 1e-6)
    )
    updater = BayesianImpactUpdater(prior, seed=seed)
    obs_history = updater.simulate_observation_campaign(n_obs=30, sigma_arcsec=0.3)

    idata = None
    if run_mcmc:
        idata = updater.full_mcmc_posterior(n_obs=20)

    updater.plot_bayesian_update(obs_history, idata, "figures/fig4_bayesian_update.png")
    elapsed_bayes = time.time() - t0

    final_prob = obs_history[-1]['impact_prob']
    log_event("phase4", "file_written", "neo_bayesian_update",
               files_written=["figures/fig4_bayesian_update.png"],
               handoff_out={"final_impact_prob": final_prob,
                             "n_observations": len(obs_history)})

    # ----------------------------------------------------------------
    # 6. Impact Energy & Damage Model
    # ----------------------------------------------------------------
    print("\n[Phase 5] Impact Energy & Damage Estimation")
    t0 = time.time()
    scenario_apophis = ImpactScenario(
        diameter_km=neo.diameter,
        density_kg_m3=neo.density,
        velocity_km_s=encounter_vel,
        entry_angle_deg=20.0,
        target_type='land',
    )
    damage_model = ImpactDamageModel(scenario_apophis)
    damage = damage_model.estimate_damage()
    plot_damage_map(damage, scenario_apophis, "figures/fig5_impact_damage.png")

    # Atmospheric entry trajectory
    entry_model = AtmosphericEntryModel(scenario_apophis)
    entry_results = entry_model.integrate()
    _plot_entry_trajectory(entry_results, "figures/fig5b_entry_trajectory.png")

    # Multi-scenario comparison
    scenarios_compare = [
        ImpactScenario(d, 2000, encounter_vel, 20.0, 'land')
        for d in [0.01, 0.05, 0.14, 0.37, 1.0, 5.0, 10.0]
    ]
    damages_compare = compare_scenarios(scenarios_compare, "figures/fig5c_scenario_comparison.png")
    elapsed_damage = time.time() - t0

    log_event("phase5", "file_written", "neo_impact_model",
               files_written=["figures/fig5_impact_damage.png",
                               "figures/fig5b_entry_trajectory.png",
                               "figures/fig5c_scenario_comparison.png"],
               handoff_out={"KE_Mt": damage.kinetic_energy_Mt,
                             "crater_km": damage.crater_diameter_km,
                             "torino": damage.torino_scale})

    # ----------------------------------------------------------------
    # 7. DART/Hera Deflection Simulation
    # ----------------------------------------------------------------
    print("\n[Phase 6] Deflection Mission Simulation (DART/Hera)")
    t0 = time.time()
    neo_phys = NEOPhysicalParams(
        diameter_km=neo.diameter,
        density_kg_m3=neo.density,
        beta_mean=2.5, beta_min=1.0, beta_max=5.0
    )
    ki_mission = KineticImpactorMission(
        spacecraft_mass_kg=570.0, impact_velocity_km_s=6.14,
        lead_time_years=10.0, mission_name="DART-like"
    )
    gt_mission = GravityTractorMission(
        spacecraft_mass_kg=1000.0, hover_distance_m=300.0,
        thrust_N=0.04, operating_time_years=5.0
    )
    deflect_results = simulate_full_deflection_campaign(
        neo_phys, ki_mission, gt_mission, "figures/fig6_deflection.png"
    )
    elapsed_deflect = time.time() - t0

    log_event("phase6", "file_written", "neo_deflection",
               files_written=["figures/fig6_deflection.png"],
               handoff_out={"dart_dv_m_s": deflect_results['dart']['mean_delta_v_m_s'],
                             "success_prob": deflect_results['dart']['success_probability']})

    # ----------------------------------------------------------------
    # 8. Summary Dashboard Figure
    # ----------------------------------------------------------------
    print("\n[Final] Generating summary dashboard")
    _plot_summary_dashboard(
        mc_results, perturb_results, yark_band, keyholes,
        obs_history, damage, deflect_results, scales,
        "figures/fig0_summary_dashboard.png"
    )

    # ----------------------------------------------------------------
    # 9. Save Numerical Results
    # ----------------------------------------------------------------
    results_summary = {
        "neo_name": "Apophis-like NEO (demo)",
        "orbital_elements": {
            "a_AU": neo.a, "e": neo.e, "inc_deg": np.degrees(neo.inc),
            "diameter_km": neo.diameter, "density_kg_m3": neo.density
        },
        "monte_carlo": {
            "n_clones": mc_results['n_clones'],
            "t_years": mc_results['t_years'],
            "impact_probability_raw": mc_results['impact_probability'],
            "n_impacts": mc_results['n_impacts'],
        },
        "yarkovsky": {
            "da_dt_AU_per_yr": float(perturb_results.get('yarkovsky_rate_au_per_yr', 0)),
            "total_drift_100yr_AU": float(perturb_results.get('yarkovsky_total_drift_au', 0)),
            "uncertainty_band_p5_p95": [float(yark_band['p5']), float(yark_band['p95'])],
        },
        "keyholes": {
            "n_found": len(keyholes),
            "top_keyhole": {
                "resonance": f"{keyholes[0].resonance_n}:{keyholes[0].resonance_m}" if keyholes else "N/A",
                "width_km": keyholes[0].half_width_km if keyholes else 0,
                "impact_prob": keyholes[0].impact_prob if keyholes else 0,
            }
        },
        "bayesian_update": {
            "prior_impact_prob": prior.impact_prob_prior,
            "posterior_impact_prob_after_30obs": final_prob,
            "uncertainty_reduction_a": obs_history[-1]['a_std'] / prior.a_std,
        },
        "impact_damage": {
            "kinetic_energy_Mt": float(damage.kinetic_energy_Mt),
            "crater_diameter_km": float(damage.crater_diameter_km),
            "blast_severe_radius_km": float(damage.blast_radius_km['severe_damage_0_6_bar']),
            "thermal_radius_km": float(damage.thermal_radius_km),
            "casualties_median": int(damage.casualties_estimate['median']),
            "torino_scale": int(damage.torino_scale),
            "palermo_scale": float(damage.palermo_scale),
        },
        "deflection": {
            "dart_delta_v_mean_m_s": float(deflect_results['dart']['mean_delta_v_m_s']),
            "dart_miss_distance_mean_km": float(deflect_results['dart']['mean_miss_km']),
            "dart_success_probability": float(deflect_results['dart']['success_probability']),
            "gravity_tractor_delta_v_m_s": float(deflect_results['gravity_tractor']['delta_v_m_s']),
            "hera_beta_improvement_factor": float(deflect_results['hera_improvement']['beta_improvement_factor']),
        },
        "torino_palermo": scales,
        "elapsed_times_s": {
            "monte_carlo": elapsed_mc,
            "perturbations": elapsed_perturb,
            "keyholes": elapsed_keyhole,
            "bayesian": elapsed_bayes,
            "damage": elapsed_damage,
            "deflection": elapsed_deflect,
        }
    }

    with open("results/neo_risk_assessment_summary.json", "w") as f:
        json.dump(results_summary, f, indent=2)

    # Detailed MC distance array
    np.save("results/mc_min_distances_au.npy", mc_results['min_distances_au'])

    # Observation campaign history CSV
    import csv
    with open("results/bayesian_observation_history.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=['obs_number', 'epoch', 'impact_prob', 'a_std', 'e_std'])
        writer.writeheader()
        writer.writerows(obs_history)

    print("\n[Results] Saved to results/neo_risk_assessment_summary.json")

    log_event("pipeline", "run_completed", "neo_pipeline",
               files_written=[
                   "results/neo_risk_assessment_summary.json",
                   "results/mc_min_distances_au.npy",
                   "results/bayesian_observation_history.csv",
               ],
               handoff_out=results_summary,
               status="ok")

    print("\n" + "="*60)
    print(" PIPELINE COMPLETE")
    print("="*60)
    print(f"  Impact probability (MC):    {mc_results['impact_probability']:.2e}")
    print(f"  Posterior (30 obs):         {final_prob:.2e}")
    print(f"  Yarkovsky drift (100 yr):   {yark_band['mean']:.2e} AU")
    print(f"  Keyholes found:             {len(keyholes)}")
    print(f"  Kinetic energy:             {damage.kinetic_energy_Mt:.1f} Mt TNT")
    print(f"  DART deflection success:    {deflect_results['dart']['success_probability']:.1%}")
    print(f"  Torino / Palermo scale:     {scales['torino_scale']} / {scales['palermo_scale']:.2f}")

    return results_summary


def _plot_entry_trajectory(entry: dict, save_path: str) -> None:
    """Plot atmospheric entry trajectory."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    t = entry['t']
    v = entry['v_trajectory_km_s']
    h = entry['h_trajectory_km']
    mask = h > 0

    axes[0].plot(t[mask], v[mask], color='crimson', linewidth=1.5)
    axes[0].set_xlabel('Time [s]')
    axes[0].set_ylabel('Velocity [km/s]')
    axes[0].set_title('Velocity During Atmospheric Entry')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(v[mask], h[mask], color='steelblue', linewidth=1.5)
    axes[1].set_xlabel('Velocity [km/s]')
    axes[1].set_ylabel('Altitude [km]')
    axes[1].set_title('Velocity-Altitude Profile')
    if entry['airburst_altitude_m']:
        ab_km = entry['airburst_altitude_m'] / 1e3
        axes[1].axhline(ab_km, color='orange', linestyle='--',
                         label=f'Airburst: {ab_km:.1f} km')
        axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Plot] Saved: {save_path}")


def _plot_summary_dashboard(mc_results, perturb_results, yark_band, keyholes,
                              obs_history, damage, deflect_results, scales,
                              save_path: str) -> None:
    """Create 6-panel summary figure."""
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    # Panel 1: MC distance distribution
    ax1 = fig.add_subplot(gs[0, 0])
    dists = mc_results['min_distances_au'] / 0.00257
    finite = dists[dists < 500]
    ax1.hist(finite, bins=40, color='steelblue', edgecolor='white', alpha=0.8)
    ax1.axvline(1.0, color='gold', linestyle='--', label='Moon')
    ax1.set_xlabel('Min Earth Distance [LD]')
    ax1.set_ylabel('Count')
    ax1.set_title(f'MC Approach Distribution\nP_impact={mc_results["impact_probability"]:.2e}')
    ax1.legend(fontsize=8)

    # Panel 2: Orbital evolution (a)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(perturb_results['times'], perturb_results['a'], color='crimson', linewidth=0.8)
    ax2.set_xlabel('Time [yr]')
    ax2.set_ylabel('Semi-major axis [AU]')
    ax2.set_title(f'Secular Evolution + Yarkovsky\nda/dt={perturb_results.get("yarkovsky_rate_au_per_yr", 0):.2e} AU/yr')

    # Panel 3: Bayesian update
    ax3 = fig.add_subplot(gs[0, 2])
    obs_nums = [h['obs_number'] for h in obs_history]
    probs = [h['impact_prob'] for h in obs_history]
    ax3.semilogy(obs_nums, probs, 'o-', color='seagreen', markersize=4)
    ax3.set_xlabel('Observation Number')
    ax3.set_ylabel('Impact Probability')
    ax3.set_title(f'Bayesian Update (30 obs)\nFinal P={probs[-1]:.2e}')
    ax3.grid(True, alpha=0.3)

    # Panel 4: Top keyholes
    ax4 = fig.add_subplot(gs[1, 0])
    if keyholes:
        top = keyholes[:8]
        labels = [f"{k.resonance_n}:{k.resonance_m}" for k in top]
        kh_probs = [k.impact_prob for k in top]
        colors = plt.cm.YlOrRd(np.linspace(0.4, 1.0, len(top)))
        ax4.barh(range(len(labels)), kh_probs, color=colors)
        ax4.set_yticks(range(len(labels)))
        ax4.set_yticklabels(labels, fontsize=9)
        ax4.set_xscale('log')
        ax4.invert_yaxis()
        ax4.set_xlabel('Impact Probability')
        ax4.set_title(f'Top Keyholes ({len(keyholes)} found)')

    # Panel 5: Damage zones
    ax5 = fig.add_subplot(gs[1, 1])
    radii_info = [
        ('Total destruction', damage.blast_radius_km['total_destruction_1_4_bar'], '#8B0000'),
        ('Severe damage',     damage.blast_radius_km['severe_damage_0_6_bar'], '#DC143C'),
        ('Thermal burns',     damage.thermal_radius_km, '#FF8C00'),
        ('Window breakage',   damage.blast_radius_km['window_breakage_0_007_bar'], '#FFD700'),
    ]
    labels_d = [r[0] for r in radii_info]
    vals_d = [r[1] for r in radii_info]
    colors_d = [r[2] for r in radii_info]
    ax5.barh(labels_d, vals_d, color=colors_d, alpha=0.8, edgecolor='white')
    ax5.set_xlabel('Radius [km]')
    ax5.set_title(f'Impact Damage Zones\nKE={damage.kinetic_energy_Mt:.1f} Mt, T={damage.torino_scale}')

    # Panel 6: Deflection comparison
    ax6 = fig.add_subplot(gs[1, 2])
    missions = ['Kinetic\nImpactor\n(DART-like)', 'Gravity\nTractor', 'Combined']
    dvs = [deflect_results['dart']['mean_delta_v_m_s'],
            deflect_results['gravity_tractor']['delta_v_m_s'],
            deflect_results['combined_dv_m_s']]
    colors_df = ['#1f77b4', '#ff7f0e', '#2ca02c']
    bars = ax6.bar(missions, dvs, color=colors_df, edgecolor='white')
    for bar, dv in zip(bars, dvs):
        ax6.text(bar.get_x() + bar.get_width()/2, dv * 1.02,
                  f'{dv:.3f}', ha='center', va='bottom', fontsize=8)
    ax6.set_ylabel('Δv Imparted [m/s]')
    ax6.set_title(f'Deflection Mission Performance\nDART P(success)={deflect_results["dart"]["success_probability"]:.1%}')

    fig.suptitle('NEO Planetary Defense Risk Assessment — Summary Dashboard',
                  fontsize=15, fontweight='bold', y=1.01)
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Plot] Saved: {save_path}")


if __name__ == "__main__":
    results = run_pipeline(
        n_clones=300,
        t_years=100.0,
        run_mcmc=False,  # Set True for full MCMC (~3 min)
        seed=42
    )
