"""
Module 6: Astrobiology — Chemical Evolution on Enceladus and Titan
Environmental parameter models and chemical evolution feasibility assessment.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

np.random.seed(202)

# --- Environmental Parameters ---
ENVIRONMENTS = {
    'Early Earth': {
        'temp_K': 353,        # ~80°C
        'pressure_atm': 1.0,
        'pH': 7.0,
        'solvent': 'H2O',
        'energy_sources': ['UV', 'Lightning', 'Volcanic'],
        'reducing': True,
        'available_elements': ['C', 'H', 'O', 'N', 'S', 'P', 'Fe'],
        'surface_gravity_ms2': 9.8,
    },
    'Enceladus Ocean': {
        'temp_K': 275,        # ~2°C subsurface ocean
        'pressure_atm': 10,   # estimated
        'pH': 9.5,            # alkaline
        'solvent': 'H2O',
        'energy_sources': ['Hydrothermal', 'Radiolysis'],
        'reducing': True,
        'available_elements': ['C', 'H', 'O', 'N', 'S'],
        'surface_gravity_ms2': 0.113,
    },
    'Enceladus Vent': {
        'temp_K': 363,        # ~90°C hydrothermal
        'pressure_atm': 50,
        'pH': 10.0,
        'solvent': 'H2O',
        'energy_sources': ['Hydrothermal', 'Serpentinization'],
        'reducing': True,
        'available_elements': ['C', 'H', 'O', 'N', 'S', 'Fe', 'Ni'],
        'surface_gravity_ms2': 0.113,
    },
    'Titan Surface': {
        'temp_K': 94,         # ~-179°C
        'pressure_atm': 1.45,
        'pH': None,           # no water
        'solvent': 'CH4/C2H6',
        'energy_sources': ['UV (weak)', 'Cosmic rays'],
        'reducing': True,
        'available_elements': ['C', 'H', 'N'],
        'surface_gravity_ms2': 1.352,
    },
    'Titan Subsurface': {
        'temp_K': 255,
        'pressure_atm': 100,
        'pH': 8.0,
        'solvent': 'H2O/NH3',
        'energy_sources': ['Radioactive decay', 'Tidal'],
        'reducing': True,
        'available_elements': ['C', 'H', 'O', 'N', 'S'],
        'surface_gravity_ms2': 1.352,
    },
}

def arrhenius_rate(T, Ea=80000, A=1e10):
    """Arrhenius rate constant. Ea in J/mol."""
    R = 8.314
    return A * np.exp(-Ea / (R * T))

def reaction_feasibility_score(env):
    """
    Compute a composite feasibility score for prebiotic chemistry.
    Factors: temperature, energy, elemental diversity, solvent polarity.
    """
    T = env['temp_K']
    score = 0.0

    # Temperature factor (Arrhenius-like, normalised)
    k_rel = arrhenius_rate(T) / arrhenius_rate(353)  # relative to early Earth
    score += min(k_rel, 2.0) * 25  # max 50 points

    # Energy source diversity
    score += len(env['energy_sources']) * 8  # max ~24

    # Elemental diversity
    essential = {'C', 'H', 'O', 'N', 'S', 'P'}
    available = set(env['available_elements'])
    score += len(available & essential) / len(essential) * 20  # max 20

    # Solvent (water = best)
    if 'H2O' in env['solvent']:
        score += 15
    elif 'NH3' in env['solvent']:
        score += 8
    else:
        score += 3  # hydrocarbon solvent

    # pH factor (near neutral to mildly alkaline is best)
    if env['pH'] is not None:
        pH_opt = 1.0 - abs(env['pH'] - 8.5) / 5.0
        score += max(0, pH_opt * 10)

    return min(score, 100)

def simulate_chemical_kinetics_comparison(environments, t_span=500):
    """Compare reaction kinetics across environments."""
    t = np.linspace(0, t_span, 1000)
    kinetics = {}

    for name, env in environments.items():
        T = env['temp_K']
        k = arrhenius_rate(T, Ea=80000)

        # Simple A -> B -> C kinetics
        k1 = k * 1e-8  # normalise
        k2 = k * 5e-9

        A = 100 * np.exp(-k1 * t)
        B = 100 * k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t)) if abs(k2 - k1) > 1e-20 \
            else 100 * k1 * t * np.exp(-k1 * t)
        C = 100 - A - B

        kinetics[name] = {'t': t.tolist(), 'A': A.tolist(), 'B': B.tolist(), 'C': C.tolist(),
                          'k1': k1, 'k2': k2}

    return kinetics

def monte_carlo_habitability(environments, n_samples=10000):
    """Monte Carlo sampling of parameter uncertainties for habitability."""
    results = {}
    for name, env in environments.items():
        scores = []
        for _ in range(n_samples):
            env_perturbed = dict(env)
            env_perturbed['temp_K'] = env['temp_K'] + np.random.normal(0, env['temp_K'] * 0.1)
            if env['pH'] is not None:
                env_perturbed['pH'] = env['pH'] + np.random.normal(0, 0.5)
            scores.append(reaction_feasibility_score(env_perturbed))
        results[name] = {
            'mean': float(np.mean(scores)),
            'std': float(np.std(scores)),
            'p95': float(np.percentile(scores, 95)),
            'p05': float(np.percentile(scores, 5)),
        }
    return results

def plot_astrobiology(environments, kinetics, mc_results):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Astrobiology: Chemical Evolution Feasibility Across Worlds', fontsize=14)

    # A) Feasibility scores
    ax = axes[0, 0]
    names = list(environments.keys())
    scores = [reaction_feasibility_score(environments[n]) for n in names]
    colors = plt.cm.viridis(np.array(scores) / 100)
    bars = ax.barh(names, scores, color=colors)
    ax.set_xlabel('Feasibility Score (0-100)')
    ax.set_title('A) Prebiotic Chemistry Feasibility')
    ax.set_xlim(0, 100)
    for bar, s in zip(bars, scores):
        ax.text(s + 1, bar.get_y() + bar.get_height()/2, f'{s:.0f}', va='center')

    # B) Kinetics comparison (product C formation)
    ax = axes[0, 1]
    for name in kinetics:
        t = kinetics[name]['t']
        C = kinetics[name]['C']
        ax.plot(t, C, label=name)
    ax.set_xlabel('Time')
    ax.set_ylabel('[Product]')
    ax.set_title('B) Product Formation Kinetics')
    ax.legend(fontsize=7)

    # C) Temperature vs rate
    ax = axes[1, 0]
    temps = np.linspace(80, 400, 200)
    rates = [arrhenius_rate(T) for T in temps]
    ax.semilogy(temps - 273.15, rates, color='#F44336')
    for name, env in environments.items():
        T = env['temp_K']
        ax.axvline(T - 273.15, alpha=0.5, linestyle='--', label=name)
        ax.annotate(name, (T - 273.15, arrhenius_rate(T)),
                    fontsize=6, rotation=45)
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Rate constant (s⁻¹)')
    ax.set_title('C) Arrhenius Rate vs Temperature')

    # D) Monte Carlo habitability
    ax = axes[1, 1]
    names_mc = list(mc_results.keys())
    means = [mc_results[n]['mean'] for n in names_mc]
    stds = [mc_results[n]['std'] for n in names_mc]
    ax.barh(names_mc, means, xerr=stds, color='#9C27B0', alpha=0.7, capsize=3)
    ax.set_xlabel('Feasibility Score (Monte Carlo)')
    ax.set_title('D) Habitability (with Uncertainty)')
    ax.set_xlim(0, 100)

    plt.tight_layout()
    plt.savefig('figures/fig9_astrobiology.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig9_astrobiology.svg', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("Running Astrobiology assessment...")

    kinetics = simulate_chemical_kinetics_comparison(ENVIRONMENTS)
    mc_results = monte_carlo_habitability(ENVIRONMENTS)

    plot_astrobiology(ENVIRONMENTS, kinetics, mc_results)

    # Summary
    results = {}
    for name, env in ENVIRONMENTS.items():
        results[name] = {
            'feasibility_score': reaction_feasibility_score(env),
            'temperature_K': env['temp_K'],
            'solvent': env['solvent'],
            'mc_mean': mc_results[name]['mean'],
            'mc_std': mc_results[name]['std'],
            'arrhenius_rate': float(arrhenius_rate(env['temp_K'])),
        }
    with open('results/astrobiology_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("Module 6 complete.")
    for name, r in results.items():
        print(f"  {name}: score={r['feasibility_score']:.0f}, MC={r['mc_mean']:.1f}±{r['mc_std']:.1f}")
