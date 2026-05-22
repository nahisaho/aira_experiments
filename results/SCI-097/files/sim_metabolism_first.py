"""
Module 3: Metabolism-First Hypothesis — Hydrothermal Vent Model
ODE + stochastic simulation of autocatalytic metabolic cycles at alkaline vents.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

np.random.seed(456)

# --- Reverse TCA cycle simplified model ---
# Species: CO2, H2, Acetate, Pyruvate, Oxaloacetate, Citrate, FeS_catalyst
SPECIES = ['CO2', 'H2', 'Acetate', 'Pyruvate', 'Oxaloacetate',
           'Citrate', 'FeS', 'Succinate', 'Fumarate', 'Malate']

def vent_ode(t, y, params):
    """ODE for simplified reverse TCA at hydrothermal vent."""
    CO2, H2, Ac, Pyr, OAA, Cit, FeS, Suc, Fum, Mal = y
    k = params

    # CO2 fixation: CO2 + H2 -> Acetate (FeS catalysed)
    r1 = k['k1'] * CO2 * H2 * FeS / (k['Km1'] + CO2)
    # Acetate -> Pyruvate
    r2 = k['k2'] * Ac * FeS / (k['Km2'] + Ac)
    # Pyruvate + CO2 -> Oxaloacetate
    r3 = k['k3'] * Pyr * CO2 / (k['Km3'] + Pyr)
    # OAA -> Citrate (condensation with Acetate)
    r4 = k['k4'] * OAA * Ac / (k['Km4'] + OAA)
    # Citrate -> Succinate + CO2 (reductive)
    r5 = k['k5'] * Cit * H2 / (k['Km5'] + Cit)
    # Succinate -> Fumarate
    r6 = k['k6'] * Suc / (k['Km6'] + Suc)
    # Fumarate -> Malate
    r7 = k['k7'] * Fum * H2 / (k['Km7'] + Fum)
    # Malate -> OAA (closing the cycle)
    r8 = k['k8'] * Mal / (k['Km8'] + Mal)

    # Influx from vent
    influx_CO2 = k['influx_CO2']
    influx_H2 = k['influx_H2']
    # Dilution
    d = k['dilution']

    dydt = [
        influx_CO2 - r1 - r3 + r5 - d * CO2,      # CO2
        influx_H2 - r1 - r5 - r7 - d * H2,         # H2
        r1 - r2 - r4 - d * Ac,                       # Acetate
        r2 - r3 - d * Pyr,                           # Pyruvate
        r3 + r8 - r4 - d * OAA,                      # Oxaloacetate
        r4 - r5 - d * Cit,                           # Citrate
        -0.001 * FeS,                                 # FeS (slow degradation)
        r5 - r6 - d * Suc,                           # Succinate
        r6 - r7 - d * Fum,                           # Fumarate
        r7 - r8 - d * Mal,                           # Malate
    ]
    return dydt

def run_vent_simulation(temp_celsius=90.0, pH=9.0, t_span=(0, 500)):
    """Run ODE for given temperature and pH (affects rate constants)."""
    temp_factor = np.exp(-30000 / (8.314 * (temp_celsius + 273.15))) / \
                  np.exp(-30000 / (8.314 * 363.15))  # normalise to 90°C
    pH_factor = 1.0 / (1 + 10**(pH - 10))  # alkaline optimum ~10

    base_k = 0.1 * temp_factor * (1 + pH_factor)

    params = {
        'k1': base_k * 2.0, 'Km1': 10.0,
        'k2': base_k * 0.5, 'Km2': 5.0,
        'k3': base_k * 0.3, 'Km3': 8.0,
        'k4': base_k * 0.4, 'Km4': 6.0,
        'k5': base_k * 0.2, 'Km5': 7.0,
        'k6': base_k * 0.6, 'Km6': 4.0,
        'k7': base_k * 0.3, 'Km7': 5.0,
        'k8': base_k * 0.5, 'Km8': 3.0,
        'influx_CO2': 50.0,
        'influx_H2': 80.0,
        'dilution': 0.01,
    }

    y0 = [100, 200, 0, 0, 0, 0, 50, 0, 0, 0]  # FeS=50 as mineral catalyst
    sol = solve_ivp(vent_ode, t_span, y0, args=(params,),
                    method='LSODA', max_step=2.0, rtol=1e-3, atol=1e-5,
                    t_eval=np.linspace(t_span[0], t_span[1], 300))
    return sol

def scan_temperature_pH():
    """Scan T and pH for metabolic cycle turnover."""
    temps = np.linspace(50, 140, 8)
    pHs = np.linspace(7, 11, 8)
    turnover = np.zeros((len(temps), len(pHs)))

    for i, T in enumerate(temps):
        for j, pH in enumerate(pHs):
            try:
                sol = run_vent_simulation(temp_celsius=T, pH=pH, t_span=(0, 50))
                if sol.success and sol.y.shape[1] > 1:
                    turnover[i, j] = np.trapz(sol.y[5], sol.t)
                else:
                    turnover[i, j] = 0.0
            except Exception:
                turnover[i, j] = 0.0

    return temps, pHs, turnover

def plot_vent_dynamics(sol):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Metabolism-First: Hydrothermal Vent rTCA Cycle', fontsize=14)

    ax = axes[0, 0]
    ax.plot(sol.t, sol.y[0], label='CO₂', color='#607D8B')
    ax.plot(sol.t, sol.y[1], label='H₂', color='#03A9F4')
    ax.set_xlabel('Time')
    ax.set_ylabel('Concentration')
    ax.set_title('A) Inorganic Substrates')
    ax.legend()

    ax = axes[0, 1]
    ax.plot(sol.t, sol.y[2], label='Acetate')
    ax.plot(sol.t, sol.y[3], label='Pyruvate')
    ax.set_xlabel('Time')
    ax.set_ylabel('Concentration')
    ax.set_title('B) C2-C3 Metabolites')
    ax.legend()

    ax = axes[1, 0]
    for idx, name in [(4, 'OAA'), (5, 'Citrate'), (7, 'Succinate'), (8, 'Fumarate'), (9, 'Malate')]:
        ax.plot(sol.t, sol.y[idx], label=name)
    ax.set_xlabel('Time')
    ax.set_ylabel('Concentration')
    ax.set_title('C) TCA Cycle Intermediates')
    ax.legend()

    ax = axes[1, 1]
    ax.plot(sol.t, sol.y[6], label='FeS Catalyst', color='#795548')
    ax.set_xlabel('Time')
    ax.set_ylabel('Concentration')
    ax.set_title('D) Mineral Catalyst (FeS)')
    ax.legend()

    plt.tight_layout()
    plt.savefig('figures/fig5_metabolism_first.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig5_metabolism_first.svg', bbox_inches='tight')
    plt.close()

def plot_turnover_heatmap(temps, pHs, turnover):
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(turnover, origin='lower', aspect='auto', cmap='inferno',
                   extent=[pHs[0], pHs[-1], temps[0], temps[-1]])
    ax.set_xlabel('pH')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title('rTCA Cycle Turnover: Temperature–pH Phase Space')
    plt.colorbar(im, ax=ax, label='Cumulative Citrate (∫[Cit]dt)')
    plt.tight_layout()
    plt.savefig('figures/fig6_vent_phase_space.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("Running Metabolism-First simulation...")
    sol = run_vent_simulation()
    plot_vent_dynamics(sol)

    print("Scanning T-pH parameter space...")
    temps, pHs, turnover = scan_temperature_pH()
    plot_turnover_heatmap(temps, pHs, turnover)

    # Optimal conditions
    opt_idx = np.unravel_index(np.argmax(turnover), turnover.shape)
    results = {
        'final_concentrations': {SPECIES[i]: float(sol.y[i, -1]) for i in range(10)},
        'optimal_temperature_C': float(temps[opt_idx[0]]),
        'optimal_pH': float(pHs[opt_idx[1]]),
        'max_turnover': float(turnover[opt_idx]),
    }
    with open('results/metabolism_first_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Module 3 complete. Optimal: T={results['optimal_temperature_C']:.0f}°C, pH={results['optimal_pH']:.1f}")
