"""
Module 5: Population Genetics Model for Local Adaptation / Evolutionary Response
Wright-Fisher model with selection under changing ocean acidification conditions.
"""

import numpy as np
import json


def wright_fisher_simulation(N_e, n_loci, generations, selection_coefficients,
                              mutation_rate=1e-5, recombination_rate=0.01,
                              initial_freq=None, env_change_func=None,
                              seed=42):
    """
    Wright-Fisher forward simulation for coral population adaptation.
    
    Parameters
    ----------
    N_e : int - Effective population size
    n_loci : int - Number of loci under selection
    generations : int - Number of generations
    selection_coefficients : array - Base selection coefficients per locus
    mutation_rate : float - Per-locus mutation rate
    recombination_rate : float - Recombination rate between adjacent loci
    initial_freq : array - Initial allele frequencies
    env_change_func : callable - Function(gen) returning environmental stress multiplier
    seed : int - Random seed
    
    Returns
    -------
    freq_history : array (generations x n_loci) - Allele frequency trajectories
    fitness_history : array - Mean population fitness over time
    """
    np.random.seed(seed)
    
    if initial_freq is None:
        initial_freq = np.random.uniform(0.01, 0.15, n_loci)
    
    freq = initial_freq.copy()
    freq_history = np.zeros((generations, n_loci))
    fitness_history = np.zeros(generations)
    Va_history = np.zeros(generations)  # additive genetic variance
    
    s = np.array(selection_coefficients)
    
    for gen in range(generations):
        # Environmental change modifies selection pressure
        if env_change_func is not None:
            env_factor = env_change_func(gen)
            s_effective = s * env_factor
        else:
            s_effective = s.copy()
        
        # Fitness calculation (multiplicative model)
        # W = product over loci of (1 + s_i * freq_i)
        mean_fitness = np.prod(1 + s_effective * freq)
        
        # Selection: adjust frequencies
        for i in range(n_loci):
            w_AA = (1 + s_effective[i])**2
            w_Aa = 1 + s_effective[i]
            w_aa = 1.0
            p = freq[i]
            q = 1 - p
            w_bar = p**2 * w_AA + 2*p*q * w_Aa + q**2 * w_aa
            if w_bar > 0:
                freq[i] = (p**2 * w_AA + p*q * w_Aa) / w_bar
        
        # Mutation
        for i in range(n_loci):
            freq[i] = freq[i] * (1 - mutation_rate) + (1 - freq[i]) * mutation_rate
        
        # Genetic drift (binomial sampling)
        for i in range(n_loci):
            freq[i] = np.random.binomial(2 * N_e, freq[i]) / (2 * N_e)
        
        freq = np.clip(freq, 0, 1)
        
        # Additive genetic variance
        Va = np.sum(2 * freq * (1-freq) * s_effective**2)
        
        freq_history[gen] = freq
        fitness_history[gen] = mean_fitness
        Va_history[gen] = Va
    
    return freq_history, fitness_history, Va_history


def coral_adaptation_model(scenario='RCP8.5', N_e=5000, n_generations=80):
    """
    Model coral adaptation to ocean acidification.
    
    Coral generation time ~5-7 years for Acropora, ~10-15 for Porites.
    By 2100, this corresponds to ~12-16 Acropora generations.
    """
    # OA-tolerance loci (hypothetical QTLs)
    # Based on heritability estimates from Matz et al. (2018)
    n_loci = 20
    
    # Selection coefficients: effect sizes of tolerance alleles
    # Larger effects for major QTLs, smaller for minor ones
    s_base = np.concatenate([
        np.array([0.08, 0.06, 0.05, 0.04, 0.04]),   # 5 major QTLs
        np.array([0.02] * 5),                          # 5 moderate QTLs
        np.array([0.01] * 10)                          # 10 minor QTLs
    ])
    
    # Environmental stress scenarios
    env_functions = {
        'RCP2.6': lambda gen: 1.0 + 0.3 * np.sin(gen / n_generations * np.pi),
        'RCP4.5': lambda gen: 1.0 + 1.5 * (gen / n_generations),
        'RCP8.5': lambda gen: 1.0 + 4.0 * (gen / n_generations)**1.5,
    }
    
    # Species-specific parameters
    species_configs = {
        'Acropora': {
            'N_e': N_e,
            'gen_time_years': 5,
            'heritability_h2': 0.25,
            'standing_variation': 0.10,
            'mutation_rate': 2e-5,
        },
        'Porites': {
            'N_e': N_e * 2,
            'gen_time_years': 12,
            'heritability_h2': 0.35,
            'standing_variation': 0.15,
            'mutation_rate': 1e-5,
        },
        'Stylophora': {
            'N_e': N_e // 2,
            'gen_time_years': 4,
            'heritability_h2': 0.20,
            'standing_variation': 0.08,
            'mutation_rate': 2e-5,
        }
    }
    
    results = {}
    for sp_name, config in species_configs.items():
        sp_results = {}
        n_gen = int(80 / config['gen_time_years'])  # generations by 2100
        
        for scen_name, env_func in env_functions.items():
            initial_freq = np.random.uniform(
                0.01, config['standing_variation'], n_loci
            )
            
            freq_hist, fit_hist, Va_hist = wright_fisher_simulation(
                N_e=config['N_e'],
                n_loci=n_loci,
                generations=n_gen,
                selection_coefficients=s_base,
                mutation_rate=config['mutation_rate'],
                initial_freq=initial_freq,
                env_change_func=env_func,
                seed=42
            )
            
            # Compute adaptation metrics
            mean_freq_initial = np.mean(freq_hist[0])
            mean_freq_final = np.mean(freq_hist[-1])
            delta_freq = mean_freq_final - mean_freq_initial
            
            # Breeder's equation: R = h² × S
            # Selection differential from environmental pressure
            S = env_func(n_gen) * np.mean(s_base)
            R = config['heritability_h2'] * S
            
            # Adaptation lag: compare rate of environmental change vs adaptation
            env_change_rate = (env_func(n_gen) - env_func(0)) / n_gen
            adaptation_rate = delta_freq / n_gen if n_gen > 0 else 0
            adaptation_lag = env_change_rate / max(adaptation_rate, 1e-10)
            
            sp_results[scen_name] = {
                'n_generations': n_gen,
                'initial_mean_freq': round(float(mean_freq_initial), 4),
                'final_mean_freq': round(float(mean_freq_final), 4),
                'delta_freq': round(float(delta_freq), 4),
                'initial_fitness': round(float(fit_hist[0]), 4),
                'final_fitness': round(float(fit_hist[-1]), 4),
                'fitness_change_pct': round(float((fit_hist[-1]/fit_hist[0] - 1)*100), 2),
                'Va_initial': round(float(Va_hist[0]), 6),
                'Va_final': round(float(Va_hist[-1]), 6),
                'breeder_response_R': round(float(R), 4),
                'adaptation_lag': round(float(adaptation_lag), 2),
                'freq_trajectory': freq_hist.mean(axis=1).tolist(),
                'fitness_trajectory': fit_hist.tolist(),
                'Va_trajectory': Va_hist.tolist(),
            }
        
        results[sp_name] = {
            'config': {k: v for k, v in config.items() if not callable(v)},
            'scenarios': sp_results
        }
    
    # Evolutionary rescue analysis
    rescue_analysis = {}
    for sp_name in species_configs:
        sp = species_configs[sp_name]
        rcp85 = results[sp_name]['scenarios']['RCP8.5']
        
        # Can adaptation keep pace with environmental change?
        can_rescue = rcp85['adaptation_lag'] < 5.0  # threshold
        rescue_analysis[sp_name] = {
            'adaptation_lag': rcp85['adaptation_lag'],
            'evolutionary_rescue_possible': can_rescue,
            'fitness_retained_pct': round(rcp85['final_fitness'] / rcp85['initial_fitness'] * 100, 1),
            'Va_retained_pct': round(rcp85['Va_final'] / max(rcp85['Va_initial'], 1e-10) * 100, 1),
        }
    
    all_results = {
        'species_results': {k: {kk: vv for kk, vv in v.items()} for k, v in results.items()},
        'evolutionary_rescue': rescue_analysis,
        'model_parameters': {
            'n_loci': n_loci,
            'selection_coefficients': s_base.tolist(),
            'n_generations_to_2100': {sp: int(80/c['gen_time_years']) 
                                       for sp, c in species_configs.items()},
        }
    }
    
    # Remove trajectories for JSON (too large), save separately
    for sp in all_results['species_results']:
        for scen in all_results['species_results'][sp]['scenarios']:
            for key in ['freq_trajectory', 'fitness_trajectory', 'Va_trajectory']:
                all_results['species_results'][sp]['scenarios'][scen].pop(key, None)
    
    with open('results/popgen_model.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    # Save trajectories separately
    traj_data = {}
    for sp_name in species_configs:
        for scen_name in env_functions:
            key = f"{sp_name}_{scen_name}".replace('.', '')
            sp_r = results[sp_name]['scenarios'][scen_name]
            # Re-run to get trajectories
            config = species_configs[sp_name]
            n_gen = int(80 / config['gen_time_years'])
            initial_freq = np.random.uniform(0.01, config['standing_variation'], n_loci)
            freq_hist, fit_hist, Va_hist = wright_fisher_simulation(
                N_e=config['N_e'], n_loci=n_loci, generations=n_gen,
                selection_coefficients=s_base, mutation_rate=config['mutation_rate'],
                initial_freq=initial_freq, env_change_func=env_functions[scen_name],
                seed=42
            )
            traj_data[f"{key}_freq"] = freq_hist.mean(axis=1)
            traj_data[f"{key}_fitness"] = fit_hist
            traj_data[f"{key}_Va"] = Va_hist
    
    np.savez('results/popgen_trajectories.npz', **traj_data)
    
    return all_results


if __name__ == '__main__':
    results = coral_adaptation_model()
    print("=== Population Genetics Model Results ===")
    print("\nEvolutionary Rescue Assessment:")
    for sp, data in results['evolutionary_rescue'].items():
        print(f"  {sp}: rescue_possible={data['evolutionary_rescue_possible']}, "
              f"fitness_retained={data['fitness_retained_pct']:.1f}%, "
              f"adaptation_lag={data['adaptation_lag']:.1f}")
    print("\nResults saved to results/popgen_model.json")
