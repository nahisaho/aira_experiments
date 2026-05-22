"""
Module 6: Fischer-Tropsch Synthesis Case Study
===============================================
Microkinetic model for CO hydrogenation on Co(0001) surface.

Mechanism (carbide mechanism):
  1. CO  + *   → CO*           (CO adsorption)
  2. H2  + 2*  → 2H*           (H2 dissociative adsorption)
  3. CO* + H*  → HCO* + *      (CO hydrogenation)
  4. HCO* + H* → CH2O* + *     (HCO hydrogenation)
  5. CH2O*     → CH2* + O*     (C-O bond scission)
  6. CH2* + H* → CH3* + *      (CH2 hydrogenation)
  7. CH3* + H* → CH4 + 2*      (Methane formation)
  8. O*  + H*  → OH* + *       (O removal step 1)
  9. OH* + H*  → H2O + 2*      (O removal step 2 → water)
  10. CH2* + CH2* → C2H4 + 2*  (Chain growth / C2 formation, simplified)
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .rate_constants import TransitionState, calculate_tst_rate, RateConstant, arrhenius_parameters
from .adsorption import (AdsorptionParameters, competitive_langmuir,
                         sticking_rate, langmuir_isotherm, temkin_isotherm)
from .lateral import (LateralInteractionParams, mean_field_interaction_energy,
                      modified_rate_constants, solve_coverage_self_consistent)
from .rds_identifier import degree_of_rate_control, energy_span_analysis
from .reactor import (ReactorConditions, PFRReactor, CSTRReactor,
                      solve_steady_state_coverages)


# ============================================================
# DFT-derived parameters for Co(0001) FT synthesis
# (Representative values from literature: Zhuo et al., JACS 2009;
#  Ojeda et al., J. Catal. 2010; van Santen et al., PCCP 2011)
# ============================================================

FT_TRANSITION_STATES = [
    TransitionState("CO_adsorption", E_activation=0.0, E_reaction=-1.30,
                    frequencies_real=[350, 380, 1800], frequency_imaginary=0, symmetry_number=1),
    TransitionState("H2_dissociation", E_activation=0.05, E_reaction=-0.50,
                    frequencies_real=[800, 900, 3100], frequency_imaginary=850, symmetry_number=2),
    TransitionState("CO_hydrogenation", E_activation=0.80, E_reaction=0.20,
                    frequencies_real=[300, 450, 500, 1200, 1600], frequency_imaginary=420, symmetry_number=1),
    TransitionState("HCO_hydrogenation", E_activation=0.55, E_reaction=-0.10,
                    frequencies_real=[280, 400, 600, 1100, 1400, 2800], frequency_imaginary=380, symmetry_number=1),
    TransitionState("C-O_scission", E_activation=1.20, E_reaction=-0.40,
                    frequencies_real=[250, 350, 500, 700], frequency_imaginary=550, symmetry_number=1),
    TransitionState("CH2_hydrogenation", E_activation=0.60, E_reaction=-0.30,
                    frequencies_real=[400, 600, 800, 1300, 2900], frequency_imaginary=450, symmetry_number=1),
    TransitionState("CH4_formation", E_activation=0.95, E_reaction=0.10,
                    frequencies_real=[500, 700, 1000, 1350, 2950, 3000], frequency_imaginary=520, symmetry_number=3),
    TransitionState("O_hydrogenation", E_activation=0.70, E_reaction=0.15,
                    frequencies_real=[350, 500, 700], frequency_imaginary=480, symmetry_number=1),
    TransitionState("OH_hydrogenation", E_activation=0.85, E_reaction=-0.20,
                    frequencies_real=[400, 600, 800, 3600], frequency_imaginary=510, symmetry_number=1),
    TransitionState("C2_coupling", E_activation=0.75, E_reaction=-0.50,
                    frequencies_real=[300, 400, 500, 600, 700, 800, 1200, 2900], frequency_imaginary=350, symmetry_number=2),
]

FT_ADSORPTION_PARAMS = {
    'CO': AdsorptionParameters('CO', delta_H_ads=-1.30, delta_S_ads=-1.5e-3, site_type='top'),
    'H':  AdsorptionParameters('H',  delta_H_ads=-0.50, delta_S_ads=-0.8e-3, site_type='hollow'),
    'O':  AdsorptionParameters('O',  delta_H_ads=-2.00, delta_S_ads=-1.2e-3, site_type='hollow'),
}

FT_LATERAL_INTERACTIONS = [
    LateralInteractionParams(('CO', 'CO'), epsilon_nn=-0.10, z_nn=6),
    LateralInteractionParams(('CO', 'H'),  epsilon_nn=-0.02, z_nn=6),
    LateralInteractionParams(('H',  'H'),  epsilon_nn=-0.01, z_nn=6),
    LateralInteractionParams(('O',  'CO'), epsilon_nn=-0.05, z_nn=6),
    LateralInteractionParams(('O',  'O'),  epsilon_nn=-0.15, z_nn=6),
]

# Surface species: CO*, H*, HCO*, CH2O*, CH2*, CH3*, O*, OH*
FT_SURFACE_SPECIES = ['CO*', 'H*', 'HCO*', 'CH2O*', 'CH2*', 'CH3*', 'O*', 'OH*']
FT_GAS_SPECIES = ['CO', 'H2', 'CH4', 'H2O', 'C2H4']


def build_stoichiometric_matrices():
    """
    Build stoichiometric matrices for gas and surface species.

    Returns
    -------
    stoich_gas : np.ndarray [n_gas x n_rxns]
    stoich_surface : np.ndarray [n_surface x n_rxns]
    """
    n_rxns = 10
    n_gas = 5    # CO, H2, CH4, H2O, C2H4
    n_surf = 8   # CO*, H*, HCO*, CH2O*, CH2*, CH3*, O*, OH*

    # Gas: [CO, H2, CH4, H2O, C2H4] x [10 reactions]
    stoich_gas = np.zeros((n_gas, n_rxns))
    stoich_gas[0, 0] = -1   # rxn1: CO consumed
    stoich_gas[1, 1] = -1   # rxn2: H2 consumed
    stoich_gas[2, 6] = +1   # rxn7: CH4 produced
    stoich_gas[3, 8] = +1   # rxn9: H2O produced
    stoich_gas[4, 9] = +1   # rxn10: C2H4 produced

    # Surface: [CO*, H*, HCO*, CH2O*, CH2*, CH3*, O*, OH*] x [10 reactions]
    stoich_surface = np.zeros((n_surf, n_rxns))
    # rxn1: CO + * → CO*
    stoich_surface[0, 0] = +1
    # rxn2: H2 + 2* → 2H*
    stoich_surface[1, 1] = +2
    # rxn3: CO* + H* → HCO* + *
    stoich_surface[0, 2] = -1; stoich_surface[1, 2] = -1; stoich_surface[2, 2] = +1
    # rxn4: HCO* + H* → CH2O* + *
    stoich_surface[2, 3] = -1; stoich_surface[1, 3] = -1; stoich_surface[3, 3] = +1
    # rxn5: CH2O* → CH2* + O*
    stoich_surface[3, 4] = -1; stoich_surface[4, 4] = +1; stoich_surface[6, 4] = +1
    # rxn6: CH2* + H* → CH3* + *
    stoich_surface[4, 5] = -1; stoich_surface[1, 5] = -1; stoich_surface[5, 5] = +1
    # rxn7: CH3* + H* → CH4 + 2*
    stoich_surface[5, 6] = -1; stoich_surface[1, 6] = -1
    # rxn8: O* + H* → OH* + *
    stoich_surface[6, 7] = -1; stoich_surface[1, 7] = -1; stoich_surface[7, 7] = +1
    # rxn9: OH* + H* → H2O + 2*
    stoich_surface[7, 8] = -1; stoich_surface[1, 8] = -1
    # rxn10: 2CH2* → C2H4 + 2*
    stoich_surface[4, 9] = -2

    return stoich_gas, stoich_surface


def ft_rate_expressions(rate_constants_list, theta, pressures, T):
    """
    Rate expressions for FT synthesis elementary steps.

    Parameters
    ----------
    rate_constants_list : list of RateConstant
    theta : np.ndarray
        Surface coverages [CO*, H*, HCO*, CH2O*, CH2*, CH3*, O*, OH*].
    pressures : dict
        Gas-phase partial pressures [bar].
    T : float
        Temperature [K].

    Returns
    -------
    np.ndarray
        Net rates for each elementary step [1/s].
    """
    n_rxns = 10
    rates = np.zeros(n_rxns)

    theta_v = max(1.0 - np.sum(theta), 1e-15)  # Vacant sites
    P_CO = pressures.get('CO', 0)
    P_H2 = pressures.get('H2', 0)

    rc = rate_constants_list

    # rxn1: CO + * → CO* (non-activated adsorption)
    rates[0] = rc[0].k_forward * P_CO * theta_v - rc[0].k_reverse * theta[0]

    # rxn2: H2 + 2* → 2H* (dissociative)
    rates[1] = rc[1].k_forward * P_H2 * theta_v**2 - rc[1].k_reverse * theta[1]**2

    # rxn3: CO* + H* → HCO* + *
    rates[2] = rc[2].k_forward * theta[0] * theta[1] - rc[2].k_reverse * theta[2] * theta_v

    # rxn4: HCO* + H* → CH2O* + *
    rates[3] = rc[3].k_forward * theta[2] * theta[1] - rc[3].k_reverse * theta[3] * theta_v

    # rxn5: CH2O* → CH2* + O*
    rates[4] = rc[4].k_forward * theta[3] - rc[4].k_reverse * theta[4] * theta[6]

    # rxn6: CH2* + H* → CH3* + *
    rates[5] = rc[5].k_forward * theta[4] * theta[1] - rc[5].k_reverse * theta[5] * theta_v

    # rxn7: CH3* + H* → CH4 + 2*
    rates[6] = rc[6].k_forward * theta[5] * theta[1] - rc[6].k_reverse * 1e-5 * theta_v**2

    # rxn8: O* + H* → OH* + *
    rates[7] = rc[7].k_forward * theta[6] * theta[1] - rc[7].k_reverse * theta[7] * theta_v

    # rxn9: OH* + H* → H2O + 2*
    rates[8] = rc[8].k_forward * theta[7] * theta[1] - rc[8].k_reverse * 1e-5 * theta_v**2

    # rxn10: 2CH2* → C2H4 + 2*
    rates[9] = rc[9].k_forward * theta[4]**2 - rc[9].k_reverse * 1e-5 * theta_v**2

    return rates


def run_ft_case_study(T: float = 500.0, P_total: float = 20.0,
                      H2_CO_ratio: float = 2.0,
                      tunneling: str = "wigner",
                      reactor_type: str = "PFR",
                      catalyst_mass: float = 0.1) -> dict:
    """
    Run complete FT synthesis case study.

    Parameters
    ----------
    T : float
        Temperature [K].
    P_total : float
        Total pressure [bar].
    H2_CO_ratio : float
        H2/CO feed ratio.
    tunneling : str
        Tunneling correction method.
    reactor_type : str
        "PFR" or "CSTR".
    catalyst_mass : float
        Catalyst mass [kg].

    Returns
    -------
    dict
        Complete results dictionary.
    """
    results = {}

    # --- Step 1: Calculate rate constants from DFT ---
    rate_constants = []
    for ts in FT_TRANSITION_STATES:
        rc = calculate_tst_rate(ts, T, tunneling)
        rate_constants.append(rc)

    results['rate_constants'] = {
        rc.label: {
            'k_forward': rc.k_forward,
            'k_reverse': rc.k_reverse,
            'K_eq': rc.K_eq,
            'E_act_forward': rc.E_act_forward,
            'tunneling_correction': rc.tunneling_correction
        } for rc in rate_constants
    }

    # --- Step 2: Adsorption analysis ---
    x_CO = 1.0 / (1.0 + H2_CO_ratio)
    x_H2 = H2_CO_ratio / (1.0 + H2_CO_ratio)
    P_CO = x_CO * P_total
    P_H2 = x_H2 * P_total

    cov_langmuir = competitive_langmuir(
        {'CO': P_CO, 'H': P_H2},
        T, FT_ADSORPTION_PARAMS
    )

    cov_temkin_CO = temkin_isotherm(P_CO, T, FT_ADSORPTION_PARAMS['CO'],
                                    alpha=0.5, delta_E=0.3)

    results['adsorption'] = {
        'langmuir_coverages': cov_langmuir,
        'temkin_CO_coverage': cov_temkin_CO,
        'P_CO': P_CO,
        'P_H2': P_H2
    }

    # --- Step 3: Self-consistent coverages with lateral interactions ---
    cov_lateral = solve_coverage_self_consistent(
        pressures={'CO': P_CO, 'H': P_H2, 'O': 0.01},
        T=T,
        adsorption_energies={'CO': -1.30, 'H': -0.50, 'O': -2.00},
        entropy_contributions={'CO': -1.5e-3, 'H': -0.8e-3, 'O': -1.2e-3},
        interactions=FT_LATERAL_INTERACTIONS
    )
    results['lateral_interaction_coverages'] = cov_lateral

    # --- Step 4: Steady-state microkinetic solution ---
    stoich_gas, stoich_surface = build_stoichiometric_matrices()
    gas_pressures = {'CO': P_CO, 'H2': P_H2, 'CH4': 0.001, 'H2O': 0.01, 'C2H4': 0.0001}

    theta_ss = solve_steady_state_coverages(
        rate_constants, gas_pressures, T,
        FT_SURFACE_SPECIES, stoich_surface,
        ft_rate_expressions
    )
    results['steady_state_coverages'] = {sp: theta_ss[i]
                                          for i, sp in enumerate(FT_SURFACE_SPECIES)}

    # --- Step 5: RDS identification ---
    def overall_rate_func(rc_list, cov, conds):
        pressures_local = {'CO': conds.get('P_CO', P_CO),
                           'H2': conds.get('P_H2', P_H2),
                           'CH4': 0.001, 'H2O': 0.01, 'C2H4': 0.0001}
        rates = ft_rate_expressions(rc_list, theta_ss, pressures_local, T)
        return rates[6]  # CH4 formation rate as overall rate

    rds_result = degree_of_rate_control(
        rate_constants, overall_rate_func, theta_ss,
        {'T': T, 'P_CO': P_CO, 'P_H2': P_H2}
    )
    results['rds_analysis'] = {
        'step_labels': rds_result.step_labels,
        'X_RC': rds_result.X_RC.tolist(),
        'rds': rds_result.rds_label,
        'rds_index': rds_result.rds_index
    }

    # Energy span analysis
    intermediate_energies = [0.0, -1.30, -0.50, -0.60, -0.70, -1.10, -1.40, -2.00, -1.85, -2.35]
    esp = energy_span_analysis(rate_constants, intermediate_energies)
    results['energy_span'] = esp

    # --- Step 6: Reactor simulation ---
    feed = {'CO': x_CO, 'H2': x_H2, 'CH4': 0.0, 'H2O': 0.0, 'C2H4': 0.0}
    conditions = ReactorConditions(
        T=T, P_total=P_total, feed_composition=feed,
        F_total=0.01, catalyst_mass=catalyst_mass,
        site_density=1.5e-5
    )

    if reactor_type == "PFR":
        reactor = PFRReactor(conditions, FT_SURFACE_SPECIES, FT_GAS_SPECIES,
                             stoich_gas, stoich_surface, ft_rate_expressions)
        reactor_result = reactor.solve(rate_constants)
    else:
        reactor = CSTRReactor(conditions, FT_SURFACE_SPECIES, FT_GAS_SPECIES,
                              stoich_gas, stoich_surface, ft_rate_expressions)
        reactor_result = reactor.solve(rate_constants)

    results['reactor'] = {
        'type': reactor_result.reactor_type,
        'conversion': reactor_result.conversion,
        'selectivities': reactor_result.selectivities,
        'TOF': reactor_result.turnover_frequency,
        'STY': reactor_result.space_time_yield,
    }
    results['reactor_result_obj'] = reactor_result

    # --- Step 7: Arrhenius analysis ---
    arrh = arrhenius_parameters(FT_TRANSITION_STATES[4], T_range=(400, 600),
                                tunneling=tunneling)
    results['arrhenius'] = {
        'A_forward': arrh['A_forward'],
        'Ea_forward': arrh['Ea_forward'],
        'temperatures': arrh['temperatures'].tolist(),
        'k_forward': arrh['k_forward'].tolist()
    }

    return results
