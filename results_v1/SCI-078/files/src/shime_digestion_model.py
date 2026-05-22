"""
SHIME (Simulator of the Human Intestinal Microbial Ecosystem) Digestion Model
==============================================================================
Multi-compartment kinetic model simulating food component digestion and absorption
through the gastrointestinal tract.

Compartments:
  1. Stomach: Acid hydrolysis, pepsin proteolysis
  2. Small intestine (duodenum/jejunum/ileum): Enzymatic digestion, bile salts, absorption
  3. Ascending colon: Proximal fermentation
  4. Transverse colon: Main fermentation
  5. Descending colon: Distal fermentation, water reabsorption
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional


@dataclass
class SHIMEParameters:
    """Parameters for the SHIME digestion model."""
    # Compartment volumes (mL)
    V_stomach: float = 200.0
    V_small_intestine: float = 300.0
    V_ascending_colon: float = 500.0
    V_transverse_colon: float = 400.0
    V_descending_colon: float = 300.0

    # Transit rates (h^-1) - gastric emptying and intestinal transit
    k_gastric_emptying: float = 0.5       # stomach → small intestine
    k_si_transit: float = 0.3             # small intestine → ascending colon
    k_ac_transit: float = 0.08            # ascending → transverse colon
    k_tc_transit: float = 0.06            # transverse → descending colon
    k_dc_transit: float = 0.04            # descending colon → excretion

    # Enzymatic digestion rates (h^-1)
    k_protein_stomach: float = 0.15       # pepsin activity
    k_protein_si: float = 0.8            # trypsin/chymotrypsin
    k_starch_si: float = 1.2             # amylase
    k_lipid_si: float = 0.6              # lipase + bile salts
    k_fiber_colon: float = 0.05          # microbial fermentation

    # Absorption rates (h^-1) in small intestine
    k_abs_amino_acids: float = 0.9
    k_abs_glucose: float = 1.5
    k_abs_fatty_acids: float = 0.7
    k_abs_vitamins: float = 0.4

    # pH values per compartment
    pH_stomach: float = 2.0
    pH_small_intestine: float = 6.5
    pH_ascending_colon: float = 5.5
    pH_transverse_colon: float = 6.2
    pH_descending_colon: float = 6.8

    # Michaelis-Menten parameters
    Km_starch: float = 10.0    # g/L
    Km_protein: float = 8.0
    Km_lipid: float = 5.0
    Km_fiber: float = 15.0


@dataclass
class FoodComposition:
    """Nutritional composition of a food item (grams per serving)."""
    name: str = "mixed_meal"
    protein: float = 30.0
    starch: float = 50.0
    simple_sugars: float = 10.0
    dietary_fiber: float = 15.0
    soluble_fiber: float = 5.0
    insoluble_fiber: float = 10.0
    lipids: float = 20.0
    polyphenols: float = 0.5
    resistant_starch: float = 5.0


def pH_activity_modifier(pH: float, pH_opt: float, pH_width: float = 1.5) -> float:
    """Gaussian pH-dependent enzyme activity modifier."""
    return np.exp(-0.5 * ((pH - pH_opt) / pH_width) ** 2)


def michaelis_menten(S: float, Vmax: float, Km: float) -> float:
    """Michaelis-Menten kinetics."""
    return Vmax * S / (Km + S) if S > 0 else 0.0


def shime_ode_system(t, y, params: SHIMEParameters):
    """
    ODE system for SHIME multi-compartment digestion.

    State variables (25 total):
    Compartment 1 - Stomach (0-4): protein, starch, lipid, fiber, polyphenols
    Compartment 2 - Small Intestine (5-9): same
    Compartment 3 - Ascending Colon (10-14): same
    Compartment 4 - Transverse Colon (15-19): same
    Compartment 5 - Descending Colon (20-24): same
    """
    p = params
    dydt = np.zeros(25)

    # Unpack state variables
    # Stomach
    prot_st, starch_st, lipid_st, fiber_st, poly_st = y[0:5]
    # Small Intestine
    prot_si, starch_si, lipid_si, fiber_si, poly_si = y[5:10]
    # Ascending Colon
    prot_ac, starch_ac, lipid_ac, fiber_ac, poly_ac = y[10:15]
    # Transverse Colon
    prot_tc, starch_tc, lipid_tc, fiber_tc, poly_tc = y[15:20]
    # Descending Colon
    prot_dc, starch_dc, lipid_dc, fiber_dc, poly_dc = y[20:25]

    # --- Stomach ---
    pepsin_activity = pH_activity_modifier(p.pH_stomach, pH_opt=2.0, pH_width=1.0)
    r_prot_st = michaelis_menten(prot_st, p.k_protein_stomach * pepsin_activity, p.Km_protein)

    dydt[0] = -r_prot_st - p.k_gastric_emptying * prot_st
    dydt[1] = -p.k_gastric_emptying * starch_st
    dydt[2] = -p.k_gastric_emptying * lipid_st
    dydt[3] = -p.k_gastric_emptying * fiber_st
    dydt[4] = -p.k_gastric_emptying * poly_st

    # --- Small Intestine ---
    si_pH_mod = pH_activity_modifier(p.pH_small_intestine, pH_opt=7.0, pH_width=1.5)
    r_prot_si = michaelis_menten(prot_si, p.k_protein_si * si_pH_mod, p.Km_protein)
    r_starch_si = michaelis_menten(starch_si, p.k_starch_si * si_pH_mod, p.Km_starch)
    r_lipid_si = michaelis_menten(lipid_si, p.k_lipid_si * si_pH_mod, p.Km_lipid)

    # Digestion removes mass; absorption removes digested products proportionally
    dydt[5] = p.k_gastric_emptying * prot_st - r_prot_si - p.k_abs_amino_acids * prot_si * 0.5 - p.k_si_transit * prot_si
    dydt[6] = p.k_gastric_emptying * starch_st - r_starch_si - p.k_abs_glucose * starch_si * 0.5 - p.k_si_transit * starch_si
    dydt[7] = p.k_gastric_emptying * lipid_st - r_lipid_si - p.k_abs_fatty_acids * lipid_si * 0.5 - p.k_si_transit * lipid_si
    dydt[8] = p.k_gastric_emptying * fiber_st - p.k_si_transit * fiber_si
    dydt[9] = p.k_gastric_emptying * poly_st - p.k_abs_vitamins * poly_si * 0.3 - p.k_si_transit * poly_si

    # --- Ascending Colon ---
    ac_ferm = pH_activity_modifier(p.pH_ascending_colon, pH_opt=5.5, pH_width=1.0)
    r_fiber_ac = michaelis_menten(fiber_ac, p.k_fiber_colon * ac_ferm * 1.5, p.Km_fiber)
    r_starch_ac = michaelis_menten(starch_ac, p.k_fiber_colon * ac_ferm * 2.0, p.Km_starch)

    dydt[10] = p.k_si_transit * prot_si - p.k_ac_transit * prot_ac * 0.3
    dydt[11] = p.k_si_transit * starch_si - r_starch_ac - p.k_ac_transit * starch_ac
    dydt[12] = p.k_si_transit * lipid_si - p.k_ac_transit * lipid_ac
    dydt[13] = p.k_si_transit * fiber_si - r_fiber_ac - p.k_ac_transit * fiber_ac
    dydt[14] = p.k_si_transit * poly_si - p.k_ac_transit * poly_ac

    # --- Transverse Colon ---
    tc_ferm = pH_activity_modifier(p.pH_transverse_colon, pH_opt=6.0, pH_width=1.2)
    r_fiber_tc = michaelis_menten(fiber_tc, p.k_fiber_colon * tc_ferm, p.Km_fiber)

    dydt[15] = p.k_ac_transit * prot_ac - p.k_tc_transit * prot_tc * 0.2
    dydt[16] = p.k_ac_transit * starch_ac - p.k_tc_transit * starch_tc
    dydt[17] = p.k_ac_transit * lipid_ac - p.k_tc_transit * lipid_tc
    dydt[18] = p.k_ac_transit * fiber_ac - r_fiber_tc - p.k_tc_transit * fiber_tc
    dydt[19] = p.k_ac_transit * poly_ac - p.k_tc_transit * poly_tc

    # --- Descending Colon ---
    dc_ferm = pH_activity_modifier(p.pH_descending_colon, pH_opt=6.5, pH_width=1.0)
    r_fiber_dc = michaelis_menten(fiber_dc, p.k_fiber_colon * dc_ferm * 0.5, p.Km_fiber)

    dydt[20] = p.k_tc_transit * prot_tc - p.k_dc_transit * prot_dc
    dydt[21] = p.k_tc_transit * starch_tc - p.k_dc_transit * starch_dc
    dydt[22] = p.k_tc_transit * lipid_tc - p.k_dc_transit * lipid_dc
    dydt[23] = p.k_tc_transit * fiber_tc - r_fiber_dc - p.k_dc_transit * fiber_dc
    dydt[24] = p.k_tc_transit * poly_tc - p.k_dc_transit * poly_dc

    return dydt


def run_shime_simulation(
    food: FoodComposition,
    params: SHIMEParameters = None,
    t_span: tuple = (0, 72),
    t_points: int = 500
) -> dict:
    """Run SHIME digestion simulation."""
    if params is None:
        params = SHIMEParameters()

    # Initial conditions: food components in stomach, rest empty
    y0 = np.zeros(25)
    y0[0] = food.protein
    y0[1] = food.starch + food.resistant_starch
    y0[2] = food.lipids
    y0[3] = food.dietary_fiber
    y0[4] = food.polyphenols

    t_eval = np.linspace(t_span[0], t_span[1], t_points)

    sol = solve_ivp(
        shime_ode_system, t_span, y0,
        args=(params,),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8, atol=1e-10,
        max_step=0.1
    )

    if not sol.success:
        raise RuntimeError(f"SHIME simulation failed: {sol.message}")

    compartment_names = ['Stomach', 'Small_Intestine', 'Ascending_Colon',
                         'Transverse_Colon', 'Descending_Colon']
    nutrient_names = ['Protein', 'Starch', 'Lipid', 'Fiber', 'Polyphenols']

    results = {
        'time': sol.t,
        'solution': sol.y,
        'compartments': compartment_names,
        'nutrients': nutrient_names,
        'food': food,
        'params': params
    }

    # Compute substrate availability for colonic microbiota
    results['colonic_substrates'] = {
        'fiber_ascending': sol.y[13],
        'fiber_transverse': sol.y[18],
        'fiber_descending': sol.y[23],
        'starch_ascending': sol.y[11],
        'protein_ascending': sol.y[10],
        'polyphenols_ascending': sol.y[14],
    }

    return results


def compute_absorption_efficiency(results: dict) -> dict:
    """Compute nutrient absorption efficiency from SHIME results."""
    food = results['food']
    sol = results['solution']

    # Estimate absorbed fraction = initial - remaining at end
    total_initial = {
        'protein': food.protein,
        'starch': food.starch + food.resistant_starch,
        'lipid': food.lipids,
        'fiber': food.dietary_fiber,
        'polyphenols': food.polyphenols
    }

    # Sum remaining across all compartments at final time
    remaining = {}
    for i, nutrient in enumerate(['protein', 'starch', 'lipid', 'fiber', 'polyphenols']):
        total_remaining = sum(sol[j * 5 + i, -1] for j in range(5))
        remaining[nutrient] = max(0, total_remaining)

    efficiency = {}
    for nutrient in total_initial:
        if total_initial[nutrient] > 0:
            absorbed = total_initial[nutrient] - remaining[nutrient]
            efficiency[nutrient] = max(0, min(1, absorbed / total_initial[nutrient]))
        else:
            efficiency[nutrient] = 0.0

    return {
        'total_initial': total_initial,
        'remaining': remaining,
        'absorption_efficiency': efficiency
    }


if __name__ == "__main__":
    food = FoodComposition(
        name="high_fiber_meal",
        protein=25.0, starch=40.0, simple_sugars=8.0,
        dietary_fiber=20.0, soluble_fiber=8.0, insoluble_fiber=12.0,
        lipids=15.0, polyphenols=1.0, resistant_starch=8.0
    )
    results = run_shime_simulation(food)
    efficiency = compute_absorption_efficiency(results)
    print("Absorption Efficiency:")
    for k, v in efficiency['absorption_efficiency'].items():
        print(f"  {k}: {v:.1%}")
