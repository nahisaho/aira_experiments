"""
Short-Chain Fatty Acid (SCFA) Production Flux Model
=====================================================
Predicts SCFA (acetate, propionate, butyrate) production rates
based on microbial community composition and substrate availability.

Uses stoichiometric coefficients from known metabolic pathways:
  - Acetate: glycolysis → pyruvate → acetyl-CoA → acetate
  - Propionate: succinate pathway / acrylate pathway
  - Butyrate: acetyl-CoA condensation → butyryl-CoA → butyrate
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple


SCFA_NAMES = ['Acetate', 'Propionate', 'Butyrate']

# Species SCFA production stoichiometries (mmol/g_substrate/h per unit biomass)
# Rows: species, Columns: [acetate, propionate, butyrate]
SCFA_STOICHIOMETRY = np.array([
    [0.8, 0.3, 0.1],   # B. thetaiotaomicron - acetate/propionate producer
    [0.3, 0.1, 0.9],   # F. prausnitzii - major butyrate producer
    [0.2, 0.1, 0.8],   # R. intestinalis - butyrate producer
    [0.7, 0.1, 0.2],   # B. longum - acetate/lactate producer
    [0.5, 0.4, 0.1],   # A. muciniphila - acetate/propionate
    [0.4, 0.1, 0.0],   # E. coli - mixed acid fermentation
    [0.3, 0.0, 0.1],   # L. rhamnosus - lactate (converted to SCFA by others)
    [0.2, 0.1, 0.3],   # C. difficile - mixed
    [0.6, 0.5, 0.1],   # P. copri - acetate/propionate (succinate pathway)
    [0.4, 0.1, 0.2],   # R. bromii - acetate, feeds butyrate producers
])


@dataclass
class SCFAParameters:
    """Parameters for SCFA flux model."""
    # SCFA absorption rates from colon (h^-1)
    k_absorption_acetate: float = 0.15
    k_absorption_propionate: float = 0.12
    k_absorption_butyrate: float = 0.08  # butyrate primarily used by colonocytes

    # Colonocyte butyrate utilization fraction
    colonocyte_butyrate_usage: float = 0.70  # 70% of butyrate used by colonocytes

    # pH inhibition of SCFA production
    pH_opt_scfa: float = 6.0
    pH_width_scfa: float = 1.0

    # Substrate conversion efficiency
    fiber_to_scfa_yield: float = 0.45    # g SCFA per g fiber fermented
    starch_to_scfa_yield: float = 0.55   # g SCFA per g resistant starch
    protein_to_scfa_yield: float = 0.20  # g SCFA per g protein (also produces BCFA)

    # Branched-chain fatty acid production from protein
    protein_to_bcfa_yield: float = 0.10  # isobutyrate, isovalerate

    # Cross-feeding lactate → butyrate conversion rate
    lactate_to_butyrate_rate: float = 0.3  # fraction of lactate converted


def compute_scfa_production_rates(
    microbial_abundances: np.ndarray,
    substrate_concentrations: np.ndarray,
    params: SCFAParameters = None,
    pH: float = 5.8
) -> Dict[str, np.ndarray]:
    """
    Compute instantaneous SCFA production rates.

    Args:
        microbial_abundances: shape (n_species,) or (n_species, n_timepoints)
        substrate_concentrations: shape (5,) or (5, n_timepoints)
            [fiber, starch, protein, simple_sugars, mucin]
        params: SCFA model parameters
        pH: colonic pH

    Returns:
        Dictionary with production rates for each SCFA
    """
    if params is None:
        params = SCFAParameters()

    # pH modifier
    pH_mod = np.exp(-0.5 * ((pH - params.pH_opt_scfa) / params.pH_width_scfa) ** 2)

    # Handle both 1D and 2D inputs
    if microbial_abundances.ndim == 1:
        microbial_abundances = microbial_abundances[:, np.newaxis]
        substrate_concentrations = substrate_concentrations[:, np.newaxis]
        squeeze = True
    else:
        squeeze = False

    n_species, n_time = microbial_abundances.shape

    # Effective substrate utilization per species
    # Total fermentable substrate
    fermentable = (
        substrate_concentrations[0] * params.fiber_to_scfa_yield +
        substrate_concentrations[1] * params.starch_to_scfa_yield +
        substrate_concentrations[2] * params.protein_to_scfa_yield +
        substrate_concentrations[3] * 0.3 +
        substrate_concentrations[4] * 0.2
    )

    # SCFA production = stoichiometry × abundance × substrate × pH_modifier
    # Shape: (3, n_time)
    scfa_rates = np.zeros((3, n_time))
    for i in range(n_species):
        for j in range(3):  # acetate, propionate, butyrate
            scfa_rates[j] += (
                SCFA_STOICHIOMETRY[i, j] *
                microbial_abundances[i] *
                fermentable *
                pH_mod * 0.001  # scaling factor
            )

    # Cross-feeding: lactate from Bifidobacterium/Lactobacillus → butyrate
    lactate_production = (
        microbial_abundances[3] * 0.5 +  # B. longum
        microbial_abundances[6] * 0.8     # L. rhamnosus
    ) * fermentable * 0.001
    butyrate_from_lactate = lactate_production * params.lactate_to_butyrate_rate
    scfa_rates[2] += butyrate_from_lactate

    # Branched-chain fatty acids from protein fermentation
    bcfa_rate = (
        substrate_concentrations[2] *
        params.protein_to_bcfa_yield *
        microbial_abundances.sum(axis=0) * 0.0001
    )

    if squeeze:
        scfa_rates = scfa_rates.squeeze()
        bcfa_rate = bcfa_rate.squeeze()

    return {
        'acetate_rate': scfa_rates[0],
        'propionate_rate': scfa_rates[1],
        'butyrate_rate': scfa_rates[2],
        'bcfa_rate': bcfa_rate,
        'total_scfa_rate': scfa_rates.sum(axis=0),
        'scfa_ratios': {
            'acetate_fraction': scfa_rates[0] / (scfa_rates.sum(axis=0) + 1e-12),
            'propionate_fraction': scfa_rates[1] / (scfa_rates.sum(axis=0) + 1e-12),
            'butyrate_fraction': scfa_rates[2] / (scfa_rates.sum(axis=0) + 1e-12),
        },
        'pH_modifier': pH_mod,
    }


def compute_scfa_accumulation(
    scfa_rates: Dict,
    params: SCFAParameters = None,
    dt: float = 1.0
) -> Dict[str, float]:
    """Compute SCFA concentrations considering absorption."""
    if params is None:
        params = SCFAParameters()

    # Simple steady-state approximation: production / absorption
    acetate_conc = scfa_rates['acetate_rate'] / params.k_absorption_acetate
    propionate_conc = scfa_rates['propionate_rate'] / params.k_absorption_propionate
    butyrate_total = scfa_rates['butyrate_rate'] / params.k_absorption_butyrate
    butyrate_colonocyte = butyrate_total * params.colonocyte_butyrate_usage
    butyrate_systemic = butyrate_total * (1 - params.colonocyte_butyrate_usage)

    return {
        'acetate_mM': acetate_conc,
        'propionate_mM': propionate_conc,
        'butyrate_total_mM': butyrate_total,
        'butyrate_colonocyte_mM': butyrate_colonocyte,
        'butyrate_systemic_mM': butyrate_systemic,
        'total_scfa_mM': acetate_conc + propionate_conc + butyrate_total,
        'acetate_propionate_ratio': acetate_conc / (propionate_conc + 1e-12),
    }


def compute_scfa_timecourse(glv_results: dict, params: SCFAParameters = None) -> dict:
    """Compute SCFA production over entire gLV simulation timecourse."""
    if params is None:
        params = SCFAParameters()

    abundances = glv_results['abundances']
    substrates = glv_results['substrates']
    time = glv_results['time']

    scfa_data = compute_scfa_production_rates(abundances, substrates, params)
    accumulation = compute_scfa_accumulation(scfa_data, params)

    return {
        'time': time,
        'production_rates': scfa_data,
        'concentrations': accumulation,
        'params': params,
    }


if __name__ == "__main__":
    # Test with sample data
    abundances = np.array([100, 80, 60, 70, 30, 20, 40, 5, 50, 45], dtype=float)
    substrates = np.array([5.0, 3.0, 2.0, 1.0, 2.0])

    rates = compute_scfa_production_rates(abundances, substrates)
    accum = compute_scfa_accumulation(rates)

    print("SCFA Production Rates (mmol/h):")
    print(f"  Acetate:    {rates['acetate_rate']:.4f}")
    print(f"  Propionate: {rates['propionate_rate']:.4f}")
    print(f"  Butyrate:   {rates['butyrate_rate']:.4f}")
    print(f"\nSCFA Concentrations (mM):")
    print(f"  Acetate:    {accum['acetate_mM']:.2f}")
    print(f"  Propionate: {accum['propionate_mM']:.2f}")
    print(f"  Butyrate:   {accum['butyrate_total_mM']:.2f}")
    print(f"  Total SCFA: {accum['total_scfa_mM']:.2f}")
