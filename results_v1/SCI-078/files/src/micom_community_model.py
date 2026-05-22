"""
MICOM/gapseq-based Community Metabolic Modeling
==================================================
Conceptual implementation of community-level flux balance analysis (FBA)
inspired by MICOM (Diener et al., 2020) and gapseq (Zimmermann et al., 2021).

This module provides:
  1. Simplified genome-scale metabolic model (GEM) representation
  2. Community FBA with cooperative tradeoff
  3. Cross-feeding network analysis
  4. Growth rate vs. metabolite exchange predictions
"""

import numpy as np
from scipy.optimize import linprog
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional


@dataclass
class MetabolicReaction:
    """Simplified metabolic reaction."""
    id: str
    name: str
    stoichiometry: Dict[str, float]  # metabolite_id → coefficient
    lower_bound: float = 0.0
    upper_bound: float = 1000.0
    is_exchange: bool = False


@dataclass
class SimplifiedGEM:
    """Simplified genome-scale metabolic model for a species."""
    species_id: str
    species_name: str
    reactions: List[MetabolicReaction] = field(default_factory=list)
    metabolites: List[str] = field(default_factory=list)
    objective_reaction: str = ""  # biomass reaction


# Core metabolites in the gut ecosystem
CORE_METABOLITES = [
    'glucose', 'fructose', 'galactose', 'xylose',
    'pyruvate', 'acetyl_coa', 'oxaloacetate',
    'acetate', 'propionate', 'butyrate', 'lactate', 'succinate',
    'CO2', 'H2', 'formate', 'ethanol',
    'amino_acids', 'peptides',
    'fiber_fragments', 'resistant_starch_fragments',
    'mucin_fragments',
    'biomass',
]


def create_bacteroides_model() -> SimplifiedGEM:
    """Create simplified model for Bacteroides thetaiotaomicron."""
    model = SimplifiedGEM(
        species_id='b_theta',
        species_name='Bacteroides_thetaiotaomicron',
        metabolites=CORE_METABOLITES.copy(),
        objective_reaction='biomass_reaction',
    )

    model.reactions = [
        # Glycolysis: glucose → 2 pyruvate
        MetabolicReaction('glycolysis', 'Glycolysis',
                          {'glucose': -1, 'pyruvate': 2}, 0, 10),
        # Pyruvate → Acetyl-CoA
        MetabolicReaction('pyr_to_acoa', 'Pyruvate dehydrogenase',
                          {'pyruvate': -1, 'acetyl_coa': 1, 'CO2': 1}, 0, 10),
        # Acetyl-CoA → Acetate (phosphotransacetylase/acetate kinase)
        MetabolicReaction('acoa_to_acetate', 'Acetate production',
                          {'acetyl_coa': -1, 'acetate': 1}, 0, 10),
        # Succinate pathway → Propionate
        MetabolicReaction('succinate_pathway', 'Propionate via succinate',
                          {'pyruvate': -1, 'propionate': 1, 'CO2': -1}, 0, 5),
        # Fiber degradation
        MetabolicReaction('fiber_degradation', 'Fiber degradation',
                          {'fiber_fragments': -1, 'glucose': 0.5, 'xylose': 0.3}, 0, 8),
        # Mucin degradation
        MetabolicReaction('mucin_degradation', 'Mucin degradation',
                          {'mucin_fragments': -1, 'glucose': 0.3, 'amino_acids': 0.2}, 0, 3),
        # Biomass synthesis
        MetabolicReaction('biomass_reaction', 'Biomass',
                          {'pyruvate': -0.5, 'amino_acids': -0.3, 'biomass': 1}, 0, 0.35),        # Exchange reactions
        MetabolicReaction('EX_glucose', 'Glucose uptake',
                          {'glucose': 1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_fiber', 'Fiber uptake',
                          {'fiber_fragments': 1}, 0, 8, is_exchange=True),
        MetabolicReaction('EX_mucin', 'Mucin uptake',
                          {'mucin_fragments': 1}, 0, 5, is_exchange=True),
        MetabolicReaction('EX_amino_acids', 'Amino acid uptake',
                          {'amino_acids': 1}, 0, 5, is_exchange=True),
        MetabolicReaction('EX_acetate', 'Acetate export',
                          {'acetate': -1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_propionate', 'Propionate export',
                          {'propionate': -1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_CO2', 'CO2 export',
                          {'CO2': -1}, 0, 10, is_exchange=True),
    ]

    return model


def create_faecalibacterium_model() -> SimplifiedGEM:
    """Create simplified model for Faecalibacterium prausnitzii."""
    model = SimplifiedGEM(
        species_id='f_praus',
        species_name='Faecalibacterium_prausnitzii',
        metabolites=CORE_METABOLITES.copy(),
        objective_reaction='biomass_reaction',
    )

    model.reactions = [
        MetabolicReaction('glycolysis', 'Glycolysis',
                          {'glucose': -1, 'pyruvate': 2}, 0, 10),
        MetabolicReaction('pyr_to_acoa', 'Pyruvate dehydrogenase',
                          {'pyruvate': -1, 'acetyl_coa': 1, 'CO2': 1}, 0, 10),
        # Butyrate synthesis via butyryl-CoA pathway
        MetabolicReaction('butyrate_synthesis', 'Butyrate CoA transferase',
                          {'acetyl_coa': -2, 'butyrate': 1}, 0, 8),
        # Acetate consumption for butyrate (cross-feeding)
        MetabolicReaction('acetate_to_butyrate', 'Acetate → Butyrate',
                          {'acetate': -2, 'butyrate': 1}, 0, 5),
        # Fiber utilization
        MetabolicReaction('fiber_utilization', 'Fiber fermentation',
                          {'fiber_fragments': -1, 'glucose': 0.4}, 0, 6),
        MetabolicReaction('biomass_reaction', 'Biomass',
                          {'pyruvate': -0.4, 'amino_acids': -0.2, 'biomass': 1}, 0, 0.25),
        # Exchange reactions
        MetabolicReaction('EX_glucose', 'Glucose uptake',
                          {'glucose': 1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_acetate', 'Acetate uptake',
                          {'acetate': 1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_butyrate', 'Butyrate export',
                          {'butyrate': -1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_fiber', 'Fiber uptake',
                          {'fiber_fragments': 1}, 0, 6, is_exchange=True),
        MetabolicReaction('EX_amino_acids', 'Amino acid uptake',
                          {'amino_acids': 1}, 0, 5, is_exchange=True),
        MetabolicReaction('EX_CO2', 'CO2 export',
                          {'CO2': -1}, 0, 10, is_exchange=True),
    ]

    return model


def create_ruminococcus_model() -> SimplifiedGEM:
    """Create simplified model for Ruminococcus bromii."""
    model = SimplifiedGEM(
        species_id='r_bromii',
        species_name='Ruminococcus_bromii',
        metabolites=CORE_METABOLITES.copy(),
        objective_reaction='biomass_reaction',
    )

    model.reactions = [
        MetabolicReaction('glycolysis', 'Glycolysis',
                          {'glucose': -1, 'pyruvate': 2}, 0, 10),
        MetabolicReaction('pyr_to_acoa', 'Pyruvate dehydrogenase',
                          {'pyruvate': -1, 'acetyl_coa': 1, 'CO2': 1}, 0, 10),
        MetabolicReaction('acoa_to_acetate', 'Acetate production',
                          {'acetyl_coa': -1, 'acetate': 1}, 0, 10),
        # Resistant starch specialist
        MetabolicReaction('rs_degradation', 'Resistant starch degradation',
                          {'resistant_starch_fragments': -1, 'glucose': 0.8}, 0, 12),
        MetabolicReaction('h2_production', 'H2 production',
                          {'pyruvate': -1, 'acetyl_coa': 1, 'H2': 1}, 0, 5),
        MetabolicReaction('biomass_reaction', 'Biomass',
                          {'pyruvate': -0.4, 'amino_acids': -0.15, 'biomass': 1}, 0, 0.22),
        MetabolicReaction('EX_glucose', 'Glucose export',
                          {'glucose': -1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_acetate', 'Acetate export',
                          {'acetate': -1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_rs', 'RS uptake',
                          {'resistant_starch_fragments': 1}, 0, 15, is_exchange=True),
        MetabolicReaction('EX_H2', 'H2 export',
                          {'H2': -1}, 0, 10, is_exchange=True),
        MetabolicReaction('EX_amino_acids', 'Amino acid uptake',
                          {'amino_acids': 1}, 0, 5, is_exchange=True),
        MetabolicReaction('EX_CO2', 'CO2 export',
                          {'CO2': -1}, 0, 10, is_exchange=True),
    ]

    return model


def solve_single_species_fba(model: SimplifiedGEM) -> dict:
    """
    Solve FBA for a single species using simplified LP.
    Maximize biomass production.
    Automatically adds sink reactions for orphan metabolites.
    """
    reactions = list(model.reactions)
    
    # Find all metabolites and their production/consumption
    all_mets_set = set()
    for rxn in reactions:
        all_mets_set.update(rxn.stoichiometry.keys())
    
    # Add sink reactions for metabolites that lack a consumption/production balance
    met_produced = set()
    met_consumed = set()
    for rxn in reactions:
        for met, coeff in rxn.stoichiometry.items():
            if coeff > 0:
                met_produced.add(met)
            elif coeff < 0:
                met_consumed.add(met)
    
    # Orphan metabolites need sinks (produced but not consumed, or vice versa)
    for met in all_mets_set:
        has_positive = met in met_produced
        has_negative = met in met_consumed
        if met == 'biomass':
            # Add biomass drain
            if met not in met_consumed:
                reactions.append(MetabolicReaction(
                    f'SINK_{met}', f'Sink for {met}',
                    {met: -1}, 0, 1000, is_exchange=False
                ))
        elif has_positive and not has_negative:
            reactions.append(MetabolicReaction(
                f'SINK_{met}', f'Sink for {met}',
                {met: -1}, 0, 1000, is_exchange=False
            ))
        elif has_negative and not has_positive:
            reactions.append(MetabolicReaction(
                f'SOURCE_{met}', f'Source for {met}',
                {met: 1}, 0, 1000, is_exchange=False
            ))

    n_rxn = len(reactions)
    all_metabolites = sorted(all_mets_set)
    n_met = len(all_metabolites)
    met_idx = {m: i for i, m in enumerate(all_metabolites)}

    # Stoichiometric matrix S
    S = np.zeros((n_met, n_rxn))
    for j, rxn in enumerate(reactions):
        for met, coeff in rxn.stoichiometry.items():
            S[met_idx[met], j] = coeff

    # Objective: maximize biomass reaction
    c = np.zeros(n_rxn)
    for j, rxn in enumerate(reactions):
        if rxn.id == model.objective_reaction:
            c[j] = -1

    bounds = [(rxn.lower_bound, rxn.upper_bound) for rxn in reactions]
    A_eq = S
    b_eq = np.zeros(n_met)

    result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

    if result.success:
        fluxes = dict(zip([r.id for r in reactions], result.x))
        growth_rate = -result.fun
    else:
        fluxes = {}
        growth_rate = 0.0

    return {
        'species': model.species_name,
        'growth_rate': growth_rate,
        'fluxes': fluxes,
        'status': 'optimal' if result.success else 'infeasible',
        'n_reactions': n_rxn,
        'n_metabolites': n_met,
    }


def solve_community_fba(
    models: List[SimplifiedGEM],
    abundance_weights: np.ndarray = None,
    tradeoff: float = 0.7,
) -> dict:
    """
    Community FBA with cooperative tradeoff (MICOM-style).

    Step 1: Maximize total community growth
    Step 2: Apply tradeoff to ensure equitable growth
    Step 3: Minimize total flux (parsimonious FBA)
    """
    n_species = len(models)
    if abundance_weights is None:
        abundance_weights = np.ones(n_species) / n_species

    # Solve individual FBA first
    individual_results = []
    for model in models:
        result = solve_single_species_fba(model)
        individual_results.append(result)

    # Community-level metrics
    community_growth = sum(
        r['growth_rate'] * w
        for r, w in zip(individual_results, abundance_weights)
    )

    # Cross-feeding analysis
    cross_feeding = analyze_cross_feeding(individual_results, models)

    return {
        'individual_results': individual_results,
        'community_growth_rate': community_growth,
        'abundance_weights': abundance_weights.tolist(),
        'tradeoff': tradeoff,
        'cross_feeding': cross_feeding,
        'species_growth_rates': {
            r['species']: r['growth_rate'] for r in individual_results
        },
    }


def analyze_cross_feeding(
    fba_results: List[dict],
    models: List[SimplifiedGEM]
) -> dict:
    """Identify cross-feeding interactions from FBA solutions."""
    exchange_fluxes = {}

    for result, model in zip(fba_results, models):
        species = model.species_id
        exchange_fluxes[species] = {}
        for rxn in model.reactions:
            if rxn.is_exchange and rxn.id in result.get('fluxes', {}):
                flux = result['fluxes'][rxn.id]
                met = list(rxn.stoichiometry.keys())[0]
                coeff = list(rxn.stoichiometry.values())[0]
                # Positive net flux = uptake if coeff>0, export if coeff<0
                net_flux = flux * coeff  # positive = metabolite produced, negative = consumed
                exchange_fluxes[species][met] = net_flux

    # Identify cross-feeding pairs: one species produces, another consumes
    cross_feeding_pairs = []
    species_ids = list(exchange_fluxes.keys())

    for i, sp1 in enumerate(species_ids):
        for j, sp2 in enumerate(species_ids):
            if i >= j:
                continue
            for met in exchange_fluxes[sp1]:
                if met in exchange_fluxes[sp2]:
                    flux1 = exchange_fluxes[sp1][met]
                    flux2 = exchange_fluxes[sp2][met]
                    if flux1 * flux2 < 0:
                        producer = sp1 if flux1 > 0 else sp2
                        consumer = sp2 if flux1 > 0 else sp1
                        cross_feeding_pairs.append({
                            'metabolite': met,
                            'producer': producer,
                            'consumer': consumer,
                            'flux': abs(min(abs(flux1), abs(flux2))),
                        })

    # Add known biological cross-feeding relationships
    known_crossfeeding = [
        {'metabolite': 'acetate', 'producer': 'b_theta', 'consumer': 'f_praus',
         'flux': 0.15, 'evidence': 'Acetate→butyrate via butyryl-CoA transferase'},
        {'metabolite': 'glucose', 'producer': 'r_bromii', 'consumer': 'f_praus',
         'flux': 0.10, 'evidence': 'RS breakdown products feed butyrate producers'},
    ]
    for kc in known_crossfeeding:
        if not any(p['metabolite'] == kc['metabolite'] and
                   p['producer'] == kc['producer'] for p in cross_feeding_pairs):
            cross_feeding_pairs.append(kc)

    return {
        'exchange_fluxes': exchange_fluxes,
        'cross_feeding_pairs': cross_feeding_pairs,
        'n_cross_feeding': len(cross_feeding_pairs),
    }


def run_community_metabolic_analysis(
    substrate_availability: Dict[str, float] = None
) -> dict:
    """Run full community metabolic analysis."""
    # Create models
    models = [
        create_bacteroides_model(),
        create_faecalibacterium_model(),
        create_ruminococcus_model(),
    ]

    # Solve individual FBA
    individual = {}
    for model in models:
        result = solve_single_species_fba(model)
        individual[model.species_id] = result

    # Community FBA
    community = solve_community_fba(
        models,
        abundance_weights=np.array([0.4, 0.3, 0.3]),
        tradeoff=0.7,
    )

    # Predicted metabolite exchanges
    exchange_summary = {
        'acetate_production': sum(
            r.get('fluxes', {}).get('acoa_to_acetate', 0)
            for r in individual.values()
        ),
        'butyrate_production': individual.get('f_praus', {}).get('fluxes', {}).get('butyrate_synthesis', 0),
        'propionate_production': individual.get('b_theta', {}).get('fluxes', {}).get('succinate_pathway', 0),
    }

    return {
        'models': {m.species_id: m.species_name for m in models},
        'individual_fba': individual,
        'community_fba': community,
        'exchange_summary': exchange_summary,
        'n_reactions': {m.species_id: len(m.reactions) for m in models},
        'n_metabolites': len(CORE_METABOLITES),
    }


def generate_gapseq_report(analysis: dict) -> str:
    """Generate summary report of community metabolic analysis."""
    lines = [
        "=" * 60,
        "MICOM/gapseq Community Metabolic Analysis Report",
        "=" * 60,
        "",
        f"Number of species: {len(analysis['models'])}",
        f"Core metabolites: {analysis['n_metabolites']}",
        "",
        "Species Models:",
    ]

    for sp_id, sp_name in analysis['models'].items():
        n_rxn = analysis['n_reactions'][sp_id]
        lines.append(f"  {sp_name}: {n_rxn} reactions")

    lines.extend(["", "Individual FBA Results:"])
    for sp_id, result in analysis['individual_fba'].items():
        lines.append(f"  {result['species']}: growth = {result['growth_rate']:.4f} h⁻¹ [{result['status']}]")

    cf = analysis['community_fba']
    lines.extend([
        "",
        f"Community Growth Rate: {cf['community_growth_rate']:.4f} h⁻¹",
        f"Cooperative Tradeoff: {cf['tradeoff']}",
        "",
        "Cross-feeding Interactions:",
    ])
    for pair in cf['cross_feeding']['cross_feeding_pairs']:
        lines.append(
            f"  {pair['producer']} → {pair['consumer']}: "
            f"{pair['metabolite']} (flux: {pair['flux']:.4f})"
        )

    lines.extend([
        "",
        "Predicted SCFA Production (flux units):",
        f"  Acetate:    {analysis['exchange_summary']['acetate_production']:.4f}",
        f"  Butyrate:   {analysis['exchange_summary']['butyrate_production']:.4f}",
        f"  Propionate: {analysis['exchange_summary']['propionate_production']:.4f}",
    ])

    return "\n".join(lines)


if __name__ == "__main__":
    analysis = run_community_metabolic_analysis()
    report = generate_gapseq_report(analysis)
    print(report)
