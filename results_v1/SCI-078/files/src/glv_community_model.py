"""
Generalized Lotka-Volterra (gLV) Community Model
==================================================
Models competitive and cooperative interactions among gut microbial species
using generalized Lotka-Volterra equations with resource-dependent growth.

dX_i/dt = X_i * (mu_i(S) + sum_j(A_ij * X_j) - delta_i)

Where:
  X_i: abundance of species i
  mu_i(S): substrate-dependent growth rate (Monod kinetics)
  A_ij: interaction matrix (competition/cooperation)
  delta_i: dilution/death rate
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


# Representative gut microbial species
SPECIES_NAMES = [
    "Bacteroides_thetaiotaomicron",
    "Faecalibacterium_prausnitzii",
    "Roseburia_intestinalis",
    "Bifidobacterium_longum",
    "Akkermansia_muciniphila",
    "Escherichia_coli",
    "Lactobacillus_rhamnosus",
    "Clostridium_difficile",
    "Prevotella_copri",
    "Ruminococcus_bromii",
]

# Short names for display
SPECIES_SHORT = [
    "B. theta.", "F. praus.", "R. intest.", "B. longum",
    "A. mucin.", "E. coli", "L. rhamn.",
    "C. diff.", "P. copri", "R. bromii"
]


@dataclass
class GLVParameters:
    """Parameters for the gLV community model."""
    n_species: int = 10

    # Maximum growth rates (h^-1)
    mu_max: np.ndarray = field(default_factory=lambda: np.array([
        0.35,  # B. thetaiotaomicron - versatile glycan degrader
        0.25,  # F. prausnitzii - butyrate producer
        0.28,  # R. intestinalis - butyrate producer
        0.30,  # B. longum - fiber fermenter
        0.15,  # A. muciniphila - mucin specialist
        0.45,  # E. coli - fast grower
        0.32,  # L. rhamnosus - lactic acid producer
        0.20,  # C. difficile - opportunistic
        0.33,  # P. copri - plant polysaccharide specialist
        0.22,  # R. bromii - resistant starch specialist
    ]))

    # Substrate affinity (Monod Ks, g/L)
    Ks: np.ndarray = field(default_factory=lambda: np.array([
        2.0, 3.0, 3.5, 2.5, 5.0, 1.0, 2.0, 4.0, 2.5, 4.0
    ]))

    # Dilution/death rates (h^-1)
    delta: np.ndarray = field(default_factory=lambda: np.array([
        0.02, 0.03, 0.03, 0.02, 0.04, 0.05, 0.03, 0.06, 0.02, 0.03
    ]))

    # Substrate utilization preferences [species x substrate_type]
    # Substrates: [fiber, starch, protein, simple_sugars, mucin]
    substrate_preference: np.ndarray = field(default_factory=lambda: np.array([
        [0.4, 0.3, 0.1, 0.1, 0.1],  # B. theta - generalist
        [0.5, 0.2, 0.1, 0.1, 0.1],  # F. praus - fiber specialist
        [0.4, 0.3, 0.1, 0.2, 0.0],  # R. intest
        [0.3, 0.2, 0.0, 0.4, 0.1],  # B. longum - oligosaccharide preference
        [0.0, 0.0, 0.0, 0.0, 1.0],  # A. mucin - mucin specialist
        [0.1, 0.1, 0.3, 0.5, 0.0],  # E. coli - simple sugars
        [0.1, 0.1, 0.1, 0.6, 0.1],  # L. rhamnosus
        [0.1, 0.1, 0.3, 0.4, 0.1],  # C. difficile
        [0.5, 0.3, 0.1, 0.1, 0.0],  # P. copri - plant polysaccharides
        [0.1, 0.7, 0.0, 0.1, 0.1],  # R. bromii - resistant starch
    ]))


def build_interaction_matrix(n_species: int = 10, seed: int = 42) -> np.ndarray:
    """
    Build species interaction matrix A_ij.
    Negative values = competition, positive = cooperation (cross-feeding).
    Diagonal = self-limitation (carrying capacity).
    """
    rng = np.random.RandomState(seed)

    A = np.zeros((n_species, n_species))

    # Self-limitation (negative diagonal)
    np.fill_diagonal(A, -np.array([
        0.001, 0.0012, 0.0011, 0.001, 0.0015,
        0.0008, 0.001, 0.0018, 0.001, 0.0013
    ]))

    # Known cross-feeding interactions (positive)
    # R. bromii → F. prausnitzii (starch breakdown products → butyrate)
    A[1, 9] = 0.0003
    # B. thetaiotaomicron → F. prausnitzii (acetate cross-feeding)
    A[1, 0] = 0.0002
    # B. longum → R. intestinalis (oligosaccharide release)
    A[2, 3] = 0.0002
    # B. theta → B. longum (glycan breakdown products)
    A[3, 0] = 0.00015

    # Competition interactions (negative)
    # E. coli vs L. rhamnosus (niche overlap)
    A[5, 6] = -0.0004
    A[6, 5] = -0.0003
    # C. difficile inhibited by many commensals
    A[7, 1] = -0.0005  # F. praus inhibits C. diff
    A[7, 6] = -0.0006  # L. rhamnosus inhibits C. diff
    A[7, 3] = -0.0004  # B. longum inhibits C. diff
    # P. copri vs B. theta (niche overlap on plant polysaccharides)
    A[8, 0] = -0.00025
    A[0, 8] = -0.0002

    # Add weak random interactions
    for i in range(n_species):
        for j in range(n_species):
            if i != j and A[i, j] == 0:
                A[i, j] = rng.normal(0, 0.00005)

    return A


def compute_growth_rate(
    species_idx: int,
    substrates: np.ndarray,
    params: GLVParameters
) -> float:
    """Compute substrate-dependent growth rate using Monod kinetics."""
    # Weighted substrate availability
    S_eff = np.dot(params.substrate_preference[species_idx], substrates)
    return params.mu_max[species_idx] * S_eff / (params.Ks[species_idx] + S_eff)


def glv_ode_system(t, y, params: GLVParameters, A: np.ndarray,
                   substrate_func=None):
    """
    gLV ODE system with resource-dependent growth.

    State: y[0:n_species] = species abundances
           y[n_species:n_species+5] = substrate concentrations
    """
    n = params.n_species
    X = np.maximum(y[:n], 0)  # species abundances
    S = np.maximum(y[n:n+5], 0)  # substrates [fiber, starch, protein, sugars, mucin]

    dXdt = np.zeros(n)
    dSdt = np.zeros(5)

    # External substrate supply (simulating continuous feeding)
    if substrate_func is not None:
        S_supply = substrate_func(t)
    else:
        S_supply = np.array([0.5, 0.3, 0.2, 0.1, 0.2])  # basal supply (g/L/h)

    # Species dynamics
    for i in range(n):
        mu_i = compute_growth_rate(i, S, params)
        interaction = np.dot(A[i], X)
        dXdt[i] = X[i] * (mu_i + interaction - params.delta[i])

        # Substrate consumption by species i
        for k in range(5):
            consumption = params.substrate_preference[i, k] * mu_i * X[i] * 0.1
            dSdt[k] -= consumption

    # Substrate dynamics: supply - consumption - dilution
    dilution_rate = 0.04  # h^-1
    dSdt += S_supply - dilution_rate * S

    dydt = np.concatenate([dXdt, dSdt])
    return dydt


def run_glv_simulation(
    initial_abundances: np.ndarray = None,
    initial_substrates: np.ndarray = None,
    params: GLVParameters = None,
    A: np.ndarray = None,
    substrate_func=None,
    t_span: tuple = (0, 720),  # 30 days in hours
    t_points: int = 1000
) -> dict:
    """Run gLV community simulation."""
    if params is None:
        params = GLVParameters()
    if A is None:
        A = build_interaction_matrix(params.n_species)

    n = params.n_species

    if initial_abundances is None:
        # Typical healthy gut composition (relative abundances scaled)
        initial_abundances = np.array([
            100, 80, 60, 70, 30, 20, 40, 5, 50, 45
        ], dtype=float)

    if initial_substrates is None:
        initial_substrates = np.array([5.0, 3.0, 2.0, 1.0, 2.0])

    y0 = np.concatenate([initial_abundances, initial_substrates])
    t_eval = np.linspace(t_span[0], t_span[1], t_points)

    sol = solve_ivp(
        glv_ode_system, t_span, y0,
        args=(params, A, substrate_func),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8, atol=1e-10,
        max_step=1.0
    )

    if not sol.success:
        raise RuntimeError(f"gLV simulation failed: {sol.message}")

    abundances = np.maximum(sol.y[:n], 0)
    substrates = np.maximum(sol.y[n:n+5], 0)

    # Compute diversity metrics
    total = abundances.sum(axis=0)
    rel_abundances = abundances / (total + 1e-12)

    # Shannon diversity
    shannon = -np.sum(
        rel_abundances * np.log(rel_abundances + 1e-12), axis=0
    )

    # Simpson diversity
    simpson = 1 - np.sum(rel_abundances ** 2, axis=0)

    return {
        'time': sol.t,
        'abundances': abundances,
        'substrates': substrates,
        'relative_abundances': rel_abundances,
        'shannon_diversity': shannon,
        'simpson_diversity': simpson,
        'species_names': SPECIES_NAMES,
        'species_short': SPECIES_SHORT,
        'interaction_matrix': A,
        'params': params,
    }


def compute_steady_state_composition(results: dict) -> dict:
    """Extract steady-state community composition."""
    # Use last 10% of simulation
    n_tail = max(1, len(results['time']) // 10)
    mean_abundance = results['abundances'][:, -n_tail:].mean(axis=1)
    total = mean_abundance.sum()
    rel_abundance = mean_abundance / (total + 1e-12)

    return {
        'species': SPECIES_NAMES,
        'mean_abundance': mean_abundance,
        'relative_abundance': rel_abundance,
        'shannon_diversity': results['shannon_diversity'][-n_tail:].mean(),
        'simpson_diversity': results['simpson_diversity'][-n_tail:].mean(),
        'dominant_species': SPECIES_NAMES[np.argmax(rel_abundance)],
    }


if __name__ == "__main__":
    results = run_glv_simulation()
    ss = compute_steady_state_composition(results)
    print("Steady-state composition:")
    for name, ra in zip(ss['species'], ss['relative_abundance']):
        print(f"  {name}: {ra:.3f}")
    print(f"Shannon diversity: {ss['shannon_diversity']:.3f}")
    print(f"Dominant species: {ss['dominant_species']}")
