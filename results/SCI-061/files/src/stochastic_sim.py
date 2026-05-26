"""
Stochastic Simulation Module — Gillespie SSA and Tau-Leaping
for synthetic gene circuit dynamics.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Callable, Tuple, Optional


@dataclass
class Reaction:
    name: str
    propensity_func: Callable  # f(state, params) -> rate
    stoichiometry: Dict[str, int]  # species -> change


@dataclass
class StochasticModel:
    species: List[str]
    reactions: List[Reaction]
    initial_state: Dict[str, int]
    parameters: Dict[str, float]


def gillespie_ssa(
    model: StochasticModel,
    t_end: float,
    seed: int = 42,
    max_steps: int = 10_000_000
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Exact Gillespie Stochastic Simulation Algorithm.
    Returns (times, states) arrays.
    """
    rng = np.random.RandomState(seed)
    n_species = len(model.species)
    sp_idx = {s: i for i, s in enumerate(model.species)}

    state = np.zeros(n_species, dtype=np.float64)
    for s, v in model.initial_state.items():
        if s in sp_idx:
            state[sp_idx[s]] = v

    stoich_matrix = np.zeros((len(model.reactions), n_species))
    for j, rxn in enumerate(model.reactions):
        for sp, delta in rxn.stoichiometry.items():
            if sp in sp_idx:
                stoich_matrix[j, sp_idx[sp]] = delta

    times = [0.0]
    states = [state.copy()]
    t = 0.0
    step = 0

    while t < t_end and step < max_steps:
        state_dict = {s: state[sp_idx[s]] for s in model.species}
        props = np.array([
            rxn.propensity_func(state_dict, model.parameters)
            for rxn in model.reactions
        ])
        props = np.maximum(props, 0.0)
        a0 = props.sum()

        if a0 == 0:
            break

        tau = rng.exponential(1.0 / a0)
        t += tau

        if t > t_end:
            break

        j = rng.choice(len(model.reactions), p=props / a0)
        state += stoich_matrix[j]
        state = np.maximum(state, 0)

        times.append(t)
        states.append(state.copy())
        step += 1

    return np.array(times), np.array(states)


def tau_leaping(
    model: StochasticModel,
    t_end: float,
    tau: float = 0.1,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Tau-leaping approximate stochastic simulation.
    """
    rng = np.random.RandomState(seed)
    n_species = len(model.species)
    sp_idx = {s: i for i, s in enumerate(model.species)}

    state = np.zeros(n_species, dtype=np.float64)
    for s, v in model.initial_state.items():
        if s in sp_idx:
            state[sp_idx[s]] = v

    stoich_matrix = np.zeros((len(model.reactions), n_species))
    for j, rxn in enumerate(model.reactions):
        for sp, delta in rxn.stoichiometry.items():
            if sp in sp_idx:
                stoich_matrix[j, sp_idx[sp]] = delta

    times = [0.0]
    states = [state.copy()]
    t = 0.0

    while t < t_end:
        state_dict = {s: state[sp_idx[s]] for s in model.species}
        props = np.array([
            rxn.propensity_func(state_dict, model.parameters)
            for rxn in model.reactions
        ])
        props = np.maximum(props, 0.0)

        # Number of firings per reaction in interval tau
        firings = rng.poisson(props * tau)
        state += stoich_matrix.T @ firings
        state = np.maximum(state, 0)

        t += tau
        times.append(t)
        states.append(state.copy())

    return np.array(times), np.array(states)


# ---- Pre-built circuit models ----

def build_toggle_switch_model(params: Optional[Dict] = None) -> StochasticModel:
    """Build stochastic model for the genetic toggle switch."""
    default_params = {
        "alpha1": 3.0,   # max expression rate of repressor 1
        "alpha2": 2.5,   # max expression rate of repressor 2
        "K1": 50.0,      # half-max constant 1
        "K2": 50.0,      # half-max constant 2
        "n1": 2.0,       # Hill coefficient 1
        "n2": 2.5,       # Hill coefficient 2
        "delta1": 0.05,  # degradation rate 1
        "delta2": 0.05,  # degradation rate 2
        "leak1": 0.01,   # basal leakage 1
        "leak2": 0.01,   # basal leakage 2
        "IPTG": 0.0,     # inducer concentration
        "aTc": 0.0,
    }
    if params:
        default_params.update(params)

    species = ["LacI", "TetR"]

    reactions = [
        Reaction(
            "LacI_production",
            lambda s, p: (p["leak1"] + p["alpha1"] *
                          (p["K1"]**p["n1"]) /
                          (p["K1"]**p["n1"] + (s["TetR"] * max(0.01, 1.0 - p["aTc"]))**p["n1"])),
            {"LacI": 1}
        ),
        Reaction(
            "TetR_production",
            lambda s, p: (p["leak2"] + p["alpha2"] *
                          (p["K2"]**p["n2"]) /
                          (p["K2"]**p["n2"] + (s["LacI"] * max(0.01, 1.0 - p["IPTG"]))**p["n2"])),
            {"TetR": 1}
        ),
        Reaction("LacI_degradation",
                 lambda s, p: p["delta1"] * s["LacI"], {"LacI": -1}),
        Reaction("TetR_degradation",
                 lambda s, p: p["delta2"] * s["TetR"], {"TetR": -1}),
    ]

    return StochasticModel(
        species=species,
        reactions=reactions,
        initial_state={"LacI": 10, "TetR": 50},
        parameters=default_params,
    )


def build_repressilator_model(params: Optional[Dict] = None) -> StochasticModel:
    """Build stochastic model for the repressilator."""
    default_params = {
        "alpha": 3.0,
        "alpha0": 0.01,
        "K": 40.0,
        "n": 2.0,
        "delta_m": 0.1,   # mRNA degradation
        "delta_p": 0.02,  # protein degradation
        "beta": 0.5,      # translation rate
    }
    if params:
        default_params.update(params)

    species = ["mRNA_lacI", "mRNA_tetR", "mRNA_cI",
               "LacI", "TetR", "cI"]

    def hill_repression(repressor_level, p):
        return (p["alpha0"] + p["alpha"] *
                (p["K"]**p["n"]) /
                (p["K"]**p["n"] + repressor_level**p["n"]))

    reactions = [
        # Transcription
        Reaction("transcribe_lacI",
                 lambda s, p: (p["alpha0"] + p["alpha"] *
                               (p["K"]**p["n"]) /
                               (p["K"]**p["n"] + s["cI"]**p["n"])),
                 {"mRNA_lacI": 1}),
        Reaction("transcribe_tetR",
                 lambda s, p: (p["alpha0"] + p["alpha"] *
                               (p["K"]**p["n"]) /
                               (p["K"]**p["n"] + s["LacI"]**p["n"])),
                 {"mRNA_tetR": 1}),
        Reaction("transcribe_cI",
                 lambda s, p: (p["alpha0"] + p["alpha"] *
                               (p["K"]**p["n"]) /
                               (p["K"]**p["n"] + s["TetR"]**p["n"])),
                 {"mRNA_cI": 1}),
        # Translation
        Reaction("translate_LacI",
                 lambda s, p: p["beta"] * s["mRNA_lacI"],
                 {"LacI": 1}),
        Reaction("translate_TetR",
                 lambda s, p: p["beta"] * s["mRNA_tetR"],
                 {"TetR": 1}),
        Reaction("translate_cI",
                 lambda s, p: p["beta"] * s["mRNA_cI"],
                 {"cI": 1}),
        # mRNA degradation
        Reaction("degrade_mRNA_lacI",
                 lambda s, p: p["delta_m"] * s["mRNA_lacI"],
                 {"mRNA_lacI": -1}),
        Reaction("degrade_mRNA_tetR",
                 lambda s, p: p["delta_m"] * s["mRNA_tetR"],
                 {"mRNA_tetR": -1}),
        Reaction("degrade_mRNA_cI",
                 lambda s, p: p["delta_m"] * s["mRNA_cI"],
                 {"mRNA_cI": -1}),
        # Protein degradation
        Reaction("degrade_LacI",
                 lambda s, p: p["delta_p"] * s["LacI"],
                 {"LacI": -1}),
        Reaction("degrade_TetR",
                 lambda s, p: p["delta_p"] * s["TetR"],
                 {"TetR": -1}),
        Reaction("degrade_cI",
                 lambda s, p: p["delta_p"] * s["cI"],
                 {"cI": -1}),
    ]

    return StochasticModel(
        species=species,
        reactions=reactions,
        initial_state={
            "mRNA_lacI": 5, "mRNA_tetR": 0, "mRNA_cI": 0,
            "LacI": 50, "TetR": 10, "cI": 10,
        },
        parameters=default_params,
    )
