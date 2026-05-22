from __future__ import annotations

import json
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, RESULTS_DIR, SEED

MISSION_JSON = RESULTS_DIR / "optimal_mission_sequence.json"
MISSION_FIG = FIGURES_DIR / "mission_optimization.png"


plt.style.use("seaborn-v0_8-whitegrid")


@dataclass
class MissionSolution:
    route: list[str]
    total_delta_v: float
    total_time_days: float
    objective: float
    score_sum: float
    targets_serviced: int


def route_metrics(
    order: list[str],
    dv_df: pd.DataFrame,
    time_df: pd.DataFrame,
    score_map: dict[str, float],
    fuel_budget_m_s: float,
    start: str = "CHASER-INITIAL",
    lambda_time: float = 1.2,
    reward_scale: float = 1400.0,
) -> MissionSolution:
    current = start
    route: list[str] = []
    total_dv = 0.0
    total_time = 0.0
    score_sum = 0.0

    for target in order:
        leg_dv = float(dv_df.loc[current, target])
        leg_time = float(time_df.loc[current, target])
        if total_dv + leg_dv > fuel_budget_m_s:
            break
        route.append(target)
        total_dv += leg_dv
        total_time += leg_time
        score_sum += float(score_map[target])
        current = target

    if not route:
        objective = 1e9
    else:
        objective = total_dv + lambda_time * total_time - reward_scale * score_sum - 180.0 * len(route)
    return MissionSolution(route=route, total_delta_v=total_dv, total_time_days=total_time, objective=objective, score_sum=score_sum, targets_serviced=len(route))


def nearest_neighbor(names: list[str], dv_df: pd.DataFrame, start: str = "CHASER-INITIAL") -> list[str]:
    unvisited = set(names)
    current = start
    route = []
    while unvisited:
        nxt = min(unvisited, key=lambda name: dv_df.loc[current, name])
        route.append(nxt)
        unvisited.remove(nxt)
        current = nxt
    return route


def two_opt(route: list[str], dv_df: pd.DataFrame, time_df: pd.DataFrame, score_map: dict[str, float], fuel_budget_m_s: float) -> list[str]:
    best = route[:]
    improved = True
    best_cost = route_metrics(best, dv_df, time_df, score_map, fuel_budget_m_s).objective
    while improved:
        improved = False
        for i in range(len(best) - 2):
            for j in range(i + 2, len(best) + 1):
                candidate = best[:i] + best[i:j][::-1] + best[j:]
                candidate_cost = route_metrics(candidate, dv_df, time_df, score_map, fuel_budget_m_s).objective
                if candidate_cost < best_cost:
                    best = candidate
                    best_cost = candidate_cost
                    improved = True
        route = best
    return best


def crossover(parent1: list[str], parent2: list[str], rng: np.random.Generator) -> list[str]:
    size = len(parent1)
    i, j = sorted(rng.choice(size, 2, replace=False))
    child = [None] * size
    child[i:j] = parent1[i:j]
    fill = [gene for gene in parent2 if gene not in child]
    fill_iter = iter(fill)
    return [gene if gene is not None else next(fill_iter) for gene in child]


def mutate(route: list[str], rng: np.random.Generator, rate: float = 0.25) -> list[str]:
    candidate = route[:]
    if rng.random() < rate:
        i, j = sorted(rng.choice(len(candidate), 2, replace=False))
        candidate[i], candidate[j] = candidate[j], candidate[i]
    return candidate


def genetic_algorithm(
    names: list[str],
    dv_df: pd.DataFrame,
    time_df: pd.DataFrame,
    score_map: dict[str, float],
    fuel_budget_m_s: float,
    generations: int = 220,
    pop_size: int = 80,
) -> tuple[MissionSolution, list[MissionSolution]]:
    rng = np.random.default_rng(SEED)
    seeds = [nearest_neighbor(names, dv_df)]
    for start_name in names:
        ordered = [start_name] + [name for name in nearest_neighbor([n for n in names if n != start_name], dv_df, start=start_name)]
        seeds.append(ordered)
    population = seeds + [list(rng.permutation(names)) for _ in range(max(pop_size - len(seeds), 0))]
    history: list[MissionSolution] = []

    best_solution = route_metrics(population[0], dv_df, time_df, score_map, fuel_budget_m_s)
    for _ in range(generations):
        scored = sorted((route_metrics(route, dv_df, time_df, score_map, fuel_budget_m_s) for route in population), key=lambda x: x.objective)
        history.extend(scored[:10])
        if scored[0].objective < best_solution.objective:
            best_solution = scored[0]
        elite_routes = [population[idx] for idx in np.argsort([route_metrics(route, dv_df, time_df, score_map, fuel_budget_m_s).objective for route in population])[:14]]
        new_population = [route[:] for route in elite_routes[:6]]
        while len(new_population) < pop_size:
            parent_indices = rng.choice(len(elite_routes), 2, replace=False)
            parent1 = elite_routes[int(parent_indices[0])]
            parent2 = elite_routes[int(parent_indices[1])]
            child = crossover(parent1[:], parent2[:], rng)
            child = mutate(child, rng)
            new_population.append(child)
        population = new_population
    return best_solution, history


def pareto_front(solutions: list[MissionSolution]) -> list[MissionSolution]:
    front = []
    for solution in solutions:
        dominated = False
        for other in solutions:
            if other is solution:
                continue
            if other.total_delta_v <= solution.total_delta_v and other.total_time_days <= solution.total_time_days and (
                other.total_delta_v < solution.total_delta_v or other.total_time_days < solution.total_time_days
            ):
                dominated = True
                break
        if not dominated:
            front.append(solution)
    unique = {(tuple(sol.route), round(sol.total_delta_v, 6), round(sol.total_time_days, 6)): sol for sol in front}
    return sorted(unique.values(), key=lambda s: (s.total_delta_v, s.total_time_days, -s.targets_serviced))


def optimize_mission(top_targets: pd.DataFrame, dv_df: pd.DataFrame, time_df: pd.DataFrame, fuel_budget_m_s: float = 2000.0) -> dict[str, object]:
    selected = top_targets.head(10).copy()
    target_names = selected["debris_id"].tolist()
    score_map = selected.set_index("debris_id")["combined_score"].to_dict()

    heuristic_orders = []
    base_order = nearest_neighbor(target_names, dv_df)
    heuristic_orders.append(base_order)
    heuristic_orders.append(sorted(target_names, key=lambda name: score_map[name], reverse=True))
    heuristic_orders.append(sorted(target_names, key=lambda name: score_map[name] / max(float(dv_df.loc["CHASER-INITIAL", name]), 1.0), reverse=True))
    for start_name in target_names:
        remaining = [name for name in target_names if name != start_name]
        heuristic_orders.append([start_name] + nearest_neighbor(remaining, dv_df, start=start_name))

    optimized_orders = [two_opt(order, dv_df, time_df, score_map, fuel_budget_m_s) for order in heuristic_orders]
    heuristic_solutions = [route_metrics(order, dv_df, time_df, score_map, fuel_budget_m_s) for order in optimized_orders]
    heuristic_solution = min(heuristic_solutions, key=lambda sol: sol.objective)
    ga_solution, history = genetic_algorithm(target_names, dv_df, time_df, score_map, fuel_budget_m_s)

    candidate_solutions = history + heuristic_solutions + [ga_solution]
    feasible = [sol for sol in candidate_solutions if sol.total_delta_v < fuel_budget_m_s and sol.targets_serviced > 0]
    best_feasible = min(feasible, key=lambda sol: sol.objective) if feasible else min(candidate_solutions, key=lambda sol: sol.objective)
    pareto = pareto_front(feasible if feasible else candidate_solutions)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    route_df = selected.set_index("debris_id").loc[best_feasible.route].reset_index() if best_feasible.route else selected.head(0)
    if not route_df.empty:
        axes[0].plot(route_df["altitude_km"], route_df["inclination_deg"], marker="o", color="#4c72b0")
        for idx, row in route_df.iterrows():
            axes[0].annotate(f"{idx+1}:{row['debris_id']}", (row["altitude_km"], row["inclination_deg"]), fontsize=8, xytext=(3, 3), textcoords="offset points")
    axes[0].set_xlabel("Altitude (km)")
    axes[0].set_ylabel("Inclination (deg)")
    axes[0].set_title("Optimized ADR mission sequence")

    axes[1].scatter(
        [sol.total_delta_v for sol in candidate_solutions],
        [sol.total_time_days for sol in candidate_solutions],
        alpha=0.2,
        color="#999999",
        label="Candidates",
    )
    axes[1].plot(
        [sol.total_delta_v for sol in pareto],
        [sol.total_time_days for sol in pareto],
        color="#dd8452",
        linewidth=2.0,
        marker="o",
        label="Pareto front",
    )
    axes[1].scatter([best_feasible.total_delta_v], [best_feasible.total_time_days], color="crimson", s=80, label="Selected mission")
    axes[1].axvline(fuel_budget_m_s, color="black", linestyle="--", linewidth=1.0, label="Fuel budget")
    axes[1].set_xlabel("Total mission ΔV (m/s)")
    axes[1].set_ylabel("Total transfer time (days)")
    axes[1].set_title("Fuel-time trade-off")
    axes[1].legend()
    fig.savefig(MISSION_FIG, dpi=300)
    plt.close(fig)

    output = {
        "selected_route": best_feasible.route,
        "targets_serviced": best_feasible.targets_serviced,
        "route_score_sum": best_feasible.score_sum,
        "total_delta_v_m_s": best_feasible.total_delta_v,
        "total_time_days": best_feasible.total_time_days,
        "objective": best_feasible.objective,
        "fuel_budget_m_s": fuel_budget_m_s,
        "feasible": best_feasible.total_delta_v < fuel_budget_m_s,
        "heuristic_route": heuristic_solution.route,
        "heuristic_delta_v_m_s": heuristic_solution.total_delta_v,
        "ga_route": ga_solution.route,
        "ga_delta_v_m_s": ga_solution.total_delta_v,
        "pareto_front": [
            {
                "route": sol.route,
                "targets_serviced": sol.targets_serviced,
                "route_score_sum": sol.score_sum,
                "total_delta_v_m_s": sol.total_delta_v,
                "total_time_days": sol.total_time_days,
                "objective": sol.objective,
            }
            for sol in pareto[:20]
        ],
    }
    MISSION_JSON.write_text(json.dumps(output, indent=2))
    return output


if __name__ == "__main__":
    from .debris_catalog import generate_debris_catalog
    from .target_selection import score_targets
    from .orbit_transition import build_delta_v_matrix

    catalog = generate_debris_catalog()
    _, top_targets = score_targets(catalog)
    dv_df, time_df = build_delta_v_matrix(top_targets)
    result = optimize_mission(top_targets, dv_df, time_df)
    print(json.dumps(result, indent=2))
