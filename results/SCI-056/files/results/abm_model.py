"""Agent-based epidemiological modeling utilities.

This module implements a reproducible, pure-Python SEIRV agent-based model
with configurable network topology, intervention hooks, and utilities for
comparing ABM trajectories with an aggregate ODE approximation.

Performance notes
-----------------
The implementation is designed for clarity and moderate-size experiments.
For pure Python, ABM studies become increasingly expensive once the product of
population size, mean degree, simulation days, and replicate count becomes
large. In practice, tens of thousands of agents with many Monte Carlo runs can
already be slow; for populations above roughly 50,000-100,000 agents, an ODE
or hybrid approach is often preferable unless fine-grained heterogeneity is
essential.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
import math
import random
from statistics import mean, pstdev
from time import perf_counter
from typing import Any, Deque, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


AGE_GROUPS: Tuple[str, ...] = ("0-19", "20-39", "40-64", "65+")
STATE_NAMES: Tuple[str, ...] = ("S", "E", "I", "R", "V")
DEFAULT_AGE_MIXING_MATRIX: Dict[str, Dict[str, float]] = {
    "0-19": {"0-19": 1.00, "20-39": 0.65, "40-64": 0.35, "65+": 0.20},
    "20-39": {"0-19": 0.65, "20-39": 1.00, "40-64": 0.60, "65+": 0.25},
    "40-64": {"0-19": 0.35, "20-39": 0.60, "40-64": 1.00, "65+": 0.55},
    "65+": {"0-19": 0.20, "20-39": 0.25, "40-64": 0.55, "65+": 0.95},
}


@dataclass
class Agent:
    """Individual agent with epidemiological state and tracing metadata.

    Parameters
    ----------
    state:
        Epidemiological state. One of ``S``, ``E``, ``I``, ``R``, or ``V``.
    age_group:
        Demographic band: ``0-19``, ``20-39``, ``40-64``, or ``65+``.
    location:
        Patch identifier for spatially structured simulations.
    compliance:
        Continuous compliance score in the interval ``[0, 1]``.
    infection_timer / incubation_timer / recovery_timer:
        State-related counters updated once per day.
    contact_history:
        Rolling list of recent contacts for simple tracing studies.
    """

    agent_id: int
    state: str
    age_group: str
    location: int
    compliance: float
    infection_timer: int = 0
    incubation_timer: int = 0
    recovery_timer: int = 0
    vaccination_efficacy: float = 0.0
    isolated_until: int = -1
    contact_history: Deque[int] = field(default_factory=lambda: deque(maxlen=20))

    def record_contact(self, other_agent_id: int) -> None:
        """Add a contact to the rolling tracing buffer."""
        self.contact_history.append(other_agent_id)

    def is_isolated(self, day: int) -> bool:
        """Return whether the agent should avoid contacts on a given day."""
        return day <= self.isolated_until

    def susceptibility(self) -> float:
        """Return relative susceptibility, accounting for vaccination."""
        if self.state == "V":
            return max(0.0, 1.0 - self.vaccination_efficacy)
        if self.state == "S":
            return 1.0
        return 0.0


class ABMSimulator:
    """Stochastic agent-based SEIRV simulator with interventions.

    Notes
    -----
    - Daily contacts are generated from an explicit contact network modulated by
      age mixing and patch-level assortativity.
    - State transitions are stochastic Bernoulli trials using ``sigma`` and
      ``gamma``.
    - This implementation is pure Python for portability and reproducibility.
      For large populations or many replicates, performance may become limiting.
    """

    def __init__(
        self,
        population_size: int,
        demographic_distribution: Optional[Mapping[str, float]] = None,
        age_mixing_matrix: Optional[Mapping[str, Mapping[str, float]]] = None,
        topology: str = "small-world",
        mean_degree: int = 10,
        beta: float = 0.045,
        sigma: float = 1 / 4,
        gamma: float = 1 / 7,
        initial_exposed: int = 10,
        initial_infected: int = 5,
        n_patches: int = 4,
        seed: Optional[int] = None,
        contact_history_size: int = 20,
        compliance_range: Tuple[float, float] = (0.5, 1.0),
    ) -> None:
        if population_size <= 0:
            raise ValueError("population_size must be positive")
        if topology not in {"random", "small-world", "scale-free"}:
            raise ValueError("topology must be 'random', 'small-world', or 'scale-free'")
        if mean_degree <= 0:
            raise ValueError("mean_degree must be positive")
        if not (0.0 <= beta <= 1.0 and 0.0 <= sigma <= 1.0 and 0.0 <= gamma <= 1.0):
            raise ValueError("beta, sigma, and gamma must be in [0, 1]")

        self.population_size = population_size
        self.topology = topology
        self.mean_degree = mean_degree
        self.beta = beta
        self.sigma = sigma
        self.gamma = gamma
        self.initial_exposed = max(0, initial_exposed)
        self.initial_infected = max(0, initial_infected)
        self.n_patches = max(1, n_patches)
        self.seed = seed
        self.rng = random.Random(seed)
        random.seed(seed)
        self.contact_history_size = contact_history_size
        self.day = 0
        self.edge_activation = 0.30
        self.same_patch_bonus = 1.15
        self.demographic_distribution = self._normalize_distribution(
            demographic_distribution
            or {"0-19": 0.24, "20-39": 0.30, "40-64": 0.30, "65+": 0.16}
        )
        self.age_mixing_matrix = self._normalize_age_mixing_matrix(
            age_mixing_matrix or DEFAULT_AGE_MIXING_MATRIX
        )
        self.lockdown_policy: Optional[Dict[str, float]] = None
        self.vaccination_policy: Optional[Dict[str, Any]] = None
        self.testing_policy: Optional[Dict[str, float]] = None
        self.daily_counts: List[Dict[str, float]] = []
        self.agents = self._initialize_population(compliance_range)
        self.network = self._build_contact_network()
        self._seed_initial_infections()
        self._record_daily_counts(
            day=0,
            new_exposed=0,
            new_infectious=0,
            new_recovered=0,
            new_vaccinated=0,
            tests_administered=0,
            isolations=0,
        )

    def _normalize_distribution(self, distribution: Mapping[str, float]) -> Dict[str, float]:
        missing = [group for group in AGE_GROUPS if group not in distribution]
        if missing:
            raise ValueError(f"demographic_distribution missing groups: {missing}")
        total = float(sum(distribution.values()))
        if total <= 0.0:
            raise ValueError("demographic_distribution must have positive total mass")
        return {group: distribution[group] / total for group in AGE_GROUPS}

    def _normalize_age_mixing_matrix(
        self,
        matrix: Mapping[str, Mapping[str, float]],
    ) -> Dict[str, Dict[str, float]]:
        normalized = {row: {col: float(matrix[row][col]) for col in AGE_GROUPS} for row in AGE_GROUPS}
        max_value = max(value for row in normalized.values() for value in row.values())
        if max_value <= 0.0:
            raise ValueError("age_mixing_matrix must contain positive values")
        return {
            row: {col: normalized[row][col] / max_value for col in AGE_GROUPS}
            for row in AGE_GROUPS
        }

    def _draw_age_group(self) -> str:
        threshold = self.rng.random()
        cumulative = 0.0
        for age_group in AGE_GROUPS:
            cumulative += self.demographic_distribution[age_group]
            if threshold <= cumulative:
                return age_group
        return AGE_GROUPS[-1]

    def _initialize_population(self, compliance_range: Tuple[float, float]) -> List[Agent]:
        low, high = compliance_range
        if not (0.0 <= low <= high <= 1.0):
            raise ValueError("compliance_range must lie within [0, 1]")
        agents: List[Agent] = []
        for agent_id in range(self.population_size):
            agents.append(
                Agent(
                    agent_id=agent_id,
                    state="S",
                    age_group=self._draw_age_group(),
                    location=self.rng.randrange(self.n_patches),
                    compliance=self.rng.uniform(low, high),
                    contact_history=deque(maxlen=self.contact_history_size),
                )
            )
        return agents

    def _seed_initial_infections(self) -> None:
        total_seeded = min(self.population_size, self.initial_exposed + self.initial_infected)
        sampled_ids = self.rng.sample(range(self.population_size), k=total_seeded)
        for agent_id in sampled_ids[: self.initial_exposed]:
            agent = self.agents[agent_id]
            agent.state = "E"
            agent.incubation_timer = 0
            agent.infection_timer = 1
        for agent_id in sampled_ids[self.initial_exposed : self.initial_exposed + self.initial_infected]:
            agent = self.agents[agent_id]
            agent.state = "I"
            agent.recovery_timer = 0
            agent.infection_timer = 1

    def _build_contact_network(self) -> List[set[int]]:
        if self.topology == "random":
            return self._build_random_network()
        if self.topology == "small-world":
            return self._build_small_world_network()
        return self._build_scale_free_network()

    def _build_random_network(self) -> List[set[int]]:
        adjacency = [set() for _ in range(self.population_size)]
        probability = min(1.0, self.mean_degree / max(1, self.population_size - 1))
        for i in range(self.population_size):
            for j in range(i + 1, self.population_size):
                if self.rng.random() < probability:
                    adjacency[i].add(j)
                    adjacency[j].add(i)
        return adjacency

    def _build_small_world_network(self) -> List[set[int]]:
        adjacency = [set() for _ in range(self.population_size)]
        degree = min(self.mean_degree, self.population_size - 1)
        if degree % 2 == 1:
            degree -= 1
        degree = max(2, degree)
        rewire_probability = 0.10
        half_degree = degree // 2

        for i in range(self.population_size):
            for step in range(1, half_degree + 1):
                j = (i + step) % self.population_size
                adjacency[i].add(j)
                adjacency[j].add(i)

        for i in range(self.population_size):
            for step in range(1, half_degree + 1):
                j = (i + step) % self.population_size
                if i < j and self.rng.random() < rewire_probability:
                    adjacency[i].discard(j)
                    adjacency[j].discard(i)
                    candidates = [
                        node
                        for node in range(self.population_size)
                        if node != i and node not in adjacency[i]
                    ]
                    if candidates:
                        new_j = self.rng.choice(candidates)
                        adjacency[i].add(new_j)
                        adjacency[new_j].add(i)
                    else:
                        adjacency[i].add(j)
                        adjacency[j].add(i)
        return adjacency

    def _build_scale_free_network(self) -> List[set[int]]:
        adjacency = [set() for _ in range(self.population_size)]
        m = max(1, min(self.mean_degree // 2, self.population_size - 1))
        initial_nodes = min(max(3, m + 1), self.population_size)

        for i in range(initial_nodes):
            for j in range(i + 1, initial_nodes):
                adjacency[i].add(j)
                adjacency[j].add(i)

        repeated_nodes: List[int] = []
        for node in range(initial_nodes):
            repeated_nodes.extend([node] * max(1, len(adjacency[node])))

        for new_node in range(initial_nodes, self.population_size):
            targets: set[int] = set()
            while len(targets) < min(m, new_node):
                if repeated_nodes:
                    targets.add(self.rng.choice(repeated_nodes))
                else:
                    targets.add(self.rng.randrange(new_node))
            for target in targets:
                adjacency[new_node].add(target)
                adjacency[target].add(new_node)
            repeated_nodes.extend(targets)
            repeated_nodes.extend([new_node] * max(1, len(targets)))
        return adjacency

    def apply_lockdown(self, reduction_factor: float, compliance_rate: float) -> None:
        """Reduce contacts among compliant agents.

        Parameters
        ----------
        reduction_factor:
            Fractional reduction among compliant agents.
        compliance_rate:
            Population-level compliance fraction.
        """
        self.lockdown_policy = {
            "reduction_factor": min(max(reduction_factor, 0.0), 1.0),
            "compliance_rate": min(max(compliance_rate, 0.0), 1.0),
        }

    def apply_vaccination(
        self,
        daily_rate: float,
        age_priority: Sequence[str],
        efficacy: float,
    ) -> None:
        """Configure a daily age-prioritized vaccination campaign."""
        self.vaccination_policy = {
            "daily_rate": max(daily_rate, 0.0),
            "age_priority": [group for group in age_priority if group in AGE_GROUPS] or list(AGE_GROUPS),
            "efficacy": min(max(efficacy, 0.0), 1.0),
        }

    def apply_testing_isolation(self, test_rate: float, isolation_compliance: float) -> None:
        """Configure daily testing and isolation."""
        self.testing_policy = {
            "test_rate": min(max(test_rate, 0.0), 1.0),
            "isolation_compliance": min(max(isolation_compliance, 0.0), 1.0),
        }

    def _daily_contact_modifiers(self, day: int) -> List[float]:
        modifiers: List[float] = []
        for agent in self.agents:
            if agent.is_isolated(day):
                modifiers.append(0.0)
                continue
            modifier = 1.0
            if self.lockdown_policy is not None:
                adherence_probability = self.lockdown_policy["compliance_rate"] * agent.compliance
                if self.rng.random() < adherence_probability:
                    modifier *= 1.0 - self.lockdown_policy["reduction_factor"]
            modifiers.append(max(0.0, modifier))
        return modifiers

    def _generate_daily_contacts(self, day: int) -> List[Tuple[int, int]]:
        modifiers = self._daily_contact_modifiers(day)
        contacts: List[Tuple[int, int]] = []
        for i, neighbors in enumerate(self.network):
            if modifiers[i] == 0.0:
                continue
            agent_i = self.agents[i]
            for j in neighbors:
                if j <= i or modifiers[j] == 0.0:
                    continue
                agent_j = self.agents[j]
                age_weight = self.age_mixing_matrix[agent_i.age_group][agent_j.age_group]
                spatial_weight = self.same_patch_bonus if agent_i.location == agent_j.location else 0.85
                activation_probability = min(
                    1.0,
                    self.edge_activation * modifiers[i] * modifiers[j] * age_weight * spatial_weight,
                )
                if self.rng.random() < activation_probability:
                    contacts.append((i, j))
                    agent_i.record_contact(j)
                    agent_j.record_contact(i)
        return contacts

    def _vaccinate_daily(self) -> int:
        if self.vaccination_policy is None:
            return 0
        daily_rate = self.vaccination_policy["daily_rate"]
        if daily_rate == 0.0:
            return 0
        target_count = int(round(daily_rate * self.population_size)) if daily_rate <= 1 else int(daily_rate)
        if target_count <= 0:
            return 0
        efficacy = self.vaccination_policy["efficacy"]
        priority_order = {group: index for index, group in enumerate(self.vaccination_policy["age_priority"])}
        eligible = [agent for agent in self.agents if agent.state == "S"]
        self.rng.shuffle(eligible)
        eligible.sort(key=lambda agent: priority_order.get(agent.age_group, len(priority_order)))
        vaccinated = 0
        for agent in eligible[:target_count]:
            agent.state = "V"
            agent.vaccination_efficacy = efficacy
            agent.infection_timer = 0
            agent.incubation_timer = 0
            agent.recovery_timer = 0
            vaccinated += 1
        return vaccinated

    def _testing_and_isolation(self, day: int) -> Tuple[int, int]:
        if self.testing_policy is None:
            return 0, 0
        test_rate = self.testing_policy["test_rate"]
        if test_rate <= 0.0:
            return 0, 0
        test_count = min(self.population_size, int(round(test_rate * self.population_size)))
        if test_count <= 0:
            return 0, 0
        sampled = self.rng.sample(self.agents, k=test_count)
        isolations = 0
        duration = max(1, math.ceil(1.0 / max(self.gamma, 1e-9)))
        for agent in sampled:
            positive = agent.state == "I" or (agent.state == "E" and agent.incubation_timer >= 1)
            if positive:
                compliance_probability = self.testing_policy["isolation_compliance"] * agent.compliance
                if self.rng.random() < compliance_probability:
                    agent.isolated_until = max(agent.isolated_until, day + duration)
                    isolations += 1
        return test_count, isolations

    def _transmission_probability(self, source: Agent, target: Agent) -> float:
        if source.state != "I":
            return 0.0
        susceptibility = target.susceptibility()
        if susceptibility <= 0.0:
            return 0.0
        behavioral_modifier = 1.0 - 0.20 * ((source.compliance + target.compliance) / 2.0)
        return max(0.0, min(1.0, self.beta * susceptibility * behavioral_modifier))

    def _apply_transmission(self, contacts: Iterable[Tuple[int, int]]) -> int:
        new_exposed: set[int] = set()
        for source_id, target_id in contacts:
            source = self.agents[source_id]
            target = self.agents[target_id]
            if target.state in {"S", "V"}:
                probability = self._transmission_probability(source, target)
                if self.rng.random() < probability:
                    new_exposed.add(target_id)
            if source.state in {"S", "V"}:
                probability = self._transmission_probability(target, source)
                if self.rng.random() < probability:
                    new_exposed.add(source_id)
        for agent_id in new_exposed:
            agent = self.agents[agent_id]
            if agent.state in {"S", "V"}:
                agent.state = "E"
                agent.infection_timer = 1
                agent.incubation_timer = 0
                agent.recovery_timer = 0
        return len(new_exposed)

    def _advance_states(self) -> Tuple[int, int]:
        became_infectious = 0
        recovered = 0
        for agent in self.agents:
            if agent.state == "E":
                agent.infection_timer += 1
                agent.incubation_timer += 1
                if self.rng.random() < self.sigma:
                    agent.state = "I"
                    agent.recovery_timer = 0
                    became_infectious += 1
            elif agent.state == "I":
                agent.infection_timer += 1
                agent.recovery_timer += 1
                if self.rng.random() < self.gamma:
                    agent.state = "R"
                    agent.isolated_until = -1
                    recovered += 1
        return became_infectious, recovered

    def _count_states(self) -> Dict[str, int]:
        counts = {state: 0 for state in STATE_NAMES}
        for agent in self.agents:
            counts[agent.state] += 1
        return counts

    def _record_daily_counts(
        self,
        day: int,
        new_exposed: int,
        new_infectious: int,
        new_recovered: int,
        new_vaccinated: int,
        tests_administered: int,
        isolations: int,
    ) -> Dict[str, float]:
        counts = self._count_states()
        entry: Dict[str, float] = {
            "day": day,
            **counts,
            "new_exposed": new_exposed,
            "new_infectious": new_infectious,
            "new_recovered": new_recovered,
            "new_vaccinated": new_vaccinated,
            "tests_administered": tests_administered,
            "isolations": isolations,
        }
        self.daily_counts.append(entry)
        return entry

    def daily_step(self) -> Dict[str, float]:
        """Advance the simulation by one day and return summary counts."""
        self.day += 1
        vaccinated = self._vaccinate_daily()
        tests_administered, isolations = self._testing_and_isolation(self.day)
        contacts = self._generate_daily_contacts(self.day)
        new_exposed = self._apply_transmission(contacts)
        new_infectious, recovered = self._advance_states()
        return self._record_daily_counts(
            day=self.day,
            new_exposed=new_exposed,
            new_infectious=new_infectious,
            new_recovered=recovered,
            new_vaccinated=vaccinated,
            tests_administered=tests_administered,
            isolations=isolations,
        )

    def simulate(self, days: int) -> List[Dict[str, float]]:
        """Run the simulation for ``days`` days and return recorded counts."""
        for _ in range(days):
            self.daily_step()
        return self.daily_counts


def recommend_model_type(
    population_size: int,
    heterogeneity_level: Any,
    stochastic_effects_important: bool,
    network_effects: bool,
    computational_budget_hours: float,
    spatial_resolution: str,
) -> Tuple[str, str]:
    """Recommend ``ODE``, ``ABM``, or ``hybrid`` with a rationale string.

    Decision rules
    --------------
    - Population < 10,000 or important stochastic effects -> ABM
    - Homogeneous mixing, large population, and no network effects -> ODE
    - Network effects critical -> ABM
    - Sub-city spatial resolution -> ABM
    - Computational budget < 1 hour with population > 1M -> ODE
    - Conflicting requirements -> hybrid
    """

    def _heterogeneity_score(value: Any) -> float:
        if isinstance(value, (int, float)):
            return float(value)
        lookup = {"low": 0.2, "medium": 0.5, "moderate": 0.5, "high": 0.9}
        return lookup.get(str(value).strip().lower(), 0.5)

    heterogeneity_score = _heterogeneity_score(heterogeneity_level)
    spatial_text = str(spatial_resolution).strip().lower()
    fine_spatial = any(
        token in spatial_text for token in ("neighborhood", "district", "patch", "ward", "sub-city", "block")
    )

    abm_reasons: List[str] = []
    ode_reasons: List[str] = []

    if population_size < 10_000 or stochastic_effects_important:
        abm_reasons.append("population scale or stochasticity favors individual-level simulation")
    if network_effects:
        abm_reasons.append("network effects such as clustering or superspreading are important")
    if fine_spatial:
        abm_reasons.append("sub-city spatial resolution is easier to represent with agents and patches")
    if heterogeneity_score >= 0.7:
        abm_reasons.append("strong heterogeneity weakens homogeneous-mixing assumptions")

    if heterogeneity_score <= 0.3 and population_size >= 100_000 and not network_effects and not stochastic_effects_important:
        ode_reasons.append("large, fairly homogeneous population supports a compartmental approximation")
    if computational_budget_hours < 1.0 and population_size > 1_000_000:
        ode_reasons.append("very large population with limited compute budget favors ODE efficiency")
    if not network_effects and heterogeneity_score <= 0.3 and not fine_spatial:
        ode_reasons.append("homogeneous mixing without explicit network structure is consistent with ODEs")

    if abm_reasons and not ode_reasons:
        return "ABM", "; ".join(abm_reasons)
    if ode_reasons and not abm_reasons:
        return "ODE", "; ".join(ode_reasons)
    if abm_reasons and ode_reasons:
        return (
            "hybrid",
            "mixed requirements suggest ODE for the large-scale background epidemic and ABM for focal high-risk areas; "
            f"ABM drivers: {'; '.join(abm_reasons)}; ODE drivers: {'; '.join(ode_reasons)}",
        )
    return "ODE", "defaulting to ODE because no strong ABM-specific requirement was identified"


def _solve_seirv_ode(params: Mapping[str, Any], days: int) -> List[Dict[str, float]]:
    population = float(params["population_size"])
    beta = float(params.get("beta", 0.045))
    sigma = float(params.get("sigma", 1 / 4))
    gamma = float(params.get("gamma", 1 / 7))
    initial_exposed = float(params.get("initial_exposed", 10))
    initial_infected = float(params.get("initial_infected", 5))
    vaccination = params.get("vaccination", {}) or {}
    lockdown = params.get("lockdown", {}) or {}
    testing = params.get("testing", {}) or {}

    vaccination_rate = float(vaccination.get("daily_rate", 0.0))
    vaccine_efficacy = float(vaccination.get("efficacy", 0.0))
    lockdown_effect = float(lockdown.get("reduction_factor", 0.0)) * float(lockdown.get("compliance_rate", 0.0))
    testing_effect = float(testing.get("test_rate", 0.0)) * float(testing.get("isolation_compliance", 0.0))
    effective_beta = beta * (1.0 - 0.75 * lockdown_effect) * (1.0 - 0.35 * testing_effect)

    s = max(0.0, population - initial_exposed - initial_infected)
    e = initial_exposed
    i = initial_infected
    r = 0.0
    v = 0.0
    trajectory = [{"day": 0, "S": s, "E": e, "I": i, "R": r, "V": v}]

    for day in range(1, days + 1):
        infectious_pressure = effective_beta * i / max(population, 1.0)
        infections_from_s = infectious_pressure * s
        infections_from_v = infectious_pressure * (1.0 - vaccine_efficacy) * v
        vaccinations = min(s, vaccination_rate * population if vaccination_rate <= 1 else vaccination_rate)
        new_infectious = sigma * e
        new_recovered = gamma * i

        s = max(0.0, s - infections_from_s - vaccinations)
        v = max(0.0, v + vaccinations - infections_from_v)
        e = max(0.0, e + infections_from_s + infections_from_v - new_infectious)
        i = max(0.0, i + new_infectious - new_recovered)
        r = max(0.0, r + new_recovered)

        total = s + e + i + r + v
        if total > 0.0 and abs(total - population) > 1e-6:
            scale = population / total
            s *= scale
            e *= scale
            i *= scale
            r *= scale
            v *= scale
        trajectory.append({"day": day, "S": s, "E": e, "I": i, "R": r, "V": v})
    return trajectory


def _extract_trajectory(records: Sequence[Mapping[str, float]], state: str) -> List[float]:
    return [float(entry[state]) for entry in records]


def _pearson_correlation(x_values: Sequence[float], y_values: Sequence[float]) -> float:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        return 0.0
    x_mean = mean(x_values)
    y_mean = mean(y_values)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    denominator_x = math.sqrt(sum((x - x_mean) ** 2 for x in x_values))
    denominator_y = math.sqrt(sum((y - y_mean) ** 2 for y in y_values))
    denominator = denominator_x * denominator_y
    return 0.0 if denominator == 0.0 else numerator / denominator


def compare_abm_ode(params: Mapping[str, Any], n_abm_runs: int = 100) -> Dict[str, Any]:
    """Compare repeated ABM runs against an aggregate ODE solution.

    Parameters
    ----------
    params:
        Simulator configuration dictionary.
    n_abm_runs:
        Number of stochastic ABM replicates to average.

    Returns
    -------
    dict
        Plot-ready trajectories, uncertainty summaries, timing, and metrics.
    """
    if n_abm_runs <= 0:
        raise ValueError("n_abm_runs must be positive")

    days = int(params.get("days", 180))
    abm_runs: List[List[Dict[str, float]]] = []
    base_seed = params.get("seed")
    seeds = [None if base_seed is None else int(base_seed) + run for run in range(n_abm_runs)]

    abm_start = perf_counter()
    for run_index in range(n_abm_runs):
        simulator = ABMSimulator(
            population_size=int(params["population_size"]),
            demographic_distribution=params.get("demographic_distribution"),
            age_mixing_matrix=params.get("age_mixing_matrix"),
            topology=str(params.get("topology", "small-world")),
            mean_degree=int(params.get("mean_degree", 10)),
            beta=float(params.get("beta", 0.045)),
            sigma=float(params.get("sigma", 1 / 4)),
            gamma=float(params.get("gamma", 1 / 7)),
            initial_exposed=int(params.get("initial_exposed", 10)),
            initial_infected=int(params.get("initial_infected", 5)),
            n_patches=int(params.get("n_patches", 4)),
            seed=seeds[run_index],
            contact_history_size=int(params.get("contact_history_size", 20)),
            compliance_range=tuple(params.get("compliance_range", (0.5, 1.0))),
        )
        vaccination = params.get("vaccination")
        if vaccination:
            simulator.apply_vaccination(
                daily_rate=float(vaccination.get("daily_rate", 0.0)),
                age_priority=vaccination.get("age_priority", AGE_GROUPS),
                efficacy=float(vaccination.get("efficacy", 0.0)),
            )
        lockdown = params.get("lockdown")
        if lockdown:
            simulator.apply_lockdown(
                reduction_factor=float(lockdown.get("reduction_factor", 0.0)),
                compliance_rate=float(lockdown.get("compliance_rate", 0.0)),
            )
        testing = params.get("testing")
        if testing:
            simulator.apply_testing_isolation(
                test_rate=float(testing.get("test_rate", 0.0)),
                isolation_compliance=float(testing.get("isolation_compliance", 0.0)),
            )
        abm_runs.append(simulator.simulate(days))
    abm_time = perf_counter() - abm_start

    ode_start = perf_counter()
    ode_trajectory = _solve_seirv_ode(params, days)
    ode_time = perf_counter() - ode_start

    state_means: Dict[str, List[float]] = {state: [] for state in STATE_NAMES}
    state_stds: Dict[str, List[float]] = {state: [] for state in STATE_NAMES}
    for state in STATE_NAMES:
        per_day_values = list(zip(*[_extract_trajectory(run, state) for run in abm_runs]))
        state_means[state] = [mean(values) for values in per_day_values]
        state_stds[state] = [pstdev(values) for values in per_day_values]

    ode_states = {state: _extract_trajectory(ode_trajectory, state) for state in STATE_NAMES}
    infected_mae = mean(abs(a - b) for a, b in zip(state_means["I"], ode_states["I"]))
    state_mae = {
        state: mean(abs(a - b) for a, b in zip(state_means[state], ode_states[state]))
        for state in STATE_NAMES
    }
    infected_correlation = _pearson_correlation(state_means["I"], ode_states["I"])

    return {
        "days": list(range(days + 1)),
        "abm_mean": state_means,
        "abm_std": state_stds,
        "ode": ode_states,
        "metrics": {
            "mean_absolute_error_infected": infected_mae,
            "mean_absolute_error_by_state": state_mae,
            "infected_trajectory_correlation": infected_correlation,
            "abm_runtime_seconds": abm_time,
            "ode_runtime_seconds": ode_time,
            "abm_to_ode_runtime_ratio": math.inf if ode_time == 0.0 else abm_time / ode_time,
            "n_abm_runs": n_abm_runs,
        },
    }


def _format_peak_infected(records: Sequence[Mapping[str, float]]) -> Tuple[int, float]:
    peak = max(records, key=lambda entry: entry["I"])
    return int(peak["day"]), float(peak["I"])


if __name__ == "__main__":
    demo_params: Dict[str, Any] = {
        "population_size": 1000,
        "days": 180,
        "topology": "small-world",
        "mean_degree": 10,
        "beta": 0.050,
        "sigma": 1 / 5,
        "gamma": 1 / 8,
        "initial_exposed": 12,
        "initial_infected": 6,
        "n_patches": 5,
        "seed": 42,
        "vaccination": {
            "daily_rate": 0.002,
            "age_priority": ["65+", "40-64", "20-39", "0-19"],
            "efficacy": 0.80,
        },
        "lockdown": {
            "reduction_factor": 0.30,
            "compliance_rate": 0.70,
        },
        "testing": {
            "test_rate": 0.03,
            "isolation_compliance": 0.80,
        },
    }

    simulator = ABMSimulator(
        population_size=demo_params["population_size"],
        topology=demo_params["topology"],
        mean_degree=demo_params["mean_degree"],
        beta=demo_params["beta"],
        sigma=demo_params["sigma"],
        gamma=demo_params["gamma"],
        initial_exposed=demo_params["initial_exposed"],
        initial_infected=demo_params["initial_infected"],
        n_patches=demo_params["n_patches"],
        seed=demo_params["seed"],
    )
    simulator.apply_vaccination(**demo_params["vaccination"])
    simulator.apply_lockdown(**demo_params["lockdown"])
    simulator.apply_testing_isolation(**demo_params["testing"])

    abm_start = perf_counter()
    records = simulator.simulate(int(demo_params["days"]))
    single_run_seconds = perf_counter() - abm_start
    peak_day, peak_infected = _format_peak_infected(records)

    comparison = compare_abm_ode(demo_params, n_abm_runs=20)
    metrics = comparison["metrics"]

    print("ABM demo summary")
    print("----------------")
    print(f"Population: {demo_params['population_size']}")
    print(f"Days simulated: {demo_params['days']}")
    print(
        f"Final counts: S={records[-1]['S']}, E={records[-1]['E']}, I={records[-1]['I']}, "
        f"R={records[-1]['R']}, V={records[-1]['V']}"
    )
    print(f"Peak infected: day {peak_day}, count {peak_infected:.1f}")
    print(f"Single ABM runtime (s): {single_run_seconds:.4f}")
    print(f"Mean infected MAE vs ODE: {metrics['mean_absolute_error_infected']:.3f}")
    print(f"Infected trajectory correlation: {metrics['infected_trajectory_correlation']:.3f}")
    print(f"ABM comparison runtime (s): {metrics['abm_runtime_seconds']:.4f}")
    print(f"ODE runtime (s): {metrics['ode_runtime_seconds']:.6f}")
    print()
    print("Model recommendations")
    print("---------------------")

    scenarios = [
        {
            "name": "Small stochastic outbreak",
            "args": (5000, "high", True, True, 4.0, "district"),
        },
        {
            "name": "National homogeneous planning",
            "args": (5_000_000, "low", False, False, 0.5, "country"),
        },
        {
            "name": "Regional mixed strategy",
            "args": (250_000, "medium", False, True, 2.0, "city"),
        },
    ]
    for scenario in scenarios:
        recommendation, rationale = recommend_model_type(*scenario["args"])
        print(f"{scenario['name']}: {recommendation} -> {rationale}")
