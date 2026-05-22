"""Self-contained infectious disease modeling orchestrator framework.

This module provides:
- ModelStructureSelector: decision framework for choosing epidemic models
- FrameworkPipeline: end-to-end orchestration pipeline
- Scenario analysis engine
- Pathogen configuration templates
- A runnable demonstration in the __main__ block

The implementation is intentionally self-contained and uses only the Python
standard library so it can run as a standalone framework demo.
"""

from __future__ import annotations

import copy
import json
import math
import random
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


PATHOGEN_LIBRARY: Dict[str, Dict[str, Any]] = {
    "covid19": {
        "display_name": "COVID-19",
        "latent_period_significant": True,
        "default_parameters": {
            "r0": 3.0,
            "latent_days": 3.5,
            "infectious_days": 6.0,
            "ifr": 0.006,
            "beta": 0.50,
            "sigma": 1.0 / 3.5,
            "gamma": 1.0 / 6.0,
            "travel_rate": 0.03,
            "waning_immunity_days": 240,
        },
        "parameter_ranges": {
            "r0": [2.0, 8.0],
            "latent_days": [2.0, 5.0],
            "infectious_days": [5.0, 10.0],
            "ifr": [0.002, 0.015],
            "serial_interval_days": [4.0, 6.5],
        },
        "notes": "Typical literature-based ranges spanning ancestral to highly transmissible variants.",
    },
    "influenza": {
        "display_name": "Influenza",
        "latent_period_significant": True,
        "default_parameters": {
            "r0": 1.5,
            "latent_days": 1.5,
            "infectious_days": 4.0,
            "ifr": 0.001,
            "beta": 0.375,
            "sigma": 1.0 / 1.5,
            "gamma": 1.0 / 4.0,
            "travel_rate": 0.04,
            "waning_immunity_days": 365,
        },
        "parameter_ranges": {
            "r0": [1.2, 2.0],
            "latent_days": [1.0, 2.0],
            "infectious_days": [3.0, 5.0],
            "ifr": [0.0001, 0.002],
            "serial_interval_days": [2.0, 4.0],
        },
        "notes": "Representative seasonal influenza ranges; pandemic strains can exceed these values.",
    },
    "measles": {
        "display_name": "Measles",
        "latent_period_significant": True,
        "default_parameters": {
            "r0": 15.0,
            "latent_days": 10.0,
            "infectious_days": 7.0,
            "ifr": 0.002,
            "beta": 2.14,
            "sigma": 1.0 / 10.0,
            "gamma": 1.0 / 7.0,
            "travel_rate": 0.01,
            "waning_immunity_days": 3650,
        },
        "parameter_ranges": {
            "r0": [12.0, 18.0],
            "latent_days": [8.0, 12.0],
            "infectious_days": [6.0, 8.0],
            "ifr": [0.001, 0.005],
            "serial_interval_days": [11.0, 14.0],
        },
        "notes": "Typical values for fully susceptible populations without high vaccine coverage.",
    },
    "ebola": {
        "display_name": "Ebola",
        "latent_period_significant": True,
        "default_parameters": {
            "r0": 1.9,
            "latent_days": 8.0,
            "infectious_days": 10.0,
            "ifr": 0.45,
            "beta": 0.19,
            "sigma": 1.0 / 8.0,
            "gamma": 1.0 / 10.0,
            "travel_rate": 0.005,
            "waning_immunity_days": 3650,
        },
        "parameter_ranges": {
            "r0": [1.5, 2.5],
            "latent_days": [5.0, 12.0],
            "infectious_days": [7.0, 14.0],
            "ifr": [0.25, 0.70],
            "serial_interval_days": [12.0, 18.0],
        },
        "notes": "Approximate outbreak-era values with large uncertainty across settings and case ascertainment.",
    },
}


def _normalize_pathogen(pathogen: Optional[str]) -> str:
    if not pathogen:
        return "covid19"
    key = str(pathogen).strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "covid": "covid19",
        "covid19": "covid19",
        "sarscov2": "covid19",
        "flu": "influenza",
    }
    return aliases.get(key, key)


def get_config_template(pathogen: str = "covid19") -> Dict[str, Any]:
    """Return a pre-filled configuration template for known pathogens."""
    key = _normalize_pathogen(pathogen)
    if key not in PATHOGEN_LIBRARY:
        supported = ", ".join(sorted(PATHOGEN_LIBRARY))
        raise ValueError(f"Unknown pathogen '{pathogen}'. Supported values: {supported}")

    spec = copy.deepcopy(PATHOGEN_LIBRARY[key])
    params = spec["default_parameters"]
    return {
        "pathogen": key,
        "pathogen_display_name": spec["display_name"],
        "description": spec["notes"],
        "population_size": 100_000,
        "simulation_days": 180,
        "time_step_days": 1.0,
        "n_age_groups": 3,
        "spatial_patches": 1,
        "initial_state": {
            "initial_exposed": 25,
            "initial_infectious": 15,
            "initial_recovered": 0,
        },
        "model_defaults": params,
        "literature_parameter_ranges": spec["parameter_ranges"],
        "default_requirements": {
            "population_size": 100_000,
            "age_groups": key in {"covid19", "influenza", "measles"},
            "spatial_patches": 1,
            "network_effects": False,
            "stochastic": False,
            "real_time": False,
            "intervention_types": ["vaccination", "testing", "contact_reduction"],
            "data_availability": "moderate",
            "computational_budget": "medium",
            "pathogen": key,
            "comparison_type": "nested",
        },
    }


def _format_number(value: Any, digits: int = 3) -> str:
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _table_string(headers: List[str], rows: Iterable[Iterable[Any]]) -> str:
    normalized_rows = [[str(cell) for cell in row] for row in rows]
    widths = [len(h) for h in headers]
    for row in normalized_rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def make_row(values: List[str]) -> str:
        return "| " + " | ".join(v.ljust(widths[i]) for i, v in enumerate(values)) + " |"

    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    lines = [sep, make_row(headers), sep]
    for row in normalized_rows:
        lines.append(make_row(row))
    lines.append(sep)
    return "\n".join(lines)


class BaseEpidemicModel:
    """Minimal deterministic/stochastic epidemic model interface."""

    model_type = "base"

    def __init__(self, params: Dict[str, Any], requirements: Dict[str, Any]):
        self.params = copy.deepcopy(params)
        self.requirements = copy.deepcopy(requirements)
        self.population = int(requirements.get("population_size", 100_000))
        self.pathogen = _normalize_pathogen(requirements.get("pathogen"))
        self.rng = random.Random(int(self.params.get("seed", 42)))

    def clone_with_params(self, params: Dict[str, Any]) -> "BaseEpidemicModel":
        return self.__class__(params=params, requirements=self.requirements)

    def run(self, days: int, dt: float = 1.0, interventions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        raise NotImplementedError

    @staticmethod
    def summarize(history: Dict[str, List[float]]) -> Dict[str, float]:
        infectious = history["I"]
        cumulative = history["cumulative_infections"]
        reff = history["Reff"]
        peak_value = max(infectious)
        peak_day = infectious.index(peak_value)
        return {
            "peak_infectious": float(peak_value),
            "peak_day": int(peak_day),
            "total_infected": float(cumulative[-1]),
            "final_reff": float(reff[-1]),
        }

    def _initial_current_params(self) -> Dict[str, Any]:
        current = copy.deepcopy(self.params)
        current.setdefault("beta", current.get("r0", 2.0) / max(current.get("infectious_days", 5.0), 1e-6))
        current.setdefault("gamma", 1.0 / max(current.get("infectious_days", 5.0), 1e-6))
        current.setdefault("sigma", 1.0 / max(current.get("latent_days", 3.0), 1e-6))
        return current

    @staticmethod
    def _apply_interventions(day: int, current_params: Dict[str, Any], interventions: Optional[List[Dict[str, Any]]]) -> None:
        if not interventions:
            return
        for intervention in interventions:
            if intervention.get("day") != day:
                continue
            changes = intervention.get("changes", {})
            for key, value in changes.items():
                if key.endswith("_multiplier"):
                    target = key[: -len("_multiplier")]
                    current_params[target] = current_params.get(target, 1.0) * float(value)
                else:
                    current_params[key] = value


class SIRModel(BaseEpidemicModel):
    model_type = "SIR"

    def run(self, days: int, dt: float = 1.0, interventions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        current = self._initial_current_params()
        initial_i = float(self.params.get("initial_infectious", 10))
        initial_r = float(self.params.get("initial_recovered", 0))
        s = float(self.population - initial_i - initial_r)
        i = initial_i
        r = initial_r
        cumulative = self.population - s
        history = {"day": [], "S": [], "E": [], "I": [], "R": [], "Reff": [], "cumulative_infections": []}

        for day in range(days + 1):
            history["day"].append(day)
            history["S"].append(s)
            history["E"].append(0.0)
            history["I"].append(i)
            history["R"].append(r)
            history["Reff"].append((current["beta"] / max(current["gamma"], 1e-9)) * (s / self.population))
            history["cumulative_infections"].append(cumulative)
            if day == days:
                break

            self._apply_interventions(day, current, interventions)
            new_inf = min(s, current["beta"] * s * i / self.population * dt)
            recoveries = min(i, current["gamma"] * i * dt)
            s -= new_inf
            i += new_inf - recoveries
            r += recoveries
            cumulative += new_inf

        return {"history": history, "summary": self.summarize(history)}


class SEIRModel(BaseEpidemicModel):
    model_type = "SEIR"

    def run(self, days: int, dt: float = 1.0, interventions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        current = self._initial_current_params()
        initial_e = float(self.params.get("initial_exposed", 10))
        initial_i = float(self.params.get("initial_infectious", 10))
        initial_r = float(self.params.get("initial_recovered", 0))
        s = float(self.population - initial_e - initial_i - initial_r)
        e = initial_e
        i = initial_i
        r = initial_r
        cumulative = self.population - s
        history = {"day": [], "S": [], "E": [], "I": [], "R": [], "Reff": [], "cumulative_infections": []}

        for day in range(days + 1):
            history["day"].append(day)
            history["S"].append(s)
            history["E"].append(e)
            history["I"].append(i)
            history["R"].append(r)
            history["Reff"].append((current["beta"] / max(current["gamma"], 1e-9)) * (s / self.population))
            history["cumulative_infections"].append(cumulative)
            if day == days:
                break

            self._apply_interventions(day, current, interventions)
            new_exp = min(s, current["beta"] * s * i / self.population * dt)
            new_inf = min(e, current["sigma"] * e * dt)
            recoveries = min(i, current["gamma"] * i * dt)
            s -= new_exp
            e += new_exp - new_inf
            i += new_inf - recoveries
            r += recoveries
            cumulative += new_exp

        return {"history": history, "summary": self.summarize(history)}


class AgeStructuredSEIRModel(BaseEpidemicModel):
    model_type = "age-SEIR"

    def run(self, days: int, dt: float = 1.0, interventions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        current = self._initial_current_params()
        n_age = int(self.requirements.get("n_age_groups", 3))
        weights = list(self.params.get("age_weights", [0.22, 0.60, 0.18]))
        if len(weights) < n_age:
            weights.extend([1.0 / n_age] * (n_age - len(weights)))
        weights = weights[:n_age]
        total_weight = sum(weights)
        weights = [w / total_weight for w in weights]
        contact = self.params.get(
            "contact_matrix",
            [[1.4, 0.9, 0.6], [0.9, 1.2, 0.8], [0.6, 0.8, 1.0]],
        )
        if len(contact) < n_age:
            contact = [[1.0 for _ in range(n_age)] for _ in range(n_age)]

        s = [self.population * w for w in weights]
        e = [0.0 for _ in range(n_age)]
        i = [0.0 for _ in range(n_age)]
        r = [0.0 for _ in range(n_age)]
        i[min(1, n_age - 1)] = float(self.params.get("initial_infectious", 12))
        e[min(1, n_age - 1)] = float(self.params.get("initial_exposed", 20))
        s[min(1, n_age - 1)] -= i[min(1, n_age - 1)] + e[min(1, n_age - 1)]
        cumulative = sum(i) + sum(e) + sum(r)
        history = {"day": [], "S": [], "E": [], "I": [], "R": [], "Reff": [], "cumulative_infections": []}

        for day in range(days + 1):
            total_i = sum(i)
            total_s = sum(s)
            history["day"].append(day)
            history["S"].append(total_s)
            history["E"].append(sum(e))
            history["I"].append(total_i)
            history["R"].append(sum(r))
            history["Reff"].append((current["beta"] / max(current["gamma"], 1e-9)) * (total_s / self.population))
            history["cumulative_infections"].append(cumulative)
            if day == days:
                break

            self._apply_interventions(day, current, interventions)
            next_s, next_e, next_i, next_r = s[:], e[:], i[:], r[:]
            for a in range(n_age):
                force = 0.0
                for b in range(n_age):
                    pop_b = max(s[b] + e[b] + i[b] + r[b], 1.0)
                    force += contact[a][b] * (i[b] / pop_b)
                new_exp = min(s[a], current["beta"] * s[a] * force * dt / max(n_age, 1))
                new_inf = min(e[a], current["sigma"] * e[a] * dt)
                recoveries = min(i[a], current["gamma"] * i[a] * dt)
                next_s[a] -= new_exp
                next_e[a] += new_exp - new_inf
                next_i[a] += new_inf - recoveries
                next_r[a] += recoveries
                cumulative += new_exp
            s, e, i, r = next_s, next_e, next_i, next_r

        return {"history": history, "summary": self.summarize(history)}


class MetapopulationSEIRModel(BaseEpidemicModel):
    model_type = "metapop-SEIR"

    def run(self, days: int, dt: float = 1.0, interventions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        current = self._initial_current_params()
        n_patches = max(int(self.requirements.get("spatial_patches", 1)), 1)
        use_age = bool(self.requirements.get("age_groups", False))
        n_age = int(self.requirements.get("n_age_groups", 3)) if use_age else 1
        age_weights = list(self.params.get("age_weights", [0.22, 0.60, 0.18]))[:n_age]
        if len(age_weights) < n_age:
            age_weights.extend([1.0 / n_age] * (n_age - len(age_weights)))
        age_weights = [w / sum(age_weights) for w in age_weights]
        contact = self.params.get(
            "contact_matrix",
            [[1.4, 0.9, 0.6], [0.9, 1.2, 0.8], [0.6, 0.8, 1.0]],
        )

        patch_pop = self.population / n_patches
        s = [[patch_pop * age_weights[a] for a in range(n_age)] for _ in range(n_patches)]
        e = [[0.0 for _ in range(n_age)] for _ in range(n_patches)]
        i = [[0.0 for _ in range(n_age)] for _ in range(n_patches)]
        r = [[0.0 for _ in range(n_age)] for _ in range(n_patches)]
        e[0][min(1, n_age - 1)] = float(self.params.get("initial_exposed", 20))
        i[0][min(1, n_age - 1)] = float(self.params.get("initial_infectious", 10))
        s[0][min(1, n_age - 1)] -= e[0][min(1, n_age - 1)] + i[0][min(1, n_age - 1)]
        cumulative = sum(sum(row) for row in e) + sum(sum(row) for row in i) + sum(sum(row) for row in r)
        history = {"day": [], "S": [], "E": [], "I": [], "R": [], "Reff": [], "cumulative_infections": []}

        for day in range(days + 1):
            total_s = sum(sum(row) for row in s)
            total_e = sum(sum(row) for row in e)
            total_i = sum(sum(row) for row in i)
            total_r = sum(sum(row) for row in r)
            history["day"].append(day)
            history["S"].append(total_s)
            history["E"].append(total_e)
            history["I"].append(total_i)
            history["R"].append(total_r)
            history["Reff"].append((current["beta"] / max(current["gamma"], 1e-9)) * (total_s / self.population))
            history["cumulative_infections"].append(cumulative)
            if day == days:
                break

            self._apply_interventions(day, current, interventions)
            next_s = copy.deepcopy(s)
            next_e = copy.deepcopy(e)
            next_i = copy.deepcopy(i)
            next_r = copy.deepcopy(r)
            travel = float(current.get("travel_rate", 0.02))
            global_i = total_i / self.population
            for p in range(n_patches):
                patch_total = sum(s[p]) + sum(e[p]) + sum(i[p]) + sum(r[p])
                local_i = sum(i[p]) / max(patch_total, 1.0)
                for a in range(n_age):
                    if use_age:
                        within_age_force = 0.0
                        for b in range(n_age):
                            patch_age_pop = s[p][b] + e[p][b] + i[p][b] + r[p][b]
                            within_age_force += contact[a][b] * (i[p][b] / max(patch_age_pop, 1.0))
                        force = (1.0 - travel) * within_age_force + travel * global_i * n_age
                    else:
                        force = (1.0 - travel) * local_i + travel * global_i

                    new_exp = min(s[p][a], current["beta"] * s[p][a] * force * dt / max(n_age, 1))
                    new_inf = min(e[p][a], current["sigma"] * e[p][a] * dt)
                    recoveries = min(i[p][a], current["gamma"] * i[p][a] * dt)
                    next_s[p][a] -= new_exp
                    next_e[p][a] += new_exp - new_inf
                    next_i[p][a] += new_inf - recoveries
                    next_r[p][a] += recoveries
                    cumulative += new_exp
            s, e, i, r = next_s, next_e, next_i, next_r

        return {"history": history, "summary": self.summarize(history)}


class AgentBasedModel(BaseEpidemicModel):
    model_type = "ABM"

    def _stochastic_event(self, n: int, p: float) -> int:
        n = max(int(round(n)), 0)
        p = min(max(float(p), 0.0), 1.0)
        if n == 0 or p == 0.0:
            return 0
        if p == 1.0:
            return n
        if n < 1000:
            return sum(1 for _ in range(n) if self.rng.random() < p)
        mean = n * p
        std = math.sqrt(max(mean * (1.0 - p), 1e-9))
        return max(0, min(n, int(round(self.rng.gauss(mean, std)))))

    def run(self, days: int, dt: float = 1.0, interventions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        current = self._initial_current_params()
        latent = PATHOGEN_LIBRARY.get(self.pathogen, PATHOGEN_LIBRARY["covid19"]).get("latent_period_significant", True)
        s = int(self.population - self.params.get("initial_exposed", 0) - self.params.get("initial_infectious", 10))
        e = int(self.params.get("initial_exposed", 20 if latent else 0))
        i = int(self.params.get("initial_infectious", 10))
        r = int(self.params.get("initial_recovered", 0))
        cumulative = self.population - s
        history = {"day": [], "S": [], "E": [], "I": [], "R": [], "Reff": [], "cumulative_infections": []}

        for day in range(days + 1):
            history["day"].append(day)
            history["S"].append(float(s))
            history["E"].append(float(e))
            history["I"].append(float(i))
            history["R"].append(float(r))
            history["Reff"].append((current["beta"] / max(current["gamma"], 1e-9)) * (s / self.population))
            history["cumulative_infections"].append(float(cumulative))
            if day == days:
                break

            self._apply_interventions(day, current, interventions)
            p_exposure = min(current["beta"] * (i / self.population) * dt, 1.0)
            p_infectious = min(current["sigma"] * dt, 1.0)
            p_recovery = min(current["gamma"] * dt, 1.0)
            new_exp = self._stochastic_event(s, p_exposure)
            new_inf = self._stochastic_event(e, p_infectious) if latent else new_exp
            recoveries = self._stochastic_event(i, p_recovery)
            s -= new_exp
            if latent:
                e += new_exp - new_inf
            i += new_inf - recoveries
            r += recoveries
            cumulative += new_exp

        return {"history": history, "summary": self.summarize(history)}


class ModelStructureSelector:
    """Framework for selecting appropriate epidemic model structure."""

    def _primary_model_type(self, requirements: Dict[str, Any]) -> str:
        req = self._normalize_requirements(requirements)
        pop = req["population_size"]
        age_groups = req["age_groups"]
        patches = req["spatial_patches"]
        network_effects = req["network_effects"]
        stochastic = req["stochastic"]
        latent = self._has_significant_latent_period(req["pathogen"])

        if network_effects or (stochastic and pop < 100_000):
            return "ABM"
        if age_groups and patches > 1:
            return "hybrid"
        if age_groups:
            return "age-SEIR"
        if patches > 1:
            return "metapop-SEIR"
        if latent:
            return "SEIR"
        return "SIR"

    def select_model(self, requirements: dict) -> dict:
        """
        Given research requirements, recommend model structure.

        Args:
            requirements: dict with keys:
                - population_size: int
                - age_groups: bool
                - spatial_patches: int
                - network_effects: bool
                - stochastic: bool
                - real_time: bool
                - intervention_types: list
                - data_availability: str
                - computational_budget: str
                - pathogen: str

        Returns:
            dict with model, estimation, and model comparison guidance.
        """
        req = self._normalize_requirements(requirements)
        rationale: List[str] = []
        warnings: List[str] = []
        extensions: List[str] = []

        pop = req["population_size"]
        age_groups = req["age_groups"]
        patches = req["spatial_patches"]
        network_effects = req["network_effects"]
        stochastic = req["stochastic"]
        real_time = req["real_time"]
        data_availability = req["data_availability"]
        budget = req["computational_budget"]
        pathogen = req["pathogen"]
        latent = self._has_significant_latent_period(pathogen)
        comparison_type = req.get("comparison_type")

        model_type = self._primary_model_type(req)
        if model_type == "ABM":
            rationale.append("Selected ABM because individual-level heterogeneity/network effects or small-population stochasticity are important.")
            extensions.extend(["Calibrate contact network structure", "Add superspreader heterogeneity", "Track household/workplace mixing"])
        elif model_type == "hybrid":
            rationale.append("Selected hybrid model because both age structure and multiple spatial patches are required.")
            rationale.append("Hybrid is implemented as a metapopulation SEIR with age stratification.")
            extensions.extend(["Patch-specific mobility inference", "Age-specific contact matrix estimation"])
        elif model_type == "age-SEIR":
            rationale.append("Selected age-SEIR because age structure is needed to capture heterogeneity in susceptibility/contact patterns.")
            extensions.extend(["Infer age-specific ascertainment", "Add age-specific severity states"])
        elif model_type == "metapop-SEIR":
            rationale.append("Selected metapop-SEIR because multiple spatial units require explicit between-patch transmission/mobility.")
            extensions.extend(["Estimate mobility matrix", "Add patch-level intervention triggers"])
        elif model_type == "SEIR":
            rationale.append("Selected SEIR because the pathogen has a meaningful latent period that should be represented explicitly.")
            extensions.extend(["Add asymptomatic compartment", "Model reporting delays"])
        else:
            rationale.append("Selected SIR as the parsimonious default when latent structure, age structure, spatial coupling, and network effects are not required.")
            extensions.extend(["Add observation model", "Add time-varying transmission"])

        if model_type == "ABM":
            estimation_method = "ABC"
            rationale.append("Selected ABC because ABM likelihoods are typically intractable.")
        elif real_time:
            estimation_method = "particle_filter"
            rationale.append("Selected particle filtering because online/real-time updating is required.")
        elif data_availability == "rich":
            estimation_method = "MCMC"
            rationale.append("Selected MCMC because rich data support full posterior inference.")
        elif data_availability == "moderate":
            estimation_method = "MCMC with informative priors"
            rationale.append("Selected MCMC with informative priors because data are moderate rather than fully identifying.")
        else:
            estimation_method = "ABC"
            rationale.append("Selected ABC because sparse data make robust likelihood specification difficult.")

        if comparison_type == "nested":
            model_selection = "WAIC/LOO-CV"
            rationale.append("Selected WAIC/LOO-CV for comparing nested or near-nested Bayesian model variants.")
        elif comparison_type == "non-nested":
            model_selection = "Bayes Factor"
            rationale.append("Selected Bayes Factor because the model comparison is non-nested.")
        elif budget == "low":
            model_selection = "AIC/BIC"
            rationale.append("Selected AIC/BIC because computational budget is low.")
        else:
            model_selection = "WAIC + LOO-CV"
            rationale.append("Selected WAIC + LOO-CV as the default robust Bayesian model comparison toolkit.")

        if model_type == "ABM" and pop >= 100_000:
            warnings.append("ABM at this population size may be computationally heavy; consider tau-leaping or hybrid approximations.")
        if model_type in {"hybrid", "metapop-SEIR", "age-SEIR"} and budget == "low":
            warnings.append("Requested structure may strain a low computational budget; consider reducing dimensionality.")
        if data_availability == "sparse" and model_type in {"hybrid", "metapop-SEIR", "age-SEIR", "ABM"}:
            warnings.append("Sparse data may not identify a high-dimensional model without strong priors or external constraints.")
        if real_time and estimation_method == "ABC":
            warnings.append("Real-time decision support may be slow with ABC; sequential Monte Carlo ABC or surrogate models may be needed.")
        if "vaccination" in req["intervention_types"] and not age_groups:
            warnings.append("Vaccination strategies are often age-targeted; consider age structure if prioritization matters.")

        return {
            "model_type": model_type,
            "estimation_method": estimation_method,
            "model_selection": model_selection,
            "rationale": rationale,
            "warnings": warnings,
            "suggested_extensions": sorted(set(extensions)),
            "candidate_models": self.candidate_models(req),
        }

    def candidate_models(self, requirements: Dict[str, Any]) -> List[str]:
        req = self._normalize_requirements(requirements)
        primary = req.get("selected_model") or self._primary_model_type(req)
        ordered = [primary]
        if primary == "ABM":
            ordered.extend(["SEIR", "SIR"])
        elif primary == "hybrid":
            ordered.extend(["metapop-SEIR", "age-SEIR", "SEIR"])
        elif primary == "age-SEIR":
            ordered.extend(["SEIR", "SIR"])
        elif primary == "metapop-SEIR":
            ordered.extend(["SEIR", "SIR"])
        elif primary == "SEIR":
            ordered.append("SIR")
        else:
            ordered.append("SEIR")
        deduped: List[str] = []
        for item in ordered:
            if item not in deduped:
                deduped.append(item)
        return deduped

    @staticmethod
    def _normalize_requirements(requirements: Dict[str, Any]) -> Dict[str, Any]:
        req = copy.deepcopy(requirements)
        req.setdefault("population_size", 100_000)
        req.setdefault("age_groups", False)
        req.setdefault("spatial_patches", 1)
        req.setdefault("network_effects", False)
        req.setdefault("stochastic", False)
        req.setdefault("real_time", False)
        req.setdefault("intervention_types", [])
        req.setdefault("data_availability", "moderate")
        req.setdefault("computational_budget", "medium")
        req.setdefault("pathogen", "covid19")
        req["pathogen"] = _normalize_pathogen(req["pathogen"])
        return req

    @staticmethod
    def _has_significant_latent_period(pathogen: str) -> bool:
        spec = PATHOGEN_LIBRARY.get(_normalize_pathogen(pathogen))
        return bool(spec and spec.get("latent_period_significant"))


def run_scenarios(model: BaseEpidemicModel, base_params: Dict[str, Any], scenarios_config: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Run a scenario set and return comparison metrics."""
    results = []
    days = int(base_params.get("simulation_days", 180))
    dt = float(base_params.get("time_step_days", 1.0))

    for scenario in scenarios_config:
        params = copy.deepcopy(base_params)
        params.update(copy.deepcopy(scenario.get("parameter_modifications", {})))
        scenario_model = model.clone_with_params(params)
        run_output = scenario_model.run(days=days, dt=dt, interventions=scenario.get("intervention_schedule", []))
        summary = run_output["summary"]
        results.append(
            {
                "scenario": scenario["name"],
                "peak_infectious": round(summary["peak_infectious"], 2),
                "total_infected": round(summary["total_infected"], 2),
                "peak_day": int(summary["peak_day"]),
                "Reff_final": round(summary["final_reff"], 3),
            }
        )
    return results


class FrameworkPipeline:
    """End-to-end modeling pipeline."""

    MODEL_MAP = {
        "SIR": SIRModel,
        "SEIR": SEIRModel,
        "age-SEIR": AgeStructuredSEIRModel,
        "metapop-SEIR": MetapopulationSEIRModel,
        "hybrid": MetapopulationSEIRModel,
        "ABM": AgentBasedModel,
    }

    def __init__(self, selector: Optional[ModelStructureSelector] = None):
        self.selector = selector or ModelStructureSelector()

    def initialize_model(self, selection: Dict[str, Any], requirements: Dict[str, Any], parameter_overrides: Optional[Dict[str, Any]] = None) -> BaseEpidemicModel:
        pathogen_template = get_config_template(requirements.get("pathogen", "covid19"))
        params = copy.deepcopy(pathogen_template["model_defaults"])
        params.update(copy.deepcopy(pathogen_template["initial_state"]))
        params.update(
            {
                "simulation_days": pathogen_template["simulation_days"],
                "time_step_days": pathogen_template["time_step_days"],
                "seed": 42,
            }
        )
        if parameter_overrides:
            params.update(copy.deepcopy(parameter_overrides))
        model_cls = self.MODEL_MAP[selection["model_type"]]
        enriched_requirements = copy.deepcopy(requirements)
        enriched_requirements.setdefault("n_age_groups", pathogen_template.get("n_age_groups", 3))
        return model_cls(params=params, requirements=enriched_requirements)

    def estimate_parameters(self, model: BaseEpidemicModel, data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Very small placeholder estimator for demonstration purposes."""
        estimated = copy.deepcopy(model.params)
        if not data:
            return estimated

        cases = data.get("cases") or []
        if len(cases) >= 2 and cases[0] > 0:
            growth_factors = []
            for prev, curr in zip(cases[:-1], cases[1:]):
                if prev > 0:
                    growth_factors.append(curr / prev)
            if growth_factors:
                avg_growth = sum(growth_factors) / len(growth_factors)
                gamma = estimated.get("gamma", 1.0 / max(estimated.get("infectious_days", 5.0), 1e-6))
                inferred_r0 = max(0.8, min(estimated.get("r0", 2.5) * avg_growth ** 0.5, 20.0))
                estimated["r0"] = inferred_r0
                estimated["beta"] = inferred_r0 * gamma

        if data.get("hospitalization_ratio") is not None:
            estimated["hospitalization_ratio"] = float(data["hospitalization_ratio"])
        return estimated

    def select_best_candidate(self, requirements: Dict[str, Any], candidates: List[str], data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run lightweight model comparison heuristics or fit proxies across candidates."""
        budget = requirements.get("computational_budget", "medium")
        data_availability = requirements.get("data_availability", "moderate")
        complexity_penalty = {"SIR": 1.0, "SEIR": 1.2, "age-SEIR": 1.6, "metapop-SEIR": 1.8, "hybrid": 2.2, "ABM": 2.4}
        support_bonus = {"rich": 1.0, "moderate": 0.6, "sparse": 0.2}
        budget_penalty = {"high": 0.0, "medium": 0.15, "low": 0.35}
        rows = []
        cases = (data or {}).get("cases") or []

        for candidate in candidates:
            temp_selection = {
                "model_type": candidate,
                "estimation_method": "ABC" if candidate == "ABM" else "MCMC",
                "model_selection": "WAIC + LOO-CV",
            }
            model = self.initialize_model(temp_selection, requirements)
            fitted_params = self.estimate_parameters(model, data)
            model = model.clone_with_params(fitted_params)
            score = support_bonus.get(data_availability, 0.6) - complexity_penalty[candidate] * budget_penalty.get(budget, 0.15)

            rmse = None
            if cases:
                sim = model.run(days=len(cases) - 1)
                predicted = sim["history"]["I"][: len(cases)]
                rmse = math.sqrt(sum((predicted[i] - cases[i]) ** 2 for i in range(len(cases))) / len(cases))
                score -= rmse / max(sum(cases) / len(cases), 1.0)
            rows.append(
                {
                    "model": candidate,
                    "score": round(score, 4),
                    "rmse": None if rmse is None else round(rmse, 3),
                    "criterion": "fit proxy + complexity penalty",
                }
            )

        best = max(rows, key=lambda row: row["score"])
        return {"best_model": best["model"], "comparison_table": rows}

    def default_scenarios(self, requirements: Dict[str, Any]) -> List[Dict[str, Any]]:
        interventions = set(requirements.get("intervention_types", []))
        scenarios = [
            {"name": "Baseline", "parameter_modifications": {}, "intervention_schedule": []},
            {
                "name": "Moderate distancing",
                "parameter_modifications": {},
                "intervention_schedule": [{"day": 20, "changes": {"beta_multiplier": 0.80}}],
            },
            {
                "name": "Strong suppression",
                "parameter_modifications": {},
                "intervention_schedule": [{"day": 15, "changes": {"beta_multiplier": 0.60}}],
            },
        ]
        if "vaccination" in interventions:
            scenarios.append(
                {
                    "name": "Vaccination rollout",
                    "parameter_modifications": {"initial_recovered": 0.08 * requirements.get("population_size", 100_000)},
                    "intervention_schedule": [{"day": 30, "changes": {"beta_multiplier": 0.85}}],
                }
            )
        if "testing" in interventions:
            scenarios.append(
                {
                    "name": "Testing and isolation",
                    "parameter_modifications": {},
                    "intervention_schedule": [{"day": 10, "changes": {"gamma_multiplier": 1.20, "beta_multiplier": 0.90}}],
                }
            )
        return scenarios

    def run(self, requirements, data=None):
        """
        1. Select model structure
        2. Initialize model with default/provided parameters
        3. If data provided: estimate parameters
        4. Run model selection if multiple candidates
        5. Generate scenario analyses
        6. Return comprehensive results
        """
        selection = self.selector.select_model(requirements)
        model = self.initialize_model(selection, requirements)
        estimated_params = self.estimate_parameters(model, data)
        model = model.clone_with_params(estimated_params)
        candidates = selection.get("candidate_models", [selection["model_type"]])
        comparison = self.select_best_candidate(requirements, candidates, data) if len(candidates) > 1 else None
        if comparison and comparison["best_model"] != selection["model_type"]:
            selection["rationale"].append(
                f"Model comparison found {comparison['best_model']} to be a simpler competitive candidate, but the selected structure was retained because it better matches the stated requirements."
            )

        scenarios = self.default_scenarios(requirements)
        scenario_results = run_scenarios(model, model.params, scenarios)
        baseline_run = model.run(days=int(model.params.get("simulation_days", 180)), dt=float(model.params.get("time_step_days", 1.0)))
        return {
            "selection": selection,
            "parameters": estimated_params,
            "baseline_summary": baseline_run["summary"],
            "model_comparison": comparison,
            "scenario_results": scenario_results,
        }


def print_framework_summary() -> None:
    """Print the complete decision framework as a formatted table."""
    print("\nMODEL TYPE DECISION FRAMEWORK")
    model_rows = [
        ["network_effects or (stochastic and pop < 100K)", "ABM", "Need heterogeneity/network-aware stochastic representation in smaller populations"],
        ["age_groups and spatial_patches > 1", "hybrid", "Requires metapopulation SEIR with age structure"],
        ["age_groups", "age-SEIR", "Need age-specific mixing, susceptibility, or intervention targeting"],
        ["spatial_patches > 1", "metapop-SEIR", "Need multi-patch transmission and mobility"],
        ["latent period significant for pathogen", "SEIR", "Explicit exposed compartment improves realism"],
        ["otherwise", "SIR", "Parsimonious baseline when extra structure is unnecessary"],
    ]
    print(_table_string(["Decision criterion", "Recommended model", "Rationale"], model_rows))

    print("\nESTIMATION METHOD FRAMEWORK")
    estimation_rows = [
        ["model_type == ABM", "ABC", "Likelihood often intractable"],
        ["real_time", "Particle Filter", "Sequential online estimation"],
        ["data_availability == rich", "MCMC", "Supports full posterior inference"],
        ["data_availability == moderate", "MCMC with informative priors", "Regularizes partially identified models"],
        ["data_availability == sparse", "ABC", "Useful when likelihood/model mismatch is hard to specify"],
    ]
    print(_table_string(["Decision criterion", "Recommended method", "Rationale"], estimation_rows))

    print("\nMODEL SELECTION FRAMEWORK")
    selection_rows = [
        ["comparing nested models", "WAIC/LOO-CV", "Robust predictive performance comparison"],
        ["comparing non-nested models", "Bayes Factor", "Direct evidence ratio across distinct structures"],
        ["computational_budget == low", "AIC/BIC", "Fast approximate comparison"],
        ["default", "WAIC + LOO-CV", "Use both when feasible for stable selection"],
    ]
    print(_table_string(["Decision criterion", "Criterion", "Rationale"], selection_rows))


def _print_scenario_table(scenario_results: List[Dict[str, Any]]) -> None:
    rows = [
        [row["scenario"], _format_number(row["peak_infectious"], 2), _format_number(row["total_infected"], 2), row["peak_day"], _format_number(row["Reff_final"], 3)]
        for row in scenario_results
    ]
    print(_table_string(["Scenario", "Peak I", "Total infected", "Peak day", "Final Reff"], rows))


if __name__ == "__main__":
    demo_requirements = {
        "population_size": 250_000,
        "age_groups": True,
        "n_age_groups": 3,
        "spatial_patches": 3,
        "network_effects": False,
        "stochastic": False,
        "real_time": True,
        "intervention_types": ["vaccination", "testing", "contact_reduction"],
        "data_availability": "moderate",
        "computational_budget": "medium",
        "pathogen": "covid19",
        "comparison_type": "nested",
    }

    print("INFECTIOUS DISEASE MODELING FRAMEWORK DEMO")
    print("=" * 60)

    selector = ModelStructureSelector()
    decision = selector.select_model(demo_requirements)
    print("\n1) MODEL SELECTION DECISION")
    print(json.dumps(decision, indent=2))

    pipeline = FrameworkPipeline(selector=selector)
    demo_data = {
        "cases": [20, 24, 31, 40, 52, 67, 84, 103, 120, 138],
        "hospitalization_ratio": 0.05,
    }
    results = pipeline.run(demo_requirements, data=demo_data)

    print("\n2) PIPELINE RESULTS")
    print(json.dumps({
        "selection": results["selection"],
        "baseline_summary": results["baseline_summary"],
        "model_comparison": results["model_comparison"],
    }, indent=2))

    print("\n3) SCENARIO COMPARISON")
    _print_scenario_table(results["scenario_results"])

    print("\n4) FRAMEWORK SUMMARY")
    print_framework_summary()

    template = get_config_template("covid19")
    print("\n5) COVID-19 CONFIG TEMPLATE")
    print(json.dumps(template, indent=2))

    output_path = Path(__file__).with_name("framework_demo_results.json")
    output_path.write_text(json.dumps({
        "requirements": demo_requirements,
        "results": results,
        "config_template": template,
    }, indent=2), encoding="utf-8")
    print(f"\nSaved demo results to: {output_path}")
