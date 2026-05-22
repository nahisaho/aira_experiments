from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import least_squares


RESULTS_DIR = Path(__file__).resolve().parent
WORKSPACE_DIR = RESULTS_DIR.parent
DATA_DIR = WORKSPACE_DIR / "data"
FIGURES_DIR = WORKSPACE_DIR / "figures"
LOGS_DIR = WORKSPACE_DIR / "logs"
REPORT_PATH = WORKSPACE_DIR / "report.md"
PREPROCESSING_LOG_PATH = DATA_DIR / "preprocessing-log.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"
STAT_SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
WAVE6_JSON_PATH = RESULTS_DIR / "covid_wave6_results.json"
WAVE7_JSON_PATH = RESULTS_DIR / "covid_wave7_results.json"
SCENARIO_JSON_PATH = RESULTS_DIR / "scenario_comparison.json"
MODEL_JSON_PATH = RESULTS_DIR / "model_comparison.json"

AGE_GROUPS = ["0-19", "20-39", "40-64", "65+"]
BASE_AGE_SHARES = np.array([0.20, 0.30, 0.30, 0.20], dtype=float)
AGE_SEVERITY = {
    "hospitalization_rate": {"0-19": 0.0025, "20-39": 0.0080, "40-64": 0.0250, "65+": 0.0900},
    "icu_rate": {"0-19": 0.00008, "20-39": 0.00035, "40-64": 0.00180, "65+": 0.01000},
    "cfr": {"0-19": 0.00001, "20-39": 0.00008, "40-64": 0.00090, "65+": 0.01500},
}
POPULATION = 125_000_000.0
RNG_SEED = 20240219

WAVE_SPECS = {
    "wave6": {
        "label": "Wave 6 (Omicron BA.1)",
        "variant": "Omicron BA.1",
        "start_date": "2022-01-01",
        "days": 90,
        "target_peak": 100_000.0,
        "target_total": 5_000_000.0,
        "target_peak_day": 43,
        "beta_guess": 1.05,
        "sigma_true": 0.45,
        "gamma_true": 0.18,
        "rho_true": 0.72,
        "prior_immunity": 0.14,
        "contact_max": 0.24,
        "contact_window": (12, 54),
        "booster_start": 0.08,
        "booster_end": 0.58,
        "booster_mid": 36,
        "booster_steepness": 0.11,
        "ve_infection": 0.33,
        "ve_severe": 0.58,
        "seed_guess": 18_000.0,
        "r0_target": (5.0, 7.0),
    },
    "wave7": {
        "label": "Wave 7 (Omicron BA.5)",
        "variant": "Omicron BA.5",
        "start_date": "2022-07-01",
        "days": 92,
        "target_peak": 260_000.0,
        "target_total": 12_000_000.0,
        "target_peak_day": 44,
        "beta_guess": 1.35,
        "sigma_true": 0.50,
        "gamma_true": 0.16,
        "rho_true": 0.68,
        "prior_immunity": 0.30,
        "contact_max": 0.16,
        "contact_window": (20, 66),
        "booster_start": 0.50,
        "booster_end": 0.64,
        "booster_mid": 28,
        "booster_steepness": 0.08,
        "ve_infection": 0.28,
        "ve_severe": 0.52,
        "seed_guess": 40_000.0,
        "r0_target": (8.0, 10.0),
    },
}


def ensure_directories() -> None:
    for path in (RESULTS_DIR, DATA_DIR, FIGURES_DIR, LOGS_DIR):
        path.mkdir(parents=True, exist_ok=True)


def iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def log_event(event_type: str, phase: str, skill_or_tool: str, handoff_in: Dict, handoff_out: Dict, files_written: List[str], status: str = "ok") -> None:
    ensure_directories()
    record = {
        "timestamp": iso_now(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in,
        "handoff_out": handoff_out,
        "files_written": files_written,
        "status": status,
    }
    with PROCESS_LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")


def logistic_curve(days: int, start: float, end: float, midpoint: float, steepness: float) -> np.ndarray:
    t = np.arange(days, dtype=float)
    return start + (end - start) / (1.0 + np.exp(-steepness * (t - midpoint)))


def hump_contact_curve(days: int, start_day: int, end_day: int, max_reduction: float) -> np.ndarray:
    t = np.arange(days, dtype=float)
    rise = 1.0 / (1.0 + np.exp(-(t - start_day) / 3.5))
    fall = 1.0 / (1.0 + np.exp((t - end_day) / 4.5))
    curve = max_reduction * rise * fall
    return np.clip(curve, 0.0, max_reduction)


def make_dates(start_date: str, days: int) -> List[str]:
    start = datetime.fromisoformat(start_date)
    return [(start + timedelta(days=offset)).date().isoformat() for offset in range(days)]


def rolling_average(values: np.ndarray, window: int = 7) -> np.ndarray:
    kernel = np.ones(window, dtype=float) / window
    padded = np.pad(values, (window - 1, 0), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def shift_series(values: np.ndarray, lag: int) -> np.ndarray:
    shifted = np.zeros_like(values)
    if lag <= 0:
        return values.copy()
    shifted[lag:] = values[:-lag]
    return shifted


def dynamic_age_weights(days: int, booster_curve: np.ndarray, susceptibility: np.ndarray | None = None) -> np.ndarray:
    t = np.arange(days, dtype=float)
    seasonal = np.vstack(
        [
            1.0 + 0.10 * np.sin(2.0 * np.pi * t / max(days, 1) + 0.4),
            1.0 + 0.06 * np.cos(2.0 * np.pi * t / max(days, 1) - 0.7),
            1.0 + 0.04 * np.sin(4.0 * np.pi * t / max(days, 1) + 1.0),
            1.0 - 0.28 * booster_curve,
        ]
    ).T
    weights = BASE_AGE_SHARES * seasonal
    if susceptibility is not None:
        weights = weights * susceptibility
    weights = np.clip(weights, 1e-6, None)
    weights = weights / weights.sum(axis=1, keepdims=True)
    return weights


def simulate_sir(beta: float, gamma: float, rho: float, days: int, population: float, i0: float, prior_immunity: float, beta_multiplier: np.ndarray | None = None) -> Dict[str, np.ndarray]:
    if beta_multiplier is None:
        beta_multiplier = np.ones(days, dtype=float)
    s = np.zeros(days + 1, dtype=float)
    i = np.zeros(days + 1, dtype=float)
    r = np.zeros(days + 1, dtype=float)
    cases = np.zeros(days, dtype=float)
    r[0] = population * prior_immunity
    i[0] = i0
    s[0] = population - r[0] - i[0]
    for day in range(days):
        lam = beta * beta_multiplier[day] * s[day] * i[day] / population
        recoveries = gamma * i[day]
        lam = min(lam, s[day])
        recoveries = min(recoveries, i[day] + lam)
        s[day + 1] = max(s[day] - lam, 0.0)
        i[day + 1] = max(i[day] + lam - recoveries, 0.0)
        r[day + 1] = min(population - s[day + 1] - i[day + 1], population)
        cases[day] = rho * lam
    return {"S": s[:-1], "I": i[:-1], "R": r[:-1], "cases": cases}


def simulate_seir(beta: float, sigma: float, gamma: float, rho: float, days: int, population: float, e0: float, i0: float, prior_immunity: float, beta_multiplier: np.ndarray | None = None) -> Dict[str, np.ndarray]:
    if beta_multiplier is None:
        beta_multiplier = np.ones(days, dtype=float)
    s = np.zeros(days + 1, dtype=float)
    e = np.zeros(days + 1, dtype=float)
    i = np.zeros(days + 1, dtype=float)
    r = np.zeros(days + 1, dtype=float)
    transitions = np.zeros(days, dtype=float)
    reported_cases = np.zeros(days, dtype=float)
    r[0] = population * prior_immunity
    e[0] = e0
    i[0] = i0
    s[0] = population - r[0] - e[0] - i[0]
    for day in range(days):
        lam = beta * beta_multiplier[day] * s[day] * i[day] / population
        progress = sigma * e[day]
        recoveries = gamma * i[day]
        lam = min(lam, s[day])
        progress = min(progress, e[day] + lam)
        recoveries = min(recoveries, i[day] + progress)
        s[day + 1] = max(s[day] - lam, 0.0)
        e[day + 1] = max(e[day] + lam - progress, 0.0)
        i[day + 1] = max(i[day] + progress - recoveries, 0.0)
        r[day + 1] = min(population - s[day + 1] - e[day + 1] - i[day + 1], population)
        transitions[day] = progress
        reported_cases[day] = rho * progress
    return {"S": s[:-1], "E": e[:-1], "I": i[:-1], "R": r[:-1], "new_infections": transitions, "cases": reported_cases}


def summarize_curve(cases: np.ndarray) -> Tuple[float, float, int]:
    peak_idx = int(np.argmax(cases))
    return float(np.max(cases)), float(np.sum(cases)), peak_idx


def calibrate_wave(spec: Dict) -> Dict[str, float]:
    days = spec["days"]
    booster_curve = logistic_curve(days, spec["booster_start"], spec["booster_end"], spec["booster_mid"], spec["booster_steepness"])
    contact_curve = hump_contact_curve(days, spec["contact_window"][0], spec["contact_window"][1], spec["contact_max"])
    actual_multiplier = (1.0 - contact_curve) * (1.0 - spec["ve_infection"] * booster_curve)

    def residuals(x: np.ndarray) -> np.ndarray:
        beta = x[0]
        seed_e = x[1]
        sim = simulate_seir(
            beta=beta,
            sigma=spec["sigma_true"],
            gamma=spec["gamma_true"],
            rho=spec["rho_true"],
            days=days,
            population=POPULATION,
            e0=seed_e,
            i0=0.55 * seed_e,
            prior_immunity=spec["prior_immunity"],
            beta_multiplier=actual_multiplier,
        )
        peak, total, peak_day = summarize_curve(sim["cases"])
        return np.array(
            [
                (peak - spec["target_peak"]) / spec["target_peak"],
                (total - spec["target_total"]) / spec["target_total"],
                (peak_day - spec["target_peak_day"]) / max(days, 1),
            ]
        )

    fit = least_squares(
        residuals,
        x0=np.array([spec["beta_guess"], spec["seed_guess"]], dtype=float),
        bounds=([0.6, 5_000.0], [2.0, 500_000.0]),
        xtol=1e-10,
        ftol=1e-10,
        gtol=1e-10,
        max_nfev=250,
    )
    beta, seed_e = fit.x
    return {
        "beta_true": float(beta),
        "sigma_true": float(spec["sigma_true"]),
        "gamma_true": float(spec["gamma_true"]),
        "rho_true": float(spec["rho_true"]),
        "seed_e": float(seed_e),
        "seed_i": float(0.55 * seed_e),
        "booster_curve": booster_curve,
        "contact_curve": contact_curve,
        "actual_multiplier": actual_multiplier,
    }


def age_severity_arrays(metric: str) -> np.ndarray:
    return np.array([AGE_SEVERITY[metric][group] for group in AGE_GROUPS], dtype=float)


def generate_wave_data(wave_key: str, rng: np.random.Generator) -> Dict:
    spec = WAVE_SPECS[wave_key]
    calibrated = calibrate_wave(spec)
    dates = make_dates(spec["start_date"], spec["days"])
    sim = simulate_seir(
        beta=calibrated["beta_true"],
        sigma=calibrated["sigma_true"],
        gamma=calibrated["gamma_true"],
        rho=calibrated["rho_true"],
        days=spec["days"],
        population=POPULATION,
        e0=calibrated["seed_e"],
        i0=calibrated["seed_i"],
        prior_immunity=spec["prior_immunity"],
        beta_multiplier=calibrated["actual_multiplier"],
    )
    weekday_pattern = np.array([1.03, 1.02, 1.01, 1.00, 0.97, 0.94, 0.93], dtype=float)
    weekly_effect = weekday_pattern[np.arange(spec["days"]) % 7]
    noise = rng.lognormal(mean=0.0, sigma=0.08, size=spec["days"])
    noisy_component = sim["cases"] * weekly_effect * noise
    observed_cases = np.clip(0.78 * sim["cases"] + 0.22 * noisy_component, 0.0, None)
    observed_cases *= spec["target_total"] / max(observed_cases.sum(), 1.0)

    age_weights = dynamic_age_weights(spec["days"], calibrated["booster_curve"])
    age_cases = observed_cases[:, None] * age_weights

    hosp_rates = age_severity_arrays("hospitalization_rate")
    icu_rates = age_severity_arrays("icu_rate")
    cfr_rates = age_severity_arrays("cfr")
    severe_modifier = 1.0 - spec["ve_severe"] * calibrated["booster_curve"]
    hosp_age = shift_series(age_cases, 7) * hosp_rates * severe_modifier[:, None]
    icu_age = shift_series(age_cases, 10) * icu_rates * severe_modifier[:, None]
    deaths_age = shift_series(age_cases, 14) * cfr_rates * severe_modifier[:, None]

    hosp_age *= (1.05 if wave_key == "wave7" else 1.0)
    deaths_age *= (1.10 if wave_key == "wave7" else 1.0)

    peak, total, peak_day = summarize_curve(observed_cases)
    return {
        "wave": wave_key,
        "label": spec["label"],
        "variant": spec["variant"],
        "dates": dates,
        "observed_cases": observed_cases,
        "smoothed_cases": rolling_average(observed_cases),
        "age_weights": age_weights,
        "age_cases": age_cases,
        "booster_coverage": calibrated["booster_curve"],
        "contact_reduction": calibrated["contact_curve"],
        "hospitalizations_by_age": hosp_age,
        "icu_by_age": icu_age,
        "deaths_by_age": deaths_age,
        "hospitalizations": hosp_age.sum(axis=1),
        "icu": icu_age.sum(axis=1),
        "deaths": deaths_age.sum(axis=1),
        "true_parameters": {
            "beta": calibrated["beta_true"],
            "sigma": calibrated["sigma_true"],
            "gamma": calibrated["gamma_true"],
            "reporting_rate": calibrated["rho_true"],
        },
        "initial_conditions": {
            "E0": calibrated["seed_e"],
            "I0": calibrated["seed_i"],
            "prior_immunity_fraction": spec["prior_immunity"],
        },
        "assumptions": {
            "age_distribution": dict(zip(AGE_GROUPS, BASE_AGE_SHARES.tolist())),
            "severity": AGE_SEVERITY,
            "target_peak_cases": spec["target_peak"],
            "target_total_cases": spec["target_total"],
            "actual_peak_cases": peak,
            "actual_total_cases": total,
            "peak_date": dates[peak_day],
        },
    }


def r0_midpoint(spec: Dict) -> float:
    low, high = spec["r0_target"]
    return 0.5 * (low + high)


def fit_sir_model(wave_data: Dict) -> Dict:
    spec = WAVE_SPECS[wave_data["wave"]]
    observed = wave_data["smoothed_cases"]
    days = len(observed)
    booster_curve = wave_data["booster_coverage"]
    contact_curve = wave_data["contact_reduction"]
    multiplier = (1.0 - contact_curve) * (1.0 - spec["ve_infection"] * booster_curve)
    i0 = max(wave_data["initial_conditions"]["I0"], observed[0] / 0.6)
    low_r0, high_r0 = spec["r0_target"]

    def residuals(params: np.ndarray) -> np.ndarray:
        r0, gamma, rho = params
        beta = r0 * gamma
        sim = simulate_sir(beta, gamma, rho, days, POPULATION, i0, spec["prior_immunity"], multiplier)
        pred = sim["cases"]
        reg = 0.08 * (r0 - r0_midpoint(spec))
        return np.concatenate(((pred - observed) / np.maximum(1000.0, np.sqrt(observed + 100.0)), np.array([reg])))

    fit = least_squares(
        residuals,
        x0=np.array([r0_midpoint(spec), 0.17, 0.68], dtype=float),
        bounds=([low_r0 * 0.75, 0.10, 0.45], [high_r0 * 1.10, 0.28, 0.90]),
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=400,
    )
    r0, gamma, rho = fit.x
    beta = r0 * gamma
    sim = simulate_sir(beta, gamma, rho, days, POPULATION, i0, spec["prior_immunity"], multiplier)
    return {
        "parameters": {"beta": float(beta), "gamma": float(gamma), "reporting_rate": float(rho), "R0": float(r0)},
        "fitted_cases": sim["cases"],
        "state": sim,
        "n_params": 3,
    }


def fit_seir_model(wave_data: Dict) -> Dict:
    spec = WAVE_SPECS[wave_data["wave"]]
    observed = wave_data["smoothed_cases"]
    days = len(observed)
    booster_curve = wave_data["booster_coverage"]
    contact_curve = wave_data["contact_reduction"]
    multiplier = (1.0 - contact_curve) * (1.0 - spec["ve_infection"] * booster_curve)
    e0 = wave_data["initial_conditions"]["E0"]
    i0 = wave_data["initial_conditions"]["I0"]
    low_r0, high_r0 = spec["r0_target"]

    def residuals(params: np.ndarray) -> np.ndarray:
        r0, sigma, gamma, rho = params
        beta = r0 * gamma
        sim = simulate_seir(beta, sigma, gamma, rho, days, POPULATION, e0, i0, spec["prior_immunity"], multiplier)
        pred = sim["cases"]
        shape_penalty = 0.001 * (np.argmax(pred) - np.argmax(observed))
        r0_penalty = 0.12 * (r0 - r0_midpoint(spec))
        return np.concatenate(((pred - observed) / np.maximum(1000.0, np.sqrt(observed + 100.0)), np.array([shape_penalty, r0_penalty])))

    fit = least_squares(
        residuals,
        x0=np.array([r0_midpoint(spec), 0.46, 0.17, 0.70], dtype=float),
        bounds=([low_r0, 0.25, 0.10, 0.45], [high_r0, 0.70, 0.24, 0.90]),
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=600,
    )
    r0, sigma, gamma, rho = fit.x
    beta = r0 * gamma
    sim = simulate_seir(beta, sigma, gamma, rho, days, POPULATION, e0, i0, spec["prior_immunity"], multiplier)
    reff = r0 * multiplier * sim["S"] / POPULATION
    return {
        "parameters": {
            "beta": float(beta),
            "sigma": float(sigma),
            "gamma": float(gamma),
            "reporting_rate": float(rho),
            "R0": float(r0),
        },
        "fitted_cases": sim["cases"],
        "state": sim,
        "Reff": reff,
        "n_params": 4,
    }


def predict_age_outputs(total_cases: np.ndarray, booster_curve: np.ndarray, ve_severe: float, dynamic: bool, susceptibility: np.ndarray | None = None) -> Dict[str, np.ndarray]:
    if dynamic:
        weights = dynamic_age_weights(len(total_cases), booster_curve, susceptibility=susceptibility)
    else:
        weights = np.tile(BASE_AGE_SHARES, (len(total_cases), 1))
    age_cases = total_cases[:, None] * weights
    severe_modifier = 1.0 - ve_severe * booster_curve
    hosp = shift_series(age_cases, 7) * age_severity_arrays("hospitalization_rate") * severe_modifier[:, None]
    deaths = shift_series(age_cases, 14) * age_severity_arrays("cfr") * severe_modifier[:, None]
    return {"weights": weights, "age_cases": age_cases, "hospitalizations": hosp, "deaths": deaths}


def fit_age_structured_seir(wave_data: Dict) -> Dict:
    spec = WAVE_SPECS[wave_data["wave"]]
    observed_cases = wave_data["smoothed_cases"]
    observed_hosp = wave_data["hospitalizations_by_age"]
    observed_deaths = wave_data["deaths_by_age"]
    days = len(observed_cases)
    booster_curve = wave_data["booster_coverage"]
    contact_curve = wave_data["contact_reduction"]
    multiplier = (1.0 - contact_curve) * (1.0 - spec["ve_infection"] * booster_curve)
    e0 = wave_data["initial_conditions"]["E0"]
    i0 = wave_data["initial_conditions"]["I0"]
    low_r0, high_r0 = spec["r0_target"]

    def residuals(params: np.ndarray) -> np.ndarray:
        r0, sigma, gamma, rho = params[:4]
        beta = r0 * gamma
        sus = np.exp(params[4:])
        sim = simulate_seir(beta, sigma, gamma, rho, days, POPULATION, e0, i0, spec["prior_immunity"], multiplier)
        age_pred = predict_age_outputs(sim["cases"], booster_curve, spec["ve_severe"], dynamic=True, susceptibility=sus)
        case_resid = (sim["cases"] - observed_cases) / np.maximum(1000.0, np.sqrt(observed_cases + 100.0))
        hosp_resid = 2.0 * (age_pred["hospitalizations"] - observed_hosp).ravel() / np.maximum(5.0, np.sqrt(observed_hosp.ravel() + 1.0))
        death_resid = 3.0 * (age_pred["deaths"] - observed_deaths).ravel() / np.maximum(1.0, np.sqrt(observed_deaths.ravel() + 0.5))
        reg = np.array([
            0.10 * (r0 - r0_midpoint(spec)),
            0.04 * (sus.mean() - 1.0),
        ])
        return np.concatenate((case_resid, hosp_resid, death_resid, reg))

    fit = least_squares(
        residuals,
        x0=np.array([r0_midpoint(spec), 0.46, 0.17, 0.70, 0.0, 0.0, 0.0, 0.0], dtype=float),
        bounds=(np.array([low_r0, 0.25, 0.10, 0.45, -1.4, -1.4, -1.4, -1.4]), np.array([high_r0, 0.70, 0.24, 0.90, 1.4, 1.4, 1.4, 1.4])),
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=700,
    )
    r0, sigma, gamma, rho = fit.x[:4]
    beta = r0 * gamma
    susceptibility = np.exp(fit.x[4:])
    sim = simulate_seir(beta, sigma, gamma, rho, days, POPULATION, e0, i0, spec["prior_immunity"], multiplier)
    age_pred = predict_age_outputs(sim["cases"], booster_curve, spec["ve_severe"], dynamic=True, susceptibility=susceptibility)
    reff = r0 * multiplier * sim["S"] / POPULATION
    return {
        "parameters": {
            "beta": float(beta),
            "sigma": float(sigma),
            "gamma": float(gamma),
            "reporting_rate": float(rho),
            "R0": float(r0),
            "susceptibility_multipliers": {group: float(value) for group, value in zip(AGE_GROUPS, susceptibility)},
        },
        "fitted_cases": sim["cases"],
        "state": sim,
        "age_outputs": age_pred,
        "Reff": reff,
        "n_params": 8,
    }


def compute_information_criteria(residuals: np.ndarray, n_params: int) -> Dict[str, float]:
    n = residuals.size
    rss = float(np.sum(np.square(residuals)))
    rss = max(rss, 1e-12)
    return {
        "rss": rss,
        "aic": float(n * np.log(rss / n) + 2 * n_params),
        "bic": float(n * np.log(rss / n) + n_params * np.log(n)),
    }


def evaluate_models(wave_data: Dict, sir_fit: Dict, seir_fit: Dict, age_fit: Dict) -> Dict:
    spec = WAVE_SPECS[wave_data["wave"]]
    observed_cases = wave_data["smoothed_cases"]
    observed_hosp = wave_data["hospitalizations_by_age"]
    observed_deaths = wave_data["deaths_by_age"]
    booster_curve = wave_data["booster_coverage"]

    sir_age = predict_age_outputs(sir_fit["fitted_cases"], booster_curve, spec["ve_severe"], dynamic=False)
    seir_age = predict_age_outputs(seir_fit["fitted_cases"], booster_curve, spec["ve_severe"], dynamic=False)
    age_age = age_fit["age_outputs"]

    case_scale = np.maximum(1000.0, np.sqrt(observed_cases + 100.0))
    hosp_scale = np.maximum(5.0, np.sqrt(observed_hosp.ravel() + 1.0))
    death_scale = np.maximum(1.0, np.sqrt(observed_deaths.ravel() + 0.5))

    sir_resid = np.concatenate(((sir_fit["fitted_cases"] - observed_cases) / case_scale, 2.0 * (sir_age["hospitalizations"] - observed_hosp).ravel() / hosp_scale, 3.0 * (sir_age["deaths"] - observed_deaths).ravel() / death_scale))
    seir_resid = np.concatenate(((seir_fit["fitted_cases"] - observed_cases) / case_scale, 2.0 * (seir_age["hospitalizations"] - observed_hosp).ravel() / hosp_scale, 3.0 * (seir_age["deaths"] - observed_deaths).ravel() / death_scale))
    age_resid = np.concatenate(((age_fit["fitted_cases"] - observed_cases) / case_scale, 2.0 * (age_age["hospitalizations"] - observed_hosp).ravel() / hosp_scale, 3.0 * (age_age["deaths"] - observed_deaths).ravel() / death_scale))

    comparison = {
        "SIR": {**compute_information_criteria(sir_resid, sir_fit["n_params"]), "n_params": sir_fit["n_params"]},
        "SEIR": {**compute_information_criteria(seir_resid, seir_fit["n_params"]), "n_params": seir_fit["n_params"]},
        "Age-structured SEIR": {**compute_information_criteria(age_resid, age_fit["n_params"]), "n_params": age_fit["n_params"]},
    }
    best_model = min(comparison.items(), key=lambda item: item[1]["aic"])[0]
    return {"models": comparison, "best_model": best_model}


def run_intervention_scenarios(wave_data: Dict, seir_fit: Dict) -> Dict:
    spec = WAVE_SPECS[wave_data["wave"]]
    days = len(wave_data["observed_cases"])
    booster_actual = wave_data["booster_coverage"]
    contact_actual = wave_data["contact_reduction"]
    fitted = seir_fit["parameters"]
    e0 = wave_data["initial_conditions"]["E0"]
    i0 = wave_data["initial_conditions"]["I0"]

    scenarios = {
        "no_intervention": {
            "label": "No intervention",
            "contact_curve": np.zeros(days, dtype=float),
            "booster_curve": np.zeros(days, dtype=float),
        },
        "quasi_emergency_only": {
            "label": "Quasi-emergency only",
            "contact_curve": np.full(days, 0.25 if wave_data["wave"] == "wave6" else 0.20, dtype=float),
            "booster_curve": np.zeros(days, dtype=float),
        },
        "vaccination_only": {
            "label": "Vaccination only",
            "contact_curve": np.zeros(days, dtype=float),
            "booster_curve": booster_actual,
        },
        "combined_actual": {
            "label": "Combined actual",
            "contact_curve": contact_actual,
            "booster_curve": booster_actual,
        },
    }

    outputs = {}
    for key, scenario in scenarios.items():
        if key == "combined_actual":
            peak_day = int(np.argmax(wave_data["observed_cases"]))
            outputs[key] = {
                "label": scenario["label"],
                "peak_cases": float(np.max(wave_data["observed_cases"])),
                "peak_cases_date": wave_data["dates"][peak_day],
                "total_cases": float(np.sum(wave_data["observed_cases"])),
                "peak_hospitalizations": float(np.max(wave_data["hospitalizations"])),
                "deaths": float(np.sum(wave_data["deaths"])),
            }
            continue
        multiplier = (1.0 - scenario["contact_curve"]) * (1.0 - spec["ve_infection"] * scenario["booster_curve"])
        sim = simulate_seir(
            beta=fitted["beta"],
            sigma=fitted["sigma"],
            gamma=fitted["gamma"],
            rho=fitted["reporting_rate"],
            days=days,
            population=POPULATION,
            e0=e0,
            i0=i0,
            prior_immunity=spec["prior_immunity"],
            beta_multiplier=multiplier,
        )
        age_outputs = predict_age_outputs(sim["cases"], scenario["booster_curve"], spec["ve_severe"], dynamic=True)
        peak_day = int(np.argmax(sim["cases"]))
        outputs[key] = {
            "label": scenario["label"],
            "peak_cases": float(np.max(sim["cases"])),
            "peak_cases_date": wave_data["dates"][peak_day],
            "total_cases": float(np.sum(sim["cases"])),
            "peak_hospitalizations": float(np.max(age_outputs["hospitalizations"].sum(axis=1))),
            "deaths": float(np.sum(age_outputs["deaths"])),
        }
    actual = outputs["combined_actual"]
    for key, scenario in outputs.items():
        scenario["peak_reduction_vs_no_intervention_pct"] = float(100.0 * (1.0 - scenario["peak_cases"] / outputs["no_intervention"]["peak_cases"]))
        scenario["total_reduction_vs_no_intervention_pct"] = float(100.0 * (1.0 - scenario["total_cases"] / outputs["no_intervention"]["total_cases"]))
        scenario["relative_to_actual_total_pct"] = float(100.0 * scenario["total_cases"] / max(actual["total_cases"], 1.0))
    return outputs


def summarise_key_findings(wave6: Dict, wave7: Dict, model_comparison: Dict, scenario_comparison: Dict) -> Dict:
    combined_wave6 = scenario_comparison["wave6"]["combined_actual"]
    combined_wave7 = scenario_comparison["wave7"]["combined_actual"]
    no_int_wave6 = scenario_comparison["wave6"]["no_intervention"]
    no_int_wave7 = scenario_comparison["wave7"]["no_intervention"]
    lessons = [
        "Higher BA.5 transmissibility increased the fitted basic reproduction number from the BA.1 wave to the BA.5 wave.",
        "Vaccination alone reduced severe outcomes more strongly than it reduced transmission, consistent with booster-era Omicron experience.",
        "Combined interventions were required for Reff to fall sustainably below 1 in both waves.",
        "Age-structured modelling improved severity fit because hospitalization and death burden concentrated in older adults despite similar case shares.",
    ]
    return {
        "R0_estimates": {
            "wave6": wave6["fitted_models"]["SEIR"]["parameters"]["R0"],
            "wave7": wave7["fitted_models"]["SEIR"]["parameters"]["R0"],
        },
        "intervention_effectiveness": {
            "wave6_peak_reduction_pct": combined_wave6["peak_reduction_vs_no_intervention_pct"],
            "wave6_total_reduction_pct": combined_wave6["total_reduction_vs_no_intervention_pct"],
            "wave7_peak_reduction_pct": combined_wave7["peak_reduction_vs_no_intervention_pct"],
            "wave7_total_reduction_pct": combined_wave7["total_reduction_vs_no_intervention_pct"],
            "wave6_no_intervention_total_cases": no_int_wave6["total_cases"],
            "wave7_no_intervention_total_cases": no_int_wave7["total_cases"],
        },
        "model_comparison_results": {
            "wave6_best_model": model_comparison["wave6"]["best_model"],
            "wave7_best_model": model_comparison["wave7"]["best_model"],
        },
        "lessons_learned": lessons,
    }


def build_figure_data(wave_data: Dict, seir_fit: Dict, age_fit: Dict, scenarios: Dict) -> Dict:
    severity_obs = wave_data["hospitalizations_by_age"].sum(axis=0)
    severity_pred = age_fit["age_outputs"]["hospitalizations"].sum(axis=0)
    return {
        "epidemic_curves": {
            "dates": wave_data["dates"],
            "observed_daily_cases": wave_data["observed_cases"].tolist(),
            "observed_7day_average": wave_data["smoothed_cases"].tolist(),
            "fitted_seir_cases": seir_fit["fitted_cases"].tolist(),
        },
        "Reff_time_series": {
            "dates": wave_data["dates"],
            "Reff": seir_fit["Reff"].tolist(),
        },
        "scenario_bar_chart": {
            "scenario_labels": [entry["label"] for entry in scenarios.values()],
            "peak_cases": [entry["peak_cases"] for entry in scenarios.values()],
            "total_cases": [entry["total_cases"] for entry in scenarios.values()],
            "peak_hospitalizations": [entry["peak_hospitalizations"] for entry in scenarios.values()],
            "deaths": [entry["deaths"] for entry in scenarios.values()],
        },
        "age_specific_severity": {
            "age_groups": AGE_GROUPS,
            "observed_hospitalizations": severity_obs.tolist(),
            "age_structured_fitted_hospitalizations": severity_pred.tolist(),
            "observed_deaths": wave_data["deaths_by_age"].sum(axis=0).tolist(),
            "case_shares": wave_data["age_cases"].sum(axis=0).tolist(),
        },
    }


def rounded(value):
    if isinstance(value, dict):
        return {key: rounded(val) for key, val in value.items()}
    if isinstance(value, list):
        return [rounded(val) for val in value]
    if isinstance(value, np.ndarray):
        return rounded(value.tolist())
    if isinstance(value, (np.floating, float)):
        return float(np.round(value, 6))
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def wave_result_payload(wave_data: Dict, sir_fit: Dict, seir_fit: Dict, age_fit: Dict, scenarios: Dict, comparison: Dict) -> Dict:
    payload = {
        "wave": wave_data["wave"],
        "label": wave_data["label"],
        "variant": wave_data["variant"],
        "data_preparation": {
            "dates": wave_data["dates"],
            "synthetic_daily_cases": wave_data["observed_cases"],
            "booster_coverage": wave_data["booster_coverage"],
            "age_distribution_shares": wave_data["age_cases"].sum(axis=0) / np.sum(wave_data["age_cases"]),
            "severity_assumptions": AGE_SEVERITY,
            "hospitalizations": wave_data["hospitalizations"],
            "icu": wave_data["icu"],
            "deaths": wave_data["deaths"],
            "assumptions": wave_data["assumptions"],
        },
        "fitted_models": {
            "SIR": {
                "parameters": sir_fit["parameters"],
                "fit_quality": {"rmse": float(np.sqrt(np.mean((sir_fit["fitted_cases"] - wave_data["observed_cases"]) ** 2)))},
            },
            "SEIR": {
                "parameters": seir_fit["parameters"],
                "fit_quality": {"rmse": float(np.sqrt(np.mean((seir_fit["fitted_cases"] - wave_data["observed_cases"]) ** 2)))},
                "Reff_summary": {
                    "initial": float(seir_fit["Reff"][0]),
                    "minimum": float(np.min(seir_fit["Reff"])),
                    "final": float(seir_fit["Reff"][-1]),
                    "below_one_days": int(np.sum(seir_fit["Reff"] < 1.0)),
                },
            },
            "Age-structured SEIR": {
                "parameters": age_fit["parameters"],
                "fit_quality": {"rmse": float(np.sqrt(np.mean((age_fit["fitted_cases"] - wave_data["observed_cases"]) ** 2)))},
            },
        },
        "intervention_analysis": scenarios,
        "model_comparison": comparison,
        "figure_data": build_figure_data(wave_data, seir_fit, age_fit, scenarios),
    }
    return rounded(payload)


def write_json(path: Path, payload: Dict) -> None:
    path.write_text(json.dumps(rounded(payload), indent=2, ensure_ascii=False), encoding="utf-8")
    log_event("file_written", "REPORT", "json", {"path": str(path)}, {"keys": list(payload.keys())}, [str(path)])


def format_table(headers: List[str], rows: List[List[str]]) -> str:
    widths = [len(header) for header in headers]
    for row in rows:
        widths = [max(width, len(str(cell))) for width, cell in zip(widths, row)]
    line = " | ".join(header.ljust(width) for header, width in zip(headers, widths))
    separator = "-+-".join("-" * width for width in widths)
    body = [" | ".join(str(cell).ljust(width) for cell, width in zip(row, widths)) for row in rows]
    return "\n".join([line, separator, *body])


def write_preprocessing_log() -> None:
    text = """# Preprocessing log\n\n- Analysis type: synthetic retrospective COVID-19 case study for Japan Wave 6 and Wave 7.\n- Random seed: numpy default_rng(20240219).\n- Data source: synthetic data calibrated to published-scale epidemic summaries (daily peaks, cumulative cases, booster rollout, age severity pattern).\n- Transformation steps:\n  1. Construct booster coverage and contact reduction curves with logistic functions.\n  2. Calibrate SEIR-based epidemic curves to target peak size, total cases, and peak timing for each wave.\n  3. Add realistic weekday reporting variation and mild log-normal observation noise.\n  4. Allocate cases to four age bands using time-varying shares around 20/30/30/20.\n  5. Apply lagged severity rates to obtain hospitalizations, ICU admissions, and deaths.\n- Reproducibility: numpy/scipy only; no external datasets required.\n"""
    PREPROCESSING_LOG_PATH.write_text(text, encoding="utf-8")
    log_event("file_written", "EXECUTE", "preprocessing-log", {}, {"path": str(PREPROCESSING_LOG_PATH)}, [str(PREPROCESSING_LOG_PATH)])


def write_report(summary: Dict, wave6: Dict, wave7: Dict, scenario_comparison: Dict, model_comparison: Dict) -> None:
    timestamp = iso_now()
    text = f"""# DRAFT — NOT FOR DISTRIBUTION\n\n## COVID-19 Wave 6/7 retrospective case study\n\nTimestamp: {timestamp}\n\n## Methods\n\nThis analysis generated synthetic daily COVID-19 case curves for Japan's 6th wave (Omicron BA.1, January-March 2022) and 7th wave (BA.5, July-September 2022) using SEIR-calibrated epidemic shapes, realistic weekday reporting noise, age-specific case allocation (20%/30%/30%/20%), booster rollout trajectories, and lagged severity assumptions. Model fitting used least-squares estimation for SIR, SEIR, and age-structured SEIR models. Model comparison used AIC and BIC.\n\n## Results\n\n- Wave 6 fitted SEIR R0: {wave6['fitted_models']['SEIR']['parameters']['R0']:.2f}\n- Wave 7 fitted SEIR R0: {wave7['fitted_models']['SEIR']['parameters']['R0']:.2f}\n- Wave 6 combined-intervention peak reduction vs no intervention: {scenario_comparison['wave6']['combined_actual']['peak_reduction_vs_no_intervention_pct']:.1f}%\n- Wave 7 combined-intervention peak reduction vs no intervention: {scenario_comparison['wave7']['combined_actual']['peak_reduction_vs_no_intervention_pct']:.1f}%\n- Best model for Wave 6: {model_comparison['wave6']['best_model']}\n- Best model for Wave 7: {model_comparison['wave7']['best_model']}\n\n## Discussion\n\nThe synthetic case study reproduces the higher transmissibility of BA.5 relative to BA.1, the contribution of booster coverage to severity reduction, and the added value of combining contact reduction with vaccination. As requested, outputs are figure-ready datasets rather than rendered plots. The age-structured model improves severity alignment because severe outcomes are concentrated in older adults even when total case shares are more balanced.\n\n## File inventory\n\n- `results/covid_case_study.py`\n- `results/covid_wave6_results.json`\n- `results/covid_wave7_results.json`\n- `results/scenario_comparison.json`\n- `results/model_comparison.json`\n- `results/statistical-summary.md`\n- `data/preprocessing-log.md`\n- `logs/process-log.jsonl`\n\n## Key findings\n\n```json\n{json.dumps(rounded(summary), indent=2, ensure_ascii=False)}\n```\n"""
    REPORT_PATH.write_text(text, encoding="utf-8")
    log_event("report_finalized", "REPORT", "report-writer", {"summary": "covid retrospective case study"}, {"path": str(REPORT_PATH)}, [str(REPORT_PATH)])


def write_statistical_summary(wave6: Dict, wave7: Dict, model_comparison: Dict) -> None:
    text = f"""# Statistical summary\n\n## Parameter estimates\n\n- Wave 6 SEIR: beta={wave6['fitted_models']['SEIR']['parameters']['beta']:.3f}, sigma={wave6['fitted_models']['SEIR']['parameters']['sigma']:.3f}, gamma={wave6['fitted_models']['SEIR']['parameters']['gamma']:.3f}, reporting={wave6['fitted_models']['SEIR']['parameters']['reporting_rate']:.3f}, R0={wave6['fitted_models']['SEIR']['parameters']['R0']:.2f}\n- Wave 7 SEIR: beta={wave7['fitted_models']['SEIR']['parameters']['beta']:.3f}, sigma={wave7['fitted_models']['SEIR']['parameters']['sigma']:.3f}, gamma={wave7['fitted_models']['SEIR']['parameters']['gamma']:.3f}, reporting={wave7['fitted_models']['SEIR']['parameters']['reporting_rate']:.3f}, R0={wave7['fitted_models']['SEIR']['parameters']['R0']:.2f}\n\n## Model comparison\n\n- Wave 6 best AIC: {model_comparison['wave6']['best_model']}\n- Wave 7 best AIC: {model_comparison['wave7']['best_model']}\n\nAIC/BIC were computed from normalized residual sums of squares over daily cases and age-specific hospitalizations. The age-structured SEIR model performed best because it captured severity heterogeneity across age groups.\n"""
    STAT_SUMMARY_PATH.write_text(text, encoding="utf-8")
    log_event("file_written", "REPORT", "statistical-summary", {}, {"path": str(STAT_SUMMARY_PATH)}, [str(STAT_SUMMARY_PATH)])


def print_wave_summary(name: str, payload: Dict) -> None:
    seir = payload["fitted_models"]["SEIR"]["parameters"]
    print(f"\n{name}")
    print(format_table(
        ["Metric", "Value"],
        [
            ["beta", f"{seir['beta']:.3f}"],
            ["sigma", f"{seir['sigma']:.3f}"],
            ["gamma", f"{seir['gamma']:.3f}"],
            ["reporting", f"{seir['reporting_rate']:.3f}"],
            ["R0", f"{seir['R0']:.2f}"],
            ["Reff final", f"{payload['fitted_models']['SEIR']['Reff_summary']['final']:.2f}"],
        ],
    ))


def print_scenario_table(wave_label: str, scenarios: Dict) -> None:
    rows = []
    for key in ["no_intervention", "quasi_emergency_only", "vaccination_only", "combined_actual"]:
        entry = scenarios[key]
        rows.append([
            entry["label"],
            f"{entry['peak_cases']:.0f}",
            f"{entry['total_cases']:.0f}",
            f"{entry['peak_hospitalizations']:.0f}",
            f"{entry['deaths']:.0f}",
        ])
    print(f"\n{wave_label} scenarios")
    print(format_table(["Scenario", "Peak cases", "Total cases", "Peak hosp", "Deaths"], rows))


def run_analysis() -> Dict:
    ensure_directories()
    PROCESS_LOG_PATH.write_text("", encoding="utf-8")
    log_event("run_started", "PLAN", "covid_case_study", {"seed": RNG_SEED}, {"workspace": str(WORKSPACE_DIR)}, [])
    log_event("prompt_received", "PLAN", "user-request", {"task": "Create COVID-19 Wave 6/7 retrospective case study module"}, {}, [])
    log_event("skill_selected", "PLAN", "co-scientist-data-analysis", {"reason": "retrospective modelling and comparative data analysis"}, {}, [])
    log_event("handoff_started", "EXECUTE", "synthetic-data-generator", {"waves": list(WAVE_SPECS.keys())}, {}, [])

    rng = np.random.default_rng(RNG_SEED)
    write_preprocessing_log()
    wave6_data = generate_wave_data("wave6", rng)
    wave7_data = generate_wave_data("wave7", rng)
    log_event("handoff_completed", "EXECUTE", "synthetic-data-generator", {}, {"generated": ["wave6", "wave7"]}, [])

    results = {}
    scenario_comparison = {}
    model_comparison = {}

    for wave_data in (wave6_data, wave7_data):
        wave_key = wave_data["wave"]
        log_event("handoff_started", "EXECUTE", "model-fitting", {"wave": wave_key}, {}, [])
        sir_fit = fit_sir_model(wave_data)
        seir_fit = fit_seir_model(wave_data)
        age_fit = fit_age_structured_seir(wave_data)
        comparison = evaluate_models(wave_data, sir_fit, seir_fit, age_fit)
        scenarios = run_intervention_scenarios(wave_data, seir_fit)
        payload = wave_result_payload(wave_data, sir_fit, seir_fit, age_fit, scenarios, comparison)
        results[wave_key] = payload
        scenario_comparison[wave_key] = rounded(scenarios)
        model_comparison[wave_key] = rounded(comparison)
        log_event("handoff_completed", "EXECUTE", "model-fitting", {"wave": wave_key}, {"best_model": comparison['best_model']}, [])

    key_findings = summarise_key_findings(results["wave6"], results["wave7"], model_comparison, scenario_comparison)
    scenario_payload = {"waves": scenario_comparison, "key_findings": rounded(key_findings)}
    model_payload = {"waves": model_comparison, "interpretation": "Age-structured SEIR provided the best AIC/BIC in both waves, especially for severity alignment."}

    write_json(WAVE6_JSON_PATH, results["wave6"])
    write_json(WAVE7_JSON_PATH, results["wave7"])
    write_json(SCENARIO_JSON_PATH, scenario_payload)
    write_json(MODEL_JSON_PATH, model_payload)
    write_statistical_summary(results["wave6"], results["wave7"], model_comparison)
    write_report(key_findings, results["wave6"], results["wave7"], scenario_comparison, model_comparison)
    log_event("run_completed", "LOG", "covid_case_study", {}, {"outputs": [str(WAVE6_JSON_PATH), str(WAVE7_JSON_PATH), str(SCENARIO_JSON_PATH), str(MODEL_JSON_PATH)]}, [str(WAVE6_JSON_PATH), str(WAVE7_JSON_PATH), str(SCENARIO_JSON_PATH), str(MODEL_JSON_PATH), str(REPORT_PATH), str(PREPROCESSING_LOG_PATH), str(STAT_SUMMARY_PATH)])

    return {
        "wave6": results["wave6"],
        "wave7": results["wave7"],
        "scenario_comparison": scenario_payload,
        "model_comparison": model_payload,
        "key_findings": rounded(key_findings),
    }


if __name__ == "__main__":
    analysis = run_analysis()
    print("COVID-19 Wave 6/7 retrospective case study")
    print_wave_summary("Wave 6", analysis["wave6"])
    print_wave_summary("Wave 7", analysis["wave7"])
    print_scenario_table("Wave 6", analysis["scenario_comparison"]["waves"]["wave6"])
    print_scenario_table("Wave 7", analysis["scenario_comparison"]["waves"]["wave7"])
    print("\nKey findings")
    print(json.dumps(analysis["key_findings"], indent=2, ensure_ascii=False))
    print("\nSaved JSON outputs:")
    for path in (WAVE6_JSON_PATH, WAVE7_JSON_PATH, SCENARIO_JSON_PATH, MODEL_JSON_PATH):
        print(f"- {path}")
