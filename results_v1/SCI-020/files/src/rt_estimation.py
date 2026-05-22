"""
Rt Estimation Module — Improved EpiEstim implementation.
Bayesian real-time estimation of effective reproduction number Rt.
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")

try:
    from statsmodels.tsa.arima.model import ARIMA
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


def gamma_serial_interval(mean: float = 4.7, sd: float = 2.9,
                           max_si: int = 30) -> np.ndarray:
    """
    Compute discretized gamma serial interval distribution.
    Returns array of length max_si+1 (index = days).
    """
    if sd <= 0 or mean <= 0:
        raise ValueError("Serial interval mean and sd must be positive.")
    shape = (mean / sd) ** 2
    scale = sd ** 2 / mean
    dist = stats.gamma(a=shape, scale=scale)
    # Discretize: P(SI = k) = F(k+0.5) - F(k-0.5)
    w = np.array([dist.cdf(k + 0.5) - dist.cdf(max(0, k - 0.5))
                  for k in range(max_si + 1)])
    w[0] = 0  # No transmission on day 0
    w /= w.sum()
    return w


def _compute_lambda(incidence: np.ndarray, w: np.ndarray, t: int) -> float:
    """Compute the force of infection (lambda) at time t."""
    max_tau = min(t, len(w) - 1)
    lam = sum(incidence[t - s] * w[s] for s in range(1, max_tau + 1))
    return max(lam, 1e-9)


class ImprovedEpiEstim:
    """
    Bayesian Rt estimator based on Cori et al. (2013) EpiEstim.
    Improved with: sliding-window credible intervals, nowcasting adjustment,
    and ARIMA-based forecasting.
    """

    def __init__(self, serial_interval_mean: float = 4.7,
                 serial_interval_sd: float = 2.9,
                 window_size: int = 7, tau: int = 7):
        self.si_mean = serial_interval_mean
        self.si_sd = serial_interval_sd
        self.window = window_size
        self.tau = tau
        self.w = gamma_serial_interval(serial_interval_mean, serial_interval_sd)

    def estimate_rt(self, incidence: np.ndarray,
                    prior_mean: float = 5.0,
                    prior_cv: float = 0.5) -> pd.DataFrame:
        """
        Estimate Rt using gamma-Poisson Bayesian conjugate model.
        Prior: Gamma(a0, b0) with mean=prior_mean, cv=prior_cv
        Posterior: Gamma(a0 + sum(I_t), b0 + sum(lambda_t))

        Returns DataFrame with columns: t, Rt_mean, Rt_median, CI_lower, CI_upper, lambda_t
        """
        n = len(incidence)
        # Prior hyperparameters
        a0 = 1.0 / (prior_cv ** 2)
        b0 = a0 / prior_mean

        results = []
        for t in range(self.tau, n):
            # Sliding window [t-tau+1, t]
            t_start = max(0, t - self.tau + 1)
            sum_I = float(np.sum(incidence[t_start:t + 1]))
            sum_lambda = sum(_compute_lambda(incidence, self.w, s)
                             for s in range(t_start, t + 1))

            # Posterior parameters
            a_post = a0 + sum_I
            b_post = b0 + sum_lambda

            # Posterior statistics
            rt_mean = a_post / b_post
            rt_median = stats.gamma.ppf(0.5, a=a_post, scale=1 / b_post)
            rt_lower = stats.gamma.ppf(0.025, a=a_post, scale=1 / b_post)
            rt_upper = stats.gamma.ppf(0.975, a=a_post, scale=1 / b_post)
            lambda_t = _compute_lambda(incidence, self.w, t)

            results.append({
                "t": t,
                "incidence": int(incidence[t]),
                "Rt_mean": round(float(rt_mean), 4),
                "Rt_median": round(float(rt_median), 4),
                "CI_lower_95": round(float(rt_lower), 4),
                "CI_upper_95": round(float(rt_upper), 4),
                "lambda_t": round(float(lambda_t), 2),
            })

        return pd.DataFrame(results)

    def detect_threshold_crossings(self, rt_df: pd.DataFrame,
                                    threshold: float = 1.0) -> pd.DataFrame:
        """Detect when Rt crosses threshold (e.g., above 1.0 = growing epidemic)."""
        df = rt_df.copy()
        df["above_threshold"] = df["Rt_mean"] > threshold
        df["crossing_up"] = (~df["above_threshold"].shift(1).fillna(False)) & df["above_threshold"]
        df["crossing_down"] = df["above_threshold"].shift(1).fillna(False) & (~df["above_threshold"])
        crossings = df[df["crossing_up"] | df["crossing_down"]][["t", "Rt_mean", "crossing_up", "crossing_down"]]
        return crossings

    def forecast_rt(self, rt_df: pd.DataFrame, days_ahead: int = 14) -> pd.DataFrame:
        """ARIMA-based Rt forecast with confidence intervals."""
        series = rt_df["Rt_mean"].values

        if HAS_STATSMODELS and len(series) >= 20:
            try:
                model = ARIMA(series, order=(2, 1, 2))
                fit = model.fit()
                forecast = fit.get_forecast(steps=days_ahead)
                fcast_mean = forecast.predicted_mean
                ci = forecast.conf_int(alpha=0.05)
                records = []
                for i in range(days_ahead):
                    records.append({
                        "t_future": int(rt_df["t"].max()) + i + 1,
                        "Rt_forecast": round(float(fcast_mean[i]), 4),
                        "CI_lower": round(float(ci.iloc[i, 0]), 4),
                        "CI_upper": round(float(ci.iloc[i, 1]), 4),
                    })
                return pd.DataFrame(records)
            except Exception:
                pass

        # Fallback: simple linear extrapolation
        last_n = series[-7:] if len(series) >= 7 else series
        slope = float(np.polyfit(range(len(last_n)), last_n, 1)[0])
        last_val = float(series[-1])
        records = []
        for i in range(days_ahead):
            predicted = max(0.1, last_val + slope * (i + 1))
            records.append({
                "t_future": int(rt_df["t"].max()) + i + 1,
                "Rt_forecast": round(predicted, 4),
                "CI_lower": round(max(0.1, predicted - 0.3), 4),
                "CI_upper": round(predicted + 0.3, 4),
            })
        return pd.DataFrame(records)

    def adjust_for_reporting_delay(self, incidence: np.ndarray,
                                    mean_delay: float = 5.0,
                                    sd_delay: float = 2.0) -> np.ndarray:
        """
        Nowcasting: adjust incidence for reporting delay using Richardson-Lucy deconvolution.
        """
        n = len(incidence)
        delay_dist = gamma_serial_interval(mean_delay, sd_delay, max_si=20)[:n]
        delay_dist = delay_dist / delay_dist.sum()

        # Simple deconvolution via iterative Richardson-Lucy
        estimate = incidence.copy().astype(float)
        for _ in range(10):
            conv = np.convolve(estimate, delay_dist)[:n]
            conv = np.maximum(conv, 1e-9)
            ratio = incidence.astype(float) / conv
            correction = np.convolve(ratio[::-1], delay_dist)[:n][::-1]
            estimate = estimate * np.maximum(correction, 0.01)

        return np.maximum(estimate, 0).astype(int)


class RtEstimator:
    """High-level interface for Rt estimation."""

    def __init__(self, config: Optional[Dict] = None):
        config = config or {}
        self.model = ImprovedEpiEstim(
            serial_interval_mean=config.get("serial_interval_mean", 4.7),
            serial_interval_sd=config.get("serial_interval_sd", 2.9),
            window_size=config.get("window_size", 7),
            tau=config.get("tau", 7),
        )

    def run(self, incidence: np.ndarray, country: str = "Global") -> Dict:
        """Full Rt estimation pipeline for a single country."""
        # Adjust for reporting delay
        adjusted = self.model.adjust_for_reporting_delay(incidence)

        # Estimate Rt
        rt_df = self.model.estimate_rt(adjusted)

        # Forecast
        forecast_df = self.model.forecast_rt(rt_df)

        # Threshold crossings
        crossings = self.model.detect_threshold_crossings(rt_df)

        # Summary stats
        latest_rt = float(rt_df["Rt_mean"].iloc[-1]) if len(rt_df) else 1.0
        max_rt = float(rt_df["Rt_mean"].max()) if len(rt_df) else 1.0
        min_rt = float(rt_df["Rt_mean"].min()) if len(rt_df) else 1.0
        ci_lower = float(rt_df["CI_lower_95"].iloc[-1]) if len(rt_df) else 0.8
        ci_upper = float(rt_df["CI_upper_95"].iloc[-1]) if len(rt_df) else 1.2

        return {
            "country": country,
            "rt_df": rt_df,
            "forecast_df": forecast_df,
            "threshold_crossings": crossings,
            "latest_rt": round(latest_rt, 3),
            "max_rt": round(max_rt, 3),
            "min_rt": round(min_rt, 3),
            "ci_lower": round(ci_lower, 3),
            "ci_upper": round(ci_upper, 3),
            "n_upward_crossings": int((crossings["crossing_up"]).sum()) if len(crossings) else 0,
            "epidemic_growing": latest_rt > 1.0,
        }


def run_rt_estimation(cases_df: Optional[pd.DataFrame] = None) -> Dict:
    """Run Rt estimation for all countries in the case data."""
    estimator = RtEstimator()

    if cases_df is not None:
        cases_df = cases_df.copy()
        cases_df["date"] = pd.to_datetime(cases_df["date"])
        countries = cases_df["country"].unique()
    else:
        # Generate synthetic incidence
        from epidemiology_integration import _generate_epidemic_curve
        countries = ["Global"]

    all_results = {}
    rt_records = []

    for country in countries:
        if cases_df is not None:
            grp = cases_df[cases_df["country"] == country].sort_values("date")
            incidence = grp["new_cases"].values
        else:
            incidence = _generate_epidemic_curve(90, seed=42)

        if len(incidence) < 15:
            continue

        res = estimator.run(incidence, country=country)
        all_results[country] = res

        # Flatten Rt estimates for CSV
        for _, row in res["rt_df"].iterrows():
            rt_records.append({
                "country": country,
                "t": row["t"],
                "Rt_mean": row["Rt_mean"],
                "CI_lower_95": row["CI_lower_95"],
                "CI_upper_95": row["CI_upper_95"],
            })

    rt_summary_df = pd.DataFrame(rt_records)
    latest_rts = {c: all_results[c]["latest_rt"] for c in all_results}
    n_growing = sum(1 for c in all_results if all_results[c]["epidemic_growing"])

    return {
        "country_results": all_results,
        "rt_summary_df": rt_summary_df,
        "latest_rts": latest_rts,
        "n_countries_growing": n_growing,
        "global_rt_mean": round(float(np.mean(list(latest_rts.values()))), 3) if latest_rts else 1.0,
    }


if __name__ == "__main__":
    results = run_rt_estimation()
    print(f"Countries growing (Rt > 1): {results['n_countries_growing']}")
    print(f"Global mean Rt: {results['global_rt_mean']}")
    print("Latest Rt by country:", results["latest_rts"])
