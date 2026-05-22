"""
Epidemiology Data Integration Module.
Integrates case counts, mobility data, and wastewater surveillance signals.
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")

try:
    from sklearn.ensemble import IsolationForest
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


def _generate_epidemic_curve(n_days: int = 90, seed: int = 42,
                              r0: float = 1.3, peak_day: int = 45) -> np.ndarray:
    """Generate a synthetic epidemic incidence curve (log-normal + noise)."""
    rng = np.random.default_rng(seed)
    t = np.arange(n_days)
    # Gaussian-like epidemic curve
    base = 500 * np.exp(-0.5 * ((t - peak_day) / 15) ** 2)
    noise = rng.normal(0, base * 0.1 + 5)
    curve = np.maximum(base + noise, 0).astype(int)
    return curve


class CaseDataIntegrator:
    """Fetches and processes case surveillance data."""

    def fetch_case_data(self, countries: List[str] = None,
                        n_days: int = 90, seed: int = 42) -> pd.DataFrame:
        """Simulate daily case counts for multiple countries."""
        countries = countries or ["Japan", "USA", "Germany", "Brazil", "India"]
        rng = np.random.default_rng(seed)
        now = datetime.now()
        records = []
        for country in countries:
            peak_day = int(rng.integers(30, 70))
            base_pop = {"Japan": 125e6, "USA": 330e6, "Germany": 83e6,
                        "Brazil": 215e6, "India": 1400e6}.get(country, 50e6)
            scale = base_pop / 1e7
            curve = _generate_epidemic_curve(n_days, seed=seed + hash(country) % 1000,
                                             peak_day=peak_day)
            curve = (curve * scale).astype(int)
            for d in range(n_days):
                records.append({
                    "date": (now - timedelta(days=n_days - d)).strftime("%Y-%m-%d"),
                    "country": country,
                    "new_cases": int(curve[d]),
                    "cumulative_cases": int(curve[:d + 1].sum()),
                    "new_deaths": int(max(0, curve[d] * rng.normal(0.01, 0.002))),
                    "hospitalizations": int(max(0, curve[d] * rng.normal(0.05, 0.005))),
                })
        return pd.DataFrame(records)

    def compute_growth_rate(self, cases_df: pd.DataFrame,
                             window: int = 7) -> pd.DataFrame:
        """Compute rolling growth rate and doubling time."""
        df = cases_df.copy()
        df["date"] = pd.to_datetime(df["date"])
        results = []
        for country, grp in df.groupby("country"):
            grp = grp.sort_values("date").copy()
            grp["rolling_cases"] = grp["new_cases"].rolling(window, min_periods=3).mean()
            grp["growth_rate"] = grp["rolling_cases"].pct_change(periods=7).fillna(0)
            grp["doubling_time"] = np.where(
                grp["growth_rate"] > 0,
                np.log(2) / np.log(1 + grp["growth_rate"]),
                np.inf
            )
            results.append(grp)
        return pd.concat(results, ignore_index=True)

    def detect_anomalies(self, cases_df: pd.DataFrame,
                          contamination: float = 0.05) -> pd.DataFrame:
        """Detect anomalous case counts using Isolation Forest."""
        df = cases_df.copy()
        df["date"] = pd.to_datetime(df["date"])
        results = []
        for country, grp in df.groupby("country"):
            grp = grp.sort_values("date").copy()
            if len(grp) < 10:
                grp["anomaly"] = False
                results.append(grp)
                continue
            features = grp[["new_cases", "new_deaths", "hospitalizations"]].fillna(0)
            if HAS_SKLEARN:
                iso = IsolationForest(contamination=contamination, random_state=42)
                grp["anomaly"] = iso.fit_predict(features) == -1
            else:
                z = (grp["new_cases"] - grp["new_cases"].mean()) / (grp["new_cases"].std() + 1e-9)
                grp["anomaly"] = z.abs() > 3
            results.append(grp)
        return pd.concat(results, ignore_index=True)


class MobilityDataAnalyzer:
    """Processes mobility data (Google/Apple-like mobility reports)."""

    def fetch_mobility_data(self, countries: List[str] = None,
                             n_days: int = 90, seed: int = 77) -> pd.DataFrame:
        """Simulate mobility index data."""
        countries = countries or ["Japan", "USA", "Germany", "Brazil", "India"]
        rng = np.random.default_rng(seed)
        now = datetime.now()
        records = []
        for country in countries:
            baseline = 100.0
            trend = rng.uniform(-30, 10, n_days)
            cumulative = baseline + np.cumsum(trend * 0.05)
            cumulative = np.clip(cumulative, 40, 120)
            for d in range(n_days):
                records.append({
                    "date": (now - timedelta(days=n_days - d)).strftime("%Y-%m-%d"),
                    "country": country,
                    "retail_recreation": round(float(cumulative[d]) + rng.normal(0, 3), 1),
                    "transit_stations": round(float(cumulative[d]) * 0.9 + rng.normal(0, 4), 1),
                    "workplaces": round(float(cumulative[d]) * 0.95 + rng.normal(0, 2), 1),
                    "residential": round(100 + (100 - cumulative[d]) * 0.3 + rng.normal(0, 1), 1),
                    "mobility_index": round(float(cumulative[d]), 2),
                })
        return pd.DataFrame(records)

    def process_mobility_data(self, mobility_df: pd.DataFrame) -> pd.DataFrame:
        """Compute mobility change relative to baseline and trend."""
        df = mobility_df.copy()
        df["date"] = pd.to_datetime(df["date"])
        df["mobility_change_pct"] = df["mobility_index"] - 100
        df["mobility_7d_avg"] = df.groupby("country")["mobility_index"].transform(
            lambda x: x.rolling(7, min_periods=1).mean()
        )
        return df


class WastewaterSurveillance:
    """Processes wastewater epidemiology signals."""

    def fetch_wastewater_data(self, sites: List[str] = None,
                               n_days: int = 90, seed: int = 55) -> pd.DataFrame:
        """Simulate wastewater viral load measurements."""
        sites = sites or ["Tokyo", "New_York", "Berlin", "Sao_Paulo", "Mumbai"]
        rng = np.random.default_rng(seed)
        now = datetime.now()
        records = []
        for site in sites:
            base_curve = _generate_epidemic_curve(n_days, seed=seed + hash(site) % 500,
                                                   peak_day=int(rng.integers(35, 65)))
            for d in range(n_days):
                viral_load = max(0, base_curve[d] * rng.lognormal(0, 0.3) * 1000)
                records.append({
                    "date": (now - timedelta(days=n_days - d)).strftime("%Y-%m-%d"),
                    "site": site,
                    "viral_load_gc_per_L": round(float(viral_load), 1),
                    "normalized_load": round(float(viral_load / 1e6), 4),
                    "sample_quality": rng.choice(["good", "acceptable", "poor"],
                                                  p=[0.7, 0.2, 0.1]),
                })
        return pd.DataFrame(records)

    def analyze_wastewater_signal(self, ww_df: pd.DataFrame,
                                   cases_df: pd.DataFrame,
                                   lag_days: int = 7) -> Dict:
        """Cross-correlate wastewater signal with reported cases."""
        results = {}
        ww_df = ww_df.copy()
        ww_df["date"] = pd.to_datetime(ww_df["date"])
        cases_df = cases_df.copy()
        cases_df["date"] = pd.to_datetime(cases_df["date"])

        # Aggregate wastewater by date
        ww_agg = ww_df.groupby("date")["viral_load_gc_per_L"].mean().reset_index()
        ww_agg.columns = ["date", "ww_signal"]

        # Aggregate cases by date
        cases_agg = cases_df.groupby("date")["new_cases"].sum().reset_index()

        merged = pd.merge(ww_agg, cases_agg, on="date").sort_values("date")
        if len(merged) > lag_days:
            ww_sig = merged["ww_signal"].values
            cases_sig = merged["new_cases"].values
            # Cross-correlation
            from numpy.fft import fft, ifft
            n = len(ww_sig)
            ww_norm = (ww_sig - ww_sig.mean()) / (ww_sig.std() + 1e-9)
            c_norm = (cases_sig - cases_sig.mean()) / (cases_sig.std() + 1e-9)
            xcorr = np.correlate(ww_norm, c_norm, mode="full")
            lags = np.arange(-n + 1, n)
            best_lag = int(lags[np.argmax(xcorr)])
            max_corr = float(xcorr.max() / n)
            results = {
                "best_lag_days": best_lag,
                "max_cross_correlation": round(max_corr, 4),
                "n_data_points": len(merged),
                "ww_leads_cases": best_lag < 0,
            }
        return results

    def compute_ww_trend(self, ww_df: pd.DataFrame) -> pd.DataFrame:
        """Compute rolling trend for wastewater signal."""
        df = ww_df.copy()
        df["date"] = pd.to_datetime(df["date"])
        df = df[df["sample_quality"] != "poor"]
        df["ww_7d_avg"] = df.groupby("site")["viral_load_gc_per_L"].transform(
            lambda x: x.rolling(7, min_periods=1).mean()
        )
        df["ww_growth_rate"] = df.groupby("site")["ww_7d_avg"].transform(
            lambda x: x.pct_change(7).fillna(0)
        )
        return df


def integrate_signals(cases_df: pd.DataFrame, mobility_df: pd.DataFrame,
                       wastewater_df: pd.DataFrame) -> Dict:
    """Compute unified epidemiological signal score."""
    # Latest growth rate per country
    cases_df["date"] = pd.to_datetime(cases_df["date"])
    latest_cases = cases_df.sort_values("date").groupby("country").last()
    avg_growth = float(latest_cases.get("growth_rate", pd.Series([0])).mean())

    # Mobility deviation from baseline
    mobility_df["date"] = pd.to_datetime(mobility_df["date"])
    latest_mob = mobility_df.sort_values("date").groupby("country").last()
    avg_mobility = float(latest_mob["mobility_index"].mean())
    mobility_signal = max(0, (avg_mobility - 80) / 40)  # 0-1 scaled

    # Wastewater trend
    wastewater_df["date"] = pd.to_datetime(wastewater_df["date"])
    ww_trend = wastewater_df.sort_values("date").groupby("site").last()
    avg_ww_growth = float(ww_trend.get("ww_growth_rate", pd.Series([0])).mean())

    # Composite score (0-100)
    growth_score = min(100, max(0, avg_growth * 200))
    ww_score = min(100, max(0, avg_ww_growth * 100))
    mobility_score = mobility_signal * 100

    composite = 0.4 * growth_score + 0.3 * ww_score + 0.3 * mobility_score

    return {
        "composite_epi_score": round(float(composite), 2),
        "avg_case_growth_rate": round(avg_growth, 4),
        "avg_mobility_index": round(avg_mobility, 2),
        "avg_ww_growth_rate": round(avg_ww_growth, 4),
        "growth_score": round(growth_score, 2),
        "ww_score": round(ww_score, 2),
        "mobility_score": round(mobility_score, 2),
    }


def run_epidemiology_pipeline(config: Optional[Dict] = None) -> Dict:
    """Run the full epidemiology integration pipeline."""
    config = config or {}
    integrator = CaseDataIntegrator()
    mobility_analyzer = MobilityDataAnalyzer()
    ww_analyzer = WastewaterSurveillance()

    cases_df = integrator.fetch_case_data()
    cases_df = integrator.compute_growth_rate(cases_df)
    cases_df = integrator.detect_anomalies(cases_df)

    mobility_df = mobility_analyzer.fetch_mobility_data()
    mobility_df = mobility_analyzer.process_mobility_data(mobility_df)

    ww_df = ww_analyzer.fetch_wastewater_data()
    ww_df = ww_analyzer.compute_ww_trend(ww_df)
    ww_correlation = ww_analyzer.analyze_wastewater_signal(ww_df, cases_df)

    signals = integrate_signals(cases_df, mobility_df, ww_df)

    n_anomalies = int(cases_df["anomaly"].sum()) if "anomaly" in cases_df.columns else 0

    return {
        "cases_df": cases_df,
        "mobility_df": mobility_df,
        "wastewater_df": ww_df,
        "ww_correlation": ww_correlation,
        "integrated_signals": signals,
        "n_anomalies_detected": n_anomalies,
    }


if __name__ == "__main__":
    results = run_epidemiology_pipeline()
    print(f"Anomalies detected: {results['n_anomalies_detected']}")
    print(f"Integrated signals: {results['integrated_signals']}")
    print(f"WW cross-correlation: {results['ww_correlation']}")
