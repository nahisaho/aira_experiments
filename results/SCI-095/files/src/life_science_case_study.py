from __future__ import annotations

import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
np.random.seed(SEED)
random.seed(SEED)

ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"
REPORT_PATH = ROOT / "report.md"
RESULTS_JSON = RESULTS_DIR / "life_science_results.json"
STATS_SUMMARY = RESULTS_DIR / "statistical-summary.md"
PREPROCESSING_LOG = DATA_DIR / "preprocessing-log.md"

PALETTE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "grey": "#4D4D4D",
}


def ensure_dirs() -> None:
    for directory in (FIGURES_DIR, RESULTS_DIR, DATA_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def bootstrap_ci(frame: pd.DataFrame, func, n_boot: int = 2000, seed: int = SEED) -> tuple[float, tuple[float, float]]:
    rng = np.random.default_rng(seed)
    estimate = float(func(frame))
    samples = np.empty(n_boot)
    for i in range(n_boot):
        sample = frame.iloc[rng.integers(0, len(frame), len(frame))]
        samples[i] = float(func(sample))
    lower, upper = np.percentile(samples, [2.5, 97.5])
    return estimate, (float(lower), float(upper))


def metric_summary(frame: pd.DataFrame, case_name: str) -> dict:
    roi_func = lambda x: x["research_value_usd"].sum() / x["sharing_cost_usd"].sum()
    acceleration_func = lambda x: x["closed_time_months"].sum() / x["open_time_months"].sum()
    elasticity_func = lambda x: stats.linregress(
        np.log1p(x["reuse_count"]), np.log1p(x["network_outputs"])
    ).slope

    roi, roi_ci = bootstrap_ci(frame, roi_func)
    acceleration, acceleration_ci = bootstrap_ci(frame, acceleration_func)
    elasticity, elasticity_ci = bootstrap_ci(frame, elasticity_func)
    spearman_rho, spearman_p = stats.spearmanr(frame["reuse_count"], frame["network_outputs"])

    return {
        "case_name": case_name,
        "datasets": int(len(frame)),
        "total_reuse_count": int(frame["reuse_count"].sum()),
        "total_network_outputs": int(frame["network_outputs"].sum()),
        "total_sharing_cost_usd": round(float(frame["sharing_cost_usd"].sum()), 2),
        "total_research_value_usd": round(float(frame["research_value_usd"].sum()), 2),
        "roi": round(float(roi), 3),
        "roi_ci_95": [round(float(roi_ci[0]), 3), round(float(roi_ci[1]), 3)],
        "acceleration_factor": round(float(acceleration), 3),
        "acceleration_factor_ci_95": [
            round(float(acceleration_ci[0]), 3),
            round(float(acceleration_ci[1]), 3),
        ],
        "network_elasticity": round(float(elasticity), 3),
        "network_elasticity_ci_95": [
            round(float(elasticity_ci[0]), 3),
            round(float(elasticity_ci[1]), 3),
        ],
        "network_spearman_rho": round(float(spearman_rho), 3),
        "network_spearman_p": float(spearman_p),
    }


def simulate_genomic_case() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    years = np.arange(2015, 2025)
    deposits_year = np.array([18, 22, 26, 31, 36, 48, 52, 56, 58, 60])
    annual_reuse = np.array([24, 32, 43, 58, 79, 142, 188, 176, 215, 264])
    records = []
    dataset_id = 1

    for year, n_deposits in zip(years, deposits_year):
        for _ in range(int(n_deposits)):
            covid_relevance = np.random.binomial(1, 0.62 if year in (2020, 2021) else 0.14)
            maturity = 1 + 0.28 * (2024 - year)
            growth = np.exp(0.13 * (year - 2015))
            covid_boost = 2.4 if covid_relevance else 1.0
            reuse_lambda = 1.8 * maturity * growth * covid_boost
            reuse_count = int(np.random.poisson(reuse_lambda) + 1)
            secondary_publications = int(np.random.poisson(0.32 * reuse_count + 0.5))
            tool_development_count = int(np.random.poisson(0.11 * reuse_count + 0.15 * (year >= 2018)))
            sharing_cost_usd = float(np.random.lognormal(np.log(42000), 0.24))
            research_value_usd = float(
                reuse_count * 82000
                + secondary_publications * 145000
                + tool_development_count * 225000
                + np.random.normal(0, 30000)
            )
            open_time_months = float(np.clip(np.random.normal(13.5 - 0.05 * reuse_count, 2.0), 4, None))
            closed_time_months = float(
                open_time_months * (1.65 + 0.07 * secondary_publications + 0.09 * tool_development_count)
            )
            records.append(
                {
                    "dataset_id": f"GEO-SRA-{dataset_id:04d}",
                    "deposit_year": int(year),
                    "reuse_count": reuse_count,
                    "secondary_publications": secondary_publications,
                    "tool_development_count": tool_development_count,
                    "sharing_cost_usd": round(sharing_cost_usd, 2),
                    "research_value_usd": round(max(research_value_usd, 60000), 2),
                    "open_time_months": round(open_time_months, 2),
                    "closed_time_months": round(closed_time_months, 2),
                    "network_outputs": secondary_publications + tool_development_count,
                    "covid_relevance": int(covid_relevance),
                }
            )
            dataset_id += 1

    frame = pd.DataFrame.from_records(records)
    yearly = pd.DataFrame(
        {
            "year": years,
            "deposits": deposits_year,
            "reuse_events": annual_reuse,
        }
    )
    summary = metric_summary(frame, "Genomic Data Sharing")
    summary["yearly_deposits"] = {str(y): int(v) for y, v in zip(years, deposits_year)}
    summary["yearly_reuse_events"] = {str(y): int(v) for y, v in zip(years, annual_reuse)}
    return frame, yearly, summary


def protein_method_probabilities(year: int) -> list[float]:
    if year <= 2017:
        return [0.58, 0.12, 0.30, 0.00]
    if year == 2018:
        return [0.50, 0.25, 0.25, 0.00]
    if year == 2019:
        return [0.44, 0.36, 0.17, 0.03]
    if year == 2020:
        return [0.40, 0.37, 0.13, 0.10]
    if year == 2021:
        return [0.31, 0.35, 0.09, 0.25]
    if year == 2022:
        return [0.26, 0.31, 0.08, 0.35]
    if year == 2023:
        return [0.23, 0.28, 0.08, 0.41]
    return [0.20, 0.25, 0.07, 0.48]


def simulate_protein_case() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    years = np.arange(2015, 2025)
    deposits_year = np.array([20, 22, 24, 26, 28, 30, 32, 34, 38, 46])
    methods = ["X-ray", "Cryo-EM", "NMR", "AlphaFold"]
    method_bonus = {"X-ray": 1.0, "Cryo-EM": 1.35, "NMR": 0.82, "AlphaFold": 1.55}
    share_cost = {"X-ray": 62000, "Cryo-EM": 86000, "NMR": 52000, "AlphaFold": 24000}
    closed_multiplier = {"X-ray": 1.40, "Cryo-EM": 1.72, "NMR": 1.18, "AlphaFold": 2.18}
    open_months = {"X-ray": 18.0, "Cryo-EM": 12.5, "NMR": 20.5, "AlphaFold": 8.8}
    records = []
    structure_id = 1

    for year, n_deposits in zip(years, deposits_year):
        year_factor = 1 + 0.08 * (year - 2015)
        for _ in range(int(n_deposits)):
            method = np.random.choice(methods, p=protein_method_probabilities(int(year)))
            reuse_lambda = 2.2 * method_bonus[method] * year_factor
            reuse_in_drug_discovery = int(np.random.poisson(reuse_lambda) + 1)
            citation_count = int(
                np.random.poisson(reuse_in_drug_discovery * 6.5 + {"X-ray": 16, "Cryo-EM": 24, "NMR": 10, "AlphaFold": 18}[method])
            )
            commercial_applications = int(
                np.random.poisson(0.18 * reuse_in_drug_discovery + {"X-ray": 0.5, "Cryo-EM": 0.8, "NMR": 0.2, "AlphaFold": 1.0}[method])
            )
            sharing_cost_usd = float(np.random.lognormal(np.log(share_cost[method]), 0.22))
            research_value_usd = float(
                reuse_in_drug_discovery * 175000
                + citation_count * 8500
                + commercial_applications * 255000
                + np.random.normal(0, 40000)
            )
            open_time_months = float(np.clip(np.random.normal(open_months[method], 2.3), 4, None))
            closed_time_months = float(open_time_months * closed_multiplier[method])
            records.append(
                {
                    "structure_id": f"PDB-CS-{structure_id:04d}",
                    "deposit_year": int(year),
                    "method_used": method,
                    "reuse_count": reuse_in_drug_discovery,
                    "reuse_in_drug_discovery": reuse_in_drug_discovery,
                    "citation_count": citation_count,
                    "commercial_applications": commercial_applications,
                    "sharing_cost_usd": round(sharing_cost_usd, 2),
                    "research_value_usd": round(max(research_value_usd, 80000), 2),
                    "open_time_months": round(open_time_months, 2),
                    "closed_time_months": round(closed_time_months, 2),
                    "network_outputs": citation_count + commercial_applications,
                }
            )
            structure_id += 1

    frame = pd.DataFrame.from_records(records)
    yearly = (
        frame.groupby(["deposit_year", "method_used"]).size().unstack(fill_value=0).reindex(columns=methods, fill_value=0)
    )
    yearly["reuse_in_drug_discovery"] = frame.groupby("deposit_year")["reuse_in_drug_discovery"].sum()
    yearly = yearly.reset_index().rename(columns={"deposit_year": "year"})
    summary = metric_summary(frame, "Protein Structure Data")
    summary["method_counts"] = {k: int(v) for k, v in frame["method_used"].value_counts().sort_index().items()}
    return frame, yearly, summary


def simulate_clinical_case() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    years = np.arange(2015, 2025)
    deposits_year = np.array([12, 14, 16, 18, 20, 22, 24, 24, 25, 25])
    records = []
    trial_id = 1

    for year, n_deposits in zip(years, deposits_year):
        maturity = 1 + 0.05 * (year - 2015)
        for _ in range(int(n_deposits)):
            delay_center = 30 - 2.25 * (year - 2015)
            sharing_delay_months = float(np.clip(np.random.normal(delay_center, 4.5), 4, 36))
            reuse_lambda = max(0.8, 7.5 - 0.17 * sharing_delay_months + 0.18 * (year - 2015)) * maturity
            reuse_count = int(np.random.poisson(reuse_lambda) + 1)
            meta_analysis_inclusion = int(np.random.binomial(1, np.clip(0.16 + 0.035 * reuse_count - 0.005 * sharing_delay_months, 0.05, 0.92)))
            regulatory_impact = int(np.clip(np.random.poisson(0.18 * reuse_count + 0.65 * meta_analysis_inclusion), 0, 5))
            patient_benefit_score = float(
                np.clip(np.random.normal(34 + 3.6 * reuse_count + 8 * meta_analysis_inclusion + 5 * regulatory_impact - 0.8 * sharing_delay_months, 8), 0, 100)
            )
            sharing_cost_usd = float(np.random.lognormal(np.log(36000), 0.20))
            research_value_usd = float(
                reuse_count * 108000
                + meta_analysis_inclusion * 190000
                + regulatory_impact * 162000
                + patient_benefit_score * 5200
                + np.random.normal(0, 28000)
            )
            open_time_months = float(np.clip(np.random.normal(15 + sharing_delay_months / 4 - patient_benefit_score / 18, 1.8), 4, None))
            closed_time_months = float(open_time_months * (1.34 + sharing_delay_months / 45))
            records.append(
                {
                    "trial_id": f"CSDR-YODA-{trial_id:04d}",
                    "deposit_year": int(year),
                    "sharing_delay_months": round(sharing_delay_months, 2),
                    "reuse_count": reuse_count,
                    "meta_analysis_inclusion": meta_analysis_inclusion,
                    "regulatory_impact": regulatory_impact,
                    "patient_benefit_score": round(patient_benefit_score, 2),
                    "sharing_cost_usd": round(sharing_cost_usd, 2),
                    "research_value_usd": round(max(research_value_usd, 75000), 2),
                    "open_time_months": round(open_time_months, 2),
                    "closed_time_months": round(closed_time_months, 2),
                    "network_outputs": meta_analysis_inclusion + regulatory_impact,
                }
            )
            trial_id += 1

    frame = pd.DataFrame.from_records(records)
    yearly = (
        frame.groupby("deposit_year")
        .agg(
            deposits=("trial_id", "size"),
            mean_sharing_delay_months=("sharing_delay_months", "mean"),
            reuse_count=("reuse_count", "sum"),
            mean_patient_benefit_score=("patient_benefit_score", "mean"),
        )
        .reset_index()
        .rename(columns={"deposit_year": "year"})
    )
    summary = metric_summary(frame, "Clinical Trial Data Sharing")
    summary["mean_sharing_delay_months"] = round(float(frame["sharing_delay_months"].mean()), 2)
    summary["mean_patient_benefit_score"] = round(float(frame["patient_benefit_score"].mean()), 2)
    return frame, yearly, summary


def make_genomic_figure(yearly: pd.DataFrame) -> str:
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    ax1.bar(yearly["year"], yearly["deposits"], color=PALETTE["blue"], alpha=0.85, label="Deposits")
    ax2.plot(yearly["year"], yearly["reuse_events"], color=PALETTE["vermillion"], marker="o", linewidth=2.5, label="Reuse events")
    ax1.axvspan(2019.5, 2021.5, color=PALETTE["yellow"], alpha=0.18)
    ax1.text(2020.5, yearly["deposits"].max() * 1.05, "COVID-19 reuse surge", ha="center", va="bottom", color=PALETTE["grey"])

    ax1.set_title("Genomic Data Sharing Growth and Reuse")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Datasets deposited")
    ax2.set_ylabel("Annual reuse events")
    ax1.set_xticks(yearly["year"])

    lines, labels = [], []
    for ax in (ax1, ax2):
        h, l = ax.get_legend_handles_labels()
        lines.extend(h)
        labels.extend(l)
    ax1.legend(lines, labels, loc="upper left")
    fig.tight_layout()

    output = FIGURES_DIR / "case_genomic_growth.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(output)


def make_protein_figure(yearly: pd.DataFrame) -> str:
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    method_order = ["X-ray", "Cryo-EM", "NMR", "AlphaFold"]
    colors = [PALETTE["blue"], PALETTE["green"], PALETTE["purple"], PALETTE["orange"]]
    bottom = np.zeros(len(yearly))

    for method, color in zip(method_order, colors):
        ax1.bar(yearly["year"], yearly[method], bottom=bottom, label=method, color=color)
        bottom += yearly[method].to_numpy()

    ax2.plot(
        yearly["year"],
        yearly["reuse_in_drug_discovery"],
        color=PALETTE["vermillion"],
        marker="o",
        linewidth=2.5,
        label="Drug discovery reuse",
    )

    ax1.set_title("Protein Structure Methods and Reuse Revolution")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Structures deposited")
    ax2.set_ylabel("Reuse in drug discovery")
    ax1.set_xticks(yearly["year"])

    ax1.annotate("Cryo-EM scaling", xy=(2019, yearly.loc[yearly["year"] == 2019, "Cryo-EM"].iloc[0] + 10), xytext=(2016.8, 32), arrowprops={"arrowstyle": "->", "color": PALETTE["green"]}, color=PALETTE["green"])
    ax1.annotate("AlphaFold era", xy=(2023, yearly.loc[yearly["year"] == 2023, "AlphaFold"].iloc[0] + 8), xytext=(2020.8, 45), arrowprops={"arrowstyle": "->", "color": PALETTE["orange"]}, color=PALETTE["orange"])

    lines, labels = [], []
    for ax in (ax1, ax2):
        h, l = ax.get_legend_handles_labels()
        lines.extend(h)
        labels.extend(l)
    ax1.legend(lines, labels, loc="upper left")
    fig.tight_layout()

    output = FIGURES_DIR / "case_protein_revolution.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(output)


def make_clinical_figure(yearly: pd.DataFrame) -> str:
    fig, ax1 = plt.subplots(figsize=(11, 6))
    ax2 = ax1.twinx()

    ax1.bar(yearly["year"], yearly["deposits"], color=PALETTE["sky"], alpha=0.9, label="Shared trial datasets")
    ax2.plot(
        yearly["year"],
        yearly["mean_sharing_delay_months"],
        color=PALETTE["vermillion"],
        marker="o",
        linewidth=2.5,
        label="Mean sharing delay",
    )

    milestones = {
        2015: "Policy launch",
        2016: "Independent review model",
        2018: "Template harmonization",
        2020: "Pandemic-era urgency",
        2023: "Patient-centric access",
    }
    for year, label in milestones.items():
        ax1.axvline(year, color=PALETTE["grey"], linestyle="--", alpha=0.35)
        ax1.text(year + 0.05, yearly["deposits"].max() + 1.1, label, rotation=90, va="bottom", ha="left", fontsize=9)

    ax1.set_title("Clinical Trial Data Sharing Milestones")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Shared trial datasets")
    ax2.set_ylabel("Mean sharing delay (months)")
    ax1.set_xticks(yearly["year"])
    ax1.set_ylim(0, yearly["deposits"].max() + 7)

    lines, labels = [], []
    for ax in (ax1, ax2):
        h, l = ax.get_legend_handles_labels()
        lines.extend(h)
        labels.extend(l)
    ax1.legend(lines, labels, loc="upper left")
    fig.tight_layout()

    output = FIGURES_DIR / "case_clinical_sharing.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(output)


def make_roi_figure(summary_rows: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(9, 6))
    names = [row["case_name"] for row in summary_rows]
    rois = [row["roi"] for row in summary_rows]
    lower_err = [row["roi"] - row["roi_ci_95"][0] for row in summary_rows]
    upper_err = [row["roi_ci_95"][1] - row["roi"] for row in summary_rows]
    colors = [PALETTE["blue"], PALETTE["green"], PALETTE["orange"]]

    bars = ax.bar(names, rois, color=colors, yerr=np.vstack([lower_err, upper_err]), capsize=6)
    ax.set_title("ROI of Open Data Sharing Across Life Science Cases")
    ax.set_ylabel("ROI (research value / sharing cost)")
    ax.set_ylim(0, max(rois) * 1.25)
    ax.tick_params(axis="x", rotation=12)

    for bar, roi in zip(bars, rois):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.06, f"{roi:.2f}", ha="center", va="bottom")

    fig.tight_layout()
    output = FIGURES_DIR / "case_roi_comparison.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(output)


def save_supporting_files(genomic: pd.DataFrame, protein: pd.DataFrame, clinical: pd.DataFrame, results: dict) -> None:
    genomic.to_csv(DATA_DIR / "genomic_case_study.csv", index=False)
    protein.to_csv(DATA_DIR / "protein_case_study.csv", index=False)
    clinical.to_csv(DATA_DIR / "clinical_trial_case_study.csv", index=False)

    PREPROCESSING_LOG.write_text(
        "# Preprocessing Log\n\n"
        "- Random seeds fixed at 42 for numpy and random.\n"
        "- Simulated deposit cohorts for 2015-2024 using case-specific growth assumptions.\n"
        "- ROI and acceleration factors were estimated with nonparametric bootstrap 95% confidence intervals (2,000 resamples).\n"
        "- Network effects were summarized with log-log reuse elasticity and Spearman association between reuse and downstream outputs.\n"
        "- No parametric hypothesis tests were used, so distributional assumptions were minimized.\n",
        encoding="utf-8",
    )

    results_payload = json.dumps(results, indent=2)
    RESULTS_JSON.write_text(results_payload + "\n", encoding="utf-8")

    summary_lines = [
        "# Statistical Summary",
        "",
        "Bootstrap-based intervals were used for ROI, acceleration factor, and network elasticity.",
        "Spearman correlation was used for network dependence to avoid normality assumptions.",
        "",
    ]
    for case_key in ("genomic_data_sharing", "protein_structure_data", "clinical_trial_data_sharing"):
        case = results["case_studies"][case_key]
        summary_lines.extend(
            [
                f"## {case['case_name']}",
                f"- ROI: {case['roi']:.3f} (95% CI {case['roi_ci_95'][0]:.3f} to {case['roi_ci_95'][1]:.3f})",
                f"- Acceleration factor: {case['acceleration_factor']:.3f} (95% CI {case['acceleration_factor_ci_95'][0]:.3f} to {case['acceleration_factor_ci_95'][1]:.3f})",
                f"- Network elasticity: {case['network_elasticity']:.3f} (95% CI {case['network_elasticity_ci_95'][0]:.3f} to {case['network_elasticity_ci_95'][1]:.3f})",
                f"- Spearman rho: {case['network_spearman_rho']:.3f}; p={case['network_spearman_p']:.3e}",
                "",
            ]
        )
    STATS_SUMMARY.write_text("\n".join(summary_lines), encoding="utf-8")

    REPORT_PATH.write_text(
        "# DRAFT — NOT FOR DISTRIBUTION\n\n"
        "## Life science open data case study\n\n"
        "This workspace contains a simulated case study of open-data impact across genomic repositories, protein structure resources, and shared clinical trial datasets."
        " Figures use English labels and colorblind-friendly palettes. Numeric summaries and confidence intervals are saved in `results/life_science_results.json` and `results/statistical-summary.md`."
        "\n\n## File inventory\n"
        "- `src/life_science_case_study.py`: simulation and plotting workflow.\n"
        "- `data/*.csv`: simulated datasets for the three case studies.\n"
        "- `figures/*.png`: publication-style figures at 300 DPI.\n"
        "- `results/life_science_results.json`: machine-readable summary.\n"
        "- `results/statistical-summary.md`: human-readable metric summary.\n",
        encoding="utf-8",
    )


def main() -> None:
    ensure_dirs()

    genomic_frame, genomic_yearly, genomic_summary = simulate_genomic_case()
    protein_frame, protein_yearly, protein_summary = simulate_protein_case()
    clinical_frame, clinical_yearly, clinical_summary = simulate_clinical_case()

    figure_paths = {
        "case_genomic_growth": make_genomic_figure(genomic_yearly),
        "case_protein_revolution": make_protein_figure(protein_yearly),
        "case_clinical_sharing": make_clinical_figure(clinical_yearly),
        "case_roi_comparison": make_roi_figure([genomic_summary, protein_summary, clinical_summary]),
    }

    results = {
        "seed": SEED,
        "methodology": {
            "roi_definition": "total estimated research value generated divided by total sharing cost",
            "acceleration_factor_definition": "estimated closed-access discovery time divided by open-data discovery time",
            "network_effect_definition": "log-log elasticity between reuse and downstream outputs with Spearman correlation",
        },
        "case_studies": {
            "genomic_data_sharing": genomic_summary,
            "protein_structure_data": protein_summary,
            "clinical_trial_data_sharing": clinical_summary,
        },
        "figures": figure_paths,
    }

    save_supporting_files(genomic_frame, protein_frame, clinical_frame, results)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
