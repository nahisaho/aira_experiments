from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, RESULTS_DIR, SEED

TARGET_SCORES_PATH = RESULTS_DIR / "target_scores.csv"
TARGET_FIG_PATH = FIGURES_DIR / "target_selection.png"


plt.style.use("seaborn-v0_8-whitegrid")


def _normalize(series: pd.Series) -> pd.Series:
    span = series.max() - series.min()
    if float(span) == 0.0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - series.min()) / span


def score_targets(catalog: pd.DataFrame, w1: float = 0.55, w2: float = 0.45) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    altitude = catalog["altitude_km"]
    area = catalog["area_m2"]
    cross_section = catalog["radar_cross_section_m2"]
    incl_rad = np.deg2rad(catalog["inclination_deg"])

    density = (
        0.9 * np.exp(-0.5 * ((altitude - 750.0) / 120.0) ** 2)
        + 1.2 * np.exp(-0.5 * ((altitude - 950.0) / 140.0) ** 2)
        + 0.6 * np.exp(-0.5 * ((altitude - 1350.0) / 180.0) ** 2)
    )
    relative_velocity_km_s = 9.5 * np.sqrt(np.sin(incl_rad / 2.0) ** 2 + 0.15 * catalog["eccentricity"] + 0.02)
    collision_proxy = density * cross_section * relative_velocity_km_s
    removal_effect = catalog["mass_kg"] * area / catalog["decay_lifetime_days"]

    scored = catalog.copy()
    scored["orbital_density_score"] = _normalize(pd.Series(density))
    scored["relative_velocity_km_s"] = relative_velocity_km_s
    scored["collision_probability_score"] = _normalize(pd.Series(collision_proxy))
    scored["removal_effect_score"] = _normalize(pd.Series(removal_effect))
    scored["combined_score"] = w1 * scored["collision_probability_score"] + w2 * scored["removal_effect_score"]
    scored["score_rank"] = scored["combined_score"].rank(ascending=False, method="first").astype(int)
    scored = scored.sort_values("combined_score", ascending=False).reset_index(drop=True)

    top10 = scored.head(10).copy()
    scored.to_csv(TARGET_SCORES_PATH, index=False)

    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    scatter = ax.scatter(
        scored["altitude_km"],
        scored["inclination_deg"],
        c=scored["combined_score"],
        s=50 + 180 * _normalize(scored["mass_kg"]),
        cmap="viridis",
        alpha=0.8,
        edgecolor="black",
        linewidth=0.3,
    )
    ax.scatter(
        top10["altitude_km"],
        top10["inclination_deg"],
        facecolors="none",
        edgecolors="crimson",
        s=220,
        linewidth=1.8,
        label="Top 10 targets",
    )
    for _, row in top10.iterrows():
        ax.annotate(row["debris_id"], (row["altitude_km"], row["inclination_deg"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Combined target score")
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Inclination (deg)")
    ax.set_title("ADR target scoring across the debris catalog")
    ax.legend(loc="upper right")
    fig.savefig(TARGET_FIG_PATH, dpi=300)
    plt.close(fig)

    summary = {
        "top_target": top10.iloc[0]["debris_id"],
        "top_score": float(top10.iloc[0]["combined_score"]),
        "mean_collision_probability_score": float(scored["collision_probability_score"].mean()),
    }
    (RESULTS_DIR / "target_selection_summary.json").write_text(json.dumps(summary, indent=2))
    return scored, top10


if __name__ == "__main__":
    from .debris_catalog import generate_debris_catalog

    catalog = generate_debris_catalog()
    scored_catalog, top_targets = score_targets(catalog)
    print(top_targets[["debris_id", "combined_score"]].to_string(index=False))
