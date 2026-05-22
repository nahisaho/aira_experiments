"""
Epigenetic Clock Data Simulator
Generates realistic DNA methylation beta-values with age-associated patterns,
tissue-specific effects, and intervention responses.
"""
import numpy as np
import pandas as pd
from pathlib import Path

SEED = 42
np.random.seed(SEED)

N_CPG_TOTAL = 2000
TISSUES = ["blood", "brain", "liver", "skin", "muscle"]


def _generate_age_coefficients(n_cpg: int) -> np.ndarray:
    coefs = np.zeros(n_cpg)
    n_active = int(n_cpg * 0.25)
    active_idx = np.random.choice(n_cpg, n_active, replace=False)
    coefs[active_idx] = np.random.normal(0, 0.02, n_active)
    strong = np.random.choice(active_idx, min(20, n_active), replace=False)
    coefs[strong] = np.random.normal(0, 0.08, len(strong))
    return coefs


def _tissue_effect(tissue: str, n_cpg: int) -> np.ndarray:
    tissue_seeds = {"blood": 1, "brain": 2, "liver": 3, "skin": 4, "muscle": 5}
    rng = np.random.RandomState(tissue_seeds.get(tissue, 0))
    return rng.normal(0, 0.05, n_cpg)


def generate_methylation_dataset(
    n_samples: int = 500,
    n_cpg: int = N_CPG_TOTAL,
    tissue: str = "blood",
    age_range: tuple = (20, 90),
    include_interventions: bool = False,
    longevity_cohort: bool = False,
) -> pd.DataFrame:
    ages = np.random.uniform(*age_range, n_samples)
    if longevity_cohort:
        ages = np.random.uniform(80, 105, n_samples)

    sex = np.random.binomial(1, 0.5, n_samples)
    coefs = _generate_age_coefficients(n_cpg)
    tissue_offset = _tissue_effect(tissue, n_cpg)
    age_transformed = np.log(ages + 1)

    beta = np.outer(age_transformed, coefs) + tissue_offset
    beta += np.random.normal(0, 0.03, (n_samples, n_cpg))
    sex_cpgs = np.random.choice(n_cpg, 50, replace=False)
    beta[:, sex_cpgs] += np.outer(sex, np.random.normal(0, 0.02, 50))
    beta = 1 / (1 + np.exp(-beta * 5 + 2.5))

    intervention = np.full(n_samples, "none")
    if include_interventions:
        n_each = n_samples // 5
        for i, intv in enumerate(["exercise", "caloric_restriction", "metformin", "rapamycin"]):
            start = (i + 1) * n_each
            end = min(start + n_each, n_samples)
            intervention[start:end] = intv
            effect_strength = {"exercise": 0.3, "caloric_restriction": 0.4,
                               "metformin": 0.35, "rapamycin": 0.5}
            effect_cpgs = np.random.choice(n_cpg, 100, replace=False)
            beta[start:end, effect_cpgs] -= (
                coefs[effect_cpgs] * effect_strength[intv] *
                age_transformed[start:end, np.newaxis] * 0.1
            )

    beta = np.clip(beta, 0.001, 0.999)
    age_accel = np.random.normal(0, 3, n_samples)
    if longevity_cohort:
        age_accel -= 5

    cpg_names = [f"cg{i:08d}" for i in range(n_cpg)]
    df = pd.DataFrame(beta, columns=cpg_names)
    df.insert(0, "chronological_age", ages)
    df.insert(1, "biological_age_offset", age_accel)
    df.insert(2, "sex", sex)
    df.insert(3, "tissue", tissue)
    df.insert(4, "intervention", intervention)
    df["true_bio_age"] = ages + age_accel
    return df


def generate_all_datasets() -> dict:
    datasets = {}
    print("Generating blood training data (n=800)...")
    datasets["blood_train"] = generate_methylation_dataset(800, tissue="blood")
    for tissue in TISSUES:
        print(f"Generating {tissue} data (n=300)...")
        datasets[tissue] = generate_methylation_dataset(300, tissue=tissue)
    print("Generating intervention data (n=500)...")
    datasets["intervention"] = generate_methylation_dataset(
        500, tissue="blood", include_interventions=True
    )
    print("Generating longevity cohort (n=200)...")
    datasets["longevity"] = generate_methylation_dataset(
        200, tissue="blood", longevity_cohort=True
    )
    return datasets


if __name__ == "__main__":
    out_dir = Path("data")
    out_dir.mkdir(exist_ok=True)
    datasets = generate_all_datasets()
    for name, df in datasets.items():
        path = out_dir / f"{name}.csv"
        df.to_csv(path, index=False)
        print(f"  Saved {path} — {df.shape}")
    print("Done.")
