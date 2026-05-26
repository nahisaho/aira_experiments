"""
Synthetic DNA methylation data generator for epigenetic clock experiments.
Simulates realistic CpG methylation patterns with tissue-specific effects,
age-related changes, and intervention effects.
"""
import numpy as np
import pandas as pd

np.random.seed(42)

N_SAMPLES = 1200
N_CPG_SITES = 500
TISSUES = ['blood', 'brain', 'liver', 'muscle', 'skin']
INTERVENTIONS = ['none', 'exercise', 'diet', 'drug']

def generate_age_related_cpgs(n_cpgs, n_samples, ages):
    """Generate CpG sites with age-correlated methylation changes."""
    beta_values = np.zeros((n_samples, n_cpgs))
    
    # Hyper-methylating CpGs (increase with age)
    n_hyper = n_cpgs // 3
    for i in range(n_hyper):
        rate = np.random.uniform(0.002, 0.008)
        base = np.random.uniform(0.1, 0.4)
        noise = np.random.normal(0, 0.03, n_samples)
        beta_values[:, i] = np.clip(base + rate * ages + noise, 0, 1)
    
    # Hypo-methylating CpGs (decrease with age)
    n_hypo = n_cpgs // 3
    for i in range(n_hyper, n_hyper + n_hypo):
        rate = np.random.uniform(-0.008, -0.002)
        base = np.random.uniform(0.6, 0.9)
        noise = np.random.normal(0, 0.03, n_samples)
        beta_values[:, i] = np.clip(base + rate * ages + noise, 0, 1)
    
    # Non-linear age-related CpGs
    n_nonlinear = n_cpgs - n_hyper - n_hypo
    for i in range(n_hyper + n_hypo, n_cpgs):
        base = np.random.uniform(0.3, 0.7)
        noise = np.random.normal(0, 0.04, n_samples)
        beta_values[:, i] = np.clip(
            base + 0.003 * ages + 0.00005 * (ages - 50)**2 + noise, 0, 1
        )
    
    return beta_values


def add_tissue_effects(beta_values, tissues):
    """Add tissue-specific methylation offsets."""
    tissue_offsets = {
        'blood': np.random.normal(0, 0.02, beta_values.shape[1]),
        'brain': np.random.normal(0.05, 0.03, beta_values.shape[1]),
        'liver': np.random.normal(-0.03, 0.025, beta_values.shape[1]),
        'muscle': np.random.normal(0.02, 0.02, beta_values.shape[1]),
        'skin': np.random.normal(-0.01, 0.015, beta_values.shape[1]),
    }
    modified = beta_values.copy()
    for i, tissue in enumerate(tissues):
        modified[i] += tissue_offsets[tissue]
    return np.clip(modified, 0, 1)


def add_intervention_effects(beta_values, ages, interventions):
    """Simulate intervention effects on methylation (age deceleration)."""
    modified = beta_values.copy()
    effect_cpgs = list(range(0, 50))  # First 50 CpGs are intervention-sensitive
    
    for i, intervention in enumerate(interventions):
        if intervention == 'exercise':
            modified[i, effect_cpgs] -= np.random.uniform(0.005, 0.015, len(effect_cpgs))
        elif intervention == 'diet':
            modified[i, effect_cpgs] -= np.random.uniform(0.003, 0.010, len(effect_cpgs))
        elif intervention == 'drug':
            modified[i, effect_cpgs] -= np.random.uniform(0.008, 0.020, len(effect_cpgs))
    
    return np.clip(modified, 0, 1)


def add_biological_age_deviation(ages):
    """Generate biological ages with deviation from chronological age."""
    bio_ages = ages + np.random.normal(0, 3.5, len(ages))
    # Add some outliers (accelerated agers)
    n_outliers = len(ages) // 20
    outlier_idx = np.random.choice(len(ages), n_outliers, replace=False)
    bio_ages[outlier_idx] += np.random.uniform(5, 15, n_outliers)
    return bio_ages


def generate_dataset():
    """Generate the full synthetic dataset."""
    ages = np.random.uniform(20, 90, N_SAMPLES)
    tissues = np.random.choice(TISSUES, N_SAMPLES)
    interventions = np.random.choice(INTERVENTIONS, N_SAMPLES, p=[0.5, 0.2, 0.15, 0.15])
    
    beta_values = generate_age_related_cpgs(N_CPG_SITES, N_SAMPLES, ages)
    beta_values = add_tissue_effects(beta_values, tissues)
    beta_values = add_intervention_effects(beta_values, ages, interventions)
    
    bio_ages = add_biological_age_deviation(ages)
    
    cpg_names = [f'cg{str(i).zfill(8)}' for i in range(N_CPG_SITES)]
    
    df_methylation = pd.DataFrame(beta_values, columns=cpg_names)
    df_metadata = pd.DataFrame({
        'chronological_age': ages,
        'biological_age': bio_ages,
        'tissue': tissues,
        'intervention': interventions,
        'age_acceleration': bio_ages - ages,
    })
    
    return df_methylation, df_metadata


if __name__ == '__main__':
    df_meth, df_meta = generate_dataset()
    df_meth.to_csv('data_methylation.csv', index=False)
    df_meta.to_csv('data_metadata.csv', index=False)
    print(f"Generated methylation data: {df_meth.shape}")
    print(f"Generated metadata: {df_meta.shape}")
    print(f"Tissues: {df_meta['tissue'].value_counts().to_dict()}")
    print(f"Interventions: {df_meta['intervention'].value_counts().to_dict()}")
