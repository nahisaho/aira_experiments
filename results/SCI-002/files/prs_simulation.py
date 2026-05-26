#!/usr/bin/env python3
"""
Cross-Ancestry PRS Transferability Simulation Framework

Simulates polygenic risk score (PRS) transferability from a European-ancestry
population (modeled on UK Biobank) to an East Asian-ancestry population
(modeled on BioBank Japan), with Bayesian LD correction, multi-ancestry
meta-analysis, and local ancestry-informed PRS adjustment.

Case study: Type 2 Diabetes (T2D)
"""

import numpy as np
import pandas as pd
from scipy import stats, linalg
from scipy.special import expit
from sklearn.metrics import roc_auc_score, r2_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
import os
import json

warnings.filterwarnings('ignore')
np.random.seed(42)

FIGURE_DIR = "figures"
os.makedirs(FIGURE_DIR, exist_ok=True)

# ============================================================
# 1. Population Genetic Parameters
# ============================================================

class PopulationParams:
    """Parameters for simulating two diverged populations."""
    def __init__(self, n_snps=500, n_causal=50, fst=0.1,
                 n_eur=10000, n_eas=5000, h2=0.5,
                 prevalence=0.10, ld_block_size=10):
        self.n_snps = n_snps
        self.n_causal = n_causal
        self.fst = fst
        self.n_eur = n_eur
        self.n_eas = n_eas
        self.h2 = h2
        self.prevalence = prevalence
        self.ld_block_size = ld_block_size


# ============================================================
# 2. LD Matrix Generation with Population-Specific Structure
# ============================================================

def generate_ld_matrix(n_snps, block_size, decay_rate=0.8, noise_scale=0.05):
    """Generate a block-diagonal LD matrix with exponential decay."""
    R = np.eye(n_snps)
    n_blocks = n_snps // block_size
    for b in range(n_blocks):
        start = b * block_size
        end = min(start + block_size, n_snps)
        for i in range(start, end):
            for j in range(i + 1, end):
                r = decay_rate ** (j - i)
                r += np.random.normal(0, noise_scale)
                r = np.clip(r, -0.99, 0.99)
                R[i, j] = r
                R[j, i] = r
    # Ensure positive semi-definite
    eigvals = np.linalg.eigvalsh(R)
    if eigvals.min() < 0:
        R += (-eigvals.min() + 1e-6) * np.eye(n_snps)
        d = np.sqrt(np.diag(R))
        R = R / np.outer(d, d)
    return R


def generate_population_ld(params):
    """Generate diverged LD matrices for EUR and EAS populations."""
    R_eur = generate_ld_matrix(params.n_snps, params.ld_block_size,
                               decay_rate=0.85, noise_scale=0.03)
    R_eas = generate_ld_matrix(params.n_snps, params.ld_block_size,
                               decay_rate=0.75, noise_scale=0.05)
    return R_eur, R_eas


# ============================================================
# 3. Allele Frequency Divergence via Fst (Balding-Nichols Model)
# ============================================================

def simulate_allele_frequencies(n_snps, fst):
    """Simulate diverged allele frequencies using the Balding-Nichols model."""
    p_anc = np.random.uniform(0.05, 0.95, n_snps)
    p_eur = np.zeros(n_snps)
    p_eas = np.zeros(n_snps)

    for i in range(n_snps):
        a = p_anc[i] * (1 - fst) / fst
        b = (1 - p_anc[i]) * (1 - fst) / fst
        if a > 0 and b > 0:
            p_eur[i] = np.clip(np.random.beta(a, b), 0.01, 0.99)
            p_eas[i] = np.clip(np.random.beta(a, b), 0.01, 0.99)
        else:
            p_eur[i] = p_anc[i]
            p_eas[i] = p_anc[i]

    return p_anc, p_eur, p_eas


# ============================================================
# 4. Genotype Simulation with LD Structure
# ============================================================

def simulate_genotypes(n_individuals, allele_freqs, R):
    """Simulate genotypes incorporating LD structure via Cholesky decomposition."""
    n_snps = len(allele_freqs)
    L = np.linalg.cholesky(R + 1e-6 * np.eye(n_snps))
    Z = np.random.randn(n_individuals, n_snps)
    Z_corr = Z @ L.T
    # Convert to genotypes using allele frequencies as thresholds
    thresholds = stats.norm.ppf(1 - allele_freqs)
    genotypes = np.zeros((n_individuals, n_snps), dtype=np.float64)
    for allele in range(2):
        Z_draw = np.random.randn(n_individuals, n_snps) @ L.T
        genotypes += (Z_draw > thresholds).astype(np.float64)
    return genotypes


# ============================================================
# 5. True Effect Sizes and Phenotype Simulation
# ============================================================

def simulate_true_effects(params):
    """Simulate true causal SNP effects with partial sharing across pops."""
    causal_idx = np.random.choice(params.n_snps, params.n_causal, replace=False)
    beta_true = np.zeros(params.n_snps)
    beta_true[causal_idx] = np.random.normal(0, np.sqrt(params.h2 / params.n_causal),
                                              params.n_causal)
    # Population-specific effect modifier (correlation ~0.8 across pops)
    beta_eur = beta_true.copy()
    beta_eas = beta_true.copy()
    modifier = np.random.normal(1.0, 0.2, params.n_snps)
    beta_eas *= modifier
    # Re-normalize to maintain h2
    beta_eas *= np.sqrt(params.h2 / params.n_causal) / (np.std(beta_eas[causal_idx]) + 1e-10)

    return causal_idx, beta_true, beta_eur, beta_eas


def simulate_phenotype_continuous(genotypes, beta, h2):
    """Simulate continuous phenotype."""
    genetic_value = genotypes @ beta
    var_g = np.var(genetic_value)
    if var_g < 1e-10:
        var_g = 1e-10
    var_e = var_g * (1 - h2) / h2
    noise = np.random.normal(0, np.sqrt(var_e), len(genotypes))
    y = genetic_value + noise
    return y


def simulate_phenotype_binary(genotypes, beta, h2, prevalence):
    """Simulate binary (disease) phenotype using liability threshold model."""
    y_cont = simulate_phenotype_continuous(genotypes, beta, h2)
    threshold = np.quantile(y_cont, 1 - prevalence)
    y_binary = (y_cont >= threshold).astype(int)
    return y_binary, y_cont


# ============================================================
# 6. GWAS: Marginal Association Testing
# ============================================================

def run_gwas(genotypes, phenotype, binary=False):
    """Run marginal GWAS (linear or logistic regression per SNP)."""
    n, p = genotypes.shape
    beta_hat = np.zeros(p)
    se_hat = np.zeros(p)
    pvalues = np.zeros(p)

    y = phenotype - np.mean(phenotype)

    for j in range(p):
        x = genotypes[:, j]
        x_centered = x - np.mean(x)
        var_x = np.var(x_centered)
        if var_x < 1e-10:
            beta_hat[j] = 0
            se_hat[j] = 1.0
            pvalues[j] = 1.0
            continue
        b = np.dot(x_centered, y) / (np.dot(x_centered, x_centered))
        residuals = y - b * x_centered
        se = np.sqrt(np.sum(residuals ** 2) / ((n - 2) * np.sum(x_centered ** 2)))
        beta_hat[j] = b
        se_hat[j] = se
        if se > 0:
            z = b / se
            pvalues[j] = 2 * stats.norm.sf(abs(z))
        else:
            pvalues[j] = 1.0

    return beta_hat, se_hat, pvalues


# ============================================================
# 7. PRS Methods
# ============================================================

def prs_clump_threshold(beta_hat, pvalues, R, p_threshold=5e-2):
    """Simple clumping and thresholding PRS."""
    selected = pvalues < p_threshold
    weights = np.zeros_like(beta_hat)
    weights[selected] = beta_hat[selected]
    return weights


def prs_standard(genotypes, weights):
    """Calculate PRS given genotypes and weights."""
    return genotypes @ weights


# ============================================================
# 8. Bayesian LD Correction (inspired by PRS-CS/PRS-CSx)
# ============================================================

def bayesian_ld_correction(beta_hat_eur, se_eur, R_eur, R_eas, n_eur,
                           phi=1e-2, n_iter=100):
    """
    Bayesian shrinkage estimator that adjusts EUR GWAS effects
    for EAS LD structure using a continuous shrinkage prior.

    beta_eas_adjusted ~ N(D_eas * R_eas^{-1} * R_eur * beta_hat_eur, ...)
    with global-local shrinkage prior.
    """
    p = len(beta_hat_eur)
    # Regularized inverse of EUR LD
    R_eur_reg = R_eur + phi * np.eye(p)
    R_eur_inv = np.linalg.inv(R_eur_reg)
    # Transformation matrix: project from EUR LD space to EAS LD space
    R_eas_reg = R_eas + phi * np.eye(p)
    # Bayesian posterior mean under Gaussian prior
    precision = n_eur * R_eur_reg + (1.0 / phi) * np.eye(p)
    precision_inv = np.linalg.inv(precision)
    # Posterior mean in EUR LD space
    posterior_mean_eur = precision_inv @ (n_eur * R_eur_reg @ beta_hat_eur)
    # Transform to EAS LD space
    transfer_matrix = R_eas_reg @ np.linalg.inv(R_eur_reg)
    beta_adjusted = transfer_matrix @ posterior_mean_eur

    return beta_adjusted


# ============================================================
# 9. Multi-Ancestry Meta-Analysis
# ============================================================

def multi_ancestry_meta_analysis(beta_eur, se_eur, beta_eas, se_eas,
                                  method='fixed'):
    """
    Fixed-effect or random-effect meta-analysis of SNP effects
    across EUR and EAS GWAS.
    """
    p = len(beta_eur)
    beta_meta = np.zeros(p)
    se_meta = np.zeros(p)

    for j in range(p):
        w_eur = 1.0 / (se_eur[j] ** 2 + 1e-10)
        w_eas = 1.0 / (se_eas[j] ** 2 + 1e-10)

        if method == 'fixed':
            beta_meta[j] = (w_eur * beta_eur[j] + w_eas * beta_eas[j]) / (w_eur + w_eas)
            se_meta[j] = np.sqrt(1.0 / (w_eur + w_eas))
        elif method == 'random':
            # DerSimonian-Laird random effects
            Q = w_eur * (beta_eur[j] - beta_meta[j]) ** 2 + \
                w_eas * (beta_eas[j] - beta_meta[j]) ** 2
            tau2 = max(0, (Q - 1) / (w_eur + w_eas - (w_eur ** 2 + w_eas ** 2) / (w_eur + w_eas)))
            w_eur_re = 1.0 / (se_eur[j] ** 2 + tau2 + 1e-10)
            w_eas_re = 1.0 / (se_eas[j] ** 2 + tau2 + 1e-10)
            beta_meta[j] = (w_eur_re * beta_eur[j] + w_eas_re * beta_eas[j]) / (w_eur_re + w_eas_re)
            se_meta[j] = np.sqrt(1.0 / (w_eur_re + w_eas_re))

    return beta_meta, se_meta


# ============================================================
# 10. Local Ancestry-Informed PRS Correction
# ============================================================

def simulate_local_ancestry(n_individuals, n_snps, prop_eur=0.0, n_segments=10):
    """
    Simulate local ancestry assignments along the genome.
    0 = EAS ancestry, 1 = EUR ancestry at each locus.
    For a purely EAS population, prop_eur ≈ 0.
    """
    segment_size = n_snps // n_segments
    lanc = np.zeros((n_individuals, n_snps), dtype=int)

    for i in range(n_individuals):
        for s in range(n_segments):
            start = s * segment_size
            end = min(start + segment_size, n_snps)
            if np.random.random() < prop_eur:
                lanc[i, start:end] = 1
            else:
                lanc[i, start:end] = 0
    return lanc


def local_ancestry_prs(genotypes, beta_eur, beta_eas_adjusted, local_ancestry):
    """
    Compute PRS using ancestry-specific weights based on local ancestry.
    At each locus, use EUR weights where local ancestry = EUR,
    and EAS-adjusted weights where local ancestry = EAS.
    """
    weights = np.where(local_ancestry == 1,
                       beta_eur[np.newaxis, :],
                       beta_eas_adjusted[np.newaxis, :])
    prs = np.sum(genotypes * weights, axis=1)
    return prs


# ============================================================
# 11. Evaluation Metrics
# ============================================================

def evaluate_prs(prs_scores, y_true, y_liability=None, binary=True):
    """Evaluate PRS performance."""
    results = {}
    if binary:
        if len(np.unique(y_true)) > 1:
            results['AUC'] = roc_auc_score(y_true, prs_scores)
        else:
            results['AUC'] = 0.5
    if y_liability is not None:
        corr = np.corrcoef(prs_scores, y_liability)[0, 1]
        results['R2_liability'] = corr ** 2
    corr_obs = np.corrcoef(prs_scores, y_true)[0, 1]
    results['R2_observed'] = corr_obs ** 2
    return results


# ============================================================
# 12. Main Simulation Pipeline
# ============================================================

def run_simulation(params, verbose=True):
    """Run the complete PRS transferability simulation."""
    if verbose:
        print(f"=== Simulation: Fst={params.fst}, N_EUR={params.n_eur}, "
              f"N_EAS={params.n_eas}, h2={params.h2} ===")

    # Generate LD matrices
    R_eur, R_eas = generate_population_ld(params)

    # Simulate allele frequencies
    _, p_eur, p_eas = simulate_allele_frequencies(params.n_snps, params.fst)

    # Simulate true effects
    causal_idx, beta_true, beta_eur, beta_eas = simulate_true_effects(params)

    # Simulate genotypes
    G_eur = simulate_genotypes(params.n_eur, p_eur, R_eur)
    G_eas = simulate_genotypes(params.n_eas, p_eas, R_eas)

    # Simulate phenotypes (binary: T2D)
    y_eur_bin, y_eur_liab = simulate_phenotype_binary(G_eur, beta_eur,
                                                       params.h2, params.prevalence)
    y_eas_bin, y_eas_liab = simulate_phenotype_binary(G_eas, beta_eas,
                                                       params.h2, params.prevalence)

    # Run GWAS in both populations
    beta_hat_eur, se_eur, pval_eur = run_gwas(G_eur, y_eur_liab)
    beta_hat_eas, se_eas, pval_eas = run_gwas(G_eas, y_eas_liab)

    # === Method 1: Direct Transfer (C+T from EUR) ===
    w_ct_eur = prs_clump_threshold(beta_hat_eur, pval_eur, R_eur)
    prs_direct = prs_standard(G_eas, w_ct_eur)
    eval_direct = evaluate_prs(prs_direct, y_eas_bin, y_eas_liab)

    # === Method 2: Target-population PRS (oracle) ===
    w_ct_eas = prs_clump_threshold(beta_hat_eas, pval_eas, R_eas)
    prs_target = prs_standard(G_eas, w_ct_eas)
    eval_target = evaluate_prs(prs_target, y_eas_bin, y_eas_liab)

    # === Method 3: Bayesian LD Correction ===
    beta_bayes = bayesian_ld_correction(beta_hat_eur, se_eur, R_eur, R_eas,
                                         params.n_eur)
    prs_bayes = prs_standard(G_eas, beta_bayes)
    eval_bayes = evaluate_prs(prs_bayes, y_eas_bin, y_eas_liab)

    # === Method 4: Multi-Ancestry Meta-Analysis ===
    beta_meta, se_meta = multi_ancestry_meta_analysis(
        beta_hat_eur, se_eur, beta_hat_eas, se_eas, method='fixed')
    prs_meta = prs_standard(G_eas, beta_meta)
    eval_meta = evaluate_prs(prs_meta, y_eas_bin, y_eas_liab)

    # === Method 5: Local Ancestry PRS ===
    lanc_eas = simulate_local_ancestry(params.n_eas, params.n_snps, prop_eur=0.05)
    prs_lanc = local_ancestry_prs(G_eas, w_ct_eur, beta_bayes, lanc_eas)
    eval_lanc = evaluate_prs(prs_lanc, y_eas_bin, y_eas_liab)

    # === Method 6: Combined (Meta + Bayesian LD + LA) ===
    beta_combined = 0.5 * beta_meta + 0.5 * beta_bayes
    prs_combined_base = prs_standard(G_eas, beta_combined)
    # LA correction on combined
    prs_combined = local_ancestry_prs(G_eas, beta_meta, beta_combined, lanc_eas)
    eval_combined = evaluate_prs(prs_combined, y_eas_bin, y_eas_liab)

    results = {
        'Direct Transfer (EUR→EAS)': eval_direct,
        'Target Pop (EAS GWAS)': eval_target,
        'Bayesian LD Correction': eval_bayes,
        'Multi-Ancestry Meta': eval_meta,
        'Local Ancestry PRS': eval_lanc,
        'Combined (Proposed)': eval_combined,
    }

    if verbose:
        print(f"  {'Method':<30} {'AUC':>8} {'R²(liab)':>10} {'R²(obs)':>10}")
        print("  " + "-" * 60)
        for method, metrics in results.items():
            print(f"  {method:<30} {metrics.get('AUC', 0):.4f}   "
                  f"{metrics.get('R2_liability', 0):.4f}     "
                  f"{metrics.get('R2_observed', 0):.4f}")

    return results, {
        'R_eur': R_eur, 'R_eas': R_eas,
        'p_eur': p_eur, 'p_eas': p_eas,
        'beta_hat_eur': beta_hat_eur, 'beta_hat_eas': beta_hat_eas,
        'pval_eur': pval_eur, 'pval_eas': pval_eas,
        'beta_true': beta_true, 'beta_eur': beta_eur, 'beta_eas': beta_eas,
        'causal_idx': causal_idx,
        'prs_direct': prs_direct, 'prs_target': prs_target,
        'prs_bayes': prs_bayes, 'prs_meta': prs_meta,
        'prs_lanc': prs_lanc, 'prs_combined': prs_combined,
        'y_eas_bin': y_eas_bin, 'y_eas_liab': y_eas_liab,
    }


# ============================================================
# 13. Parameter Sweep Experiments
# ============================================================

def sweep_fst(fst_values=[0.01, 0.05, 0.1, 0.15, 0.2]):
    """Sweep over Fst values to assess impact of population divergence."""
    all_results = {}
    for fst in fst_values:
        params = PopulationParams(fst=fst)
        results, _ = run_simulation(params, verbose=False)
        all_results[fst] = results
    return all_results


def sweep_sample_size(n_eas_values=[500, 1000, 2000, 5000, 10000]):
    """Sweep over target population sample sizes."""
    all_results = {}
    for n_eas in n_eas_values:
        params = PopulationParams(n_eas=n_eas)
        results, _ = run_simulation(params, verbose=False)
        all_results[n_eas] = results
    return all_results


def sweep_heritability(h2_values=[0.1, 0.2, 0.3, 0.5, 0.7]):
    """Sweep over heritability values."""
    all_results = {}
    for h2 in h2_values:
        params = PopulationParams(h2=h2)
        results, _ = run_simulation(params, verbose=False)
        all_results[h2] = results
    return all_results


# ============================================================
# 14. Visualization Functions
# ============================================================

def plot_method_comparison(results, data, save_path=None):
    """Bar plot comparing all methods on AUC and R²."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    methods = list(results.keys())
    aucs = [results[m].get('AUC', 0) for m in methods]
    r2s = [results[m].get('R2_liability', 0) for m in methods]

    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']

    axes[0].barh(methods, aucs, color=colors, edgecolor='black', linewidth=0.5)
    axes[0].set_xlabel('AUC', fontsize=12)
    axes[0].set_title('Prediction Accuracy (AUC)', fontsize=13, fontweight='bold')
    axes[0].set_xlim(0.45, max(aucs) + 0.05)
    axes[0].axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)

    axes[1].barh(methods, r2s, color=colors, edgecolor='black', linewidth=0.5)
    axes[1].set_xlabel('R² (liability scale)', fontsize=12)
    axes[1].set_title('Variance Explained (R²)', fontsize=13, fontweight='bold')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_prs_distribution(data, save_path=None):
    """PRS distribution by case/control status for each method."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()

    prs_keys = ['prs_direct', 'prs_target', 'prs_bayes',
                'prs_meta', 'prs_lanc', 'prs_combined']
    titles = ['Direct Transfer', 'Target Pop', 'Bayesian LD',
              'Meta-Analysis', 'Local Ancestry', 'Combined (Proposed)']
    colors_case = '#e74c3c'
    colors_ctrl = '#3498db'

    y = data['y_eas_bin']
    for i, (key, title) in enumerate(zip(prs_keys, titles)):
        prs = data[key]
        prs_cases = prs[y == 1]
        prs_controls = prs[y == 0]

        axes[i].hist(prs_controls, bins=30, alpha=0.6, color=colors_ctrl,
                     label='Controls', density=True)
        axes[i].hist(prs_cases, bins=30, alpha=0.6, color=colors_case,
                     label='Cases', density=True)
        axes[i].set_title(title, fontsize=11, fontweight='bold')
        axes[i].legend(fontsize=8)
        axes[i].set_xlabel('PRS')
        axes[i].set_ylabel('Density')

    plt.suptitle('PRS Distribution by Disease Status (T2D, EAS Population)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_fst_sweep(fst_results, save_path=None):
    """Plot method performance across Fst values."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fst_vals = sorted(fst_results.keys())
    methods = list(fst_results[fst_vals[0]].keys())
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    markers = ['o', 's', '^', 'D', 'v', 'P']

    for i, method in enumerate(methods):
        aucs = [fst_results[f][method].get('AUC', 0.5) for f in fst_vals]
        r2s = [fst_results[f][method].get('R2_liability', 0) for f in fst_vals]
        axes[0].plot(fst_vals, aucs, marker=markers[i], color=colors[i],
                     label=method, linewidth=2, markersize=6)
        axes[1].plot(fst_vals, r2s, marker=markers[i], color=colors[i],
                     label=method, linewidth=2, markersize=6)

    axes[0].set_xlabel('Fst (Population Divergence)', fontsize=12)
    axes[0].set_ylabel('AUC', fontsize=12)
    axes[0].set_title('AUC vs Population Divergence', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=7, loc='lower left')
    axes[0].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    axes[1].set_xlabel('Fst (Population Divergence)', fontsize=12)
    axes[1].set_ylabel('R² (liability scale)', fontsize=12)
    axes[1].set_title('R² vs Population Divergence', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=7, loc='lower left')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_sample_size_sweep(ss_results, save_path=None):
    """Plot method performance across target sample sizes."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ss_vals = sorted(ss_results.keys())
    methods = list(ss_results[ss_vals[0]].keys())
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    markers = ['o', 's', '^', 'D', 'v', 'P']

    for i, method in enumerate(methods):
        aucs = [ss_results[s][method].get('AUC', 0.5) for s in ss_vals]
        r2s = [ss_results[s][method].get('R2_liability', 0) for s in ss_vals]
        axes[0].plot(ss_vals, aucs, marker=markers[i], color=colors[i],
                     label=method, linewidth=2, markersize=6)
        axes[1].plot(ss_vals, r2s, marker=markers[i], color=colors[i],
                     label=method, linewidth=2, markersize=6)

    axes[0].set_xlabel('Target Population Sample Size (N_EAS)', fontsize=12)
    axes[0].set_ylabel('AUC', fontsize=12)
    axes[0].set_title('AUC vs Sample Size', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=7)

    axes[1].set_xlabel('Target Population Sample Size (N_EAS)', fontsize=12)
    axes[1].set_ylabel('R² (liability scale)', fontsize=12)
    axes[1].set_title('R² vs Sample Size', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_ld_comparison(R_eur, R_eas, save_path=None):
    """Visualize LD structure differences between populations."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    n_show = min(50, R_eur.shape[0])

    im0 = axes[0].imshow(R_eur[:n_show, :n_show], cmap='RdBu_r', vmin=-1, vmax=1)
    axes[0].set_title('EUR LD Matrix', fontsize=12, fontweight='bold')
    plt.colorbar(im0, ax=axes[0], fraction=0.046)

    im1 = axes[1].imshow(R_eas[:n_show, :n_show], cmap='RdBu_r', vmin=-1, vmax=1)
    axes[1].set_title('EAS LD Matrix', fontsize=12, fontweight='bold')
    plt.colorbar(im1, ax=axes[1], fraction=0.046)

    diff = R_eur[:n_show, :n_show] - R_eas[:n_show, :n_show]
    im2 = axes[2].imshow(diff, cmap='RdBu_r', vmin=-0.5, vmax=0.5)
    axes[2].set_title('LD Difference (EUR - EAS)', fontsize=12, fontweight='bold')
    plt.colorbar(im2, ax=axes[2], fraction=0.046)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_effect_size_comparison(data, save_path=None):
    """Compare estimated effect sizes across populations."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    beta_eur = data['beta_hat_eur']
    beta_eas = data['beta_hat_eas']
    causal = data['causal_idx']

    mask_causal = np.zeros(len(beta_eur), dtype=bool)
    mask_causal[causal] = True

    # EUR vs EAS effect sizes
    axes[0].scatter(beta_eur[~mask_causal], beta_eas[~mask_causal],
                    alpha=0.3, s=10, color='gray', label='Non-causal')
    axes[0].scatter(beta_eur[mask_causal], beta_eas[mask_causal],
                    alpha=0.8, s=30, color='red', label='Causal')
    lim = max(abs(beta_eur).max(), abs(beta_eas).max()) * 1.1
    axes[0].plot([-lim, lim], [-lim, lim], 'k--', alpha=0.5)
    axes[0].set_xlabel('EUR Effect Size (β̂)', fontsize=11)
    axes[0].set_ylabel('EAS Effect Size (β̂)', fontsize=11)
    axes[0].set_title('Effect Size Comparison', fontsize=12, fontweight='bold')
    axes[0].legend(fontsize=9)
    r = np.corrcoef(beta_eur, beta_eas)[0, 1]
    axes[0].text(0.05, 0.95, f'r = {r:.3f}', transform=axes[0].transAxes,
                 fontsize=11, va='top')

    # True vs estimated
    beta_true = data['beta_true']
    axes[1].scatter(beta_true[~mask_causal], beta_eur[~mask_causal],
                    alpha=0.3, s=10, color='gray', label='Non-causal')
    axes[1].scatter(beta_true[mask_causal], beta_eur[mask_causal],
                    alpha=0.8, s=30, color='blue', label='Causal')
    axes[1].plot([-lim, lim], [-lim, lim], 'k--', alpha=0.5)
    axes[1].set_xlabel('True Effect Size (β)', fontsize=11)
    axes[1].set_ylabel('EUR Estimated Effect (β̂)', fontsize=11)
    axes[1].set_title('True vs Estimated Effects (EUR)', fontsize=12, fontweight='bold')
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_heritability_sweep(h2_results, save_path=None):
    """Plot method performance across heritability values."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    h2_vals = sorted(h2_results.keys())
    methods = list(h2_results[h2_vals[0]].keys())
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    markers = ['o', 's', '^', 'D', 'v', 'P']

    for i, method in enumerate(methods):
        aucs = [h2_results[h][method].get('AUC', 0.5) for h in h2_vals]
        r2s = [h2_results[h][method].get('R2_liability', 0) for h in h2_vals]
        axes[0].plot(h2_vals, aucs, marker=markers[i], color=colors[i],
                     label=method, linewidth=2, markersize=6)
        axes[1].plot(h2_vals, r2s, marker=markers[i], color=colors[i],
                     label=method, linewidth=2, markersize=6)

    axes[0].set_xlabel('Heritability (h²)', fontsize=12)
    axes[0].set_ylabel('AUC', fontsize=12)
    axes[0].set_title('AUC vs Heritability', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=7)
    axes[0].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    axes[1].set_xlabel('Heritability (h²)', fontsize=12)
    axes[1].set_ylabel('R² (liability scale)', fontsize=12)
    axes[1].set_title('R² vs Heritability', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_allele_freq_divergence(p_eur, p_eas, fst, save_path=None):
    """Plot allele frequency divergence between populations."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.scatter(p_eur, p_eas, alpha=0.4, s=15, color='#2c3e50')
    ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=1.5)
    ax.set_xlabel('EUR Allele Frequency', fontsize=12)
    ax.set_ylabel('EAS Allele Frequency', fontsize=12)
    ax.set_title(f'Allele Frequency Divergence (Fst={fst})',
                 fontsize=13, fontweight='bold')
    r = np.corrcoef(p_eur, p_eas)[0, 1]
    ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes,
            fontsize=12, va='top')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# 15. Run All Experiments
# ============================================================

def main():
    print("=" * 70)
    print("Cross-Ancestry PRS Transferability Simulation")
    print("Case Study: Type 2 Diabetes (EUR → EAS)")
    print("=" * 70)

    # --- Main simulation ---
    print("\n[1/5] Running main simulation (Fst=0.1, N_EUR=10000, N_EAS=5000)...")
    params = PopulationParams()
    results_main, data_main = run_simulation(params)

    # --- Figures for main simulation ---
    print("\n[2/5] Generating figures...")
    plot_method_comparison(results_main, data_main,
                          os.path.join(FIGURE_DIR, 'method_comparison.png'))
    plot_prs_distribution(data_main,
                          os.path.join(FIGURE_DIR, 'prs_distributions.png'))
    plot_ld_comparison(data_main['R_eur'], data_main['R_eas'],
                       os.path.join(FIGURE_DIR, 'ld_comparison.png'))
    plot_effect_size_comparison(data_main,
                                os.path.join(FIGURE_DIR, 'effect_sizes.png'))
    plot_allele_freq_divergence(data_main['p_eur'], data_main['p_eas'], 0.1,
                                os.path.join(FIGURE_DIR, 'allele_freq_divergence.png'))

    # --- Parameter sweeps ---
    print("\n[3/5] Fst sweep...")
    fst_results = sweep_fst()
    plot_fst_sweep(fst_results, os.path.join(FIGURE_DIR, 'fst_sweep.png'))

    print("\n[4/5] Sample size sweep...")
    ss_results = sweep_sample_size()
    plot_sample_size_sweep(ss_results,
                           os.path.join(FIGURE_DIR, 'sample_size_sweep.png'))

    print("\n[5/5] Heritability sweep...")
    h2_results = sweep_heritability()
    plot_heritability_sweep(h2_results,
                            os.path.join(FIGURE_DIR, 'heritability_sweep.png'))

    # --- Summary table ---
    print("\n" + "=" * 70)
    print("SUMMARY OF RESULTS")
    print("=" * 70)
    summary_data = []
    for method, metrics in results_main.items():
        summary_data.append({
            'Method': method,
            'AUC': f"{metrics.get('AUC', 0):.4f}",
            'R2_liability': f"{metrics.get('R2_liability', 0):.4f}",
            'R2_observed': f"{metrics.get('R2_observed', 0):.4f}",
        })
    df_summary = pd.DataFrame(summary_data)
    print(df_summary.to_string(index=False))

    # Save results to JSON
    results_json = {}
    results_json['main'] = {m: {k: float(v) for k, v in metrics.items()}
                            for m, metrics in results_main.items()}
    results_json['fst_sweep'] = {
        str(fst): {m: {k: float(v) for k, v in metrics.items()}
                   for m, metrics in res.items()}
        for fst, res in fst_results.items()
    }
    results_json['sample_size_sweep'] = {
        str(ss): {m: {k: float(v) for k, v in metrics.items()}
                  for m, metrics in res.items()}
        for ss, res in ss_results.items()
    }
    results_json['heritability_sweep'] = {
        str(h2): {m: {k: float(v) for k, v in metrics.items()}
                  for m, metrics in res.items()}
        for h2, res in h2_results.items()
    }

    with open('simulation_results.json', 'w') as f:
        json.dump(results_json, f, indent=2)

    print(f"\nResults saved to simulation_results.json")
    print(f"Figures saved to {FIGURE_DIR}/")
    print("Done!")

    return results_json


if __name__ == '__main__':
    results = main()
