#!/usr/bin/env python3
"""
PRS Transferability Framework: Improving Cross-Ethnic Polygenic Risk Score Transfer
From UK Biobank (European) to BioBank Japan (Japanese)

Implements:
1. PRS transfer problem formulation
2. Bayesian LD-aware correction
3. Multi-ethnic meta-analysis SNP effect re-estimation
4. Local ancestry-informed PRS correction
5. Simulation experiments (true effects, Fst, sample sizes)
6. Type 2 diabetes case study
"""

import numpy as np
import pandas as pd
from scipy import stats, linalg
from scipy.optimize import minimize
from sklearn.metrics import roc_auc_score, r2_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# 1. Population Genetics Simulation Engine
# ============================================================

class PopulationSimulator:
    """Simulate genotype data for two populations with controlled Fst."""
    
    def __init__(self, n_snps, n_causal, fst, n_eur, n_eas, heritability=0.5):
        self.n_snps = n_snps
        self.n_causal = n_causal
        self.fst = fst
        self.n_eur = n_eur
        self.n_eas = n_eas
        self.h2 = heritability
        self.causal_idx = np.sort(np.random.choice(n_snps, n_causal, replace=False))
        
    def _simulate_allele_freqs(self):
        """Simulate ancestral and population-specific allele frequencies using Balding-Nichols model."""
        p_anc = np.random.uniform(0.05, 0.95, self.n_snps)
        p_eur = np.zeros(self.n_snps)
        p_eas = np.zeros(self.n_snps)
        
        for i in range(self.n_snps):
            a = p_anc[i] * (1 - self.fst) / self.fst
            b = (1 - p_anc[i]) * (1 - self.fst) / self.fst
            if a > 0 and b > 0:
                p_eur[i] = np.clip(np.random.beta(a, b), 0.01, 0.99)
                p_eas[i] = np.clip(np.random.beta(a, b), 0.01, 0.99)
            else:
                p_eur[i] = p_anc[i]
                p_eas[i] = p_anc[i]
        
        return p_anc, p_eur, p_eas
    
    def _simulate_ld_matrix(self, n_snps, decay_rate=0.1, block_size=50):
        """Simulate block-diagonal LD matrix with exponential decay."""
        R = np.eye(n_snps)
        n_blocks = n_snps // block_size + 1
        
        for block in range(n_blocks):
            start = block * block_size
            end = min((block + 1) * block_size, n_snps)
            for i in range(start, end):
                for j in range(i + 1, end):
                    r = np.exp(-decay_rate * abs(i - j)) * np.random.uniform(0.5, 1.0)
                    R[i, j] = r
                    R[j, i] = r
        
        # Ensure positive definiteness
        eigvals = np.linalg.eigvalsh(R)
        if eigvals.min() < 1e-6:
            R += (abs(eigvals.min()) + 1e-4) * np.eye(n_snps)
            D = np.sqrt(np.diag(1.0 / np.diag(R)))
            R = D @ R @ D
        
        return R
    
    def _simulate_genotypes(self, n_samples, allele_freqs, ld_matrix):
        """Simulate genotypes with LD structure using Cholesky decomposition."""
        L = np.linalg.cholesky(ld_matrix)
        Z = np.random.randn(n_samples, self.n_snps)
        correlated_Z = Z @ L.T
        
        # Convert to genotypes using allele frequencies as thresholds
        genotypes = np.zeros((n_samples, self.n_snps), dtype=np.float64)
        for j in range(self.n_snps):
            p = allele_freqs[j]
            thresholds = stats.norm.ppf([p**2, p**2 + 2*p*(1-p)])
            genotypes[:, j] = np.where(correlated_Z[:, j] < thresholds[0], 0,
                                       np.where(correlated_Z[:, j] < thresholds[1], 1, 2))
        
        return genotypes
    
    def _simulate_true_effects(self):
        """Simulate true causal effect sizes with possible population-shared and specific effects."""
        beta_true = np.zeros(self.n_snps)
        # Shared causal effects (most effects are shared)
        beta_true[self.causal_idx] = np.random.normal(0, np.sqrt(self.h2 / self.n_causal), self.n_causal)
        return beta_true
    
    def simulate(self):
        """Run complete simulation."""
        p_anc, p_eur, p_eas = self._simulate_allele_freqs()
        
        # Different LD structures for each population
        ld_eur = self._simulate_ld_matrix(self.n_snps, decay_rate=0.08)
        ld_eas = self._simulate_ld_matrix(self.n_snps, decay_rate=0.12)
        
        geno_eur = self._simulate_genotypes(self.n_eur, p_eur, ld_eur)
        geno_eas = self._simulate_genotypes(self.n_eas, p_eas, ld_eas)
        
        # Standardize genotypes
        for j in range(self.n_snps):
            for geno, p in [(geno_eur, p_eur), (geno_eas, p_eas)]:
                mu = 2 * p[j]
                sd = np.sqrt(2 * p[j] * (1 - p[j]))
                if sd > 0:
                    geno[:, j] = (geno[:, j] - mu) / sd
        
        beta_true = self._simulate_true_effects()
        
        # Generate phenotypes
        g_eur = geno_eur @ beta_true
        g_eas = geno_eas @ beta_true
        
        var_g_eur = np.var(g_eur) if np.var(g_eur) > 0 else 1.0
        var_g_eas = np.var(g_eas) if np.var(g_eas) > 0 else 1.0
        
        noise_var_eur = var_g_eur * (1 - self.h2) / self.h2 if self.h2 > 0 else 1.0
        noise_var_eas = var_g_eas * (1 - self.h2) / self.h2 if self.h2 > 0 else 1.0
        
        y_eur = g_eur + np.random.normal(0, np.sqrt(noise_var_eur), self.n_eur)
        y_eas = g_eas + np.random.normal(0, np.sqrt(noise_var_eas), self.n_eas)
        
        return {
            'geno_eur': geno_eur, 'geno_eas': geno_eas,
            'y_eur': y_eur, 'y_eas': y_eas,
            'beta_true': beta_true,
            'p_eur': p_eur, 'p_eas': p_eas, 'p_anc': p_anc,
            'ld_eur': ld_eur, 'ld_eas': ld_eas,
            'causal_idx': self.causal_idx
        }


# ============================================================
# 2. PRS Methods
# ============================================================

class StandardPRS:
    """Standard PRS: direct transfer of EUR GWAS effects."""
    
    def __init__(self):
        self.name = "Standard PRS (EUR→EAS direct)"
    
    def fit(self, geno_eur, y_eur, **kwargs):
        """Estimate marginal effects from EUR GWAS."""
        n, p = geno_eur.shape
        self.beta_hat = np.zeros(p)
        self.se = np.zeros(p)
        self.pvals = np.ones(p)
        
        for j in range(p):
            x = geno_eur[:, j]
            if np.std(x) == 0:
                continue
            slope, intercept, r, pval, se = stats.linregress(x, y_eur)
            self.beta_hat[j] = slope
            self.se[j] = se
            self.pvals[j] = pval
        
        return self
    
    def predict(self, geno_target, p_threshold=5e-2):
        """Predict PRS in target population."""
        mask = self.pvals < p_threshold
        return geno_target[:, mask] @ self.beta_hat[mask]


class BayesianLDCorrectedPRS:
    """Bayesian PRS with LD correction across populations.
    
    Model: beta_target ~ N(C * beta_source, sigma^2 * (R_target)^{-1})
    where C is LD-ratio correction matrix.
    """
    
    def __init__(self, prior_variance=0.01):
        self.name = "Bayesian LD-Corrected PRS"
        self.prior_var = prior_variance
    
    def fit(self, geno_eur, y_eur, ld_eur, ld_eas, **kwargs):
        """Fit Bayesian model with LD correction."""
        n, p = geno_eur.shape
        
        # Step 1: Get marginal EUR effects
        beta_marginal = np.zeros(p)
        for j in range(p):
            x = geno_eur[:, j]
            if np.std(x) > 0:
                beta_marginal[j] = np.dot(x, y_eur) / np.dot(x, x)
        
        # Step 2: LD correction - estimate joint effects
        # R_eur * beta_joint = beta_marginal (approximately)
        reg_eur = ld_eur + 0.1 * np.eye(p)
        beta_joint_eur = np.linalg.solve(reg_eur, beta_marginal)
        
        # Step 3: Bayesian shrinkage with target LD
        # Posterior: beta_target = (R_eas + (1/prior_var)*I)^{-1} * R_eas * beta_joint_eur
        precision = ld_eas + (1.0 / self.prior_var) * np.eye(p)
        self.beta_corrected = np.linalg.solve(precision, ld_eas @ beta_joint_eur)
        
        # Estimate posterior variance
        self.posterior_var = np.diag(np.linalg.inv(precision)) * np.var(y_eur - geno_eur @ beta_joint_eur)
        
        return self
    
    def predict(self, geno_target, **kwargs):
        return geno_target @ self.beta_corrected


class MultiEthnicMetaAnalysisPRS:
    """Multi-ethnic meta-analysis with RE2 random-effects model.
    
    beta_meta = sum(w_k * beta_k) / sum(w_k)
    where w_k = 1 / (se_k^2 + tau^2), tau^2 estimated via DerSimonian-Laird.
    """
    
    def __init__(self, lambda_shrink=0.5):
        self.name = "Multi-Ethnic Meta-Analysis PRS"
        self.lambda_shrink = lambda_shrink
    
    def fit(self, geno_eur, y_eur, geno_eas_train=None, y_eas_train=None, **kwargs):
        """Fit using multi-ethnic meta-analysis."""
        n_eur, p = geno_eur.shape
        
        # EUR GWAS
        beta_eur = np.zeros(p)
        se_eur = np.zeros(p)
        for j in range(p):
            x = geno_eur[:, j]
            if np.std(x) > 0:
                slope, _, _, _, se = stats.linregress(x, y_eur)
                beta_eur[j] = slope
                se_eur[j] = se if se > 0 else 1.0
            else:
                se_eur[j] = 1.0
        
        # EAS GWAS (smaller sample)
        beta_eas = np.zeros(p)
        se_eas = np.zeros(p)
        if geno_eas_train is not None:
            for j in range(p):
                x = geno_eas_train[:, j]
                if np.std(x) > 0:
                    slope, _, _, _, se = stats.linregress(x, y_eas_train)
                    beta_eas[j] = slope
                    se_eas[j] = se if se > 0 else 1.0
                else:
                    se_eas[j] = 1.0
        else:
            beta_eas = beta_eur * 0.8 + np.random.normal(0, 0.01, p)
            se_eas = se_eur * 2.0
        
        # DerSimonian-Laird random-effects meta-analysis per SNP
        self.beta_meta = np.zeros(p)
        for j in range(p):
            betas = np.array([beta_eur[j], beta_eas[j]])
            ses = np.array([se_eur[j], se_eas[j]])
            variances = ses**2
            
            # Fixed-effect estimate
            w_fe = 1.0 / variances
            beta_fe = np.sum(w_fe * betas) / np.sum(w_fe)
            
            # Cochran's Q
            Q = np.sum(w_fe * (betas - beta_fe)**2)
            df = len(betas) - 1
            
            # Tau^2 (between-study variance)
            C = np.sum(w_fe) - np.sum(w_fe**2) / np.sum(w_fe)
            tau2 = max(0, (Q - df) / C) if C > 0 else 0
            
            # Random-effects weights
            w_re = 1.0 / (variances + tau2)
            self.beta_meta[j] = np.sum(w_re * betas) / np.sum(w_re)
        
        return self
    
    def predict(self, geno_target, **kwargs):
        return geno_target @ self.beta_meta


class LocalAncestryPRS:
    """PRS with local ancestry correction.
    
    PRS_corrected(i) = sum_j beta_adj(j) * G(i,j)
    where beta_adj(j) = alpha * beta_eur(j) + (1-alpha) * beta_eas(j)
    and alpha depends on local ancestry at locus j.
    """
    
    def __init__(self):
        self.name = "Local Ancestry-Corrected PRS"
    
    def _estimate_local_ancestry(self, geno_target, p_eur, p_eas):
        """Estimate local ancestry proportions per SNP using allele frequency similarity."""
        n, p_snps = geno_target.shape
        ancestry = np.zeros((n, p_snps))
        
        for j in range(p_snps):
            freq_obs = np.mean(geno_target[:, j] > 0)
            d_eur = abs(freq_obs - p_eur[j])
            d_eas = abs(freq_obs - p_eas[j])
            total = d_eur + d_eas
            if total > 0:
                # Proportion of EAS ancestry (closer to EAS freq = higher EAS ancestry)
                ancestry[:, j] = d_eur / total
            else:
                ancestry[:, j] = 0.5
        
        return ancestry
    
    def fit(self, geno_eur, y_eur, geno_eas_train=None, y_eas_train=None,
            p_eur=None, p_eas=None, **kwargs):
        """Fit ancestry-corrected model."""
        n, p = geno_eur.shape
        
        # EUR effects
        self.beta_eur = np.zeros(p)
        for j in range(p):
            x = geno_eur[:, j]
            if np.std(x) > 0:
                self.beta_eur[j] = np.dot(x, y_eur) / np.dot(x, x)
        
        # EAS effects
        self.beta_eas = np.zeros(p)
        if geno_eas_train is not None:
            for j in range(p):
                x = geno_eas_train[:, j]
                if np.std(x) > 0:
                    self.beta_eas[j] = np.dot(x, y_eas_train) / np.dot(x, x)
        else:
            self.beta_eas = self.beta_eur * 0.7
        
        self.p_eur = p_eur
        self.p_eas = p_eas
        
        return self
    
    def predict(self, geno_target, **kwargs):
        n, p = geno_target.shape
        ancestry = self._estimate_local_ancestry(geno_target, self.p_eur, self.p_eas)
        
        # Ancestry-weighted effect sizes
        prs = np.zeros(n)
        for i in range(n):
            beta_adj = ancestry[i] * self.beta_eas + (1 - ancestry[i]) * self.beta_eur
            prs[i] = np.dot(geno_target[i], beta_adj)
        
        return prs


class PenalizedTransferPRS:
    """Penalized regression transfer learning PRS.
    
    Minimizes: ||y_target - X_target * beta||^2 + lambda1 * ||beta||^2 + lambda2 * ||beta - beta_source||^2
    """
    
    def __init__(self, lambda1=0.1, lambda2=1.0):
        self.name = "Penalized Transfer PRS"
        self.lambda1 = lambda1
        self.lambda2 = lambda2
    
    def fit(self, geno_eur, y_eur, geno_eas_train=None, y_eas_train=None, **kwargs):
        """Fit penalized transfer model."""
        n_eur, p = geno_eur.shape
        
        # Source estimates
        self.beta_source = np.zeros(p)
        for j in range(p):
            x = geno_eur[:, j]
            if np.std(x) > 0:
                self.beta_source[j] = np.dot(x, y_eur) / np.dot(x, x)
        
        if geno_eas_train is not None and y_eas_train is not None:
            n_eas = geno_eas_train.shape[0]
            # Closed-form: beta = (X'X + lambda1*I + lambda2*I)^{-1} (X'y + lambda2*beta_source)
            XtX = geno_eas_train.T @ geno_eas_train
            Xty = geno_eas_train.T @ y_eas_train
            A = XtX + (self.lambda1 + self.lambda2) * np.eye(p)
            b = Xty + self.lambda2 * self.beta_source
            self.beta_transfer = np.linalg.solve(A, b)
        else:
            self.beta_transfer = self.beta_source
        
        return self
    
    def predict(self, geno_target, **kwargs):
        return geno_target @ self.beta_transfer


# ============================================================
# 3. Evaluation Metrics
# ============================================================

def evaluate_prs(y_true, prs_scores, is_binary=False):
    """Evaluate PRS performance."""
    results = {}
    
    # R-squared (partial)
    correlation = np.corrcoef(y_true, prs_scores)[0, 1]
    results['r2'] = correlation**2 if not np.isnan(correlation) else 0.0
    results['correlation'] = correlation if not np.isnan(correlation) else 0.0
    
    # Regression slope
    if np.std(prs_scores) > 0:
        slope, intercept, _, pval, _ = stats.linregress(prs_scores, y_true)
        results['slope'] = slope
        results['p_value'] = pval
    else:
        results['slope'] = 0.0
        results['p_value'] = 1.0
    
    # AUC for binary outcomes
    if is_binary:
        try:
            results['auc'] = roc_auc_score(y_true, prs_scores)
        except:
            results['auc'] = 0.5
    
    # Decile analysis
    if np.std(prs_scores) > 0:
        deciles = pd.qcut(prs_scores, 10, labels=False, duplicates='drop')
        top_decile = y_true[deciles == deciles.max()]
        bottom_decile = y_true[deciles == deciles.min()]
        results['top_bottom_ratio'] = np.mean(top_decile) / np.mean(bottom_decile) if np.mean(bottom_decile) != 0 else np.nan
    
    return results


# ============================================================
# 4. Main Simulation Experiments
# ============================================================

def run_experiment(n_snps=200, n_causal=20, fst=0.1, n_eur=5000, n_eas=1000,
                   h2=0.5, eas_train_frac=0.3, seed=42):
    """Run a single simulation experiment."""
    np.random.seed(seed)
    
    sim = PopulationSimulator(n_snps, n_causal, fst, n_eur, n_eas, h2)
    data = sim.simulate()
    
    # Split EAS into train/test
    n_eas_train = int(n_eas * eas_train_frac)
    idx = np.random.permutation(n_eas)
    eas_train_idx = idx[:n_eas_train]
    eas_test_idx = idx[n_eas_train:]
    
    geno_eas_train = data['geno_eas'][eas_train_idx]
    y_eas_train = data['y_eas'][eas_train_idx]
    geno_eas_test = data['geno_eas'][eas_test_idx]
    y_eas_test = data['y_eas'][eas_test_idx]
    
    methods = {
        'Standard': StandardPRS(),
        'Bayesian LD': BayesianLDCorrectedPRS(prior_variance=0.01),
        'Meta-Analysis': MultiEthnicMetaAnalysisPRS(),
        'Local Ancestry': LocalAncestryPRS(),
        'Penalized Transfer': PenalizedTransferPRS(lambda1=0.1, lambda2=1.0),
    }
    
    results = {}
    for name, method in methods.items():
        fit_kwargs = {
            'geno_eur': data['geno_eur'], 'y_eur': data['y_eur'],
            'ld_eur': data['ld_eur'], 'ld_eas': data['ld_eas'],
            'geno_eas_train': geno_eas_train, 'y_eas_train': y_eas_train,
            'p_eur': data['p_eur'], 'p_eas': data['p_eas'],
        }
        method.fit(**fit_kwargs)
        prs = method.predict(geno_eas_test)
        metrics = evaluate_prs(y_eas_test, prs)
        metrics['method'] = name
        results[name] = metrics
    
    # EUR within-population performance (oracle)
    std = StandardPRS()
    std.fit(data['geno_eur'][:n_eur//2], data['y_eur'][:n_eur//2])
    prs_eur = std.predict(data['geno_eur'][n_eur//2:])
    eur_metrics = evaluate_prs(data['y_eur'][n_eur//2:], prs_eur)
    eur_metrics['method'] = 'EUR Within-Pop'
    results['EUR Within-Pop'] = eur_metrics
    
    return results, data


def run_fst_sweep():
    """Sweep over Fst values to study transferability decay."""
    fst_values = [0.01, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]
    all_results = []
    
    for fst in fst_values:
        print(f"  Fst = {fst:.2f}...")
        res, _ = run_experiment(fst=fst, seed=42)
        for method, metrics in res.items():
            metrics['fst'] = fst
            all_results.append(metrics)
    
    return pd.DataFrame(all_results)


def run_sample_size_sweep():
    """Sweep over EAS training sample sizes."""
    n_eas_values = [200, 500, 1000, 2000, 5000]
    all_results = []
    
    for n_eas in n_eas_values:
        print(f"  N_EAS = {n_eas}...")
        res, _ = run_experiment(n_eas=n_eas, eas_train_frac=0.3, seed=42)
        for method, metrics in res.items():
            metrics['n_eas'] = n_eas
            all_results.append(metrics)
    
    return pd.DataFrame(all_results)


def run_heritability_sweep():
    """Sweep over heritability values."""
    h2_values = [0.1, 0.2, 0.3, 0.5, 0.7]
    all_results = []
    
    for h2 in h2_values:
        print(f"  h2 = {h2:.1f}...")
        res, _ = run_experiment(h2=h2, seed=42)
        for method, metrics in res.items():
            metrics['h2'] = h2
            all_results.append(metrics)
    
    return pd.DataFrame(all_results)


# ============================================================
# 5. Type 2 Diabetes Case Study
# ============================================================

def t2d_case_study():
    """Simulate T2D-like PRS transfer scenario."""
    print("\n=== Type 2 Diabetes Case Study ===")
    
    # T2D parameters: ~400 known loci, h2_SNP ~0.20, Fst_EUR-EAS ~0.11
    n_snps = 300
    n_causal = 40
    h2 = 0.20
    fst = 0.11
    n_eur = 8000  # UK Biobank-scale (reduced for simulation)
    n_eas = 2000  # BioBank Japan-scale (reduced)
    
    np.random.seed(123)
    sim = PopulationSimulator(n_snps, n_causal, fst, n_eur, n_eas, h2)
    data = sim.simulate()
    
    # Convert to binary trait (T2D prevalence ~10% EUR, ~12% EAS)
    threshold_eur = np.percentile(data['y_eur'], 90)
    threshold_eas = np.percentile(data['y_eas'], 88)
    y_eur_binary = (data['y_eur'] > threshold_eur).astype(int)
    y_eas_binary = (data['y_eas'] > threshold_eas).astype(int)
    
    # Split EAS
    n_eas_train = 600
    idx = np.random.permutation(n_eas)
    eas_train_idx = idx[:n_eas_train]
    eas_test_idx = idx[n_eas_train:]
    
    methods = {
        'Standard': StandardPRS(),
        'Bayesian LD': BayesianLDCorrectedPRS(prior_variance=0.01),
        'Meta-Analysis': MultiEthnicMetaAnalysisPRS(),
        'Local Ancestry': LocalAncestryPRS(),
        'Penalized Transfer': PenalizedTransferPRS(lambda1=0.1, lambda2=1.0),
    }
    
    results = {}
    prs_dict = {}
    for name, method in methods.items():
        fit_kwargs = {
            'geno_eur': data['geno_eur'], 'y_eur': data['y_eur'],
            'ld_eur': data['ld_eur'], 'ld_eas': data['ld_eas'],
            'geno_eas_train': data['geno_eas'][eas_train_idx],
            'y_eas_train': data['y_eas'][eas_train_idx],
            'p_eur': data['p_eur'], 'p_eas': data['p_eas'],
        }
        method.fit(**fit_kwargs)
        prs = method.predict(data['geno_eas'][eas_test_idx])
        metrics = evaluate_prs(y_eas_binary[eas_test_idx], prs, is_binary=True)
        metrics['method'] = name
        results[name] = metrics
        prs_dict[name] = prs
    
    return results, prs_dict, y_eas_binary[eas_test_idx], data


# ============================================================
# 6. Visualization
# ============================================================

def plot_fst_sweep(df, save_path='figures/fst_sweep.png'):
    """Plot R² vs Fst for all methods."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    methods_to_plot = [m for m in df['method'].unique() if m != 'EUR Within-Pop']
    colors = sns.color_palette('viridis', len(methods_to_plot))
    markers = ['o', 's', 'D', '^', 'v']
    
    for i, method in enumerate(methods_to_plot):
        subset = df[df['method'] == method]
        ax.plot(subset['fst'], subset['r2'], marker=markers[i % len(markers)],
                color=colors[i], label=method, linewidth=2, markersize=8)
    
    # EUR within-pop reference
    eur = df[df['method'] == 'EUR Within-Pop']
    if len(eur) > 0:
        ax.axhline(y=eur['r2'].mean(), color='gray', linestyle='--',
                    label='EUR Within-Pop (mean)', alpha=0.7)
    
    ax.set_xlabel('Population Differentiation (Fst)', fontsize=13)
    ax.set_ylabel('Prediction R²', fontsize=13)
    ax.set_title('PRS Transferability vs Population Differentiation', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_sample_size_sweep(df, save_path='figures/sample_size_sweep.png'):
    """Plot R² vs EAS sample size."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    methods_to_plot = [m for m in df['method'].unique() if m != 'EUR Within-Pop']
    colors = sns.color_palette('cividis', len(methods_to_plot))
    markers = ['o', 's', 'D', '^', 'v']
    
    for i, method in enumerate(methods_to_plot):
        subset = df[df['method'] == method]
        ax.plot(subset['n_eas'], subset['r2'], marker=markers[i % len(markers)],
                color=colors[i], label=method, linewidth=2, markersize=8)
    
    ax.set_xlabel('EAS Sample Size', fontsize=13)
    ax.set_ylabel('Prediction R²', fontsize=13)
    ax.set_title('PRS Performance vs Target Population Sample Size', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_heritability_sweep(df, save_path='figures/heritability_sweep.png'):
    """Plot R² vs heritability."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    methods_to_plot = [m for m in df['method'].unique() if m != 'EUR Within-Pop']
    colors = sns.color_palette('viridis', len(methods_to_plot))
    markers = ['o', 's', 'D', '^', 'v']
    
    for i, method in enumerate(methods_to_plot):
        subset = df[df['method'] == method]
        ax.plot(subset['h2'], subset['r2'], marker=markers[i % len(markers)],
                color=colors[i], label=method, linewidth=2, markersize=8)
    
    ax.set_xlabel('Heritability (h²)', fontsize=13)
    ax.set_ylabel('Prediction R²', fontsize=13)
    ax.set_title('PRS Performance vs Trait Heritability', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_method_comparison(results, save_path='figures/method_comparison.png'):
    """Bar plot comparing methods."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    methods = list(results.keys())
    r2_vals = [results[m]['r2'] for m in methods]
    corr_vals = [results[m]['correlation'] for m in methods]
    
    colors = sns.color_palette('viridis', len(methods))
    
    axes[0].barh(methods, r2_vals, color=colors)
    axes[0].set_xlabel('R²', fontsize=13)
    axes[0].set_title('Prediction R² by Method', fontsize=14)
    axes[0].grid(True, alpha=0.3, axis='x')
    
    axes[1].barh(methods, corr_vals, color=colors)
    axes[1].set_xlabel('Correlation', fontsize=13)
    axes[1].set_title('Correlation by Method', fontsize=14)
    axes[1].grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_t2d_results(results, prs_dict, y_true, save_path='figures/t2d_case_study.png'):
    """Plot T2D case study results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # AUC comparison
    methods = list(results.keys())
    aucs = [results[m].get('auc', 0.5) for m in methods]
    colors = sns.color_palette('viridis', len(methods))
    
    axes[0].barh(methods, aucs, color=colors)
    axes[0].set_xlabel('AUC', fontsize=13)
    axes[0].set_title('T2D Prediction AUC', fontsize=14)
    axes[0].axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='Random')
    axes[0].grid(True, alpha=0.3, axis='x')
    axes[0].legend()
    
    # PRS distribution by case/control for best method
    best_method = max(results, key=lambda m: results[m].get('auc', 0))
    prs = prs_dict[best_method]
    cases = prs[y_true == 1]
    controls = prs[y_true == 0]
    
    axes[1].hist(controls, bins=30, alpha=0.6, label='Controls', color='steelblue', density=True)
    axes[1].hist(cases, bins=30, alpha=0.6, label='Cases', color='coral', density=True)
    axes[1].set_xlabel('PRS', fontsize=13)
    axes[1].set_ylabel('Density', fontsize=13)
    axes[1].set_title(f'PRS Distribution ({best_method})', fontsize=14)
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)
    
    # Odds ratio by PRS decile
    if np.std(prs) > 0:
        try:
            deciles = pd.qcut(prs, 10, labels=range(1, 11), duplicates='drop')
            or_per_decile = []
            decile_labels = []
            for d in sorted(np.unique(deciles)):
                mask = deciles == d
                n_case = np.sum(y_true[mask] == 1)
                n_ctrl = np.sum(y_true[mask] == 0)
                # Reference: median decile
                ref_mask = deciles == np.median(np.unique(deciles))
                ref_case = np.sum(y_true[ref_mask] == 1)
                ref_ctrl = np.sum(y_true[ref_mask] == 0)
                
                if n_ctrl > 0 and ref_case > 0 and ref_ctrl > 0:
                    OR = (n_case / n_ctrl) / (ref_case / ref_ctrl)
                    or_per_decile.append(OR)
                    decile_labels.append(str(d))
            
            axes[2].bar(decile_labels, or_per_decile, color=sns.color_palette('viridis', len(decile_labels)))
            axes[2].axhline(y=1.0, color='red', linestyle='--', alpha=0.5)
            axes[2].set_xlabel('PRS Decile', fontsize=13)
            axes[2].set_ylabel('Odds Ratio (vs median)', fontsize=13)
            axes[2].set_title('T2D Risk by PRS Decile', fontsize=14)
            axes[2].grid(True, alpha=0.3, axis='y')
        except Exception:
            axes[2].text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=axes[2].transAxes)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_allele_freq_comparison(data, save_path='figures/allele_freq_comparison.png'):
    """Plot allele frequency comparison between populations."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    axes[0].scatter(data['p_eur'], data['p_eas'], alpha=0.5, s=20, c='teal')
    axes[0].plot([0, 1], [0, 1], 'r--', alpha=0.5)
    axes[0].set_xlabel('EUR Allele Frequency', fontsize=13)
    axes[0].set_ylabel('EAS Allele Frequency', fontsize=13)
    axes[0].set_title('Allele Frequency Comparison', fontsize=14)
    axes[0].grid(True, alpha=0.3)
    
    # LD comparison (first 50x50 block)
    n_show = min(50, data['ld_eur'].shape[0])
    im = axes[1].imshow(data['ld_eur'][:n_show, :n_show] - data['ld_eas'][:n_show, :n_show],
                        cmap='RdBu_r', vmin=-0.5, vmax=0.5, aspect='auto')
    axes[1].set_title('LD Difference (EUR - EAS)', fontsize=14)
    axes[1].set_xlabel('SNP Index', fontsize=13)
    axes[1].set_ylabel('SNP Index', fontsize=13)
    plt.colorbar(im, ax=axes[1], label='LD Difference')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_effect_size_comparison(data, save_path='figures/effect_size_analysis.png'):
    """Plot effect size analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    beta = data['beta_true']
    causal = data['causal_idx']
    
    # Effect size distribution
    axes[0].hist(beta[beta != 0], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('True Effect Size', fontsize=13)
    axes[0].set_ylabel('Count', fontsize=13)
    axes[0].set_title('Distribution of Causal Effect Sizes', fontsize=14)
    axes[0].axvline(x=0, color='red', linestyle='--', alpha=0.5)
    axes[0].grid(True, alpha=0.3)
    
    # Manhattan-style plot
    positions = np.arange(len(beta))
    colors = ['gray'] * len(beta)
    for idx in causal:
        colors[idx] = 'red'
    
    axes[1].scatter(positions, np.abs(beta), c=colors, s=15, alpha=0.7)
    axes[1].set_xlabel('SNP Position', fontsize=13)
    axes[1].set_ylabel('|Effect Size|', fontsize=13)
    axes[1].set_title('Effect Size Manhattan Plot (Red = Causal)', fontsize=14)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


# ============================================================
# 7. Main Execution
# ============================================================

def main():
    print("=" * 60)
    print("PRS Cross-Ethnic Transferability Simulation Framework")
    print("=" * 60)
    
    # --- Experiment 1: Baseline comparison ---
    print("\n[1/5] Running baseline experiment...")
    base_results, base_data = run_experiment(
        n_snps=200, n_causal=20, fst=0.10, n_eur=5000, n_eas=1000, h2=0.5
    )
    
    print("\nBaseline Results (Fst=0.10, N_EUR=5000, N_EAS=1000, h²=0.5):")
    print(f"{'Method':<30} {'R²':>8} {'Corr':>8} {'Slope':>8}")
    print("-" * 56)
    for name, m in base_results.items():
        print(f"{name:<30} {m['r2']:8.4f} {m['correlation']:8.4f} {m['slope']:8.4f}")
    
    # --- Experiment 2: Fst sweep ---
    print("\n[2/5] Running Fst sweep...")
    fst_df = run_fst_sweep()
    
    # --- Experiment 3: Sample size sweep ---
    print("\n[3/5] Running sample size sweep...")
    ss_df = run_sample_size_sweep()
    
    # --- Experiment 4: Heritability sweep ---
    print("\n[4/5] Running heritability sweep...")
    h2_df = run_heritability_sweep()
    
    # --- Experiment 5: T2D case study ---
    print("\n[5/5] Running T2D case study...")
    t2d_results, t2d_prs, t2d_y, t2d_data = t2d_case_study()
    
    print("\nT2D Case Study Results:")
    print(f"{'Method':<30} {'AUC':>8} {'R²':>8}")
    print("-" * 48)
    for name, m in t2d_results.items():
        print(f"{name:<30} {m.get('auc', 'N/A'):>8.4f} {m['r2']:8.4f}")
    
    # --- Generate plots ---
    print("\n[Generating Figures]")
    plot_fst_sweep(fst_df)
    plot_sample_size_sweep(ss_df)
    plot_heritability_sweep(h2_df)
    plot_method_comparison(base_results)
    plot_t2d_results(t2d_results, t2d_prs, t2d_y)
    plot_allele_freq_comparison(base_data)
    plot_effect_size_comparison(base_data)
    
    # --- Save numeric results ---
    base_df = pd.DataFrame(base_results).T
    base_df.to_csv('results/baseline_results.csv', index=True)
    fst_df.to_csv('results/fst_sweep_results.csv', index=False)
    ss_df.to_csv('results/sample_size_sweep_results.csv', index=False)
    h2_df.to_csv('results/heritability_sweep_results.csv', index=False)
    t2d_df = pd.DataFrame(t2d_results).T
    t2d_df.to_csv('results/t2d_case_study_results.csv', index=True)
    
    print("\n✓ All results saved to results/")
    print("✓ All figures saved to figures/")
    
    return {
        'base_results': base_results,
        'fst_df': fst_df,
        'ss_df': ss_df,
        'h2_df': h2_df,
        't2d_results': t2d_results,
    }


if __name__ == '__main__':
    results = main()
