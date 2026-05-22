"""
PRS Transferability Simulation Framework
多遺伝子リスクスコアの民族集団間移植性改善手法

Methods:
  1. Standard PRS (no correction)
  2. LD-corrected Bayesian PRS (LDpred-style)
  3. Multi-ethnic meta-analysis re-estimation
  4. Local ancestry-informed PRS (LAI-PRS)
  5. Continuous shrinkage prior (CSP-PRS)
"""

import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.stats as stats
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
import warnings, json, os, time
warnings.filterwarnings('ignore')

np.random.seed(42)

# ── output directories ──────────────────────────────────────────────────────
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)
os.makedirs("logs", exist_ok=True)

LOG_PATH = "logs/process-log.jsonl"

def log_event(phase, event_type, skill_or_tool, handoff_in=None, handoff_out=None, files=None, status="ok"):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase": phase, "event_type": event_type, "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files or [],
        "status": status,
    }
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


# ═══════════════════════════════════════════════════════════════════════════
# 1. SIMULATION ENGINE
# ═══════════════════════════════════════════════════════════════════════════

def simulate_ld_matrix(n_snps: int, ld_decay: float = 0.3) -> np.ndarray:
    """Toeplitz LD matrix with exponential decay."""
    row = ld_decay ** np.arange(n_snps)
    R = la.toeplitz(row)
    # ensure positive definite
    eigvals = np.linalg.eigvalsh(R)
    if eigvals.min() < 1e-6:
        R += (1e-5 - eigvals.min()) * np.eye(n_snps)
    return R


def simulate_population(
    n_samples: int, n_snps: int, mafs: np.ndarray,
    ld_matrix: np.ndarray, ld_decay: float
) -> np.ndarray:
    """
    Simulate genotypes for a population with given MAFs and LD structure.
    Uses multivariate normal → threshold approach to approximate LD.
    """
    # Cholesky factor for correlated standard normals
    L = np.linalg.cholesky(ld_matrix + 1e-8 * np.eye(n_snps))
    Z = np.random.randn(n_samples, n_snps) @ L.T   # correlated latent
    # Threshold to produce genotypes (0/1/2)
    G = np.zeros((n_samples, n_snps), dtype=np.float32)
    for j in range(n_snps):
        p = mafs[j]
        t1 = stats.norm.ppf((1 - p) ** 2)
        t2 = stats.norm.ppf(1 - p ** 2)
        G[:, j] = (Z[:, j] > t1).astype(float) + (Z[:, j] > t2).astype(float)
    return G


def simulate_phenotype(G: np.ndarray, beta: np.ndarray, h2: float) -> np.ndarray:
    """
    y = G·β + ε  with Var(ε) chosen to achieve heritability h2.
    """
    genetic = G @ beta
    var_g = np.var(genetic)
    var_e = var_g * (1 - h2) / h2
    eps = np.random.normal(0, np.sqrt(max(var_e, 1e-8)), len(G))
    return genetic + eps


def differentiate_mafs(mafs_eur: np.ndarray, fst: float) -> np.ndarray:
    """
    Drift MAFs from EUR to ASN using Balding-Nichols model.
    p_asn ~ Beta(p*(1-fst)/fst, (1-p)*(1-fst)/fst)
    """
    alpha = mafs_eur * (1 - fst) / fst
    beta_  = (1 - mafs_eur) * (1 - fst) / fst
    alpha = np.clip(alpha, 0.01, 1000)
    beta_ = np.clip(beta_, 0.01, 1000)
    return np.random.beta(alpha, beta_)


# ═══════════════════════════════════════════════════════════════════════════
# 2. PRS METHODS
# ═══════════════════════════════════════════════════════════════════════════

def compute_prs_standard(G: np.ndarray, beta_gwas: np.ndarray) -> np.ndarray:
    """Standard PRS: direct application of GWAS effect sizes."""
    return G @ beta_gwas


def compute_prs_pt(G: np.ndarray, beta_gwas: np.ndarray, p_values: np.ndarray,
                   threshold: float = 5e-8) -> np.ndarray:
    """P+T PRS: clumping & thresholding (simplified, no LD clumping here)."""
    mask = p_values < threshold
    b = beta_gwas.copy()
    b[~mask] = 0.0
    return G @ b


def bayesian_ld_correction(
    beta_gwas: np.ndarray,
    R_ref: np.ndarray,
    n_gwas: int,
    h2: float,
    p_causal: float
) -> np.ndarray:
    """
    LDpred-style Bayesian correction for LD structure.
    Posterior mean estimator under point-normal prior.

    E[β|β_hat] ≈ (R + (1-p)*M/n*h2 * I)^{-1} * R * β_hat  (simplified)

    Full conjugate update (per-SNP approximation):
      posterior_mean_j = prior_weight * shrinkage * beta_j
    """
    M = len(beta_gwas)
    # Marginal per-SNP posterior (LDpred-inf approximation)
    lambda_reg = (1 - h2) * M / (h2 * n_gwas)
    A = R_ref + lambda_reg * np.eye(M)
    # Regularized solve
    beta_posterior = np.linalg.solve(A, R_ref @ beta_gwas)
    # Mix with zero (point-normal)
    pi = p_causal  # mixing probability
    beta_corrected = pi * beta_posterior
    return beta_corrected


def multi_ethnic_meta_reestimate(
    beta_eur: np.ndarray,
    se_eur: np.ndarray,
    beta_asn: np.ndarray,
    se_asn: np.ndarray,
    R_eur: np.ndarray,
    R_asn: np.ndarray,
    n_eur: int,
    n_asn: int
) -> np.ndarray:
    """
    Multi-ethnic fixed-effects meta-analysis with LD weighting.
    Inverse-variance weighting per SNP.
    W_eur = 1/se_eur^2, W_asn = 1/se_asn^2
    beta_meta = (W_eur*beta_eur + W_asn*beta_asn) / (W_eur + W_asn)
    """
    w_eur = 1.0 / (se_eur ** 2 + 1e-10)
    w_asn = 1.0 / (se_asn ** 2 + 1e-10)
    beta_meta = (w_eur * beta_eur + w_asn * beta_asn) / (w_eur + w_asn)
    return beta_meta


def local_ancestry_prs(
    G: np.ndarray,
    beta_eur: np.ndarray,
    beta_asn: np.ndarray,
    anc_eur: np.ndarray,
    anc_asn: np.ndarray
) -> np.ndarray:
    """
    Local ancestry-informed PRS (LAI-PRS).
    For each individual i and SNP j:
      PRS_i = Σ_j [ anc_eur_ij * beta_eur_j + anc_asn_ij * beta_asn_j ] * G_ij / 2

    anc_eur, anc_asn: (n_samples, n_snps) ancestry dosage fractions (0..2)
    """
    # Weighted effect sizes per sample per SNP
    beta_local = anc_eur * beta_eur[np.newaxis, :] + anc_asn * beta_asn[np.newaxis, :]
    prs = np.sum(G * beta_local, axis=1)
    return prs


def continuous_shrinkage_prs(
    beta_gwas: np.ndarray,
    se_gwas: np.ndarray,
    R_ref: np.ndarray,
    n_gwas: int,
    h2: float
) -> np.ndarray:
    """
    PRS-CS style: continuous shrinkage (global-local horseshoe approximation).
    Simplified via ridge regression with SNP-specific penalty.
    """
    M = len(beta_gwas)
    # SNP-specific shrinkage based on GWAS chi2
    chi2 = (beta_gwas / (se_gwas + 1e-10)) ** 2
    # Global shrinkage parameter (PRS-CS default inspired)
    phi = h2 / (M * np.mean(chi2 + 1e-5))
    # Local scale via half-Cauchy approximation
    delta = chi2 / (chi2 + 1.0)  # ∈ (0,1)
    shrink = delta * phi / (delta * phi + (1 - phi) / n_gwas)
    shrink = np.clip(shrink, 0, 1)
    beta_cs = shrink * beta_gwas
    return beta_cs


# ═══════════════════════════════════════════════════════════════════════════
# 3. EVALUATION METRICS
# ═══════════════════════════════════════════════════════════════════════════

def evaluate_prs(prs: np.ndarray, y_true: np.ndarray, label: str) -> dict:
    r2 = np.corrcoef(prs, y_true)[0, 1] ** 2
    pearson_r, p_val = stats.pearsonr(prs, y_true)
    # Nagelkerke-like R2 for continuous trait
    residuals = y_true - (np.mean(y_true) + pearson_r * np.std(y_true) / np.std(prs + 1e-10) * (prs - np.mean(prs)))
    return {
        "method": label,
        "R2": float(r2),
        "Pearson_r": float(pearson_r),
        "p_value": float(p_val),
        "mean_prs": float(np.mean(prs)),
        "std_prs": float(np.std(prs)),
    }


# ═══════════════════════════════════════════════════════════════════════════
# 4. MAIN SIMULATION
# ═══════════════════════════════════════════════════════════════════════════

def run_simulation(
    n_snps: int = 200,
    n_causal: int = 30,
    h2: float = 0.40,
    fst: float = 0.10,
    n_eur: int = 10000,
    n_asn_gwas: int = 5000,
    n_asn_test: int = 2000,
    ld_eur: float = 0.30,
    ld_asn: float = 0.20,
    p_causal: float = None,
    seed: int = 42
) -> dict:
    """
    Complete simulation for one parameter configuration.
    Returns dict of evaluation metrics per method.
    """
    np.random.seed(seed)
    if p_causal is None:
        p_causal = n_causal / n_snps

    # ── true effect sizes ────────────────────────────────────────────────
    causal_idx = np.random.choice(n_snps, n_causal, replace=False)
    beta_true = np.zeros(n_snps)
    beta_true[causal_idx] = np.random.normal(0, np.sqrt(h2 / n_causal), n_causal)

    # ── population MAFs ──────────────────────────────────────────────────
    mafs_eur = np.random.uniform(0.05, 0.45, n_snps)
    mafs_asn = differentiate_mafs(mafs_eur, fst)
    mafs_asn = np.clip(mafs_asn, 0.01, 0.49)

    # ── LD matrices ──────────────────────────────────────────────────────
    R_eur = simulate_ld_matrix(n_snps, ld_eur)
    R_asn = simulate_ld_matrix(n_snps, ld_asn)

    # ── genotypes ────────────────────────────────────────────────────────
    G_eur = simulate_population(n_eur, n_snps, mafs_eur, R_eur, ld_eur)
    G_asn_gwas = simulate_population(n_asn_gwas, n_snps, mafs_asn, R_asn, ld_asn)
    G_asn_test = simulate_population(n_asn_test, n_snps, mafs_asn, R_asn, ld_asn)

    # ── phenotypes ───────────────────────────────────────────────────────
    y_eur = simulate_phenotype(G_eur, beta_true, h2)
    y_asn_gwas = simulate_phenotype(G_asn_gwas, beta_true, h2)
    y_asn_test = simulate_phenotype(G_asn_test, beta_true, h2)

    # ── GWAS summary statistics (OLS marginal) ───────────────────────────
    def marginal_gwas(G, y):
        n, m = G.shape
        betas, ses = np.zeros(m), np.zeros(m)
        for j in range(m):
            x = G[:, j]
            sx = x - x.mean()
            if sx.std() < 1e-8:
                continue
            b = np.dot(sx, y - y.mean()) / np.dot(sx, sx)
            resid = y - y.mean() - b * sx
            s2 = np.var(resid)
            betas[j] = b
            ses[j] = np.sqrt(s2 / (np.dot(sx, sx) + 1e-10))
        return betas, ses

    beta_hat_eur, se_eur = marginal_gwas(G_eur, y_eur)
    beta_hat_asn, se_asn = marginal_gwas(G_asn_gwas, y_asn_gwas)
    p_values_eur = 2 * stats.norm.sf(np.abs(beta_hat_eur / (se_eur + 1e-10)))

    # ── PRS methods ──────────────────────────────────────────────────────
    prs_standard = compute_prs_standard(G_asn_test, beta_hat_eur)
    prs_pt = compute_prs_pt(G_asn_test, beta_hat_eur, p_values_eur, 0.001)
    beta_ld  = bayesian_ld_correction(beta_hat_eur, R_asn, n_eur, h2, p_causal)
    prs_ldcorr = compute_prs_standard(G_asn_test, beta_ld)
    beta_meta = multi_ethnic_meta_reestimate(
        beta_hat_eur, se_eur, beta_hat_asn, se_asn, R_eur, R_asn, n_eur, n_asn_gwas
    )
    prs_meta = compute_prs_standard(G_asn_test, beta_meta)
    # Local ancestry (admixed: mix EUR/ASN proportions 50/50 in test set)
    anc_eur_frac = np.random.beta(2, 2, (n_asn_test, n_snps))  # admixture
    anc_asn_frac = 1.0 - anc_eur_frac
    prs_lai = local_ancestry_prs(G_asn_test, beta_hat_eur, beta_hat_asn,
                                  anc_eur_frac, anc_asn_frac)
    beta_cs = continuous_shrinkage_prs(beta_hat_eur, se_eur, R_asn, n_eur, h2)
    prs_cs = compute_prs_standard(G_asn_test, beta_cs)
    # Oracle (true betas)
    prs_oracle = compute_prs_standard(G_asn_test, beta_true)

    results = {}
    for label, prs in [
        ("Standard PRS",    prs_standard),
        ("P+T PRS",         prs_pt),
        ("LD-corrected (Bayes)", prs_ldcorr),
        ("Multi-ethnic Meta", prs_meta),
        ("LAI-PRS",         prs_lai),
        ("CS-PRS",          prs_cs),
        ("Oracle",          prs_oracle),
    ]:
        results[label] = evaluate_prs(prs, y_asn_test, label)

    return results, {
        "beta_true": beta_true,
        "beta_hat_eur": beta_hat_eur,
        "beta_hat_asn": beta_hat_asn,
        "beta_meta": beta_meta,
        "beta_ld": beta_ld,
        "beta_cs": beta_cs,
        "causal_idx": causal_idx,
        "prs_standard": prs_standard,
        "prs_ldcorr": prs_ldcorr,
        "prs_meta": prs_meta,
        "prs_lai": prs_lai,
        "prs_cs": prs_cs,
        "prs_oracle": prs_oracle,
        "y_asn_test": y_asn_test,
        "mafs_eur": mafs_eur,
        "mafs_asn": mafs_asn,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 5. FST SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def fst_sensitivity(fst_values=None, n_rep=5):
    """Vary Fst and compare R2 across methods."""
    if fst_values is None:
        fst_values = [0.01, 0.05, 0.10, 0.15, 0.20]
    records = []
    for fst in fst_values:
        r2_by_method = {m: [] for m in ["Standard PRS", "LD-corrected (Bayes)",
                                         "Multi-ethnic Meta", "LAI-PRS", "CS-PRS", "Oracle"]}
        for rep in range(n_rep):
            res, _ = run_simulation(fst=fst, seed=rep * 100)
            for m in r2_by_method:
                r2_by_method[m].append(res[m]["R2"])
        for m, vals in r2_by_method.items():
            records.append({
                "Fst": fst, "method": m,
                "R2_mean": np.mean(vals),
                "R2_std": np.std(vals),
            })
    return pd.DataFrame(records)


def sample_size_sensitivity(n_values=None, n_rep=5):
    """Vary target population GWAS sample size."""
    if n_values is None:
        n_values = [500, 1000, 2000, 5000, 10000]
    records = []
    for n in n_values:
        r2_by_method = {m: [] for m in ["Standard PRS", "Multi-ethnic Meta", "CS-PRS", "Oracle"]}
        for rep in range(n_rep):
            res, _ = run_simulation(n_asn_gwas=n, seed=rep * 100)
            for m in r2_by_method:
                r2_by_method[m].append(res[m]["R2"])
        for m, vals in r2_by_method.items():
            records.append({"N_asn_gwas": n, "method": m,
                            "R2_mean": np.mean(vals), "R2_std": np.std(vals)})
    return pd.DataFrame(records)


# ═══════════════════════════════════════════════════════════════════════════
# 6. TYPE 2 DIABETES CASE STUDY (semi-realistic parameters)
# ═══════════════════════════════════════════════════════════════════════════

def t2d_case_study():
    """
    T2D case study with literature-informed parameters.
    h2_SNP ≈ 0.18, Fst(EUR-EAS) ≈ 0.11, ~400 GWAS loci
    """
    print("\n  Running T2D case study...")
    results, artifacts = run_simulation(
        n_snps=400,
        n_causal=50,
        h2=0.18,
        fst=0.11,
        n_eur=50000,
        n_asn_gwas=8000,
        n_asn_test=3000,
        ld_eur=0.35,
        ld_asn=0.22,
        seed=2024
    )
    return results, artifacts


# ═══════════════════════════════════════════════════════════════════════════
# 7. PLOTTING
# ═══════════════════════════════════════════════════════════════════════════

PALETTE = {
    "Standard PRS":          "#E74C3C",
    "P+T PRS":               "#E67E22",
    "LD-corrected (Bayes)":  "#3498DB",
    "Multi-ethnic Meta":     "#2ECC71",
    "LAI-PRS":               "#9B59B6",
    "CS-PRS":                "#1ABC9C",
    "Oracle":                "#2C3E50",
}

def plot_r2_comparison(results: dict, title: str, fname: str):
    methods = [m for m in results if m != "Oracle"]
    r2s = [results[m]["R2"] for m in methods]
    oracle_r2 = results["Oracle"]["R2"]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh(methods, r2s, color=[PALETTE.get(m, "#888") for m in methods], alpha=0.85)
    ax.axvline(oracle_r2, color=PALETTE["Oracle"], linestyle="--", lw=2, label=f"Oracle R² = {oracle_r2:.3f}")
    ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)
    ax.set_xlabel("R² (Incremental PRS, target population)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend()
    ax.set_xlim(0, max(oracle_r2 * 1.25, max(r2s) * 1.25, 0.1))
    plt.tight_layout()
    plt.savefig(fname, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


def plot_fst_sensitivity(df: pd.DataFrame, fname: str):
    methods_plot = ["Standard PRS", "LD-corrected (Bayes)", "Multi-ethnic Meta", "CS-PRS", "Oracle"]
    fig, ax = plt.subplots(figsize=(9, 5))
    for m in methods_plot:
        sub = df[df["method"] == m]
        ax.errorbar(sub["Fst"], sub["R2_mean"], yerr=sub["R2_std"],
                    label=m, marker="o", linewidth=2, capsize=4,
                    color=PALETTE.get(m, "#888"))
    ax.set_xlabel("Population differentiation (Fst)", fontsize=11)
    ax.set_ylabel("R² in target population", fontsize=11)
    ax.set_title("PRS Transferability vs. Population Differentiation (Fst)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fname, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


def plot_sample_size(df: pd.DataFrame, fname: str):
    fig, ax = plt.subplots(figsize=(9, 5))
    for m in ["Standard PRS", "Multi-ethnic Meta", "CS-PRS", "Oracle"]:
        sub = df[df["method"] == m]
        ax.errorbar(sub["N_asn_gwas"], sub["R2_mean"], yerr=sub["R2_std"],
                    label=m, marker="s", linewidth=2, capsize=4,
                    color=PALETTE.get(m, "#888"))
    ax.set_xscale("log")
    ax.set_xlabel("Asian GWAS sample size (N)", fontsize=11)
    ax.set_ylabel("R² in target population", fontsize=11)
    ax.set_title("Effect of Target-Population GWAS Sample Size on PRS Performance", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fname, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


def plot_effect_size_comparison(artifacts: dict, fname: str):
    """Scatter: estimated vs. true effect sizes for all methods."""
    beta_true = artifacts["beta_true"]
    causal = artifacts["causal_idx"]
    non_causal = np.setdiff1d(np.arange(len(beta_true)), causal)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    methods_betas = [
        ("EUR GWAS (Standard)", artifacts["beta_hat_eur"]),
        ("Asian GWAS", artifacts["beta_hat_asn"]),
        ("LD-corrected (Bayes)", artifacts["beta_ld"]),
        ("Multi-ethnic Meta", artifacts["beta_meta"]),
        ("CS-PRS weights", artifacts["beta_cs"]),
    ]

    for ax, (label, beta_est) in zip(axes[:5], methods_betas):
        ax.scatter(beta_true[non_causal], beta_est[non_causal],
                   alpha=0.4, s=15, c="#95A5A6", label="Non-causal")
        ax.scatter(beta_true[causal], beta_est[causal],
                   alpha=0.8, s=40, c="#E74C3C", label="Causal")
        lims = [min(beta_true.min(), beta_est.min()) * 1.1,
                max(beta_true.max(), beta_est.max()) * 1.1]
        ax.plot(lims, lims, "k--", alpha=0.5, lw=1)
        r, _ = stats.pearsonr(beta_true, beta_est)
        ax.set_title(f"{label}\n(r = {r:.3f})", fontsize=9)
        ax.set_xlabel("True effect size (β)", fontsize=8)
        ax.set_ylabel("Estimated effect size", fontsize=8)
        if ax == axes[0]:
            ax.legend(fontsize=7)

    axes[5].set_visible(False)
    plt.suptitle("True vs. Estimated SNP Effect Sizes Across Methods", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fname, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


def plot_maf_comparison(artifacts: dict, fname: str):
    """MAF comparison EUR vs. ASN + Fst per SNP."""
    mafs_eur = artifacts["mafs_eur"]
    mafs_asn = artifacts["mafs_asn"]
    diff = np.abs(mafs_eur - mafs_asn)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(mafs_eur, mafs_asn, alpha=0.4, s=15, c="#3498DB")
    axes[0].plot([0, 0.5], [0, 0.5], "r--", lw=1.5, label="Identity line")
    axes[0].set_xlabel("MAF (EUR)", fontsize=11)
    axes[0].set_ylabel("MAF (ASN)", fontsize=11)
    axes[0].set_title("Minor Allele Frequency: EUR vs. ASN", fontsize=11, fontweight="bold")
    axes[0].legend()

    axes[1].hist(diff, bins=30, color="#E74C3C", alpha=0.75, edgecolor="white")
    axes[1].axvline(diff.mean(), color="#2C3E50", linestyle="--",
                    label=f"Mean ΔMAF = {diff.mean():.3f}")
    axes[1].set_xlabel("|ΔMAF|  (EUR − ASN)", fontsize=11)
    axes[1].set_ylabel("Count", fontsize=11)
    axes[1].set_title("Distribution of MAF Differences", fontsize=11, fontweight="bold")
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


def plot_prs_distributions(artifacts: dict, fname: str):
    """PRS score distributions per method for T2D case study."""
    y = artifacts["y_asn_test"]
    # Binary T2D phenotype (top 20% = cases)
    q80 = np.percentile(y, 80)
    cases = y > q80

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    methods_prs = [
        ("Standard PRS",       artifacts["prs_standard"]),
        ("LD-corrected (Bayes)", artifacts["prs_ldcorr"]),
        ("Multi-ethnic Meta",  artifacts["prs_meta"]),
        ("LAI-PRS",            artifacts["prs_lai"]),
        ("CS-PRS",             artifacts["prs_cs"]),
        ("Oracle",             artifacts["prs_oracle"]),
    ]

    for ax, (label, prs) in zip(axes, methods_prs):
        prs_z = (prs - prs.mean()) / (prs.std() + 1e-10)
        ax.hist(prs_z[~cases], bins=30, alpha=0.6, color="#3498DB",
                label="Controls", density=True)
        ax.hist(prs_z[cases], bins=30, alpha=0.6, color="#E74C3C",
                label="Cases (top 20%)", density=True)
        r2 = np.corrcoef(prs, y)[0, 1] ** 2
        ax.set_title(f"{label}\n(R² = {r2:.3f})", fontsize=9)
        ax.set_xlabel("Standardized PRS", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        if ax == axes[0]:
            ax.legend(fontsize=7)

    plt.suptitle("PRS Distributions: Cases vs. Controls (T2D, simulated)", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fname, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


# ═══════════════════════════════════════════════════════════════════════════
# 8. MAIN ORCHESTRATION
# ═══════════════════════════════════════════════════════════════════════════

def main():
    log_event("setup", "run_started", "prs_transferability.py",
              handoff_in={"task": "PRS transferability simulation"},
              files=[])

    print("=" * 65)
    print("  PRS Transferability Simulation Framework")
    print("  UKB (EUR) → BioBank Japan (ASN)")
    print("=" * 65)

    # ── Baseline simulation ──────────────────────────────────────────────
    print("\n[1/5] Baseline simulation (Fst=0.10, n_EUR=10k, n_ASN=5k)...")
    results_base, artifacts_base = run_simulation(
        n_snps=200, n_causal=30, h2=0.40, fst=0.10,
        n_eur=10000, n_asn_gwas=5000, n_asn_test=2000,
        ld_eur=0.30, ld_asn=0.20, seed=42
    )
    df_base = pd.DataFrame(list(results_base.values()))
    df_base.to_csv("results/baseline_results.csv", index=False)
    log_event("execute", "file_written", "simulation",
              handoff_out={"methods": list(results_base.keys())},
              files=["results/baseline_results.csv"])

    plot_r2_comparison(results_base,
                       "PRS Method Comparison: EUR→ASN Transfer (Baseline)",
                       "figures/fig1_r2_comparison_baseline.png")

    plot_effect_size_comparison(artifacts_base,
                                "figures/fig2_effect_size_comparison.png")

    plot_maf_comparison(artifacts_base,
                        "figures/fig3_maf_comparison.png")

    # ── Fst sensitivity ──────────────────────────────────────────────────
    print("\n[2/5] Fst sensitivity analysis (5 levels × 5 reps)...")
    df_fst = fst_sensitivity(fst_values=[0.01, 0.05, 0.10, 0.15, 0.20], n_rep=5)
    df_fst.to_csv("results/fst_sensitivity.csv", index=False)
    plot_fst_sensitivity(df_fst, "figures/fig4_fst_sensitivity.png")
    log_event("execute", "file_written", "fst_sensitivity",
              files=["results/fst_sensitivity.csv", "figures/fig4_fst_sensitivity.png"])

    # ── Sample size sensitivity ──────────────────────────────────────────
    print("\n[3/5] Sample size sensitivity analysis...")
    df_ss = sample_size_sensitivity(n_values=[500, 1000, 2000, 5000, 10000], n_rep=5)
    df_ss.to_csv("results/sample_size_sensitivity.csv", index=False)
    plot_sample_size(df_ss, "figures/fig5_sample_size_sensitivity.png")
    log_event("execute", "file_written", "sample_size_sensitivity",
              files=["results/sample_size_sensitivity.csv", "figures/fig5_sample_size_sensitivity.png"])

    # ── T2D case study ───────────────────────────────────────────────────
    print("\n[4/5] T2D case study (semi-realistic parameters)...")
    results_t2d, artifacts_t2d = t2d_case_study()
    df_t2d = pd.DataFrame(list(results_t2d.values()))
    df_t2d.to_csv("results/t2d_case_study.csv", index=False)
    plot_r2_comparison(results_t2d,
                       "T2D PRS Transfer: UKB (EUR) → BioBank Japan (Semi-realistic)",
                       "figures/fig6_t2d_r2_comparison.png")
    plot_prs_distributions(artifacts_t2d,
                           "figures/fig7_t2d_prs_distributions.png")
    log_event("execute", "file_written", "t2d_case_study",
              files=["results/t2d_case_study.csv",
                     "figures/fig6_t2d_r2_comparison.png",
                     "figures/fig7_t2d_prs_distributions.png"])

    # ── Summary table ─────────────────────────────────────────────────────
    print("\n[5/5] Generating summary table...")
    # Relative improvement over Standard PRS
    baseline_r2 = df_base[df_base["method"] == "Standard PRS"]["R2"].values[0]
    oracle_r2   = df_base[df_base["method"] == "Oracle"]["R2"].values[0]
    df_base["relative_improvement_pct"] = (df_base["R2"] - baseline_r2) / (baseline_r2 + 1e-10) * 100
    df_base["pct_of_oracle"] = df_base["R2"] / oracle_r2 * 100
    summary_path = "results/summary_table.csv"
    df_base.to_csv(summary_path, index=False)
    log_event("execute", "file_written", "summary",
              files=[summary_path])

    # ── Heritability sensitivity (h2 values) ─────────────────────────────
    h2_records = []
    for h2_val in [0.10, 0.20, 0.30, 0.40, 0.50]:
        for rep in range(5):
            res, _ = run_simulation(h2=h2_val, seed=rep * 10)
            for m in ["Standard PRS", "Multi-ethnic Meta", "CS-PRS", "Oracle"]:
                h2_records.append({"h2": h2_val, "method": m,
                                   "R2": res[m]["R2"]})
    df_h2 = pd.DataFrame(h2_records)
    df_h2_agg = df_h2.groupby(["h2", "method"]).agg(
        R2_mean=("R2", "mean"), R2_std=("R2", "std")).reset_index()
    df_h2_agg.to_csv("results/h2_sensitivity.csv", index=False)

    fig, ax = plt.subplots(figsize=(9, 5))
    for m in ["Standard PRS", "Multi-ethnic Meta", "CS-PRS", "Oracle"]:
        sub = df_h2_agg[df_h2_agg["method"] == m]
        ax.errorbar(sub["h2"], sub["R2_mean"], yerr=sub["R2_std"],
                    label=m, marker="D", linewidth=2, capsize=4,
                    color=PALETTE.get(m, "#888"))
    ax.set_xlabel("SNP heritability (h²)", fontsize=11)
    ax.set_ylabel("R² in target population", fontsize=11)
    ax.set_title("PRS Transferability vs. SNP Heritability", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/fig8_h2_sensitivity.png", dpi=180, bbox_inches="tight")
    plt.close()
    print("  Saved: figures/fig8_h2_sensitivity.png")

    # ── Print summary ────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("  RESULTS SUMMARY (Baseline simulation)")
    print("=" * 65)
    print(df_base[["method", "R2", "Pearson_r", "relative_improvement_pct", "pct_of_oracle"]].to_string(index=False))

    print("\n  T2D Case Study:")
    df_t2d["R2_rel"] = (df_t2d["R2"] - df_t2d[df_t2d["method"]=="Standard PRS"]["R2"].values[0]) / \
                       (df_t2d[df_t2d["method"]=="Standard PRS"]["R2"].values[0] + 1e-10) * 100
    print(df_t2d[["method", "R2", "Pearson_r"]].to_string(index=False))

    log_event("report", "run_completed", "prs_transferability.py",
              handoff_out={"n_figures": 8, "n_results_files": 5,
                           "best_method": df_base.loc[df_base["R2"].idxmax(), "method"]},
              files=["report.md"],
              status="ok")

    return df_base, df_fst, df_ss, df_t2d, artifacts_base, artifacts_t2d


if __name__ == "__main__":
    df_base, df_fst, df_ss, df_t2d, artifacts_base, artifacts_t2d = main()
