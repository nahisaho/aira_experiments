"""
Allosteric Transcription Factor-Based Biosensor Rational Design Framework
=========================================================================
Computational framework integrating:
  1. Ligand binding pocket analysis & docking simulation
  2. Allosteric communication pathway analysis
  3. Extended Hill equation dose-response modeling
  4. Variant library computational design
  5. Dynamic range optimization
  6. Environmental pollutant detection (heavy metals / organic solvents)

Author: GitHub Copilot (claude-sonnet-4.6)
Date  : 2026-05-31
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.optimize import curve_fit, minimize
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings('ignore')

# ── reproducibility ──────────────────────────────────────────────────────────
SEED = 42
np.random.seed(SEED)

# output paths
FIG_DIR = "figures"
DATA_DIR = "data/raw"

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Ligand Binding Pocket Analysis & Docking Simulation
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION 1: Ligand Binding Pocket Analysis & Docking Simulation")
print("=" * 70)

# Simulate docking scores for representative aTF systems
# Based on known aTF families: MerR (Hg), ArsR (As/Cd), CadC (Cd/Pb),
#   CueR (Cu/Ag), BenM (benzoate/adipic acid), TtgR (flavonoids/drugs)
tf_systems = {
    "MerR_Hg2+":  {"pocket_vol": 285, "hydrophobic_ratio": 0.31, "polar_contacts": 4, "docking_score": -8.7},
    "ArsR_As3+":  {"pocket_vol": 198, "hydrophobic_ratio": 0.22, "polar_contacts": 5, "docking_score": -7.9},
    "CadC_Cd2+":  {"pocket_vol": 312, "hydrophobic_ratio": 0.28, "polar_contacts": 6, "docking_score": -9.1},
    "CueR_Cu+":   {"pocket_vol": 178, "hydrophobic_ratio": 0.35, "polar_contacts": 3, "docking_score": -7.2},
    "BenM_adipate":{"pocket_vol": 425, "hydrophobic_ratio": 0.48, "polar_contacts": 8, "docking_score": -10.3},
    "TtgR_naringenin":{"pocket_vol": 390, "hydrophobic_ratio": 0.61, "polar_contacts": 5, "docking_score": -9.8},
    "HucR_urate": {"pocket_vol": 345, "hydrophobic_ratio": 0.40, "polar_contacts": 7, "docking_score": -10.1},
    "SRTF1_progesterone":{"pocket_vol": 510, "hydrophobic_ratio": 0.72, "polar_contacts": 4, "docking_score": -11.2},
}

df_pocket = pd.DataFrame(tf_systems).T.reset_index()
df_pocket.columns = ["TF_System","Pocket_Volume_A3","Hydrophobic_Ratio","Polar_Contacts","Docking_Score_kcal_mol"]
df_pocket.to_csv(f"{DATA_DIR}/pocket_analysis.csv", index=False)

print("\nLigand Binding Pocket Properties:")
print(df_pocket.to_string(index=False))

# Add noise for simulated replicate measurements (n=3)
docking_scores = df_pocket["Docking_Score_kcal_mol"].values
docking_noise = np.random.normal(0, 0.15, (3, len(docking_scores)))
docking_replicates = docking_scores + docking_noise  # shape (3, 8)
docking_mean = docking_replicates.mean(axis=0)
docking_std  = docking_replicates.std(axis=0)

print(f"\nMean docking score (all systems): {docking_mean.mean():.2f} ± {docking_std.mean():.2f} kcal/mol")

# Pearson r between pocket volume and docking score
r_pv_ds, p_pv_ds = pearsonr(df_pocket["Pocket_Volume_A3"].values, docking_mean)
print(f"Pearson r (pocket volume vs docking score): r={r_pv_ds:.3f}, p={p_pv_ds:.4f}")

r_hyd_ds, p_hyd_ds = pearsonr(df_pocket["Hydrophobic_Ratio"].values, docking_mean)
print(f"Pearson r (hydrophobic ratio vs docking score): r={r_hyd_ds:.3f}, p={p_hyd_ds:.4f}")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Allosteric Communication Pathway Analysis (MD-proxy)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SECTION 2: Allosteric Communication Pathway Analysis")
print("=" * 70)

# Simulate a covariance/mutual information matrix for allosteric network
# Representing residue-residue dynamic coupling in MerR (196 residue dimer)
N_RES = 40  # representative subset of key residues
np.random.seed(SEED)

# Ligand-binding domain: residues 0-14
# Allosteric linker:     residues 15-24
# DNA-binding domain:    residues 25-39
domain_labels = (["LBD"] * 15 + ["linker"] * 10 + ["DBD"] * 15)

# Create mutual information matrix with structured allosteric coupling
MI_matrix = np.zeros((N_RES, N_RES))
for i in range(N_RES):
    for j in range(N_RES):
        if i == j:
            MI_matrix[i, j] = 1.0
        else:
            # intra-domain coupling is stronger than inter-domain
            same_domain = (i < 15 and j < 15) or (15 <= i < 25 and 15 <= j < 25) or (i >= 25 and j >= 25)
            base = 0.65 if same_domain else 0.25
            # allosteric pathway: LBD → linker → DBD
            if (i < 15 and 15 <= j < 25) or (15 <= i < 25 and j < 15):
                base = 0.55
            elif (15 <= i < 25 and j >= 25) or (i >= 25 and 15 <= j < 25):
                base = 0.50
            elif i < 15 and j >= 25:
                base = 0.35
            MI_matrix[i, j] = base + np.random.normal(0, 0.06)

MI_matrix = (MI_matrix + MI_matrix.T) / 2
np.fill_diagonal(MI_matrix, 1.0)
MI_matrix = np.clip(MI_matrix, 0, 1)

# Identify top allosteric pathway residues (LBD -> DBD coupling)
lbd_indices  = list(range(0, 15))
dbd_indices  = list(range(25, 40))
cross_coupling = MI_matrix[np.ix_(lbd_indices, dbd_indices)]
top_lbd  = lbd_indices[np.argmax(cross_coupling.max(axis=1))]
top_dbd  = dbd_indices[np.argmax(cross_coupling.max(axis=0))]
max_cross = cross_coupling.max()

print(f"\nAllosteric coupling analysis (N={N_RES} residues)")
print(f"Peak LBD-DBD mutual information: {max_cross:.3f}")
print(f"  Key LBD residue index: {top_lbd+1}  |  Key DBD residue index: {top_dbd+1}")
print(f"Mean LBD-DBD coupling: {cross_coupling.mean():.3f} ± {cross_coupling.std():.3f}")
print(f"Mean intra-LBD coupling: {MI_matrix[:15,:15][np.triu_indices(15,1)].mean():.3f}")
print(f"Mean intra-DBD coupling: {MI_matrix[25:,25:][np.triu_indices(15,1)].mean():.3f}")
print(f"Mean linker coupling: {MI_matrix[15:25,15:25][np.triu_indices(10,1)].mean():.3f}")

# Save MI data
np.save(f"{DATA_DIR}/mi_matrix.npy", MI_matrix)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Extended Hill Equation Dose-Response Modeling
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SECTION 3: Extended Hill Equation Dose-Response Modeling")
print("=" * 70)

def hill_extended(L, Vmin, Vmax, K_d, n, K_coop):
    """Extended Hill equation with cooperativity term.

    Parameters
    ----------
    L        : ligand concentration (nM)
    Vmin     : basal reporter output (AU)
    Vmax     : maximum reporter output (AU)
    K_d      : apparent dissociation constant (nM)
    n        : Hill coefficient (cooperativity)
    K_coop   : cooperative binding constant (higher-order, nM)
    """
    term1 = (L / K_d)**n
    term2 = (L / K_coop)**(n + 1)
    return Vmin + (Vmax - Vmin) * (term1 + term2) / (1 + term1 + term2)

def hill_basic(L, Vmin, Vmax, K_d, n):
    """Standard Hill equation (4-parameter logistic)."""
    return Vmin + (Vmax - Vmin) * (L/K_d)**n / (1 + (L/K_d)**n)

# ── Simulate dose-response data for 6 target analytes ──────────────────────
analytes = {
    "Hg(II)":   {"K_d": 2.5,   "n": 1.8, "K_coop": 150,  "Vmin": 0.05, "Vmax": 1.0},
    "Cd(II)":   {"K_d": 12.0,  "n": 1.5, "K_coop": 800,  "Vmin": 0.08, "Vmax": 0.85},
    "As(III)":  {"K_d": 8.0,   "n": 1.3, "K_coop": 500,  "Vmin": 0.04, "Vmax": 0.92},
    "Cu(II)":   {"K_d": 5.0,   "n": 2.1, "K_coop": 200,  "Vmin": 0.06, "Vmax": 0.78},
    "Toluene":  {"K_d": 45.0,  "n": 1.1, "K_coop": 2000, "Vmin": 0.03, "Vmax": 0.65},
    "Benzene":  {"K_d": 80.0,  "n": 1.0, "K_coop": 4000, "Vmin": 0.02, "Vmax": 0.60},
}

conc_range = np.logspace(-2, 4, 300)  # 0.01 to 10000 nM
dr_results = {}
fit_params_all = {}

for analyte, params in analytes.items():
    true_output = hill_extended(conc_range, **params)
    # Add experimental noise (SNR ~20)
    noise_sigma = 0.025
    np.random.seed(SEED + list(analytes.keys()).index(analyte))
    noise = np.random.normal(0, noise_sigma, len(conc_range))
    measured = np.clip(true_output + noise, 0, 1.2)

    # Fit standard Hill equation
    try:
        popt, pcov = curve_fit(
            hill_basic, conc_range, measured,
            p0=[0.05, 0.9, params["K_d"], params["n"]],
            bounds=([0, 0, 0.01, 0.5], [0.5, 1.5, 1e4, 5]),
            maxfev=5000
        )
        perr = np.sqrt(np.diag(pcov))
        fitted_output = hill_basic(conc_range, *popt)
        residuals = measured - fitted_output
        r2 = r2_score(measured, fitted_output)
        fit_params_all[analyte] = {
            "Vmin": popt[0], "Vmax": popt[1], "K_d_fit": popt[2], "n_fit": popt[3],
            "K_d_err": perr[2], "n_err": perr[3], "R2": r2,
            "dynamic_range": popt[1] / popt[0] if popt[0] > 0 else np.inf,
            "LOD_nM": popt[2] * ((0.1 * (popt[1] - popt[0]) + popt[0] - popt[0]) /
                                 (popt[1] - popt[0]))**(1/popt[3]) if popt[3]>0 else np.nan
        }
    except Exception as e:
        print(f"  Fit failed for {analyte}: {e}")
        fit_params_all[analyte] = {}

    dr_results[analyte] = {"conc": conc_range, "measured": measured, "true": true_output}

df_fit = pd.DataFrame(fit_params_all).T.reset_index()
df_fit.columns = ["Analyte"] + list(df_fit.columns[1:])
df_fit.to_csv(f"{DATA_DIR}/dose_response_params.csv", index=False)

print("\nDose-Response Fitting Results (Standard Hill Equation):")
for a, p in fit_params_all.items():
    if p:
        dr = p.get('dynamic_range', float('nan'))
        dr_str = f"{dr:.1f}" if np.isfinite(dr) else ">100"
        print(f"  {a:12s}: Kd={p['K_d_fit']:.2f}±{p['K_d_err']:.2f} nM, "
              f"n={p['n_fit']:.2f}±{p['n_err']:.2f}, R²={p['R2']:.4f}, "
              f"DR={dr_str}x")

# Dynamic range (Vmax/Vmin) statistics
dr_values = [p["dynamic_range"] for p in fit_params_all.values() if p and np.isfinite(p.get("dynamic_range", np.inf))]
print(f"\nDynamic range: mean={np.mean(dr_values):.1f}x, "
      f"max={np.max(dr_values):.1f}x, min={np.min(dr_values):.1f}x")
print(f"Hg(II) LOD (10% activation): {fit_params_all['Hg(II)']['K_d_fit'] * 0.1**(1/fit_params_all['Hg(II)']['n_fit']):.3f} nM")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Variant Library Computational Design
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SECTION 4: Variant Library Computational Design")
print("=" * 70)

# Simulate a variant library (1000 mutants) with in silico predicted features
N_VARIANTS = 1000
np.random.seed(SEED)

# Amino acid substitution features (one-hot + physicochemical)
features = {
    "delta_volume":       np.random.normal(0, 15, N_VARIANTS),      # Å³
    "delta_hydrophobicity": np.random.normal(0, 0.8, N_VARIANTS),
    "charge_change":      np.random.choice([-2,-1,0,1,2], N_VARIANTS, p=[0.05,0.15,0.6,0.15,0.05]),
    "contact_count":      np.random.randint(2, 12, N_VARIANTS),
    "B_factor":           np.random.exponential(20, N_VARIANTS),     # flexibility
    "SASA_change":        np.random.normal(0, 25, N_VARIANTS),       # Å²
    "rosetta_ddG":        np.random.normal(0.5, 2.0, N_VARIANTS),    # REU
    "conservation_score": np.random.beta(3, 2, N_VARIANTS),          # 0-1
    "polarity_index":     np.random.uniform(0, 1, N_VARIANTS),
    "is_binding_site":    np.random.choice([0, 1], N_VARIANTS, p=[0.75, 0.25]),
    "is_allosteric_path": np.random.choice([0, 1], N_VARIANTS, p=[0.60, 0.40]),
}
df_variants = pd.DataFrame(features)

# Target: binding affinity change (ΔΔG_binding in kcal/mol, negative = improved)
# Physics-informed synthetic target with realistic structure
delta_G = (
    0.12 * df_variants["rosetta_ddG"]
    - 0.05 * df_variants["delta_hydrophobicity"] * df_variants["is_binding_site"]
    + 0.008 * df_variants["delta_volume"] * df_variants["is_allosteric_path"]
    - 0.15 * df_variants["conservation_score"]
    - 0.04 * df_variants["contact_count"] * df_variants["is_binding_site"]
    + 0.002 * df_variants["B_factor"]
    + np.random.normal(0, 0.35, N_VARIANTS)   # experimental noise
)
df_variants["delta_deltaG"] = delta_G

df_variants.to_csv(f"{DATA_DIR}/variant_library.csv", index=False)

# ── Machine-learning model for binding affinity prediction ─────────────────
feature_cols = [c for c in df_variants.columns if c != "delta_deltaG"]
X = df_variants[feature_cols].values
y = df_variants["delta_deltaG"].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

models = {
    "Random Forest":          RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1),
    "Gradient Boosting":      GradientBoostingRegressor(n_estimators=200, random_state=SEED),
}

kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
cv_results = {}

for name, model in models.items():
    cv_r2  = cross_val_score(model, X_scaled, y, cv=kf, scoring="r2")
    cv_mse = cross_val_score(model, X_scaled, y, cv=kf, scoring="neg_mean_squared_error")
    cv_results[name] = {
        "R2_mean": cv_r2.mean(), "R2_std": cv_r2.std(),
        "RMSE_mean": np.sqrt(-cv_mse).mean(), "RMSE_std": np.sqrt(-cv_mse).std()
    }
    print(f"  {name}: R²={cv_r2.mean():.3f}±{cv_r2.std():.3f}, "
          f"RMSE={np.sqrt(-cv_mse).mean():.3f}±{np.sqrt(-cv_mse).std():.3f} kcal/mol")

# Feature importance (RF)
rf_model = RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1)
rf_model.fit(X_scaled, y)
importances = rf_model.feature_importances_
feat_imp = pd.Series(importances, index=feature_cols).sort_values(ascending=False)
print("\nTop 5 features (RF importance):")
print(feat_imp.head(5).to_string())

# Select top-10% predicted improvers (ΔΔG < 0)
rf_pred = rf_model.predict(X_scaled)
top_mask = rf_pred < np.percentile(rf_pred, 10)
n_top = top_mask.sum()
print(f"\nTop 10% predicted improvers: N={n_top} variants")
print(f"Predicted ΔΔG_binding: {rf_pred[top_mask].mean():.3f}±{rf_pred[top_mask].std():.3f} kcal/mol")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Dynamic Range Optimization
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SECTION 5: Dynamic Range Optimization")
print("=" * 70)

# Model how promoter architecture affects dynamic range
def dynamic_range_model(params, alpha_basal, beta_max, sigma_RBS, n_op):
    """
    Model for reporter output dynamic range as function of:
      params:
        alpha_basal  : basal transcription (leakiness)
        beta_max     : max induced transcription
        sigma_RBS    : RBS strength (translation efficiency)
        n_op         : operator copy number
    """
    output_min = alpha_basal * sigma_RBS
    output_max = beta_max * sigma_RBS * n_op
    return output_max / output_min

# Parameter sweep for optimization
alpha_vals   = np.linspace(0.001, 0.10, 30)   # basal rate
beta_vals    = np.linspace(0.5, 5.0, 30)      # max induction
sigma_vals   = np.array([0.5, 1.0, 2.0, 4.0]) # RBS strength
n_op_vals    = np.array([1, 2, 3, 4])         # operator copies

# Fixed sigma=1.0, n_op=2 → sweep alpha × beta
DR_grid = np.zeros((30, 30))
for i, a in enumerate(alpha_vals):
    for j, b in enumerate(beta_vals):
        DR_grid[i, j] = b / a  # simplified ratio

# Find optimum
opt_i, opt_j = np.unravel_index(DR_grid.argmax(), DR_grid.shape)
print(f"Optimal alpha_basal: {alpha_vals[opt_i]:.4f} (minimized leakiness)")
print(f"Optimal beta_max:    {beta_vals[opt_j]:.2f}")
print(f"Theoretical max dynamic range: {DR_grid.max():.0f}x")

# Effect of operator copy number on dynamic range
dr_vs_nop = []
for n_op in n_op_vals:
    for sigma in sigma_vals:
        dr = dynamic_range_model(None, 0.01, 2.0, sigma, n_op)
        dr_vs_nop.append({"n_operators": n_op, "RBS_strength": sigma, "Dynamic_Range": dr})

df_dr = pd.DataFrame(dr_vs_nop)
print("\nDynamic range vs operator copies (beta_max=2.0, alpha=0.01):")
pivot = df_dr.pivot(index="n_operators", columns="RBS_strength", values="Dynamic_Range")
print(pivot.to_string())

# Optimized vs WT comparison
wt_dr  = dynamic_range_model(None, 0.05, 1.0, 1.0, 1)
opt_dr = dynamic_range_model(None, 0.005, 4.0, 2.0, 3)
fold_improvement = opt_dr / wt_dr
print(f"\nWT dynamic range:  {wt_dr:.0f}x")
print(f"Optimized dynamic range: {opt_dr:.0f}x  (fold improvement: {fold_improvement:.1f}x)")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Heavy Metal / Organic Solvent Detection Panel
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SECTION 6: Environmental Pollutant Detection Performance")
print("=" * 70)

# Simulated performance metrics for a panel of biosensors
pollutants = [
    # analyte, LOD(nM), LOQ(nM), linear_range_low(nM), linear_range_high(nM),
    # sensitivity(AU/nM), specificity_score(0-1), reg_threshold(nM), DR
    ("Hg(II)",   0.12, 0.45,  0.5,   50,   0.082, 0.95, 1.0,   18),
    ("Cd(II)",   0.65, 2.1,   2.0,  200,   0.041, 0.91, 3.0,   11),
    ("As(III)",  0.38, 1.2,   1.5,  120,   0.063, 0.89, 6.7,   14),
    ("Cu(II)",   0.29, 0.95,  1.0,   80,   0.074, 0.93, 13.0,  13),
    ("Pb(II)",   0.82, 2.8,   3.0,  250,   0.035, 0.87, 0.1,    9),
    ("Cr(VI)",   1.50, 4.9,   5.0,  500,   0.022, 0.84, 0.19,   8),
    ("Toluene",  45.0, 145,  150,  5000,  0.008, 0.78, 5700,   6),
    ("Benzene",  82.0, 265,  300,  9000,  0.005, 0.75,10000,   5),
    ("Xylene",   38.0, 122,  130,  4000,  0.009, 0.80, 1000,   7),
]

cols = ["Analyte","LOD_nM","LOQ_nM","Linear_Low_nM","Linear_High_nM",
        "Sensitivity_AU_nM","Specificity","Reg_Threshold_nM","Dynamic_Range_x"]
df_perf = pd.DataFrame(pollutants, columns=cols)
df_perf["Exceeds_Threshold"] = df_perf["LOD_nM"] < df_perf["Reg_Threshold_nM"]
df_perf["LOD_vs_threshold_ratio"] = df_perf["Reg_Threshold_nM"] / df_perf["LOD_nM"]
df_perf.to_csv(f"{DATA_DIR}/detection_performance.csv", index=False)

print("\nDetection Performance Summary:")
print(df_perf[["Analyte","LOD_nM","Sensitivity_AU_nM","Specificity",
               "Dynamic_Range_x","Exceeds_Threshold"]].to_string(index=False))

heavy_metal_mask = df_perf["Analyte"].str.contains(r"\(")
print(f"\nHeavy metal sensors: LOD range = {df_perf[heavy_metal_mask]['LOD_nM'].min():.2f} – "
      f"{df_perf[heavy_metal_mask]['LOD_nM'].max():.2f} nM")
print(f"Organic solvent sensors: LOD range = {df_perf[~heavy_metal_mask]['LOD_nM'].min():.1f} – "
      f"{df_perf[~heavy_metal_mask]['LOD_nM'].max():.1f} nM")
print(f"All sensors below regulatory threshold: {df_perf['Exceeds_Threshold'].all()}")
print(f"Mean specificity score: {df_perf['Specificity'].mean():.3f} ± {df_perf['Specificity'].std():.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE GENERATION
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

sns.set_theme(style="whitegrid", font_scale=1.0)
palette = sns.color_palette("tab10")

# ── Figure 1: Docking scores & pocket properties ───────────────────────────
fig1, axes = plt.subplots(1, 3, figsize=(15, 5))
fig1.suptitle("Figure 1: Ligand Binding Pocket Analysis", fontsize=13, fontweight="bold")

ax = axes[0]
colors = [palette[i] for i in range(len(df_pocket))]
bars = ax.bar(range(len(df_pocket)), docking_mean, yerr=docking_std,
              color=colors, capsize=4, alpha=0.85)
ax.set_xticks(range(len(df_pocket)))
ax.set_xticklabels([s.split("_")[0] for s in df_pocket["TF_System"]], rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Docking Score (kcal/mol)")
ax.set_title("(A) Docking Scores")
ax.axhline(docking_mean.mean(), color="red", linestyle="--", alpha=0.5, label=f"Mean: {docking_mean.mean():.1f}")
ax.legend(fontsize=8)

ax = axes[1]
ax.scatter(df_pocket["Pocket_Volume_A3"], docking_mean,
           c=colors, s=80, edgecolors="k", linewidth=0.5)
for i, row in df_pocket.iterrows():
    ax.annotate(row["TF_System"].split("_")[0], (row["Pocket_Volume_A3"], docking_mean[i]),
                fontsize=7, xytext=(3, 3), textcoords="offset points")
m, b = np.polyfit(df_pocket["Pocket_Volume_A3"], docking_mean, 1)
x_line = np.linspace(df_pocket["Pocket_Volume_A3"].min(), df_pocket["Pocket_Volume_A3"].max(), 100)
ax.plot(x_line, m * x_line + b, "r--", alpha=0.6, label=f"r={r_pv_ds:.2f}, p={p_pv_ds:.3f}")
ax.set_xlabel("Pocket Volume (Å³)")
ax.set_ylabel("Docking Score (kcal/mol)")
ax.set_title("(B) Volume vs Docking")
ax.legend(fontsize=8)

ax = axes[2]
ax.scatter(df_pocket["Hydrophobic_Ratio"], docking_mean,
           c=colors, s=80, edgecolors="k", linewidth=0.5)
for i, row in df_pocket.iterrows():
    ax.annotate(row["TF_System"].split("_")[0], (row["Hydrophobic_Ratio"], docking_mean[i]),
                fontsize=7, xytext=(3, 3), textcoords="offset points")
m2, b2 = np.polyfit(df_pocket["Hydrophobic_Ratio"], docking_mean, 1)
x2 = np.linspace(df_pocket["Hydrophobic_Ratio"].min(), df_pocket["Hydrophobic_Ratio"].max(), 100)
ax.plot(x2, m2 * x2 + b2, "r--", alpha=0.6, label=f"r={r_hyd_ds:.2f}, p={p_hyd_ds:.3f}")
ax.set_xlabel("Hydrophobic Ratio")
ax.set_ylabel("Docking Score (kcal/mol)")
ax.set_title("(C) Hydrophobicity vs Docking")
ax.legend(fontsize=8)

plt.tight_layout()
fig1.savefig(f"{FIG_DIR}/fig1_docking_analysis.png", dpi=150, bbox_inches="tight")
plt.close(fig1)
print("Saved fig1_docking_analysis.png")

# ── Figure 2: Allosteric MI heatmap ────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
fig2.suptitle("Figure 2: Allosteric Communication Pathway Analysis", fontsize=13, fontweight="bold")

ax = axes2[0]
domain_colors = {"LBD": "#e74c3c", "linker": "#f39c12", "DBD": "#2980b9"}
cmap = sns.diverging_palette(220, 10, as_cmap=True)
im = sns.heatmap(MI_matrix, ax=ax, cmap="YlOrRd", vmin=0, vmax=1,
                 xticklabels=False, yticklabels=False, cbar_kws={"label": "Mutual Information"})
# Domain boundary lines
ax.axhline(15, color="white", linewidth=1.5)
ax.axhline(25, color="white", linewidth=1.5)
ax.axvline(15, color="white", linewidth=1.5)
ax.axvline(25, color="white", linewidth=1.5)
ax.set_title("(A) Residue-Residue Mutual Information", fontsize=10)
ax.set_xlabel("Residue Index")
ax.set_ylabel("Residue Index")
# Domain labels
ax.text(7, -1, "LBD", ha="center", color="#e74c3c", fontsize=9, fontweight="bold")
ax.text(20, -1, "Linker", ha="center", color="#f39c12", fontsize=9, fontweight="bold")
ax.text(32, -1, "DBD", ha="center", color="#2980b9", fontsize=9, fontweight="bold")

ax2 = axes2[1]
# Average MI profile for each residue (sum of coupling)
coupling_profile = MI_matrix.sum(axis=1) - 1  # subtract self
x_res = np.arange(N_RES)
domain_color_list = [domain_colors[d] for d in domain_labels]
bars2 = ax2.bar(x_res, coupling_profile, color=domain_color_list, alpha=0.8)
ax2.axvline(15, color="gray", linestyle="--", alpha=0.5)
ax2.axvline(25, color="gray", linestyle="--", alpha=0.5)
ax2.set_xlabel("Residue Index")
ax2.set_ylabel("Total Coupling Score")
ax2.set_title("(B) Per-residue Allosteric Coupling", fontsize=10)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor="#e74c3c", label="LBD"),
                   Patch(facecolor="#f39c12", label="Linker"),
                   Patch(facecolor="#2980b9", label="DBD")]
ax2.legend(handles=legend_elements, fontsize=8)

plt.tight_layout()
fig2.savefig(f"{FIG_DIR}/fig2_allosteric_network.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
print("Saved fig2_allosteric_network.png")

# ── Figure 3: Dose-Response Curves ─────────────────────────────────────────
fig3, axes3 = plt.subplots(2, 3, figsize=(15, 10))
fig3.suptitle("Figure 3: Extended Hill Equation Dose-Response Modeling", fontsize=13, fontweight="bold")

analyte_list = list(analytes.keys())
for idx, (analyte, ax) in enumerate(zip(analyte_list, axes3.flatten())):
    data = dr_results[analyte]
    params = analytes[analyte]
    fp = fit_params_all.get(analyte, {})

    ax.semilogx(data["conc"], data["true"], "b-", linewidth=2, label="True (extended Hill)", alpha=0.8)
    ax.semilogx(data["conc"], data["measured"], ".", markersize=1.5, color="gray", alpha=0.4, label="Simulated data")

    if fp:
        fitted = hill_basic(data["conc"], fp["Vmin"], fp["Vmax"], fp["K_d_fit"], fp["n_fit"])
        ax.semilogx(data["conc"], fitted, "r--", linewidth=1.5,
                    label=f"Fit: Kd={fp['K_d_fit']:.1f} nM\nn={fp['n_fit']:.2f}, R²={fp['R2']:.3f}")
        # Mark LOD
        lod_y = fp["Vmin"] + 0.1 * (fp["Vmax"] - fp["Vmin"])
        ax.axhline(lod_y, color="orange", linestyle=":", alpha=0.7)
        ax.axvline(fp["K_d_fit"], color="green", linestyle=":", alpha=0.7)

    ax.set_xlabel("[Analyte] (nM)")
    ax.set_ylabel("Reporter Output (AU)")
    ax.set_title(f"({chr(65+idx)}) {analyte}", fontsize=10)
    ax.legend(fontsize=7)
    ax.set_xlim(0.01, 1e4)
    ax.set_ylim(-0.05, 1.25)

plt.tight_layout()
fig3.savefig(f"{FIG_DIR}/fig3_dose_response.png", dpi=150, bbox_inches="tight")
plt.close(fig3)
print("Saved fig3_dose_response.png")

# ── Figure 4: Variant Library ML Analysis ─────────────────────────────────
fig4, axes4 = plt.subplots(1, 3, figsize=(15, 5))
fig4.suptitle("Figure 4: Variant Library Computational Design", fontsize=13, fontweight="bold")

ax = axes4[0]
ax.hist(y, bins=50, color="steelblue", alpha=0.7, edgecolor="white", linewidth=0.5)
ax.axvline(0, color="red", linestyle="--", linewidth=1.5, label="ΔΔG=0 (WT)")
ax.axvline(np.percentile(y, 10), color="orange", linestyle="--", linewidth=1.5,
           label=f"10th pctile = {np.percentile(y, 10):.2f}")
ax.set_xlabel("ΔΔG_binding (kcal/mol)")
ax.set_ylabel("Count")
ax.set_title("(A) Variant ΔΔG Distribution")
ax.legend(fontsize=8)

ax = axes4[1]
ax.barh(feat_imp.index[:8][::-1], feat_imp.values[:8][::-1],
        color=palette[:8], alpha=0.85, edgecolor="k", linewidth=0.4)
ax.set_xlabel("Feature Importance (RF)")
ax.set_title("(B) Feature Importance")

ax = axes4[2]
ax.scatter(y, rf_pred, alpha=0.3, s=8, c="steelblue", edgecolors="none")
lim = [min(y.min(), rf_pred.min()), max(y.max(), rf_pred.max())]
ax.plot(lim, lim, "r--", linewidth=1.5, label="y=x")
r2_all = r2_score(y, rf_pred)
ax.set_xlabel("True ΔΔG (kcal/mol)")
ax.set_ylabel("Predicted ΔΔG (kcal/mol)")
ax.set_title(f"(C) RF Predictions (train R²={r2_all:.3f})")
ax.legend(fontsize=8)

plt.tight_layout()
fig4.savefig(f"{FIG_DIR}/fig4_variant_library.png", dpi=150, bbox_inches="tight")
plt.close(fig4)
print("Saved fig4_variant_library.png")

# ── Figure 5: Dynamic Range Optimization ──────────────────────────────────
fig5, axes5 = plt.subplots(1, 2, figsize=(12, 5))
fig5.suptitle("Figure 5: Dynamic Range Optimization", fontsize=13, fontweight="bold")

ax = axes5[0]
im5 = ax.contourf(beta_vals, alpha_vals, DR_grid, levels=20, cmap="plasma")
plt.colorbar(im5, ax=ax, label="Dynamic Range (x)")
ax.set_xlabel("β_max (max induction)")
ax.set_ylabel("α_basal (leakiness)")
ax.set_title("(A) DR Landscape (α × β sweep)")
ax.axhline(alpha_vals[opt_i], color="white", linestyle="--", linewidth=1)
ax.axvline(beta_vals[opt_j], color="white", linestyle="--", linewidth=1)

ax = axes5[1]
df_dr_pivot = df_dr.pivot(index="n_operators", columns="RBS_strength", values="Dynamic_Range")
sns.heatmap(df_dr_pivot, ax=ax, cmap="YlGnBu", annot=True, fmt=".0f",
            cbar_kws={"label": "Dynamic Range (x)"})
ax.set_xlabel("RBS Strength")
ax.set_ylabel("Operator Copy Number")
ax.set_title("(B) DR vs n_operators × RBS")

plt.tight_layout()
fig5.savefig(f"{FIG_DIR}/fig5_dynamic_range.png", dpi=150, bbox_inches="tight")
plt.close(fig5)
print("Saved fig5_dynamic_range.png")

# ── Figure 6: Detection Panel ───────────────────────────────────────────────
fig6, axes6 = plt.subplots(1, 2, figsize=(14, 6))
fig6.suptitle("Figure 6: Environmental Pollutant Detection Panel", fontsize=13, fontweight="bold")

ax = axes6[0]
analyte_names = df_perf["Analyte"].values
lods = df_perf["LOD_nM"].values
thresholds = df_perf["Reg_Threshold_nM"].values
x_pos = np.arange(len(analyte_names))
bars_lod = ax.bar(x_pos - 0.2, lods, width=0.4, color="steelblue", alpha=0.85, label="LOD (nM)")
bars_thr = ax.bar(x_pos + 0.2, thresholds, width=0.4, color="tomato", alpha=0.75, label="Regulatory limit (nM)")
ax.set_yscale("log")
ax.set_xticks(x_pos)
ax.set_xticklabels(analyte_names, rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Concentration (nM, log scale)")
ax.set_title("(A) LOD vs Regulatory Limits")
ax.legend(fontsize=9)

ax = axes6[1]
ax.scatter(df_perf["Dynamic_Range_x"], df_perf["Specificity"],
           s=100, c=range(len(df_perf)), cmap="tab10", edgecolors="k", linewidth=0.5, zorder=5)
for i, row in df_perf.iterrows():
    ax.annotate(row["Analyte"], (row["Dynamic_Range_x"], row["Specificity"]),
                fontsize=8, xytext=(3, 3), textcoords="offset points")
r_dr_sp, p_dr_sp = pearsonr(df_perf["Dynamic_Range_x"], df_perf["Specificity"])
ax.set_xlabel("Dynamic Range (x)")
ax.set_ylabel("Specificity Score")
ax.set_title(f"(B) DR vs Specificity (r={r_dr_sp:.2f}, p={p_dr_sp:.3f})")

plt.tight_layout()
fig6.savefig(f"{FIG_DIR}/fig6_detection_panel.png", dpi=150, bbox_inches="tight")
plt.close(fig6)
print("Saved fig6_detection_panel.png")

# ─────────────────────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)
print(f"\n[Cell 1] Docking scores: {docking_mean.mean():.2f} ± {docking_std.mean():.2f} kcal/mol (N=8 systems)")
print(f"[Cell 1] Pocket volume – docking correlation: r={r_pv_ds:.3f}, p={p_pv_ds:.4f}")
print(f"[Cell 1] Hydrophobicity – docking correlation: r={r_hyd_ds:.3f}, p={p_hyd_ds:.4f}")
print(f"\n[Cell 2] Peak LBD-DBD MI: {max_cross:.3f}")
print(f"[Cell 2] Mean LBD-DBD coupling: {cross_coupling.mean():.3f} ± {cross_coupling.std():.3f}")

print("\n[Cell 3] Dose-Response Hill Equation Fits:")
for a, p in fit_params_all.items():
    if p:
        print(f"  {a}: Kd={p['K_d_fit']:.2f} nM, n={p['n_fit']:.2f}, R²={p['R2']:.4f}")

print(f"\n[Cell 4] ML Variant Design:")
for name, res in cv_results.items():
    print(f"  {name}: 5-fold CV R²={res['R2_mean']:.3f}±{res['R2_std']:.3f}, "
          f"RMSE={res['RMSE_mean']:.3f}±{res['RMSE_std']:.3f} kcal/mol")

print(f"\n[Cell 5] Dynamic range: WT={wt_dr:.0f}x → Optimized={opt_dr:.0f}x ({fold_improvement:.1f}x improvement)")
print(f"\n[Cell 6] Detection: {df_perf['Exceeds_Threshold'].sum()}/{len(df_perf)} sensors below regulatory thresholds")
print(f"[Cell 6] Mean specificity: {df_perf['Specificity'].mean():.3f} ± {df_perf['Specificity'].std():.3f}")
print(f"[Cell 6] Dynamic range – specificity correlation: r={r_dr_sp:.3f}, p={p_dr_sp:.3f}")

print("\n✓ All analyses complete. Results saved in data/raw/ and figures/")
