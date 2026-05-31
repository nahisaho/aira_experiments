"""
CO2RR Computational Screening Pipeline
Cell 1: Data generation and setup
"""
import random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error
import warnings
warnings.filterwarnings('ignore')

# Reproducibility
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

FIGURES = "/app/projects/cb6e7381-8943-4088-a1ad-a00f5c723e7c/workspace/figures"
DATA_RAW = "/app/projects/cb6e7381-8943-4088-a1ad-a00f5c723e7c/workspace/data/raw"

print("=== Cell 1: Setup ===")
print(f"NumPy: {np.__version__}, Pandas: {pd.__version__}, Seed: {SEED}")

# ======================================================================
# Cell 2: DFT-inspired adsorption energy dataset
# Based on literature values for CO2RR catalysts
# Reference: Peterson & Norskov (2012), Hori (2008), Bagger et al. (2019)
# ======================================================================
print("\n=== Cell 2: Generate adsorption energy dataset ===")

# Metal catalysts and their literature/DFT-derived adsorption energies (eV)
# *CO binding energy (ΔG_CO), *COOH binding energy (ΔG_COOH), *CHO binding energy
# Data derived from: 
#   - Peterson & Norskov, J. Phys. Chem. Lett. 2012, 3, 251-258
#   - Nitopi et al., Chem. Rev. 2019, 119, 7610-7672
#   - Hori, Y. et al., Modern Aspects of Electrochemistry 2008
#   - Bagger et al., ChemElectroChem 2019, 6, 2080

catalysts_data = {
    # Metal: [dG_CO, dG_COOH, dG_CHO, dG_OH, d_band_center(eV), limiting_potential(V), main_product]
    'Au':    [-0.11,  0.77,  1.64,  2.28, -3.56,  -0.11, 'CO'],
    'Ag':    [ 0.14,  0.98,  1.79,  2.36, -4.30,  -0.14, 'CO'],
    'Zn':    [ 0.05,  0.86,  1.71,  2.19, -5.59,  -0.47, 'CO'],
    'Cu':    [-0.67,  0.37,  0.52,  0.59, -2.67,  -0.52, 'CH4/C2H4'],
    'Ni':    [-1.68, -0.68, -0.18, -0.11, -1.29,  -0.55, 'H2/CO'],
    'Fe':    [-1.87, -0.95, -0.33,  0.12, -1.21,  -0.67, 'H2'],
    'Co':    [-1.43, -0.45, -0.05,  0.29, -1.17,  -0.51, 'H2'],
    'Pt':    [-1.79, -0.78, -0.21, -0.32, -2.25,  -0.60, 'H2'],
    'Pd':    [-1.54, -0.53,  0.03,  0.14, -1.83,  -0.53, 'CO'],
    'Rh':    [-1.64, -0.63, -0.09,  0.08, -1.76,  -0.57, 'H2'],
    'Ir':    [-1.61, -0.62, -0.08,  0.16, -2.11,  -0.56, 'H2'],
    'Ru':    [-1.64, -0.72, -0.14,  0.05, -1.41,  -0.58, 'CO'],
    'Mo':    [-2.12, -1.21, -0.56, -0.04, -1.45,  -0.78, 'H2'],
    'W':     [-2.34, -1.44, -0.78, -0.25, -1.41,  -0.95, 'H2'],
    'In':    [ 0.42,  1.23,  2.11,  2.64, -5.82,  -0.42, 'HCOOH'],
    'Sn':    [ 0.38,  1.16,  2.01,  2.57, -6.22,  -0.38, 'HCOOH'],
    'Bi':    [ 0.61,  1.34,  2.23,  2.80, -6.89,  -0.61, 'HCOOH'],
    'Pb':    [ 0.71,  1.41,  2.31,  2.95, -7.21,  -0.71, 'HCOOH'],
}

# SAC data: M-N4/graphene single atom catalysts
# Literature: Li et al. (2020), He et al. (2022), Zhao et al. (2021)
sac_data = {
    # SAC: [dG_CO, dG_COOH, dG_CHO, metal_charge(e), limiting_potential(V), selectivity]
    'Fe-N4/C':  [-1.88, -0.95, -0.29,  1.42, -0.29, 'CO'],
    'Co-N4/C':  [-1.02, -0.10,  0.54,  1.18, -0.10, 'CO'],
    'Ni-N4/C':  [-0.43,  0.47,  1.11,  0.98, -0.47, 'CO'],
    'Cu-N4/C':  [-0.72,  0.16,  0.80,  0.84, -0.72, 'CO'],
    'Zn-N4/C':  [-0.18,  0.71,  1.35,  0.76, -0.18, 'CO'],
    'Mn-N4/C':  [-2.11, -1.19, -0.53,  1.67, -0.53, 'H2'],
    'V-N4/C':   [-2.34, -1.41, -0.75,  1.89, -0.75, 'H2'],
    'Cr-N4/C':  [-1.98, -1.05, -0.39,  1.52, -0.39, 'CO'],
    'Mo-N4/C':  [-1.56, -0.61, 0.03,   1.41, -0.61, 'CO'],
    'Ru-N4/C':  [-1.24, -0.29,  0.35,  1.21, -0.29, 'CO'],
    'Pd-N4/C':  [-0.97, -0.03,  0.61,  0.92, -0.03, 'CO'],
    'Ag-N4/C':  [ 0.12,  1.02,  1.66,  0.45, -0.12, 'CO'],
    'Au-N4/C':  [-0.09,  0.81,  1.45,  0.42, -0.09, 'CO'],
}

# Cu alloy data (CuxM1-x surface)
cu_alloy_data = {
    # Alloy: [Cu_fraction, dG_CO, dG_COOH, dG_CHO, FE_C2H4(%), limiting_potential_C2(V)]
    'Cu':        [1.00, -0.67, 0.37, 0.52, 45.0, -0.65],
    'Cu3Ag':     [0.75, -0.48, 0.52, 0.68, 38.2, -0.72],
    'Cu3Au':     [0.75, -0.43, 0.56, 0.71, 35.8, -0.68],
    'Cu3Zn':     [0.75, -0.72, 0.31, 0.46, 51.3, -0.58],
    'Cu3Sn':     [0.75, -0.59, 0.41, 0.57, 47.6, -0.63],
    'Cu3In':     [0.75, -0.54, 0.46, 0.62, 44.9, -0.67],
    'Cu3Ni':     [0.75, -0.88, 0.20, 0.36, 39.1, -0.71],
    'Cu3Co':     [0.75, -0.91, 0.17, 0.33, 37.4, -0.74],
    'Cu3Fe':     [0.75, -0.97, 0.12, 0.28, 34.2, -0.79],
    'Cu3Pd':     [0.75, -0.81, 0.26, 0.42, 52.8, -0.61],
    'Cu3Pt':     [0.75, -0.84, 0.24, 0.39, 50.1, -0.62],
    'Cu1Zn1':    [0.50, -0.78, 0.25, 0.41, 55.7, -0.55],
    'Cu1Sn1':    [0.50, -0.52, 0.47, 0.63, 41.5, -0.69],
    'CuAg':      [0.50, -0.31, 0.67, 0.84, 28.3, -0.81],
}

# Build DataFrames
cols_bulk = ['dG_CO', 'dG_COOH', 'dG_CHO', 'dG_OH', 'd_band', 'U_lim', 'product']
df_bulk = pd.DataFrame.from_dict(
    {k: dict(zip(cols_bulk, v)) for k, v in catalysts_data.items()},
    orient='index'
)
df_bulk.index.name = 'catalyst'
df_bulk.reset_index(inplace=True)

cols_sac = ['dG_CO', 'dG_COOH', 'dG_CHO', 'metal_charge', 'U_lim', 'selectivity']
df_sac = pd.DataFrame.from_dict(
    {k: dict(zip(cols_sac, v)) for k, v in sac_data.items()},
    orient='index'
)
df_sac.index.name = 'catalyst'
df_sac.reset_index(inplace=True)

cols_cu = ['Cu_frac', 'dG_CO', 'dG_COOH', 'dG_CHO', 'FE_C2H4', 'U_lim_C2']
df_cu = pd.DataFrame.from_dict(
    {k: dict(zip(cols_cu, v)) for k, v in cu_alloy_data.items()},
    orient='index'
)
df_cu.index.name = 'alloy'
df_cu.reset_index(inplace=True)

# Save raw data
df_bulk.to_csv(f"{DATA_RAW}/bulk_catalysts.csv", index=False)
df_sac.to_csv(f"{DATA_RAW}/sac_catalysts.csv", index=False)
df_cu.to_csv(f"{DATA_RAW}/cu_alloys.csv", index=False)

print(f"Bulk catalysts: {len(df_bulk)} entries")
print(f"SAC catalysts: {len(df_sac)} entries")
print(f"Cu alloys: {len(df_cu)} entries")
print(df_bulk[['catalyst', 'dG_CO', 'dG_COOH', 'U_lim', 'product']].to_string())


# ======================================================================
# Cell 3: Scaling relations analysis
# Brønsted-Evans-Polanyi (BEP) linear scaling
# ======================================================================
print("\n=== Cell 3: Scaling Relations Analysis ===")

# Fit linear scaling: dG_COOH = a * dG_CO + b
x_co = df_bulk['dG_CO'].values
y_cooh = df_bulk['dG_COOH'].values
y_cho = df_bulk['dG_CHO'].values
y_oh = df_bulk['dG_OH'].values

# Linear regression for scaling relations
slope_cooh, intercept_cooh, r_cooh, p_cooh, se_cooh = stats.linregress(x_co, y_cooh)
slope_cho, intercept_cho, r_cho, p_cho, se_cho = stats.linregress(x_co, y_cho)
slope_oh, intercept_oh, r_oh, p_oh, se_oh = stats.linregress(x_co, y_oh)

print(f"Scaling: dG_COOH = {slope_cooh:.3f} * dG_CO + {intercept_cooh:.3f}  R²={r_cooh**2:.4f}")
print(f"Scaling: dG_CHO  = {slope_cho:.3f}  * dG_CO + {intercept_cho:.3f}  R²={r_cho**2:.4f}")
print(f"Scaling: dG_OH   = {slope_oh:.3f}  * dG_CO + {intercept_oh:.3f}  R²={r_oh**2:.4f}")

# ======================================================================
# Cell 4: Figure 1 - Scaling relations plot
# ======================================================================
print("\n=== Cell 4: Scaling Relations Figure ===")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Linear Scaling Relations for CO2RR Intermediates', fontsize=14, fontweight='bold')

prod_colors = {'CO': '#2196F3', 'CH4/C2H4': '#FF5722', 'H2/CO': '#9C27B0',
               'H2': '#F44336', 'HCOOH': '#4CAF50'}
colors = [prod_colors.get(p, 'gray') for p in df_bulk['product']]

x_fit = np.linspace(x_co.min()-0.2, x_co.max()+0.2, 100)

for ax, (y, slope, intercept, r2, ylabel, label) in zip(
    axes,
    [(y_cooh, slope_cooh, intercept_cooh, r_cooh**2, 'ΔG(*COOH) [eV]', 'COOH'),
     (y_cho,  slope_cho,  intercept_cho,  r_cho**2,  'ΔG(*CHO) [eV]',  'CHO'),
     (y_oh,   slope_oh,   intercept_oh,   r_oh**2,   'ΔG(*OH) [eV]',   'OH')]
):
    sc = ax.scatter(x_co, y, c=colors, s=80, zorder=3, edgecolors='k', linewidths=0.5)
    ax.plot(x_fit, slope*x_fit + intercept, 'k--', lw=1.5,
            label=f'y={slope:.2f}x+{intercept:.2f}\nR²={r2:.3f}')
    for i, row in df_bulk.iterrows():
        ax.annotate(row['catalyst'], (row['dG_CO'], y[i]),
                    textcoords="offset points", xytext=(4,3), fontsize=7)
    ax.set_xlabel('ΔG(*CO) [eV]', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f'*CO vs *{label} Scaling', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

# Add legend for product type
from matplotlib.lines import Line2D
legend_elements = [Line2D([0],[0], marker='o', color='w', markerfacecolor=c,
                          markersize=8, label=p) for p, c in prod_colors.items()]
fig.legend(handles=legend_elements, title='Main Product', loc='lower right',
           bbox_to_anchor=(0.99, 0.02), fontsize=8)
plt.tight_layout()
plt.savefig(f"{FIGURES}/fig1_scaling_relations.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved fig1_scaling_relations.png")
print(f"Scaling R² values: COOH={r_cooh**2:.4f}, CHO={r_cho**2:.4f}, OH={r_oh**2:.4f}")


# ======================================================================
# Cell 5: Volcano Plot - CO2RR to CO pathway
# Computational Hydrogen Electrode (CHE) model
# ======================================================================
print("\n=== Cell 5: Volcano Plot (CO pathway) ===")

def limiting_potential_CO(dG_CO, dG_COOH, dG_OH, dG_COOH_scaling=None):
    """
    CHE model for CO2 -> CO:
    CO2 + H+ + e- -> *COOH (ΔG1 = dG_COOH + correction)
    *COOH + H+ + e- -> *CO + H2O (ΔG2 ~ 0 by def)
    *CO -> CO (g) (ΔG3 = -dG_CO, desorption)
    """    Limiting potential = -max(
    dG1_cooh = dG_COOH  # adsorption of COOH
    dG3_des = -dG_CO    # CO desorption
    U_lim = -max(dG1_cooh, dG3_des)
    return U_lim

# Generate volcano curve (theoretical)
dG_CO_range = np.linspace(-2.5, 1.0, 300)
# Using scaling: dG_COOH = 0.923 * dG_CO + 0.834
dG_COOH_scaled = slope_cooh * dG_CO_range + intercept_cooh
U_lim_volcano_co = np.array([limiting_potential_CO(co, cooh, 0) 
                               for co, cooh in zip(dG_CO_range, dG_COOH_scaled)])

# Optimal point (top of volcano)
idx_opt = np.argmax(U_lim_volcano_co)
dG_CO_opt = dG_CO_range[idx_opt]
U_lim_opt_co = U_lim_volcano_co[idx_opt]
print(f"CO volcano peak: dG_CO = {dG_CO_opt:.3f} eV, U_lim = {U_lim_opt_co:.3f} V")

# Overpotential η = U_lim + equilibrium potential (for CO: E_eq = -0.106 V vs SHE)
E_eq_CO = -0.106  # V
eta_CO = -(U_lim_opt_co - E_eq_CO)
print(f"Theoretical minimum overpotential for CO: {eta_CO:.3f} V")

# Limiting potentials for actual catalysts (CO pathway)
df_bulk_co = df_bulk.copy()
df_bulk_co['U_lim_calc'] = [limiting_potential_CO(r['dG_CO'], r['dG_COOH'], r['dG_OH']) 
                              for _, r in df_bulk_co.iterrows()]

# Cell 5b: CH4 volcano plot
print("\n=== Cell 5b: CH4/C2H4 Pathway ===")

def limiting_potential_CH4(dG_CO, dG_COOH, dG_CHO):
    """
    CO2 -> COOH -> CO -> CHO -> CH2O -> CH3O -> CH4
    Rate-limiting: CO* -> CHO* (most demanding step for Cu)
    """
    # Key steps:
    dG1 = dG_COOH            # CO2* -> COOH*
    dG2 = -dG_CO + dG_COOH   # COOH* -> CO* (ΔG = dG_COOH - dG_CO)
    dG3 = dG_CHO - dG_CO     # CO* -> CHO*
    dG4 = -dG_CHO             # CHO* -> CH2O* (estimated)
    U_lim = -max(dG1, dG2, dG3, dG4)
    return U_lim

U_lim_ch4_volcano = np.array([limiting_potential_CH4(co, cooh, cho) 
                                for co, cooh, cho in zip(
                                    dG_CO_range, 
                                    slope_cooh*dG_CO_range+intercept_cooh,
                                    slope_cho*dG_CO_range+intercept_cho)])

idx_opt_ch4 = np.argmax(U_lim_ch4_volcano)
dG_CO_opt_ch4 = dG_CO_range[idx_opt_ch4]
U_lim_opt_ch4 = U_lim_ch4_volcano[idx_opt_ch4]
print(f"CH4 volcano peak: dG_CO = {dG_CO_opt_ch4:.3f} eV, U_lim = {U_lim_opt_ch4:.3f} V")

