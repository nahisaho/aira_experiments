"""
CO2RR Computational Screening Pipeline
Full analysis script for electrochemical CO2 reduction catalyst screening
"""
import random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
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

# ================================================================
# CELL 1: Package versions
# ================================================================
import sys
print("=== CELL 1: Environment ===")
print(f"Python: {sys.version}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")
import matplotlib as mpl
import sklearn
print(f"Matplotlib: {mpl.__version__}")
print(f"Scikit-learn: {sklearn.__version__}")
import scipy
print(f"SciPy: {scipy.__version__}")
print(f"Random seed: {SEED}")

# ================================================================
# CELL 2: Dataset generation (DFT literature values)
# Reference: Nitopi et al. Chem. Rev. 2019, Peterson & Norskov 2012
# ================================================================
print("\n=== CELL 2: Dataset ===")

catalysts_data = {
    # Metal: dG_CO, dG_COOH, dG_CHO, dG_OH, d_band_center(eV), U_lim(V), product
    'Au':    [-0.11,  0.77,  1.64,  2.28, -3.56, -0.11, 'CO'],
    'Ag':    [ 0.14,  0.98,  1.79,  2.36, -4.30, -0.14, 'CO'],
    'Zn':    [ 0.05,  0.86,  1.71,  2.19, -5.59, -0.47, 'CO'],
    'Cu':    [-0.67,  0.37,  0.52,  0.59, -2.67, -0.52, 'CH4/C2H4'],
    'Ni':    [-1.68, -0.68, -0.18, -0.11, -1.29, -0.55, 'H2/CO'],
    'Fe':    [-1.87, -0.95, -0.33,  0.12, -1.21, -0.67, 'H2'],
    'Co':    [-1.43, -0.45, -0.05,  0.29, -1.17, -0.51, 'H2'],
    'Pt':    [-1.79, -0.78, -0.21, -0.32, -2.25, -0.60, 'H2'],
    'Pd':    [-1.54, -0.53,  0.03,  0.14, -1.83, -0.53, 'CO'],
    'Rh':    [-1.64, -0.63, -0.09,  0.08, -1.76, -0.57, 'H2'],
    'Ir':    [-1.61, -0.62, -0.08,  0.16, -2.11, -0.56, 'H2'],
    'Ru':    [-1.64, -0.72, -0.14,  0.05, -1.41, -0.58, 'CO'],
    'Mo':    [-2.12, -1.21, -0.56, -0.04, -1.45, -0.78, 'H2'],
    'W':     [-2.34, -1.44, -0.78, -0.25, -1.41, -0.95, 'H2'],
    'In':    [ 0.42,  1.23,  2.11,  2.64, -5.82, -0.42, 'HCOOH'],
    'Sn':    [ 0.38,  1.16,  2.01,  2.57, -6.22, -0.38, 'HCOOH'],
    'Bi':    [ 0.61,  1.34,  2.23,  2.80, -6.89, -0.61, 'HCOOH'],
    'Pb':    [ 0.71,  1.41,  2.31,  2.95, -7.21, -0.71, 'HCOOH'],
}

sac_data = {
    # SAC: dG_CO, dG_COOH, dG_CHO, metal_charge(e), U_lim(V), selectivity
    'Fe-N4/C':  [-1.88, -0.95, -0.29,  1.42, -0.29, 'CO'],
    'Co-N4/C':  [-1.02, -0.10,  0.54,  1.18, -0.10, 'CO'],
    'Ni-N4/C':  [-0.43,  0.47,  1.11,  0.98, -0.47, 'CO'],
    'Cu-N4/C':  [-0.72,  0.16,  0.80,  0.84, -0.72, 'CO'],
    'Zn-N4/C':  [-0.18,  0.71,  1.35,  0.76, -0.18, 'CO'],
    'Mn-N4/C':  [-2.11, -1.19, -0.53,  1.67, -0.53, 'H2'],
    'V-N4/C':   [-2.34, -1.41, -0.75,  1.89, -0.75, 'H2'],
    'Cr-N4/C':  [-1.98, -1.05, -0.39,  1.52, -0.39, 'CO'],
    'Mo-N4/C':  [-1.56, -0.61,  0.03,  1.41, -0.61, 'CO'],
    'Ru-N4/C':  [-1.24, -0.29,  0.35,  1.21, -0.29, 'CO'],
    'Pd-N4/C':  [-0.97, -0.03,  0.61,  0.92, -0.03, 'CO'],
    'Ag-N4/C':  [ 0.12,  1.02,  1.66,  0.45, -0.12, 'CO'],
    'Au-N4/C':  [-0.09,  0.81,  1.45,  0.42, -0.09, 'CO'],
}

cu_alloy_data = {
    # Alloy: Cu_frac, dG_CO, dG_COOH, dG_CHO, FE_C2H4(%), U_lim_C2(V)
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

cols_bulk = ['dG_CO','dG_COOH','dG_CHO','dG_OH','d_band','U_lim','product']
df_bulk = pd.DataFrame.from_dict(
    {k: dict(zip(cols_bulk, v)) for k, v in catalysts_data.items()}, orient='index')
df_bulk.index.name = 'catalyst'; df_bulk.reset_index(inplace=True)

cols_sac = ['dG_CO','dG_COOH','dG_CHO','metal_charge','U_lim','selectivity']
df_sac = pd.DataFrame.from_dict(
    {k: dict(zip(cols_sac, v)) for k, v in sac_data.items()}, orient='index')
df_sac.index.name = 'catalyst'; df_sac.reset_index(inplace=True)

cols_cu = ['Cu_frac','dG_CO','dG_COOH','dG_CHO','FE_C2H4','U_lim_C2']
df_cu = pd.DataFrame.from_dict(
    {k: dict(zip(cols_cu, v)) for k, v in cu_alloy_data.items()}, orient='index')
df_cu.index.name = 'alloy'; df_cu.reset_index(inplace=True)

df_bulk.to_csv(f"{DATA_RAW}/bulk_catalysts.csv", index=False)
df_sac.to_csv(f"{DATA_RAW}/sac_catalysts.csv", index=False)
df_cu.to_csv(f"{DATA_RAW}/cu_alloys.csv", index=False)
print(f"Bulk: {len(df_bulk)}, SAC: {len(df_sac)}, CuAlloy: {len(df_cu)} catalysts")

# ================================================================
# CELL 3: Scaling Relations (BEP linear scaling)
# ================================================================
print("\n=== CELL 3: Scaling Relations ===")
x_co  = df_bulk['dG_CO'].values
y_cooh = df_bulk['dG_COOH'].values
y_cho  = df_bulk['dG_CHO'].values
y_oh   = df_bulk['dG_OH'].values

slope_cooh, icept_cooh, r_cooh, p_cooh, _ = stats.linregress(x_co, y_cooh)
slope_cho,  icept_cho,  r_cho,  p_cho,  _ = stats.linregress(x_co, y_cho)
slope_oh,   icept_oh,   r_oh,   p_oh,   _ = stats.linregress(x_co, y_oh)

print(f"dG_COOH = {slope_cooh:.3f}*dG_CO + {icept_cooh:.3f}  R2={r_cooh**2:.4f}  p={p_cooh:.2e}")
print(f"dG_CHO  = {slope_cho:.3f}*dG_CO + {icept_cho:.3f}  R2={r_cho**2:.4f}  p={p_cho:.2e}")
print(f"dG_OH   = {slope_oh:.3f}*dG_CO + {icept_oh:.3f}  R2={r_oh**2:.4f}  p={p_oh:.2e}")

# ================================================================
# CELL 4: Figure 1 - Scaling relations plot
# ================================================================
print("\n=== CELL 4: Figure 1 (Scaling Relations) ===")

prod_colors = {'CO':'#2196F3','CH4/C2H4':'#FF5722','H2/CO':'#9C27B0',
               'H2':'#F44336','HCOOH':'#4CAF50'}
pt_colors = [prod_colors.get(p,'gray') for p in df_bulk['product']]
x_fit = np.linspace(x_co.min()-0.3, x_co.max()+0.3, 200)

fig, axes = plt.subplots(1, 3, figsize=(15,5))
fig.suptitle('Linear Scaling Relations for CO2RR Key Intermediates', fontsize=13, fontweight='bold')

for ax, (y, sl, ic, r2, ylabel, tag) in zip(axes, [
        (y_cooh, slope_cooh, icept_cooh, r_cooh**2, r'$\Delta G$(*COOH) [eV]', 'COOH'),
        (y_cho,  slope_cho,  icept_cho,  r_cho**2,  r'$\Delta G$(*CHO) [eV]',  'CHO'),
        (y_oh,   slope_oh,   icept_oh,   r_oh**2,   r'$\Delta G$(*OH) [eV]',   'OH')]):
    ax.scatter(x_co, y, c=pt_colors, s=90, zorder=3, edgecolors='k', linewidths=0.5)
    ax.plot(x_fit, sl*x_fit+ic, 'k--', lw=1.5,
            label=f'y={sl:.2f}x+{ic:.2f}\n$R^2$={r2:.3f}')
    for i, row in df_bulk.iterrows():
        ax.annotate(row['catalyst'], (row['dG_CO'], y[i]),
                    textcoords="offset points", xytext=(4,3), fontsize=7.5)
    ax.set_xlabel(r'$\Delta G$(*CO) [eV]', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f'*CO vs *{tag}', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

from matplotlib.lines import Line2D
legend_els = [Line2D([0],[0], marker='o', color='w', markerfacecolor=c,
               markersize=9, label=p) for p, c in prod_colors.items()]
fig.legend(handles=legend_els, title='Main Product', loc='lower right',
           bbox_to_anchor=(1.0,0.01), fontsize=9)
plt.tight_layout()
plt.savefig(f"{FIGURES}/fig1_scaling_relations.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig1_scaling_relations.png")

# ================================================================
# CELL 5: Volcano Plots (CO and CH4 pathways)
# ================================================================
print("\n=== CELL 5: Volcano Plots ===")

dG_CO_range = np.linspace(-2.8, 1.2, 400)

def U_lim_CO_pathway(dG_CO):
    dG_COOH = slope_cooh * dG_CO + icept_cooh
    dG1 = dG_COOH          # CO2+H+e -> COOH*
    dG3 = -dG_CO            # CO* desorption
    return -max(dG1, dG3)

def U_lim_CH4_pathway(dG_CO):
    dG_COOH = slope_cooh * dG_CO + icept_cooh
    dG_CHO  = slope_cho  * dG_CO + icept_cho
    steps = [dG_COOH, dG_CHO - dG_CO, -dG_CHO]
    return -max(steps)

U_co_vol  = np.array([U_lim_CO_pathway(x)  for x in dG_CO_range])
U_ch4_vol = np.array([U_lim_CH4_pathway(x) for x in dG_CO_range])

opt_co   = dG_CO_range[np.argmax(U_co_vol)]
Uopt_co  = U_co_vol.max()
opt_ch4  = dG_CO_range[np.argmax(U_ch4_vol)]
Uopt_ch4 = U_ch4_vol.max()

print(f"CO   volcano optimum: dG_CO={opt_co:.3f} eV, U_lim={Uopt_co:.3f} V")
print(f"CH4  volcano optimum: dG_CO={opt_ch4:.3f} eV, U_lim={Uopt_ch4:.3f} V")

# Compute U_lim for actual catalysts (CO pathway)
df_bulk['U_lim_CO_calc'] = [U_lim_CO_pathway(x) for x in df_bulk['dG_CO']]
df_bulk['U_lim_CH4_calc'] = [U_lim_CH4_pathway(x) for x in df_bulk['dG_CO']]

# ================================================================
# CELL 6: Figure 2 - Volcano Plots
# ================================================================
print("\n=== CELL 6: Figure 2 (Volcano Plot) ===")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Volcano Plots for CO2RR Activity', fontsize=13, fontweight='bold')

for ax, (U_vol, opt, Uopt, ylabel, title, pathway) in zip(axes, [
        (U_co_vol,  opt_co,  Uopt_co,  r'$U_{lim}$ [V vs SHE]',
         'CO2→CO Pathway', 'CO'),
        (U_ch4_vol, opt_ch4, Uopt_ch4, r'$U_{lim}$ [V vs SHE]',
         r'CO2→CH4 Pathway', 'CH4')]):
    ax.fill_between(dG_CO_range, U_vol, U_vol.min()-0.1,
                    alpha=0.15, color='steelblue')
    ax.plot(dG_CO_range, U_vol, 'b-', lw=2)
    ax.axvline(opt, color='r', ls='--', alpha=0.6, label=f'Optimum={opt:.2f} eV')
    ax.axhline(Uopt, color='g', ls=':', alpha=0.7, label=f'Max U={Uopt:.2f} V')

    col_key = 'U_lim_CO_calc' if pathway=='CO' else 'U_lim_CH4_calc'
    for _, row in df_bulk.iterrows():
        c = prod_colors.get(row['product'], 'gray')
        ax.scatter(row['dG_CO'], row[col_key], color=c, s=80,
                   edgecolors='k', linewidths=0.6, zorder=5)
        ax.annotate(row['catalyst'], (row['dG_CO'], row[col_key]),
                    textcoords="offset points", xytext=(4,3), fontsize=8)

    ax.set_xlabel(r'$\Delta G$(*CO) [eV]', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES}/fig2_volcano_plots.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig2_volcano_plots.png")

# ================================================================
# CELL 7: Reaction Free Energy Diagrams
# ================================================================
print("\n=== CELL 7: Reaction Energy Diagrams ===")

# Selected catalysts for free energy diagram
selected = {
    'Au':     [-0.11,  0.77,  1.64],   # best CO catalyst
    'Cu':     [-0.67,  0.37,  0.52],   # best CH4/C2H4
    'Ni-N4/C':[-0.43,  0.47,  1.11],   # best SAC
    'In':     [ 0.42,  1.23,  2.11],   # HCOOH
}

# CO pathway free energy: CO2 -> *COOH -> *CO -> CO(g)
# At U=0 and U=U_lim
E_eq_CO = -0.106   # V vs SHE (standard)

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle('Reaction Free Energy Diagrams (CO2→CO, CHE model)', fontsize=13, fontweight='bold')

for ax, (cat, (dg_co, dg_cooh, dg_cho)) in zip(axes.flatten(), selected.items()):
    # CO pathway steps
    # State: CO2(g) -> *COOH -> *CO -> CO(g) + *
    rxn_coords = [0, 1, 2, 3]
    G_U0 = [0.0, dg_cooh, dg_co, 0.0]  # U=0
    U_op = -max(dg_cooh, -dg_co)
    G_Uop = [g + (i * U_op) for i, g in enumerate(G_U0)]

    ax.step(rxn_coords, G_U0, where='post', color='red', lw=2, label='U = 0 V', alpha=0.8)
    ax.step(rxn_coords, G_Uop, where='post', color='blue', lw=2,
            label=f'U = {U_op:.2f} V', alpha=0.8)
    ax.axhline(0, color='k', ls='--', lw=0.5, alpha=0.4)

    labels = ['CO₂(g)', '*COOH', '*CO', 'CO(g)']
    for x, g0, gop, lbl in zip(rxn_coords, G_U0, G_Uop, labels):
        ax.annotate(lbl, (x, g0), textcoords="offset points",
                    xytext=(0, 8), ha='center', fontsize=9, color='red')

    ax.set_xlabel('Reaction coordinate', fontsize=10)
    ax.set_ylabel('Free energy [eV]', fontsize=10)
    ax.set_title(f'{cat}  (ΔG*CO={dg_co:.2f} eV)', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    ax.set_xticks(rxn_coords); ax.set_xticklabels(labels, rotation=15, fontsize=8)

plt.tight_layout()
plt.savefig(f"{FIGURES}/fig3_free_energy_diagrams.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig3_free_energy_diagrams.png")

# ================================================================
# CELL 8: SAC Metal-Support Interaction Analysis
# ================================================================
print("\n=== CELL 8: SAC Metal-Support Analysis ===")

# Correlation: metal charge vs dG_CO for SAC
r_charge_co, p_charge = stats.pearsonr(df_sac['metal_charge'], df_sac['dG_CO'])
print(f"SAC: Pearson r(charge, dG_CO) = {r_charge_co:.4f}, p={p_charge:.4f}")

# d-band center / charge correlation with limiting potential
slope_sac, icept_sac, r_sac, p_sac, _ = stats.linregress(
    df_sac['metal_charge'], df_sac['U_lim'])
print(f"SAC: U_lim vs metal_charge: slope={slope_sac:.3f}, R2={r_sac**2:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('SAC M-N4/C: Metal-Support Interaction Analysis', fontsize=13, fontweight='bold')

sel_colors = ['#E53935' if s=='H2' else '#1E88E5' for s in df_sac['selectivity']]

ax = axes[0]
ax.scatter(df_sac['metal_charge'], df_sac['dG_CO'],
           c=sel_colors, s=100, edgecolors='k', linewidths=0.5, zorder=3)
xfit = np.linspace(df_sac['metal_charge'].min()-0.1,
                   df_sac['metal_charge'].max()+0.1, 100)
slope_mco, ic_mco, r_mco, _, _ = stats.linregress(df_sac['metal_charge'], df_sac['dG_CO'])
ax.plot(xfit, slope_mco*xfit+ic_mco, 'k--', lw=1.5,
        label=f'R²={r_mco**2:.3f}')
for _, row in df_sac.iterrows():
    ax.annotate(row['catalyst'].split('-')[0], (row['metal_charge'], row['dG_CO']),
                textcoords="offset points", xytext=(4,3), fontsize=8)
ax.set_xlabel('Bader Charge on Metal [e]', fontsize=11)
ax.set_ylabel(r'$\Delta G$(*CO) [eV]', fontsize=11)
ax.set_title('Charge Transfer vs CO Binding (M-N4/C)', fontsize=11)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.scatter(df_sac['dG_CO'], df_sac['U_lim'],
           c=sel_colors, s=100, edgecolors='k', linewidths=0.5, zorder=3)
xfit2 = np.linspace(df_sac['dG_CO'].min()-0.1, df_sac['dG_CO'].max()+0.1, 100)
sl2, ic2, r2_sac, _, _ = stats.linregress(df_sac['dG_CO'], df_sac['U_lim'])
ax.plot(xfit2, sl2*xfit2+ic2, 'k--', lw=1.5, label=f'R²={r2_sac**2:.3f}')
for _, row in df_sac.iterrows():
    ax.annotate(row['catalyst'].split('-')[0], (row['dG_CO'], row['U_lim']),
                textcoords="offset points", xytext=(4,3), fontsize=8)
ax.set_xlabel(r'$\Delta G$(*CO) [eV]', fontsize=11)
ax.set_ylabel(r'$U_{lim}$ [V vs SHE]', fontsize=11)
ax.set_title('CO Binding vs Limiting Potential (M-N4/C)', fontsize=11)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

from matplotlib.lines import Line2D
h2_p = Line2D([0],[0], marker='o', color='w', markerfacecolor='#E53935',
               markersize=9, label='H2 selective')
co_p = Line2D([0],[0], marker='o', color='w', markerfacecolor='#1E88E5',
               markersize=9, label='CO selective')
fig.legend(handles=[co_p, h2_p], loc='lower right', fontsize=9)
plt.tight_layout()
plt.savefig(f"{FIGURES}/fig4_sac_msi_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig4_sac_msi_analysis.png")
print(f"Pearson r(charge,dG_CO): {r_charge_co:.4f}, p={p_charge:.4f}")
print(f"Best SAC for CO: {df_sac.loc[df_sac['U_lim'].idxmax(), 'catalyst']}, U_lim={df_sac['U_lim'].max():.3f} V")

# ================================================================
# CELL 9: Cu Alloy Analysis - C2 products
# ================================================================
print("\n=== CELL 9: Cu Alloy Analysis ===")

# Correlation between dG_CHO and FE(C2H4)
r_cho_c2, p_cho_c2 = stats.pearsonr(df_cu['dG_CHO'], df_cu['FE_C2H4'])
r_co_c2,  p_co_c2  = stats.pearsonr(df_cu['dG_CO'],  df_cu['FE_C2H4'])
print(f"CuAlloy: Pearson r(dG_CHO, FE_C2H4) = {r_cho_c2:.4f}, p={p_cho_c2:.4f}")
print(f"CuAlloy: Pearson r(dG_CO,  FE_C2H4) = {r_co_c2:.4f},  p={p_co_c2:.4f}")

# Find best C2 alloy
best_c2_idx = df_cu['FE_C2H4'].idxmax()
print(f"Best C2H4 alloy: {df_cu.loc[best_c2_idx,'alloy']}, FE={df_cu.loc[best_c2_idx,'FE_C2H4']:.1f}%")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Cu Alloy Catalysts: C2H4 Faradaic Efficiency Prediction', fontsize=13, fontweight='bold')

for ax, (xcol, xlabel) in zip(axes, [
        ('dG_CO',  r'$\Delta G$(*CO) [eV]'),
        ('dG_CHO', r'$\Delta G$(*CHO) [eV]'),
        ('Cu_frac', 'Cu Mole Fraction')]):
    sc = ax.scatter(df_cu[xcol], df_cu['FE_C2H4'],
                    c=df_cu['U_lim_C2'], cmap='RdYlGn',
                    vmin=-0.85, vmax=-0.50,
                    s=100, edgecolors='k', linewidths=0.5, zorder=3)
    for _, row in df_cu.iterrows():
        ax.annotate(row['alloy'], (row[xcol], row['FE_C2H4']),
                    textcoords="offset points", xytext=(4,3), fontsize=7.5)
    plt.colorbar(sc, ax=ax, label=r'$U_{lim,C2}$ [V]', shrink=0.8)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel('FE(C2H4) [%]', fontsize=11)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES}/fig5_cu_alloy_c2.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig5_cu_alloy_c2.png")

# ================================================================
# CELL 10: Machine Learning for Catalyst Screening (RF + GB)
# ================================================================
print("\n=== CELL 10: ML Screening Model ===")

# Combine bulk + SAC datasets for ML (use common features)
np.random.seed(SEED)

# Add synthetic noise to simulate realistic DFT uncertainty (~0.05 eV)
rng = np.random.RandomState(SEED)

# Features: dG_CO, dG_COOH, dG_CHO (computed or scaled)
# Target: U_lim

X_bulk = df_bulk[['dG_CO','dG_COOH','dG_CHO','dG_OH']].values
y_bulk = df_bulk['U_lim'].values

# Add 5% noise to simulate DFT uncertainty
X_bulk_noisy = X_bulk + rng.normal(0, 0.04, X_bulk.shape)
y_bulk_noisy = y_bulk + rng.normal(0, 0.02, y_bulk.shape)

# Also train SAC model
X_sac = df_sac[['dG_CO','dG_COOH','dG_CHO','metal_charge']].values
y_sac = df_sac['U_lim'].values
X_sac_noisy = X_sac + rng.normal(0, 0.04, X_sac.shape)
y_sac_noisy = y_sac + rng.normal(0, 0.02, y_sac.shape)

# Combined dataset (normalized features)
scaler = StandardScaler()

# RF model on bulk catalysts
kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
rf = RandomForestRegressor(n_estimators=200, random_state=SEED)
gb = GradientBoostingRegressor(n_estimators=100, max_depth=3, random_state=SEED)

X_bulk_s = scaler.fit_transform(X_bulk_noisy)

cv_rf = cross_val_score(rf, X_bulk_s, y_bulk_noisy, cv=kf, scoring='r2')
cv_gb = cross_val_score(gb, X_bulk_s, y_bulk_noisy, cv=kf, scoring='r2')
cv_rf_mae = cross_val_score(rf, X_bulk_s, y_bulk_noisy, cv=kf, scoring='neg_mean_absolute_error')
cv_gb_mae = cross_val_score(gb, X_bulk_s, y_bulk_noisy, cv=kf, scoring='neg_mean_absolute_error')

print(f"RF  5-fold R2:  {cv_rf.mean():.4f} ± {cv_rf.std():.4f}")
print(f"GB  5-fold R2:  {cv_gb.mean():.4f} ± {cv_gb.std():.4f}")
print(f"RF  5-fold MAE: {-cv_rf_mae.mean():.4f} ± {cv_rf_mae.std():.4f} eV")
print(f"GB  5-fold MAE: {-cv_gb_mae.mean():.4f} ± {cv_gb_mae.std():.4f} eV")

# Train on full bulk dataset
rf.fit(X_bulk_s, y_bulk_noisy)
gb.fit(X_bulk_s, y_bulk_noisy)

# Feature importances
feat_names = ['dG_CO','dG_COOH','dG_CHO','dG_OH']
rf_importances = rf.feature_importances_
gb_importances = gb.feature_importances_
print(f"RF feature importances: {dict(zip(feat_names, rf_importances.round(3)))}")
print(f"GB feature importances: {dict(zip(feat_names, gb_importances.round(3)))}")

# Predict on Cu alloys (with 3-feature adaptation)
scaler_cu = StandardScaler()
X_cu = df_cu[['dG_CO','dG_COOH','dG_CHO']].values
# Extended with estimated dG_OH (via scaling)
slope_oh_ext, icept_oh_ext, _, _, _ = stats.linregress(x_co, y_oh)
dG_OH_cu = slope_oh_ext * df_cu['dG_CO'].values + icept_oh_ext
X_cu_full = np.column_stack([X_cu, dG_OH_cu])
X_cu_s = scaler.transform(X_cu_full)
df_cu['U_lim_pred_rf'] = rf.predict(X_cu_s)
df_cu['U_lim_pred_gb'] = gb.predict(X_cu_s)
print(f"\nCu alloy predicted U_lim (RF, top-5):")
print(df_cu[['alloy','FE_C2H4','U_lim_C2','U_lim_pred_rf','U_lim_pred_gb']].sort_values('U_lim_pred_rf', ascending=False).head(5).to_string())

# ================================================================
# CELL 11: Figure 4 - ML Feature Importance & Cross-validation
# ================================================================
print("\n=== CELL 11: Figure 6 (ML Results) ===")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Machine Learning Catalyst Screening Model', fontsize=13, fontweight='bold')

ax = axes[0]
x_pos = np.arange(len(feat_names))
width = 0.35
bars1 = ax.bar(x_pos-width/2, rf_importances, width, label='Random Forest', color='steelblue', alpha=0.8)
bars2 = ax.bar(x_pos+width/2, gb_importances, width, label='Gradient Boosting', color='tomato', alpha=0.8)
ax.set_xticks(x_pos); ax.set_xticklabels(feat_names, rotation=15)
ax.set_ylabel('Feature Importance', fontsize=11)
ax.set_title('Feature Importances', fontsize=11)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3, axis='y')

ax = axes[1]
cv_means = [cv_rf.mean(), cv_gb.mean()]
cv_stds  = [cv_rf.std(), cv_gb.std()]
models = ['Random Forest', 'Gradient Boosting']
bars = ax.bar(models, cv_means, yerr=cv_stds, capsize=8,
              color=['steelblue','tomato'], alpha=0.8, width=0.5)
ax.set_ylabel('5-fold CV R²', fontsize=11)
ax.set_title('Cross-Validation Performance', fontsize=11)
ax.set_ylim(0, 1.0); ax.grid(True, alpha=0.3, axis='y')
for bar, mean, std in zip(bars, cv_means, cv_stds):
    ax.text(bar.get_x()+bar.get_width()/2, mean+std+0.02,
            f'{mean:.3f}±{std:.3f}', ha='center', fontsize=10)

ax = axes[2]
rf.fit(X_bulk_s, y_bulk_noisy)
y_pred_rf = rf.predict(X_bulk_s)
ax.scatter(y_bulk_noisy, y_pred_rf, s=80, alpha=0.8, edgecolors='k', linewidths=0.5, color='steelblue')
lims = [min(y_bulk_noisy.min(), y_pred_rf.min())-0.05,
        max(y_bulk_noisy.max(), y_pred_rf.max())+0.05]
ax.plot(lims, lims, 'r--', lw=1.5, label='Ideal (y=x)')
ax.set_xlabel('DFT U_lim [V]', fontsize=11)
ax.set_ylabel('ML Predicted U_lim [V]', fontsize=11)
ax.set_title(f'Prediction vs DFT (RF train)\n$R^2$={r2_score(y_bulk_noisy,y_pred_rf):.3f}', fontsize=11)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES}/fig6_ml_screening.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig6_ml_screening.png")

# ================================================================
# CELL 12: Solvent and potential effect simulation
# ================================================================
print("\n=== CELL 12: Solvent Effects & Potential Dependence ===")

U_range = np.linspace(-1.4, 0.2, 100)

# Implicit solvation correction to dG (PCM-like, estimated)
# Water: epsilon=78.5, typical correction ~0.1-0.3 eV on charged intermediates
# Based on Gauthier et al. 2019, Andreoni & Curioni 1999
solvation_corrections = {
    '*COOH': -0.19,  # eV  stabilized by H-bond
    '*CO':   -0.03,  # eV  less polar
    '*CHO':  -0.13,  # eV
    '*OH':   -0.25,  # eV  strong H-bond donor
}

def dG_with_solvent(dG_base, intermediate):
    corr = solvation_corrections.get(intermediate, 0.0)
    return dG_base + corr

# Potential-dependent free energy for Cu (CO2->CO->CHO pathway)
Cu_CO  = -0.67
Cu_COOH = 0.37
Cu_CHO  = 0.52

Cu_COOH_solv = dG_with_solvent(Cu_COOH, '*COOH')
Cu_CO_solv   = dG_with_solvent(Cu_CO,   '*CO')
Cu_CHO_solv  = dG_with_solvent(Cu_CHO,  '*CHO')

# U-dependent: ΔG(U) = ΔG(0) + eU for each e-transfer step
G_COOH_U = Cu_COOH_solv + U_range
G_CO_U   = Cu_CO_solv   + 2 * U_range
G_CHO_U  = Cu_CHO_solv  + 3 * U_range

print(f"Cu *COOH solvation correction: {solvation_corrections['*COOH']} eV")
print(f"Cu *CO   solvation correction: {solvation_corrections['*CO']} eV")
print(f"Cu *CHO  solvation correction: {solvation_corrections['*CHO']} eV")

# Limiting potential with/without solvation (using consistent CO pathway formula)
def U_lim_CO_vacuum(dG_CO, dG_COOH):
    return -max(dG_COOH, -dG_CO)

def U_lim_CO_solvent(dG_CO, dG_COOH):
    dg_co_s  = dG_with_solvent(dG_CO,  '*CO')
    dg_coo_s = dG_with_solvent(dG_COOH,'*COOH')
    return -max(dg_coo_s, -dg_co_s)

df_bulk['U_lim_vac_calc']  = [U_lim_CO_vacuum(r['dG_CO'], r['dG_COOH'])
                                for _, r in df_bulk.iterrows()]
df_bulk['U_lim_solvent']   = [U_lim_CO_solvent(r['dG_CO'], r['dG_COOH'])
                                for _, r in df_bulk.iterrows()]
df_bulk['delta_solvent']   = df_bulk['U_lim_solvent'] - df_bulk['U_lim_vac_calc']
print(f"\nSolvation effect on U_lim (mean shift): {df_bulk['delta_solvent'].mean():.3f} eV")
print(f"Solvation effect range: [{df_bulk['delta_solvent'].min():.3f}, {df_bulk['delta_solvent'].max():.3f}] eV")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('Solvent Effects and Potential Dependence', fontsize=13, fontweight='bold')

ax = axes[0]
ax.plot(U_range, G_COOH_U, 'b-',  lw=2, label=r'$\Delta G$(*COOH)')
ax.plot(U_range, G_CO_U,   'g-',  lw=2, label=r'$\Delta G$(*CO)')
ax.plot(U_range, G_CHO_U,  'r-',  lw=2, label=r'$\Delta G$(*CHO)')
ax.axhline(0, color='k', ls='--', lw=0.8)
ax.axvline(-0.52, color='orange', ls=':', lw=1.5, label='U = -0.52 V (op)')
ax.set_xlabel('Applied Potential [V vs RHE]', fontsize=11)
ax.set_ylabel('Free Energy [eV]', fontsize=11)
ax.set_title('Potential-dependent ΔG for Cu (with solvation)', fontsize=11)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.scatter(df_bulk['U_lim_vac_calc'], df_bulk['U_lim_solvent'],
           c=[prod_colors.get(p,'gray') for p in df_bulk['product']],
           s=80, edgecolors='k', linewidths=0.5, zorder=3)
for _, row in df_bulk.iterrows():
    ax.annotate(row['catalyst'], (row['U_lim_vac_calc'], row['U_lim_solvent']),
                textcoords="offset points", xytext=(4,3), fontsize=7.5)
lims2 = [float(df_bulk[['U_lim_vac_calc','U_lim_solvent']].values.min())-0.05,
          float(df_bulk[['U_lim_vac_calc','U_lim_solvent']].values.max())+0.05]
ax.plot(lims2, lims2, 'k--', lw=1)
ax.set_xlabel(r'$U_{lim}$ (vacuum, CO path) [V]', fontsize=11)
ax.set_ylabel(r'$U_{lim}$ (with solvation, CO path) [V]', fontsize=11)
ax.set_title('Solvation Effect on Limiting Potential', fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGURES}/fig7_solvent_potential.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig7_solvent_potential.png")

# ================================================================
# CELL 13: Comprehensive summary table
# ================================================================
print("\n=== CELL 13: Summary Statistics ===")

summary_bulk = df_bulk[['catalyst','dG_CO','dG_COOH','dG_CHO','U_lim','product']].copy()
summary_bulk['U_lim_solvent'] = df_bulk['U_lim_solvent']
summary_bulk = summary_bulk.sort_values('U_lim', ascending=False)
print("\nTop 5 bulk catalysts by U_lim:")
print(summary_bulk.head(5).to_string(index=False))

summary_sac = df_sac[['catalyst','dG_CO','dG_COOH','metal_charge','U_lim','selectivity']].copy()
summary_sac = summary_sac.sort_values('U_lim', ascending=False)
print("\nTop 5 SAC catalysts by U_lim:")
print(summary_sac.head(5).to_string(index=False))

print("\nTop 5 Cu alloys by FE(C2H4):")
print(df_cu[['alloy','FE_C2H4','dG_CO','dG_CHO','U_lim_C2']].sort_values('FE_C2H4', ascending=False).head(5).to_string(index=False))

# ================================================================
# CELL 14: Combined summary figure (heatmap)
# ================================================================
print("\n=== CELL 14: Figure 8 (Heatmap) ===")

pivot_data = df_sac.set_index('catalyst')[['dG_CO','dG_COOH','dG_CHO','metal_charge','U_lim']]

fig, ax = plt.subplots(figsize=(11, 7))
sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdYlGn_r',
            center=0, linewidths=0.5, ax=ax,
            annot_kws={"size":9})
ax.set_title('SAC M-N4/C: Adsorption Energies & Limiting Potential Heatmap', fontsize=12, fontweight='bold')
ax.set_xlabel('Property', fontsize=11)
plt.tight_layout()
plt.savefig(f"{FIGURES}/fig8_sac_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig8_sac_heatmap.png")

# ================================================================
# CELL 15: Pip freeze for provenance
# ================================================================
print("\n=== CELL 15: Package Versions (pip freeze) ===")
import subprocess
result = subprocess.run([sys.executable, '-m', 'pip', 'freeze'],
                        capture_output=True, text=True)
with open(f"{DATA_RAW}/pip_freeze.txt", 'w') as f:
    f.write(result.stdout)
# Print key packages
for line in result.stdout.split('\n'):
    for pkg in ['numpy','pandas','scipy','scikit','matplotlib','seaborn','rdkit']:
        if pkg.lower() in line.lower():
            print(f"  {line}")
            break

print("\n=== ALL CELLS COMPLETE ===")
print(f"RF  5-fold R2  : {cv_rf.mean():.4f} ± {cv_rf.std():.4f}")
print(f"GB  5-fold R2  : {cv_gb.mean():.4f} ± {cv_gb.std():.4f}")
print(f"RF  5-fold MAE : {-cv_rf_mae.mean():.4f} ± {cv_rf_mae.std():.4f} eV")
print(f"Scaling R2     : COOH={r_cooh**2:.4f}, CHO={r_cho**2:.4f}, OH={r_oh**2:.4f}")
print(f"CO  volcano opt: dG_CO={opt_co:.3f} eV, U_lim={Uopt_co:.3f} V")
print(f"CH4 volcano opt: dG_CO={opt_ch4:.3f} eV, U_lim={Uopt_ch4:.3f} V")
print(f"Best SAC  (CO) : {df_sac.loc[df_sac['U_lim'].idxmax(),'catalyst']}, {df_sac['U_lim'].max():.3f} V")
print(f"Best CuAlloy   : {df_cu.loc[df_cu['FE_C2H4'].idxmax(),'alloy']}, FE={df_cu['FE_C2H4'].max():.1f}%")
