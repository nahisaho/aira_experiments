"""
BIM-Integrated Environmental Performance Simulation System
Jupyter-equivalent analysis script — Cell 0 through Cell 8
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error
import warnings, json, random, os
warnings.filterwarnings('ignore')

# ── Cell 0: dirs & seeds ─────────────────────────────────────────────
np.random.seed(42)
random.seed(42)
os.makedirs("figures", exist_ok=True)
os.makedirs("data/raw", exist_ok=True)
print("[Cell 0] Seeds fixed (np=42, random=42), directories created.")

# ════════════════════════════════════════════════════════════
# Cell 1 : IFC Geometry Model — synthetic ZEB office building
# ════════════════════════════════════════════
print("\n[Cell 1] IFC Geometry Model")

ifc_data = {
    "building": "ZEB_Office_Tokyo",
    "floors": 5,
    "floor_area_m2": 1200,         # per floor
    "total_area_m2": 6000,
    "orientation_deg": 15,          # east of south15
    "glazing_ratio_south": 0.45,
    "glazing_ratio_north": 0.25,
    "glazing_ratio_east": 0.30,
    "glazing_ratio_west": 0.30,
    "wall_U_value": 0.25,           # W/m²K
    "roof_U_value": 0.15,
    "window_U_value": 1.20,
    "SHGC": 0.30,
    "infiltration_ACH": 0.10,
    "occupancy_density": 0.10,      # persons/m²
    "lighting_W_m2": 8.0,
    "equipment_W_m2": 15.0,
    "location": "Tokyo, Japan",
    "climate_zone": "Cfa (humid subtropical)",
    "HDD": 1340,                    # heating degree days
    "CDD": 1060,                    # cooling degree days
}

# Save IFC mock data
with open("data/raw/ifc_building_data.json", "w") as f:
    json.dump(ifc_data, f, indent=2)

df_ifc = pd.DataFrame([ifc_data])
print(df_ifc[["building","total_area_m2","wall_U_value","glazing_ratio_south","climate_zone"]].to_string())

# ═══════════════════════════════════════════
# Cell 2 : Thermal Load Simulation (EnergyPlus-proxy)
# ═══════════════════════════════════
print("\n[Cell 2] Thermal Load Simulation")

months = np.arange(1, 13)
month_names = ['Jan','Feb','Mar','Apr','May','Jun',
               'Jul','Aug','Sep','Oct','Nov','Dec']
T_out = np.array([-3, -1, 5, 12, 18, 23, 28, 30, 24, 17, 10, 2])   # Tokyo avg °C
T_setpoint_cool = 26.0
T_setpoint_heat = 22.0

# Internal gains [W] = (lighting + equipment + occupancy) * area
internal_gain_W = (ifc_data["lighting_W_m2"] + ifc_data["equipment_W_m2"]
                   + ifc_data["occupancy_density"] * 100) * ifc_data["total_area_m2"]

# Solar gains (simplified) [W] – south facade dominant
solar_irrad = np.array([160,190,260,340,380,360,330,320,280,240,175,145])  # W/m²
south_facade_area = 24 * 5 * ifc_data["floors"]  # width × height × floors
solar_gain_W = solar_irrad * south_facade_area * ifc_data["SHGC"] * ifc_data["glazing_ratio_south"]

# Transmission loss [W] – UA * ΔT
UA_total = (ifc_data["wall_U_value"] * 1800 +
            ifc_data["window_U_value"] * 900 +
            ifc_data["roof_U_value"] * 1200)   # approximate envelope UA
delta_T_cool = np.maximum(T_out - T_setpoint_cool, 0)
delta_T_heat = np.maximum(T_setpoint_heat - T_out, 0)

# Monthly energy demand [kWh] — 720 h/month
cooling_kWh = ((solar_gain_W + internal_gain_W + UA_total * delta_T_cool)
               / 3.5 * 720 / 1000)           # COP=3.5
heating_kWh = (UA_total * delta_T_heat
               / 4.0 * 720 / 1000)           # COP=4.0 heat pump

# PV generation (roof + facade BIPV)
pv_capacity_kWp = 480   # kWp
pv_efficiency = 0.20
# PV area: capacity [kWp] / (std irradiance 1 kW/m² × efficiency) → m²
pv_area_m2 = pv_capacity_kWp / pv_efficiency          # = 2400 m² (roof + BIPV facade)
# Monthly generation: irradiance[W/m²] × area[m²] × η × hours / 1000 → kWh
pv_gen_kWh = solar_irrad * pv_area_m2 * pv_efficiency * 720 / 1000

# Other building loads
lighting_kWh = np.full(12, ifc_data["lighting_W_m2"] * ifc_data["total_area_m2"] * 720 / 1000)
equip_kWh    = np.full(12, ifc_data["equipment_W_m2"] * ifc_data["total_area_m2"] * 720 / 1000)

total_demand_kWh  = cooling_kWh + heating_kWh + lighting_kWh + equip_kWh
net_energy_kWh    = total_demand_kWh - pv_gen_kWh
annual_demand     = total_demand_kWh.sum()
annual_pv         = pv_gen_kWh.sum()
annual_net        = net_energy_kWh.sum()
EUI               = annual_demand / ifc_data["total_area_m2"]  # kWh/m²/yr

print(f"  Annual demand : {annual_demand:,.0f} kWh/yr")
print(f"  PV generation : {annual_pv:,.0f} kWh/yr")
print(f"  Net energy    : {annual_net:,.0f} kWh/yr")
print(f"  EUI           : {EUI:.1f} kWh/m²/yr")

df_thermal = pd.DataFrame({
    "month": month_names,
    "T_out_C": T_out,
    "cooling_kWh": cooling_kWh.round(0),
    "heating_kWh": heating_kWh.round(0),
    "lighting_kWh": lighting_kWh.round(0),
    "equip_kWh"   : equip_kWh.round(0),
    "pv_gen_kWh"  : pv_gen_kWh.round(0),
    "net_kWh"     : net_energy_kWh.round(0),
})
df_thermal.to_csv("data/raw/monthly_energy.csv", index=False)
print(df_thermal[["month","cooling_kWh","heating_kWh","pv_gen_kWh","net_kWh"]].to_string(index=False))

# ══════════════════
# Cell 3 : CFD Natural Ventilation (simplified analytical model)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[Cell 3] CFD Natural Ventilation")

# Cross-ventilation: Q = Cd * A * sqrt(2*ΔP/ρ)
Cd = 0.65          # discharge coefficient
rho = 1.2          # kg/m³
wind_speeds = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])

# Opening area per floor (windward + leeward)
A_open = 24 * 2.5 * 0.30    # width × height × WWR = 18 m²
H_stack = 3.0                # floor height
T_indoor = 26.0              # °C indoor
T_outdoor_summer = 30.0      # °C outdoor (worst case)

# Wind-driven pressure difference
Cp_windward = 0.8
Cp_leeward  = -0.5
dCp = Cp_windward - Cp_leeward  # = 1.3

Q_wind = Cd * A_open * wind_speeds * np.sqrt(dCp)
# Buoyancy-driven (stack effect)
g = 9.81
dT = T_outdoor_summer - T_indoor
Q_stack = Cd * A_open * np.sqrt(2 * g * H_stack * abs(dT) / (T_indoor + 273))

# Air changes per hour
room_volume = ifc_data["floor_area_m2"] * 3.0  # m³ per floor
ACH_wind  = Q_wind  * 3600 / room_volume
ACH_stack = Q_stack * 3600 / room_volume

# Thermal comfort (ASHRAE 55 adaptive)
PMV_cooling = np.where(ACH_wind >= 4, 'Comfortable', 
              np.where(ACH_wind >= 2, 'Slightly warm', 'Hot'))

print(f"  Stack-effect Q = {Q_stack:.2f} m³/s  (ACH = {ACH_stack:.2f})")
print(f"  At 3 m/s wind  : Q = {Q_wind[5]:.2f} m³/s  (ACH = {ACH_wind[5]:.2f})")
print(f"  Comfort threshold 4 ACH reached at wind ≥ {wind_speeds[ACH_wind>=4][0]:.1f} m/s")

df_cfd = pd.DataFrame({
    "wind_speed_m_s": wind_speeds,
    "Q_wind_m3_s"   : Q_wind.round(3),
    "ACH_wind"      : ACH_wind.round(2),
    "comfort"       : PMV_cooling,
})
df_cfd.to_csv("data/raw/cfd_ventilation.csv", index=False)
print(df_cfd.to_string(index=False))

# ═══════════════════════════════════════════════════════════════
# Cell 4 : Daylight Simulation (Radiance/Honeybee proxy)
# ════════════════════════════════════════════════════
print("\n[Cell 4] Daylight Simulation")

# Spatial Daylight Autonomy (sDA) and Annual Sunlight Exposure (ASE)
np.random.seed(42)
n_points = 200   # sensor grid points

# Illuminance distribution (lognormal, typical office)
illum_mean = 350   # lux target
illum_std  = 180
illuminance = np.random.lognormal(np.log(illum_mean) - 0.5*np.log(1+(illum_std/illum_mean)**2),
                                   np.log(1+(illum_std/illum_mean)**2)**0.5, n_points)
illuminance = np.clip(illuminance, 50, 5000)

sDA300 = (illuminance >= 300).mean() * 100   # % points ≥ 300 lux (LEED target: ≥55%)
ASE1000 = (illuminance >= 1000).mean() * 100  # % points ≥ 1000 lux (< 10% acceptable)
cDA = np.percentile(illuminance, 50)           # median illuminance (Continuous DA)

# Unified Glare Rating proxy
UGR = 16.5 + 5.0 * np.log10(illum_mean / 500 + 0.01)  # simplified
DGP = 0.18 + 0.015 * np.log1p(illum_mean / 100)         # Daylight Glare Probability

print(f"  sDA(300 lux)  = {sDA300:.1f}%   (LEED target >=55%: {'PASS' if sDA300>=55 else 'FAIL'})")
print(f"  ASE(1000 lux) = {ASE1000:.1f}%  (target <10%:    {'PASS' if ASE1000<10 else 'FAIL'})")
print(f"  Median illuminance = {cDA:.0f} lux")
print(f"  UGR ≈ {UGR:.1f}  (target <19: {'PASS' if UGR<19 else 'FAIL'})")
print(f"  DGP ≈ {DGP:.3f} (target <0.35: {'PASS' if DGP<0.35 else 'FAIL'})")

df_daylight = pd.DataFrame({"sensor_illuminance_lux": illuminance.round(1)})
df_daylight.to_csv("data/raw/daylight_sensor_grid.csv", index=False)

# ═══════════════════════
# Cell 5 : ML model — ZEB Design Parameter Optimization
# ═════════════
print("\n[Cell 5] ML ZEB Optimization")

np.random.seed(42)
N = 500   # design variants

# Features: window U-value, SHGC, WWR-south, wall U-value, infiltration, PV capacity
X_raw = np.column_stack([
    np.random.uniform(0.8, 3.0, N),    # window_U
    np.random.uniform(0.15, 0.60, N),  # SHGC
    np.random.uniform(0.15, 0.55, N),  # WWR_south
    np.random.uniform(0.10, 0.60, N),  # wall_U
    np.random.uniform(0.05, 0.30, N),  # infiltration_ACH
    np.random.uniform(200, 700, N),    # PV_kWp
])
feat_names = ["window_U","SHGC","WWR_south","wall_U","infiltration","PV_kWp"]

# Target: annual net energy [kWh/m²] — physics-based synthetic label
window_U, SHGC, WWR_s, wall_U, infil, pv_kWp = X_raw.T
EUI_base = 75
y = (EUI_base
     + 8 * window_U
     + 5 * wall_U
     - 12 * SHGC * WWR_s          # solar gain reduction benefit
     + 15 * infil
     - 0.05 * pv_kWp               # PV offsets
     + np.random.normal(0, 3, N))  # noise

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_raw)

rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
gb = GradientBoostingRegressor(n_estimators=200, learning_rate=0.05, random_state=42)

cv = KFold(n_splits=5, shuffle=True, random_state=42)
rf_cv  = cross_val_score(rf, X_scaled, y, cv=cv, scoring='r2')
gb_cv  = cross_val_score(gb, X_scaled, y, cv=cv, scoring='r2')
rf_mae = -cross_val_score(rf, X_scaled, y, cv=cv, scoring='neg_mean_absolute_error')
gb_mae = -cross_val_score(gb, X_scaled, y, cv=cv, scoring='neg_mean_absolute_error')

rf.fit(X_scaled, y)
gb.fit(X_scaled, y)

print(f"  RandomForest  R² = {rf_cv.mean():.3f} ± {rf_cv.std():.3f}  MAE = {rf_mae.mean():.2f} ± {rf_mae.std():.2f} kWh/m²")
print(f"  GradBoost     R² = {gb_cv.mean():.3f} ± {gb_cv.std():.3f}  MAE = {gb_mae.mean():.2f} ± {gb_mae.std():.2f} kWh/m²")

# Feature importance
fi = pd.Series(rf.feature_importances_, index=feat_names).sort_values(ascending=False)
print("\n  Feature importances (RandomForest):")
print(fi.to_string())

df_ml = pd.DataFrame(X_raw, columns=feat_names)
df_ml["net_EUI"] = y
df_ml.to_csv("data/raw/zeb_design_variants.csv", index=False)

# ═══════════════════════════════════════════════════════════════════════════════
# Cell 6 : Statistical Analysis
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[Cell 6] Statistical Analysis")

# Correlation between design params and EUI
corr_matrix = pd.DataFrame(X_raw, columns=feat_names)
corr_matrix["EUI"] = y
corr_vals = corr_matrix.corr()["EUI"].drop("EUI")
print("  Pearson correlation with net EUI:")
print(corr_vals.round(3).to_string())

# One-sample t-test: Is mean net EUI significantly below 50 kWh/m²/yr (NZEB threshold)?
t_stat, p_val = stats.ttest_1samp(y, 50)
print(f"\n  One-sample t-test vs 50 kWh/m²/yr: t={t_stat:.3f}, p={p_val:.4f}")
ci_low, ci_high = stats.t.interval(0.95, len(y)-1, loc=np.mean(y), scale=stats.sem(y))
print(f"  Mean EUI = {y.mean():.2f} ± {y.std():.2f}  95%CI [{ci_low:.2f}, {ci_high:.2f}]")

# ═══════════════════════════════════════════════════════════════════════════════
# Cell 7 : Figures
# ═════════════════════
print("\n[Cell 7] Generating figures…")
sns.set_theme(style="whitegrid", palette="muted")

# ── Fig 1: Monthly energy balance ────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
x = np.arange(12)
w = 0.22
ax.bar(x - 1.5*w, df_thermal.cooling_kWh/1000, w, label="Cooling", color="#e74c3c")
ax.bar(x - 0.5*w, df_thermal.heating_kWh/1000, w, label="Heating", color="#3498db")
ax.bar(x + 0.5*w, df_thermal.lighting_kWh/1000 + df_thermal.equip_kWh/1000, w,
       label="Lighting+Equip", color="#95a5a6")
ax.bar(x + 1.5*w, -df_thermal.pv_gen_kWh/1000, w, label="PV Generation", color="#2ecc71")
ax.axhline(0, color='black', linewidth=0.8)
ax.set_xticks(x); ax.set_xticklabels(month_names, rotation=45, ha='right')
ax.set_xlabel("Month"); ax.set_ylabel("Energy [MWh]")
ax.set_title(f"Monthly Energy Balance\n(Annual net: {annual_net/1000:+.0f} MWh/yr)")
ax.legend(fontsize=8)

# ── Net line energy ───────────────────
ax2 = axes[1]
ax2.plot(month_names, df_thermal.net_kWh/1000, 'o-', color="#e67e22", linewidth=2)
ax2.axhline(0, color='green', linestyle='--', label='ZEB target (net=0)')
ax2.fill_between(range(12), df_thermal.net_kWh/1000, 0,
                 where=(df_thermal.net_kWh > 0), alpha=0.3, color='red', label='Net import')
ax2.fill_between(range(12), df_thermal.net_kWh/1000, 0,
                 where=(df_thermal.net_kWh <= 0), alpha=0.3, color='green', label='Net export')
ax2.set_xticks(range(12)); ax2.set_xticklabels(month_names, rotation=45, ha='right')
ax2.set_xlabel("Month"); ax2.set_ylabel("Net Energy [MWh]")
ax2.set_title("Monthly Net Energy Balance (+ = import, − = export)")
ax2.legend(fontsize=8)
plt.tight_layout()
plt.savefig("figures/fig1_energy_balance.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figures/fig1_energy_balance.png")

# ── Fig 2: CFD ventilation ────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.plot(wind_speeds, Q_wind, 'b-o', label='Wind-driven Q')
ax.axhline(Q_stack, color='r', linestyle='--', label=f'Stack Q = {Q_stack:.2f} m³/s')
ax.set_xlabel("Wind Speed [m/s]"); ax.set_ylabel("Airflow Rate [m³/s]")
ax.set_title("Cross-ventilation Airflow vs Wind Speed")
ax.legend()

ax2 = axes[1]
ax2.plot(wind_speeds, ACH_wind, 'g-s', label='ACH (wind)')
ax2.axhline(ACH_stack, color='orange', linestyle='--', label=f'ACH stack = {ACH_stack:.2f}')
ax2.axhline(4, color='red', linestyle=':', label='ASHRAE comfort threshold (4 ACH)')
ax2.set_xlabel("Wind Speed [m/s]"); ax2.set_ylabel("Air Changes per Hour [ACH]")
ax2.set_title("ACH vs Wind Speed — Comfort Assessment")
ax2.legend()
plt.tight_layout()
plt.savefig("figures/fig2_cfd_ventilation.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figures/fig2_cfd_ventilation.png")

# ── Fig 3: Daylight distribution ─────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.hist(illuminance, bins=40, color='gold', edgecolor='darkorange', alpha=0.8)
ax.axvline(300, color='green', linestyle='--', linewidth=1.5, label='sDA threshold (300 lux)')
ax.axvline(1000, color='red', linestyle='--', linewidth=1.5, label='ASE threshold (1000 lux)')
ax.axvline(cDA, color='blue', linestyle=':', linewidth=1.5, label=f'Median {cDA:.0f} lux')
ax.set_xlabel("Illuminance [lux]"); ax.set_ylabel("Frequency")
ax.set_title(f"Daylight Illuminance Distribution\nsDA={sDA300:.1f}%  ASE={ASE1000:.1f}%")
ax.legend(fontsize=8)

# Sensor heatmap (10×20 grid)
ax2 = axes[1]
grid = illuminance[:200].reshape(10, 20)
im = ax2.imshow(grid, cmap='YlOrRd', aspect='auto', vmin=0, vmax=2000)
plt.colorbar(im, ax=ax2, label='Illuminance [lux]')
ax2.set_xlabel("East-West sensor position"); ax2.set_ylabel("North-South sensor position")
ax2.set_title("Daylight Sensor Grid (Radiance proxy)")
plt.tight_layout()
plt.savefig("figures/fig3_daylight.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figures/fig3_daylight.png")

# ── Fig 4: ML results ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
# Predicted vs actual
y_pred_rf = rf.predict(X_scaled)
ax = axes[0]
ax.scatter(y, y_pred_rf, alpha=0.4, s=20, color='steelblue')
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'r--')
ax.set_xlabel("Actual EUI [kWh/m²/yr]"); ax.set_ylabel("Predicted EUI")
ax.set_title(f"RF: R²={r2_score(y,y_pred_rf):.3f}  MAE={mean_absolute_error(y,y_pred_rf):.2f}")

# Feature importance
ax2 = axes[1]
fi_sorted = fi.sort_values()
fi_sorted.plot(kind='barh', ax=ax2, color='coral')
ax2.set_xlabel("Importance"); ax2.set_title("Feature Importance (RF)")

# CV scores
ax3 = axes[2]
cv_df = pd.DataFrame({"RandomForest": rf_cv, "GradBoost": gb_cv})
cv_df.boxplot(ax=ax3)
ax3.set_ylabel("R² Score")
ax3.set_title(f"5-fold CV R²\nRF={rf_cv.mean():.3f}±{rf_cv.std():.3f}  GB={gb_cv.mean():.3f}±{gb_cv.std():.3f}")
ax3.axhline(0.9, color='green', linestyle='--', linewidth=0.8)
plt.tight_layout()
plt.savefig("figures/fig4_ml_results.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figures/fig4_ml_results.png")

# ── Fig 5: Integrated ZEB Dashboard ──────────────────────────────────────────
fig = plt.figure(figsize=(15, 9))
fig.suptitle("BIM-Integrated Environmental Performance Dashboard — ZEB Office Tokyo",
             fontsize=14, fontweight='bold')

# KPI gauges as colored bars
kpis = {
    "EUI\n[kWh/m²]":     (EUI, 55, 150, "↓ lower better"),
    "sDA\n[%≥300lux]":   (sDA300, 55, 100, "↑ higher better"),
    "ASE\n[%≥1000lux]":  (ASE1000, 0, 10, "↓ lower better"),
    "Net Energy\n[MWh]": (annual_net/1000, -500, 500, "→ target 0"),
    "ACH@3m/s":          (ACH_wind[5], 4, 12, "↑ higher better"),
}
ax_kpi = fig.add_subplot(2, 3, 1)
names = list(kpis.keys())
vals  = [kpis[k][0] for k in names]
lo    = [kpis[k][1] for k in names]
hi    = [kpis[k][2] for k in names]
colors = ["#e74c3c" if v > h else "#2ecc71"
          for v, l, h in zip(vals, lo, hi)]
# Override: EUI and ASE lower is better
if vals[0] <= 55: colors[0] = "#2ecc71"
if vals[2] <= 10: colors[2] = "#2ecc71"
bars = ax_kpi.barh(names, vals, color=colors, edgecolor='gray')
for bar, val in zip(bars, vals):
    ax_kpi.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                f"{val:.1f}", va='center', fontsize=9)
ax_kpi.set_xlabel("Value")
ax_kpi.set_title("KPI Summary")

# Energy pie
ax_pie = fig.add_subplot(2, 3, 2)
pie_vals = [cooling_kWh.sum(), heating_kWh.sum(), lighting_kWh.sum(), equip_kWh.sum()]
pie_lbls = ["Cooling", "Heating", "Lighting", "Equipment"]
ax_pie.pie(pie_vals, labels=pie_lbls, autopct='%1.1f%%',
           colors=["#e74c3c","#3498db","#f39c12","#95a5a6"])
ax_pie.set_title("Annual Energy Breakdown")

# Monthly net
ax_net = fig.add_subplot(2, 3, 3)
colors_bar = ['#e74c3c' if v > 0 else '#2ecc71' for v in df_thermal.net_kWh]
ax_net.bar(month_names, df_thermal.net_kWh/1000, color=colors_bar)
ax_net.axhline(0, color='black', linewidth=0.8)
ax_net.set_xticklabels(month_names, rotation=45, ha='right', fontsize=7)
ax_net.set_ylabel("Net MWh"); ax_net.set_title("Monthly Net Energy")

# Ventilation
ax_vent = fig.add_subplot(2, 3, 4)
ax_vent.plot(wind_speeds, ACH_wind, 'b-o', markersize=5)
ax_vent.axhline(4, color='red', linestyle='--', label='Comfort threshold')
ax_vent.set_xlabel("Wind [m/s]"); ax_vent.set_ylabel("ACH")
ax_vent.set_title("Ventilation Performance")
ax_vent.legend(fontsize=8)

# Daylight histogram (compact)
ax_dl = fig.add_subplot(2, 3, 5)
ax_dl.hist(illuminance, bins=30, color='gold', edgecolor='darkorange', alpha=0.8)
ax_dl.axvline(300, color='green', linestyle='--', linewidth=1.2)
ax_dl.set_xlabel("Illuminance [lux]"); ax_dl.set_ylabel("Count")
ax_dl.set_title(f"Daylighting  sDA={sDA300:.1f}%")

# ML feature importance
ax_ml = fig.add_subplot(2, 3, 6)
fi.sort_values().plot(kind='barh', ax=ax_ml, color='steelblue')
ax_ml.set_xlabel("Importance"); ax_ml.set_title(f"ML Feature Importance\nRF R²={rf_cv.mean():.3f}")

plt.tight_layout()
plt.savefig("figures/fig5_dashboard.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figures/fig5_dashboard.png")

# ═════════════════════════════════════════════════════════════════════
# Cell 8 : pip freeze (environment record)
# ════════════════════════════════════════════════════════════════════════
import subprocess
result = subprocess.run(["pip", "freeze"], capture_output=True, text=True)
with open("data/raw/pip_freeze.txt", "w") as f:
    f.write(result.stdout)
print("\n[Cell 8] pip freeze saved to data/raw/pip_freeze.txt")
key_pkgs = [l for l in result.stdout.splitlines()
            if any(p in l.lower() for p in ["numpy","pandas","scikit","scipy","matplotlib","seaborn"])]
print("\n  Key packages:")
for p in key_pkgs:
    print(f"    {p}")

print("\n✅ All cells completed successfully.")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary dict for paper
# ════════════════════
summary = {
    "annual_demand_kWh": float(annual_demand),
    "annual_pv_kWh": float(annual_pv),
    "annual_net_kWh": float(annual_net),
    "EUI_kWh_m2_yr": float(EUI),
    "Q_stack_m3_s": float(Q_stack),
    "ACH_stack": float(ACH_stack),
    "ACH_at_3ms": float(ACH_wind[5]),
    "sDA300_pct": float(sDA300),
    "ASE1000_pct": float(ASE1000),
    "median_lux": float(cDA),
    "UGR": float(UGR),
    "DGP": float(DGP),
    "RF_R2_mean": float(rf_cv.mean()),
    "RF_R2_std": float(rf_cv.std()),
    "RF_MAE_mean": float(rf_mae.mean()),
    "RF_MAE_std": float(rf_mae.std()),
    "GB_R2_mean": float(gb_cv.mean()),
    "GB_R2_std": float(gb_cv.std()),
    "GB_MAE_mean": float(gb_mae.mean()),
    "GB_MAE_std": float(gb_mae.std()),
    "t_stat": float(t_stat),
    "p_val": float(p_val),
    "EUI_mean": float(y.mean()),
    "EUI_std": float(y.std()),
    "EUI_CI_low": float(ci_low),
    "EUI_CI_high": float(ci_high),
    "feature_importance": fi.to_dict(),
    "corr_EUI": corr_vals.to_dict(),
}
with open("data/raw/results_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nResults summary saved to data/raw/results_summary.json")
import pprint; pprint.pprint(summary)
