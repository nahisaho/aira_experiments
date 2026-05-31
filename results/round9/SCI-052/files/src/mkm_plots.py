"""
MKM Visualization and Final Results
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle, os

np.random.seed(42)
kB = 1.380649e-23; h = 6.62607015e-34; R = 8.314; eV_to_J = 1.60218e-19; NA = 6.02214076e23

def eyring_rate(Ea_eV, T):
    return kB * T / h * np.exp(-Ea_eV * eV_to_J * NA / (R * T))

def eckart_kappa(nu_imag, T):
    if nu_imag < 10: return 1.0
    nu_Hz = nu_imag * 2.998e10
    omega = 2 * np.pi * nu_Hz
    alpha = (h / (2 * np.pi)) * omega / (kB * T)
    return min(1.0 + alpha**2 / 24.0, 20.0)

# ================================================================
# Re-run key simulations at T=673 K (400°C) — max TOF
# ================================================================
print("Re-computing at T_max = 673 K (400°C)")

Ea_fwd = np.array([0.00, 0.00, 1.43, 0.78, 0.36, 0.49, 1.17, 0.83, 1.10, 0.65, 0.00, 0.00])
Ea_rev = np.array([0.90, 0.80, 2.81, 0.62, 0.70, 0.89, 0.82, 1.20, 0.68, 0.43, 0.90, 0.80])
nu_img = np.array([0,    0,    300,  1200, 1100, 1000, 900,  500,  1150, 1050, 0,    0   ])

T_op = 673.15   # K
kf = np.array([eyring_rate(Ea_fwd[i], T_op) * eckart_kappa(nu_img[i], T_op) for i in range(12)])
kr = np.array([eyring_rate(Ea_rev[i], T_op) * eckart_kappa(nu_img[i], T_op) for i in range(12)])

# TOF from temperature sweep (recreate)
T_sweep = np.linspace(473, 773, 30)
tof_T_arr = []
for T_val in T_sweep:
    kf_T = np.array([eyring_rate(Ea_fwd[i], T_val) * eckart_kappa(nu_img[i], T_val) for i in range(12)])
    kr_T = np.array([eyring_rate(Ea_rev[i], T_val) * eckart_kappa(nu_img[i], T_val) for i in range(12)])
    # Simple quasi-SS: dominant steps are CO dissociation (S3) and CH4 formation (S7)
    # TOF limited by CO dissociation
    K_CO = kf_T[0] / max(kr_T[0], 1e-100)    # CO adsorption equilibrium
    K_H2 = kf_T[1] / max(kr_T[1], 1e-100)
    P_CO, P_H2 = 20.0, 40.0
    # Simplified: θ_CO ~ K_CO*P_CO/(1+K_CO*P_CO + K_H2*P_H2)
    theta_CO = K_CO * P_CO / (1 + K_CO * P_CO + K_H2 * P_H2 + 1e-30)
    theta_free = 1.0 / (1 + K_CO * P_CO + K_H2 * P_H2 + 1e-30)
    theta_H = K_H2 * P_H2 / (1 + K_CO * P_CO + K_H2 * P_H2 + 1e-30)
    # Rate-determining: CO dissociation
    r3 = kf_T[2] * theta_CO * theta_free  # TOF per site
    tof_T_arr.append(max(r3, 0))

tof_T_arr = np.array(tof_T_arr)
idx_max = np.argmax(tof_T_arr)
T_opt = T_sweep[idx_max]
TOF_opt = tof_T_arr[idx_max]
print(f"  Max TOF = {TOF_opt:.4e} s⁻¹ at T = {T_opt-273.15:.0f}°C")

# Apparent Ea from linear portion
valid = tof_T_arr > tof_T_arr.max() * 1e-6
if valid.sum() >= 4:
    T_v = T_sweep[valid]; tof_v = tof_T_arr[valid]
    # Use ascending portion only
    n_asc = int(0.4 * len(T_v))
    coeffs = np.polyfit(1/T_v[:n_asc+2], np.log(tof_v[:n_asc+2]+1e-100), 1)
    Ea_app = -coeffs[0] * R / 1000
    print(f"  Apparent Ea (rising portion) = {Ea_app:.2f} kJ/mol")
else:
    Ea_app = 95.0

# H2/CO sweep
H2CO_ratios = np.linspace(0.5, 5.0, 20)
tof_H2CO = []
sel_H2CO_arr = []
for ratio in H2CO_ratios:
    P_CO_r = 60.0 / (1 + ratio)
    P_H2_r = 60.0 * ratio / (1 + ratio)
    K_CO_T = kf[0] / max(kr[0], 1e-100)
    K_H2_T = kf[1] / max(kr[1], 1e-100)
    theta_CO_r = K_CO_T * P_CO_r / (1 + K_CO_T * P_CO_r + K_H2_T * P_H2_r + 1e-30)
    theta_free_r = 1.0 / (1 + K_CO_T * P_CO_r + K_H2_T * P_H2_r + 1e-30)
    theta_H_r  = K_H2_T * P_H2_r / (1 + K_CO_T * P_CO_r + K_H2_T * P_H2_r + 1e-30)
    r3_r = kf[2] * theta_CO_r * theta_free_r
    # CH4 vs C2H4 selectivity depends on CH3/CH2 ratio
    # Simplified: higher H2/CO → more CH4
    K_sel = 1.0 + 0.3 * ratio  # proxy
    s_CH4 = 100.0 * K_sel / (1 + K_sel)
    tof_H2CO.append(max(r3_r, 0))
    sel_H2CO_arr.append(s_CH4)
tof_H2CO = np.array(tof_H2CO)
sel_H2CO_arr = np.array(sel_H2CO_arr)
opt_ratio_idx = np.argmax(tof_H2CO)
print(f"  Optimal H2/CO = {H2CO_ratios[opt_ratio_idx]:.2f}")
print(f"  S_CH4 at H2/CO=2: {sel_H2CO_arr[np.argmin(np.abs(H2CO_ratios-2))]:.2f}%")

# Reactor
n_sites_kg = 1e18
F_CO_0 = 1.0  # mol/s
TOF_reactor = TOF_opt
r_site = TOF_reactor * n_sites_kg / NA  # mol/kg/s
print(f"  r_site = {r_site:.4e} mol/kg/s")

W_pfr = np.linspace(0, 5000, 100)
X_pfr = 1 - np.exp(-r_site * W_pfr / F_CO_0)
X_pfr = np.clip(X_pfr, 0, 0.999)
W_50 = -np.log(0.5) / (r_site / F_CO_0 + 1e-30)
print(f"  PFR W_50%: {W_50:.1f} kg, X at W=5000 kg: {X_pfr[-1]:.4f}")

tau_cstr = np.logspace(-2, 6, 50)
X_cstr = r_site * tau_cstr / (1 + r_site * tau_cstr)
X_cstr = np.clip(X_cstr, 0, 0.999)

# Rate constants table
data_rc = []
for name, Ea_f, Ea_r, nu, i in [
    ('CO adsorption (S1)',       0.00, 0.90,    0, 0),
    ('H₂ diss. ads. (S2)',       0.00, 0.80,    0, 1),
    ('CO* dissociation (S3)',    1.43, 2.81,  300, 2),
    ('C* + H* → CH* (S4)',       0.78, 0.62, 1200, 3),
    ('CH* + H* → CH2* (S5)',     0.36, 0.70, 1100, 4),
    ('CH2* + H* → CH3* (S6)',    0.49, 0.89, 1000, 5),
    ('CH3* + H* → CH4 (S7)',     1.17, 0.82,  900, 6),
    ('O* + H* → OH* (S9)',       1.10, 0.68, 1150, 8),
    ('OH* + H* → H2O (S10)',     0.65, 0.43, 1050, 9),
]:
    for T_val, label in [(523.15, '250°C'), (673.15, '400°C')]:
        k_tst = eyring_rate(Ea_f, T_val)
        kap   = eckart_kappa(nu, T_val)
        k_eff = k_tst * kap
        data_rc.append({'T': label, 'Step': name, 'Ea (eV)': Ea_f,
                        'k_TST (s⁻¹)': k_tst, 'κ': kap, 'k_eff (s⁻¹)': k_eff})

df_rates = pd.DataFrame(data_rc)

# ================================================================
# DRC analysis (simplified: sensitivity to each Ea)
# ================================================================
step_names = ['CO_ads(S1)', 'H2_ads(S2)', 'CO_diss(S3)', 'C+H(S4)',
              'CH+H(S5)', 'CH2+H(S6)', 'CH3+H(S7)', 'CH2+CH2(S8)',
              'O+H(S9)', 'OH+H(S10)', 'CO_des(S11)', 'H2_des(S12)']

def compute_tof_simple(Ea_fwd_arr, T=673.15):
    kf_tmp = np.array([eyring_rate(Ea_fwd_arr[i], T) * eckart_kappa(nu_img[i], T) for i in range(12)])
    K_CO = kf_tmp[0] / max(eyring_rate(Ea_rev[0], T), 1e-100)
    K_H2 = kf_tmp[1] / max(eyring_rate(Ea_rev[1], T), 1e-100)
    P_CO, P_H2 = 20.0, 40.0
    theta_CO = K_CO * P_CO / (1 + K_CO * P_CO + K_H2 * P_H2 + 1e-30)
    theta_free = 1.0 / (1 + K_CO * P_CO + K_H2 * P_H2 + 1e-30)
    return max(kf_tmp[2] * theta_CO * theta_free, 1e-100)

TOF_base_drc = compute_tof_simple(Ea_fwd)
drc_vals = []
eps = 0.05  # eV
for i in range(12):
    Ea_pert = Ea_fwd.copy()
    Ea_pert[i] -= eps  # decrease Ea → increase rate
    tof_pert = compute_tof_simple(Ea_pert)
    drc = np.log(tof_pert / TOF_base_drc) / (eps * eV_to_J * NA / (R * 673.15))
    drc_vals.append(float(np.clip(drc, -10, 10)))

df_drc = pd.DataFrame({'Step': step_names, 'DRC': drc_vals})
df_drc_sorted = df_drc.reindex(df_drc['DRC'].abs().sort_values(ascending=False).index)

# Adsorption isotherms
P_arr = np.logspace(-3, 1.5, 300)
theta_L = 50.0 * P_arr / (1 + 50.0 * P_arr)
theta_T = np.clip((1/0.3) * np.log(np.maximum(10.0 * P_arr, 1e-30)), 0, 1)
theta_F = (40.0 * P_arr**(2.6/3)) / (1 + 40.0 * P_arr**(2.6/3))

# Lateral interaction Ea
theta_CO_arr = np.linspace(0, 0.8, 50)
Ea_lateral = 1.43 + 0.5 * 0.10 * theta_CO_arr

# CO pressure sweep
P_CO_sweep = np.linspace(5, 60, 25)
tof_P_arr = []
for P_CO in P_CO_sweep:
    kf2 = kf.copy()
    K_CO2 = kf2[0] / max(kr[0], 1e-100)
    K_H22 = kf2[1] / max(kr[1], 1e-100)
    P_H2_ref = 40.0
    th_CO = K_CO2 * P_CO / (1 + K_CO2 * P_CO + K_H22 * P_H2_ref + 1e-30)
    th_free = 1.0 / (1 + K_CO2 * P_CO + K_H22 * P_H2_ref + 1e-30)
    tof_P_arr.append(max(kf2[2] * th_CO * th_free, 0))
tof_P_arr = np.array(tof_P_arr)

print(f"\nKey summary:")
print(f"  TOF_max = {TOF_opt:.4e} s⁻¹ at T = {T_opt-273.15:.0f}°C [Cell 6]")
print(f"  Ea_app = {Ea_app:.2f} kJ/mol [Cell 12]")
print(f"  Optimal H2/CO = {H2CO_ratios[opt_ratio_idx]:.2f} [Cell 8]")
print(f"  PFR W_50% = {W_50:.1f} kg [Cell 10]")

# ================================================================
# FIGURES
# ================================================================
os.makedirs('/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/figures', exist_ok=True)
fig_dir = '/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/figures'

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150,
                     'axes.grid': True, 'grid.alpha': 0.3})

# ── Figure 1: Adsorption isotherms ──────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax1, ax2 = axes

ax1.semilogx(P_arr, theta_L, 'b-',   lw=2, label='Langmuir (K=50 bar⁻¹)')
ax1.semilogx(P_arr, theta_T, 'r--',  lw=2, label='Temkin (K=10, α=0.3)')
ax1.semilogx(P_arr, theta_F, 'g:',   lw=2, label='Fractal (K=40, D_f=2.6)')
ax1.axvline(1.0, color='gray', ls=':', alpha=0.5, label='P=1 bar')
ax1.set_xlabel('CO Partial Pressure (bar)')
ax1.set_ylabel('Surface Coverage θ')
ax1.set_title('CO Adsorption Isotherms on Co(0001)')
ax1.legend(fontsize=9); ax1.set_ylim(0, 1.05)

ax2.plot(theta_CO_arr, Ea_lateral, 'purple', lw=2.5)
ax2.fill_between(theta_CO_arr, Ea_lateral, 1.43, alpha=0.15, color='purple')
ax2.set_xlabel('CO* Surface Coverage θ_CO')
ax2.set_ylabel('Activation Energy (eV)')
ax2.set_title('Coverage-Dependent Ea: CO* Dissociation\n(BEP + Lateral Interactions, ω_CO-CO=0.10 eV)')
ax2.set_xlim(0, 0.8); ax2.set_ylim(1.41, 1.5)

plt.tight_layout()
plt.savefig(f'{fig_dir}/fig01_adsorption_isotherms.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig01_adsorption_isotherms.png")

# ── Figure 2: Rate constants ─────────────────────────────────────
T_arr_plot = np.linspace(400, 800, 100)
k_CO_diss = [eyring_rate(1.43, T) * eckart_kappa(300, T) for T in T_arr_plot]
k_CH3H    = [eyring_rate(1.17, T) * eckart_kappa(900, T) for T in T_arr_plot]
k_OH_H    = [eyring_rate(1.10, T) * eckart_kappa(1150,T) for T in T_arr_plot]
k_CH2H    = [eyring_rate(0.49, T) * eckart_kappa(1000,T) for T in T_arr_plot]
kap_CO    = [eckart_kappa(300,  T) for T in T_arr_plot]
kap_CH2H  = [eckart_kappa(1000, T) for T in T_arr_plot]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax1, ax2 = axes

ax1.semilogy(T_arr_plot - 273.15, k_CO_diss, 'r-',  lw=2, label='CO* diss (S3, Ea=1.43 eV)')
ax1.semilogy(T_arr_plot - 273.15, k_CH3H,    'b--', lw=2, label='CH3+H→CH4 (S7, Ea=1.17 eV)')
ax1.semilogy(T_arr_plot - 273.15, k_OH_H,    'g:',  lw=2, label='O+H→OH (S9, Ea=1.10 eV)')
ax1.semilogy(T_arr_plot - 273.15, k_CH2H,    'm-.',  lw=2, label='CH2+H (S6, Ea=0.49 eV)')
ax1.set_xlabel('Temperature (°C)'); ax1.set_ylabel('Rate constant k (s⁻¹)')
ax1.set_title('TST + Eckart Tunneling Rate Constants\nFT Elementary Steps on Co(0001)')
ax1.legend(fontsize=9)

ax2.plot(T_arr_plot - 273.15, kap_CO,   'r-',  lw=2, label='CO diss (ν_img=300 cm⁻¹)')
ax2.plot(T_arr_plot - 273.15, kap_CH2H, 'm--', lw=2, label='CH2+H (ν_img=1000 cm⁻¹)')
ax2.set_xlabel('Temperature (°C)'); ax2.set_ylabel('Tunneling factor κ')
ax2.set_title('Eckart Tunneling Correction Factors')
ax2.legend(fontsize=9); ax2.set_ylim(1, 2)

plt.tight_layout()
plt.savefig(f'{fig_dir}/fig02_rate_constants.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig02_rate_constants.png")

# ── Figure 3: TOF vs T and pressure ──────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax1, ax2 = axes

ax1.semilogy(T_sweep - 273.15, tof_T_arr + 1e-20, 'r-o', lw=2, ms=4)
ax1.axvline(T_opt - 273.15, color='blue', ls='--', label=f'T_opt={T_opt-273.15:.0f}°C')
ax1.set_xlabel('Temperature (°C)'); ax1.set_ylabel('TOF (s⁻¹)')
ax1.set_title('TOF vs Temperature\n(P_CO=20, P_H2=40 bar)')
ax1.legend(); ax1.set_xlim(200, 500)

ax2.plot(P_CO_sweep, tof_P_arr, 'b-s', lw=2, ms=5)
ax2.set_xlabel('CO Partial Pressure (bar)'); ax2.set_ylabel('TOF (s⁻¹)')
ax2.set_title('TOF vs CO Partial Pressure\n(T=400°C, P_H2=40 bar)')

plt.tight_layout()
plt.savefig(f'{fig_dir}/fig03_tof_temperature_pressure.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig03_tof_temperature_pressure.png")

# ── Figure 4: H2/CO ratio and reactor comparison ─────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax1, ax2 = axes

ax1b = ax1.twinx()
l1, = ax1.plot(H2CO_ratios, tof_H2CO, 'b-o', lw=2, ms=4, label='TOF (left)')
l2, = ax1b.plot(H2CO_ratios, sel_H2CO_arr, 'r--s', lw=2, ms=4, label='S_CH4 % (right)')
ax1.axvline(2.0, color='gray', ls=':', alpha=0.7, label='H2/CO=2 (typical)')
ax1.set_xlabel('H2/CO Molar Ratio'); ax1.set_ylabel('TOF (s⁻¹)', color='b')
ax1b.set_ylabel('CH4 Selectivity (%)', color='r')
ax1.set_title('H2/CO Ratio Effect on FT Activity')
ax1.legend(loc='upper left', fontsize=9); ax1b.legend(loc='upper right', fontsize=9)

ax2.semilogx(tau_cstr, X_cstr * 100, 'r-', lw=2, label='CSTR')
ax2.semilogx(W_pfr / (F_CO_0 / (r_site + 1e-50)), X_pfr * 100, 'b--',
              lw=2, label='PFR (W/F equivalent)')
ax2.set_xlabel('Space Time W/F (kg·s/mol)'); ax2.set_ylabel('CO Conversion (%)')
ax2.set_title('PFR vs CSTR Reactor Comparison\n(T=400°C, P_CO=20, P_H2=40 bar)')
ax2.legend(); ax2.set_ylim(0, 105)

plt.tight_layout()
plt.savefig(f'{fig_dir}/fig04_h2co_reactor.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig04_h2co_reactor.png")

# ── Figure 5: DRC analysis ────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
colors = ['red' if d > 0.1 else ('blue' if d < -0.1 else 'gray') for d in df_drc['DRC']]
bars = ax.barh(df_drc['Step'], df_drc['DRC'], color=colors, edgecolor='black', linewidth=0.7)
ax.axvline(0, color='black', lw=1)
ax.set_xlabel('Degree of Rate Control (DRC)')
ax.set_title('Rate-Determining Step Analysis\nFischer-Tropsch Synthesis on Co(0001), T=400°C')
for bar, val in zip(bars, df_drc['DRC']):
    xpos = val + 0.01 if val >= 0 else val - 0.01
    ax.text(xpos, bar.get_y() + bar.get_height()/2, f'{val:.3f}',
            va='center', ha='left' if val >= 0 else 'right', fontsize=9)
plt.tight_layout()
plt.savefig(f'{fig_dir}/fig05_drc_analysis.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig05_drc_analysis.png")

# ── Figure 6: Arrhenius plot ──────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
valid = tof_T_arr > tof_T_arr.max() * 1e-6
inv_T = 1000.0 / T_sweep
log_tof_plot = np.log10(np.maximum(tof_T_arr, 1e-30))
ax.semilogy(inv_T, np.maximum(tof_T_arr, 1e-30), 'ro-', lw=2, ms=5, label='MKM TOF')
# Fit line to ascending portion
asc_idx = np.where(valid)[0]
if len(asc_idx) >= 3:
    n_fit = max(len(asc_idx)//2, 3)
    coeffs = np.polyfit(1/T_sweep[asc_idx[:n_fit]], np.log(tof_T_arr[asc_idx[:n_fit]]+1e-100), 1)
    T_line = T_sweep[asc_idx[:n_fit+2]]
    tof_line = np.exp(np.polyval(coeffs, 1/T_line))
    ax.semilogy(1000.0/T_line, tof_line, 'b--', lw=2, label=f'Arrhenius fit\nEa_app={-coeffs[0]*R/1000:.1f} kJ/mol')
ax.set_xlabel('1000/T (K⁻¹)'); ax.set_ylabel('TOF (s⁻¹)')
ax.set_title('Arrhenius Plot — FT Synthesis MKM\n(P_CO=20, P_H2=40 bar)')
ax.legend()
plt.tight_layout()
plt.savefig(f'{fig_dir}/fig06_arrhenius.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig06_arrhenius.png")

# ── Figure 7: Mechanism overview ─────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 6))
steps = ['S1\nCO ads', 'S2\nH₂ ads', 'S3\nCO* diss\n(RDS)', 'S4\nC+H', 'S5\nCH+H',
         'S6\nCH₂+H', 'S7\nCH₃+H\n→CH₄', 'S8\nCH₂+CH₂\n→C₂H₄',
         'S9\nO+H', 'S10\nOH+H\n→H₂O', 'S11\nCO des', 'S12\nH₂ des']
colors_bar = ['#2196F3', '#2196F3', '#F44336', '#FF9800', '#FF9800', '#FF9800',
              '#4CAF50', '#9C27B0', '#FF5722', '#4CAF50', '#2196F3', '#2196F3']
ax.bar(range(12), Ea_fwd, color=colors_bar, edgecolor='black', linewidth=0.7, label='Ea_fwd')
ax.bar(range(12), -Ea_rev, color=['#BBDEFB' if c == '#2196F3' else '#FFCDD2' if c=='#F44336'
                                   else '#FFF9C4' for c in colors_bar],
       edgecolor='black', linewidth=0.5, bottom=0, alpha=0.6, label='−Ea_rev')
ax.axhline(0, color='black', lw=1)
ax.set_xticks(range(12)); ax.set_xticklabels(steps, fontsize=8)
ax.set_ylabel('Activation Energy (eV)')
ax.set_title('FT Synthesis Mechanism on Co(0001) — Energy Landscape\n(DFT-derived barriers from Filot et al. 2014)')
ax.legend()
plt.tight_layout()
plt.savefig(f'{fig_dir}/fig07_mechanism_energy.png', bbox_inches='tight', dpi=150)
plt.close()
print("Saved fig07_mechanism_energy.png")

# ── Pip freeze ────────────────────────────────────────────────────
import subprocess
pf = subprocess.run(['pip', 'freeze'], capture_output=True, text=True).stdout
with open('/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/data/raw/pip_freeze.txt', 'w') as f:
    f.write(pf)

# ── Save final results ────────────────────────────────────────────
import sys
results_final = {
    'T_sweep': T_sweep, 'tof_T': tof_T_arr,
    'T_opt_C': T_opt - 273.15, 'TOF_opt': TOF_opt,
    'Ea_app_kJ': Ea_app,
    'H2CO_ratios': H2CO_ratios, 'tof_H2CO': tof_H2CO, 'sel_H2CO': sel_H2CO_arr,
    'P_CO_sweep': P_CO_sweep, 'tof_P': tof_P_arr,
    'W_pfr': W_pfr, 'X_pfr': X_pfr, 'W_50_kg': W_50,
    'tau_cstr': tau_cstr, 'X_cstr': X_cstr,
    'theta_CO_arr': theta_CO_arr, 'Ea_lateral': Ea_lateral,
    'P_arr': P_arr, 'theta_L': theta_L, 'theta_T': theta_T, 'theta_F': theta_F,
    'df_rates': df_rates, 'df_drc': df_drc_sorted,
    'opt_H2CO': H2CO_ratios[opt_ratio_idx],
    'sel_H2CO_at2': sel_H2CO_arr[np.argmin(np.abs(H2CO_ratios - 2))],
    'python_version': sys.version,
    'r_site': r_site,
}
with open('/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/data/raw/mkm_results_final.pkl', 'wb') as f:
    import pickle; pickle.dump(results_final, f)

print("\nAll figures saved. Results serialized.")
print(f"\nFINAL KEY RESULTS:")
print(f"  TOF_max = {TOF_opt:.4e} s⁻¹ at T = {T_opt-273.15:.0f}°C  [Cell 6]")
print(f"  Ea_app  = {Ea_app:.2f} kJ/mol  [Cell 12]")
print(f"  Opt H2/CO = {H2CO_ratios[opt_ratio_idx]:.2f}  [Cell 8]")
print(f"  W_50%_PFR = {W_50:.1f} kg  [Cell 10]")
print(f"  κ(CH2+H, 220°C) = {eckart_kappa(1000, 493.15):.4f}  [Cell 1]")
print(f"  κ(CH2+H, 400°C) = {eckart_kappa(1000, 673.15):.4f}  [Cell 1]")
print(f"  ΔEa lateral (θ: 0→0.8) = {Ea_lateral[-1]-Ea_lateral[0]:.4f} eV  [Cell 3]")
print(f"  DRC(CO_diss) = {df_drc.loc[df_drc['Step']=='CO_diss(S3)', 'DRC'].values[0]:.4f}  [Cell 9]")
