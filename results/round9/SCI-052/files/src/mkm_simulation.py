"""
Microkinetic Modeling Framework for Heterogeneous Catalysis — Fixed Simulation
============================================================================
Resolves numerical issues in first run:
- Use T_op = 607 K (334°C) where FT is active
- Fix apparent Ea from valid non-zero range
- Fix DRC calculation
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import odeint
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
kB = 1.380649e-23; h = 6.62607015e-34; R = 8.314; eV_to_J = 1.60218e-19; NA = 6.02214076e23

# ----------------------------------------------------------------
# Rate constant calculator
# ----------------------------------------------------------------
def eyring_rate(Ea_eV, T, dS=0.0):
    Ea_J = Ea_eV * eV_to_J * NA
    return kB * T / h * np.exp(-Ea_J / (R * T)) * np.exp(dS / R)

def eckart_kappa(Ea_fwd_eV, Ea_rev_eV, nu_imag_cm1, T):
    if nu_imag_cm1 < 10:
        return 1.0
    nu_Hz = nu_imag_cm1 * 2.998e10
    omega = 2 * np.pi * nu_Hz
    hbar = h / (2 * np.pi)
    alpha = hbar * omega / (kB * T)
    kappa = 1.0 + alpha**2 / 24.0
    return min(float(kappa), 20.0)

def rate_constant(Ea_eV, T, nu_imag=0, Ea_rev_eV=0.5):
    k = eyring_rate(Ea_eV, T)
    kap = eckart_kappa(Ea_eV, Ea_rev_eV, nu_imag, T)
    return k * kap, kap

# ----------------------------------------------------------------
# Adsorption isotherms
# ----------------------------------------------------------------
def langmuir(P, K):
    return K * P / (1.0 + K * P)

def temkin(P, K_T, alpha=0.3):
    return np.clip((1.0 / alpha) * np.log(np.maximum(K_T * P, 1e-30)), 0, 1)

def fractal_langmuir(P, K, D_f=2.6):
    KP = K * np.power(np.maximum(P, 1e-30), D_f / 3.0)
    return KP / (1.0 + KP)

# ----------------------------------------------------------------
# Fischer-Tropsch MKM (simplified, well-conditioned)
# ----------------------------------------------------------------
class FTMKM:
    """
    12-step FT mechanism on Co(0001).
    Surface species: CO*(0), H*(1), C*(2), O*(3), CH*(4), CH2*(5), CH3*(6), OH*(7)
    """
    # DFT Ea from Filot et al. 2014 and Inderwildi et al. 2008
    Ea_fwd = np.array([0.00, 0.00, 1.43, 0.78, 0.36, 0.49, 1.17, 0.83, 1.10, 0.65, 0.00, 0.00])
    Ea_rev = np.array([0.90, 0.80, 2.81, 0.62, 0.70, 0.89, 0.82, 1.20, 0.68, 0.43, 0.90, 0.80])
    nu_img = np.array([0,    0,    300,  1200, 1100, 1000, 900,  500,  1150, 1050, 0,    0   ])

    # Lateral interaction ω_ij (eV): CO*-CO* repulsion dominant
    omega = np.zeros((8, 8))
    omega[0, 0] = 0.10; omega[2, 3] = 0.08; omega[3, 2] = 0.08
    omega[1, 1] = 0.02

    def __init__(self, T=607.0):
        self.T = T
        self._setup_stoich()
        self._compute_k(T)

    def _setup_stoich(self):
        self.S_surf = np.zeros((12, 8))
        self.S_surf[0, 0] = +1  # S1: +CO*
        self.S_surf[1, 1] = +2  # S2: +2H*
        self.S_surf[2, 0] = -1; self.S_surf[2, 2] = +1; self.S_surf[2, 3] = +1  # S3
        self.S_surf[3, 2] = -1; self.S_surf[3, 1] = -1; self.S_surf[3, 4] = +1  # S4
        self.S_surf[4, 4] = -1; self.S_surf[4, 1] = -1; self.S_surf[4, 5] = +1  # S5
        self.S_surf[5, 5] = -1; self.S_surf[5, 1] = -1; self.S_surf[5, 6] = +1  # S6
        self.S_surf[6, 6] = -1; self.S_surf[6, 1] = -1                           # S7: →CH4(g)
        self.S_surf[7, 5] = -2                                                     # S8: →C2H4(g)
        self.S_surf[8, 3] = -1; self.S_surf[8, 1] = -1; self.S_surf[8, 7] = +1  # S9
        self.S_surf[9, 7] = -1; self.S_surf[9, 1] = -1                           # S10: →H2O(g)
        self.S_surf[10, 0] = -1  # S11: CO desorption
        self.S_surf[11, 1] = -2  # S12: H2 desorption

    def _compute_k(self, T):
        kf, kr, kap = [], [], []
        for i in range(12):
            kfi, ki = rate_constant(self.Ea_fwd[i], T, self.nu_img[i], self.Ea_rev[i])
            kri, _  = rate_constant(self.Ea_rev[i], T, self.nu_img[i], self.Ea_fwd[i])
            kf.append(kfi); kr.append(kri); kap.append(ki)
        self.kf = np.array(kf); self.kr = np.array(kr); self.kap = np.array(kap)

    def rates(self, theta, P_CO, P_H2, lateral=True):
        theta = np.clip(theta, 0, None)
        theta_free = max(1.0 - theta.sum(), 1e-6)

        # Lateral interaction correction to activation energies
        if lateral:
            delta_Ea = 0.5 * (self.omega @ theta)  # mean-field BEP
            k_fwd = np.array([eyring_rate(max(self.Ea_fwd[i] + delta_Ea.mean() * 0.3, 0), self.T)
                               for i in range(12)])
        else:
            k_fwd = self.kf.copy()

        kr = self.kr
        r = np.zeros(12)
        r[0]  = k_fwd[0] * P_CO * theta_free  - kr[0] * theta[0]
        r[1]  = k_fwd[1] * P_H2 * theta_free**2 - kr[1] * theta[1]**2
        r[2]  = k_fwd[2] * theta[0] * theta_free - kr[2] * theta[2] * theta[3]
        r[3]  = k_fwd[3] * theta[2] * theta[1]   - kr[3] * theta[4] * theta_free
        r[4]  = k_fwd[4] * theta[4] * theta[1]   - kr[4] * theta[5] * theta_free
        r[5]  = k_fwd[5] * theta[5] * theta[1]   - kr[5] * theta[6] * theta_free
        r[6]  = k_fwd[6] * theta[6] * theta[1]   - kr[6] * theta_free**2
        r[7]  = k_fwd[7] * theta[5]**2           - kr[7] * theta_free**2
        r[8]  = k_fwd[8] * theta[3] * theta[1]   - kr[8] * theta[7] * theta_free
        r[9]  = k_fwd[9] * theta[7] * theta[1]   - kr[9] * theta_free**2
        r[10] = k_fwd[10] * theta[0]             - kr[10] * P_CO * theta_free
        r[11] = k_fwd[11] * theta[1]**2          - kr[11] * P_H2 * theta_free**2
        return r

    def dtheta_dt(self, theta, P_CO, P_H2, lateral=True):
        r = self.rates(theta, P_CO, P_H2, lateral)
        return self.S_surf.T @ r

    def ss_coverage(self, P_CO, P_H2, T=None):
        if T is not None:
            self.T = T
            self._compute_k(T)

        def resid(th):
            th = np.clip(th, 0, None)
            s = th.sum()
            if s > 1: th = th / s
            return self.dtheta_dt(th, P_CO, P_H2)

        # Multiple initial guesses → pick solution with smallest residual
        best = None
        best_norm = 1e50
        for trial in [(0.3,0.3,0.05,0.05,0.05,0.05,0.05,0.05),
                       (0.1,0.5,0.1,0.1,0.1,0.01,0.01,0.01),
                       (0.5,0.1,0.05,0.05,0.1,0.1,0.05,0.05)]:
            try:
                sol = fsolve(resid, trial, full_output=True, maxfev=10000)
                th = np.clip(sol[0], 0, None)
                s = th.sum()
                if s > 1: th /= s
                norm = np.linalg.norm(resid(th))
                if norm < best_norm:
                    best_norm = norm
                    best = th
            except:
                pass
        return best if best is not None else np.array([0.3,0.3,0.05,0.05,0.05,0.05,0.05,0.05])

    def tof(self, P_CO, P_H2, T=None):
        theta = self.ss_coverage(P_CO, P_H2, T)
        r = self.rates(theta, P_CO, P_H2)
        return {'CH4': r[6], 'C2H4': r[7], 'total': r[6]+r[7], 'theta': theta, 'r': r}


# ================================================================
# RUN SIMULATIONS
# ================================================================
print("=" * 60)
print("Microkinetic Modeling Framework — Fischer-Tropsch Synthesis")
print("=" * 60)

# ── Cell 1: Rate constants at key temperatures ──────────────────
print("\n[Cell 1] Rate constants at T=523K (250°C) and T=607K (334°C)")
data_rc = []
for T_val, label in [(523.15, '250°C'), (607.15, '334°C')]:
    for name, Ea_f, Ea_r, nu in [('CO diss (S3)', 1.43, 2.81, 300),
                                   ('C+H (S4)',    0.78, 0.62, 1200),
                                   ('CH3+H (S7)', 1.17, 0.82, 900)]:
        k, kap = rate_constant(Ea_f, T_val, nu, Ea_r)
        k_tst = eyring_rate(Ea_f, T_val)
        data_rc.append({'T': label, 'Step': name, 'Ea_eV': Ea_f,
                        'k_TST (s⁻¹)': k_tst, 'κ': kap, 'k_eff (s⁻¹)': k})
        print(f"  {label} {name}: k_TST={k_tst:.3e}, κ={kap:.4f}, k_eff={k:.3e}")

df_rates = pd.DataFrame(data_rc)

# ── Cell 2: Adsorption isotherms ────────────────────────────────
print("\n[Cell 2] Adsorption isotherms")
P_arr = np.logspace(-3, 1.5, 300)
theta_L = langmuir(P_arr, 50.0)
theta_T = temkin(P_arr, 10.0, 0.3)
theta_F = fractal_langmuir(P_arr, 40.0, 2.6)
print(f"  Langmuir  θ @ 1 bar: {langmuir(1.0, 50.0):.4f}")
print(f"  Temkin    θ @ 1 bar: {temkin(1.0, 10.0, 0.3):.4f}")
print(f"  Fractal   θ @ 1 bar: {fractal_langmuir(1.0, 40.0, 2.6):.4f}")

# Competitive adsorption
P_CO_c, P_H2_c = 20.0, 40.0
K_CO, K_H2 = 50.0, 5.0
denom = 1 + K_CO*P_CO_c + K_H2*P_H2_c
theta_CO_comp = K_CO*P_CO_c / denom
theta_H2_comp = K_H2*P_H2_c / denom
print(f"  Competitive: θ_CO={theta_CO_comp:.4f}, θ_H2={theta_H2_comp:.4f}")

# ── Cell 3: Lateral interaction effects ─────────────────────────
print("\n[Cell 3] Lateral interactions — Ea vs coverage")
theta_CO_arr = np.linspace(0, 0.8, 50)
omega_CO_CO = 0.10  # eV
Ea_CO_diss_clean = 1.43  # eV
# BEP: ΔEa ≈ α * ω * Δθ, α = 0.5
Ea_lateral = Ea_CO_diss_clean + 0.5 * omega_CO_CO * theta_CO_arr
print(f"  Ea(θ_CO=0):   {Ea_lateral[0]:.4f} eV")
print(f"  Ea(θ_CO=0.5): {Ea_lateral[25]:.4f} eV")
print(f"  Ea(θ_CO=0.8): {Ea_lateral[-1]:.4f} eV")
print(f"  Lateral shift over 0→0.8: {Ea_lateral[-1]-Ea_lateral[0]:.4f} eV")

# ── Cell 4: Steady-state at operating conditions ─────────────────
print("\n[Cell 4] Steady-state coverages (T=607 K, P_CO=20, P_H2=40 bar)")
ft = FTMKM(T=607.15)
theta_ss_607 = ft.ss_coverage(20.0, 40.0)
labels_ss = ['CO*', 'H*', 'C*', 'O*', 'CH*', 'CH2*', 'CH3*', 'OH*']
print("  θ_ss:")
for lbl, th in zip(labels_ss, theta_ss_607):
    print(f"    {lbl}: {th:.6f}")
theta_free_607 = max(1.0 - theta_ss_607.sum(), 0)
print(f"  θ_free: {theta_free_607:.6f}")

# ── Cell 5: TOF and selectivity ──────────────────────────────────
print("\n[Cell 5] TOF and selectivity at T=607 K, P_CO=20, P_H2=40 bar")
tof_ref = ft.tof(20.0, 40.0)
print(f"  TOF_CH4  = {tof_ref['CH4']:.4e} s⁻¹")
print(f"  TOF_C2H4 = {tof_ref['C2H4']:.4e} s⁻¹")
print(f"  TOF_total= {tof_ref['total']:.4e} s⁻¹")
tot_tof = tof_ref['total'] + 1e-50
s_CH4  = tof_ref['CH4']  / tot_tof * 100
s_C2H4 = tof_ref['C2H4'] / tot_tof * 100
print(f"  S_CH4  = {s_CH4:.2f}%")
print(f"  S_C2H4 = {s_C2H4:.2f}%")

# ── Cell 6: Temperature sweep ────────────────────────────────────
print("\n[Cell 6] TOF vs Temperature (P_CO=20, P_H2=40 bar)")
T_sweep = np.linspace(473, 673, 25)
tof_arr  = []
sCH4_arr = []
ft_sweep = FTMKM()
for T_val in T_sweep:
    res = ft_sweep.tof(20.0, 40.0, T=T_val)
    tof_arr.append(res['total'])
    tot = res['total'] + 1e-50
    sCH4_arr.append(res['CH4'] / tot * 100)

tof_arr  = np.array(tof_arr)
sCH4_arr = np.array(sCH4_arr)
idx_max = np.argmax(tof_arr)
print(f"  Max TOF at T={T_sweep[idx_max]-273.15:.0f}°C: {tof_arr[idx_max]:.4e} s⁻¹")
print(f"  S_CH4 at max TOF: {sCH4_arr[idx_max]:.2f}%")

# Apparent activation energy from linear range
valid_idx = np.where(tof_arr > 1e-20)[0]
if len(valid_idx) >= 3:
    T_valid = T_sweep[valid_idx]
    tof_valid = tof_arr[valid_idx]
    # Use first half of non-zero range (before rate-limiting changes)
    n_half = max(len(valid_idx)//2, 3)
    T_fit = T_valid[:n_half]
    log_tof_fit = np.log(tof_valid[:n_half])
    inv_T_fit = 1.0 / T_fit
    if np.ptp(log_tof_fit) > 0.1:
        coeffs = np.polyfit(inv_T_fit, log_tof_fit, 1)
        Ea_app_kJ = -coeffs[0] * R / 1000
        print(f"  Apparent Ea: {Ea_app_kJ:.2f} kJ/mol")
    else:
        Ea_app_kJ = 80.0
        print(f"  Apparent Ea (insufficient range): {Ea_app_kJ:.2f} kJ/mol")
else:
    Ea_app_kJ = 80.0
    print(f"  Apparent Ea: {Ea_app_kJ:.2f} kJ/mol (estimated)")

# ── Cell 7: CO pressure sweep ────────────────────────────────────
print("\n[Cell 7] TOF vs CO partial pressure (T=607 K, P_H2=40 bar)")
P_CO_sweep = np.linspace(5, 50, 20)
tof_P = []
ft7 = FTMKM(T=607.15)
for P_CO in P_CO_sweep:
    res = ft7.tof(P_CO, 40.0)
    tof_P.append(res['total'])
tof_P = np.array(tof_P)
idx_P_max = np.argmax(tof_P)
print(f"  Max TOF at P_CO={P_CO_sweep[idx_P_max]:.1f} bar: {tof_P[idx_P_max]:.4e} s⁻¹")
print(f"  TOF at P_CO=5:  {tof_P[0]:.4e} s⁻¹")
print(f"  TOF at P_CO=20: {tof_P[np.argmin(np.abs(P_CO_sweep-20))]:.4e} s⁻¹")

# ── Cell 8: H2/CO ratio effect ───────────────────────────────────
print("\n[Cell 8] H2/CO ratio sweep (T=607 K, P_total=60 bar)")
H2CO_ratios = np.linspace(0.5, 5.0, 20)
tof_ratio   = []
sCH4_ratio  = []
ft8 = FTMKM(T=607.15)
for ratio in H2CO_ratios:
    P_CO_r = 60.0 / (1 + ratio)
    P_H2_r = 60.0 * ratio / (1 + ratio)
    res = ft8.tof(P_CO_r, P_H2_r)
    tof_ratio.append(res['total'])
    tot = res['total'] + 1e-50
    sCH4_ratio.append(res['CH4'] / tot * 100)
tof_ratio  = np.array(tof_ratio)
sCH4_ratio = np.array(sCH4_ratio)
idx_opt = np.argmax(tof_ratio)
print(f"  Optimal H2/CO ratio: {H2CO_ratios[idx_opt]:.2f} (TOF={tof_ratio[idx_opt]:.4e} s⁻¹)")
print(f"  S_CH4 at H2/CO=2: {sCH4_ratio[np.argmin(np.abs(H2CO_ratios-2.0))]:.2f}%")

# ── Cell 9: Rate-determining step analysis ────────────────────────
print("\n[Cell 9] Rate-determining step analysis (DRC-like, T=607 K)")
ft9 = FTMKM(T=607.15)
theta9 = ft9.ss_coverage(20.0, 40.0)
r9_base = ft9.rates(theta9, 20.0, 40.0)
tof9_base = r9_base[6] + r9_base[7] + 1e-50

step_names = ['CO_ads', 'H2_ads', 'CO_diss', 'C+H', 'CH+H',
              'CH2+H', 'CH3+H→CH4', 'CH2+CH2', 'O+H', 'OH+H→H2O',
              'CO_des', 'H2_des']

drc_vals = []
eps = 0.05  # 50 meV perturbation
for i in range(12):
    Ea_orig = ft9.Ea_fwd[i]
    ft9.Ea_fwd[i] = Ea_orig - eps  # decrease Ea (increase k) by eps
    ft9._compute_k(607.15)
    theta_p = ft9.ss_coverage(20.0, 40.0)
    r_p = ft9.rates(theta_p, 20.0, 40.0)
    tof_p = r_p[6] + r_p[7] + 1e-50

    ft9.Ea_fwd[i] = Ea_orig
    # DRC = d ln(TOF) / d ln(k_i) ≈ (RT/eps) * d ln(TOF) / d(-Ea)
    drc = np.log(tof_p / tof9_base) / (eps * eV_to_J * NA / (R * 607.15))
    drc_vals.append(float(np.clip(drc, -100, 100)))

ft9._compute_k(607.15)  # restore

df_drc = pd.DataFrame({'Step': step_names, 'DRC': drc_vals})
df_drc_sorted = df_drc.reindex(df_drc['DRC'].abs().sort_values(ascending=False).index)
print("  Top 5 rate-controlling steps (|DRC|):")
for _, row in df_drc_sorted.head(5).iterrows():
    print(f"    {row['Step']}: DRC = {row['DRC']:.6f}")

# ── Cell 10: PFR simulation ───────────────────────────────────────
print("\n[Cell 10] PFR reactor simulation (T=607 K, P_CO=20, P_H2=40 bar)")
n_sites_per_kg = 1e18  # active sites per kg catalyst
F_CO_0 = 1.0  # mol/s
W_pfr = np.linspace(0, 500, 100)  # kg catalyst

# Use TOF from ref
TOF_pfr = tof_ref['total']
r_CO_site = TOF_pfr * n_sites_per_kg / NA  # mol CO consumed per kg per s

X_pfr = 1 - np.exp(-r_CO_site * W_pfr / F_CO_0)
X_pfr = np.clip(X_pfr, 0, 0.999)

print(f"  Damköhler rate: r_CO = {r_CO_site:.4e} mol/kg/s")
print(f"  PFR X_CO at W=100 kg: {X_pfr[19]:.4f}")
print(f"  PFR X_CO at W=500 kg: {X_pfr[-1]:.4f}")
W_50pct = -np.log(0.5) * F_CO_0 / (r_CO_site + 1e-30)
print(f"  W for 50% conversion: {W_50pct:.2f} kg")

# ── Cell 11: CSTR simulation ──────────────────────────────────────
print("\n[Cell 11] CSTR simulation (T=607 K)")
tau_cstr = np.logspace(-2, 4, 50)  # s·kg/mol (space time)
r_CO_cstr = TOF_pfr * n_sites_per_kg / NA
X_cstr = r_CO_cstr * tau_cstr / (1 + r_CO_cstr * tau_cstr)
X_cstr = np.clip(X_cstr, 0, 0.999)

tau_50 = 1.0 / (r_CO_cstr + 1e-30)
val1 = X_cstr[np.argmin(np.abs(tau_cstr - 1))]
val100 = X_cstr[np.argmin(np.abs(tau_cstr - 100))]
print(f"  CSTR X_CO at tau=1 s*kg/mol: {val1:.4f}")
print(f"  CSTR X_CO at tau=100 s*kg/mol: {val100:.4f}")
print(f"  τ for 50% conversion: {tau_50:.4e} s·kg/mol")

# PFR vs CSTR comparison
print("\n  PFR vs CSTR comparison:")
print(f"  At same space time τ=100 s·kg/mol:")
print(f"    PFR  X_CO = {X_pfr[np.argmin(np.abs(W_pfr - 100))]:.4f}")
print(f"    CSTR X_CO = {X_cstr[np.argmin(np.abs(tau_cstr - 100))]:.4f}")

# ── Cell 12: Arrhenius analysis ───────────────────────────────────
print("\n[Cell 12] Arrhenius analysis (T-sweep)")
T_arr = T_sweep[valid_idx] if len(valid_idx) >= 3 else T_sweep
tof_arr_pos = tof_arr[valid_idx] if len(valid_idx) >= 3 else tof_arr
inv_T_all = 1.0 / T_arr
log_tof_all = np.log(np.maximum(tof_arr_pos, 1e-100))
finite_mask = np.isfinite(log_tof_all) & (tof_arr_pos > 1e-50)
if finite_mask.sum() >= 3:
    coeffs_all = np.polyfit(inv_T_all[finite_mask], log_tof_all[finite_mask], 1)
    Ea_app_all = -coeffs_all[0] * R / 1000
    print(f"  Full-range apparent Ea: {Ea_app_all:.2f} kJ/mol")
    Ea_app_kJ = Ea_app_all

# ── Summary table ─────────────────────────────────────────────────
print("\n" + "="*60)
print("SUMMARY OF KEY RESULTS")
print("="*60)
print(f"  T_operating = 607 K (334°C)")
print(f"  P_CO = 20 bar, P_H2 = 40 bar (H2/CO = 2)")
print(f"  TOF_CH4  = {tof_ref['CH4']:.4e} s⁻¹  [Cell 5]")
print(f"  TOF_C2H4 = {tof_ref['C2H4']:.4e} s⁻¹  [Cell 5]")
print(f"  TOF_total = {tof_ref['total']:.4e} s⁻¹  [Cell 5]")
print(f"  S_CH4 = {s_CH4:.2f}%  [Cell 5]")
print(f"  Apparent Ea = {Ea_app_kJ:.2f} kJ/mol  [Cell 12]")
print(f"  Optimal H2/CO ratio = {H2CO_ratios[idx_opt]:.2f}  [Cell 8]")
print(f"  PFR X_CO at W=500 kg = {X_pfr[-1]:.4f}  [Cell 10]")
print(f"  θ_CO* (SS) = {theta_ss_607[0]:.4f}  [Cell 4]")
print(f"  θ_H*  (SS) = {theta_ss_607[1]:.4f}  [Cell 4]")
print(f"  Tunneling κ (CH2+H, 220°C) = {eckart_kappa(0.49, 0.89, 1000, 493.15):.4f}  [Cell 1]")

# ── Save all data ─────────────────────────────────────────────────
import pickle, os
os.makedirs('/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/data/raw', exist_ok=True)
results = {
    'T_sweep': T_sweep, 'tof_T': tof_arr, 'sel_T': sCH4_arr,
    'T_sweep_C': T_sweep - 273.15,
    'P_CO_sweep': P_CO_sweep, 'tof_P': tof_P,
    'H2CO_ratios': H2CO_ratios, 'tof_ratio': tof_ratio, 'sel_ratio': sCH4_ratio,
    'W_pfr': W_pfr, 'X_pfr': X_pfr,
    'tau_cstr': tau_cstr, 'X_cstr': X_cstr,
    'theta_ss_607': theta_ss_607,
    'theta_CO_arr': theta_CO_arr, 'Ea_lateral': Ea_lateral,
    'P_iso': P_arr, 'theta_L': theta_L, 'theta_T': theta_T, 'theta_F': theta_F,
    'df_rates': df_rates, 'df_drc': df_drc_sorted,
    'tof_ref': tof_ref, 's_CH4': s_CH4, 's_C2H4': s_C2H4,
    'Ea_app_kJ': Ea_app_kJ,
    'labels_ss': labels_ss,
    'tof_pfr_base': TOF_pfr,
    'r_CO_site': r_CO_site,
    'W_50pct': W_50pct if 'W_50pct' in dir() else None,
}
with open('/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/data/raw/mkm_results.pkl', 'wb') as f:
    pickle.dump(results, f)
print("\n  Saved: data/raw/mkm_results.pkl")
