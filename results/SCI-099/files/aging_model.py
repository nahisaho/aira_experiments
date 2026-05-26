#!/usr/bin/env python3
"""
Integrated Mathematical Model of Aging: ODE-based simulation framework.

This module implements a multi-hallmark aging model integrating:
- Hallmarks of aging interaction network (telomere, epigenetic, mitochondrial, senescence, etc.)
- Reliability theory (damage accumulation) with antagonistic pleiotropy
- Senolytic therapy modeling
- Caloric restriction / Rapamycin / NAD+ precursor interventions
- Interspecies lifespan scaling
- Combination intervention optimization
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json

sns.set_theme(style="whitegrid", font_scale=1.1)
FIGURE_DIR = "figures"
os.makedirs(FIGURE_DIR, exist_ok=True)

# =============================================================================
# 1. INTEGRATED HALLMARKS OF AGING ODE MODEL
# =============================================================================

class AgingHallmarksModel:
    """
    ODE system modeling 8 key hallmarks of aging and their interactions.
    
    State variables:
        T: Telomere integrity (1 = full, 0 = depleted)
        E: Epigenetic integrity (1 = youthful, 0 = fully altered)
        M: Mitochondrial function (1 = full, 0 = dysfunctional)
        S: Senescent cell fraction (0 = none, 1 = fully senescent)
        P: Proteostasis capacity (1 = full, 0 = collapsed)
        N: Nutrient sensing regulation (1 = optimal, 0 = deregulated)
        I: Inflammation level (0 = none, 1 = maximal chronic inflammation)
        SC: Stem cell function (1 = full, 0 = exhausted)
    """
    
    def __init__(self, params=None):
        self.default_params = {
            # Intrinsic decay rates
            'alpha_T': 0.008,   # Telomere shortening rate
            'alpha_E': 0.006,   # Epigenetic drift rate
            'alpha_M': 0.007,   # Mitochondrial decline rate
            'alpha_S': 0.010,   # Senescence accumulation rate
            'alpha_P': 0.005,   # Proteostasis decline rate
            'alpha_N': 0.004,   # Nutrient sensing deregulation rate
            'alpha_I': 0.009,   # Inflammation increase rate
            'alpha_SC': 0.006,  # Stem cell exhaustion rate
            
            # Cross-hallmark interaction strengths
            'beta_TM': 0.003,   # Telomere → Mitochondria coupling
            'beta_TS': 0.004,   # Telomere → Senescence coupling
            'beta_MS': 0.005,   # Mitochondria → Senescence coupling
            'beta_MI': 0.004,   # Mitochondria → Inflammation coupling
            'beta_SI': 0.006,   # Senescence → Inflammation (SASP)
            'beta_IE': 0.003,   # Inflammation → Epigenetic coupling
            'beta_IP': 0.003,   # Inflammation → Proteostasis coupling
            'beta_EP': 0.002,   # Epigenetic → Proteostasis coupling
            'beta_NSC': 0.003,  # Nutrient sensing → Stem cell coupling
            'beta_ISC': 0.004,  # Inflammation → Stem cell coupling
            'beta_SSC': 0.003,  # Senescence → Stem cell coupling
            'beta_NM': 0.002,   # Nutrient sensing → Mitochondria coupling
            
            # Senolytic parameters
            'gamma_sen': 0.0,   # Senolytic clearance rate
            
            # Intervention parameters
            'cr_factor': 1.0,   # Caloric restriction factor (1 = no CR)
            'rapa_factor': 1.0, # Rapamycin factor (1 = no drug)
            'nad_factor': 1.0,  # NAD+ precursor factor (1 = no supplement)
        }
        if params:
            self.default_params.update(params)
        self.params = self.default_params
    
    def ode_system(self, t, y):
        T, E, M, S, P, N, I, SC = y
        p = self.params
        
        # Intervention modulation
        cr = p['cr_factor']
        rapa = p['rapa_factor']
        nad = p['nad_factor']
        
        # Telomere integrity decline (accelerated by inflammation)
        dT = -p['alpha_T'] * T - p['beta_TM'] * (1 - M) * T
        
        # Epigenetic integrity decline (accelerated by inflammation)
        dE = -p['alpha_E'] * E - p['beta_IE'] * I * E
        
        # Mitochondrial function decline
        # NAD+ precursors slow mitochondrial decline; CR helps via nutrient sensing
        mito_protection = (1.0 / nad) * (1.0 / cr)
        dM = -p['alpha_M'] * M * mito_protection - p['beta_NM'] * (1 - N) * M
        
        # Senescent cell accumulation
        # Driven by telomere shortening, mitochondrial dysfunction
        # Senolytics clear senescent cells
        dS = (p['alpha_S'] * (1 - S) * ((1 - T) + (1 - M)) / 2.0
              + p['beta_TS'] * (1 - T) * (1 - S)
              + p['beta_MS'] * (1 - M) * (1 - S)
              - p['gamma_sen'] * S)
        
        # Proteostasis decline (accelerated by inflammation, epigenetic drift)
        dP = -p['alpha_P'] * P - p['beta_IP'] * I * P - p['beta_EP'] * (1 - E) * P
        
        # Nutrient sensing deregulation (Rapamycin/CR improve mTOR regulation)
        nutrient_protection = (1.0 / rapa) * (1.0 / cr)
        dN = -p['alpha_N'] * N * nutrient_protection
        
        # Inflammation (SASP from senescence, mitochondrial ROS)
        dI = (p['alpha_I'] * (1 - I) * 0.3
              + p['beta_SI'] * S * (1 - I)
              + p['beta_MI'] * (1 - M) * (1 - I)
              - 0.02 * I * (1.0 / cr))  # CR reduces inflammation
        
        # Stem cell exhaustion
        dSC = (-p['alpha_SC'] * SC
               - p['beta_NSC'] * (1 - N) * SC
               - p['beta_ISC'] * I * SC
               - p['beta_SSC'] * S * SC)
        
        return [dT, dE, dM, dS, dP, dN, dI, dSC]
    
    def simulate(self, t_span=(0, 100), y0=None, t_eval=None):
        if y0 is None:
            y0 = [0.95, 0.95, 0.95, 0.02, 0.95, 0.95, 0.05, 0.95]
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 1000)
        
        sol = solve_ivp(self.ode_system, t_span, y0, t_eval=t_eval,
                        method='RK45', max_step=0.5,
                        rtol=1e-8, atol=1e-10)
        return sol
    
    def compute_health_index(self, sol):
        """Weighted composite health index from all hallmarks."""
        weights = np.array([0.12, 0.10, 0.15, -0.18, 0.13, 0.10, -0.12, 0.10])
        # Normalize senescence and inflammation (higher = worse)
        y = sol.y.copy()
        health = np.zeros(len(sol.t))
        for i, w in enumerate(weights):
            health += w * y[i]
        # Normalize to [0, 1]
        health = (health - health.min()) / (health.max() - health.min() + 1e-12)
        return health
    
    def compute_mortality_rate(self, health):
        """Gompertz-like mortality rate from health index."""
        mu_0 = 0.0001
        gamma = 8.0
        return mu_0 * np.exp(gamma * (1.0 - health))
    
    def compute_survival(self, t, mortality_rate):
        """Survival curve from mortality rate."""
        from scipy.integrate import cumulative_trapezoid
        cum_hazard = cumulative_trapezoid(mortality_rate, t, initial=0)
        return np.exp(-cum_hazard)
    
    def compute_median_lifespan(self, t, survival):
        """Find median lifespan from survival curve."""
        idx = np.searchsorted(-survival, -0.5)
        if idx < len(t):
            return t[idx]
        return t[-1]


# =============================================================================
# 2. RELIABILITY THEORY + ANTAGONISTIC PLEIOTROPY
# =============================================================================

class ReliabilityAgingModel:
    """
    Reliability theory model with antagonistic pleiotropy integration.
    
    Models organism as a system of n redundant components, each with
    failure rate lambda. AP genes provide early-life benefit b but
    increase late-life failure rate by factor delta.
    """
    
    def __init__(self, n_components=500, base_failure_rate=0.001,
                 n_ap_genes=20, ap_early_benefit=0.3, ap_late_penalty=0.5,
                 ap_onset_age=40):
        self.n = n_components
        self.lam = base_failure_rate
        self.n_ap = n_ap_genes
        self.b = ap_early_benefit
        self.delta = ap_late_penalty
        self.ap_onset = ap_onset_age
    
    def failure_rate(self, t):
        """Time-dependent component failure rate with AP effects."""
        base = self.lam
        if t > self.ap_onset:
            ap_effect = self.delta * (1 - np.exp(-(t - self.ap_onset) / 20))
            base *= (1 + ap_effect * self.n_ap / self.n)
        return base
    
    def system_reliability(self, t_array):
        """System reliability R(t) for redundant parallel system."""
        R = np.zeros(len(t_array))
        for i, t in enumerate(t_array):
            lam_t = self.failure_rate(t)
            # Component reliability
            r_comp = np.exp(-lam_t * t)
            # System reliability (1 - (1-r)^n for parallel system)
            R[i] = 1.0 - (1.0 - r_comp) ** self.n
        return R
    
    def mortality_rate_from_reliability(self, t_array, R):
        """Compute mortality rate mu(t) = -R'(t)/R(t) using Gompertz approximation."""
        dt = np.diff(t_array)
        dR = np.diff(R)
        mu = -dR / (R[:-1] * dt + 1e-15)
        mu = np.clip(mu, 1e-8, 10)
        # Smooth with moving average
        kernel = 5
        if len(mu) > kernel:
            mu_smooth = np.convolve(mu, np.ones(kernel)/kernel, mode='same')
            mu = mu_smooth
        return mu
    
    def early_life_fitness(self, t):
        """Fitness benefit from AP genes in early life."""
        return 1.0 + self.b * self.n_ap / self.n * np.exp(-t / 30)


# =============================================================================
# 3. SENOLYTIC THERAPY MODEL
# =============================================================================

class SenolyticModel:
    """
    ODE model for senolytic therapy effects.
    
    dN/dt = r*N*(1 - (N+S)/K) - alpha*N
    dS/dt = alpha*N - beta*S - gamma(t)*S
    dD/dt = delta*S  (cumulative SASP-driven damage)
    
    gamma(t) = senolytic clearance rate (pulsed or continuous)
    """
    
    def __init__(self, r=0.05, K=1.0, alpha=0.008, beta=0.002,
                 delta=0.01, gamma_base=0.0, dose_schedule=None):
        self.r = r
        self.K = K
        self.alpha = alpha
        self.beta = beta
        self.delta = delta
        self.gamma_base = gamma_base
        self.dose_schedule = dose_schedule or []  # [(start, end, dose), ...]
    
    def gamma(self, t):
        """Time-dependent senolytic clearance rate."""
        g = self.gamma_base
        for start, end, dose in self.dose_schedule:
            if start <= t <= end:
                g += dose
        return g
    
    def ode(self, t, y):
        N, S, D = y
        g = self.gamma(t)
        dN = self.r * N * (1 - (N + S) / self.K) - self.alpha * N
        dS = self.alpha * N - self.beta * S - g * S
        dD = self.delta * S
        return [dN, dS, dD]
    
    def simulate(self, t_span=(0, 100), y0=None, t_eval=None):
        if y0 is None:
            y0 = [0.9, 0.05, 0.0]
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 2000)
        sol = solve_ivp(self.ode, t_span, y0, t_eval=t_eval,
                        method='RK45', max_step=0.2, rtol=1e-8, atol=1e-10)
        return sol


# =============================================================================
# 4. INTERVENTION PATHWAY MODEL (mTOR/AMPK/SIRT1/NAD+)
# =============================================================================

class InterventionPathwayModel:
    """
    ODE model of mTOR/AMPK/SIRT1/NAD+ signaling with interventions.
    
    State variables:
        mTOR: mTOR activity (0-1)
        AMPK: AMPK activity (0-1)
        SIRT1: SIRT1 activity (0-1)
        NAD: NAD+ level (0-1)
        autophagy: Autophagy level (0-1)
        damage: Accumulated damage (0+)
    """
    
    def __init__(self, cr_level=0.0, rapa_dose=0.0, nad_supplement=0.0):
        self.cr = cr_level         # 0-1, 0=no CR, 1=40% CR
        self.rapa = rapa_dose      # 0-1, rapamycin dose
        self.nad_sup = nad_supplement  # 0-1, NAD+ precursor dose
    
    def ode(self, t, y):
        mTOR, AMPK, SIRT1, NAD, autophagy, damage = y
        
        # mTOR dynamics: activated by nutrients, inhibited by AMPK and rapamycin
        nutrient_signal = 1.0 - self.cr * 0.4
        dmTOR = 0.5 * (nutrient_signal - mTOR) - 0.3 * AMPK * mTOR - 0.8 * self.rapa * mTOR
        
        # AMPK: activated by energy stress (CR), inhibited by mTOR
        dAMPK = 0.4 * (self.cr * 0.5 + 0.3 - AMPK) - 0.2 * mTOR * AMPK + 0.1 * (1 - AMPK)
        
        # SIRT1: dependent on NAD+, activated by CR
        dSIRT1 = 0.3 * NAD * (1 - SIRT1) + 0.2 * self.cr * (1 - SIRT1) - 0.1 * SIRT1
        
        # NAD+: declines with age/damage, boosted by supplement
        dNAD = -0.005 * damage * NAD + 0.3 * self.nad_sup * (1 - NAD) - 0.02 * NAD
        
        # Autophagy: inhibited by mTOR, activated by AMPK and SIRT1
        dautophagy = 0.3 * AMPK * (1 - autophagy) + 0.2 * SIRT1 * (1 - autophagy) \
                     - 0.5 * mTOR * autophagy - 0.05 * autophagy
        
        # Damage accumulation: reduced by autophagy and SIRT1-mediated repair
        ddamage = 0.02 * (1 + 0.5 * mTOR) - 0.03 * autophagy - 0.02 * SIRT1
        ddamage = max(ddamage, -0.01 * damage)  # Can't have negative damage rate below repair
        
        return [dmTOR, dAMPK, dSIRT1, dNAD, dautophagy, ddamage]
    
    def simulate(self, t_span=(0, 80), y0=None, t_eval=None):
        if y0 is None:
            y0 = [0.7, 0.3, 0.5, 0.8, 0.3, 0.0]
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 1000)
        sol = solve_ivp(self.ode, t_span, y0, t_eval=t_eval,
                        method='RK45', max_step=0.2, rtol=1e-8, atol=1e-10)
        return sol


# =============================================================================
# 5. INTERSPECIES LIFESPAN SCALING MODEL
# =============================================================================

class InterspeciesModel:
    """
    Allometric scaling model for lifespan across species.
    
    L = a * M^b * R^c * D^d
    
    where M = body mass, R = metabolic rate scaling, D = DNA repair capacity
    """
    
    # Species data: (name, body_mass_kg, max_lifespan_years, relative_metabolic_rate, relative_dna_repair)
    SPECIES_DATA = [
        ("Mouse", 0.03, 4, 7.0, 0.3),
        ("Rat", 0.3, 5, 5.5, 0.35),
        ("Rabbit", 2.0, 12, 3.5, 0.4),
        ("Dog", 15.0, 20, 2.5, 0.5),
        ("Cat", 4.0, 25, 3.0, 0.45),
        ("Naked Mole Rat", 0.035, 32, 4.0, 0.9),
        ("Bat (Myotis)", 0.008, 40, 8.0, 0.85),
        ("Human", 70.0, 122, 1.0, 1.0),
        ("Elephant", 5000.0, 70, 0.4, 0.8),
        ("Bowhead Whale", 75000.0, 211, 0.15, 0.95),
        ("Greenland Shark", 400.0, 400, 0.05, 0.7),
        ("Tortoise", 250.0, 190, 0.1, 0.65),
    ]
    
    @staticmethod
    def fit_allometric_model():
        """Fit allometric scaling parameters using least squares."""
        data = InterspeciesModel.SPECIES_DATA
        names = [d[0] for d in data]
        M = np.array([d[1] for d in data])
        L = np.array([d[2] for d in data])
        R = np.array([d[3] for d in data])
        D = np.array([d[4] for d in data])
        
        # Log-linear fit: log(L) = log(a) + b*log(M) + c*log(R) + d*log(D)
        X = np.column_stack([np.ones(len(M)), np.log(M), np.log(R), np.log(D)])
        y = np.log(L)
        coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
        
        a = np.exp(coeffs[0])
        b, c, d_coeff = coeffs[1], coeffs[2], coeffs[3]
        
        L_pred = a * M**b * R**c * D**d_coeff
        r_squared = 1 - np.sum((L - L_pred)**2) / np.sum((L - np.mean(L))**2)
        
        return {
            'a': a, 'b': b, 'c': c, 'd': d_coeff,
            'names': names, 'M': M, 'L_actual': L, 'L_predicted': L_pred,
            'R_squared': r_squared
        }


# =============================================================================
# 6. COMBINATION INTERVENTION OPTIMIZATION
# =============================================================================

class InterventionOptimizer:
    """
    Optimize combination of anti-aging interventions using differential evolution.
    
    Decision variables: [cr_level, rapa_dose, nad_dose, sen_frequency, sen_dose]
    Objective: maximize healthspan (time until health index < 0.5)
    """
    
    def __init__(self):
        self.bounds = [
            (0.0, 0.6),   # CR level (0-60%)
            (0.0, 1.0),   # Rapamycin dose
            (0.0, 1.0),   # NAD+ supplement dose
            (0.0, 0.1),   # Senolytic clearance rate
        ]
    
    def objective(self, x):
        cr, rapa, nad, sen = x
        
        model = AgingHallmarksModel(params={
            'cr_factor': 1.0 + cr,
            'rapa_factor': 1.0 + rapa * 0.5,
            'nad_factor': 1.0 + nad * 0.3,
            'gamma_sen': sen,
        })
        
        sol = model.simulate(t_span=(0, 120))
        health = model.compute_health_index(sol)
        
        # Healthspan: time until health < 0.5
        below_threshold = np.where(health < 0.5)[0]
        if len(below_threshold) > 0:
            healthspan = sol.t[below_threshold[0]]
        else:
            healthspan = sol.t[-1]
        
        # Penalty for side effects (high doses)
        side_effect_penalty = 0.1 * (cr**2 + rapa**2 + nad**2 + (sen * 10)**2)
        
        # Maximize healthspan minus side effects (negative for minimization)
        return -(healthspan - side_effect_penalty)
    
    def optimize(self, maxiter=50):
        result = differential_evolution(
            self.objective, self.bounds,
            maxiter=maxiter, seed=42, tol=1e-4,
            popsize=15, mutation=(0.5, 1.5), recombination=0.7
        )
        return {
            'cr_level': result.x[0],
            'rapa_dose': result.x[1],
            'nad_dose': result.x[2],
            'sen_rate': result.x[3],
            'healthspan': -result.fun,
            'success': result.success,
        }


# =============================================================================
# SIMULATION AND VISUALIZATION
# =============================================================================

def run_all_simulations():
    results = {}
    
    # ---- Simulation 1: Baseline aging hallmarks ----
    print("Running Simulation 1: Baseline Aging Hallmarks...")
    model_base = AgingHallmarksModel()
    sol_base = model_base.simulate()
    health_base = model_base.compute_health_index(sol_base)
    mortality_base = model_base.compute_mortality_rate(health_base)
    survival_base = model_base.compute_survival(sol_base.t, mortality_base)
    
    labels = ['Telomere', 'Epigenetic', 'Mitochondrial', 'Senescence',
              'Proteostasis', 'Nutrient Sensing', 'Inflammation', 'Stem Cell']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel A: All hallmarks over time
    ax = axes[0, 0]
    colors = sns.color_palette("husl", 8)
    for i in range(8):
        ax.plot(sol_base.t, sol_base.y[i], label=labels[i], color=colors[i], linewidth=1.5)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('State Value')
    ax.set_title('A) Hallmarks of Aging Dynamics')
    ax.legend(fontsize=7, ncol=2)
    ax.set_xlim(0, 100)
    
    # Panel B: Health index
    ax = axes[0, 1]
    ax.plot(sol_base.t, health_base, 'k-', linewidth=2)
    ax.axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='Healthspan threshold')
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Health Index')
    ax.set_title('B) Composite Health Index')
    ax.legend()
    
    # Panel C: Mortality rate
    ax = axes[1, 0]
    ax.semilogy(sol_base.t, mortality_base, 'b-', linewidth=1.5)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Mortality Rate (log scale)')
    ax.set_title('C) Gompertz Mortality Rate')
    
    # Panel D: Survival curve
    ax = axes[1, 1]
    ax.plot(sol_base.t, survival_base, 'g-', linewidth=2)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Survival Probability')
    ax.set_title('D) Survival Curve')
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig1_hallmarks_baseline.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig1_hallmarks_baseline.png")
    
    # ---- Simulation 2: Hallmark interaction network heatmap ----
    print("Running Simulation 2: Interaction Network...")
    interaction_matrix = np.array([
        [0,    0,    0.003, 0.004, 0,     0,     0,     0],
        [0,    0,    0,     0,     0.002, 0,     0,     0],
        [0,    0,    0,     0.005, 0,     0,     0.004, 0],
        [0,    0,    0,     0,     0,     0,     0.006, 0.003],
        [0,    0,    0,     0,     0,     0,     0,     0],
        [0,    0,    0.002, 0,     0,     0,     0,     0.003],
        [0,    0.003, 0,    0,     0.003, 0,     0,     0.004],
        [0,    0,    0,     0,     0,     0,     0,     0],
    ])
    
    fig, ax = plt.subplots(figsize=(9, 7))
    mask = interaction_matrix == 0
    sns.heatmap(interaction_matrix, annot=True, fmt='.3f', cmap='YlOrRd',
                xticklabels=labels, yticklabels=labels, mask=mask,
                linewidths=0.5, ax=ax, cbar_kws={'label': 'Coupling Strength'})
    ax.set_title('Hallmark Interaction Network (Coupling Strengths)')
    ax.set_xlabel('Target Hallmark')
    ax.set_ylabel('Source Hallmark')
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig2_interaction_network.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig2_interaction_network.png")
    
    # ---- Simulation 3: Reliability Theory + Antagonistic Pleiotropy ----
    print("Running Simulation 3: Reliability Theory...")
    t_array = np.linspace(0, 100, 2000)
    
    # Compare with and without AP
    rel_no_ap = ReliabilityAgingModel(n_ap_genes=0)
    rel_with_ap = ReliabilityAgingModel(n_ap_genes=20)
    rel_high_ap = ReliabilityAgingModel(n_ap_genes=50)
    
    R_no_ap = rel_no_ap.system_reliability(t_array)
    R_with_ap = rel_with_ap.system_reliability(t_array)
    R_high_ap = rel_high_ap.system_reliability(t_array)
    
    mu_no_ap = rel_no_ap.mortality_rate_from_reliability(t_array, R_no_ap)
    mu_with_ap = rel_with_ap.mortality_rate_from_reliability(t_array, R_with_ap)
    mu_high_ap = rel_high_ap.mortality_rate_from_reliability(t_array, R_high_ap)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    ax = axes[0]
    ax.plot(t_array, R_no_ap, label='No AP genes', linewidth=2)
    ax.plot(t_array, R_with_ap, label='20 AP genes', linewidth=2)
    ax.plot(t_array, R_high_ap, label='50 AP genes', linewidth=2)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('System Reliability')
    ax.set_title('A) Reliability Curves')
    ax.legend()
    
    ax = axes[1]
    ax.semilogy(t_array[:-1], mu_no_ap, label='No AP', linewidth=1.5)
    ax.semilogy(t_array[:-1], mu_with_ap, label='20 AP', linewidth=1.5)
    ax.semilogy(t_array[:-1], mu_high_ap, label='50 AP', linewidth=1.5)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Mortality Rate (log)')
    ax.set_title('B) Mortality Rate (Gompertz-like)')
    ax.legend()
    
    ax = axes[2]
    fitness_t = np.linspace(0, 100, 500)
    f1 = rel_no_ap.early_life_fitness(fitness_t)
    f2 = rel_with_ap.early_life_fitness(fitness_t)
    f3 = rel_high_ap.early_life_fitness(fitness_t)
    ax.plot(fitness_t, f1, label='No AP', linewidth=2)
    ax.plot(fitness_t, f2, label='20 AP', linewidth=2)
    ax.plot(fitness_t, f3, label='50 AP', linewidth=2)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Fitness')
    ax.set_title('C) Early-Life Fitness Benefit')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig3_reliability_theory.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig3_reliability_theory.png")
    
    # ---- Simulation 4: Senolytic Therapy ----
    print("Running Simulation 4: Senolytic Therapy...")
    
    # No treatment
    sen_none = SenolyticModel()
    sol_sen_none = sen_none.simulate()
    
    # Continuous senolytic
    sen_cont = SenolyticModel(gamma_base=0.05)
    sol_sen_cont = sen_cont.simulate()
    
    # Intermittent senolytic (every 5 years for 1 year)
    schedule = [(i, i+1, 0.15) for i in range(30, 100, 5)]
    sen_pulse = SenolyticModel(dose_schedule=schedule)
    sol_sen_pulse = sen_pulse.simulate()
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    ax = axes[0]
    ax.plot(sol_sen_none.t, sol_sen_none.y[0], label='Normal cells (no tx)', linewidth=1.5)
    ax.plot(sol_sen_cont.t, sol_sen_cont.y[0], '--', label='Normal cells (continuous)', linewidth=1.5)
    ax.plot(sol_sen_pulse.t, sol_sen_pulse.y[0], ':', label='Normal cells (pulsed)', linewidth=1.5)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Normal Cell Fraction')
    ax.set_title('A) Normal Cell Population')
    ax.legend(fontsize=8)
    
    ax = axes[1]
    ax.plot(sol_sen_none.t, sol_sen_none.y[1], label='No treatment', linewidth=2, color='red')
    ax.plot(sol_sen_cont.t, sol_sen_cont.y[1], '--', label='Continuous', linewidth=2, color='blue')
    ax.plot(sol_sen_pulse.t, sol_sen_pulse.y[1], ':', label='Pulsed', linewidth=2, color='green')
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Senescent Cell Fraction')
    ax.set_title('B) Senescent Cell Burden')
    ax.legend()
    
    ax = axes[2]
    ax.plot(sol_sen_none.t, sol_sen_none.y[2], label='No treatment', linewidth=2, color='red')
    ax.plot(sol_sen_cont.t, sol_sen_cont.y[2], '--', label='Continuous', linewidth=2, color='blue')
    ax.plot(sol_sen_pulse.t, sol_sen_pulse.y[2], ':', label='Pulsed', linewidth=2, color='green')
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Cumulative SASP Damage')
    ax.set_title('C) Cumulative Damage')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig4_senolytic_therapy.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig4_senolytic_therapy.png")
    
    # ---- Simulation 5: Intervention Pathways ----
    print("Running Simulation 5: Intervention Pathways...")
    
    scenarios = {
        'No intervention': (0.0, 0.0, 0.0),
        'CR (30%)': (0.75, 0.0, 0.0),
        'Rapamycin': (0.0, 0.7, 0.0),
        'NAD+ precursor': (0.0, 0.0, 0.7),
        'Combined': (0.5, 0.5, 0.5),
    }
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    var_names = ['mTOR', 'AMPK', 'SIRT1', 'NAD+', 'Autophagy', 'Damage']
    
    for name, (cr, rapa, nad) in scenarios.items():
        model = InterventionPathwayModel(cr_level=cr, rapa_dose=rapa, nad_supplement=nad)
        sol = model.simulate()
        for idx, ax in enumerate(axes.flat):
            ax.plot(sol.t, sol.y[idx], label=name, linewidth=1.5)
    
    for idx, ax in enumerate(axes.flat):
        ax.set_xlabel('Time (years)')
        ax.set_ylabel(var_names[idx])
        ax.set_title(var_names[idx])
        if idx == 0:
            ax.legend(fontsize=7, loc='best')
    
    plt.suptitle('Intervention Pathway Dynamics: mTOR/AMPK/SIRT1/NAD+', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig5_intervention_pathways.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig5_intervention_pathways.png")
    
    # ---- Simulation 6: Interspecies Lifespan Scaling ----
    print("Running Simulation 6: Interspecies Scaling...")
    
    allometric = InterspeciesModel.fit_allometric_model()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax = axes[0]
    ax.scatter(allometric['M'], allometric['L_actual'], s=100, c='blue', zorder=5, label='Observed')
    ax.scatter(allometric['M'], allometric['L_predicted'], s=60, marker='^', c='red', zorder=5, label='Predicted')
    for i, name in enumerate(allometric['names']):
        ax.annotate(name, (allometric['M'][i], allometric['L_actual'][i]),
                    fontsize=7, ha='left', va='bottom')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Body Mass (kg)')
    ax.set_ylabel('Maximum Lifespan (years)')
    ax.set_title(f'A) Allometric Scaling (R² = {allometric["R_squared"]:.3f})')
    ax.legend()
    
    ax = axes[1]
    residuals = allometric['L_actual'] - allometric['L_predicted']
    colors_bar = ['green' if r > 0 else 'red' for r in residuals]
    ax.barh(allometric['names'], residuals, color=colors_bar, alpha=0.7)
    ax.set_xlabel('Residual (Actual - Predicted) years')
    ax.set_title('B) Prediction Residuals')
    ax.axvline(x=0, color='k', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig6_interspecies_scaling.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig6_interspecies_scaling.png")
    
    results['allometric'] = {
        'a': float(allometric['a']),
        'b': float(allometric['b']),
        'c': float(allometric['c']),
        'd': float(allometric['d']),
        'R_squared': float(allometric['R_squared'])
    }
    
    # ---- Simulation 7: Combination Optimization ----
    print("Running Simulation 7: Combination Optimization...")
    
    optimizer = InterventionOptimizer()
    opt_result = optimizer.optimize(maxiter=30)
    results['optimization'] = opt_result
    print(f"  Optimal: CR={opt_result['cr_level']:.3f}, Rapa={opt_result['rapa_dose']:.3f}, "
          f"NAD={opt_result['nad_dose']:.3f}, Sen={opt_result['sen_rate']:.4f}")
    print(f"  Predicted healthspan: {opt_result['healthspan']:.1f} years")
    
    # Compare baseline vs optimal
    model_opt = AgingHallmarksModel(params={
        'cr_factor': 1.0 + opt_result['cr_level'],
        'rapa_factor': 1.0 + opt_result['rapa_dose'] * 0.5,
        'nad_factor': 1.0 + opt_result['nad_dose'] * 0.3,
        'gamma_sen': opt_result['sen_rate'],
    })
    sol_opt = model_opt.simulate(t_span=(0, 120))
    health_opt = model_opt.compute_health_index(sol_opt)
    survival_opt = model_opt.compute_survival(sol_opt.t,
                                               model_opt.compute_mortality_rate(health_opt))
    
    # Individual interventions for comparison
    interventions_single = {
        'CR only': {'cr_factor': 1.4},
        'Rapamycin only': {'rapa_factor': 1.3},
        'NAD+ only': {'nad_factor': 1.2},
        'Senolytic only': {'gamma_sen': 0.05},
    }
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax = axes[0]
    ax.plot(sol_base.t, health_base, 'k-', label='No intervention', linewidth=2)
    for name, params_i in interventions_single.items():
        m = AgingHallmarksModel(params=params_i)
        s = m.simulate(t_span=(0, 120))
        h = m.compute_health_index(s)
        ax.plot(s.t, h, '--', label=name, linewidth=1.5)
    ax.plot(sol_opt.t, health_opt, 'r-', label='Optimized combination', linewidth=2.5)
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Health Index')
    ax.set_title('A) Health Index: Single vs. Combined Interventions')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 120)
    
    ax = axes[1]
    ax.plot(sol_base.t, survival_base, 'k-', label='No intervention', linewidth=2)
    for name, params_i in interventions_single.items():
        m = AgingHallmarksModel(params=params_i)
        s = m.simulate(t_span=(0, 120))
        h = m.compute_health_index(s)
        surv = m.compute_survival(s.t, m.compute_mortality_rate(h))
        ax.plot(s.t, surv, '--', label=name, linewidth=1.5)
    ax.plot(sol_opt.t, survival_opt, 'r-', label='Optimized combination', linewidth=2.5)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Survival Probability')
    ax.set_title('B) Survival Curves: Single vs. Combined Interventions')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig7_optimization_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig7_optimization_comparison.png")
    
    # ---- Summary figure: Sensitivity analysis ----
    print("Running Simulation 8: Sensitivity Analysis...")
    
    param_ranges = {
        'CR level': np.linspace(0, 0.6, 20),
        'Rapamycin dose': np.linspace(0, 1.0, 20),
        'NAD+ dose': np.linspace(0, 1.0, 20),
        'Senolytic rate': np.linspace(0, 0.1, 20),
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    for idx, (pname, prange) in enumerate(param_ranges.items()):
        ax = axes.flat[idx]
        healthspans = []
        for val in prange:
            if pname == 'CR level':
                params = {'cr_factor': 1.0 + val}
            elif pname == 'Rapamycin dose':
                params = {'rapa_factor': 1.0 + val * 0.5}
            elif pname == 'NAD+ dose':
                params = {'nad_factor': 1.0 + val * 0.3}
            else:
                params = {'gamma_sen': val}
            
            m = AgingHallmarksModel(params=params)
            s = m.simulate(t_span=(0, 120))
            h = m.compute_health_index(s)
            below = np.where(h < 0.5)[0]
            hs = s.t[below[0]] if len(below) > 0 else 120
            healthspans.append(hs)
        
        ax.plot(prange, healthspans, 'b-o', markersize=4, linewidth=2)
        ax.set_xlabel(pname)
        ax.set_ylabel('Healthspan (years)')
        ax.set_title(f'Sensitivity: {pname}')
    
    plt.suptitle('Single-Parameter Sensitivity Analysis', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/fig8_sensitivity_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig8_sensitivity_analysis.png")
    
    # Save numerical results
    print("\n=== KEY RESULTS ===")
    print(f"Allometric model R²: {results['allometric']['R_squared']:.3f}")
    print(f"Allometric coefficients: a={results['allometric']['a']:.3f}, "
          f"b={results['allometric']['b']:.3f}, c={results['allometric']['c']:.3f}, "
          f"d={results['allometric']['d']:.3f}")
    print(f"Optimal intervention: CR={results['optimization']['cr_level']:.3f}, "
          f"Rapa={results['optimization']['rapa_dose']:.3f}, "
          f"NAD={results['optimization']['nad_dose']:.3f}, "
          f"Sen={results['optimization']['sen_rate']:.4f}")
    print(f"Optimized healthspan: {results['optimization']['healthspan']:.1f} years")
    
    with open('simulation_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    return results


if __name__ == '__main__':
    results = run_all_simulations()
    print("\nAll simulations completed successfully!")
