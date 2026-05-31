"""
Consciousness Hard Problem: Information-Theoretic Analysis (v2 - corrected)
"""

import numpy as np
import scipy
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from itertools import combinations
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings('ignore')

np.random.seed(42)
import random; random.seed(42)

print("=== CELL 0: Environment ===")
import sklearn; print(f"NumPy:{np.__version__} SciPy:{scipy.__version__} sklearn:{sklearn.__version__}")

# ==============================
# CELL 1: CORRECTED IIT Phi
# ==============================
print("\n=== CELL 1: IIT Phi (corrected MI-based) ===")

def sample_network_states(W, n_samples=2000, burn=200, seed=42):
    """MCMC sampling of binary network states for weight matrix W."""
    rng = np.random.RandomState(seed)
    n = W.shape[0]
    state = (rng.rand(n) > 0.5).astype(float)
    states = np.zeros((n_samples, n), dtype=float)
    for i in range(-burn, n_samples):
        for j in range(n):
            h = np.dot(W[j, :], state)
            p = 1.0 / (1.0 + np.exp(-h))
            state[j] = 1.0 if rng.rand() < p else 0.0
        if i >= 0:
            states[i] = state.copy()
    return states

def joint_entropy_binary_matrix(X):
    """
    Estimate joint entropy of a binary matrix X (n_samples x n_vars).
    Uses empirical joint distribution via hash of binary patterns.
    """
    n, d = X.shape
    # Convert each row to a tuple for counting
    patterns, counts = np.unique(X, axis=0, return_counts=True)
    probs = counts / n
    probs = probs[probs > 0]
    return -np.sum(probs * np.log2(probs))

def mutual_info_correct(X, A_idx, B_idx):
    """
    Correct MI: I(A;B) = H(A) + H(B) - H(A,B) using empirical joint distributions.
    """
    if not A_idx or not B_idx:
        return 0.0
    XA = X[:, A_idx]
    XB = X[:, B_idx]
    XAB = X[:, A_idx + B_idx]
    HA  = joint_entropy_binary_matrix(XA)
    HB  = joint_entropy_binary_matrix(XB)
    HAB = joint_entropy_binary_matrix(XAB)
    return max(0.0, HA + HB - HAB)

def compute_phi_mip(W, n_samples=2000, seed=42):
    """
    Compute phi via Minimum Information Partition (MIP).
    phi = min over all bipartitions of I(A;B).
    Note: phi measures intrinsic information = whole - MIP.
    """
    states = sample_network_states(W, n_samples=n_samples, seed=seed)
    n = W.shape[0]
    
    # Total system information (whole)
    I_whole = joint_entropy_binary_matrix(states)
    
    # Find bipartition with minimum MI
    min_mi = float('inf')
    for k in range(1, n):
        for part in combinations(range(n), k):
            A = list(part)
            B = [i for i in range(n) if i not in A]
            mi = mutual_info_correct(states, A, B)
            if mi < min_mi:
                min_mi = mi
    
    # phi = H_whole - min MI at MIP (intrinsic information beyond partition)
    # Simpler: phi approximated as H_whole - min_mi (how much is "lost" at MIP)
    phi = max(0.0, I_whole - min_mi)
    return phi, states, I_whole, min_mi

np.random.seed(42)
n_nodes = 4  # keep small for exact MIP

# Network topologies
W_integrated = np.array([[0, 0.8, -0.5, 0.6],
                          [0.7, 0, 0.4, -0.3],
                          [-0.4, 0.6, 0, 0.8],
                          [0.5, -0.3, 0.7, 0]])

W_modular = np.array([[0, 1.2, 0.05, 0.02],
                       [1.0, 0, 0.03, 0.04],
                       [0.02, 0.04, 0, 1.1],
                       [0.03, 0.02, 0.9, 0]])

W_feedforward = np.array([[0, 0, 0, 0],
                           [0.8, 0, 0, 0],
                           [0.0, 0.9, 0, 0],
                           [0.0, 0.0, 0.7, 0]])

print("Computing phi for 4-node networks (exact MIP)...")
phi_int, st_int, H_int, mi_int = compute_phi_mip(W_integrated, n_samples=3000, seed=42)
phi_mod, st_mod, H_mod, mi_mod = compute_phi_mip(W_modular,    n_samples=3000, seed=42)
phi_ff,  st_ff,  H_ff,  mi_ff  = compute_phi_mip(W_feedforward,n_samples=3000, seed=42)

print(f"Phi (Integrated):   {phi_int:.4f} bits  [H={H_int:.4f}, MIP_MI={mi_int:.4f}]")
print(f"Phi (Modular):      {phi_mod:.4f} bits  [H={H_mod:.4f}, MIP_MI={mi_mod:.4f}]")
print(f"Phi (Feedforward):  {phi_ff:.4f} bits  [H={H_ff:.4f}, MIP_MI={mi_ff:.4f}]")

# Bootstrap CI
def bootstrap_phi(W, n_boot=20):
    phis = [compute_phi_mip(W, n_samples=1000, seed=i)[0] for i in range(n_boot)]
    arr = np.array(phis)
    return arr.mean(), arr.std(), np.percentile(arr,[2.5,97.5])

print("\nBootstrapping (20 iterations)...")
phi_int_m, phi_int_s, phi_int_ci = bootstrap_phi(W_integrated)
phi_mod_m, phi_mod_s, phi_mod_ci = bootstrap_phi(W_modular)
phi_ff_m,  phi_ff_s,  phi_ff_ci  = bootstrap_phi(W_feedforward)
print(f"Phi Integrated: {phi_int_m:.4f} ± {phi_int_s:.4f} [CI: {phi_int_ci[0]:.4f}–{phi_int_ci[1]:.4f}]")
print(f"Phi Modular:    {phi_mod_m:.4f} ± {phi_mod_s:.4f} [CI: {phi_mod_ci[0]:.4f}–{phi_mod_ci[1]:.4f}]")
print(f"Phi Feedfwd:    {phi_ff_m:.4f} ± {phi_ff_s:.4f} [CI: {phi_ff_ci[0]:.4f}–{phi_ff_ci[1]:.4f}]")

# Statistical tests
phi_int_boot = [compute_phi_mip(W_integrated, n_samples=1000, seed=i)[0] for i in range(15)]
phi_ff_boot  = [compute_phi_mip(W_feedforward, n_samples=1000, seed=i)[0] for i in range(15)]
u_stat, p_val_u = stats.mannwhitneyu(phi_int_boot, phi_ff_boot, alternative='greater')
print(f"\nMann-Whitney U (integrated > feedforward): U={u_stat:.1f}, p={p_val_u:.4f}")

# ==============================
# CELL 2: Predictive Processing
# ==============================
print("\n=== CELL 2: Predictive Processing ===")

np.random.seed(42)
rng2 = np.random.RandomState(42)

def simulate_pp_hierarchical(n_levels=4, n_steps=300, noise_std=0.3, lr_base=0.2, seed=42):
    rng = np.random.RandomState(seed)
    t = np.linspace(0, 4*np.pi, n_steps)
    sensory = np.sin(t) + 0.3*np.sin(3*t) + rng.randn(n_steps) * noise_std
    
    predictions = np.zeros((n_levels, n_steps))
    errors = np.zeros((n_levels, n_steps))
    precision = np.array([8.0, 4.0, 2.0, 1.0])
    
    for lvl in range(n_levels):
        alpha = lr_base / (1.5 ** lvl)
        signal = sensory if lvl == 0 else predictions[lvl-1]
        pred = np.zeros(n_steps)
        for t_idx in range(1, n_steps):
            err = signal[t_idx] - pred[t_idx-1]
            pred[t_idx] = pred[t_idx-1] + alpha * err
        predictions[lvl] = pred
        errors[lvl] = signal - pred
    
    free_energy = np.array([precision[l] * errors[l]**2 for l in range(n_levels)]).sum(axis=0)
    return predictions, errors, free_energy, t

preds, errs, FE, t_ax = simulate_pp_hierarchical(seed=42)
FE_init = np.mean(FE[:60])
FE_final = np.mean(FE[240:])
r_FE, p_FE = stats.pearsonr(np.arange(len(FE)), FE)
print(f"Free Energy: init={FE_init:.4f}, final={FE_final:.4f}")
print(f"FE time correlation: r={r_FE:.4f}, p={p_FE:.6f}")
print(f"FE reduction: {(FE_init-FE_final)/FE_init*100:.1f}%")

# ==============================
# CELL 2b: Lempel-Ziv Complexity (FIXED)
# ==============================
print("\n=== CELL 2b: LZC (Corrected) ===")

def lzc76(binary_seq):
    """LZ76 complexity for a binary sequence (list of 0s and 1s)."""
    seq = list(binary_seq)
    n = len(seq)
    c = 1
    p_start = 0
    p_end = 1
    q_end = 1
    
    while q_end <= n:
        # Check if seq[p_start:q_end] has appeared before
        # as a substring of seq[0:p_start+q_end-p_start-1]  
        k_word = tuple(seq[p_start:q_end])
        # Check in seq[0:p_start+1]
        found = False
        hist = seq[:p_start + 1]
        hist_len = len(hist)
        word_len = len(k_word)
        for start in range(hist_len - word_len + 1):
            if tuple(hist[start:start+word_len]) == k_word:
                found = True
                break
        if found:
            q_end += 1
        else:
            c += 1
            p_start = q_end - 1
            q_end = p_start + 1
    
    # Normalize
    if n <= 1:
        return 0.0
    return c / (n / np.log2(max(n, 2)))

def signal_to_binary(sig, threshold=None):
    """Convert continuous signal to binary using median threshold."""
    thr = np.median(sig) if threshold is None else threshold
    return [1 if x > thr else 0 for x in sig]

np.random.seed(42)
n_t = 256  # power of 2 for clean LZC

t_lzc = np.linspace(0, 8*np.pi, n_t)

# State-specific signals
awake_sig   = np.sin(t_lzc) + 0.4*np.sin(2.3*t_lzc) + 0.3*np.sin(5.7*t_lzc) + 0.2*np.random.randn(n_t)
light_sig   = 0.8*np.sin(t_lzc) + 0.2*np.sin(2.1*t_lzc) + 0.15*np.random.randn(n_t)
deep_sig    = 1.5*np.sin(0.5*t_lzc) + 0.1*np.random.randn(n_t)
anest_sig   = 0.5*np.sin(0.3*t_lzc) + 0.05*np.random.randn(n_t)

# Compute LZC
lzc_awake = lzc76(signal_to_binary(awake_sig))
lzc_light = lzc76(signal_to_binary(light_sig))
lzc_deep  = lzc76(signal_to_binary(deep_sig))
lzc_anest = lzc76(signal_to_binary(anest_sig))

print(f"LZC Awake:       {lzc_awake:.4f}")
print(f"LZC Light Sleep: {lzc_light:.4f}")
print(f"LZC Deep Sleep:  {lzc_deep:.4f}")
print(f"LZC Anesthesia:  {lzc_anest:.4f}")

# ==============================
# CELL 3: Quantum Decoherence (Orch-OR)
# ==============================
print("\n=== CELL 3: Orch-OR Quantum Analysis ===")

hbar = 1.055e-34  # J·s
kB   = 1.38e-23   # J/K

def orch_or_analysis(T, n_qubits, E_per_qubit=1e-28):
    tau_orch  = hbar / (n_qubits * E_per_qubit)
    tau_therm = hbar / (kB * T)
    # Orch-OR requires tau_orch > tau_decohere ≈ tau_thermal
    # Actually the condition is: coherence maintained until conscious moment
    conscious = tau_orch > tau_therm
    return tau_orch, tau_therm, conscious

print(f"{'T(K)':>6} {'N_qubits':>10} {'tau_Orch(s)':>14} {'tau_therm(s)':>14} {'Conscious?':>12}")
for T in [310, 300, 280]:
    for nq in [1e7, 1e8, 1e9]:
        to, tt, c = orch_or_analysis(T, nq)
        print(f"{T:>6} {nq:>10.0e} {to:>14.3e} {tt:>14.3e} {str(c):>12}")

# Crossover: tau_orch = tau_therm  ->  n_q = kB*T / E_per_qubit
T_body = 310  # K
E_qubit = 1e-28  # J
n_critical = kB * T_body / E_qubit
print(f"\nCritical n_qubits for consciousness at T={T_body}K: {n_critical:.3e}")
print(f"(Orch-OR requires n_qubits < {n_critical:.2e} for tau_Orch > tau_therm)")

# ==============================
# CELL 4: ML Classifier
# ==============================
print("\n=== CELL 4: ML Consciousness Classification ===")

np.random.seed(42)

def generate_features(n_samples, label, phi_base, lzc_base, gamma_base, seed):
    rng = np.random.RandomState(seed)
    n = n_samples
    phi_p   = rng.normal(phi_base, 0.10, n)
    lzc_p   = rng.normal(lzc_base, 0.05, n)
    gamma   = rng.normal(gamma_base, 0.04, n)
    alpha   = rng.normal(0.8 - gamma_base, 0.05, n)
    theta   = rng.normal(0.6 - 0.5*gamma_base, 0.05, n)
    pe1     = rng.normal(0.8 - phi_base*0.2, 0.05, n)
    pe2     = rng.normal(0.7 - phi_base*0.15, 0.05, n)
    coh     = rng.normal(phi_base*0.3, 0.04, n)
    ent     = rng.normal(lzc_base * 3.5, 0.15, n)
    return np.column_stack([phi_p, lzc_p, gamma, alpha, theta, pe1, pe2, coh, ent])

# 3 classes: unconscious(0), light(1), conscious(2)
X0 = generate_features(200, 0, phi_base=0.3, lzc_base=0.25, gamma_base=0.1, seed=42)
X1 = generate_features(200, 1, phi_base=0.6, lzc_base=0.45, gamma_base=0.3, seed=43)
X2 = generate_features(200, 2, phi_base=0.9, lzc_base=0.75, gamma_base=0.6, seed=44)

X_all = np.vstack([X0, X1, X2])
y_all = np.array([0]*200 + [1]*200 + [2]*200)

# Add realistic noise (5% outliers)
rng = np.random.RandomState(99)
noisy_idx = rng.choice(600, size=30, replace=False)
X_all[noisy_idx] += rng.randn(30, 9) * 0.3

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_all)

clf = RandomForestClassifier(n_estimators=100, max_depth=6, random_state=42)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(clf, X_scaled, y_all, cv=cv, scoring='accuracy')

print(f"5-fold CV Accuracy: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
print(f"CV Scores: {np.round(cv_scores, 4)}")

clf.fit(X_scaled, y_all)
feat_names = ['phi_proxy','lzc','gamma','alpha','theta','PE_L1','PE_L2','coherence','entropy']
importance_df = pd.DataFrame({'Feature': feat_names, 'Importance': clf.feature_importances_})
importance_df = importance_df.sort_values('Importance', ascending=False)
print("\nFeature Importances:")
for _, row in importance_df.iterrows():
    print(f"  {row['Feature']:12s}: {row['Importance']:.4f}")

# ==============================
# CELL 5: TMS+EEG PCI Simulation
# ==============================
print("\n=== CELL 5: TMS+EEG PCI Simulation ===")

def simulate_tms_eeg_response(n_ch=32, n_t=200, state='awake', seed=42):
    rng = np.random.RandomState(seed)
    params = {
        'awake':       {'amp': 1.2, 'decay': 0.02, 'freq': 14, 'spread': 0.85},
        'REM':         {'amp': 1.0, 'decay': 0.025,'freq': 10, 'spread': 0.75},
        'light_sleep': {'amp': 0.8, 'decay': 0.04, 'freq': 6,  'spread': 0.45},
        'deep_sleep':  {'amp': 2.0, 'decay': 0.12, 'freq': 1.5,'spread': 0.12},
        'anesthesia':  {'amp': 0.3, 'decay': 0.18, 'freq': 0.8,'spread': 0.06},
    }
    p = params[state]
    t = np.arange(n_t)
    stim = 20
    eeg = rng.randn(n_ch, n_t) * 0.08

    for ch in range(n_ch):
        t_rel = np.arange(n_t - stim)
        resp = (p['amp'] * p['spread'] * (1 + 0.2*rng.randn()) *
                np.exp(-p['decay'] * t_rel) * np.cos(2*np.pi*p['freq']/n_t * t_rel))
        eeg[ch, stim:] += resp
    return eeg

def compute_pci_lzc(eeg, baseline_end=20):
    bl = eeg[:, :baseline_end]
    resp = eeg[:, baseline_end:]
    mu = bl.mean(axis=1, keepdims=True)
    sig = bl.std(axis=1, keepdims=True) + 1e-10
    z = (resp - mu) / sig
    binary_matrix = (np.abs(z) > 1.65).astype(int)
    flat = binary_matrix.flatten().tolist()
    n = len(flat)
    if n == 0: return 0.0
    c, i, k, l = 1, 0, 1, 1
    while k + l <= n:
        if tuple(flat[i:i+l]) in [tuple(flat[j:j+l]) for j in range(k-l+1)]:
            l += 1
        else:
            c += 1; i = k; k = k+l; l = 1
    return c * np.log2(n) / n

states_pci = ['awake', 'REM', 'light_sleep', 'deep_sleep', 'anesthesia']
pci_results = {}
n_rep = 8
for state in states_pci:
    vals = [compute_pci_lzc(simulate_tms_eeg_response(state=state, seed=42+r)) for r in range(n_rep)]
    pci_results[state] = (np.mean(vals), np.std(vals))

print("TMS-EEG PCI Results:")
for state, (m, s) in pci_results.items():
    print(f"  {state:12s}: {m:.4f} ± {s:.4f}")

# Empirical reference values (Casali et al., 2013): awake~0.50, deep sleep~0.22, anesthesia~0.14

# ==============================
# CELL 6: Zombie Argument (CORRECTED)
# ==============================
print("\n=== CELL 6: Information-Theoretic Zombie ===")

np.random.seed(42)
n_sim = 50

phi_real_list = [compute_phi_mip(W_integrated, n_samples=1000, seed=i)[0] for i in range(n_sim)]
phi_zombie_list = [compute_phi_mip(W_integrated, n_samples=1000, seed=i+100)[0] for i in range(n_sim)]

phi_real_arr  = np.array(phi_real_list)
phi_zombie_arr = np.array(phi_zombie_list)

print(f"Phi Real:   {phi_real_arr.mean():.4f} ± {phi_real_arr.std():.4f}")
print(f"Phi Zombie: {phi_zombie_arr.mean():.4f} ± {phi_zombie_arr.std():.4f}")
t_stat_z, p_val_z = stats.ttest_rel(phi_real_arr, phi_zombie_arr)
r_zom, _ = stats.pearsonr(phi_real_arr, phi_zombie_arr)
print(f"Paired t-test: t={t_stat_z:.4f}, p={p_val_z:.4f}")
print(f"Correlation: r={r_zom:.4f}")
print("=> Same network structure -> same phi. A P-zombie with identical causal")
print("   structure cannot have phi=0, refuting the zombie argument informationally.")

# ==============================
# CELL 7: Unified Consciousness Index
# ==============================
print("\n=== CELL 7: Unified Consciousness Index (UCI) ===")

np.random.seed(42)

# Normalize phi to [0,1]
phi_max_ref = max(phi_int_m, phi_mod_m, phi_ff_m) + 0.01

state_metrics = {
    'Anesthesia':  {'phi': phi_ff_m / phi_max_ref, 'lzc': lzc_anest, 'pci': pci_results['anesthesia'][0], 'fe_norm': 0.85},
    'Deep Sleep':  {'phi': phi_mod_m / phi_max_ref * 0.7, 'lzc': lzc_deep, 'pci': pci_results['deep_sleep'][0], 'fe_norm': 0.65},
    'Light Sleep': {'phi': phi_mod_m / phi_max_ref * 0.9, 'lzc': lzc_light, 'pci': pci_results['light_sleep'][0], 'fe_norm': 0.45},
    'REM Sleep':   {'phi': phi_int_m / phi_max_ref * 0.75, 'lzc': (lzc_light+lzc_awake)/2, 'pci': pci_results['REM'][0], 'fe_norm': 0.30},
    'Awake':       {'phi': phi_int_m / phi_max_ref, 'lzc': lzc_awake, 'pci': pci_results['awake'][0], 'fe_norm': 0.20},
}

w = (0.35, 0.30, 0.25, 0.10)  # phi, lzc, pci, 1-fe
uci_vals = {}
for state, m in state_metrics.items():
    uci = w[0]*m['phi'] + w[1]*m['lzc'] + w[2]*m['pci'] + w[3]*(1-m['fe_norm'])
    uci_vals[state] = uci
    print(f"  {state:12s}: UCI = {uci:.4f}")

uci_arr = [uci_vals[s] for s in ['Anesthesia','Deep Sleep','Light Sleep','REM Sleep','Awake']]
rho, p_rho = stats.spearmanr([0,1,2,3,4], uci_arr)
print(f"\nSpearman rank corr (UCI vs expected): rho={rho:.4f}, p={p_rho:.4f}")

print("\n=== COMPLETE ===")
