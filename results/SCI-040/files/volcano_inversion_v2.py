"""
Volcanic Crustal Deformation Inversion Framework
=================================================
Bayesian inversion for 3D magma supply system structure.
Implements Mogi, spheroid, and FEM-based forward models with
joint GNSS+InSAR+gravity inversion and Kalman filter tracking.
Uses emcee for MCMC sampling.
"""

import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import emcee
import warnings
import os
import json
from multiprocessing import Pool

warnings.filterwarnings('ignore')
np.random.seed(42)

FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# 1. Forward Models
# ============================================================

def mogi_forward(x_obs, y_obs, xs, ys, d, dV, nu=0.25):
    """Mogi point pressure source in elastic half-space."""
    dx = x_obs - xs
    dy = y_obs - ys
    R = np.sqrt(dx**2 + dy**2 + d**2)
    C = (1 - nu) * dV / np.pi
    ux = C * dx / R**3
    uy = C * dy / R**3
    uz = C * d / R**3
    rho = 2700.0
    G_grav = 6.674e-11
    dg = -2 * G_grav * rho * dV / R + 0.3086e-5 * uz
    return ux, uy, uz, dg


def spheroid_forward(x_obs, y_obs, xs, ys, d, a, b, dP, mu=3e10, nu=0.25):
    """Prolate spheroid pressure source (Yang et al., 1988 approx)."""
    dx = x_obs - xs
    dy = y_obs - ys
    r = np.sqrt(dx**2 + dy**2)
    R = np.sqrt(r**2 + d**2)
    A = a / b if a > b else b / a
    dV = (4.0/3.0) * np.pi * a * b**2 * dP / (mu * (1 + (A-1)*0.5))
    f_geom = 1.0 + 0.5 * (A - 1) * (d**2 / R**2)
    C = (1 - nu) * dV / np.pi
    ux = C * dx / R**3 * f_geom
    uy = C * dy / R**3 * f_geom
    uz = C * d / R**3 * f_geom
    rho = 2700.0
    G_grav = 6.674e-11
    dg = -2 * G_grav * rho * dV / R + 0.3086e-5 * uz
    return ux, uy, uz, dg


def fem_forward_simplified(x_obs, y_obs, xs, ys, d, dV,
                           mu=3e10, nu=0.25, eta=1e19, t=0):
    """Simplified FEM-like model with Maxwell viscoelastic relaxation."""
    ux_e, uy_e, uz_e, dg_e = mogi_forward(x_obs, y_obs, xs, ys, d, dV, nu)
    tau = eta / mu
    if t > 0 and tau > 0:
        relax = 1.0 + (1.0/(2*(1-nu))) * (1 - np.exp(-t/tau))
        ux = ux_e * relax
        uy = uy_e * relax
        uz = uz_e * relax
    else:
        ux, uy, uz = ux_e, uy_e, uz_e
    rho = 2700.0
    G_grav = 6.674e-11
    dg = -2 * G_grav * rho * dV / np.sqrt((x_obs-xs)**2+(y_obs-ys)**2+d**2) \
         + 0.3086e-5 * uz
    return ux, uy, uz, dg


# ============================================================
# 2. Synthetic Data Generation
# ============================================================

def generate_station_network(n_gnss=15, n_insar=400, n_grav=8, extent=20000):
    rng = np.random.RandomState(42)
    gnss_x = rng.uniform(-extent, extent, n_gnss)
    gnss_y = rng.uniform(-extent, extent, n_gnss)
    n_side = int(np.sqrt(n_insar))
    ix = np.linspace(-extent, extent, n_side)
    iy = np.linspace(-extent, extent, n_side)
    insar_x, insar_y = np.meshgrid(ix, iy)
    insar_x, insar_y = insar_x.ravel(), insar_y.ravel()
    grav_x = rng.uniform(-extent*0.5, extent*0.5, n_grav)
    grav_y = rng.uniform(-extent*0.5, extent*0.5, n_grav)
    return {'gnss': (gnss_x, gnss_y), 'insar': (insar_x, insar_y),
            'gravity': (grav_x, grav_y)}


def generate_synthetic_data(true_params, stations, model='mogi',
                            sigma_gnss=0.003, sigma_insar=0.005,
                            sigma_grav=5e-8):
    rng = np.random.RandomState(123)
    data = {}
    for dtype, (xo, yo) in stations.items():
        if model == 'mogi':
            ux, uy, uz, dg = mogi_forward(xo, yo, **true_params)
        elif model == 'spheroid':
            ux, uy, uz, dg = spheroid_forward(xo, yo, **true_params)
        elif model == 'fem':
            ux, uy, uz, dg = fem_forward_simplified(xo, yo, **true_params)
        n = len(xo)
        if dtype == 'gnss':
            data[dtype] = {'ux': ux + rng.normal(0, sigma_gnss, n),
                           'uy': uy + rng.normal(0, sigma_gnss, n),
                           'uz': uz + rng.normal(0, sigma_gnss, n),
                           'sigma': sigma_gnss}
        elif dtype == 'insar':
            los = -0.07*ux + 0.39*uy + 0.92*uz
            data[dtype] = {'los': los + rng.normal(0, sigma_insar, n),
                           'sigma': sigma_insar}
        elif dtype == 'gravity':
            data[dtype] = {'dg': dg + rng.normal(0, sigma_grav, n),
                           'sigma': sigma_grav}
    return data


# ============================================================
# 3. Bayesian MCMC Inversion (emcee)
# ============================================================

def log_prior_mogi(theta):
    xs, ys, d, dV = theta
    if not (-15000 < xs < 15000): return -np.inf
    if not (-15000 < ys < 15000): return -np.inf
    if not (500 < d < 20000): return -np.inf
    if not (1e4 < dV < 5e7): return -np.inf
    # Gaussian priors
    lp = -0.5 * ((xs/5000)**2 + (ys/5000)**2)
    lp += -0.5 * ((dV - 5e6)/(3e6))**2
    return lp


def log_likelihood_mogi(theta, stations, data):
    xs, ys, d, dV = theta
    ll = 0.0
    # GNSS
    gx, gy = stations['gnss']
    ux, uy, uz, _ = mogi_forward(gx, gy, xs, ys, d, dV)
    sig = data['gnss']['sigma']
    ll += -0.5 * np.sum(((data['gnss']['ux'] - ux)/sig)**2)
    ll += -0.5 * np.sum(((data['gnss']['uy'] - uy)/sig)**2)
    ll += -0.5 * np.sum(((data['gnss']['uz'] - uz)/sig)**2)
    # InSAR
    ix, iy = stations['insar']
    uxi, uyi, uzi, _ = mogi_forward(ix, iy, xs, ys, d, dV)
    los_pred = -0.07*uxi + 0.39*uyi + 0.92*uzi
    sig_i = data['insar']['sigma']
    ll += -0.5 * np.sum(((data['insar']['los'] - los_pred)/sig_i)**2)
    # Gravity
    grx, gry = stations['gravity']
    _, _, _, dg_pred = mogi_forward(grx, gry, xs, ys, d, dV)
    sig_g = data['gravity']['sigma']
    ll += -0.5 * np.sum(((data['gravity']['dg'] - dg_pred)/sig_g)**2)
    return ll


def log_posterior_mogi(theta, stations, data):
    lp = log_prior_mogi(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_mogi(theta, stations, data)


def run_mcmc_mogi(stations, data, n_walkers=32, n_steps=2000, n_burn=500):
    ndim = 4
    # Initialize near true values with scatter
    p0_center = np.array([0.0, 0.0, 4000.0, 4e6])
    p0 = p0_center + np.random.randn(n_walkers, ndim) * np.array([500, 500, 500, 5e5])
    # Ensure valid initial points
    for i in range(n_walkers):
        p0[i, 2] = max(600, p0[i, 2])
        p0[i, 3] = max(2e4, p0[i, 3])

    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_posterior_mogi,
                                     args=(stations, data))
    sampler.run_mcmc(p0, n_steps, progress=True)
    samples = sampler.get_chain(discard=n_burn, flat=True)
    chain = sampler.get_chain(discard=n_burn)
    return samples, chain, sampler


def log_prior_spheroid(theta):
    xs, ys, d, a, b, dP = theta
    if not (-15000 < xs < 15000): return -np.inf
    if not (-15000 < ys < 15000): return -np.inf
    if not (500 < d < 20000): return -np.inf
    if not (200 < a < 8000): return -np.inf
    if not (100 < b < 5000): return -np.inf
    if not (1e5 < dP < 5e8): return -np.inf
    return 0.0


def log_likelihood_spheroid(theta, stations, data):
    xs, ys, d, a, b, dP = theta
    ll = 0.0
    gx, gy = stations['gnss']
    ux, uy, uz, _ = spheroid_forward(gx, gy, xs, ys, d, a, b, dP)
    sig = data['gnss']['sigma']
    ll += -0.5 * np.sum(((data['gnss']['ux'] - ux)/sig)**2)
    ll += -0.5 * np.sum(((data['gnss']['uy'] - uy)/sig)**2)
    ll += -0.5 * np.sum(((data['gnss']['uz'] - uz)/sig)**2)
    ix, iy = stations['insar']
    uxi, uyi, uzi, _ = spheroid_forward(ix, iy, xs, ys, d, a, b, dP)
    los_pred = -0.07*uxi + 0.39*uyi + 0.92*uzi
    sig_i = data['insar']['sigma']
    ll += -0.5 * np.sum(((data['insar']['los'] - los_pred)/sig_i)**2)
    return ll


def log_posterior_spheroid(theta, stations, data):
    lp = log_prior_spheroid(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_spheroid(theta, stations, data)


def run_mcmc_spheroid(stations, data, n_walkers=48, n_steps=2000, n_burn=500):
    ndim = 6
    p0_center = np.array([0.0, 0.0, 4000.0, 1500.0, 800.0, 1e7])
    p0 = p0_center + np.random.randn(n_walkers, ndim) * \
         np.array([500, 500, 500, 200, 100, 2e6])
    for i in range(n_walkers):
        p0[i, 2] = max(600, p0[i, 2])
        p0[i, 3] = max(300, p0[i, 3])
        p0[i, 4] = max(150, p0[i, 4])
        p0[i, 5] = max(2e5, p0[i, 5])

    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_posterior_spheroid,
                                     args=(stations, data))
    sampler.run_mcmc(p0, n_steps, progress=True)
    samples = sampler.get_chain(discard=n_burn, flat=True)
    chain = sampler.get_chain(discard=n_burn)
    return samples, chain, sampler


# ============================================================
# 4. Kalman Filter for Time-Varying Source
# ============================================================

class VolcanicKalmanFilter:
    def __init__(self, x0, P0, Q, R_gnss, stations):
        self.x = np.array(x0, dtype=float)
        self.P = np.array(P0, dtype=float)
        self.Q = np.array(Q, dtype=float)
        self.R_gnss = R_gnss
        self.stations = stations
        self.history = {'x': [], 'P': [], 'innovation': []}

    def predict(self, dt=1.0, dV_rate=0.0):
        F = np.eye(4)
        self.x[3] += dV_rate * dt
        self.P = F @ self.P @ F.T + self.Q * dt

    def _jacobian(self, x_obs, y_obs, state):
        eps = [10.0, 10.0, 10.0, 1e3]
        n = len(x_obs)
        H = np.zeros((3*n, 4))
        for j in range(4):
            sp, sm = state.copy(), state.copy()
            sp[j] += eps[j]; sm[j] -= eps[j]
            uxp, uyp, uzp, _ = mogi_forward(x_obs, y_obs, *sp)
            uxm, uym, uzm, _ = mogi_forward(x_obs, y_obs, *sm)
            H[:n, j] = (uxp - uxm) / (2*eps[j])
            H[n:2*n, j] = (uyp - uym) / (2*eps[j])
            H[2*n:, j] = (uzp - uzm) / (2*eps[j])
        return H

    def update_gnss(self, obs_ux, obs_uy, obs_uz):
        gnss_x, gnss_y = self.stations['gnss']
        n = len(gnss_x)
        ux, uy, uz, _ = mogi_forward(gnss_x, gnss_y, *self.x)
        z_pred = np.concatenate([ux, uy, uz])
        z_obs = np.concatenate([obs_ux, obs_uy, obs_uz])
        H = self._jacobian(gnss_x, gnss_y, self.x)
        R = np.eye(3*n) * self.R_gnss
        innov = z_obs - z_pred
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.solve(S, np.eye(S.shape[0]))
        self.x = self.x + K @ innov
        self.P = (np.eye(4) - K @ H) @ self.P
        self.history['x'].append(self.x.copy())
        self.history['P'].append(self.P.copy())
        self.history['innovation'].append(np.mean(np.abs(innov)))


# ============================================================
# 5. Viscoelastic Correction
# ============================================================

def compute_viscoelastic_correction(t_obs, d_elastic, mu=3e10,
                                     eta=1e19, nu=0.25):
    tau = eta / mu
    factor = 1.0 + (1.0/(2*(1-nu))) * (1.0 - np.exp(-t_obs/tau))
    return d_elastic * factor


# ============================================================
# 6. Main Experiment
# ============================================================

def run_full_experiment():
    results = {}
    print("=" * 60)
    print("Volcanic Crustal Deformation Inversion Framework")
    print("=" * 60)

    true_mogi = {'xs': 500, 'ys': -300, 'd': 5000, 'dV': 5e6}
    true_spheroid = {'xs': 500, 'ys': -300, 'd': 5000,
                     'a': 2000, 'b': 1000, 'dP': 1e7}

    stations = generate_station_network(n_gnss=15, n_insar=400, n_grav=8)

    # -------------------------------------------------------
    # Exp 1: Forward Model Comparison
    # -------------------------------------------------------
    print("\n[1] Forward Model Comparison...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    insar_x, insar_y = stations['insar']
    n_side = int(np.sqrt(len(insar_x)))

    models_fwd = {
        'Mogi': mogi_forward(insar_x, insar_y, 500, -300, 5000, 5e6),
        'Spheroid': spheroid_forward(insar_x, insar_y, 500, -300, 5000, 2000, 1000, 1e7),
        'FEM+Viscoelastic': fem_forward_simplified(insar_x, insar_y, 500, -300, 5000, 5e6, t=3.15e7)
    }

    for ax, (name, (ux, uy, uz, dg)) in zip(axes, models_fwd.items()):
        uz_grid = uz.reshape(n_side, n_side)
        im = ax.imshow(uz_grid * 1000, extent=[-20, 20, -20, 20],
                      cmap='RdBu_r', origin='lower')
        ax.set_title(f'{name}\nmax uz = {np.max(uz)*1000:.2f} mm')
        ax.set_xlabel('X (km)'); ax.set_ylabel('Y (km)')
        plt.colorbar(im, ax=ax, label='Uz (mm)')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/model_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    results['model_comparison'] = {
        name: {'max_uz_mm': float(np.max(uz)*1000),
               'max_dg_uGal': float(np.max(np.abs(dg))*1e8)}
        for name, (ux, uy, uz, dg) in models_fwd.items()
    }
    print("  Done.")

    # -------------------------------------------------------
    # Exp 2: Bayesian MCMC Inversion (Mogi)
    # -------------------------------------------------------
    print("\n[2] Bayesian MCMC Inversion (Mogi)...")
    data_mogi = generate_synthetic_data(true_mogi, stations, model='mogi')
    samples_mogi, chain_mogi, sampler_mogi = run_mcmc_mogi(
        stations, data_mogi, n_walkers=32, n_steps=2000, n_burn=500)

    param_names = ['xs', 'ys', 'd', 'dV']
    post_mean = {p: float(np.mean(samples_mogi[:, i]))
                 for i, p in enumerate(param_names)}
    post_std = {p: float(np.std(samples_mogi[:, i]))
                for i, p in enumerate(param_names)}
    post_median = {p: float(np.median(samples_mogi[:, i]))
                   for i, p in enumerate(param_names)}

    results['mogi_inversion'] = {
        'posterior_mean': post_mean, 'posterior_std': post_std,
        'posterior_median': post_median, 'true_values': true_mogi,
        'n_samples': len(samples_mogi)
    }

    print(f"  True:  xs={true_mogi['xs']}, ys={true_mogi['ys']}, d={true_mogi['d']}, dV={true_mogi['dV']:.2e}")
    print(f"  Mean:  xs={post_mean['xs']:.0f}±{post_std['xs']:.0f}, ys={post_mean['ys']:.0f}±{post_std['ys']:.0f}, "
          f"d={post_mean['d']:.0f}±{post_std['d']:.0f}, dV={post_mean['dV']:.2e}±{post_std['dV']:.2e}")

    # Trace plot
    fig, axes = plt.subplots(4, 2, figsize=(14, 10))
    for i, pname in enumerate(param_names):
        axes[i, 0].plot(chain_mogi[:, :, i], alpha=0.3, lw=0.5)
        axes[i, 0].set_ylabel(pname)
        axes[i, 0].axhline(list(true_mogi.values())[i], color='r', ls='--', lw=2)
        axes[i, 1].hist(samples_mogi[:, i], bins=50, density=True, alpha=0.7, color='steelblue')
        axes[i, 1].axvline(list(true_mogi.values())[i], color='r', ls='--', lw=2, label='True')
        axes[i, 1].axvline(post_mean[pname], color='k', ls='-', lw=1.5, label='Mean')
        axes[i, 1].legend(fontsize=8)
    axes[0, 0].set_title('Chain'); axes[0, 1].set_title('Posterior')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/mcmc_trace_mogi.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Corner plot
    fig, axes = plt.subplots(4, 4, figsize=(12, 12))
    labels = ['xs (m)', 'ys (m)', 'd (m)', 'dV (m³)']
    true_vals = list(true_mogi.values())
    for i in range(4):
        for j in range(4):
            ax = axes[i, j]
            if i == j:
                ax.hist(samples_mogi[:, i], bins=40, density=True, color='steelblue', alpha=0.7)
                ax.axvline(true_vals[i], color='r', ls='--', lw=1.5)
            elif i > j:
                ax.scatter(samples_mogi[::5, j], samples_mogi[::5, i],
                          s=1, alpha=0.1, color='steelblue')
                ax.axhline(true_vals[i], color='r', ls='--', lw=0.5, alpha=0.5)
                ax.axvline(true_vals[j], color='r', ls='--', lw=0.5, alpha=0.5)
            else:
                ax.axis('off')
            if i == 3: ax.set_xlabel(labels[j], fontsize=8)
            if j == 0: ax.set_ylabel(labels[i], fontsize=8)
            ax.tick_params(labelsize=6)
    plt.suptitle('Mogi Source Posterior Distribution', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/posterior_corner_mogi.png', dpi=150, bbox_inches='tight')
    plt.close()

    # -------------------------------------------------------
    # Exp 3: Joint Inversion Data Fit
    # -------------------------------------------------------
    print("\n[3] Joint Inversion Data Fit...")
    gnss_x, gnss_y = stations['gnss']
    ux_pred, uy_pred, uz_pred, _ = mogi_forward(
        gnss_x, gnss_y, post_mean['xs'], post_mean['ys'],
        post_mean['d'], post_mean['dV'])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    comps = [('East (Ux)', data_mogi['gnss']['ux'], ux_pred),
             ('North (Uy)', data_mogi['gnss']['uy'], uy_pred),
             ('Vertical (Uz)', data_mogi['gnss']['uz'], uz_pred)]
    residuals = {}
    for ax, (comp, obs, pred) in zip(axes, comps):
        ax.scatter(obs*1000, pred*1000, c='steelblue', edgecolor='k', s=60)
        lim = max(np.max(np.abs(obs)), np.max(np.abs(pred))) * 1000 * 1.2
        ax.plot([-lim, lim], [-lim, lim], 'r--', lw=1.5)
        ax.set_xlabel('Observed (mm)'); ax.set_ylabel('Predicted (mm)')
        ax.set_title(f'GNSS {comp}')
        rms = np.sqrt(np.mean((obs - pred)**2)) * 1000
        ax.text(0.05, 0.95, f'RMS = {rms:.3f} mm', transform=ax.transAxes, va='top')
        residuals[comp] = rms
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/joint_inversion_fit.png', dpi=150, bbox_inches='tight')
    plt.close()

    # InSAR fit
    insar_x, insar_y = stations['insar']
    ux_i, uy_i, uz_i, _ = mogi_forward(insar_x, insar_y,
        post_mean['xs'], post_mean['ys'], post_mean['d'], post_mean['dV'])
    los_pred = -0.07*ux_i + 0.39*uy_i + 0.92*uz_i
    los_obs = data_mogi['insar']['los']

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    n_side = int(np.sqrt(len(insar_x)))
    vmax = max(np.max(np.abs(los_obs)), np.max(np.abs(los_pred))) * 1000
    for ax, (title, vals) in zip(axes, [('Observed LOS', los_obs),
                                         ('Predicted LOS', los_pred),
                                         ('Residual', los_obs - los_pred)]):
        grid = vals.reshape(n_side, n_side) * 1000
        if title == 'Residual':
            vm = np.max(np.abs(grid)) * 1.1
            im = ax.imshow(grid, extent=[-20,20,-20,20], cmap='RdBu_r',
                          origin='lower', vmin=-vm, vmax=vm)
        else:
            im = ax.imshow(grid, extent=[-20,20,-20,20], cmap='RdBu_r',
                          origin='lower', vmin=-vmax, vmax=vmax)
        ax.set_title(title); ax.set_xlabel('X (km)'); ax.set_ylabel('Y (km)')
        plt.colorbar(im, ax=ax, label='mm')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/insar_fit.png', dpi=150, bbox_inches='tight')
    plt.close()

    insar_rms = np.sqrt(np.mean((los_obs - los_pred)**2)) * 1000
    residuals['InSAR_LOS'] = insar_rms
    results['residuals_mm'] = {k: float(v) for k, v in residuals.items()}
    print(f"  GNSS RMS: E={residuals['East (Ux)']:.3f}, N={residuals['North (Uy)']:.3f}, U={residuals['Vertical (Uz)']:.3f} mm")
    print(f"  InSAR LOS RMS: {insar_rms:.3f} mm")

    # -------------------------------------------------------
    # Exp 4: Kalman Filter Time-Varying Source
    # -------------------------------------------------------
    print("\n[4] Kalman Filter for Time-Varying Source...")
    n_epochs = 50
    rng_kf = np.random.RandomState(99)
    true_dV_history = 5e6 + np.cumsum(rng_kf.normal(1e5, 5e4, n_epochs))

    x0 = [0, 0, 4000, 4e6]
    P0 = np.diag([1e6, 1e6, 1e6, 1e12])
    Q = np.diag([100, 100, 100, 1e8])

    kf = VolcanicKalmanFilter(x0, P0, Q, R_gnss=(0.003)**2, stations=stations)
    kf_true = [500, -300, 5000]
    gnss_x, gnss_y = stations['gnss']

    for epoch in range(n_epochs):
        kf.predict(dt=1.0, dV_rate=1e5)
        true_dV = true_dV_history[epoch]
        ux_t, uy_t, uz_t, _ = mogi_forward(gnss_x, gnss_y,
                                             kf_true[0], kf_true[1], kf_true[2], true_dV)
        noise = rng_kf.normal(0, 0.003, len(gnss_x))
        kf.update_gnss(ux_t + noise, uy_t + noise, uz_t + noise)

    kf_states = np.array(kf.history['x'])
    kf_covs = np.array(kf.history['P'])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    params_kf = ['xs (m)', 'ys (m)', 'Depth (m)', 'Volume Change (m³)']
    true_kf = [500, -300, 5000, None]
    for i, (ax, pname) in enumerate(zip(axes.flat, params_kf)):
        est = kf_states[:, i]
        std = np.sqrt(kf_covs[:, i, i])
        epochs = np.arange(n_epochs)
        ax.plot(epochs, est, 'b-', lw=2, label='KF Estimate')
        ax.fill_between(epochs, est - 2*std, est + 2*std,
                        alpha=0.2, color='blue', label='±2σ')
        if i == 3:
            ax.plot(epochs, true_dV_history, 'r--', lw=1.5, label='True')
        elif true_kf[i] is not None:
            ax.axhline(true_kf[i], color='r', ls='--', lw=1.5, label='True')
        ax.set_xlabel('Epoch'); ax.set_ylabel(pname)
        ax.set_title(pname); ax.legend(loc='best')
    plt.suptitle('Kalman Filter Source Tracking', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/kalman_filter_tracking.png', dpi=150, bbox_inches='tight')
    plt.close()

    results['kalman_filter'] = {
        'final_state': kf_states[-1].tolist(),
        'final_std': np.sqrt(np.diag(kf_covs[-1])).tolist(),
        'mean_innovation': float(np.mean(kf.history['innovation'])),
        'dV_rmse': float(np.sqrt(np.mean((kf_states[:, 3] - true_dV_history)**2)))
    }
    print(f"  KF final: xs={kf_states[-1,0]:.0f}, ys={kf_states[-1,1]:.0f}, "
          f"d={kf_states[-1,2]:.0f}, dV={kf_states[-1,3]:.2e}")
    print(f"  dV RMSE: {results['kalman_filter']['dV_rmse']:.2e} m³")

    # -------------------------------------------------------
    # Exp 5: Viscoelastic Correction
    # -------------------------------------------------------
    print("\n[5] Viscoelastic Correction...")
    t_years = np.linspace(0.01, 10, 100)
    t_seconds = t_years * 3.15e7
    uz_elastic = mogi_forward(np.array([0.0]), np.array([0.0]),
                               500, -300, 5000, 5e6)[2][0]
    viscosities = [1e18, 1e19, 1e20]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    ax.axhline(uz_elastic*1000, color='k', ls='--', lw=1, label='Elastic')
    for eta in viscosities:
        uz_ve = compute_viscoelastic_correction(t_seconds, uz_elastic, eta=eta)
        ax.plot(t_years, uz_ve*1000, lw=2, label=f'η = {eta:.0e} Pa·s')
    ax.set_xlabel('Time (years)'); ax.set_ylabel('Vertical Displacement (mm)')
    ax.set_title('Viscoelastic Relaxation'); ax.legend(); ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    for eta in viscosities:
        uz_ve = compute_viscoelastic_correction(t_seconds, uz_elastic, eta=eta)
        ax2.plot(t_years, uz_ve/uz_elastic, lw=2, label=f'η = {eta:.0e} Pa·s')
    ax2.set_xlabel('Time (years)'); ax2.set_ylabel('Amplification Factor')
    ax2.set_title('Displacement Amplification Ratio'); ax2.legend(); ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/viscoelastic_correction.png', dpi=150, bbox_inches='tight')
    plt.close()

    ve_ratios = {}
    for eta in viscosities:
        uz_ve_10yr = compute_viscoelastic_correction(10*3.15e7, uz_elastic, eta=eta)
        ve_ratios[f'eta_{eta:.0e}'] = float(uz_ve_10yr / uz_elastic)
    results['viscoelastic_ratios_10yr'] = ve_ratios
    print(f"  Amplification at 10 yr: {ve_ratios}")

    # -------------------------------------------------------
    # Exp 6: Case Studies
    # -------------------------------------------------------
    print("\n[6] Case Study: Sakurajima & Aso...")
    sakurajima_params = {'xs': 0, 'ys': 0, 'd': 3000, 'dV': 2e6}
    aso_params = {'xs': 0, 'ys': 0, 'd': 6000, 'dV': 8e6}
    case_stations = generate_station_network(n_gnss=20, n_insar=625, n_grav=10,
                                              extent=15000)

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    for row, (vname, vparams) in enumerate([('Sakurajima', sakurajima_params),
                                              ('Aso', aso_params)]):
        case_data = generate_synthetic_data(vparams, case_stations, model='mogi')
        ci_x, ci_y = case_stations['insar']
        cg_x, cg_y = case_stations['gnss']
        cn = int(np.sqrt(len(ci_x)))

        los = case_data['insar']['los']
        im = axes[row, 0].imshow(los.reshape(cn, cn)*1000, extent=[-15,15,-15,15],
                                  cmap='RdBu_r', origin='lower')
        axes[row, 0].set_title(f'{vname}: InSAR LOS (mm)')
        plt.colorbar(im, ax=axes[row, 0])

        ux_c = case_data['gnss']['ux']; uy_c = case_data['gnss']['uy']
        uz_c = case_data['gnss']['uz']
        axes[row, 1].quiver(cg_x/1000, cg_y/1000, ux_c*1000, uy_c*1000,
                            scale=0.5, color='blue')
        sc = axes[row, 1].scatter(cg_x/1000, cg_y/1000, c=uz_c*1000,
                                   cmap='RdBu_r', s=80, edgecolor='k', zorder=5)
        axes[row, 1].set_title(f'{vname}: GNSS (horiz=arrows, vert=color mm)')
        plt.colorbar(sc, ax=axes[row, 1])
        axes[row, 1].set_xlim(-15, 15); axes[row, 1].set_ylim(-15, 15)

        dg_c = case_data['gravity']['dg']
        gr_x, gr_y = case_stations['gravity']
        sc2 = axes[row, 2].scatter(gr_x/1000, gr_y/1000, c=dg_c*1e8,
                                    cmap='viridis', s=120, edgecolor='k')
        axes[row, 2].set_title(f'{vname}: Gravity Change (µGal)')
        plt.colorbar(sc2, ax=axes[row, 2])
        axes[row, 2].set_xlim(-15, 15); axes[row, 2].set_ylim(-15, 15)

    for ax in axes.flat:
        ax.set_xlabel('X (km)'); ax.set_ylabel('Y (km)')
    plt.suptitle('Synthetic Case Studies: Sakurajima & Aso', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/case_study_data.png', dpi=150, bbox_inches='tight')
    plt.close()

    # MCMC for case studies
    for vname, vparams in [('Sakurajima', sakurajima_params),
                           ('Aso', aso_params)]:
        print(f"  MCMC for {vname}...")
        vdata = generate_synthetic_data(vparams, case_stations, model='mogi')
        vsamples, _, _ = run_mcmc_mogi(case_stations, vdata,
                                        n_walkers=32, n_steps=1500, n_burn=400)
        vpost = {p: float(np.mean(vsamples[:, i])) for i, p in enumerate(param_names)}
        vstd = {p: float(np.std(vsamples[:, i])) for i, p in enumerate(param_names)}
        results[f'{vname.lower()}_inversion'] = {
            'true': vparams, 'posterior_mean': vpost, 'posterior_std': vstd}
        print(f"    d={vpost['d']:.0f}±{vstd['d']:.0f} m, dV={vpost['dV']:.2e}±{vstd['dV']:.2e} m³")

    # -------------------------------------------------------
    # Summary Figure
    # -------------------------------------------------------
    print("\n[7] Summary figure...")
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    mn = list(results['model_comparison'].keys())
    muz = [results['model_comparison'][m]['max_uz_mm'] for m in mn]
    ax1.bar(mn, muz, color=['#2196F3', '#FF9800', '#4CAF50'])
    ax1.set_ylabel('Max Uz (mm)'); ax1.set_title('Forward Model Comparison')
    ax1.tick_params(axis='x', rotation=15)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(samples_mogi[:, 2], bins=40, density=True, alpha=0.7, color='steelblue')
    ax2.axvline(5000, color='r', ls='--', lw=2, label='True')
    ax2.set_xlabel('Depth (m)'); ax2.set_ylabel('Density')
    ax2.set_title('Mogi Depth Posterior'); ax2.legend()

    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(np.arange(n_epochs), kf_states[:, 3], 'b-', lw=2, label='KF')
    ax3.plot(np.arange(n_epochs), true_dV_history, 'r--', lw=1.5, label='True')
    ax3.set_xlabel('Epoch'); ax3.set_ylabel('dV (m³)')
    ax3.set_title('Volume Change Tracking'); ax3.legend()

    ax4 = fig.add_subplot(gs[1, 0])
    for eta in viscosities:
        uz_ve = compute_viscoelastic_correction(t_seconds, uz_elastic, eta=eta)
        ax4.plot(t_years, uz_ve/uz_elastic, lw=2, label=f'η={eta:.0e}')
    ax4.set_xlabel('Time (years)'); ax4.set_ylabel('Amplification')
    ax4.set_title('Viscoelastic Amplification'); ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    ax5 = fig.add_subplot(gs[1, 1:])
    x_pos = np.arange(4)
    width = 0.2
    for idx, (vn, color) in enumerate([('sakurajima', '#E53935'), ('aso', '#1E88E5')]):
        vt = results[f'{vn}_inversion']['true']
        ve = results[f'{vn}_inversion']['posterior_mean']
        vs = results[f'{vn}_inversion']['posterior_std']
        true_v = [vt['xs'], vt['ys'], vt['d'], vt['dV']/1e6]
        est_v = [ve['xs'], ve['ys'], ve['d'], ve['dV']/1e6]
        err_v = [vs['xs'], vs['ys'], vs['d'], vs['dV']/1e6]
        ax5.bar(x_pos - width + idx*width, true_v, width, alpha=0.5,
                color=color, label=f'{vn.title()} True')
        ax5.bar(x_pos + idx*width, est_v, width, yerr=err_v,
                alpha=0.8, color=color, label=f'{vn.title()} Est.', capsize=3)
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels(['xs (m)', 'ys (m)', 'd (m)', 'dV (×10⁶ m³)'])
    ax5.set_title('Case Study Parameters'); ax5.legend(fontsize=7, ncol=2)

    plt.suptitle('Volcanic Deformation Inversion: Summary', fontsize=16, fontweight='bold')
    plt.savefig(f'{FIGURES_DIR}/summary.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save results
    with open('results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print("\n" + "=" * 60)
    print("All experiments completed!")
    print(f"Figures saved to {FIGURES_DIR}/")
    print("=" * 60)
    return results


if __name__ == '__main__':
    results = run_full_experiment()
