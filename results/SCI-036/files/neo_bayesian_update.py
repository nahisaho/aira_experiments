"""
Bayesian Collision Probability Update Module
Sequential Bayesian inference for NEO impact probability using astrometric observations.
"""

import numpy as np
import pymc as pm
import arviz as az
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from dataclasses import dataclass, field
from typing import Optional
import warnings
warnings.filterwarnings('ignore')


@dataclass
class AstrometricObservation:
    """Single astrometric observation of a NEO."""
    epoch: float            # Julian Date
    ra_deg: float           # Right Ascension [deg]
    dec_deg: float          # Declination [deg]
    sigma_arcsec: float     # 1-sigma positional uncertainty [arcsec]
    observatory_code: str = "500"   # MPC observatory code


@dataclass
class PriorBelief:
    """Prior distribution over NEO orbital elements."""
    a_mean: float; a_std: float     # Semi-major axis [AU]
    e_mean: float; e_std: float     # Eccentricity
    impact_prob_prior: float = 1e-4  # Initial impact probability


class BayesianImpactUpdater:
    """
    Sequential Bayesian updating of NEO impact probability.

    Model: P(impact | obs) ∝ P(obs | orbital params) × P(orbital params)

    Uses a combination of:
    - Analytic Gaussian update for orbital elements (Kalman-like)
    - Full MCMC for posterior over impact indicator
    """

    def __init__(self, prior: PriorBelief, seed: int = 42):
        self.prior = prior
        self.seed = seed
        self.observations: list[AstrometricObservation] = []
        self.posterior_history: list[dict] = []
        self.rng = np.random.default_rng(seed)

        # Current posterior state (Gaussian approximation)
        self.post_a_mean = prior.a_mean
        self.post_a_std = prior.a_std
        self.post_e_mean = prior.e_mean
        self.post_e_std = prior.e_std
        self.post_impact_prob = prior.impact_prob_prior

    def _likelihood_ratio(self, obs: AstrometricObservation,
                           a_trial: float, e_trial: float) -> float:
        """
        Simplified observation likelihood: how well does a trial orbit
        predict the observed RA/Dec?
        In full implementation this requires numerical orbit fitting.
        Here we model the residual as Gaussian.
        """
        # Simulate predicted RA/Dec from trial orbital elements
        # (simplified: use perturbation from reference orbit)
        da = a_trial - self.prior.a_mean
        de = e_trial - self.prior.e_mean
        # Orbital position sensitivity (rough, in arcsec / AU)
        dra_pred = da * 500.0 + de * 200.0  # arcsec
        ddec_pred = da * 300.0 - de * 150.0
        sigma = obs.sigma_arcsec
        chi2 = (dra_pred**2 + ddec_pred**2) / sigma**2
        return np.exp(-0.5 * chi2)

    def kalman_update(self, obs: AstrometricObservation) -> dict:
        """
        Kalman-filter style sequential update of orbital elements posterior.
        Bayesian update: posterior ∝ likelihood × prior (Gaussian × Gaussian = Gaussian).
        """
        # Observation variance from astrometry
        sigma_obs = obs.sigma_arcsec / 500.0   # crude scaling to AU units

        # Update a (semi-major axis)
        K_a = self.post_a_std**2 / (self.post_a_std**2 + sigma_obs**2)
        innovation_a = self.rng.normal(0, obs.sigma_arcsec / 1000)   # simulated residual
        self.post_a_mean += K_a * innovation_a
        self.post_a_std = np.sqrt((1 - K_a) * self.post_a_std**2)

        # Update e (eccentricity)
        K_e = self.post_e_std**2 / (self.post_e_std**2 + (sigma_obs * 0.3)**2)
        innovation_e = self.rng.normal(0, obs.sigma_arcsec / 5000)
        self.post_e_mean = np.clip(self.post_e_mean + K_e * innovation_e, 0, 0.99)
        self.post_e_std = np.sqrt((1 - K_e) * self.post_e_std**2)

        # Update impact probability via Bayes rule
        # P(impact | obs) = P(obs | impact) × P(impact) / P(obs)
        # P(obs | impact) / P(obs | no impact) encodes whether obs moves
        # the orbit solution closer to or further from impact corridor
        # Simplified: if uncertainty decreases, so can impact probability
        sigma_shrink_ratio = self.post_a_std / self.prior.a_std
        # Impact probability scales roughly with orbital uncertainty volume
        self.post_impact_prob *= sigma_shrink_ratio**2
        self.post_impact_prob = np.clip(self.post_impact_prob, 1e-12, 1.0)

        state = {
            'epoch': obs.epoch,
            'a_mean': self.post_a_mean, 'a_std': self.post_a_std,
            'e_mean': self.post_e_mean, 'e_std': self.post_e_std,
            'impact_prob': self.post_impact_prob,
        }
        self.observations.append(obs)
        self.posterior_history.append(state)
        return state

    def full_mcmc_posterior(self, n_obs: int = 20) -> az.InferenceData:
        """
        Full PyMC MCMC posterior over orbital elements and impact flag.
        """
        print(f"[MCMC] Running full Bayesian posterior with {n_obs} observations")
        # Generate synthetic observations for demonstration
        true_a = self.prior.a_mean + self.rng.normal(0, self.prior.a_std * 0.3)
        true_e = self.prior.e_mean + self.rng.normal(0, self.prior.e_std * 0.3)
        epochs = np.linspace(2459000, 2459000 + 365 * 2, n_obs)
        obs_ra = true_a * 50 + true_e * 30 + self.rng.normal(0, 0.5, n_obs)
        obs_dec = true_a * 20 - true_e * 15 + self.rng.normal(0, 0.5, n_obs)
        obs_sigma = np.full(n_obs, 0.5)   # arcsec

        with pm.Model() as model:
            # Priors on orbital elements
            a = pm.Normal('a', mu=self.prior.a_mean, sigma=self.prior.a_std)
            e = pm.TruncatedNormal('e', mu=self.prior.e_mean, sigma=self.prior.e_std,
                                    lower=0.0, upper=0.99)

            # Forward model (linear approximation)
            ra_pred = a * 50 + e * 30
            dec_pred = a * 20 - e * 15

            # Likelihood
            pm.Normal('ra_obs', mu=ra_pred, sigma=obs_sigma, observed=obs_ra)
            pm.Normal('dec_obs', mu=dec_pred, sigma=obs_sigma, observed=obs_dec)

            # Derived impact probability: logistic function of orbital deviation
            a_dev = pm.math.abs_(a - 1.0) + pm.math.abs_(e - 0.1)
            log_impact = -5.0 - 10.0 * a_dev
            p_impact = pm.Deterministic('p_impact', pm.math.invlogit(log_impact))

            # Sample
            idata = pm.sample(500, tune=300, chains=2, progressbar=False,
                               random_seed=self.seed,
                               target_accept=0.9, return_inferencedata=True)

        print(f"[MCMC] Sampling complete")
        return idata

    def simulate_observation_campaign(self, n_obs: int = 30,
                                       sigma_arcsec: float = 0.3,
                                       base_epoch: float = 2459000.0) -> list[dict]:
        """
        Simulate a sequential observation campaign and track Bayesian updates.
        """
        print(f"[Bayesian Update] Simulating {n_obs}-observation campaign")
        history = []
        for i in range(n_obs):
            epoch = base_epoch + i * 30.0  # 1 obs per month
            obs = AstrometricObservation(
                epoch=epoch,
                ra_deg=self.rng.normal(180.0, sigma_arcsec / 3600),
                dec_deg=self.rng.normal(15.0, sigma_arcsec / 3600),
                sigma_arcsec=sigma_arcsec * self.rng.uniform(0.8, 1.2)
            )
            state = self.kalman_update(obs)
            history.append({
                'obs_number': i + 1,
                'epoch': epoch,
                'impact_prob': state['impact_prob'],
                'a_std': state['a_std'],
                'e_std': state['e_std'],
            })
        return history

    def plot_bayesian_update(self, history: list[dict], idata: Optional[az.InferenceData],
                              save_path: str) -> None:
        """Visualize Bayesian probability evolution and MCMC posterior."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        obs_nums = [h['obs_number'] for h in history]
        probs = [h['impact_prob'] for h in history]
        a_stds = [h['a_std'] for h in history]
        e_stds = [h['e_std'] for h in history]

        # Impact probability evolution
        axes[0, 0].semilogy(obs_nums, probs, 'o-', color='crimson', markersize=5)
        axes[0, 0].axhline(self.prior.impact_prob_prior, color='gray', linestyle='--',
                            label=f'Prior: {self.prior.impact_prob_prior:.1e}')
        axes[0, 0].set_xlabel('Observation Number')
        axes[0, 0].set_ylabel('Impact Probability')
        axes[0, 0].set_title('Bayesian Update: Impact Probability vs. Observations')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        # Uncertainty reduction
        axes[0, 1].plot(obs_nums, np.array(a_stds) / self.prior.a_std,
                         '-', color='steelblue', label='a uncertainty')
        axes[0, 1].plot(obs_nums, np.array(e_stds) / self.prior.e_std,
                         '-', color='seagreen', label='e uncertainty')
        axes[0, 1].set_xlabel('Observation Number')
        axes[0, 1].set_ylabel('Uncertainty / Prior σ')
        axes[0, 1].set_title('Orbital Uncertainty Reduction')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        # MCMC posterior
        if idata is not None:
            try:
                a_post = idata.posterior['a'].values.flatten()
                e_post = idata.posterior['e'].values.flatten()
                axes[1, 0].hist2d(a_post, e_post, bins=40, cmap='Blues')
                axes[1, 0].axvline(self.prior.a_mean, color='red', linestyle='--', label='Prior mean')
                axes[1, 0].set_xlabel('Semi-major axis [AU]')
                axes[1, 0].set_ylabel('Eccentricity')
                axes[1, 0].set_title('MCMC Posterior: Orbital Elements')
                axes[1, 0].legend()
            except Exception:
                axes[1, 0].text(0.5, 0.5, 'MCMC posterior\nnot available',
                                 ha='center', va='center', transform=axes[1, 0].transAxes)

            try:
                p_imp = idata.posterior['p_impact'].values.flatten()
                axes[1, 1].hist(np.log10(p_imp + 1e-15), bins=40, color='tomato',
                                 edgecolor='white', density=True)
                axes[1, 1].set_xlabel('log₁₀(Impact Probability)')
                axes[1, 1].set_ylabel('Density')
                axes[1, 1].set_title('MCMC Posterior: Impact Probability Distribution')
            except Exception:
                axes[1, 1].text(0.5, 0.5, 'Impact probability\nposterior not available',
                                 ha='center', va='center', transform=axes[1, 1].transAxes)
        else:
            for ax in [axes[1, 0], axes[1, 1]]:
                ax.text(0.5, 0.5, 'MCMC skipped\n(fast mode)', ha='center', va='center',
                         transform=ax.transAxes, fontsize=12)

        plt.suptitle('Bayesian Collision Probability Assessment', fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"[Plot] Saved: {save_path}")
