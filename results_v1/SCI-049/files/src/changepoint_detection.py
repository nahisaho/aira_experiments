"""
Module 1: Change Point Detection (PELT / BOCPD)
- PELT (Pruned Exact Linear Time) for offline batch detection
- BOCPD (Bayesian Online Changepoint Detection) for streaming
"""
import numpy as np
from scipy import stats


class PELTDetector:
    """Pruned Exact Linear Time change point detection using ruptures-compatible logic."""

    def __init__(self, model="rbf", penalty=None, min_size=2):
        self.model = model
        self.penalty = penalty
        self.min_size = min_size
        self.breakpoints_ = []

    def fit_predict(self, signal: np.ndarray, pen: float = None) -> dict:
        import ruptures
        penalty = pen or self.penalty or np.log(len(signal)) * np.var(signal)
        algo = ruptures.Pelt(model=self.model, min_size=self.min_size).fit(signal)
        self.breakpoints_ = algo.predict(pen=penalty)
        segments = self._compute_segments(signal)
        return {
            "breakpoints": self.breakpoints_,
            "n_changepoints": len(self.breakpoints_) - 1,
            "segments": segments,
            "penalty_used": penalty,
        }

    def _compute_segments(self, signal):
        segs = []
        prev = 0
        for bp in self.breakpoints_:
            seg = signal[prev:bp]
            segs.append({"start": prev, "end": bp, "mean": float(np.mean(seg)),
                         "std": float(np.std(seg)), "length": len(seg)})
            prev = bp
        return segs


class BOCPDDetector:
    """Bayesian Online Changepoint Detection for streaming data."""

    def __init__(self, hazard_rate=1/250, mu0=0, kappa0=1, alpha0=1, beta0=1):
        self.hazard = 1.0 / hazard_rate if hazard_rate < 1 else hazard_rate
        self.mu0 = mu0
        self.kappa0 = kappa0
        self.alpha0 = alpha0
        self.beta0 = beta0
        self.run_length_probs = []
        self.changepoint_probs = []

    def detect(self, data: np.ndarray, threshold: float = 0.5) -> dict:
        T = len(data)
        R = np.zeros((T + 1, T + 1))
        R[0, 0] = 1.0

        mu_params = np.array([self.mu0])
        kappa_params = np.array([self.kappa0])
        alpha_params = np.array([self.alpha0])
        beta_params = np.array([self.beta0])

        changepoints = []
        cp_probs = np.zeros(T)

        for t in range(T):
            x = data[t]
            pred_probs = self._student_t_pdf(
                x, mu_params, kappa_params, alpha_params, beta_params
            )

            H = 1.0 / self.hazard
            R[1:t+2, t+1] = R[:t+1, t] * pred_probs * (1 - H)
            R[0, t+1] = np.sum(R[:t+1, t] * pred_probs * H)

            evidence = np.sum(R[:t+2, t+1])
            if evidence > 0:
                R[:t+2, t+1] /= evidence

            cp_probs[t] = R[0, t+1]
            if cp_probs[t] > threshold:
                changepoints.append(t)

            # Update sufficient statistics
            new_mu = np.append([self.mu0],
                               (kappa_params * mu_params + x) / (kappa_params + 1))
            new_kappa = np.append([self.kappa0], kappa_params + 1)
            new_alpha = np.append([self.alpha0], alpha_params + 0.5)
            new_beta = np.append([self.beta0],
                                 beta_params + 0.5 * kappa_params * (x - mu_params)**2 / (kappa_params + 1))
            mu_params = new_mu
            kappa_params = new_kappa
            alpha_params = new_alpha
            beta_params = new_beta

        self.changepoint_probs = cp_probs
        return {
            "changepoints": changepoints,
            "n_changepoints": len(changepoints),
            "changepoint_probabilities": cp_probs,
            "threshold": threshold,
        }

    @staticmethod
    def _student_t_pdf(x, mu, kappa, alpha, beta):
        df = 2 * alpha
        scale = np.sqrt(beta * (kappa + 1) / (alpha * kappa))
        scale = np.maximum(scale, 1e-10)
        return stats.t.pdf(x, df=df, loc=mu, scale=scale)
