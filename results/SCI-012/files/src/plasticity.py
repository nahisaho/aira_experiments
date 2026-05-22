"""
Synaptic plasticity rules:
  - Spike-Timing-Dependent Plasticity (STDP)
  - Homeostatic (synaptic scaling) plasticity
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# STDP
# ---------------------------------------------------------------------------

@dataclass
class STDPParams:
    A_plus:   float = 0.005    # LTP amplitude
    A_minus:  float = 0.005    # LTD amplitude
    tau_plus: float = 20.0     # pre-before-post time constant [ms]
    tau_minus:float = 20.0     # post-before-pre time constant [ms]
    w_min:    float = 0.0      # minimum synaptic weight
    w_max:    float = 1.0      # maximum synaptic weight
    # Nearest-neighbour STDP (Song et al. 2000 variant)


class STDPSynapse:
    """
    Pair-wise STDP rule with eligibility traces.

    Maintains per-synapse weight and per-neuron eligibility traces.
    """

    def __init__(self, n_pre: int, n_post: int,
                 connectivity: Optional[np.ndarray] = None,
                 params: Optional[STDPParams] = None,
                 rng: Optional[np.random.Generator] = None):
        self.p = params or STDPParams()
        self.rng = rng or np.random.default_rng(42)

        self.n_pre  = n_pre
        self.n_post = n_post

        if connectivity is None:
            # sparse random: each post receives from ~10% of pre
            self.W = (self.rng.random((n_post, n_pre)) < 0.1).astype(float)
            self.W *= self.rng.uniform(0.1, 0.5, (n_post, n_pre))
        else:
            self.W = connectivity.copy().astype(float)

        # Eligibility traces
        self.x_pre  = np.zeros(n_pre)   # pre-synaptic trace
        self.x_post = np.zeros(n_post)  # post-synaptic trace

    def update(self, pre_spikes: np.ndarray, post_spikes: np.ndarray,
               dt: float) -> None:
        """
        Update weights given spike vectors for current time step.

        Args:
            pre_spikes:  bool array (n_pre,), True if pre neuron fired
            post_spikes: bool array (n_post,), True if post neuron fired
            dt: time step [ms]
        """
        p = self.p

        # Decay traces
        self.x_pre  *= np.exp(-dt / p.tau_plus)
        self.x_post *= np.exp(-dt / p.tau_minus)

        # LTD: post fires, depress synapses from all pre that had recent activity
        if post_spikes.any():
            # shape: (n_post_fired, n_pre)
            fired_post = np.where(post_spikes)[0]
            self.W[np.ix_(fired_post, np.arange(self.n_pre))] -= (
                p.A_minus * self.x_pre[np.newaxis, :]
            )

        # Update pre traces THEN do LTP (post fires first in time already handled)
        if pre_spikes.any():
            self.x_pre[pre_spikes] += 1.0

        # LTP: pre fires, potentiate synapses to all post that had recent activity
        if pre_spikes.any():
            fired_pre = np.where(pre_spikes)[0]
            self.W[np.ix_(np.arange(self.n_post), fired_pre)] += (
                p.A_plus * self.x_post[:, np.newaxis]
            )

        # Update post traces
        if post_spikes.any():
            self.x_post[post_spikes] += 1.0

        # Clip weights
        np.clip(self.W, p.w_min, p.w_max, out=self.W)

    def get_input_current(self, pre_spikes: np.ndarray) -> np.ndarray:
        """
        Compute synaptic input to each post neuron from pre spike vector.

        Returns shape (n_post,) in units of weight.
        """
        return self.W @ pre_spikes.astype(float)

    def weight_stats(self) -> dict:
        mask = self.W > 0
        return {
            "mean":   float(self.W[mask].mean()) if mask.any() else 0.0,
            "std":    float(self.W[mask].std())  if mask.any() else 0.0,
            "min":    float(self.W[mask].min())  if mask.any() else 0.0,
            "max":    float(self.W[mask].max())  if mask.any() else 0.0,
            "n_syn":  int(mask.sum()),
        }


# ---------------------------------------------------------------------------
# Homeostatic (synaptic scaling) plasticity
# ---------------------------------------------------------------------------

@dataclass
class HomeostaticParams:
    target_rate:  float = 5.0     # target firing rate [Hz]
    tau_scaling:  float = 1e4     # scaling time constant [ms]  (slow)
    eta:          float = 0.001   # scaling learning rate
    w_min:        float = 0.0
    w_max:        float = 2.0
    # Turrigiano et al. 1998 model


class HomeostaticScaling:
    """
    Multiplicative synaptic scaling that nudges population firing rate
    toward target_rate.
    """

    def __init__(self, n_neurons: int, params: Optional[HomeostaticParams] = None):
        self.p = params or HomeostaticParams()
        self.n = n_neurons
        self.rate_estimate = np.zeros(n_neurons)   # running firing rate estimate [Hz]
        self._tau_rate = 1000.0                    # rate estimation time constant [ms]

    def update(self, spikes: np.ndarray, W: np.ndarray, dt: float) -> np.ndarray:
        """
        Update weight matrix W in-place using homeostatic scaling.

        Args:
            spikes: bool array (n_neurons,)
            W:      weight matrix (n_post, n_pre) – will be scaled in-place
            dt:     time step [ms]

        Returns:
            Updated W
        """
        p = self.p

        # Update firing rate estimate (exponential moving average)
        alpha = dt / self._tau_rate
        self.rate_estimate += alpha * (spikes.astype(float) / (dt * 1e-3) - self.rate_estimate)

        # Scale factor for each post neuron
        error = p.target_rate - self.rate_estimate          # shape (n,)
        scale = 1.0 + p.eta * (dt / p.tau_scaling) * error # per neuron

        # Multiplicative scaling: rows of W correspond to post neurons
        W *= scale[:, np.newaxis]
        np.clip(W, p.w_min, p.w_max, out=W)
        return W

    def rate_stats(self) -> dict:
        return {
            "mean_rate": float(self.rate_estimate.mean()),
            "std_rate":  float(self.rate_estimate.std()),
            "max_rate":  float(self.rate_estimate.max()),
        }


# ---------------------------------------------------------------------------
# Combined STDP + Homeostatic simulation helper
# ---------------------------------------------------------------------------

def run_plasticity_demo(n_pre: int = 100, n_post: int = 100,
                        T_ms: float = 5000.0, dt: float = 1.0,
                        input_rate: float = 10.0,
                        rng: Optional[np.random.Generator] = None) -> dict:
    """
    Simulate Poisson pre-synaptic population driving post-synaptic neurons
    through STDP + homeostatic synapses.

    Returns weight evolution and firing rate history.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    stdp   = STDPSynapse(n_pre, n_post, rng=rng)
    homeo  = HomeostaticScaling(n_post)
    p_fire = input_rate * dt * 1e-3   # spike probability per step

    T = int(T_ms / dt)
    weight_history  = []  # sampled every 500 steps
    rate_history    = []
    post_V = np.full(n_post, -65.0)
    post_spikes_prev = np.zeros(n_post, dtype=bool)

    for t in range(T):
        pre_spikes  = rng.random(n_pre)  < p_fire
        # Simple LIF-like post: integrate and fire
        I_syn = stdp.get_input_current(pre_spikes) * 5.0  # scale
        post_V += (-0.1*(post_V + 65.0) + I_syn) * dt / 10.0
        post_spikes = post_V >= -55.0
        post_V[post_spikes] = -65.0  # reset

        stdp.update(pre_spikes, post_spikes, dt)
        homeo.update(post_spikes, stdp.W, dt)

        if t % 500 == 0:
            ws = stdp.weight_stats()
            rs = homeo.rate_stats()
            weight_history.append({"t_ms": t*dt, **ws})
            rate_history.append({"t_ms": t*dt, **rs})

    return {
        "weight_history": weight_history,
        "rate_history":   rate_history,
        "final_W":        stdp.W.copy(),
        "final_rates":    homeo.rate_estimate.copy(),
    }
