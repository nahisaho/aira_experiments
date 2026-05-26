"""
Synaptic Plasticity: STDP and Homeostatic Plasticity implementations.
"""
import numpy as np


class STDP:
    """Spike-Timing-Dependent Plasticity with exponential kernel."""
    
    def __init__(self, N_pre, N_post, A_plus=0.01, A_minus=0.012,
                 tau_plus=20.0, tau_minus=20.0, w_max=1.0, w_min=0.0):
        self.N_pre = N_pre
        self.N_post = N_post
        self.A_plus = A_plus
        self.A_minus = A_minus
        self.tau_plus = tau_plus
        self.tau_minus = tau_minus
        self.w_max = w_max
        self.w_min = w_min
        self.W = np.random.uniform(0.1, 0.5, (N_pre, N_post))
        self.trace_pre = np.zeros(N_pre)
        self.trace_post = np.zeros(N_post)
    
    def update(self, dt, pre_spikes, post_spikes):
        self.trace_pre *= np.exp(-dt / self.tau_plus)
        self.trace_post *= np.exp(-dt / self.tau_minus)
        
        pre_idx = np.where(pre_spikes > 0)[0]
        post_idx = np.where(post_spikes > 0)[0]
        
        # LTP: post fires after pre
        if len(post_idx) > 0:
            for j in post_idx:
                self.W[:, j] += self.A_plus * self.trace_pre
        
        # LTD: pre fires after post
        if len(pre_idx) > 0:
            for i in pre_idx:
                self.W[i, :] -= self.A_minus * self.trace_post
        
        self.trace_pre[pre_idx] += 1.0
        self.trace_post[post_idx] += 1.0
        
        np.clip(self.W, self.w_min, self.w_max, out=self.W)
        return self.W.copy()


class HomeostaticPlasticity:
    """Homeostatic synaptic scaling to maintain target firing rate."""
    
    def __init__(self, N, target_rate=5.0, tau_homeo=10000.0, eta=0.001):
        self.N = N
        self.target_rate = target_rate  # Hz
        self.tau_homeo = tau_homeo  # ms
        self.eta = eta
        self.rate_estimate = target_rate * np.ones(N)
        self.scaling_factor = np.ones(N)
    
    def update(self, dt, spike_counts, weights):
        self.rate_estimate += dt/self.tau_homeo * (
            spike_counts * 1000.0/dt - self.rate_estimate
        )
        error = self.target_rate - self.rate_estimate
        self.scaling_factor += self.eta * error * dt/self.tau_homeo
        self.scaling_factor = np.clip(self.scaling_factor, 0.1, 10.0)
        
        scaled_weights = weights * self.scaling_factor[np.newaxis, :]
        return scaled_weights, self.rate_estimate.copy()


def simulate_stdp_learning(N_pre=100, N_post=50, T=5000, dt=0.5):
    """Simulate STDP learning with pre/post spike trains."""
    stdp = STDP(N_pre, N_post)
    homeo = HomeostaticPlasticity(N_post)
    
    steps = int(T / dt)
    weight_history = []
    rate_history = []
    
    for t in range(steps):
        pre_rate = 10 + 5*np.sin(2*np.pi*t*dt/1000)
        pre_spikes = (np.random.rand(N_pre) < pre_rate*dt/1000).astype(float)
        
        I_syn = stdp.W.T @ pre_spikes
        post_spikes = (np.random.rand(N_post) < (I_syn/N_pre*50)*dt/1000).astype(float)
        
        stdp.update(dt, pre_spikes, post_spikes)
        scaled_W, rates = homeo.update(dt, post_spikes, stdp.W)
        
        if t % 200 == 0:
            weight_history.append(stdp.W.mean())
            rate_history.append(rates.mean())
    
    return {
        'weight_history': np.array(weight_history),
        'rate_history': np.array(rate_history),
        'final_weights': stdp.W.copy(),
    }
