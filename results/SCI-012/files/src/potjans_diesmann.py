"""
Potjans-Diesmann Cortical Microcircuit Model reimplementation.
Simplified version with 8 populations (4 layers × 2 types).
"""
import numpy as np


# Original Potjans-Diesmann (2014) parameters
LAYER_NAMES = ['L2/3', 'L4', 'L5', 'L6']
POP_NAMES = ['E', 'I']

# Full-scale neuron numbers per population
FULL_SCALE_N = {
    'L2/3E': 20683, 'L2/3I': 5834,
    'L4E': 21915, 'L4I': 5479,
    'L5E': 4850, 'L5I': 1065,
    'L6E': 14395, 'L6I': 2948,
}

# Connection probabilities (Table 5, Potjans & Diesmann 2014)
CONN_PROBS = np.array([
    [0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0.0, 0.0076, 0.0],
    [0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0.0, 0.0042, 0.0],
    [0.0077, 0.0059, 0.0497, 0.135, 0.0067, 0.0003, 0.0453, 0.0],
    [0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0.0, 0.1057, 0.0],
    [0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.0],
    [0.0548, 0.0269, 0.0257, 0.0022, 0.06, 0.3158, 0.0086, 0.0],
    [0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252],
    [0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443],
])

# Mean synaptic weights (pA)
W_EXC = 87.8  # excitatory
W_INH = -351.2  # inhibitory (4× excitatory)

# Background input rate
BG_RATE = 8.0  # Hz per synapse
N_EXT = 1600  # external synapses


class PotjansDiesmannCircuit:
    """Simplified Potjans-Diesmann cortical microcircuit."""
    
    def __init__(self, scale=0.1):
        self.scale = scale
        self.pop_names = [f"{l}{t}" for l in LAYER_NAMES for t in POP_NAMES]
        self.N_pops = {name: max(int(FULL_SCALE_N[name] * scale), 10) 
                       for name in self.pop_names}
        self.N_total = sum(self.N_pops.values())
        
        # Build population indices
        self.pop_slices = {}
        idx = 0
        for name in self.pop_names:
            n = self.N_pops[name]
            self.pop_slices[name] = slice(idx, idx + n)
            idx += n
        
        # Initialize Izhikevich neurons
        self.v = -65.0 * np.ones(self.N_total)
        self.u = -13.0 * np.ones(self.N_total)
        
        # Set excitatory/inhibitory parameters
        self.a = np.full(self.N_total, 0.02)
        self.b = np.full(self.N_total, 0.2)
        self.c = np.full(self.N_total, -65.0)
        self.d = np.full(self.N_total, 8.0)
        
        for name in self.pop_names:
            if 'I' in name:
                s = self.pop_slices[name]
                self.a[s] = 0.1
                self.d[s] = 2.0
        
        # Build sparse connectivity
        self._build_connectivity()
    
    def _build_connectivity(self):
        """Build weight matrix based on connection probabilities."""
        self.W = np.zeros((self.N_total, self.N_total))
        
        for i, pre_name in enumerate(self.pop_names):
            for j, post_name in enumerate(self.pop_names):
                p = CONN_PROBS[j, i]
                if p < 1e-6:
                    continue
                
                pre_s = self.pop_slices[pre_name]
                post_s = self.pop_slices[post_name]
                n_pre = self.N_pops[pre_name]
                n_post = self.N_pops[post_name]
                
                mask = np.random.rand(n_pre, n_post) < p
                w = W_EXC / 1000.0 if 'E' in pre_name else W_INH / 1000.0
                self.W[pre_s, post_s] = mask * w
    
    def step(self, dt, external_input=None):
        fired = self.v >= 30
        self.v[fired] = self.c[fired]
        self.u[fired] += self.d[fired]
        
        I_syn = self.W.T @ fired.astype(float)
        I_bg = np.random.randn(self.N_total) * 2.0 + 5.0
        I_ext = external_input if external_input is not None else np.zeros(self.N_total)
        I_total = I_syn + I_bg + I_ext
        
        self.v += dt * (0.04*self.v**2 + 5*self.v + 140 - self.u + I_total)
        self.u += dt * self.a * (self.b*self.v - self.u)
        
        return fired.astype(float)
    
    def simulate(self, T=1000, dt=0.5, stim_pop=None, stim_start=200, stim_end=400, stim_amp=10):
        steps = int(T / dt)
        spike_trains = {name: [] for name in self.pop_names}
        pop_rates = {name: np.zeros(steps) for name in self.pop_names}
        
        for t in range(steps):
            ext = np.zeros(self.N_total)
            if stim_pop and stim_start <= t*dt < stim_end:
                ext[self.pop_slices[stim_pop]] = stim_amp
            
            spikes = self.step(dt, ext)
            
            for name in self.pop_names:
                s = self.pop_slices[name]
                n_spikes = spikes[s].sum()
                pop_rates[name][t] = n_spikes / self.N_pops[name] * 1000/dt
                
                spike_idx = np.where(spikes[s] > 0)[0]
                for idx in spike_idx:
                    spike_trains[name].append((t*dt, idx))
        
        return {
            'pop_rates': pop_rates,
            'spike_trains': spike_trains,
            'T': T, 'dt': dt,
        }
