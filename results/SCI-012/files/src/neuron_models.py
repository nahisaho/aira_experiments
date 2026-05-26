"""
Neuron Models: Hodgkin-Huxley, Izhikevich, and AdEx implementations.
Provides comparative simulation and benchmarking.
"""
import numpy as np
import time


class HodgkinHuxley:
    """Hodgkin-Huxley neuron model with Na+, K+, and leak channels."""
    
    def __init__(self, N=1):
        self.N = N
        self.C_m = 1.0  # uF/cm^2
        self.g_Na = 120.0; self.g_K = 36.0; self.g_L = 0.3
        self.E_Na = 50.0; self.E_K = -77.0; self.E_L = -54.387
        self.V = -65.0 * np.ones(N)
        self.m = 0.05 * np.ones(N)
        self.h = 0.6 * np.ones(N)
        self.n = 0.32 * np.ones(N)
    
    def _alpha_m(self, V): return 0.1*(V+40)/(1 - np.exp(-(V+40)/10) + 1e-12)
    def _beta_m(self, V): return 4.0*np.exp(-(V+65)/18)
    def _alpha_h(self, V): return 0.07*np.exp(-(V+65)/20)
    def _beta_h(self, V): return 1.0/(1+np.exp(-(V+35)/10))
    def _alpha_n(self, V): return 0.01*(V+55)/(1-np.exp(-(V+55)/10)+1e-12)
    def _beta_n(self, V): return 0.125*np.exp(-(V+65)/80)
    
    def step(self, dt, I_ext):
        am, bm = self._alpha_m(self.V), self._beta_m(self.V)
        ah, bh = self._alpha_h(self.V), self._beta_h(self.V)
        an, bn = self._alpha_n(self.V), self._beta_n(self.V)
        
        self.m += dt * (am*(1-self.m) - bm*self.m)
        self.h += dt * (ah*(1-self.h) - bh*self.h)
        self.n += dt * (an*(1-self.n) - bn*self.n)
        
        I_Na = self.g_Na * self.m**3 * self.h * (self.V - self.E_Na)
        I_K = self.g_K * self.n**4 * (self.V - self.E_K)
        I_L = self.g_L * (self.V - self.E_L)
        
        self.V += dt/self.C_m * (I_ext - I_Na - I_K - I_L)
        spikes = (self.V > 0).astype(float)
        return self.V.copy(), spikes


class Izhikevich:
    """Izhikevich neuron model (2003)."""
    
    def __init__(self, N=1, mode='RS'):
        self.N = N
        params = {
            'RS': (0.02, 0.2, -65, 8),
            'IB': (0.02, 0.2, -55, 4),
            'CH': (0.02, 0.2, -50, 2),
            'FS': (0.1, 0.2, -65, 2),
            'LTS': (0.02, 0.25, -65, 2),
        }
        self.a, self.b, self.c, self.d = params.get(mode, params['RS'])
        self.v = self.c * np.ones(N)
        self.u = self.b * self.v
    
    def step(self, dt, I_ext):
        fired = self.v >= 30
        self.v[fired] = self.c
        self.u[fired] += self.d
        
        self.v += dt * (0.04*self.v**2 + 5*self.v + 140 - self.u + I_ext)
        self.u += dt * self.a * (self.b*self.v - self.u)
        
        spikes = fired.astype(float)
        return self.v.copy(), spikes


class AdEx:
    """Adaptive Exponential Integrate-and-Fire model."""
    
    def __init__(self, N=1):
        self.N = N
        self.C = 281.0  # pF
        self.gL = 30.0  # nS
        self.EL = -70.6  # mV
        self.VT = -50.4  # mV
        self.DeltaT = 2.0  # mV
        self.Vr = -70.6  # mV (reset)
        self.Vpeak = 20.0  # mV
        self.a = 4.0  # nS
        self.b = 80.5  # pA
        self.tau_w = 144.0  # ms
        
        self.V = self.EL * np.ones(N)
        self.w = 0.0 * np.ones(N)
    
    def step(self, dt, I_ext):
        exp_term = self.DeltaT * np.exp((self.V - self.VT) / self.DeltaT)
        dV = (-self.gL*(self.V - self.EL) + self.gL*exp_term - self.w + I_ext) / self.C
        dw = (self.a*(self.V - self.EL) - self.w) / self.tau_w
        
        self.V += dt * dV
        self.w += dt * dw
        
        fired = self.V >= self.Vpeak
        spikes = fired.astype(float)
        self.V[fired] = self.Vr
        self.w[fired] += self.b
        
        return self.V.copy(), spikes


def benchmark_models(N=1000, T=1000, dt=0.1):
    """Benchmark all three neuron models."""
    steps = int(T / dt)
    results = {}
    
    for name, ModelClass in [('HH', HodgkinHuxley), ('Izhikevich', Izhikevich), ('AdEx', AdEx)]:
        model = ModelClass(N=N)
        spike_counts = np.zeros(N)
        V_trace = []
        
        t0 = time.time()
        for t in range(steps):
            I_ext = np.random.randn(N) * 5 + 10
            V, spikes = model.step(dt, I_ext)
            spike_counts += spikes
            if t % 100 == 0:
                V_trace.append(V[0])
        
        elapsed = time.time() - t0
        results[name] = {
            'time': elapsed,
            'mean_rate': np.mean(spike_counts) / (T/1000),
            'std_rate': np.std(spike_counts) / (T/1000),
            'V_trace': np.array(V_trace),
            'spike_counts': spike_counts,
        }
    
    return results
