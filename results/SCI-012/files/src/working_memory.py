"""
Working Memory Task: Delayed match-to-sample using SNN.
Models persistent activity and compares with experimental data.
"""
import numpy as np
from scipy import sparse


class WorkingMemoryNetwork:
    """SNN model for delayed match-to-sample working memory task."""
    
    def __init__(self, N_exc=800, N_inh=200, N_selective=4, f=0.15):
        self.N_exc = N_exc
        self.N_inh = N_inh
        self.N = N_exc + N_inh
        self.N_selective = N_selective
        self.f = f  # fraction of E neurons in each selective pool
        self.pool_size = int(N_exc * f)
        self.non_selective = N_exc - N_selective * self.pool_size
        
        # Initialize Izhikevich neurons
        self.v = -65.0 * np.ones(self.N)
        self.u = -13.0 * np.ones(self.N)
        
        # Excitatory params
        self.a = np.full(self.N, 0.02)
        self.b = np.full(self.N, 0.2)
        self.c = np.full(self.N, -65.0)
        self.d = np.full(self.N, 8.0)
        
        # Inhibitory params
        self.a[N_exc:] = 0.1
        self.d[N_exc:] = 2.0
        
        self._build_connectivity()
        # Convert to sparse for fast matmul
        self.W_sparse = sparse.csr_matrix(self.W)
    
    def _build_connectivity(self):
        """Build structured connectivity for working memory."""
        self.W = np.zeros((self.N, self.N))
        
        w_plus = 1.7  # potentiated weight within selective pools
        w_minus = 1.0 - self.f*(w_plus - 1.0)/(1.0 - self.f)
        
        # E→E connections
        for i in range(self.N_selective):
            pool_start = i * self.pool_size
            pool_end = pool_start + self.pool_size
            
            # Within-pool (strong)
            mask = np.random.rand(self.pool_size, self.pool_size) < 0.2
            self.W[pool_start:pool_end, pool_start:pool_end] = mask * w_plus * 0.05
            
            # Between pools (weak)
            for j in range(self.N_selective):
                if i != j:
                    other_start = j * self.pool_size
                    other_end = other_start + self.pool_size
                    mask = np.random.rand(self.pool_size, self.pool_size) < 0.2
                    self.W[pool_start:pool_end, other_start:other_end] = mask * w_minus * 0.05
        
        # E→I and I→E connections
        mask_ei = np.random.rand(self.N_exc, self.N_inh) < 0.2
        self.W[:self.N_exc, self.N_exc:] = mask_ei.astype(float) * 0.04
        
        mask_ie = np.random.rand(self.N_inh, self.N_exc) < 0.2
        self.W[self.N_exc:, :self.N_exc] = -mask_ie.astype(float) * 0.15
        
        # I→I connections
        mask_ii = np.random.rand(self.N_inh, self.N_inh) < 0.2
        self.W[self.N_exc:, self.N_exc:] = -mask_ii.astype(float) * 0.10
    
    def run_trial(self, stim_pool=0, T=2000, dt=0.5,
                  stim_start=300, stim_end=500, stim_amp=15,
                  probe_start=1500, probe_end=1700):
        """Run a single delayed match-to-sample trial."""
        steps = int(T / dt)
        
        pool_start = stim_pool * self.pool_size
        pool_end = pool_start + self.pool_size
        
        # Reset
        self.v = -65.0 * np.ones(self.N)
        self.u = -13.0 * np.ones(self.N)
        
        # Recording
        pool_rates = np.zeros((self.N_selective, steps))
        inh_rates = np.zeros(steps)
        
        for t in range(steps):
            time_ms = t * dt
            
            fired = self.v >= 30
            self.v[fired] = self.c[fired]
            self.u[fired] += self.d[fired]
            
            I_syn = self.W_sparse.T @ fired.astype(float) * 100
            I_bg = np.random.randn(self.N) * 3.0 + 5.0
            I_stim = np.zeros(self.N)
            
            # Sample stimulus
            if stim_start <= time_ms < stim_end:
                I_stim[pool_start:pool_end] = stim_amp
            
            # Probe stimulus
            if probe_start <= time_ms < probe_end:
                I_stim[pool_start:pool_end] = stim_amp * 0.5
            
            I_total = I_syn + I_bg + I_stim
            
            # Euler integration with voltage clamping for stability
            dv = 0.04*self.v**2 + 5*self.v + 140 - self.u + I_total
            self.v += dt * np.clip(dv, -500, 500)
            self.v = np.clip(self.v, -100, 40)
            self.u += dt * self.a * (self.b*self.v - self.u)
            
            for p in range(self.N_selective):
                ps = p * self.pool_size
                pe = ps + self.pool_size
                pool_rates[p, t] = fired[ps:pe].sum() / self.pool_size * 1000/dt
            
            inh_rates[t] = fired[self.N_exc:].sum() / self.N_inh * 1000/dt
        
        times = np.arange(steps) * dt
        return {
            'times': times,
            'pool_rates': pool_rates,
            'inh_rates': inh_rates,
            'stim_pool': stim_pool,
            'T': T,
        }


def run_working_memory_experiment(n_trials=5):
    """Run multiple trials of working memory task."""
    net = WorkingMemoryNetwork(N_exc=400, N_inh=100, N_selective=4, f=0.15)
    
    results = []
    for trial in range(n_trials):
        stim_pool = trial % net.N_selective
        result = net.run_trial(stim_pool=stim_pool)
        results.append(result)
    
    # Compute metrics
    delay_rates = []
    baseline_rates = []
    stim_rates = []
    
    for r in results:
        sp = r['stim_pool']
        dt = r['times'][1] - r['times'][0]
        
        baseline_mask = r['times'] < 300
        stim_mask = (r['times'] >= 300) & (r['times'] < 500)
        delay_mask = (r['times'] >= 600) & (r['times'] < 1500)
        
        baseline_rates.append(r['pool_rates'][sp, baseline_mask].mean())
        stim_rates.append(r['pool_rates'][sp, stim_mask].mean())
        delay_rates.append(r['pool_rates'][sp, delay_mask].mean())
    
    return {
        'trials': results,
        'baseline_rate': np.mean(baseline_rates),
        'stim_rate': np.mean(stim_rates),
        'delay_rate': np.mean(delay_rates),
        'baseline_std': np.std(baseline_rates),
        'stim_std': np.std(stim_rates),
        'delay_std': np.std(delay_rates),
    }
