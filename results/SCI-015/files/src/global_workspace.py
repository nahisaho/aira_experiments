"""
Global Workspace Theory (GWT) Computational Model

Implements a simplified Global Workspace architecture for simulating
conscious access and information broadcasting.
"""

import numpy as np
import networkx as nx


class GlobalWorkspaceModel:
    """
    Computational model of Global Neuronal Workspace Theory.
    
    Architecture:
    - Multiple specialized processors (sensory, motor, memory, etc.)
    - A central workspace that broadcasts information globally
    - Competition and ignition dynamics
    """
    
    def __init__(self, n_processors=6, workspace_size=10, processor_size=20):
        self.n_processors = n_processors
        self.workspace_size = workspace_size
        self.processor_size = processor_size
        
        self.workspace = np.zeros(workspace_size)
        self.processors = [np.zeros(processor_size) for _ in range(n_processors)]
        
        self.W_to_workspace = [np.random.randn(workspace_size, processor_size) * 0.1 
                                for _ in range(n_processors)]
        self.W_from_workspace = [np.random.randn(processor_size, workspace_size) * 0.1 
                                  for _ in range(n_processors)]
        
        self.ignition_threshold = 0.5
        self.decay_rate = 0.1
        self.broadcast_gain = 1.5
        
    def sigmoid(self, x, gain=1.0):
        return 1 / (1 + np.exp(-gain * x))
    
    def step(self, inputs, noise_level=0.05):
        """Single timestep of GWT dynamics."""
        # Bottom-up: processors send to workspace
        workspace_input = np.zeros(self.workspace_size)
        for i, proc in enumerate(self.processors):
            workspace_input += self.W_to_workspace[i] @ proc
        
        workspace_input /= self.n_processors
        
        # Ignition check - use absolute max of input
        input_strength = np.max(np.abs(workspace_input))
        ignited = input_strength > self.ignition_threshold
        
        if ignited:
            self.workspace = self.sigmoid(workspace_input * self.broadcast_gain, gain=3.0)
        else:
            self.workspace = self.workspace * (1 - self.decay_rate) + workspace_input * 0.05
        
        # Top-down: workspace broadcasts to all processors
        for i in range(self.n_processors):
            external = inputs[i] if i < len(inputs) else np.zeros(self.processor_size)
            
            if ignited:
                broadcast = self.W_from_workspace[i] @ self.workspace
                self.processors[i] = self.sigmoid(
                    external + broadcast + np.random.randn(self.processor_size) * noise_level
                )
            else:
                self.processors[i] = self.sigmoid(
                    external * 0.5 + self.processors[i] * 0.5 + 
                    np.random.randn(self.processor_size) * noise_level
                )
        
        return ignited
    
    def run_simulation(self, stimulus_sequence, n_timesteps=200):
        """
        Run full simulation with stimulus sequence.
        
        Returns history of workspace activity, ignition events, and processor states.
        """
        workspace_history = np.zeros((n_timesteps, self.workspace_size))
        processor_history = np.zeros((n_timesteps, self.n_processors, self.processor_size))
        ignition_history = np.zeros(n_timesteps, dtype=bool)
        
        for t in range(n_timesteps):
            if t < len(stimulus_sequence):
                inputs = stimulus_sequence[t]
            else:
                inputs = [np.zeros(self.processor_size)] * self.n_processors
            
            ignited = self.step(inputs)
            
            workspace_history[t] = self.workspace
            for i, proc in enumerate(self.processors):
                processor_history[t, i] = proc
            ignition_history[t] = ignited
        
        return workspace_history, processor_history, ignition_history


def compute_workspace_metrics(workspace_history, processor_history, ignition_history):
    """Compute information-theoretic metrics for workspace dynamics."""
    metrics = {}
    
    # Ignition rate
    metrics['ignition_rate'] = np.mean(ignition_history)
    
    # Workspace entropy over time
    workspace_entropy = []
    for t in range(len(workspace_history)):
        p = workspace_history[t]
        p_norm = np.abs(p) / (np.sum(np.abs(p)) + 1e-10)
        p_norm = p_norm[p_norm > 0]
        workspace_entropy.append(-np.sum(p_norm * np.log2(p_norm + 1e-10)))
    metrics['mean_workspace_entropy'] = np.mean(workspace_entropy)
    
    # Inter-processor synchrony
    n_proc = processor_history.shape[1]
    sync_values = []
    for t in range(len(processor_history)):
        proc_means = [np.mean(processor_history[t, i]) for i in range(n_proc)]
        if np.std(proc_means) > 0:
            sync_values.append(1 - np.std(proc_means) / np.mean(np.abs(proc_means) + 1e-10))
    metrics['mean_synchrony'] = np.mean(sync_values) if sync_values else 0
    
    # Information integration proxy
    all_activity = np.concatenate([workspace_history] + 
                                   [processor_history[:, i] for i in range(n_proc)], axis=1)
    cov = np.cov(all_activity.T)
    eigenvalues = np.linalg.eigvalsh(cov)
    eigenvalues = eigenvalues[eigenvalues > 1e-10]
    metrics['integration_proxy'] = np.sum(np.log(eigenvalues)) if len(eigenvalues) > 0 else 0
    
    return metrics


def compare_gwt_iit(n_trials=10):
    """Compare GWT workspace metrics with IIT-inspired measures across conditions."""
    conditions = {
        'conscious': {'threshold': 0.3, 'gain': 2.0, 'noise': 0.05},
        'subliminal': {'threshold': 0.8, 'gain': 0.5, 'noise': 0.1},
        'anesthesia': {'threshold': 0.9, 'gain': 0.3, 'noise': 0.2},
    }
    
    results = {}
    
    for cond_name, params in conditions.items():
        trial_metrics = []
        
        for trial in range(n_trials):
            np.random.seed(trial * 42 + hash(cond_name) % 1000)
            
            model = GlobalWorkspaceModel()
            model.ignition_threshold = params['threshold']
            model.broadcast_gain = params['gain']
            
            stimulus = []
            for t in range(200):
                if 50 <= t <= 70:
                    inputs = [np.random.rand(model.processor_size) * 0.8 
                             for _ in range(model.n_processors)]
                else:
                    inputs = [np.random.rand(model.processor_size) * 0.1 
                             for _ in range(model.n_processors)]
                stimulus.append(inputs)
            
            wh, ph, ih = model.run_simulation(stimulus, n_timesteps=200)
            metrics = compute_workspace_metrics(wh, ph, ih)
            trial_metrics.append(metrics)
        
        avg_metrics = {}
        for key in trial_metrics[0]:
            vals = [m[key] for m in trial_metrics]
            avg_metrics[key] = {'mean': np.mean(vals), 'std': np.std(vals)}
        
        results[cond_name] = avg_metrics
    
    return results
