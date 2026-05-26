#!/usr/bin/env python3
"""
Main experiment runner for SNN simulation framework.
Generates all figures and results for report.md and paper.md.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import json
import time

from src.neuron_models import HodgkinHuxley, Izhikevich, AdEx, benchmark_models
from src.plasticity import STDP, HomeostaticPlasticity, simulate_stdp_learning
from src.gpu_architecture import GPUBlockSimulator, benchmark_scaling, estimate_gpu_performance
from src.potjans_diesmann import PotjansDiesmannCircuit
from src.analysis_tools import (compute_firing_rates, compute_cv_isi,
                                 compute_phase_synchrony, compute_transfer_entropy,
                                 compute_power_spectrum, compute_fano_factor)
from src.working_memory import WorkingMemoryNetwork, run_working_memory_experiment

np.random.seed(42)
os.makedirs('figures', exist_ok=True)
all_results = {}

# ============================================================
# Experiment 1: Neuron Model Comparison
# ============================================================
print("=" * 60)
print("Experiment 1: Neuron Model Comparison")
print("=" * 60)

fig, axes = plt.subplots(3, 2, figsize=(14, 12))
fig.suptitle('Neuron Model Comparison: HH vs Izhikevich vs AdEx', fontsize=14, fontweight='bold')

models_info = {
    'HH': {'class': HodgkinHuxley, 'dt': 0.01, 'I_range': np.arange(0, 25, 1)},
    'Izhikevich': {'class': Izhikevich, 'dt': 0.1, 'I_range': np.arange(0, 25, 1)},
    'AdEx': {'class': AdEx, 'dt': 0.1, 'I_range': np.arange(0, 600, 25)},
}

for row, (name, info) in enumerate(models_info.items()):
    # Voltage trace
    model = info['class'](N=1)
    T = 500
    dt = info['dt']
    steps = int(T / dt)
    V_trace = np.zeros(steps)
    t_arr = np.arange(steps) * dt

    I_val = info['I_range'][len(info['I_range'])//2]
    for t in range(steps):
        I_ext = np.array([I_val])
        V, spikes = model.step(dt, I_ext)
        V_trace[t] = V[0]

    axes[row, 0].plot(t_arr, V_trace, 'k', linewidth=0.8)
    axes[row, 0].set_title(f'{name}: Voltage Trace (I={I_val})')
    axes[row, 0].set_xlabel('Time (ms)')
    axes[row, 0].set_ylabel('V (mV)')
    axes[row, 0].set_xlim(0, T)

    # F-I curve
    fi_rates = []
    for I_val in info['I_range']:
        model_fi = info['class'](N=1)
        spike_count = 0
        T_fi = 1000
        steps_fi = int(T_fi / dt)
        for t in range(steps_fi):
            V, spikes = model_fi.step(dt, np.array([I_val]))
            spike_count += spikes[0]
        fi_rates.append(spike_count / (T_fi / 1000))

    axes[row, 1].plot(info['I_range'], fi_rates, 'b-o', markersize=3)
    axes[row, 1].set_title(f'{name}: F-I Curve')
    axes[row, 1].set_xlabel('Input Current')
    axes[row, 1].set_ylabel('Firing Rate (Hz)')

plt.tight_layout()
plt.savefig('figures/neuron_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/neuron_comparison.png")

# Benchmark
print("  Running benchmark (N=1000, T=1000ms)...")
bench_results = benchmark_models(N=1000, T=1000, dt=0.1)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Neuron Model Benchmark (N=1000, T=1000ms)', fontsize=13, fontweight='bold')

names = list(bench_results.keys())
times = [bench_results[n]['time'] for n in names]
rates = [bench_results[n]['mean_rate'] for n in names]
rate_stds = [bench_results[n]['std_rate'] for n in names]

colors = ['#2196F3', '#4CAF50', '#FF9800']
axes[0].bar(names, times, color=colors)
axes[0].set_ylabel('Computation Time (s)')
axes[0].set_title('Simulation Speed')

axes[1].bar(names, rates, yerr=rate_stds, color=colors, capsize=5)
axes[1].set_ylabel('Firing Rate (Hz)')
axes[1].set_title('Mean Firing Rate')

# Spike count distributions
for i, n in enumerate(names):
    axes[2].hist(bench_results[n]['spike_counts'], bins=30, alpha=0.6, label=n, color=colors[i])
axes[2].set_xlabel('Spike Count')
axes[2].set_ylabel('Neuron Count')
axes[2].set_title('Spike Count Distribution')
axes[2].legend()

plt.tight_layout()
plt.savefig('figures/benchmark_results.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/benchmark_results.png")

all_results['benchmark'] = {
    n: {'time_s': bench_results[n]['time'],
        'mean_rate_hz': bench_results[n]['mean_rate'],
        'std_rate_hz': bench_results[n]['std_rate']}
    for n in names
}

# ============================================================
# Experiment 2: STDP and Homeostatic Plasticity
# ============================================================
print("\n" + "=" * 60)
print("Experiment 2: STDP and Homeostatic Plasticity")
print("=" * 60)

stdp_results = simulate_stdp_learning(N_pre=100, N_post=50, T=5000, dt=0.5)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Synaptic Plasticity: STDP + Homeostatic Scaling', fontsize=14, fontweight='bold')

axes[0, 0].plot(stdp_results['weight_history'], 'b-', linewidth=1.5)
axes[0, 0].set_xlabel('Time Step (×100)')
axes[0, 0].set_ylabel('Mean Weight')
axes[0, 0].set_title('STDP Weight Evolution')

axes[0, 1].plot(stdp_results['rate_history'], 'r-', linewidth=1.5)
axes[0, 1].axhline(y=5.0, color='k', linestyle='--', label='Target Rate')
axes[0, 1].set_xlabel('Time Step (×100)')
axes[0, 1].set_ylabel('Firing Rate (Hz)')
axes[0, 1].set_title('Homeostatic Rate Regulation')
axes[0, 1].legend()

im = axes[1, 0].imshow(stdp_results['final_weights'], aspect='auto', cmap='hot')
axes[1, 0].set_xlabel('Post-synaptic')
axes[1, 0].set_ylabel('Pre-synaptic')
axes[1, 0].set_title('Final Weight Matrix')
plt.colorbar(im, ax=axes[1, 0])

axes[1, 1].hist(stdp_results['final_weights'].flatten(), bins=50, color='purple', alpha=0.7)
axes[1, 1].set_xlabel('Weight')
axes[1, 1].set_ylabel('Count')
axes[1, 1].set_title('Weight Distribution')

plt.tight_layout()
plt.savefig('figures/plasticity_results.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/plasticity_results.png")

all_results['plasticity'] = {
    'final_mean_weight': float(stdp_results['final_weights'].mean()),
    'final_std_weight': float(stdp_results['final_weights'].std()),
    'final_mean_rate': float(stdp_results['rate_history'][-1]),
}

# ============================================================
# Experiment 3: GPU Scaling Analysis
# ============================================================
print("\n" + "=" * 60)
print("Experiment 3: GPU Scaling Analysis")
print("=" * 60)

sizes = [1000, 5000, 10000, 50000]
scaling_results = benchmark_scaling(sizes=sizes, T=50, dt=0.5)

gpu_est = estimate_gpu_performance()

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('GPU Parallel Scaling Analysis', fontsize=14, fontweight='bold')

ns = list(scaling_results.keys())
ts = [scaling_results[n]['time_s'] for n in ns]
throughput = [scaling_results[n]['neurons_per_s'] for n in ns]

axes[0].loglog(ns, ts, 'bo-', linewidth=2, markersize=8)
axes[0].set_xlabel('Network Size (neurons)')
axes[0].set_ylabel('Simulation Time (s)')
axes[0].set_title('Scaling: Time vs Network Size')
axes[0].grid(True, alpha=0.3)

axes[1].semilogx(ns, throughput, 'rs-', linewidth=2, markersize=8)
axes[1].set_xlabel('Network Size (neurons)')
axes[1].set_ylabel('Throughput (neuron·steps/s)')
axes[1].set_title('Computational Throughput')
axes[1].grid(True, alpha=0.3)

# GPU memory estimation
mem_sizes = [10000, 100000, 500000, 1000000, 5000000]
neuron_mem = [n * 48 / 1e9 for n in mem_sizes]
synapse_mem = [n * 1000 * 8 / 1e9 for n in mem_sizes]
total_mem = [nm + sm for nm, sm in zip(neuron_mem, synapse_mem)]

axes[2].semilogy(mem_sizes, neuron_mem, 'b--', label='Neuron State')
axes[2].semilogy(mem_sizes, synapse_mem, 'r--', label='Synaptic Weights')
axes[2].semilogy(mem_sizes, total_mem, 'k-', linewidth=2, label='Total')
axes[2].axhline(y=80, color='green', linestyle=':', label='A100 80GB')
axes[2].axhline(y=40, color='orange', linestyle=':', label='A6000 48GB')
axes[2].set_xlabel('Network Size (neurons)')
axes[2].set_ylabel('Memory (GB)')
axes[2].set_title('GPU Memory Requirements')
axes[2].legend(fontsize=8)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/gpu_scaling.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/gpu_scaling.png")

all_results['gpu'] = {
    'scaling': {str(n): {'time_s': scaling_results[n]['time_s'],
                          'throughput': scaling_results[n]['neurons_per_s']} 
                for n in ns},
    'gpu_estimate': gpu_est,
}

# ============================================================
# Experiment 4: Potjans-Diesmann Cortical Microcircuit
# ============================================================
print("\n" + "=" * 60)
print("Experiment 4: Potjans-Diesmann Cortical Microcircuit")
print("=" * 60)

print("  Building circuit (scale=0.1)...")
circuit = PotjansDiesmannCircuit(scale=0.1)
print(f"  Total neurons: {circuit.N_total}")

print("  Simulating 1000ms...")
pd_results = circuit.simulate(T=1000, dt=0.5, stim_pop='L4E', stim_start=300, stim_end=500, stim_amp=15)

fig = plt.figure(figsize=(16, 14))
gs = GridSpec(4, 2, figure=fig)
fig.suptitle('Potjans-Diesmann Cortical Microcircuit (10% scale)', fontsize=14, fontweight='bold')

pop_names = circuit.pop_names
colors_pd = plt.cm.tab10(np.linspace(0, 1, len(pop_names)))
time_axis = np.arange(len(pd_results['pop_rates'][pop_names[0]])) * pd_results['dt']

# Population firing rates
for i, (layer, pops) in enumerate(zip(['L2/3', 'L4', 'L5', 'L6'],
                                       [('L2/3E', 'L2/3I'), ('L4E', 'L4I'),
                                        ('L5E', 'L5I'), ('L6E', 'L6I')])):
    ax = fig.add_subplot(gs[i, 0])
    for pop_name in pops:
        rate = pd_results['pop_rates'][pop_name]
        # Smooth with moving average
        kernel = np.ones(20) / 20
        smoothed = np.convolve(rate, kernel, mode='same')
        color = 'blue' if 'E' in pop_name else 'red'
        ax.plot(time_axis, smoothed, color=color, alpha=0.8, label=pop_name, linewidth=1)
    ax.set_ylabel('Rate (Hz)')
    ax.set_title(f'{layer}')
    ax.legend(loc='upper right', fontsize=8)
    ax.axvspan(300, 500, alpha=0.15, color='yellow', label='Stim')
    if i == 3:
        ax.set_xlabel('Time (ms)')

# Raster plots
for i, (layer, pops) in enumerate(zip(['L2/3', 'L4', 'L5', 'L6'],
                                       [('L2/3E', 'L2/3I'), ('L4E', 'L4I'),
                                        ('L5E', 'L5I'), ('L6E', 'L6I')])):
    ax = fig.add_subplot(gs[i, 1])
    offset = 0
    for pop_name in pops:
        spikes = pd_results['spike_trains'][pop_name]
        if spikes:
            ts, ns = zip(*spikes[:2000])
            color = 'blue' if 'E' in pop_name else 'red'
            ax.scatter(ts, np.array(ns) + offset, s=0.3, c=color, alpha=0.5)
        offset += circuit.N_pops[pop_name]
    ax.set_title(f'{layer} Raster')
    ax.set_ylabel('Neuron ID')
    if i == 3:
        ax.set_xlabel('Time (ms)')
    ax.axvspan(300, 500, alpha=0.15, color='yellow')

plt.tight_layout()
plt.savefig('figures/potjans_diesmann.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/potjans_diesmann.png")

# Compute statistics
pd_stats = {}
for pop_name in pop_names:
    rate = pd_results['pop_rates'][pop_name]
    spikes = pd_results['spike_trains'][pop_name]
    n_neurons = circuit.N_pops[pop_name]
    
    mean_rate = len(spikes) / n_neurons / (pd_results['T'] / 1000)
    cv_isi, cv_std = compute_cv_isi(spikes, n_neurons, pd_results['T'])
    
    pd_stats[pop_name] = {
        'mean_rate_hz': float(mean_rate),
        'cv_isi': float(cv_isi),
        'n_neurons': n_neurons,
    }
    print(f"  {pop_name}: rate={mean_rate:.1f} Hz, CV_ISI={cv_isi:.2f}")

all_results['potjans_diesmann'] = pd_stats

# ============================================================
# Experiment 5: Analysis Tools Demonstration
# ============================================================
print("\n" + "=" * 60)
print("Experiment 5: Analysis Tools")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Neural Signal Analysis', fontsize=14, fontweight='bold')

# Phase synchrony between L4E and L2/3E
rate_L4E = pd_results['pop_rates']['L4E']
rate_L23E = pd_results['pop_rates']['L2/3E']
sync = compute_phase_synchrony(rate_L4E, rate_L23E, pd_results['dt'])
print(f"  Phase synchrony L4E-L2/3E: PLV={sync['PLV']:.3f}, PLI={sync['PLI']:.3f}")

# Power spectrum
freqs_L4, psd_L4 = compute_power_spectrum(rate_L4E, pd_results['dt'])
freqs_L23, psd_L23 = compute_power_spectrum(rate_L23E, pd_results['dt'])

axes[0, 0].semilogy(freqs_L4, psd_L4, 'b-', label='L4E')
axes[0, 0].semilogy(freqs_L23, psd_L23, 'r-', label='L2/3E')
axes[0, 0].set_xlabel('Frequency (Hz)')
axes[0, 0].set_ylabel('Power')
axes[0, 0].set_title('Power Spectral Density')
axes[0, 0].legend()
axes[0, 0].set_xlim(0, 200)

# Transfer entropy
te_values = []
pop_pairs = [('L4E', 'L2/3E'), ('L2/3E', 'L5E'), ('L5E', 'L6E'), ('L6E', 'L4E')]
for src, tgt in pop_pairs:
    te = compute_transfer_entropy(pd_results['pop_rates'][src][:500],
                                   pd_results['pop_rates'][tgt][:500])
    te_values.append(te)
    print(f"  TE({src}→{tgt}) = {te:.4f}")

pair_labels = [f'{s}→{t}' for s, t in pop_pairs]
axes[0, 1].bar(pair_labels, te_values, color=['#2196F3', '#4CAF50', '#FF9800', '#9C27B0'])
axes[0, 1].set_ylabel('Transfer Entropy (bits)')
axes[0, 1].set_title('Information Transfer Between Layers')
axes[0, 1].tick_params(axis='x', rotation=30)

# Fano factors
fano_values = []
for pop_name in ['L2/3E', 'L4E', 'L5E', 'L6E']:
    ff = compute_fano_factor(pd_results['spike_trains'][pop_name],
                              circuit.N_pops[pop_name], pd_results['T'])
    fano_values.append(ff)

axes[1, 0].bar(['L2/3E', 'L4E', 'L5E', 'L6E'], fano_values, color='steelblue')
axes[1, 0].axhline(y=1.0, color='red', linestyle='--', label='Poisson')
axes[1, 0].set_ylabel('Fano Factor')
axes[1, 0].set_title('Spike Count Variability')
axes[1, 0].legend()

# Phase locking values between all E populations
e_pops = ['L2/3E', 'L4E', 'L5E', 'L6E']
plv_matrix = np.zeros((4, 4))
for i, p1 in enumerate(e_pops):
    for j, p2 in enumerate(e_pops):
        if i != j:
            s = compute_phase_synchrony(pd_results['pop_rates'][p1],
                                         pd_results['pop_rates'][p2],
                                         pd_results['dt'])
            plv_matrix[i, j] = s['PLV']
        else:
            plv_matrix[i, j] = 1.0

im = axes[1, 1].imshow(plv_matrix, cmap='YlOrRd', vmin=0, vmax=1)
axes[1, 1].set_xticks(range(4))
axes[1, 1].set_yticks(range(4))
axes[1, 1].set_xticklabels(e_pops)
axes[1, 1].set_yticklabels(e_pops)
axes[1, 1].set_title('Phase Locking Value Matrix')
plt.colorbar(im, ax=axes[1, 1])

for i in range(4):
    for j in range(4):
        axes[1, 1].text(j, i, f'{plv_matrix[i,j]:.2f}', ha='center', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('figures/analysis_tools.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/analysis_tools.png")

all_results['analysis'] = {
    'phase_sync_L4E_L23E': sync,
    'transfer_entropy': dict(zip(pair_labels, te_values)),
    'fano_factors': dict(zip(e_pops, fano_values)),
}

# ============================================================
# Experiment 6: Working Memory Task
# ============================================================
print("\n" + "=" * 60)
print("Experiment 6: Working Memory Task")
print("=" * 60)

print("  Running working memory experiment (10 trials)...")
wm_results = run_working_memory_experiment(n_trials=5)

print(f"  Baseline rate: {wm_results['baseline_rate']:.1f} ± {wm_results['baseline_std']:.1f} Hz")
print(f"  Stimulus rate: {wm_results['stim_rate']:.1f} ± {wm_results['stim_std']:.1f} Hz")
print(f"  Delay rate: {wm_results['delay_rate']:.1f} ± {wm_results['delay_std']:.1f} Hz")

fig = plt.figure(figsize=(16, 12))
gs = GridSpec(3, 2, figure=fig)
fig.suptitle('Working Memory: Delayed Match-to-Sample Task', fontsize=14, fontweight='bold')

# Single trial example
trial = wm_results['trials'][0]
ax1 = fig.add_subplot(gs[0, :])
for p in range(trial['pool_rates'].shape[0]):
    kernel = np.ones(40) / 40
    smoothed = np.convolve(trial['pool_rates'][p], kernel, mode='same')
    label = f'Pool {p}' + (' (stim)' if p == trial['stim_pool'] else '')
    lw = 2 if p == trial['stim_pool'] else 1
    ax1.plot(trial['times'], smoothed, linewidth=lw, label=label, alpha=0.8)

ax1.axvspan(300, 500, alpha=0.2, color='green', label='Sample')
ax1.axvspan(1500, 1700, alpha=0.2, color='orange', label='Probe')
ax1.axvspan(500, 1500, alpha=0.1, color='gray', label='Delay')
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Firing Rate (Hz)')
ax1.set_title('Single Trial: Pool-Specific Firing Rates')
ax1.legend(loc='upper right', fontsize=8)

# Mean rates across trials
ax2 = fig.add_subplot(gs[1, 0])
phases = ['Baseline', 'Stimulus', 'Delay']
means = [wm_results['baseline_rate'], wm_results['stim_rate'], wm_results['delay_rate']]
stds = [wm_results['baseline_std'], wm_results['stim_std'], wm_results['delay_std']]
bars = ax2.bar(phases, means, yerr=stds, color=['gray', 'green', 'orange'], capsize=8)
ax2.set_ylabel('Firing Rate (Hz)')
ax2.set_title('Mean Firing Rates Across Task Phases')

# Comparison with experimental data
ax3 = fig.add_subplot(gs[1, 1])
exp_rates = {'Baseline': 5.0, 'Stimulus': 25.0, 'Delay': 12.0}
model_rates = {'Baseline': wm_results['baseline_rate'],
               'Stimulus': wm_results['stim_rate'],
               'Delay': wm_results['delay_rate']}

x = np.arange(3)
width = 0.35
ax3.bar(x - width/2, list(exp_rates.values()), width, label='Experimental (ref)', color='steelblue')
ax3.bar(x + width/2, list(model_rates.values()), width, label='Model', color='coral')
ax3.set_xticks(x)
ax3.set_xticklabels(phases)
ax3.set_ylabel('Firing Rate (Hz)')
ax3.set_title('Model vs Experimental Comparison')
ax3.legend()

# Inhibitory activity
ax4 = fig.add_subplot(gs[2, 0])
kernel = np.ones(40) / 40
inh_smooth = np.convolve(trial['inh_rates'], kernel, mode='same')
ax4.plot(trial['times'], inh_smooth, 'r-', linewidth=1.5)
ax4.axvspan(300, 500, alpha=0.2, color='green')
ax4.axvspan(1500, 1700, alpha=0.2, color='orange')
ax4.set_xlabel('Time (ms)')
ax4.set_ylabel('Firing Rate (Hz)')
ax4.set_title('Inhibitory Population Activity')

# Trial variability
ax5 = fig.add_subplot(gs[2, 1])
delay_rates_per_trial = []
for r in wm_results['trials']:
    sp = r['stim_pool']
    delay_mask = (r['times'] >= 600) & (r['times'] < 1500)
    delay_rates_per_trial.append(r['pool_rates'][sp, delay_mask].mean())

ax5.bar(range(len(delay_rates_per_trial)), delay_rates_per_trial, color='orange', alpha=0.7)
ax5.axhline(y=np.mean(delay_rates_per_trial), color='red', linestyle='--', label='Mean')
ax5.set_xlabel('Trial')
ax5.set_ylabel('Delay Period Rate (Hz)')
ax5.set_title('Trial-by-Trial Variability')
ax5.legend()

plt.tight_layout()
plt.savefig('figures/working_memory.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/working_memory.png")

all_results['working_memory'] = {
    'baseline_rate': float(wm_results['baseline_rate']),
    'stim_rate': float(wm_results['stim_rate']),
    'delay_rate': float(wm_results['delay_rate']),
    'baseline_std': float(wm_results['baseline_std']),
    'stim_std': float(wm_results['stim_std']),
    'delay_std': float(wm_results['delay_std']),
}

# ============================================================
# Architecture Overview Figure
# ============================================================
print("\n" + "=" * 60)
print("Creating Architecture Overview")
print("=" * 60)

fig, ax = plt.subplots(1, 1, figsize=(14, 8))
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('SNN Simulation Framework Architecture', fontsize=16, fontweight='bold', pad=20)

# Draw modules
modules = [
    (1, 6, 3, 1.2, 'Neuron Models\n(HH/Izh/AdEx)', '#E3F2FD'),
    (5, 6, 3, 1.2, 'Synaptic Plasticity\n(STDP/Homeo)', '#E8F5E9'),
    (9, 6, 3, 1.2, 'GPU Engine\n(CUDA Blocks)', '#FFF3E0'),
    (1, 3.5, 3, 1.2, 'Cortical Circuit\n(Potjans-Diesmann)', '#F3E5F5'),
    (5, 3.5, 3, 1.2, 'Analysis Tools\n(Rate/Sync/TE)', '#FFEBEE'),
    (9, 3.5, 3, 1.2, 'Working Memory\n(DMTS Task)', '#E0F7FA'),
    (3, 1, 6, 1.2, 'Unified API & Visualization', '#FFF9C4'),
]

for x, y, w, h, label, color in modules:
    rect = plt.Rectangle((x, y), w, h, linewidth=2, edgecolor='black',
                           facecolor=color, zorder=2)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, label, ha='center', va='center',
            fontsize=10, fontweight='bold', zorder=3)

# Arrows
arrow_style = dict(arrowstyle='->', color='gray', lw=1.5)
from matplotlib.patches import FancyArrowPatch
arrows = [
    ((2.5, 6), (2.5, 4.7)),
    ((6.5, 6), (6.5, 4.7)),
    ((10.5, 6), (10.5, 4.7)),
    ((4, 4.1), (5, 4.1)),
    ((8, 4.1), (9, 4.1)),
    ((6, 3.5), (6, 2.2)),
]
for (x1, y1), (x2, y2) in arrows:
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='gray', lw=2))

plt.savefig('figures/architecture.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/architecture.png")

# Save all results
with open('results.json', 'w') as f:
    def convert(o):
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        return str(o)
    json.dump(all_results, f, indent=2, default=convert)

print("\n" + "=" * 60)
print("All experiments complete!")
print("=" * 60)
print(f"Results saved to: results.json")
print(f"Figures saved to: figures/")
