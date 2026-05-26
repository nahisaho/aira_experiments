#!/usr/bin/env python3
"""
Main Experiment Runner: Neural Correlates of Consciousness (NCC) 
Information-Theoretic Analysis Framework

Runs all experiments and generates figures for the paper.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import json

from src.iit_phi import (compute_phi_approximate, generate_network, 
                          compute_tpm, effective_information,
                          compute_stochastic_interaction)
from src.pci_simulation import (compute_pci_across_conditions, 
                                 simulate_neural_mass_model, 
                                 apply_tms_perturbation, compute_pci)
from src.consciousness_classifier import (simulate_doc_dataset, 
                                           classify_consciousness_states,
                                           extract_consciousness_features)
from src.global_workspace import compare_gwt_iit, GlobalWorkspaceModel, compute_workspace_metrics

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

all_results = {}

# =============================================================================
# Experiment 1: IIT Φ Computation Across Network Architectures
# =============================================================================
print("=" * 60)
print("Experiment 1: IIT Φ Computation")
print("=" * 60)

network_types = ['integrated', 'modular', 'feedforward', 'disconnected']
n_nodes_range = [3, 4, 5]
phi_results = {}

for net_type in network_types:
    phi_results[net_type] = {}
    for n_nodes in n_nodes_range:
        np.random.seed(42)
        W = generate_network(n_nodes, net_type)
        phi, mip, data, phi_si = compute_phi_approximate(W, n_samples=500, noise_level=0.1)
        phi_results[net_type][n_nodes] = {'phi': phi, 'phi_si': phi_si, 'mip': str(mip)}
        print(f"  {net_type} (n={n_nodes}): Φ_G = {phi:.4f}, Φ_SI = {phi_si:.4f}")

all_results['experiment1_phi'] = {
    net: {str(n): v for n, v in vals.items()} 
    for net, vals in phi_results.items()
}

# Figure 1: Φ across network types
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Bar plot
x_labels = network_types
for idx, n in enumerate(n_nodes_range):
    phi_vals = [phi_results[nt][n]['phi'] for nt in network_types]
    x_pos = np.arange(len(network_types)) + idx * 0.25
    axes[0].bar(x_pos, phi_vals, width=0.2, label=f'n={n}', alpha=0.8)

axes[0].set_xlabel('Network Architecture')
axes[0].set_ylabel('Φ (Integrated Information)')
axes[0].set_title('Integrated Information (Φ) by Network Type')
axes[0].set_xticks(np.arange(len(network_types)) + 0.25)
axes[0].set_xticklabels(network_types, rotation=15)
axes[0].legend()
axes[0].grid(axis='y', alpha=0.3)

# Connectivity matrices visualization
np.random.seed(42)
for idx, net_type in enumerate(['integrated', 'modular']):
    W = generate_network(5, net_type)
    ax_inset = fig.add_axes([0.55 + idx * 0.22, 0.55, 0.18, 0.35])
    im = ax_inset.imshow(W, cmap='RdBu_r', vmin=-1, vmax=1)
    ax_inset.set_title(f'{net_type.capitalize()}', fontsize=9)
    ax_inset.set_xticks([])
    ax_inset.set_yticks([])

axes[1].axis('off')
axes[1].set_title('Example Connectivity Matrices')

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig1_phi_network_types.png'))
plt.close()
print("  → Saved fig1_phi_network_types.png")

# Figure 2: Φ scaling analysis
fig, ax = plt.subplots(figsize=(8, 5))
for net_type in network_types:
    phi_vals = [phi_results[net_type][n]['phi'] for n in n_nodes_range]
    ax.plot(n_nodes_range, phi_vals, 'o-', label=net_type, linewidth=2, markersize=8)

ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Φ (Integrated Information)')
ax.set_title('Scaling of Φ with System Size')
ax.legend()
ax.grid(alpha=0.3)
plt.savefig(os.path.join(FIGURES_DIR, 'fig2_phi_scaling.png'))
plt.close()
print("  → Saved fig2_phi_scaling.png")

# =============================================================================
# Experiment 2: Anesthesia Simulation - Consciousness Level Estimation
# =============================================================================
print("\n" + "=" * 60)
print("Experiment 2: Consciousness Level Estimation (Anesthesia)")
print("=" * 60)

conditions = ['awake', 'light_sedation', 'deep_anesthesia']
anesthesia_results = {}

for cond in conditions:
    np.random.seed(42)
    E, I = simulate_neural_mass_model(16, 1000, consciousness_level=cond)
    
    features = extract_consciousness_features(E, fs=256)
    anesthesia_results[cond] = features
    
    print(f"  {cond}:")
    print(f"    Shannon H = {features['shannon_entropy_mean']:.3f} ± {features['shannon_entropy_std']:.3f}")
    print(f"    Spectral H = {features['spectral_entropy_mean']:.3f} ± {features['spectral_entropy_std']:.3f}")
    print(f"    LZC = {features['lzc_mean']:.3f} ± {features['lzc_std']:.3f}")
    print(f"    Connectivity = {features['mean_connectivity']:.3f}")

all_results['experiment2_anesthesia'] = anesthesia_results

# Figure 3: Information-theoretic metrics across anesthesia levels
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

metrics_to_plot = [
    ('shannon_entropy_mean', 'shannon_entropy_std', 'Shannon Entropy'),
    ('spectral_entropy_mean', 'spectral_entropy_std', 'Spectral Entropy'),
    ('lzc_mean', 'lzc_std', 'Lempel-Ziv Complexity'),
    ('permutation_entropy_mean', 'permutation_entropy_std', 'Permutation Entropy'),
]

colors_cond = {'awake': '#2ecc71', 'light_sedation': '#f39c12', 'deep_anesthesia': '#e74c3c'}

for idx, (mean_key, std_key, title) in enumerate(metrics_to_plot):
    ax = axes[idx // 2, idx % 2]
    means = [anesthesia_results[c][mean_key] for c in conditions]
    stds = [anesthesia_results[c][std_key] for c in conditions]
    colors = [colors_cond[c] for c in conditions]
    
    bars = ax.bar(range(len(conditions)), means, yerr=stds, 
                  color=colors, alpha=0.8, capsize=5, edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(['Awake', 'Light\nSedation', 'Deep\nAnesthesia'])
    ax.set_ylabel(title)
    ax.set_title(title)
    ax.grid(axis='y', alpha=0.3)

plt.suptitle('Information-Theoretic Metrics Across Anesthesia Levels', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig3_anesthesia_metrics.png'))
plt.close()
print("  → Saved fig3_anesthesia_metrics.png")

# =============================================================================
# Experiment 3: PCI Simulation
# =============================================================================
print("\n" + "=" * 60)
print("Experiment 3: Perturbational Complexity Index (PCI)")
print("=" * 60)

pci_results = compute_pci_across_conditions(n_channels=16, n_timepoints=500, n_trials=8)

for cond, vals in pci_results.items():
    print(f"  {cond}: PCI = {vals['mean']:.4f} ± {vals['std']:.4f}")

all_results['experiment3_pci'] = {
    k: {'mean': v['mean'], 'std': v['std']} for k, v in pci_results.items()
}

# Figure 4: PCI across consciousness states
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

cond_labels = ['Awake', 'Light\nSedation', 'Deep\nAnesthesia', 'Vegetative', 'MCS']
cond_keys = ['awake', 'light_sedation', 'deep_anesthesia', 'vegetative', 'mcs']
pci_means = [pci_results[k]['mean'] for k in cond_keys]
pci_stds = [pci_results[k]['std'] for k in cond_keys]
colors_pci = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad', '#3498db']

axes[0].bar(range(len(cond_labels)), pci_means, yerr=pci_stds,
            color=colors_pci, alpha=0.85, capsize=5, edgecolor='black', linewidth=0.5)
axes[0].set_xticks(range(len(cond_labels)))
axes[0].set_xticklabels(cond_labels)
axes[0].set_ylabel('PCI Value')
axes[0].set_title('Perturbational Complexity Index (PCI)')
axes[0].grid(axis='y', alpha=0.3)

# Box plot of individual trials
pci_data = []
pci_labels = []
for k, label in zip(cond_keys, cond_labels):
    pci_data.extend(pci_results[k]['values'])
    pci_labels.extend([label.replace('\n', ' ')] * len(pci_results[k]['values']))

import pandas as pd
try:
    df_pci = pd.DataFrame({'PCI': pci_data, 'Condition': pci_labels})
    sns.boxplot(data=df_pci, x='Condition', y='PCI', ax=axes[1], palette=colors_pci)
except:
    axes[1].boxplot([pci_results[k]['values'] for k in cond_keys], labels=[l.replace('\n',' ') for l in cond_labels])

axes[1].set_title('PCI Distribution Across Conditions')
axes[1].set_ylabel('PCI Value')
axes[1].grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig4_pci_conditions.png'))
plt.close()
print("  → Saved fig4_pci_conditions.png")

# Figure 5: Spatiotemporal response example
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
example_conditions = ['awake', 'deep_anesthesia']

for idx, cond in enumerate(example_conditions):
    np.random.seed(42)
    E, I = simulate_neural_mass_model(16, 500, consciousness_level=cond)
    E_stim = apply_tms_perturbation(E, I, 8, 100)
    response = E_stim - E
    
    axes[0, idx].imshow(E_stim, aspect='auto', cmap='hot', interpolation='nearest')
    axes[0, idx].axvline(x=100, color='cyan', linestyle='--', linewidth=1.5, label='TMS pulse')
    axes[0, idx].set_xlabel('Time (samples)')
    axes[0, idx].set_ylabel('Channel')
    axes[0, idx].set_title(f'EEG Response: {cond.replace("_", " ").title()}')
    axes[0, idx].legend(fontsize=8)
    
    _, binary = compute_pci(E_stim, pre_stim_samples=100)
    axes[1, idx].imshow(binary, aspect='auto', cmap='binary', interpolation='nearest')
    axes[1, idx].set_xlabel('Time (post-stimulus)')
    axes[1, idx].set_ylabel('Channel')
    axes[1, idx].set_title(f'Binary Significant Response: {cond.replace("_", " ").title()}')

plt.suptitle('TMS-EEG Response Patterns', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig5_tms_response_patterns.png'))
plt.close()
print("  → Saved fig5_tms_response_patterns.png")

# =============================================================================
# Experiment 4: Global Workspace Theory Integration
# =============================================================================
print("\n" + "=" * 60)
print("Experiment 4: Global Workspace Theory (GWT)")
print("=" * 60)

gwt_results = compare_gwt_iit(n_trials=8)

for cond, metrics in gwt_results.items():
    print(f"  {cond}:")
    for key, val in metrics.items():
        print(f"    {key}: {val['mean']:.4f} ± {val['std']:.4f}")

all_results['experiment4_gwt'] = {
    cond: {k: v['mean'] for k, v in metrics.items()} 
    for cond, metrics in gwt_results.items()
}

# Figure 6: GWT metrics comparison
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

gwt_conditions = list(gwt_results.keys())
gwt_colors = ['#2ecc71', '#f39c12', '#e74c3c']

metric_keys = ['ignition_rate', 'mean_workspace_entropy', 'mean_synchrony']
metric_titles = ['Ignition Rate', 'Workspace Entropy', 'Inter-processor Synchrony']

for idx, (mkey, mtitle) in enumerate(zip(metric_keys, metric_titles)):
    means = [gwt_results[c][mkey]['mean'] for c in gwt_conditions]
    stds = [gwt_results[c][mkey]['std'] for c in gwt_conditions]
    
    axes[idx].bar(range(len(gwt_conditions)), means, yerr=stds,
                  color=gwt_colors, alpha=0.85, capsize=5, edgecolor='black', linewidth=0.5)
    axes[idx].set_xticks(range(len(gwt_conditions)))
    axes[idx].set_xticklabels([c.capitalize() for c in gwt_conditions])
    axes[idx].set_ylabel(mtitle)
    axes[idx].set_title(mtitle)
    axes[idx].grid(axis='y', alpha=0.3)

plt.suptitle('Global Workspace Theory Metrics', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig6_gwt_metrics.png'))
plt.close()
print("  → Saved fig6_gwt_metrics.png")

# Figure 7: GWT workspace dynamics
fig, axes = plt.subplots(2, 1, figsize=(14, 8))

for idx, (cond, params) in enumerate([('conscious', {'threshold': 0.3, 'gain': 2.0}),
                                        ('anesthesia', {'threshold': 0.9, 'gain': 0.3})]):
    np.random.seed(42)
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
    
    ax = axes[idx]
    im = ax.imshow(wh.T, aspect='auto', cmap='viridis', interpolation='nearest')
    ax.axvline(x=50, color='red', linestyle='--', alpha=0.7, label='Stimulus onset')
    ax.axvline(x=70, color='red', linestyle=':', alpha=0.7, label='Stimulus offset')
    
    ignition_times = np.where(ih)[0]
    for it in ignition_times:
        ax.axvline(x=it, color='yellow', alpha=0.1, linewidth=0.5)
    
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Workspace Unit')
    ax.set_title(f'Workspace Dynamics: {cond.capitalize()} (Ignition rate: {np.mean(ih):.2f})')
    ax.legend(fontsize=8)
    plt.colorbar(im, ax=ax, label='Activation')

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig7_workspace_dynamics.png'))
plt.close()
print("  → Saved fig7_workspace_dynamics.png")

# =============================================================================
# Experiment 5: DoC Classification (VS/UWS vs MCS vs Healthy)
# =============================================================================
print("\n" + "=" * 60)
print("Experiment 5: Disorders of Consciousness Classification")
print("=" * 60)

X, y, feature_names = simulate_doc_dataset(n_subjects_per_class=25, n_channels=8, n_timepoints=500)
classification_results = classify_consciousness_states(X, y, feature_names)

print(f"  SVM Accuracy: {classification_results['SVM']['mean_accuracy']:.3f} ± {classification_results['SVM']['std_accuracy']:.3f}")
print(f"  RF Accuracy: {classification_results['RandomForest']['mean_accuracy']:.3f} ± {classification_results['RandomForest']['std_accuracy']:.3f}")
print(f"  Confusion Matrix:\n{classification_results['confusion_matrix']}")

all_results['experiment5_classification'] = {
    'svm_accuracy': classification_results['SVM']['mean_accuracy'],
    'svm_std': classification_results['SVM']['std_accuracy'],
    'rf_accuracy': classification_results['RandomForest']['mean_accuracy'],
    'rf_std': classification_results['RandomForest']['std_accuracy'],
}

# Figure 8: Classification results
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Confusion matrix
cm = classification_results['confusion_matrix']
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0],
            xticklabels=['VS/UWS', 'MCS', 'Healthy'],
            yticklabels=['VS/UWS', 'MCS', 'Healthy'])
axes[0].set_xlabel('Predicted')
axes[0].set_ylabel('True')
axes[0].set_title('Confusion Matrix (SVM)')

# Classifier comparison
classifiers = ['SVM', 'RandomForest']
accs = [classification_results[c]['mean_accuracy'] for c in classifiers]
stds = [classification_results[c]['std_accuracy'] for c in classifiers]
axes[1].bar(range(len(classifiers)), accs, yerr=stds, 
            color=['#3498db', '#2ecc71'], alpha=0.85, capsize=5, edgecolor='black')
axes[1].set_xticks(range(len(classifiers)))
axes[1].set_xticklabels(['SVM (RBF)', 'Random Forest'])
axes[1].set_ylabel('5-Fold CV Accuracy')
axes[1].set_title('Classifier Comparison')
axes[1].set_ylim(0, 1.1)
axes[1].grid(axis='y', alpha=0.3)

# Feature importance
if 'feature_importance' in classification_results:
    imp = classification_results['feature_importance']
    sorted_imp = sorted(imp.items(), key=lambda x: x[1], reverse=True)[:10]
    names, values = zip(*sorted_imp)
    short_names = [n.replace('_', '\n')[:20] for n in names]
    axes[2].barh(range(len(short_names)), values, color='#9b59b6', alpha=0.8)
    axes[2].set_yticks(range(len(short_names)))
    axes[2].set_yticklabels(short_names, fontsize=8)
    axes[2].set_xlabel('Feature Importance')
    axes[2].set_title('Top Features (Random Forest)')
    axes[2].grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig8_doc_classification.png'))
plt.close()
print("  → Saved fig8_doc_classification.png")

# =============================================================================
# Experiment 6: Artificial Systems Consciousness Assessment
# =============================================================================
print("\n" + "=" * 60)
print("Experiment 6: Artificial Systems Consciousness Criteria")
print("=" * 60)

artificial_systems = {
    'Feedforward NN': {'type': 'feedforward', 'n_nodes': 5},
    'Recurrent NN': {'type': 'integrated', 'n_nodes': 5},
    'Modular NN': {'type': 'modular', 'n_nodes': 4},
    'Disconnected': {'type': 'disconnected', 'n_nodes': 5},
}

ai_results = {}
for name, config in artificial_systems.items():
    np.random.seed(42)
    W = generate_network(config['n_nodes'], config['type'])
    phi, mip, data, phi_si = compute_phi_approximate(W, n_samples=500, noise_level=0.1)
    
    ai_results[name] = {
        'phi': phi,
        'phi_si': phi_si,
        'n_nodes': config['n_nodes'],
        'type': config['type'],
    }
    print(f"  {name}: Φ_G = {phi:.4f}, Φ_SI = {phi_si:.4f}")

all_results['experiment6_artificial'] = ai_results

# Figure 9: Artificial systems consciousness assessment
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sys_names = list(ai_results.keys())
phi_vals = [ai_results[n]['phi'] for n in sys_names]
ei_vals = [ai_results[n]['phi_si'] for n in sys_names]

colors_ai = ['#e74c3c', '#2ecc71', '#f39c12', '#95a5a6']
axes[0].bar(range(len(sys_names)), phi_vals, color=colors_ai, alpha=0.85, 
            edgecolor='black', linewidth=0.5)
axes[0].set_xticks(range(len(sys_names)))
axes[0].set_xticklabels(sys_names, rotation=15)
axes[0].set_ylabel('Φ (Integrated Information)')
axes[0].set_title('Integrated Information in Artificial Systems')
axes[0].axhline(y=0.1, color='red', linestyle='--', alpha=0.5, label='Consciousness threshold (hypothetical)')
axes[0].legend()
axes[0].grid(axis='y', alpha=0.3)

axes[1].scatter(phi_vals, ei_vals, c=colors_ai, s=200, edgecolors='black', linewidth=1, zorder=5)
for i, name in enumerate(sys_names):
    axes[1].annotate(name, (phi_vals[i], ei_vals[i]), textcoords="offset points",
                     xytext=(10, 5), fontsize=9)
axes[1].set_xlabel('Φ_G (Geometric Integrated Information)')
axes[1].set_ylabel('Φ_SI (Stochastic Interaction)')
axes[1].set_title('Φ_G vs Φ_SI')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig9_artificial_systems.png'))
plt.close()
print("  → Saved fig9_artificial_systems.png")

# =============================================================================
# Summary Figure
# =============================================================================
print("\n" + "=" * 60)
print("Generating Summary Figure")
print("=" * 60)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Panel A: Φ by network type
for idx, n in enumerate(n_nodes_range):
    phi_v = [phi_results[nt][n]['phi'] for nt in network_types]
    x_pos = np.arange(len(network_types)) + idx * 0.25
    axes[0, 0].bar(x_pos, phi_v, width=0.2, label=f'n={n}', alpha=0.8)
axes[0, 0].set_title('A. Φ by Network Type')
axes[0, 0].set_xticks(np.arange(len(network_types)) + 0.25)
axes[0, 0].set_xticklabels(network_types, rotation=15, fontsize=8)
axes[0, 0].legend(fontsize=8)
axes[0, 0].set_ylabel('Φ')

# Panel B: Anesthesia entropy
for i, metric in enumerate(['shannon_entropy_mean', 'spectral_entropy_mean']):
    vals = [anesthesia_results[c][metric] for c in conditions]
    axes[0, 1].bar(np.arange(len(conditions)) + i*0.3, vals, width=0.25, 
                   label=metric.replace('_mean','').replace('_',' ').title(), alpha=0.8)
axes[0, 1].set_title('B. Entropy Under Anesthesia')
axes[0, 1].set_xticks(np.arange(len(conditions)) + 0.15)
axes[0, 1].set_xticklabels(['Awake', 'Light Sed.', 'Deep Anes.'], fontsize=8)
axes[0, 1].legend(fontsize=7)
axes[0, 1].set_ylabel('Entropy')

# Panel C: PCI
pci_m = [pci_results[k]['mean'] for k in cond_keys]
pci_s = [pci_results[k]['std'] for k in cond_keys]
axes[0, 2].bar(range(len(cond_keys)), pci_m, yerr=pci_s, color=colors_pci, alpha=0.85, capsize=3)
axes[0, 2].set_title('C. PCI Across States')
axes[0, 2].set_xticks(range(len(cond_keys)))
axes[0, 2].set_xticklabels(['Awake', 'Light', 'Deep', 'VS', 'MCS'], fontsize=8)
axes[0, 2].set_ylabel('PCI')

# Panel D: GWT ignition
for i, mkey in enumerate(['ignition_rate', 'mean_workspace_entropy']):
    vals = [gwt_results[c][mkey]['mean'] for c in gwt_conditions]
    axes[1, 0].bar(np.arange(len(gwt_conditions)) + i*0.3, vals, width=0.25,
                   label=mkey.replace('_', ' ').title(), alpha=0.8)
axes[1, 0].set_title('D. GWT Metrics')
axes[1, 0].set_xticks(np.arange(len(gwt_conditions)) + 0.15)
axes[1, 0].set_xticklabels([c.capitalize() for c in gwt_conditions], fontsize=8)
axes[1, 0].legend(fontsize=7)

# Panel E: Classification
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1, 1],
            xticklabels=['VS', 'MCS', 'H'], yticklabels=['VS', 'MCS', 'H'])
axes[1, 1].set_title('E. DoC Classification')

# Panel F: AI systems
axes[1, 2].bar(range(len(sys_names)), phi_vals, color=colors_ai, alpha=0.85)
axes[1, 2].set_xticks(range(len(sys_names)))
axes[1, 2].set_xticklabels(sys_names, rotation=15, fontsize=8)
axes[1, 2].set_title('F. AI Systems Φ')
axes[1, 2].set_ylabel('Φ')

plt.suptitle('Neural Correlates of Consciousness: Information-Theoretic Analysis', fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig10_summary.png'))
plt.close()
print("  → Saved fig10_summary.png")

# Save all results
with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results.json'), 'w') as f:
    json.dump(all_results, f, indent=2, default=str)

print("\n" + "=" * 60)
print("ALL EXPERIMENTS COMPLETED SUCCESSFULLY")
print("=" * 60)
