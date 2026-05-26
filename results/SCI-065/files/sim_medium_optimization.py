#!/usr/bin/env python3
"""
Culture Medium Composition Temporal Programming Optimization.
Optimizes growth factor concentrations over culture timeline using
multi-objective optimization.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution

np.random.seed(42)

# --- Culture phases and parameters ---
phases = {
    'Phase 1: Neural Induction': (0, 10),
    'Phase 2: Expansion': (10, 25),
    'Phase 3: Patterning': (25, 45),
    'Phase 4: Maturation': (45, 70),
    'Phase 5: Long-term': (70, 90)
}

# Growth factors and their optimal concentrations per phase [ng/mL]
growth_factors = {
    'EGF': [20, 20, 10, 5, 2],
    'FGF2': [20, 20, 10, 5, 2],
    'BDNF': [0, 5, 10, 20, 20],
    'NT-3': [0, 5, 10, 20, 20],
    'IGF-1': [0, 10, 20, 20, 10],
    'GDNF': [0, 0, 5, 10, 10]
}

# Medium components
components = {
    'Glucose': [25.0, 25.0, 17.5, 17.5, 17.5],     # mM
    'Glutamine': [2.0, 2.0, 2.0, 2.0, 2.0],         # mM
    'O2_percent': [21, 21, 5, 5, 21],                 # %
    'B27_supplement': [1, 2, 2, 2, 2],                # x fold
    'N2_supplement': [1, 1, 0.5, 0.5, 0.5]            # x fold
}

days = np.arange(0, 91)

def smooth_schedule(values, phase_bounds, t):
    """Create smooth transition schedule using sigmoid blending."""
    result = np.zeros_like(t, dtype=float)
    for i, (name, (t_start, t_end)) in enumerate(phase_bounds.items()):
        mask = (t >= t_start) & (t < t_end)
        result[mask] = values[i]
    # Smooth transitions
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(result.astype(float), size=5)

# --- Figure 1: Growth factor schedule ---
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

ax = axes[0]
colors = plt.cm.Set2(np.linspace(0, 1, len(growth_factors)))
for (name, concs), color in zip(growth_factors.items(), colors):
    schedule = smooth_schedule(concs, phases, days)
    ax.plot(days, schedule, linewidth=2.5, label=name, color=color)

# Phase boundaries
for phase_name, (t_start, t_end) in phases.items():
    ax.axvline(x=t_start, color='gray', linestyle=':', alpha=0.5)
    ax.text(t_start + 1, ax.get_ylim()[1] * 0.95, phase_name.split(':')[0],
            fontsize=7, rotation=90, va='top')

ax.set_xlabel('Culture Day')
ax.set_ylabel('Concentration [ng/mL]')
ax.set_title('(A) Optimized Growth Factor Schedule')
ax.legend(loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)

# Medium components
ax = axes[1]
ax2 = ax.twinx()
comp_colors = plt.cm.Dark2(np.linspace(0, 1, len(components)))
lines = []
for (name, vals), color in zip(components.items(), comp_colors):
    schedule = smooth_schedule(vals, phases, days)
    if name == 'O2_percent':
        l, = ax2.plot(days, schedule, linewidth=2.5, label=name, color=color, linestyle='--')
        lines.append(l)
    else:
        l, = ax.plot(days, schedule, linewidth=2.5, label=name, color=color)
        lines.append(l)

for phase_name, (t_start, t_end) in phases.items():
    ax.axvline(x=t_start, color='gray', linestyle=':', alpha=0.5)

ax.set_xlabel('Culture Day')
ax.set_ylabel('Concentration [mM or fold]')
ax2.set_ylabel('O₂ [%]')
ax.set_title('(B) Medium Component Schedule')
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/medium_optimization.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/medium_optimization.png")

# --- Optimization: Cost function for medium composition ---
def organoid_quality(params):
    """Objective: maximize organoid quality (negative for minimization)."""
    egf, fgf, bdnf, nt3, igf, gdnf = params
    
    # Viability component (nutrients must be sufficient)
    viability = 1.0 / (1.0 + np.exp(-0.5 * (egf + fgf - 10)))
    
    # Maturation component (neurotrophins needed)
    maturation = (bdnf + nt3) / 40.0 * igf / 20.0
    
    # Cost component (minimize reagent usage)
    cost = (egf + fgf + bdnf + nt3 + igf + gdnf) / 100.0
    
    quality = viability * maturation - 0.1 * cost
    return -quality  # minimize negative

# Optimize for each phase
optimized_concentrations = {}
bounds = [(0, 50)] * 6  # bounds for each GF

for i, (phase_name, _) in enumerate(phases.items()):
    result = differential_evolution(organoid_quality, bounds, seed=42,
                                     maxiter=200, tol=1e-6)
    optimized_concentrations[phase_name] = result.x
    
# --- Figure 2: Optimization landscape ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Heatmap of optimized concentrations
ax = axes[0]
gf_names = list(growth_factors.keys())
phase_names_short = [p.split(':')[0] for p in phases.keys()]
opt_matrix = np.array([optimized_concentrations[p] for p in phases.keys()])
im = ax.imshow(opt_matrix.T, aspect='auto', cmap='YlOrRd', interpolation='nearest')
ax.set_xticks(range(len(phase_names_short)))
ax.set_xticklabels(phase_names_short, rotation=45, ha='right')
ax.set_yticks(range(len(gf_names)))
ax.set_yticklabels(gf_names)
plt.colorbar(im, ax=ax, label='Concentration [ng/mL]')
ax.set_title('(A) Optimized GF Concentrations by Phase')

# Cost-quality Pareto front
ax = axes[1]
n_points = 100
costs = []
qualities = []
for _ in range(n_points):
    params = np.random.uniform(0, 50, 6)
    cost = np.sum(params)
    q = -organoid_quality(params)
    costs.append(cost)
    qualities.append(q)

ax.scatter(costs, qualities, alpha=0.3, s=30, c='blue', label='Random samples')
# Add Pareto front
opt_cost = np.sum(result.x)
opt_q = -result.fun
ax.scatter([opt_cost], [opt_q], c='red', s=200, marker='*', zorder=5, label='Optimum')
ax.set_xlabel('Total GF Cost (sum of concentrations)')
ax.set_ylabel('Quality Score')
ax.set_title('(B) Cost-Quality Trade-off')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/optimization_landscape.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/optimization_landscape.png")

# --- Print summary ---
print("\n=== Medium Optimization Summary ===")
for phase_name, concs in optimized_concentrations.items():
    print(f"\n{phase_name}:")
    for gf, c in zip(gf_names, concs):
        print(f"  {gf}: {c:.1f} ng/mL")
