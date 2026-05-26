#!/usr/bin/env python3
"""
Biomarker Monitoring Strategy for Brain Organoid Maturation Assessment.
Models temporal expression of key neural biomarkers and monitoring protocols.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(42)

# --- Biomarker temporal models ---
days = np.arange(0, 91)

def logistic(t, L, k, t0):
    """Logistic growth model."""
    return L / (1 + np.exp(-k * (t - t0)))

def gaussian_pulse(t, peak_day, width, amplitude):
    """Transient expression peak."""
    return amplitude * np.exp(-((t - peak_day)**2) / (2 * width**2))

# Neural progenitor markers (decrease with maturation)
SOX2 = 1.0 * np.exp(-0.03 * days) + 0.05 * np.random.randn(len(days)) * 0.1
PAX6 = gaussian_pulse(days, 15, 10, 0.9) + 0.05 * np.random.randn(len(days)) * 0.05
NESTIN = 0.8 * np.exp(-0.025 * days) + gaussian_pulse(days, 10, 8, 0.3)

# Neuronal markers (increase with maturation)
TUJ1 = logistic(days, 0.8, 0.1, 25) + np.random.randn(len(days)) * 0.02
MAP2 = logistic(days, 0.85, 0.08, 35) + np.random.randn(len(days)) * 0.02
NEUN = logistic(days, 0.75, 0.07, 45) + np.random.randn(len(days)) * 0.02

# Synaptic markers (late maturation)
SYN1 = logistic(days, 0.7, 0.06, 55) + np.random.randn(len(days)) * 0.02
PSD95 = logistic(days, 0.65, 0.05, 60) + np.random.randn(len(days)) * 0.02

# Glial markers
GFAP = logistic(days, 0.6, 0.05, 50) + np.random.randn(len(days)) * 0.02
OLIG2 = logistic(days, 0.4, 0.04, 55) + np.random.randn(len(days)) * 0.015

# Metabolic indicators
lactate = 2.0 + 3.0 * logistic(days, 1, 0.05, 30) + np.random.randn(len(days)) * 0.2
glucose_consumed = 5.5 - 3.0 * logistic(days, 1, 0.04, 25) + np.random.randn(len(days)) * 0.15
LDH_release = 0.1 + gaussian_pulse(days, 5, 5, 0.3) + 0.05 * np.exp(0.02 * days) * 0.1

# Clip to valid ranges
for arr in [SOX2, PAX6, NESTIN, TUJ1, MAP2, NEUN, SYN1, PSD95, GFAP, OLIG2]:
    arr[:] = np.clip(arr, 0, 1)

# --- Figure 1: Biomarker temporal profiles ---
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Progenitor markers
ax = axes[0, 0]
ax.plot(days, SOX2, 'b-', linewidth=2, label='SOX2')
ax.plot(days, PAX6, 'r-', linewidth=2, label='PAX6')
ax.plot(days, NESTIN, 'g-', linewidth=2, label='Nestin')
ax.fill_between(days, 0, SOX2, alpha=0.1, color='blue')
ax.set_xlabel('Culture Day')
ax.set_ylabel('Expression Level (normalized)')
ax.set_title('(A) Neural Progenitor Markers')
ax.legend()
ax.grid(True, alpha=0.3)
# Phase annotations
for phase, (d1, d2) in [('Induction', (0, 10)), ('Expansion', (10, 25)),
                          ('Patterning', (25, 45)), ('Maturation', (45, 90))]:
    ax.axvspan(d1, d2, alpha=0.05, color='gray')
    ax.text((d1+d2)/2, 1.05, phase, ha='center', fontsize=7, style='italic')

# Neuronal markers
ax = axes[0, 1]
ax.plot(days, TUJ1, 'm-', linewidth=2, label='TUJ1 (β-III tubulin)')
ax.plot(days, MAP2, 'c-', linewidth=2, label='MAP2')
ax.plot(days, NEUN, 'orange', linewidth=2, label='NeuN')
ax.set_xlabel('Culture Day')
ax.set_ylabel('Expression Level (normalized)')
ax.set_title('(B) Neuronal Maturation Markers')
ax.legend()
ax.grid(True, alpha=0.3)

# Synaptic & Glial markers
ax = axes[1, 0]
ax.plot(days, SYN1, 'r-', linewidth=2, label='Synapsin-1')
ax.plot(days, PSD95, 'b-', linewidth=2, label='PSD-95')
ax.plot(days, GFAP, 'g--', linewidth=2, label='GFAP (astrocytes)')
ax.plot(days, OLIG2, 'm--', linewidth=2, label='OLIG2 (oligodendrocytes)')
ax.set_xlabel('Culture Day')
ax.set_ylabel('Expression Level (normalized)')
ax.set_title('(C) Synaptic & Glial Markers')
ax.legend()
ax.grid(True, alpha=0.3)

# Metabolic indicators
ax = axes[1, 1]
ax2 = ax.twinx()
l1, = ax.plot(days, glucose_consumed, 'b-', linewidth=2, label='Glucose')
l2, = ax.plot(days, lactate, 'r-', linewidth=2, label='Lactate')
l3, = ax2.plot(days, LDH_release, 'k--', linewidth=2, label='LDH (damage)')
ax.set_xlabel('Culture Day')
ax.set_ylabel('Concentration [mM]')
ax2.set_ylabel('LDH Activity [AU]')
ax.set_title('(D) Metabolic Indicators')
lines = [l1, l2, l3]
ax.legend(lines, [l.get_label() for l in lines])
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/biomarker_monitoring.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/biomarker_monitoring.png")

# --- Figure 2: Monitoring strategy ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Sampling schedule
ax = axes[0]
monitoring_methods = {
    'Immunostaining': {'days': [7, 14, 30, 60, 90], 'invasive': True},
    'qPCR': {'days': [7, 14, 21, 30, 45, 60, 75, 90], 'invasive': True},
    'MEA Recording': {'days': list(range(30, 91, 7)), 'invasive': False},
    'Medium Metabolites': {'days': list(range(0, 91, 3)), 'invasive': False},
    'Morphometry': {'days': list(range(0, 91, 1)), 'invasive': False}
}

y_pos = 0
for method, info in monitoring_methods.items():
    color = 'red' if info['invasive'] else 'green'
    ax.scatter(info['days'], [y_pos] * len(info['days']),
              c=color, s=50, zorder=5, alpha=0.7)
    ax.text(-5, y_pos, method, ha='right', va='center', fontsize=9)
    y_pos += 1

ax.set_xlabel('Culture Day')
ax.set_title('(A) Monitoring Schedule')
ax.set_yticks([])
ax.set_xlim(-30, 95)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='green', label='Non-invasive'),
                   Patch(color='red', label='Invasive')], loc='lower right')
ax.grid(True, alpha=0.3, axis='x')

# Composite maturation score
ax = axes[1]
maturation_score = (
    (1 - SOX2) * 0.15 +
    TUJ1 * 0.15 +
    MAP2 * 0.20 +
    NEUN * 0.15 +
    SYN1 * 0.15 +
    PSD95 * 0.10 +
    GFAP * 0.10
)

ax.plot(days, maturation_score, 'k-', linewidth=3, label='Composite Score')
ax.fill_between(days, maturation_score - 0.05, maturation_score + 0.05,
                alpha=0.2, color='blue', label='±1 SD')
# Milestones
milestones = [(10, 'Neural rosettes'), (25, 'Cortical layering'),
              (45, 'Synaptic activity'), (70, 'Network maturation')]
for day, label in milestones:
    idx = day
    ax.annotate(label, (day, maturation_score[idx]),
               xytext=(day + 5, maturation_score[idx] + 0.1),
               arrowprops=dict(arrowstyle='->', color='red'),
               fontsize=8, color='red')
ax.set_xlabel('Culture Day')
ax.set_ylabel('Maturation Score')
ax.set_title('(B) Composite Maturation Score Over Time')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/monitoring_strategy.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/monitoring_strategy.png")

# --- Print summary ---
print("\n=== Biomarker Monitoring Summary ===")
print(f"Day 30 maturation score: {maturation_score[30]:.3f}")
print(f"Day 60 maturation score: {maturation_score[60]:.3f}")
print(f"Day 90 maturation score: {maturation_score[90]:.3f}")
print(f"50% maturation reached at day: {np.argmin(np.abs(maturation_score - 0.5))}")
print(f"80% maturation reached at day: {np.argmin(np.abs(maturation_score - 0.8))}")
