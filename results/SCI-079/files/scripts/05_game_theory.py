"""
Module 5: Pathogen-host coevolution game theory analysis
Gene-for-gene model, replicator dynamics, trench warfare
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

w, c_R, d_a, d, c_a, c_v = 0.8, 0.15, 0.3, 0.7, 0.1, 0.05

host_payoff = np.array([[w-c_R, -c_R], [-d_a, -d]])
pathogen_payoff = np.array([[-w, d_a-c_a], [0, d-c_v]])

def replicator(state, t, A_h, A_p):
    p, q = np.clip(state, 0.001, 0.999)
    f_R = A_h[0,0]*q + A_h[0,1]*(1-q)
    f_r = A_h[1,0]*q + A_h[1,1]*(1-q)
    f_Avr = A_p[0,0]*p + A_p[0,1]*(1-p)
    f_avr = A_p[1,0]*p + A_p[1,1]*(1-p)
    dp = p*(f_R - p*f_R - (1-p)*f_r)
    dq = q*(f_Avr - q*f_Avr - (1-q)*f_avr)
    return [dp, dq]

t = np.linspace(0, 200, 5000)
starts = [(0.1,0.9),(0.9,0.1),(0.5,0.5),(0.3,0.7),(0.7,0.3),(0.2,0.2),(0.8,0.8),(0.4,0.6),(0.6,0.4)]
trajectories = [odeint(replicator, s, t, args=(host_payoff, pathogen_payoff)) for s in starts]

# Nash equilibrium
dq = host_payoff[1,0]-host_payoff[1,1]-host_payoff[0,0]+host_payoff[0,1]
q_star = (host_payoff[0,1]-host_payoff[1,1])/dq if abs(dq)>1e-10 else 0.5
dp = pathogen_payoff[1,0]-pathogen_payoff[1,1]-pathogen_payoff[0,0]+pathogen_payoff[0,1]
p_star = (pathogen_payoff[0,1]-pathogen_payoff[1,1])/dp if abs(dp)>1e-10 else 0.5

# Trench warfare with mutation
def replicator_mut(state, t, A_h, A_p, mu=0.01):
    p, q = np.clip(state, 0.001, 0.999)
    f_R = A_h[0,0]*q + A_h[0,1]*(1-q)
    f_r = A_h[1,0]*q + A_h[1,1]*(1-q)
    f_Avr = A_p[0,0]*p + A_p[0,1]*(1-p)
    f_avr = A_p[1,0]*p + A_p[1,1]*(1-p)
    dp = p*(f_R - p*f_R-(1-p)*f_r) + mu*(0.5-p)
    dq = q*(f_Avr - q*f_Avr-(1-q)*f_avr) + mu*(0.5-q)
    return [dp, dq]

t_long = np.linspace(0, 500, 10000)
sol_mut = odeint(replicator_mut, [0.3, 0.7], t_long, args=(host_payoff, pathogen_payoff, 0.01))

# Multi-locus fitness landscape
def multi_locus_fitness(p1, q1):
    p_def = p1*q1; hf = p_def*w - p1*c_R; return hf

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

ax = axes[0, 0]
for traj in trajectories:
    ax.plot(traj[:,0], traj[:,1], '-', alpha=0.7, lw=1.5)
    ax.plot(traj[0,0], traj[0,1], 'go', markersize=5)
if 0<=p_star<=1 and 0<=q_star<=1:
    ax.plot(p_star, q_star, 'k*', markersize=15, label=f'Eq ({p_star:.2f},{q_star:.2f})')
ax.set_xlabel('Freq(R) Host'); ax.set_ylabel('Freq(Avr) Pathogen')
ax.set_title('Coevolutionary Phase Portrait'); ax.legend(); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.grid(True, alpha=0.3)

ax = axes[0, 1]
sol_main = trajectories[2]
ax.plot(t, sol_main[:,0], 'b-', label='Freq(R)', lw=2)
ax.plot(t, sol_main[:,1], 'r-', label='Freq(Avr)', lw=2)
ax.set_xlabel('Generation'); ax.set_ylabel('Frequency'); ax.set_title('Gene-for-Gene Dynamics'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[0, 2]; ax.axis('off')
txt = f"HOST PAYOFF\n       Avr     avr\nR  [{host_payoff[0,0]:+.2f}  {host_payoff[0,1]:+.2f}]\n"
txt += f"r  [{host_payoff[1,0]:+.2f}  {host_payoff[1,1]:+.2f}]\n\n"
txt += f"PATHOGEN PAYOFF\n       R       r\nAvr [{pathogen_payoff[0,0]:+.2f}  {pathogen_payoff[0,1]:+.2f}]\n"
txt += f"avr [{pathogen_payoff[1,0]:+.2f}  {pathogen_payoff[1,1]:+.2f}]\n\nNash: p*={p_star:.3f}, q*={q_star:.3f}"
ax.text(0.1, 0.5, txt, family='monospace', fontsize=11, va='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow'))
ax.set_title('Payoff Matrices')

ax = axes[1, 0]
ax.plot(t_long, sol_mut[:,0], 'b-', label='Freq(R)', lw=1.5)
ax.plot(t_long, sol_mut[:,1], 'r-', label='Freq(Avr)', lw=1.5)
ax.set_xlabel('Generation'); ax.set_ylabel('Frequency'); ax.set_title('Trench Warfare (μ=0.01)'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1, 1]
pr, qr = np.meshgrid(np.linspace(0,1,50), np.linspace(0,1,50))
Hfit = np.vectorize(multi_locus_fitness)(pr, qr)
im = ax.contourf(pr, qr, Hfit, levels=20, cmap='RdYlGn')
plt.colorbar(im, ax=ax, label='Host Fitness')
ax.set_xlabel('Freq(R)'); ax.set_ylabel('Freq(Avr)'); ax.set_title('Host Fitness Landscape')

ax = axes[1, 2]
costs = np.linspace(0.01, 0.4, 30)
eq_Avr = []
for c in costs:
    hp = np.array([[w-c, -c], [-d_a, -d]])
    dq2 = hp[1,0]-hp[1,1]-hp[0,0]+hp[0,1]
    eq_Avr.append(np.clip((hp[0,1]-hp[1,1])/dq2, 0, 1) if abs(dq2)>1e-10 else 0.5)
ax.plot(costs, eq_Avr, 'r-s', label='Eq Freq(Avr)', markersize=4, lw=2)
ax.set_xlabel('Cost of R gene'); ax.set_ylabel('Equilibrium Frequency')
ax.set_title('R Gene Cost vs Equilibrium'); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/05_game_theory.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/05_game_theory.svg', bbox_inches='tight')
plt.close()

results = {
    'host_payoff': host_payoff.tolist(), 'pathogen_payoff': pathogen_payoff.tolist(),
    'nash_equilibrium': {'p_star_R': float(p_star), 'q_star_Avr': float(q_star)},
    'parameters': {'w': w, 'c_R': c_R, 'd_a': d_a, 'd': d, 'c_a': c_a, 'c_v': c_v},
    'trench_warfare_final': {'R': float(sol_mut[-1,0]), 'Avr': float(sol_mut[-1,1])}
}
with open('results/05_game_theory.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Module 5 completed.")
print(f"Nash: p*(R)={p_star:.4f}, q*(Avr)={q_star:.4f}")
