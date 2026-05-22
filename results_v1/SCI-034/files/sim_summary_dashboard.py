"""
Summary figure: QKD Network Overview Dashboard
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json

np.random.seed(42)

# Load results
with open('results/bb84_e91_results.json') as f:
    r_bb84 = json.load(f)
with open('results/quantum_repeater_results.json') as f:
    r_rep = json.load(f)
with open('results/distillation_results.json') as f:
    r_dist = json.load(f)
with open('results/routing_results.json') as f:
    r_route = json.load(f)
with open('results/decoherence_results.json') as f:
    r_deco = json.load(f)
with open('results/tokyo_casestudy_results.json') as f:
    r_tokyo = json.load(f)

fig = plt.figure(figsize=(18, 12))
fig.patch.set_facecolor('#f8f9fa')

gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.5, wspace=0.4)

# Title
fig.text(0.5, 0.97, 'Quantum Internet: QKD Network Protocol Design — Summary Dashboard',
         ha='center', va='top', fontsize=16, fontweight='bold', color='#2c3e50')
fig.text(0.5, 0.94, 'Tokyo QKD Network Scale Case Study | BB84/E91 · Quantum Repeaters · Entanglement Distillation · Fidelity-Aware Routing',
         ha='center', va='top', fontsize=10, color='#7f8c8d')

# ─── Panel 1: BB84 key rate vs block size ────────────────────────────────────
ax1 = fig.add_subplot(gs[0, :2])
n_vals = np.logspace(3, 8, 60).astype(int)

def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

def bb84_rate(n, q, eps=1e-8):
    leak_ec = n * h2(q) + np.log2(1/eps) * np.sqrt(n)
    delta = np.sqrt((2 * np.log(1/eps) + np.log(2*n+1)) / (2*n))
    q_p = min(q + delta, 0.5)
    l = n * (1 - h2(q_p)) - leak_ec - 2*np.log2(1/eps)
    return max(l/n, 0)

qbers = [0.01, 0.03, 0.05, 0.08, 0.11]
colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(qbers)))
for q, c in zip(qbers, colors):
    rates = [bb84_rate(n, q) for n in n_vals]
    ax1.semilogx(n_vals, rates, color=c, lw=2, label=f'QBER={q}')
    ax1.axhline(1-h2(q), color=c, ls='--', alpha=0.35, lw=1)

ax1.set_xlabel('Block size n', fontsize=10)
ax1.set_ylabel('Key rate (bits/sifted bit)', fontsize=10)
ax1.set_title('BB84 Finite-Key Rate vs Block Size\n(dashed = asymptotic limit, ε=10⁻⁸)', fontsize=11)
ax1.legend(fontsize=8, ncol=2)
ax1.set_ylim(0, 1.05)
ax1.grid(True, alpha=0.25)
ax1.set_facecolor('#fdfefe')

# ─── Panel 2: Repeater rate vs segments ─────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2])
chain = r_rep['chain_analysis']
n_segs = [c['n_segments'] for c in chain]
rates_r = [c['chain_rate_Hz'] for c in chain]
fids_r = [c['fidelity_estimate'] for c in chain]

ax2_twin = ax2.twinx()
ax2.semilogy(n_segs, rates_r, 'bo-', lw=2, markersize=7, label='Rate (Hz)')
ax2_twin.plot(n_segs, fids_r, 'r^--', lw=2, markersize=7, label='Fidelity')
ax2.set_xlabel('Repeater segments', fontsize=10)
ax2.set_ylabel('E2E rate (Hz)', fontsize=10, color='blue')
ax2_twin.set_ylabel('Fidelity', fontsize=10, color='red')
ax2.set_title('Repeater Chain\n(250 km)', fontsize=11)
ax2.tick_params(axis='y', labelcolor='blue')
ax2_twin.tick_params(axis='y', labelcolor='red')
ax2.grid(True, alpha=0.25)
ax2.set_facecolor('#fdfefe')

# ─── Panel 3: Memory requirements ───────────────────────────────────────────
ax3 = fig.add_subplot(gs[0, 3])
mem_reqs = [c['min_memory_per_node'] for c in chain]
t_cohs = [c['t_coh_required_ms'] for c in chain]
ax3.bar(range(len(n_segs)), np.log10([max(m, 1) for m in mem_reqs]),
        color=plt.cm.viridis(np.linspace(0.2, 0.9, len(n_segs))), alpha=0.85)
ax3.set_xticks(range(len(n_segs)))
ax3.set_xticklabels([str(n) for n in n_segs], fontsize=8)
ax3.set_xlabel('Segments', fontsize=10)
ax3.set_ylabel('log₁₀(Min memories/node)', fontsize=9)
ax3.set_title('Memory per Node\n(log scale)', fontsize=11)
ax3.grid(True, alpha=0.25, axis='y')
ax3.set_facecolor('#fdfefe')

# ─── Panel 4: Distillation convergence ──────────────────────────────────────
ax4 = fig.add_subplot(gs[1, :2])

def dejmps_step(F):
    p = F**2 + ((1-F)/3)**2 * 1 + 2*F*(1-F)/3
    p_num = F**2
    q_num = ((1-F)/3)**2
    p_suc = p_num + q_num + 2*F*(1-F)/3
    F_out = (p_num + q_num) / p_suc
    return F_out, p_suc

F0_vals_plot = [0.55, 0.65, 0.75, 0.85, 0.92]
c_plot = plt.cm.cool(np.linspace(0.1, 0.9, len(F0_vals_plot)))
for F0, c in zip(F0_vals_plot, c_plot):
    fids = [F0]
    for _ in range(12):
        if fids[-1] >= 0.999:
            break
        Fn, _ = dejmps_step(fids[-1])
        if Fn <= fids[-1]:
            break
        fids.append(Fn)
    ax4.plot(range(len(fids)), fids, 'o-', color=c, lw=2, markersize=6, label=f'F₀={F0}')

ax4.axhline(0.99, color='red', ls='--', lw=1.5, alpha=0.7, label='Target F=0.99')
ax4.axhline(2/3, color='gray', ls=':', lw=1.5, label='Distillation threshold')
ax4.set_xlabel('Distillation round', fontsize=10)
ax4.set_ylabel('Fidelity', fontsize=10)
ax4.set_title('DEJMPS Entanglement Distillation Convergence', fontsize=11)
ax4.legend(fontsize=8, ncol=2)
ax4.set_ylim(0.5, 1.02)
ax4.grid(True, alpha=0.25)
ax4.set_facecolor('#fdfefe')

# ─── Panel 5: Channel model ──────────────────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 2:])
distances = np.linspace(1, 200, 200)

def fiber_model(L):
    eta = 10**(-0.2*L/10) * 0.85
    R_s = 1e6 * 0.1 * eta
    R_d = 100
    qber = R_d/(R_s+R_d) + 0.01
    def h(p): return 0.0 if p<=0 or p>=1 else -p*np.log2(p)-(1-p)*np.log2(1-p)
    rate = max(1e3*0.1*eta*(1-h(qber)-h(qber)), 0)
    return eta, qber, rate

etas, qbers_ch, rates_ch = zip(*[fiber_model(L) for L in distances])
etas, qbers_ch, rates_ch = np.array(etas), np.array(qbers_ch)*100, np.array(rates_ch)

ax5_t = ax5.twinx()
l1, = ax5.semilogy(distances, np.maximum(rates_ch, 1e-4), 'b-', lw=2.5, label='Key rate (kHz)')
l2, = ax5_t.plot(distances, qbers_ch, 'r--', lw=2, label='QBER (%)')
ax5_t.axhline(11, color='orange', ls=':', lw=2, label='QBER=11% limit')
ax5.set_xlabel('Fiber distance (km)', fontsize=10)
ax5.set_ylabel('Secret key rate (kHz)', fontsize=10, color='blue')
ax5_t.set_ylabel('QBER (%)', fontsize=10, color='red')
ax5.set_title('BB84 Channel Performance vs Distance\n(1 GHz source, 85% detector, 100 dark counts/s)', fontsize=11)
ax5.tick_params(axis='y', labelcolor='blue')
ax5_t.tick_params(axis='y', labelcolor='red')
lines = [l1, l2]
ax5.legend(lines, [l.get_label() for l in lines], fontsize=9)
ax5.grid(True, alpha=0.25)
ax5.set_facecolor('#fdfefe')

# ─── Panel 6: Tokyo network summary table ────────────────────────────────────
ax6 = fig.add_subplot(gs[2, :])
ax6.axis('off')

table_data = [
    ['Metric', 'Value', 'Protocol/Model', 'Notes'],
    ['BB84 finite key rate (n=10⁶, QBER=3%, ε=10⁻⁸)', '55.98%', 'Scarani-Renner bound', 'vs 80.6% asymptotic'],
    ['BB84 minimum block size (QBER=3%)', '7,038 bits', 'Finite-key analysis', 'for positive key rate >1%'],
    ['E91 key rate (pure Bell, n=10⁶)', '97.0%', 'Acín et al. DI bound', 'CHSH S=2√2'],
    ['Optimal repeater segments (250 km)', '16 segments', 'Sangouard model', 'rate=28.5 Hz, F=0.67'],
    ['Memory requirement (16 segments)', '8 qubits/node', 'Geometric model', 'T_coh ≥ 0.1 ms needed'],
    ['DEJMPS rounds to F=0.99 (F₀=0.75)', '9 rounds', 'DEJMPS recurrence', '1575× pair overhead'],
    ['Max QKD distance (α=0.2 dB/km)', '98.7 km', 'BB84+dark counts', '100 dark counts/s, μ=0.1'],
    ['Tokyo network mean E2E fidelity', '37.3% ± 20.1%', 'Fidelity-optimal routing', '10 nodes, 15 links, 279 km'],
    ['All 15 links Eve-detectable', '100%', 'Intercept-resend model', 'QBER jump >11% threshold'],
    ['Bell pair fidelity (rare-earth, 1 ms)', '97.3%', 'T₁=100 ms, T₂=10 ms', 'Suitable for QR memory'],
]

table = ax6.table(cellText=table_data[1:], colLabels=table_data[0],
                   cellLoc='left', loc='center',
                   colWidths=[0.35, 0.15, 0.22, 0.28])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.6)

# Style header
for j in range(4):
    table[0, j].set_facecolor('#2c3e50')
    table[0, j].set_text_props(color='white', fontweight='bold')

for i in range(1, len(table_data)):
    for j in range(4):
        bg = '#eaf2ff' if i % 2 == 0 else 'white'
        table[i, j].set_facecolor(bg)

ax6.set_title('Key Results Summary', fontsize=12, fontweight='bold', pad=10)

plt.savefig('figures/qkd_network_summary_dashboard.png', dpi=200, bbox_inches='tight',
            facecolor='#f8f9fa')
plt.savefig('figures/qkd_network_summary_dashboard.svg', bbox_inches='tight',
            facecolor='#f8f9fa')
plt.close()

print("Summary dashboard saved.")
