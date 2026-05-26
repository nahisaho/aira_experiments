#!/usr/bin/env python3
"""
Plant Immunity Signaling Model: PTI and ETI Pathway Simulations
Integrates receptor binding, MAPK cascade, SA/JA crosstalk,
WRKY/TGA network, coevolution game theory, and rice blast case study.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import networkx as nx
import os
import warnings
warnings.filterwarnings('ignore')

FIGDIR = 'figures'
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 150
})

# ============================================================
# 1. Receptor-Level Ligand Binding & Signal Initiation Model
# ============================================================
def run_receptor_model():
    print("=== 1. Receptor-Ligand Binding Model ===")

    # Parameters
    kon = 0.1      # association rate (nM^-1 min^-1)
    koff = 0.01    # dissociation rate (min^-1)
    kact = 0.05    # receptor activation rate (min^-1)
    kdeact = 0.02  # receptor deactivation rate (min^-1)
    kint = 0.005   # receptor internalization rate (min^-1)
    ksynth = 0.1   # receptor synthesis rate (nM/min)
    kdeg = 0.001   # receptor degradation rate (min^-1)

    # [R_free, R_bound, R_active, R_internal, Signal]
    def odes(t, y, PAMP_conc):
        Rf, Rb, Ra, Ri, S = y
        binding = kon * Rf * PAMP_conc - koff * Rb
        activation = kact * Rb - kdeact * Ra
        internalization = kint * Ra
        synthesis = ksynth - kdeg * Rf

        dRf = -binding + synthesis
        dRb = binding - activation
        dRa = activation - internalization
        dRi = internalization - 0.01 * Ri
        dS = 0.1 * Ra - 0.02 * S
        return [dRf, dRb, dRa, dRi, dS]

    t_span = (0, 120)
    t_eval = np.linspace(0, 120, 500)
    y0 = [10.0, 0, 0, 0, 0]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) Time course at different PAMP concentrations
    pamp_concs = [0.5, 1.0, 5.0, 10.0, 50.0]
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(pamp_concs)))
    for pc, c in zip(pamp_concs, colors):
        sol = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(pc,), method='RK45')
        axes[0, 0].plot(sol.t, sol.y[4], color=c, label=f'PAMP={pc} nM')
    axes[0, 0].set_xlabel('Time (min)')
    axes[0, 0].set_ylabel('Signal Intensity (a.u.)')
    axes[0, 0].set_title('(a) Signal Output vs PAMP Concentration')
    axes[0, 0].legend(fontsize=8)

    # (b) Receptor state dynamics at PAMP=5 nM
    sol = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(5.0,), method='RK45')
    labels = ['Free Receptor', 'Bound Receptor', 'Active Receptor', 'Internalized']
    for i, lbl in enumerate(labels):
        axes[0, 1].plot(sol.t, sol.y[i], label=lbl)
    axes[0, 1].set_xlabel('Time (min)')
    axes[0, 1].set_ylabel('Concentration (nM)')
    axes[0, 1].set_title('(b) Receptor State Dynamics (PAMP=5 nM)')
    axes[0, 1].legend(fontsize=8)

    # (c) Dose-response curve
    pamp_range = np.logspace(-2, 2, 50)
    max_signals = []
    ec50_signals = []
    for pc in pamp_range:
        sol = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(pc,), method='RK45')
        max_signals.append(np.max(sol.y[4]))
    max_signals = np.array(max_signals)
    axes[1, 0].semilogx(pamp_range, max_signals / max_signals.max(), 'b-', linewidth=2)
    axes[1, 0].axhline(0.5, color='r', linestyle='--', alpha=0.5, label='EC50')
    axes[1, 0].set_xlabel('PAMP Concentration (nM)')
    axes[1, 0].set_ylabel('Normalized Max Signal')
    axes[1, 0].set_title('(c) Dose-Response Curve')
    axes[1, 0].legend()

    # (d) PTI vs ETI signal comparison
    # ETI: stronger, sustained signal via NLR
    def odes_eti(t, y):
        Rf, Rb, Ra, Ri, S = y
        eff_conc = 2.0
        binding = 0.2 * Rf * eff_conc - 0.005 * Rb
        activation = 0.1 * Rb - 0.01 * Ra
        internalization = 0.001 * Ra
        synthesis = ksynth - kdeg * Rf
        dRf = -binding + synthesis
        dRb = binding - activation
        dRa = activation - internalization
        dRi = internalization - 0.01 * Ri
        dS = 0.2 * Ra - 0.01 * S
        return [dRf, dRb, dRa, dRi, dS]

    sol_pti = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(5.0,), method='RK45')
    sol_eti = solve_ivp(odes_eti, t_span, y0, t_eval=t_eval, method='RK45')
    axes[1, 1].plot(sol_pti.t, sol_pti.y[4], label='PTI (FLS2-flg22)', color='blue')
    axes[1, 1].plot(sol_eti.t, sol_eti.y[4], label='ETI (NLR-Effector)', color='red')
    axes[1, 1].set_xlabel('Time (min)')
    axes[1, 1].set_ylabel('Signal Intensity (a.u.)')
    axes[1, 1].set_title('(d) PTI vs ETI Signal Dynamics')
    axes[1, 1].legend()

    plt.suptitle('Figure 1: Receptor-Level Ligand Binding and Signal Initiation', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig1_receptor_binding.png')
    plt.close()
    print("  Saved fig1_receptor_binding.png")

    return {
        'ec50_approx': pamp_range[np.argmin(np.abs(max_signals/max_signals.max() - 0.5))],
        'pti_max_signal': np.max(sol_pti.y[4]),
        'eti_max_signal': np.max(sol_eti.y[4]),
    }


# ============================================================
# 2. MAPK Cascade Dynamics Simulation
# ============================================================
def run_mapk_cascade():
    print("=== 2. MAPK Cascade Dynamics ===")

    # 3-tier MAPK cascade: MAPKKK -> MAPKK -> MAPK
    # Each with phosphorylation/dephosphorylation
    def mapk_odes(t, y, params):
        M3, M3p, M2, M2p, M1, M1p = y
        S = params['signal']
        # Stimulus function
        stim = S * np.exp(-0.01 * t) if params.get('transient', False) else S

        k1, k2 = params['k1'], params['k2']   # MAPKKK activation/deactivation
        k3, k4 = params['k3'], params['k4']   # MAPKK
        k5, k6 = params['k5'], params['k6']   # MAPK
        Km = params['Km']

        # Michaelis-Menten kinetics
        v1 = k1 * stim * M3 / (Km + M3)
        v2 = k2 * M3p / (Km + M3p)
        v3 = k3 * M3p * M2 / (Km + M2)
        v4 = k4 * M2p / (Km + M2p)
        v5 = k5 * M2p * M1 / (Km + M1)
        v6 = k6 * M1p / (Km + M1p)

        dM3 = -v1 + v2
        dM3p = v1 - v2
        dM2 = -v3 + v4
        dM2p = v3 - v4
        dM1 = -v5 + v6
        dM1p = v5 - v6

        return [dM3, dM3p, dM2, dM2p, dM1, dM1p]

    base_params = {
        'signal': 1.0, 'k1': 0.5, 'k2': 0.1, 'k3': 0.3, 'k4': 0.1,
        'k5': 0.3, 'k6': 0.1, 'Km': 0.5, 'transient': False
    }

    t_span = (0, 100)
    t_eval = np.linspace(0, 100, 500)
    y0 = [1.0, 0, 1.0, 0, 1.0, 0]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) Cascade activation dynamics
    sol = solve_ivp(mapk_odes, t_span, y0, t_eval=t_eval, args=(base_params,), method='RK45')
    axes[0, 0].plot(sol.t, sol.y[1], label='MAPKKK*', linewidth=2)
    axes[0, 0].plot(sol.t, sol.y[3], label='MAPKK*', linewidth=2)
    axes[0, 0].plot(sol.t, sol.y[5], label='MAPK*', linewidth=2)
    axes[0, 0].set_xlabel('Time (min)')
    axes[0, 0].set_ylabel('Active Fraction')
    axes[0, 0].set_title('(a) MAPK Cascade Activation Dynamics')
    axes[0, 0].legend()

    # (b) Signal amplification (ultrasensitivity)
    signals = np.logspace(-2, 1, 30)
    ss_mapk = []
    for s in signals:
        p = base_params.copy()
        p['signal'] = s
        sol = solve_ivp(mapk_odes, (0, 200), y0, t_eval=[200], args=(p,), method='RK45')
        ss_mapk.append(sol.y[5, -1])
    ss_mapk = np.array(ss_mapk)
    axes[0, 1].semilogx(signals, ss_mapk / (ss_mapk.max() + 1e-10), 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Input Signal')
    axes[0, 1].set_ylabel('Normalized MAPK* (Steady State)')
    axes[0, 1].set_title('(b) Ultrasensitivity (Hill-like Response)')

    # Calculate Hill coefficient
    norm_mapk = ss_mapk / (ss_mapk.max() + 1e-10)
    try:
        idx10 = np.argmin(np.abs(norm_mapk - 0.1))
        idx90 = np.argmin(np.abs(norm_mapk - 0.9))
        hill_coeff = np.log(81) / np.log(signals[idx90] / signals[idx10])
    except:
        hill_coeff = 'N/A'
    axes[0, 1].text(0.05, 0.9, f'Hill coeff ≈ {hill_coeff:.2f}', transform=axes[0, 1].transAxes)

    # (c) Transient vs sustained signal
    p_trans = base_params.copy()
    p_trans['transient'] = True
    p_sust = base_params.copy()
    p_sust['transient'] = False
    sol_t = solve_ivp(mapk_odes, t_span, y0, t_eval=t_eval, args=(p_trans,), method='RK45')
    sol_s = solve_ivp(mapk_odes, t_span, y0, t_eval=t_eval, args=(p_sust,), method='RK45')
    axes[1, 0].plot(sol_t.t, sol_t.y[5], label='Transient (PTI-like)', linewidth=2)
    axes[1, 0].plot(sol_s.t, sol_s.y[5], label='Sustained (ETI-like)', linewidth=2)
    axes[1, 0].set_xlabel('Time (min)')
    axes[1, 0].set_ylabel('Active MAPK')
    axes[1, 0].set_title('(c) Transient vs Sustained MAPK Activation')
    axes[1, 0].legend()

    # (d) Sensitivity analysis
    param_names = ['k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'Km']
    sensitivities = []
    base_output = ss_mapk[-1]
    for pn in param_names:
        p = base_params.copy()
        p[pn] = base_params[pn] * 1.1
        sol = solve_ivp(mapk_odes, (0, 200), y0, t_eval=[200], args=(p,), method='RK45')
        delta = (sol.y[5, -1] - base_output) / (0.1 * base_params[pn] + 1e-10) * base_params[pn] / (base_output + 1e-10)
        sensitivities.append(delta)

    colors_sens = ['red' if s > 0 else 'blue' for s in sensitivities]
    axes[1, 1].barh(param_names, sensitivities, color=colors_sens, alpha=0.7)
    axes[1, 1].set_xlabel('Sensitivity Coefficient')
    axes[1, 1].set_title('(d) Parameter Sensitivity Analysis')
    axes[1, 1].axvline(0, color='black', linewidth=0.5)

    plt.suptitle('Figure 2: MAPK Cascade Dynamics Simulation', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig2_mapk_cascade.png')
    plt.close()
    print("  Saved fig2_mapk_cascade.png")

    return {'hill_coefficient': hill_coeff, 'sensitivities': dict(zip(param_names, sensitivities))}


# ============================================================
# 3. SA/JA Crosstalk Model
# ============================================================
def run_sa_ja_crosstalk():
    print("=== 3. SA/JA Crosstalk Model ===")

    def crosstalk_odes(t, y, params):
        SA, JA, NPR1, JAZ, PR, PDF = y
        pathogen_signal = params.get('pathogen', 0)
        herbivore_signal = params.get('herbivore', 0)

        # SA synthesis (induced by biotrophic pathogen)
        sa_synth = params['v_sa'] * pathogen_signal / (params['K_sa'] + pathogen_signal)
        sa_deg = params['d_sa'] * SA

        # JA synthesis (induced by necrotrophic pathogen / herbivore)
        ja_synth = params['v_ja'] * herbivore_signal / (params['K_ja'] + herbivore_signal)
        ja_deg = params['d_ja'] * JA

        # NPR1 activation by SA
        npr1_act = params['k_npr1'] * SA / (params['Kn'] + SA) - params['d_npr1'] * NPR1

        # JAZ degradation by JA (SCF^COI1 pathway)
        jaz_prod = params['v_jaz'] - params['k_jaz'] * JA * JAZ / (params['Kj'] + JAZ)
        jaz_deg = params['d_jaz'] * JAZ

        # SA-JA antagonism
        sa_inhibit_ja = params['alpha'] * SA / (params['Ka'] + SA)  # SA suppresses JA
        ja_inhibit_sa = params['beta'] * JA / (params['Kb'] + JA)   # JA suppresses SA

        # PR gene expression (SA-dependent, via NPR1)
        pr_expr = params['v_pr'] * NPR1 / (params['Kp'] + NPR1) - params['d_pr'] * PR
        # PDF1.2 expression (JA-dependent, suppressed by SA)
        pdf_expr = params['v_pdf'] * (1 - JAZ / (params['Kz'] + JAZ)) * (1 - sa_inhibit_ja) - params['d_pdf'] * PDF

        dSA = sa_synth - sa_deg - ja_inhibit_sa * SA
        dJA = ja_synth - ja_deg - sa_inhibit_ja * JA
        dNPR1 = npr1_act
        dJAZ = jaz_prod - jaz_deg
        dPR = pr_expr
        dPDF = pdf_expr

        return [dSA, dJA, dNPR1, dJAZ, dPR, dPDF]

    base_params = {
        'v_sa': 5.0, 'K_sa': 1.0, 'd_sa': 0.1,
        'v_ja': 5.0, 'K_ja': 1.0, 'd_ja': 0.1,
        'k_npr1': 0.5, 'Kn': 1.0, 'd_npr1': 0.05,
        'v_jaz': 0.5, 'k_jaz': 0.3, 'Kj': 0.5, 'd_jaz': 0.05,
        'alpha': 0.8, 'Ka': 2.0,
        'beta': 0.3, 'Kb': 2.0,
        'v_pr': 2.0, 'Kp': 0.5, 'd_pr': 0.05,
        'v_pdf': 2.0, 'Kz': 0.5, 'd_pdf': 0.05,
    }

    t_span = (0, 200)
    t_eval = np.linspace(0, 200, 1000)
    y0 = [0.1, 0.1, 0, 1.0, 0, 0]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # (a) Biotrophic pathogen (SA pathway)
    p = base_params.copy()
    p['pathogen'] = 5.0
    p['herbivore'] = 0
    sol = solve_ivp(crosstalk_odes, t_span, y0, t_eval=t_eval, args=(p,), method='RK45')
    axes[0, 0].plot(sol.t, sol.y[0], 'r-', label='SA', linewidth=2)
    axes[0, 0].plot(sol.t, sol.y[1], 'b-', label='JA', linewidth=2)
    axes[0, 0].set_title('(a) Biotroph Attack: SA/JA')
    axes[0, 0].set_xlabel('Time (min)')
    axes[0, 0].set_ylabel('Concentration (a.u.)')
    axes[0, 0].legend()

    # (b) Necrotrophic / herbivore (JA pathway)
    p2 = base_params.copy()
    p2['pathogen'] = 0
    p2['herbivore'] = 5.0
    sol2 = solve_ivp(crosstalk_odes, t_span, y0, t_eval=t_eval, args=(p2,), method='RK45')
    axes[0, 1].plot(sol2.t, sol2.y[0], 'r-', label='SA', linewidth=2)
    axes[0, 1].plot(sol2.t, sol2.y[1], 'b-', label='JA', linewidth=2)
    axes[0, 1].set_title('(b) Necrotroph Attack: SA/JA')
    axes[0, 1].set_xlabel('Time (min)')
    axes[0, 1].legend()

    # (c) Combined attack
    p3 = base_params.copy()
    p3['pathogen'] = 5.0
    p3['herbivore'] = 5.0
    sol3 = solve_ivp(crosstalk_odes, t_span, y0, t_eval=t_eval, args=(p3,), method='RK45')
    axes[0, 2].plot(sol3.t, sol3.y[0], 'r-', label='SA', linewidth=2)
    axes[0, 2].plot(sol3.t, sol3.y[1], 'b-', label='JA', linewidth=2)
    axes[0, 2].set_title('(c) Combined Attack: SA/JA')
    axes[0, 2].set_xlabel('Time (min)')
    axes[0, 2].legend()

    # (d) Defense gene expression comparison
    axes[1, 0].plot(sol.t, sol.y[4], 'r-', label='PR1 (Biotroph)', linewidth=2)
    axes[1, 0].plot(sol2.t, sol2.y[4], 'r--', label='PR1 (Necrotroph)', linewidth=1)
    axes[1, 0].plot(sol.t, sol.y[5], 'b-', label='PDF1.2 (Biotroph)', linewidth=2)
    axes[1, 0].plot(sol2.t, sol2.y[5], 'b--', label='PDF1.2 (Necrotroph)', linewidth=1)
    axes[1, 0].set_title('(d) Defense Gene Expression')
    axes[1, 0].set_xlabel('Time (min)')
    axes[1, 0].set_ylabel('Expression Level')
    axes[1, 0].legend(fontsize=7)

    # (e) NPR1 and JAZ dynamics
    axes[1, 1].plot(sol.t, sol.y[2], 'g-', label='NPR1 (Biotroph)', linewidth=2)
    axes[1, 1].plot(sol2.t, sol2.y[2], 'g--', label='NPR1 (Necrotroph)', linewidth=1)
    axes[1, 1].plot(sol.t, sol.y[3], 'm-', label='JAZ (Biotroph)', linewidth=2)
    axes[1, 1].plot(sol2.t, sol2.y[3], 'm--', label='JAZ (Necrotroph)', linewidth=1)
    axes[1, 1].set_title('(e) NPR1 & JAZ Dynamics')
    axes[1, 1].set_xlabel('Time (min)')
    axes[1, 1].legend(fontsize=7)

    # (f) Phase portrait SA vs JA
    # Different alpha values
    for alpha_val, color, lbl in [(0.2, 'green', 'Weak ant.'), (0.8, 'orange', 'Moderate'), (1.5, 'red', 'Strong ant.')]:
        p4 = base_params.copy()
        p4['pathogen'] = 5.0
        p4['herbivore'] = 5.0
        p4['alpha'] = alpha_val
        sol4 = solve_ivp(crosstalk_odes, t_span, y0, t_eval=t_eval, args=(p4,), method='RK45')
        axes[1, 2].plot(sol4.y[0], sol4.y[1], color=color, label=f'α={alpha_val} ({lbl})', linewidth=1.5)
        axes[1, 2].plot(sol4.y[0, -1], sol4.y[1, -1], 'o', color=color, markersize=8)
    axes[1, 2].set_xlabel('SA Concentration')
    axes[1, 2].set_ylabel('JA Concentration')
    axes[1, 2].set_title('(f) SA-JA Phase Portrait')
    axes[1, 2].legend(fontsize=7)

    plt.suptitle('Figure 3: Salicylic Acid / Jasmonic Acid Pathway Crosstalk', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig3_sa_ja_crosstalk.png')
    plt.close()
    print("  Saved fig3_sa_ja_crosstalk.png")

    return {
        'biotroph_SA_ss': sol.y[0, -1],
        'biotroph_JA_ss': sol.y[1, -1],
        'necrotroph_SA_ss': sol2.y[0, -1],
        'necrotroph_JA_ss': sol2.y[1, -1],
    }


# ============================================================
# 4. Transcription Factor Network (WRKY/TGA)
# ============================================================
def run_tf_network():
    print("=== 4. WRKY/TGA Transcription Factor Network ===")

    # Build regulatory network
    G = nx.DiGraph()

    # Nodes: TFs and target genes
    tfs = ['WRKY18', 'WRKY33', 'WRKY40', 'WRKY53', 'WRKY70',
           'TGA2', 'TGA3', 'TGA5', 'NPR1', 'MYC2']
    targets = ['PR1', 'PR2', 'PR5', 'PDF1.2', 'VSP2', 'LOX2',
               'PAD4', 'EDS1', 'ICS1', 'FMO1']

    for n in tfs:
        G.add_node(n, node_type='TF')
    for n in targets:
        G.add_node(n, node_type='target')

    # Regulatory edges (activation=+1, repression=-1)
    edges = [
        ('NPR1', 'TGA2', 1), ('NPR1', 'TGA3', 1), ('NPR1', 'TGA5', 1),
        ('TGA2', 'PR1', 1), ('TGA3', 'PR1', 1), ('TGA5', 'PR2', 1),
        ('WRKY18', 'PR1', -1), ('WRKY18', 'WRKY40', 1),
        ('WRKY33', 'PAD4', 1), ('WRKY33', 'PDF1.2', 1),
        ('WRKY40', 'WRKY33', -1), ('WRKY40', 'PR1', -1),
        ('WRKY53', 'PR1', 1), ('WRKY53', 'PAD4', 1),
        ('WRKY70', 'PR1', 1), ('WRKY70', 'PR2', 1), ('WRKY70', 'PDF1.2', -1),
        ('WRKY70', 'VSP2', -1),
        ('MYC2', 'VSP2', 1), ('MYC2', 'LOX2', 1), ('MYC2', 'PDF1.2', -1),
        ('TGA2', 'WRKY70', 1),
        ('WRKY53', 'EDS1', 1), ('WRKY33', 'ICS1', -1),
        ('NPR1', 'WRKY70', 1), ('TGA2', 'FMO1', 1),
    ]

    for u, v, w in edges:
        G.add_edge(u, v, weight=w)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # (a) Network visualization
    ax = axes[0, 0]
    pos = nx.spring_layout(G, seed=42, k=2)
    tf_nodes = [n for n in G.nodes if G.nodes[n].get('node_type') == 'TF']
    tgt_nodes = [n for n in G.nodes if G.nodes[n].get('node_type') == 'target']

    nx.draw_networkx_nodes(G, pos, nodelist=tf_nodes, node_color='lightcoral',
                           node_size=600, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=tgt_nodes, node_color='lightblue',
                           node_size=400, ax=ax)
    act_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] > 0]
    rep_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] < 0]
    nx.draw_networkx_edges(G, pos, edgelist=act_edges, edge_color='green',
                           arrows=True, ax=ax, alpha=0.7, width=1.5)
    nx.draw_networkx_edges(G, pos, edgelist=rep_edges, edge_color='red',
                           arrows=True, ax=ax, alpha=0.7, width=1.5, style='dashed')
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
    ax.set_title('(a) WRKY/TGA Regulatory Network')

    # (b) Network centrality
    ax = axes[0, 1]
    betw = nx.betweenness_centrality(G)
    deg = dict(G.degree())
    nodes_sorted = sorted(G.nodes, key=lambda x: betw[x], reverse=True)[:12]
    x_pos = np.arange(len(nodes_sorted))
    ax.bar(x_pos - 0.15, [betw[n] for n in nodes_sorted], 0.3, label='Betweenness', color='coral')
    ax.bar(x_pos + 0.15, [deg[n] / max(deg.values()) for n in nodes_sorted], 0.3,
           label='Degree (norm)', color='skyblue')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(nodes_sorted, rotation=45, ha='right', fontsize=7)
    ax.set_title('(b) Network Centrality Measures')
    ax.legend(fontsize=8)

    # (c) Boolean simulation of defense gene expression
    ax = axes[1, 0]
    # Simplified Boolean model: SA-on, JA-off scenario vs SA-off, JA-on
    scenarios = {
        'SA-pathway ON': {'NPR1': 1, 'TGA2': 1, 'TGA3': 1, 'TGA5': 1,
                          'WRKY70': 1, 'WRKY53': 1, 'WRKY33': 0, 'WRKY18': 0,
                          'WRKY40': 0, 'MYC2': 0},
        'JA-pathway ON': {'NPR1': 0, 'TGA2': 0, 'TGA3': 0, 'TGA5': 0,
                          'WRKY70': 0, 'WRKY53': 0, 'WRKY33': 1, 'WRKY18': 1,
                          'WRKY40': 1, 'MYC2': 1},
        'Both ON':       {'NPR1': 1, 'TGA2': 1, 'TGA3': 1, 'TGA5': 1,
                          'WRKY70': 1, 'WRKY53': 1, 'WRKY33': 1, 'WRKY18': 1,
                          'WRKY40': 1, 'MYC2': 1},
    }

    gene_activity = {}
    for scenario_name, tf_states in scenarios.items():
        gene_scores = {}
        for target in targets:
            score = 0
            for tf in tfs:
                if G.has_edge(tf, target):
                    w = G[tf][target]['weight']
                    score += w * tf_states.get(tf, 0)
            gene_scores[target] = max(0, score)
        gene_activity[scenario_name] = gene_scores

    x_pos = np.arange(len(targets))
    width = 0.25
    for i, (sname, scores) in enumerate(gene_activity.items()):
        ax.bar(x_pos + i * width, [scores[g] for g in targets], width, label=sname)
    ax.set_xticks(x_pos + width)
    ax.set_xticklabels(targets, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Predicted Activity Score')
    ax.set_title('(c) Boolean Network: Defense Gene Prediction')
    ax.legend(fontsize=7)

    # (d) Motif analysis
    ax = axes[1, 1]
    # Count feedforward loops and feedback loops
    ffl_count = 0
    fbl_count = 0
    for n1 in G.nodes:
        for n2 in G.successors(n1):
            for n3 in G.successors(n2):
                if G.has_edge(n1, n3):
                    ffl_count += 1
            if G.has_edge(n2, n1):
                fbl_count += 1

    motifs = {'Feed-Forward\nLoops': ffl_count, 'Feedback\nLoops': fbl_count,
              'Activating\nEdges': len(act_edges), 'Repressing\nEdges': len(rep_edges)}
    ax.bar(motifs.keys(), motifs.values(), color=['#4CAF50', '#F44336', '#2196F3', '#FF9800'])
    ax.set_title('(d) Network Motif Analysis')
    ax.set_ylabel('Count')

    plt.suptitle('Figure 4: WRKY/TGA Transcription Factor Network', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig4_tf_network.png')
    plt.close()
    print("  Saved fig4_tf_network.png")

    return {
        'nodes': len(G.nodes), 'edges': len(G.edges),
        'ffl_count': ffl_count, 'fbl_count': fbl_count,
        'top_centrality': nodes_sorted[:3],
        'gene_activity': gene_activity,
    }


# ============================================================
# 5. Pathogen-Host Coevolution Game Theory
# ============================================================
def run_game_theory():
    print("=== 5. Game Theory Analysis ===")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # (a) Basic gene-for-gene payoff matrix
    # Host strategies: R (resistant), S (susceptible)
    # Pathogen strategies: Avr (avirulent), vir (virulent)
    cost_R = 0.1    # fitness cost of resistance
    cost_vir = 0.15 # fitness cost of virulence
    benefit = 1.0   # benefit of successful infection

    payoff_host = np.array([
        [1 - cost_R, 1 - cost_R],           # R vs [Avr, vir]
        [1 - benefit, 1]                     # S vs [Avr, vir]
    ])
    # Corrected: R-gene detects Avr → no infection
    payoff_host[0, 0] = 1 - cost_R   # R vs Avr: resistance works
    payoff_host[0, 1] = 1 - cost_R - benefit * 0.3  # R vs vir: partial infection
    payoff_host[1, 0] = 1 - benefit  # S vs Avr: fully infected
    payoff_host[1, 1] = 1 - benefit  # S vs vir: fully infected

    payoff_pathogen = np.array([
        [0, benefit - cost_vir],        # Avr vs [R, S]
        [benefit, benefit]               # vir vs [R, S] -- but Avr detected by R
    ])
    # Corrected
    payoff_pathogen[0, 0] = 0           # Avr vs R: detected, no success
    payoff_pathogen[0, 1] = benefit     # Avr vs S: full success, no cost
    payoff_pathogen[1, 0] = benefit - cost_vir  # vir vs R: evade but costly
    payoff_pathogen[1, 1] = benefit - cost_vir  # vir vs S: success but unnecessary cost

    ax = axes[0, 0]
    im = ax.imshow(payoff_host, cmap='RdYlGn', aspect='auto')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Avr', 'Vir'])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['R-gene', 'Susceptible'])
    ax.set_xlabel('Pathogen Strategy')
    ax.set_ylabel('Host Strategy')
    ax.set_title('(a) Host Payoff Matrix')
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f'{payoff_host[i, j]:.2f}', ha='center', va='center', fontsize=12)
    plt.colorbar(im, ax=ax)

    ax = axes[0, 1]
    im2 = ax.imshow(payoff_pathogen, cmap='RdYlGn', aspect='auto')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['R-gene Host', 'Susceptible Host'])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Avirulent', 'Virulent'])
    ax.set_xlabel('Host Strategy')
    ax.set_ylabel('Pathogen Strategy')
    ax.set_title('(b) Pathogen Payoff Matrix')
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f'{payoff_pathogen[i, j]:.2f}', ha='center', va='center', fontsize=12)
    plt.colorbar(im2, ax=ax)

    # (c) Replicator dynamics
    def replicator(t, y, payoffs_h, payoffs_p):
        p = y[0]  # freq of R in host
        q = y[1]  # freq of Avr in pathogen
        p = np.clip(p, 0.001, 0.999)
        q = np.clip(q, 0.001, 0.999)

        # Host fitness
        f_R = payoffs_h[0, 0] * q + payoffs_h[0, 1] * (1 - q)
        f_S = payoffs_h[1, 0] * q + payoffs_h[1, 1] * (1 - q)
        f_avg_h = p * f_R + (1 - p) * f_S

        # Pathogen fitness
        f_Avr = payoffs_p[0, 0] * p + payoffs_p[0, 1] * (1 - p)
        f_Vir = payoffs_p[1, 0] * p + payoffs_p[1, 1] * (1 - p)
        f_avg_p = q * f_Avr + (1 - q) * f_Vir

        dp = p * (f_R - f_avg_h)
        dq = q * (f_Avr - f_avg_p)
        return [dp, dq]

    t_span = (0, 200)
    t_eval = np.linspace(0, 200, 2000)

    ax = axes[0, 2]
    for p0, q0 in [(0.1, 0.9), (0.5, 0.5), (0.9, 0.1), (0.3, 0.7), (0.7, 0.3)]:
        sol = solve_ivp(replicator, t_span, [p0, q0], t_eval=t_eval,
                        args=(payoff_host, payoff_pathogen), method='RK45')
        ax.plot(sol.t, sol.y[0], '-', label=f'Host R (p₀={p0})')
    ax.set_xlabel('Generations')
    ax.set_ylabel('Frequency of R-gene hosts')
    ax.set_title('(c) Replicator Dynamics: Host')
    ax.legend(fontsize=7)

    # (d) Phase portrait
    ax = axes[1, 0]
    for p0 in np.arange(0.1, 1.0, 0.2):
        for q0 in np.arange(0.1, 1.0, 0.2):
            sol = solve_ivp(replicator, t_span, [p0, q0], t_eval=t_eval,
                            args=(payoff_host, payoff_pathogen), method='RK45')
            ax.plot(sol.y[0], sol.y[1], 'b-', alpha=0.3, linewidth=0.8)
            ax.plot(sol.y[0, 0], sol.y[1, 0], 'go', markersize=4)
            ax.plot(sol.y[0, -1], sol.y[1, -1], 'rs', markersize=4)
    ax.set_xlabel('Freq. R-gene (Host)')
    ax.set_ylabel('Freq. Avirulent (Pathogen)')
    ax.set_title('(d) Coevolutionary Phase Portrait')

    # (e) Arms race escalation model
    ax = axes[1, 1]
    generations = 100
    host_defense = np.zeros(generations)
    path_attack = np.zeros(generations)
    host_defense[0] = 1.0
    path_attack[0] = 0.5
    mu_h, mu_p = 0.1, 0.12  # mutation rates
    c_h, c_p = 0.05, 0.04   # costs

    for g in range(1, generations):
        if path_attack[g-1] > host_defense[g-1]:
            host_defense[g] = host_defense[g-1] + mu_h * (1 + 0.1 * np.random.randn())
        else:
            host_defense[g] = host_defense[g-1] - c_h + 0.02 * np.random.randn()
        if host_defense[g-1] > path_attack[g-1]:
            path_attack[g] = path_attack[g-1] + mu_p * (1 + 0.1 * np.random.randn())
        else:
            path_attack[g] = path_attack[g-1] - c_p + 0.02 * np.random.randn()

    ax.plot(range(generations), host_defense, 'b-', label='Host Defense Level', linewidth=2)
    ax.plot(range(generations), path_attack, 'r-', label='Pathogen Attack Level', linewidth=2)
    ax.fill_between(range(generations), host_defense, path_attack,
                     where=host_defense > path_attack, alpha=0.2, color='blue', label='Host advantage')
    ax.fill_between(range(generations), host_defense, path_attack,
                     where=path_attack > host_defense, alpha=0.2, color='red', label='Pathogen advantage')
    ax.set_xlabel('Generations')
    ax.set_ylabel('Trait Level')
    ax.set_title('(e) Arms Race Escalation')
    ax.legend(fontsize=7)

    # (f) ESS analysis with varying cost
    ax = axes[1, 2]
    costs = np.linspace(0.01, 0.5, 50)
    ess_host = []
    ess_path = []
    for c in costs:
        ph = np.array([
            [1 - c, 1 - c - benefit * 0.3],
            [1 - benefit, 1 - benefit]
        ])
        pp = np.array([
            [0, benefit],
            [benefit - c, benefit - c]
        ])
        # Find interior equilibrium if exists
        # p* where pathogen is indifferent
        denom = (pp[0,0] - pp[1,0]) - (pp[0,1] - pp[1,1])
        if abs(denom) > 1e-10:
            p_star = (pp[0,1] - pp[1,1]) / denom
        else:
            p_star = 0.5
        denom2 = (ph[0,0] - ph[0,1]) - (ph[1,0] - ph[1,1])
        if abs(denom2) > 1e-10:
            q_star = (ph[1,1] - ph[0,1]) / denom2
        else:
            q_star = 0.5
        ess_host.append(np.clip(p_star, 0, 1))
        ess_path.append(np.clip(q_star, 0, 1))

    ax.plot(costs, ess_host, 'b-', label='ESS: Host R freq.', linewidth=2)
    ax.plot(costs, ess_path, 'r-', label='ESS: Pathogen Avr freq.', linewidth=2)
    ax.set_xlabel('Cost of Resistance/Virulence')
    ax.set_ylabel('ESS Frequency')
    ax.set_title('(f) ESS vs Cost Parameter')
    ax.legend()

    plt.suptitle('Figure 5: Pathogen-Host Coevolution Game Theory', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig5_game_theory.png')
    plt.close()
    print("  Saved fig5_game_theory.png")

    return {'payoff_host': payoff_host.tolist(), 'payoff_pathogen': payoff_pathogen.tolist()}


# ============================================================
# 6. Rice Blast Resistance Case Study
# ============================================================
def run_rice_blast():
    print("=== 6. Rice Blast Case Study ===")

    # Integrated model for rice (Oryza sativa) - M. oryzae interaction
    # Key R genes: Pi-ta, Pi-b, Pi-k, Piz-t, Pi9
    # Key pathways: OsCEBiP/OsCERK1 (PTI), NLR (ETI)

    def rice_blast_model(t, y, params):
        # State variables
        chitin = y[0]      # fungal chitin (PAMP signal)
        CEBiP = y[1]       # OsCEBiP receptor
        CERK1 = y[2]       # OsCERK1 active complex
        MAPK_rice = y[3]   # OsMPK3/6
        SA_rice = y[4]     # SA in rice
        JA_rice = y[5]     # JA in rice
        WRKY45 = y[6]      # OsWRKY45 (SA pathway TF)
        WRKY13 = y[7]      # OsWRKY13 (regulates SA/JA balance)
        PR1a = y[8]        # Defense PR gene
        Pita = y[9]        # Pi-ta NLR protein level
        ROS = y[10]        # Reactive oxygen species
        HR = y[11]         # Hypersensitive response

        effector = params.get('effector', 0)
        avr_pita = params.get('avr_pita', 0)  # AvrPi-ta effector
        has_pita = params.get('has_pita', 1)   # Whether rice has Pi-ta

        # PTI signaling
        chitin_prod = params['chitin_rate'] * np.exp(-0.02 * t)
        cebip_bind = 0.1 * chitin * CEBiP
        cerk1_act = 0.08 * cebip_bind - 0.02 * CERK1
        mapk_act = 0.15 * CERK1 - 0.05 * MAPK_rice

        # ETI signaling (Pi-ta recognizes AvrPi-ta)
        eti_signal = has_pita * 0.5 * Pita * avr_pita / (1.0 + avr_pita)

        # Hormone synthesis
        sa_prod = 0.3 * MAPK_rice + 0.5 * eti_signal - 0.05 * SA_rice
        ja_prod = 0.1 * MAPK_rice - 0.03 * JA_rice
        # WRKY13 promotes SA, suppresses JA
        sa_prod += 0.2 * WRKY13
        ja_prod -= 0.15 * WRKY13

        # TF activation
        wrky45_act = 0.2 * SA_rice / (1.0 + SA_rice) - 0.05 * WRKY45
        wrky13_act = 0.15 * MAPK_rice / (0.5 + MAPK_rice) + 0.1 * eti_signal - 0.04 * WRKY13

        # Defense gene expression
        pr1a_expr = 0.3 * WRKY45 + 0.2 * SA_rice / (1.0 + SA_rice) - 0.03 * PR1a

        # ROS burst
        ros_prod = 0.2 * CERK1 + 0.4 * eti_signal - 0.1 * ROS

        # HR (only with strong ETI)
        hr_prog = 0.1 * eti_signal * ROS / (1.0 + ROS) - 0.01 * HR

        dchitin = chitin_prod - 0.05 * chitin
        dCEBiP = -cebip_bind + 0.05  # recycling
        dCERK1 = cerk1_act
        dMAPK = mapk_act + 0.3 * eti_signal
        dSA = sa_prod
        dJA = ja_prod
        dWRKY45 = wrky45_act
        dWRKY13 = wrky13_act
        dPR1a = pr1a_expr
        dPita = -0.001 * Pita  # slow degradation
        dROS = ros_prod
        dHR = hr_prog

        return [dchitin, dCEBiP, dCERK1, dMAPK, dSA, dJA,
                dWRKY45, dWRKY13, dPR1a, dPita, dROS, dHR]

    t_span = (0, 120)
    t_eval = np.linspace(0, 120, 600)
    y0 = [0, 1.0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # (a) Compatible interaction (no Pi-ta, virulent strain)
    params_compat = {'chitin_rate': 2.0, 'effector': 1, 'avr_pita': 1.0, 'has_pita': 0}
    sol_c = solve_ivp(rice_blast_model, t_span, y0, t_eval=t_eval,
                      args=(params_compat,), method='RK45')
    axes[0, 0].plot(sol_c.t, sol_c.y[3], label='MAPK', linewidth=2)
    axes[0, 0].plot(sol_c.t, sol_c.y[4], label='SA', linewidth=2)
    axes[0, 0].plot(sol_c.t, sol_c.y[10], label='ROS', linewidth=2)
    axes[0, 0].plot(sol_c.t, sol_c.y[11], label='HR', linewidth=2)
    axes[0, 0].set_title('(a) Compatible (No Pi-ta)')
    axes[0, 0].set_xlabel('Time (hpi)')
    axes[0, 0].set_ylabel('Level (a.u.)')
    axes[0, 0].legend(fontsize=7)

    # (b) Incompatible interaction (Pi-ta present, AvrPi-ta)
    params_incompat = {'chitin_rate': 2.0, 'effector': 1, 'avr_pita': 1.0, 'has_pita': 1}
    sol_i = solve_ivp(rice_blast_model, t_span, y0, t_eval=t_eval,
                      args=(params_incompat,), method='RK45')
    axes[0, 1].plot(sol_i.t, sol_i.y[3], label='MAPK', linewidth=2)
    axes[0, 1].plot(sol_i.t, sol_i.y[4], label='SA', linewidth=2)
    axes[0, 1].plot(sol_i.t, sol_i.y[10], label='ROS', linewidth=2)
    axes[0, 1].plot(sol_i.t, sol_i.y[11], label='HR', linewidth=2)
    axes[0, 1].set_title('(b) Incompatible (Pi-ta + AvrPi-ta)')
    axes[0, 1].set_xlabel('Time (hpi)')
    axes[0, 1].legend(fontsize=7)

    # (c) Defense gene comparison
    axes[0, 2].plot(sol_c.t, sol_c.y[8], 'b--', label='PR1a (Compatible)', linewidth=2)
    axes[0, 2].plot(sol_i.t, sol_i.y[8], 'r-', label='PR1a (Incompatible)', linewidth=2)
    axes[0, 2].plot(sol_c.t, sol_c.y[6], 'g--', label='WRKY45 (Compatible)', linewidth=1.5)
    axes[0, 2].plot(sol_i.t, sol_i.y[6], 'g-', label='WRKY45 (Incompatible)', linewidth=1.5)
    axes[0, 2].set_title('(c) Defense Gene Expression')
    axes[0, 2].set_xlabel('Time (hpi)')
    axes[0, 2].set_ylabel('Expression Level')
    axes[0, 2].legend(fontsize=7)

    # (d) R-gene dose effect
    doses = [0, 0.25, 0.5, 1.0, 2.0]
    max_hr = []
    max_ros = []
    max_pr = []
    for dose in doses:
        y0_d = y0.copy()
        y0_d[9] = dose
        sol_d = solve_ivp(rice_blast_model, t_span, y0_d, t_eval=t_eval,
                          args=(params_incompat,), method='RK45')
        max_hr.append(np.max(sol_d.y[11]))
        max_ros.append(np.max(sol_d.y[10]))
        max_pr.append(np.max(sol_d.y[8]))

    ax = axes[1, 0]
    ax.plot(doses, max_hr, 'rs-', label='Max HR', linewidth=2)
    ax.plot(doses, max_ros, 'bo-', label='Max ROS', linewidth=2)
    ax.plot(doses, max_pr, 'g^-', label='Max PR1a', linewidth=2)
    ax.set_xlabel('Pi-ta Expression Level')
    ax.set_ylabel('Max Response')
    ax.set_title('(d) R-gene Dosage Effect')
    ax.legend()

    # (e) Multi-R-gene stacking simulation
    ax = axes[1, 1]
    # Simulate adding multiple R genes
    r_gene_counts = [0, 1, 2, 3, 4, 5]
    resistance_scores = []
    for n_genes in r_gene_counts:
        params_stack = params_incompat.copy()
        # Each R gene contributes additively to ETI signal
        y0_stack = y0.copy()
        y0_stack[9] = 0.5 * n_genes  # cumulative R gene effect
        sol_s = solve_ivp(rice_blast_model, t_span, y0_stack, t_eval=t_eval,
                          args=(params_stack,), method='RK45')
        # Resistance score: HR + PR expression
        score = np.max(sol_s.y[11]) + np.max(sol_s.y[8])
        resistance_scores.append(score)

    ax.bar(r_gene_counts, resistance_scores, color='forestgreen', alpha=0.7)
    ax.set_xlabel('Number of Stacked R Genes')
    ax.set_ylabel('Resistance Score')
    ax.set_title('(e) R-gene Stacking Effect')

    # (f) Temporal dynamics of PTI+ETI integration
    ax = axes[1, 2]
    # Show how PTI primes ETI
    ax.fill_between(sol_i.t, 0, sol_i.y[2], alpha=0.3, color='blue', label='PTI (CERK1)')
    ax.fill_between(sol_i.t, 0, sol_i.y[11], alpha=0.3, color='red', label='ETI (HR)')
    ax.plot(sol_i.t, sol_i.y[3], 'k-', label='MAPK (shared)', linewidth=2)
    ax.set_xlabel('Time (hpi)')
    ax.set_ylabel('Signal Intensity')
    ax.set_title('(f) PTI-ETI Integration')
    ax.legend(fontsize=7)

    plt.suptitle('Figure 6: Rice Blast Resistance Case Study (Oryza sativa - M. oryzae)',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig6_rice_blast.png')
    plt.close()
    print("  Saved fig6_rice_blast.png")

    return {
        'compatible_max_HR': np.max(sol_c.y[11]),
        'incompatible_max_HR': np.max(sol_i.y[11]),
        'compatible_max_PR1a': np.max(sol_c.y[8]),
        'incompatible_max_PR1a': np.max(sol_i.y[8]),
        'resistance_scores': dict(zip(r_gene_counts, resistance_scores)),
    }


# ============================================================
# 7. SBML/COPASI-compatible pathway specification
# ============================================================
def generate_copasi_spec():
    print("=== 7. COPASI/CellDesigner Pathway Specification ===")

    spec = """<?xml version="1.0" encoding="UTF-8"?>
<!-- SBML Model: Plant PTI-ETI Signaling Network -->
<!-- Compatible with COPASI and CellDesigner -->
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="PlantImmunityModel" name="PTI-ETI Signaling Network">
    <listOfCompartments>
      <compartment id="apoplast" name="Apoplast" size="1" constant="true"/>
      <compartment id="cytoplasm" name="Cytoplasm" size="1" constant="true"/>
      <compartment id="nucleus" name="Nucleus" size="1" constant="true"/>
    </listOfCompartments>

    <listOfSpecies>
      <!-- Receptors -->
      <species id="FLS2" name="FLS2" compartment="cytoplasm" initialConcentration="10"/>
      <species id="flg22" name="flg22" compartment="apoplast" initialConcentration="5"/>
      <species id="FLS2_flg22" name="FLS2:flg22 Complex" compartment="cytoplasm" initialConcentration="0"/>
      <species id="BAK1" name="BAK1" compartment="cytoplasm" initialConcentration="5"/>
      <species id="FLS2_BAK1" name="FLS2:BAK1 Active Complex" compartment="cytoplasm" initialConcentration="0"/>

      <!-- NLR (ETI) -->
      <species id="NLR" name="NLR Receptor" compartment="cytoplasm" initialConcentration="5"/>
      <species id="Effector" name="Pathogen Effector" compartment="cytoplasm" initialConcentration="0"/>
      <species id="NLR_Eff" name="NLR:Effector Complex" compartment="cytoplasm" initialConcentration="0"/>

      <!-- MAPK Cascade -->
      <species id="MAPKKK" name="MAPKKK (MEKK1)" compartment="cytoplasm" initialConcentration="1"/>
      <species id="MAPKKK_P" name="MAPKKK-P" compartment="cytoplasm" initialConcentration="0"/>
      <species id="MAPKK" name="MAPKK (MKK4/5)" compartment="cytoplasm" initialConcentration="1"/>
      <species id="MAPKK_P" name="MAPKK-P" compartment="cytoplasm" initialConcentration="0"/>
      <species id="MAPK" name="MAPK (MPK3/6)" compartment="cytoplasm" initialConcentration="1"/>
      <species id="MAPK_P" name="MAPK-P" compartment="cytoplasm" initialConcentration="0"/>

      <!-- Hormones -->
      <species id="SA" name="Salicylic Acid" compartment="cytoplasm" initialConcentration="0"/>
      <species id="JA" name="Jasmonic Acid" compartment="cytoplasm" initialConcentration="0"/>
      <species id="NPR1_inactive" name="NPR1 (inactive)" compartment="cytoplasm" initialConcentration="5"/>
      <species id="NPR1_active" name="NPR1 (active)" compartment="nucleus" initialConcentration="0"/>
      <species id="JAZ" name="JAZ Repressor" compartment="nucleus" initialConcentration="5"/>

      <!-- Transcription Factors -->
      <species id="WRKY33" name="WRKY33" compartment="nucleus" initialConcentration="0"/>
      <species id="WRKY70" name="WRKY70" compartment="nucleus" initialConcentration="0"/>
      <species id="TGA2" name="TGA2" compartment="nucleus" initialConcentration="1"/>
      <species id="MYC2" name="MYC2" compartment="nucleus" initialConcentration="0"/>

      <!-- Defense Outputs -->
      <species id="PR1" name="PR1" compartment="cytoplasm" initialConcentration="0"/>
      <species id="PDF12" name="PDF1.2" compartment="cytoplasm" initialConcentration="0"/>
      <species id="ROS" name="Reactive Oxygen Species" compartment="apoplast" initialConcentration="0"/>
      <species id="Callose" name="Callose Deposit" compartment="cytoplasm" initialConcentration="0"/>
    </listOfSpecies>

    <listOfReactions>
      <!-- R1: FLS2 + flg22 -> FLS2:flg22 -->
      <reaction id="R1" name="PAMP binding" reversible="false">
        <listOfReactants>
          <speciesReference species="FLS2" stoichiometry="1"/>
          <speciesReference species="flg22" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="FLS2_flg22" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><times/><cn>0.1</cn><ci>FLS2</ci><ci>flg22</ci></apply>
        </math></kineticLaw>
      </reaction>

      <!-- R2: FLS2:flg22 + BAK1 -> FLS2:BAK1 -->
      <reaction id="R2" name="Co-receptor recruitment" reversible="false">
        <listOfReactants>
          <speciesReference species="FLS2_flg22" stoichiometry="1"/>
          <speciesReference species="BAK1" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="FLS2_BAK1" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><times/><cn>0.05</cn><ci>FLS2_flg22</ci><ci>BAK1</ci></apply>
        </math></kineticLaw>
      </reaction>

      <!-- R3-R8: MAPK Cascade (phosphorylation/dephosphorylation) -->
      <reaction id="R3" name="MAPKKK activation" reversible="false">
        <listOfReactants><speciesReference species="MAPKKK" stoichiometry="1"/></listOfReactants>
        <listOfProducts><speciesReference species="MAPKKK_P" stoichiometry="1"/></listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><divide/>
            <apply><times/><cn>0.5</cn><ci>FLS2_BAK1</ci><ci>MAPKKK</ci></apply>
            <apply><plus/><cn>0.5</cn><ci>MAPKKK</ci></apply>
          </apply>
        </math></kineticLaw>
      </reaction>

      <!-- R9: SA synthesis -->
      <reaction id="R9" name="SA biosynthesis" reversible="false">
        <listOfProducts><speciesReference species="SA" stoichiometry="1"/></listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><times/><cn>0.3</cn><ci>MAPK_P</ci></apply>
        </math></kineticLaw>
      </reaction>

      <!-- R10: NPR1 activation by SA -->
      <reaction id="R10" name="NPR1 activation" reversible="false">
        <listOfReactants><speciesReference species="NPR1_inactive" stoichiometry="1"/></listOfReactants>
        <listOfProducts><speciesReference species="NPR1_active" stoichiometry="1"/></listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><divide/>
            <apply><times/><cn>0.2</cn><ci>SA</ci><ci>NPR1_inactive</ci></apply>
            <apply><plus/><cn>1.0</cn><ci>SA</ci></apply>
          </apply>
        </math></kineticLaw>
      </reaction>

      <!-- R11: PR1 expression -->
      <reaction id="R11" name="PR1 transcription" reversible="false">
        <listOfProducts><speciesReference species="PR1" stoichiometry="1"/></listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><divide/>
            <apply><times/><cn>1.0</cn><ci>NPR1_active</ci><ci>TGA2</ci></apply>
            <apply><plus/><cn>0.5</cn><ci>NPR1_active</ci></apply>
          </apply>
        </math></kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
"""
    with open('pti_eti_model.sbml', 'w') as f:
        f.write(spec)
    print("  Saved pti_eti_model.sbml")
    return True


# ============================================================
# Main Execution
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Plant PTI-ETI Signaling Model - Full Simulation Suite")
    print("=" * 60)

    results = {}
    results['receptor'] = run_receptor_model()
    results['mapk'] = run_mapk_cascade()
    results['sa_ja'] = run_sa_ja_crosstalk()
    results['tf_network'] = run_tf_network()
    results['game_theory'] = run_game_theory()
    results['rice_blast'] = run_rice_blast()
    generate_copasi_spec()

    print("\n" + "=" * 60)
    print("Summary of Key Results:")
    print("=" * 60)
    print(f"Receptor Model:")
    print(f"  EC50 ≈ {results['receptor']['ec50_approx']:.2f} nM")
    print(f"  PTI max signal: {results['receptor']['pti_max_signal']:.2f}")
    print(f"  ETI max signal: {results['receptor']['eti_max_signal']:.2f}")
    print(f"  ETI/PTI ratio: {results['receptor']['eti_max_signal']/results['receptor']['pti_max_signal']:.2f}")
    print(f"\nMAPK Cascade:")
    print(f"  Hill coefficient: {results['mapk']['hill_coefficient']:.2f}")
    print(f"\nSA/JA Crosstalk:")
    print(f"  Biotroph SS: SA={results['sa_ja']['biotroph_SA_ss']:.2f}, JA={results['sa_ja']['biotroph_JA_ss']:.2f}")
    print(f"  Necrotroph SS: SA={results['sa_ja']['necrotroph_SA_ss']:.2f}, JA={results['sa_ja']['necrotroph_JA_ss']:.2f}")
    print(f"\nTF Network:")
    print(f"  Nodes: {results['tf_network']['nodes']}, Edges: {results['tf_network']['edges']}")
    print(f"  FFLs: {results['tf_network']['ffl_count']}, FBLs: {results['tf_network']['fbl_count']}")
    print(f"  Top centrality: {results['tf_network']['top_centrality']}")
    print(f"\nRice Blast:")
    print(f"  Compatible HR: {results['rice_blast']['compatible_max_HR']:.4f}")
    print(f"  Incompatible HR: {results['rice_blast']['incompatible_max_HR']:.4f}")
    print(f"  HR ratio (incompat/compat): {results['rice_blast']['incompatible_max_HR']/(results['rice_blast']['compatible_max_HR']+1e-10):.2f}")
    print(f"  R-gene stacking scores: {results['rice_blast']['resistance_scores']}")
    print("\nAll simulations completed successfully!")
    print(f"Figures saved to: {FIGDIR}/")
    print(f"SBML model saved to: pti_eti_model.sbml")
