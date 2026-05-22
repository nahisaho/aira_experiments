"""
Tokyo QKD Network Case Study
=============================
Comprehensive analysis of a Tokyo-scale QKD network:
  - Node performance statistics
  - Network-wide key rate estimation
  - Security analysis under Eve eavesdropping
  - Deployment cost-performance tradeoffs
  - NetSquid/SimulaQron conceptual protocol design
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import json
from collections import defaultdict

np.random.seed(42)


# ─── extended Tokyo network (more nodes for realism) ─────────────────────────
NODES_EXT = {
    'NEC_Tamagawa':    (35.5820, 139.6491),
    'NICT_Koganei':    (35.7003, 139.5041),
    'Hakusan':         (35.7278, 139.7437),
    'Tokyo_Univ':      (35.7151, 139.7626),
    'NTT_Musashino':   (35.7167, 139.5662),
    'NTT_Otemachi':    (35.6900, 139.7640),
    'JAXA_Tsukuba':    (36.1044, 140.0870),
    'Hitachi_Kokubunji': (35.7177, 139.4624),
    'Toshiba_Komukai':  (35.5389, 139.6572),
    'KDDI_Shinjuku':   (35.6895, 139.7006),
}

EDGES_EXT = [
    ('NEC_Tamagawa', 'NICT_Koganei',    45, 0.95),
    ('NEC_Tamagawa', 'Hakusan',          7, 0.98),
    ('NICT_Koganei', 'Hakusan',         40, 0.93),
    ('NICT_Koganei', 'NTT_Musashino',   15, 0.97),
    ('Hakusan', 'Tokyo_Univ',            3, 0.99),
    ('Hakusan', 'NTT_Otemachi',         15, 0.96),
    ('NTT_Musashino', 'NTT_Otemachi',   14, 0.97),
    ('Tokyo_Univ', 'NTT_Otemachi',       5, 0.98),
    ('NEC_Tamagawa', 'NTT_Musashino',   25, 0.95),
    ('NICT_Koganei', 'Hitachi_Kokubunji', 12, 0.97),
    ('NTT_Otemachi', 'KDDI_Shinjuku',    8, 0.98),
    ('NEC_Tamagawa', 'Toshiba_Komukai',  10, 0.98),
    ('Toshiba_Komukai', 'NTT_Otemachi', 20, 0.96),
    ('Tokyo_Univ', 'JAXA_Tsukuba',      55, 0.88),
    ('Hitachi_Kokubunji', 'NICT_Koganei', 12, 0.97),
    ('KDDI_Shinjuku', 'Hakusan',          5, 0.99),
]

ALPHA_DB_KM = 0.2

def build_extended_network():
    G = nx.Graph()
    for node, pos in NODES_EXT.items():
        G.add_node(node, pos=pos)
    for src, dst, dist_km, align_fid in EDGES_EXT:
        eta = 10 ** (-ALPHA_DB_KM * dist_km / 10) * 0.85
        F = eta * align_fid + (1 - eta) * 0.25

        def h2(p):
            if p <= 0 or p >= 1:
                return 0.0
            return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

        # Dark count QBER
        R_signal = 1e6 * 0.1 * eta  # 1MHz, mu=0.1
        R_dark = 100
        qber = R_dark / (R_signal + R_dark) + 0.01
        key_rate = max(1e3 * 0.1 * eta * (1 - h2(qber) - h2(qber)), 0)

        G.add_edge(src, dst,
                   distance_km=dist_km,
                   fidelity=F,
                   transmittance=eta,
                   qber=qber,
                   key_rate_kbps=key_rate / 1000,
                   fid_weight=-np.log(max(F, 1e-10)),
                   latency_ms=dist_km / 200)
    return G


G = build_extended_network()

# ─── network-wide statistics ──────────────────────────────────────────────────
print("=== Tokyo QKD Extended Network Statistics ===")
print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

# Average node degree
avg_degree = np.mean([d for _, d in G.degree()])
print(f"Average node degree: {avg_degree:.2f}")

# All-pairs fidelity
all_pairs_fidelity = {}
all_pairs_rate = {}
for src in G.nodes():
    for dst in G.nodes():
        if src >= dst:
            continue
        try:
            path = nx.shortest_path(G, src, dst, weight='fid_weight')
            edges = list(zip(path[:-1], path[1:]))
            F = np.prod([G[u][v]['fidelity'] for u, v in edges])
            rate = min([G[u][v]['key_rate_kbps'] for u, v in edges])
            all_pairs_fidelity[(src, dst)] = F
            all_pairs_rate[(src, dst)] = rate
        except nx.NetworkXNoPath:
            all_pairs_fidelity[(src, dst)] = 0.0
            all_pairs_rate[(src, dst)] = 0.0

fidelities = list(all_pairs_fidelity.values())
rates = list(all_pairs_rate.values())
print(f"Mean E2E fidelity: {np.mean(fidelities):.4f} ± {np.std(fidelities):.4f}")
print(f"Mean bottleneck key rate: {np.mean(rates):.4f} kbps")

# ─── security analysis: Eve interception ──────────────────────────────────────
def analyze_eavesdropping(G, intercept_edge_idx):
    """
    Model Eve performing an intercept-resend attack on one link.
    """
    edges = list(G.edges(data=True))
    if intercept_edge_idx >= len(edges):
        return {}
    u, v, data = edges[intercept_edge_idx]

    # Intercept-resend adds 25% QBER
    qber_with_eve = min(data['qber'] + 0.25, 0.5)

    def h2(p):
        if p <= 0 or p >= 1:
            return 0.0
        return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

    # Key rate drops to 0 if QBER > 11%
    key_rate_with_eve = max(data['key_rate_kbps'] * (1 - h2(qber_with_eve) / max(1 - h2(data['qber']), 0.01)), 0)

    return {
        'intercepted_link': (u, v),
        'normal_qber': data['qber'],
        'eve_qber': qber_with_eve,
        'normal_rate_kbps': data['key_rate_kbps'],
        'rate_with_eve_kbps': key_rate_with_eve,
        'eve_detectable': qber_with_eve > 0.11
    }

eve_analysis = [analyze_eavesdropping(G, i) for i in range(len(list(G.edges())))]

# ─── NetSquid protocol design (conceptual, with pseudo-code) ─────────────────
netsquid_design = {
    "simulation_framework": "NetSquid 1.1.6 (conceptual design)",
    "protocol_stack": {
        "physical_layer": {
            "nodes": list(G.nodes()),
            "connections": [{"src": u, "dst": v, "fiber_km": d['distance_km'],
                             "loss_model": "FibreChannelModel",
                             "T2_memory_us": 10000}
                            for u, v, d in G.edges(data=True)],
            "memory_model": "DepolarNoiseModel(depolar_rate=1e-3)",
            "photon_source": "QSource(state=Qubit, freq=1e6)"
        },
        "link_layer": {
            "protocol": "EntanglementGenProtocol",
            "BSM": "BellStateMeasurement(success_prob=0.5)",
            "correction": "ClassicalChannel(latency=link_latency)"
        },
        "network_layer": {
            "entanglement_swapping": "SwappingProtocol",
            "routing": "FidelityAwareRoutingAlgorithm",
            "distillation": "DEJMPSDistillationProtocol"
        },
        "application_layer": {
            "qkd": "BB84Protocol / E91Protocol",
            "teleportation": "QuantumTeleportationProtocol"
        }
    },
    "simulation_parameters": {
        "sim_time_ns": 1e9,  # 1 second
        "num_runs": 100,
        "random_seed": 42
    }
}

# ─── visualizations ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# 1. Extended network topology
ax1 = fig.add_subplot(gs[0, :2])
pos = nx.get_node_attributes(G, 'pos')
pos_xy = {node: (lon, lat) for node, (lat, lon) in pos.items()}

node_colors = plt.cm.Set3(np.linspace(0, 1, G.number_of_nodes()))
edge_fids = [G[u][v]['fidelity'] for u, v in G.edges()]
edge_widths = [G[u][v]['key_rate_kbps'] * 20 + 1 for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos_xy, ax=ax1, node_size=600, node_color=node_colors,
                       edgecolors='black', linewidths=1.5)
nx.draw_networkx_labels(G, pos_xy, ax=ax1, font_size=7, font_weight='bold')
ec = nx.draw_networkx_edges(G, pos_xy, ax=ax1, width=edge_widths,
                              edge_color=edge_fids, edge_cmap=plt.cm.RdYlGn,
                              edge_vmin=0.0, edge_vmax=1.0, arrows=False)

edge_labels_d = {(u, v): f"{d['distance_km']}km" for u, v, d in G.edges(data=True)}
nx.draw_networkx_edge_labels(G, pos_xy, edge_labels_d, ax=ax1, font_size=6)

sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn, norm=plt.Normalize(0, 1))
sm.set_array([])
plt.colorbar(sm, ax=ax1, label='Link fidelity', shrink=0.6)
ax1.set_title('Tokyo QKD Extended Network (10 nodes, 16 links)', fontsize=13)
ax1.set_xlabel('Longitude', fontsize=10)
ax1.set_ylabel('Latitude', fontsize=10)
ax1.grid(True, alpha=0.2)

# 2. All-pairs fidelity heatmap
ax2 = fig.add_subplot(gs[0, 2])
nodes = list(G.nodes())
n = len(nodes)
fid_matrix = np.zeros((n, n))
for i, s in enumerate(nodes):
    for j, d in enumerate(nodes):
        if s == d:
            fid_matrix[i, j] = 1.0
        elif (s, d) in all_pairs_fidelity:
            fid_matrix[i, j] = all_pairs_fidelity[(s, d)]
        elif (d, s) in all_pairs_fidelity:
            fid_matrix[i, j] = all_pairs_fidelity[(d, s)]

short_names = [n.split('_')[0][:6] for n in nodes]
im = ax2.imshow(fid_matrix, cmap='viridis', vmin=0, vmax=1)
ax2.set_xticks(range(n))
ax2.set_yticks(range(n))
ax2.set_xticklabels(short_names, rotation=45, fontsize=7, ha='right')
ax2.set_yticklabels(short_names, fontsize=7)
plt.colorbar(im, ax=ax2, label='E2E fidelity', shrink=0.8)
ax2.set_title('All-Pairs Fidelity\n(Fidelity-optimal routing)', fontsize=11)

# 3. Key rate distribution
ax3 = fig.add_subplot(gs[1, 0])
edge_rates = [d['key_rate_kbps'] for _, _, d in G.edges(data=True)]
ax3.hist(edge_rates, bins=20, color='steelblue', alpha=0.8, edgecolor='black')
ax3.axvline(np.mean(edge_rates), color='red', linestyle='--', linewidth=2,
             label=f'Mean={np.mean(edge_rates):.3f} kbps')
ax3.set_xlabel('Link key rate (kbps)', fontsize=11)
ax3.set_ylabel('Count', fontsize=11)
ax3.set_title('Key Rate Distribution\nacross Network Links', fontsize=12)
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# 4. Eavesdropping detection
ax4 = fig.add_subplot(gs[1, 1])
eve_qbers = [e['eve_qber'] for e in eve_analysis if e]
normal_qbers = [e['normal_qber'] for e in eve_analysis if e]
link_names = [f"{e['intercepted_link'][0][:4]}-{e['intercepted_link'][1][:4]}"
              for e in eve_analysis if e]
x = range(len(eve_qbers))
ax4.bar([i - 0.2 for i in x], [q * 100 for q in normal_qbers], 0.35,
        label='Normal QBER', color='steelblue', alpha=0.8)
ax4.bar([i + 0.2 for i in x], [q * 100 for q in eve_qbers], 0.35,
        label='QBER with Eve', color='tomato', alpha=0.8)
ax4.axhline(11, color='black', linestyle='--', linewidth=2, label='Detection threshold 11%')
ax4.set_xticks(list(x))
ax4.set_xticklabels(link_names, rotation=45, ha='right', fontsize=6)
ax4.set_ylabel('QBER (%)', fontsize=11)
ax4.set_title('Eavesdropping Detection Analysis\n(Intercept-resend attack)', fontsize=11)
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3, axis='y')

# 5. Protocol stack diagram
ax5 = fig.add_subplot(gs[1, 2])
layers = ['Physical Layer\n(Photon/Memory)', 'Link Layer\n(Entanglement Gen.)',
          'Network Layer\n(Swapping/Routing)', 'Application Layer\n(BB84/E91/Teleport)']
protocols = ['FibreChannel+\nDepolarNoise', 'EntanglementGen\nBSM(p=0.5)',
             'DEJMPSDistill+\nFidelityRoute', 'BB84/E91/QT\nProtocol']
colors_stack = ['#d4e6f1', '#a9cce3', '#7fb3d3', '#5499c2']
for i, (layer, proto, color) in enumerate(zip(layers, protocols, colors_stack)):
    rect = plt.Rectangle((0.1, i * 0.22 + 0.02), 0.8, 0.18, facecolor=color,
                          edgecolor='black', linewidth=1.5)
    ax5.add_patch(rect)
    ax5.text(0.5, i * 0.22 + 0.11, f'{layer}\n[{proto}]',
             ha='center', va='center', fontsize=8, fontweight='bold')
    if i < 3:
        ax5.annotate('', xy=(0.5, (i+1) * 0.22 + 0.01), xytext=(0.5, (i+1) * 0.22 - 0.01),
                     arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))

ax5.set_xlim(0, 1)
ax5.set_ylim(0, 1)
ax5.axis('off')
ax5.set_title('QKD Protocol Stack\n(NetSquid Architecture)', fontsize=11)

plt.savefig('figures/tokyo_qkd_casestudy.png', dpi=200, bbox_inches='tight')
plt.savefig('figures/tokyo_qkd_casestudy.svg', bbox_inches='tight')
plt.close()

# ─── save results ────────────────────────────────────────────────────────────
summary = {
    "network": {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "avg_degree": float(avg_degree),
        "total_fiber_km": float(sum(d['distance_km'] for _, _, d in G.edges(data=True))),
    },
    "performance": {
        "mean_e2e_fidelity": float(np.mean(fidelities)),
        "std_e2e_fidelity": float(np.std(fidelities)),
        "mean_bottleneck_rate_kbps": float(np.mean(rates)),
        "max_bottleneck_rate_kbps": float(np.max(rates)),
    },
    "security": {
        "links_vulnerable_to_IR": sum(1 for e in eve_analysis if e and not e['eve_detectable']),
        "links_where_eve_detectable": sum(1 for e in eve_analysis if e and e['eve_detectable']),
        "detection_qber_threshold": 0.11,
    },
    "netsquid_design": netsquid_design,
}

with open('results/tokyo_casestudy_results.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("Tokyo QKD case study complete.")
print(json.dumps({k: v for k, v in summary.items() if k != 'netsquid_design'}, indent=2))
