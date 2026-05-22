"""
Quantum Network Routing: Fidelity-Aware Path Selection
=======================================================
Implements:
  - Fidelity-weighted Dijkstra (maximize end-to-end fidelity)
  - Bandwidth-fidelity Pareto routing
  - Tokyo QKD network topology
Reference: Van Meter et al. ACM SIGCOMM (2013), Caleffi (2017)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import json
from itertools import islice

np.random.seed(42)


# ─── Tokyo QKD network topology ───────────────────────────────────────────────
# Based on the Tokyo QKD Network (2010-2015) with 6 nodes
# Refs: Sasaki et al. Optics Express 19, 10387 (2011)
TOKYO_NODES = {
    'NEC':         (35.720, 139.742),   # NEC Tamagawa
    'Koganei':     (35.700, 139.504),   # NICT Koganei
    'Hakusan':     (35.728, 139.744),   # JGN Hakusan
    'Tokyo_Univ':  (35.715, 139.762),   # Tokyo University
    'NTT':         (35.677, 139.753),   # NTT Musashino
    'Otemachi':    (35.690, 139.764),   # NTT Otemachi
}

# Distances in km (approximate fiber lengths)
TOKYO_EDGES = [
    ('NEC', 'Koganei',    45, 'DPS-QKD'),
    ('NEC', 'Hakusan',     7, 'BB84'),
    ('Koganei', 'Hakusan', 40, 'BB84'),
    ('Hakusan', 'Tokyo_Univ', 3, 'BB84'),
    ('Hakusan', 'NTT',    15, 'COW'),
    ('NTT', 'Otemachi',   14, 'DPS-QKD'),
    ('Tokyo_Univ', 'Otemachi', 5, 'BB84'),
    ('NEC', 'NTT',        25, 'BB84'),
    ('Koganei', 'NTT',    30, 'E91'),
]


def build_tokyo_network():
    """Build Tokyo QKD network graph with fidelity and capacity attributes."""
    G = nx.Graph()

    for node, pos in TOKYO_NODES.items():
        G.add_node(node, pos=pos)

    for src, dst, dist_km, proto in TOKYO_EDGES:
        # Fiber transmittance
        eta = 10 ** (-0.2 * dist_km / 10)
        # Fidelity model: F = eta * 0.98 + (1-eta)*0.25 (depolarizing noise)
        F = eta * 0.97 + (1 - eta) * 0.25
        # Key rate (simplified, proportional to transmittance)
        key_rate_kbps = max(eta * 100 - 0.1, 0)  # kbps
        # Weight for fidelity routing: -log(F) so path maximizes F
        fid_weight = -np.log(max(F, 1e-10))
        G.add_edge(src, dst,
                   distance_km=dist_km,
                   protocol=proto,
                   fidelity=F,
                   key_rate_kbps=key_rate_kbps,
                   fid_weight=fid_weight,
                   latency_ms=dist_km / 200)  # 200 km/ms fiber speed
    return G


# ─── routing algorithms ───────────────────────────────────────────────────────
def fidelity_dijkstra(G, source, target):
    """
    Find path maximizing end-to-end fidelity (= minimizing sum of -log(F_i)).
    """
    try:
        path = nx.shortest_path(G, source=source, target=target, weight='fid_weight')
        # Compute end-to-end fidelity (product of edge fidelities)
        edges = list(zip(path[:-1], path[1:]))
        fidelity = 1.0
        total_dist = 0
        for u, v in edges:
            fidelity *= G[u][v]['fidelity']
            total_dist += G[u][v]['distance_km']
        return path, fidelity, total_dist
    except nx.NetworkXNoPath:
        return None, 0.0, 0.0


def k_shortest_paths(G, source, target, k=5, weight='fid_weight'):
    """K-shortest paths using Yen's algorithm."""
    try:
        paths = list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))
        results = []
        for path in paths:
            edges = list(zip(path[:-1], path[1:]))
            F = 1.0
            dist = 0
            rate = float('inf')
            latency = 0
            for u, v in edges:
                F *= G[u][v]['fidelity']
                dist += G[u][v]['distance_km']
                rate = min(rate, G[u][v]['key_rate_kbps'])
                latency += G[u][v]['latency_ms']
            results.append({
                'path': path,
                'fidelity': F,
                'distance_km': dist,
                'bottleneck_rate_kbps': rate,
                'latency_ms': latency,
                'hops': len(path) - 1
            })
        return results
    except Exception:
        return []


def bandwidth_aware_routing(G, source, target, fidelity_threshold=0.5):
    """
    Route maximizing key rate subject to fidelity threshold.
    """
    G_filtered = nx.Graph(
        [(u, v, d) for u, v, d in G.edges(data=True) if d['fidelity'] >= fidelity_threshold]
    )
    G_filtered.add_nodes_from(G.nodes(data=True))

    if not nx.has_path(G_filtered, source, target):
        return None, 0, 0, 0

    # Use inverse key_rate as weight to maximize rate
    for u, v, d in G_filtered.edges(data=True):
        G_filtered[u][v]['rate_weight'] = 1.0 / max(d['key_rate_kbps'], 1e-6)

    try:
        path = nx.shortest_path(G_filtered, source, target, weight='rate_weight')
        edges = list(zip(path[:-1], path[1:]))
        F = 1.0
        bottleneck_rate = float('inf')
        for u, v in edges:
            F *= G_filtered[u][v]['fidelity']
            bottleneck_rate = min(bottleneck_rate, G_filtered[u][v]['key_rate_kbps'])
        return path, F, bottleneck_rate, len(path) - 1
    except Exception:
        return None, 0, 0, 0


# ─── build network and run routing ───────────────────────────────────────────
G = build_tokyo_network()

print("Tokyo QKD Network Routing Analysis")
print("=" * 50)

routing_results = []
node_pairs = [
    ('Koganei', 'Otemachi'),
    ('NEC', 'Otemachi'),
    ('NEC', 'Tokyo_Univ'),
    ('Koganei', 'Tokyo_Univ'),
]

for src, dst in node_pairs:
    path, F, dist = fidelity_dijkstra(G, src, dst)
    k_paths = k_shortest_paths(G, src, dst, k=3)
    bw_path, bw_F, bw_rate, bw_hops = bandwidth_aware_routing(G, src, dst, fidelity_threshold=0.3)

    result = {
        'source': src,
        'destination': dst,
        'fidelity_optimal': {
            'path': path,
            'fidelity': F,
            'distance_km': dist,
        },
        'k_shortest_paths': k_paths[:3],
        'bandwidth_optimal': {
            'path': bw_path,
            'fidelity': bw_F,
            'bottleneck_rate_kbps': bw_rate,
            'hops': bw_hops
        }
    }
    routing_results.append(result)

    print(f"\n{src} → {dst}")
    print(f"  Fidelity-optimal: {' → '.join(path)} | F={F:.4f}, D={dist} km")
    if k_paths:
        print(f"  Top-3 alternative paths:")
        for i, p in enumerate(k_paths[:3]):
            print(f"    [{i+1}] {' → '.join(p['path'])} | F={p['fidelity']:.4f}, D={p['distance_km']} km, "
                  f"rate={p['bottleneck_rate_kbps']:.2f} kbps")

# ─── visualize Tokyo network ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

ax = axes[0]
pos = nx.get_node_attributes(G, 'pos')
# Swap lat/lon for plotting
pos_xy = {node: (lon, lat) for node, (lat, lon) in pos.items()}

edge_colors = [G[u][v]['fidelity'] for u, v in G.edges()]
edge_widths = [G[u][v]['key_rate_kbps'] / 5 + 1 for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos_xy, ax=ax, node_size=700, node_color='lightblue',
                       edgecolors='darkblue', linewidths=2)
nx.draw_networkx_labels(G, pos_xy, ax=ax, font_size=9, font_weight='bold')
edges_drawn = nx.draw_networkx_edges(G, pos_xy, ax=ax, width=edge_widths,
                                      edge_color=edge_colors, edge_cmap=plt.cm.RdYlGn,
                                      edge_vmin=0.3, edge_vmax=1.0, arrows=False)

# Fidelity labels on edges
edge_labels = {(u, v): f"F={G[u][v]['fidelity']:.2f}" for u, v in G.edges()}
nx.draw_networkx_edge_labels(G, pos_xy, edge_labels, ax=ax, font_size=7, alpha=0.8)

sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn, norm=plt.Normalize(0.3, 1.0))
sm.set_array([])
plt.colorbar(sm, ax=ax, label='Link fidelity', shrink=0.7)

ax.set_title('Tokyo QKD Network Topology\n(Edge color = fidelity, width = key rate)', fontsize=12)
ax.set_xlabel('Longitude', fontsize=10)
ax.set_ylabel('Latitude', fontsize=10)
ax.grid(True, alpha=0.2)

# Highlight optimal path: Koganei → Otemachi
opt_path, opt_F, _ = fidelity_dijkstra(G, 'Koganei', 'Otemachi')
if opt_path:
    path_edges = list(zip(opt_path[:-1], opt_path[1:]))
    nx.draw_networkx_edges(G, pos_xy, edgelist=path_edges, ax=ax, width=5,
                           edge_color='red', alpha=0.7, style='dashed', arrows=False)
ax.annotate('━ Optimal path: Koganei→Otemachi', xy=(0.02, 0.02),
            xycoords='axes fraction', color='red', fontsize=9)

# 2. Routing comparison bar chart
ax = axes[1]
pairs_labels = [f"{r['source'][:3]}→{r['destination'][:3]}" for r in routing_results]
fid_vals = [r['fidelity_optimal']['fidelity'] for r in routing_results]
bw_vals = [r['bandwidth_optimal']['fidelity'] for r in routing_results]
x = np.arange(len(pairs_labels))
w = 0.35

bars1 = ax.bar(x - w/2, fid_vals, w, label='Fidelity-optimal routing', color='steelblue', alpha=0.8)
bars2 = ax.bar(x + w/2, bw_vals, w, label='Bandwidth-optimal routing', color='tomato', alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(pairs_labels, fontsize=10)
ax.set_ylabel('End-to-end fidelity', fontsize=11)
ax.set_title('Routing Strategy Comparison\n(End-to-end fidelity per node pair)', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim(0, 1)

for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=8)
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig('figures/quantum_routing.png', dpi=200, bbox_inches='tight')
plt.savefig('figures/quantum_routing.svg', bbox_inches='tight')
plt.close()

# ─── Pareto frontier (fidelity vs rate) ──────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
for r in routing_results:
    paths = r['k_shortest_paths']
    if paths:
        fids_p = [p['fidelity'] for p in paths]
        rates_p = [p['bottleneck_rate_kbps'] for p in paths]
        label = f"{r['source'][:4]}→{r['destination'][:4]}"
        ax.scatter(rates_p, fids_p, s=100, label=label, zorder=5)
        for i, (rr, ff) in enumerate(zip(rates_p, fids_p)):
            ax.annotate(f"path{i+1}", (rr, ff), textcoords='offset points',
                        xytext=(5, 5), fontsize=7)

ax.set_xlabel('Bottleneck key rate (kbps)', fontsize=12)
ax.set_ylabel('End-to-end fidelity', fontsize=12)
ax.set_title('Pareto Frontier: Fidelity vs Key Rate\n(Tokyo QKD Network)', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/routing_pareto.png', dpi=200, bbox_inches='tight')
plt.close()

# ─── save results ────────────────────────────────────────────────────────────
with open('results/routing_results.json', 'w') as f:
    json.dump(routing_results, f, indent=2)

print("\nRouting analysis complete. Results saved.")
