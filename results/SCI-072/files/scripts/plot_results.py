#!/usr/bin/env python3
"""Generate visualization figures for MAPF benchmark results."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150, 'savefig.bbox': 'tight'})
FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

def read_csv(path):
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)

# ============================================================
# Figure 1: Scalability - Runtime vs Agents
# ============================================================
def plot_scalability():
    data = read_csv('benchmarks/scalability.csv')
    algos = {'CBS': [], 'EECBS': [], 'PP': []}
    for row in data:
        if row['mapType'] != 'random': continue
        alg = row['algorithm']
        if alg not in algos: continue
        agents = int(row['agents'])
        solved = int(row['solved'])
        rt = float(row['runtime_ms']) if solved else None
        algos[alg].append((agents, rt))

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = {'CBS': '#e74c3c', 'EECBS': '#3498db', 'PP': '#2ecc71'}
    markers = {'CBS': 'o', 'EECBS': 's', 'PP': '^'}

    for alg, pts in algos.items():
        pts.sort()
        x_solved = [p[0] for p in pts if p[1] is not None and p[1] > 0]
        y_solved = [p[1] for p in pts if p[1] is not None and p[1] > 0]
        x_fail = [p[0] for p in pts if p[1] is None or p[1] <= 0]

        ax.plot(x_solved, y_solved, f'-{markers[alg]}', color=colors[alg],
                label=alg, markersize=8, linewidth=2)
        if x_fail:
            ax.scatter(x_fail, [30000]*len(x_fail), marker='x', color=colors[alg],
                      s=100, zorder=5)

    ax.set_yscale('log')
    ax.set_xlabel('Number of Agents', fontsize=13)
    ax.set_ylabel('Runtime (ms, log scale)', fontsize=13)
    ax.set_title('Algorithm Scalability: Runtime vs. Number of Agents', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.1, 100000)
    ax.axhline(y=30000, color='gray', linestyle='--', alpha=0.5, label='Timeout (30s)')
    plt.savefig(f'{FIGURES_DIR}/scalability_runtime.png')
    plt.close()
    print("  Created scalability_runtime.png")

# ============================================================
# Figure 2: Solution Quality (Cost) vs Agents
# ============================================================
def plot_solution_quality():
    data = read_csv('benchmarks/scalability.csv')
    algos = {'CBS': [], 'EECBS': [], 'PP': []}
    for row in data:
        if row['mapType'] != 'random': continue
        alg = row['algorithm']
        if alg not in algos: continue
        if int(row['solved']):
            algos[alg].append((int(row['agents']), int(row['totalCost'])))

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = {'CBS': '#e74c3c', 'EECBS': '#3498db', 'PP': '#2ecc71'}
    for alg, pts in algos.items():
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], '-o', color=colors[alg],
                label=alg, markersize=7, linewidth=2)

    ax.set_xlabel('Number of Agents', fontsize=13)
    ax.set_ylabel('Total Path Cost (Sum of Costs)', fontsize=13)
    ax.set_title('Solution Quality: Total Cost vs. Number of Agents', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.savefig(f'{FIGURES_DIR}/solution_quality.png')
    plt.close()
    print("  Created solution_quality.png")

# ============================================================
# Figure 3: Map Type Comparison
# ============================================================
def plot_map_comparison():
    data = read_csv('benchmarks/scalability.csv')
    map_types = ['empty', 'random10', 'random20', 'warehouse']
    algos_list = ['CBS', 'EECBS', 'PP']

    results = {}
    for row in data:
        if row['mapType'] not in map_types: continue
        key = (row['algorithm'], row['mapType'])
        if int(row['solved']):
            results[key] = float(row['runtime_ms'])
        else:
            results[key] = None

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(map_types))
    width = 0.25
    colors = ['#e74c3c', '#3498db', '#2ecc71']

    for i, alg in enumerate(algos_list):
        vals = []
        for mt in map_types:
            v = results.get((alg, mt))
            vals.append(v if v and v > 0 else 0.01)
        ax.bar(x + i*width, vals, width, label=alg, color=colors[i], alpha=0.85)

    ax.set_yscale('log')
    ax.set_xlabel('Map Type', fontsize=13)
    ax.set_ylabel('Runtime (ms, log scale)', fontsize=13)
    ax.set_title('Runtime Comparison Across Map Types (20 agents)', fontsize=14)
    ax.set_xticks(x + width)
    ax.set_xticklabels(['Empty', 'Random 10%', 'Random 20%', 'Warehouse'])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    plt.savefig(f'{FIGURES_DIR}/map_comparison.png')
    plt.close()
    print("  Created map_comparison.png")

# ============================================================
# Figure 4: Suboptimality Analysis
# ============================================================
def plot_suboptimality():
    data = read_csv('benchmarks/suboptimality.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    agent_groups = {}
    for row in data:
        a = int(row['agents'])
        if a not in agent_groups: agent_groups[a] = []
        agent_groups[a].append(row)

    colors = {5: '#e74c3c', 10: '#3498db', 15: '#2ecc71', 20: '#9b59b6'}

    for a, rows in sorted(agent_groups.items()):
        ws = [float(r['w']) for r in rows if float(r['totalCost']) > 0]
        rts = [float(r['runtime_ms']) for r in rows if float(r['totalCost']) > 0]
        ratios = [float(r['ratio']) for r in rows if float(r['ratio']) > 0]
        ws_r = [float(r['w']) for r in rows if float(r['ratio']) > 0]

        if ws and rts:
            ax1.plot(ws, rts, '-o', color=colors[a], label=f'{a} agents', markersize=7)
        if ws_r and ratios:
            ax2.plot(ws_r, ratios, '-s', color=colors[a], label=f'{a} agents', markersize=7)

    ax1.set_xlabel('Suboptimality Bound (w)', fontsize=12)
    ax1.set_ylabel('Runtime (ms)', fontsize=12)
    ax1.set_title('EECBS: Runtime vs. w', fontsize=13)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('Suboptimality Bound (w)', fontsize=12)
    ax2.set_ylabel('Cost Ratio (actual/optimal)', fontsize=12)
    ax2.set_title('EECBS: Solution Quality vs. w', fontsize=13)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0.95, 1.15)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/suboptimality_analysis.png')
    plt.close()
    print("  Created suboptimality_analysis.png")

# ============================================================
# Figure 5: Lifelong MAPF
# ============================================================
def plot_lifelong():
    data = read_csv('benchmarks/lifelong.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    agents = [int(r['agents']) for r in data]
    completed = [int(r['completed']) for r in data]
    tasks = [int(r['tasks']) for r in data]
    service = [float(r['avgServiceTime']) for r in data]
    runtime = [float(r['runtime_ms']) for r in data]
    completion_rate = [c/t*100 for c, t in zip(completed, tasks)]

    ax1.bar(range(len(agents)), completion_rate, color='#3498db', alpha=0.85)
    ax1.set_xticks(range(len(agents)))
    ax1.set_xticklabels(agents)
    ax1.set_xlabel('Number of Agents', fontsize=12)
    ax1.set_ylabel('Task Completion Rate (%)', fontsize=12)
    ax1.set_title('Lifelong MAPF: Task Completion Rate', fontsize=13)
    ax1.grid(True, alpha=0.3, axis='y')

    ax2.plot(agents, service, '-o', color='#e74c3c', markersize=8, linewidth=2)
    ax2.set_xlabel('Number of Agents', fontsize=12)
    ax2.set_ylabel('Average Service Time (steps/task)', fontsize=12)
    ax2.set_title('Lifelong MAPF: Service Efficiency', fontsize=13)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/lifelong_mapf.png')
    plt.close()
    print("  Created lifelong_mapf.png")

# ============================================================
# Figure 6: Distributed MAPF - Communication Impact
# ============================================================
def plot_distributed():
    data = read_csv('benchmarks/distributed.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Conflicts vs comm radius for different drop rates (50 agents)
    for drop in [0.0, 0.1, 0.3]:
        subset = [r for r in data if int(r['agents']) == 50 and float(r['dropRate']) == drop]
        radii = [float(r['commRadius']) for r in subset]
        conflicts = [int(r['conflicts']) for r in subset]
        ax1.plot(radii, conflicts, '-o', label=f'Drop={drop}', markersize=8, linewidth=2)

    ax1.set_xlabel('Communication Radius', fontsize=12)
    ax1.set_ylabel('Number of Conflicts', fontsize=12)
    ax1.set_title('Distributed MAPF: Conflicts vs. Comm. Radius\n(50 agents)', fontsize=13)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Avg messages vs agents
    for radius in [5.0, 10.0, 20.0, 50.0]:
        subset = [r for r in data if float(r['commRadius']) == radius and float(r['dropRate']) == 0.0]
        ag = [int(r['agents']) for r in subset]
        msgs = [float(r['avgMsgs']) for r in subset]
        ax2.plot(ag, msgs, '-s', label=f'r={int(radius)}', markersize=7, linewidth=2)

    ax2.set_xlabel('Number of Agents', fontsize=12)
    ax2.set_ylabel('Avg Messages per Step', fontsize=12)
    ax2.set_title('Communication Load vs. Agents', fontsize=13)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/distributed_mapf.png')
    plt.close()
    print("  Created distributed_mapf.png")

# ============================================================
# Figure 7: Warehouse Large-Scale
# ============================================================
def plot_warehouse():
    data = read_csv('benchmarks/warehouse_large.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    for ms in [32, 64]:
        subset = [r for r in data if int(r['mapSize']) == ms and int(r['solved'])]
        agents = [int(r['agents']) for r in subset]
        runtime = [float(r['runtime_ms']) for r in subset]
        cost = [int(r['totalCost']) for r in subset]

        ax1.plot(agents, runtime, '-o', label=f'{ms}×{ms} map', markersize=8, linewidth=2)
        ax2.plot(agents, cost, '-s', label=f'{ms}×{ms} map', markersize=8, linewidth=2)

    ax1.set_xlabel('Number of Agents', fontsize=12)
    ax1.set_ylabel('Runtime (ms)', fontsize=12)
    ax1.set_title('Warehouse Scale: PP Runtime', fontsize=13)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('Number of Agents', fontsize=12)
    ax2.set_ylabel('Total Path Cost', fontsize=12)
    ax2.set_title('Warehouse Scale: Solution Cost', fontsize=13)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/warehouse_scale.png')
    plt.close()
    print("  Created warehouse_scale.png")

# ============================================================
# Figure 8: Algorithm Overview Heatmap
# ============================================================
def plot_overview_heatmap():
    fig, ax = plt.subplots(figsize=(9, 5))

    algorithms = ['CBS', 'EECBS\n(w=1.5)', 'PP', 'LaCAM']
    metrics = ['Optimality', 'Speed\n(small)', 'Speed\n(large)', 'Completeness', 'Scalability']

    # Scores 0-5
    scores = np.array([
        [5, 3, 1, 5, 1],  # CBS
        [4, 4, 2, 4, 2],  # EECBS
        [2, 5, 5, 3, 5],  # PP
        [1, 5, 4, 2, 4],  # LaCAM
    ])

    im = ax.imshow(scores, cmap='RdYlGn', aspect='auto', vmin=0, vmax=5)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, fontsize=11)
    ax.set_yticks(range(len(algorithms)))
    ax.set_yticklabels(algorithms, fontsize=11)

    for i in range(len(algorithms)):
        for j in range(len(metrics)):
            ax.text(j, i, str(scores[i, j]), ha='center', va='center',
                   fontsize=14, fontweight='bold', color='black')

    ax.set_title('Algorithm Comparison Overview (1=Low, 5=High)', fontsize=14)
    plt.colorbar(im, ax=ax, shrink=0.8)
    plt.savefig(f'{FIGURES_DIR}/algorithm_overview.png')
    plt.close()
    print("  Created algorithm_overview.png")

if __name__ == '__main__':
    print("Generating MAPF benchmark figures...")
    plot_scalability()
    plot_solution_quality()
    plot_map_comparison()
    plot_suboptimality()
    plot_lifelong()
    plot_distributed()
    plot_warehouse()
    plot_overview_heatmap()
    print("All figures generated successfully!")
