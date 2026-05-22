"""
Module 5: Membrane Self-Organization and Protocell Formation
Agent-based model of amphiphile self-assembly into vesicles / protocells.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

np.random.seed(101)

class Amphiphile:
    __slots__ = ['x', 'y', 'angle', 'cluster_id']
    def __init__(self, x, y, angle):
        self.x = x
        self.y = y
        self.angle = angle
        self.cluster_id = -1

def simulate_protocell(
    n_amphiphiles=500,
    box_size=50.0,
    n_steps=1000,
    interaction_radius=3.0,
    alignment_strength=0.3,
    attraction_strength=0.5,
    noise=0.3,
    critical_aggregation=8,
    fatty_acid_synthesis_rate=0.002,
    encapsulation_prob=0.05,
):
    """
    2D self-assembly simulation of amphiphiles forming micelles/vesicles.
    Uses Vicsek-like alignment + short-range attraction.
    """
    # Initialise
    amps = []
    for _ in range(n_amphiphiles):
        a = Amphiphile(
            np.random.uniform(0, box_size),
            np.random.uniform(0, box_size),
            np.random.uniform(0, 2 * np.pi)
        )
        amps.append(a)

    history = {
        'step': [], 'n_clusters': [], 'max_cluster': [],
        'n_vesicles': [], 'n_amphiphiles': [],
        'mean_cluster_size': [],
    }
    snapshots = []

    for step in range(n_steps):
        n = len(amps)
        xs = np.array([a.x for a in amps])
        ys = np.array([a.y for a in amps])
        angles = np.array([a.angle for a in amps])

        # Pairwise distances (periodic BC)
        dx = xs[:, None] - xs[None, :]
        dy = ys[:, None] - ys[None, :]
        dx = dx - box_size * np.round(dx / box_size)
        dy = dy - box_size * np.round(dy / box_size)
        dist = np.sqrt(dx**2 + dy**2)

        # Find neighbours
        neighbours = dist < interaction_radius
        np.fill_diagonal(neighbours, False)

        # Alignment + attraction
        new_angles = np.copy(angles)
        new_xs = np.copy(xs)
        new_ys = np.copy(ys)

        for i in range(n):
            nb_idx = np.where(neighbours[i])[0]
            if len(nb_idx) > 0:
                # Alignment
                mean_angle = np.arctan2(
                    np.mean(np.sin(angles[nb_idx])),
                    np.mean(np.cos(angles[nb_idx]))
                )
                new_angles[i] = (1 - alignment_strength) * angles[i] + \
                                alignment_strength * mean_angle

                # Attraction toward centre of neighbours
                cx = np.mean(dx[i, nb_idx])
                cy = np.mean(dy[i, nb_idx])
                d = np.sqrt(cx**2 + cy**2) + 1e-10
                new_xs[i] += attraction_strength * cx / d
                new_ys[i] += attraction_strength * cy / d

            # Noise + movement
            new_angles[i] += noise * np.random.uniform(-np.pi, np.pi)
            speed = 0.5
            new_xs[i] += speed * np.cos(new_angles[i])
            new_ys[i] += speed * np.sin(new_angles[i])

        # Periodic BC
        new_xs %= box_size
        new_ys %= box_size

        for i in range(n):
            amps[i].x = new_xs[i]
            amps[i].y = new_ys[i]
            amps[i].angle = new_angles[i]

        # New amphiphile synthesis
        n_new = np.random.poisson(fatty_acid_synthesis_rate * n)
        for _ in range(n_new):
            amps.append(Amphiphile(
                np.random.uniform(0, box_size),
                np.random.uniform(0, box_size),
                np.random.uniform(0, 2 * np.pi)
            ))

        # Clustering (simple distance-based)
        if step % 20 == 0:
            clusters = cluster_amphiphiles(amps, interaction_radius * 1.5, box_size)
            cluster_sizes = [len(c) for c in clusters if len(c) >= 2]
            n_vesicles = sum(1 for s in cluster_sizes if s >= critical_aggregation)

            history['step'].append(step)
            history['n_clusters'].append(len(cluster_sizes))
            history['max_cluster'].append(max(cluster_sizes) if cluster_sizes else 0)
            history['n_vesicles'].append(n_vesicles)
            history['n_amphiphiles'].append(len(amps))
            history['mean_cluster_size'].append(
                np.mean(cluster_sizes) if cluster_sizes else 0)

            if step in [0, n_steps // 4, n_steps // 2, n_steps - 1]:
                snapshots.append({
                    'step': step,
                    'xs': [a.x for a in amps],
                    'ys': [a.y for a in amps],
                    'cluster_ids': [a.cluster_id for a in amps],
                })

    return history, snapshots

def cluster_amphiphiles(amps, radius, box_size):
    """Simple single-linkage clustering."""
    n = len(amps)
    visited = [False] * n
    clusters = []

    xs = np.array([a.x for a in amps])
    ys = np.array([a.y for a in amps])

    for i in range(n):
        if visited[i]:
            continue
        cluster = []
        stack = [i]
        while stack:
            j = stack.pop()
            if visited[j]:
                continue
            visited[j] = True
            cluster.append(j)
            amps[j].cluster_id = len(clusters)
            # Find neighbours
            dx = xs - xs[j]
            dy = ys - ys[j]
            dx = dx - box_size * np.round(dx / box_size)
            dy = dy - box_size * np.round(dy / box_size)
            dist = np.sqrt(dx**2 + dy**2)
            for k in range(n):
                if not visited[k] and dist[k] < radius:
                    stack.append(k)
        clusters.append(cluster)
    return clusters

def plot_protocell(history, snapshots, box_size=50.0):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Protocell Formation: Amphiphile Self-Assembly', fontsize=14)

    # Snapshots
    for idx, (ax, snap) in enumerate(zip(axes[0], snapshots[:3])):
        cids = np.array(snap['cluster_ids'])
        ax.scatter(snap['xs'], snap['ys'], c=cids % 20, cmap='tab20',
                   s=3, alpha=0.6)
        ax.set_xlim(0, box_size)
        ax.set_ylim(0, box_size)
        ax.set_title(f'Step {snap["step"]}')
        ax.set_aspect('equal')

    ax = axes[1, 0]
    ax.plot(history['step'], history['n_vesicles'], color='#E91E63')
    ax.set_xlabel('Step')
    ax.set_ylabel('Count')
    ax.set_title('D) Vesicle Count (≥8 amphiphiles)')

    ax = axes[1, 1]
    ax.plot(history['step'], history['max_cluster'], color='#2196F3')
    ax.set_xlabel('Step')
    ax.set_ylabel('Size')
    ax.set_title('E) Largest Cluster Size')

    ax = axes[1, 2]
    ax.plot(history['step'], history['n_amphiphiles'], color='#4CAF50')
    ax.set_xlabel('Step')
    ax.set_ylabel('Count')
    ax.set_title('F) Total Amphiphile Population')

    plt.tight_layout()
    plt.savefig('figures/fig8_protocell.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/fig8_protocell.svg', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("Running Protocell simulation...")
    history, snapshots = simulate_protocell()

    plot_protocell(history, snapshots)

    results = {
        'final_vesicles': history['n_vesicles'][-1] if history['n_vesicles'] else 0,
        'final_max_cluster': history['max_cluster'][-1] if history['max_cluster'] else 0,
        'final_amphiphiles': history['n_amphiphiles'][-1] if history['n_amphiphiles'] else 0,
        'final_mean_cluster_size': history['mean_cluster_size'][-1] if history['mean_cluster_size'] else 0,
    }
    with open('results/protocell_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Module 5 complete. Vesicles={results['final_vesicles']}, "
          f"Max cluster={results['final_max_cluster']}")
