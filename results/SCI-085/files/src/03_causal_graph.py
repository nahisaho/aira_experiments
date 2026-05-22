"""
Module 3: Causal Graph Estimation from Perturbation Effects
============================================================
- Perturbation → gene regulatory network inference
- PC algorithm (constraint-based) for causal DAG
- Bootstrap stability assessment
- Network visualization
"""

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations
from scipy.stats import pearsonr, spearmanr
from scipy.stats import norm
import os
import json
import warnings
warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)


def build_perturbation_effect_matrix(adata, de_path="results/02_de_results.csv",
                                      n_top_genes=50):
    """
    Build perturbation × gene effect matrix from DE results.
    Each row = perturbation, each column = gene, value = log2FC.
    """
    de_df = pd.read_csv(de_path)

    # Select top DE genes across all perturbations
    top_genes = (de_df[de_df["significant"]]
                 .groupby("gene")["p_adj"]
                 .min()
                 .nsmallest(n_top_genes)
                 .index.tolist())

    # Build effect matrix
    perts = de_df["perturbation"].unique()
    effect_matrix = pd.DataFrame(0.0, index=perts, columns=top_genes)

    for _, row in de_df[de_df["gene"].isin(top_genes)].iterrows():
        effect_matrix.loc[row["perturbation"], row["gene"]] = row["log2FC"]

    return effect_matrix


def partial_correlation(X, i, j, S):
    """
    Compute partial correlation between variables i and j given set S.
    Uses recursive formula for computational efficiency.
    """
    if len(S) == 0:
        r, _ = pearsonr(X[:, i], X[:, j])
        return r

    # Residualize i and j on S
    from sklearn.linear_model import LinearRegression
    reg_i = LinearRegression().fit(X[:, list(S)], X[:, i])
    reg_j = LinearRegression().fit(X[:, list(S)], X[:, j])
    res_i = X[:, i] - reg_i.predict(X[:, list(S)])
    res_j = X[:, j] - reg_j.predict(X[:, list(S)])

    if np.std(res_i) < 1e-10 or np.std(res_j) < 1e-10:
        return 0.0

    r, _ = pearsonr(res_i, res_j)
    return r


def pc_algorithm(data, alpha=0.05, max_cond_size=3):
    """
    PC algorithm for causal DAG estimation.
    data: np.array (samples × variables)
    Returns: adjacency matrix and separation sets.
    """
    n_vars = data.shape[1]
    n_samples = data.shape[0]

    # Start with complete undirected graph
    adj = np.ones((n_vars, n_vars), dtype=bool)
    np.fill_diagonal(adj, False)
    sep_sets = {(i, j): set() for i in range(n_vars) for j in range(n_vars)}

    # Edge removal by conditional independence testing
    for cond_size in range(max_cond_size + 1):
        for i in range(n_vars):
            for j in range(i + 1, n_vars):
                if not adj[i, j]:
                    continue

                # Get neighbors of i (excluding j)
                neighbors_i = [k for k in range(n_vars) if adj[i, k] and k != j]

                if len(neighbors_i) < cond_size:
                    continue

                # Test all conditioning sets of given size
                from itertools import combinations as combins
                for S in combins(neighbors_i, cond_size):
                    S = set(S)
                    pcor = partial_correlation(data, i, j, S)

                    # Fisher's z-test
                    z = 0.5 * np.log((1 + pcor + 1e-10) / (1 - pcor + 1e-10))
                    z_stat = np.sqrt(n_samples - len(S) - 3) * abs(z)
                    p_value = 2 * (1 - norm.cdf(z_stat))

                    if p_value > alpha:
                        adj[i, j] = False
                        adj[j, i] = False
                        sep_sets[(i, j)] = S
                        sep_sets[(j, i)] = S
                        break

    return adj, sep_sets


def orient_edges(adj, sep_sets, n_vars):
    """
    Orient edges using Meek's rules (v-structure detection).
    Returns directed adjacency matrix.
    """
    directed = adj.copy().astype(int)

    # Rule 1: V-structure orientation
    # If i - k - j and k not in sep(i,j), orient as i → k ← j
    for k in range(n_vars):
        neighbors = [n for n in range(n_vars) if adj[k, n]]
        for i, j in combinations(neighbors, 2):
            if not adj[i, j] and not adj[j, i]:  # i and j not adjacent
                if k not in sep_sets.get((i, j), set()):
                    # Orient i → k ← j
                    directed[k, i] = 0  # remove k → i
                    directed[k, j] = 0  # remove k → j

    return directed


def bootstrap_stability(data, n_bootstrap=100, alpha=0.05):
    """
    Assess edge stability via bootstrap resampling.
    Returns edge frequency matrix.
    """
    n_samples, n_vars = data.shape
    edge_freq = np.zeros((n_vars, n_vars))

    for b in range(n_bootstrap):
        boot_idx = np.random.choice(n_samples, size=n_samples, replace=True)
        boot_data = data[boot_idx]

        try:
            adj, sep_sets = pc_algorithm(boot_data, alpha=alpha, max_cond_size=2)
            directed = orient_edges(adj, sep_sets, n_vars)
            edge_freq += (directed > 0).astype(float)
        except Exception:
            continue

    edge_freq /= n_bootstrap
    return edge_freq


def plot_causal_graph(adj_matrix, var_names, edge_freq=None, out_dir="figures"):
    """Visualize causal graph."""
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    # 1. Causal DAG
    ax = axes[0]
    G = nx.DiGraph()
    for i, name in enumerate(var_names):
        G.add_node(name)
    for i in range(len(var_names)):
        for j in range(len(var_names)):
            if adj_matrix[i, j] > 0 and i != j:
                G.add_edge(var_names[i], var_names[j])

    if len(G.edges()) > 0:
        pos = nx.spring_layout(G, seed=SEED, k=2)
        nx.draw(G, pos, ax=ax, with_labels=True, node_color="lightblue",
                node_size=600, font_size=7, arrows=True,
                arrowsize=15, edge_color="gray", width=1.5)
    ax.set_title(f"Causal DAG (PC Algorithm)\n{len(G.nodes())} nodes, {len(G.edges())} edges")

    # 2. Edge stability heatmap
    ax = axes[1]
    if edge_freq is not None:
        n_show = min(30, len(var_names))
        im = ax.imshow(edge_freq[:n_show, :n_show], cmap="YlOrRd", aspect="auto",
                       vmin=0, vmax=1)
        plt.colorbar(im, ax=ax, label="Bootstrap Frequency")
        ax.set_xticks(range(n_show))
        ax.set_yticks(range(n_show))
        ax.set_xticklabels(var_names[:n_show], rotation=90, fontsize=6)
        ax.set_yticklabels(var_names[:n_show], fontsize=6)
    ax.set_title("Edge Stability (Bootstrap)")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "03_causal_graph.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, "03_causal_graph.svg"), bbox_inches="tight")
    plt.close()
    print(f"✓ Causal graph plots saved")


def run_causal_pipeline(adata_path="data/perturbseq_processed.h5ad"):
    """Main causal inference pipeline."""
    adata = ad.read_h5ad(adata_path)
    print(f"Loaded: {adata.shape}")

    # Build effect matrix
    print("Building perturbation effect matrix...")
    effect_matrix = build_perturbation_effect_matrix(adata, n_top_genes=30)
    print(f"  Effect matrix: {effect_matrix.shape}")

    # Run PC algorithm
    print("Running PC algorithm for causal graph estimation...")
    data = effect_matrix.values
    adj, sep_sets = pc_algorithm(data, alpha=0.1, max_cond_size=2)
    n_vars = data.shape[1]
    directed = orient_edges(adj, sep_sets, n_vars)

    # Bootstrap stability
    print("Assessing edge stability (bootstrap)...")
    edge_freq = bootstrap_stability(data, n_bootstrap=50, alpha=0.1)

    # Plot
    var_names = effect_matrix.columns.tolist()
    plot_causal_graph(directed, var_names, edge_freq)

    # Save results
    os.makedirs("results", exist_ok=True)

    G = nx.DiGraph()
    for i, name in enumerate(var_names):
        G.add_node(name)
    for i in range(len(var_names)):
        for j in range(len(var_names)):
            if directed[i, j] > 0 and i != j:
                G.add_edge(var_names[i], var_names[j],
                           stability=float(edge_freq[i, j]))

    edge_list = [(u, v, d.get("stability", 0)) for u, v, d in G.edges(data=True)]
    edge_df = pd.DataFrame(edge_list, columns=["source", "target", "bootstrap_freq"])
    edge_df.to_csv("results/03_causal_edges.csv", index=False)

    summary = {
        "n_genes_in_graph": len(var_names),
        "n_edges": len(G.edges()),
        "n_nodes_with_edges": len([n for n in G.nodes() if G.degree(n) > 0]),
        "mean_stability": float(edge_freq[edge_freq > 0].mean()) if edge_freq.any() else 0,
        "stable_edges_gt50": int((edge_freq > 0.5).sum()),
        "graph_density": float(nx.density(G)) if len(G.nodes()) > 1 else 0,
    }
    with open("results/03_causal_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Save adjacency
    pd.DataFrame(directed, index=var_names, columns=var_names).to_csv(
        "results/03_adjacency_matrix.csv")
    pd.DataFrame(edge_freq, index=var_names, columns=var_names).to_csv(
        "results/03_edge_stability.csv")

    print(f"✓ Causal graph complete:")
    print(f"  Nodes: {summary['n_genes_in_graph']}, Edges: {summary['n_edges']}")
    print(f"  Stable edges (>50%): {summary['stable_edges_gt50']}")
    print(f"  Graph density: {summary['graph_density']:.3f}")

    return G, summary


if __name__ == "__main__":
    run_causal_pipeline()
