"""
Module 1: Research Network Structure Analysis
Generates and analyzes co-authorship and citation networks.
"""

import networkx as nx
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
import os

np.random.seed(42)


def generate_coauthorship_network(n_researchers=200, n_communities=5):
    """Generate a co-authorship network using stochastic block model."""
    sizes = [n_researchers // n_communities] * n_communities
    sizes[-1] += n_researchers - sum(sizes)

    p_matrix = np.full((n_communities, n_communities), 0.02)
    np.fill_diagonal(p_matrix, 0.15)

    G = nx.stochastic_block_model(sizes, p_matrix.tolist(), seed=42)

    fields = ["Physics", "Biology", "CS", "Chemistry", "Medicine"]
    genders = ["M", "F", "Other"]
    regions = ["Asia", "Europe", "NorthAmerica", "LatinAmerica", "Africa"]

    for node in G.nodes():
        block = G.nodes[node].get("block", 0)
        G.nodes[node]["field"] = fields[block % len(fields)]
        G.nodes[node]["gender"] = np.random.choice(genders, p=[0.60, 0.35, 0.05])
        G.nodes[node]["region"] = np.random.choice(regions, p=[0.30, 0.30, 0.25, 0.10, 0.05])
        G.nodes[node]["career_stage"] = np.random.choice(
            ["early", "mid", "senior"], p=[0.40, 0.35, 0.25]
        )
        G.nodes[node]["h_index"] = max(0, int(np.random.lognormal(2.0, 0.8)))
        G.nodes[node]["productivity"] = max(1, int(np.random.lognormal(1.5, 0.7)))

    return G


def generate_citation_network(coauth_G):
    """Generate a directed citation network based on co-authorship."""
    C = nx.DiGraph()
    nodes = list(coauth_G.nodes())
    C.add_nodes_from(nodes)

    for u in nodes:
        n_citations = max(1, int(np.random.exponential(3)))
        potential = [v for v in nodes if v != u]
        weights = []
        for v in potential:
            w = 1.0
            if coauth_G.has_edge(u, v):
                w += 5.0
            w += coauth_G.nodes[v]["h_index"] * 0.1
            if coauth_G.nodes[v]["field"] == coauth_G.nodes[u]["field"]:
                w += 2.0
            weights.append(w)
        weights = np.array(weights)
        weights /= weights.sum()

        targets = np.random.choice(potential, size=min(n_citations, len(potential)),
                                   replace=False, p=weights)
        for t in targets:
            C.add_edge(u, t)

    return C


def analyze_network(G, label="coauthorship"):
    """Compute network metrics."""
    metrics = {}
    metrics["n_nodes"] = G.number_of_nodes()
    metrics["n_edges"] = G.number_of_edges()

    if G.is_directed():
        metrics["density"] = nx.density(G)
        in_deg = dict(G.in_degree())
        metrics["avg_in_degree"] = np.mean(list(in_deg.values()))
        pr = nx.pagerank(G)
        metrics["top_pagerank"] = sorted(pr.items(), key=lambda x: -x[1])[:10]
        metrics["pagerank_gini"] = _gini(list(pr.values()))
    else:
        metrics["density"] = nx.density(G)
        deg = dict(G.degree())
        metrics["avg_degree"] = np.mean(list(deg.values()))
        metrics["clustering_coefficient"] = nx.average_clustering(G)

        communities = nx.community.greedy_modularity_communities(G)
        metrics["n_communities"] = len(communities)
        metrics["modularity"] = nx.community.modularity(G, communities)

        bc = nx.betweenness_centrality(G)
        metrics["top_betweenness"] = sorted(bc.items(), key=lambda x: -x[1])[:10]
        metrics["betweenness_gini"] = _gini(list(bc.values()))

    return metrics


def _gini(values):
    """Compute Gini coefficient."""
    values = np.array(sorted(values))
    n = len(values)
    if n == 0 or values.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * values) - (n + 1) * np.sum(values)) / (n * np.sum(values))


def plot_network(G, filename, title="Network"):
    """Plot network with community coloring."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    if not G.is_directed():
        communities = nx.community.greedy_modularity_communities(G)
        color_map = {}
        cmap = plt.cm.viridis
        for i, comm in enumerate(communities):
            for node in comm:
                color_map[node] = cmap(i / max(1, len(communities) - 1))
        colors = [color_map.get(n, (0.5, 0.5, 0.5, 1.0)) for n in G.nodes()]
    else:
        colors = "steelblue"

    pos = nx.spring_layout(G, seed=42, k=0.3)
    nx.draw_networkx_nodes(G, pos, node_size=20, node_color=colors, alpha=0.7, ax=ax)
    nx.draw_networkx_edges(G, pos, alpha=0.1, arrows=G.is_directed(), ax=ax)
    ax.set_title(title, fontsize=14)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_degree_distribution(G, filename, title="Degree Distribution"):
    """Plot degree distribution (log-log)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if G.is_directed():
        degrees = [d for _, d in G.in_degree()]
        ax.set_xlabel("In-Degree", fontsize=12)
    else:
        degrees = [d for _, d in G.degree()]
        ax.set_xlabel("Degree", fontsize=12)

    unique, counts = np.unique(degrees, return_counts=True)
    ax.loglog(unique + 1, counts, "o", markersize=5, alpha=0.7)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_network_analysis(output_dir="figures", results_dir="results"):
    """Run full network analysis pipeline."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    coauth_G = generate_coauthorship_network(n_researchers=200)
    citation_G = generate_citation_network(coauth_G)

    coauth_metrics = analyze_network(coauth_G, "coauthorship")
    citation_metrics = analyze_network(citation_G, "citation")

    plot_network(coauth_G, f"{output_dir}/fig1_coauthorship_network.png",
                 "Co-authorship Network (SBM, N=200)")
    plot_network(citation_G, f"{output_dir}/fig2_citation_network.png",
                 "Citation Network (N=200)")
    plot_degree_distribution(coauth_G, f"{output_dir}/fig3_coauth_degree_dist.png",
                             "Co-authorship Degree Distribution")
    plot_degree_distribution(citation_G, f"{output_dir}/fig4_citation_degree_dist.png",
                             "Citation In-Degree Distribution")

    serializable_coauth = {k: v for k, v in coauth_metrics.items()
                           if not isinstance(v, list)}
    serializable_citation = {k: v for k, v in citation_metrics.items()
                             if not isinstance(v, list)}

    results = {"coauthorship": serializable_coauth, "citation": serializable_citation}

    with open(f"{results_dir}/network_metrics.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print("Network analysis complete.")
    print(f"  Co-authorship: {coauth_metrics['n_nodes']} nodes, {coauth_metrics['n_edges']} edges")
    print(f"  Modularity: {coauth_metrics.get('modularity', 'N/A'):.3f}")
    print(f"  Betweenness Gini: {coauth_metrics.get('betweenness_gini', 'N/A'):.3f}")
    print(f"  Citation: {citation_metrics['n_nodes']} nodes, {citation_metrics['n_edges']} edges")
    print(f"  PageRank Gini: {citation_metrics.get('pagerank_gini', 'N/A'):.3f}")

    return coauth_G, citation_G, results


if __name__ == "__main__":
    run_network_analysis()
