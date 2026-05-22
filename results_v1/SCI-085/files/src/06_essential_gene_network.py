"""
Module 6: Essential Gene Network Estimation — Case Study
=========================================================
- Identify essential gene candidates from perturbation fitness effects
- Build essentiality network from gene-gene interactions
- Community detection for functional modules
- Comparison with known essential gene databases
"""

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from scipy.stats import zscore
from sklearn.metrics.pairwise import cosine_similarity
import os
import json
import warnings
warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)

# Known essential gene pathways (simulated reference)
ESSENTIAL_PATHWAYS = {
    "DNA_replication": ["Gene_0", "Gene_1", "Gene_2", "Gene_3", "Gene_4"],
    "transcription": ["Gene_10", "Gene_11", "Gene_12", "Gene_13", "Gene_14"],
    "translation": ["Gene_20", "Gene_21", "Gene_22", "Gene_23", "Gene_24"],
    "proteasome": ["Gene_30", "Gene_31", "Gene_32", "Gene_33", "Gene_34"],
    "metabolism": ["Gene_40", "Gene_41", "Gene_42", "Gene_43", "Gene_44"],
}


def compute_fitness_effects(adata, control_key="non-targeting"):
    """
    Estimate fitness effects from perturbation data.
    Proxy: cells with essential gene perturbations have fewer UMIs / distinct genes.
    """
    if "log_normalized" in adata.layers:
        X = adata.layers["log_normalized"]
    else:
        X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    ctrl_mask = adata.obs["perturbation"] == control_key
    ctrl_mean = X[ctrl_mask].mean(axis=0)
    ctrl_total_umi = np.asarray(adata[ctrl_mask].X.sum(axis=1)).flatten().mean()

    fitness_scores = {}
    perturbation_effects = {}

    for pert in adata.obs["perturbation"].unique():
        if pert == control_key or "|" in pert:
            continue
        mask = adata.obs["perturbation"] == pert
        if mask.sum() < 5:
            continue

        pert_expr = X[mask]
        pert_total_umi = np.asarray(adata[mask].X.sum(axis=1)).flatten().mean()

        # Fitness proxy: UMI ratio (lower = more essential)
        fitness_ratio = pert_total_umi / (ctrl_total_umi + 1e-10)

        # Expression effect magnitude
        effect = pert_expr.mean(axis=0) - ctrl_mean
        effect_magnitude = np.linalg.norm(effect)

        # DE gene count
        n_de_genes = np.sum(np.abs(effect) > np.std(effect) * 1.5)

        fitness_scores[pert] = {
            "fitness_ratio": float(fitness_ratio),
            "effect_magnitude": float(effect_magnitude),
            "n_de_genes": int(n_de_genes),
            "n_cells": int(mask.sum()),
            "mean_umi": float(pert_total_umi),
        }
        perturbation_effects[pert] = effect

    return fitness_scores, perturbation_effects


def build_essentiality_network(perturbation_effects, fitness_scores,
                                sim_threshold=0.3, min_effect=0.5):
    """
    Build gene interaction network based on perturbation effect similarity.
    Edges connect perturbations with correlated transcriptional responses.
    """
    perts = list(perturbation_effects.keys())
    effects_matrix = np.array([perturbation_effects[p] for p in perts])

    # Cosine similarity between perturbation effects
    sim = cosine_similarity(effects_matrix)

    # Build graph
    G = nx.Graph()
    for i, p in enumerate(perts):
        fs = fitness_scores[p]
        essentiality_score = (1 - fs["fitness_ratio"]) * fs["effect_magnitude"]
        G.add_node(p,
                   fitness_ratio=fs["fitness_ratio"],
                   effect_magnitude=fs["effect_magnitude"],
                   essentiality_score=essentiality_score,
                   n_de_genes=fs["n_de_genes"])

    for i in range(len(perts)):
        for j in range(i + 1, len(perts)):
            if abs(sim[i, j]) > sim_threshold:
                G.add_edge(perts[i], perts[j], weight=float(abs(sim[i, j])))

    return G, sim


def community_detection(G):
    """Louvain-style community detection using greedy modularity."""
    if len(G.edges()) == 0:
        return {n: 0 for n in G.nodes()}

    communities = nx.community.greedy_modularity_communities(G)
    community_map = {}
    for idx, comm in enumerate(communities):
        for node in comm:
            community_map[node] = idx
    return community_map


def annotate_essential_genes(G, fitness_scores, z_threshold=-1.5):
    """
    Classify genes as essential based on fitness score z-score.
    """
    ratios = [fitness_scores[n]["fitness_ratio"] for n in G.nodes()]
    z_scores = zscore(ratios) if len(ratios) > 2 else np.zeros(len(ratios))

    essential_genes = []
    for i, node in enumerate(G.nodes()):
        if z_scores[i] < z_threshold:
            essential_genes.append(node)
            G.nodes[node]["essential"] = True
        else:
            G.nodes[node]["essential"] = False

    return essential_genes


def pathway_enrichment(essential_genes, pathways=ESSENTIAL_PATHWAYS):
    """
    Simple overlap-based pathway enrichment for essential genes.
    """
    # Map perturbation names back to gene names
    ess_gene_names = set()
    for eg in essential_genes:
        gene_name = eg.replace("_guide", "").replace("gene_", "Gene_")
        ess_gene_names.add(gene_name)

    enrichment = {}
    for pathway, genes in pathways.items():
        overlap = ess_gene_names.intersection(genes)
        enrichment[pathway] = {
            "overlap": list(overlap),
            "n_overlap": len(overlap),
            "pathway_size": len(genes),
            "fraction": len(overlap) / len(genes) if genes else 0,
        }
    return enrichment


def plot_essential_network(G, community_map, fitness_scores, enrichment, out_dir="figures"):
    """Visualize essential gene network."""
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(18, 16))

    # 1. Network with communities
    ax = axes[0, 0]
    if len(G.nodes()) > 0:
        pos = nx.spring_layout(G, seed=SEED, k=1.5)

        # Node colors by community
        n_communities = len(set(community_map.values()))
        cmap = plt.cm.Set3(np.linspace(0, 1, max(n_communities, 1)))
        node_colors = [cmap[community_map.get(n, 0)] for n in G.nodes()]

        # Node sizes by essentiality score
        node_sizes = [max(G.nodes[n].get("essentiality_score", 1) * 100, 50) for n in G.nodes()]

        # Essential gene markers
        essential_nodes = [n for n in G.nodes() if G.nodes[n].get("essential", False)]
        non_essential = [n for n in G.nodes() if not G.nodes[n].get("essential", False)]

        nx.draw_networkx_nodes(G, pos, nodelist=non_essential, node_color=[
            cmap[community_map.get(n, 0)] for n in non_essential],
            node_size=[max(G.nodes[n].get("essentiality_score", 1) * 100, 50) for n in non_essential],
            alpha=0.6, ax=ax)

        if essential_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=essential_nodes, node_color="red",
                                   node_size=[max(G.nodes[n].get("essentiality_score", 1) * 200, 100)
                                              for n in essential_nodes],
                                   alpha=0.9, edgecolors="black", linewidths=2, ax=ax)

        nx.draw_networkx_edges(G, pos, alpha=0.2, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=5, ax=ax)

    ax.set_title(f"Essential Gene Network\n{len(G.nodes())} genes, {len(G.edges())} interactions")

    # 2. Fitness score distribution
    ax = axes[0, 1]
    ratios = [fitness_scores[n]["fitness_ratio"] for n in fitness_scores]
    ax.hist(ratios, bins=20, color="steelblue", edgecolor="white", alpha=0.7)
    ax.axvline(np.mean(ratios) - 1.5 * np.std(ratios), color="red", linestyle="--",
               label="Essential threshold")
    ax.set_xlabel("Fitness Ratio (perturbed / control)")
    ax.set_ylabel("N Perturbations")
    ax.set_title("Fitness Score Distribution")
    ax.legend()

    # 3. Pathway enrichment
    ax = axes[1, 0]
    if enrichment:
        pathways_sorted = sorted(enrichment.items(), key=lambda x: x[1]["fraction"], reverse=True)
        names = [p[0] for p in pathways_sorted]
        fracs = [p[1]["fraction"] for p in pathways_sorted]
        colors_bar = ["crimson" if f > 0 else "grey" for f in fracs]
        ax.barh(names, fracs, color=colors_bar, edgecolor="white")
        ax.set_xlabel("Overlap Fraction")
        ax.set_title("Pathway Enrichment of Essential Genes")

    # 4. Effect magnitude vs fitness
    ax = axes[1, 1]
    for pert, fs in fitness_scores.items():
        color = "red" if any(pert == eg for eg in
                              [n for n in G.nodes() if G.nodes[n].get("essential", False)]) else "steelblue"
        ax.scatter(fs["fitness_ratio"], fs["effect_magnitude"], c=color, s=50, alpha=0.6)
    ax.set_xlabel("Fitness Ratio")
    ax.set_ylabel("Transcriptional Effect Magnitude")
    ax.set_title("Fitness vs Transcriptional Effect")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "06_essential_network.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, "06_essential_network.svg"), bbox_inches="tight")
    plt.close()
    print(f"✓ Essential gene network plots saved")


def run_essential_pipeline(adata_path="data/perturbseq_with_latent.h5ad"):
    """Main essential gene network pipeline."""
    adata = ad.read_h5ad(adata_path)
    print(f"Loaded: {adata.shape}")

    # Compute fitness effects
    print("Computing fitness effects...")
    fitness_scores, pert_effects = compute_fitness_effects(adata)
    print(f"  Perturbations analyzed: {len(fitness_scores)}")

    # Build network
    print("Building essentiality network...")
    G, sim = build_essentiality_network(pert_effects, fitness_scores)
    print(f"  Network: {len(G.nodes())} nodes, {len(G.edges())} edges")

    # Community detection
    print("Detecting communities...")
    community_map = community_detection(G)
    n_communities = len(set(community_map.values()))
    print(f"  Communities: {n_communities}")

    # Essential gene classification
    print("Classifying essential genes...")
    essential_genes = annotate_essential_genes(G, fitness_scores)
    print(f"  Essential genes: {len(essential_genes)}")

    # Pathway enrichment
    enrichment = pathway_enrichment(essential_genes)

    # Plot
    plot_essential_network(G, community_map, fitness_scores, enrichment)

    # Save results
    os.makedirs("results", exist_ok=True)

    # Save fitness scores
    pd.DataFrame(fitness_scores).T.to_csv("results/06_fitness_scores.csv")

    # Save network
    edge_list = [(u, v, d["weight"]) for u, v, d in G.edges(data=True)]
    pd.DataFrame(edge_list, columns=["source", "target", "weight"]).to_csv(
        "results/06_network_edges.csv", index=False)

    summary = {
        "n_perturbations": len(fitness_scores),
        "n_network_nodes": len(G.nodes()),
        "n_network_edges": len(G.edges()),
        "n_communities": n_communities,
        "n_essential_genes": len(essential_genes),
        "essential_genes": essential_genes,
        "community_sizes": {str(c): sum(1 for v in community_map.values() if v == c)
                            for c in set(community_map.values())},
        "pathway_enrichment": enrichment,
        "mean_fitness_ratio": float(np.mean([f["fitness_ratio"] for f in fitness_scores.values()])),
        "network_density": float(nx.density(G)) if len(G.nodes()) > 1 else 0,
        "network_clustering_coeff": float(nx.average_clustering(G)) if len(G.nodes()) > 1 else 0,
    }
    with open("results/06_essential_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✓ Essential gene network complete:")
    print(f"  Essential genes: {len(essential_genes)}")
    print(f"  Communities: {n_communities}")
    print(f"  Network density: {summary['network_density']:.3f}")
    print(f"  Clustering coeff: {summary['network_clustering_coeff']:.3f}")

    return G, summary


if __name__ == "__main__":
    run_essential_pipeline()
